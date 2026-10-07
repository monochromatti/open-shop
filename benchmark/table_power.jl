module TablePower
using JuMP
import OpenSHOP
export add_table_power!

const _FRACTIONS=(0.25,0.5,0.75,1.0)
_value(x,assigned)=x isa Number ? Float64(x) : JuMP.value(v->assigned[v],x)
_fixed_off(u)=u isa Number ? u==0 : is_fixed(u)&&fix_value(u)==0

function _support_cache()
    cache=IdDict{OpenSHOP.TurbineTable,Dict{NTuple{7,Float64},Float64}}()
    evaluations=Ref(0);hits=Ref(0);seconds=Ref(0.0)
    function support(table,qlo,qhi,hlo,hhi,em,a,b=0.0)
        entries=get!(cache,table) do
            Dict{NTuple{7,Float64},Float64}()
        end
        key=Float64.((qlo,qhi,hlo,hhi,em,a,b))
        if haskey(entries,key)
            hits[]+=1
            return entries[key]
        end
        evaluations[]+=1
        began=time_ns()
        result=OpenSHOP._power_support(table,qlo,qhi,hlo,hhi,em,a,b)
        seconds[]+=(time_ns()-began)/1e9
        entries[key]=result
    end
    support,evaluations,hits,seconds
end

# A nodal envelope alone misses concave interpolation between knots. Bound the
# residual against each original coordinate chord on its physical on portion,
# then raise both endpoints by at least that cell's certified excess.
function _head_coefficients(table,qbox,hbox,nodes,em,a,support)
    base=[support(table,qbox...,h,h,em,a) for h in nodes]
    correction=zeros(length(nodes))
    for j in 1:(length(nodes)-1)
        lo=max(nodes[j],hbox[1]);hi=min(nodes[j+1],hbox[2])
        lo>hi && continue
        beta=(base[j+1]-base[j])/(nodes[j+1]-nodes[j])
        gamma=base[j]-beta*nodes[j]
        upper=support(table,qbox...,lo,hi,em,a,beta)
        delta=max(0.0,upper-gamma)+1e-9*max(1.0,abs(upper),abs(gamma))
        correction[j]=max(correction[j],delta)
        correction[j+1]=max(correction[j+1],delta)
    end
    base.+correction
end

function _discharge_coefficients(table,qbox,hbox,nodes,em,b,support)
    base=[support(table,q,q,hbox...,em,0.0,b) for q in nodes]
    correction=zeros(length(nodes))
    for j in 1:(length(nodes)-1)
        lo=max(nodes[j],qbox[1]);hi=min(nodes[j+1],qbox[2])
        lo>hi && continue
        alpha=(base[j+1]-base[j])/(nodes[j+1]-nodes[j])
        gamma=base[j]-alpha*nodes[j]
        upper=support(table,lo,hi,hbox...,em,alpha,b)
        delta=max(0.0,upper-gamma)+1e-9*max(1.0,abs(upper),abs(gamma))
        correction[j]=max(correction[j],delta)
        correction[j+1]=max(correction[j+1],delta)
    end
    base.+correction
end

function _intercept!(m,coefficients,weights,u,gate,name,records)
    c=gate ? coefficients : max.(0.0,coefficients)
    # The coordinate weights sum to one, including singleton coordinates.
    if all(==(first(c)),c)
        return gate ? first(c)*u : first(c)
    end
    expression=sum(c[j]*weights[j] for j in eachindex(c))
    !gate && return expression
    known=u isa Number ? u : is_fixed(u) ? fix_value(u) : nothing
    known!==nothing && return known*expression
    lo,hi=extrema(c)
    z=@variable(m,lower_bound=min(0.0,lo),upper_bound=max(0.0,hi),base_name="$(name)_on")
    @constraint(m,z>=lo*u)
    @constraint(m,z<=hi*u)
    @constraint(m,z>=expression-hi*(1-u))
    @constraint(m,z<=expression-lo*(1-u))
    push!(records,(z=z,u=u,expression=expression))
    z
end

function _lift_gates!(records,assigned)
    for record in records
        assigned[record.z]=_value(record.u,assigned)*_value(record.expression,assigned)
    end
    assigned
end

"""Add certified table-coordinate power supports without new integer variables.

`axes` selects head, discharge, or both existing SOS2 coordinates. Ungated
intercepts are nonnegative for off validity. Gated intercepts retain their signs
and use at most one binary-product auxiliary per supporting row.
"""
function add_table_power!(b,c;axes=:head,gate=false)
    axes in (:head,:discharge,:both) || throw(ArgumentError("unsupported table-power axis"))
    began=time_ns()
    m=b.m
    haskey(m.ext,:table_power_profile) && throw(ArgumentError("table power has already been added"))
    tensors=get(m.ext,:global_tensor_turbines,nothing)
    hulls=get(m.ext,:global_power_hulls,nothing)
    records=Any[]
    before_variables=num_variables(m)
    before_constraints=num_constraints(m;count_variable_in_set_constraints=false)
    support,evaluations,hits,certificate_seconds=_support_cache()
    skipped=Dict("analytic"=>0,"fixed_off"=>0,"electrical"=>0,"graph_or_hull"=>0,"nonfinite"=>0)
    units=0;head_rows=0;discharge_rows=0
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices)
        g.turbine_table===nothing && (skipped["analytic"]+=1;continue)
        u=b.u[i,t]
        _fixed_off(u) && (skipped["fixed_off"]+=1;continue)
        plo=max(0.0,OpenSHOP.opinterval(c,g.name,:pmin,t,g.pmin))
        phi=min(g.pmax,OpenSHOP.opinterval(c,g.name,:pmax,t,g.pmax))
        electrical=if g.generator_efficiency_curve===nothing
            [1.0]
        elseif isfinite(plo)&&isfinite(phi)&&plo<=phi
            curve=g.generator_efficiency_curve
            [OpenSHOP.table_value(curve,p;extrapolation=:linear) for p in
                OpenSHOP._global_tensor_nodes(curve.x,plo,phi)]
        else
            Float64[]
        end
        if isempty(electrical)||!all(isfinite,electrical)||minimum(electrical)<0
            skipped["electrical"]+=1;continue
        end
        tensor=tensors===nothing ? nothing : get(tensors,"turbine_$(i)_$(t)",nothing)
        hull=hulls===nothing ? nothing : get(hulls,"power_hull_$(i)_$(t)",nothing)
        if tensor===nothing || hull===nothing
            skipped["graph_or_hull"]+=1;continue
        end
        qbox=(first(hull.axes[1]),last(hull.axes[1]))
        hbox=(first(hull.axes[2]),last(hull.axes[2]))
        em=maximum(electrical)
        added=false
        for axis in (axes==:both ? (:head,:discharge) : (axes,))
            coordinate=axis==:head ? tensor.hcoordinate : tensor.qcoordinate
            slopes=[0.00981*em*(axis==:head ? hbox[2] : qbox[2])*f for f in _FRACTIONS]
            coefficients=[axis==:head ?
                _head_coefficients(g.turbine_table,qbox,hbox,coordinate.nodes,em,a,support) :
                _discharge_coefficients(g.turbine_table,qbox,hbox,coordinate.nodes,em,a,support)
                for a in slopes]
            if !all(v->all(isfinite,v),coefficients) || !all(isfinite,slopes)
                skipped["nonfinite"]+=1;continue
            end
            for k in eachindex(slopes)
                n="table_power_$(axis)_$(i)_$(t)_$(k)"
                intercept=_intercept!(m,coefficients[k],coordinate.weights,u,gate,n,records)
                term=slopes[k]*(axis==:head ? b.GQ[i,t] : hull.onhead)
                @constraint(m,(b.P[i,t]-term-intercept)/40<=0)
                axis==:head ? (head_rows+=1) : (discharge_rows+=1)
            end
            added=true
        end
        added && (units+=1)
    end
    m.ext[:table_power_gates]=records
    !isempty(records) && push!(get!(m.ext,:experiment_start_lifters,Any[]),
        assigned->_lift_gates!(records,assigned))
    profile=Dict("axes"=>string(axes),"gate"=>gate,"units_added"=>units,
        "head_rows"=>head_rows,"discharge_rows"=>discharge_rows,
        "variables_added"=>num_variables(m)-before_variables,
        "constraints_added"=>num_constraints(m;count_variable_in_set_constraints=false)-before_constraints,
        "gates_added"=>length(records),"certificate_evaluations"=>evaluations[],
        "certificate_cache_hits"=>hits[],"certificate_seconds"=>certificate_seconds[],
        "setup_seconds"=>(time_ns()-began)/1e9,"skipped"=>skipped)
    m.ext[:table_power_profile]=profile
    profile
end
end
