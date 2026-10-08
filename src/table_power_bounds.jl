struct _TablePowerGateRecord
    z::VariableRef
    u::VariableRef
    expression::AffExpr
end

struct _TablePowerCoordinateRecord
    onweights::Vector{VariableRef}
    u::VariableRef
    weights::Vector{VariableRef}
end

# Private experiment policies preserve the same certified coefficient vectors.
function _table_power_onweights!(m,coordinate,u,moment,name,records;affine=false)
    known=u isa Number ? u : is_fixed(u) ? fix_value(u) : nothing
    known!==nothing && return known.*coordinate.weights
    length(coordinate.nodes)==1 && return [u]
    if affine && first(coordinate.nodes)==0.0 && all(>(0.0),coordinate.nodes[2:end])
        # Off discharge is zero, hence its weight vector is uniquely (1,0,...).
        @constraint(m,sum(coordinate.weights[2:end])<=u)
        return vcat(coordinate.weights[1]+u-1,coordinate.weights[2:end])
    end
    n=length(coordinate.nodes)
    z=@variable(m,[1:n],lower_bound=0,upper_bound=1,base_name="$(name)_onweight")
    @constraint(m,[j=1:n],z[j]<=coordinate.weights[j])
    @constraint(m,sum(z)==u)
    # Centering gives the same moment with less cancellation at large heads.
    @constraint(m,sum((coordinate.nodes[j]-coordinate.nodes[1])*z[j] for j in 2:n)==moment-coordinate.nodes[1]*u)
    push!(records,_TablePowerCoordinateRecord(z,u,coordinate.weights))
    z
end

const _TABLE_POWER_FRACTIONS=(0.25,0.5,0.75,1.0)
_table_power_value(x,assigned)=x isa Number ? Float64(x) : JuMP.value(v->assigned[v],x)
_table_power_fixed_off(u)=u isa Number ? u==0 : is_fixed(u)&&fix_value(u)==0

function _table_power_support_cache()
    cache=IdDict{TurbineTable,Dict{NTuple{7,Float64},Float64}}()
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
        result=_power_support(table,qlo,qhi,hlo,hhi,em,a,b)
        seconds[]+=(time_ns()-began)/1e9
        entries[key]=result
    end
    support,evaluations,hits,seconds
end

# A nodal envelope alone misses concave interpolation between knots. Bound the
# residual against each original coordinate chord on its physical on portion,
# then raise both endpoints by at least that cell's certified excess.
function _table_power_head_coefficients(table,qbox,hbox,nodes,em,a,support)
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

function _table_power_discharge_coefficients(table,qbox,hbox,nodes,em,b,support)
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

function _table_power_intercept!(m,c,weights,u,name,records)
    # The coordinate weights sum to one, including singleton coordinates.
    if all(==(first(c)),c)
        return first(c)*u
    end
    expression=sum(c[j]*weights[j] for j in eachindex(c))
    known=u isa Number ? u : is_fixed(u) ? fix_value(u) : nothing
    known!==nothing && return known*expression
    lo,hi=extrema(c)
    z=@variable(m,lower_bound=min(0.0,lo),upper_bound=max(0.0,hi),base_name="$(name)_on")
    @constraint(m,z>=lo*u)
    @constraint(m,z<=hi*u)
    @constraint(m,z>=expression-hi*(1-u))
    @constraint(m,z<=expression-lo*(1-u))
    push!(records,_TablePowerGateRecord(z,u,expression))
    z
end

function _lift_table_power_bounds!(m,assigned)
    records=get(m.ext,:global_table_power_gates,nothing)
    for record in something(records,_TablePowerGateRecord[])::Vector{_TablePowerGateRecord}
        assigned[record.z]=_table_power_value(record.u,assigned)*_table_power_value(record.expression,assigned)
    end
    for record in get(m.ext,:global_table_power_coordinates,_TablePowerCoordinateRecord[])::Vector{_TablePowerCoordinateRecord}
        for j in eachindex(record.onweights)
            assigned[record.onweights[j]]=_table_power_value(record.u,assigned)*_table_power_value(record.weights[j],assigned)
        end
    end
    assigned
end

"""Add gated power supports on both existing turbine-table SOS2 coordinates.

Nodal intercepts enclose exact turbine power on the conditional on-state box.
A guarded polynomial certificate raises adjacent intercepts to cover the
interpolation residual. Binary-product rows gate signed intercepts off, while
retaining the original table graphs and nonlinear power equations.
"""
function _add_table_power_bounds!(b,c;policy=:baseline)
    policy in (:baseline,:head_shared,:discharge_affine,:shared) || throw(ArgumentError("unknown table power policy"))
    began=time_ns()
    m=b.m
    haskey(m.ext,:global_table_power_bounds_profile) &&
        throw(ArgumentError("table power bounds have already been added"))
    tensors=get(m.ext,:global_tensor_turbines,nothing)
    hulls=get(m.ext,:global_power_hulls,nothing)
    records=_TablePowerGateRecord[]
    coordinate_records=_TablePowerCoordinateRecord[]
    support_records=Any[]
    before_variables=num_variables(m)
    before_constraints=num_constraints(m;count_variable_in_set_constraints=false)
    support,evaluations,hits,certificate_seconds=_table_power_support_cache()
    skipped=Dict("analytic"=>0,"fixed_off"=>0,"electrical"=>0,"graph_or_hull"=>0,"nonfinite"=>0)
    units=0;head_rows=0;discharge_rows=0
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices)
        g.turbine_table===nothing && (skipped["analytic"]+=1;continue)
        u=b.u[i,t]
        _table_power_fixed_off(u) && (skipped["fixed_off"]+=1;continue)
        plo=max(0.0,opinterval(c,g.name,:pmin,t,g.pmin))
        phi=min(g.pmax,opinterval(c,g.name,:pmax,t,g.pmax))
        electrical=if g.generator_efficiency_curve===nothing
            [1.0]
        elseif isfinite(plo)&&isfinite(phi)&&plo<=phi
            curve=g.generator_efficiency_curve
            [table_value(curve,p;extrapolation=:linear) for p in
                _global_tensor_nodes(curve.x,plo,phi)]
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
        for axis in (:head,:discharge)
            coordinate=axis==:head ? tensor.hcoordinate : tensor.qcoordinate
            slopes=[0.00981*em*(axis==:head ? hbox[2] : qbox[2])*f for f in _TABLE_POWER_FRACTIONS]
            coefficients=[axis==:head ?
                _table_power_head_coefficients(g.turbine_table,qbox,hbox,coordinate.nodes,em,a,support) :
                _table_power_discharge_coefficients(g.turbine_table,qbox,hbox,coordinate.nodes,em,a,support)
                for a in slopes]
            if !all(v->all(isfinite,v),coefficients) || !all(isfinite,slopes)
                skipped["nonfinite"]+=1;continue
            end
            share=axis==:head ? policy in (:head_shared,:shared) : policy in (:discharge_affine,:shared)
            onweights=share ? _table_power_onweights!(m,coordinate,u,
                axis==:head ? hull.onhead : b.GQ[i,t],"table_power_$(axis)_$(i)_$(t)",coordinate_records;
                affine=axis==:discharge) : nothing
            push!(support_records,(unit=i,interval=t,axis=axis,slopes=slopes,coefficients=coefficients,
                onweights=onweights,qbox=qbox,hbox=hbox,em=em,table=g.turbine_table,nodes=coordinate.nodes,
                term=axis==:head ? b.GQ[i,t] : hull.onhead,power=b.P[i,t],u=u))
            for k in eachindex(slopes)
                n="table_power_$(axis)_$(i)_$(t)_$(k)"
                term=slopes[k]*(axis==:head ? b.GQ[i,t] : hull.onhead)
                intercept=share ? sum(coefficients[k][j]*onweights[j] for j in eachindex(onweights)) :
                    _table_power_intercept!(m,coefficients[k],coordinate.weights,u,n,records)
                @constraint(m,(b.P[i,t]-term-intercept)/40<=0)
                axis==:head ? (head_rows+=1) : (discharge_rows+=1)
            end
            added=true
        end
        added && (units+=1)
    end
    m.ext[:global_table_power_gates]=records
    m.ext[:global_table_power_coordinates]=coordinate_records
    m.ext[:global_table_power_supports]=support_records
    profile=Dict("policy"=>string(policy),"axes"=>"both","gate"=>true,"units_added"=>units,
        "head_rows"=>head_rows,"discharge_rows"=>discharge_rows,
        "variables_added"=>num_variables(m)-before_variables,
        "constraints_added"=>num_constraints(m;count_variable_in_set_constraints=false)-before_constraints,
        "gates_added"=>length(records),"shared_coordinates"=>length(coordinate_records),"certificate_evaluations"=>evaluations[],
        "certificate_cache_hits"=>hits[],"certificate_seconds"=>certificate_seconds[],
        "setup_seconds"=>(time_ns()-began)/1e9,"skipped"=>skipped)
    m.ext[:global_table_power_bounds_profile]=profile
    profile
end
