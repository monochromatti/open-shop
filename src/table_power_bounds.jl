const _TablePowerTerm = Union{Float64,VariableRef,AffExpr}

struct _TablePowerCoordinateRecord
    onweights::Vector{VariableRef}
    u::VariableRef
    weights::Vector{VariableRef}
end

struct _TablePowerSupportCoordinate
    unit::Int
    interval::Int
    axis::Symbol
    table::TurbineTable
    qbox::NTuple{2,Float64}
    hbox::NTuple{2,Float64}
    electrical_max::Float64
    nodes::Vector{Float64}
    onweights::Vector{_TablePowerTerm}
    term::_TablePowerTerm
    power::_TablePowerTerm
    commitment::_TablePowerTerm
    discharge::_TablePowerTerm
    onhead::_TablePowerTerm
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

function _table_power_onweights!(m,coordinate,u,moment,name,records;affine=false)
    known=u isa Number ? u : is_fixed(u) ? fix_value(u) : nothing
    if known!==nothing
        return known==0 ? _TablePowerTerm[0.0 for _ in coordinate.weights] :
            _TablePowerTerm[known==1 ? w : known*w for w in coordinate.weights]
    end
    length(coordinate.nodes)==1 && return _TablePowerTerm[u]
    if affine && first(coordinate.nodes)==0.0 && all(>(0.0),coordinate.nodes[2:end])
        # Zero discharge has the unique coordinate vector (1,0,...), so its
        # on-state weights are affine and need no additional variables.
        @constraint(m,sum(coordinate.weights[2:end])<=u)
        return _TablePowerTerm[coordinate.weights[1]+u-1;coordinate.weights[2:end]]
    end
    n=length(coordinate.nodes)
    z=@variable(m,[1:n],lower_bound=0,upper_bound=1,base_name="$(name)_onweight")
    @constraint(m,[j=1:n],z[j]<=coordinate.weights[j])
    @constraint(m,sum(z)==u)
    @constraint(m,sum((coordinate.nodes[j]-coordinate.nodes[1])*z[j] for j in 2:n)==
        moment-coordinate.nodes[1]*u)
    push!(records,_TablePowerCoordinateRecord(z,u,coordinate.weights))
    _TablePowerTerm[z...]
end

function _lift_table_power_bounds!(m,assigned)
    records=get(m.ext,:global_table_power_coordinates,nothing)
    records===nothing && return assigned
    for record in records::Vector{_TablePowerCoordinateRecord}
        state=_table_power_value(record.u,assigned)
        for j in eachindex(record.onweights)
            assigned[record.onweights[j]]=state*_table_power_value(record.weights[j],assigned)
        end
    end
    assigned
end

"""Add power supports using shared on-state turbine-table coordinates.

Nodal intercepts enclose exact turbine power on the conditional on-state box.
A guarded polynomial certificate raises adjacent intercepts to cover the
interpolation residual. Shared weights sum to commitment and retain the
existing on-state moment, preserving signed intercepts and off-state heads.
The original table graphs and nonlinear power equations remain intact.
"""
function _add_table_power_bounds!(b,c)
    began=time_ns()
    m=b.m
    haskey(m.ext,:global_table_power_bounds_profile) &&
        throw(ArgumentError("table power bounds have already been added"))
    tensors=get(m.ext,:global_tensor_turbines,nothing)
    hulls=get(m.ext,:global_power_hulls,nothing)
    records=_TablePowerCoordinateRecord[]
    supports=_TablePowerSupportCoordinate[]
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
            onweights=_table_power_onweights!(m,coordinate,u,
                axis==:head ? hull.onhead : b.GQ[i,t],"table_power_$(axis)_$(i)_$(t)",records;
                affine=axis==:discharge)
            term=axis==:head ? b.GQ[i,t] : hull.onhead
            push!(supports,_TablePowerSupportCoordinate(i,t,axis,g.turbine_table,qbox,hbox,em,
                coordinate.nodes,onweights,term,b.P[i,t],u isa Number ? Float64(u) : u,
                b.GQ[i,t],hull.onhead))
            for k in eachindex(slopes)
                intercept=sum(coefficients[k][j]*onweights[j] for j in eachindex(onweights))
                @constraint(m,(b.P[i,t]-slopes[k]*term-intercept)/40<=0)
                axis==:head ? (head_rows+=1) : (discharge_rows+=1)
            end
            added=true
        end
        added && (units+=1)
    end
    m.ext[:global_table_power_coordinates]=records
    m.ext[:global_table_power_supports]=supports
    profile=Dict("axes"=>"both","on_state_coordinates"=>"shared","units_added"=>units,
        "head_rows"=>head_rows,"discharge_rows"=>discharge_rows,
        "variables_added"=>num_variables(m)-before_variables,
        "constraints_added"=>num_constraints(m;count_variable_in_set_constraints=false)-before_constraints,
        "shared_coordinates"=>length(records),"certificate_evaluations"=>evaluations[],
        "certificate_cache_hits"=>hits[],"certificate_seconds"=>certificate_seconds[],
        "setup_seconds"=>(time_ns()-began)/1e9,"skipped"=>skipped)
    m.ext[:global_table_power_bounds_profile]=profile
    profile
end
