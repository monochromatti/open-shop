# River costs act on upstream interval-average flow. Motion is total variation
# of that sequence; start/stop costs use the disjunction q=0 or q>=q_min.
const RIVER_COST_ATTRIBUTES = (:flow_change_cost,:flow_start_cost,:flow_stop_cost)
river_cost_enabled(c,r,attributes=RIVER_COST_ATTRIBUTES)=
    any(a->getfield(r,a)>0,attributes) ||
    any(z->z.object==r.name && z.attribute in attributes && any(>(0),z.values),c.operations)

function validate_river_operations(c)
    for r in c.system.rivers
        all(a->isfinite(getfield(r,a)) && getfield(r,a)>=0,RIVER_COST_ATTRIBUTES) ||
            throw(ArgumentError("$(r.name): river costs must be finite and nonnegative"))
        has_transitions=river_cost_enabled(c,r,(:flow_start_cost,:flow_stop_cost))
        if r.min_allowed_flow===nothing
            has_transitions && throw(ArgumentError("$(r.name): river start/stop costs require min_allowed_flow"))
            continue
        end
        !is_channel(r) && isfinite(r.min_allowed_flow) && 0<r.min_allowed_flow<=r.capacity ||
            throw(ArgumentError("$(r.name): flow commitment needs a positive minimum within routed river capacity"))
        initial,_=river_initial_output(r)
        initial!==nothing && (initial==0 || initial>=r.min_allowed_flow) ||
            throw(ArgumentError("$(r.name): flow commitment needs zero or on-state initial_release/history"))
    end
    nothing
end

# Midpoint classification is unambiguous on the two allowed domains. The
# independent gap residual below rejects values between them.
river_flow_state(r,q)=r.min_allowed_flow===nothing ? 0 : Int(q>=r.min_allowed_flow/2)
function river_on_schedule(c,release)
    [river_flow_state(r,release[i,t]) for (i,r) in enumerate(c.system.rivers),t in eachindex(c.prices)]
end

function fixed_river_on(c,on,warm)
    if on===nothing
        return warm===nothing ? [Int(r.min_allowed_flow!==nothing) for r in c.system.rivers,t in eachindex(c.prices)] :
            river_on_schedule(c,warm.river_release)
    end
    size(on)==(length(c.system.rivers),length(c.prices)) && all(x->x in (0,1),on) ||
        throw(ArgumentError("river_on must be a binary [river,interval] matrix"))
    Matrix{Int}(on)
end

function river_operating_cost(c,release)
    change,start,stop=0.0,0.0,0.0
    for (i,r) in enumerate(c.system.rivers)
        river_cost_enabled(c,r) || continue
        initial,_=river_initial_output(r)
        for t in eachindex(c.prices)
            previous=t==1 ? initial : release[i,t-1]
            previous===nothing && continue
            q=release[i,t]
            change+=opvalue(c,r.name,:flow_change_cost,c.grid[t],r.flow_change_cost)*abs(q-previous)
            if r.min_allowed_flow!==nothing
                old,new=river_flow_state(r,previous),river_flow_state(r,q)
                start+=opvalue(c,r.name,:flow_start_cost,c.grid[t],r.flow_start_cost)*max(0,new-old)
                stop+=opvalue(c,r.name,:flow_stop_cost,c.grid[t],r.flow_stop_cost)*max(0,old-new)
            end
        end
    end
    (;change,start,stop,total=change+start+stop)
end

function constrain_river_operations!(m,c,release;fixed=nothing,relaxed=false)
    D,T=size(release)
    on=Matrix{Union{Nothing,VariableRef,Float64}}(nothing,D,T)
    starts=similar(on);fill!(starts,nothing)
    stops=similar(on);fill!(stops,nothing)
    changes=Matrix{Union{Nothing,VariableRef}}(nothing,D,T)
    cost=AffExpr(0.0)
    for (i,r) in enumerate(c.system.rivers)
        initial,_=river_initial_output(r)
        if r.min_allowed_flow!==nothing
            for t in 1:T
                state=fixed===nothing ? @variable(m,lower_bound=0.0,upper_bound=1.0,
                    base_name="river_on[$i,$t]") : Float64(fixed[i,t])
                state isa VariableRef && !relaxed && set_binary(state)
                on[i,t]=state
                @constraint(m,release[i,t]>=r.min_allowed_flow*state)
                @constraint(m,release[i,t]<=opinterval(c,r.name,:capacity,t,r.capacity)*state)
                prior=t==1 ? river_flow_state(r,initial) : on[i,t-1]
                if fixed===nothing
                    su=@variable(m,lower_bound=0.0,upper_bound=1.0,base_name="river_start[$i,$t]")
                    sd=@variable(m,lower_bound=0.0,upper_bound=1.0,base_name="river_stop[$i,$t]")
                    @constraint(m,state-prior==su-sd)
                    @constraint(m,su<=state)
                    @constraint(m,su<=1-prior)
                    @constraint(m,sd<=prior)
                    @constraint(m,sd<=1-state)
                else
                    su,sd=max(0,state-prior),max(0,prior-state)
                end
                starts[i,t],stops[i,t]=su,sd
                cost+=opvalue(c,r.name,:flow_start_cost,c.grid[t],r.flow_start_cost)*su+
                      opvalue(c,r.name,:flow_stop_cost,c.grid[t],r.flow_stop_cost)*sd
            end
        end
        river_cost_enabled(c,r,(:flow_change_cost,)) || continue
        for t in 1:T
            previous=t==1 ? initial : release[i,t-1]
            previous===nothing && continue
            coefficient=opvalue(c,r.name,:flow_change_cost,c.grid[t],r.flow_change_cost)
            coefficient==0 && continue
            delta=@variable(m,lower_bound=0.0,upper_bound=r.capacity+max(r.capacity,initial===nothing ? 0.0 : abs(initial)),base_name="river_change[$i,$t]")
            @constraint(m,delta>=release[i,t]-previous)
            @constraint(m,delta>=previous-release[i,t])
            changes[i,t]=delta
            cost+=coefficient*delta
        end
    end
    (;on,starts,stops,changes,cost)
end

function river_operation_residuals(c,release)
    result=Dict{String,Float64}()
    for (i,r) in enumerate(c.system.rivers)
        r.min_allowed_flow===nothing && continue
        result["river_flow_state_$(r.name)"]=maximum(q->min(abs(q),max(0.0,r.min_allowed_flow-q)),view(release,i,:);init=0.0)
    end
    result
end

function river_operation_start_values(operations,c,release)
    pairs=Pair{VariableRef,Float64}[]
    for (i,r) in enumerate(c.system.rivers),t in eachindex(c.prices)
        initial,_=river_initial_output(r)
        previous=t==1 ? initial : release[i,t-1]
        state=river_flow_state(r,release[i,t])
        prior=previous===nothing ? 0 : river_flow_state(r,previous)
        for (variable,value) in ((operations.on[i,t],state),
            (operations.starts[i,t],max(0,state-prior)),(operations.stops[i,t],max(0,prior-state)),
            (operations.changes[i,t],previous===nothing ? 0.0 : abs(release[i,t]-previous)))
            variable isa VariableRef && push!(pairs,variable=>Float64(value))
        end
    end
    pairs
end
