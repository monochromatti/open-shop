function _commitment_transitions!(m, grid, u, su, sd, initial_on, minup, mindown)
    for t in eachindex(u)
        prev=t==1 ? initial_on : u[t-1]
        @constraint(m,u[t]-prev==su[t]-sd[t])
        @constraint(m,su[t]<=u[t])
        @constraint(m,su[t]<=1-prev)
        @constraint(m,sd[t]<=prev)
        @constraint(m,sd[t]<=1-u[t])
        # Binary state and these inequalities determine exact transitions.
        for k in t:length(u)
            grid[k]-grid[t]<minup-1e-8 && @constraint(m,u[k]>=su[t])
            grid[k]-grid[t]<mindown-1e-8 && @constraint(m,u[k]<=1-sd[t])
        end
    end
end

function _plant_commitment_states!(m,c,u,start,bounds;relaxed=false)
    P,T=length(c.system.plants),length(c.prices)
    plant_u=Matrix{Union{Nothing,VariableRef}}(nothing,P,T)
    plant_su=similar(plant_u)
    plant_sd=similar(plant_u)
    fill!(plant_su,nothing)
    fill!(plant_sd,nothing)
    for (i,p) in enumerate(c.system.plants)
        members=bounds.members[i]
        if p.minup==0 && p.mindown==0
            # On requires a member; off has already propagated to unit bounds.
            # Legacy plants with no aggregate controls add no state variables.
            for t in 1:T
                bounds.plant_lower[i,t]==1 && @constraint(m,
                    sum(u[j,t] for j in members;init=0.0)>=1)
            end
            continue
        end
        states=[@variable(m,lower_bound=0.0,upper_bound=1.0,
            base_name="plant_on[$i,$t]") for t in 1:T]
        starts=[@variable(m,lower_bound=0.0,upper_bound=1.0,
            base_name="plant_start[$i,$t]") for t in 1:T]
        stops=[@variable(m,lower_bound=0.0,upper_bound=1.0,
            base_name="plant_stop[$i,$t]") for t in 1:T]
        for t in 1:T
            relaxed || set_binary(states[t])
            for j in members
                @constraint(m,states[t]>=u[j,t])
            end
            @constraint(m,states[t]<=sum(u[j,t] for j in members;init=0.0))
            lo,hi=bounds.plant_lower[i,t],bounds.plant_upper[i,t]
            lo==hi && fix(states[t],lo;force=true)
            state=any(start[j,t]==1 for j in members) ? 1 : 0
            set_start_value(states[t],state)
            prior=t==1 ? bounds.plant_initial_on[i] :
                (any(start[j,t-1]==1 for j in members) ? 1 : 0)
            set_start_value(starts[t],max(0,state-prior))
            set_start_value(stops[t],max(0,prior-state))
        end
        _commitment_transitions!(m,c.grid,states,starts,stops,
            bounds.plant_initial_on[i],p.minup,p.mindown)
        plant_u[i,:]=states
        plant_su[i,:]=starts
        plant_sd[i,:]=stops
    end
    (;plant_u,plant_su,plant_sd)
end

function _joint_states!(m,c;relaxed=false,fixed_u=nothing,free_mask=nothing,incumbent=nothing)
    G,T=length(c.system.generators),length(c.prices)
    bounds=_commitment_bounds(c)
    start=incumbent===nothing ? ones(Int,G,T) : copy(incumbent["u"])
    size(start)==(G,T) || throw(ArgumentError("incorrect incumbent shape"))
    if fixed_u!==nothing
        _commitment_admissible(c,fixed_u;bounds) || throw(ArgumentError("invalid fixed commitment"))
        start=copy(fixed_u)
    end
    if free_mask!==nothing
        incumbent===nothing && throw(ArgumentError("a restricted search needs an incumbent"))
        size(free_mask)==(G,T) || throw(ArgumentError("incorrect free-mask shape"))
        all(x->x==true || x==false,free_mask) || throw(ArgumentError("free-mask must be boolean"))
        _commitment_admissible(c,start;bounds) || throw(ArgumentError("invalid incumbent commitment"))
    end
    @variable(m,0<=u[1:G,1:T]<=1)
    @variable(m,0<=su[1:G,1:T]<=1)
    @variable(m,0<=sd[1:G,1:T]<=1)
    for (j,g) in enumerate(c.system.generators)
        for t in 1:T
            relaxed || set_binary(u[j,t])
            lo,hi=bounds.unit_lower[j,t],bounds.unit_upper[j,t]
            if fixed_u!==nothing || (free_mask!==nothing && !free_mask[j,t])
                lo=hi=Int(start[j,t])
            end
            lo==hi && fix(u[j,t],lo;force=true)
            start[j,t]=clamp(start[j,t],lo,hi)
            set_start_value(u[j,t],start[j,t])
            prior=t==1 ? g.initial_on : start[j,t-1]
            set_start_value(su[j,t],max(0,start[j,t]-prior))
            set_start_value(sd[j,t],max(0,prior-start[j,t]))
        end
        _commitment_transitions!(m,c.grid,view(u,j,:),view(su,j,:),view(sd,j,:),
            g.initial_on,g.minup,g.mindown)
    end
    plant_states=_plant_commitment_states!(m,c,u,start,bounds;relaxed)
    (;u,su,sd,start,plant_states...,bounds)
end
