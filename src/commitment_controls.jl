# Fixed operating inputs are normalized once into unit and plant state bounds.
# Every scheduler and independent commitment check uses the same interpretation.

function _plant_initial_commitment(c, plant, members)
    state=any(c.system.generators[j].initial_on==1 for j in members) ? 1 : 0
    plant.initial_on===nothing || plant.initial_on==state ||
        throw(ArgumentError("plant $(plant.name) initial state disagrees with its units"))
    known_age=state==1 ? maximum((c.system.generators[j].initial_age
        for j in members if c.system.generators[j].initial_on==1);init=0.0) :
        (isempty(members) ? 0.0 : minimum(c.system.generators[j].initial_age for j in members))
    if plant.initial_age===nothing
        plant.minup==0 && plant.mindown==0 || throw(ArgumentError(
            "plant $(plant.name) requires initial_age when minimum dwell is active"))
        # An on plant can have stayed on through earlier unit handovers: the
        # longest currently-on unit age is a known lower bound, not full history.
        age=known_age
    else
        age=plant.initial_age
        consistent=state==1 ? age+1e-8>=known_age : abs(age-known_age)<=1e-8
        consistent || throw(ArgumentError(
            "plant $(plant.name) initial age disagrees with its unit histories"))
    end
    (;state,age)
end

function _commitment_operating_bounds(c, object, t)
    lower,upper=0,1
    forced=opinterval(c,object,:forced_on,t,-1.0)
    if forced>=0
        lower=max(lower,Int(forced))
        upper=min(upper,Int(forced))
    end
    if opinterval(c,object,:maintenance,t,0.0)==1
        upper=0
    end
    for attribute in (:power,:discharge)
        schedule=optional_operation(c,object,attribute,t)
        schedule===nothing && continue
        if schedule==0
            upper=0
        else
            lower=1
        end
    end
    (lower,upper)
end

"""Compile hard commitment restrictions without precedence-based overwrites.

Maintenance, forced states, zero/positive schedules and residual initial dwell
are intersected. Incompatible inputs raise an error, rather than weakening one
restriction. Plant state means at least one member unit is on. Off restrictions
propagate to every member; on restrictions retain the choice of member units.
"""
function _commitment_bounds(c)
    s=c.system
    G,P,T=length(s.generators),length(s.plants),length(c.prices)
    unit_lower=zeros(Int,G,T)
    unit_upper=ones(Int,G,T)
    plant_lower=zeros(Int,P,T)
    plant_upper=ones(Int,P,T)
    members=[findall(g->g.plant==p.name,s.generators) for p in s.plants]
    plant_initial_on=zeros(Int,P)
    plant_initial_age=zeros(P)
    for (j,g) in enumerate(s.generators),t in 1:T
        lo,hi=_commitment_operating_bounds(c,g.name,t)
        residual=(g.initial_on==1 ? g.minup : g.mindown)-g.initial_age
        if c.grid[t]-first(c.grid)<residual-1e-8
            lo=max(lo,g.initial_on)
            hi=min(hi,g.initial_on)
        end
        unit_lower[j,t],unit_upper[j,t]=lo,hi
    end
    for (i,p) in enumerate(s.plants)
        initial=_plant_initial_commitment(c,p,members[i])
        plant_initial_on[i],plant_initial_age[i]=initial.state,initial.age
        for t in 1:T
            lo,hi=_commitment_operating_bounds(c,p.name,t)
            residual=(initial.state==1 ? p.minup : p.mindown)-initial.age
            if c.grid[t]-first(c.grid)<residual-1e-8
                lo=max(lo,initial.state)
                hi=min(hi,initial.state)
            end
            if hi==0
                for j in members[i]
                    unit_upper[j,t]=0
                end
            end
            any(unit_lower[j,t]==1 for j in members[i]) && (lo=1)
            any(unit_upper[j,t]==1 for j in members[i]) || (hi=0)
            lo<=hi || throw(ArgumentError(
                "conflicting commitment controls for plant $(p.name) at $(c.grid[t])"))
            plant_lower[i,t],plant_upper[i,t]=lo,hi
        end
    end
    for (j,g) in enumerate(s.generators),t in 1:T
        unit_lower[j,t]<=unit_upper[j,t] || throw(ArgumentError(
            "conflicting commitment controls for unit $(g.name) at $(c.grid[t])"))
    end
    (;unit_lower,unit_upper,plant_lower,plant_upper,members,
        plant_initial_on,plant_initial_age)
end

function _commitment_dwell_admissible(grid, states, initial_on, initial_age, minup, mindown)
    state=initial_on
    since=first(grid)-initial_age
    for t in eachindex(states)
        if states[t]!=state
            grid[t]-since+1e-8 >= (state==1 ? minup : mindown) || return false
            state=states[t]
            since=grid[t]
        end
    end
    true
end

function _commitment_admissible(c,u;bounds=nothing)
    size(u)==(length(c.system.generators),length(c.prices)) || return false
    all(x->x==0 || x==1,u) || return false
    if bounds===nothing
        bounds=try
            _commitment_bounds(c)
        catch error
            error isa ArgumentError || rethrow()
            return false
        end
    end
    all(bounds.unit_lower .<= u .<= bounds.unit_upper) || return false
    for (j,g) in enumerate(c.system.generators)
        _commitment_dwell_admissible(c.grid,view(u,j,:),g.initial_on,g.initial_age,
            g.minup,g.mindown) || return false
    end
    for (i,p) in enumerate(c.system.plants)
        states=[any(u[j,t]==1 for j in bounds.members[i]) ? 1 : 0
            for t in eachindex(c.prices)]
        all(bounds.plant_lower[i,:] .<= states .<= bounds.plant_upper[i,:]) || return false
        _commitment_dwell_admissible(c.grid,states,bounds.plant_initial_on[i],
            bounds.plant_initial_age[i],p.minup,p.mindown) || return false
    end
    true
end

"""Independent aggregate plant state, elapsed age and remaining dwell at an edge."""
function _plant_commitment_history(c,u,time=last(c.grid))
    edge=findfirst(==(time),c.grid)
    edge===nothing && throw(ArgumentError("plant history time must be a grid edge"))
    history=NamedTuple[]
    for p in c.system.plants
        members=findall(g->g.plant==p.name,c.system.generators)
        initial=_plant_initial_commitment(c,p,members)
        state=initial.state
        since=first(c.grid)-initial.age
        for t in 1:(edge-1)
            current=any(u[j,t]==1 for j in members) ? 1 : 0
            if current!=state
                state=current
                since=c.grid[t]
            end
        end
        age=Float64(time-since)
        residual=max(0.0,(state==1 ? p.minup : p.mindown)-age)
        push!(history,(plant=p.name,state=state,age=age,residual=residual))
    end
    history
end
