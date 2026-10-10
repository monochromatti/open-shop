# Ordinary operating constraints shared by nonlinear dispatch and its proposal.
# Optional inputs are resolved during construction; the solver receives algebraic rows.

const RAMP_ATTRIBUTES = (
    :ramp_up, :ramp_down, :discharge_ramp_up, :discharge_ramp_down,
    :volume_ramp_up, :volume_ramp_down, :level_ramp_up, :level_ramp_down,
)

function ramp_default(object, attribute)
    value=getfield(object, attribute)
    object isa Plant && attribute in (:ramp_up, :ramp_down) && value===nothing ?
        object.ramp : value
end

"""Integrate a rate limit over physical hours. An unrestricted segment disables
the comparison; series before their first knot retain the static default.
"""
function ramp_allowance(c, object, attribute, a, b)
    i=findfirst(z->z.object==object.name && z.attribute==attribute, c.operations)
    knots=i===nothing ? Float64[] : c.operations[i].times
    edges=vcat(a, [t for t in knots if a<t<b], b)
    total=0.0
    for k in 1:(length(edges)-1)
        rate=optional_opvalue(c, object.name, attribute, edges[k], ramp_default(object, attribute))
        rate===nothing && return nothing
        total+=rate*(edges[k+1]-edges[k])
    end
    isfinite(total) || throw(ArgumentError("nonfinite integrated ramp for $(object.name).$attribute"))
    total
end

function interval_midpoints(grid, t, historical_hours)
    b=grid[t]+(grid[t+1]-grid[t])/2
    a=t==1 ? first(grid)-historical_hours/2 : grid[t-1]+(grid[t]-grid[t-1])/2
    a, b
end

function unit_transitions(c, u; transitions=nothing)
    transitions!==nothing && return transitions.su, transitions.sd
    G, T=size(u)
    su=zeros(G, T)
    sd=zeros(G, T)
    for (j, g) in enumerate(c.system.generators), t in 1:T
        previous=t==1 ? g.initial_on : u[j, t-1]
        su[j, t]=max(0, u[j, t]-previous)
        sd[j, t]=max(0, previous-u[j, t])
    end
    su, sd
end

function transition_costs(c, u; transitions=nothing)
    su, sd=unit_transitions(c, u; transitions)
    sum(opinterval(c, g.name, :startup, t, g.startup)*su[j,t]+
        opinterval(c, g.name, :shutdown, t, g.shutdown)*sd[j,t]
        for (j,g) in enumerate(c.system.generators), t in eachindex(c.prices); init=0.0)
end

function power_transition_allowances(c, members, t, su, sd)
    up=sum(opinterval(c, c.system.generators[j].name, :pmin, t,
        c.system.generators[j].pmin)*su[j,t] for j in members; init=0.0)
    down=sum((t==1 ? opvalue(c,c.system.generators[j].name,:pmin,first(c.grid),
        c.system.generators[j].pmin;side=:left) :
        opinterval(c, c.system.generators[j].name, :pmin, t-1,
            c.system.generators[j].pmin))*sd[j,t] for j in members; init=0.0)
    up, down
end

function constrain_interval_ramps!(m, c, object, values, initial, historical_hours,
    attributes; su=nothing, sd=nothing, members=Int[], margin=0.0)
    for t in eachindex(c.prices)
        previous=t==1 ? initial : values[t-1]
        previous===nothing && continue
        a,b=interval_midpoints(c.grid,t,historical_hours)
        up=ramp_allowance(c,object,attributes[1],a,b)
        down=ramp_allowance(c,object,attributes[2],a,b)
        jump_up,jump_down=su===nothing ? (0.0,0.0) : power_transition_allowances(c,members,t,su,sd)
        change=values[t]-previous
        up===nothing || @constraint(m, change<=max(0.0,up-margin)+jump_up)
        down===nothing || @constraint(m, -change<=max(0.0,down-margin)+jump_down)
    end
    nothing
end

function constrain_schedules!(m,c,object,p,q; proposal=false)
    for t in eachindex(c.prices)
        scheduled_power=optional_operation(c,object.name,:power,t)
        scheduled_discharge=optional_operation(c,object.name,:discharge,t)
        # The proposal's P≈alpha*Q cannot enforce a physical power schedule exactly.
        # Its shared commitment rules retain positive/zero schedule implications.
        proposal || scheduled_power===nothing || @constraint(m,p[t]==scheduled_power)
        scheduled_discharge===nothing || @constraint(m,q[t]==scheduled_discharge)
    end
end

function plant_initial_output(c,plant,members,attribute)
    field=attribute==:power ? :initial_power : :initial_discharge
    explicit=getfield(plant,field)
    explicit!==nothing && return explicit,plant.initial_interval_hours
    attributes=attribute==:power ? (:ramp_up,:ramp_down) :
        (:discharge_ramp_up,:discharge_ramp_down)
    enabled=any(ramp_default(plant,a)!==nothing for a in attributes) ||
        any(z.object==plant.name && z.attribute in attributes for z in c.operations)
    enabled || return nothing,plant.initial_interval_hours
    isempty(members) && return nothing,plant.initial_interval_hours
    values=[getfield(c.system.generators[j],field) for j in members]
    any(isnothing,values) && return nothing,plant.initial_interval_hours
    durations=[c.system.generators[j].initial_interval_hours for j in members]
    all(==(first(durations)),durations) || throw(ArgumentError(
        "plant $(plant.name) needs aggregate historical $attribute when unit history windows differ"))
    sum(values),first(durations)
end

function river_initial_output(r)
    r.initial_release!==nothing && return r.initial_release,r.initial_interval_hours
    isempty(r.history_release) && return nothing,1.0
    last(r.history_release),r.history_grid[end]-r.history_grid[end-1]
end

function constrain_dispatch_operations!(m,c,u,P,GQ,RQ;
    transitions=nothing,margin=0.0,proposal=false)
    su,sd=unit_transitions(c,u;transitions)
    for (i,g) in enumerate(c.system.generators)
        constrain_schedules!(m,c,g,view(P,i,:),view(GQ,i,:);proposal)
        constrain_interval_ramps!(m,c,g,view(P,i,:),g.initial_power,g.initial_interval_hours,
            (:ramp_up,:ramp_down);su,sd,members=[i],margin)
        constrain_interval_ramps!(m,c,g,view(GQ,i,:),g.initial_discharge,g.initial_interval_hours,
            (:discharge_ramp_up,:discharge_ramp_down))
    end
    for p in c.system.plants
        members=findall(g->g.plant==p.name,c.system.generators)
        total_p=[sum(P[j,t] for j in members;init=0.0) for t in eachindex(c.prices)]
        total_q=[sum(GQ[j,t] for j in members;init=0.0) for t in eachindex(c.prices)]
        for t in eachindex(c.prices)
            @constraint(m,total_p[t]<=max(0.0,opinterval(c,p.name,:pmax,t,p.pmax)-margin))
            pmin=opinterval(c,p.name,:pmin,t,p.pmin)
            qmin=opinterval(c,p.name,:qmin,t,p.qmin)
            for j in members
                pmin>0 && @constraint(m,total_p[t]>=pmin*u[j,t])
                qmin>0 && @constraint(m,total_q[t]>=qmin*u[j,t])
            end
            qmax=optional_operation(c,p.name,:qmax,t;default=p.qmax)
            qmax===nothing || @constraint(m,total_q[t]<=qmax)
        end
        constrain_schedules!(m,c,p,total_p,total_q;proposal)
        initial_p,p_hours=plant_initial_output(c,p,members,:power)
        initial_q,q_hours=plant_initial_output(c,p,members,:discharge)
        constrain_interval_ramps!(m,c,p,total_p,initial_p,p_hours,
            (:ramp_up,:ramp_down);su,sd,members,margin)
        constrain_interval_ramps!(m,c,p,total_q,initial_q,q_hours,
            (:discharge_ramp_up,:discharge_ramp_down))
    end
    for (i,r) in enumerate(c.system.rivers)
        initial,duration=river_initial_output(r)
        constrain_interval_ramps!(m,c,r,view(RQ,i,:),initial,duration,(:ramp_up,:ramp_down))
    end
    nothing
end

function constrain_reservoir_ramps!(m,c,V;level=nothing)
    for (i,r) in enumerate(c.system.reservoirs)
        levels=nothing
        needs_level=any(ramp_allowance(c,r,attribute,c.grid[t],c.grid[t+1])!==nothing
            for attribute in (:level_ramp_up,:level_ramp_down), t in eachindex(c.prices))
        needs_level && level!==nothing && (levels=[level(i,t) for t in eachindex(c.grid)])
        for t in eachindex(c.prices), (attributes,values) in (
            ((:volume_ramp_up,:volume_ramp_down),view(V,i,:)),
            ((:level_ramp_up,:level_ramp_down),levels),
        )
            values===nothing && continue
            change=values[t+1]-values[t]
            up=ramp_allowance(c,r,attributes[1],c.grid[t],c.grid[t+1])
            down=ramp_allowance(c,r,attributes[2],c.grid[t],c.grid[t+1])
            scale=attributes[1]==:volume_ramp_up ? r.vmax : 250.0
            up===nothing || @constraint(m,change/scale<=up/scale)
            down===nothing || @constraint(m,-change/scale<=down/scale)
        end
    end
    nothing
end

function interval_ramp_violation(c,object,values,initial,historical_hours,attributes;
    su=nothing,sd=nothing,members=Int[])
    violation=0.0
    for t in eachindex(c.prices)
        previous=t==1 ? initial : values[t-1]
        previous===nothing && continue
        a,b=interval_midpoints(c.grid,t,historical_hours)
        up=ramp_allowance(c,object,attributes[1],a,b)
        down=ramp_allowance(c,object,attributes[2],a,b)
        ju,jd=su===nothing ? (0.0,0.0) : power_transition_allowances(c,members,t,su,sd)
        up===nothing || (violation=max(violation,values[t]-previous-up-ju))
        down===nothing || (violation=max(violation,previous-values[t]-down-jd))
    end
    violation
end

function schedule_violation(c,object,p,q)
    violation=0.0
    for t in eachindex(c.prices), (attribute,values) in ((:power,p),(:discharge,q))
        required=optional_operation(c,object.name,attribute,t)
        required===nothing || (violation=max(violation,abs(values[t]-required)))
    end
    violation
end

"""Numerical residuals for operating rows, evaluated from physical schedule values."""
function operating_residuals(c,x)
    residuals=Dict{String,Float64}()
    u,P,GQ,RQ,V=x["u"],x["power"],x["generator_q"],x["river_release"],x["V"]
    su,sd=unit_transitions(c,u)
    for (i,g) in enumerate(c.system.generators)
        p,q=view(P,i,:),view(GQ,i,:)
        residuals["schedule_$(g.name)"]=schedule_violation(c,g,p,q)
        residuals["unit_ramp_$(g.name)"]=interval_ramp_violation(c,g,p,g.initial_power,
            g.initial_interval_hours,(:ramp_up,:ramp_down);su,sd,members=[i])
        residuals["discharge_ramp_$(g.name)"]=interval_ramp_violation(c,g,q,g.initial_discharge,
            g.initial_interval_hours,(:discharge_ramp_up,:discharge_ramp_down))
    end
    for p in c.system.plants
        members=findall(g->g.plant==p.name,c.system.generators)
        tp=vec(sum(P[members,:];dims=1))
        tq=vec(sum(GQ[members,:];dims=1))
        capacity=0.0
        for t in eachindex(c.prices)
            on=any(u[j,t]==1 for j in members)
            capacity=max(capacity,tp[t]-opinterval(c,p.name,:pmax,t,p.pmax),
                opinterval(c,p.name,:pmin,t,p.pmin)*on-tp[t],
                opinterval(c,p.name,:qmin,t,p.qmin)*on-tq[t])
            maximum_q=optional_operation(c,p.name,:qmax,t;default=p.qmax)
            maximum_q===nothing || (capacity=max(capacity,tq[t]-maximum_q))
        end
        residuals["plant_capacity_$(p.name)"]=capacity
        residuals["schedule_$(p.name)"]=schedule_violation(c,p,tp,tq)
        initial_p,p_hours=plant_initial_output(c,p,members,:power)
        initial_q,q_hours=plant_initial_output(c,p,members,:discharge)
        residuals["plant_ramp_$(p.name)"]=interval_ramp_violation(c,p,tp,initial_p,
            p_hours,(:ramp_up,:ramp_down);su,sd,members)
        residuals["discharge_ramp_$(p.name)"]=interval_ramp_violation(c,p,tq,initial_q,
            q_hours,(:discharge_ramp_up,:discharge_ramp_down))
    end
    for (i,r) in enumerate(c.system.rivers)
        initial,duration=river_initial_output(r)
        residuals["river_ramp_$(r.name)"]=interval_ramp_violation(c,r,view(RQ,i,:),
            initial,duration,(:ramp_up,:ramp_down))
    end
    for (i,r) in enumerate(c.system.reservoirs)
        levels=head.(Ref(r),view(V,i,:))
        for (label,attributes,values) in (
            ("volume_ramp",(:volume_ramp_up,:volume_ramp_down),view(V,i,:)),
            ("level_ramp",(:level_ramp_up,:level_ramp_down),levels),
        )
            violation=0.0
            for t in eachindex(c.prices)
                up=ramp_allowance(c,r,attributes[1],c.grid[t],c.grid[t+1])
                down=ramp_allowance(c,r,attributes[2],c.grid[t],c.grid[t+1])
                up===nothing || (violation=max(violation,values[t+1]-values[t]-up))
                down===nothing || (violation=max(violation,values[t]-values[t+1]-down))
            end
            residuals["$(label)_$(r.name)"]=violation
        end
    end
    residuals
end

"""Finer replay checks interval operating rows on their original windows."""
function without_interval_controls(c)
    strip(object)=begin
        changes=Dict{Symbol,Any}(attribute=>nothing for attribute in RAMP_ATTRIBUTES
            if attribute in fieldnames(typeof(object)))
        object isa Plant && (changes[:ramp]=nothing)
        _river_replace(object;changes...)
    end
    system=_river_replace(c.system;
        reservoirs=strip.(c.system.reservoirs),
        generators=strip.(c.system.generators),
        plants=strip.(c.system.plants),
        rivers=strip.(c.system.rivers))
    operations=[z for z in c.operations if z.attribute ∉ RAMP_ATTRIBUTES &&
        z.attribute ∉ (:power,:discharge)]
    _river_replace(c;system,operations)
end

function operating_window_values(c,z,grid,indices,u)
    dt=diff(grid)
    original_dt=diff(c.grid)
    average(values)=[sum(values[i,k]*dt[k] for k in eachindex(dt) if indices[k]==t)/original_dt[t]
        for i in axes(values,1), t in eachindex(c.prices)]
    vertices=[something(findfirst(==(time),grid)) for time in c.grid]
    Dict("u"=>u,"power"=>average(z["power"]),
        "generator_q"=>average(z["generator_q"]),
        "river_release"=>average(z["river_release"]),"V"=>z["V"][:,vertices])
end
