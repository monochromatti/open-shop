# Reflect every field so restarts retain all declared constitutive and operating
# data. This is a state update; it does not construct a simplified watercourse.
function _restart_replace(x; kwargs...)
    values=Dict(f=>getfield(x, f) for f in fieldnames(typeof(x)))
    merge!(values, Dict(kwargs))
    typeof(x)(; values...)
end

"""Build a complete scheduling continuation at an interior time-grid boundary.

The source schedule must pass independent validation. Reservoir states, unit
state/dwell age and previous plant interval power are carried into the new case.
Each river keeps original historical and executed release cohorts, without
rebinning. Absolute times are retained. Declared arrival averaging windows must
also have a boundary at the restart time; partial-window accounting is rejected.
Flow-observation windows and reservoir level ranges can cross the restart;
their original observed intervals and storage vertices are retained as history.

The resulting objective measures incremental wealth from the restart state, so
its constant stock/transit reference differs from the original full horizon.
"""
function restart_case(c::ScheduleCase, x, time::Real)
    isfinite(time) || throw(ArgumentError("restart time must be finite"))
    edge=findfirst(==(time), c.grid)
    edge!==nothing && 1<edge<length(c.grid) || throw(
        ArgumentError("restart time must be a strict interior scheduling-grid boundary"),
    )
    for r in c.system.rivers
        windows=isempty(r.arrival_window_grid) ? c.grid : r.arrival_window_grid
        time in windows || throw(
            ArgumentError(
                "restart crosses declared arrival averaging window for $(r.name)",
            ),
        )
    end
    for rule in c.soft_constraints
        rule.end_time===nothing || !(rule.time<time<rule.end_time) ||
            throw(ArgumentError("restart crosses soft assessment window $(rule.name)"))
    end
    audit=validate(c, x)
    audit["valid"] || throw(
        ArgumentError(
            "restart source schedule fails independent validation: $(join(audit["errors"], "; "))",
        ),
    )
    x=_schedule_input(x; case=c)
    s=c.system
    previous=edge-1
    reservoirs=[begin
        needs_history=!isempty(r.level_windows) || !isempty(r.history_grid)
        first_new=isempty(r.history_grid) ? 1 : 2
        history_grid=needs_history ? vcat(r.history_grid,c.grid[first_new:edge]) : Float64[]
        history_storage=needs_history ? vcat(r.history_storage,
            vec(x.storage[i,first_new:edge])) : Float64[]
        _restart_replace(r;v0=Float64(x.storage[i,edge]),history_grid,history_storage,
            level_windows=[w for w in r.level_windows if w.end_time>time])
    end for (i,r) in enumerate(s.reservoirs)]
    generators=Generator[]
    for (i, g) in enumerate(s.generators)
        state=g.initial_on
        since=first(c.grid)-g.initial_age
        for t in 1:previous
            if x.commitment[i, t]!=state
                state=Int(x.commitment[i, t])
                since=c.grid[t]
            end
        end
        push!(
            generators,
            _restart_replace(g; initial_on = state, initial_age = Float64(time-since),
                initial_power=Float64(x.power[i,previous]),
                initial_discharge=Float64(x.generator_q[i,previous]),
                initial_interval_hours=Float64(c.grid[edge]-c.grid[edge-1])),
        )
    end
    opening=_physical_tunnel_openings(c,get(x,"tunnel_opening",nothing))
    tunnels=[_restart_replace(e;initial_opening=Float64(opening[i,previous]))
        for (i,e) in enumerate(s.tunnels)]
    switch_history=tunnel_switch_history(c,opening,time)
    tunnel_interlocks=[_restart_replace(rule;
        initial_forward_wait=history.forward_wait,initial_backward_wait=history.backward_wait)
        for (rule,history) in zip(s.tunnel_interlocks,switch_history)]
    plant_history=_plant_commitment_history(c,x.commitment,time)
    plants=[
        _restart_replace(
            p;
            initial_power = sum(
                x.power[i, previous] for
                (i, g) in enumerate(s.generators) if g.plant==p.name;
                init = 0.0,
            ),
            initial_discharge=sum(x.generator_q[i,previous] for
                (i,g) in enumerate(s.generators) if g.plant==p.name;init=0.0),
            initial_on=plant_history[j].state,
            initial_age=plant_history[j].age,
            initial_interval_hours = Float64(c.grid[edge]-c.grid[edge - 1]),
        ) for (j,p) in enumerate(s.plants)
    ]
    routed=routed_schedule(c,x.river_release;generator_q=x.generator_q,tunnel_q=x.tunnel_q)
    rivers=River[]
    for (i, r) in enumerate(s.rivers)
        windows=isempty(r.arrival_window_grid) ? Float64[] :
                Float64.(
            r.arrival_window_grid[findfirst(==(time), r.arrival_window_grid):end],
        )
        if is_channel(r)
            history_grid=[Float64(time)-1.0,Float64(time)]
            history_release=[0.0]
        else
            cohorts=routed.cohorts[i]
            executed=release_history(cohorts.grid, cohorts.release, time)
            history_grid=vcat(r.history_grid, executed.grid[2:end])
            history_release=vcat(r.history_release, executed.release)
        end
        push!(
            rivers,
            _restart_replace(
                r;
                history_grid,
                history_release,
                arrival_window_grid = windows,
                initial_release=Float64(x.river_release[i,previous]),
                initial_interval_hours=Float64(c.grid[edge]-c.grid[edge-1]),
            ),
        )
    end
    flow_requirements=[begin
        observed=observed_flow(c,rule,x.generator_q,x.river_release,x.arrival_volume,x.tunnel_q)
        counted=sum(0.0036*(c.grid[t+1]-c.grid[t])*observed[t] for t in 1:previous;init=0.0)
        needs_history=!isempty(rule.volume_windows) || !isempty(rule.rolling_ramps) ||
            !isempty(rule.period_ramps) || !isempty(rule.history_grid)
        first_new=isempty(rule.history_grid) ? 1 : 2
        history_grid=needs_history ? vcat(rule.history_grid,
            c.grid[first_new:edge]) : Float64[]
        history_flow=needs_history ? vcat(rule.history_flow,observed[1:previous]) : Float64[]
        _restart_replace(rule;
            initial_volume=rule.min_volume===nothing && rule.max_volume===nothing ?
                rule.initial_volume : rule.initial_volume+counted,
            history_grid,history_flow,
            volume_windows=[w for w in rule.volume_windows if w.end_time>time])
    end for rule in c.flow_requirements]
    system=_restart_replace(s; reservoirs, generators, plants, rivers, tunnels,tunnel_interlocks)
    result=ScheduleCase(
        name = c.name*"_restart_"*string(time),
        system = system,
        grid = Float64.(c.grid[edge:end]),
        prices = Float64.(c.prices[edge:end]),
        operations = c.operations,
        flow_requirements = flow_requirements,
        soft_constraints = [rule for rule in c.soft_constraints
            if (rule.end_time===nothing ? rule.time : rule.end_time)>time],
        reserve_groups = c.reserve_groups,
        routing_grid = isempty(c.routing_grid) ? Float64[] : filter(>=(time),c.routing_grid),
    )
    validate_inputs(result)
    result
end
