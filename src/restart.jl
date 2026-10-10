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
    audit=validate(c, x)
    audit["valid"] || throw(
        ArgumentError(
            "restart source schedule fails independent validation: $(join(audit["errors"], "; "))",
        ),
    )
    s=c.system
    previous=edge-1
    reservoirs=[
        _restart_replace(r; v0 = Float64(x["V"][i, edge])) for
        (i, r) in enumerate(s.reservoirs)
    ]
    generators=Generator[]
    for (i, g) in enumerate(s.generators)
        state=g.initial_on
        since=first(c.grid)-g.initial_age
        for t in 1:previous
            if x["u"][i, t]!=state
                state=Int(x["u"][i, t])
                since=c.grid[t]
            end
        end
        push!(
            generators,
            _restart_replace(g; initial_on = state, initial_age = Float64(time-since),
                initial_power=Float64(x["power"][i,previous]),
                initial_discharge=Float64(x["generator_q"][i,previous]),
                initial_interval_hours=Float64(c.grid[edge]-c.grid[edge-1])),
        )
    end
    plant_history=_plant_commitment_history(c,x["u"],time)
    plants=[
        _restart_replace(
            p;
            initial_power = sum(
                x["power"][i, previous] for
                (i, g) in enumerate(s.generators) if g.plant==p.name;
                init = 0.0,
            ),
            initial_discharge=sum(x["generator_q"][i,previous] for
                (i,g) in enumerate(s.generators) if g.plant==p.name;init=0.0),
            initial_on=plant_history[j].state,
            initial_age=plant_history[j].age,
            initial_interval_hours = Float64(c.grid[edge]-c.grid[edge - 1]),
        ) for (j,p) in enumerate(s.plants)
    ]
    exact=all(r.deterministic_delay!==nothing for r in s.rivers)
    routed=exact ? route_network_exact(c, x["river_release"];generator_q=x["generator_q"],tunnel_q=x["tunnel_q"]) : nothing
    rivers=River[]
    for (i, r) in enumerate(s.rivers)
        windows=isempty(r.arrival_window_grid) ? Float64[] :
                Float64.(
            r.arrival_window_grid[findfirst(==(time), r.arrival_window_grid):end],
        )
        if exact
            cohorts=routed["release_cohorts"][i]
            executed=release_history(cohorts.grid, cohorts.release, time)
            history_grid=vcat(r.history_grid, executed.grid[2:end])
            history_release=vcat(r.history_release, executed.release)
        else
            history_grid=vcat(r.history_grid, c.grid[2:edge])
            history_release=vcat(
                r.history_release,
                Float64.(x["river_release"][i, 1:previous]),
            )
        end
        push!(
            rivers,
            _restart_replace(
                r;
                history_grid,
                history_release,
                arrival_window_grid = windows,
                initial_release=Float64(x["river_release"][i,previous]),
                initial_interval_hours=Float64(c.grid[edge]-c.grid[edge-1]),
            ),
        )
    end
    system=_restart_replace(s; reservoirs, generators, plants, rivers)
    result=ScheduleCase(
        name = c.name*"_restart_"*string(time),
        system = system,
        grid = Float64.(c.grid[edge:end]),
        prices = Float64.(c.prices[edge:end]),
        operations = c.operations,
        flow_requirements = c.flow_requirements,
    )
    validate_inputs(result)
    result
end
