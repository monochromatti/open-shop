"""Independently audit equations, operating limits, routing, and chronological replay.

Residuals use physical units except objective (relative error). Terminal minimum
up/down obligations are reported for continuation rather than silently discarded.
"""
function validate(c::ScheduleCase, result; tolerance = 1e-4, transport = nothing)
    errors=String[]
    residuals=Dict{String,Float64}()
    obligations=Any[]
    finish() = Dict(
        "valid"=>isempty(errors),
        "max_residual"=>maximum(values(residuals); init = 0.0),
        "residuals"=>residuals,
        "errors"=>errors,
        "terminal_commitment_obligations"=>obligations,
    )
    record(name, x) = begin
        v=maximum(abs, x; init = 0.0)
        residuals[name]=v
        (!isfinite(v) || v>tolerance) &&
            push!(errors, "$name residual $v exceeds $tolerance")
    end
    bound(name, x) = record(name, max.(x, 0.0))
    s=c.system
    T=length(c.prices)
    dt=diff(c.grid)
    R=length(s.reservoirs)
    G=length(s.generators)
    D=length(s.rivers)
    ix=nodeindex(s)
    dimensions=Dict(
        "V"=>(R, T+1),
        "H"=>(length(nodes(s)), T),
        "tunnel_q"=>(length(s.tunnels), T),
        "generator_q"=>(G, T),
        "power"=>(G, T),
        "river_release"=>(D, T),
        "gate"=>(D, T),
        "arrival_volume"=>(D, T),
        "u"=>(G, T),
        "terminal_transit"=>(D,),
    )
    for (key, shape) in dimensions
        if !haskey(result, key) || size(result[key])!=shape
            push!(errors, "missing or malformed $key")
            return finish()
        end
        if !all(isfinite, result[key])
            push!(errors, "nonfinite $key")
            return finish()
        end
    end
    if !haskey(result, "objective") || !isfinite(result["objective"])
        push!(errors, "missing or nonfinite objective")
        return finish()
    end
    V=result["V"]
    H=result["H"]
    Q=result["tunnel_q"]
    GQ=result["generator_q"]
    P=result["power"]
    RQ=result["river_release"]
    A=result["arrival_volume"]
    ga=result["gate"]
    u=result["u"]
    admissible(c, u) ||
        push!(errors, "invalid binary commitment or minimum up/down history")
    record("initial_storage", V[:, 1]-[r.v0 for r in s.reservoirs])
    for (i, r) in enumerate(s.reservoirs)
        bound(
            "storage_$(r.name)",
            vcat(
                [storage_bounds(c, r, c.grid[t])[1]-V[i, t] for t in 1:(T + 1)],
                [V[i, t]-storage_bounds(c, r, c.grid[t])[2] for t in 1:(T + 1)],
            ),
        )
        record(
            "head_storage_$(r.name)",
            [H[i, t]-head(r, (V[i, t]+V[i, t + 1])/2) for t in 1:T],
        )
    end
    for (i, j) in enumerate(s.junctions)
        bound("junction_head_$(j.name)", vcat(j.hmin .- H[R + i, :], H[R + i, :] .- j.hmax))
    end
    for b in s.boundaries
        record("boundary_head_$(b.name)", H[ix[b.name], :] .- b.head)
    end
    for (i, e) in enumerate(s.tunnels)
        bound(
            "tunnel_capacity_$(e.name)",
            [abs(Q[i, t])-opinterval(c, e.name, :capacity, t, e.capacity) for t in 1:T],
        )
        record(
            "tunnel_loss_$(e.name)",
            [
                opinterval(c, e.name, :opening, t, e.opening)==0 ? Q[i, t] :
                opinterval(c, e.name, :opening, t, e.opening)*(
                    H[ix[e.source], t]-H[ix[e.target], t]
                )-e.resistance*Q[i, t]*abs(Q[i, t]) for t in 1:T
            ],
        )
    end
    for (i, g) in enumerate(s.generators)
        pl=plantof(s, g)
        hd=H[ix[pl.source], :]-outlet_head.(Ref(pl), H[ix[pl.target], :])-[
            tailwater(
                pl,
                sum(GQ[j, t] for (j, z) in enumerate(s.generators) if z.plant==g.plant),
            ) for t in 1:T
        ]
        bound(
            "unit_bounds_$(g.name)",
            vcat(
                [opinterval(c, g.name, :qmin, t, g.qmin)*u[i, t]-GQ[i, t] for t in 1:T],
                [GQ[i, t]-opinterval(c, g.name, :qmax, t, g.qmax)*u[i, t] for t in 1:T],
                [opinterval(c, g.name, :pmin, t, g.pmin)*u[i, t]-P[i, t] for t in 1:T],
                [P[i, t]-opinterval(c, g.name, :pmax, t, g.pmax)*u[i, t] for t in 1:T],
            ),
        )
        record(
            "unit_conversion_$(g.name)",
            [P[i, t]-power(g, GQ[i, t], hd[t]) for t in 1:T],
        )
        on=findall(==(1), u[i, :])
        eta=[efficiency(g, GQ[i, t], hd[t]) for t in on]
        bound(
            "unit_head_efficiency_$(g.name)",
            vcat(g.hmin .- hd[on], hd[on] .- g.hmax, g.min_efficiency .- eta, eta .- 1),
        )
        if g.turbine_table!==nothing
            bound(
                "turbine_envelope_$(g.name)",
                vcat(
                    [
                        turbine_qmin(g.turbine_table, hd[t]; extrapolation = :linear)-GQ[
                            i,
                            t,
                        ] for t in on
                    ],
                    [
                        GQ[i, t]-turbine_qmax(
                            g.turbine_table,
                            hd[t];
                            extrapolation = :linear,
                        ) for t in on
                    ],
                ),
            )
        end
        state=g.initial_on
        since=c.grid[1]-g.initial_age
        for t in 1:T
            if u[i, t]!=state
                state=u[i, t]
                since=c.grid[t]
            end
        end
        left=max(0.0, (state==1 ? g.minup : g.mindown)-(last(c.grid)-since))
        left>0 && push!(
            obligations,
            Dict("generator"=>string(g.name), "state"=>state, "remaining_hours"=>left),
        )
    end
    for pl in s.plants
        ids=findall(g->g.plant==pl.name, s.generators)
        total=vec(sum(P[ids, :]; dims = 1))
        bound(
            "plant_capacity_$(pl.name)",
            [total[t]-opinterval(c, pl.name, :pmax, t, pl.pmax) for t in 1:T],
        )
        violations=Float64[]
        if pl.initial_power!==nothing
            base=pl.ramp*(pl.initial_interval_hours+dt[1])/2
            su=sum(
                opinterval(c, s.generators[j].name, :pmin, 1, s.generators[j].pmin)*max(
                    0,
                    u[j, 1]-s.generators[j].initial_on,
                ) for j in ids
            )
            sd=sum(
                opinterval(c, s.generators[j].name, :pmin, 1, s.generators[j].pmin)*max(
                    0,
                    s.generators[j].initial_on-u[j, 1],
                ) for j in ids
            )
            push!(
                violations,
                total[1]-pl.initial_power-base-su,
                pl.initial_power-total[1]-base-sd,
            )
        end
        for t in 2:T
            base=pl.ramp*(dt[t]+dt[t - 1])/2
            su=sum(
                opinterval(c, s.generators[j].name, :pmin, t, s.generators[j].pmin)*max(
                    0,
                    u[j, t]-u[j, t - 1],
                ) for j in ids
            )
            sd=sum(
                opinterval(c, s.generators[j].name, :pmin, t, s.generators[j].pmin)*max(
                    0,
                    u[j, t - 1]-u[j, t],
                ) for j in ids
            )
            push!(violations, total[t]-total[t - 1]-base-su, total[t - 1]-total[t]-base-sd)
        end
        bound("plant_ramp_$(pl.name)", violations)
    end
    bound("gate_bounds", vcat(vec(-ga), vec(ga .- 1)))
    exact=all(r.deterministic_delay!==nothing for r in s.rivers)
    nd=exact ? _transport_data(c, transport) : nothing
    !exact &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed validation",
            ),
        )
    for (d, r) in enumerate(s.rivers)
        bound(
            "gate_min_$(r.name)",
            [opinterval(c, r.name, :gate_min, t, r.gate_min)-ga[d, t] for t in 1:T],
        )
        bound(
            "gate_max_$(r.name)",
            [ga[d, t]-opinterval(c, r.name, :gate_max, t, 1.0) for t in 1:T],
        )
        bound(
            "river_capacity_$(r.name)",
            vcat(
                -RQ[d, :],
                [RQ[d, t]-opinterval(c, r.name, :capacity, t, r.capacity) for t in 1:T],
            ),
        )
        windows=isempty(r.arrival_window_grid) ? c.grid : r.arrival_window_grid
        averaged=exact ?
                 [
            _transport_value(
                _pulse_expression(nd.arrival_pulses[d], windows[k], windows[k + 1]),
                RQ,
            ) for k in 1:(length(windows) - 1)
        ] :
                 windows==c.grid ? vec(A[d, :]) :
                 route_volumes(r, windows, c.grid, vec(RQ[d, :]))+route_volumes(
            r,
            windows,
            r.history_grid,
            r.history_release,
        )
        bound(
            "environmental_arrival_$(r.name)",
            [
                opaverage(
                    c,
                    r.name,
                    :min_arrival,
                    windows[k],
                    windows[k + 1],
                    r.min_arrival,
                )-averaged[k]/(0.0036*(windows[k + 1]-windows[k])) for
                k in 1:(length(windows) - 1)
            ],
        )
        if exact
            for kind in (:arrival, :release)
                knots=sort!(unique!(vcat(deterministic_knots(nd, d; kind), c.grid)))
                deviations=Float64[]
                for time in knots,
                    side in (
                        time==first(c.grid) ? (:right,) :
                        time==last(c.grid) ? (:left,) : (:left, :right)
                    )

                    rate=_transport_value(
                        deterministic_point_data(nd, d, time; side, kind),
                        RQ,
                    )
                    kind==:release && push!(
                        deviations,
                        rate-opvalue(c, r.name, :capacity, time, r.capacity; side),
                    )
                    kind==:arrival &&
                        r.arrival_policy==:pointwise &&
                        push!(
                            deviations,
                            opvalue(c, r.name, :min_arrival, time, r.min_arrival; side)-rate,
                        )
                end
                bound("exact_$(kind)_$(r.name)", deviations)
            end
        elseif r.arrival_policy==:pointwise
            bound(
                "pointwise_arrival_$(r.name)",
                [arrival_requirement_violation(c, r, c.grid, vec(RQ[d, :]))],
            )
        end
        if r.law!=:junction
            if !r.allow_dry && r.discharge_curve===nothing && r.law!=:controlled
                bound("river_law_domain_$(r.name)", r.crest .- H[ix[r.source], :])
            end
            expected=[
                river_law_value(
                    r.law==:controlled ?
                    _river_replace(
                        r;
                        capacity = opinterval(c, r.name, :capacity, t, r.capacity),
                    ) : r,
                    H[ix[r.source], t],
                    ga[d, t],
                ) for t in 1:T
            ]
            record("river_law_$(r.name)", RQ[d, :]-expected)
            r.law==:weir && record("weir_gate_$(r.name)", ga[d, :] .- 1)
        else
            record("junction_gate_$(r.name)", ga[d, :])
        end
    end
    requirements=release_requirements(c, RQ)
    for rule in c.flow_requirements
        gs, rs=flow_requirement_indices(c, rule)
        bound("flow_requirement_$(rule.name)", [
            opinterval(c, rule.name, :min_flow, t, rule.min_flow)-
            opinterval(c, rule.name, :inflow, t, rule.inflow)-
            sum(GQ[i, t] for i in gs; init=0.0)-
            sum(RQ[i, t] for i in rs; init=0.0) for t in 1:T
        ])
    end
    bound("hard_release_requirements", requirements.hard)
    if any(opinterval(c, r.name, :release_penalty, t, 0.0)>0 for r in s.rivers, t in 1:T)
        if !haskey(result, "shortfall_release") ||
           size(result["shortfall_release"])!=(D, T) ||
           !haskey(result, "release_penalty_cost")
            push!(errors, "missing soft release accounting")
        else
            record("release_shortfall", result["shortfall_release"]-requirements.shortfall)
            record(
                "release_penalty_cost",
                [
                    (result["release_penalty_cost"]-requirements.cost)/max(
                        1.0,
                        abs(requirements.cost),
                    ),
                ],
            )
        end
    end
    for j in s.river_junctions
        outgoing=only(findall(r->r.source==j.name, s.rivers))
        incoming=findall(r->r.target==j.name, s.rivers)
        record(
            "river_junction_$(j.name)",
            0.0036 .* dt .* RQ[outgoing, :]-vec(sum(A[incoming, :]; dims = 1)),
        )
    end
    for name in nodes(s)[1:(R + length(s.junctions))]
        net=zeros(T)
        for (d, e) in enumerate(s.tunnels)
            net .+= ((e.target==name)-(e.source==name)) .* Q[d, :]
        end
        for (d, g) in enumerate(s.generators)
            p=plantof(s, g)
            net .+= ((p.target==name)-(p.source==name)) .* GQ[d, :]
        end
        for (d, r) in enumerate(s.rivers)
            r.target==name && (net .+= A[d, :] ./ (0.0036 .* dt))
            r.source==name && (net .-= RQ[d, :])
        end
        i=ix[name]
        record(
            "balance_$name",
            i<=R ?
            diff(V[i, :])-0.0036 .* dt .* (
                net .+ [
                    opinterval(c, s.reservoirs[i].name, :inflow, t, s.reservoirs[i].inflow)
                    for t in 1:T
                ]
            ) : net,
        )
    end
    try
        routed=exact ? route_network_exact(c, RQ; transport = nd) : route_network(c, RQ)
        record("routing_release", routed["release"]-RQ)
        record("routing_arrival", routed["arrival_volume"]-A)
        record("terminal_transit", routed["transit"][:, end]-result["terminal_transit"])
        W=routed["transit"]
        external=zeros(T)
        for b in s.boundaries
            for (d, e) in enumerate(s.tunnels)
                external .+= ((e.target==b.name)-(e.source==b.name)) .* Q[d, :]
            end
            for (d, g) in enumerate(s.generators)
                p=plantof(s, g)
                external .+= ((p.target==b.name)-(p.source==b.name)) .* GQ[d, :]
            end
            for (d, r) in enumerate(s.rivers)
                r.target==b.name && (external .+= A[d, :] ./ (0.0036 .* dt))
            end
        end
        initial=sum(V[:, 1])+sum(W[:, 1])
        inflow=[
            sum(opinterval(c, r.name, :inflow, t, r.inflow) for r in s.reservoirs) for
            t in 1:T
        ]
        # One column reduction and prefix sum instead of allocating and summing
        # every time prefix separately (quadratic work in horizon length).
        inventory=vec(sum(V; dims = 1)) .+ vec(sum(W; dims = 1))
        exchange=cumsum(dt .* (external .- inflow))
        record("global_water", (@view inventory[2:end]) .- initial .+ 0.0036 .* exchange)
        startup=sum(
            g.startup*max(0, u[i, t]-(t==1 ? g.initial_on : u[i, t - 1])) for
            (i, g) in enumerate(s.generators), t in 1:T;
            init = 0.0,
        )
        shutdown=sum(
            g.shutdown*max(0, (t==1 ? g.initial_on : u[i, t - 1])-u[i, t]) for
            (i, g) in enumerate(s.generators), t in 1:T;
            init = 0.0,
        )
        objective=sum(c.prices[t]*dt[t]*P[i, t] for i in 1:G, t in 1:T; init = 0.0)-startup-shutdown +
                  sum(
                      r.water_value*(V[i, end]-r.v0) for (i, r) in enumerate(s.reservoirs)
                  ) +
                  sum(
                      r.water_value*(W[d, end]-W[d, 1]) for (d, r) in enumerate(s.rivers);
                      init = 0.0,
                  )
        objective-=requirements.cost
        record(
            "objective_relative",
            [(objective-result["objective"])/max(1.0, abs(objective))],
        )
    catch e
        push!(errors, "independent routing audit failed: $(sprint(showerror,e))")
    end
    try
        replay=simulate(c, GQ, ga; transport = nd)
        replay["converged"] || push!(
            errors,
            "independent forward replay did not converge in the river-law domain",
        )
        record("replay_residual", [replay["residual"]])
        for key in ["V", "H", "tunnel_q", "power", "river_release", "arrival_volume"]
            record("replay_$key", replay[key]-result[key])
        end
    catch e
        push!(errors, "independent forward replay failed: $(sprint(showerror,e))")
    end
    finish()
end
