function admissible(c, u)
    size(u)==(length(c.system.generators), length(c.prices)) || return false
    all(x->x==0 || x==1, u) || return false
    for (j, g) in enumerate(c.system.generators)
        state=g.initial_on
        since=c.grid[1]-g.initial_age
        for t in axes(u, 2)
            forced=opinterval(c, g.name, :forced_on, t, -1.0)
            forced>=0 && u[j, t]!=forced && return false
            if u[j, t]!=state
                c.grid[t]-since+1e-8 >= (state==1 ? g.minup : g.mindown) || return false
                state=u[j, t]
                since=c.grid[t]
            end
        end
    end
    true
end
function _build_dispatch(
    c;
    u = ones(Int, length(c.system.generators), length(c.prices)),
    solver = "Ipopt",
    seed = 1,
    warm = nothing,
    arrival_margin = 0.0,
    operational_margin = 0.0,
    time_limit = 45.0,
    feasibility_only = false,
    joint = false,
    relaxed = false,
    fixed_u = nothing,
    free_mask = nothing,
    transport = nothing,
)
    joint && throw(
        ArgumentError(
            "local dispatch requires fixed commitment; use the global solver for joint scheduling",
        ),
    )
    starttime=time()
    isfinite(arrival_margin) && arrival_margin>=0 ||
        throw(ArgumentError("invalid arrival margin"))
    isfinite(operational_margin) && operational_margin>=0 ||
        throw(ArgumentError("invalid operational margin"))
    isfinite(time_limit) && time_limit>0 ||
        throw(ArgumentError("time limit must be finite and positive"))
    validate_inputs(c)
    joint || admissible(c, u)||error("Invalid commitment")
    s=c.system
    node_names=nodes(s)
    T=length(c.prices)
    R=length(s.reservoirs)
    N=length(node_names)
    E=length(s.tunnels)
    G=length(s.generators)
    D=length(s.rivers)
    ix=nodeindex(s)
    dt=diff(c.grid)
    generator_plants=[plantof(s, g) for g in s.generators]
    plant_generators=Dict(
        plant.name=>findall(g->g.plant==plant.name, s.generators) for plant in s.plants
    )
    tunnel_incidence=[
        [
            (j, Int(e.target==name)-Int(e.source==name)) for
            (j, e) in enumerate(s.tunnels) if e.target==name || e.source==name
        ] for name in node_names
    ]
    generator_incidence=[
        [
            (j, Int(plant.target==name)-Int(plant.source==name)) for
            (j, plant) in enumerate(generator_plants) if
            plant.target==name || plant.source==name
        ] for name in node_names
    ]
    river_incidence=[
        [
            (j, r.target==name, r.source==name) for
            (j, r) in enumerate(s.rivers) if r.target==name || r.source==name
        ] for name in node_names
    ]
    river_junction_incidence=[
        (
            only(findall(r->r.source==j.name, s.rivers)),
            findall(r->r.target==j.name, s.rivers),
        ) for j in s.river_junctions
    ]
    m=Model()
    set_silent(m)
    ustart=copy(u)
    transitions=nothing
    if joint
        states=_joint_states!(m, c; relaxed, fixed_u, free_mask, incumbent = warm)
        u=states.u
        transitions=states
        ustart=states.start
    end
    @variable(m, v[1:R, 1:(T + 1)])
    @variable(m, h[1:N, 1:T])
    @variable(m, q[1:E, 1:T])
    @variable(m, gq[1:G, 1:T]>=0)
    @variable(m, p[1:G, 1:T]>=0)
    @variable(m, rq[1:D, 1:T]>=0)
    @variable(m, 0<=a[1:D, 1:T]<=1)
    @variable(m, shortfall_release[1:D, 1:T]>=0)
    level_ops=[
        r.level_curve===nothing ? nothing :
        table_operator(m, r.level_curve; name = Symbol("level_", i)) for
        (i, r) in enumerate(s.reservoirs)
    ]
    turbine_ops=[
        g.turbine_table===nothing ? nothing :
        turbine_operator(m, g.turbine_table; name = Symbol("turbine_", i)) for
        (i, g) in enumerate(s.generators)
    ]
    electrical_ops=[
        g.generator_efficiency_curve===nothing ? nothing :
        table_operator(m, g.generator_efficiency_curve; name = Symbol("electrical_", i)) for
        (i, g) in enumerate(s.generators)
    ]
    tailwater_ops=Dict(
        p.name=>(
            p.tailwater_curve===nothing ? nothing :
            table_operator(m, p.tailwater_curve; name = Symbol("tailwater_", p.name))
        ) for p in s.plants
    )
    river_ops=[
        r.discharge_curve===nothing ? nothing :
        table_operator(m, r.discharge_curve; name = Symbol("river_law_", i)) for
        (i, r) in enumerate(s.rivers)
    ]
    turbine_flow_ops=[
        if !joint && g.turbine_table!==nothing && any(t->u[i, t]!=0, 1:T)
            table=g.turbine_table
            (
                table_operator(
                    m,
                    TableCurve(table.heads, table.qmin);
                    name = Symbol("qlo_", i),
                ),
                table_operator(
                    m,
                    TableCurve(table.heads, table.qmax);
                    name = Symbol("qhi_", i),
                ),
            )
        else
            nothing
        end for (i, g) in enumerate(s.generators)
    ]
    V=[s.reservoirs[i].vmax*v[i, t] for i in 1:R, t in 1:(T + 1)]
    H=250 .* h
    Q=50 .* q
    GQ=50 .* gq
    P=40 .* p
    RQ=100 .* rq
    constrain_flow_requirements!(m, c, GQ, RQ)
    for (i, r) in enumerate(s.reservoirs)
        for t in 1:(T + 1)
            lo, hi=storage_bounds(c, r, c.grid[t])
            set_lower_bound(v[i, t], lo/r.vmax)
            set_upper_bound(v[i, t], hi/r.vmax)
            set_start_value(v[i, t], r.v0/r.vmax)
        end
        fix(v[i, 1], r.v0/r.vmax; force = true)
        for t in 1:T
            set_lower_bound(h[i, t], head(r, r.vmin)/250)
            set_upper_bound(h[i, t], head(r, r.vmax)/250)
            mid=(V[i, t]+V[i, t + 1])/2
            level=level_ops[i]===nothing ? r.z0+r.slope*mid+r.curvature*mid^2 :
                  level_ops[i](mid)
            @constraint(m, (H[i, t]-level)/250==0)
        end
    end
    for (i, j) in enumerate(s.junctions), t in 1:T
        set_lower_bound(h[R + i, t], j.hmin/250)
        set_upper_bound(h[R + i, t], j.hmax/250)
        set_start_value(h[R + i, t], (j.hmin+j.hmax)/500)
    end
    for (i, b) in enumerate(s.boundaries), t in 1:T
        fix(h[R + length(s.junctions) + i, t], b.head/250; force = true)
    end
    for (i, e) in enumerate(s.tunnels), t in 1:T
        cap=opinterval(c, e.name, :capacity, t, e.capacity)
        opening=opinterval(c, e.name, :opening, t, e.opening)
        set_lower_bound(q[i, t], -cap/50)
        set_upper_bound(q[i, t], cap/50)
        if opening==0
            fix(q[i, t], 0; force = true)
        else
            @constraint(
                m,
                (
                    opening*(H[ix[e.source], t]-H[ix[e.target], t])-e.resistance*Q[i, t]*abs(
                        Q[i, t],
                    )
                )/100==0
            )
        end
    end
    outlet_ops=Dict(
        plant.name=>begin
            floor=plant.outlet_head_floor
            receiver=h[ix[plant.target], 1]
            lo=250*(is_fixed(receiver) ? fix_value(receiver) : lower_bound(receiver))
            hi=250*(is_fixed(receiver) ? fix_value(receiver) : upper_bound(receiver))
            floor===nothing || lo>=floor ? nothing :
            hi<=floor ? floor :
            table_operator(
                m,
                TableCurve([lo, floor, hi], [floor, floor, hi]);
                name = Symbol("outlet_head_", plant.name),
            )
        end for plant in s.plants
    )
    for (i, g) in enumerate(s.generators), t in 1:T
        plant=generator_plants[i]
        totalq=sum(GQ[j, t] for j in plant_generators[plant.name])
        tail=tailwater_ops[plant.name]===nothing ? 0.0 : tailwater_ops[plant.name](totalq)
        op=outlet_ops[plant.name]
        receiver=op===nothing ? H[ix[plant.target], t] :
                 op isa Real ? op : op(H[ix[plant.target], t])
        hd=H[ix[plant.source], t]-receiver-tail
        eta=turbine_ops[i]===nothing ?
            g.efficiency-g.qcurvature*((GQ[i, t]-g.qbest)/g.qbest)^2-g.hcurvature*(
            (hd-g.hbest)/g.hbest
        )^2 : turbine_ops[i](GQ[i, t], hd)
        electrical=electrical_ops[i]===nothing ? 1.0 : electrical_ops[i](P[i, t])
        qmin=opinterval(c, g.name, :qmin, t, g.qmin)
        qmax=opinterval(c, g.name, :qmax, t, g.qmax)
        pmin=opinterval(c, g.name, :pmin, t, g.pmin)
        pmax=opinterval(c, g.name, :pmax, t, g.pmax)
        begin
            set_lower_bound(gq[i, t], qmin*u[i, t]/50)
            set_upper_bound(gq[i, t], qmax*u[i, t]/50)
            set_lower_bound(p[i, t], (pmin+operational_margin)*u[i, t]/40)
            set_upper_bound(p[i, t], max(0.0, pmax-operational_margin)*u[i, t]/40)
            if u[i, t]==0
                fix(gq[i, t], 0; force = true)
                fix(p[i, t], 0; force = true)
            else
                @constraint(m, (P[i, t]-0.00981*GQ[i, t]*hd*eta*electrical)/40==0)
                @constraint(m, g.hmin<=hd<=g.hmax)
                @constraint(m, g.min_efficiency<=eta<=1.0)
                if g.turbine_table!==nothing
                    qlo, qhi=turbine_flow_ops[i]
                    flow_margin=operational_margin/(
                        0.00981*g.hbest*max(g.min_efficiency, g.efficiency, eps(Float64))
                    )
                    @constraint(m, GQ[i, t]>=qlo(hd)+flow_margin)
                    @constraint(m, GQ[i, t]<=qhi(hd)-flow_margin)
                end
            end
        end
    end
    for plant in s.plants
        ids=plant_generators[plant.name]
        for t in 1:T
            @constraint(
                m,
                sum(P[j, t] for j in ids)<=max(
                    0.0,
                    opinterval(c, plant.name, :pmax, t, plant.pmax)-operational_margin,
                )
            )
            if t==1 && plant.initial_power!==nothing
                ramp=max(
                    0.0,
                    plant.ramp*(plant.initial_interval_hours+dt[t])/2-operational_margin,
                )
                su=sum(
                    opinterval(c, s.generators[j].name, :pmin, t, s.generators[j].pmin)*(
                        joint ? transitions.su[j, t] :
                        max(0, u[j, t]-s.generators[j].initial_on)
                    ) for j in ids
                )
                sd=sum(
                    opinterval(c, s.generators[j].name, :pmin, t, s.generators[j].pmin)*(
                        joint ? transitions.sd[j, t] :
                        max(0, s.generators[j].initial_on-u[j, t])
                    ) for j in ids
                )
                @constraint(m, sum(P[j, t] for j in ids)-plant.initial_power>=-ramp-sd)
                @constraint(m, sum(P[j, t] for j in ids)-plant.initial_power<=ramp+su)
            elseif t>1
                ramp=max(0.0, plant.ramp*(dt[t - 1]+dt[t])/2-operational_margin)
                su=sum(
                    opinterval(c, s.generators[j].name, :pmin, t, s.generators[j].pmin)*(
                        joint ? transitions.su[j, t] : max(0, u[j, t]-u[j, t - 1])
                    ) for j in ids
                )
                sd=sum(
                    opinterval(c, s.generators[j].name, :pmin, t, s.generators[j].pmin)*(
                        joint ? transitions.sd[j, t] : max(0, u[j, t - 1]-u[j, t])
                    ) for j in ids
                )
                @constraint(m, sum(P[j, t]-P[j, t - 1] for j in ids)>=-ramp-sd)
                @constraint(m, sum(P[j, t]-P[j, t - 1] for j in ids)<=ramp+su)
            end
        end
    end
    exact=all(r.deterministic_delay!==nothing for r in s.rivers)
    nd=exact ? _transport_data(c, transport) : nothing
    !exact &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed dispatch",
            ),
        )
    rd=exact ? nothing : routing_data(c)
    arrivals=Matrix{Any}(undef, D, T)
    terminal=Any[]
    for (i, r) in enumerate(s.rivers)
        B=exact ? nothing : rd["B"][i]
        for t in 1:T
            cap=opinterval(c, r.name, :capacity, t, r.capacity)
            set_upper_bound(rq[i, t], cap/100)
            set_lower_bound(a[i, t], opinterval(c, r.name, :gate_min, t, r.gate_min))
            set_upper_bound(a[i, t], opinterval(c, r.name, :gate_max, t, 1.0))
            requirement=opinterval(c, r.name, :min_release, t, 0.0)
            penalty=opinterval(c, r.name, :release_penalty, t, 0.0)
            if penalty>0
                @constraint(
                    m,
                    (shortfall_release[i, t]-0.0036*dt[t]*(requirement-RQ[i, t]))/0.3>=0
                )
            else
                fix(shortfall_release[i, t], 0.0; force = true)
                requirement>0 && @constraint(m, RQ[i, t]>=requirement)
            end
            if exact
                arrivals[i, t]=nd.arrival_history[i, t]+sum(
                    v*RQ[d, k] for (d, k, v) in nd.arrival_terms[i, t];
                    init = 0.0,
                )
            else
                # Skip structural zeros to preserve banded routing sparsity.
                ks=[k for k in 1:T if B[t, k, 1]!=0 || B[t, k, 2]!=0]
                arrivals[i, t]=rd["history_arrival"][i, t]+sum(
                    0.0036*dt[k]*RQ[i, k]*(
                        B[t, k, 1]+RQ[i, k]/r.capacity*(B[t, k, 2]-B[t, k, 1])
                    ) for k in ks;
                    init = 0.0,
                )
            end
            if r.arrival_policy==:pointwise || isempty(r.arrival_window_grid)
                @constraint(
                    m,
                    arrivals[i, t]>=0.0036*dt[t]*opinterval(
                        c,
                        r.name,
                        :min_arrival,
                        t,
                        r.min_arrival,
                    )
                )
            end
            if r.law==:junction
                fix(a[i, t], 0.0; force = true)
            elseif river_ops[i]!==nothing
                r.law==:weir && fix(a[i, t], 1.0; force = true)
                level=H[ix[r.source], t]
                @constraint(m, first(r.discharge_curve.x)<=level<=last(r.discharge_curve.x))
                @constraint(m, (RQ[i, t]-a[i, t]*river_ops[i](level))/100==0)
            elseif r.law==:controlled
                @constraint(m, RQ[i, t]==cap*a[i, t])
            elseif r.law in (:orifice, :weir)
                r.law==:weir && fix(a[i, t], 1.0; force = true)
                hd=H[ix[r.source], t]-r.crest
                wet=r.allow_dry ? @expression(m, max(hd, 0.0)) : hd
                if r.law==:orifice
                    @constraint(m, (RQ[i, t]-r.coefficient*a[i, t]*sqrt(wet))/100==0)
                else
                    @constraint(m, (RQ[i, t]-r.coefficient*wet^1.5)/100==0)
                end
            else
                error("Unknown river law")
            end
        end
        if r.arrival_policy==:interval_average && !isempty(r.arrival_window_grid)
            windows=r.arrival_window_grid
            K=exact ? nothing : transfer_coefficients(r, windows, c.grid)
            hist=exact ? nothing :
                 route_volumes(r, windows, r.history_grid, r.history_release)
            for j in 1:(length(windows) - 1)
                if exact
                    ex=_pulse_expression(nd.arrival_pulses[i], windows[j], windows[j + 1])
                    volume=ex.history+sum(v*RQ[d, k] for (d, k, v) in ex.terms; init = 0.0)
                else
                    ks=[k for k in 1:T if K[j, k, 1]!=0 || K[j, k, 2]!=0]
                    volume=hist[j]+sum(
                        0.0036*dt[k]*RQ[i, k]*(
                            K[j, k, 1]+RQ[i, k]/r.capacity*(K[j, k, 2]-K[j, k, 1])
                        ) for k in ks;
                        init = 0.0,
                    )
                end
                @constraint(
                    m,
                    volume>=0.0036*(windows[j + 1]-windows[j])*opaverage(
                        c,
                        r.name,
                        :min_arrival,
                        windows[j],
                        windows[j + 1],
                        r.min_arrival,
                    )
                )
            end
        end
        if r.arrival_policy==:pointwise
            knots=exact ?
                  sort!(unique!(vcat(deterministic_knots(nd, i), [z for z in c.grid]))) :
                  operational_arrival_knots(c, r, c.grid)
            for time in knots
                sides=time==first(c.grid) ? (:right,) :
                      time==last(c.grid) ? (:left,) : (:left, :right)
                for side in sides
                    if exact
                        ex=deterministic_point_data(nd, i, time; side)
                        rate=ex.history+sum(
                            v*RQ[d, k] for (d, k, v) in ex.terms;
                            init = 0.0,
                        )
                    else
                        K=point_coefficients(r, c.grid, time; side = side)
                        history=point_arrival(r, c.grid, zeros(T), time; side = side)
                        ks=[k for k in 1:T if K[k, 1]!=0 || K[k, 2]!=0]
                        rate=history+sum(
                            RQ[i, k]*(K[k, 1]+RQ[i, k]/r.capacity*(K[k, 2]-K[k, 1])) for
                            k in ks;
                            init = 0.0,
                        )
                    end
                    @constraint(
                        m,
                        rate>=opvalue(c, r.name, :min_arrival, time, r.min_arrival; side)+arrival_margin
                    )
                end
            end
        end
        if exact
            push!(
                terminal,
                nd.terminal_history[i]+sum(
                    v*RQ[d, k] for (d, k, v) in nd.terminal_terms[i];
                    init = 0.0,
                ),
            )
            # Intermediate capacity constrains actual arrival-shaped releases, not only their averages.
            for time in
                sort!(unique!(vcat(deterministic_knots(nd, i; kind = :release), c.grid)))
                for side in (
                    time==first(c.grid) ? (:right,) :
                    time==last(c.grid) ? (:left,) : (:left, :right)
                )
                    ex=deterministic_point_data(nd, i, time; side, kind = :release)
                    rate=ex.history+sum(v*RQ[d, k] for (d, k, v) in ex.terms; init = 0.0)
                    @constraint(
                        m,
                        rate<=opvalue(c, r.name, :capacity, time, r.capacity; side)
                    )
                end
            end
        else
            push!(
                terminal,
                rd["history_terminal"][i]+sum(
                    0.0036*dt[k]*RQ[i, k]*(
                        1-sum(B[:, k, 1])-RQ[i, k]/r.capacity*sum(B[:, k, 2]-B[:, k, 1])
                    ) for k in 1:T
                ),
            )
        end
    end
    for (outgoing, incoming) in river_junction_incidence, t in 1:T
        @constraint(
            m,
            (0.0036*dt[t]*RQ[outgoing, t]-sum(arrivals[k, t] for k in incoming))/0.3==0
        )
    end
    for i in 1:(R + length(s.junctions)), t in 1:T
        net=sum(sign*Q[j, t] for (j, sign) in tunnel_incidence[i]; init = 0.0)
        net+=sum(sign*GQ[j, t] for (j, sign) in generator_incidence[i]; init = 0.0)
        net+=sum(
            (incoming ? arrivals[j, t]/(0.0036*dt[t]) : 0)-(outgoing ? RQ[j, t] : 0) for
            (j, incoming, outgoing) in river_incidence[i];
            init = 0.0,
        )
        if i<=R
            @constraint(
                m,
                (
                    V[i, t + 1]-V[i, t]-0.0036*dt[t]*(
                        net+opinterval(
                            c,
                            s.reservoirs[i].name,
                            :inflow,
                            t,
                            s.reservoirs[i].inflow,
                        )
                    )
                )/s.reservoirs[i].vmax==0
            )
        else
            @constraint(m, net/50==0)
        end
    end
    startup=sum(
        g.startup*(
            joint ? transitions.su[j, t] :
            max(0, u[j, t]-(t==1 ? g.initial_on : u[j, t - 1]))
        ) for (j, g) in enumerate(s.generators), t in 1:T
    )
    shutdown=sum(
        g.shutdown*(
            joint ? transitions.sd[j, t] :
            max(0, (t==1 ? g.initial_on : u[j, t - 1])-u[j, t])
        ) for (j, g) in enumerate(s.generators), t in 1:T
    )
    history_initial=exact ? nd.initial_transit : rd["history_initial"]
    obj=sum(c.prices[t]*dt[t]*P[j, t] for j in 1:G, t in 1:T)-startup-shutdown+sum(
        r.water_value*(V[i, T + 1]-r.v0) for (i, r) in enumerate(s.reservoirs)
    )+sum(
        r.water_value*(terminal[i]-history_initial[i]) for (i, r) in enumerate(s.rivers);
        init = 0.0,
    )
    penalty_cost=sum(
        opinterval(c, r.name, :release_penalty, t, 0.0)*shortfall_release[i, t] for
        (i, r) in enumerate(s.rivers), t in 1:T;
        init = 0.0,
    )
    obj-=penalty_cost
    @objective(m, Max, feasibility_only ? 0.0 : obj/10000)
    rng=MersenneTwister(seed)
    for i in 1:R, t in 1:T
        set_start_value(h[i, t], head(s.reservoirs[i], s.reservoirs[i].v0)/250)
    end
    for (i, g) in enumerate(s.generators), t in 1:T
        guess=clamp((seed==1 ? 0.6 : 0.45+0.3rand(rng))*g.qmax, g.qmin, g.qmax)*ustart[i, t]
        set_start_value(gq[i, t], guess/50)
        set_start_value(p[i, t], power(g, guess, g.hbest)/40)
    end
    for (i, r) in enumerate(s.rivers), t in 1:T
        set_start_value(rq[i, t], mean(r.history_release)/100)
        r.law==:orifice && set_start_value(a[i, t], 0.5)
    end
    if warm!==nothing && haskey(warm, "V")
        for i in 1:R, t in 1:(T + 1)
            set_start_value(v[i, t], warm["V"][i, t]/s.reservoirs[i].vmax)
        end
        for i in 1:N, t in 1:T
            set_start_value(h[i, t], warm["H"][i, t]/250)
        end
        for i in 1:E, t in 1:T
            set_start_value(q[i, t], warm["tunnel_q"][i, t]/50)
        end
        for (i, g) in enumerate(s.generators), t in 1:T
            set_start_value(
                gq[i, t],
                clamp(warm["generator_q"][i, t], g.qmin*ustart[i, t], g.qmax*ustart[i, t])/50,
            )
            set_start_value(
                p[i, t],
                clamp(warm["power"][i, t], g.pmin*ustart[i, t], g.pmax*ustart[i, t])/40,
            )
        end
        for i in 1:D, t in 1:T
            set_start_value(rq[i, t], warm["river_release"][i, t]/100)
            set_start_value(a[i, t], warm["gate"][i, t])
        end
    end
    (;
        m,
        obj,
        V,
        H,
        Q,
        GQ,
        P,
        RQ,
        a,
        arrivals,
        terminal,
        u,
        transitions,
        starttime,
        shortfall_release,
        penalty_cost,
        transport = nd,
    )
end

function _dispatch_values(b)
    Dict{String,Any}(
        "objective"=>value(b.obj),
        "V"=>value.(b.V),
        "H"=>value.(b.H),
        "tunnel_q"=>value.(b.Q),
        "generator_q"=>value.(b.GQ),
        "power"=>value.(b.P),
        "river_release"=>value.(b.RQ),
        "gate"=>value.(b.a),
        "arrival_volume"=>value.(b.arrivals),
        "terminal_transit"=>value.(b.terminal),
        "shortfall_release"=>value.(b.shortfall_release),
        "release_penalty_cost"=>value(b.penalty_cost),
    )
end

function _local_optimizer!(m, solver, time_limit)
    solver == "Ipopt" || throw(ArgumentError("Unsupported local solver $solver"))
    set_optimizer(m, Ipopt.Optimizer)
    set_silent(m)
    if solver=="Ipopt"
        set_optimizer_attribute(m, "tol", 1e-8)
        set_optimizer_attribute(m, "bound_relax_factor", 0.0)
        set_optimizer_attribute(m, "max_iter", 1000)
        set_optimizer_attribute(m, "max_cpu_time", max(0.001, time_limit))
        set_optimizer_attribute(m, "max_wall_time", max(0.001, time_limit))
    end
end

function solve_case(
    c;
    u = ones(Int, length(c.system.generators), length(c.prices)),
    solver = "Ipopt",
    seed = 1,
    warm = nothing,
    arrival_margin = 0.0,
    operational_margin = 0.0,
    time_limit = 45.0,
    feasibility_only = false,
    transport = nothing,
)
    solver=string(solver)
    solver == "Ipopt" || throw(ArgumentError("Unsupported local solver $solver"))
    b=_build_dispatch(
        c;
        u,
        solver,
        seed,
        warm,
        arrival_margin,
        operational_margin,
        time_limit,
        feasibility_only,
        transport,
    )
    construction_seconds=time()-b.starttime
    remaining=time_limit-construction_seconds
    result=Dict{String,Any}(
        "case"=>c.name,
        "solver"=>solver,
        "seed"=>seed,
        "status"=>"CONSTRUCTION_BUDGET_EXHAUSTED",
        "seconds"=>0.0,
        "feasibility_only"=>feasibility_only,
        "u"=>u,
        "variable_count"=>num_variables(b.m),
        "construction_seconds"=>construction_seconds,
        "arrival_margin"=>arrival_margin,
        "operational_margin"=>operational_margin,
    )
    if remaining<=0
        result["total_seconds"]=time()-b.starttime
        return result
    end
    _local_optimizer!(b.m, solver, remaining)
    construction_seconds=time()-b.starttime
    result["construction_seconds"]=construction_seconds
    if time_limit-construction_seconds<=0
        result["total_seconds"]=time()-b.starttime
        return result
    end
    elapsed=@elapsed optimize!(b.m)
    result["status"]=string(termination_status(b.m))
    result["seconds"]=elapsed
    if has_values(b.m)
        merge!(result, _dispatch_values(b))
        audit_start=time()
        result["validation"]=validate(c, result; transport = b.transport)
        result["validation_seconds"]=time()-audit_start
        if has_operational_data(c) && !result["validation"]["valid"]
            try
                repaired=dispatch_from_controls(
                    c,
                    u,
                    result["generator_q"],
                    result["gate"];
                    transport = b.transport,
                )
                result["forward_reconstruction_audit"]=repaired["validation"]
                if repaired["validation"]["valid"]
                    result["raw_solver_validation"]=result["validation"]
                    result["raw_solver_objective"]=result["objective"]
                    result["forward_reconstructed"]=true
                    for key in (
                        "V",
                        "H",
                        "tunnel_q",
                        "power",
                        "river_release",
                        "arrival_volume",
                        "terminal_transit",
                        "objective",
                        "shortfall_release",
                        "release_penalty_cost",
                        "validation",
                    )
                        result[key]=repaired[key]
                    end
                end
            catch e
                result["forward_reconstruction_error"]=sprint(showerror, e)
            end
        end
    end
    result["total_seconds"]=time()-b.starttime
    result
end
