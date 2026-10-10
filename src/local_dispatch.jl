admissible(c, u) = _commitment_admissible(c, u)
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
    transport = nothing,
)
    starttime=time()
    isfinite(arrival_margin) && arrival_margin>=0 ||
        throw(ArgumentError("invalid arrival margin"))
    isfinite(operational_margin) && operational_margin>=0 ||
        throw(ArgumentError("invalid operational margin"))
    isfinite(time_limit) && time_limit>0 ||
        throw(ArgumentError("time limit must be finite and positive"))
    validate_inputs(c)
    admissible(c, u) || error("Invalid commitment")
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
            (j, water_incidence(e,name)) for
            (j, e) in enumerate(s.tunnels) if e.target==name || e.source==name
        ] for name in node_names
    ]
    generator_incidence=[
        [
            (j, water_incidence(plant,name)) for
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
    m=Model()
    set_silent(m)
    ustart=copy(u)
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
        if g.turbine_table!==nothing && any(t->u[i, t]!=0, 1:T)
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
    level_at_vertex=(i,t)->begin
        r=s.reservoirs[i]
        volume=V[i,t]
        level_ops[i]===nothing ? r.z0+r.slope*volume+r.curvature*volume^2 : level_ops[i](volume)
    end
    constrain_reservoir_ramps!(m,c,V; level=level_at_vertex)
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
        set_lower_bound(q[i, t], (e.discharge_river===nothing ? -cap : 0.0)/50)
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
    constrain_dispatch_operations!(m,c,u,P,GQ,RQ; margin=operational_margin)
    outlet_bounds = _dispatch_outlet_bounds(c)
    injections=river_injections(c,RQ,GQ,Q)
    routing = _dispatch_transport_expressions(m,c, RQ; transport, injections)
    (; exact, nd, rd, arrivals, terminal) = routing
    transfer=(i,k,b)->river_transfer_expression(m,s.rivers[i],RQ[i,k],b)
    for (i, r) in enumerate(s.rivers)
        for t in 1:T
            cap=outlet_bounds.capacity[i, t]
            set_upper_bound(rq[i, t], cap/100)
            set_lower_bound(rq[i,t],opinterval(c,r.name,:inflow,t,r.inflow)/100)
            set_lower_bound(a[i, t], outlet_bounds.gate_lower[i, t])
            set_upper_bound(a[i, t], outlet_bounds.gate_upper[i, t])
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
            natural=opinterval(c,r.name,:inflow,t,r.inflow)
            outlet_release=RQ[i,t]-natural
            if r.law==:junction
                fix(a[i, t], 0.0; force = true)
            elseif river_ops[i]!==nothing
                r.law==:weir && fix(a[i, t], 1.0; force = true)
                level=H[ix[r.source], t]
                @constraint(m, first(r.discharge_curve.x)<=level<=last(r.discharge_curve.x))
                @constraint(m, (outlet_release-a[i, t]*river_ops[i](level))/100==0)
            elseif r.law==:controlled
                @constraint(m, outlet_release==(cap-natural)*a[i, t])
            elseif r.law in (:orifice, :weir)
                r.law==:weir && fix(a[i, t], 1.0; force = true)
                hd=H[ix[r.source], t]-r.crest
                wet=r.allow_dry ? @expression(m, max(hd, 0.0)) : hd
                if r.law==:orifice
                    @constraint(m, (outlet_release-r.coefficient*a[i, t]*sqrt(wet))/100==0)
                else
                    @constraint(m, (outlet_release-r.coefficient*wet^1.5)/100==0)
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
                    volume=ex.history+sum(v*injections[d, k] for (d, k, v) in ex.terms; init = 0.0)
                else
                    ks=[k for k in 1:T if any(!iszero,view(K,j,k,:))]
                    volume=hist[j]+sum(
                        0.0036*dt[k]*transfer(i,k,collect(view(K,j,k,1:river_curve_count(r)))) for k in ks;
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
                            v*injections[d, k] for (d, k, v) in ex.terms;
                            init = 0.0,
                        )
                    else
                        K=point_coefficients(r, c.grid, time; side = side)
                        history=point_arrival(r, c.grid, zeros(T), time; side = side)
                        ks=[k for k in 1:T if any(!iszero,view(K,k,:))]
                        rate=history+sum(
                            transfer(i,k,collect(view(K,k,1:river_curve_count(r)))) for
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
            # Intermediate capacity constrains actual arrival-shaped releases, not only their averages.
            for time in
                sort!(unique!(vcat(deterministic_knots(nd, i; kind = :release), c.grid)))
                for side in (
                    time==first(c.grid) ? (:right,) :
                    time==last(c.grid) ? (:left,) : (:left, :right)
                )
                    ex=deterministic_point_data(nd, i, time; side, kind = :release)
                    rate=ex.history+sum(v*injections[d, k] for (d, k, v) in ex.terms; init = 0.0)
                    @constraint(
                        m,
                        rate<=opvalue(c, r.name, :capacity, time, r.capacity; side)
                    )
                end
            end
        end
    end
    constrain_river_sources!(m,c,RQ,arrivals,injections)
    for i in 1:(R + length(s.junctions)), t in 1:T
        net=sum(sign*Q[j, t] for (j, sign) in tunnel_incidence[i]; init = 0.0)
        net+=sum(sign*GQ[j, t] for (j, sign) in generator_incidence[i]; init = 0.0)
        net+=sum(
            (incoming ? arrivals[j, t]/(0.0036*dt[t]) : 0)-(outgoing ? RQ[j, t]-opinterval(c,s.rivers[j].name,:inflow,t,s.rivers[j].inflow) : 0) for
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
    operating_cost=transition_costs(c,u)
    history_initial=exact ? nd.initial_transit : rd["history_initial"]
    obj=sum(c.prices[t]*dt[t]*P[j, t] for j in 1:G, t in 1:T;init=0.0)-operating_cost+sum(
        r.water_value*(V[i, T + 1]-r.v0) for (i, r) in enumerate(s.reservoirs);
        init = 0.0,
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
    warm, warm_start=_dispatch_start(c, u, warm)
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
        starttime,
        shortfall_release,
        penalty_cost,
        transport = nd,
        warm_start,
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

function _repair_dispatch!(c, result; transport = nothing)
    began=time()
    try
        repaired, correction=_reconstruct_candidate(c, result; transport)
        result["forward_reconstruction_audit"]=repaired["validation"]
        result["forward_reconstruction_correction"]=correction
        if repaired["validation"]["valid"]
            result["raw_solver_validation"]=result["validation"]
            result["raw_solver_objective"]=result["objective"]
            result["forward_reconstructed"]=true
            for key in (
                "generator_q", "gate", "V", "H", "tunnel_q", "power",
                "river_release", "arrival_volume", "terminal_transit", "objective",
                "shortfall_release", "release_penalty_cost", "validation",
            )
                result[key]=repaired[key]
            end
        end
    catch error
        result["forward_reconstruction_error"]=sprint(showerror, error)
    end
    result["forward_reconstruction_seconds"]=time()-began
    result
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
        "grid"=>copy(c.grid),
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
        "warm_start"=>b.warm_start,
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
    result["iterations"]=MOI.get(b.m, MOI.BarrierIterations())
    if has_values(b.m)
        merge!(result, _dispatch_values(b))
        audit_start=time()
        result["validation"]=validate(c, result; transport = b.transport)
        result["validation_seconds"]=time()-audit_start
        result["validation"]["valid"] || _repair_dispatch!(c, result; transport = b.transport)
    end
    result["total_seconds"]=time()-b.starttime
    result
end
