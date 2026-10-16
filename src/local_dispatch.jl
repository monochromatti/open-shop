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
    tunnel_opening = nothing,
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
    topology = HydraulicTopology(s)
    node_names = topology.node_names
    T=length(c.prices)
    R=length(s.reservoirs)
    N=length(node_names)
    E=length(s.tunnels)
    G=length(s.generators)
    D=length(s.rivers)
    ix = topology.node_index
    dt=diff(c.grid)
    generator_plants = s.plants[topology.generator_plant]
    plant_generators = Dict(p.name => topology.plant_generators[i] for (i, p) in enumerate(s.plants))
    warm, warm_start=_dispatch_start(c, u, warm)
    m=Model()
    set_silent(m)
    openings=tunnel_opening===nothing ? tunnel_opening_schedule(c) : tunnel_opening
    gates=_tunnel_gate_states!(m,c;fixed=openings)
    constrain_tunnel_switching!(m,c,gates.open,u)
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
    mean_level_ops=[r.level_curve===nothing ? nothing :
        _mean_head_operator(m,r;name=Symbol("mean_level_",i)) for
        (i,r) in enumerate(s.reservoirs)]
    turbine_ops=[
        g.turbine_table===nothing ? nothing :
        turbine_operator(m, g.turbine_table; name = Symbol("turbine_", i)) for
        (i, g) in enumerate(s.generators)
    ]
    characteristic_ops=[g.characteristic===nothing ? nothing : (
        JuMP.add_nonlinear_operator(m,1,h->pump_discharge(g.characteristic,h;extrapolation=:linear);
            name=Symbol("pump_flow_",i)),
        JuMP.add_nonlinear_operator(m,1,h->pump_efficiency(g.characteristic,h;extrapolation=:linear);
            name=Symbol("pump_eta_",i))) for (i,g) in enumerate(s.generators)]
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
    channel_ops=[is_channel(r) ? level_flow_operator(m,channel_curve(r);name=Symbol("channel_",i)) : nothing
        for (i,r) in enumerate(s.rivers)]
    turbine_flow_ops=[
        if g.turbine_table!==nothing && any(t->u[i, t]!=0, 1:T)
            table=g.turbine_table
            (
                turbine_flow_operator(
                    m,
                    table, table.qmin;
                    name = Symbol("qlo_", i),
                ),
                turbine_flow_operator(
                    m,
                    table, table.qmax;
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
            lo, hi=soft_storage_bounds(c, r, c.grid[t])
            set_lower_bound(v[i, t], lo/r.vmax)
            set_upper_bound(v[i, t], hi/r.vmax)
            set_start_value(v[i, t], r.v0/r.vmax)
        end
        fix(v[i, 1], r.v0/r.vmax; force = true)
        for t in 1:T
            set_lower_bound(h[i, t], head(r, r.vmin)/250)
            set_upper_bound(h[i, t], head(r, r.vmax)/250)
            level=mean_level_ops[i]===nothing ? mean_head(r,V[i,t],V[i,t+1]) :
                  mean_level_ops[i](V[i,t],V[i,t+1])
            @constraint(m, (H[i, t]-level)/250==0)
        end
    end
    level_at_vertex=(i,t)->begin
        r=s.reservoirs[i]
        volume=V[i,t]
        level_ops[i]===nothing ? r.z0+r.slope*volume+r.curvature*volume^2 : level_ops[i](volume)
    end
    constrain_reservoir_ramps!(m,c,V; level=level_at_vertex)
    soft=constrain_soft_constraints!(m,c,V,P,GQ;level=level_at_vertex)
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
        opening=gates.opening[i,t]
        factor, resistance=tunnel_coefficients(e, opening)
        set_lower_bound(q[i, t], (e.discharge_river===nothing ? -cap : 0.0)/50)
        set_upper_bound(q[i, t], cap/50)
        if opening==0
            fix(q[i, t], 0; force = true)
        else
            @constraint(
                m,
                (
                    factor*(H[ix[e.source], t]-H[ix[e.target], t])-resistance*Q[i, t]*abs(
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
        hd=unit_power_sign(g)*(H[ix[plant.source], t]-receiver-tail)
        eta=characteristic_ops[i]!==nothing ? characteristic_ops[i][2](hd) : turbine_ops[i]===nothing ?
            g.efficiency-g.qcurvature*((GQ[i, t]-g.qbest)/g.qbest)^2-g.hcurvature*(
            (hd-g.hbest)/g.hbest
        )^2 : turbine_ops[i](GQ[i, t], hd)
        electrical=electrical_ops[i]===nothing ? 1.0 : electrical_ops[i](P[i, t])
        qmin=opinterval(c, g.name, :qmin, t, g.qmin)
        qmax=opinterval(c, g.name, :qmax, t, g.qmax)
        pmin=opinterval(c, g.name, :pmin, t, g.pmin)
        pmax=opinterval(c, g.name, :pmax, t, g.pmax)
        conversion=g.mode==:pump ?
            P[i,t]*eta*g.motor_efficiency-0.00981*GQ[i,t]*hd :
            P[i,t]-0.00981*GQ[i,t]*hd*eta*electrical
        begin
            set_lower_bound(gq[i, t], qmin*u[i, t]/50)
            set_upper_bound(gq[i, t], qmax*u[i, t]/50)
            set_lower_bound(p[i, t], (pmin+operational_margin)*u[i, t]/40)
            set_upper_bound(p[i, t], max(0.0, pmax-operational_margin)*u[i, t]/40)
            if u[i, t]==0
                fix(gq[i, t], 0; force = true)
                fix(p[i, t], 0; force = true)
            else
                if characteristic_ops[i]!==nothing
                    @constraint(m,(GQ[i,t]-characteristic_ops[i][1](hd))/50==0)
                end
                @constraint(m, conversion/40==0)
                @constraint(m, g.hmin<=hd<=g.hmax)
                @constraint(m, g.min_efficiency<=eta<=1.0)
                if g.turbine_table!==nothing
                    qlo, qhi=turbine_flow_ops[i]
                    flow_margin=operational_margin/unit_reference_power_slope(g)
                    @constraint(m, GQ[i, t]>=qlo(hd)+flow_margin)
                    @constraint(m, GQ[i, t]<=qhi(hd)-flow_margin)
                end
            end
        end
    end
    reserves=_local_reserve_points!(m,c,u,H,P,GQ,turbine_ops,electrical_ops,tailwater_ops,outlet_ops,turbine_flow_ops)
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
            set_lower_bound(rq[i,t],river_flow_lower(r,cap,opinterval(c,r.name,:inflow,t,r.inflow))/100)
            set_lower_bound(a[i, t], outlet_bounds.gate_lower[i, t])
            set_upper_bound(a[i, t], outlet_bounds.gate_upper[i, t])
            if is_channel(r)
                fix(a[i,t],0.0;force=true)
                fix(shortfall_release[i,t],0.0;force=true)
                upstream,downstream=H[ix[r.source],t],H[ix[r.target],t]
                flow=channel_ops[i](upstream-downstream,channel_reference(r,upstream,downstream))
                @constraint(m,(RQ[i,t]-flow)/100==0)
                continue
            end
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
                source=only(x for x in s.reservoirs if x.name==r.source)
                upper=sqrt(max(0.0,head(source,source.vmax)-r.crest))
                root=@variable(m,lower_bound=0.0,upper_bound=upper,
                    base_name="outlet_root_$(i)_$(t)")
                initial=warm===nothing ? head(source,source.v0) : warm.head[ix[r.source],t]
                set_start_value(root,sqrt(max(0.0,initial-r.crest)))
                # The squared-root equality avoids singular sqrt/power
                # derivatives without redundant complementarity inequalities.
                if upper==0
                    fix(root,0.0;force=true)
                else
                    wet=r.allow_dry ? @expression(m,max(hd,0.0)) : hd
                    @constraint(m,root^2==wet)
                end
                if r.law==:orifice
                    @constraint(m, (outlet_release-r.coefficient*a[i, t]*root)/100==0)
                else
                    @constraint(m, (outlet_release-r.coefficient*root^3)/100==0)
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
    constrain_water_balance!(m, c, topology, V, Q, GQ, RQ, arrivals)
    history_initial=exact ? nd.initial_transit : rd["history_initial"]
    obj,penalty_cost=dispatch_objective(m,c,P,V,terminal,u,shortfall_release;
        initial_transit=history_initial, gate_cost=gates.cost, soft_cost=soft.cost)
    reserves===nothing || (obj+=reserves.value)
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
            set_start_value(v[i, t], warm.storage[i, t]/s.reservoirs[i].vmax)
        end
        for i in 1:N, t in 1:T
            set_start_value(h[i, t], warm.head[i, t]/250)
        end
        for i in 1:E, t in 1:T
            set_start_value(q[i, t], warm.tunnel_q[i, t]/50)
        end
        for (i, g) in enumerate(s.generators), t in 1:T
            set_start_value(
                gq[i, t],
                clamp(warm.generator_q[i, t], g.qmin*ustart[i, t], g.qmax*ustart[i, t])/50,
            )
            set_start_value(
                p[i, t],
                clamp(warm.power[i, t], g.pmin*ustart[i, t], g.pmax*ustart[i, t])/40,
            )
        end
        for i in 1:D, t in 1:T
            set_start_value(rq[i, t], warm.river_release[i, t]/100)
            set_start_value(a[i, t], warm.gate[i, t])
        end
    end
    _start_reserves!(reserves,c,warm)
    (;
        reserves,
        case=c,
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
        soft,
        transport = nd,
        warm_start,
        gates,
    )
end

function _dispatch_values(b)
    c=b.case
    release=value.(b.RQ)
    opening=_tunnel_gate_solution(b.gates)
    raw=Schedule(copy(c.grid), round.(Int,value.(b.u)), value.(b.V), value.(b.H),
        value.(b.Q), Matrix{Float64}(opening), value.(b.GQ), value.(b.P), value.(b.a),
        release, Matrix{Float64}(value.(b.arrivals)), Vector{Float64}(value.(b.terminal)), value.(b.shortfall_release),
        value(b.obj), value(b.gates.cost), value(b.penalty_cost),
        Dict{String,Any}("case"=>c.name,
            "fractional_tunnel_opening"=>value.(b.gates.opening));reserves=_reserve_values(c,b.reserves))
    with_soft_assessment(c,raw;included_cost=value(b.soft.cost))
end

function _repair_dispatch!(c, solution::Schedule, report; transport=nothing)
    began=time()
    try
        repaired, correction=_reconstruct_candidate(c,solution;transport)
        report["forward_reconstruction_audit"]=repaired.diagnostics["validation"]
        report["forward_reconstruction_correction"]=correction
        if repaired.diagnostics["validation"]["valid"]
            report["raw_solver_validation"]=solution.diagnostics["validation"]
            report["raw_solver_objective"]=solution.objective
            report["forward_reconstructed"]=true
            solution=repaired
        end
    catch error
        report["forward_reconstruction_error"]=sprint(showerror,error)
    end
    report["forward_reconstruction_seconds"]=time()-began
    solution
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

"""
    solve_case(case; u=ones(Int, length(case.system.generators), length(case.prices)),
               tunnel_opening=nothing, time_limit=45.0, warm=nothing, ...)

Optimize continuous dispatch with Ipopt for fixed commitment `u[unit, interval]`
and fixed tunnel gate positions. Gate positions default to the case's prescribed
schedule; this routine never searches discrete states. Use [`solve`](@ref) for
the global mixed-integer problem, or [`schedule_case`](@ref) to prepare a feasible
schedule by combining MILP proposals with nonlinear dispatch.

Returns a [`SolveResult`](@ref). Its `solution` is a [`Schedule`](@ref), or
`nothing` if no candidate exists. `accepted` means the candidate passed the
original-grid physical and operating checks; use [`solve_verified`](@ref) for
finer replay. A local optimum is not a global certificate.
"""
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
    tunnel_opening = nothing,
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
        tunnel_opening,
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
        return SolveResult(result)
    end
    _local_optimizer!(b.m, solver, remaining)
    construction_seconds=time()-b.starttime
    result["construction_seconds"]=construction_seconds
    if time_limit-construction_seconds<=0
        result["total_seconds"]=time()-b.starttime
        return SolveResult(result)
    end
    elapsed=@elapsed optimize!(b.m)
    result["status"]=string(termination_status(b.m))
    result["seconds"]=elapsed
    result["iterations"]=MOI.get(b.m, MOI.BarrierIterations())
    solution=nothing
    if has_values(b.m)
        solution=_dispatch_values(b)
        audit_start=time()
        solution.diagnostics["validation"]=validate(c,solution; transport=b.transport)
        result["validation_seconds"]=time()-audit_start
        if !solution.diagnostics["validation"]["valid"]
            solution=_repair_dispatch!(c,solution,result; transport=b.transport)
        end
        result["accepted"]=solution.diagnostics["validation"]["valid"]
    end
    result["total_seconds"]=time()-b.starttime
    SolveResult(result;solution)
end
