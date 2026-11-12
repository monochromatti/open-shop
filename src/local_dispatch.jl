admissible(c, u) = _commitment_admissible(c, u)

"""Intersect a variable domain with its physical bounds, including shared cells."""
function intersect_bounds!(variable, lower, upper)
    if is_fixed(variable)
        lower<=fix_value(variable)<=upper || throw(ArgumentError("incompatible variable bounds"))
        return
    end
    lo=has_lower_bound(variable) ? max(lower,lower_bound(variable)) : lower
    hi=has_upper_bound(variable) ? min(upper,upper_bound(variable)) : upper
    lo<=hi || throw(ArgumentError("incompatible variable bounds"))
    if lo==hi
        fix(variable,lo;force=true)
    else
        set_lower_bound(variable,lo)
        set_upper_bound(variable,hi)
    end
    nothing
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
    transport = nothing,
    tunnel_opening = nothing,
    river_on = nothing,
    model = Model(),
    prefix = "",
    shared_controls = nothing,
    control_index = nothing,
    integral_grid = c.grid,
    grid = c.grid,
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
    operating_case=c
    mesh=TimeDiscretization(c.grid,grid)
    operating_u=u
    u=u[:,mesh.interval]
    c=hydraulic_case(c,mesh.grid)
    s=c.system
    topology = HydraulicTopology(s)
    T=length(c.prices)
    R=length(s.reservoirs)
    D=length(s.rivers)
    ix = topology.node_index
    dt=diff(c.grid)
    generator_plants = s.plants[topology.generator_plant]
    plant_generators = Dict(p.name => topology.plant_generators[i] for (i, p) in enumerate(s.plants))
    warm,warm_start=_dispatch_start(operating_case,operating_u,warm;grid=mesh.grid)
    m=model
    set_silent(m)
    openings=tunnel_opening===nothing ? tunnel_opening_schedule(operating_case) : tunnel_opening
    gates=_tunnel_gate_states!(m,operating_case;fixed=openings)
    constrain_tunnel_switching!(m,operating_case,gates.open,operating_u)
    ustart=copy(u)
    variables=DispatchVariables(m,s,mesh;prefix,shared_controls,
        control_index=control_index===nothing ? mesh.interval : control_index)
    (;v,h,q,gq,p,rq,a,V,H,Q,GQ,P,RQ)=variables
    shortfall_release=@variable(m,[1:D,1:T],lower_bound=0,base_name=prefix*"shortfall_release")
    level_ops=[
        r.level_curve===nothing ? nothing :
        table_operator(m, r.level_curve; name = Symbol(prefix,"level_", i)) for
        (i, r) in enumerate(s.reservoirs)
    ]
    mean_level_ops=[r.level_curve===nothing ? nothing :
        _mean_head_operator(m,r;name=Symbol(prefix,"mean_level_",i)) for
        (i,r) in enumerate(s.reservoirs)]
    turbine_ops=[
        g.turbine_table===nothing ? nothing :
        turbine_operator(m, g.turbine_table; name = Symbol(prefix,"turbine_", i)) for
        (i, g) in enumerate(s.generators)
    ]
    characteristic_ops=[g.characteristic===nothing ? nothing : (
        JuMP.add_nonlinear_operator(m,1,h->pump_discharge(g.characteristic,h;extrapolation=:linear);
            name=Symbol(prefix,"pump_flow_",i)),
        JuMP.add_nonlinear_operator(m,1,h->pump_efficiency(g.characteristic,h;extrapolation=:linear);
            name=Symbol(prefix,"pump_eta_",i))) for (i,g) in enumerate(s.generators)]
    electrical_ops=[
        g.generator_efficiency_curve===nothing ? nothing :
        table_operator(m, g.generator_efficiency_curve; name = Symbol(prefix,"electrical_", i)) for
        (i, g) in enumerate(s.generators)
    ]
    tailwater_ops=Dict(
        p.name=>(
            p.tailwater_curve===nothing ? nothing :
            table_operator(m, p.tailwater_curve; name = Symbol(prefix,"tailwater_", p.name))
        ) for p in s.plants
    )
    river_ops=[
        r.discharge_curve===nothing ? nothing :
        table_operator(m, r.discharge_curve; name = Symbol(prefix,"river_law_", i)) for
        (i, r) in enumerate(s.rivers)
    ]
    channel_ops=[is_channel(r) ? level_flow_operator(m,channel_curve(r);name=Symbol(prefix,"channel_",i)) : nothing
        for (i,r) in enumerate(s.rivers)]
    turbine_flow_ops=[
        if g.turbine_table!==nothing && any(t->u[i, t]!=0, 1:T)
            table=g.turbine_table
            (
                turbine_flow_operator(
                    m,
                    table, table.qmin;
                    name = Symbol(prefix,"qlo_", i),
                ),
                turbine_flow_operator(
                    m,
                    table, table.qmax;
                    name = Symbol(prefix,"qhi_", i),
                ),
            )
        else
            nothing
        end for (i, g) in enumerate(s.generators)
    ]
    for (i, r) in enumerate(s.reservoirs)
        for t in 1:(T + 1)
            lo, hi=soft_storage_bounds(c, r, c.grid[t])
            set_lower_bound(v[i, t], lo/r.vmax)
            set_upper_bound(v[i, t], hi/r.vmax)
            set_start_value(v[i, t], r.v0/r.vmax)
        end
        fix(v[i, 1], r.v0/r.vmax; force = true)
        for t in 1:T
            set_lower_bound(h[i, t], head(r, r.vmin)/DISPATCH_HEAD_SCALE)
            set_upper_bound(h[i, t], head(r, r.vmax)/DISPATCH_HEAD_SCALE)
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
    operating_V=V[:,mesh.vertices]
    operating_level=(i,t)->level_at_vertex(i,mesh.vertices[t])
    operating_head=(i,t)->begin
        r=s.reservoirs[i]
        op=mean_level_ops[i]
        op===nothing ? mean_head(r,operating_V[i,t],operating_V[i,t+1]) :
            op(operating_V[i,t],operating_V[i,t+1])
    end
    constrain_reservoir_ramps!(m,operating_case,operating_V;level=operating_level)
    level_windows=constrain_level_windows!(m,c,level_at_vertex)
    constrain_gross_heads!(m,operating_case,operating_u,operating_head,operating_level)
    operating_P=interval_means(mesh,P)
    operating_GQ=interval_means(mesh,GQ)
    soft=constrain_soft_constraints!(m,operating_case,operating_V,operating_P,operating_GQ;
        level=operating_level)
    for (i, j) in enumerate(s.junctions), t in 1:T
        set_lower_bound(h[R + i, t], j.hmin/DISPATCH_HEAD_SCALE)
        set_upper_bound(h[R + i, t], j.hmax/DISPATCH_HEAD_SCALE)
        set_start_value(h[R + i, t], (j.hmin+j.hmax)/(2DISPATCH_HEAD_SCALE))
    end
    for (i, b) in enumerate(s.boundaries), t in 1:T
        fix(h[R + length(s.junctions) + i, t], b.head/DISPATCH_HEAD_SCALE; force = true)
    end
    for (i, e) in enumerate(s.tunnels), t in 1:T
        cap=opinterval(c, e.name, :capacity, t, e.capacity)
        opening=gates.opening[i,mesh.interval[t]]
        factor, resistance=tunnel_coefficients(e, opening)
        set_lower_bound(q[i, t], (e.discharge_river===nothing ? -cap : 0.0)/DISPATCH_FLOW_SCALE)
        set_upper_bound(q[i, t], cap/DISPATCH_FLOW_SCALE)
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
            lo=DISPATCH_HEAD_SCALE*(is_fixed(receiver) ? fix_value(receiver) : lower_bound(receiver))
            hi=DISPATCH_HEAD_SCALE*(is_fixed(receiver) ? fix_value(receiver) : upper_bound(receiver))
            floor===nothing || lo>=floor ? nothing :
            hi<=floor ? floor :
            table_operator(
                m,
                TableCurve([lo, floor, hi], [floor, floor, hi]);
                name = Symbol(prefix,"outlet_head_", plant.name),
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
            intersect_bounds!(gq[i,t],qmin*u[i,t]/DISPATCH_FLOW_SCALE,qmax*u[i,t]/DISPATCH_FLOW_SCALE)
            plo,phi=inset_bounds(pmin,pmax,operational_margin)
            set_lower_bound(p[i, t], plo*u[i, t]/DISPATCH_POWER_SCALE)
            set_upper_bound(p[i, t], phi*u[i, t]/DISPATCH_POWER_SCALE)
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
                    @constraint(m, GQ[i, t]>=qlo(hd))
                    @constraint(m, GQ[i, t]<=qhi(hd))
                end
            end
        end
    end
    reserves=_local_reserve_points!(m,c,u,H,P,GQ,turbine_ops,electrical_ops,tailwater_ops,outlet_ops,turbine_flow_ops;operating_case)
    outlet_bounds=(;
        capacity=[opinterval(c,r.name,:capacity,t,r.capacity) for r in s.rivers,t in 1:T],
        gate_lower=[opinterval(c,r.name,:gate_min,t,r.gate_min) for r in s.rivers,t in 1:T],
        gate_upper=[opinterval(c,r.name,:gate_max,t,1.0) for r in s.rivers,t in 1:T])
    injections=river_injections(c,RQ,GQ,Q)
    routing = RiverTransport(m,c, RQ; transport, injections,integral_grid)
    (;arrivals,terminal)=routing
    for (i, r) in enumerate(s.rivers)
        for t in 1:T
            cap=outlet_bounds.capacity[i, t]
            set_upper_bound(rq[i, t], cap/DISPATCH_RELEASE_SCALE)
            set_lower_bound(rq[i,t],river_flow_lower(r,cap,opinterval(c,r.name,:inflow,t,r.inflow))/DISPATCH_RELEASE_SCALE)
            intersect_bounds!(a[i,t],outlet_bounds.gate_lower[i,t],outlet_bounds.gate_upper[i,t])
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
            for j in 1:(length(windows)-1)
                volume=routing_volume_expression(m,c,routing,i,windows[j],windows[j+1])
                @constraint(m,volume>=0.0036*(windows[j+1]-windows[j])*opaverage(
                    c,r.name,:min_arrival,windows[j],windows[j+1],r.min_arrival))
            end
        end
        if r.arrival_policy==:pointwise
            for time in routing_arrival_knots(c,routing,i)
                sides=time==first(c.grid) ? (:right,) :
                    time==last(c.grid) ? (:left,) : (:left,:right)
                for side in sides
                    rate=routing_point_expression(m,c,routing,i,time;side)
                    @constraint(m,rate>=opvalue(c,r.name,:min_arrival,time,r.min_arrival;side)+arrival_margin)
                end
            end
        end
    end
    constrain_river_sources!(m,c,RQ,arrivals,injections)
    windows=operating_window_values(mesh,(;power=P,generator_q=GQ,river_release=RQ,
        storage=V,tunnel_q=Q,arrival_volume=arrivals),operating_u)
    constrain_plant_capacities!(m,c,u,P,GQ;margin=operational_margin)
    constrain_dispatch_operations!(m,operating_case,operating_u,windows.power,
        windows.generator_q,windows.river_release;margin=operational_margin,capacities=false)
    flow_integrals=constrain_flow_requirements!(m,operating_case,windows.generator_q,
        windows.river_release,windows.arrival_volume,windows.tunnel_q)
    operating_warm=warm===nothing ? nothing : operating_window_values(operating_case,warm)
    river_operations=constrain_river_operations!(m,operating_case,windows.river_release;
        fixed=fixed_river_on(operating_case,river_on,operating_warm))
    constrain_water_balance!(m, c, topology, V, Q, GQ, RQ, arrivals)
    constrain_storage_paths!(m, c, topology, V, Q, GQ, RQ, routing)
    history_initial=routing.initial_transit
    obj,penalty_cost=dispatch_objective(m,c,P,V,terminal,u,shortfall_release;
        initial_transit=history_initial,gate_cost=gates.cost,soft_cost=soft.cost,river_cost=river_operations.cost,
        operating_case,operating_u)
    reserves===nothing || (obj+=reserves.value)
    @objective(m, Max, feasibility_only ? 0.0 : obj/10000)
    initialize_dispatch!(variables,c,ustart;warm,seed)
    if warm!==nothing
        for (variable,value) in river_operation_start_values(river_operations,operating_case,operating_warm.river_release)
            set_start_value(variable,value)
        end
    end
    _start_reserves!(reserves,c,warm)
    start_transport!(routing)
    if warm!==nothing
        for (variable,value) in window_start_values(operating_case,operating_warm,flow_integrals,[])
            set_start_value(variable,value)
        end
        for (variable,value) in window_start_values(c,warm,[],level_windows)
            set_start_value(variable,value)
        end
    end
    (;
        variables,reserves,
        case=c,
        operating_case,mesh,operating_u,
        m,
        obj,
        V,
        H,
        Q,
        GQ,
        P,
        RQ,
        u,
        starttime,
        shortfall_release,
        penalty_cost,
        soft,
        transport = routing.compiled,
        routing,
        warm_start,
        gates,
        river_operations,
        flow_integrals,level_windows,
    )
end

# A second hydraulic mesh shares the requested controls. Its state variables are
# witnesses of the finer equations, not additional scheduling decisions. The
# objective and returned schedule still belong to the original model.
function constrain_replay!(b,grid;warm=nothing,river_on=nothing)
    objective=objective_function(b.m)
    parent=TimeDiscretization(b.case.grid,grid)
    refined=replay_case(b.case,grid)
    operating=_river_replace(b.operating_case;routing_grid=refined.routing_grid)
    f=_build_dispatch(operating;grid,u=b.operating_u,model=b.m,prefix="replay_",
        warm,tunnel_opening=Matrix{Float64}(b.gates.opening),river_on,
        shared_controls=b.variables,control_index=parent.interval)
    if b.reserves!==nothing
        base_cells=first.(b.mesh.cells)
        fine_cells=first.(f.mesh.cells)
        for (fine,base) in ((f.reserves.up,b.reserves.up),(f.reserves.down,b.reserves.down))
            for i in axes(base,1),k in axes(base,2),j in eachindex(base_cells)
                @constraint(b.m,fine[i,k,fine_cells[j]]==base[i,k,base_cells[j]])
            end
        end
    end
    set_objective_function(b.m,objective)
    nothing
end

function _repair_dispatch!(c, solution::Schedule, report; transport=nothing,deadline=Inf)
    began=time()
    try
        repaired, correction=_reconstruct_candidate(c,solution;transport,deadline)
        report["forward_reconstruction_audit"]=repaired.diagnostics["validation"]
        report["forward_reconstruction_correction"]=correction
        if repaired.diagnostics["validation"]["valid"]
            report["raw_solver_validation"]=solution.diagnostics["validation"]
            report["raw_solver_objective"]=solution.objective
            report["forward_reconstructed"]=true
            solution=repaired
        end
    catch error
        error isa Union{ArgumentError,DomainError} || rethrow()
        report["forward_reconstruction_error"]=sprint(showerror,error)
    end
    report["forward_reconstruction_seconds"]=time()-began
    solution
end

function _local_optimizer!(m, solver, time_limit;feasibility_only=false)
    solver == "Ipopt" || throw(ArgumentError("Unsupported local solver $solver"))
    set_optimizer(m, Ipopt.Optimizer)
    set_silent(m)
    if solver=="Ipopt"
        set_optimizer_attribute(m, "tol", 1e-8)
        set_optimizer_attribute(m, "bound_relax_factor", 0.0)
        set_optimizer_attribute(m, "max_iter", 1000)
        set_optimizer_attribute(m, "max_cpu_time", max(0.001, time_limit))
        set_optimizer_attribute(m, "max_wall_time", max(0.001, time_limit))
        if feasibility_only
            # With constant objective, any feasible point is optimal. Do not
            # wait for barrier multipliers; check unscaled rows and bounds.
            lower=Float64[];upper=Float64[];violation=Float64[]
            callback=function(mode,iteration,objective,primal,args...)
                mode==0 && primal<=1e-8 || return true
                native=unsafe_backend(m).inner
                resize!(lower,native.n);resize!(upper,native.n);resize!(violation,native.m)
                Ipopt.GetIpoptCurrentViolations(native,false,native.n,lower,upper,
                    C_NULL,C_NULL,C_NULL,native.m,violation,C_NULL)
                max(maximum(lower;init=0.0),maximum(upper;init=0.0),
                    maximum(violation;init=0.0))>1e-8
            end
            MOI.set(m,Ipopt.CallbackFunction(),callback)
        end
    end
end

"""
    solve_case(case; u=ones(Int, length(case.system.generators), length(case.prices)),
               tunnel_opening=nothing, river_on=nothing, time_limit=45.0, warm=nothing, ...)

Optimize continuous dispatch with Ipopt for fixed commitment `u[unit, interval]`
and fixed tunnel gate positions. Enabled river states are fixed by `river_on`,
or inferred from a warm schedule; otherwise they default to on.
Gate positions default to the case's prescribed
schedule; this routine never searches discrete states. Use [`solve`](@ref) for
the global mixed-integer problem, or [`schedule_case`](@ref) to prepare a feasible
schedule by combining MILP proposals with nonlinear dispatch.

Returns a [`SolveResult`](@ref). Its `solution` is a [`Schedule`](@ref), or
`nothing` if no candidate exists. `accepted` means the candidate passed the
physical equations and original operating-window checks. `grid=case.grid`
optionally supplies a finer primary hydraulic mesh; controls retain their
original intervals and the returned schedule records the physical grid.
`replay_grid` optionally adds finer
hydraulic equations sharing these controls; it does not change the returned
grid or objective. Use [`solve_verified`](@ref) for
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
    river_on = nothing,
    replay_grid = nothing,
    grid = c.grid,
    deadline = Inf,
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
        river_on,
        grid,
    )
    replay_grid===nothing || constrain_replay!(b,replay_grid;warm=_dispatch_start(c,u,warm;grid)[1],river_on)
    construction_seconds=time()-b.starttime
    remaining=min(time_limit-construction_seconds,deadline-time())
    result=Dict{String,Any}(
        "case"=>c.name,
        "problem"=>c,
        "grid"=>copy(b.mesh.grid),
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
    _local_optimizer!(b.m, solver, remaining;feasibility_only)
    construction_seconds=time()-b.starttime
    result["construction_seconds"]=construction_seconds
    if time_limit-construction_seconds<=0 || time()>=deadline
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
        solution=_repair_dispatch!(c,solution,result; transport=b.transport,deadline)
        result["accepted"]=get(result,"forward_reconstructed",false)
    end
    result["total_seconds"]=time()-b.starttime
    SolveResult(result;solution)
end
