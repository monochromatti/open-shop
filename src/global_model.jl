# Finite exact polynomial reformulation of the dispatch model.
# Deterministic reaches preserve rectangular cohorts exactly. Distributed,
# flow-dependent reaches use the existing exact quadratic transfer coefficients
# on the scheduling grid; downstream junction mixing is a grid approximation.
# A global bound for that formulation does not bound continuous-time mixing.
function _global_analytic_eta_bounds(g, qlo, qhi, hlo, hhi)
    function contribution(coef, best, lo, hi)
        vals=[-coef*((lo-best)/best)^2, -coef*((hi-best)/best)^2]
        lo<=best<=hi && push!(vals, 0.0)
        extrema(vals)
    end
    a, b=contribution(g.qcurvature, g.qbest, qlo, qhi)
    c, d=contribution(g.hcurvature, g.hbest, hlo, hhi)
    (g.efficiency+a+c, g.efficiency+b+d)
end

function _build_global_dispatch(
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
    bound_start=time()
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
    domains=_global_capacity_bounds(c, nd, rd)
    bounds_seconds=time()-bound_start
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
    node_bounds=Dict{Tuple{Symbol,Int},Tuple{Float64,Float64}}()
    for (i, r) in enumerate(s.reservoirs), t in 1:T
        node_bounds[(r.name, t)]=(domains.hlo[i, t], domains.hhi[i, t])
    end
    for j in s.junctions, t in 1:T
        node_bounds[(j.name, t)]=(j.hmin, j.hmax)
    end
    for z in s.boundaries, t in 1:T
        node_bounds[(z.name, t)]=(z.head, z.head)
    end
    V=[s.reservoirs[i].vmax*v[i, t] for i in 1:R, t in 1:(T + 1)]
    H=250 .* h
    Q=50 .* q
    GQ=50 .* gq
    P=40 .* p
    RQ=100 .* rq
    constrain_flow_requirements!(m, c, GQ, RQ)
    for (i, r) in enumerate(s.reservoirs)
        for t in 1:(T + 1)
            lo, hi=domains.lower[i, t], domains.upper[i, t]
            set_lower_bound(v[i, t], lo/r.vmax)
            set_upper_bound(v[i, t], hi/r.vmax)
            set_start_value(v[i, t], r.v0/r.vmax)
        end
        fix(v[i, 1], r.v0/r.vmax; force = true)
        for t in 1:T
            set_lower_bound(h[i, t], node_bounds[(r.name, t)][1]/250)
            set_upper_bound(h[i, t], node_bounds[(r.name, t)][2]/250)
            mid=(V[i, t]+V[i, t + 1])/2
            midlo=(domains.lower[i, t]+domains.lower[i, t + 1])/2
            midhi=(domains.upper[i, t]+domains.upper[i, t + 1])/2
            level=r.level_curve===nothing ? r.z0+r.slope*mid+r.curvature*mid^2 :
                  _global_tensor_table!(
                m,
                r.level_curve,
                mid,
                midlo,
                midhi;
                name = Symbol("level_", i, "_", t),
            )
            @constraint(m, (H[i, t]-level)/250==0)
        end
    end
    level_at_vertex=(i,t)->begin
        r=s.reservoirs[i]
        volume=V[i,t]
        r.level_curve===nothing ? r.z0+r.slope*volume+r.curvature*volume^2 :
            _global_tensor_table!(m,r.level_curve,volume,domains.lower[i,t],domains.upper[i,t];
                name=Symbol("level_vertex_",i,"_",t))
    end
    constrain_reservoir_ramps!(m,c,V; level=level_at_vertex)
    for (i, j) in enumerate(s.junctions), t in 1:T
        set_lower_bound(h[R + i, t], node_bounds[(j.name,t)][1]/250)
        set_upper_bound(h[R + i, t], node_bounds[(j.name,t)][2]/250)
        set_start_value(h[R + i, t], (j.hmin+j.hmax)/500)
    end
    for (i, b) in enumerate(s.boundaries), t in 1:T
        fix(h[R + length(s.junctions) + i, t], b.head/250; force = true)
    end
    for (i, e) in enumerate(s.tunnels), t in 1:T
        cap=opinterval(c, e.name, :capacity, t, e.capacity)
        opening=opinterval(c, e.name, :opening, t, e.opening)
        qlo,qhi=-cap,cap
        set_lower_bound(q[i, t], (e.discharge_river===nothing ? qlo : max(0.0,qlo))/50)
        set_upper_bound(q[i, t], qhi/50)
        if opening==0
            fix(q[i, t], 0; force = true)
        elseif e.discharge_river!==nothing
            @constraint(m,(opening*(H[ix[e.source],t]-H[ix[e.target],t])-e.resistance*Q[i,t]^2)/100==0)
        else
            qp=@variable(
                m,
                lower_bound=0,
                upper_bound=cap,
                base_name="tunnel_positive_$(i)_$(t)"
            )
            qm=@variable(
                m,
                lower_bound=0,
                upper_bound=cap,
                base_name="tunnel_negative_$(i)_$(t)"
            )
            direction=@variable(m, binary=true, base_name="tunnel_direction_$(i)_$(t)")
            @constraint(m, qp<=cap*direction)
            @constraint(m, qm<=cap*(1-direction))
            @constraint(m, Q[i, t]==qp-qm)
            @constraint(
                m,
                (opening*(H[ix[e.source], t]-H[ix[e.target], t])-e.resistance*(qp^2-qm^2))/100==0
            )
        end
    end
    shared_heads=Dict{Tuple{Symbol,Int},VariableRef}()
    for (i, g) in enumerate(s.generators), t in 1:T
        plant=generator_plants[i]
        key=(plant.name, t)
        if haskey(shared_heads, key)
            hd=shared_heads[key]
            hlo=lower_bound(hd)
            hhi=upper_bound(hd)
        else
            totalq=sum(GQ[j, t] for j in plant_generators[plant.name])
            aggregate_max=sum(
                s.generators[j].qmax for j in plant_generators[plant.name];
                init = 0.0,
            )
            tail=plant.tailwater_curve===nothing ? 0.0 :
                 _global_tensor_table!(
                m,
                plant.tailwater_curve,
                totalq,
                0.0,
                aggregate_max;
                name = Symbol("tailwater_", i, "_", t),
            )
            tail_lo=plant.tailwater_curve===nothing ? 0.0 : minimum(plant.tailwater_curve.y)
            tail_hi=plant.tailwater_curve===nothing ? 0.0 : maximum(plant.tailwater_curve.y)
            src=node_bounds[(plant.source, t)]
            dst=node_bounds[(plant.target, t)]
            receiver=H[ix[plant.target], t]
            floor=plant.outlet_head_floor
            if floor!==nothing
                receiver=dst[2]<=floor ? floor :
                         dst[1]>=floor ? receiver :
                         _global_tensor_table!(
                    m,
                    TableCurve([dst[1], floor, dst[2]], [floor, floor, dst[2]]),
                    receiver,
                    dst[1],
                    dst[2];
                    name = Symbol("outlet_head_", i, "_", t),
                )
                dst=(max(dst[1], floor), max(dst[2], floor))
            end
            hlo=src[1]-dst[2]-tail_hi
            hhi=src[2]-dst[1]-tail_lo
            hd=@variable(
                m,
                lower_bound=hlo,
                upper_bound=hhi,
                base_name="net_head_$(i)_$(t)"
            )
            @constraint(m, hd==H[ix[plant.source], t]-receiver-tail)
            shared_heads[key]=hd
        end
        flowmax=g.qmax
        powmax=g.pmax
        # Fixed states need only their physical branch; retaining a redundant
        # off/on disjunction creates degenerate table equations in presolve.
        known_state=joint ? (is_fixed(u[i,t]) ? fix_value(u[i,t]) : nothing) : u[i,t]
        flowmin=known_state==1 ? opinterval(c,g.name,:qmin,t,g.qmin) : 0.0
        known_state==0 && (flowmax=0.0)
        eta=if g.turbine_table===nothing
            emin, emax=_global_analytic_eta_bounds(g, flowmin, flowmax, hlo, hhi)
            z=@variable(m, lower_bound=emin, upper_bound=emax, base_name="eta_$(i)_$(t)")
            @constraint(
                m,
                z==g.efficiency-g.qcurvature*((GQ[i, t]-g.qbest)/g.qbest)^2-g.hcurvature*(
                    (hd-g.hbest)/g.hbest
                )^2
            )
            z
        else
            _global_tensor_turbine!(m,g.turbine_table,GQ[i,t],hd,flowmin,flowmax,hlo,hhi;
                name=Symbol("turbine_",i,"_",t))
        end
        electrical=g.generator_efficiency_curve===nothing ? 1.0 :
                   _global_tensor_table!(
            m,
            g.generator_efficiency_curve,
            P[i, t],
            0.0,
            powmax;
            name = Symbol("electrical_", i, "_", t),
        )
        qmin=opinterval(c, g.name, :qmin, t, g.qmin)
        qmax=opinterval(c, g.name, :qmax, t, g.qmax)
        pmin=opinterval(c, g.name, :pmin, t, g.pmin)
        pmax=opinterval(c, g.name, :pmax, t, g.pmax)
        effective=GQ[i,t]*eta
        if joint
            set_upper_bound(gq[i, t], flowmax/50)
            set_upper_bound(p[i, t], powmax/40)
            if (fixed_u!==nothing && fixed_u[i, t]==0) ||
               (free_mask!==nothing && !free_mask[i, t] && ustart[i, t]==0)
                fix(gq[i, t], 0.0; force = true)
                fix(p[i, t], 0.0; force = true)
            end
            @constraint(m, GQ[i, t]>=qmin*u[i, t])
            @constraint(m, GQ[i, t]<=qmax*u[i, t])
            @constraint(m, P[i, t]>=(pmin+operational_margin)*u[i, t])
            @constraint(m, P[i, t]<=max(0.0, pmax-operational_margin)*u[i, t])
            @constraint(m, (P[i, t]-0.00981*effective*hd*electrical)/40==0)
            ranges=(
                head_min = hlo,
                head_max = hhi,
                efficiency_min = lower_bound(eta),
                efficiency_max = upper_bound(eta),
            )
            @constraint(m, hd>=g.hmin-max(0.0, g.hmin-ranges.head_min)*(1-u[i, t]))
            @constraint(m, hd<=g.hmax+max(0.0, ranges.head_max-g.hmax)*(1-u[i, t]))
            @constraint(
                m,
                eta>=g.min_efficiency-max(0.0, g.min_efficiency-ranges.efficiency_min)*(
                    1-u[i, t]
                )
            )
            @constraint(m, eta<=1.0 + max(0.0, ranges.efficiency_max-1.0)*(1-u[i, t]))
        else
            set_lower_bound(gq[i, t], qmin*u[i, t]/50)
            set_upper_bound(gq[i, t], qmax*u[i, t]/50)
            set_lower_bound(p[i, t], (pmin+operational_margin)*u[i, t]/40)
            set_upper_bound(p[i, t], max(0.0, pmax-operational_margin)*u[i, t]/40)
            if u[i, t]==0
                fix(gq[i, t], 0; force = true)
                fix(p[i, t], 0; force = true)
            else
                @constraint(m, (P[i, t]-0.00981*effective*hd*electrical)/40==0)
                @constraint(m, g.hmin<=hd<=g.hmax)
                @constraint(m, g.min_efficiency<=eta<=1.0)

            end
        end
        if g.turbine_table!==nothing
            table=g.turbine_table
            qlo=_global_tensor_table!(m,TableCurve(table.heads,table.qmin),hd,hlo,hhi;
                name=Symbol("qlo_",i,"_",t))
            qhi=_global_tensor_table!(m,TableCurve(table.heads,table.qmax),hd,hlo,hhi;
                name=Symbol("qhi_",i,"_",t))
            flow_margin=operational_margin/(
                0.00981*g.hbest*max(g.min_efficiency, g.efficiency, eps(Float64))
            )
            if joint
                @constraint(
                    m,
                    GQ[i, t]>=qlo+flow_margin-max(0.0, upper_bound(qlo)+flow_margin)*(
                        1-u[i, t]
                    )
                )
                @constraint(
                    m,
                    GQ[i, t]<=qhi-flow_margin+max(0.0, g.qmax-lower_bound(qhi)+flow_margin)*(
                        1-u[i, t]
                    )
                )
            elseif u[i, t]!=0
                @constraint(m, GQ[i, t]>=qlo+flow_margin)
                @constraint(m, GQ[i, t]<=qhi-flow_margin)
            end
        end
    end
    constrain_dispatch_operations!(m,c,u,P,GQ,RQ; margin=operational_margin, transitions)
    injections=river_injections(c,RQ,GQ,Q)
    coordinates=exact ? nothing : [
        river_reference_coordinate!(m,r,RQ[i,t];name=Symbol("routing_",i,"_",t))
        for (i,r) in enumerate(s.rivers),t in 1:T
    ]
    transfer=(i,k,b)->river_transfer_expression(m,s.rivers[i],RQ[i,k],b; coordinate=coordinates[i,k])
    arrivals=Matrix{Any}(undef, D, T)
    terminal=Any[]
    for (i, r) in enumerate(s.rivers)
        B=exact ? nothing : rd["B"][i]
        for t in 1:T
            cap=opinterval(c, r.name, :capacity, t, r.capacity)
            set_upper_bound(rq[i, t], cap/100)
            set_lower_bound(rq[i,t],opinterval(c,r.name,:inflow,t,r.inflow)/100)
            set_lower_bound(a[i, t], opinterval(c, r.name, :gate_min, t, r.gate_min))
            set_upper_bound(a[i, t], opinterval(c, r.name, :gate_max, t, 1.0))
            requirement=opinterval(c, r.name, :min_release, t, 0.0)
            penalty=opinterval(c, r.name, :release_penalty, t, 0.0)
            set_upper_bound(shortfall_release[i, t], 0.0036*dt[t]*max(0.0, requirement))
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
                    v*injections[d, k] for (d, k, v) in nd.arrival_terms[i, t];
                    init = 0.0,
                )
            else
                # Skip structural zeros to preserve banded routing sparsity.
                ks=[k for k in 1:T if any(!iszero,view(B,t,k,:))]
                arrivals[i, t]=rd["history_arrival"][i, t]+sum(
                    0.0036*dt[k]*transfer(i,k,collect(view(B,t,k,:))) for k in ks;
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
            natural=opinterval(c,r.name,:inflow,t,r.inflow)
            outlet_release=RQ[i,t]-natural
            if r.law==:junction
                fix(a[i, t], 0.0; force = true)
            elseif r.discharge_curve!==nothing
                r.law==:weir && fix(a[i, t], 1.0; force = true)
                level=H[ix[r.source], t]
                @constraint(m, first(r.discharge_curve.x)<=level<=last(r.discharge_curve.x))
                lawlo,lawhi=first(r.discharge_curve.x),last(r.discharge_curve.x)
                lawlo<=lawhi || throw(ArgumentError("river law outside reachable domain for $(r.name)"))
                discharge=_global_tensor_table!(
                    m,
                    r.discharge_curve,
                    level,
                    lawlo,
                    lawhi;
                    name = Symbol("river_law_", i, "_", t),
                )
                @constraint(m, (outlet_release-a[i, t]*discharge)/100==0)
            elseif r.law==:controlled
                @constraint(m, outlet_release==(cap-natural)*a[i, t])
            elseif r.law in (:orifice, :weir)
                r.law==:weir && fix(a[i, t], 1.0; force = true)
                hd=H[ix[r.source], t]-r.crest
                lo, hi=node_bounds[(r.source, t)]
                lo-=r.crest
                hi-=r.crest
                hi>=0 ||
                    r.allow_dry ||
                    throw(ArgumentError("river $(r.name) has no wet domain"))
                wet=@variable(
                    m,
                    lower_bound=0.0,
                    upper_bound=max(0.0, hi),
                    base_name="wet_$(i)_$(t)"
                )
                if !r.allow_dry || lo>=0
                    @constraint(m, wet==hd)
                elseif hi<=0
                    fix(wet, 0.0; force = true)
                else
                    z=@variable(m, binary=true, base_name="wet_branch_$(i)_$(t)")
                    @constraint(m, wet>=hd)
                    @constraint(m, wet<=hi*z)
                    @constraint(m, wet<=hd-lo*(1-z))
                end
                root=@variable(
                    m,
                    lower_bound=0.0,
                    upper_bound=sqrt(max(0.0, hi)),
                    base_name="wet_root_$(i)_$(t)"
                )
                @constraint(m, root^2==wet)
                if r.law==:orifice
                    @constraint(m, (outlet_release-r.coefficient*a[i, t]*root)/100==0)
                else
                    @constraint(m, (outlet_release-r.coefficient*wet*root)/100==0)
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
            push!(
                terminal,
                nd.terminal_history[i]+sum(
                    v*injections[d, k] for (d, k, v) in nd.terminal_terms[i];
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
                    rate=ex.history+sum(v*injections[d, k] for (d, k, v) in ex.terms; init = 0.0)
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
                    0.0036*dt[k]*transfer(i,k,[1-sum(view(B,:,k,l)) for l in 1:river_curve_count(r)]) for k in 1:T
                ),
            )
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
    operating_cost=transition_costs(c,u; transitions)
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
    b=(;
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
        domains,
        bounds_seconds,
        shared_heads,
        node_bounds,
    )
    _add_power_bounds!(b,c)
    joint && _add_power_hull!(b,c)
    joint && _add_table_power_bounds!(b,c)
    for variable in all_variables(m)
        bounded=is_binary(variable) || (
            is_fixed(variable) ? isfinite(fix_value(variable)) :
            has_lower_bound(variable) &&
            has_upper_bound(variable) &&
            isfinite(lower_bound(variable)) &&
            isfinite(upper_bound(variable))
        )
        bounded || throw(
            ArgumentError("global solver requires finite bounds for $(name(variable))"),
        )
    end
    b
end
