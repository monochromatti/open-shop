"""A linear surrogate proposes commitments; its objective/bound is NOT a bound
for the nonlinear scheduling problem. Every proposed commitment needs an NLP.
"""
function propose_commitment(
    c;
    water_scale = 1.0,
    reference = nothing,
    time_limit = 20.0,
    incumbent = nothing,
    radius = nothing,
    deadline = Inf,
    transport = nothing,
)
    started=time()
    validate_inputs(c)
    isfinite(time_limit) && time_limit>0 ||
        throw(ArgumentError("Proposal time limit must be finite and positive"))
    isnan(deadline) && throw(ArgumentError("Invalid proposal deadline"))
    s=c.system
    G=length(s.generators)
    T=length(c.prices)
    R=length(s.reservoirs)
    E=length(s.tunnels)
    D=length(s.rivers)
    N=length(nodes(s))
    ix=nodeindex(s)
    dt=diff(c.grid)
    exact=all(r.deterministic_delay!==nothing for r in s.rivers)
    data=exact ? nothing : routing_data(c)
    nd=exact ? _transport_data(c, transport) : nothing
    !exact &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed proposal",
            ),
        )
    m=Model(HiGHS.Optimizer)
    set_silent(m)
    set_optimizer_attribute(m, "mip_rel_gap", 0.002)
    states=_joint_states!(m,c)
    u=states.u
    @variable(m, V[1:R, 1:(T + 1)])
    @variable(m, Q[1:E, 1:T])
    @variable(m, gq[1:G, 1:T]>=0)
    @variable(m, rq[1:D, 1:T]>=0)
    @variable(m, release_shortfall[1:D, 1:T]>=0)
    constrain_flow_requirements!(m, c, gq, rq)
    alpha=zeros(G, T)
    for (j, g) in enumerate(s.generators), t in 1:T
        h=g.hbest
        if reference!==nothing && haskey(reference, "H")
            p=plantof(s, g)
            qtotal=haskey(reference, "generator_q") ?
                   sum(
                reference["generator_q"][k, t] for
                (k, z) in enumerate(s.generators) if z.plant==g.plant
            ) : 0.0
            h=reference["H"][ix[p.source], t]-outlet_head(
                p,
                reference["H"][ix[p.target], t],
            )-tailwater(p, qtotal)
            h=clamp(h, g.hmin, g.hmax)
        end
        qref=clamp(g.qbest, g.qmin, g.qmax)
        alpha[j, t]=power(g, qref, h)/qref
        qlo=opvalue(c, g.name, :qmin, c.grid[t], g.qmin)
        qhi=opvalue(c, g.name, :qmax, c.grid[t], g.qmax)
        if g.turbine_table!==nothing
            qlo=max(qlo, turbine_qmin(g.turbine_table, h; extrapolation = :linear))
            qhi=min(qhi, turbine_qmax(g.turbine_table, h; extrapolation = :linear))
        end
        plo=opvalue(c, g.name, :pmin, c.grid[t], g.pmin)
        phi=opvalue(c, g.name, :pmax, c.grid[t], g.pmax)
        @constraint(m, gq[j, t]>=qlo*u[j, t])
        @constraint(m, gq[j, t]<=qhi*u[j, t])
        @constraint(m, alpha[j, t]*gq[j, t]>=plo*u[j, t])
        @constraint(m, alpha[j, t]*gq[j, t]<=phi*u[j, t])
    end
    P=[alpha[j, t]*gq[j, t] for j in 1:G, t in 1:T]
    constrain_dispatch_operations!(m,c,u,P,gq,rq; transitions=states, proposal=true)
    constrain_reservoir_ramps!(m,c,V) # Only linear storage ramps belong to the proposal.
    for (i, r) in enumerate(s.reservoirs)
        for t in 1:(T + 1)
            lo, hi=storage_bounds(c, r, c.grid[t])
            set_lower_bound(V[i, t], lo)
            set_upper_bound(V[i, t], hi)
        end
        fix(V[i, 1], r.v0; force = true)
    end
    for (j, e) in enumerate(s.tunnels), t in 1:T
        cap=opvalue(c, e.name, :opening, c.grid[t], e.opening)==0 ? 0.0 :
            opvalue(c, e.name, :capacity, c.grid[t], e.capacity)
        set_lower_bound(Q[j, t], e.discharge_river===nothing ? -cap : 0.0)
        set_upper_bound(Q[j, t], cap)
    end
    injections=river_injections(c,rq,gq,Q)
    arrivals=Matrix{Any}(undef, D, T)
    terminal=Any[]
    for (d, r) in enumerate(s.rivers)
        K=nothing
        if !exact
            B=data["B"][d]
            qref=reference!==nothing && haskey(reference,"river_release") ?
                clamp.(reference["river_release"][d,:],0.0,r.capacity) : fill(r.capacity/2,T)
            K=[r.deterministic_delay!==nothing ? B[t,k,1] :
                sum(w*B[t,k,l] for (l,w) in RiverRouting.blend(r.curves,qref[k]))
                for t in 1:T,k in 1:T]
        end
        for t in 1:T
            set_upper_bound(rq[d, t], opvalue(c, r.name, :capacity, c.grid[t], r.capacity))
            set_lower_bound(rq[d,t],opinterval(c,r.name,:inflow,t,r.inflow))
            minimum_release=opvalue(c, r.name, :min_release, c.grid[t], 0.0)
            penalty=opvalue(c, r.name, :release_penalty, c.grid[t], 0.0)
            penalty>0 || fix(release_shortfall[d, t], 0; force = true)
            @constraint(m, rq[d, t]+release_shortfall[d, t]>=minimum_release)
            arrivals[d, t]=exact ?
                           nd.arrival_history[d, t]+sum(
                v*injections[i, k] for (i, k, v) in nd.arrival_terms[d, t];
                init = 0.0,
            ) :
                           data["history_arrival"][d, t]+sum(
                0.0036*dt[k]*K[t, k]*rq[d, k] for k in 1:T
            )
            @constraint(
                m,
                arrivals[d, t]>=0.0036*dt[t]*opvalue(
                    c,
                    r.name,
                    :min_arrival,
                    c.grid[t],
                    r.min_arrival,
                )
            )
        end
        push!(
            terminal,
            exact ?
            nd.terminal_history[d]+sum(
                v*injections[i, k] for (i, k, v) in nd.terminal_terms[d];
                init = 0.0,
            ) :
            data["history_terminal"][d]+sum(
                0.0036*dt[k]*(1-sum(K[:, k]))*rq[d, k] for k in 1:T
            ),
        )
    end
    constrain_river_sources!(m,c,rq,arrivals,injections)
    # Hydraulic node continuity is retained, but pressure-loss equations and
    # nonlinear outlet laws are deliberately absent from this proposal model.
    for i in 1:(R + length(s.junctions)), t in 1:T
        name=nodes(s)[i]
        net=sum(
            water_incidence(e,name)*Q[j, t] for
            (j, e) in enumerate(s.tunnels);
            init = 0.0,
        )
        for (j, g) in enumerate(s.generators)
            p=plantof(s, g)
            net+=water_incidence(p,name)*gq[j, t]
        end
        net+=sum(
            (r.target==name ? arrivals[d, t]/(0.0036*dt[t]) : 0)-(
                r.source==name ? rq[d, t]-opinterval(c,r.name,:inflow,t,r.inflow) : 0
            ) for (d, r) in enumerate(s.rivers);
            init = 0.0,
        )
        if i<=R
            @constraint(
                m,
                V[i, t + 1]-V[i, t]==0.0036*dt[t]*(
                    net+opaverage(
                        c,
                        s.reservoirs[i].name,
                        :inflow,
                        c.grid[t],
                        c.grid[t + 1],
                        s.reservoirs[i].inflow,
                    )
                )
            )
        else
            @constraint(m, net==0)
        end
    end
    if incumbent!==nothing
        size(incumbent)==(G, T) || error("Incorrect incumbent shape")
        all(x->x==0||x==1, incumbent) || error("Nonbinary incumbent")
        for j in 1:G, t in 1:T
            set_start_value(u[j, t], incumbent[j, t])
        end
        if radius!==nothing
            @constraint(
                m,
                sum(incumbent[j, t]==1 ? 1-u[j, t] : u[j, t] for j in 1:G, t in 1:T)<=radius
            )
        end
    end
    valuechange=sum(
        r.water_value*(V[i, T + 1]-r.v0) for (i, r) in enumerate(s.reservoirs);
        init = 0.0,
    ) + sum(
        r.water_value*(
            terminal[d]-(exact ? nd.initial_transit[d] : data["history_initial"][d])
        ) for (d, r) in enumerate(s.rivers);
        init = 0.0,
    )
    @objective(
        m,
        Max,
        sum(
            c.prices[t]*dt[t]*P[j, t] for j in 1:G, t in 1:T
        )-transition_costs(c,u; transitions=states)+water_scale*valuechange-sum(
            0.0036*dt[t]*opvalue(c, r.name, :release_penalty, c.grid[t], 0.0)*release_shortfall[
                d,
                t,
            ] for (d, r) in enumerate(s.rivers), t in 1:T;
            init = 0.0,
        )
    )
    construction=time()-started
    remaining=min(time_limit-construction, deadline-time())
    remaining<=0 && return Dict{String,Any}(
        "status"=>"CONSTRUCTION_BUDGET_EXHAUSTED",
        "seconds"=>0.0,
        "construction_seconds"=>construction,
        "total_seconds"=>time()-started,
        "water_scale"=>water_scale,
        "nonlinear_bound"=>false,
    )
    set_time_limit_sec(m, remaining)
    elapsed=@elapsed optimize!(m)
    result=Dict{String,Any}(
        "status"=>string(termination_status(m)),
        "seconds"=>elapsed,
        "construction_seconds"=>construction,
        "water_scale"=>water_scale,
        "nonlinear_bound"=>false,
    )
    if has_values(m)
        result["u"]=round.(Int, value.(u))
        result["generator_q"]=value.(gq)
        result["river_release"]=value.(rq)
        result["proposal_objective"]=objective_value(m)
    end
    result["total_seconds"]=time()-started
    result
end
