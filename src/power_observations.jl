"""Pressure-dependent unit conversion and operating limits, in physical units."""
function constrain_unit_physics!(m,s,u,variables,operators)
    (;H,GQ,P,gq,p)=variables
    (;turbine_ops,characteristic_ops,electrical_ops,tailwater_ops,outlet_ops,turbine_flow_ops)=operators
    topology=HydraulicTopology(s)
    ix=topology.node_index
    generator_plants=s.plants[topology.generator_plant]
    plant_generators=Dict(p.name=>topology.plant_generators[i] for (i,p) in enumerate(s.plants))
    heads=Matrix{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,length(s.generators),size(GQ,2))
    for (i, g) in enumerate(s.generators), t in axes(GQ,2)
        plant=generator_plants[i]
        totalq=sum(GQ[j, t] for j in plant_generators[plant.name])
        tail=tailwater_ops[plant.name]===nothing ? 0.0 : tailwater_ops[plant.name](totalq)
        op=outlet_ops[plant.name]
        receiver=op===nothing ? H[ix[plant.target], t] :
                 op isa Real ? op : op(H[ix[plant.target], t])
        hd=unit_power_sign(g)*(H[ix[plant.source], t]-receiver-tail)
        heads[i,t]=hd
        eta=characteristic_ops[i]!==nothing ? characteristic_ops[i][2](hd) : turbine_ops[i]===nothing ?
            g.efficiency-g.qcurvature*((GQ[i, t]-g.qbest)/g.qbest)^2-g.hcurvature*(
            (hd-g.hbest)/g.hbest
        )^2 : turbine_ops[i](GQ[i, t], hd)
        electrical=electrical_ops[i]===nothing ? 1.0 : electrical_ops[i](P[i, t])
        conversion=g.mode==:pump ?
            P[i,t]*eta*g.motor_efficiency-0.00981*GQ[i,t]*hd :
            P[i,t]-0.00981*GQ[i,t]*hd*eta*electrical
        begin
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
    heads
end

"""Evaluate pressure and unit limits at both ends of each hydraulic interval.

These are algebraic states at prescribed storage, not additional integration
steps. Adjustable unit flows and tunnel openings retain their control times.
Each side of a control change has its own pressure state. Mean-step power
continues to determine energy and the objective.
"""
function constrain_power_observations!(b;operational_margin=0.0,prefix="",observation_times=Float64[])
    c,s,m=b.case,b.case.system,b.m
    topology=HydraulicTopology(s)
    nr,nj,ng,ne=length(s.reservoirs),length(s.junctions),length(s.generators),length(s.tunnels)
    all(time->time isa Real && isfinite(time) && first(c.grid)<=time<=last(c.grid),observation_times) ||
        throw(ArgumentError("equipment observation times must lie within the scheduling horizon"))
    ng==0 && ne==0 && return nothing
    periods=Int[];times=Float64[];cells=UnitRange{Int}[]
    for t in eachindex(c.prices)
        points=sort!(unique!(vcat(c.grid[t],
            filter(time->c.grid[t]<time<c.grid[t+1],observation_times),c.grid[t+1])))
        push!(cells,(length(times)+1):(length(times)+length(points)))
        append!(times,points)
        append!(periods,fill(t,length(points)))
    end
    K=length(periods)
    h=@variable(m,[1:nj,1:K],base_name=prefix*"point_head")
    q=@variable(m,[1:ne,1:K],base_name=prefix*"point_tunnel_q")
    # Unit power bounds are imposed over the head span below. Fixing these
    # observation powers too would duplicate equality restrictions when the
    # unit has identical minimum and maximum power.
    p=@variable(m,[1:ng,1:K],lower_bound=0.0,base_name=prefix*"point_power")
    gq=Matrix{VariableRef}(undef,ng,K)
    for k in 1:K,(i,g) in enumerate(s.generators)
        gq[i,k]=g.characteristic===nothing ? b.variables.gq[i,periods[k]] :
            @variable(m,lower_bound=0.0,base_name=prefix*"point_pump_q[$i,$k]")
        set_start_value(p[i,k],something(start_value(b.variables.p[i,periods[k]]),0.0))
        g.characteristic===nothing ||
            set_start_value(gq[i,k],something(start_value(b.variables.gq[i,periods[k]]),0.0))
    end
    H=Matrix{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,length(topology.node_names),K)
    storage=[reservoir_volume_expression(m,c,topology,b.V,b.Q,b.GQ,b.RQ,b.routing,i,time)
        for i in 1:nr,time in times]
    for k in 1:K
        t=periods[k]
        for (i,r) in enumerate(s.reservoirs)
            V=storage[i,k]
            op=b.unit_operators.level_ops[i]
            H[i,k]=op===nothing ? r.z0+r.slope*V+r.curvature*V^2 : op(V)
        end
        for (i,j) in enumerate(s.junctions)
            set_lower_bound(h[i,k],j.hmin/DISPATCH_HEAD_SCALE)
            set_upper_bound(h[i,k],j.hmax/DISPATCH_HEAD_SCALE)
            set_start_value(h[i,k],something(start_value(b.variables.h[nr+i,t]),0.0))
            H[nr+i,k]=DISPATCH_HEAD_SCALE*h[i,k]
        end
        for (i,boundary) in enumerate(s.boundaries)
            H[nr+nj+i,k]=AffExpr(boundary.head)
        end
        for (i,e) in enumerate(s.tunnels)
            cap=opinterval(c,e.name,:capacity,t,e.capacity)/DISPATCH_FLOW_SCALE
            set_lower_bound(q[i,k],e.discharge_river===nothing ? -cap : 0.0)
            set_upper_bound(q[i,k],cap)
            set_start_value(q[i,k],something(start_value(b.variables.q[i,t]),0.0))
            opening=b.gates.opening[i,b.mesh.interval[t]]
            if opening==0
                fix(q[i,k],0.0;force=true)
            else
                factor,resistance=tunnel_coefficients(e,opening)
                source,target=topology.tunnel_nodes[i]
                flow=DISPATCH_FLOW_SCALE*q[i,k]
                @constraint(m,(factor*(H[source,k]-H[target,k])-resistance*flow*abs(flow))/100==0)
            end
        end
    end
    Q=DISPATCH_FLOW_SCALE.*q
    GQ=DISPATCH_FLOW_SCALE.*gq
    P=DISPATCH_POWER_SCALE.*p
    # River transport enters reservoir balances, not junction pressure
    # equations. Its prescribed-storage pressure problem needs no river states.
    pressure_system=_river_replace(s;rivers=River[],river_junctions=RiverJunction[])
    for k in 1:K
        t=periods[k]
        volume=[value(v->something(start_value(v),0.0),storage[i,k]) for i in 1:nr]
        discharge=[DISPATCH_FLOW_SCALE*something(start_value(b.variables.gq[i,t]),0.0) for i in 1:ng]
        reference=[DISPATCH_HEAD_SCALE*something(start_value(b.variables.h[i,t]),0.0)
            for i in eachindex(topology.node_names)]
        point=try
            forward_step(pressure_system,volume,0.0,discharge,Float64[],Float64[];
                tunnel_opening=Float64.(b.gates.opening[:,b.mesh.interval[t]]),
                commitment=b.u[:,t],head_reference=reference)
        catch error
            error isa Union{ArgumentError,DomainError} || rethrow()
            nothing
        end
        point===nothing && continue
        point.converged || continue
        for i in 1:nj
            set_start_value(h[i,k],point.head[nr+i]/DISPATCH_HEAD_SCALE)
        end
        set_start_value.(q[:,k],point.tunnel_q./DISPATCH_FLOW_SCALE)
        set_start_value.(p[:,k],point.power./DISPATCH_POWER_SCALE)
        for (i,g) in enumerate(s.generators)
            g.characteristic===nothing || set_start_value(gq[i,k],point.generator_q[i]/DISPATCH_FLOW_SCALE)
        end
    end
    for k in 1:K,i in (nr+1):(nr+nj)
        net=sum(sign*Q[j,k] for (j,sign) in topology.tunnel_incidence[i];init=0.0)+
            sum(sign*GQ[j,k] for (j,sign) in topology.generator_incidence[i];init=0.0)
        @constraint(m,net/50==0)
    end
    u=b.u[:,periods]
    heads=constrain_unit_physics!(m,s,u,(;H,GQ,P,gq,p),b.unit_operators)
    # Endpoint inequalities are smooth even when a fixed unit output makes
    # both sides active. Do not turn these observation powers into fixed vars.
    for (i,g) in enumerate(s.generators),(k,t) in enumerate(periods)
        b.u[i,t]==0 && continue
        lo,hi=inset_bounds(opinterval(c,g.name,:pmin,t,g.pmin),
            opinterval(c,g.name,:pmax,t,g.pmax),operational_margin)
        @constraint(m,p[i,k]>=lo/DISPATCH_POWER_SCALE)
        @constraint(m,p[i,k]<=hi/DISPATCH_POWER_SCALE)
        if g.characteristic!==nothing
            @constraint(m,gq[i,k]>=opinterval(c,g.name,:qmin,t,g.qmin)/DISPATCH_FLOW_SCALE)
            @constraint(m,gq[i,k]<=opinterval(c,g.name,:qmax,t,g.qmax)/DISPATCH_FLOW_SCALE)
        end
    end
    increasing=power_increases_with_head.(s.generators)
    eta_bounds=unit_efficiency_bounds.(s.generators)
    monotone_flow=discharge_envelope_monotone.(s.generators)
    # A continuous head trajectory attains every head between its endpoint
    # values. Enforce the equipment limits over that whole span, including
    # interior efficiency-table extrema, without adding integration steps.
    for (i,g) in enumerate(s.generators),t in eachindex(c.prices),k in first(cells[t]):last(cells[t])-1
        b.u[i,t]==0 && continue
        pmin,pmax=inset_bounds(opinterval(c,g.name,:pmin,t,g.pmin),
            opinterval(c,g.name,:pmax,t,g.pmax),operational_margin)
        qmin=opinterval(c,g.name,:qmin,t,g.qmin)
        qmax=opinterval(c,g.name,:qmax,t,g.qmax)
        for (component,scale) in ((:upper,DISPATCH_POWER_SCALE),(:lower,DISPATCH_POWER_SCALE),
            (:domain,1.0),(:flow,DISPATCH_FLOW_SCALE))
            component in (:upper,:lower) && increasing[i] && continue
            component==:domain && first(eta_bounds[i])>=g.min_efficiency && last(eta_bounds[i])<=1.0 && continue
            component==:flow && monotone_flow[i] && continue
            envelope=autodiff_operator(m,3,
                (q,a,b)->getproperty(unit_head_envelope(g,q,a,b;pmin,pmax,qmin,qmax),component)/scale;
                name=Symbol(prefix,"unit_head_span_",component,"_",i,"_",k))
            @constraint(m,envelope(b.GQ[i,t],heads[i,k],heads[i,k+1])<=0)
        end
    end
    constrain_plant_capacities!(m,c,u,P,GQ;intervals=periods,margin=operational_margin)
    nothing
end

"""Operating residuals from instantaneous pressure states at hydraulic vertices.

Both limits at control changes use their own held controls. Reservoir volumes
are prescribed; this calculation advances no time and changes no trajectory.
It also checks each unit over the head span between its two observations.
This remains a finite observation check, not a bound on every continuous-time
head excursion or interior peak of total plant power.
"""
function point_operating_residuals(c,x;path_data,tolerance=1e-4)
    s=c.system
    ng,ne=length(s.generators),length(s.tunnels)
    unit_power=zeros(ng);unit_q=zeros(ng);unit_domain=zeros(ng)
    plant_power=zeros(length(s.plants));plant_q=similar(plant_power);fill!(plant_q,0.0)
    tunnel_q=zeros(ne);pressure_head=zeros(length(s.junctions))
    times=Float64[]
    convergence=0.0
    ng==0 && ne==0 && return (;unit_power,unit_q,unit_domain,plant_power,plant_q,tunnel_q,pressure_head,times,convergence)
    pressure_system=_river_replace(s;rivers=River[],river_junctions=RiverJunction[])
    topology=HydraulicTopology(s)
    nr=length(s.reservoirs)
    observed_heads=Vector{Matrix{Float64}}(undef,length(c.prices))
    observations=Vector{Vector{Float64}}(undef,length(c.prices))
    for t in eachindex(c.prices)
        observations[t]=sort!(unique!(vcat(c.grid[t],
            filter(time->c.grid[t]<time<c.grid[t+1],path_data.crossings),c.grid[t+1])))
        observed_heads[t]=fill(NaN,ng,length(observations[t]))
    end
    for t in eachindex(c.prices),(k,time) in enumerate(observations[t])
        vertex=searchsortedfirst(c.grid,time)
        volume=vertex<=length(c.grid) && c.grid[vertex]==time ? x.storage[:,vertex] :
            [path(time) for path in path_data.paths]
        point=try
            forward_step(pressure_system,volume,0.0,
                view(x.generator_q,:,t),Float64[],Float64[];
                tunnel_opening=view(x.tunnel_opening,:,t),commitment=view(x.commitment,:,t),
                head_reference=view(x.head,:,t))
        catch error
            error isa Union{ArgumentError,DomainError} || rethrow()
            nothing
        end
        if point===nothing || !point.converged
            convergence=Inf
            push!(times,time)
            continue
        end
        observed_heads[t][:,k].=point.generator_head
        violation=0.0
        for (i,g) in enumerate(s.generators)
            u=x.commitment[i,t]
            q,p=point.generator_q[i],point.power[i]
            vp=max(0.0,opinterval(c,g.name,:pmin,t,g.pmin)*u-p,
                p-opinterval(c,g.name,:pmax,t,g.pmax)*u)
            vq=max(0.0,opinterval(c,g.name,:qmin,t,g.qmin)*u-q,
                q-opinterval(c,g.name,:qmax,t,g.qmax)*u)
            domain=0.0
            if u==1
                h=point.generator_head[i]
                eta=efficiency(g,q,h)
                domain=max(0.0,g.hmin-h,h-g.hmax,g.min_efficiency-eta,eta-1.0)
                if g.turbine_table!==nothing
                    domain=max(domain,turbine_qmin(g.turbine_table,h;extrapolation=:linear)-q,
                        q-turbine_qmax(g.turbine_table,h;extrapolation=:linear))
                end
            end
            unit_power[i]=max(unit_power[i],vp)
            unit_q[i]=max(unit_q[i],vq)
            unit_domain[i]=max(unit_domain[i],domain)
            violation=max(violation,vp,vq,domain)
        end
        for (i,p) in enumerate(s.plants)
            members=topology.plant_generators[i]
            on=any(x.commitment[j,t]==1 for j in members)
            power_sum=sum(point.power[j] for j in members;init=0.0)
            flow_sum=sum(point.generator_q[j] for j in members;init=0.0)
            vp=max(0.0,power_sum-opinterval(c,p.name,:pmax,t,p.pmax),
                opinterval(c,p.name,:pmin,t,p.pmin)*on-power_sum)
            upper=optional_operation(c,p.name,:qmax,t;default=p.qmax)
            vq=max(0.0,opinterval(c,p.name,:qmin,t,p.qmin)*on-flow_sum,
                upper===nothing ? 0.0 : flow_sum-upper)
            plant_power[i]=max(plant_power[i],vp)
            plant_q[i]=max(plant_q[i],vq)
            violation=max(violation,vp,vq)
        end
        for (i,e) in enumerate(s.tunnels)
            excess=max(0.0,abs(point.tunnel_q[i])-opinterval(c,e.name,:capacity,t,e.capacity))
            tunnel_q[i]=max(tunnel_q[i],excess)
            violation=max(violation,excess)
        end
        for (i,j) in enumerate(s.junctions)
            h=point.head[nr+i]
            excess=max(0.0,j.hmin-h,h-j.hmax)
            pressure_head[i]=max(pressure_head[i],excess)
            violation=max(violation,excess)
        end
        violation>tolerance && push!(times,time)
    end
    for t in eachindex(c.prices),(i,g) in enumerate(s.generators),k in 1:length(observations[t])-1
        x.commitment[i,t]==0 && continue
        a,b=observed_heads[t][i,k],observed_heads[t][i,k+1]
        isfinite(a) && isfinite(b) || continue
        pmin=opinterval(c,g.name,:pmin,t,g.pmin)
        pmax=opinterval(c,g.name,:pmax,t,g.pmax)
        qmin=opinterval(c,g.name,:qmin,t,g.qmin)
        qmax=opinterval(c,g.name,:qmax,t,g.qmax)
        span=unit_head_envelope(g,x.generator_q[i,t],a,b;pmin,pmax,qmin,qmax)
        # The transformed inequalities are equivalent to power bounds whenever
        # efficiency is positive. Evaluate actual MW at their extremal heads;
        # never report shaft-work residuals as electrical power.
        for h in (span.at_upper,span.at_lower)
            q=g.characteristic===nothing ? x.generator_q[i,t] :
                pump_discharge(g.characteristic,h;extrapolation=:linear)
            eta=efficiency(g,q,h)
            if eta>0
                p=power(g,q,h)
                unit_power[i]=max(unit_power[i],pmin-p,p-pmax)
            end
        end
        unit_q[i]=max(unit_q[i],span.flow)
        unit_domain[i]=max(unit_domain[i],span.domain)
        if max(unit_power[i],unit_q[i],unit_domain[i])>tolerance
            push!(times,observations[t][k],observations[t][k+1])
        end
    end
    (;unit_power,unit_q,unit_domain,plant_power,plant_q,tunnel_q,pressure_head,
        times=sort!(unique!(times)),convergence)
end
