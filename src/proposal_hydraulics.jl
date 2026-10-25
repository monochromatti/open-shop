# A gate proposal must distinguish the heads supplied by alternative intakes.
# These affine head/storage and secant loss equations are a scheduling surrogate,
# never a physical audit or a bound on the nonlinear problem. They are used only
# when the case contains selectable gates. Prescribed-gate proposals instead
# enclose omitted pressure laws through their physical flow domains.
function _proposal_reference(c, point)
    if point isa Schedule || point isa SolveResult
        point=_schedule_input(point;case=c)
        return (;point.storage,point.head,point.tunnel_q,point.generator_q,point.river_release)
    end
    # A linearization point may prescribe only some physical quantities.
    # Missing components retain the surrogate's existing nominal defaults.
    matrix(key,rows,columns) = point===nothing || !haskey(point,key) ? nothing : begin
        values=_schedule_matrix(point[key],columns)
        size(values)==(rows,columns) || throw(DimensionMismatch("proposal reference $key"))
        values
    end
    s,T=c.system,length(c.prices)
    (;storage=matrix("V",length(s.reservoirs),T+1),
        head=matrix("H",length(nodes(s)),T),
        tunnel_q=matrix("tunnel_q",length(s.tunnels),T),
        generator_q=matrix("generator_q",length(s.generators),T),
        river_release=matrix("river_release",length(s.rivers),T))
end

function _proposal_gate_hydraulics!(m,c,V,Q,gq,u,alpha,reference_heads,gates,reference)
    s=c.system
    T=length(c.prices)
    names=nodes(s);ix=nodeindex(s)
    H=Matrix{Any}(undef,length(names),T)
    bounds=Matrix{Tuple{Float64,Float64}}(undef,length(names),T)
    reference_head=fill(NaN,length(names),T)
    for (i,r) in enumerate(s.reservoirs),t in 1:T
        a,b=reference.storage===nothing ? (r.v0,r.v0) :
            (reference.storage[i,t],reference.storage[i,t+1])
        da,db=_mean_head_partials(r,a,b)
        reference_head[i,t]=mean_head(r,a,b)
        H[i,t]=reference_head[i,t]+da*(V[i,t]-a)+db*(V[i,t+1]-b)
        lo1,hi1=storage_bounds(c,r,c.grid[t]);lo2,hi2=storage_bounds(c,r,c.grid[t+1])
        bounds[i,t]=extrema(reference_head[i,t]+da*(x-a)+db*(y-b)
            for x in (lo1,hi1),y in (lo2,hi2))
    end
    for (j,node) in enumerate(s.junctions),t in 1:T
        i=length(s.reservoirs)+j
        H[i,t]=@variable(m,lower_bound=node.hmin,upper_bound=node.hmax,
            base_name="proposal_head_$(i)_$(t)")
        bounds[i,t]=(node.hmin,node.hmax)
        reference.head!==nothing &&
            (reference_head[i,t]=reference.head[i,t])
    end
    for node in s.boundaries,t in 1:T
        H[ix[node.name],t]=node.head
        bounds[ix[node.name],t]=(node.head,node.head)
        reference_head[ix[node.name],t]=node.head
    end
    total_nominal=sum(clamp(g.qbest,g.qmin,g.qmax) for g in s.generators;init=0.0)
    for (i,e) in enumerate(s.tunnels),t in 1:T
        positions=tunnel_allowed_positions(c,e,t)
        all(iszero,positions) && continue
        src,dst=ix[e.source],ix[e.target]
        dh=H[src,t]-H[dst,t]
        dlo,dhi=bounds[src,t][1]-bounds[dst,t][2],bounds[src,t][2]-bounds[dst,t][1]
        onhead=_bounded_binary_product!(m,dh,gates.open[i,t],dlo,dhi;name="proposal_gate_head_$(i)_$(t)")
        cap=opinterval(c,e.name,:capacity,t,e.capacity)
        direct_nominal=sum(clamp(g.qbest,g.qmin,g.qmax) for g in s.generators
            if plantof(s,g).source in (e.source,e.target);init=0.0)
        nominal=direct_nominal>0 ? direct_nominal : total_nominal/max(1,length(s.tunnels))
        anchor=reference.tunnel_q===nothing ? nominal :
            abs(reference.tunnel_q[i,t])
        if (reference.tunnel_q===nothing || anchor==0) &&
            isfinite(reference_head[src,t]) && isfinite(reference_head[dst,t])
            # A passive reservoir/boundary link has a physical reference flow.
            # A unit's nominal flow can otherwise make its secant law infeasible
            # even when the pressure difference permits flow below capacity.
            position=length(positions)==1 ? only(positions) :
                e.opening>0 ? e.opening : last(positions)
            hydraulic=sqrt(abs(reference_head[src,t]-reference_head[dst,t])/
                tunnel_effective_resistance(e,position))
            hydraulic>0 && (anchor=hydraulic)
        end
        anchor=clamp(anchor>0 ? anchor : (nominal>0 ? nominal : cap),0.0,cap)
        if length(positions)==1
            loss=tunnel_effective_resistance(e,only(positions))*anchor*Q[i,t]
        else
            weights=gates.weights[(i,t)]
            loss=AffExpr(0.0)
            positive=findall(>(0),positions)
            for k in positive
                modeflow=length(positive)==1 ? Q[i,t] :
                    _bounded_binary_product!(m,Q[i,t],weights[k],e.discharge_river===nothing ? -cap : 0.0,
                        cap;name="proposal_gate_flow_$(i)_$(t)_$(k)")
                add_to_expression!(loss,tunnel_effective_resistance(e,positions[k])*anchor,modeflow)
            end
        end
        @constraint(m,onhead==loss)
    end
    P=Matrix{Any}(undef,length(s.generators),T)
    for (i,g) in enumerate(s.generators),t in 1:T
        plant=plantof(s,g)
        src,dst=ix[plant.source],ix[plant.target]
        receiver=H[dst,t]
        rlo,rhi=bounds[dst,t]
        if plant.outlet_head_floor!==nothing
            # A first-order approximation selects the active outlet-head branch.
            ref=reference.head===nothing ? (rlo+rhi)/2 : reference.head[dst,t]
            if ref<plant.outlet_head_floor
                receiver=plant.outlet_head_floor
                rlo=rhi=plant.outlet_head_floor
            end
        end
        qref=g.characteristic===nothing ? clamp(g.qbest,g.qmin,g.qmax) :
            pump_discharge(g.characteristic,reference_heads[i,t];extrapolation=:linear)
        tailref=reference.generator_q===nothing ?
            sum(clamp(z.qbest,z.qmin,z.qmax) for z in s.generators
            if z.plant==plant.name;init=0.0) :
            sum(reference.generator_q[j,t] for (j,z) in enumerate(s.generators) if z.plant==plant.name;init=0.0)
        tail=tailwater(plant,tailref)
        nethead=unit_power_sign(g)*(H[src,t]-receiver-tail)
        hlo,hhi=bounds[src,t][1]-rhi-tail,bounds[src,t][2]-rlo-tail
        g.mode==:pump && ((hlo,hhi)=(-hhi,-hlo))
        onhead=_bounded_binary_product!(m,nethead,u[i,t],hlo,hhi;name="proposal_unit_head_$(i)_$(t)")
        href=reference_heads[i,t]
        if g.characteristic!==nothing
            @constraint(m,gq[i,t]==qref*u[i,t]+_pump_slope(g.characteristic,href)*(onhead-href*u[i,t]))
        end
        head_gain=g.mode==:pump ? 0.00981*qref/(
            efficiency(g,qref,href)*g.motor_efficiency) :
            0.00981*qref*efficiency(g,qref,href)*generator_efficiency(g,power(g,qref,href))
        P[i,t]=alpha[i,t]*gq[i,t]+head_gain*(onhead-href*u[i,t])
    end
    P
end
"""Enclose passive flows over the permitted heads and gate positions.

Reservoir interval heads are volume averages, so endpoint storage domains
enclose them. No trajectory or proposed discrete state enters these bounds.
Flow-dependent delay laws do not change these local hydraulic domains.
"""
function _proposal_flow_bounds(c)
    s=c.system
    T=length(c.prices)
    ix=nodeindex(s)
    lower=zeros(length(nodes(s)),T)
    upper=similar(lower)
    for (i,r) in enumerate(s.reservoirs),t in 1:T
        a,b=soft_storage_bounds(c,r,c.grid[t])
        d,e=soft_storage_bounds(c,r,c.grid[t+1])
        lower[i,t],upper[i,t]=_global_level_range(r,min(a,d),max(b,e))
    end
    for (j,node) in enumerate(s.junctions)
        lower[length(s.reservoirs)+j,:].=node.hmin
        upper[length(s.reservoirs)+j,:].=node.hmax
    end
    for (j,node) in enumerate(s.boundaries)
        i=length(s.reservoirs)+length(s.junctions)+j
        lower[i,:].=node.head
        upper[i,:].=node.head
    end
    tunnel_lower=zeros(length(s.tunnels),T)
    tunnel_upper=similar(tunnel_lower)
    for (i,e) in enumerate(s.tunnels),t in 1:T
        a,b=ix[e.source],ix[e.target]
        lo,hi=lower[a,t]-upper[b,t],upper[a,t]-lower[b,t]
        tunnel_lower[i,t],tunnel_upper[i,t]=extrema(
            tunnel_flow(e,h,p) for p in tunnel_allowed_positions(c,e,t),h in (lo,hi))
    end
    river_lower=zeros(length(s.rivers),T)
    river_upper=similar(river_lower)
    for (i,r) in enumerate(s.rivers),t in 1:T
        cap=opinterval(c,r.name,:capacity,t,r.capacity)
        natural=opinterval(c,r.name,:inflow,t,r.inflow)
        river_lower[i,t],river_upper[i,t]=river_flow_lower(r,cap,natural),cap
        if is_channel(r)
            a,b=ix[r.source],ix[r.target]
            lo,hi=lower[a,t]-upper[b,t],upper[a,t]-lower[b,t]
            ref=r.delta_level_ref_up_flow_curve===nothing ? b : a
            river_lower[i,t],river_upper[i,t]=level_flow_bounds(
                channel_curve(r),lo,hi,lower[ref,t],upper[ref,t])
        elseif r.law!=:junction
            a=ix[r.source]
            gl=r.law==:weir ? 1.0 : opinterval(c,r.name,:gate_min,t,r.gate_min)
            gh=r.law==:weir ? 1.0 : opinterval(c,r.name,:gate_max,t,1.0)
            river_lower[i,t],river_upper[i,t]=extrema(
                natural+river_law_value(r,h,g;inflow=natural,capacity=cap)
                for h in (lower[a,t],upper[a,t]),g in (gl,gh))
        end
    end
    (;tunnel_lower,tunnel_upper,river_lower,river_upper)
end
