# A gate proposal must distinguish the heads supplied by alternative intakes.
# These affine head/storage and secant loss equations are a scheduling surrogate,
# never a physical audit or a bound on the nonlinear problem. They are used only
# when the case contains selectable gates; the existing prescribed-gate proposal
# keeps its original construction.
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
        vref=reference.storage===nothing ? r.v0 :
            (reference.storage[i,t]+reference.storage[i,t+1])/2
        slope=r.level_curve===nothing ? r.slope+2r.curvature*vref :
            table_slope(r.level_curve,vref;extrapolation=:linear)
        intercept=head(r,vref)-slope*vref
        reference_head[i,t]=head(r,vref)
        H[i,t]=intercept+slope*(V[i,t]+V[i,t+1])/2
        lo1,hi1=storage_bounds(c,r,c.grid[t]);lo2,hi2=storage_bounds(c,r,c.grid[t+1])
        bounds[i,t]=extrema((intercept+slope*(lo1+lo2)/2,intercept+slope*(hi1+hi2)/2))
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
        qref=clamp(g.qbest,g.qmin,g.qmax)
        tailref=reference.generator_q===nothing ?
            sum(clamp(z.qbest,z.qmin,z.qmax) for z in s.generators
            if z.plant==plant.name;init=0.0) :
            sum(reference.generator_q[j,t] for (j,z) in enumerate(s.generators) if z.plant==plant.name;init=0.0)
        tail=tailwater(plant,tailref)
        nethead=H[src,t]-receiver-tail
        hlo,hhi=bounds[src,t][1]-rhi-tail,bounds[src,t][2]-rlo-tail
        onhead=_bounded_binary_product!(m,nethead,u[i,t],hlo,hhi;name="proposal_unit_head_$(i)_$(t)")
        href=reference_heads[i,t]
        head_gain=0.00981*qref*efficiency(g,qref,href)*
            generator_efficiency(g,power(g,qref,href))
        P[i,t]=alpha[i,t]*gq[i,t]+head_gain*(onhead-href*u[i,t])
    end
    P
end
