# Conservative interval domains for the exact scheduling-grid equations.
# Release coefficients are bounded independently: no assumed commitment,
# tunnel direction, or local trajectory is used to remove feasible states.
function _global_quadratic_range(a, b, lo, hi)
    vals=(a*lo+b*lo^2, a*hi+b*hi^2)
    lower, upper=extrema(vals)
    if b!=0
        x=-a/(2b)
        if lo<x<hi
            y=a*x+b*x^2
            lower=min(lower, y)
            upper=max(upper, y)
        end
    end
    lower, upper
end
function _global_level_range(r, lo, hi)
    if r.level_curve===nothing
        a, b=_global_quadratic_range(r.slope, r.curvature, lo, hi)
        return r.z0+a, r.z0+b
    end
    # TableCurve is piecewise linear, including its linear extrapolation.
    vals=[head(r, lo), head(r, hi)]
    for x in r.level_curve.x
        lo<x<hi && push!(vals, head(r, x))
    end
    extrema(vals)
end
function _global_capacity_bounds(c, nd, rd)
    s=c.system
    T=length(c.prices)
    R=length(s.reservoirs)
    dt=diff(c.grid)
    release_cap=[
        opinterval(c, r.name, :capacity, t, r.capacity) for r in s.rivers, t in 1:T
    ]
    al=zeros(length(s.rivers), T)
    au=similar(al)
    for (i, r) in enumerate(s.rivers), t in 1:T
        if nd!==nothing
            lo=hi=nd.arrival_history[i, t]
            for (d, k, v) in nd.arrival_terms[i, t]
                z=v*release_cap[d, k]
                lo+=min(0.0, z)
                hi+=max(0.0, z)
            end
        else
            lo=hi=rd["history_arrival"][i, t]
            B=rd["B"][i]
            for k in 1:T
                a=0.0036*dt[k]*B[t, k, 1]
                b=0.0036*dt[k]*(B[t, k, 2]-B[t, k, 1])/r.capacity
                l, h=_global_quadratic_range(a, b, 0.0, release_cap[i, k])
                lo+=l
                hi+=h
            end
        end
        al[i, t]=lo
        au[i, t]=hi
    end
    lower=zeros(R, T+1)
    upper=similar(lower)
    delta_lo=zeros(R, T)
    delta_hi=similar(delta_lo)
    plants=[plantof(s, g) for g in s.generators]
    for (i, r) in enumerate(s.reservoirs)
        tunnels=[e for e in s.tunnels if e.source==r.name || e.target==r.name]
        gen_out=findall(p->p.source==r.name, plants)
        gen_in=findall(p->p.target==r.name, plants)
        river_out=findall(reach->reach.source==r.name, s.rivers)
        river_in=findall(reach->reach.target==r.name, s.rivers)
        for t in 1:(T + 1)
            lower[i, t], upper[i, t]=storage_bounds(c, r, c.grid[t])
        end
        lower[i, 1]=upper[i, 1]=r.v0
        for t in 1:T
            lo=hi=opinterval(c, r.name, :inflow, t, r.inflow)
            for e in tunnels
                cap=opinterval(c, e.name, :opening, t, e.opening)==0 ? 0.0 :
                    opinterval(c, e.name, :capacity, t, e.capacity)
                lo-=cap
                hi+=cap
            end
            for j in gen_out
                g=s.generators[j]
                lo-=opinterval(c, g.name, :qmax, t, g.qmax)
            end
            for j in gen_in
                g=s.generators[j]
                hi+=opinterval(c, g.name, :qmax, t, g.qmax)
            end
            for j in river_out
                lo-=release_cap[j, t]
            end
            for j in river_in
                lo+=al[j, t]/(0.0036*dt[t])
                hi+=au[j, t]/(0.0036*dt[t])
            end
            # Outward slack protects the enclosure against floating arithmetic.
            delta_lo[i, t]=0.0036*dt[t]*lo-1e-7
            delta_hi[i, t]=0.0036*dt[t]*hi+1e-7
        end
        for t in 1:T
            lower[i, t + 1]=max(lower[i, t + 1], lower[i, t]+delta_lo[i, t])
            upper[i, t + 1]=min(upper[i, t + 1], upper[i, t]+delta_hi[i, t])
        end
        for t in T:-1:2
            lower[i, t]=max(lower[i, t], lower[i, t + 1]-delta_hi[i, t])
            upper[i, t]=min(upper[i, t], upper[i, t + 1]-delta_lo[i, t])
        end
        any(lower[i, :] .> upper[i, :]) &&
            throw(ArgumentError("unreachable storage limits for $(r.name)"))
    end
    hlo=zeros(R, T)
    hhi=similar(hlo)
    for (i, r) in enumerate(s.reservoirs), t in 1:T
        hlo[i, t], hhi[i, t]=_global_level_range(
            r,
            (lower[i, t]+lower[i, t + 1])/2,
            (upper[i, t]+upper[i, t + 1])/2,
        )
        hlo[i, t]-=1e-7
        hhi[i, t]+=1e-7
    end
    (; lower, upper, hlo, hhi, arrival_lower = al, arrival_upper = au)
end

# Each intersection is an enclosure, not a numerical feasibility projection.
# Absolute outward slack is kept in physical units; existing declared bounds
# (including exact zero-flow bounds) are never relaxed by that slack.
const _GLOBAL_DOMAIN_SLACK = 1e-7

function _global_intersect!(lower, upper, i, t, lo, hi; slack = _GLOBAL_DOMAIN_SLACK)
    a=max(lower[i, t], lo-slack)
    b=min(upper[i, t], hi+slack)
    a<=b || throw(ArgumentError("inconsistent reachable interval at ($i,$t)"))
    change=max(a-lower[i, t], upper[i, t]-b)
    lower[i, t], upper[i, t]=a, b
    change
end

function _global_head_intersect!(bounds, key, lo, hi)
    oldlo, oldhi=bounds[key]
    a=max(oldlo, lo-_GLOBAL_DOMAIN_SLACK)
    b=min(oldhi, hi+_GLOBAL_DOMAIN_SLACK)
    a<=b || throw(ArgumentError("inconsistent hydraulic head bounds for $key"))
    bounds[key]=(a, b)
    max(a-oldlo, oldhi-b)
end

_global_signed_root(x) = x==0 ? 0.0 : copysign(sqrt(abs(x)), x)
_global_signed_loss(x) = x*abs(x)
_global_signed_interval(sign, lo, hi) = sign>0 ? (lo, hi) : (-hi, -lo)

function _global_transport_bounds!(al, au, c, nd, rd, rl, ru)
    dt=diff(c.grid)
    T=length(c.prices)
    change=0.0
    for (i, r) in enumerate(c.system.rivers), t in 1:T
        lo=hi=nd===nothing ? rd["history_arrival"][i, t] : nd.arrival_history[i, t]
        if nd!==nothing
            for (d, k, coefficient) in nd.arrival_terms[i, t]
                a=coefficient*rl[d, k]
                b=coefficient*ru[d, k]
                lo+=min(a, b)
                hi+=max(a, b)
            end
        else
            B=rd["B"][i]
            for k in 1:T
                a=0.0036*dt[k]*B[t, k, 1]
                b=0.0036*dt[k]*(B[t, k, 2]-B[t, k, 1])/r.capacity
                l, h=_global_quadratic_range(a, b, rl[i, k], ru[i, k])
                lo+=l
                hi+=h
            end
        end
        change=max(change, _global_intersect!(al, au, i, t, max(0.0, lo), hi))
    end
    change
end

function _global_river_law_bounds!(rl, ru, c, heads, al, au, confluences)
    dt=diff(c.grid)
    T=length(c.prices)
    change=0.0
    for (i, r) in enumerate(c.system.rivers), t in 1:T
        cap=opinterval(c, r.name, :capacity, t, r.capacity)
        gate_lo=opinterval(c, r.name, :gate_min, t, r.gate_min)
        gate_hi=opinterval(c, r.name, :gate_max, t, 1.0)
        if r.law==:junction
            incoming=confluences[r.source]
            lo=sum(al[j, t] for j in incoming)/(0.0036*dt[t])
            hi=sum(au[j, t] for j in incoming)/(0.0036*dt[t])
        elseif r.law==:controlled && r.discharge_curve===nothing
            lo, hi=cap*gate_lo, cap*gate_hi
        else
            head_lo, head_hi=heads[(r.source, t)]
            if r.discharge_curve!==nothing
                # Input validation requires this physical law to be monotone.
                flow_lo=table_value(r.discharge_curve, head_lo; extrapolation = :linear)
                flow_hi=table_value(r.discharge_curve, head_hi; extrapolation = :linear)
            else
                wet_lo=max(0.0, head_lo-r.crest)
                wet_hi=max(0.0, head_hi-r.crest)
                exponent=r.law==:orifice ? 0.5 : 1.5
                flow_lo=r.coefficient*wet_lo^exponent
                flow_hi=r.coefficient*wet_hi^exponent
            end
            r.law==:weir && (gate_lo=gate_hi=1.0)
            lo=max(0.0, flow_lo)*gate_lo
            hi=max(0.0, flow_hi)*gate_hi
        end
        change=max(change, _global_intersect!(rl, ru, i, t, lo, hi))
    end
    change
end

function _global_tunnel_bounds!(ql, qu, c, heads, junction_incidence, gl, gu)
    T=length(c.prices)
    change=0.0
    for (i, e) in enumerate(c.system.tunnels), t in 1:T
        opening=opinterval(c, e.name, :opening, t, e.opening)
        opening==0 && continue # A closed tunnel imposes no head relationship.
        src=(e.source, t)
        dst=(e.target, t)
        alo, ahi=heads[src]
        blo, bhi=heads[dst]
        lo=_global_signed_root(opening*(alo-bhi)/e.resistance)
        hi=_global_signed_root(opening*(ahi-blo)/e.resistance)
        change=max(change, _global_intersect!(ql, qu, i, t, lo, hi))
        loss_lo=e.resistance*_global_signed_loss(ql[i, t])/opening
        loss_hi=e.resistance*_global_signed_loss(qu[i, t])/opening
        change=max(change, _global_head_intersect!(heads, src, blo+loss_lo, bhi+loss_hi))
        alo, ahi=heads[src]
        change=max(change, _global_head_intersect!(heads, dst, alo-loss_hi, ahi-loss_lo))
    end
    # Hydraulic junctions have no storage or river releases. The signed sum of
    # their tunnel flows and generating discharges is exactly zero. In addition
    # to interval sums, a one-sided cone proves a zero bound without erasing it
    # with rounding slack. This detects intake chains without assuming that
    # arbitrary reservoir-to-reservoir tunnels run forward.
    for (tunnels, generators) in junction_incidence, t in 1:T
        for (i, sign) in tunnels
            lo=hi=0.0
            all_nonnegative=true
            all_nonpositive=true
            for (j, coefficient) in tunnels
                j==i && continue
                a, b=_global_signed_interval(coefficient, ql[j, t], qu[j, t])
                lo+=a
                hi+=b
                all_nonnegative &= a>=0
                all_nonpositive &= b<=0
            end
            for (j, coefficient) in generators
                a, b=_global_signed_interval(coefficient, gl[j, t], gu[j, t])
                lo+=a
                hi+=b
                all_nonnegative &= a>=0
                all_nonpositive &= b<=0
            end
            a, b=sign>0 ? (-hi, -lo) : (lo, hi)
            change=max(change, _global_intersect!(ql, qu, i, t, a, b))
            nonnegative=sign>0 ? all_nonpositive : all_nonnegative
            nonpositive=sign>0 ? all_nonnegative : all_nonpositive
            if nonnegative && ql[i, t]<0
                change=max(change, -ql[i, t])
                ql[i, t]=0.0
            end
            if nonpositive && qu[i, t]>0
                change=max(change, qu[i, t])
                qu[i, t]=0.0
            end
            ql[i, t]<=qu[i, t] ||
                throw(ArgumentError("inconsistent tunnel continuity for ($i,$t)"))
        end
    end
    change
end

function _global_storage_bounds!(lower, upper, c, ql, qu, gl, gu, rl, ru, al, au, incidence)
    dt=diff(c.grid)
    T=length(c.prices)
    change=0.0
    for (i, r) in enumerate(c.system.reservoirs)
        tunnels, generators, outgoing, incoming=incidence[i]
        dl=zeros(T)
        du=similar(dl)
        for t in 1:T
            lo=hi=opinterval(c, r.name, :inflow, t, r.inflow)
            for (j, sign) in tunnels
                a, b=_global_signed_interval(sign, ql[j, t], qu[j, t])
                lo+=a
                hi+=b
            end
            for (j, sign) in generators
                a, b=_global_signed_interval(sign, gl[j, t], gu[j, t])
                lo+=a
                hi+=b
            end
            for j in outgoing
                lo-=ru[j, t]
                hi-=rl[j, t]
            end
            for j in incoming
                lo+=al[j, t]/(0.0036*dt[t])
                hi+=au[j, t]/(0.0036*dt[t])
            end
            dl[t]=0.0036*dt[t]*lo-_GLOBAL_DOMAIN_SLACK
            du[t]=0.0036*dt[t]*hi+_GLOBAL_DOMAIN_SLACK
        end
        for t in 1:T
            change=max(change, _global_intersect!(lower, upper, i, t + 1,
                lower[i, t]+dl[t], upper[i, t]+du[t]; slack = 0.0))
        end
        for t in T:-1:2
            change=max(change, _global_intersect!(lower, upper, i, t,
                lower[i, t + 1]-du[t], upper[i, t + 1]-dl[t]; slack = 0.0))
        end
    end
    change
end

"""Enclose the full feasible set using conservation and monotone hydraulic laws.

`tightened=false` reproduces the original capacity-only domains. Tightened
domains use no commitment candidate, tunnel-direction guess or local trajectory.
Finite propagation passes may leave slack; they cannot remove feasible states.
"""
function _global_reachable_bounds(c, nd, rd; tightened = true)
    old=_global_capacity_bounds(c, nd, rd)
    tightened || return old
    s=c.system
    T=length(c.prices)
    R=length(s.reservoirs)
    E=length(s.tunnels)
    G=length(s.generators)
    D=length(s.rivers)
    lower, upper=copy(old.lower), copy(old.upper)
    ql=zeros(E, T)
    qu=similar(ql)
    gl=zeros(G, T)
    gu=similar(gl)
    rl=zeros(D, T)
    ru=similar(rl)
    for (i, e) in enumerate(s.tunnels), t in 1:T
        cap=opinterval(c, e.name, :opening, t, e.opening)==0 ? 0.0 :
            opinterval(c, e.name, :capacity, t, e.capacity)
        ql[i, t], qu[i, t]=-cap, cap
    end
    for (i, g) in enumerate(s.generators), t in 1:T
        forced=opinterval(c, g.name, :forced_on, t, -1.0)
        gl[i, t]=forced==1 ? opinterval(c, g.name, :qmin, t, g.qmin) : 0.0
        gu[i, t]=forced==0 ? 0.0 : opinterval(c, g.name, :qmax, t, g.qmax)
    end
    for (i, r) in enumerate(s.rivers), t in 1:T
        penalty=opinterval(c, r.name, :release_penalty, t, 0.0)
        rl[i, t]=penalty>0 ? 0.0 : opinterval(c, r.name, :min_release, t, 0.0)
        ru[i, t]=opinterval(c, r.name, :capacity, t, r.capacity)
        rl[i, t]<=ru[i, t] || throw(ArgumentError("infeasible river release limits"))
    end
    al, au=copy(old.arrival_lower), copy(old.arrival_upper)
    heads=Dict{Tuple{Symbol,Int},Tuple{Float64,Float64}}()
    for (i, r) in enumerate(s.reservoirs), t in 1:T
        heads[(r.name, t)]=(old.hlo[i, t], old.hhi[i, t])
    end
    for j in s.junctions, t in 1:T
        heads[(j.name, t)]=(j.hmin, j.hmax)
    end
    for b in s.boundaries, t in 1:T
        heads[(b.name, t)]=(b.head, b.head)
    end
    plants=[plantof(s, g) for g in s.generators]
    function hydraulic_incidence(name)
        tunnels=[(i, Int(e.target==name)-Int(e.source==name))
            for (i, e) in enumerate(s.tunnels) if e.target==name || e.source==name]
        generators=[(i, Int(p.target==name)-Int(p.source==name))
            for (i, p) in enumerate(plants) if p.target!=p.source &&
            (p.target==name || p.source==name)]
        tunnels, generators
    end
    junction_incidence=[hydraulic_incidence(j.name) for j in s.junctions]
    incidence=[(hydraulic_incidence(r.name)...,
        findall(reach->reach.source==r.name, s.rivers),
        findall(reach->reach.target==r.name, s.rivers)) for r in s.reservoirs]
    confluences=Dict(j.name=>findall(r->r.target==j.name, s.rivers)
        for j in s.river_junctions)
    passes=0
    # Bounded preprocessing, even for a loop with asymptotic propagation.
    for pass in 1:24
        passes=pass
        change=_global_tunnel_bounds!(ql, qu, c, heads, junction_incidence, gl, gu)
        change=max(change, _global_river_law_bounds!(rl, ru, c, heads, al, au, confluences))
        change=max(change, _global_transport_bounds!(al, au, c, nd, rd, rl, ru))
        change=max(change, _global_storage_bounds!(lower, upper, c,
            ql, qu, gl, gu, rl, ru, al, au, incidence))
        for (i, r) in enumerate(s.reservoirs), t in 1:T
            a, b=_global_level_range(r,
                (lower[i, t]+lower[i, t + 1])/2,
                (upper[i, t]+upper[i, t + 1])/2)
            change=max(change, _global_head_intersect!(heads, (r.name, t), a, b))
        end
        change<=1e-8 && break
    end
    hlo=zeros(R, T)
    hhi=similar(hlo)
    for (i, r) in enumerate(s.reservoirs), t in 1:T
        hlo[i, t], hhi[i, t]=heads[(r.name, t)]
    end
    (; lower, upper, hlo, hhi, arrival_lower = al, arrival_upper = au,
        tunnel_lower = ql, tunnel_upper = qu, generator_lower = gl,
        generator_upper = gu, release_lower = rl, release_upper = ru,
        node_head_bounds = heads, tightening_passes = passes)
end
