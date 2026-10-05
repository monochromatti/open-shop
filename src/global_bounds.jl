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
function _global_reachable_bounds(c, nd, rd)
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
