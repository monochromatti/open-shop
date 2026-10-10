# Instantaneous arrival laws for the same original release cohorts used by
# RiverRouting. Arrival rates are m³/s; volume conversion does not enter here.
function delay_cdf(c::RiverRouting.DelayCurve, x::Real)
    sum(
        w * clamp((x-a)/(b-a), 0.0, 1.0) for
        (a, b, w) in zip(c.edges[1:(end - 1)], c.edges[2:end], c.weights)
    )
end

"""Instantaneous arrival coefficients, cohort × reference curve.

Each distributed reference curve has a column. A single curve or deterministic
delay retains two identical columns for compatibility. A release cohort [a,b) contributes
q*(F(t-a)-F(t-b)). For an atom, `side` chooses the one-sided value at a jump.
"""
function point_coefficients(r::River, release_grid, time::Real; side::Symbol = :right)
    RiverRouting.check_grid(release_grid)
    isfinite(time) || throw(ArgumentError("arrival query time must be finite"))
    side in (:left, :right) || throw(ArgumentError("side must be :left or :right"))
    n = length(release_grid)-1
    if r.deterministic_delay !== nothing
        d=r.deterministic_delay
        isfinite(d) && d>=0 ||
            throw(ArgumentError("deterministic delay must be finite and nonnegative"))
        K=zeros(n, 2)
        for k in 1:n
            a, b=release_grid[k]+d, release_grid[k + 1]+d
            K[k, :] .= side==:right ? (a<=time<b) : (a<time<=b)
        end
        return K
    end
    RiverRouting.check_curves(r.curves)
    K=zeros(n, max(2, length(r.curves)))
    for k in 1:n, j in eachindex(r.curves)
        a, b=release_grid[k], release_grid[k + 1]
        K[k, j]=clamp(
            delay_cdf(r.curves[j], time-a)-delay_cdf(r.curves[j], time-b),
            0.0,
            1.0,
        )
    end
    length(r.curves)==1 && (K[:, 2].=K[:, 1])
    K
end

function cohort_point_arrival(r::River, grid, q, time; side = :right)
    RiverRouting.check_releases(grid, q)
    K=point_coefficients(r, grid, time; side)
    total=0.0
    for k in eachindex(q)
        if r.deterministic_delay!==nothing
            total+=q[k]*K[k, 1]
        else
            for (j, w) in RiverRouting.blend(r.curves, q[k])
                total+=q[k]*w*K[k, j]
            end
        end
    end
    total
end

"""Exact instantaneous arrivals from historical and scheduled cohorts."""
function point_arrival(r::River, current_grid, q, time::Real; side::Symbol = :right)
    cohort_point_arrival(r, current_grid, q, time; side) +
    cohort_point_arrival(r, r.history_grid, r.history_release, time; side)
end

"""All breakpoints of the piecewise-linear instantaneous arrival rate.

Original cohort edges and all delay-bin edges are retained. Thus evaluating
both sides of these knots is sufficient for an exact minimum for prescribed
piecewise-constant releases, including flow-dependent mixtures selected at
release time. This statement applies to a reach's supplied cohorts; a confluence
rebinned onto scheduling intervals remains an approximation to continuous
upstream mixing.
"""
function arrival_knots(r::River, grid)
    RiverRouting.check_grid(grid)
    RiverRouting.check_grid(r.history_grid)
    delays =
        r.deterministic_delay===nothing ? unique(vcat((c.edges for c in r.curves)...)) :
        [r.deterministic_delay]
    isempty(delays) && throw(ArgumentError("a delay distribution is required"))
    all(isfinite, delays) && all(delays .>= 0) ||
        throw(ArgumentError("delays must be finite and nonnegative"))
    lo, hi=first(grid), last(grid)
    sort!(
        unique!(
            vcat(
                [lo, hi],
                [t+d for t in vcat(r.history_grid, grid) for d in delays if lo<=t+d<=hi],
            ),
        ),
    )
end

"""Minimum instantaneous arrival over the closed scheduling horizon.

The horizon start uses its right limit and horizon end its left limit, avoiding
an artificial post-horizon zero caused by the finite release schedule.
"""
function minimum_arrival(r::River, grid, q)
    best=Inf
    at=first(grid)
    bestside=:right
    for t in arrival_knots(r, grid)
        sides=t==first(grid) ? (:right,) : t==last(grid) ? (:left,) : (:left, :right)
        for side in sides
            value=point_arrival(r, grid, q, t; side)
            if value<best
                best=value
                at=t
                bestside=side
            end
        end
    end
    Dict("value"=>best, "time"=>at, "side"=>bestside)
end

"""Clip original release cohorts to a restart time, retaining their flows.

This is restart history, not arrival bucket averages. The returned edges are
absolute times; the caller may translate all edges by the same restart origin.
A restart preceding every release returns one edge and no releases.
"""
function release_history(grid, q, time::Real)
    RiverRouting.check_grid(grid)
    RiverRouting.check_releases(grid, q)
    isfinite(time) || throw(ArgumentError("restart time must be finite"))
    time<=first(grid) && return (grid = [Float64(time)], release = Float64[])
    k=min(searchsortedlast(grid, time), length(q))
    edges=Float64.(grid[1:(k + 1)])
    edges[end]=min(Float64(time), last(grid))
    if edges[end]==edges[end - 1]
        pop!(edges)
        k-=1
    end
    (grid = edges, release = Float64.(q[1:k]))
end
