module RiverRouting

export DelayCurve,
    coefficients,
    fixed_coefficients,
    route,
    remaining_volume,
    mean_flows,
    merge_arrivals,
    chain_route

const FLOW_HOUR_TO_MM3 = 0.0036

"""Finite uniform delay bins, in hours. Weights are bin probabilities.

The SHOP representation with one weight per edge is accepted only when its last
weight is zero. Internally there is one weight per interval. Reference flow is
in m³/s. This module deliberately does not silently normalize invalid inputs.
"""
struct DelayCurve
    reference_flow::Float64
    edges::Vector{Float64}
    weights::Vector{Float64}
    function DelayCurve(
        reference_flow::Real,
        edges::AbstractVector{<:Real},
        weights::AbstractVector{<:Real},
    )
        r, e, w = Float64(reference_flow), Float64.(edges), Float64.(weights)
        isfinite(r) && r >= 0 ||
            throw(ArgumentError("reference flow must be finite and nonnegative"))
        length(e) >= 2 && all(isfinite, e) && all(diff(e) .> 0) && first(e) >= 0 || throw(
            ArgumentError(
                "delay edges must be finite, nonnegative, and strictly increasing",
            ),
        )
        if length(w) == length(e)
            last(w) == 0 || throw(ArgumentError("SHOP trailing edge weight must be zero"))
            pop!(w)
        end
        length(w) == length(e)-1 && all(isfinite, w) && all(w .>= 0) ||
            throw(ArgumentError("one finite nonnegative weight is required per delay bin"))
        isapprox(sum(w), 1.0; atol = 1e-12, rtol = 1e-12) ||
            throw(ArgumentError("delay weights must sum to one"))
        new(r, e, w)
    end
end

function check_grid(grid)
    length(grid) >= 2 && all(isfinite, grid) && all(diff(grid) .> 0) || throw(
        ArgumentError("time grids need at least two finite, strictly increasing edges"),
    )
    nothing
end

function check_curves(curves)
    !isempty(curves) || throw(ArgumentError("at least one delay curve is required"))
    all(diff([c.reference_flow for c in curves]) .> 0) || throw(ArgumentError("reference flows must be strictly increasing"))
    nothing
end

# Stable difference of CDF antiderivatives G(hi)-G(lo), where
# G(x)=sum(w*((x-a)_+²-(x-b)_+²)/(2*(b-a))). Integrating ramp and
# saturated portions separately avoids catastrophic cancellation far downstream.
function integral_cdf(c::DelayCurve, lo::Real, hi::Real)
    hi >= lo || throw(ArgumentError("reversed integration interval"))
    total = 0.0
    for k in eachindex(c.weights)
        a, b, w = c.edges[k], c.edges[k + 1], c.weights[k]
        l, h = max(lo, a), min(hi, b)
        ramp = h > l ? (h-l)*((l-a)+(h-a))/(2*(b-a)) : 0.0
        tail = max(0.0, hi-max(lo, b))
        total += w*(ramp+tail)
    end
    total
end

"""Exact transfer fractions [arrival interval, release interval, reference curve].

Releases are constant within each release interval. Arrival/release time grids
need not coincide, and may include negative times for historical cohorts.
"""
function coefficients(arrival_grid, release_grid, curves)
    check_grid(arrival_grid)
    check_grid(release_grid)
    check_curves(curves)
    B = zeros(length(arrival_grid)-1, length(release_grid)-1, length(curves))
    for l in eachindex(curves),
        k in 1:(length(release_grid) - 1),
        t in 1:(length(arrival_grid) - 1)

        a, b = arrival_grid[t], arrival_grid[t + 1]
        s, e = release_grid[k], release_grid[k + 1]
        # Impossible/full-support cases also protect distant-time subtraction.
        if b <= s+first(curves[l].edges) || a >= e+last(curves[l].edges)
            B[t, k, l] = 0.0
        else
            x =
                (integral_cdf(curves[l], b-e, b-s) - integral_cdf(curves[l], a-e, a-s))/(
                    e-s
                )
            -1e-10 <= x <= 1+1e-10 || error("invalid transfer fraction: $x")
            B[t, k, l] = clamp(x, 0.0, 1.0)
        end
    end
    B
end

"""Exact uniform-cohort fractions for a deterministic nonnegative delay."""
function fixed_coefficients(arrival_grid, release_grid, delay::Real)
    check_grid(arrival_grid)
    check_grid(release_grid)
    isfinite(delay) && delay >= 0 ||
        throw(ArgumentError("delay must be finite and nonnegative"))
    [
        max(
            0.0,
            min(arrival_grid[t + 1], release_grid[k + 1]+delay) -
            max(arrival_grid[t], release_grid[k]+delay),
        )/(release_grid[k + 1]-release_grid[k]) for
        t in 1:(length(arrival_grid) - 1), k in 1:(length(release_grid) - 1)
    ]
end

function check_releases(release_grid, releases)
    length(releases) == length(release_grid)-1 &&
    all(isfinite, releases) &&
    all(releases .>= 0) ||
        throw(ArgumentError("one finite nonnegative mean release is required per interval"))
end

# One curve is flow-independent. Multiple curves use a convex linear blend of
# distributions at the contemporaneous release flow, with no extrapolation.
function blend(curves, q)
    length(curves) == 1 && return [(1, 1.0)]
    refs = [c.reference_flow for c in curves]
    first(refs) <= q <= last(refs) ||
        throw(ArgumentError("release flow outside reference range"))
    q == last(refs) && return [(length(refs), 1.0)]
    j = searchsortedlast(refs, q)
    w = (q-refs[j])/(refs[j + 1]-refs[j])
    [(j, 1-w), (j+1, w)]
end

"""Route original release cohorts; returns arrival volumes in Mm³.

Flow-dependent distributions are selected at release time. The entire original
cohort must be retained for restart on a finer grid: arrival bucket averages
generally cannot reconstruct its within-bucket arrival shape.
"""
function route(arrival_grid, release_grid, releases, curves)
    check_releases(release_grid, releases)
    B = coefficients(arrival_grid, release_grid, curves)
    volumes = zeros(length(arrival_grid)-1)
    for k in eachindex(releases), (l, w) in blend(curves, releases[k])
        volumes .+=
            (FLOW_HOUR_TO_MM3*(release_grid[k + 1]-release_grid[k])*releases[k]*w) .*
            B[:, k, l]
    end
    volumes
end

"""Already released but not yet arrived volume (Mm³) at `time`.

Future, unreleased parts of cohorts are excluded. This is physical water in
transit, rather than the total amount scheduled to arrive after `time`.
"""
function remaining_volume(time::Real, release_grid, releases, curves)
    isfinite(time) || throw(ArgumentError("time must be finite"))
    check_grid(release_grid)
    check_releases(release_grid, releases)
    check_curves(curves)
    volume = 0.0
    for k in eachindex(releases)
        mix = blend(curves, releases[k])
        a, b = release_grid[k], min(time, release_grid[k + 1])
        b <= a && continue
        for (l, w) in mix
            arrived_time = integral_cdf(curves[l], time-b, time-a)
            volume += FLOW_HOUR_TO_MM3*releases[k]*w*max(0.0, (b-a)-arrived_time)
        end
    end
    volume
end

"""Convert interval volumes in Mm³ to interval-average flows in m³/s."""
function mean_flows(grid, volumes)
    check_grid(grid)
    length(volumes) == length(grid)-1 && all(isfinite, volumes) && all(volumes .>= 0) ||
        throw(ArgumentError("one finite nonnegative volume is required per interval"))
    volumes ./ (FLOW_HOUR_TO_MM3 .* diff(grid))
end

"""Merge disjoint incoming water streams on the same grid, adding each once.

The caller owns topology: passing the same physical stream twice is invalid.
"""
function merge_arrivals(arrivals::AbstractVector...)
    !isempty(arrivals) || throw(ArgumentError("at least one arrival stream is required"))
    n = length(first(arrivals))
    all(v -> length(v)==n && all(isfinite, v) && all(v .>= 0), arrivals) ||
        throw(ArgumentError("arrival streams must share a grid and be nonnegative"))
    reduce(+, arrivals)
end

"""Route a downstream reach after interval-average mixing at a junction.

This is a deliberate discretization: incoming volume is redistributed uniformly
within each junction interval before entering the next reach. It conserves mass
but can spread pulses or move arrivals within an interval. Refine the junction
grid to study this error; never use this operation to reconstruct restart history.
"""
function chain_route(arrival_grid, junction_grid, incoming_volumes, curves)
    route(arrival_grid, junction_grid, mean_flows(junction_grid, incoming_volumes), curves)
end

end # module
