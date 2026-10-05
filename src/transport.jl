# Exact deterministic-delay routing retains rectangular release cohorts through
# every confluence. Labels identify source decision variables; label (0,0) is
# historical water. No intermediate scheduling-bin averaging is performed.
struct _TransportPulse
    start::Float64
    stop::Float64
    coefficient::Float64
    source::Int
    period::Int
end

const _TransportTerm=Tuple{Int,Int,Float64}

# A borrowed, read-only range into the compiled coefficient buffer. Keeping
# this concrete avoids copying a coefficient vector for every point query.
struct _TransportTerms <: AbstractVector{_TransportTerm}
    data::Vector{_TransportTerm}
    first::Int
    last::Int
end
Base.size(x::_TransportTerms) = (max(0, x.last-x.first+1),)
Base.IndexStyle(::Type{_TransportTerms}) = IndexLinear()
function Base.getindex(x::_TransportTerms, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds x.data[x.first + i - 1]
end

struct _TransportPointIndex
    times::Vector{Float64}
    offsets::Vector{Int}
    terms::Vector{_TransportTerm}
    history::Vector{Float64}
end

struct _TransportSignature
    reservoirs::Vector{Symbol}
    names::Vector{Symbol}
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    delays::Vector{Float64}
    history_grid::Vector{Vector{Float64}}
    history_release::Vector{Vector{Float64}}
end

"""Owned transport coefficients for one source grid and evaluation grid.

Use concrete fields in numerical kernels. String indexing preserves the v0.2
read interface. Point expressions borrow read-only coefficient ranges. Reuse
is explicit and checked against transport input snapshots; there is no global
cache keyed by a mutable case. Treat the compiled buffers as read-only.
"""
struct CompiledTransport
    arrival_terms::Matrix{Vector{_TransportTerm}}
    arrival_history::Matrix{Float64}
    release_terms::Matrix{Vector{_TransportTerm}}
    release_history::Matrix{Float64}
    terminal_terms::Vector{Vector{_TransportTerm}}
    terminal_history::Vector{Float64}
    initial_transit::Vector{Float64}
    release_pulses::Vector{Vector{_TransportPulse}}
    arrival_pulses::Vector{Vector{_TransportPulse}}
    grid::Vector{Float64}
    source_grid::Vector{Float64}
    order::Vector{Int}
    arrival_index::Vector{_TransportPointIndex}
    release_index::Vector{_TransportPointIndex}
    signature::_TransportSignature
end
const _TRANSPORT_KEYS=(
    "arrival_terms",
    "arrival_history",
    "release_terms",
    "release_history",
    "terminal_terms",
    "terminal_history",
    "initial_transit",
    "release_pulses",
    "arrival_pulses",
    "grid",
    "source_grid",
    "order",
)
Base.keys(::CompiledTransport) = _TRANSPORT_KEYS
Base.haskey(::CompiledTransport, key) = key in _TRANSPORT_KEYS
function Base.getindex(data::CompiledTransport, key::AbstractString)
    haskey(data, key) || throw(KeyError(key))
    getproperty(data, Symbol(key))
end

function _transport_signature(c)
    s=c.system
    _TransportSignature(
        [r.name for r in s.reservoirs],
        [r.name for r in s.rivers],
        [r.source for r in s.rivers],
        [r.target for r in s.rivers],
        Float64[r.deterministic_delay for r in s.rivers],
        [copy(r.history_grid) for r in s.rivers],
        [copy(r.history_release) for r in s.rivers],
    )
end

function _transport_data(c, data::CompiledTransport; arrival_grid = c.grid)
    # Reuse must retain the topology/history/capacity checks of fresh compile.
    # Capacity is a live restriction, not a coefficient cache dependency.
    river_order(c)
    s=c.system
    sig=data.signature
    data.source_grid==c.grid && data.grid==arrival_grid ||
        throw(ArgumentError("compiled transport grid differs from case"))
    length(sig.reservoirs)==length(s.reservoirs) && length(sig.names)==length(s.rivers) ||
        throw(ArgumentError("compiled transport topology differs from case"))
    for (i, r) in enumerate(s.reservoirs)
        sig.reservoirs[i]==r.name ||
            throw(ArgumentError("compiled transport reservoir order differs from case"))
    end
    for (i, r) in enumerate(s.rivers)
        sig.names[i]==r.name &&
        sig.sources[i]==r.source &&
        sig.targets[i]==r.target &&
        sig.delays[i]==r.deterministic_delay &&
        sig.history_grid[i]==r.history_grid &&
        sig.history_release[i]==r.history_release ||
            throw(ArgumentError("compiled transport inputs changed for $(r.name)"))
    end
    data
end
_transport_data(c, ::Nothing; arrival_grid = c.grid) =
    deterministic_network_data(c; arrival_grid)

# Group exact event times, then recompute the active expression in original
# pulse order. This preserves history summation (including large cancellation)
# rather than accumulating floating-point start/stop updates indefinitely.
function _transport_point_index(pulses::Vector{_TransportPulse})
    events=Tuple{Float64,Int,Bool}[]
    sizehint!(events, 2length(pulses))
    for (i, p) in enumerate(pulses)
        push!(events, (p.start, i, true))
        push!(events, (p.stop, i, false))
    end
    sort!(events; by = first)
    times=Float64[]
    offsets=Int[1]
    terms=_TransportTerm[]
    history=Float64[]
    sizehint!(times, length(events))
    sizehint!(offsets, length(events)+1)
    sizehint!(history, length(events))
    active=BitSet()
    coefficients=Dict{Tuple{Int,Int},Float64}()
    j=1
    while j<=length(events)
        time=events[j][1]
        while j<=length(events) && events[j][1]==time
            _, i, starts=events[j]
            starts ? push!(active, i) : delete!(active, i)
            j+=1
        end
        empty!(coefficients)
        historical=0.0
        for i in active
            p=pulses[i]
            if p.source==0
                historical+=p.coefficient
            else
                key=(p.source, p.period)
                coefficients[key]=get(coefficients, key, 0.0)+p.coefficient
            end
        end
        for (key, value) in sort!(collect(coefficients); by = first)
            push!(terms, (key[1], key[2], value))
        end
        push!(times, time)
        push!(history, historical)
        push!(offsets, length(terms)+1)
    end
    _TransportPointIndex(times, offsets, terms, history)
end

function _pulse_clip(p::_TransportPulse, a, b)
    lo, hi=max(p.start, a), min(p.stop, b)
    hi>lo ? _TransportPulse(lo, hi, p.coefficient, p.source, p.period) : nothing
end
_pulse_shift(p::_TransportPulse, d) =
    _TransportPulse(p.start+d, p.stop+d, p.coefficient, p.source, p.period)

function _pulse_expression(pulses, a, b)
    history=0.0
    terms=Dict{Tuple{Int,Int},Float64}()
    for p in pulses
        v=0.0036*p.coefficient*max(0.0, min(b, p.stop)-max(a, p.start))
        v==0 && continue
        if p.source==0
            history+=v
        else
            key=(p.source, p.period)
            terms[key]=get(terms, key, 0.0)+v
        end
    end
    (
        terms = [(k[1], k[2], v) for (k, v) in sort!(collect(terms); by = first)],
        history = history,
    )
end

"""Sparse exact arrival/release expressions for a deterministic river network.

Each `arrival_terms[d,t]` / `release_terms[d,t]` is a vector of
`(source_river, original_period, coefficient)` triples. Its dot product with
root releases in m³/s plus the corresponding `*_history[d,t]` gives volume in
Mm³. Non-root release rows are therefore derived quantities. Terminal terms
have the same units. Historical cohorts independently declared on each reach
represent water already entering that reach before horizon start; upstream
historical arrivals propagate downstream only after the horizon starts. This
clipping avoids counting pre-horizon water twice.

The evaluation grid is independent of the source decision grid. Pointwise
expressions and exact event knots are available with deterministic_point_data.
"""
function deterministic_network_data(c::ScheduleCase; arrival_grid = c.grid)
    order=river_order(c)
    RiverRouting.check_grid(arrival_grid)
    first(arrival_grid)==first(c.grid) && last(arrival_grid)==last(c.grid) ||
        throw(ArgumentError("transport evaluation grid must span the scheduling horizon"))
    s=c.system
    D=length(s.rivers)
    T=length(arrival_grid)-1
    lo, hi=first(c.grid), last(c.grid)
    all(r->r.deterministic_delay!==nothing, s.rivers) || throw(
        ArgumentError(
            "exact network operator requires deterministic delays on every reach",
        ),
    )
    reservoirs=Set(r.name for r in s.reservoirs)
    releases=[_TransportPulse[] for _ in 1:D]
    arrivals=[_TransportPulse[] for _ in 1:D]
    initial=zeros(D)
    for d in order
        r=s.rivers[d]
        delay=r.deterministic_delay
        isfinite(delay) && delay>=0 || throw(ArgumentError("invalid deterministic delay"))
        if r.source in reservoirs
            append!(
                releases[d],
                [
                    _TransportPulse(c.grid[k], c.grid[k + 1], 1.0, d, k) for
                    k in 1:(length(c.grid) - 1)
                ],
            )
        else
            for u in findall(e->e.target==r.source, s.rivers), p in arrivals[u]
                clipped=_pulse_clip(p, lo, hi)
                clipped===nothing || push!(releases[d], clipped)
            end
        end
        append!(arrivals[d], [_pulse_shift(p, delay) for p in releases[d]])
        for k in eachindex(r.history_release)
            p=_TransportPulse(
                r.history_grid[k],
                r.history_grid[k + 1],
                r.history_release[k],
                0,
                0,
            )
            push!(arrivals[d], _pulse_shift(p, delay))
        end
        initial[d]=remaining_volume(r, lo, r.history_grid, r.history_release)
    end
    AT=[Tuple{Int,Int,Float64}[] for _ in 1:D, _ in 1:T]
    RT=[Tuple{Int,Int,Float64}[] for _ in 1:D, _ in 1:T]
    AH=zeros(D, T)
    RH=zeros(D, T)
    TT=[Tuple{Int,Int,Float64}[] for _ in 1:D]
    TH=zeros(D)
    for d in 1:D
        for t in 1:T
            a=_pulse_expression(arrivals[d], arrival_grid[t], arrival_grid[t + 1])
            r=_pulse_expression(releases[d], arrival_grid[t], arrival_grid[t + 1])
            AT[d, t]=a.terms
            AH[d, t]=a.history
            RT[d, t]=r.terms
            RH[d, t]=r.history
        end
        a=_pulse_expression(arrivals[d], lo, hi)
        r=_pulse_expression(releases[d], lo, hi)
        coefficients=Dict{Tuple{Int,Int},Float64}()
        for (i, k, v) in r.terms
            coefficients[(i, k)]=get(coefficients, (i, k), 0.0)+v
        end
        for (i, k, v) in a.terms
            coefficients[(i, k)]=get(coefficients, (i, k), 0.0)-v
        end
        TT[d]=[
            (k[1], k[2], v) for
            (k, v) in sort!(collect(coefficients); by = first) if abs(v)>1e-15
        ]
        TH[d]=initial[d]+r.history-a.history
    end
    CompiledTransport(
        AT,
        AH,
        RT,
        RH,
        TT,
        TH,
        initial,
        releases,
        arrivals,
        Float64.(arrival_grid),
        copy(c.grid),
        order,
        [_transport_point_index(p) for p in arrivals],
        [_transport_point_index(p) for p in releases],
        _transport_signature(c),
    )
end

"""Exact one-sided rate expression (m³/s) and event knots for fixed delays."""
function deterministic_point_data(
    data::CompiledTransport,
    d,
    time;
    side = :right,
    kind = :arrival,
)
    side in (:left, :right) || throw(ArgumentError("side must be left or right"))
    kind in (:arrival, :release) || throw(ArgumentError("kind must be arrival or release"))
    isfinite(time) || throw(ArgumentError("transport query time must be finite"))
    index=kind==:arrival ? data.arrival_index[d] : data.release_index[d]
    k=side==:right ? searchsortedlast(index.times, time) :
      searchsortedfirst(index.times, time)-1
    k==0 ? (terms = _TransportTerms(index.terms, 1, 0), history = 0.0) :
    (
        terms = _TransportTerms(index.terms, index.offsets[k], index.offsets[k + 1]-1),
        history = index.history[k],
    )
end
function deterministic_knots(data::CompiledTransport, d; kind = :arrival)
    kind in (:arrival, :release) || throw(ArgumentError("kind must be arrival or release"))
    lo, hi=first(data.grid), last(data.grid)
    times=(kind==:arrival ? data.arrival_index[d] : data.release_index[d]).times
    a=searchsortedlast(times, lo)+1
    b=searchsortedfirst(times, hi)-1
    result=Vector{Float64}(undef, max(0, b-a+1)+2)
    result[1]=lo
    result[end]=hi
    copyto!(result, 2, times, a, max(0, b-a+1))
    result
end
_transport_value(x, q) = x.history+sum(v*q[d, k] for (d, k, v) in x.terms; init = 0.0)

"""Evaluate exact deterministic network routing, including restart cohorts."""
function route_network_exact(
    c::ScheduleCase,
    source_release::AbstractMatrix;
    grid = c.grid,
    transport = nothing,
)
    size(source_release)==(length(c.system.rivers), length(c.grid)-1) ||
        throw(DimensionMismatch("source releases"))
    all(isfinite, source_release) && all(>=(0), source_release) ||
        throw(ArgumentError("source releases must be finite and nonnegative"))
    data=_transport_data(c, transport; arrival_grid = grid)
    D=length(c.system.rivers)
    T=length(grid)-1
    A=zeros(D, T)
    R=zeros(D, T)
    W=zeros(D, T+1)
    cohorts=Any[]
    W[:, 1].=data.initial_transit
    for t in 1:T, d in 1:D
        A[d, t]=data.arrival_history[d, t]+sum(
            v*source_release[i, k] for (i, k, v) in data.arrival_terms[d, t];
            init = 0.0,
        )
        volume=data.release_history[d, t]+sum(
            v*source_release[i, k] for (i, k, v) in data.release_terms[d, t];
            init = 0.0,
        )
        R[d, t]=volume/(0.0036*(grid[t + 1]-grid[t]))
        W[d, t + 1]=W[d, t]+volume-A[d, t]
    end
    for d in 1:D
        knots=deterministic_knots(data, d; kind = :release)
        q=[
            _transport_value(
                deterministic_point_data(data, d, t; kind = :release),
                source_release,
            ) for t in knots[1:(end - 1)]
        ]
        all(x->x<=c.system.rivers[d].capacity+1e-10, q) || throw(
            ArgumentError(
                "instantaneous release exceeds capacity of $(c.system.rivers[d].name)",
            ),
        )
        push!(cohorts, (grid = knots, release = q))
    end
    Dict{String,Any}(
        "release"=>R,
        "arrival_volume"=>A,
        "transit"=>W,
        "order"=>data.order,
        "release_cohorts"=>cohorts,
        "exact"=>true,
        "converged"=>true,
        "transport_grid"=>Float64.(grid),
        "estimated_cumulative_error_Mm3"=>0.0,
        "estimated_transit_error_Mm3"=>0.0,
    )
end

function _transport_subdivide(grid, factor)
    vcat(
        [
            [grid[t]+(grid[t + 1]-grid[t])*j/factor for j in 0:(factor - 1)] for
            t in 1:(length(grid) - 1)
        ]...,
        last(grid),
    )
end

function _transport_evaluate_previous(c, previous, previous_grid, newgrid)
    reservoirs=Set(r.name for r in c.system.reservoirs)
    A=zeros(length(c.system.rivers), length(newgrid)-1)
    W=zeros(length(c.system.rivers), length(newgrid))
    for (d, r) in enumerate(c.system.rivers)
        oldgrid=r.source in reservoirs ? c.grid : previous_grid
        q=r.source in reservoirs ? previous["source_release"][d, :] :
          previous["release"][d, :]
        A[d, :]=route_volumes(r, newgrid, oldgrid, q)+route_volumes(
            r,
            newgrid,
            r.history_grid,
            r.history_release,
        )
        W[d, :]=[
            remaining_volume(r, t, oldgrid, q)+remaining_volume(
                r,
                t,
                r.history_grid,
                r.history_release,
            ) for t in newgrid
        ]
    end
    A, W
end

"""Converged transport audit on a grid independent of scheduling decisions.

All deterministic networks are evaluated exactly. General flow-dependent
networks refine junction mixing, retaining original reservoir release cohorts.
The estimate compares cumulative arrivals and physical transit on successive
grids; it is a numerical convergence check, not a rigorous continuum bound.
`converged=false` is explicit when the refinement limit is exhausted.
"""
function route_network_controlled(
    c::ScheduleCase,
    q::AbstractMatrix;
    grid = c.grid,
    absolute_tolerance = 1e-5,
    max_refinements = 6,
    max_internal_intervals = 2048,
    transport = nothing,
)
    isfinite(absolute_tolerance) && absolute_tolerance>0 ||
        throw(ArgumentError("positive transport tolerance required"))
    max_refinements>=1 ||
        throw(ArgumentError("at least one refinement comparison required"))
    max_internal_intervals>=1 ||
        throw(ArgumentError("positive internal interval limit required"))
    all(r->r.deterministic_delay!==nothing, c.system.rivers) &&
        return route_network_exact(c, q; grid, transport)
    transport===nothing || throw(
        ArgumentError("compiled deterministic transport cannot serve distributed routing"),
    )
    RiverRouting.check_grid(grid)
    first(grid)==first(c.grid) && last(grid)==last(c.grid) ||
        throw(ArgumentError("transport grid must span horizon"))
    current_grid=sort!(unique!(vcat(Float64.(grid), c.grid)))
    previous=route_network(c, q; grid = current_grid)
    previous["source_release"]=q
    error=Inf
    transiterror=Inf
    converged=false
    rounds=0
    for iteration in 1:max_refinements
        2(length(current_grid)-1)>max_internal_intervals && break
        newgrid=_transport_subdivide(current_grid, 2)
        current=route_network(c, q; grid = newgrid)
        oldA, oldW=_transport_evaluate_previous(c, previous, current_grid, newgrid)
        error=maximum(abs, cumsum(current["arrival_volume"]-oldA; dims = 2); init = 0.0)
        transiterror=maximum(abs, current["transit"]-oldW; init = 0.0)
        current["source_release"]=q
        previous=current
        current_grid=newgrid
        rounds=iteration
        if max(error, transiterror)<=absolute_tolerance
            converged=true
            break
        end
    end
    D=length(c.system.rivers)
    T=length(grid)-1
    A=zeros(D, T)
    R=zeros(D, T)
    for d in 1:D
        A[d, :]=0.0036 .* diff(grid) .* source_averages(
            grid,
            current_grid,
            previous["arrival_volume"][d, :] ./ (0.0036 .* diff(current_grid)),
        )
        R[d, :]=source_averages(grid, current_grid, previous["release"][d, :])
    end
    indices=[only(findall(==(t), current_grid)) for t in grid]
    Dict{String,Any}(
        "release"=>R,
        "arrival_volume"=>A,
        "transit"=>previous["transit"][:, indices],
        "order"=>previous["order"],
        "exact"=>false,
        "converged"=>converged,
        "transport_grid"=>current_grid,
        "refinements"=>rounds,
        "estimated_cumulative_error_Mm3"=>error,
        "estimated_transit_error_Mm3"=>transiterror,
        "termination"=>converged ? "tolerance" :
                       rounds==max_refinements ? "refinement_limit" : "interval_limit",
        "release_cohorts"=>[
            (
                grid = r.source in Set(x.name for x in c.system.reservoirs) ? copy(c.grid) :
                       current_grid,
                release = r.source in Set(x.name for x in c.system.reservoirs) ?
                          collect(q[d, :]) : collect(previous["release"][d, :]),
            ) for (d, r) in enumerate(c.system.rivers)
        ],
    )
end

"""Compare a stored schedule with independently controlled network transport."""
function transport_audit(
    c::ScheduleCase,
    x;
    absolute_tolerance = 1e-5,
    max_refinements = 6,
    schedule_tolerance = 2e-3,
    max_internal_intervals = 2048,
    transport = nothing,
)
    routed=route_network_controlled(
        c,
        x["river_release"];
        absolute_tolerance,
        max_refinements,
        max_internal_intervals,
        transport,
    )
    discrepancy=maximum(
        abs,
        cumsum(routed["arrival_volume"]-x["arrival_volume"]; dims = 2);
        init = 0.0,
    )
    terminal=maximum(abs, routed["transit"][:, end]-x["terminal_transit"]; init = 0.0)
    Dict(
        "accepted"=>routed["converged"] &&
                    max(discrepancy, terminal)+max(
            routed["estimated_cumulative_error_Mm3"],
            routed["estimated_transit_error_Mm3"],
        )<=schedule_tolerance,
        "converged"=>routed["converged"],
        "exact"=>routed["exact"],
        "cumulative_schedule_discrepancy_Mm3"=>discrepancy,
        "terminal_schedule_discrepancy_Mm3"=>terminal,
        "estimated_cumulative_error_Mm3"=>routed["estimated_cumulative_error_Mm3"],
        "estimated_transit_error_Mm3"=>routed["estimated_transit_error_Mm3"],
        "internal_intervals"=>length(routed["transport_grid"])-1,
        "scope"=>"transport with prescribed source releases; coupled hydraulic replay remains a separate check",
    )
end
