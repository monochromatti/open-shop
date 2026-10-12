"""Data for one implicit midpoint hydraulic step, with prescribed controls.

The unknown vector is `[Vnew; junction_head; tunnel_q; outlet_release]`.
Storage is in Mm³, heads in metres, rates in m³/s, and duration in hours.
`network` is either a same-interval river transfer matrix or `nothing` for
flow-dependent response curves. The type parameter specializes that choice.
"""
struct HydraulicStepData{M, S}
    system::HydroSystem
    topology::HydraulicTopology
    previous_storage::Vector{Float64}
    duration::Float64
    generator_q::Vector{Float64}
    gate::Vector{Float64}
    opening::Vector{Float64}
    previous_arrival::Vector{Float64}
    transfer::Matrix{Float64}
    network::M
    sources::S
    active_rivers::Vector{Int}
    river_order::Vector{Int}
    fixed_inflow::Vector{Float64}
end

function HydraulicStepData(s, Vprev, dt, generator_q, gate, arrival, opening, transfer, network)
    previous_storage = Vector{Float64}(Vprev)
    duration = Float64(dt)
    discharge = Vector{Float64}(generator_q)
    river_gate = Vector{Float64}(gate)
    tunnel_opening = Vector{Float64}(opening)
    previous_arrival = Vector{Float64}(arrival)
    current_transfer = Matrix{Float64}(transfer)
    current_network = network === nothing ? nothing : Float64.(network)
    isfinite(duration) && duration > 0 ||
        throw(ArgumentError("duration must remain finite and positive in Float64"))
    all(
        A -> all(isfinite, A),
        (
            previous_storage,
            discharge,
            river_gate,
            tunnel_opening,
            previous_arrival,
            current_transfer,
        ),
    ) && (current_network === nothing || all(isfinite, current_network)) ||
        throw(ArgumentError("hydraulic inputs must remain finite in Float64"))
    # Exact zero identifies closed gates and inactive generating units. Reject
    # underflow instead of changing that physical state during normalization.
    for (original, normalized) in
        ((generator_q, discharge), (gate, river_gate), (opening, tunnel_opening))
        all(i -> iszero(original[i]) == iszero(normalized[i]), eachindex(normalized)) ||
            throw(ArgumentError("nonzero hydraulic controls must remain nonzero in Float64"))
    end
    topology = HydraulicTopology(s)
    nr, nj = length(s.reservoirs), length(s.junctions)
    fixed = zeros(nr + nj)
    for (i, reservoir) in enumerate(s.reservoirs)
        fixed[i] = reservoir.inflow
    end
    for i in eachindex(fixed)
        for (j, sign) in topology.generator_incidence[i]
            fixed[i] += sign * discharge[j]
        end
    end
    active = findall(r -> r.law != :junction, s.rivers)
    all(d -> 1 <= topology.river_source[d] <= nr, active) ||
        throw(ArgumentError("operational river law requires reservoir source"))
    HydraulicStepData(
        s,
        topology,
        previous_storage,
        duration,
        discharge,
        river_gate,
        tunnel_opening,
        previous_arrival,
        current_transfer,
        current_network,
        river_source_indices(s),
        active,
        _forward_river_order(s),
        fixed,
    )
end

"""Scratch arrays for residual evaluation. No solver or result metadata."""
struct HydraulicWorkspace{T}
    head::Vector{T}
    net_inflow::Vector{T}
    release::Vector{T}
    arrival::Vector{T}
    injection::Vector{T}
end

function HydraulicWorkspace(data::HydraulicStepData, ::Type{T} = Float64) where {T}
    s = data.system
    nn = length(data.topology.node_names)
    nr, nj, nd = length(s.reservoirs), length(s.junctions), length(s.rivers)
    HydraulicWorkspace(zeros(T, nn), zeros(T, nr + nj), zeros(T, nd), zeros(T, nd), zeros(T, nd))
end

function midpoint_heads!(h, x, data::HydraulicStepData)
    s = data.system
    nr, nj = length(s.reservoirs), length(s.junctions)
    for (i, reservoir) in enumerate(s.reservoirs)
        h[i] = head(reservoir, (data.previous_storage[i] + x[i]) / 2)
    end
    for i in (nr + 1):(nr + nj)
        h[i] = x[i]
    end
    for (i, boundary) in enumerate(s.boundaries)
        h[nr + nj + i] = boundary.head
    end
    h
end

function outlet_flow(r::River, h, gate, source, target)
    is_channel(r) ? channel_flow(r, h[source], h[target]) : river_law_value(r, h[source], gate)
end

# Clamping only extends trial-point residuals. Accepted releases are checked
# against the original routing domain; no accepted physical state is clipped.
function _current_arrival_fraction(r, q, transfer, d)
    (_forward_atom(r) || length(r.curves) == 1) && return transfer[d, 1]
    curves = r.curves
    q = clamp(q, first(curves).reference_flow, last(curves).reference_flow)
    j = RiverRouting.reference_segment(curves, q; reference = c -> c.reference_flow)
    w = (q - curves[j].reference_flow) / (curves[j + 1].reference_flow - curves[j].reference_flow)
    (1 - w) * transfer[d, j] + w * transfer[d, j + 1]
end

function river_rates!(work, x, data::HydraulicStepData)
    s, topology = data.system, data.topology
    nr, nj, ne = length(s.reservoirs), length(s.junctions), length(s.tunnels)
    release, arrival, injection = work.release, work.arrival, work.injection
    fill!(release, 0)
    copyto!(arrival, data.previous_arrival)
    for (k, d) in enumerate(data.active_rivers)
        release[d] = x[nr + nj + ne + k]
    end
    for (d, river) in enumerate(s.rivers)
        injection[d] = if data.sources.reservoirs[d]
            release[d]
        else
            river.inflow +
            sum(data.generator_q[j] for j in data.sources.generators[d]; init = 0.0) +
            sum(x[nr + nj + j] for j in data.sources.tunnels[d]; init = 0.0)
        end
    end
    if data.network !== nothing
        mul!(arrival, data.network, injection, 1, 1)
    end
    for d in data.river_order
        river = s.rivers[d]
        if river.law == :junction
            release[d] =
                injection[d] + sum(arrival[j] for j in topology.river_upstream[d]; init = 0.0)
        end
        if data.network === nothing
            arrival[d] +=
                is_channel(river) ? release[d] :
                max(0.0, release[d]) *
                _current_arrival_fraction(river, release[d], data.transfer, d)
        end
    end
    work
end

"""Evaluate storage balance, junction balance, tunnel losses and outlet laws.

For a reservoir, `(Vnew - Vprev) / (0.0036Δt) = net_inflow`.
For a junction, `net_inflow = 0`. Heads use midpoint storage. Every physical
balance is retained, including those in floating pressure components.
The scratch arrays and `residual` are overwritten; `x` and `data` are read-only.
"""
function hydraulic_residual!(residual, x, data::HydraulicStepData, work::HydraulicWorkspace)
    s, topology = data.system, data.topology
    nr, nj, ne = length(s.reservoirs), length(s.junctions), length(s.tunnels)
    base, nq = nr + nj, nr + nj + ne
    h, net = work.head, work.net_inflow
    midpoint_heads!(h, x, data)
    river_rates!(work, x, data)
    for i in eachindex(net)
        net[i] =
            sum(sign * x[base + j] for (j, sign) in topology.tunnel_incidence[i]; init = 0.0) +
            data.fixed_inflow[i]
    end
    for (d, river) in enumerate(s.rivers)
        source, target = topology.river_source[d], topology.river_target[d]
        1 <= source <= base && (net[source] -= work.release[d] - river.inflow)
        1 <= target <= base && (net[target] += work.arrival[d])
    end
    for i in 1:nr
        residual[i] = (x[i] - data.previous_storage[i]) / (0.0036 * data.duration) - net[i]
    end
    for i in (nr + 1):base
        residual[i] = -net[i]
    end
    for (j, tunnel) in enumerate(s.tunnels)
        source, target = topology.tunnel_nodes[j]
        residual[base + j] =
            tunnel_loss_residual(tunnel, h[source] - h[target], x[base + j], data.opening[j])
    end
    for (k, d) in enumerate(data.active_rivers)
        residual[nq + k] =
            x[nq + k] - s.rivers[d].inflow - outlet_flow(
                s.rivers[d],
                h,
                data.gate[d],
                topology.river_source[d],
                topology.river_target[d],
            )
    end
    residual
end

"""Damped Newton iteration with central differences and reusable residual buffers.

One head column per pressure nullspace is fixed. All residual rows remain, so
Julia's rectangular `\\` solve uses least squares for those components. This
cannot hide incompatible water injections by discarding a conservation row.
"""
function solve_hydraulics!(x, data::HydraulicStepData, work; tolerance, max_iterations)
    floating = _floating_pressure_components(
        data.system,
        data.opening,
        data.generator_q;
        ix = data.topology.node_index,
    )
    _normalize_pressure_references!(x, data.system, floating)
    columns = isempty(floating) ? collect(eachindex(x)) : setdiff(eachindex(x), first.(floating))
    residual, plus, minus, trial, step = (similar(x) for _ in 1:5)
    jacobian = Matrix{Float64}(undef, length(x), length(columns))
    iterations = 0
    for k in 1:max_iterations
        hydraulic_residual!(residual, x, data, work)
        maximum(abs, residual; init = 0.0) <= tolerance && break
        copyto!(trial, x)
        for (column, j) in enumerate(columns)
            ε = 1e-5 * max(abs(x[j]), 1.0)
            trial[j] = x[j] + ε
            hydraulic_residual!(plus, trial, data, work)
            trial[j] = x[j] - ε
            hydraulic_residual!(minus, trial, data, work)
            for i in eachindex(residual)
                jacobian[i, column] = (plus[i] - minus[i]) / (2ε)
            end
            trial[j] = x[j]
        end
        direction = try
            -(jacobian \ residual)
        catch error
            error isa Union{SingularException, LAPACKException, RankDeficientException} ||
                rethrow()
            -pinv(jacobian) * residual
        end
        all(isfinite, direction) || (direction = -pinv(jacobian) * residual)
        fill!(step, 0)
        step[columns] .= direction
        merit = sum(abs2, residual)
        accepted, scale = false, 1.0
        for _ in 1:45
            @. trial = x + scale * step
            hydraulic_residual!(plus, trial, data, work)
            if all(isfinite, plus) && sum(abs2, plus) < merit
                copyto!(x, trial)
                accepted = true
                break
            end
            scale /= 2
        end
        iterations = k
        accepted || break
    end
    _normalize_pressure_references!(x, data.system, floating)
    iterations
end
