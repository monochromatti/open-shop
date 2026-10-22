"""Independent discrete-gradient hydraulic step. Rates are m³/s; storage is Mm³.

`river_arrival` contains known history/prior-cohort arrivals. `current_transfer`
contains the fractions of each current release cohort arriving in this same
interval, one column per reference curve. River junction mixing and current
arrivals are evaluated inside Newton's residual, coupled to storage and heads.
For deterministic transport, `current_network` maps local injections to arrivals;
upstream water entering a confluence is represented by those same source pulses.
Operational bounds are reported, never imposed by clipping physical states.
For a fixed-speed pump, optional `commitment` selects operation and its flow
entry is a Newton guess. Without commitment, a positive entry means on and zero
means off. The returned `generator_q` contains solved flows.
State and control inputs are normalized to Float64. Values that overflow, or
nonzero controls that underflow to zero, throw `ArgumentError`. Sparse
`current_network` matrices retain their sparse storage. All rate/state inputs
are vectors, and transfer inputs are matrices.
Optional `head_reference[node]` supplies a Newton head guess and selects the
constant pressure reference in floating junction components. It imposes no
additional head equation and cannot override a hydraulic law or head bound.
Newton uses automatic derivatives. Water-balance convergence includes a
two-ULP allowance for stored-volume rounding, converted to m³/s; other equations
use `tolerance` directly. `residual` retains the unadjusted maximum equation
error, so it can exceed `tolerance` even when `converged` is true.
"""
function forward_step(
    s::HydroSystem,
    Vprev,
    dt,
    generator_q,
    gate,
    river_arrival;
    tunnel_opening = nothing,
    current_transfer = zeros(length(s.rivers), 2),
    current_network = nothing,
    initial_guess = nothing,
    tolerance = 1e-9,
    max_iterations = 100,
    commitment = nothing,
    head_reference = nothing,
)
    validate_unit_modes(s)
    validate_delay_curves(s)
    validate_river_allocations(s)
    validate_curve_data(s)
    nr, nj, ne, nd=length(s.reservoirs), length(s.junctions), length(s.tunnels), length(s.rivers)
    ng=length(s.generators)
    ta=tunnel_opening===nothing ? [t.opening for t in s.tunnels] : tunnel_opening
    Vprev isa AbstractVector && length(Vprev)==nr &&
    generator_q isa AbstractVector && length(generator_q)==ng &&
    gate isa AbstractVector && length(gate)==nd &&
    river_arrival isa AbstractVector && length(river_arrival)==nd &&
    ta isa AbstractVector &&
    length(ta)==ne || throw(DimensionMismatch("forward inputs"))
    current_transfer isa AbstractMatrix && size(current_transfer, 1)==nd &&
    size(current_transfer, 2)>=maximum((length(r.curves) for r in s.rivers); init = 1) ||
        throw(DimensionMismatch("current transfer fractions"))
    isfinite(dt) && dt>0 || throw(ArgumentError("finite positive duration required"))
    all(>=(0), generator_q) &&
    all(x->0<=x<=1, gate) &&
    all(x->0<=x<=1, ta) &&
    all(>=(0), river_arrival) &&
    all(x->0<=x<=1, current_transfer) ||
        throw(ArgumentError("invalid controls or arrival fractions"))
    current_network === nothing ||
        (current_network isa AbstractMatrix && size(current_network) == (nd, nd)) ||
        throw(ArgumentError("current network must be a river-by-river matrix"))
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("positive finite tolerance required"))
    max_iterations isa Integer && max_iterations >= 0 ||
        throw(ArgumentError("nonnegative iteration limit required"))
    data = HydraulicStepData(
        s,
        Vprev,
        dt,
        generator_q,
        gate,
        river_arrival,
        ta,
        current_transfer,
        current_network,
        ;commitment,
    )
    reference=head_reference===nothing ? nothing : Vector{Float64}(head_reference)
    reference===nothing || (length(reference)==length(nodes(s)) && all(isfinite,reference)) ||
        throw(ArgumentError("head reference must contain one finite value per hydraulic node"))
    _forward_step(data; initial_guess, tolerance, max_iterations,head_reference=reference)
end

function _forward_step(data::HydraulicStepData; initial_guess, tolerance, max_iterations,
    head_reference=nothing)
    s, topology = data.system, data.topology
    nr, nj, ne = length(s.reservoirs), length(s.junctions), length(s.tunnels)
    base, nq = nr + nj, nr + nj + ne
    active = data.active_rivers
    offset=nq+length(active)
    nv = offset+length(data.fixed_pumps)
    ix = topology.node_index
    Vprev, generator_q, gate, ta = data.previous_storage, data.generator_q, data.gate, data.opening
    work = HydraulicWorkspace(data)
    x=zeros(nv)
    for (k,i) in enumerate(data.fixed_pumps)
        x[offset+k]=generator_q[i]
    end
    x[1:nr].=Vprev
    known=vcat(
        [head(r, Vprev[i]) for (i, r) in enumerate(s.reservoirs)],
        [b.head for b in s.boundaries],
    )
    x[(nr + 1):base].=isempty(known) ? 0.0 : sum(known)/length(known)
    head_reference===nothing || (x[(nr+1):base].=head_reference[(nr+1):base])
    h = interval_heads!(work.head, x, data)
    for (k, t) in enumerate(s.tunnels)
        dh=h[ix[t.source]]-h[ix[t.target]]
        x[base + k]=tunnel_flow(t, dh, ta[k])
    end
    for (k, d) in enumerate(active)
        x[nq + k]=data.river_inflow[d]+outlet_flow(
            s.rivers[d],
            h,
            gate[d],
            topology.river_source[d],
            topology.river_target[d],
            inflow = data.river_inflow[d],
            capacity = data.river_capacity[d],
        )
    end
    if initial_guess!==nothing
        if initial_guess isa HydraulicStep
            x[1:nr].=initial_guess.storage
            x[(nr + 1):base].=initial_guess.head[(nr + 1):base]
            x[(base + 1):nq].=initial_guess.tunnel_q
            x[(nq + 1):offset].=initial_guess.river_release[active]
            x[(offset+1):end].=initial_guess.generator_q[data.fixed_pumps]
        elseif initial_guess isa AbstractDict
            x[1:nr].=initial_guess["Vnew"]
            x[(nr + 1):base].=initial_guess["H"][(nr + 1):base]
            x[(base + 1):nq].=initial_guess["tunnel_q"]
            x[(nq + 1):offset].=initial_guess["river_release"][active]
            haskey(initial_guess,"generator_q") &&
                (x[(offset+1):end].=initial_guess["generator_q"][data.fixed_pumps])
        else
            length(initial_guess)==nv || throw(DimensionMismatch("initial guess"))
            x.=initial_guess
        end
    end
    # q|q| has zero derivative at q=0. A zero or opposed flow guess can land
    # there during Newton. Seed those guesses from the passive pressure law;
    # retain warm flows that already have the pressure difference's direction.
    h = interval_heads!(work.head, x, data)
    for (j, tunnel) in enumerate(s.tunnels)
        source, target = topology.tunnel_nodes[j]
        difference = h[source] - h[target]
        if ta[j] > 0 && sign(x[base + j]) != sign(difference)
            x[base + j] = tunnel_flow(tunnel, difference, ta[j])
        end
    end
    iterations = solve_hydraulics!(x, data, work; tolerance, max_iterations,head_reference)
    h = interval_heads!(work.head, x, data)
    for k in eachindex(s.tunnels)
        ta[k]==0 && (x[base + k]=0.0)
    end
    # Exact law evaluation removes Newton roundoff from a shut outlet. Its
    # storage consequence is included in the recomputed convergence residual.
    for (k, d) in enumerate(active)
        x[nq + k]=data.river_inflow[d]+outlet_flow(
            s.rivers[d],
            h,
            gate[d],
            topology.river_source[d],
            topology.river_target[d],
            inflow = data.river_inflow[d],
            capacity = data.river_capacity[d],
        )
    end
    river_rates!(work, x, data)
    generator_q=work.generator_q
    release, arrival, injection = work.release, work.arrival, work.injection
    gh=[
        begin
            p=topology.generator_plant[i]
            plant=s.plants[p]
            unit_head(s.generators[i],plant,h[ix[plant.source]],h[ix[plant.target]],
                sum(generator_q[j] for j in topology.plant_generators[p]; init = 0.0))
        end for i in eachindex(s.generators)
    ]
    pw=[power(g, generator_q[i], gh[i]) for (i, g) in enumerate(s.generators)]
    exchange=zeros(length(s.boundaries))
    function boundary!(node, q)
        j=get(ix, node, 0)-base
        j>0 && (exchange[j]+=q)
    end
    for (i, t) in enumerate(s.tunnels)
        boundary!(t.source, -x[base + i])
        t.discharge_river===nothing && boundary!(t.target, x[base + i])
    end
    for i in eachindex(s.generators)
        p=s.plants[topology.generator_plant[i]]
        boundary!(p.source, -generator_q[i])
        p.discharge_river===nothing && boundary!(p.target, generator_q[i])
    end
    for (i, r) in enumerate(s.rivers)
        boundary!(r.source, -(release[i]-data.river_inflow[i]))
        boundary!(r.target, arrival[i])
    end
    residual_vector = hydraulic_residual!(similar(x), x, data, work)
    residual = maximum(abs, residual_vector; init = 0.0)
    converged = hydraulic_converged(residual_vector, x, data, tolerance)
    domain=all(isfinite, vcat(x, h, pw, release, arrival)) &&
           all(
               _forward_domain(r, release[d];
                   capacity = r.law in (:controlled,:delta) ? data.river_capacity[d] : r.capacity)
               for (d, r) in enumerate(s.rivers)
           ) &&
           all(
               t.discharge_river===nothing || x[base + i]>=-1e-10 for (i, t) in enumerate(s.tunnels)
           )
    violations=HydraulicBoundViolations(
        maximum(
            (max(r.vmin-x[i], x[i]-r.vmax, 0.0) for (i, r) in enumerate(s.reservoirs));
            init = 0.0,
        ),
        maximum(
            (
                max(
                    abs(x[base + i])-data.tunnel_capacity[i],
                    r.discharge_river===nothing ? 0.0 : -x[base + i],
                    0.0,
                ) for (i, r) in enumerate(s.tunnels)
            );
            init = 0.0,
        ),
        maximum(
            (
                max(
                    river_flow_lower(r, data.river_capacity[d])-release[d],
                    release[d]-data.river_capacity[d],
                    0.0,
                ) for
                (d, r) in enumerate(s.rivers)
            );
            init = 0.0,
        ),
    )
    HydraulicStep(
        x[1:nr], h, x[(base + 1):nq], copy(ta), release, injection, arrival,
        gh, copy(generator_q), pw, residual, converged && domain,
        iterations, exchange, domain, violations,
        all(
            is_channel(s.rivers[d]) ||
                s.rivers[d].law==:controlled ||
                s.rivers[d].allow_dry ||
                s.rivers[d].discharge_curve!==nothing ||
                h[ix[s.rivers[d].source]]>=s.rivers[d].crest for d in active
        ),
    )
end

# Topological reach evaluation also accepts boundary targets used by standalone
# integration tests. River graph/input validation belongs to the public model.
function _forward_river_order(s)
    validate_river_allocations(s)
    remaining=Set(eachindex(s.rivers))
    order=Int[]
    while !isempty(remaining)
        ready=sort!([
            d for d in remaining if !any(
                !is_channel(s.rivers[d]) &&
                    !is_channel(s.rivers[k]) &&
                    k in remaining &&
                    s.rivers[k].target==s.rivers[d].source for k in eachindex(s.rivers)
            )
        ])
        isempty(ready) && throw(ArgumentError("river graph must be acyclic"))
        for d in ready
            push!(order, d)
            delete!(remaining, d)
        end
    end
    order
end
_forward_atom(r) = hasproperty(r, :deterministic_delay) && r.deterministic_delay!==nothing
function _forward_weights(r, q)
    _forward_atom(r) && return [(1, 1.0)]
    length(r.curves)==1 && return [(1, 1.0)]
    refs=[c.reference_flow for c in r.curves]
    RiverRouting.blend(r.curves, clamp(q, first(refs), last(refs)))
end
function _forward_domain(r, q; capacity=r.capacity)
    river_flow_valid(r,q;capacity) || return false
    _forward_atom(r) ||
        length(r.curves)==1 ||
        first(r.curves).reference_flow-1e-10<=q<=last(r.curves).reference_flow+1e-9
end
function _forward_coefficients(grid, r)
    _forward_atom(r) || return RiverRouting.coefficients(grid, grid, r.curves)
    a=RiverRouting.fixed_coefficients(grid, grid, r.deterministic_delay)
    reshape(a, size(a, 1), size(a, 2), 1)
end
function _forward_history(grid, r)
    _forward_atom(r) || return (
        RiverRouting.route(grid, r.history_grid, r.history_release, r.curves),
        RiverRouting.remaining_volume(first(grid), r.history_grid, r.history_release, r.curves),
    )
    RiverRouting.check_releases(r.history_grid, r.history_release)
    a=RiverRouting.fixed_coefficients(grid, r.history_grid, r.deterministic_delay)*(
        0.0036 .* diff(r.history_grid) .* r.history_release
    )
    time=first(grid)
    w=sum(
        0.0036*r.history_release[k]*(
            max(0.0, min(time, r.history_grid[k + 1])-r.history_grid[k])-max(
                0.0,
                min(time-r.deterministic_delay, r.history_grid[k + 1])-r.history_grid[k],
            )
        ) for k in eachindex(r.history_release)
    )
    a, w
end

"""
    simulate(case; generator_q, gate, tunnel_opening=nothing, grid=case.grid)
    simulate(case, generator_q, gate; kwargs...)

Independent chronological replay with coupled same-step river arrivals.
`generator_q[unit, interval]` is discharge in m³/s; `gate[river,
interval]` is a river outlet opening in [0,1]. `tunnel_opening[tunnel, interval]`
is required when tunnels have selectable positions. Input arrays are read-only.

Returns a [`HydraulicTrajectory`](@ref) with `storage[reservoir, vertex]` in Mm³,
volume-averaged `head[node, interval]` in metres, `power[unit, interval]` in MW, river
releases in m³/s, arrivals/transit in Mm³, and convergence diagnostics. A converged physical
trajectory can violate operating limits; use [`dispatch_from_controls`](@ref)
to reconstruct and validate a complete schedule. Malformed inputs throw an
`ArgumentError` or `DimensionMismatch`. Optional `commitment[unit,interval]`
selects fixed-speed pump operation independently of its discharge guess. Without
it, positive flow entries mean on. Fixed-speed pump entries are guesses and
their returned discharges follow the characteristic at each replay interval head;
other unit discharges remain prescribed and must be zero when their supplied
commitment is off. A committed adjustable unit may have zero prescribed flow.
Optional `head_reference[node, control_interval]` supplies Newton guesses and
retains the chosen reference in floating pressure components. It does not
replace hydraulic equations or operating bounds; validate the returned schedule.

Source outlet laws, volume-averaged heads and zero/short-delay transfers are solved
simultaneously during each forward step. Deterministic routing preserves each
local injection's pulses through confluences. Distributed junction mixing
conserves volume but rebins arriving water uniformly on the replay grid; refine
the grid to assess this timing approximation. Each cohort contributes once.
"""
function simulate(
    c::ScheduleCase,
    generator_q,
    gate;
    grid = c.grid,
    transport = nothing,
    tunnel_opening = nothing,
    commitment = nothing,
    head_reference = nothing,
)
    s=c.system
    validate_unit_modes(s)
    validate_soft_constraints(c)
    validate_delay_curves(s)
    validate_river_allocations(s)
    validate_curve_data(s)
    ng, nd=length(s.generators), length(s.rivers)
    T=length(grid)-1
    reference=if head_reference===nothing
        nothing
    else
        head_reference isa AbstractMatrix && size(head_reference)==(length(nodes(s)),length(c.prices)) ||
            throw(DimensionMismatch("head reference must be node by control interval"))
        converted=Matrix{Float64}(head_reference)
        all(isfinite,converted) || throw(ArgumentError("head reference must remain finite in Float64"))
        converted
    end
    size(generator_q)==(ng, length(c.grid)-1) && size(gate)==(nd, length(c.grid)-1) ||
        throw(DimensionMismatch("schedule dimensions"))
    state=commitment===nothing ? generator_q.>0 : commitment
    size(state)==size(generator_q) && all(x->x==0||x==1,state) ||
        throw(ArgumentError("binary unit commitment dimensions required"))
    commitment_exclusions_admissible(c,state) ||
        throw(ArgumentError("unit discharges violate commitment exclusion"))
    all(>=(0), generator_q) && all(x -> 0<=x<=1, gate) ||
        throw(ArgumentError("invalid generator discharges or river gates"))
    RiverRouting.check_grid(grid)
    first(grid)==first(c.grid) && last(grid)==last(c.grid) ||
        throw(ArgumentError("replay grid must span original horizon"))
    all(t->any(isapprox(t, x; atol = 1e-12) for x in grid), c.grid) ||
        throw(ArgumentError("replay grid must contain original control change times"))
    indices=[min(searchsortedlast(c.grid, (grid[t]+grid[t + 1])/2), length(c.grid)-1) for t in 1:T]
    gq=Float64.(generator_q[:, indices])
    all(isfinite,gq) && all(iszero(generator_q[i,indices[t]])==iszero(gq[i,t])
        for i in 1:ng,t in 1:T) || throw(ArgumentError(
            "unit flows must remain finite in Float64 without changing zero states"))
    ga=gate[:, indices]
    opening=_physical_tunnel_openings(c,tunnel_opening)
    validate_tunnel_openings(c,opening,state)
    ta=opening[:,indices]
    nr=length(s.reservoirs)
    topology=HydraulicTopology(s)
    sources=river_source_indices(s)
    active_rivers=findall(r -> r.law != :junction, s.rivers)
    river_order=_forward_river_order(s)
    V=zeros(nr, T+1)
    V[:, 1].=[r.v0 for r in s.reservoirs]
    H=zeros(length(topology.node_names), T)
    Q=zeros(length(s.tunnels), T)
    P=zeros(ng, T)
    R=zeros(nd, T)
    L=zeros(nd, T)
    A=zeros(nd, T)
    W=zeros(nd, T+1)
    out=zeros(T)
    steps=HydraulicStep[]
    reservoir_inflow=zeros(nr)
    river_inflow=zeros(nd)
    tunnel_capacity=zeros(length(s.tunnels))
    river_capacity=zeros(nd)
    exact=all(r.deterministic_delay!==nothing for r in s.rivers)
    nc_case=_river_replace(c; grid = Float64.(grid), prices = zeros(T))
    ndat=exact ? _transport_data(nc_case, transport) : nothing
    !exact &&
        transport!==nothing &&
        throw(ArgumentError("compiled deterministic transport cannot serve distributed replay"))
    B=exact ? Array{Float64, 3}[] : [_forward_coefficients(grid, r) for r in s.rivers]
    # Retain reference-curve column dimensions even when the exact network
    # matrix supplies the current arrivals instead of reference fractions.
    nc=exact ? maximum((length(r.curves) for r in s.rivers); init = 1) :
       maximum((size(b, 3) for b in B); init = 1)
    transfer=zeros(nd, nc)
    network=exact ? zeros(nd, nd) : nothing
    sizehint!(steps, T)
    if exact
        W[:, 1].=ndat.initial_transit
    else
        for (d, r) in enumerate(s.rivers)
            ah, wh=_forward_history(grid, r)
            A[d, :].=ah
            W[d, 1]=wh
        end
    end
    for t in 1:T
        dt=grid[t + 1]-grid[t]
        if !exact
            for d in 1:nd, l in axes(B[d], 3)
                transfer[d, l]=B[d][t, t, l]
            end
        end
        if exact
            fill!(network, 0.0)
            for d in 1:nd
                A[d, t]=ndat.arrival_history[d, t]
                for (i, k, v) in ndat.arrival_terms[d, t]
                    if k<t
                        A[d, t]+=v*L[i, k]
                    elseif k==t
                        network[d, i]+=v/(0.0036dt)
                    else
                        error("noncausal river transport")
                    end
                end
            end
        end
        for (i, r) in enumerate(s.reservoirs)
            reservoir_inflow[i]=opaverage(c, r.name, :inflow, grid[t], grid[t + 1], r.inflow)
        end
        for (i, e) in enumerate(s.tunnels)
            tunnel_capacity[i]=opaverage(c, e.name, :capacity, grid[t], grid[t + 1], e.capacity)
        end
        for (d, r) in enumerate(s.rivers)
            river_inflow[d]=opaverage(c, r.name, :inflow, grid[t], grid[t + 1], r.inflow)
            river_capacity[d]=opaverage(c, r.name, :capacity, grid[t], grid[t + 1], r.capacity)
        end
        data=HydraulicStepData(
            s,
            @view(V[:, t]),
            dt,
            @view(gq[:, t]),
            @view(ga[:, t]),
            A[:, t]/(0.0036dt),
            @view(ta[:, t]),
            transfer,
            network;
            topology, sources, active_rivers, river_order,
            reservoir_inflow, river_inflow, tunnel_capacity, river_capacity,
            commitment=@view(state[:,indices[t]]),
        )
        z=_forward_step(
            data;
            initial_guess = isempty(steps) ? nothing : last(steps),
            tolerance = 1e-9,
            max_iterations = 100,
            head_reference = reference===nothing ? nothing : view(reference,:,indices[t]),
        )
        push!(steps, z)
        V[:, t + 1].=z.storage
        H[:, t].=z.head
        Q[:, t].=z.tunnel_q
        P[:, t].=z.power
        gq[:,t].=z.generator_q
        R[:, t].=z.river_release
        L[:, t].=z.river_injection
        out[t]=sum(z.boundary_net_inflow)
        A[:, t].=0.0036dt .* z.river_arrival
        for (d, r) in enumerate(s.rivers)
            if !exact
                for (l, w) in _forward_weights(r, R[d, t])
                    A[d, (t + 1):end].+=0.0036dt*R[d, t]*w .* B[d][(t + 1):end, t, l]
                end
            end
            W[d, t + 1]=W[d, t]+0.0036dt*R[d, t]-A[d, t]
        end
    end
    HydraulicTrajectory(
        collect(grid), V, H, Q, ta, gq, ga, P, R, L, A, W, out,
        all(z.converged && z.law_domain_valid for z in steps),
        maximum((z.residual for z in steps); init = 0.0), steps,
    )
end

# Named controls avoid mixing two physically different discharge/gate arrays.
function simulate(
    c::ScheduleCase;
    generator_q,
    gate,
    grid = c.grid,
    transport = nothing,
    tunnel_opening = nothing,
    commitment = nothing,
    head_reference = nothing,
)
    simulate(c, generator_q, gate; grid, transport, tunnel_opening,commitment,head_reference)
end
