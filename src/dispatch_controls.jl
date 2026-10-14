# Physical controls are explicit for every optimized tunnel. Legacy schedules
# without gate decisions retain prescribed operating-series/static defaults.
function _tunnel_opening_input(c, opening)
    if opening===nothing
        any(e->!isempty(e.gate_positions), c.system.tunnels) && throw(
            ArgumentError("optimized tunnel gates require a complete tunnel_opening schedule"),
        )
        return tunnel_opening_schedule(c)
    end
    opening isa AbstractMatrix && size(opening)==(length(c.system.tunnels), length(c.prices)) ||
        throw(DimensionMismatch("tunnel opening schedule dimensions"))
    all(x->x isa Real && isfinite(x) && 0<=x<=1, opening) ||
        throw(ArgumentError("tunnel openings must be finite real values in [0,1]"))
    converted=Matrix{Float64}(opening)
    all(isfinite,converted) || throw(ArgumentError("tunnel openings must remain finite in Float64"))
    for k in eachindex(converted,opening)
        iszero(converted[k])==iszero(opening[k]) ||
            throw(ArgumentError("nonzero tunnel openings must remain nonzero in Float64"))
    end
    converted
end
function _physical_tunnel_openings(c, opening)
    A=_tunnel_opening_input(c, opening)
    validate_tunnel_openings(c, A)
    # Prescribed interval averages may differ by arithmetic roundoff after
    # grid refinement. Use the prescribed value; selectable modes stay exact.
    for (i, e) in enumerate(c.system.tunnels)
        isempty(e.gate_positions) || continue
        for t in axes(A, 2)
            A[i, t]=opinterval(c, e.name, :opening, t, e.opening)
        end
    end
    validate_tunnel_openings(c, A; tolerance = 0)
    A
end

"""
    dispatch_from_controls(case; commitment, generator_q, gate, tunnel_opening=nothing)
    dispatch_from_controls(case, commitment, generator_q, gate; kwargs...)

Reconstruct and validate a physical schedule from prescribed controls.
`commitment[unit, interval]` is 0 or 1, `generator_q[unit, interval]` is m³/s,
and `gate[river, interval]` is an outlet opening in [0,1]. Supply an explicit
`tunnel_opening[tunnel, interval]` for selectable tunnel gates. Optional `reserves`
contains capacity controls; their unit/plant flow witnesses are inferred at the
scheduled hydraulic heads and never enter nominal water balances.

Returns a [`Schedule`](@ref) with physical arrays, `commitment` and `objective`.
Inspect `schedule.diagnostics["validation"]["valid"]` before using it. Invalid commitment
or nonconvergent integration throws `ArgumentError`. No bound, gate or state
is clipped, and the input arrays are not mutated. This function does not optimize.
"""
function dispatch_from_controls(
    c,
    u,
    q,
    gate;
    status = "FORWARD_RECONSTRUCTED",
    transport = nothing,
    tunnel_opening = nothing,
    reserves = nothing,
)
    admissible(c, u) || throw(ArgumentError("inadmissible commitment"))
    opening=_physical_tunnel_openings(c,tunnel_opening)
    validate_tunnel_openings(c,opening,u)
    nd=all(r->r.deterministic_delay!==nothing, c.system.rivers) ? _transport_data(c, transport) :
       nothing
    nd===nothing &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed reconstruction",
            ),
        )
    z=simulate(c, q, gate; transport = nd, tunnel_opening=opening)
    z.converged || throw(ArgumentError("forward reconstruction did not converge"))
    s=c.system
    dt=diff(c.grid)
    T=length(dt)
    req=release_requirements(c, z.river_release)
    operating_cost=transition_costs(c, u)
    adjustment_cost=tunnel_adjustment_cost(c, z.tunnel_opening)
    reserves=reserves===nothing ? _zero_reserves(c,q) :
        reserve_witnesses(c,_reserve_input(reserves,T),u,z)
    objective=sum(c.prices[t]*dt[t]*unit_power_sign(s.generators[j])*z.power[j,t]
        for j in eachindex(s.generators), t in 1:T; init=0.0) -
        operating_cost-req.cost-adjustment_cost +
        sum(r.water_value*(z.storage[i,end]-r.v0)
            for (i,r) in enumerate(s.reservoirs); init=0.0) +
        sum(r.water_value*(z.transit[i,end]-z.transit[i,1])
            for (i,r) in enumerate(s.rivers); init=0.0) + reserve_economics(c,reserves).value
    x=Schedule(copy(c.grid), copy(u), z.storage, z.head, z.tunnel_q,
        z.tunnel_opening, Matrix{Float64}(q), z.power, Matrix{Float64}(gate),
        z.river_release, z.arrival_volume, z.transit[:,end], req.shortfall,
        objective, adjustment_cost, req.cost,
        Dict{String,Any}("case"=>c.name,"status"=>status,"solver"=>"independent_forward");reserves)
    x=with_soft_assessment(c,x)
    x.diagnostics["validation"]=validate(c, x; transport=nd)
    x
end

"""Correct bounded solver roundoff, then independently rebuild all physical values."""
function _reconstruct_candidate(c, raw; transport = nothing)
    raw=_schedule_input(raw; case=c)
    admissible(c, raw.commitment) || throw(ArgumentError("inadmissible commitment"))
    q, gate = copy(raw.generator_q), copy(raw.gate)
    all(isfinite, q) && all(isfinite, gate) || throw(ArgumentError("nonfinite controls"))
    all(x -> x >= -1e-6, q) || throw(ArgumentError("negative discharge"))
    all(x -> -1e-8 <= x <= 1 + 1e-8, gate) || throw(ArgumentError("gate outside bounds"))
    for k in eachindex(q, raw.commitment)
        raw.commitment[k]==0 || continue
        abs(q[k])<=1e-6 || throw(ArgumentError("nonzero discharge at an off unit"))
        q[k]=0.0
    end
    q = max.(q, 0.0)
    gate = clamp.(gate, 0.0, 1.0)
    raw_opening=_tunnel_opening_input(c, raw.tunnel_opening)
    ta=_physical_tunnel_openings(c, raw_opening)
    correction = (
        flow = maximum(abs, q .- raw.generator_q; init = 0.0),
        gate = maximum(abs, gate .- raw.gate; init = 0.0),
        tunnel_opening = maximum(
            abs,
            ta .- get(raw, "fractional_tunnel_opening", raw_opening);
            init = 0.0,
        ),
    )
    x = dispatch_from_controls(c, raw.commitment, q, gate; transport, tunnel_opening = ta, reserves=raw.reserves)
    x, correction
end

"""Use a warm start only when its hydraulic arrays match the requested grid."""
function _dispatch_start(c, u, warm)
    diagnostics=Dict{String, Any}("used"=>false)
    warm===nothing && return nothing, diagnostics
    if haskey(warm, "grid") && warm["grid"]!=c.grid
        diagnostics["ignored_reason"]="warm control grid differs from case"
        return nothing, diagnostics
    end
    haskey(warm, "V") || return nothing, diagnostics
    T=length(c.prices)
    for (key, shape) in (
        "generator_q"=>size(u),
        "power"=>size(u),
        "gate"=>(length(c.system.rivers), T),
        "river_release"=>(length(c.system.rivers), T),
        "V"=>(length(c.system.reservoirs), T+1),
        "H"=>(length(nodes(c.system)), T),
        "tunnel_q"=>(length(c.system.tunnels), T),
    )
        if !haskey(warm, key) || size(warm[key])!=shape || !all(isfinite, warm[key])
            diagnostics["ignored_reason"]="warm $key is missing, nonfinite or has an incompatible shape"
            return nothing, diagnostics
        end
    end
    opening=try
        _physical_tunnel_openings(c, get(warm, "tunnel_opening", nothing))
    catch e
        diagnostics["ignored_reason"]="warm tunnel_opening: "*sprint(showerror, e)
        return nothing, diagnostics
    end
    normalized=try
        warm isa Schedule && warm.tunnel_opening==opening ? warm : _schedule_input(warm;case=c)
    catch error
        diagnostics["ignored_reason"]="warm schedule: "*sprint(showerror,error)
        return nothing,diagnostics
    end
    diagnostics["used"]=true
    normalized, diagnostics
end

function dispatch_from_controls(
    c::ScheduleCase;
    commitment,
    generator_q,
    gate,
    status = "FORWARD_RECONSTRUCTED",
    transport = nothing,
    tunnel_opening = nothing,
    reserves = nothing,
)
    dispatch_from_controls(c, commitment, generator_q, gate; status, transport, tunnel_opening,reserves)
end
