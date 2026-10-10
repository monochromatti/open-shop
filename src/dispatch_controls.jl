"""Build a physical schedule from prescribed controls using independent integration.
No bound, gate or state is clipped. This is reconstruction, not optimization.
"""
function dispatch_from_controls(
    c,
    u,
    q,
    gate;
    status = "FORWARD_RECONSTRUCTED",
    transport = nothing,
)
    admissible(c, u) || throw(ArgumentError("inadmissible commitment"))
    nd=all(r->r.deterministic_delay!==nothing, c.system.rivers) ?
       _transport_data(c, transport) : nothing
    nd===nothing &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed reconstruction",
            ),
        )
    z=simulate(c, q, gate; transport = nd)
    z["converged"] || throw(ArgumentError("forward reconstruction did not converge"))
    s=c.system
    dt=diff(c.grid)
    T=length(dt)
    x=Dict{String,Any}(
        "case"=>c.name,
        "grid"=>copy(c.grid),
        "status"=>status,
        "solver"=>"independent_forward",
        "u"=>copy(u),
        "generator_q"=>copy(q),
        "gate"=>copy(gate),
        "V"=>z["V"],
        "H"=>z["H"],
        "tunnel_q"=>z["tunnel_q"],
        "power"=>z["power"],
        "river_release"=>z["river_release"],
        "arrival_volume"=>z["arrival_volume"],
        "terminal_transit"=>z["transit"][:, end],
    )
    req=release_requirements(c, x["river_release"])
    x["shortfall_release"]=req.shortfall
    x["release_penalty_cost"]=req.cost
    operating_cost=transition_costs(c,u)
    x["objective"]=sum(
                       c.prices[t]*dt[t]*x["power"][j, t] for
                       j in eachindex(s.generators), t in 1:T;
                       init = 0.0,
                   )-operating_cost-req.cost +
                   sum(
                       r.water_value*(x["V"][i, end]-r.v0) for
                       (i, r) in enumerate(s.reservoirs);
                       init = 0.0,
                   ) +
                   sum(
                       r.water_value*(z["transit"][i, end]-z["transit"][i, 1]) for
                       (i, r) in enumerate(s.rivers);
                       init = 0.0,
                   )
    x["validation"]=validate(c, x; transport = nd)
    x
end

"""Correct bounded solver roundoff, then independently rebuild all physical values."""
function _reconstruct_candidate(c, raw; transport = nothing)
    admissible(c, raw["u"]) || throw(ArgumentError("inadmissible commitment"))
    q, gate = copy(raw["generator_q"]), copy(raw["gate"])
    all(isfinite, q) && all(isfinite, gate) || throw(ArgumentError("nonfinite controls"))
    all(x -> x >= -1e-6, q) || throw(ArgumentError("negative discharge"))
    all(x -> -1e-8 <= x <= 1 + 1e-8, gate) || throw(ArgumentError("gate outside bounds"))
    for k in eachindex(q, raw["u"])
        raw["u"][k]==0 || continue
        abs(q[k])<=1e-6 || throw(ArgumentError("nonzero discharge at an off unit"))
        q[k]=0.0
    end
    q = max.(q, 0.0)
    gate = clamp.(gate, 0.0, 1.0)
    correction = (
        flow = maximum(abs, q .- raw["generator_q"]; init = 0.0),
        gate = maximum(abs, gate .- raw["gate"]; init = 0.0),
    )
    x = dispatch_from_controls(c, raw["u"], q, gate; transport)
    x, correction
end

"""Use a warm start only when its hydraulic arrays match the requested grid."""
function _dispatch_start(c, u, warm)
    diagnostics=Dict{String,Any}("used"=>false)
    warm===nothing && return nothing, diagnostics
    if haskey(warm, "grid") && warm["grid"]!=c.grid
        diagnostics["ignored_reason"]="warm control grid differs from case"
        return nothing, diagnostics
    end
    haskey(warm, "V") || return nothing, diagnostics
    T=length(c.prices)
    for (key, shape) in (
        "generator_q"=>size(u), "power"=>size(u),
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
    diagnostics["used"]=true
    warm, diagnostics
end
