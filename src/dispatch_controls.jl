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
    startup=sum(
        g.startup*max(0, u[j, t]-(t==1 ? g.initial_on : u[j, t - 1])) for
        (j, g) in enumerate(s.generators), t in 1:T;
        init = 0.0,
    )
    shutdown=sum(
        g.shutdown*max(0, (t==1 ? g.initial_on : u[j, t - 1])-u[j, t]) for
        (j, g) in enumerate(s.generators), t in 1:T;
        init = 0.0,
    )
    x["objective"]=sum(
                       c.prices[t]*dt[t]*x["power"][j, t] for
                       j in eachindex(s.generators), t in 1:T;
                       init = 0.0,
                   )-startup-shutdown-req.cost +
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
