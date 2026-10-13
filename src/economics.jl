"""Economic value of dispatch, in physical currency units.

Both optimization representations maximize revenue plus changes in stored and
in-transit water value, less unit transitions, release penalties and gate motion.
Physical replay evaluates these quantities independently.
"""
function dispatch_objective(model, c, P, V, terminal, u, shortfall;
    initial_transit, gate_cost, transitions=nothing)
    dt=diff(c.grid)
    s=c.system
    operating_cost=transition_costs(c,u;transitions)
    penalty=@expression(model,
        sum(opinterval(c,r.name,:release_penalty,t,0.0)*shortfall[i,t]
            for (i,r) in enumerate(s.rivers), t in eachindex(dt)))
    objective=@expression(model,
        sum(c.prices[t]*dt[t]*P[j,t] for j in axes(P,1), t in eachindex(dt)) +
        sum(r.water_value*(V[i,end]-r.v0) for (i,r) in enumerate(s.reservoirs)) +
        sum(r.water_value*(terminal[i]-initial_transit[i]) for (i,r) in enumerate(s.rivers)) -
        operating_cost-penalty-gate_cost)
    objective, penalty
end
