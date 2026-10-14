"""Economic value of dispatch, in physical currency units.

Both optimization representations maximize revenue plus changes in stored and
in-transit water value, less unit transitions, named penalties and gate motion.
Physical replay evaluates these quantities independently.
"""
function dispatch_objective(model, c, P, V, terminal, u, shortfall;
    initial_transit, gate_cost, transitions=nothing, soft_cost=0.0)
    dt=diff(c.grid)
    s=c.system
    operating_cost=transition_costs(c,u;transitions)
    penalty=@expression(model,
        sum(opinterval(c,r.name,:release_penalty,t,0.0)*shortfall[i,t]
            for (i,r) in enumerate(s.rivers), t in eachindex(dt)))
    objective=@expression(model,
        sum(c.prices[t]*dt[t]*unit_power_sign(s.generators[j])*P[j,t]
            for j in axes(P,1), t in eachindex(dt)) +
        sum(r.water_value*(V[i,end]-r.v0) for (i,r) in enumerate(s.reservoirs)) +
        sum(r.water_value*(terminal[i]-initial_transit[i]) for (i,r) in enumerate(s.rivers)) -
        operating_cost-penalty-gate_cost-soft_cost)
    objective, penalty
end

"""
    net_power(case, power)
    net_power(case, schedule_or_trajectory)

Net electrical exchange in MW for each interval: production minus consumption.
Rows of `power` are nonnegative magnitudes in `case.system.generators` order,
the common axis used for generating and pumping modes.
"""
function net_power(c::ScheduleCase,P::AbstractMatrix)
    size(P,1)==length(c.system.generators) || throw(DimensionMismatch("unit power rows"))
    [sum(unit_power_sign(g)*P[j,t] for (j,g) in enumerate(c.system.generators);init=0.0)
        for t in axes(P,2)]
end
net_power(c::ScheduleCase,x::Union{Schedule,HydraulicTrajectory}) = net_power(c,x.power)
