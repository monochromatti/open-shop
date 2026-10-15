"""
    energy_conversion_factor(head, efficiency)

Electricity in MWh per Mm³ at a supplied reference net head (metres) and
combined turbine/electrical efficiency (fraction). The constant is consistent
with `P = 0.00981*q*h*η` MW and water volumes in Mm³.
This is a supplied-reference calculation, not a marginal value from a solve.
"""
function energy_conversion_factor(h::Real,η::Real)
    isfinite(h) && h>=0 && isfinite(η) && 0<η<=1 ||
        throw(ArgumentError("reference head must be finite and nonnegative; efficiency must lie in (0,1]"))
    factor=(0.00981/0.0036)*h*η
    isfinite(factor) || throw(ArgumentError("energy conversion factor overflow"))
    factor
end

"""
    volume_water_value(energy_value, conversion_factor; downstream_value=0)

Convert currency/MWh to currency/Mm³ using a supplied MWh/Mm³ factor. Add
`downstream_value` only when it describes additional future generation after
the conversion represented by this factor. All values must use one currency.
No head, downstream path, exchange rate or optimization dual is inferred.
"""
function volume_water_value(energy_value::Real,factor::Real;downstream_value::Real=0)
    all(isfinite,(energy_value,factor,downstream_value)) && factor>=0 ||
        throw(ArgumentError("water values must be finite and conversion factor nonnegative"))
    value=energy_value*factor+downstream_value
    isfinite(value) || throw(ArgumentError("volume water value overflow"))
    value
end

"""
    objective_components(case, schedule)

Recompute the economic components in the case's common currency. Reservoir
and river values are changes from initial storage and historical transit,
respectively. Recomputed costs use physical controls, not reported cost fields.
Use `validate` separately to establish schedule feasibility.
"""
function objective_components(c::ScheduleCase,x)
    x=_schedule_input(x;case=c)
    x.grid==c.grid || throw(ArgumentError("objective accounting requires the case control grid"))
    s=c.system;dt=diff(c.grid)
    initial=[remaining_volume(r,first(c.grid),r.history_grid,r.history_release)
        for r in s.rivers]
    exchange=net_power(c,x)
    energy_revenue=sum(c.prices[t]*dt[t]*exchange[t] for t in eachindex(dt);init=0.0)
    reservoir_value=sum(r.water_value*(x.storage[i,end]-r.v0)
        for (i,r) in enumerate(s.reservoirs);init=0.0)
    transit_value=sum(r.water_value*(x.terminal_transit[i]-initial[i])
        for (i,r) in enumerate(s.rivers);init=0.0)
    transition_cost=transition_costs(c,x.commitment)
    gate_cost=tunnel_adjustment_cost(c,x.tunnel_opening)
    release_penalty=release_requirements(c,x.river_release).cost
    soft_penalty=soft_constraint_assessment(c,x).cost
    reserve=reserve_economics(c,x.reserves)
    reserve_revenue,reserve_penalty=reserve.revenue,reserve.cost
    total=energy_revenue+reservoir_value+transit_value+reserve_revenue-
        transition_cost-gate_cost-release_penalty-soft_penalty-reserve_penalty
    (;energy_revenue,reservoir_value,transit_value,reserve_revenue,
        transition_cost,gate_cost,release_penalty,soft_penalty,reserve_penalty,total)
end
