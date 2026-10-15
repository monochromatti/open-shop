# Effective mathematical origins, not a history of constructor keyword spelling.
function penalty_series_origin(c,object,attribute,a,b;default)
    i=findfirst(z->z.object==object && z.attribute==attribute,c.operations)
    i===nothing && return default
    start=first(c.operations[i].times)
    start<=a ? :operational_series : start<b ? :default_and_series : default
end
soft_penalty_origin(bound,allowance) = bound===nothing ? :absent :
    allowance==0 ? :hard_limit : :explicit_penalty
soft_coefficient_unit(quantity) = quantity==:volume ? :currency_per_Mm3 :
    quantity==:level ? :currency_per_metre : quantity==:power ? :currency_per_MWh :
    :currency_per_Mm3

"""
    penalty_settings(case)

Effective operating and terminal penalty equations, without inspecting a result.
`soft` retains each original `SoftConstraint`, its physical unit, coefficient
unit, cost factor, activity, and each side's mathematical origin. `release` and
`reserve` give effective interval requirements, allowed shortfalls, coefficients
and their origins. A cost equals coefficient × cost_factor × physical shortfall;
soft rules have separate lower and upper terms. Inactive soft rules accrue zero.

The `defaults` tuple states OpenSHOP defaults. Positive soft allowances require
explicit positive costs; zero allowances remain hard. Normalized unused zero
soft coefficients do not retain whether their keyword was omitted or supplied.
`nothing` reserve shortfall cost is hard. These are OpenSHOP policies, not SHOP
default penalties. Currency is the same base currency as `case.prices`.
"""
function penalty_settings(c::ScheduleCase)
    RiverRouting.check_grid(c.grid)
    length(c.prices)==length(c.grid)-1 || throw(DimensionMismatch("prices and grid"))
    validate_soft_constraints(c)
    validate_reserve_inputs(c)
    validate_operations(c)
    defaults=(soft_allowance=0.0,soft_coefficient=0.0,river_minimum=0.0,
        river_release_coefficient=0.0,reserve_shortfall_coefficient=nothing)
    soft=[(;rule,active=soft_constraint_active(c,rule),unit=soft_slack_unit(rule.quantity),
        coefficient_unit=soft_coefficient_unit(rule.quantity),cost_factor=soft_cost_factor(rule),
        lower_origin=soft_penalty_origin(rule.lower,rule.max_lower_slack),
        upper_origin=soft_penalty_origin(rule.upper,rule.max_upper_slack))
        for rule in c.soft_constraints]
    release=[begin
        a,b=c.grid[t],c.grid[t+1]
        minimum=opinterval(c,r.name,:min_release,t,0.0)
        coefficient=opinterval(c,r.name,:release_penalty,t,0.0)
        (;object=r.name,time=a,end_time=b,minimum,coefficient,
            minimum_origin=penalty_series_origin(c,r.name,:min_release,a,b;default=:default_zero),
            coefficient_origin=penalty_series_origin(c,r.name,:release_penalty,a,b;default=:default_zero),
            max_shortfall=coefficient>0 ? minimum : 0.0,unit=:m3_per_s,
            coefficient_unit=:currency_per_Mm3,cost_factor=0.0036*(b-a))
    end for r in c.system.rivers for t in eachindex(c.prices) if !is_channel(r)]
    reserve=[begin
        a,b=c.grid[t],c.grid[t+1]
        up_requirement=opinterval(c,rule.name,:up_requirement,t,rule.up_requirement)
        down_requirement=opinterval(c,rule.name,:down_requirement,t,rule.down_requirement)
        coefficient=rule.shortfall_cost
        (;rule,time=a,end_time=b,up_requirement,down_requirement,coefficient,
            up_requirement_origin=penalty_series_origin(c,rule.name,:up_requirement,a,b;default=:group_value),
            down_requirement_origin=penalty_series_origin(c,rule.name,:down_requirement,a,b;default=:group_value),
            coefficient_origin=coefficient===nothing ? :default_hard : :explicit_penalty,
            max_up_shortfall=coefficient===nothing ? 0.0 : up_requirement,
            max_down_shortfall=coefficient===nothing ? 0.0 : down_requirement,
            unit=:MW,coefficient_unit=:currency_per_MWh,cost_factor=b-a)
    end for rule in c.reserve_groups for t in eachindex(c.prices)]
    (;defaults,soft,release,reserve)
end
