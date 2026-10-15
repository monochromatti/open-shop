# Slacks change named operating restrictions, never conservation or a physical law.
soft_slack_unit(quantity) = quantity==:volume ? :Mm3 : quantity==:level ? :metres :
    quantity==:power ? :MW : :m3_per_s
soft_cost_factor(rule) = rule.end_time===nothing ? 1.0 :
    (rule.end_time-rule.time)*(rule.quantity==:discharge ? 0.0036 : 1.0)

function soft_constraint_active(c,rule)
    lo,hi=first(c.grid),last(c.grid)
    rule.end_time===nothing && return lo<=rule.time<=hi
    rule.end_time<=lo || rule.time>=hi ? false : begin
        lo<=rule.time && rule.end_time<=hi ||
            throw(ArgumentError("soft constraint $(rule.name) crosses a horizon boundary"))
        true
    end
end

function validate_soft_constraints(c)
    isempty(c.soft_constraints) && return nothing
    RiverRouting.check_grid(c.grid)
    s=c.system
    objects=Dict(x.name=>x for field in fieldnames(HydroSystem) for x in getfield(s,field))
    names=Set(vcat(collect(keys(objects)),[x.name for x in c.flow_requirements]))
    for rule in c.soft_constraints
        rule.name in names && throw(ArgumentError("duplicate soft constraint name $(rule.name)"))
        push!(names,rule.name)
        object=get(objects,rule.object,nothing)
        state=rule.quantity in (:volume,:level)
        rate=rule.quantity in (:power,:discharge)
        state && object isa Reservoir || rate && object isa Union{Generator,Plant} ||
            throw(ArgumentError("$(rule.name): unsupported soft constraint quantity or object"))
        isfinite(rule.time) || throw(ArgumentError("soft constraint time must be finite"))
        state ? rule.end_time===nothing || throw(ArgumentError("state constraints use one assessment time")) :
            rule.end_time!==nothing && isfinite(rule.end_time) && rule.end_time>rule.time &&
            isfinite(rule.end_time-rule.time) ||
                throw(ArgumentError("rate constraints need a finite positive assessment window"))
        rule.lower!==nothing || rule.upper!==nothing ||
            throw(ArgumentError("soft constraint needs a lower or upper limit"))
        for (bound,allowance,cost) in ((rule.lower,rule.max_lower_slack,rule.lower_cost),
            (rule.upper,rule.max_upper_slack,rule.upper_cost))
            isfinite(allowance) && allowance>=0 && isfinite(cost) && cost>=0 ||
                throw(ArgumentError("slack allowances and costs must be finite and nonnegative"))
            if bound===nothing
                allowance==0 && cost==0 || throw(ArgumentError("absent limits cannot have slack or cost"))
            else
                isfinite(bound) && (rule.quantity==:level || bound>=0) ||
                    throw(ArgumentError("soft limits must be finite and in the quantity's physical range"))
                allowance==0 || cost>0 || throw(ArgumentError("permitted slack needs a positive penalty cost"))
            end
        end
        rule.lower===nothing || rule.upper===nothing || rule.lower<=rule.upper ||
            throw(ArgumentError("soft constraint lower limit exceeds upper limit"))
        factor=soft_cost_factor(rule)
        all(isfinite(factor*price) && isfinite(factor*price*allowance)
            for (price,allowance) in ((rule.lower_cost,rule.max_lower_slack),
                (rule.upper_cost,rule.max_upper_slack))) ||
            throw(ArgumentError("soft penalty coefficients and bounded costs must be finite"))
        if soft_constraint_active(c,rule)
            rule.time in c.grid && (state || rule.end_time in c.grid) ||
                throw(ArgumentError("grid must include soft assessment times"))
        end
    end
    nothing
end

# A bounded slack implies this relaxed domain, not the nominal target domain.
# Legacy operating bounds and static reservoir safety limits remain hard.
function soft_storage_bounds(c,reservoir,time)
    lo,hi=storage_bounds(c,reservoir,time)
    for rule in c.soft_constraints
        rule.object==reservoir.name && rule.quantity==:volume && rule.time==time || continue
        rule.lower===nothing || (lo=max(lo,rule.lower-rule.max_lower_slack))
        rule.upper===nothing || (hi=min(hi,rule.upper+rule.max_upper_slack))
    end
    lo,hi
end

function soft_constraint_value(c,rule,V,P,GQ;level)
    if rule.end_time===nothing
        i=findfirst(r->r.name==rule.object,c.system.reservoirs)
        k=findfirst(==(rule.time),c.grid)
        return rule.quantity==:volume ? V[i,k] : level(i,k)
    end
    unit=findfirst(g->g.name==rule.object,c.system.generators)
    members=unit===nothing ? findall(g->g.plant==rule.object,c.system.generators) : [unit]
    rates=rule.quantity==:power ? P : GQ
    sum((c.grid[k+1]-c.grid[k])*sum(rates[j,k] for j in members;init=0.0)
        for k in eachindex(c.prices) if rule.time<=c.grid[k] && c.grid[k+1]<=rule.end_time;
        init=0.0)/(rule.end_time-rule.time)
end

"""Canonical minimal slacks, explicit positive costs, and excess beyond allowances.

Reservoir quantities are evaluated at their original state vertices. Unit or
plant outputs are averaged over the original assessment windows, even when a
schedule is finer than those windows. The returned costs are always subtracted
from economic value; positive permitted deviations are not physical failures.
"""
function soft_constraint_assessment(c,x)
    validate_soft_constraints(c)
    x=_schedule_input(x;case=c)
    x.grid==c.grid || throw(ArgumentError("soft assessment schedule grid differs from case"))
    slacks=SoftConstraintSlack[]
    excess=Float64[]
    cost=0.0
    level=(i,k)->head(c.system.reservoirs[i],x.storage[i,k])
    for rule in c.soft_constraints
        active=soft_constraint_active(c,rule)
        value=active ? soft_constraint_value(c,rule,x.storage,x.power,x.generator_q;level) : 0.0
        lower=active && rule.lower!==nothing ? max(0.0,rule.lower-value) : 0.0
        upper=active && rule.upper!==nothing ? max(0.0,value-rule.upper) : 0.0
        penalty=soft_cost_factor(rule)*(rule.lower_cost*lower+rule.upper_cost*upper)
        push!(slacks,SoftConstraintSlack(rule.name,lower,upper,soft_slack_unit(rule.quantity),penalty))
        push!(excess,max(0.0,lower-rule.max_lower_slack,upper-rule.max_upper_slack))
        cost+=penalty
    end
    (;slacks,cost,excess)
end

function constrain_soft_constraints!(m,c,V,P,GQ;level)
    slacks=Tuple{Any,Any}[]
    cost=AffExpr(0.0)
    for (i,rule) in enumerate(c.soft_constraints)
        if !soft_constraint_active(c,rule)
            push!(slacks,(0.0,0.0))
            continue
        end
        lower=rule.max_lower_slack==0 ? 0.0 : @variable(m,
            lower_bound=0,upper_bound=rule.max_lower_slack,base_name="soft_lower_$(i)")
        upper=rule.max_upper_slack==0 ? 0.0 : @variable(m,
            lower_bound=0,upper_bound=rule.max_upper_slack,base_name="soft_upper_$(i)")
        value=soft_constraint_value(c,rule,V,P,GQ;level)
        rule.lower===nothing || @constraint(m,value+lower>=rule.lower)
        rule.upper===nothing || @constraint(m,value-upper<=rule.upper)
        factor=soft_cost_factor(rule)
        lower isa Real || add_to_expression!(cost,factor*rule.lower_cost,lower)
        upper isa Real || add_to_expression!(cost,factor*rule.upper_cost,upper)
        push!(slacks,(lower,upper))
    end
    (;slacks,cost)
end

function with_soft_assessment(c,x::Schedule;included_cost=0.0)
    isempty(c.soft_constraints) && return x
    assessment=soft_constraint_assessment(c,x)
    Schedule(x.grid,x.commitment,x.storage,x.head,x.tunnel_q,x.tunnel_opening,x.generator_q,
        x.power,x.gate,x.river_release,x.arrival_volume,x.terminal_transit,x.shortfall_release,
        x.objective+included_cost-assessment.cost,x.gate_adjustment_cost,x.release_penalty_cost,
        assessment.slacks,assessment.cost,x.diagnostics)
end
