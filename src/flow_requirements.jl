"""Observed interval rates. Arrivals are supplied as volumes in Mm³;
tunnel rates retain their source-to-target sign. No additional routing occurs.
"""
function observed_flow(c,rule,generator_q,release,arrivals,tunnel_q)
    gs=[something(findfirst(g->g.name==name,c.system.generators)) for name in rule.generators]
    rs=[something(findfirst(r->r.name==name,c.system.rivers)) for name in rule.rivers]
    ds=[something(findfirst(r->r.name==name,c.system.rivers)) for name in rule.river_arrivals]
    es=[something(findfirst(e->e.name==name,c.system.tunnels)) for name in rule.tunnels]
    [sum(generator_q[i,t] for i in gs;init=0.0)+
     sum(release[i,t] for i in rs;init=0.0)+
     sum(arrivals[i,t] for i in ds;init=0.0)/(0.0036*(c.grid[t+1]-c.grid[t]))+
     sum(tunnel_q[i,t] for i in es;init=0.0)+
     opinterval(c,rule.name,:inflow,t,rule.inflow) for t in eachindex(c.prices)]
end

function validate_flow_requirements(c)
    for rule in c.flow_requirements
        all(v->v===nothing || isfinite(v),(rule.inflow,rule.min_flow,rule.max_flow,
            rule.min_volume,rule.max_volume,rule.initial_volume)) && rule.inflow>=0 ||
            throw(ArgumentError("flow observation bounds must be finite"))
        rule.min_volume===nothing || rule.max_volume===nothing || rule.min_volume<=rule.max_volume ||
            throw(ArgumentError("inconsistent volume allowance for $(rule.name)"))
        for (names,objects) in ((rule.generators,c.system.generators),
            (rule.rivers,c.system.rivers),(rule.river_arrivals,c.system.rivers),
            (rule.tunnels,c.system.tunnels))
            length(unique(names))==length(names) ||
                throw(ArgumentError("duplicate flow observation contribution"))
            all(name->any(o->o.name==name,objects),names) ||
                throw(ArgumentError("unknown flow observation contribution"))
        end
        for t in eachindex(c.prices)
            lo=opinterval(c,rule.name,:min_flow,t,rule.min_flow)
            hi=optional_operation(c,rule.name,:max_flow,t;default=rule.max_flow)
            hi===nothing || lo<=hi || throw(ArgumentError("inconsistent flow bounds for $(rule.name)"))
        end
    end
    nothing
end

function constrain_flow_requirements!(m,c,generator_q,release,arrivals,tunnel_q)
    for rule in c.flow_requirements
        observed=observed_flow(c,rule,generator_q,release,arrivals,tunnel_q)
        for t in eachindex(c.prices)
            lo=opinterval(c,rule.name,:min_flow,t,rule.min_flow)
            hi=optional_operation(c,rule.name,:max_flow,t;default=rule.max_flow)
            if hi==lo
                @constraint(m,observed[t]==lo)
            else
                @constraint(m,observed[t]>=lo)
                hi===nothing || @constraint(m,observed[t]<=hi)
            end
        end
        volume=rule.initial_volume+sum(0.0036*(c.grid[t+1]-c.grid[t])*observed[t]
            for t in eachindex(c.prices);init=0.0)
        if rule.min_volume!==nothing && rule.min_volume==rule.max_volume
            @constraint(m,volume==rule.min_volume)
        else
            rule.min_volume===nothing || @constraint(m,volume>=rule.min_volume)
            rule.max_volume===nothing || @constraint(m,volume<=rule.max_volume)
        end
    end
    nothing
end

function flow_requirement_residuals(c,x)
    result=Dict{String,Float64}()
    for rule in c.flow_requirements
        observed=observed_flow(c,rule,x.generator_q,x.river_release,x.arrival_volume,x.tunnel_q)
        violation=0.0
        for t in eachindex(c.prices)
            violation=max(violation,opinterval(c,rule.name,:min_flow,t,rule.min_flow)-observed[t])
            hi=optional_operation(c,rule.name,:max_flow,t;default=rule.max_flow)
            hi===nothing || (violation=max(violation,observed[t]-hi))
        end
        result["flow_requirement_$(rule.name)"]=violation
        volume=rule.initial_volume+sum(0.0036*(c.grid[t+1]-c.grid[t])*observed[t]
            for t in eachindex(c.prices);init=0.0)
        result["flow_volume_$(rule.name)"]=max(0.0,
            rule.min_volume===nothing ? 0.0 : rule.min_volume-volume,
            rule.max_volume===nothing ? 0.0 : volume-rule.max_volume)
    end
    result
end
