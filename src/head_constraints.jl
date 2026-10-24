# Gross head uses named reservoir surfaces or fixed boundary levels, before
# tunnel and machine losses. Startup is assessed at the control edge; running
# limits use the same volume-averaged reservoir levels as dispatch.
const GROSS_HEAD_ATTRIBUTES = (:gross_head_min, :gross_head_max, :startup_gross_head_min)

function gross_head_references(s, p)
    names=(something(p.gross_head_source,p.source),something(p.gross_head_target,p.target))
    map(names) do name
        i=findfirst(r->r.name==name,s.reservoirs)
        i!==nothing && return (;reservoir=i,level=0.0)
        b=findfirst(b->b.name==name,s.boundaries)
        b!==nothing && return (;reservoir=0,level=s.boundaries[b].head)
        throw(ArgumentError("$(p.name): gross head needs explicit reservoir or boundary references"))
    end
end

function validate_gross_heads(c)
    for p in c.system.plants
        if p.gross_head_source!==nothing || p.gross_head_target!==nothing
            gross_head_references(c.system,p)
        end
    end
    for g in c.system.generators
        enabled=any(a->getfield(g,a)!==nothing,GROSS_HEAD_ATTRIBUTES) ||
            any(z->z.object==g.name && z.attribute in GROSS_HEAD_ATTRIBUTES,c.operations)
        enabled || continue
        gross_head_references(c.system,plantof(c.system,g))
        for a in GROSS_HEAD_ATTRIBUTES
            value=getfield(g,a)
            value===nothing || isfinite(value) && value>=0 ||
                throw(ArgumentError("$(g.name).$a must be finite and nonnegative"))
        end
        for t in eachindex(c.prices)
            lo=optional_operation(c,g.name,:gross_head_min,t;default=g.gross_head_min)
            hi=optional_operation(c,g.name,:gross_head_max,t;default=g.gross_head_max)
            lo===nothing || hi===nothing || lo<=hi ||
                throw(ArgumentError("inconsistent gross-head limits for $(g.name)"))
        end
    end
    nothing
end

reference_level(ref,level,t)=ref.reservoir==0 ? ref.level : level(ref.reservoir,t)
function gross_head_bounds(s,g,refs)
    bounds(ref)=ref.reservoir==0 ? (ref.level,ref.level) :
        (head(s.reservoirs[ref.reservoir],s.reservoirs[ref.reservoir].vmin),
         head(s.reservoirs[ref.reservoir],s.reservoirs[ref.reservoir].vmax))
    a,b=bounds.(refs)
    unit_power_sign(g)>0 ? (a[1]-b[2],a[2]-b[1]) : (b[1]-a[2],b[2]-a[1])
end

function constrain_gross_heads!(m,c,u,interval_level,vertex_level;transitions=nothing)
    for (i,g) in enumerate(c.system.generators)
        any(a->getfield(g,a)!==nothing,GROSS_HEAD_ATTRIBUTES) ||
            any(z->z.object==g.name && z.attribute in GROSS_HEAD_ATTRIBUTES,c.operations) || continue
        refs=gross_head_references(c.system,plantof(c.system,g))
        lower,upper=gross_head_bounds(c.system,g,refs)
        for t in eachindex(c.prices)
            gross=unit_power_sign(g)*(reference_level(refs[1],interval_level,t)-
                reference_level(refs[2],interval_level,t))
            lo=optional_operation(c,g.name,:gross_head_min,t;default=g.gross_head_min)
            hi=optional_operation(c,g.name,:gross_head_max,t;default=g.gross_head_max)
            lo===nothing || @constraint(m,gross/250>=(lo-max(0.0,lo-lower)*(1-u[i,t]))/250)
            hi===nothing || @constraint(m,gross/250<=(hi+max(0.0,upper-hi)*(1-u[i,t]))/250)
            start_min=optional_operation(c,g.name,:startup_gross_head_min,t;default=g.startup_gross_head_min)
            start_min===nothing && continue
            # u_t-u_(t-1) is one exactly at a startup; no new binary is needed.
            prior=t==1 ? g.initial_on : u[i,t-1]
            start=transitions===nothing ? u[i,t]-prior : transitions.su[i,t]
            at_start=unit_power_sign(g)*(reference_level(refs[1],vertex_level,t)-
                reference_level(refs[2],vertex_level,t))
            @constraint(m,at_start/250>=(start_min-max(0.0,start_min-lower)*(1-start))/250)
        end
    end
    nothing
end

function gross_head_residuals(c,x)
    result=Dict{String,Float64}()
    for (i,g) in enumerate(c.system.generators)
        any(a->getfield(g,a)!==nothing,GROSS_HEAD_ATTRIBUTES) ||
            any(z->z.object==g.name && z.attribute in GROSS_HEAD_ATTRIBUTES,c.operations) || continue
        refs=gross_head_references(c.system,plantof(c.system,g))
        interval_level=(j,t)->mean_head(c.system.reservoirs[j],x.storage[j,t],x.storage[j,t+1])
        vertex_level=(j,t)->head(c.system.reservoirs[j],x.storage[j,t])
        violation=0.0
        for t in eachindex(c.prices)
            x.commitment[i,t]==1 || continue
            gross=unit_power_sign(g)*(reference_level(refs[1],interval_level,t)-reference_level(refs[2],interval_level,t))
            lo=optional_operation(c,g.name,:gross_head_min,t;default=g.gross_head_min)
            hi=optional_operation(c,g.name,:gross_head_max,t;default=g.gross_head_max)
            lo===nothing || (violation=max(violation,lo-gross))
            hi===nothing || (violation=max(violation,gross-hi))
            prior=t==1 ? g.initial_on : x.commitment[i,t-1]
            startup=optional_operation(c,g.name,:startup_gross_head_min,t;default=g.startup_gross_head_min)
            if prior==0 && startup!==nothing
                gross=unit_power_sign(g)*(reference_level(refs[1],vertex_level,t)-reference_level(refs[2],vertex_level,t))
                violation=max(violation,startup-gross)
            end
        end
        result["gross_head_$(g.name)"]=violation
    end
    result
end
