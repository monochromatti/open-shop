const OPERATION_ATTRIBUTES=Dict(
    Reservoir=>Set((:inflow, :vmin, :vmax, :volume_ramp_up, :volume_ramp_down,
        :level_ramp_up, :level_ramp_down)),
    Generator=>Set((:qmin, :qmax, :pmin, :pmax, :forced_on, :maintenance,
        :power, :discharge, :startup, :shutdown, :ramp_up, :ramp_down,
        :discharge_ramp_up, :discharge_ramp_down)),
    Plant=>Set((:pmin, :pmax, :qmin, :qmax, :forced_on, :maintenance,
        :power, :discharge, :ramp_up, :ramp_down,
        :discharge_ramp_up, :discharge_ramp_down)),
    Tunnel=>Set((:capacity, :opening)),
    River=>Set((
        :inflow,
        :capacity,
        :gate_min,
        :gate_max,
        :min_arrival,
        :min_release,
        :release_penalty,
        :ramp_up,
        :ramp_down,
    )),
    FlowRequirement=>Set((:inflow, :min_flow)),
)

"""Optional interval input. Missing schedules or ramp limits remain disabled."""
function optional_opvalue(c, object, attribute, time, default=nothing)
    i=findfirst(z->z.object==object && z.attribute==attribute, c.operations)
    i===nothing && return default
    z=c.operations[i]
    k=searchsortedlast(z.times, time)
    k==0 ? default : z.values[k]
end
optional_operation(c, object, attribute, t; grid=c.grid, default=nothing) =
    optional_opvalue(c, object, attribute, grid[t], default)

"""Piecewise-constant operating data at absolute physical time. Before a series'
first knot the object's static default applies. Internal knots have both limits.
"""
function opvalue(c, object, attribute, time, default; side = :right)
    side in (:left, :right) ||
        throw(ArgumentError("operating-data side must be left or right"))
    isfinite(time) || throw(ArgumentError("operating-data time must be finite"))
    series=findfirst(z->z.object==object && z.attribute==attribute, c.operations)
    series===nothing && return Float64(default)
    z=c.operations[series]
    k=side==:right ? searchsortedlast(z.times, time) : searchsortedfirst(z.times, time)-1
    k<=0 ? Float64(default) : z.values[k]
end

function opaverage(c, object, attribute, a, b, default)
    b>a || throw(ArgumentError("positive operating-data averaging window required"))
    isfinite(a) && isfinite(b) ||
        throw(ArgumentError("operating-data window must be finite"))
    series=findfirst(z->z.object==object && z.attribute==attribute, c.operations)
    series===nothing && return Float64(default)
    z=c.operations[series]
    k=searchsortedlast(z.times, a)
    value=k==0 ? Float64(default) : z.values[k]
    total=0.0
    previous=a
    k+=1
    while k<=length(z.times) && z.times[k]<b
        time=z.times[k]
        total+=(time-previous)*value
        previous=time
        value=z.values[k]
        k+=1
    end
    (total+(b-previous)*value)/(b-a)
end
opinterval(c, object, attribute, t, default; grid = c.grid) =
    opaverage(c, object, attribute, grid[t], grid[t + 1], default)

"""Continuous storage satisfies both sides of an internal bound change."""
function opvertex(c, object, attribute, time, default; lower = true, grid = c.grid)
    time==first(grid) && return opvalue(c, object, attribute, time, default)
    time==last(grid) && return opvalue(c, object, attribute, time, default; side = :left)
    left=opvalue(c, object, attribute, time, default; side = :left)
    right=opvalue(c, object, attribute, time, default)
    lower ? max(left, right) : min(left, right)
end
storage_bounds(c, r, time; grid = c.grid) = (
    opvertex(c, r.name, :vmin, time, r.vmin; grid),
    opvertex(c, r.name, :vmax, time, r.vmax; lower = false, grid),
)

function validate_operations(c)
    objects=Dict(
        x.name=>x for xs in (
            c.system.reservoirs,
            c.system.generators,
            c.system.plants,
            c.system.tunnels,
            c.system.rivers,
            c.flow_requirements,
        ) for x in xs
    )
    seen=Set{Tuple{Symbol,Symbol}}()
    for z in c.operations
        haskey(objects, z.object) ||
            throw(ArgumentError("unknown operational object $(z.object)"))
        z.attribute in OPERATION_ATTRIBUTES[typeof(objects[z.object])] || throw(
            ArgumentError("unsupported operational attribute $(z.object).$(z.attribute)"),
        )
        key=(z.object, z.attribute)
        key in seen && throw(ArgumentError("duplicate operational series $key"))
        push!(seen, key)
        length(z.times)==length(z.values)>0 &&
        all(isfinite, z.times) &&
        all(isfinite, z.values) &&
        all(diff(z.times) .> 0) ||
            throw(ArgumentError("operational series needs ordered finite times and values"))
        for t in z.times
            first(c.grid)<t<last(c.grid) &&
                !any(isapprox(t, k; atol = 1e-12) for k in c.grid) &&
                throw(ArgumentError("scheduling grid misses operating change at $t"))
        end
        obj=objects[z.object]
        if z.attribute in (:vmin, :qmin, :pmin)
            all(v->v>=getfield(obj, z.attribute), z.values) ||
                throw(ArgumentError("operating minimum may only tighten its static bound"))
        elseif z.attribute in (:vmax, :qmax, :pmax, :capacity)
            ceiling=getfield(obj, z.attribute)
            if obj isa Plant && z.attribute==:qmax && ceiling===nothing
                ceiling=sum(g.qmax for g in c.system.generators if g.plant==obj.name; init=0.0)
            end
            all(v->v<=ceiling, z.values) ||
                throw(ArgumentError("operating maximum may only tighten its static bound"))
        end
        if z.attribute==:forced_on
            all(x->x in (-1.0, 0.0, 1.0), z.values) ||
                throw(ArgumentError("forced_on must be -1 (free), 0 (off), or 1 (on)"))
        elseif z.attribute==:maintenance
            all(x->x in (0.0, 1.0), z.values) ||
                throw(ArgumentError("maintenance must be 0 (available) or 1 (unavailable)"))
        elseif z.attribute in (:opening, :gate_min, :gate_max)
            all(x->0<=x<=1, z.values) ||
                throw(ArgumentError("opening and gate bounds must lie in [0,1]"))
        else
            all(>=(0), z.values) ||
                throw(ArgumentError("operational bounds/rates/costs must be nonnegative"))
        end
    end
    for r in c.system.reservoirs, t in c.grid
        opvertex(c, r.name, :vmin, t, r.vmin)<=opvertex(
            c,
            r.name,
            :vmax,
            t,
            r.vmax;
            lower = false,
        ) || throw(ArgumentError("inconsistent storage bounds for $(r.name)"))
    end
    for r in c.system.reservoirs
        lo, hi=storage_bounds(c, r, first(c.grid))
        lo<=r.v0<=hi || throw(ArgumentError("initial storage outside operative bounds"))
    end
    for g in c.system.generators, t in eachindex(c.prices)
        opinterval(c, g.name, :qmin, t, g.qmin)<=opinterval(c, g.name, :qmax, t, g.qmax) ||
            throw(ArgumentError("inconsistent unit flow bounds"))
        opinterval(c, g.name, :pmin, t, g.pmin)<=opinterval(c, g.name, :pmax, t, g.pmax) ||
            throw(ArgumentError("inconsistent unit power bounds"))
        for (attribute, lo, hi) in (
            (:power, opinterval(c, g.name, :pmin, t, g.pmin), opinterval(c, g.name, :pmax, t, g.pmax)),
            (:discharge, opinterval(c, g.name, :qmin, t, g.qmin), opinterval(c, g.name, :qmax, t, g.qmax)),
        )
            scheduled=optional_operation(c, g.name, attribute, t)
            scheduled===nothing || scheduled==0 || lo<=scheduled<=hi ||
                throw(ArgumentError("$(g.name).$attribute schedule outside running bounds"))
        end
    end
    for p in c.system.plants, t in eachindex(c.prices)
        qmax=p.qmax===nothing ? sum(g.qmax for g in c.system.generators if g.plant==p.name; init=0.0) : p.qmax
        for (attribute, lo, hi) in (
            (:power, opinterval(c, p.name, :pmin, t, p.pmin), opinterval(c, p.name, :pmax, t, p.pmax)),
            (:discharge, opinterval(c, p.name, :qmin, t, p.qmin), opinterval(c, p.name, :qmax, t, qmax)),
        )
            lo<=hi || throw(ArgumentError("inconsistent plant operating bounds"))
            scheduled=optional_operation(c, p.name, attribute, t)
            scheduled===nothing || scheduled==0 || lo<=scheduled<=hi ||
                throw(ArgumentError("$(p.name).$attribute schedule outside running bounds"))
        end
    end
    for r in c.system.rivers, t in eachindex(c.prices)
        opinterval(c, r.name, :gate_min, t, r.gate_min)<=opinterval(
            c,
            r.name,
            :gate_max,
            t,
            1.0,
        ) || throw(ArgumentError("inconsistent gate bounds"))
        r.law==:weir &&
            opinterval(c, r.name, :gate_max, t, 1.0)<1 &&
            throw(ArgumentError("uncontrolled weir gate must be one"))
        r.law==:junction &&
            opinterval(c, r.name, :gate_min, t, r.gate_min)>0 &&
            throw(ArgumentError("inferred confluence has no gate"))
        opinterval(c,r.name,:inflow,t,r.inflow)<=opinterval(c,r.name,:capacity,t,r.capacity) ||
            throw(ArgumentError("river inflow exceeds capacity of $(r.name)"))
        opinterval(c, r.name, :capacity, t, r.capacity)<=r.capacity+1e-12 || throw(
            ArgumentError(
                "river operating capacity cannot exceed routing calibration capacity",
            ),
        )
    end
    _commitment_bounds(c) # Intersect operational states with initial dwell obligations.
    nothing
end

function flow_requirement_indices(c, requirement)
    generators=Dict(g.name=>i for (i, g) in enumerate(c.system.generators))
    rivers=Dict(r.name=>i for (i, r) in enumerate(c.system.rivers))
    ([generators[n] for n in requirement.generators], [rivers[n] for n in requirement.rivers])
end

"""Add the same linear operating observation to local, proposal and global models."""
function constrain_flow_requirements!(m, c, generator_q, release)
    for rule in c.flow_requirements
        gs, rs=flow_requirement_indices(c, rule)
        for t in eachindex(c.prices)
            observed=sum(generator_q[i, t] for i in gs; init=0.0)+
                     sum(release[i, t] for i in rs; init=0.0)+
                     opinterval(c, rule.name, :inflow, t, rule.inflow)
            @constraint(m, observed>=opinterval(c, rule.name, :min_flow, t, rule.min_flow))
        end
    end
    nothing
end

function release_requirements(c, release; grid = c.grid)
    D=length(c.system.rivers)
    T=length(grid)-1
    size(release)==(D, T) || throw(DimensionMismatch("release requirements"))
    shortfall=zeros(D, T)
    cost=0.0
    hard=zeros(D, T)
    for (d, r) in enumerate(c.system.rivers), t in 1:T
        requirement=opinterval(c, r.name, :min_release, t, 0.0; grid)
        penalty=opinterval(c, r.name, :release_penalty, t, 0.0; grid)
        deficit=max(0.0, requirement-release[d, t])
        if penalty>0
            shortfall[d, t]=0.0036*(grid[t + 1]-grid[t])*deficit
            cost+=penalty*shortfall[d, t]
        else
            hard[d, t]=deficit
        end
    end
    (; shortfall, cost, hard)
end

function operational_arrival_knots(c, r, grid)
    sort!(
        unique!(
            vcat(
                arrival_knots(r, grid),
                [
                    t for z in c.operations if z.object==r.name && z.attribute==:min_arrival
                    for t in z.times if first(grid)<=t<=last(grid)
                ],
            ),
        ),
    )
end

function arrival_requirement_violation(c, r, grid, q)
    maximum(
        (
            opvalue(c, r.name, :min_arrival, t, r.min_arrival; side)-point_arrival(
                r,
                grid,
                q,
                t;
                side,
            ) for t in operational_arrival_knots(c, r, grid) for side in
            (t==first(grid) ? (:right,) : t==last(grid) ? (:left,) : (:left, :right))
        );
        init = 0.0,
    )
end

# Snapshot physical exogenous inputs for an independent hydraulic step.
# Routing calibration capacity is never changed by an operating restriction.
function operational_system(c, a, b)
    isempty(c.operations) && return c.system
    s=c.system
    reservoirs=[
        _river_replace(r; inflow = opaverage(c, r.name, :inflow, a, b, r.inflow)) for
        r in s.reservoirs
    ]
    tunnels=[
        _river_replace(e; opening = opaverage(c, e.name, :opening, a, b, e.opening)) for
        e in s.tunnels
    ]
    rivers=[
        _river_replace(r; inflow=opaverage(c,r.name,:inflow,a,b,r.inflow),
            capacity=r.law==:controlled ? opaverage(c,r.name,:capacity,a,b,r.capacity) : r.capacity)
        for r in s.rivers
    ]
    _river_replace(s; reservoirs, tunnels, rivers)
end
function river_law_value(r, h, a)
    r.discharge_curve!==nothing &&
        return a*table_value(r.discharge_curve, h; extrapolation = :linear)
    r.law==:controlled && return max(0.0,r.capacity-r.inflow)*a
    r.law==:orifice && return r.coefficient*a*sqrt(max(h-r.crest, 0.0))
    r.law==:weir && return r.coefficient*max(h-r.crest, 0.0)^1.5
    throw(ArgumentError("unsupported river law"))
end
