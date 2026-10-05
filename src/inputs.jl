function validate_inputs(c)
    s=c.system
    for objects in (
        s.reservoirs,
        s.junctions,
        s.boundaries,
        s.tunnels,
        s.plants,
        s.generators,
        s.river_junctions,
        s.rivers,
    )
        length(unique(x.name for x in objects))==length(objects) ||
            error("Duplicate object name")
        for x in objects, f in fieldnames(typeof(x))
            value=getfield(x, f)
            value isa Real && !isfinite(value) && error("Nonfinite attribute $(x.name).$f")
        end
    end
    objectnames=[
        x.name for xs in (
            s.reservoirs,
            s.junctions,
            s.boundaries,
            s.tunnels,
            s.plants,
            s.generators,
            s.river_junctions,
            s.rivers,
        ) for x in xs
    ]
    length(unique(objectnames))==length(objectnames) ||
        error("Names must be globally unique")
    ns=nodes(s)
    allnames=vcat(ns, [j.name for j in s.river_junctions])
    length(unique(allnames))==length(allnames)||error("Duplicate node")
    length(unique([g.name for g in s.generators]))==length(s.generators)||error("Duplicate generator")
    length(unique([p.name for p in s.plants]))==length(s.plants)||error("Duplicate plant")
    length(c.prices)==length(c.grid)-1 && all(isfinite, c.prices)||error("Prices")
    RiverRouting.check_grid(c.grid)
    all(diff(c.grid) .> 0)||error("Grid")
    for r in s.reservoirs
        r.vmin<=r.v0<=r.vmax && r.vmin>=0 && r.vmax>r.vmin || error("Storage bounds")
        r.level_curve!==nothing ||
            min(r.slope+2r.curvature*r.vmin, r.slope+2r.curvature*r.vmax)>0||error(
                "Nonmonotone head storage",
            )
    end
    for e in s.tunnels
        e.source in ns && e.target in ns && e.source!=e.target || error("Tunnel connection")
        e.resistance>0 && e.capacity>=0 && 0<=e.opening<=1 || error("Tunnel limits")
    end
    for p in s.plants
        p.source in ns && p.target in ns || error("Plant connection")
        p.pmax>0 && p.ramp>=0 && p.initial_interval_hours>0 || error("Plant limits")
        p.initial_power===nothing ||
            0<=p.initial_power<=p.pmax ||
            error("Initial plant production")
        p.outlet_head_floor===nothing ||
            isfinite(p.outlet_head_floor) ||
            error("Nonfinite plant outlet head floor")
    end
    for j in s.junctions
        j.hmin<=j.hmax || error("Junction limits")
    end
    for g in s.generators
        plantof(s, g)
        0<g.qmin<g.qmax && 0<g.pmin<g.pmax || error("Unit limits")
        g.qbest>0 && g.hbest>0 && g.hmin<=g.hmax && g.startup>=0 && g.shutdown>=0 ||
            error("Unit operating domain")
        0<=g.min_efficiency<=1 || error("Unit minimum efficiency must lie in [0,1]")
        g.initial_on in [0, 1] && g.initial_age>=0 && g.minup>=0 && g.mindown>=0 ||
            error("Initial commitment")
    end
    for r in s.rivers
        r.capacity>0 && r.min_arrival>=0 && 0<=r.gate_min<=1 || error("River limits")
        r.law in (:orifice, :weir) &&
            r.discharge_curve===nothing &&
            r.coefficient<=0 &&
            error("River law coefficient")
        r.law==:junction &&
            !any(x.name==r.source for x in s.river_junctions) &&
            error("Junction river source must be declared river junction")
        r.law!=:junction &&
            !any(x.name==r.source for x in s.reservoirs) &&
            error("Operational river source must be reservoir")
        r.arrival_policy in (:interval_average, :pointwise) ||
            error("Unknown environmental averaging policy")
        if !isempty(r.arrival_window_grid)
            RiverRouting.check_grid(r.arrival_window_grid)
            first(r.arrival_window_grid)==first(c.grid) &&
            last(r.arrival_window_grid)==last(c.grid) ||
                error("Arrival windows must span scheduling horizon")
        end
        if r.deterministic_delay===nothing
            length(r.curves)==2 &&
            r.curves[1].reference_flow==0 &&
            r.curves[2].reference_flow==r.capacity || error("Delay curves")
        else
            isfinite(r.deterministic_delay) && r.deterministic_delay>=0 ||
                error("Invalid deterministic delay")
        end
        r.law in (:orifice, :weir, :junction, :controlled) || error("Unknown river law")
        if r.law in (:orifice, :weir) && !r.allow_dry && r.discharge_curve===nothing
            src=only(x for x in s.reservoirs if x.name==r.source)
            head(src, src.vmin)>r.crest ||
                error("Only always-wet outlet laws currently supported")
        end
    end
    validate_curve_data(c.system)
    validate_operations(c)
    # Exact transport does not need dense transfer tensors for input checks.
    all(r->r.deterministic_delay!==nothing, s.rivers) ? river_order(c) : routing_data(c)
    true
end
jsonready(x::AbstractMatrix) = [jsonready(collect(row)) for row in eachrow(x)]
jsonready(x::AbstractArray) = [jsonready(v) for v in x]
jsonready(x::AbstractDict) = Dict(string(k)=>jsonready(v) for (k, v) in x)
jsonready(x::AbstractFloat) = isfinite(x) ? x : nothing
jsonready(x::Symbol) = string(x)
jsonready(x) = x
function writejson(path, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.pretty(io, jsonready(obj))
    end
end

"""Refine decision intervals while preserving operating times and averaging windows."""
function with_grid(c::ScheduleCase, newgrid)
    grid=Float64.(collect(newgrid))
    RiverRouting.check_grid(grid)
    first(grid)==first(c.grid) && last(grid)==last(c.grid) ||
        throw(ArgumentError("new grid must span the same horizon"))
    all(t in grid for t in c.grid) ||
        throw(ArgumentError("refinement must retain every original control boundary"))
    prices=[
        c.prices[clamp(
            searchsortedlast(c.grid, (grid[t]+grid[t + 1])/2),
            1,
            length(c.prices),
        )] for t in 1:(length(grid) - 1)
    ]
    rivers=[
        r.arrival_policy==:interval_average && isempty(r.arrival_window_grid) ?
        _river_replace(r; arrival_window_grid = copy(c.grid)) : r for r in c.system.rivers
    ]
    _river_replace(c; grid, prices, system = _river_replace(c.system; rivers))
end
