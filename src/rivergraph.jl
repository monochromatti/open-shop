const _INFERRED_RIVER_PREFIX = "__river_merge__"

# Keep this independent of restart/serialization so graph compilation remains
# available during input construction.
function _river_replace(x; kwargs...)
    values=Dict(f=>getfield(x, f) for f in fieldnames(typeof(x)))
    merge!(values, Dict(kwargs))
    typeof(x)(; values...)
end

"""Compile direct reach-to-reach targets to private zero-storage graph nodes.

An upstream river may name another river as its target; plant/tunnel outfalls
name their discharge_river. The receiving reach uses source=:auto; incoming
arrivals, directed outfalls and natural inflow are merged by conservation. Explicit
RiverJunction nodes remain supported. Inferred names are reserved and
collisions are rejected. This function does not infer a diversion allocation.
"""
function normalize_river_connections(s::HydroSystem)
    names=Set(r.name for r in s.rivers)
    length(names)==length(s.rivers) || throw(ArgumentError("duplicate reach names"))
    physical=Set(vcat(nodes(s), [j.name for j in s.river_junctions]))
    isempty(intersect(names, physical)) ||
        throw(ArgumentError("river names and node names must be disjoint"))
    destinations=Set(r.target for r in s.rivers if r.target in names)
    for x in Iterators.flatten((s.plants,s.tunnels))
        x.discharge_river===nothing && continue
        x.discharge_river in names || throw(ArgumentError("unknown discharge river for $(x.name)"))
        push!(destinations,x.discharge_river)
    end
    # An inflow-only reach also has an internal zero-storage source.
    union!(destinations, (r.name for r in s.rivers if r.source==:auto))
    isempty(destinations) && return s
    inferred=Dict(n=>Symbol(_INFERRED_RIVER_PREFIX, string(n)) for n in destinations)
    existing=Set(j.name for j in s.river_junctions if startswith(string(j.name),_INFERRED_RIVER_PREFIX))
    reused=Set(inferred[r.name] for r in s.rivers if r.name in destinations && r.source==inferred[r.name] && r.source in existing)
    isempty(intersect(setdiff(Set(values(inferred)),reused), union(names, physical))) || throw(
        ArgumentError(
            "reserved inferred river-junction name collides with an input object",
        ),
    )
    result=River[]
    for r in s.rivers
        source=r.source
        if r.name in destinations
            (source==:auto || source==inferred[r.name] && source in reused) || throw(
                ArgumentError(
                    "river $(r.name) receives direct discharge connections and must omit source or set source=auto",
                ),
            )
            r.law==:junction || throw(
                ArgumentError(
                    "a directly connected downstream river must use the junction law",
                ),
            )
            source=inferred[r.name]
        end
        push!(result, _river_replace(r; source, target = get(inferred, r.target, r.target)))
    end
    js=vcat(
        s.river_junctions,
        [
            RiverJunction(name = inferred[n]) for
            n in sort!(collect(destinations); by = string) if !(inferred[n] in reused)
        ],
    )
    _river_replace(s; rivers = result, river_junctions = js)
end

"""Return the public direct-connection form, retaining explicit junction nodes."""
function public_river_connections(s::HydroSystem)
    inferred=Dict{Symbol,Symbol}()
    for j in s.river_junctions
        startswith(string(j.name), _INFERRED_RIVER_PREFIX) || continue
        out=findall(r->r.source==j.name, s.rivers)
        length(out)==1 || throw(ArgumentError("invalid inferred river junction $(j.name)"))
        inferred[j.name]=s.rivers[only(out)].name
    end
    rivers=[
        _river_replace(
            r;
            source = haskey(inferred, r.source) ? :auto : r.source,
            target = get(inferred, r.target, r.target),
        ) for r in s.rivers
    ]
    (
        rivers = rivers,
        river_junctions = [j for j in s.river_junctions if !haskey(inferred, j.name)],
    )
end

"""Validate a directed acyclic river graph and return a reach evaluation order.

Zero-storage junctions merge all incoming reaches into exactly one outgoing
reach. Reservoir releases remain independent decisions.
"""
function river_order(c::ScheduleCase)
    s=c.system
    storages=Set(r.name for r in s.reservoirs)
    junctions=Set(j.name for j in s.river_junctions)
    length(storages)==length(s.reservoirs) ||
        throw(ArgumentError("duplicate reservoir names"))
    length(junctions)==length(s.river_junctions) ||
        throw(ArgumentError("duplicate river junction names"))
    isempty(intersect(storages, junctions)) ||
        throw(ArgumentError("river junction and reservoir names overlap"))
    length(unique(r.name for r in s.rivers))==length(s.rivers) ||
        throw(ArgumentError("duplicate reach names"))
    sources=union(storages, junctions)
    valid=union(sources, Set(b.name for b in s.boundaries))
    incoming=Dict(n=>Int[] for n in valid)
    outgoing=Dict(n=>Int[] for n in valid)
    for (d, r) in enumerate(s.rivers)
        r.source in sources && r.target in valid || throw(
            ArgumentError("reach $(r.name) has an unknown endpoint or invalid source"),
        )
        r.source!=r.target || throw(ArgumentError("reach $(r.name) is a self-loop"))
        push!(incoming[r.target], d)
        push!(outgoing[r.source], d)
        isfinite(r.capacity) && r.capacity>0 ||
            throw(ArgumentError("invalid reach capacity"))
        r.deterministic_delay!==nothing || RiverRouting.check_curves(r.curves; capacity=r.capacity)
        RiverRouting.check_grid(r.history_grid)
        RiverRouting.check_releases(r.history_grid, r.history_release)
        last(r.history_grid)==first(c.grid) ||
            throw(ArgumentError("river history must end at the horizon start"))
        all(r.history_release .<= r.capacity) ||
            throw(ArgumentError("historical release exceeds curve domain"))
    end
    for j in junctions
        length(outgoing[j])==1 || throw(ArgumentError("junction $j needs exactly one outgoing reach"))
        d=only(outgoing[j])
        river=s.rivers[d]
        supplied=!isempty(incoming[j]) ||
            any(x.discharge_river==river.name for x in Iterators.flatten((s.plants,s.tunnels))) ||
            river.inflow>0 || any(z.object==river.name && z.attribute==:inflow && any(>(0),z.values) for z in c.operations)
        supplied || throw(
            ArgumentError(
                "junction $j needs an upstream water supply",
            ),
        )
    end
    indegree=Dict(n=>length(incoming[n]) for n in valid)
    queue=sort!([n for n in valid if indegree[n]==0]; by = string)
    order=Int[]
    visited=0
    while !isempty(queue)
        n=popfirst!(queue)
        visited+=1
        for d in outgoing[n]
            push!(order, d)
            target=s.rivers[d].target
            indegree[target]-=1
            indegree[target]==0 && push!(queue, target)
        end
    end
    visited==length(valid) || throw(ArgumentError("river graph must be acyclic"))
    order
end

"""Transfer fractions for optimization, and independent pre-horizon arrivals."""
function routing_data(c::ScheduleCase)
    RiverRouting.check_grid(c.grid)
    river_order(c)
    D=length(c.system.rivers)
    T=length(c.grid)-1
    B=Vector{Array{Float64,3}}(undef, D)
    history_arrival=zeros(D, T)
    history_initial=zeros(D)
    history_terminal=zeros(D)
    for (d, r) in enumerate(c.system.rivers)
        B[d]=transfer_coefficients(r, c.grid, c.grid)
        history_arrival[d, :]=route_volumes(r, c.grid, r.history_grid, r.history_release)
        history_initial[d]=remaining_volume(
            r,
            first(c.grid),
            r.history_grid,
            r.history_release,
        )
        history_terminal[d]=remaining_volume(
            r,
            last(c.grid),
            r.history_grid,
            r.history_release,
        )
    end
    Dict(
        "B"=>B,
        "history_arrival"=>history_arrival,
        "history_initial"=>history_initial,
        "history_terminal"=>history_terminal,
    )
end

# Interval averages of a source schedule, preserving its original cohort grid.
function source_averages(grid, original_grid, q)
    [
        sum(
            q[k]*max(
                0.0,
                min(grid[t + 1], original_grid[k + 1])-max(grid[t], original_grid[k]),
            ) for k in eachindex(q)
        )/(grid[t + 1]-grid[t]) for t in 1:(length(grid) - 1)
    ]
end

"""Route all reaches in topological order, independently of model coefficients.

`source_release` is D×T on the original case grid. Reservoir-source rows are
prescribed total reach releases, including natural inflow. Other local injections
come from `generator_q`, `tunnel_q` and natural inflow. Junction release adds the
upstream arrival volumes divided by interval duration. Optional `grid` preserves
original reservoir cohorts, but
junction mixing rebins arrivals: its timing is a convergent grid approximation,
not an exact continuous-time convolution. Transit includes historical cohorts.
"""
function route_network(c::ScheduleCase, source_release::AbstractMatrix; grid = c.grid,
    generator_q=nothing, tunnel_q=nothing, injections=nothing)
    size(source_release)==(length(c.system.rivers),length(c.prices)) ||
        throw(ArgumentError("source release matrix has incorrect dimensions"))
    source_release=injections===nothing ? river_injections(c,source_release,generator_q,tunnel_q) : injections
    order=river_order(c)
    RiverRouting.check_grid(grid)
    first(grid)==first(c.grid) && last(grid)==last(c.grid) ||
        throw(ArgumentError("routing grid must span the original horizon"))
    s=c.system
    D=length(s.rivers)
    T=length(grid)-1
    size(source_release)==(D, length(c.grid)-1) ||
        throw(ArgumentError("source release matrix has incorrect dimensions"))
    reservoirs=Set(r.name for r in s.reservoirs)
    release=zeros(D, T)
    arrival=zeros(D, T)
    transit=zeros(D, T+1)
    for d in order
        r=s.rivers[d]
        if r.source in reservoirs
            q=collect(source_release[d, :])
            cohort_grid=c.grid
            RiverRouting.check_releases(cohort_grid, q)
            release[d, :]=source_averages(grid, cohort_grid, q)
        else
            upstream=findall(e->e.target==r.source, s.rivers)
            q=source_averages(grid,c.grid,view(source_release,d,:)) +
              vec(sum(arrival[upstream, :]; dims = 1)) ./ (0.0036 .* diff(grid))
            # Floating-point conservation can put a capacity-bound flow a few
            # ulps outside its reference range. Do not hide physical excess.
            q=[x<=r.capacity+1e-10 ? min(x, r.capacity) : x for x in q]
            release[d, :]=q
            cohort_grid=grid
        end
        all(q .<= r.capacity) ||
            throw(ArgumentError("release in $(r.name) exceeds routing curve domain"))
        arrival[d, :]=route_volumes(r, grid, cohort_grid, q) +
                      route_volumes(r, grid, r.history_grid, r.history_release)
        transit[d, :]=[
            remaining_volume(r, t, cohort_grid, q) +
            remaining_volume(r, t, r.history_grid, r.history_release) for t in grid
        ]
    end
    Dict("release"=>release, "arrival_volume"=>arrival, "transit"=>transit, "order"=>order)
end

# Fixed atoms and distributed bins share the same conservation interface.
function transfer_coefficients(r, arrival_grid, release_grid)
    r.deterministic_delay===nothing &&
        return RiverRouting.coefficients(arrival_grid, release_grid, r.curves)
    B=RiverRouting.fixed_coefficients(arrival_grid, release_grid, r.deterministic_delay)
    repeat(reshape(B, size(B)..., 1), 1, 1, 2)
end
function route_volumes(r, arrival_grid, release_grid, q)
    r.deterministic_delay===nothing &&
        return RiverRouting.route(arrival_grid, release_grid, q, r.curves)
    RiverRouting.check_releases(release_grid, q)
    RiverRouting.fixed_coefficients(arrival_grid, release_grid, r.deterministic_delay)*(
        0.0036 .* diff(release_grid) .* q
    )
end
function remaining_volume(r, time, grid, q)
    r.deterministic_delay===nothing &&
        return RiverRouting.remaining_volume(time, grid, q, r.curves)
    RiverRouting.check_releases(grid, q)
    sum(
        0.0036*q[k]*(
            max(0.0, min(grid[k + 1], time)-grid[k])-max(
                0.0,
                min(grid[k + 1], time-r.deterministic_delay)-grid[k],
            )
        ) for k in eachindex(q);
        init = 0.0,
    )
end
