"""Integer connectivity in model order, shared by numerical and JuMP equations.

Incidence is positive into a node. A plant or tunnel outfall contributes water
at its river destination; its hydraulic target still supplies the outlet head.
River source/target index zero denotes a river confluence, not a hydraulic node.
This contains connectivity only, so changing operating limits needs no new graph.
"""
struct HydraulicTopology
    node_names::Vector{Symbol}
    node_index::Dict{Symbol, Int}
    tunnel_nodes::Vector{Tuple{Int, Int}}
    plant_nodes::Vector{Tuple{Int, Int}}
    generator_plant::Vector{Int}
    plant_generators::Vector{Vector{Int}}
    tunnel_incidence::Vector{Vector{Tuple{Int, Int}}}
    generator_incidence::Vector{Vector{Tuple{Int, Int}}}
    river_incidence::Vector{Vector{Tuple{Int, Bool, Bool}}}
    river_source::Vector{Int}
    river_target::Vector{Int}
    river_upstream::Vector{Vector{Int}}
end

function HydraulicTopology(s::HydroSystem)
    node_names = nodes(s)
    node_index = Dict(name => i for (i, name) in enumerate(node_names))
    plant_index = Dict(p.name => i for (i, p) in enumerate(s.plants))
    tunnel_nodes = [(node_index[e.source], node_index[e.target]) for e in s.tunnels]
    plant_nodes = [(node_index[p.source], node_index[p.target]) for p in s.plants]
    generator_plant = [plant_index[g.plant] for g in s.generators]
    plant_generators = [Int[] for _ in s.plants]
    tunnel_incidence = [Tuple{Int, Int}[] for _ in node_names]
    generator_incidence = [Tuple{Int, Int}[] for _ in node_names]
    river_incidence = [Tuple{Int, Bool, Bool}[] for _ in node_names]
    for (j, e) in enumerate(s.tunnels)
        source, target = tunnel_nodes[j]
        push!(tunnel_incidence[source], (j, water_incidence(e, e.source)))
        source == target || push!(tunnel_incidence[target], (j, water_incidence(e, e.target)))
    end
    for (j, p) in enumerate(generator_plant)
        push!(plant_generators[p], j)
        plant = s.plants[p]
        source, target = plant_nodes[p]
        push!(generator_incidence[source], (j, water_incidence(plant, plant.source)))
        source == target ||
            push!(generator_incidence[target], (j, water_incidence(plant, plant.target)))
    end
    river_source = [get(node_index, r.source, 0) for r in s.rivers]
    river_target = [get(node_index, r.target, 0) for r in s.rivers]
    incoming = Dict{Symbol, Vector{Int}}()
    for (j, r) in enumerate(s.rivers)
        push!(get!(Vector{Int}, incoming, r.target), j)
        source, target = river_source[j], river_target[j]
        source == 0 || push!(river_incidence[source], (j, source == target, true))
        target == 0 || source == target || push!(river_incidence[target], (j, true, false))
    end
    river_upstream = [get(incoming, r.source, Int[]) for r in s.rivers]
    HydraulicTopology(
        node_names,
        node_index,
        tunnel_nodes,
        plant_nodes,
        generator_plant,
        plant_generators,
        tunnel_incidence,
        generator_incidence,
        river_incidence,
        river_source,
        river_target,
        river_upstream,
    )
end

"""Storage evolution and zero-storage junction conservation, in physical units.

`V` is Mm³, `Q`, `GQ`, and `RQ` are m³/s, and `arrivals` is interval volume in
Mm³. Both optimization formulations use these same linear balance equations.
"""
function constrain_water_balance!(model, c, topology, V, Q, GQ, RQ, arrivals)
    s = c.system
    dt = diff(c.grid)
    nr = length(s.reservoirs)
    for i in 1:(nr + length(s.junctions)), t in eachindex(dt)
        net = sum(sign * Q[j, t] for (j, sign) in topology.tunnel_incidence[i]; init = 0.0)
        net += sum(sign * GQ[j, t] for (j, sign) in topology.generator_incidence[i]; init = 0.0)
        net += sum(
            (incoming ? arrivals[j, t] / (0.0036 * dt[t]) : 0.0) - (
                outgoing ?
                RQ[j, t] - opinterval(c, s.rivers[j].name, :inflow, t, s.rivers[j].inflow) :
                0.0
            ) for (j, incoming, outgoing) in topology.river_incidence[i];
            init = 0.0,
        )
        if i <= nr
            reservoir = s.reservoirs[i]
            inflow = opinterval(c, reservoir.name, :inflow, t, reservoir.inflow)
            @constraint(
                model,
                (V[i, t + 1] - V[i, t] - 0.0036 * dt[t] * (net + inflow)) / reservoir.vmax == 0
            )
        else
            @constraint(model, net / 50 == 0)
        end
    end
    nothing
end

"""Bound cumulative storage at river transport vertices inside a control interval.

The hydraulic quadrature holds tunnel, turbine and outlet rates over the interval.
Routed arrivals need not be constant. Their integrated volume determines storage
at each additional transport vertex. These rows guide dispatch; independent
hydraulic and transport refinement still determine numerical acceptance.
"""
function constrain_storage_paths!(model, c, topology, V, Q, GQ, RQ, routing)
    for (i, reservoir) in enumerate(c.system.reservoirs)
        incoming = [j for (j, enters, _) in topology.river_incidence[i] if enters]
        isempty(incoming) && continue
        times = routing.exact ? sort!(unique!(vcat(
            c.grid, [routing_arrival_knots(c, routing, j) for j in incoming]...,
        ))) : routing.rd.grid
        for time in times
            time in c.grid && continue
            t = searchsortedlast(c.grid, time)
            start = c.grid[t]
            duration = time - start
            net = sum(sign * Q[j, t] for (j, sign) in topology.tunnel_incidence[i]; init = 0.0)
            net += sum(sign * GQ[j, t] for (j, sign) in topology.generator_incidence[i]; init = 0.0)
            net -= sum(RQ[j, t] for (j, _, leaves) in topology.river_incidence[i] if leaves; init = 0.0)
            natural = opaverage(c, reservoir.name, :inflow, start, time, reservoir.inflow)
            natural += sum(opaverage(c, c.system.rivers[j].name, :inflow,
                start, time, c.system.rivers[j].inflow)
                for (j, _, leaves) in topology.river_incidence[i] if leaves; init = 0.0)
            arrivals = sum(routing_volume_expression(model, c, routing, j, start, time)
                for j in incoming; init = 0.0)
            volume = V[i, t] + 0.0036 * duration * (net + natural) + arrivals
            lower, upper = soft_storage_bounds(c, reservoir, time)
            @constraint(model, volume / reservoir.vmax >= lower / reservoir.vmax)
            @constraint(model, volume / reservoir.vmax <= upper / reservoir.vmax)
        end
    end
    nothing
end

"""Worst storage-bound excess and its time for each reservoir.

Between arrival and operating-data knots the net rate is affine, so storage is
quadratic. Check its endpoints and its stationary point without time stepping.
Each control boundary restarts from the supplied storage; ordinary balance rows
separately check agreement with the next supplied vertex. This checks the routed
cohort representation, not an exact continuous river-network solution.
"""
function storage_path_residuals(c,x,routed)
    s=c.system
    topology=HydraulicTopology(s)
    excess=zeros(length(s.reservoirs))
    times=fill(first(c.grid),length(s.reservoirs))
    for (i,r) in enumerate(s.reservoirs)
        incoming=[d for (d,enters,_) in topology.river_incidence[i] if enters]
        outgoing=[d for (d,_,leaves) in topology.river_incidence[i] if leaves]
        objects=Set([r.name;[s.rivers[d].name for d in outgoing]])
        changes=[t for z in c.operations if z.object in objects for t in z.times
            if first(c.grid)<t<last(c.grid)]
        knots=sort!(unique!(vcat(c.grid,changes,
            [arrival_knots(s.rivers[d],routed.cohorts[d].grid) for d in incoming]...)))
        left=zeros(length(knots));right=similar(left);fill!(right,0.0)
        for d in incoming
            river,cohort=s.rivers[d],routed.cohorts[d]
            for (side,rates) in ((:left,left),(:right,right))
                rates .+= cohort_arrival_rates(river,cohort.grid,cohort.release,knots;side)
                rates .+= cohort_arrival_rates(river,river.history_grid,river.history_release,knots;side)
            end
        end
        volume=x.storage[i,1]
        violation(v,t)=begin
            lower,upper=storage_bounds(c,r,t)
            max(lower-v,v-upper,0.0)
        end
        update_maximum!(value,time)=begin
            residual=violation(value,time)
            if residual>excess[i] || isnan(residual)
                excess[i]=residual
                times[i]=time
            end
        end
        for k in 1:length(knots)-1
            a,b=knots[k],knots[k+1]
            t=searchsortedlast(c.grid,a)
            a==c.grid[t] && (volume=x.storage[i,t])
            net=sum(sign*x.tunnel_q[j,t] for (j,sign) in topology.tunnel_incidence[i];init=0.0)
            net+=sum(sign*x.generator_q[j,t] for (j,sign) in topology.generator_incidence[i];init=0.0)
            net-=sum(x.river_release[d,t] for d in outgoing;init=0.0)
            net+=opvalue(c,r.name,:inflow,a,r.inflow)
            net+=sum(opvalue(c,s.rivers[d].name,:inflow,a,s.rivers[d].inflow)
                for d in outgoing;init=0.0)
            before,after=net+right[k],net+left[k+1]
            duration=b-a
            finish=volume+0.0036duration*(before+after)/2
            update_maximum!(volume,a)
            update_maximum!(finish,b)
            change=after-before
            if change!=0
                fraction=-before/change
                if 0<fraction<1
                    stationary=volume+0.0036duration*fraction*(before+fraction*change/2)
                    update_maximum!(stationary,a+duration*fraction)
                end
            end
            volume=finish
        end
    end
    (;excess,times)
end
