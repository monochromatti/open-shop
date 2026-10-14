# Distance to the allowed adjacent support; tiny nonadjacent weights remain
# visible to the same numerical tolerance used for the other constraints.
function _sos2_residual(x,weights)
    length(x)<=2 && return 0.0
    ordered=x[sortperm(weights)]
    minimum(maximum(abs(ordered[j]) for j in eachindex(ordered) if j!=i && j!=i+1;init=0.0)
        for i in 1:(length(ordered)-1))
end

"""Lift a physical incumbent into all algebraic variables and audit feasibility."""
function _lift_start!(b, c, warm; tolerance = 1e-7)
    assigned=Dict{VariableRef,Float64}()
    variables=all_variables(b.m)
    for v in variables
        set_start_value(v, nothing)
    end
    function put(variable, x)
        variable===nothing && return
        y=Float64(x)
        isfinite(y) || error("nonfinite lift for $(name(variable))")
        assigned[variable]=y
        set_start_value(variable, y)
    end
    function table(graph, curve::TableCurve, x)
        table(graph,z->table_value(curve,z;extrapolation=:linear),x)
    end
    table(graph, ::Nothing, x) = nothing
    function table(graph, evaluate, x)
        graph===nothing && return
        weights=_global_tensor_coordinate_weights(graph.coordinate.nodes,x)
        for (variable,value) in zip(graph.coordinate.weights,weights)
            variable isa VariableRef && put(variable,value)
        end
        # Evaluate the physical curve independently of the interpolation rows.
        put(graph.y,evaluate(x))
    end
    try
        validate(c, warm)["valid"] || error("seed fails original physical audit")
        physical=_schedule_input(warm;case=c)
        s=c.system
        T=length(c.prices)
        ix=nodeindex(s)
        openings=physical.tunnel_opening
        for (variable,value) in _tunnel_gate_values(b.gates,c,openings)
            put(variable,value)
        end
        for variable in variables
            is_fixed(variable) && put(variable,fix_value(variable))
        end
        for (i,r) in enumerate(s.reservoirs), t in 1:(T+1)
            put(b.v[i,t],physical.storage[i,t]/r.vmax)
            table(b.graphs.reservoir_vertices[i,t],r.level_curve,physical.storage[i,t])
        end
        for i in axes(physical.head,1), t in 1:T
            put(b.h[i,t],physical.head[i,t]/250)
        end
        for (i,e) in enumerate(s.tunnels), t in 1:T
            q=physical.tunnel_q[i,t]
            put(b.q[i,t],q/50)
            graph=b.graphs.tunnel_losses[i,t]
            graph===nothing && continue
            put(graph.positive,max(q,0.0))
            put(graph.negative,max(-q,0.0))
            put(graph.signed_square,q*abs(q)/2500)
            put(graph.open_head,openings[i,t]>0 ?
                physical.head[ix[e.source],t]-physical.head[ix[e.target],t] : 0.0)
            # Fixed direction branches remain fixed even at zero flow.
            direction=graph.direction
            choice=direction!==nothing && is_fixed(direction) ? fix_value(direction) :
                (openings[i,t]>0 && q>=0 ? 1.0 : 0.0)
            put(direction,choice)
        end
        for (i,g) in enumerate(s.generators), t in 1:T
            q,p,u=physical.generator_q[i,t],physical.power[i,t],physical.commitment[i,t]
            put(b.gq[i,t],q/50)
            put(b.p[i,t],p/40)
            if b.transitions!==nothing
                put(b.transitions.u[i,t],u)
                prior=t==1 ? g.initial_on : physical.commitment[i,t-1]
                put(b.transitions.su[i,t],max(0,u-prior))
                put(b.transitions.sd[i,t],max(0,prior-u))
            end
            plant=plantof(s,g)
            total=sum(physical.generator_q[j,t] for (j,z) in enumerate(s.generators)
                if z.plant==plant.name)
            graph=b.graphs.generator_graphs[i,t]
            receiver=physical.head[ix[plant.target],t]
            h=physical.head[ix[plant.source],t]-outlet_head(plant,receiver)-tailwater(plant,total)
            put(graph.head.head,h)
            if graph.head.outlet!==nothing
                floor=plant.outlet_head_floor
                nodes=graph.head.outlet.coordinate.nodes
                table(graph.head.outlet,TableCurve(nodes,max.(nodes,floor)),receiver)
            end
            table(graph.head.tailwater,plant.tailwater_curve,total)
            if graph.turbine===nothing
                put(graph.efficiency,efficiency(g,q,h))
            else
                for (variable,value) in _global_tensor_turbine_values(graph.turbine,q,h)
                    put(variable,value)
                end
                put(graph.efficiency,turbine_efficiency(g.turbine_table,q,h;extrapolation=:linear))
                table(graph.flow_min,z->turbine_qmin(g.turbine_table,z;extrapolation=:linear),h)
                table(graph.flow_max,z->turbine_qmax(g.turbine_table,z;extrapolation=:linear),h)
            end
            table(graph.electrical,g.generator_efficiency_curve,p)
        end
        if b.transitions!==nothing
            states=b.transitions
            for (i,plant) in enumerate(s.plants), t in 1:T
                states.plant_u[i,t]===nothing && continue
                members=states.bounds.members[i]
                u=any(physical.commitment[j,t]==1 for j in members) ? 1.0 : 0.0
                prior=t==1 ? states.bounds.plant_initial_on[i] :
                    (any(physical.commitment[j,t-1]==1 for j in members) ? 1.0 : 0.0)
                put(states.plant_u[i,t],u)
                put(states.plant_su[i,t],max(0,u-prior))
                put(states.plant_sd[i,t],max(0,prior-u))
            end
        end
        for (i,r) in enumerate(s.reservoirs), t in 1:T
            table(b.graphs.reservoir_levels[i,t],r.level_curve,
                (physical.storage[i,t]+physical.storage[i,t+1])/2)
        end
        shortfall=release_requirements(c,physical.river_release).shortfall
        for (i,r) in enumerate(s.rivers), t in 1:T
            put(b.rq[i,t],physical.river_release[i,t]/100)
            put(b.a[i,t],physical.gate[i,t])
            put(b.shortfall_release[i,t],shortfall[i,t])
            graph=b.graphs.river_graphs[i,t]
            if is_channel(r)
                upstream,downstream=physical.head[ix[r.source],t],physical.head[ix[r.target],t]
                delta=upstream-downstream
                put(graph.delta,delta)
                for (variable,value) in _global_level_flow_values(graph.channel,delta,
                    channel_reference(r,upstream,downstream))
                    put(variable,value)
                end
            elseif r.law!=:junction
                h=physical.head[ix[r.source],t]
                table(graph.discharge,r.discharge_curve,h)
                wet=max(h-r.crest,0.0)
                put(graph.wet,wet)
                put(graph.root,sqrt(wet))
                put(graph.branch,h>=r.crest ? 1.0 : 0.0)
            end
        end
        _lift_distributed_references!(b.m, assigned)
        _lift_power_hulls!(b.m, assigned)
        _lift_table_power_bounds!(b.m, assigned)
        for (variable, val) in assigned
            set_start_value(variable, val)
        end
        missing=[name(v) for v in variables if !haskey(assigned, v)]
        isempty(missing) ||
            return Dict("valid"=>false, "missing"=>missing, "assigned"=>length(assigned))
        residual=0.0
        worst=""
        count=0
        violation(x, set) =
            set isa MOI.EqualTo ? abs(x-set.value) :
            set isa MOI.LessThan ? max(0.0, x-set.upper) :
            set isa MOI.GreaterThan ? max(0.0, set.lower-x) :
            set isa MOI.Interval ? max(0.0, set.lower-x, x-set.upper) :
            set isa MOI.ZeroOne ? max(abs(x-round(x)), max(0.0, -x, x-1)) :
            set isa MOI.Integer ? abs(x-round(x)) :
            error("unsupported audit set $(typeof(set))")
        for ref in all_constraints(b.m; include_variable_in_set_constraints = true)
            obj=constraint_object(ref)
            e=if obj.set isa MOI.SOS2
                x=[JuMP.value(v->assigned[v],z) for z in obj.func]
                _sos2_residual(x,obj.set.weights)
            else
                x=obj.func isa Number ? obj.func : JuMP.value(v->assigned[v], obj.func)
                violation(x,obj.set)
            end
            isfinite(e) || error("nonfinite constraint residual")
            count+=1
            if e>residual
                residual=e
                worst=string(ref)
            end
        end
        Dict(
            "valid"=>residual<=tolerance,
            "assigned"=>length(assigned),
            "variables"=>length(variables),
            "constraints"=>count,
            "max_residual"=>residual,
            "worst_constraint"=>worst,
            "tolerance"=>tolerance,
            "objective"=>JuMP.value(v->assigned[v], b.obj),
        )
    catch error
        Dict(
            "valid"=>false,
            "assigned"=>length(assigned),
            "error"=>sprint(showerror, error),
        )
    end
end
