"""Subdivide integration intervals while retaining every control-change time."""
function refined_grid(grid; factor = 2)
    factor isa Integer && factor>=1 ||
        throw(ArgumentError("positive integer refinement required"))
    fine=Vector{Float64}(undef,factor*(length(grid)-1)+1)
    for t in 1:(length(grid)-1), j in 0:(factor-1)
        fine[factor*(t-1)+j+1]=grid[t]+(grid[t+1]-grid[t])*j/factor
    end
    fine[end]=last(grid)
    fine
end

"""Refine hydraulics and distributed river mixing, including a supplied river mesh."""
function replay_case(c, grid)
    fine=with_grid(c, grid)
    all(r->r.deterministic_delay!==nothing,c.system.rivers) && return fine
    routing_grid=sort!(unique!(vcat(refined_grid(river_routing_grid(c)),grid)))
    _river_replace(fine; routing_grid)
end

storage_convergence_limits(c,absolute,relative) = [min(absolute,isinf(relative) ? Inf :
    1e-9+relative*(r.vmax-r.vmin)) for r in c.system.reservoirs]

"""
    replay_audit(case, schedule; grid=refined_grid(schedule.grid),
                 storage_convergence=2e-3, storage_relative_convergence=1e-3,
                 level_convergence=1e-2, objective_convergence=1e-3)

Check finer physical replay and independent river transport. Operating data
retain physical times; plant ramp limits retain their original decision windows.

At common vertices, each reservoir's storage difference must be at most
`min(storage_convergence, 1e-9 + storage_relative_convergence * (vmax-vmin))`
in Mm³, and its endpoint level difference must be at most `level_convergence`
in metres. The one-litre absolute allowance keeps the comparison meaningful for
fixed or very small storage ranges. `storage_convergence` also bounds terminal
transit differences.
`objective_convergence` is relative to `max(1, abs(schedule.objective))`.
These are mesh convergence tests, not bounds on continuous-time error.
"""
replay_audit(c::ScheduleCase, x; kwargs...) = _replay_report(_replay(c, x; kwargs...))

_replay_report(check) = Dict{String,Any}(string(k)=>v for (k,v) in pairs(check))

# Numerical replay returns named quantities. Diagnostic labels never determine
# a subsequent solve's margins or its physical acceptance.
function _replay(
    c::ScheduleCase,
    x;
    grid = nothing,
    tolerance = 1e-4,
    storage_convergence = 2e-3,
    storage_relative_convergence = 1e-3,
    level_convergence = 1e-2,
    objective_convergence = 1e-3,
    transport = nothing,
)
    all(x -> x isa Real && !isnan(x) && x >= 0,
        (storage_convergence, storage_relative_convergence, level_convergence,
         objective_convergence)) || throw(ArgumentError("nonnegative replay tolerances required"))
    began=time()
    errors=String[]
    violations=Dict{String,Float64}()
    shortfalls=(arrival=0.0,power=0.0,discharge=0.0)
    original=validate(c,x;transport)
    if !original["valid"]
        push!(errors,"original equation audit failed")
        return (;valid=false,errors,violations,original_audit=original,shortfalls,seconds=time()-began)
    end
    x=_schedule_input(x;case=c)
    grid=grid===nothing ? refined_grid(x.grid) : grid
    nd=all(r->r.deterministic_delay!==nothing,c.system.rivers) ?
        _transport_data(hydraulic_case(c,x.grid),transport) : nothing
    nd===nothing && transport!==nothing && throw(ArgumentError(
        "compiled deterministic transport cannot serve distributed acceptance"))
    length(grid)>length(x.grid) ||
        push!(errors,"acceptance requires a strictly finer replay grid")
    if !isempty(errors)
        return (; valid=false, errors, violations, original_audit=original,
            shortfalls, seconds=time()-began)
    end
    try
        mesh=TimeDiscretization(c.grid,x.grid)
        firstcell=first.(mesh.cells)
        original_u=x.commitment[:,firstcell]
        fine_mesh=TimeDiscretization(x.grid,grid)
        cc=replay_case(hydraulic_case(c,x.grid),grid)
        replay_input=_river_replace(c;routing_grid=cc.routing_grid)
        opening=_physical_tunnel_openings(c,x.tunnel_opening[:,firstcell])
        z=dispatch_from_controls(replay_input,original_u,x.generator_q[:,firstcell],
            x.gate[:,firstcell];grid,tunnel_opening=opening,
            reserves=_reserve_slice(x.reserves,firstcell),head_reference=interval_means(mesh,x.head))
        audit=z["validation"]
        merge!(violations, audit["residuals"])
        append!(errors, audit["errors"])
        fine_shortfalls=audit["shortfalls"]
        for (i,g) in enumerate(c.system.generators)
            violations["unit_power_$(g.name)"]=fine_shortfalls.unit_power[i]
        end
        transport_check=transport_audit(
            hydraulic_case(c,x.grid),
            x;
            schedule_tolerance = storage_convergence,
            absolute_tolerance = isfinite(storage_convergence) ? storage_convergence/10 :
                                 2e-4,
            transport = nd,
            tolerance,
        )
        for (name,residual) in transport_check["violations"]
            violations[name]=max(get(violations,name,0.0),residual)
        end
        arrival=max.(fine_shortfalls.arrival,transport_check["arrival_shortfall"])
        for (i,r) in enumerate(c.system.rivers)
            violations["environmental_$(r.name)"]=arrival[i]
        end
        shortfalls=(
            arrival=maximum((arrival[i] for (i,r) in enumerate(c.system.rivers)
                if r.arrival_policy==:pointwise);init=0.0),
            power=max(maximum(fine_shortfalls.unit_power;init=0.0),
                maximum(fine_shortfalls.plant_power;init=0.0)),
            discharge=max(maximum(fine_shortfalls.unit_discharge;init=0.0),
                maximum(fine_shortfalls.plant_discharge;init=0.0)),
        )
        transport_check["accepted"] ||
            push!(errors, "transport convergence, operating limits or schedule discrepancy fail acceptance")
        physical=isempty(errors)
        edges=fine_mesh.vertices
        storage_errors=[maximum(abs, z.storage[i, edges]-x.storage[i, :]; init=0.0)
            for i in eachindex(c.system.reservoirs)]
        storage_limits=storage_convergence_limits(c,storage_convergence,storage_relative_convergence)
        level_errors=[maximum(abs(head(r,z.storage[i,edges[t]])-head(r,x.storage[i,t]))
            for t in eachindex(x.grid); init=0.0) for (i,r) in enumerate(c.system.reservoirs)]
        storage_error=maximum(storage_errors; init=0.0)
        level_error=maximum(level_errors; init=0.0)
        transit_error=maximum(abs, z.terminal_transit-x.terminal_transit; init = 0.0)
        replayed_objective=z.objective
        objective_error=replayed_objective-x.objective
        numerical=all(storage_errors .<= storage_limits) && level_error<=level_convergence &&
                  transit_error<=storage_convergence &&
                  abs(objective_error)<=objective_convergence*max(1.0, abs(x.objective))
        numerical || push!(
            errors,
            "finer replay exceeded storage, level, transit or objective convergence tolerance",
        )
        (;
            valid=isempty(errors),
            physically_valid=physical,
            numerically_converged=numerical,
            errors, violations, shortfalls,
            original_audit=original,
            transport_audit=transport_check,
            replay_grid=collect(grid),
            replay_intervals=length(grid)-1,
            max_water_residual_Mm3=get(violations, "global_water", Inf),
            max_storage_difference_Mm3=storage_error,
            storage_differences_Mm3=storage_errors,
            storage_convergence_tolerances_Mm3=storage_limits,
            max_level_difference_m=level_error,
            level_differences_m=level_errors,
            terminal_transit_difference_Mm3=transit_error,
            replayed_objective,
            objective_difference=objective_error,
            storage_convergence_tolerance_Mm3=storage_convergence,
            relative_storage_convergence_tolerance=storage_relative_convergence,
            level_convergence_tolerance_m=level_convergence,
            relative_objective_convergence_tolerance=objective_convergence,
            seconds=time()-began,
            certificate="numerical finer-grid and transport acceptance only",
        )
    catch e
        e isa Union{ArgumentError,DomainError} || rethrow()
        push!(errors, sprint(showerror, e))
        (; valid=false, errors, violations, original_audit=original,
            shortfalls, seconds=time()-began)
    end
end

"""Solve a fixed commitment and check an independently integrated finer trajectory.

After a failed replay, refine distributed transport entry cohorts for the current
controls. If transport already meets its volume tolerance, add finer hydraulic equations,
sharing the original unit discharges, outlet openings and reserve allocations.
Fixed-speed pump discharges remain hydraulic states. Operating schedules and ramps
retain their original windows. No measured residual becomes a blanket margin.
The returned case records the transport mesh. Operating decisions retain their
requested intervals; a supplied resolved initial schedule keeps its physical
grid. The latest physical trajectory starts the next solve.
`refine_grid` is retained for source compatibility; both values now retain controls.
`operational_margin` is an explicit optional power inset, never an automatic retry rule.
"""
solve_verified(c::ScheduleCase; kwargs...) = _verified_report(_verified_dispatch(c; kwargs...))

function _verified_report(result)
    Dict{String,Any}(
        "accepted"=>result.accepted,
        "problem"=>result.case,
        "case"=>result.case,
        "solution"=>result.solution,
        "audit"=>_replay_report(result.audit),
        "attempts"=>result.attempts,
        "total_seconds"=>result.total_seconds,
        "global_optimality_proven"=>false,
    )
end

function _verified_dispatch(
    c::ScheduleCase;
    u = ones(Int, length(c.system.generators), length(c.prices)),
    initial = nothing,
    max_refinements = 2,
    solver = "Ipopt",
    seed = 1,
    replay_factor = 2,
    time_limit = 45.0,
    feasibility_only = false,
    warm = nothing,
    tunnel_opening = nothing,
    river_on = nothing,
    operational_margin = 0.0,
    deadline = Inf,
    refine_grid = true,
    storage_convergence = 2e-3,
    storage_relative_convergence = 1e-3,
    level_convergence = 1e-2,
    objective_convergence = 1e-3,
)
    isfinite(time_limit) && time_limit>0 ||
        throw(ArgumentError("positive finite per-solve time limit required"))
    isnan(deadline) && throw(ArgumentError("invalid dispatch deadline"))
    replay_factor isa Integer && replay_factor>=2 ||
        throw(ArgumentError("acceptance requires replay_factor >= 2"))
    max_refinements isa Integer && max_refinements>=0 ||
        throw(ArgumentError("nonnegative integer refinement budget required"))
    refine_grid isa Bool || throw(ArgumentError("refine_grid must be boolean"))
    all(x->x isa Real && !isnan(x) && x>=0,
        (storage_convergence,storage_relative_convergence,level_convergence,objective_convergence)) ||
        throw(ArgumentError("nonnegative replay tolerances required"))
    initial_schedule=initial===nothing ? nothing : _schedule_input(initial;case=c)
    initial_windows=initial===nothing ? nothing : operating_window_values(c,initial_schedule)
    initial!==nothing && initial_windows.commitment!=u &&
        throw(ArgumentError("initial schedule commitment differs from requested commitment"))
    if river_on!==nothing
        river_on=fixed_river_on(c,river_on,nothing)
        if initial!==nothing
            initial_states=river_on_schedule(c,initial_windows.river_release)
            all(r.min_allowed_flow===nothing || initial_states[i,:]==river_on[i,:]
                for (i,r) in enumerate(c.system.rivers)) ||
                throw(ArgumentError("initial candidate river states differ from requested states"))
        end
    end
    isfinite(operational_margin) && operational_margin >= 0 ||
        throw(ArgumentError("nonnegative finite operational_margin required"))
    openings=tunnel_opening===nothing ? nothing :
        _physical_tunnel_openings(c,tunnel_opening)
    if initial!==nothing && openings!==nothing
        _physical_tunnel_openings(c,initial_schedule.tunnel_opening[:,first.(TimeDiscretization(c.grid,initial_schedule.grid).cells)])==openings ||
            throw(ArgumentError("initial candidate tunnel openings differ from requested openings"))
    end
    began=time()
    attempts=Any[]
    integration_grid=nothing
    solution_case=c
    storage_limits=storage_convergence_limits(c,storage_convergence,storage_relative_convergence)
    x=initial
    audit=(valid=false,errors=["dispatch deadline exhausted"],shortfalls=(arrival=0.0,power=0.0,discharge=0.0))
    for iteration in 0:max_refinements
        remaining=min(time_limit,deadline-time())
        iteration==0 && initial!==nothing || remaining>0 || break
        x=iteration==0 && initial!==nothing ? initial :
          solve_case(
            c;
            u,
            solver = solver,
            seed = seed,
            warm = warm,
            operational_margin,
            replay_grid=integration_grid,
            time_limit = remaining,
            feasibility_only = feasibility_only,
            tunnel_opening = openings,
            river_on,
        )
        dispatch=x
        solution_case=c
        x=dispatch isa SolveResult ? dispatch.solution : dispatch
        audit=x!==nothing && haskey(x,"objective") ?
              _replay(c, x; grid = integration_grid===nothing ? refined_grid(get(x,"grid",c.grid);factor=replay_factor) : integration_grid,
                  storage_convergence, storage_relative_convergence, level_convergence, objective_convergence) :
              (valid=false,errors=["NLP candidate has no objective"],shortfalls=(arrival=0.0,power=0.0,discharge=0.0))
        original=get(audit, :original_audit, Dict{String,Any}())
        physical=get(original, "valid", false)
        push!(
            attempts,
            Dict(
                "iteration"=>iteration,
                "intervals"=>length(c.prices),
                "integration_intervals"=>integration_grid===nothing ? length(c.prices) : length(integration_grid)-1,
                "transport_intervals"=>length(river_routing_grid(c))-1,
                "operational_margin_MW"=>operational_margin,
                "status"=>dispatch["status"],
                "objective"=>(x===nothing ? nothing : x["objective"]),
                "model_valid"=>physical,
                "model_errors"=>get(original, "errors", String[]),
                "replay"=>_replay_report(audit),
                "solve_seconds"=>get(dispatch, "total_seconds", 0.0),
                "construction_seconds"=>get(dispatch, "construction_seconds", 0.0),
                "optimization_seconds"=>get(dispatch, "seconds", 0.0),
                "validation_seconds"=>get(dispatch, "validation_seconds", 0.0),
                "forward_reconstruction_seconds"=>get(dispatch, "forward_reconstruction_seconds", 0.0),
                "forward_reconstructed"=>get(dispatch, "forward_reconstructed", false),
                "iterations"=>get(dispatch, "iterations", nothing),
                "warm_start"=>get(dispatch, "warm_start", nothing),
            ),
        )
        x!==nothing && physical && (x=_schedule_input(x;case=c))
        audit.valid && break
        iteration==max_refinements && break
        x===nothing && break
        # Reserve the preceding build's measured cost before starting a retry.
        # Construction can still overrun; deadlines are checked between phases.
        time()+get(dispatch,"construction_seconds",0.0)<deadline || break
        try
            start,_=_reconstruct_candidate(c,x)
            grid=refine_transport_grid(c,start,storage_limits,storage_convergence;deadline)
            if grid!=river_routing_grid(c)
                next=_river_replace(c;routing_grid=grid)
                warm,_=_reconstruct_candidate(next,start)
                c=next
            else
                integration_grid=integration_grid===nothing ? refined_grid(c.grid;factor=replay_factor) :
                    refined_grid(integration_grid;factor=replay_factor)
                warm=start
            end
        catch error
            error isa Union{ArgumentError,DomainError} || rethrow()
            break
        end
    end
    (;accepted=audit.valid,case=solution_case,solution=x,audit,attempts,total_seconds=time()-began)
end
