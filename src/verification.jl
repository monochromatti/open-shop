"""Subdivide integration intervals while retaining every control-change time."""
function refined_grid(grid; factor = 2, existing = grid)
    factor isa Integer && factor>=1 ||
        throw(ArgumentError("positive integer refinement required"))
    RiverRouting.check_grid(grid)
    fine=Float64[]
    sizehint!(fine,factor*(length(grid)-1)+1)
    for t in 1:(length(grid)-1)
        a,b=grid[t],grid[t+1]
        push!(fine,a)
        for j in 1:(factor-1)
            vertex=refinement_vertex(a,b,j/factor,existing)
            a<vertex<b && push!(fine,vertex)
        end
    end
    push!(fine,last(grid))
    sort!(unique!(fine))
end

"""Refine hydraulics and distributed river mixing, including a supplied river mesh."""
function replay_case(c, grid)
    fine=with_grid(c, grid)
    all(r->r.deterministic_delay!==nothing,c.system.rivers) && return fine
    existing=sort!(unique!(vcat(river_routing_grid(c),grid)))
    routing_grid=sort!(unique!(vcat(refined_grid(river_routing_grid(c);existing),grid)))
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

"""Numerical replay quantities. Unevaluated checks remain unsuccessful.

Storage and transit differences are Mm³, level differences are metres and
objective differences are currency units. Constraint times are hours. Reports
are retained for presentation; mesh selection reads the numerical fields.
"""
Base.@kwdef struct ReplayResult
    valid::Bool = false
    original_valid::Bool = false
    physically_valid::Bool = false
    numerically_converged::Bool = false
    errors::Vector{String} = String[]
    violations::Dict{String,Float64} = Dict{String,Float64}()
    shortfalls::@NamedTuple{arrival::Float64,power::Float64,discharge::Float64} =
        (arrival=0.0,power=0.0,discharge=0.0)
    original_audit::Dict{String,Any} = Dict{String,Any}()
    transport_audit::Union{Nothing,Dict{String,Any}} = nothing
    replay_grid::Vector{Float64} = Float64[]
    constraint_times::Vector{Float64} = Float64[]
    storage_differences_Mm3::Vector{Float64} = Float64[]
    storage_convergence_tolerances_Mm3::Vector{Float64} = Float64[]
    level_differences_m::Vector{Float64} = Float64[]
    terminal_transit_difference_Mm3::Float64 = Inf
    replayed_objective::Union{Nothing,Float64} = nothing
    objective_difference::Float64 = Inf
    storage_convergence_tolerance_Mm3::Float64 = 2e-3
    relative_storage_convergence_tolerance::Float64 = 1e-3
    level_convergence_tolerance_m::Float64 = 1e-2
    relative_objective_convergence_tolerance::Float64 = 1e-3
    estimated_error_factor::Float64 = 1.0
    seconds::Float64 = 0.0
end

function _replay_report(check::ReplayResult)
    report=Dict{String,Any}(string(k)=>getfield(check,k) for k in fieldnames(ReplayResult))
    # A failed probe has no finer trajectory or numerical differences to report.
    if !isempty(check.replay_grid)
        report["replay_intervals"]=length(check.replay_grid)-1
        report["max_water_residual_Mm3"]=get(check.violations,"global_water",Inf)
        report["max_storage_difference_Mm3"]=maximum(check.storage_differences_Mm3;init=0.0)
        report["max_level_difference_m"]=maximum(check.level_differences_m;init=0.0)
        report["certificate"]="numerical finer-grid and transport acceptance only"
    else
        delete!(report,"replay_grid")
    end
    report
end

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
    error_factor = 1.0,
    require_original = true,
)
    all(x -> x isa Real && !isnan(x) && x >= 0,
        (storage_convergence, storage_relative_convergence, level_convergence,
         objective_convergence)) || throw(ArgumentError("nonnegative replay tolerances required"))
    isfinite(error_factor) && error_factor>=1 ||
        throw(ArgumentError("finite error factor >= 1 required"))
    began=time()
    errors=String[]
    violations=Dict{String,Float64}()
    shortfalls=(arrival=0.0,power=0.0,discharge=0.0)
    limits=(;storage_convergence_tolerance_Mm3=storage_convergence,
        relative_storage_convergence_tolerance=storage_relative_convergence,
        level_convergence_tolerance_m=level_convergence,
        relative_objective_convergence_tolerance=objective_convergence,
        estimated_error_factor=error_factor)
    original=validate(c,x;transport)
    if !original["valid"]
        push!(errors,"original equation audit failed")
        require_original && return ReplayResult(;errors,violations,original_audit=original,
            shortfalls,seconds=time()-began,limits...)
    end
    x=_schedule_input(x;case=c)
    grid=grid===nothing ? refined_grid(x.grid;existing=river_routing_grid(c;physical_grid=x.grid)) : grid
    nd=all(r->r.deterministic_delay!==nothing,c.system.rivers) ?
        _transport_data(hydraulic_case(c,x.grid),transport) : nothing
    nd===nothing && transport!==nothing && throw(ArgumentError(
        "compiled deterministic transport cannot serve distributed acceptance"))
    if length(grid)<=length(x.grid)
        push!(errors,"acceptance requires a strictly finer replay grid")
        return ReplayResult(;errors,violations,original_audit=original,original_valid=original["valid"],
            shortfalls,seconds=time()-began,limits...)
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
        physical=audit["valid"] && transport_check["operating_valid"]
        edges=fine_mesh.vertices
        storage_errors=[maximum(abs, z.storage[i, edges]-x.storage[i, :]; init=0.0)
            for i in eachindex(c.system.reservoirs)]
        storage_limits=storage_convergence_limits(c,storage_convergence,storage_relative_convergence)
        level_errors=[maximum(abs(head(r,z.storage[i,edges[t]])-head(r,x.storage[i,t]))
            for t in eachindex(x.grid); init=0.0) for (i,r) in enumerate(c.system.reservoirs)]
        level_error=maximum(level_errors; init=0.0)
        transit_error=maximum(abs, z.terminal_transit-x.terminal_transit; init = 0.0)
        replayed_objective=z.objective
        objective_error=replayed_objective-x.objective
        numerical=all(error_factor .*storage_errors .<= storage_limits) &&
                  error_factor*level_error<=level_convergence &&
                  error_factor*transit_error<=storage_convergence &&
                  error_factor*abs(objective_error)<=objective_convergence*max(1.0,abs(x.objective))
        transport_error=max(transport_check["cumulative_schedule_discrepancy_Mm3"],
            transport_check["terminal_schedule_discrepancy_Mm3"])+
            max(transport_check["estimated_cumulative_error_Mm3"],
                transport_check["estimated_transit_error_Mm3"])
        numerical &= transport_check["converged"] && transport_error<=storage_convergence
        numerical || push!(
            errors,
            "finer replay exceeded storage, level, transit or objective convergence tolerance",
        )
        ReplayResult(;
            valid=isempty(errors),
            original_valid=original["valid"],
            physically_valid=physical,
            numerically_converged=numerical,
            errors, violations, shortfalls,
            original_audit=original,
            transport_audit=transport_check,
            replay_grid=collect(grid),
            constraint_times=audit["storage_path_times"],
            storage_differences_Mm3=storage_errors,
            storage_convergence_tolerances_Mm3=storage_limits,
            level_differences_m=level_errors,
            terminal_transit_difference_Mm3=transit_error,
            replayed_objective,
            objective_difference=objective_error,
            seconds=time()-began,
            limits...,
        )
    catch e
        e isa Union{ArgumentError,DomainError} || rethrow()
        push!(errors, sprint(showerror, e))
        ReplayResult(;errors,violations,original_audit=original,original_valid=original["valid"],
            shortfalls,seconds=time()-began,limits...)
    end
end

"""Solve a fixed commitment with adaptive physical integration and finer replay.

Reconstruct and resolve starting controls before optimizing revenue. Local
step-doubling defects select hydraulic cells to bisect; coupled replay checks
the complete trajectory after every control change. Controls and operating
observations keep their original intervals. Finer witness equations enforce
restrictions that fail replay even after integration accuracy is resolved,
including analytically located storage extrema. After a mesh or constraint
change, restore feasibility before optimizing revenue. The best independently
accepted schedule survives later solve failures or budget exhaustion.

`initial` is checked as supplied and retains the audit-only contract; `warm`
is reconstructed under current inputs. `max_refinements=48` bounds mesh updates.
`time_limit` bounds each NLP; `deadline` bounds shared work between phases.
With `refine_grid=false`, keep the original hydraulic model and optionally add
finer witness constraints. Global preparation uses that finite-model scope.
`operational_margin` is an explicit power inset, never an automatic retry rule.
"""
solve_verified(c::ScheduleCase; kwargs...) = _verified_report(_verified_dispatch(c; kwargs...))

function _dispatch_reports(attempts)
    [Dict{String,Any}(string(k)=>(k==:replay ? _replay_report(v) : v)
        for (k,v) in pairs(attempt)) for attempt in attempts]
end

function _verified_report(result)
    SolveResult(Dict{String,Any}(
        "status"=>result.accepted ? "FEASIBLE" : "NO_ACCEPTED_SOLUTION",
        "accepted"=>result.accepted,
        "problem"=>result.case,
        "solution"=>result.solution,
        "audit"=>_replay_report(result.audit),
        "attempts"=>_dispatch_reports(result.attempts),
        "total_seconds"=>result.total_seconds,
        "feasible_lower_bound"=>result.accepted ? result.solution.objective : nothing,
        "acceptance_scope"=>"physical equations, operating limits and finer replay",
        "global_optimality_proven"=>false,
    ))
end

"""Resolve fixed controls; return a physical failure for the next NLP to repair."""
function resolve_trajectory(c,x;max_refinements,deadline,refine_grid,replay_factor,storage_limits,reconstructed,
    storage_convergence,storage_relative_convergence,level_convergence,objective_convergence)
    steps=Any[]
    refinements=0
    while true
        existing=river_routing_grid(c;physical_grid=x.grid)
        probe=refined_grid(x.grid;factor=replay_factor,existing)
        audit=_replay(c,x;grid=probe,storage_convergence,
            storage_relative_convergence,level_convergence,objective_convergence,
            error_factor=refine_grid ? replay_factor^2/(replay_factor^2-1) : 1.0,
            require_original=!reconstructed)
        push!(steps,(;intervals=length(x.grid)-1,objective=x.objective,accepted=audit.valid))
        original_valid=audit.original_valid
        if audit.valid || !reconstructed && !original_valid ||
           refinements==max_refinements || time()>=deadline
            return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
        end
        physical=hydraulic_case(c,x.grid)
        if audit.numerically_converged
            # Accuracy is resolved. A binding restriction must be enforced on
            # the checked mesh; moving the primary mesh alone can land on the
            # new bound and fail its next replay in exactly the same way.
            return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
        else
            routing_grid=try
                refine_transport_grid(physical,x,storage_limits,storage_convergence;deadline)
            catch error
                error isa Union{ArgumentError,DomainError} || rethrow()
                refine_grid || return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
                # A physical-domain failure of the finer probe is not a valid
                # trajectory. Build that finer problem and repair its controls.
                grid=probe
                next=_river_replace(c;routing_grid=replay_case(physical,grid).routing_grid)
                return (;case=next,solution=x,audit,grid,refinements=refinements+1,steps)
            end
            next=_river_replace(c;routing_grid)
            grid=x.grid
            if routing_grid==river_routing_grid(physical)
                refine_grid || return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
                indicator=hydraulic_defects(c,x;storage_limits,level_limit=level_convergence,
                    transport_limit=isfinite(storage_convergence) ? storage_convergence/10 : 2e-4,deadline)
                indicator===nothing && return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
                grid=bulk_refined_grid(x.grid,indicator;existing)
                grid==x.grid && (grid=probe)
            end
        end
        grid==x.grid && next.routing_grid==c.routing_grid &&
            return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
        time()<deadline || return (;case=c,solution=x,audit,grid=x.grid,refinements,steps)
        try
            x=controlled_trajectory(next,x;grid)
        catch error
            error isa Union{ArgumentError,DomainError} || rethrow()
            return (;case=next,solution=x,audit,grid,refinements=refinements+1,steps)
        end
        c=next
        refinements+=1
        reconstructed=true
    end
end

_verified_dispatch(c::ScheduleCase;kwargs...) = _verified_dispatch(c,solve_case;kwargs...)

function _verified_dispatch(
    c::ScheduleCase, dispatch_solver;
    u = ones(Int, length(c.system.generators), length(c.prices)),
    initial = nothing,
    max_refinements = 48,
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
    validate_inputs(c)
    admissible(c,u) || throw(ArgumentError("inadmissible commitment"))
    initial_schedule=initial===nothing ? nothing : _schedule_input(initial;case=c)
    initial_mesh=initial===nothing ? nothing : TimeDiscretization(c.grid,initial_schedule.grid)
    if initial!==nothing
        refine_grid || initial_schedule.grid==c.grid || throw(ArgumentError(
            "fixed-grid initial schedule must use the case's hydraulic grid"))
        size(initial_schedule.commitment)==(length(c.system.generators),length(initial_mesh.dt)) ||
            throw(ArgumentError("initial schedule commitment dimensions differ from its grid"))
        initial_schedule.commitment[:,first.(initial_mesh.cells)]==u ||
            throw(ArgumentError("initial schedule commitment differs from requested commitment"))
    end
    if river_on!==nothing
        river_on=fixed_river_on(c,river_on,nothing)
        if initial!==nothing
            initial_states=river_on_schedule(c,interval_means(initial_mesh,initial_schedule.river_release))
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
    storage_limits=storage_convergence_limits(c,storage_convergence,storage_relative_convergence)
    tolerances=(;storage_convergence,storage_relative_convergence,level_convergence,objective_convergence)
    best=nothing
    latest=nothing
    latest_audit=ReplayResult(errors=["dispatch deadline exhausted"])
    grid=copy(c.grid)
    replay_grid=nothing
    refinements=0
    # A supplied initial is checked as written. An unchecked start is first
    # reintegrated under the requested controls and current forecast.
    if initial!==nothing
        latest=initial_schedule
        grid=copy(latest.grid)
    elseif warm!==nothing && time()<deadline
        warm_grid=refine_grid ? get(warm,"grid",c.grid) : c.grid
        start,_=_dispatch_start(c,u,warm;grid=warm_grid)
        if start!==nothing
            try
                requested_opening=openings===nothing ? tunnel_opening_schedule(c) : openings
                latest=controlled_trajectory(c,start;u,opening=requested_opening,grid=warm_grid)
                grid=copy(latest.grid)
                if river_on!==nothing
                    states=river_on_schedule(c,operating_window_values(c,latest).river_release)
                    all(r.min_allowed_flow===nothing || states[i,:]==river_on[i,:]
                        for (i,r) in enumerate(c.system.rivers)) || (latest=nothing)
                end
            catch error
                error isa Union{ArgumentError,DomainError} || rethrow()
            end
        end
    end
    feasibility=feasibility_only || latest===nothing
    dispatch=nothing
    while true
        if latest===nothing
            remaining=min(time_limit,deadline-time())
            remaining>0 || break
            dispatch!==nothing && remaining<=get(dispatch,"construction_seconds",0.0) && break
            dispatch=dispatch_solver(c;u,solver,seed,warm,operational_margin,grid,replay_grid,
                time_limit=remaining,feasibility_only=feasibility,tunnel_opening=openings,river_on)
            latest=dispatch.solution
            latest===nothing && break
            if !dispatch.accepted
                try
                    latest,_=_reconstruct_candidate(c,latest)
                catch error
                    error isa Union{ArgumentError,DomainError} || rethrow()
                    break
                end
            end
        end
        resolution=resolve_trajectory(c,latest;max_refinements=max_refinements-refinements,
            deadline,refine_grid,replay_factor,storage_limits,reconstructed=initial===nothing,tolerances...)
        c=resolution.case
        latest=resolution.solution
        latest_audit=resolution.audit
        refinements+=resolution.refinements
        grid=copy(resolution.grid)
        replay_grid===nothing || (replay_grid=sort!(unique!(vcat(replay_grid,grid))))
        original=latest_audit.original_audit
        push!(attempts,(;
            iteration=length(attempts),intervals=length(c.prices),
            integration_intervals=length(latest.grid)-1,
            transport_intervals=length(river_routing_grid(c;physical_grid=grid))-1,
            operational_margin_MW=operational_margin,
            status=dispatch===nothing ? get(latest,"status","FORWARD_RECONSTRUCTED") : dispatch.status,
            objective=latest.objective,model_valid=latest_audit.original_valid,
            model_errors=get(original,"errors",String[]),replay=latest_audit,
            integration_steps=resolution.steps,
            solve_seconds=dispatch===nothing ? 0.0 : get(dispatch,"total_seconds",0.0),
            construction_seconds=dispatch===nothing ? 0.0 : get(dispatch,"construction_seconds",0.0),
            optimization_seconds=dispatch===nothing ? 0.0 : get(dispatch,"seconds",0.0),
            iterations=dispatch===nothing ? nothing : get(dispatch,"iterations",nothing)))
        if latest_audit.valid && (best===nothing || latest.objective>best.solution.objective)
            best=(;case=c,solution=latest,audit=latest_audit)
        end
        # Initial schedules retain the audit-only contract. A feasible numerical
        # start or feasibility solve is followed by one revenue solve.
        initial!==nothing && latest_audit.valid && break
        initial!==nothing && !latest_audit.original_valid && resolution.refinements==0 && break
        time()<deadline || break
        if latest_audit.valid
            feasibility_only && break
            dispatch!==nothing && !feasibility && break
            feasibility=false
        elseif grid!=latest.grid
            feasibility=true
        elseif !latest_audit.original_valid &&
               (!latest_audit.numerically_converged || refinements==max_refinements)
            # Exhausting the mesh budget still permits repair of restrictions
            # already represented by the primary model.
            dispatch!==nothing && resolution.refinements==0 && break
            feasibility=true
        else
            refinements<max_refinements || break
            if !refine_grid || latest_audit.numerically_converged
                parent=replay_grid===nothing ? grid : replay_grid
                replay_grid=refined_grid(parent;factor=replay_factor,
                    existing=river_routing_grid(c;physical_grid=parent))
                # A quadratic storage path can violate a bound between mesh
                # vertices. Enforce its observed extremum in the next solve.
                replay_grid=sort!(unique!(vcat(replay_grid,
                    latest_audit.constraint_times)))
                refinements+=1
                feasibility=true
            else
                break
            end
        end
        initial=nothing
        warm=latest
        latest=nothing
    end
    (;accepted=best!==nothing,case=best===nothing ? c : best.case,
        solution=best===nothing ? latest : best.solution,
        audit=best===nothing ? latest_audit : best.audit,attempts,total_seconds=time()-began)
end
