"""Subdivide integration intervals while retaining every control-change time."""
function refined_grid(grid; factor = 2)
    factor isa Integer && factor>=1 ||
        throw(ArgumentError("positive integer refinement required"))
    vcat(
        [
            collect(range(grid[t], grid[t + 1]; length = factor+1))[1:(end - 1)] for
            t in 1:(length(grid) - 1)
        ]...,
        last(grid),
    )
end

"""Accept a schedule only after finer physical replay and an independent transport check.
Operating data retain physical times; plant ramp limits retain original decision windows.
"""
replay_audit(c::ScheduleCase, x; kwargs...) = _replay_report(_replay(c, x; kwargs...))

_replay_report(check) = Dict{String,Any}(string(k)=>v for (k,v) in pairs(check))

# Numerical replay returns named quantities. Diagnostic labels never determine
# a subsequent solve's margins or its physical acceptance.
function _replay(
    c::ScheduleCase,
    x;
    grid = refined_grid(c.grid),
    tolerance = 1e-4,
    storage_convergence = 2e-3,
    objective_convergence = 1e-3,
    transport = nothing,
)
    began=time()
    errors=String[]
    violations=Dict{String,Float64}()
    shortfalls=(arrival=0.0,power=0.0,discharge=0.0)
    nd=all(r->r.deterministic_delay!==nothing, c.system.rivers) ?
       _transport_data(c, transport) : nothing
    nd===nothing &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed acceptance",
            ),
        )
    original=validate(c, x; transport = nd)
    original["valid"] || push!(errors, "original equation audit failed")
    length(grid)>length(c.grid) ||
        push!(errors, "acceptance requires a strictly finer replay grid")
    if !isempty(errors)
        return (; valid=false, errors, violations, original_audit=original,
            shortfalls, seconds=time()-began)
    end
    x=_schedule_input(x; case=c)
    try
        cc=with_grid(c, grid)
        fine=without_interval_controls(cc)
        indices=[
            clamp(searchsortedlast(c.grid, (grid[t]+grid[t + 1])/2), 1, length(c.prices))
            for t in 1:(length(grid) - 1)
        ]
        opening=_physical_tunnel_openings(c,get(x,"tunnel_opening",nothing))
        z=dispatch_from_controls(
            fine,
            x.commitment[:, indices],
            x.generator_q[:, indices],
            x.gate[:, indices];
            tunnel_opening=opening[:,indices],
            reserves=_reserve_slice(x.reserves,indices),
            head_reference=x.head[:,indices],
        )
        audit=z["validation"]
        merge!(violations, audit["residuals"])
        append!(errors, audit["errors"])
        # Reserve response is measured from each operating point, so its rates
        # remain physical limits when ordinary scheduling ramps are stripped.
        reserve_audit=validate_reserves(cc,z;tolerance)
        merge!(violations,Dict("reserve_"*key=>value for (key,value) in reserve_audit["residuals"]))
        append!(errors,reserve_audit["errors"])
        windows=operating_window_values(c,z,grid,indices,x.commitment)
        plants=plant_operating_residuals(c,windows)
        for (name,residual) in operating_residuals(c,windows;plants)
            violations[name]=max(get(violations,name,0.0),residual)
            residual>tolerance && push!(errors,"$name fails finer replay: $residual")
        end
        fine_shortfalls=audit["shortfalls"]
        for (i,g) in enumerate(c.system.generators)
            violations["unit_power_$(g.name)"]=fine_shortfalls.unit_power[i]
        end
        transport_check=transport_audit(
            c,
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
                maximum(fine_shortfalls.plant_power;init=0.0),
                maximum(plants.power_capacity;init=0.0),maximum(plants.power_ramp;init=0.0)),
            discharge=max(maximum(fine_shortfalls.unit_discharge;init=0.0),
                maximum(fine_shortfalls.plant_discharge;init=0.0),
                maximum(plants.discharge_capacity;init=0.0),maximum(plants.discharge_ramp;init=0.0)),
        )
        transport_check["accepted"] ||
            push!(errors, "transport convergence, operating limits or schedule discrepancy fail acceptance")
        physical=isempty(errors)
        edges=[searchsortedfirst(grid, t) for t in c.grid]
        storage_error=maximum(abs, z.storage[:, edges]-x.storage; init = 0.0)
        transit_error=maximum(abs, z.terminal_transit-x.terminal_transit; init = 0.0)
        replayed_objective=z.objective-river_operating_cost(c,windows.river_release).total
        objective_error=replayed_objective-x.objective
        numerical=storage_error<=storage_convergence &&
                  transit_error<=storage_convergence &&
                  abs(objective_error)<=objective_convergence*max(1.0, abs(x.objective))
        numerical || push!(
            errors,
            "finer replay exceeded storage/transit or objective convergence tolerance",
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
            terminal_transit_difference_Mm3=transit_error,
            replayed_objective,
            objective_difference=objective_error,
            storage_convergence_tolerance_Mm3=storage_convergence,
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

After a failed replay, add the finer hydraulic equations to the nonlinear model,
sharing the original unit discharges, outlet openings and reserve allocations.
Fixed-speed pump discharges remain hydraulic states. Operating schedules and ramps
retain their original windows. No measured residual becomes a blanket margin.
The returned case, schedule and objective retain the requested grid.
`refine_grid` is retained for source compatibility; both values now retain controls.
`operational_margin` is an explicit optional power inset, never an automatic retry rule.
"""
solve_verified(c::ScheduleCase; kwargs...) = _verified_report(_verified_dispatch(c; kwargs...))

function _verified_report(result)
    Dict{String,Any}(
        "accepted"=>result.accepted,
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
)
    isfinite(time_limit) && time_limit>0 ||
        throw(ArgumentError("positive finite per-solve time limit required"))
    isnan(deadline) && throw(ArgumentError("invalid dispatch deadline"))
    replay_factor isa Integer && replay_factor>=2 ||
        throw(ArgumentError("acceptance requires replay_factor >= 2"))
    max_refinements isa Integer && max_refinements>=0 ||
        throw(ArgumentError("nonnegative integer refinement budget required"))
    refine_grid isa Bool || throw(ArgumentError("refine_grid must be boolean"))
    initial!==nothing &&
        (!haskey(initial, "u") || initial["u"]!=u) &&
        throw(
            ArgumentError("initial candidate commitment differs from requested commitment"),
        )
    initial!==nothing && haskey(initial, "grid") && initial["grid"]!=c.grid &&
        throw(ArgumentError("initial candidate control grid differs from requested grid"))
    if river_on!==nothing
        river_on=fixed_river_on(c,river_on,nothing)
        if initial!==nothing
            initial_states=river_on_schedule(c,_schedule_input(initial;case=c).river_release)
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
        _physical_tunnel_openings(c,get(initial,"tunnel_opening",nothing))==openings ||
            throw(ArgumentError("initial candidate tunnel openings differ from requested openings"))
    end
    began=time()
    attempts=Any[]
    integration_grid=nothing
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
        x=dispatch isa SolveResult ? dispatch.solution : dispatch
        audit=x!==nothing && haskey(x,"objective") ?
              _replay(c, x; grid = integration_grid===nothing ? refined_grid(c.grid;factor=replay_factor) : integration_grid) :
              (valid=false,errors=["NLP candidate has no objective"],shortfalls=(arrival=0.0,power=0.0,discharge=0.0))
        original=get(audit, :original_audit, Dict{String,Any}())
        physical=get(original, "valid", false)
        push!(
            attempts,
            Dict(
                "iteration"=>iteration,
                "intervals"=>length(c.prices),
                "integration_intervals"=>integration_grid===nothing ? length(c.prices) : length(integration_grid)-1,
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
        physical || break # A failed local solve is not an infeasibility proof.
        # Constrain the mesh that exposed the error. A subsequent failure may
        # justify a finer state mesh, but never new controls or guessed buffers.
        integration_grid=integration_grid===nothing ? refined_grid(c.grid;factor=replay_factor) :
            refined_grid(integration_grid;factor=replay_factor)
        warm=x
    end
    (;accepted=audit.valid,case=c,solution=x,audit,attempts,total_seconds=time()-began)
end
