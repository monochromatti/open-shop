"""Refine every control interval without moving a control-change time."""
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
function replay_audit(
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
        return Dict(
            "valid"=>false,
            "errors"=>errors,
            "violations"=>violations,
            "original_audit"=>original,
            "seconds"=>time()-began,
        )
    end
    try
        cc=with_grid(c, grid)
        fine=without_interval_controls(cc)
        indices=[
            clamp(searchsortedlast(c.grid, (grid[t]+grid[t + 1])/2), 1, length(c.prices))
            for t in 1:(length(grid) - 1)
        ]
        z=dispatch_from_controls(
            fine,
            x["u"][:, indices],
            x["generator_q"][:, indices],
            x["gate"][:, indices],
        )
        audit=z["validation"]
        merge!(violations, audit["residuals"])
        append!(errors, audit["errors"])
        windows=operating_window_values(c,z,grid,indices,x["u"])
        for (name,residual) in operating_residuals(c,windows)
            violations[name]=max(get(violations,name,0.0),residual)
            residual>tolerance && push!(errors,"$name fails finer replay: $residual")
        end
        for g in c.system.generators
            violations["unit_power_$(g.name)"]=max(
                get(violations, "unit_bounds_$(g.name)", 0.0),
                0.00981*g.hbest*max(g.min_efficiency, g.efficiency, eps(Float64))*get(
                    violations,
                    "turbine_envelope_$(g.name)",
                    0.0,
                ),
            )
        end
        for r in c.system.rivers
            violations["environmental_$(r.name)"]=max(
                get(violations, "environmental_arrival_$(r.name)", 0.0),
                get(violations, "pointwise_arrival_$(r.name)", 0.0),
                get(violations, "exact_arrival_$(r.name)", 0.0),
            )
        end
        transport_check=transport_audit(
            c,
            x;
            schedule_tolerance = storage_convergence,
            absolute_tolerance = isfinite(storage_convergence) ? storage_convergence/10 :
                                 2e-4,
            transport = nd,
        )
        transport_check["accepted"] ||
            push!(errors, "transport convergence or schedule discrepancy exceeds tolerance")
        physical=isempty(errors)
        edges=[searchsortedfirst(grid, t) for t in c.grid]
        storage_error=maximum(abs, z["V"][:, edges]-x["V"]; init = 0.0)
        transit_error=maximum(abs, z["terminal_transit"]-x["terminal_transit"]; init = 0.0)
        objective_error=z["objective"]-x["objective"]
        numerical=storage_error<=storage_convergence &&
                  transit_error<=storage_convergence &&
                  abs(objective_error)<=objective_convergence*max(1.0, abs(x["objective"]))
        numerical || push!(
            errors,
            "finer replay exceeded storage/transit or objective convergence tolerance",
        )
        Dict(
            "valid"=>isempty(errors),
            "physically_valid"=>physical,
            "numerically_converged"=>numerical,
            "errors"=>errors,
            "violations"=>violations,
            "original_audit"=>original,
            "transport_audit"=>transport_check,
            "replay_grid"=>collect(grid),
            "replay_intervals"=>length(grid)-1,
            "max_water_residual_Mm3"=>get(violations, "global_water", Inf),
            "max_storage_difference_Mm3"=>storage_error,
            "terminal_transit_difference_Mm3"=>transit_error,
            "replayed_objective"=>z["objective"],
            "objective_difference"=>objective_error,
            "storage_convergence_tolerance_Mm3"=>storage_convergence,
            "relative_objective_convergence_tolerance"=>objective_convergence,
            "seconds"=>time()-began,
            "certificate"=>"numerical finer-grid and transport acceptance only",
        )
    catch e
        push!(errors, sprint(showerror, e))
        Dict(
            "valid"=>false,
            "errors"=>errors,
            "violations"=>violations,
            "original_audit"=>original,
            "seconds"=>time()-began,
        )
    end
end

"""Solve a fixed commitment, refine/reoptimize after failed independent replay.

The original physical environmental minimum is unchanged. An explicit arrival
margin can be introduced after a failed replay to cover the measured numerical
shortfall. `operational_margin` requests initial power/envelope tightening
when optimizing a candidate; a supplied `initial` is audited as supplied.
All attempts, margins and actual solution grids are returned.
"""
function solve_verified(
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
    operational_margin = 0.0,
)
    isfinite(time_limit) && time_limit>0 ||
        throw(ArgumentError("positive finite per-solve time limit required"))
    replay_factor isa Integer && replay_factor>=2 ||
        throw(ArgumentError("acceptance requires replay_factor >= 2"))
    max_refinements>=0 || throw(ArgumentError("negative refinement budget"))
    initial!==nothing &&
        (!haskey(initial, "u") || initial["u"]!=u) &&
        throw(
            ArgumentError("initial candidate commitment differs from requested commitment"),
        )
    initial!==nothing && haskey(initial, "grid") && initial["grid"]!=c.grid &&
        throw(ArgumentError("initial candidate control grid differs from requested grid"))
    isfinite(operational_margin) && operational_margin >= 0 ||
        throw(ArgumentError("nonnegative finite operational_margin required"))
    began=time()
    attempts=Any[]
    cc=c
    uu=copy(u)
    margin=0.0
    power_margin=operational_margin
    x=initial
    audit=nothing
    for iteration in 0:max_refinements
        x=iteration==0 && initial!==nothing ? initial :
          solve_case(
            cc;
            u = uu,
            solver = solver,
            seed = seed,
            warm = warm,
            arrival_margin = margin,
            operational_margin = power_margin,
            time_limit = time_limit,
            feasibility_only = feasibility_only,
        )
        audit=haskey(x, "objective") ?
              replay_audit(cc, x; grid = refined_grid(cc.grid; factor = replay_factor)) :
              Dict("valid"=>false, "errors"=>["NLP candidate has no objective"])
        original=get(audit, "original_audit", Dict{String,Any}())
        physical=get(original, "valid", false)
        push!(
            attempts,
            Dict(
                "iteration"=>iteration,
                "intervals"=>length(cc.prices),
                "margin_m3s"=>margin,
                "operational_margin_MW"=>power_margin,
                "status"=>x["status"],
                "objective"=>get(x, "objective", nothing),
                "model_valid"=>physical,
                "model_errors"=>get(original, "errors", String[]),
                "replay"=>audit,
                "solve_seconds"=>get(x, "total_seconds", 0.0),
                "construction_seconds"=>get(x, "construction_seconds", 0.0),
                "optimization_seconds"=>get(x, "seconds", 0.0),
                "validation_seconds"=>get(x, "validation_seconds", 0.0),
                "forward_reconstruction_seconds"=>get(x, "forward_reconstruction_seconds", 0.0),
                "forward_reconstructed"=>get(x, "forward_reconstructed", false),
                "iterations"=>get(x, "iterations", nothing),
                "warm_start"=>get(x, "warm_start", nothing),
            ),
        )
        get(audit, "valid", false) && break
        iteration==max_refinements && break
        physical || break # A failed local solve is not an infeasibility proof.
        shortfall=maximum(
            (
                get(get(audit, "violations", Dict()), "environmental_$(r.name)", 0.0) for
                r in cc.system.rivers if r.arrival_policy==:pointwise
            );
            init = 0.0,
        )
        shortfall>1e-4 && (margin=max(margin, shortfall+0.02))
        power_shortfall=maximum(
            (
                value for (name, value) in get(audit, "violations", Dict()) if
                startswith(name, "plant_") || startswith(name, "unit_power_")
            );
            init = 0.0,
        )
        power_shortfall>1e-4 && (power_margin=max(power_margin, 1.5power_shortfall+0.002))
        nextgrid=refined_grid(cc.grid)
        idx=[
            min(
                searchsortedlast(cc.grid, (nextgrid[t]+nextgrid[t + 1])/2),
                length(cc.prices),
            ) for t in 1:(length(nextgrid) - 1)
        ]
        uu=uu[:, idx]
        cc=with_grid(cc, nextgrid)
        warm=try
            dispatch_from_controls(cc, uu, x["generator_q"][:, idx], x["gate"][:, idx])
        catch
            nothing
        end
    end
    Dict(
        "accepted"=>get(audit, "valid", false),
        "case"=>cc,
        "solution"=>x,
        "audit"=>audit,
        "attempts"=>attempts,
        "total_seconds"=>time()-began,
        "global_optimality_proven"=>false,
    )
end
