include("starts.jl")

function _remove_constant_constraints!(model)
    count = 0
    for ref in all_constraints(model; include_variable_in_set_constraints = false)
        obj = constraint_object(ref)
        f = obj.func
        x = if f isa Number
            Float64(f)
        elseif f isa GenericAffExpr && isempty(f.terms)
            Float64(f.constant)
        elseif f isa GenericQuadExpr && isempty(f.terms) && isempty(f.aff.terms)
            Float64(f.aff.constant)
        else
            nothing
        end
        x === nothing && continue
        set = obj.set
        feasible =
            set isa MOI.EqualTo ? x == set.value :
            set isa MOI.LessThan ? x <= set.upper :
            set isa MOI.GreaterThan ? x >= set.lower :
            set isa MOI.Interval ? set.lower <= x <= set.upper : false
        feasible || throw(ArgumentError("infeasible constant constraint"))
        delete(model, ref)
        count += 1
    end
    count
end

function _reconstruct_candidate(c, raw; transport = nothing)
    admissible(c, raw["u"]) || throw(ArgumentError("inadmissible commitment"))
    q, gate = copy(raw["generator_q"]), copy(raw["gate"])
    all(isfinite, q) && all(isfinite, gate) || throw(ArgumentError("nonfinite controls"))
    all(x -> x >= -1e-6, q) || throw(ArgumentError("negative discharge"))
    all(x -> -1e-8 <= x <= 1 + 1e-8, gate) || throw(ArgumentError("gate outside bounds"))
    q = max.(q, 0.0)
    gate = clamp.(gate, 0.0, 1.0)
    correction = (
        flow = maximum(abs, q .- raw["generator_q"]; init = 0.0),
        gate = maximum(abs, gate .- raw["gate"]; init = 0.0),
    )
    x = dispatch_from_controls(c, raw["u"], q, gate; transport)
    x, correction
end

function _set_incumbent!(result, candidate, source, relative_gap, absolute_gap)
    result["solution"] = candidate
    result["incumbent_source"] = source
    result["validation"] = candidate["validation"]
    result["objective"] = candidate["objective"]
    result["feasible_lower_bound"] = candidate["objective"]
    result["relative_gap"] = nothing
    result["absolute_gap"] = nothing
    result["global_certificate"] = false
    upper, lower = result["global_bound"], candidate["objective"]
    if upper !== nothing && upper >= lower - 1e-6
        gap = max(0.0, upper - lower)
        result["absolute_gap"] = gap
        result["relative_gap"] = gap / max(1.0, abs(lower))
        result["global_certificate"] =
            gap <= absolute_gap || result["relative_gap"] <= relative_gap
    end
    result
end

"""
    solve(case; time_limit=60.0, relative_gap=1e-3, absolute_gap=0.0, initial=nothing, fixed_u=nothing)

Optimize generation and binary unit commitment with native SCIP. The objective
is revenue minus transition costs and release penalties, plus changes in stored
water value. `initial` may supply a physical schedule on the same control grid.

`solution` is independently reconstructed and equation-validated. `accepted`
also requires finer chronological replay. `global_certificate` requires an
enclosing solver upper bound within either requested gap. Its scope is the
declared discrete hydraulic model, under numerical solver tolerances. Supplying
`fixed_u` makes the certificate conditional on that commitment matrix.

The time allowance includes construction and optimization; extraction/replay
can overrun it and are included in `total_seconds`. Julia compilation is also
included if this is the first call in a process.
"""
function solve(
    c::ScheduleCase;
    time_limit = 60.0,
    relative_gap = 1e-3,
    absolute_gap = 0.0,
    initial = nothing,
    fixed_u = nothing,
    replay = true,
)
    isfinite(time_limit) && time_limit > 0 ||
        throw(ArgumentError("positive finite time_limit required"))
    isfinite(relative_gap) && 0 <= relative_gap < 1 ||
        throw(ArgumentError("relative_gap must lie in [0,1)"))
    isfinite(absolute_gap) && absolute_gap >= 0 ||
        throw(ArgumentError("nonnegative finite absolute_gap required"))
    began = time()
    result = Dict{String,Any}(
        "case" => c.name,
        "solver" => "SCIP",
        "status" => "CONSTRUCTION_BUDGET_EXHAUSTED",
        "accepted" => false,
        "global_certificate" => false,
        "solution" => nothing,
        "feasible_lower_bound" => nothing,
        "global_bound" => nothing,
        "relative_gap" => nothing,
        "absolute_gap" => nothing,
        "commitment_fixed" => fixed_u !== nothing,
        "certificate_scope" =>
            fixed_u === nothing ?
            "declared discrete midpoint model; numerical global bound" :
            "declared discrete midpoint model with fixed commitment; numerical global bound",
        "transport_exact" =>
            all(r -> r.deterministic_delay !== nothing, c.system.rivers),
        "requested_relative_gap" => relative_gap,
        "requested_absolute_gap" => absolute_gap,
    )
    best = nothing
    if initial !== nothing
        initial_audit = validate(c, initial)
        initial_audit["valid"] ||
            throw(ArgumentError("initial schedule fails the equation audit"))
        best, correction = _reconstruct_candidate(c, initial)
        best["validation"]["valid"] ||
            throw(ArgumentError("initial controls fail reconstruction audit"))
        result["initial_objective"] = best["objective"]
    end
    b = _build_global_dispatch(c; joint = true, warm = best, fixed_u)
    result["removed_constant_constraints"] = _remove_constant_constraints!(b.m)
    if best !== nothing
        start_audit = _lift_start!(b, c, best)
        result["start_audit"] = start_audit
        # A partial or invalid lift must not be submitted as a feasible start.
        if !start_audit["valid"]
            for variable in all_variables(b.m)
                set_start_value(variable, nothing)
            end
            best = nothing
            result["initial_rejected_by_model"] = true
        end
    end
    validated_initial = best
    result["variable_count"] = num_variables(b.m)
    result["constraint_count"] =
        num_constraints(b.m; count_variable_in_set_constraints = false)
    result["construction_seconds"] = time() - began
    remaining = time_limit - result["construction_seconds"]
    if remaining > 0
        set_optimizer(
            b.m,
            optimizer_with_attributes(
                SCIP.Optimizer,
                "display/verblevel" => 0,
                "limits/time" => remaining,
                "limits/gap" => relative_gap,
                "limits/absgap" => absolute_gap / 10000,
                "numerics/feastol" => 1e-8,
                "parallel/maxnthreads" => 1,
            ),
        )
        result["solve_seconds"] = @elapsed try
            optimize!(b.m)
        catch error
            result["solver_error"] = sprint(showerror, error)
        end
        result["status"] = string(termination_status(b.m))
        result["primal_status"] = string(primal_status(b.m))
        upper = try
            objective_bound(b.m) * 10000
        catch
            NaN
        end
        result["global_bound"] = isfinite(upper) ? upper : nothing
        if result_count(b.m) > 0
            try
                raw = _dispatch_values(b)
                fractional = value.(b.u)
                raw["fractional_commitment"] = fractional
                raw["u"] = round.(Int, fractional)
                error = maximum(abs, fractional .- raw["u"]; init = 0.0)
                raw["validation"] = validate(c, raw; transport = b.transport)
                result["raw_solver_solution"] = raw
                result["raw_integrality_error"] = error
                if error <= 1e-6 && (fixed_u === nothing || raw["u"] == fixed_u)
                    x, correction = _reconstruct_candidate(c, raw; transport = b.transport)
                    result["control_boundary_correction"] =
                        Dict("flow" => correction.flow, "gate" => correction.gate)
                    if x["validation"]["valid"] &&
                       (best === nothing || x["objective"] > best["objective"])
                        best = x
                        result["incumbent_source"] = "solver"
                    end
                end
            catch error
                result["candidate_error"] = sprint(showerror, error)
            end
        end
    end
    if best !== nothing
        source = get(result, "incumbent_source", "initial")
        _set_incumbent!(result, best, source, relative_gap, absolute_gap)
        if replay
            result["replay_audit"] = replay_audit(c, best; transport = b.transport)
            result["accepted"] = result["replay_audit"]["valid"]
            if !result["accepted"] && source == "solver" && validated_initial !== nothing
                initial_replay = replay_audit(c, validated_initial; transport = b.transport)
                if initial_replay["valid"]
                    result["discrete_candidate"] = Dict(
                        "solution" => best,
                        "objective" => result["objective"],
                        "relative_gap" => result["relative_gap"],
                        "absolute_gap" => result["absolute_gap"],
                        "global_certificate" => result["global_certificate"],
                        "replay_audit" => result["replay_audit"],
                    )
                    _set_incumbent!(
                        result,
                        validated_initial,
                        "initial",
                        relative_gap,
                        absolute_gap,
                    )
                    result["replay_audit"] = initial_replay
                    result["accepted"] = true
                end
            end
        end
    end
    result["total_seconds"] = time() - began
    result["budget_overrun_seconds"] = max(0.0, result["total_seconds"] - time_limit)
    result
end

export solve
