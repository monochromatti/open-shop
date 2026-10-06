include("starts.jl")

# Native SCIP statistics do not supply the audited incumbent. A conservative
# unresolved-root guard may withhold a certificate. Bound values use physical
# objective units; the optional logfile retains the progress trajectory.
function _scip_bound_value(optimizer, value; scale=10000.)
    isfinite(value) && !Bool(SCIP.SCIPisInfinity(optimizer, abs(value))) ? value*scale : nothing
end

function _scip_diagnostics(model)
    result=Dict{String,Any}("source"=>"SCIP native statistics; objective bounds in physical units")
    try
        optimizer=JuMP.unsafe_backend(model)
        optimizer isa SCIP.Optimizer || return merge(result,Dict("available"=>false))
        stage=SCIP.SCIPgetStage(optimizer)
        result["stage"]=string(stage)
        result["raw_status"]=string(SCIP.SCIPgetStatus(optimizer))
        # Native bound/statistic APIs are queried only after transformation.
        if Int(stage) in 3:10
            result["available"]=true
            result["nodes"]=MOI.get(optimizer,MOI.NodeCount())
            result["total_nodes"]=SCIP.SCIPgetNTotalNodes(optimizer)
            result["lp_iterations"]=SCIP.SCIPgetNLPIterations(optimizer)
            result["solutions_stored"]=SCIP.SCIPgetNSols(optimizer)
            result["solutions_found"]=SCIP.SCIPgetNSolsFound(optimizer)
            result["native_solve_seconds"]=SCIP.SCIPgetSolvingTime(optimizer)
            result["final_upper_bound"]=_scip_bound_value(optimizer,SCIP.SCIPgetDualbound(optimizer))
            result["native_incumbent_objective"]=_scip_bound_value(optimizer,SCIP.SCIPgetPrimalbound(optimizer))
            if Int(stage)>=5
                root=SCIP.SCIPgetDualboundRoot(optimizer)
                result["root_upper_bound"]=_scip_bound_value(optimizer,root)
                first=SCIP.SCIPgetFirstLPDualboundRoot(optimizer)
                result["first_root_lp_upper_bound"]=_scip_bound_value(optimizer,first)
            end
            result["root_bound_note"]="SCIP original-problem root bound; unavailable/sentinel is null, not a proof"
        else
            result["available"]=false
        end
    catch error
        result["error"]=sprint(showerror,error)
    end
    result
end


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

function _unresolved_root_certificate(diagnostics, status)
    status=="OPTIMAL" && get(diagnostics,"available",false) &&
    get(diagnostics,"lp_iterations",0)>0 &&
    haskey(diagnostics,"first_root_lp_upper_bound") &&
    diagnostics["first_root_lp_upper_bound"]===nothing
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

Optimize generation and binary unit commitment with native SCIP.
`formulation` selects `:tensor` (default), `:baseline`, `:domains` or `:tightened`; all retain the
same physical equations. `diagnostics_path` optionally writes a native SCIP
progress log. Returned `scip_diagnostics` are observational statistics, not
independent feasibility or certificate evidence. The objective
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
    formulation = :tensor,
    diagnostics_path = nothing,
)
    isfinite(time_limit) && time_limit > 0 ||
        throw(ArgumentError("positive finite time_limit required"))
    isfinite(relative_gap) && 0 <= relative_gap < 1 ||
        throw(ArgumentError("relative_gap must lie in [0,1)"))
    isfinite(absolute_gap) && absolute_gap >= 0 ||
        throw(ArgumentError("nonnegative finite absolute_gap required"))
    base_formulation=endswith(string(formulation),"_flow") ? Symbol(chop(string(formulation);tail=5)) : formulation
    base_formulation in (:baseline, :domains, :tightened, :tensor, :cartesian_ranges, :cartesian_cuts, :cartesian_refined, :tensor_pruned, :tensor_quadratic, :tensor_refined) || throw(ArgumentError("unknown formulation"))
    diagnostics_path!==nothing && (diagnostics_path=abspath(String(diagnostics_path)))
    began = time()
    result = Dict{String,Any}(
        "case" => c.name,
        "solver" => "SCIP",
        "formulation" => string(formulation),
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
    b = _build_global_dispatch(c; joint = true, warm = best, fixed_u, formulation)
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
    variables=all_variables(b.m)
    binary_names=[name(v) for v in variables if is_binary(v)]
    result["model_profile"]=Dict(
        "binary_vars"=>length(binary_names),
        "sos2_constraints"=>count(ref->constraint_object(ref).set isa MOI.SOS2,all_constraints(b.m;include_variable_in_set_constraints=false)),
        "tensor_axes"=>length(get(b.m.ext,:global_tensor_coordinates,[])),
        "tensor_axis_knots"=>sum(length(x.nodes) for x in get(b.m.ext,:global_tensor_coordinates,[]);init=0),
        "tunnel_direction_binaries"=>count(n->startswith(n,"tunnel_direction_"),binary_names),
        "turbine_cell_binaries"=>count(n->startswith(n,"turbine_") && occursin("_cell[",n),binary_names),
        "river_table_cell_binaries"=>count(n->startswith(n,"river_law_") && occursin("_cell[",n),binary_names),
        "bounds_seconds"=>b.bounds_seconds,
        "tightening_passes"=>b.domains!==nothing && hasproperty(b.domains,:tightening_passes) ? b.domains.tightening_passes : 0,
        "shared_head_entries"=>length(b.shared_heads),
        "head_variables_saved"=>length(c.system.generators)*length(c.prices)-length(b.shared_heads))
    result["variable_count"] = num_variables(b.m)
    result["constraint_count"] =
        num_constraints(b.m; count_variable_in_set_constraints = false)
    result["construction_seconds"] = time() - began
    remaining = time_limit - result["construction_seconds"]
    if remaining > 0
        factory=SCIP.Optimizer
        if diagnostics_path!==nothing
            # The builder's cached MOI.Silent would otherwise override verbosity.
            unset_silent(b.m)
            mkpath(dirname(diagnostics_path))
            isfile(diagnostics_path) && rm(diagnostics_path)
            result["diagnostics_path"]=diagnostics_path
        end
        set_optimizer(
            b.m,
            optimizer_with_attributes(
                factory,
                "display/verblevel" => (diagnostics_path===nothing ? 0 : 4),
                "display/freq" => 100,
                "limits/time" => remaining,
                "limits/gap" => relative_gap,
                "limits/absgap" => absolute_gap / 10000,
                "numerics/feastol" => 1e-8,
                "parallel/maxnthreads" => 1,
            ),
        )
        result["solve_seconds"] = @elapsed try
            # Copy resets SCIP's native instance; attach before installing its log.
            JuMP.MOI.Utilities.attach_optimizer(JuMP.backend(b.m))
            remaining_after_copy=max(0.,time_limit-(time()-began))
            set_optimizer_attribute(b.m,"limits/time",remaining_after_copy)
            if diagnostics_path!==nothing
                SCIP.SCIPsetMessagehdlrLogfile(JuMP.unsafe_backend(b.m),diagnostics_path)
            end
            optimize!(b.m)
        catch error
            result["solver_error"] = sprint(showerror, error)
        end
        result["status"] = string(termination_status(b.m))
        result["primal_status"] = string(primal_status(b.m))
        result["scip_diagnostics"]=_scip_diagnostics(b.m)
        if diagnostics_path!==nothing
            try
                SCIP.SCIPsetMessagehdlrLogfile(JuMP.unsafe_backend(b.m),C_NULL)
            catch error
                result["diagnostics_close_error"]=sprint(showerror,error)
            end
        end
        upper = try
            _scip_bound_value(JuMP.unsafe_backend(b.m),objective_bound(b.m))
        catch
            nothing
        end
        if _unresolved_root_certificate(result["scip_diagnostics"],result["status"])
            result["bound_rejection"]="SCIP terminated without a finite first root LP bound despite LP iterations; native bound is retained only in diagnostics"
            upper=nothing
        end
        result["global_bound"] = upper
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
                sos2_error=maximum((_sos2_residual(JuMP.value.(constraint_object(ref).func),constraint_object(ref).set.weights)
                    for ref in all_constraints(b.m;include_variable_in_set_constraints=false)
                    if constraint_object(ref).set isa MOI.SOS2);init=0.0)
                result["raw_sos2_residual"]=sos2_error
                if error <= 1e-6 && sos2_error <= 1e-6 && (fixed_u === nothing || raw["u"] == fixed_u)
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
