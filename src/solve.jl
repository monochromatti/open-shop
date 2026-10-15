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
            result["presolving_seconds"]=SCIP.SCIPgetPresolvingTime(optimizer)
            result["first_lp_seconds"]=SCIP.SCIPgetFirstLPTime(optimizer)
            result["last_node_separation_rounds"]=SCIP.SCIPgetNSepaRounds(optimizer)
            result["cuts_applied"]=SCIP.SCIPgetNCutsApplied(optimizer)
            result["final_upper_bound"]=_scip_bound_value(optimizer,SCIP.SCIPgetDualbound(optimizer))
            result["native_incumbent_objective"]=_scip_bound_value(optimizer,SCIP.SCIPgetPrimalbound(optimizer))
            result["native_incumbent_available"]=result["solutions_found"]>0 &&
                result["native_incumbent_objective"]!==nothing
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


# Native JSON contains timing and work counters for LPs and individual plugins.
# Capture only when diagnostics were requested; it does not affect certificates.
function _write_scip_statistics(model, path)
    file=ccall(:fopen, Ptr{Cvoid}, (Cstring, Cstring), path, "w")
    file==C_NULL && error("cannot open SCIP statistics file: $path")
    try
        SCIP.@SCIP_CALL SCIP.SCIPprintStatisticsJson(JuMP.unsafe_backend(model), file)
    finally
        ccall(:fclose, Cint, (Ptr{Cvoid},), file)
    end
    JSON3.read(read(path, String), Dict{String,Any})
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

function _set_incumbent!(result, candidate, source, relative_gap, absolute_gap)
    result["solution"] = candidate
    result["incumbent_source"] = source
    result["validation"] = candidate["validation"]
    result["objective"] = candidate.objective
    result["feasible_lower_bound"] = candidate.objective
    result["relative_gap"] = nothing
    result["absolute_gap"] = nothing
    result["global_certificate"] = false
    upper, lower = result["global_bound"], candidate.objective
    if upper !== nothing && upper < lower - 1e-6
        result["rejected_global_bound"] = upper
        result["bound_rejection"] = "solver upper bound below independently reconstructed objective"
        result["global_bound"] = nothing
        upper = nothing
    end
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
    solve(case; time_limit=60.0, relative_gap=1e-3, absolute_gap=0.0,
          initial=nothing, fixed_u=nothing, fixed_tunnel_opening=nothing,
          preparation_time_limit=min(60.0, 0.8time_limit), operational_margin=0.0)

Optimize generation, unit commitment and eligible tunnel gates with native SCIP.
Without a supplied initial schedule or fixed discrete controls, bounded MILP
proposals and Ipopt dispatch first prepare a replay-valid start. These choices
do not restrict SCIP's search. Preparation receives at most 80% of the total
allowance, capped at 60 seconds; unused time remains available to SCIP.
`preparation_time_limit=0` skips this accelerator. `operational_margin` tightens
only preparation; the global model retains the original physical restrictions.
`diagnostics_path` optionally writes a native SCIP
progress log and an accompanying `.statistics.json` file with native plugin
and LP counters. Returned `scip_diagnostics` and `scip_statistics` are
observational statistics, not independent feasibility or certificate evidence.
The objective
is revenue minus transition and gate movement costs and release penalties, plus
changes in stored water value. `initial` may supply a physical schedule on the same control grid.

`solution` is independently reconstructed and equation-validated. `accepted`
also requires finer chronological replay. `global_certificate` requires an
enclosing solver upper bound within either requested gap. Its scope is the
declared discrete hydraulic model, under numerical solver tolerances. Supplying
`fixed_u` makes the certificate conditional on that commitment matrix;
`fixed_tunnel_opening` makes it conditional on the supplied gate positions.

The time allowance includes preparation, construction and optimization; extraction/replay
can overrun it and are included in `total_seconds`. Julia compilation is also
counted after entry timing begins; compilation needed to enter the public
function can precede that clock. A known exhausted allowance prevents starting
global construction. A validated preparation/initial schedule can still be
returned, with no global bound or gap when search did not run.
"""
solve(c::ScheduleCase;kwargs...) = _solve(c,time,schedule_case,time();kwargs...)

# Clock and preparation are explicit dependencies of this orchestration routine,
# allowing deterministic deadline tests without sleeping or replacing solvers.
function _solve(
    c::ScheduleCase, clock, prepare, began=clock();
    time_limit = 60.0,
    relative_gap = 1e-3,
    absolute_gap = 0.0,
    initial = nothing,
    fixed_u = nothing,
    fixed_tunnel_opening = nothing,
    replay = true,
    diagnostics_path = nothing,
    preparation_time_limit = min(60.0, 0.8time_limit),
    operational_margin = 0.0,
)
    isfinite(time_limit) && time_limit > 0 ||
        throw(ArgumentError("positive finite time_limit required"))
    isfinite(relative_gap) && 0 <= relative_gap < 1 ||
        throw(ArgumentError("relative_gap must lie in [0,1)"))
    isfinite(absolute_gap) && absolute_gap >= 0 ||
        throw(ArgumentError("nonnegative finite absolute_gap required"))
    isfinite(preparation_time_limit) && 0<=preparation_time_limit<=time_limit ||
        throw(ArgumentError("preparation_time_limit must lie in [0,time_limit]"))
    isfinite(operational_margin) && operational_margin>=0 ||
        throw(ArgumentError("nonnegative finite operational_margin required"))
    diagnostics_path!==nothing && (diagnostics_path=abspath(String(diagnostics_path)))
    result = Dict{String,Any}(
        "case" => c.name,
        "solver" => "SCIP",
        "status" => "CONSTRUCTION_BUDGET_EXHAUSTED",
        "accepted" => false,
        "global_certificate" => false,
        "solution" => nothing,
        "feasible_lower_bound" => nothing,
        "global_bound" => nothing,
        "power_cut_statistics" => nothing,
        "relative_gap" => nothing,
        "absolute_gap" => nothing,
        "commitment_fixed" => fixed_u !== nothing,
        "tunnel_opening_fixed" => fixed_tunnel_opening !== nothing,
        "certificate_scope" =>
            fixed_u === nothing && fixed_tunnel_opening===nothing ?
            "declared discrete midpoint model; numerical global bound" :
            "declared discrete midpoint model with "*join(
                vcat(fixed_u===nothing ? String[] : ["fixed commitment"],
                    fixed_tunnel_opening===nothing ? String[] : ["fixed tunnel openings"])," and ")*
                "; numerical global bound",
        "transport_exact" =>
            all(r -> r.deterministic_delay !== nothing, c.system.rivers),
        "requested_relative_gap" => relative_gap,
        "requested_absolute_gap" => absolute_gap,
        "preparation_seconds" => 0.0,
        "preparation_accepted" => false,
        "prepared_objective" => nothing,
        "construction_seconds" => 0.0,
        "solve_seconds" => 0.0,
    )
    preparation_budget=min(preparation_time_limit,max(0.0,time_limit-(clock()-began)))
    if initial===nothing && fixed_u===nothing && fixed_tunnel_opening===nothing &&
       preparation_budget>0
        prepared=prepare(c;time_limit=preparation_budget,
            proposal_time_limit=min(5.0,preparation_budget/4),
            nlp_time_limit=min(20.0,preparation_budget),max_refinements=1,
            refine_grid=false,operational_margin)
        result["preparation_seconds"]=clock()-began
        result["preparation_accepted"]=prepared["accepted"]
        if prepared["accepted"]
            initial=prepared["solution"]
            result["prepared_objective"]=initial.objective
        end
    end
    reconstruction_started=clock()
    best = nothing
    if initial !== nothing
        initial_audit = validate(c, initial)
        initial_audit["valid"] ||
            throw(ArgumentError("initial schedule fails the equation audit"))
        best, correction = _reconstruct_candidate(c, initial)
        best["validation"]["valid"] ||
            throw(ArgumentError("initial controls fail reconstruction audit"))
        result["initial_objective"] = best.objective
    end
    result["initial_reconstruction_seconds"]=clock()-reconstruction_started
    if clock()-began>=time_limit
        # Construction normally validates these inputs. Preserve that contract
        # when the deadline prevents entering it, including fixed-state scope.
        validate_inputs(c)
        fixed_u===nothing || admissible(c,fixed_u) ||
            throw(ArgumentError("invalid fixed commitment"))
        opening=fixed_tunnel_opening===nothing ? nothing :
            _physical_tunnel_openings(c,fixed_tunnel_opening)
        opening===nothing || fixed_u===nothing || validate_tunnel_openings(c,opening,fixed_u)
        if best!==nothing && ((fixed_u!==nothing && best.commitment!=fixed_u) ||
            (opening!==nothing && best.tunnel_opening!=opening))
            best=nothing
            result["initial_rejected_by_fixed_controls"]=true
        end
        phase=result["preparation_seconds"]>=time_limit ? "preparation" : "initial_reconstruction"
        result["budget_exhausted_phase"]=phase
        result["status"]=phase=="preparation" ? "PREPARATION_BUDGET_EXHAUSTED" :
            "INITIALIZATION_BUDGET_EXHAUSTED"
        if best!==nothing
            source=result["preparation_accepted"] ? "preparation" : "initial"
            _set_incumbent!(result,best,source,relative_gap,absolute_gap)
            if replay
                result["replay_audit"]=replay_audit(c,best)
                result["accepted"]=result["replay_audit"]["valid"]
            end
        end
        result["total_seconds"]=clock()-began
        result["budget_overrun_seconds"]=max(0.0,result["total_seconds"]-time_limit)
        return SolveResult(result)
    end
    construction_started=clock()
    b = _build_global_dispatch(c; joint = true, warm = best, fixed_u,tunnel_opening=fixed_tunnel_opening)
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
        "tunnel_gate_binaries"=>count(n->startswith(n,"gate_select_"),binary_names),
        "turbine_cell_binaries"=>count(n->startswith(n,"turbine_") && occursin("_cell[",n),binary_names),
        "river_table_cell_binaries"=>count(n->startswith(n,"river_law_") && occursin("_cell[",n),binary_names),
        "bounds_seconds"=>b.bounds_seconds,
        "shared_head_entries"=>length(b.shared_heads),
        "power_hull_units"=>length(get(b.m.ext,:global_power_hulls,Dict())),
        "table_power_units"=>get(get(b.m.ext,:global_table_power_bounds_profile,Dict()),"units_added",0),
        "head_variables_saved"=>length(c.system.generators)*length(c.prices)-length(b.shared_heads))
    result["variable_count"] = num_variables(b.m)
    result["constraint_count"] =
        num_constraints(b.m; count_variable_in_set_constraints = false)
    result["construction_seconds"] = clock() - construction_started
    remaining = time_limit - (clock()-began)
    remaining>0 || (result["budget_exhausted_phase"]="global_construction")
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
                "propagating/obbt/createbilinineqs" => false,
                "propagating/obbt/itlimitfactor" => 1.0,
            ),
        )
        power_cuts=nothing
        result["solve_seconds"] = @elapsed try
            # Copy resets SCIP's native instance; attach before installing its log.
            JuMP.MOI.Utilities.attach_optimizer(JuMP.backend(b.m))
            if best!==nothing
                result["native_start_audit"]=_complete_scip_start!(b.m)
                if !result["native_start_audit"]["valid"]
                    _clear_scip_start!(b.m)
                end
            end
            if !isempty(get(b.m.ext,:global_table_power_supports,_TablePowerSupportCoordinate[]))
                power_cuts=_install_power_cuts!(b,c)
            end
            remaining_after_copy=max(0.,time_limit-(clock()-began))
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
        if power_cuts!==nothing
            try
                result["power_cut_statistics"]=_power_cut_statistics(power_cuts)
            catch error
                result["power_cut_statistics_error"]=sprint(showerror,error)
            end
        end
        result["scip_diagnostics"]=_scip_diagnostics(b.m)
        if diagnostics_path!==nothing
            try
                statistics_path=diagnostics_path*".statistics.json"
                result["scip_statistics"]=_write_scip_statistics(b.m, statistics_path)
                result["statistics_path"]=statistics_path
            catch error
                result["statistics_error"]=sprint(showerror,error)
            end
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
        _reject_power_cut_bound!(result)
        if result_count(b.m) > 0
            try
                raw = _dispatch_values(b)
                fractional = value.(b.u)
                raw.diagnostics["fractional_commitment"] = fractional
                error = maximum(abs, fractional .- raw.commitment; init = 0.0)
                gate_integrality=maximum((abs(value(v)-round(value(v))) for v in variables if is_binary(v));init=0.0)
                raw.diagnostics["validation"] = validate(c, raw; transport = b.transport)
                result["raw_solver_solution"] = raw
                result["raw_integrality_error"] = error
                result["raw_binary_integrality_error"] = gate_integrality
                sos2_error=maximum((_sos2_residual(JuMP.value.(constraint_object(ref).func),constraint_object(ref).set.weights)
                    for ref in all_constraints(b.m;include_variable_in_set_constraints=false)
                    if constraint_object(ref).set isa MOI.SOS2);init=0.0)
                result["raw_sos2_residual"]=sos2_error
                if error <= 1e-6 && gate_integrality<=1e-6 && sos2_error <= 1e-6 && (fixed_u === nothing || raw.commitment == fixed_u)
                    x, correction = _reconstruct_candidate(c, raw; transport = b.transport)
                    result["control_boundary_correction"] =
                        Dict("flow" => correction.flow, "gate" => correction.gate,
                            "tunnel_opening"=>correction.tunnel_opening)
                    if x["validation"]["valid"] &&
                       (best === nothing || x.objective > best.objective)
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
    result["total_seconds"] = clock() - began
    result["budget_overrun_seconds"] = max(0.0, result["total_seconds"] - time_limit)
    SolveResult(result)
end

export solve
