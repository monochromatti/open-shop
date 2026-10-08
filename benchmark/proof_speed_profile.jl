# Experimental proof policies. Cases and audited controls are frozen per job.
# Uploaded summaries contain scalar measurements, not model inputs or controls.
include("diagnose.jl")
using SCIP, JuMP, TOML
isfile(joinpath(@__DIR__, "targeted_power_supports.jl")) && include("targeted_power_supports.jl")

const PROOF_PROFILES = Dict{String,NamedTuple{(:policy, :parameters),Tuple{Symbol,Vector{Pair{String,Any}}}}}(
    "baseline" => (policy=:baseline, parameters=Pair{String,Any}[]),
    "no_bilin" => (policy=:baseline, parameters=Pair{String,Any}["propagating/obbt/createbilinineqs"=>false]),
    "filter" => (policy=:baseline, parameters=Pair{String,Any}["propagating/obbt/applyfilterrounds"=>true]),
    "no_bilin_filter" => (policy=:baseline, parameters=Pair{String,Any}[
        "propagating/obbt/createbilinineqs"=>false, "propagating/obbt/applyfilterrounds"=>true]),
    "head_shared" => (policy=:head_shared, parameters=Pair{String,Any}[]),
    "discharge_affine" => (policy=:discharge_affine, parameters=Pair{String,Any}[]),
    "shared" => (policy=:shared, parameters=Pair{String,Any}[]),
    "shared_cuts" => (policy=:shared, parameters=Pair{String,Any}[]),
    "shared_cuts_no_bilin" => (policy=:shared, parameters=Pair{String,Any}["propagating/obbt/createbilinineqs"=>false]),
    "shared_cuts_wide" => (policy=:shared, parameters=Pair{String,Any}[]),
    "shared_cuts_wide_no_bilin" => (policy=:shared, parameters=Pair{String,Any}["propagating/obbt/createbilinineqs"=>false]),
    "shared_no_bilin" => (policy=:shared, parameters=Pair{String,Any}["propagating/obbt/createbilinineqs"=>false]),
    "shared_filter" => (policy=:shared, parameters=Pair{String,Any}["propagating/obbt/applyfilterrounds"=>true]),
    "shared_no_bilin_filter" => (policy=:shared, parameters=Pair{String,Any}[
        "propagating/obbt/createbilinineqs"=>false, "propagating/obbt/applyfilterrounds"=>true]))

const TARGETED_PROFILE_OPTIONS = Dict(
    "shared_cuts_wide" => (max_cuts=128,max_rounds=8,max_cuts_per_round=32,
        max_coordinates_per_round=32,max_coordinate_checks=256),
    "shared_cuts_wide_no_bilin" => (max_cuts=128,max_rounds=8,max_cuts_per_round=32,
        max_coordinates_per_round=32,max_coordinate_checks=256))

const PROOF_GAP_THRESHOLDS = (0.10, 0.08, 0.07, 0.06, 0.05, 0.01, 0.0001)

mutable struct ProofBounds <: SCIP.AbstractEventhdlr
    optimizer::SCIP.Optimizer
    started_ns::UInt64
    audited_lower::Float64
    points::Vector{NamedTuple{(:seconds, :upper, :native_seconds, :run),Tuple{Float64,Float64,Float64,Int}}}
    calls::Int
    unavailable::Int
    monotonic::Bool
    maximum_increase::Float64
    capture_seconds::Float64
    errors::Vector{String}
end

function SCIP.eventinit(e::ProofBounds)
    SCIP.catch_event(e.optimizer.inner, SCIP.SCIP_EVENTTYPE_DUALBOUNDIMPROVED, e)
end

function SCIP.eventexit(e::ProofBounds)
    SCIP.drop_event(e.optimizer.inner, SCIP.SCIP_EVENTTYPE_DUALBOUNDIMPROVED, e)
end

function SCIP.eventexec(e::ProofBounds)
    began_ns = time_ns()
    e.calls += 1
    try
        # This is SCIP's global bound in original objective units. An LP value
        # observed during probing, diving or root diagnostics is never used.
        upper = OpenSHOP._scip_bound_value(e.optimizer, SCIP.SCIPgetDualbound(e.optimizer))
        if upper === nothing
            e.unavailable += 1
            return
        end
        upper >= e.audited_lower - 1e-6 || error("event upper bound below frozen audited schedule")
        if !isempty(e.points)
            increase = upper - last(e.points).upper
            tolerance = 1e-8 * max(1.0, abs(upper), abs(last(e.points).upper))
            if increase > tolerance
                e.monotonic = false
                e.maximum_increase = max(e.maximum_increase, increase)
            end
        end
        push!(e.points, (seconds=Float64(time_ns()-e.started_ns)/1e9,
            upper=Float64(upper), native_seconds=Float64(SCIP.SCIPgetSolvingTime(e.optimizer)),
            run=Int(SCIP.SCIPgetNRuns(e.optimizer))))
    catch exception
        push!(e.errors, sprint(showerror, exception))
    finally
        e.capture_seconds += Float64(time_ns()-began_ns)/1e9
    end
    nothing
end

proof_gap(upper, lower) = max(0.0, upper-lower) / max(1.0, abs(lower))

"""First observed threshold attainment; an unreached threshold is right censored."""
function proof_thresholds(points, lower, observed_seconds; capture_enabled=true)
    [begin
        index = findfirst(p -> proof_gap(p.upper, lower) <= target, points)
        Dict("relative_gap"=>target, "gap_percent"=>100target,
            "attained"=>index !== nothing,
            "first_observed_seconds"=>index === nothing ? nothing : points[index].seconds,
            "right_censored"=>index === nothing,
            "observation_end_seconds"=>observed_seconds,
            "trajectory_capture_enabled"=>capture_enabled,
            "timing_scope"=>capture_enabled ? "native global-bound events, including construction; final accepted bound added at native completion" :
                "final accepted bound only; threshold time is an upper bound on attainment")
    end for target in PROOF_GAP_THRESHOLDS]
end

function proof_runtime()
    # The checked-in Manifest fixes package versions; its hash identifies the
    # complete dependency graph without relying on the compatible Project range.
    manifest = joinpath(@__DIR__, "..", "Manifest.toml")
    dependencies = TOML.parsefile(manifest)["deps"]
    packages = Dict(package=>only(dependencies[package])["version"] for package in
        ("JuMP", "MathOptInterface", "SCIP", "SCIP_jll", "SCIP_PaPILO_jll", "SoPlex_jll", "HiGHS", "Ipopt")
        if haskey(dependencies, package))
    Dict("manifest_sha256"=>bytes2hex(sha256(read(manifest))),
        "locked_packages"=>packages,
        "native_scip_version"=>join((SCIP.SCIPmajorVersion(), SCIP.SCIPminorVersion(), SCIP.SCIPtechVersion()), '.'),
        "native_lp_solver"=>unsafe_string(SCIP.SCIPlpiGetSolverName()),
        "julia_version"=>string(VERSION), "threads"=>Threads.nthreads(),
        "blas_threads"=>OpenSHOP.LinearAlgebra.BLAS.get_num_threads(),
        "cpu_name"=>Sys.CPU_NAME, "kernel"=>string(Sys.KERNEL),
        "source_sha256"=>bytes2hex(sha256(join(read(p, String) for p in
            sort(filter(p->endswith(p, ".jl"), readdir(joinpath(@__DIR__, "..", "src"); join=true)))))),
        "experiment_sha256"=>bytes2hex(sha256(join(read(joinpath(@__DIR__, p), String) for p in
            ("proof_speed_profile.jl", "diagnose.jl", "targeted_power_supports.jl") if isfile(joinpath(@__DIR__, p))))))
end

function proof_profile(input, output; seconds=120.0, repeats=1,
        profiles=["baseline", "no_bilin", "filter", "no_bilin_filter", "head_shared", "discharge_affine", "shared"],
        capture=true, commitment="free")
    isfinite(seconds) && seconds > 0 || error("positive finite allowance required")
    repeats isa Integer && repeats >= 1 || error("positive integer repeats required")
    commitment in ("free", "fixed") || error("unknown commitment scope")
    !isempty(profiles) && all(p->haskey(PROOF_PROFILES,p), profiles) || error("unknown experimental profile")
    length(unique(profiles)) == length(profiles) || error("duplicate experimental profile")
    mkpath(output)
    frozen_case = joinpath(output, "case.json")
    benchmark_freeze(frozen_case, case_dict(readcase(input)))
    c = readcase(frozen_case)
    seedfile = joinpath(output, "seed.json")
    preparation_seconds = 0.0
    if !isfile(seedfile)
        preparation_seconds = @elapsed preparation = schedule_case(c;
            proposal_time_limit=5.0, nlp_time_limit=20.0, max_refinements=0, operational_margin=0.1)
        preparation["accepted"] || error("no audited common seed")
        seed = preparation["solution"]
        q = copy(seed["generator_q"])
        q[seed["u"].==0] .= 0.0
        seed = dispatch_from_controls(c, seed["u"], q, seed["gate"])
        seed["validation"]["valid"] || error("reconstructed seed invalid")
        benchmark_freeze(seedfile, seed)
    end
    initial = benchmark_seed(seedfile, c)
    initial !== nothing && initial["validation"]["valid"] || error("frozen seed fails physical audit")
    # Include only physical controls in this hash; model-specific lift variables
    # cannot make two nominally identical runs use different schedules.
    controlsfile = joinpath(output, "seed-controls.json")
    benchmark_freeze(controlsfile, Dict(k=>initial[k] for k in ("u", "generator_q", "gate")))
    audited_lower = Float64(initial["objective"])
    fixed_u = commitment == "fixed" ? copy(initial["u"]) : nothing
    metadata = merge(proof_runtime(), Dict("case"=>c.name,
        "case_sha256"=>bytes2hex(sha256(read(frozen_case))),
        "seed_sha256"=>bytes2hex(sha256(read(seedfile))),
        "seed_controls_sha256"=>bytes2hex(sha256(read(controlsfile))),
        "seed_objective"=>audited_lower, "threshold_audited_lower_bound"=>audited_lower,
        "preparation_seconds_excluded"=>preparation_seconds, "commitment"=>commitment,
        "capture_enabled"=>capture, "allowance_seconds"=>seconds,
        "scope"=>commitment == "free" ? "whole declared discrete model" : "declared discrete model with frozen commitment",
        "threshold_comparison"=>"all profiles use the same frozen audited lower bound; no rounded display or local LP bounds",
        "schema_version"=>1))
    # Warm common compilation before policy-specific compilation. Neither warmup
    # produces a replacement seed or contributes to a measured time.
    common_warm_seconds = @elapsed begin
        warm = OpenSHOP._solve(c; initial, fixed_u, time_limit=15.0,
            diagnostics_path=joinpath(output, "warm-common.log"))
        if warm["status"] == "CONSTRUCTION_BUDGET_EXHAUSTED"
            OpenSHOP._solve(c; initial, fixed_u, time_limit=10.0,
                diagnostics_path=joinpath(output, "warm-common-retry.log"))
        end
    end
    metadata["common_warmup_elapsed_seconds_excluded"] = common_warm_seconds
    rows = Any[]
    warm_seconds = Dict{String,Float64}()
    for repeat in 1:repeats, profile in (isodd(repeat) ? profiles : reverse(profiles))
        config = PROOF_PROFILES[profile]
        row = merge(copy(metadata), Dict("profile"=>profile, "repeat"=>repeat,
            "table_power_policy"=>string(config.policy), "parameter_changes"=>Dict(config.parameters)))
        try
            graph = Ref{Any}(nothing)
            event = Ref{Union{Nothing,ProofBounds}}(nothing)
            separator = Ref{Any}(nothing)
            started_ns = Ref(UInt64(0))
            setup = function(b)
                graph[] = b
                for (key, val) in config.parameters
                    set_optimizer_attribute(b.m, key, val)
                end
                if startswith(profile, "shared_cuts")
                    isdefined(@__MODULE__, :install_targeted_supports) || error("targeted supports installer unavailable")
                    separator[] = install_targeted_supports(b, c; get(TARGETED_PROFILE_OPTIONS,profile,(;))...)
                    Base.precompile(SCIP.exec_lp, (typeof(separator[]),))
                end
                if capture
                    optimizer = unsafe_backend(b.m)
                    observer = ProofBounds(optimizer, started_ns[], audited_lower,
                        NamedTuple{(:seconds,:upper,:native_seconds,:run),Tuple{Float64,Float64,Float64,Int}}[],
                        0, 0, true, 0.0, 0.0, String[])
                    Base.precompile(SCIP.eventexec, (typeof(observer),))
                    SCIP.include_event_handler(optimizer.inner, observer; desc="Global proof-bound timing")
                    event[] = observer
                end
            end
            if repeat == 1
                warm_seconds[profile] = @elapsed begin
                    started_ns[] = time_ns()
                    warmed = OpenSHOP._solve(c; initial, fixed_u, time_limit=10.0,
                        table_power_policy=config.policy, optimizer_setup=setup,
                        diagnostics_path=joinpath(output, "warm-$(profile).log"))
                    if warmed["status"] == "CONSTRUCTION_BUDGET_EXHAUSTED"
                        started_ns[] = time_ns()
                        OpenSHOP._solve(c; initial, fixed_u, time_limit=10.0,
                            table_power_policy=config.policy, optimizer_setup=setup,
                            diagnostics_path=joinpath(output, "warm-$(profile)-retry.log"))
                    end
                end
            end
            graph[] = nothing
            event[] = nothing
            separator[] = nothing
            logpath = joinpath(output, "$(profile)-$(commitment)-$(repeat).log")
            started_ns[] = time_ns()
            result = OpenSHOP._solve(c; initial, fixed_u, time_limit=seconds, relative_gap=1e-4,
                table_power_policy=config.policy, optimizer_setup=setup, diagnostics_path=logpath)
            row["harness_seconds"] = Float64(time_ns()-started_ns[])/1e9
            for key in ("status", "accepted", "global_certificate", "feasible_lower_bound", "global_bound",
                    "relative_gap", "absolute_gap", "construction_seconds", "solve_seconds", "total_seconds",
                    "variable_count", "constraint_count", "model_profile", "start_audit", "incumbent_source",
                    "scip_diagnostics", "scip_statistics", "solver_error", "candidate_error", "statistics_error",
                    "diagnostics_close_error", "bound_rejection", "budget_overrun_seconds", "requested_relative_gap",
                    "certificate_scope", "raw_integrality_error", "raw_sos2_residual")
                row[key] = get(result, key, nothing)
            end
            row["profile_warmup_elapsed_seconds_excluded"] = get(warm_seconds, profile, 0.0)
            row["profile_warmup_nominal_allowance_seconds_per_attempt_excluded"] = repeat == 1 ? 10.0 : 0.0
            row["progress"] = scip_progress(logpath)
            points = event[] === nothing ?
                NamedTuple{(:seconds,:upper,:native_seconds,:run),Tuple{Float64,Float64,Float64,Int}}[] : copy(event[].points)
            proof_end = max(get(result, "construction_seconds", 0.0) + get(result, "solve_seconds", 0.0),
                maximum(p->p.seconds, points; init=0.0))
            upper = get(result, "global_bound", nothing)
            if upper !== nothing
                upper >= audited_lower - 1e-6 || error("final upper bound below frozen audited schedule")
                # Presolve can finish without an improved-bound event. Its final
                # accepted bound still supplies a conservative attainment time.
                push!(points, (seconds=Float64(proof_end), upper=Float64(upper),
                    native_seconds=Float64(get(get(result, "scip_diagnostics", Dict()), "native_solve_seconds", 0.0)),
                    run=event[] === nothing || isempty(event[].points) ? 0 : last(event[].points).run))
            end
            sort!(points; by=p->p.seconds)
            row["global_bound_trajectory"] = points
            row["proof_observation_end_seconds"] = proof_end
            row["gap_against_frozen_audited_objective"] = upper === nothing ? nothing : proof_gap(upper, audited_lower)
            row["time_to_gap"] = proof_thresholds(points, audited_lower, proof_end; capture_enabled=capture)
            row["bound_event_calls"] = event[] === nothing ? 0 : event[].calls
            row["bound_event_unavailable"] = event[] === nothing ? 0 : event[].unavailable
            row["bound_event_capture_seconds"] = event[] === nothing ? 0.0 : event[].capture_seconds
            row["bound_event_capture_fraction"] = row["bound_event_capture_seconds"] / max(proof_end, eps())
            row["bound_event_monotonic"] = event[] === nothing ? nothing : event[].monotonic
            row["bound_event_maximum_increase"] = event[] === nothing ? 0.0 : event[].maximum_increase
            row["bound_event_errors"] = event[] === nothing ? String[] : event[].errors
            row["table_power_profile"] = graph[] === nothing ? nothing :
                get(graph[].m.ext, :global_table_power_bounds_profile, nothing)
            row["targeted_supports"] = separator[] === nothing ? nothing : targeted_supports_statistics(separator[])
            get(result, "accepted", false) && get(get(result, "start_audit", Dict()), "valid", false) || error("unaccepted schedule/start")
            upper !== nothing || error("no usable accepted global upper bound")
            get(result, "scip_statistics", nothing) !== nothing || error("native plugin statistics unavailable")
            for key in ("solver_error", "statistics_error", "diagnostics_close_error", "bound_rejection")
                get(result, key, nothing) === nothing || error("$(key): $(result[key])")
            end
            isempty(row["bound_event_errors"]) || error("native bound observer failed")
            sepstats=row["targeted_supports"]
            sepstats===nothing || isempty(sepstats["errors"]) || error("targeted support separator failed")
        catch exception
            row["experiment_error"] = sprint(showerror, exception)
            haskey(row, "status") || (row["status"] = "EXPERIMENT_ERROR")
            haskey(row, "accepted") || (row["accepted"] = false)
        end
        push!(rows, row)
        writejson(joinpath(output, "summary.json"), rows)
        println(now(), " ", profile, " ", commitment, " ", row["status"],
            " gap=", get(row, "gap_against_frozen_audited_objective", nothing),
            haskey(row, "experiment_error") ? " error=" * row["experiment_error"] : "")
        flush(stdout)
    end
    # Every accepted schedule is feasible in the same case. A bound that excludes
    # any one of them cannot be used even when it encloses its own incumbent.
    accepted = [r for r in rows if get(r, "accepted", false) && get(r, "feasible_lower_bound", nothing) !== nothing]
    if !isempty(accepted)
        pooled_lower = maximum(r["feasible_lower_bound"] for r in accepted)
        for row in rows
            row["matched_job_best_accepted_objective"] = pooled_lower
            points = get(row, "global_bound_trajectory", [])
            if any(p->p.upper < pooled_lower-1e-6, points)
                row["experiment_error"] = "native bound trajectory excludes an independently audited matched schedule"
                row["rejected_global_bound"] = get(row, "global_bound", nothing)
                row["global_bound"] = nothing
                row["global_certificate"] = false
                row["relative_gap"] = nothing
                row["time_to_gap"] = nothing
            else
                row["time_to_gap_against_matched_job_best_objective"] = proof_thresholds(points,
                    pooled_lower, get(row, "proof_observation_end_seconds", 0.0); capture_enabled=capture)
            end
        end
        writejson(joinpath(output, "summary.json"), rows)
    end
    any(r->haskey(r, "experiment_error"), rows) && error("one or more proof experiments failed; records retained")
    rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    proof_profile(abspath(ARGS[1]), abspath(ARGS[2]); seconds=parse(Float64, ARGS[3]),
        repeats=parse(Int, ARGS[4]), profiles=split(ARGS[5], ','),
        capture=length(ARGS)<6 || ARGS[6]=="true", commitment=length(ARGS)<7 ? "free" : ARGS[7])
end
