# Measure schedule generation with --project pointing at the checkout under test.
# Each process measures one implementation; compare matching input/control hashes.
using OpenSHOP, JSON3, SHA, Dates, LinearAlgebra

const HELP = """
Usage: julia --project=CHECKOUT benchmark/schedule.jl --label LABEL --output DIRECTORY \
    --case CASE.json [--seed SEED.json] [--case CASE.json ...]

Options:
  --proposal-seconds N    Per-MILP allowance including construction (default 5).
  --nlp-seconds N         Per-NLP allowance including construction (default 15).
  --refinements N         Additional refined-grid solves (default 0).
  --repeats N             Sequential measured repetitions (default 2).
  --operational-margin N  Candidate power/envelope margin in MW (default 0.1).
  --warmup-case PATH      Otherwise use the first case's first two intervals.
  --warmup-repeats N      Excluded schedule_case calls (default 2).
  --warmup-nlp-seconds N  Per-NLP warmup allowance (default 10).
  --warmup-proposal-seconds N  Per-MILP warmup allowance (default 2).
  --no-warmup            Diagnostic cold-process run; identify it as such.

--seed applies to the immediately preceding --case. Seeds may contain a full
schedule, a result with a solution field, or just u/generator_q/gate controls.
Cold loading, seed reconstruction and warmup are reported separately. Per-solve
allowances are not an end-to-end preparation deadline. The runner never invokes
SCIP or changes solver settings. Full local results accompany scalar summaries.
"""

function options(args)
    cfg = Dict{Symbol,Any}(
        :label => "candidate", :output => joinpath(@__DIR__, "..", "results", "scheduling"),
        :proposal => 5.0, :nlp => 15.0, :refinements => 0, :repeats => 2,
        :margin => 0.1, :warmup_case => nothing, :warmup_repeats => 2,
        :warmup_nlp => 10.0, :warmup_proposal => 2.0, :warmup => true,
    )
    cases = NamedTuple{(:path, :seed),Tuple{String,Union{Nothing,String}}}[]
    floats = Dict("--proposal-seconds" => :proposal, "--nlp-seconds" => :nlp,
        "--operational-margin" => :margin, "--warmup-nlp-seconds" => :warmup_nlp,
        "--warmup-proposal-seconds" => :warmup_proposal)
    integers = Dict("--refinements" => :refinements, "--repeats" => :repeats,
        "--warmup-repeats" => :warmup_repeats)
    strings = Dict("--label" => :label, "--output" => :output,
        "--warmup-case" => :warmup_case)
    i = 1
    while i <= length(args)
        flag = args[i]
        if flag in ("--help", "-h")
            print(HELP)
            return nothing
        elseif flag == "--no-warmup"
            cfg[:warmup] = false
            i += 1
            continue
        end
        i < length(args) || error("Missing value after $flag")
        value = args[i + 1]
        if flag == "--case"
            push!(cases, (path = abspath(value), seed = nothing))
        elseif flag == "--seed"
            isempty(cases) && error("--seed requires a preceding --case")
            cases[end].seed === nothing || error("Duplicate seed for $(cases[end].path)")
            cases[end] = (path = cases[end].path, seed = abspath(value))
        elseif haskey(floats, flag)
            cfg[floats[flag]] = parse(Float64, value)
        elseif haskey(integers, flag)
            cfg[integers[flag]] = parse(Int, value)
        elseif haskey(strings, flag)
            cfg[strings[flag]] = value
        else
            error("Unknown option $flag; use --help")
        end
        i += 2
    end
    isempty(cases) && error("At least one --case is required")
    length(unique(entry.path for entry in cases)) == length(cases) ||
        error("Repeat seeded/unseeded comparisons under separate labels, not duplicate --case paths")
    for key in (:proposal, :nlp, :warmup_nlp, :warmup_proposal)
        isfinite(cfg[key]) && cfg[key] > 0 || error("Positive finite $key required")
    end
    isfinite(cfg[:margin]) && cfg[:margin] >= 0 || error("Invalid operational margin")
    cfg[:refinements] >= 0 || error("Negative refinement count")
    cfg[:repeats] >= 1 && cfg[:warmup_repeats] >= 1 || error("Positive repetitions required")
    occursin(r"^[A-Za-z0-9_-]+$", cfg[:label]) || error("Label must be a simple directory name")
    cfg[:output] = abspath(cfg[:output])
    cfg, cases
end

function canonical(x)
    if x isa AbstractDict
        ordered = sort!(collect(keys(x)); by = string)
        return NamedTuple{Tuple(Symbol(string(k)) for k in ordered)}(
            Tuple(canonical(x[k]) for k in ordered))
    elseif x isa AbstractVector
        return canonical.(x)
    end
    x
end
json_hash(x) = bytes2hex(sha256(JSON3.write(canonical(OpenSHOP.jsonready(x)))))
file_hash(path) = isfile(path) ? bytes2hex(sha256(read(path))) : nothing
# POSIX process clock includes CPU time used by native solver worker threads.
process_cpu_seconds() = Float64(ccall(:clock, Clong, ())) / 1_000_000

function provenance()
    source_root = dirname(dirname(pathof(OpenSHOP)))
    source_files = sort!([joinpath(dir, name) for (dir, _, names) in
        walkdir(joinpath(source_root, "src")) for name in names if endswith(name, ".jl")])
    source = join((relpath(p, source_root) * "\0" * read(p, String) for p in source_files), "\0")
    project = Base.active_project()
    Dict("source_root" => source_root, "source_sha256" => bytes2hex(sha256(source)),
        "project_path" => project, "project_sha256" => project === nothing ? nothing : file_hash(project),
        "manifest_sha256" => project === nothing ? nothing : file_hash(joinpath(dirname(project), "Manifest.toml")),
        "runner_sha256" => file_hash(@__FILE__), "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "blas_configuration" => sprint(show, BLAS.lbt_get_config()),
        "cpu_name" => Sys.CPU_NAME, "kernel" => string(Sys.KERNEL),
        "machine" => Sys.MACHINE,
        "thread_environment" => Dict(k => get(ENV, k, nothing) for k in
            ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_CANCELLATION", "OMP_PROC_BIND")))
end

function matrix_rows(data, columns)
    isempty(data) && return zeros(0, columns)
    all(length(row) == columns for row in data) || error("Seed control columns differ from case grid")
    permutedims(hcat([Float64.(row) for row in data]...))
end

function prepare_seed(c, path)
    path === nothing && return nothing
    data = JSON3.read(read(path, String), Dict{String,Any})
    data = get(data, "solution", data)
    data === nothing && error("Seed has no solution")
    T = length(c.prices)
    u = matrix_rows(data["u"], T)
    all(isinteger, u) || error("Seed commitment is fractional")
    raw = Dict("u" => Int.(u), "generator_q" => matrix_rows(data["generator_q"], T),
        "gate" => matrix_rows(data["gate"], T))
    physical, correction = OpenSHOP._reconstruct_candidate(c, raw)
    physical["validation"]["valid"] || error("Seed controls fail independent equation validation")
    physical["benchmark_control_correction"] = correction
    physical
end

function audit_summary(audit)
    audit === nothing && return Dict{String,Any}()
    keys = ("valid", "physically_valid", "numerically_converged", "seconds", "errors",
        "violations", "replay_intervals", "max_water_residual_Mm3", "max_storage_difference_Mm3",
        "terminal_transit_difference_Mm3", "replayed_objective", "objective_difference",
        "storage_convergence_tolerance_Mm3", "relative_objective_convergence_tolerance")
    Dict(k => audit[k] for k in keys if haskey(audit, k))
end

function attempts_summary(result)
    map(get(result, "attempts", Any[])) do attempt
        row = Dict{String,Any}(string(k) => v for (k, v) in attempt if k != "steps")
        haskey(attempt, "u") && (row["commitment_sha256"] = json_hash(attempt["u"]))
        row["steps"] = map(get(attempt, "steps", Any[])) do step
            s = Dict{String,Any}(string(k) => v for (k, v) in step if k != "replay")
            s["replay"] = audit_summary(get(step, "replay", Dict()))
            s
        end
        row
    end
end

function acceptance_events(result, total, initial)
    attempts = get(result, "attempts", Any[])
    proposal_seconds = sum(get(p, "total_seconds", get(p, "seconds", 0.0)) for
        p in get(result, "proposals", Any[]); init = 0.0)
    measured_stages = proposal_seconds + sum(get(a, "seconds", 0.0) for a in attempts; init = 0.0)
    overhead = max(0.0, total - measured_stages)
    prefix = 0.0
    proposals_counted = false
    events = Any[]
    for (index, attempt) in enumerate(attempts)
        stage = get(attempt, "stage", "")
        supplied = stage in ("initial", "supplied_initial", "seed")
        if stage != "feasibility" && !supplied && !proposals_counted
            prefix += proposal_seconds
            proposals_counted = true
        end
        prefix += get(attempt, "seconds", 0.0)
        supplied && continue
        get(attempt, "accepted", false) || continue
        steps = get(attempt, "steps", Any[])
        objective = isempty(steps) ? nothing : get(last(steps), "objective", nothing)
        exact = get(attempt, "elapsed_seconds", nothing)
        push!(events, Dict("attempt" => index, "stage" => get(attempt, "stage", nothing),
            "objective" => objective, "elapsed_seconds" => exact,
            "elapsed_interval_seconds" => exact === nothing ? [prefix, prefix + overhead] : [exact, exact],
            "timing_scope" => exact === nothing ? "inferred from ordered stages plus unassigned overhead" :
                "recorded since schedule_case start"))
    end
    first_new = isempty(events) ? nothing : first(events)
    # A supplied audited seed is already available before the measured call.
    # Its reacceptance time is not present in baseline attempts and is not invented.
    first_available = initial === nothing ? first_new : Dict(
        "objective" => initial["objective"], "elapsed_seconds" => 0.0,
        "timing_scope" => "supplied independently validated seed; reconstruction excluded")
    best = initial === nothing ? -Inf : initial["objective"]
    improvements = Any[]
    for event in events
        objective = event["objective"]
        objective === nothing && continue
        threshold = isfinite(best) ? max(1e-6, 1e-10 * abs(best)) : 0.0
        if objective > best + threshold
            push!(improvements, event)
            best = objective
        end
    end
    (events = events, first_new = first_new, first_available = first_available,
        improvements = improvements,
        proposal_seconds = proposal_seconds, unassigned_seconds = overhead,
        stage_accounting_exceeds_total = measured_stages > total + 1e-6)
end

function proposal_summary(result)
    map(enumerate(get(result, "proposals", Any[]))) do pair
        index, p = pair
        row = Dict{String,Any}("index" => index)
        for k in ("status", "seconds", "construction_seconds", "total_seconds", "water_scale", "proposal_objective")
            haskey(p, k) && (row[k] = p[k])
        end
        if haskey(p, "u")
            row["u"] = p["u"]
            row["commitment_sha256"] = json_hash(p["u"])
        end
        row["objective_scope"] = "MILP surrogate; not an accepted physical objective or global bound"
        row
    end
end

function short_warmup_case(c)
    n = min(2, length(c.prices))
    grid = copy(c.grid[1:(n + 1)])
    # with_grid only refines an unchanged horizon. Build a separate warmup case
    # and shorten explicit observation windows as well as the decision horizon.
    rivers = map(c.system.rivers) do r
        isempty(r.arrival_window_grid) && return r
        windows = sort!(unique!(vcat(first(grid),
            [x for x in r.arrival_window_grid if first(grid) < x < last(grid)], last(grid))))
        OpenSHOP._river_replace(r; arrival_window_grid = windows)
    end
    ScheduleCase(; name = c.name * "_warmup", system = OpenSHOP._river_replace(c.system; rivers),
        grid, prices = copy(c.prices[1:n]), operations = c.operations, flow_requirements = c.flow_requirements)
end

function run_benchmark(cfg, cases)
    root = joinpath(cfg[:output], cfg[:label])
    mkpath(root)
    metadata = provenance()
    OpenSHOP.writejson(joinpath(root, "configuration.json"),
        Dict("configuration" => cfg, "provenance" => metadata, "cases" => cases,
            "started_utc" => string(now(UTC)), "scope" => "feasible schedule preparation only; no global proof"))
    warmups = Any[]
    if cfg[:warmup]
        warm_path = cfg[:warmup_case] === nothing ? first(cases).path : abspath(cfg[:warmup_case])
        c = readcase(warm_path)
        warm = cfg[:warmup_case] === nothing ? short_warmup_case(c) : c
        for repetition in 1:cfg[:warmup_repeats]
            started = time()
            x = schedule_case(warm; proposal_time_limit = cfg[:warmup_proposal],
                nlp_time_limit = cfg[:warmup_nlp], max_refinements = 0, operational_margin = cfg[:margin])
            push!(warmups, Dict("repeat" => repetition, "input_sha256" => file_hash(warm_path),
                "case_sha256" => json_hash(case_dict(warm)), "intervals" => length(warm.prices),
                "elapsed_seconds_excluded" => time() - started, "accepted" => x["accepted"],
                "attempts" => attempts_summary(x), "proposals" => proposal_summary(x)))
            OpenSHOP.writejson(joinpath(root, "warmup.json"), warmups)
        end
    end
    rows = Any[]
    for entry in cases
        loading_started = time()
        c = readcase(entry.path)
        loading_seconds = time() - loading_started
        input_hash = file_hash(entry.path)
        case_hash = json_hash(case_dict(c))
        slug = replace(c.name, r"[^A-Za-z0-9_-]" => "_") * "-" * first(case_hash, 12)
        folder = joinpath(root, slug)
        mkpath(folder)
        write(joinpath(folder, "input.json"), read(entry.path))
        seed_started = time()
        common_initial = prepare_seed(c, entry.seed)
        seed_seconds = time() - seed_started
        controls = common_initial === nothing ? nothing : Dict(k => common_initial[k] for k in ("u", "generator_q", "gate"))
        OpenSHOP.writejson(joinpath(folder, "seed-controls.json"), controls)
        for repetition in 1:cfg[:repeats]
            initial = common_initial === nothing ? nothing : deepcopy(common_initial)
            row = merge(copy(metadata), Dict{String,Any}(
                "label" => cfg[:label], "case" => c.name, "case_id" => slug,
                "input_path" => entry.path, "input_sha256" => input_hash,
                "case_sha256" => case_hash, "seed_input_sha256" => entry.seed === nothing ? nothing : file_hash(entry.seed),
                "seed_controls_sha256" => json_hash(controls), "repeat" => repetition,
                "warmup_enabled" => cfg[:warmup], "cold_loading_seconds_excluded" => loading_seconds,
                "seed_preparation_seconds_excluded" => seed_seconds,
                "seed_objective" => initial === nothing ? nothing : initial["objective"],
                "proposal_allowance_seconds" => cfg[:proposal], "nlp_allowance_seconds" => cfg[:nlp],
                "max_refinements" => cfg[:refinements], "operational_margin_MW" => cfg[:margin]))
            started = time()
            cpu_started = process_cpu_seconds()
            try
                measurement = @timed schedule_case(c; initial, proposal_time_limit = cfg[:proposal],
                    nlp_time_limit = cfg[:nlp], max_refinements = cfg[:refinements], operational_margin = cfg[:margin])
                result = measurement.value
                elapsed = time() - started
                total = max(elapsed, get(result, "seconds", elapsed))
                acceptance = acceptance_events(result, total, initial)
                solution = get(result, "solution", nothing)
                accepted = get(result, "accepted", false)
                merge!(row, Dict("status" => accepted ? "ACCEPTED" : "NO_ACCEPTED_SCHEDULE",
                    "accepted" => accepted, "accepted_objective" => accepted ? solution["objective"] : nothing,
                    "harness_seconds" => elapsed, "schedule_seconds" => get(result, "seconds", nothing),
                    "process_cpu_seconds" => process_cpu_seconds() - cpu_started,
                    "allocated_bytes" => measurement.bytes, "gc_seconds" => measurement.gctime,
                    "allocation_scope" => "schedule_case only; summary and result serialization excluded",
                    "compilation_seconds" => hasproperty(measurement, :compile_time) ? measurement.compile_time : nothing,
                    "recompilation_seconds" => hasproperty(measurement, :recompile_time) ? measurement.recompile_time : nothing,
                    "proposal_total_seconds" => acceptance.proposal_seconds,
                    "unassigned_overhead_seconds" => acceptance.unassigned_seconds,
                    "stage_accounting_exceeds_total" => acceptance.stage_accounting_exceeds_total,
                    "first_available_schedule" => acceptance.first_available,
                    "first_new_accepted_schedule" => acceptance.first_new,
                    "first_improving_schedule" => isempty(acceptance.improvements) ? nothing : first(acceptance.improvements),
                    "incumbent_improvement_events" => acceptance.improvements,
                    "acceptance_events" => acceptance.events, "attempts" => attempts_summary(result),
                    "dispatch_attempt_count" => count(a -> !(get(a, "stage", "") in
                        ("feasibility", "initial", "supplied_initial", "seed")), get(result, "attempts", Any[])),
                    "proposals" => proposal_summary(result), "audit" => audit_summary(get(result, "audit", Dict())),
                    "selected_commitment_sha256" => solution === nothing ? nothing : json_hash(solution["u"]),
                    "selected_u" => solution === nothing ? nothing : solution["u"],
                    "solution_grid_sha256" => get(result, "case", nothing) === nothing ? nothing : json_hash(result["case"]["grid"])))
                OpenSHOP.writejson(joinpath(folder, "repeat-$(lpad(string(repetition), 2, '0')).json"), result)
            catch error
                row["status"] = "BENCHMARK_ERROR"
                row["accepted"] = false
                row["error"] = sprint(showerror, error)
                row["harness_seconds"] = time() - started
                row["process_cpu_seconds"] = process_cpu_seconds() - cpu_started
            end
            push!(rows, row)
            OpenSHOP.writejson(joinpath(root, "summary.json"), rows)
            println(now(UTC), " ", cfg[:label], " ", slug, " repeat ", repetition,
                " ", row["status"], " ", round(row["harness_seconds"]; digits = 3), " s")
            flush(stdout)
        end
    end
    rows
end

parsed = options(ARGS)
if parsed !== nothing
    rows=run_benchmark(parsed...)
    all(row->row["accepted"], rows) ||
        error("benchmark lacks an accepted schedule; inspect the recorded attempts and audits")
end
