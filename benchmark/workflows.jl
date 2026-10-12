# Compare complete workflows from identical inputs, without a supplied schedule.
include("diagnose.jl")
include("tunnel-gates.jl")
using HiGHS, SCIP

function workflow_trial(c, method; budget, output, preparation_fraction=0.4,
    highs_single_thread=true)
    started=time()
    prepared=nothing
    first_validated=nothing
    initial=nothing
    if method!=:scip
        factory=method==:highs_ipopt ?
            (highs_single_thread ? optimizer_with_attributes(HiGHS.Optimizer,"threads"=>1) : HiGHS.Optimizer) :
            optimizer_with_attributes(SCIP.Optimizer,"parallel/maxnthreads"=>1,
                "numerics/feastol"=>1e-8)
        preparation_budget=min(60.0,preparation_fraction*budget)
        prepared=schedule_case(c;time_limit=preparation_budget,
            proposal_time_limit=min(5.0,preparation_budget/4),
            nlp_time_limit=min(20.0,preparation_budget),max_refinements=1,
            refine_grid=false,operational_margin=0.0,proposal_optimizer=factory)
        if prepared["accepted"]
            initial=prepared["solution"]
            first_validated=minimum(a["elapsed_seconds"] for a in prepared["attempts"]
                if get(a,"accepted",false))
        end
        writejson(joinpath(output,"preparation.json"),prepared)
    end
    preparation_seconds=time()-started
    remaining=budget-preparation_seconds
    remaining>0 || error("preparation consumed the entire workflow allowance")
    result=solve(c;initial,time_limit=remaining,preparation_time_limit=0,
        diagnostics_path=joinpath(output,"scip.log"))
    elapsed=time()-started
    writejson(joinpath(output,"result.json"),result)
    result["accepted"] && first_validated===nothing && (first_validated=elapsed)
    x=result["solution"]
    first_native=get(get(get(result,"scip_statistics",Dict()),"solution",Dict()),
        "first_solution_time",nothing)
    first_native!==nothing && (!isfinite(first_native) || first_native>=1e19) &&
        (first_native=nothing)
    return Dict(
        "case"=>c.name,"method"=>string(method),"budget_seconds"=>budget,
        "preparation_fraction"=>method==:scip ? 0.0 : preparation_fraction,
        "proposal_thread_policy"=>method==:highs_ipopt && !highs_single_thread ? "HiGHS default" : "one thread",
        "preparation_seconds"=>preparation_seconds,
        "preparation_accepted"=>prepared===nothing ? false : prepared["accepted"],
        "global_start_valid"=>get(get(result,"start_audit",Dict()),"valid",nothing),
        "prepared_objective"=>initial===nothing ? nothing : initial["objective"],
        "first_validated_observed_seconds"=>first_validated,
        "global_allowance_seconds"=>remaining,
        "total_seconds"=>elapsed,"budget_overrun_seconds"=>max(0.0,elapsed-budget),
        "accepted"=>result["accepted"],"status"=>result["status"],
        "objective"=>result["accepted"] ? result["feasible_lower_bound"] : nothing,
        "upper_bound"=>result["global_bound"],"relative_gap"=>result["relative_gap"],
        "global_certificate"=>result["global_certificate"],
        "incumbent_source"=>get(result,"incumbent_source",nothing),
        "variables"=>result["variable_count"],"constraints"=>result["constraint_count"],
        "native_statistics"=>result["scip_diagnostics"],
        "native_first_solution_seconds"=>first_native,
        "native_nlp_solvers"=>get(get(get(result,"scip_statistics",Dict()),
            "nlpi",Dict()),"nlp_solvers",nothing),
        "solver_error"=>get(result,"solver_error",nothing),
        "candidate_error"=>get(result,"candidate_error",nothing),
        "replay_errors"=>get(get(result,"replay_audit",Dict()),"errors",[]),
        "replay_storage_difference_Mm3"=>get(get(result,"replay_audit",Dict()),
            "max_storage_difference_Mm3",nothing),
        "on_intervals"=>x===nothing ? nothing : sum(x["u"]),
        "preparation_attempts"=>prepared===nothing ? [] : [
            merge(Dict(k=>get(a,k,nothing) for k in ("stage","accepted","seconds","elapsed_seconds")),
                Dict("objective"=>get(a,"accepted",false) && !isempty(get(a,"steps",[])) ?
                    last(a["steps"])["objective"] : nothing))
            for a in prepared["attempts"]],
        "proposals"=>prepared===nothing ? [] : [
            Dict(k=>get(p,k,nothing) for k in ("water_scale","status","total_seconds"))
            for p in prepared["proposals"]])
end

function workflow_quality!(records)
    for label in unique(r["label"] for r in records)
        group=filter(r->r["label"]==label,records)
        values=[r["objective"] for r in group if r["objective"]!==nothing]
        if isempty(values)
            for r in group
                r["best_observed_objective"]=nothing
                r["within_0_1_percent_of_best_observed_seconds"]=nothing
            end
            continue
        end
        best=maximum(values)
        threshold=best-0.001max(1.0,abs(best))
        for r in group
            times=[a["elapsed_seconds"] for a in r["preparation_attempts"]
                if get(a,"objective",nothing)!==nothing && a["objective"]>=threshold]
            r["objective"]!==nothing && r["objective"]>=threshold && push!(times,r["total_seconds"])
            r["best_observed_objective"]=best
            r["within_0_1_percent_of_best_observed_seconds"]=isempty(times) ? nothing : minimum(times)
        end
    end
    records
end

function workflow_cases(labels)
    root=get(ENV,"OPENSHOP_TOKKE_CASES",joinpath(@__DIR__,"..","examples","tokke_vinje","generated"))
    map(labels) do label
        c=label=="distributed" ? readcase(joinpath(@__DIR__,"cases","distributed-rivers.json")) :
          label=="tables" ? readcase(joinpath(@__DIR__,"cases","turbine-tables.json")) :
          startswith(label,"gates") ? gate_network_case(parse(Int,label[6:end])) :
          startswith(label,"tokke") ? readcase(joinpath(root,"tokke_vinje_$(label[6:end])h.json")) :
          error("unknown case: $label")
        (label,c)
    end
end

function benchmark_workflows(output,labels;budget=30.0,repeats=1,
    methods=(:scip,:highs_ipopt,:scip_ipopt),preparation_fraction=0.4)
    isfinite(budget) && budget>0 || throw(ArgumentError("positive finite budget required"))
    repeats isa Integer && repeats>=1 || throw(ArgumentError("positive repeats required"))
    0<preparation_fraction<1 || throw(ArgumentError("preparation fraction must lie in (0,1)"))
    !isempty(methods) && all(m->m in (:scip,:highs_ipopt,:scip_ipopt),methods) ||
        throw(ArgumentError("unknown or empty workflow selection"))
    mkpath(output)
    source_hash=bytes2hex(sha256(join(read.(sort(filter(p->endswith(p,".jl"),
        readdir(joinpath(@__DIR__,"..","src");join=true))),String))))
    # Compile all three paths before measuring. Warm-up schedules are discarded.
    warmcase=gate_network_case(2)
    warmups=Any[]
    for method in methods
        println("Warming $method");flush(stdout)
        folder=joinpath(output,"warmup",string(method));mkpath(folder)
        push!(warmups,workflow_trial(warmcase,method;budget=300.0,output=folder,preparation_fraction))
    end
    records=Any[]
    summary=Dict("records"=>records,"warmups"=>warmups,"source_sha256"=>source_hash,
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>OpenSHOP.BLAS.get_num_threads(),"cpu"=>Sys.CPU_NAME,
        "timing_scope"=>"shared construction/preparation/SCIP allowance; final audits may overrun; loading and warmup excluded",
        "first_validated_scope"=>"first independent acceptance observed by the workflow; SCIP incumbents are audited only after optimize! returns",
        "quality_timing_scope"=>"post-hoc acceptance within 0.1% of the best delivered objective on this input in this campaign; not an optimality proof")
    for (label,c) in workflow_cases(labels)
        frozen=joinpath(output,label,"case.json");benchmark_freeze(frozen,case_dict(c))
        case_hash=bytes2hex(sha256(read(frozen)))
        # Each network exercises different numerical kernels. Warm its complete
        # paths too; no controls or models are reused by the measured trials.
        for method in methods
            folder=joinpath(output,label,"warmup",string(method));mkpath(folder)
            println("Warming $label $method");flush(stdout)
            workflow_trial(c,method;budget=min(30.0,budget),output=folder,preparation_fraction)
        end
        for repeat in 1:repeats
            # Rotate execution order to avoid always favoring the first method.
            for method in circshift(collect(methods),repeat-1)
                folder=joinpath(output,label,"repeat-$repeat",string(method));mkpath(folder)
                println("Starting $label $method repetition $repeat ($budget seconds)");flush(stdout)
                record=workflow_trial(c,method;budget,output=folder,preparation_fraction)
                record["label"]=label;record["repeat"]=repeat;record["case_sha256"]=case_hash
                push!(records,record)
                writejson(joinpath(output,"summary.json"),summary)
                println(JSON3.write(record));flush(stdout)
            end
        end
    end
    workflow_quality!(records)
    writejson(joinpath(output,"summary.json"),summary)
    records
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)>=4 || error("usage: workflows.jl OUTPUT SECONDS REPEATS CASE...")
    benchmark_workflows(abspath(ARGS[1]),ARGS[4:end];budget=parse(Float64,ARGS[2]),
        repeats=parse(Int,ARGS[3]))
end
