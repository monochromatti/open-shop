# Confirm the proposal backend on the same formulation and dispatch algorithm.
# Loading, compilation and finer replay are measured separately from preparation.
include("workflows.jl")

function proposal_confirmation_cases()
    bases = Dict(workflow_cases(["distributed", "tables", "gates6", "gates24", "tokke6", "tokke24"]))
    gated = bases["gates24"]
    tokke = bases["tokke24"]
    volatile = gate_benchmark_copy(gated;
        name="gated network: alternating prices",
        prices=[isodd(t) ? -30.0 : 120.0 for t in eachindex(gated.prices)])
    wet = gate_benchmark_copy(gated; name="gated network: wet initial conditions",
        system=gate_benchmark_copy(gated.system; reservoirs=[
            gate_benchmark_copy(r; v0=r.v0+0.1*(r.vmax-r.v0), inflow=r.inflow+2.0)
            for r in gated.system.reservoirs]))
    restricted = gate_benchmark_copy(gated; name="gated network: longer dwell times",
        system=gate_benchmark_copy(gated.system; generators=[
            gate_benchmark_copy(g; minup=max(g.minup,4.0), mindown=max(g.mindown,3.0),
                startup=3g.startup) for g in gated.system.generators]))
    tokke_prices = gate_benchmark_copy(tokke; name="Tokke–Vinje: alternating prices",
        prices=[isodd(t) ? 0.5p : 1.5p for (t,p) in enumerate(tokke.prices)])
    tokke_wet = gate_benchmark_copy(tokke; name="Tokke–Vinje: wetter initial conditions",
        system=gate_benchmark_copy(tokke.system; reservoirs=[
            gate_benchmark_copy(r; v0=r.v0+0.02*(r.vmax-r.v0), inflow=1.2r.inflow)
            for r in tokke.system.reservoirs]))
    cases = vcat(collect(pairs(bases)), ["gates48"=>gate_network_case(48),
        "gates_volatile"=>volatile, "gates_wet"=>wet, "gates_restricted"=>restricted,
        "tokke_prices"=>tokke_prices, "tokke_wet"=>tokke_wet])
    sort!(cases; by=first)
    for (_,c) in cases
        OpenSHOP.validate_inputs(c)
    end
    cases
end

function proposal_factory(method)
    if method == :highs
        return optimizer_with_attributes(HiGHS.Optimizer,"threads"=>1)
    elseif method == :scip
        return optimizer_with_attributes(SCIP.Optimizer,"parallel/maxnthreads"=>1,
            "numerics/feastol"=>1e-8)
    end
    throw(ArgumentError("unknown proposal method: $method"))
end

function preparation_trial(c,method;budget,output)
    started=time()
    p=schedule_case(c;time_limit=budget,proposal_time_limit=min(5.0,budget/4),
        nlp_time_limit=min(20.0,budget),max_refinements=1,refine_grid=false,
        operational_margin=0.0,proposal_optimizer=proposal_factory(method))
    elapsed=time()-started
    writejson(joinpath(output,"preparation.json"),p)
    events=[Dict("seconds"=>a["elapsed_seconds"],"objective"=>last(a["steps"])["objective"])
        for a in p["attempts"] if get(a,"accepted",false) && !isempty(get(a,"steps",[]))]
    replay=nothing;replay_seconds=0.0;components=nothing
    if p["accepted"]
        x=p["solution"]
        replay_seconds=@elapsed replay=replay_audit(c,x;grid=OpenSHOP.refined_grid(c.grid;factor=4))
        dt=diff(c.grid)
        revenue=sum(c.prices[t]*dt[t]*x["power"][j,t]
            for j in eachindex(c.system.generators),t in eachindex(c.prices);init=0.0)
        components=Dict("electricity_revenue"=>revenue,
            "reservoir_value_change"=>sum(r.water_value*(x["V"][i,end]-r.v0)
                for (i,r) in enumerate(c.system.reservoirs);init=0.0),
            "transition_cost"=>OpenSHOP.transition_costs(c,x["u"]),
            "gate_cost"=>x["gate_adjustment_cost"],"release_penalty"=>x["release_penalty_cost"])
        components["river_value_change"]=x["objective"]-revenue-
            components["reservoir_value_change"]+components["transition_cost"]+
            components["gate_cost"]+components["release_penalty"]
    end
    Dict("method"=>string(method),"budget_seconds"=>budget,"total_seconds"=>elapsed,
        "accepted"=>p["accepted"],"objective"=>p["accepted"] ? p["solution"]["objective"] : nothing,
        "events"=>events,"objective_components"=>components,
        "proposal_seconds"=>sum(get(z,"total_seconds",0.0) for z in p["proposals"];init=0.0),
        "proposal_count"=>length(p["proposals"]),
        "proposal_statuses"=>[z["status"] for z in p["proposals"]],
        "dispatch_seconds"=>sum(get(a,"seconds",0.0) for a in p["attempts"];init=0.0),
        "dispatch_count"=>length(p["attempts"]),
        "attempt_errors"=>[a["error"] for a in p["attempts"] if haskey(a,"error")],
        "fourfold_replay"=>replay===nothing ? false : replay["valid"],
        "fourfold_errors"=>replay===nothing ? [] : replay["errors"],
        "fourfold_seconds"=>replay_seconds,
        "fourfold_storage_difference_Mm3"=>replay===nothing ? nothing : get(replay,"max_storage_difference_Mm3",nothing))
end

function preparation_quality!(records)
    # Use one reference across allowances so a poor short-run fallback is not
    # counted as useful merely because neither method had time to improve it.
    for label in unique(r["label"] for r in records)
        group=filter(r->r["label"]==label,records)
        values=[r["objective"] for r in group if r["objective"]!==nothing && r["fourfold_replay"]]
        isempty(values) && continue
        best=maximum(values)
        for r in group
            r["best_observed_objective"]=best
            for fraction in (0.001,0.01)
                threshold=best-fraction*max(1.0,abs(best))
                times=[e["seconds"] for e in r["events"] if e["objective"]>=threshold]
                r["time_within_$(fraction)_of_best"]=isempty(times) ? nothing : minimum(times)
            end
        end
    end
    records
end

function confirm_proposal_solvers(output;repeats=3)
    repeats isa Integer && repeats>=1 || throw(ArgumentError("positive repeats required"))
    mkpath(output)
    records=Any[]
    source_hash=bytes2hex(sha256(join(read.(sort(filter(p->endswith(p,".jl"),
        readdir(joinpath(@__DIR__,"..","src");join=true))),String))))
    summary=Dict("records"=>records,"source_sha256"=>source_hash,
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>OpenSHOP.BLAS.get_num_threads(),"cpu"=>Sys.CPU_NAME,
        "timing_scope"=>"preparation only; loading and warmups excluded; fourfold replay outside measured preparation",
        "quality_scope"=>"first independently accepted preparation event within 0.1% or 1% of best fourfold-valid delivered objective on the same input across all measured budgets; empirical, not optimality")
    for method in (:highs,:scip)
        println("Compiling $method");flush(stdout)
        folder=joinpath(output,"warmup",string(method));mkpath(folder)
        preparation_trial(gate_network_case(2),method;budget=120.0,output=folder)
    end
    for (label,c) in proposal_confirmation_cases()
        frozen=joinpath(output,label,"case.json");benchmark_freeze(frozen,case_dict(c))
        input_hash=bytes2hex(sha256(read(frozen)))
        # Exercise every network's kernels without retaining its schedule.
        for method in (:highs,:scip)
            folder=joinpath(output,label,"warmup",string(method));mkpath(folder)
            println("Warming $label $method");flush(stdout)
            preparation_trial(c,method;budget=30.0,output=folder)
        end
        budgets=label in ("gates24","gates48","tokke6","tokke24") ? (24.0,60.0) : (60.0,)
        for budget in budgets, repeat in 1:repeats
            for method in circshift([:highs,:scip],repeat-1)
                folder=joinpath(output,label,"budget-$(Int(budget))","repeat-$repeat",string(method));mkpath(folder)
                println("Starting $label $method $budget seconds repetition $repeat");flush(stdout)
                record=preparation_trial(c,method;budget,output=folder)
                merge!(record,Dict("label"=>label,"repeat"=>repeat,"case_sha256"=>input_hash))
                push!(records,record)
                preparation_quality!(records)
                writejson(joinpath(output,"summary.json"),summary)
                println(JSON3.write(record));flush(stdout)
            end
        end
    end
    summary
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: proposal-solvers.jl OUTPUT")
    confirm_proposal_solvers(abspath(ARGS[1]))
end
