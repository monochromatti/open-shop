module Qualification

using OpenSHOP, JSON3, Random, LinearAlgebra
include("fuzz_worker.jl")
include("../examples/rolling_horizon.jl")
const updated = FuzzWorker.FuzzCases.updated

# Inputs are saved before reconstruction. Reference controls are used only by
# the generous-budget comparison, never by a short cold solve.
function input_data(request)
    if haskey(request,"input")
        raw=JSON3.read(read(request["input"],String),Dict{String,Any})
        return haskey(raw,"controls") ? FuzzWorker.load_input(request["input"]) :
            (;case=case_from_dict(raw),u=nothing)
    end
    seed=parse(Int,request["seed"])
    family=request["family"]
    if family!="integrated"
        return FuzzWorker.FuzzCases.generate_case(seed,family,request["size"])
    end
    rng=Xoshiro(seed)
    data=basin_case(;hours=get(request,"hours",48),branches=get(request,"branches",4),
        pond_volume=get(request,"pond_volume",0.002),distributed=get(request,"distributed",true))
    c=data.case
    reservoirs=[updated(r;water_value=r.water_value*(0.8+0.4rand(rng))) for r in c.system.reservoirs]
    tunnels=[updated(e;resistance=e.resistance*(0.7+0.6rand(rng))) for e in c.system.tunnels]
    operations=[updated(z;values=z.values.*(0.7+0.6rand(rng))) for z in c.operations]
    prices=c.prices.*(0.8+0.4rand(rng)).+10randn(rng,length(c.prices))
    merge(data,(;case=updated(c;system=updated(c.system;reservoirs,tunnels),operations,prices)))
end

function compare(rows,reference)
    best=reference===nothing ? nothing : reference.objective
    [Dict("label"=>row["label"],"method"=>get(row,"method","schedule_case"),
        "budget"=>get(row,"budget",nothing),"accepted"=>row["accepted"],
        "qualified"=>row["qualified"],"seconds"=>row["seconds"],
        "objective"=>row["objective"],"reference_objective"=>best,
        "value_shortfall"=>!row["qualified"] || best===nothing ? nothing : max(0.,best-row["objective"]),
        "relative_shortfall"=>!row["qualified"] || best===nothing ? nothing :
            max(0.,best-row["objective"])/max(1.,abs(best))) for row in rows]
end

# An operating revision changes future inputs, never the executed history.
function revise(c,event)
    grid=c.grid[1:end-1]
    operations=copy(c.operations)
    replace_series(z)=begin
        filter!(old->(old.object,old.attribute)!=(z.object,z.attribute),operations)
        push!(operations,z)
    end
    if event=="forecast"
        for r in c.system.reservoirs
            values=[0.7OpenSHOP.opvalue(c,r.name,:inflow,t,r.inflow) for t in grid]
            replace_series(OperationalSeries(;object=r.name,attribute=:inflow,times=grid,values))
        end
        return updated(c;operations,prices=reverse(c.prices).+20.0,name=c.name*" forecast revision")
    end
    event in ("outage","restriction") || throw(ArgumentError("unknown recovery event: $event"))
    # Do not create contradictory input by overriding a prescribed discharge,
    # forced commitment or power schedule. Such units are not outage test sites.
    eligible=[g for g in c.system.generators if g.mode==:generate &&
        !any(z->z.object==g.name && z.attribute in (:discharge,:power,:forced_on),c.operations)]
    isempty(eligible) && throw(ArgumentError("no unconstrained generating unit for $event"))
    g=first(eligible)
    a=c.grid[min(3,length(c.grid)-1)]
    b=min(last(c.grid),a+max(4.0,g.mindown))
    attribute=event=="outage" ? :maintenance : :qmax
    default=event=="outage" ? 0.0 : g.qmax
    values=[begin
        original=OpenSHOP.opvalue(c,g.name,attribute,t,default)
        a<=t<b ? (event=="outage" ? 1.0 : min(original,max(g.qmin,0.6g.qmax))) : original
    end for t in grid]
    replace_series(OperationalSeries(;object=g.name,attribute,times=grid,values))
    updated(c;operations,name=c.name*" $event")
end

function run(request,directory)
    BLAS.set_num_threads(1)
    budgets=Float64.(get(request,"budgets",[30.,90.]))
    reference_budget=Float64(get(request,"reference_budget",300.))
    repetitions=Int(get(request,"repetitions",2))
    all(b->isfinite(b) && b>0,budgets) && !isempty(budgets) && length(unique(budgets))==length(budgets) &&
        isfinite(reference_budget) && reference_budget>=maximum(budgets) && repetitions>0 ||
        throw(ArgumentError("positive budgets/repetitions and reference_budget >= all short budgets required"))
    report=Dict{String,Any}("family"=>request["family"],"size"=>request["size"],
        "seed"=>request["seed"],"findings"=>Any[],"runs"=>Any[],"recovery"=>Any[],
        "scope"=>"feasible-start scheduling, independent of global bound closure",
        "versions"=>Dict("Julia"=>string(VERSION),"OpenSHOP"=>string(pkgversion(OpenSHOP)),
            "JuMP"=>string(pkgversion(OpenSHOP.JuMP)),"SCIP.jl"=>string(pkgversion(FuzzWorker.SCIP)),
            "Ipopt.jl"=>string(pkgversion(FuzzWorker.Ipopt))))
    checkpoint()=FuzzWorker.checkpoint(joinpath(directory,"partial.json"),report)
    progress(phase)=FuzzWorker.checkpoint(joinpath(directory,"progress.json"),Dict("phase"=>phase))
    function measure(c,label,budget;initial=nothing,warm=nothing,global_search=false,generous=false)
        progress(label)
        println("BEGIN $label budget=$budget intervals=$(length(c.prices))");flush(stdout)
        row=Dict{String,Any}("label"=>label,"budget"=>budget,"accepted"=>false,
            "qualified"=>false,"objective"=>nothing,"method"=>global_search ? "solve" : "schedule_case")
        x=nothing
        started=time()
        try
            timed=@timed if global_search
                solve(c;initial,time_limit=budget,diagnostics_path=joinpath(directory,label*".scip.log"))
            else
                schedule_case(c;initial,warm,time_limit=budget,proposal_time_limit=min(5.,budget),
                    nlp_time_limit=generous ? budget : min(45.,budget),max_refinements=1)
            end
            result=timed.value
            writejson(joinpath(directory,label*".json"),result)
            row["seconds"]=timed.time
            row["compile_seconds"]=timed.compile_time
            row["recompile_seconds"]=timed.recompile_time
            row["bytes"]=timed.bytes
            row["accepted"]=result["accepted"]
            if global_search
                row["status"]=result.status
                row["global_bound"]=result.global_bound
                row["global_certificate"]=result.global_certificate
            end
            attempts=global_search ? Dict() : result
            row["attempts"]=[Dict(k=>v for (k,v) in a if k ∉ ("u","tunnel_opening"))
                for a in get(attempts,"attempts",Any[])]
            row["first_solver_accepted_seconds"]=nothing
            for a in row["attempts"]
                if get(a,"accepted",false)
                    row["first_solver_accepted_seconds"]=get(a,"elapsed_seconds",nothing)
                    break
                end
            end
            if result["accepted"]
                candidate=result["solution"]
                row["objective"]=candidate.objective
                checking=@timed FuzzWorker.check_schedule!(report,c,candidate,label,4)
                checks=checking.value
                row["campaign_check_seconds"]=checking.time
                row["checks"]=checks
                row["qualified"]=checks["audit"]["valid"] && checks["replay"]["valid"] &&
                    checks["water_error_Mm3"]<=1e-6 &&
                    !any(f->f["stage"]==label && f["severity"]=="error",report["findings"])
                row["qualified"] && (x=candidate)
            end
        catch error
            row["exception"]=sprint(showerror,error,catch_backtrace())
            FuzzWorker.finding!(report,"exception",label,row["exception"])
        end
        row["seconds"]=get(row,"seconds",time()-started)
        println("END $label accepted=$(row["accepted"]) qualified=$(row["qualified"]) seconds=$(row["seconds"])");flush(stdout)
        push!(report["runs"],row)
        checkpoint()
        GC.gc()
        row,x
    end
    bestof(x,y)=x===nothing ? y : y===nothing || x.objective>=y.objective ? x : y
    try
        progress("generation")
        data=input_data(request)
        c=case_from_dict(case_dict(data.case))
        writejson(joinpath(directory,"input.json"),data.u===nothing ? case_dict(c) :
            FuzzWorker.input_dict(merge(data,(;case=c))))
        report["features"]=FuzzWorker.topology_features(c)
        witness=nothing
        if data.u!==nothing
            progress("reference_controls")
            witness=FuzzWorker.reconstruct(merge(data,(;case=c)))
            validate(c,witness)["valid"] || error("generated reference violates original equations or restrictions")
            if get(request,"restricted",false)
                c=FuzzWorker.FuzzCases.restrict_case(parse(Int,request["seed"]),c,witness)
                writejson(joinpath(directory,"input.json"),FuzzWorker.input_dict(merge(data,(;case=c))))
                witness=FuzzWorker.reconstruct(merge(data,(;case=c)))
            end
            checks=FuzzWorker.check_schedule!(report,c,witness,"reference_controls",4)
            report["reference_controls"]=checks
            checks["audit"]["valid"] && checks["replay"]["valid"] &&
                checks["water_error_Mm3"]<=1e-6 &&
                !any(f->f["stage"]=="reference_controls" && f["severity"]=="error",report["findings"]) ||
                (witness=nothing)
        end
        # Compile both public paths separately from the measured case. Remaining
        # specialization is reported per call, rather than assumed to be zero.
        progress("compilation")
        warmdata=FuzzWorker.FuzzCases.generate_case(1,"cascade","small")
        warmref=FuzzWorker.reconstruct(warmdata)
        compilation=@elapsed begin
            schedule_case(warmdata.case;time_limit=60.,nlp_time_limit=20.,max_refinements=1)
            solve(warmdata.case;initial=warmref,time_limit=5.)
        end
        report["warmup_seconds"]=compilation
        report["features"]=FuzzWorker.topology_features(c)
        checkpoint()
        rows=Any[]
        best=witness
        # Reverse the budget order on alternate repetitions to expose order and
        # residual compilation effects. No result is passed to another cold run.
        for repetition in 1:repetitions, budget in (isodd(repetition) ? budgets : reverse(budgets))
            row,x=measure(c,"cold-$(budget)-r$repetition",budget)
            push!(rows,row)
            best=bestof(best,x)
        end
        for budget in Float64.(get(request,"public_budgets",Float64[]))
            isfinite(budget) && budget>0 || throw(ArgumentError("positive finite public budget required"))
            row,x=measure(c,"public-$(budget)",budget;global_search=true)
            push!(rows,row)
            best=bestof(best,x)
        end
        _,x=measure(c,"reference-preparation",reference_budget;initial=best,generous=true)
        best=bestof(best,x)
        if get(request,"global_reference",false)
            _,x=measure(c,"reference-global",reference_budget;initial=best,global_search=true)
            best=bestof(best,x)
        end
        report["comparisons"]=compare(rows,best)
        best!==nothing && writejson(joinpath(directory,"reference.json"),best)
        report["reference_objective"]=best===nothing ? nothing : best.objective
        for (revision,event) in enumerate(get(request,"recovery",String[]))
            if best===nothing || length(c.prices)<4
                push!(report["recovery"],Dict("event"=>event,"status"=>"no executable source schedule"))
                break
            end
            edge=clamp(searchsortedfirst(c.grid,first(c.grid)+get(request,"advance",6.)),2,length(c.grid)-2)
            tail=restart_case(c,best,c.grid[edge])
            continued=dispatch_from_controls(tail,best.commitment[:,edge:end],best.generator_q[:,edge:end],
                best.gate[:,edge:end];tunnel_opening=best.tunnel_opening[:,edge:end],head_reference=best.head[:,edge:end])
            validate(tail,continued)["valid"] &&
                maximum(abs,continued.storage-best.storage[:,edge:end];init=0.)<=1e-6 ||
                error("unchanged continuation differs from executed schedule")
            next=revise(tail,event)
            writejson(joinpath(directory,"recovery-$revision-input.json"),case_dict(next))
            guess=shifted_start(best,next.grid)
            budget=maximum(budgets)
            paired=Any[]
            revised_best=nothing
            # The old trajectory is a numerical guess, never a feasible fallback
            # after changed input. Both solves see the same executed state.
            for mode in (isodd(revision) ? ("cold","warm") : ("warm","cold"))
                row,x=measure(next,"recovery-$revision-$mode",budget;warm=mode=="warm" ? guess : nothing)
                push!(paired,row)
                revised_best=bestof(revised_best,x)
            end
            short_success=revised_best!==nothing
            _,x=measure(next,"recovery-$revision-reference",reference_budget;
                initial=revised_best,warm=revised_best===nothing ? guess : nothing,generous=true)
            revised_best=bestof(revised_best,x)
            push!(report["recovery"],Dict("event"=>event,"start"=>first(next.grid),
                "end"=>last(next.grid),"short_budget_recovered"=>short_success,
                "reference_recovered"=>revised_best!==nothing,"comparisons"=>compare(paired,revised_best)))
            c,best=next,revised_best
            checkpoint()
        end
    catch error
        FuzzWorker.finding!(report,"exception","campaign",sprint(showerror,error,catch_backtrace()))
    end
    report["outcome"]=any(f->f["severity"]=="error",report["findings"]) ? "fail" :
        any(!r["qualified"] for r in report["runs"]) || !isempty(report["findings"]) ? "review" : "pass"
    FuzzWorker.checkpoint(joinpath(directory,"result.json"),report)
    report
end

function main(args)
    length(args)==1 || error("usage: qualify.jl REQUEST.json")
    path=abspath(only(args))
    request=JSON3.read(read(path,String),Dict{String,Any})
    run(request,dirname(path))
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    Qualification.main(ARGS)
end
