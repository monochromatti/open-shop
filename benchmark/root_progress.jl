# Observations at ordinary root LP boundaries. Local LP objectives are diagnostic
# only; ProofBounds independently records SCIP's globally valid upper bound.
mutable struct RootOBBTBudget
    seconds::Float64
    disabled::Bool
    disabled_seconds::Union{Nothing,Float64}
    disabled_calls::Union{Nothing,Int}
end

mutable struct RootProgress <: SCIP.AbstractEventhdlr
    optimizer::SCIP.Optimizer
    started_ns::UInt64
    kind::Symbol
    mask::UInt64
    budget::RootOBBTBudget
    points::Vector{Dict{String,Any}}
    calls::Int
    skipped::Int
    dropped::Int
    capture_seconds::Float64
    errors::Vector{String}
end

function SCIP.eventinit(e::RootProgress)
    SCIP.catch_event(e.optimizer.inner,e.mask,e)
end
function SCIP.eventexit(e::RootProgress)
    SCIP.drop_event(e.optimizer.inner,e.mask,e)
end

function root_progress_reliable(o)
    node=SCIP.SCIPgetFocusNode(o)
    node!=C_NULL && SCIP.SCIPnodeGetDepth(node)==0 &&
        !Bool(SCIP.SCIPinProbing(o)) && !Bool(SCIP.SCIPinDive(o)) &&
        !Bool(SCIP.SCIPinRepropagation(o)) &&
        SCIP.SCIPgetLPSolstat(o)==SCIP.SCIP_LPSOLSTAT_OPTIMAL &&
        Bool(SCIP.SCIPisLPPrimalReliable(o)) && Bool(SCIP.SCIPisLPDualReliable(o)) &&
        Bool(SCIP.SCIPisLPRelax(o))
end

function SCIP.eventexec(e::RootProgress)
    began=time_ns();e.calls+=1
    try
        o=e.optimizer
        if e.kind==:branch
            if !isempty(e.points)
                e.skipped+=1;return
            end
            node=SCIP.SCIPgetFocusNode(o)
            if node==C_NULL || SCIP.SCIPnodeGetDepth(node)!=0 ||
                    Bool(SCIP.SCIPinProbing(o)) || Bool(SCIP.SCIPinDive(o)) ||
                    Bool(SCIP.SCIPinRepropagation(o))
                e.skipped+=1;return
            end
            push!(e.points,Dict{String,Any}("event"=>"first_root_branch",
                "seconds"=>Float64(time_ns()-e.started_ns)/1e9,
                "native_seconds"=>Float64(SCIP.SCIPgetSolvingTime(o)),
                "run"=>Int(SCIP.SCIPgetNRuns(o)),"depth"=>0,
                "global_upper_bound"=>OpenSHOP._scip_bound_value(o,SCIP.SCIPgetDualbound(o))))
            return
        end
        if !root_progress_reliable(o)
            e.skipped+=1
            return
        end
        prop=SCIP.SCIPfindProp(o,"obbt")
        prop!=C_NULL || error("native OBBT propagator unavailable")
        seconds=Float64(SCIP.SCIPpropGetTime(prop))
        calls=Int(SCIP.SCIPpropGetNCalls(prop))
        if !e.budget.disabled && calls>=1 && seconds>=e.budget.seconds
            # A boundary stop, not an interrupt: one invocation can overrun.
            SCIP.@SCIP_CALL SCIP.SCIPsetIntParam(o,"propagating/obbt/freq",-1)
            SCIP.SCIPpropGetFreq(prop)==-1 || error("OBBT frequency change was not applied")
            e.budget.disabled=true
            e.budget.disabled_seconds=seconds
            e.budget.disabled_calls=calls
        end
        length(e.points)>=4096 && (e.dropped+=1;return)
        sepa=SCIP.SCIPfindSepa(o,"table_power")
        push!(e.points,Dict{String,Any}(
            "event"=>string(e.kind),"seconds"=>Float64(time_ns()-e.started_ns)/1e9,
            "native_seconds"=>Float64(SCIP.SCIPgetSolvingTime(o)),
            "run"=>Int(SCIP.SCIPgetNRuns(o)),"depth"=>0,
            "global_upper_bound"=>OpenSHOP._scip_bound_value(o,SCIP.SCIPgetDualbound(o)),
            "local_root_lp_objective"=>OpenSHOP._scip_bound_value(o,
                SCIP.SCIPretransformObj(o,SCIP.SCIPgetLPObjval(o))),
            "lp_iterations"=>Int(SCIP.SCIPgetNLPIterations(o)),
            "root_lp_iterations"=>Int(SCIP.SCIPgetNRootLPIterations(o)),
            "lps"=>Int(SCIP.SCIPgetNLPs(o)),
            "obbt_seconds"=>seconds,"obbt_calls"=>calls,
            "obbt_domain_reductions"=>Int(SCIP.SCIPpropGetNDomredsFound(prop)),
            "obbt_frequency"=>Int(SCIP.SCIPpropGetFreq(prop)),
            "power_cut_seconds"=>sepa==C_NULL ? 0.0 : Float64(SCIP.SCIPsepaGetTime(sepa)),
            "power_cut_calls"=>sepa==C_NULL ? 0 : Int(SCIP.SCIPsepaGetNCalls(sepa))))
    catch exception
        push!(e.errors,sprint(showerror,exception))
    finally
        e.capture_seconds+=Float64(time_ns()-began)/1e9
    end
    nothing
end

function install_root_progress!(optimizer,started_ns;obbt_cap_seconds=Inf)
    obbt_cap_seconds>=0 && !isnan(obbt_cap_seconds) || error("invalid cumulative OBBT allowance")
    budget=RootOBBTBudget(Float64(obbt_cap_seconds),false,nothing,nothing)
    observers=RootProgress[]
    # SCIP.jl hides the native event pointer/type; one handler per mask keeps
    # FIRSTLPSOLVED and LPSOLVED observations distinguishable.
    for (kind,mask) in ((:first_lp,SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED),
                       (:lp,SCIP.SCIP_EVENTTYPE_LPSOLVED),
                       (:branch,SCIP.SCIP_EVENTTYPE_NODEBRANCHED))
        e=RootProgress(optimizer,started_ns,kind,mask,budget,Dict{String,Any}[],0,0,0,0.0,String[])
        Base.precompile(SCIP.eventexec,(typeof(e),))
        SCIP.include_event_handler(optimizer.inner,e;desc="Ordinary root LP progress $(kind)")
        push!(observers,e)
    end
    observers
end

function root_progress_statistics(observers)
    points=sort!([p for e in observers if e.kind!=:branch for p in e.points];by=p->p["seconds"])
    budget=isempty(observers) ? nothing : first(observers).budget
    branches=[p for e in observers if e.kind==:branch for p in e.points]
    Dict("points"=>points,"first_root_branch"=>isempty(branches) ? nothing : first(branches),"calls"=>sum(e.calls for e in observers;init=0),
        "skipped"=>sum(e.skipped for e in observers;init=0),
        "dropped"=>sum(e.dropped for e in observers;init=0),
        "capture_seconds"=>sum(e.capture_seconds for e in observers;init=0.0),
        "errors"=>[err for e in observers for err in e.errors],
        "obbt_cap_seconds"=>budget===nothing || !isfinite(budget.seconds) ? nothing : budget.seconds,
        "obbt_disabled_at_boundary"=>budget===nothing ? false : budget.disabled,
        "obbt_seconds_at_disable"=>budget===nothing ? nothing : budget.disabled_seconds,
        "obbt_calls_at_disable"=>budget===nothing ? nothing : budget.disabled_calls,
        "scope"=>"optimal reliable root LPs outside probing/diving/repropagation; local objective is diagnostic only; cumulative plugin timers overlap LP timers; domain-reduction counts exclude probing changes and generalized bounds; cap is checked after completed OBBT calls and may overrun")
end
