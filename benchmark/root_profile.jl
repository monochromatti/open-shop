# Diagnose and screen interventions within the single production SCIP model.
# Raw inputs, controls and root solutions stay local; summary.json is scalar.
include("diagnose.jl")
using SCIP, JuMP
include("power_envelopes.jl")

mutable struct RootLP <: SCIP.AbstractEventhdlr
    optimizer::SCIP.Optimizer
    indices::Vector{JuMP.MOI.VariableIndex}
    first::Vector{Float64}
    last::Vector{Float64}
    first_seconds::Float64
    last_seconds::Float64
end
function SCIP.eventinit(e::RootLP)
    SCIP.catch_event(e.optimizer.inner,SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP.SCIP_EVENTTYPE_LPSOLVED,e)
end
function SCIP.eventexit(e::RootLP)
    SCIP.drop_event(e.optimizer.inner,SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP.SCIP_EVENTTYPE_LPSOLVED,e)
end
function SCIP.eventexec(e::RootLP)
    o=e.optimizer
    node=SCIP.SCIPgetFocusNode(o)
    node==C_NULL && return
    SCIP.SCIPnodeGetDepth(node)==0 || return
    Bool(SCIP.SCIPinProbing(o)) && return
    SCIP.SCIPgetLPSolstat(o)==SCIP.SCIP_LPSOLSTAT_OPTIMAL || return
    vals=SCIP.sol_values(o,e.indices)
    if isempty(e.first)
        e.first=vals;e.first_seconds=SCIP.SCIPgetSolvingTime(o)
    else
        e.last=vals;e.last_seconds=SCIP.SCIPgetSolvingTime(o)
    end
end

function relaxation_summary(b,c,variables,vals)
    isempty(vals) && return nothing
    assigned=Dict(zip(variables,vals));at(x)=x isa Number ? x : value(v->assigned[v],x)
    frac=maximum(abs(at(x)-round(at(x))) for x in b.u;init=0.)
    sos=maximum(OpenSHOP._sos2_residual(at.(axis.weights),collect(1.:length(axis.weights)))
        for axis in get(b.m.ext,:global_tensor_coordinates,[]);init=0.)
    generated=0.;physical=0.;worst=0.
    level_error=maximum(abs(at(b.H[i,t])-OpenSHOP.head(r,
        (at(b.V[i,t])+at(b.V[i,t+1]))/2)) for (i,r) in enumerate(c.system.reservoirs),
        t in eachindex(c.prices);init=0.)
    ix=OpenSHOP.nodeindex(c.system)
    tunnel_error=maximum(abs(OpenSHOP.opinterval(c,e.name,:opening,t,e.opening)*
        (at(b.H[ix[e.source],t])-at(b.H[ix[e.target],t]))-
        e.resistance*at(b.Q[i,t])*abs(at(b.Q[i,t])))
        for (i,e) in enumerate(c.system.tunnels),t in eachindex(c.prices);init=0.)
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices)
        q=at(b.GQ[i,t]);p=at(b.P[i,t])
        h=at(b.shared_heads[(OpenSHOP.plantof(c.system,g).name,t)])
        eta=g.turbine_table===nothing ? OpenSHOP.efficiency(g,q,h) :
            OpenSHOP.turbine_efficiency(g.turbine_table,q,h;extrapolation=:linear)
        electrical=g.generator_efficiency_curve===nothing ? 1. :
            OpenSHOP.table_value(g.generator_efficiency_curve,p;extrapolation=:linear)
        actual=0.00981*q*h*eta*electrical
        generated+=p;physical+=actual;worst=max(worst,abs(p-actual))
    end
    Dict("fractional_commitment_max"=>frac,"sos2_residual_max"=>sos,
        "sum_reported_power_mw"=>generated,"sum_power_from_actual_curves_mw"=>physical,
        "power_equation_residual_max_mw"=>worst,
        "reservoir_level_residual_max_m"=>level_error,
        "tunnel_head_loss_residual_max_m"=>tunnel_error,
        "scope"=>"diagnostic LP solution; not a feasible schedule or a certificate")
end

const ROOT_PROFILES=Dict{String,Vector{Pair{String,Any}}}(
    "baseline"=>[],
    "root20"=>["separating/maxroundsroot"=>20],
    "no_obbt"=>["propagating/obbt/freq"=>-1],
    "lean_obbt"=>["propagating/obbt/itlimitfactor"=>1.0,"propagating/obbt/minitlimit"=>1000],
    "root20_no_obbt"=>["separating/maxroundsroot"=>20,"propagating/obbt/freq"=>-1],
    "power_envelopes"=>[],
)

function root_profile(input,output;seconds=120.,repeats=1,
        profiles=["baseline"],capture=true,commitment="free")
    mkpath(output)
    c=readcase(input)
    benchmark_freeze(joinpath(output,"case.json"),case_dict(c))
    seedfile=joinpath(output,"seed.json")
    preparation_seconds=0.
    if !isfile(seedfile)
        preparation_seconds=@elapsed prep=schedule_case(c;proposal_time_limit=5.,nlp_time_limit=20.,max_refinements=0,operational_margin=.1)
        prep["accepted"] || error("no audited common seed")
        seed=prep["solution"]
        q=copy(seed["generator_q"]);q[seed["u"].==0].=0.
        seed=dispatch_from_controls(c,seed["u"],q,seed["gate"])
        seed["validation"]["valid"] || error("reconstructed seed invalid")
        benchmark_freeze(seedfile,seed)
    end
    initial=benchmark_seed(seedfile,c)
    fixed_u=commitment=="fixed" ? initial["u"] : nothing
    # Warm setup, solver paths and extraction before matched measurement.
    warm=OpenSHOP._solve(c;initial,fixed_u,time_limit=15.)
    warm["status"]=="CONSTRUCTION_BUDGET_EXHAUSTED" && OpenSHOP._solve(c;initial,fixed_u,time_limit=10.)
    rows=Any[]
    for repeat in 1:repeats,profile in (isodd(repeat) ? profiles : reverse(profiles))
        haskey(ROOT_PROFILES,profile) || error("unknown experimental profile")
        graph=Ref{Any}(nothing);event=Ref{Any}(nothing);vars=Ref{Any}(nothing)
        setup=function(b)
            graph[]=b
            for (key,val) in ROOT_PROFILES[profile]
                set_optimizer_attribute(b.m,key,val)
            end
            if capture
                vars[]=all_variables(b.m)
                o=unsafe_backend(b.m)
                e=RootLP(o,optimizer_index.(vars[]),Float64[],Float64[],0.,0.)
                SCIP.include_event_handler(o.inner,e;desc="Root LP diagnostic snapshots")
                event[]=e
            end
        end
        transform=profile=="power_envelopes" ? add_power_envelopes! : nothing
        if repeat==1
            OpenSHOP._solve(c;initial,fixed_u,time_limit=10.,
                optimizer_setup=setup,model_transform=transform)
            graph[]=nothing;event[]=nothing;vars[]=nothing
        end
        logpath=joinpath(output,"$(profile)-$(commitment)-$(repeat).log")
        r=OpenSHOP._solve(c;initial,fixed_u,time_limit=seconds,relative_gap=1e-4,
            diagnostics_path=logpath,optimizer_setup=setup,model_transform=transform)
        row=Dict(k=>get(r,k,nothing) for k in ("status","accepted","global_certificate",
            "feasible_lower_bound","global_bound","relative_gap","construction_seconds",
            "solve_seconds","total_seconds","variable_count","constraint_count","scip_diagnostics",
            "scip_statistics","start_audit","bound_rejection","solver_error","diagnostics_close_error"))
        row["profile"]=profile;row["repeat"]=repeat;row["commitment"]=commitment
        row["seed_objective"]=initial["objective"];row["seed_sha256"]=bytes2hex(sha256(read(seedfile)))
        row["case_sha256"]=bytes2hex(sha256(read(joinpath(output,"case.json"))))
        row["preparation_seconds_excluded"]=preparation_seconds
        row["profile_warmup_seconds_excluded"]=repeat==1 ? 10. : 0.
        row["source_sha256"]=bytes2hex(sha256(join(read(p,String) for p in
            sort(filter(p->endswith(p,".jl"),readdir(joinpath(@__DIR__,"..","src");join=true))))))
        row["julia_version"]=string(VERSION)
        row["cpu_name"]=Sys.CPU_NAME;row["kernel"]=string(Sys.KERNEL)
        row["progress"]=scip_progress(logpath)
        if capture && event[]!==nothing && graph[]!==nothing
            e=event[];b=graph[]
            row["first_root_lp"]=relaxation_summary(b,c,vars[],e.first)
            row["first_root_lp_snapshot_seconds"]=e.first_seconds
            if get(r["scip_diagnostics"],"total_nodes",0)==1 && !Bool(SCIP.SCIPinProbing(e.optimizer)) &&
                    SCIP.SCIPgetLPSolstat(e.optimizer)==SCIP.SCIP_LPSOLSTAT_OPTIMAL
                e.last=SCIP.sol_values(e.optimizer,e.indices)
                e.last_seconds=SCIP.SCIPgetSolvingTime(e.optimizer)
            end
            row["last_root_lp"]=relaxation_summary(b,c,vars[],e.last)
            row["last_root_lp_snapshot_seconds"]=e.last_seconds
        end
        push!(rows,row)
        writejson(joinpath(output,"summary.json"),rows)
        println(now()," ",profile," ",commitment," ",r["status"]," gap=",r["relative_gap"]);flush(stdout)
        get(r,"accepted",false) && get(get(r,"start_audit",Dict()),"valid",false) || error("unaccepted schedule/start")
        get(r,"bound_rejection",nothing)===nothing || error("rejected upper bound")
    end
    rows
end

if abspath(PROGRAM_FILE)==@__FILE__
    root_profile(abspath(ARGS[1]),abspath(ARGS[2]);seconds=parse(Float64,ARGS[3]),
        repeats=parse(Int,ARGS[4]),profiles=split(ARGS[5],','),
        capture=length(ARGS)<6 || ARGS[6]=="true",commitment=length(ARGS)<7 ? "free" : ARGS[7])
end
