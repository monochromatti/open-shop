# Isolated relaxation experiments; never a public solver policy.
# Raw inputs, controls and root solutions stay local; summary.json is scalar.
include("diagnose.jl")
using SCIP, JuMP
include("table_consistency.jl")
using .TableConsistency
include("product_hull.jl")
using .ProductHull
include("network_energy.jl")

mutable struct RootLP <: SCIP.AbstractEventhdlr
    optimizer::SCIP.Optimizer
    indices::Vector{JuMP.MOI.VariableIndex}
    variables::Vector{VariableRef}
    objective::Any
    first_state::Dict{String,Any}
    last_state::Dict{String,Any}
    first::Vector{Float64}
    last::Vector{Float64}
    first_seconds::Float64
    last_seconds::Float64
    snapshots::Int
    capture_seconds::Float64
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
    (Bool(SCIP.SCIPinProbing(o)) || Bool(SCIP.SCIPinDive(o))) && return
    SCIP.SCIPgetLPSolstat(o)==SCIP.SCIP_LPSOLSTAT_OPTIMAL || return
    began=time()
    vals=SCIP.sol_values(o,e.indices)
    e.snapshots+=1
    assigned=Dict(zip(e.variables,vals))
    physical=JuMP.value(v->assigned[v],e.objective)
    native=10000*SCIP.SCIPgetSolOrigObj(o,C_NULL)
    affine=e.objective isa Union{Number,VariableRef,GenericAffExpr}
    state=Dict{String,Any}("values_finite"=>all(isfinite,vals),
        "lp_is_relaxation"=>Bool(SCIP.SCIPisLPRelax(o)),
        "lp_primal_reliable"=>Bool(SCIP.SCIPisLPPrimalReliable(o)),
        "scip_run"=>SCIP.SCIPgetNRuns(o),
        "native_original_lp_objective"=>native,
        "captured_objective_expression"=>physical,
        "objective_affine"=>affine,
        "objective_difference"=>physical-native,
        "objective_consistent"=>affine ? abs(physical-native)<=1e-7*max(1.,abs(native)) : nothing,
        "objective_check_scope"=>affine ? "same original affine revenue objective" : "nonlinear objective also relaxed; equality is not expected")
    if isempty(e.first)
        e.first=vals;e.first_seconds=SCIP.SCIPgetSolvingTime(o);e.first_state=state
    else
        e.last=vals;e.last_seconds=SCIP.SCIPgetSolvingTime(o);e.last_state=state
    end
    e.capture_seconds+=time()-began
end

function relaxation_summary(b,c,variables,vals)
    isempty(vals) && return nothing
    assigned=Dict(zip(variables,vals));at(x)=x isa Number ? x : value(v->assigned[v],x)
    frac=maximum(abs(at(x)-round(at(x))) for x in b.u;init=0.)
    sos=maximum(OpenSHOP._sos2_residual(at.(axis.weights),collect(1.:length(axis.weights)))
        for axis in get(b.m.ext,:global_tensor_coordinates,[]);init=0.)
    generated=0.;physical=0.;worst=0.
    decomposition=Any[]
    byname=Dict(name(v)=>v for v in variables)
    priceweighted=Dict("multiplication"=>0.,"turbine_table"=>0.,"electrical_table"=>0.)
    absoluteweighted=copy(priceweighted)
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
        etavar=at(byname[g.turbine_table===nothing ? "eta_$(i)_$(t)" : "turbine_$(i)_$(t)"])
        evar=g.generator_efficiency_curve===nothing ? 1. : at(byname["electrical_$(i)_$(t)"])
        parts=(multiplication=p-0.00981*q*h*etavar*evar,
            turbine_table=0.00981*q*h*(etavar-eta)*evar,
            electrical_table=0.00981*q*h*eta*(evar-electrical))
        multiplier=c.prices[t]*(c.grid[t+1]-c.grid[t])
        for key in keys(parts)
            priceweighted[string(key)]+=multiplier*parts[key]
            absoluteweighted[string(key)]+=abs(multiplier*parts[key])
        end
        push!(decomposition,Dict("unit"=>string(g.name),"interval"=>t,
            "reported_power_mw"=>p,"turbine_efficiency_residual"=>etavar-eta,
            "electrical_efficiency_residual"=>evar-electrical,
            "multiplication_residual_mw"=>parts.multiplication,
            "table_contribution_mw"=>parts.turbine_table,
            "electrical_contribution_mw"=>parts.electrical_table,
            "original_residual_mw"=>p-actual,
            "decomposition_error_mw"=>sum(parts)-(p-actual),
            "priceweighted_residual"=>multiplier*(p-actual)))
    end
    Dict("fractional_commitment_max"=>frac,"sos2_residual_max"=>sos,
        "sum_reported_power_mw"=>generated,"sum_power_from_actual_curves_mw"=>physical,
        "power_equation_residual_max_mw"=>worst,
        "power_residual_priceweighted_components"=>priceweighted,
        "power_residual_absolute_priceweighted_components"=>absoluteweighted,
        "power_residual_decomposition_error_max_mw"=>maximum(abs(d["decomposition_error_mw"]) for d in decomposition;init=0.),
        "largest_power_relaxation_errors"=>first(sort!(decomposition;by=d->abs(d["priceweighted_residual"]),rev=true),min(12,length(decomposition))),
        "network_energy"=>network_energy_diagnostic(b,c,at),
        "reservoir_level_residual_max_m"=>level_error,
        "tunnel_head_loss_residual_max_m"=>tunnel_error,
        "scope"=>"diagnostic LP solution; not a feasible schedule or a certificate")
end

const ROOT_PROFILES=Dict{String,Vector{Pair{String,Any}}}(
    "baseline"=>[],
    "rlt10"=>["separating/rlt/maxunknownterms"=>10],
    "rlt20"=>["separating/rlt/maxunknownterms"=>20],
    "rlt_visible"=>["separating/rlt/onlyoriginal"=>false],
    "rlt_visible20"=>["separating/rlt/maxunknownterms"=>20,"separating/rlt/onlyoriginal"=>false],
    "table_bilinear"=>[],
    "table_all"=>[],
    "energy"=>[],
    "product_hull"=>[],
    "product_both"=>[],
    "product_energy"=>[],
)

function profile_transform(profile)
    profile=="table_bilinear" && return (b,c)->add_table_consistency!(b,c;mode=:bilinear)
    profile=="table_all" && return (b,c)->add_table_consistency!(b,c;mode=:all)
    profile=="energy" && return add_network_energy!
    profile=="product_hull" && return add_product_hull!
    profile=="product_both" && return (b,c)->add_product_hull!(b,c;lower_power=true)
    profile=="product_energy" && return (b,c)->begin add_product_hull!(b,c;lower_power=true);add_network_energy!(b,c) end
    nothing
end

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
    common_warm_seconds=@elapsed begin
        warm=OpenSHOP._solve(c;initial,fixed_u,time_limit=15.,diagnostics_path=joinpath(output,"warm-common.log"))
        warm["status"]=="CONSTRUCTION_BUDGET_EXHAUSTED" && OpenSHOP._solve(c;initial,fixed_u,time_limit=10.,diagnostics_path=joinpath(output,"warm-common-retry.log"))
    end
    rows=Any[]
    for repeat in 1:repeats,profile in (isodd(repeat) ? profiles : reverse(profiles))
        haskey(ROOT_PROFILES,profile) || error("unknown experimental profile")
        try
        graph=Ref{Any}(nothing);event=Ref{Any}(nothing);vars=Ref{Any}(nothing)
        setup=function(b)
            graph[]=b
            for (key,val) in ROOT_PROFILES[profile]
                set_optimizer_attribute(b.m,key,val)
            end
            if capture
                vars[]=all_variables(b.m)
                o=unsafe_backend(b.m)
                e=RootLP(o,optimizer_index.(vars[]),vars[],b.obj,Dict{String,Any}(),Dict{String,Any}(),Float64[],Float64[],0.,0.,0,0.)
                SCIP.include_event_handler(o.inner,e;desc="Root LP diagnostic snapshots")
                event[]=e
            end
        end
        transform=profile_transform(profile)
        profile_warm_seconds=0.
        if repeat==1
            profile_warm_seconds=@elapsed OpenSHOP._solve(c;initial,fixed_u,time_limit=10.,
                diagnostics_path=joinpath(output,"warm-$(profile).log"),
                optimizer_setup=setup,model_transform=transform)
            graph[]=nothing;event[]=nothing;vars[]=nothing
        end
        logpath=joinpath(output,"$(profile)-$(commitment)-$(repeat).log")
        r=OpenSHOP._solve(c;initial,fixed_u,time_limit=seconds,relative_gap=1e-4,
            diagnostics_path=logpath,optimizer_setup=setup,model_transform=transform)
        row=Dict(k=>get(r,k,nothing) for k in ("status","accepted","global_certificate",
            "feasible_lower_bound","global_bound","relative_gap","construction_seconds",
            "solve_seconds","total_seconds","variable_count","constraint_count","scip_diagnostics",
            "scip_statistics","start_audit","bound_rejection","solver_error","statistics_error","diagnostics_close_error","budget_overrun_seconds","model_profile","requested_relative_gap","certificate_scope","incumbent_source","candidate_error","raw_integrality_error","raw_sos2_residual"))
        row["profile"]=profile;row["repeat"]=repeat;row["commitment"]=commitment
        row["seed_objective"]=initial["objective"];row["seed_sha256"]=bytes2hex(sha256(read(seedfile)))
        row["case_sha256"]=bytes2hex(sha256(read(joinpath(output,"case.json"))))
        row["preparation_seconds_excluded"]=preparation_seconds
        row["profile_warmup_allowance_seconds_excluded"]=repeat==1 ? 10. : 0.
        row["profile_warmup_elapsed_seconds_excluded"]=profile_warm_seconds
        row["common_warmup_elapsed_seconds_excluded"]=common_warm_seconds
        row["capture_enabled"]=capture
        row["parameter_changes"]=Dict(ROOT_PROFILES[profile])
        row["source_sha256"]=bytes2hex(sha256(join(read(p,String) for p in
            sort(filter(p->endswith(p,".jl"),readdir(joinpath(@__DIR__,"..","src");join=true))))))
        row["model_additions"]=Dict(string(k)=>graph[].m.ext[k] for k in (:table_consistency_profile,:network_energy_profile,:product_hull_profile) if graph[]!==nothing && haskey(graph[].m.ext,k))
        row["allowance_seconds"]=seconds
        row["julia_version"]=string(VERSION)
        row["experiment_sha256"]=bytes2hex(sha256(join(read(joinpath(@__DIR__,p),String)
            for p in ("coupled_profile.jl","table_consistency.jl","network_energy.jl","product_hull.jl"))))
        row["cpu_name"]=Sys.CPU_NAME;row["kernel"]=string(Sys.KERNEL)
        row["progress"]=scip_progress(logpath)
        if capture && event[]!==nothing && graph[]!==nothing
            e=event[];b=graph[]
            row["first_root_lp"]=relaxation_summary(b,c,vars[],e.first)
            row["first_root_lp_snapshot_seconds"]=e.first_seconds
            row["first_root_lp_state"]=e.first_state
            row["last_root_lp"]=relaxation_summary(b,c,vars[],isempty(e.last) ? e.first : e.last)
            row["root_snapshot_scope"]="Only optimal root LP solved events outside probing/diving; no postsolve LP reread"
            row["last_root_lp_snapshot_seconds"]=isempty(e.last) ? e.first_seconds : e.last_seconds
            row["last_root_lp_state"]=isempty(e.last) ? e.first_state : e.last_state
            row["root_snapshot_count"]=e.snapshots
            row["root_snapshot_capture_seconds"]=e.capture_seconds
        end
        push!(rows,row)
        writejson(joinpath(output,"summary.json"),rows)
        println(now()," ",profile," ",commitment," ",r["status"]," gap=",r["relative_gap"]);flush(stdout)
        get(r,"accepted",false) && get(get(r,"start_audit",Dict()),"valid",false) || error("unaccepted schedule/start")
        get(r,"solver_error",nothing)===nothing || error("native solve failed")
        get(r,"statistics_error",nothing)===nothing || error("native statistics unavailable")
        get(r,"bound_rejection",nothing)===nothing || error("rejected upper bound")
        catch error
            saved=findlast(r->r["profile"]==profile && r["repeat"]==repeat,rows)
            if saved===nothing
                push!(rows,Dict("profile"=>profile,"repeat"=>repeat,"commitment"=>commitment,
                    "status"=>"EXPERIMENT_ERROR","accepted"=>false,
                    "seed_objective"=>initial["objective"],"allowance_seconds"=>seconds))
                saved=length(rows)
            end
            rows[saved]["experiment_error"]=sprint(showerror,error)
            writejson(joinpath(output,"summary.json"),rows)
            println(now()," ",profile," EXPERIMENT_ERROR ",sprint(showerror,error));flush(stdout)
        end
    end
    any(r->haskey(r,"experiment_error"),rows) && error("one or more experiment profiles failed; records retained")
    rows
end

if abspath(PROGRAM_FILE)==@__FILE__
    root_profile(abspath(ARGS[1]),abspath(ARGS[2]);seconds=parse(Float64,ARGS[3]),
        repeats=parse(Int,ARGS[4]),profiles=split(ARGS[5],','),
        capture=length(ARGS)<6 || ARGS[6]=="true",commitment=length(ARGS)<7 ? "free" : ARGS[7])
end
