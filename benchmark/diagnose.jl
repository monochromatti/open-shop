using OpenSHOP, JSON3, SHA, Dates

"""Parse displayed SCIP progress rows; rounded log values never certify a gap."""
function scip_progress(path; objective_scale=10000.)
    rows=Any[];columns=String[]
    isfile(path) || return Dict("available"=>false,"rows"=>rows)
    function number(x)
        clean=strip(x)
        lowercase(clean) in ("--","-","cutoff","infeasible","unbounded") && return nothing
        value=tryparse(Float64,replace(clean,"%"=>""))
        value===nothing || !isfinite(value) || abs(value)>=1e19 ? nothing : value
    end
    for line in eachline(path)
        parts=strip.(split(line,'|'))
        if any(x->lowercase(x)=="dualbound",parts) && any(x->lowercase(x)=="primalbound",parts)
            columns=lowercase.(parts);continue
        end
        isempty(columns) && continue
        length(parts)==length(columns) || continue
        tm=match(r"([0-9]+(?:\.[0-9]+)?)s",parts[1])
        tm===nothing && continue
        row=Dict{String,Any}("seconds"=>parse(Float64,tm.captures[1]),"line"=>line)
        for (label,key) in (("node","nodes"),("left","nodes_left"),("dualbound","displayed_upper_bound"),("primalbound","displayed_incumbent_objective"),("gap","displayed_gap_percent"))
            index=findfirst(==(label),columns);index===nothing && continue
            value=number(parts[index]);row[key]=value===nothing ? nothing : (label in ("dualbound","primalbound") ? value*objective_scale : value)
        end
        push!(rows,row)
    end
    Dict("available"=>!isempty(rows),"rows"=>rows,"scope"=>"native SCIP display trajectory, rounded values only; final certificate uses solver bound and independent physical audit")
end

# Sorted object keys make frozen case/seed hashes reproducible across processes.
function benchmark_freeze(path,data)
    function canonical(x)
        if x isa AbstractDict
            keys_sorted=sort!(collect(string.(keys(x))))
            NamedTuple{Tuple(Symbol.(keys_sorted))}(Tuple(canonical(x[key]) for key in keys_sorted))
        elseif x isa AbstractVector
            canonical.(x)
        else
            x
        end
    end
    ready=canonical(OpenSHOP.jsonready(data))
    mkpath(dirname(path))
    open(path,"w") do io
        JSON3.pretty(io,ready)
    end
end

function benchmark_seed(path,c)
    isfile(path) || return nothing
    data=JSON3.read(read(path,String))
    data===nothing && return nothing
    get(data,"solution",data)===nothing && return nothing
    data=get(data,"solution",data)
    matrix(key)=isempty(data[key]) ? zeros(0,length(c.prices)) : permutedims(hcat([Float64.(row) for row in data[key]]...))
    fractional=matrix("u")
    all(isinteger,fractional) || error("frozen commitment is not integral")
    u=Int.(fractional);q=matrix("generator_q")
    gate=isempty(data["gate"]) ? zeros(0,length(c.prices)) : matrix("gate")
    dispatch_from_controls(c,u,q,gate)
end

function paired_benchmark(case_paths;output,time_limit=60.,repeats=1,formulations=(:baseline,:tightened),commitments=(:free,:fixed),
        operational_margin=.1,relative_gap=1e-3,seed_directory=nothing,warmup=true,probe_time_limit=20.)
    repeats isa Integer && repeats>=1 || throw(ArgumentError("positive repeats required"))
    mkpath(output);records=Any[]
    source_hash=bytes2hex(sha256(join([read(p,String) for p in sort(filter(p->endswith(p,".jl"),readdir(joinpath(@__DIR__,"..","src");join=true)))])))
    for case_path in case_paths
        label=replace(splitext(basename(case_path))[1],r"[^a-zA-Z0-9_-]"=>"_")
        folder=joinpath(output,label);mkpath(folder)
        try
            c=readcase(case_path)
            frozen=joinpath(folder,"case.json");benchmark_freeze(frozen,case_dict(c));c=readcase(frozen)
            case_hash=bytes2hex(sha256(read(frozen)))
            frozen_seed=joinpath(folder,"seed.json")
            supplied=seed_directory===nothing ? frozen_seed : joinpath(seed_directory,label,"seed.json")
            preparation_seconds=0.;initial=nothing
            if isfile(supplied)
                initial=benchmark_seed(supplied,c)
                initial===nothing || initial["validation"]["valid"] || error("frozen seed fails the same case audit")
                benchmark_freeze(frozen_seed,initial)
            else
                began=time()
                preparation=schedule_case(c;proposal_time_limit=5.,nlp_time_limit=20.,max_refinements=0,operational_margin)
                preparation_seconds=time()-began
                writejson(joinpath(folder,"preparation.json"),preparation)
                initial=preparation["accepted"] ? preparation["solution"] : nothing
                benchmark_freeze(frozen_seed,initial)
            end
            seed_hash=bytes2hex(sha256(read(frozen_seed)))
            control_path=joinpath(folder,"seed-controls.json")
            benchmark_freeze(control_path,initial===nothing ? nothing : Dict(k=>initial[k] for k in ("u","generator_q","gate")))
            control_hash=bytes2hex(sha256(read(control_path)))
            reference_lower=initial===nothing ? nothing : initial["objective"]
            probe_seconds=0.;probe_valid=false
            if initial!==nothing && probe_time_limit>0
                probe_started=time()
                probe=solve_case(c;u=initial["u"],warm=initial,time_limit=probe_time_limit)
                probe_seconds=time()-probe_started
                writejson(joinpath(folder,"discrete-probe.json"),probe)
                probe_valid=get(get(probe,"validation",Dict()),"valid",false)
                probe_valid && (reference_lower=max(reference_lower,probe["objective"]))
            end
            metadata=Dict("case"=>c.name,"case_sha256"=>case_hash,"seed_sha256"=>seed_hash,
                "seed_controls_sha256"=>control_hash,"source_sha256"=>source_hash,"input_path"=>abspath(case_path),
                "diagnostic_trace_enabled"=>true,"warmup_enabled"=>warmup,
                "julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),
                "preparation_seconds_excluded"=>preparation_seconds,"probe_seconds_excluded"=>probe_seconds,
                "probe_discrete_valid"=>probe_valid,"known_discrete_lower_bound"=>reference_lower,
                "seed_objective"=>initial===nothing ? nothing : initial["objective"],
                "seed_shared_across_all_variants"=>true,"target_relative_gap"=>relative_gap,
                "global_allowance_seconds"=>time_limit,
                "scope"=>"same frozen discrete equations and audited seed; fixed commitment restricts feasible set")
            writejson(joinpath(folder,"metadata.json"),metadata)
            for formulation in formulations
                if warmup
                    try
                        result=solve(c;initial,time_limit=2.,relative_gap,formulation)
                        writejson(joinpath(folder,"warmup-$(formulation).json"),result)
                    catch error
                        writejson(joinpath(folder,"warmup-$(formulation).json"),Dict("error"=>sprint(showerror,error)))
                    end
                end
            end
            # Alternate order across repetitions to reduce systematic warm-cache bias.
            variants=[(f,mode) for f in formulations for mode in commitments]
            for repetition in 1:repeats
                for (formulation,mode) in (isodd(repetition) ? variants : reverse(variants))
                    row=merge(copy(metadata),Dict("formulation"=>string(formulation),"commitment"=>string(mode),"repeat"=>repetition))
                    stem="$(formulation)-$(mode)-$(lpad(string(repetition),2,'0'))"
                    if mode==:fixed && initial===nothing
                        row["status"]="SKIPPED_NO_COMMON_FEASIBLE_SEED"
                    else
                        fixed_u=mode==:fixed ? copy(initial["u"]) : nothing
                        log_path=joinpath(folder,stem*".log")
                        began=time()
                        try
                            result=solve(c;initial,fixed_u,time_limit,relative_gap,formulation,diagnostics_path=log_path)
                            row["harness_seconds"]=time()-began
                            upper=get(result,"global_bound",nothing)
                            consistent=upper===nothing || reference_lower===nothing || upper>=reference_lower-1e-4
                            row["bound_consistent_with_known_schedule"]=consistent
                            if !consistent
                                result["rejected_global_bound"]=upper
                                result["bound_rejection"]="solver upper bound below independently validated fixed-commitment probe"
                                result["global_bound"]=nothing
                                result["global_certificate"]=false
                                result["relative_gap"]=nothing
                                result["absolute_gap"]=nothing
                            end
                            for key in ("status","accepted","global_certificate","feasible_lower_bound","global_bound","absolute_gap","relative_gap","construction_seconds","solve_seconds","total_seconds","budget_overrun_seconds","variable_count","constraint_count","model_profile","incumbent_source","start_audit","scip_diagnostics","solver_error","candidate_error","bound_rejection")
                                row[key]=get(result,key,nothing)
                            end
                            row["progress"]=scip_progress(log_path)
                            writejson(joinpath(folder,stem*".json"),result)
                        catch error
                            row["status"]="BENCHMARK_ERROR";row["error"]=sprint(showerror,error)
                            row["harness_seconds"]=time()-began
                        end
                    end
                    push!(records,row);writejson(joinpath(output,"runs.json"),records)
                    println(now()," ",label," ",stem," ",get(row,"status","UNKNOWN"));flush(stdout)
                end
            end
        catch error
            push!(records,Dict("input_path"=>abspath(case_path),"status"=>"CASE_OR_SEED_ERROR","error"=>sprint(showerror,error)))
            writejson(joinpath(output,"runs.json"),records)
        end
    end
    records
end
