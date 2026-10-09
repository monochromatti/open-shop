# Profile fixed-commitment dispatch with identical starting values.
using OpenSHOP, JuMP, Ipopt, JSON3, LinearAlgebra, SHA, TOML, Libdl
const MOI=OpenSHOP.MOI

const HELP="""
Usage: julia --project=. benchmark/linear-systems.jl FIXTURE.json OPTIONS.toml OUTPUT [REPEATS]

The fixture contains a case dictionary, a fixed binary commitment u and a complete
warm schedule. Freeze it once for all comparisons. OPTIONS is a dictionary of
Ipopt attributes (an empty file uses OpenSHOP defaults). Two warmup calls are excluded;
REPEATS defaults to 3. Native timing logs, callbacks and fresh physical/replay
audits are saved under OUTPUT. Runs are serial with one BLAS thread.

Set OpenMP variables before Julia starts when testing SPRAL. CPU/wall allowances
and iteration caps are diagnostic settings; they do not alter physical tolerances.
"""

function matrix_rows(x,columns)
    isempty(x) && return zeros(0,columns)
    all(length(row)==columns for row in x) || error("Fixture matrix has incompatible columns")
    permutedims(hcat([Float64.(row) for row in x]...))
end
function load_fixture(path)
    data=JSON3.read(read(path,String),Dict{String,Any})
    c=OpenSHOP.case_from_dict(data["case"])
    T=length(c.prices)
    u=matrix_rows(data["u"],T)
    all(x->x in (0,1),u) || error("Fixture commitment must be binary")
    OpenSHOP.admissible(c,u) || error("Fixture commitment is inadmissible")
    warm=data["warm"]
    for key in ("V","H","tunnel_q","generator_q","power","river_release","gate","arrival_volume","u","shortfall_release")
        haskey(warm,key) && (warm[key]=matrix_rows(warm[key],T+(key=="V")))
    end
    _,status=OpenSHOP._dispatch_start(c,u,warm)
    get(status,"used",false) || error("Fixture warm start is incompatible: $status")
    replay_audit(c,warm)["valid"] || error("Fixture reference failed physical/replay acceptance")
    c,Int.(u),warm
end
function native_timers(log)
    Dict(m.captures[1]=>parse(Float64,m.captures[2]) for m in eachmatch(
        r"(?m)^\s*(\w+)\.+:.*?wall:\s*([0-9.]+)",log))
end
function source_hash()
    root=dirname(dirname(pathof(OpenSHOP)))
    files=sort!([joinpath(dir,name) for (dir,_,names) in walkdir(joinpath(root,"src"))
                 for name in names if endswith(name,".jl")])
    bytes2hex(sha256(join((relpath(p,root)*"\0"*read(p,String) for p in files),"\0")))
end
function measure(c,u,warm,attributes,folder,index;warmup=false)
    began=time()
    b=OpenSHOP._build_dispatch(c;u,warm,time_limit=30.0,operational_margin=0.1)
    OpenSHOP._local_optimizer!(b.m,"Ipopt",30.0)
    for (key,value) in attributes
        set_optimizer_attribute(b.m,key,value)
    end
    # Match Ipopt.jl's zero-within-bounds fallback before its interior push.
    starts=map(all_variables(b.m)) do v
        value=start_value(v)
        value!==nothing && return value
        lo=is_fixed(v) ? fix_value(v) : has_lower_bound(v) ? lower_bound(v) : -Inf
        hi=is_fixed(v) ? fix_value(v) : has_upper_bound(v) ? upper_bound(v) : Inf
        clamp(0.0,lo,hi)
    end
    start_hash=bytes2hex(sha256(reinterpret(UInt8,starts)))
    construction=time()-began
    trace=NamedTuple[]
    sizehint!(trace,1001)
    MOI.set(b.m,Ipopt.CallbackFunction(),function(mode,iteration,objective,primal,dual,mu,step,regularization,alpha_dual,alpha_primal,line_searches)
        push!(trace,(;mode,iteration,objective,primal,dual,mu,step,regularization,
                     alpha_dual,alpha_primal,line_searches))
        true
    end)
    unset_silent(b.m)
    set_optimizer_attribute(b.m,"print_level",4)
    set_optimizer_attribute(b.m,"print_timing_statistics","yes")
    prefix=joinpath(folder,(warmup ? "warmup-" : "repeat-")*lpad(string(index),2,'0'))
    measurement=nothing
    cpu=0.0
    open(prefix*".log","w") do io
        redirect_stdout(io) do
            cpu_start=Float64(ccall(:clock,Clong,()))/1e6
            measurement=@timed optimize!(b.m)
            cpu=Float64(ccall(:clock,Clong,()))/1e6-cpu_start
        end
    end
    result=Dict{String,Any}("case"=>c.name,"repeat"=>index,"warmup"=>warmup,
        "options"=>attributes,"start_sha256"=>start_hash,"variables"=>num_variables(b.m),
        "status"=>string(termination_status(b.m)),"iterations"=>MOI.get(b.m,MOI.BarrierIterations()),
        "construction_seconds"=>construction,"optimize_seconds"=>measurement.time,
        "process_cpu_seconds"=>cpu,"julia_allocated_bytes"=>measurement.bytes,"gc_seconds"=>measurement.gctime,
        "compilation_seconds"=>measurement.compile_time,"trace"=>trace,"accepted"=>false)
    result["native_timers"]=native_timers(read(prefix*".log",String))
    acceptance_started=time()
    if has_values(b.m)
        x=OpenSHOP._dispatch_values(b)
        x["u"]=u
        x["grid"]=copy(c.grid)
        x["validation"]=validate(c,x;transport=b.transport)
        result["raw_validation"]=x["validation"]
        result["raw_objective"]=x["objective"]
        x["validation"]["valid"] || OpenSHOP._repair_dispatch!(c,x;transport=b.transport)
        result["solution"]=x
        if x["validation"]["valid"]
            audit=replay_audit(c,x)
            result["audit"]=audit
            result["accepted"]=audit["valid"]
            audit["valid"] && (result["accepted_objective"]=x["objective"])
        end
    end
    result["acceptance_seconds"]=time()-acceptance_started
    result["total_seconds"]=time()-began
    writejson(prefix*".json",result)
    println(warmup ? "warmup" : "repeat",' ',index,' ',result["status"]," accepted=",
            result["accepted"]," seconds=",round(measurement.time;digits=3))
    flush(stdout)
    result
end
function main(args)
    if isempty(args) || args[1] in ("--help","-h")
        print(HELP)
        return
    end
    length(args) in (3,4) || error(HELP)
    fixture,optionfile,folder=abspath.(args[1:3])
    repeats=length(args)==4 ? parse(Int,args[4]) : 3
    repeats>=1 || error("Positive repetition count required")
    BLAS.set_num_threads(1)
    c,u,warm=load_fixture(fixture)
    attributes=TOML.parsefile(optionfile)
    isdir(folder) && !isempty(readdir(folder)) && error("Use a fresh output directory")
    mkpath(folder)
    project=Base.active_project()
    manifest=joinpath(dirname(project),"Manifest.toml")
    native=Libdl.dlpath(Ipopt.libipopt)
    writejson(joinpath(folder,"configuration.json"),Dict(
        "fixture_sha256"=>bytes2hex(sha256(read(fixture))),"options"=>attributes,
        "base_settings"=>Dict("tol"=>1e-8,"bound_relax_factor"=>0.0,"max_iter"=>1000,
            "max_cpu_time"=>30.0,"max_wall_time"=>30.0,"operational_margin_MW"=>0.1),
        "source_sha256"=>source_hash(),"manifest_sha256"=>isfile(manifest) ? bytes2hex(sha256(read(manifest))) : nothing,
        "runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),"julia_version"=>string(VERSION),
        "package_versions"=>Dict("Ipopt"=>string(Base.pkgversion(Ipopt)),
            "JuMP"=>string(Base.pkgversion(JuMP)),"MathOptInterface"=>string(Base.pkgversion(MOI))),
        "native_ipopt_library"=>native,"native_ipopt_sha256"=>bytes2hex(sha256(read(native))),
        "cpu_name"=>Sys.CPU_NAME,"kernel"=>string(Sys.KERNEL),"machine"=>Sys.MACHINE,
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "blas_configuration"=>sprint(show,BLAS.lbt_get_config()),
        "thread_environment"=>Dict(k=>get(ENV,k,nothing) for k in
            ("OMP_NUM_THREADS","OMP_CANCELLATION","OMP_PROC_BIND","OPENBLAS_NUM_THREADS"))))
    for index in 1:2
        measure(c,u,warm,attributes,folder,index;warmup=true)
    end
    rows=[measure(c,u,warm,attributes,folder,index) for index in 1:repeats]
    all(row->row["accepted"],rows) || error("Rejected dispatch; inspect recorded audits")
end
main(ARGS)
