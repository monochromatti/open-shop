# Repeat the production solver on frozen inputs and one audited initial schedule.
include("diagnose.jl")
paths=isempty(ARGS) ? [joinpath(@__DIR__,"cases",name*".json") for name in
    ("distributed-rivers","turbine-tables")] : [abspath(ARGS[1])]
output=length(ARGS)<2 ? joinpath(@__DIR__,"..","results","benchmarks") : abspath(ARGS[2])
seconds=length(ARGS)<3 ? 60.0 : parse(Float64,ARGS[3])
repeats=length(ARGS)<4 ? 3 : parse(Int,ARGS[4])
rows=benchmark_cases(paths;output,time_limit=seconds,repeats,diagnostics=false)
# Scalar summaries can be shared without bundling imported upstream data.
for row in rows
    pop!(row,"input_path",nothing);pop!(row,"progress",nothing)
    row["cpu_name"]=Sys.CPU_NAME;row["kernel"]=string(Sys.KERNEL)
end
writejson(joinpath(output,"summary.json"),rows)
all(row->get(row,"status","")!="BENCHMARK_ERROR" && get(row,"accepted",false) &&
    get(get(row,"start_audit",Dict()),"valid",false),rows) ||
    error("benchmark lacks an accepted schedule or complete initial lift; inspect summary.json")
