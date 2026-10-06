# Representation screening; losing graphs are removed before production selection.
include("diagnose.jl")
length(ARGS)>=3 || error("usage: contest.jl CASE OUTPUT screen|synthetic [SECONDS]")
input,output=abspath.(ARGS[1:2])
profile=ARGS[3]
variants=(:baseline,:cartesian_ranges,:cartesian_cuts,:cartesian_refined,
    :tensor,:tensor_pruned,:tensor_quadratic,:tensor_refined)
seconds=length(ARGS)>=4 ? parse(Float64,ARGS[4]) : (profile=="synthetic" ? 30.0 : 60.0)
allrows=Any[]
function publish_summary(rows)
    for row in rows
        pop!(row,"progress",nothing);pop!(row,"input_path",nothing)
        row["cpu_name"]=Sys.CPU_NAME;row["kernel"]=string(Sys.KERNEL)
    end
    append!(allrows,rows)
    writejson(joinpath(output,"summary.json"),allrows)
end
if profile=="screen"
    publish_summary(paired_benchmark([input];output=joinpath(output,"screen"),time_limit=seconds,
        formulations=variants,commitments=(:free,:fixed),diagnostics=false))
elseif profile=="synthetic"
    for target in (1e-3,1e-6), native_start in (true,false)
        folder=joinpath(output,"gap-$(target)-start-$(native_start)")
        publish_summary(paired_benchmark([input];output=folder,time_limit=seconds,repeats=3,
            formulations=variants,commitments=(:free,),diagnostics=false,relative_gap=target,
            native_start,probe_time_limit=0.))
    end
else
    error("unknown screening profile")
end
