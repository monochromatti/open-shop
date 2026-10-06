# Exact table-graph comparison with one frozen case and common audited controls.
include("diagnose.jl")
length(ARGS)>=2 || error("usage: tables.jl CASE OUTPUT [SEED_DIRECTORY] [SECONDS] [REPEATS] [trace|silent]")
input=abspath(ARGS[1])
output=abspath(ARGS[2])
seed_directory=length(ARGS)>=3 && ARGS[3]!="-" ? abspath(ARGS[3]) : nothing
seconds=length(ARGS)>=4 ? parse(Float64,ARGS[4]) : 60.0
repeats=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 2
mode=length(ARGS)>=6 ? ARGS[6] : "trace"
mode in ("trace","silent") || error("logging mode must be trace or silent")
paired_benchmark([input];output,seed_directory,time_limit=seconds,repeats,diagnostics=mode=="trace",
    formulations=(:baseline,:tensor),commitments=(:free,:fixed))
