# Hosted comparisons keep machine sleep out of bounded solver measurements.
include("diagnose.jl")
length(ARGS)==3 || error("usage: hosted_tables.jl CASE OUTPUT SECONDS")
input,output=abspath.(ARGS[1:2])
seconds=parse(Float64,ARGS[3])
rows=paired_benchmark([input];output,time_limit=seconds,repeats=1,
    formulations=(:baseline,:tensor),commitments=(:free,:fixed),diagnostics=false)
# Publish only scalar diagnostics, never imported cases or generated schedules.
for row in rows
    pop!(row,"progress",nothing)
    pop!(row,"input_path",nothing)
    row["cpu_name"]=Sys.CPU_NAME
    row["kernel"]=string(Sys.KERNEL)
end
writejson(joinpath(output,"summary.json"),rows)
all(row->get(row,"status","")!="BENCHMARK_ERROR" &&
    get(row,"accepted",false) && get(get(row,"start_audit",Dict()),"valid",false),rows) ||
    error("hosted comparison lacks an audited common start or accepted schedule; inspect scalar summary")
