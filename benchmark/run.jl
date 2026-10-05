using OpenSHOP, Dates

output = isempty(ARGS) ? joinpath(@__DIR__, "..", "results", "synthetic") : abspath(ARGS[1])
mkpath(output)
records = Any[]
for name in ("distributed-rivers", "turbine-tables")
    case = readcase(joinpath(@__DIR__, "cases", name * ".json"))
    # Warm the complete solve/audit path before measuring fresh instances.
    warmup = solve(case; time_limit = 60.0, relative_gap = 1e-3)
    writejson(joinpath(output, name * "-warmup.json"), warmup)
    for repetition in 1:3
        GC.gc()
        result = solve(case; time_limit = 30.0, relative_gap = 1e-3)
        writejson(joinpath(output, name * "-$(repetition).json"), result)
        record = Dict(
            k => get(result, k, nothing) for k in (
                "case",
                "status",
                "formulation",
                "accepted",
                "global_certificate",
                "objective",
                "global_bound",
                "relative_gap",
                "total_seconds",
                "construction_seconds",
                "solve_seconds",
                "variable_count",
                "constraint_count",
            )
        )
        record["repetition"] = repetition
        push!(records, record)
        writejson(joinpath(output, "runs.json"), records)
        println(now(), " ", record)
        flush(stdout)
    end
end
