using OpenSHOP, Dates

input =
    isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "tokke_vinje", "generated") :
    abspath(ARGS[1])
output =
    length(ARGS) < 2 ? joinpath(@__DIR__, "..", "results", "tokke_vinje") : abspath(ARGS[2])
operational_margin = length(ARGS) < 3 ? 0.1 : parse(Float64, ARGS[3])
mkpath(output)
records = Any[]

for hours in (2, 6, 24)
    case = readcase(joinpath(input, "tokke_vinje_$(hours)h.json"))
    OpenSHOP.validate_inputs(case)
    println(now(), " preparing ", case.name)
    flush(stdout)
    preparation = schedule_case(
        case;
        proposal_time_limit = 5.0,
        nlp_time_limit = 20.0,
        max_refinements = 0,
        operational_margin,
    )
    writejson(joinpath(output, "$(hours)h-preparation.json"), preparation)
    initial = preparation["accepted"] ? preparation["solution"] : nothing
    # Compile this model family and its audit paths outside scored runs.
    warmup = solve(case; initial, time_limit = 3.0, relative_gap = 1e-3)
    writejson(joinpath(output, "$(hours)h-warmup.json"), warmup)
    limit = hours == 2 ? 60.0 : 120.0
    println(now(), " native SCIP ", case.name, " allowance=", limit)
    flush(stdout)
    result = solve(case; initial, time_limit = limit, relative_gap = 1e-3)
    writejson(joinpath(output, "$(hours)h-solve.json"), result)
    row = Dict(
        k => get(result, k, nothing) for k in (
            "case",
            "status",
            "accepted",
            "global_certificate",
            "objective",
            "feasible_lower_bound",
            "global_bound",
            "relative_gap",
            "absolute_gap",
            "total_seconds",
            "solve_seconds",
            "construction_seconds",
            "budget_overrun_seconds",
            "variable_count",
            "constraint_count",
            "incumbent_source",
            "start_audit",
            "candidate_error",
            "solver_error",
        )
    )
    row["preparation_operational_margin_MW"] = operational_margin
    row["preparation_seconds"] = preparation["seconds"]
    row["preparation_accepted"] = preparation["accepted"]
    row["global_allowance_seconds"] = limit
    row["target_gap"] = 1e-3
    row["unit_count"] = length(case.system.generators)
    row["reservoir_count"] = length(case.system.reservoirs)
    row["tunnel_count"] = length(case.system.tunnels)
    row["river_count"] = length(case.system.rivers)
    row["replay_errors"] = get(get(result, "replay_audit", Dict()), "errors", String[])
    push!(records, row)
    writejson(joinpath(output, "runs.json"), records)
    println(row)
    flush(stdout)
end
