# Warm chronological replay measurements, independent of optimization timing.
# Run this same file against two revisions to compare the residual implementation.
using OpenSHOP, JuMP, Statistics
include(joinpath(@__DIR__, "..", "examples", "mixed_source_network.jl"))
using .MixedSourceNetworkExample: mixed_source_fixture

function measure_hydraulics(fixture; repetitions = 9)
    c = fixture.case
    simulate(c, fixture.q, fixture.gates)
    measurements = [@timed simulate(c, fixture.q, fixture.gates) for _ in 1:repetitions]
    trajectory = last(measurements).value
    trajectory["converged"] || error("hydraulic replay did not converge")
    local_model = OpenSHOP._build_dispatch(c; u = fixture.u)
    global_model = OpenSHOP._build_global_dispatch(c; u = fixture.u)
    Dict(
        "simulation" => trajectory,
        "seconds" => median(m.time for m in measurements),
        "bytes" => median(m.bytes for m in measurements),
        "local_variables" => num_variables(local_model.m),
        "local_constraints" =>
            num_constraints(local_model.m; count_variable_in_set_constraints = true),
        "global_variables" => num_variables(global_model.m),
        "global_constraints" =>
            num_constraints(global_model.m; count_variable_in_set_constraints = true),
    )
end

function hydraulics_benchmark()
    Dict(
        "deterministic" => measure_hydraulics(mixed_source_fixture()),
        "distributed" => measure_hydraulics(mixed_source_fixture(; distributed = true)),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: julia --project=. benchmark/hydraulics.jl output.json")
    writejson(only(ARGS), hydraulics_benchmark())
end
