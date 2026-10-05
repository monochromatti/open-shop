include(joinpath(@__DIR__,"..","benchmark","diagnose.jl"))

@testset "Native SCIP diagnostic log parser" begin
    mktemp() do path,io
        write(io," time | node | left | dualbound | primalbound | gap\n")
        write(io,"p 0.1s | 1 | 0 | 1.234500e+00 | 1.230000e+00 | 0.37%\n")
        write(io," 0.2s | 2 | 1 | 1.234000e+00 | -- | --\n")
        close(io)
        trace=scip_progress(path)
        @test trace["available"]
        @test length(trace["rows"])==2
        @test trace["rows"][1]["displayed_upper_bound"]≈12345.
        @test trace["rows"][1]["displayed_incumbent_objective"]≈12300.
        @test trace["rows"][2]["displayed_incumbent_objective"]===nothing
        @test trace["rows"][1]["nodes"]==1
    end
end

@testset "Paired benchmark preserves solver failures" begin
    fixture = analytic_global_fixture(:on)
    c = fixture.case
    seed = dispatch_from_controls(c, ones(Int, 1, 1),
        fill(fixture.discharge, 1, 1), zeros(0, 1))
    mktempdir() do directory
        input = joinpath(directory, "case.json")
        writejson(input, case_dict(c))
        seeds = joinpath(directory, "seeds", "case")
        mkpath(seeds)
        benchmark_freeze(joinpath(seeds, "seed.json"), seed)
        output = joinpath(directory, "output")
        records = paired_benchmark([input]; output,
            seed_directory = dirname(seeds), formulations = (:invalid,),
            commitments = (:free,), warmup = false, probe_time_limit = 0.0)
        @test length(records) == 1
        @test records[1]["status"] == "BENCHMARK_ERROR"
        @test occursin("formulation must", records[1]["error"])
    end
end
