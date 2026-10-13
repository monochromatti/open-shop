@testset "Native SCIP JSON statistics" begin
    m=Model(OpenSHOP.SCIP.Optimizer)
    set_silent(m)
    @variable(m,0<=x<=1)
    @objective(m,Max,x)
    optimize!(m)
    mktempdir() do directory
        stats=OpenSHOP._write_scip_statistics(m,joinpath(directory,"statistics.json"))
        @test stats["status"]["status"]=="optimal solution found"
        @test stats["timing"]["solving_time"]>=0
        @test haskey(stats["propagator"]["plugins"],"obbt")
        diagnostics=OpenSHOP._scip_diagnostics(m)
        @test diagnostics["available"]
        @test diagnostics["presolving_seconds"]>=0
        @test diagnostics["final_upper_bound"]≈10000.
    end
end
