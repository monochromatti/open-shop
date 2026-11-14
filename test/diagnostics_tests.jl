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
        diagnostics=OpenSHOP._scip_statistics(m)
        @test diagnostics["available"]
        @test diagnostics["presolving_seconds"]>=0
        @test diagnostics["final_upper_bound"]≈10000.
    end
end

@testset "Logging retains solver bounds and complete SCIP statistics" begin
    fixture=analytic_global_fixture(:on)
    c=fixture.case
    initial=dispatch_from_controls(c,ones(Int,1,1),fill(fixture.discharge,1,1),zeros(0,1))
    mktempdir() do directory
        path=joinpath(directory,"solve.log")
        result=solve(c;initial,time_limit=30.0,preparation_time_limit=0.0,
            diagnostics_path=path)
        @test result.accepted
        @test validate(c,result.solution)["valid"]
        @test result.objective≈fixture.objective atol=1e-3
        @test result.global_certificate
        @test result["solver_statistics"]["available"]
        @test isfinite(result["solver_statistics"]["final_upper_bound"])
        @test result["scip_statistics"]["timing"]["solving_time"]>=0
        saved=OpenSHOP.JSON3.read(read(path*".statistics.json",String),Dict{String,Any})
        @test saved["status"]==result["scip_statistics"]["status"]
    end
end
