include("../scripts/fuzz_worker.jl")

@testset "Synthetic campaign inputs have reconstructible feasible references" begin
    for (seed,family) in enumerate(FuzzWorker.FuzzCases.FAMILIES)
        @testset "$family" begin
            data = FuzzWorker.FuzzCases.generate_case(seed,family,"small")
            @test OpenSHOP.validate_inputs(data.case)
            reference = FuzzWorker.reconstruct(data)
            @test validate(data.case,reference)["valid"]
            c = FuzzWorker.FuzzCases.restrict_case(seed,data.case,reference)
            @test OpenSHOP.validate_inputs(c)
            data = merge(data,(;case=c))
            reference = FuzzWorker.reconstruct(data)
            @test validate(c,reference)["valid"]
            trajectory = simulate(c,reference.generator_q,reference.gate;
                commitment=reference.commitment,tunnel_opening=reference.tunnel_opening)
            @test FuzzWorker.water_error(c,trajectory)<1e-8
            trajectory.storage[1,end] += 0.001
            @test FuzzWorker.water_error(c,trajectory)>0.0009
            # The saved physical input, not a particular RNG implementation,
            # reproduces the case and its original operating controls.
            mktempdir() do directory
                path = joinpath(directory,"input.json")
                writejson(path,FuzzWorker.input_dict(data))
                restored = FuzzWorker.load_input(path)
                y = FuzzWorker.reconstruct(restored)
                @test validate(restored.case,y)["valid"]
                @test y.storage≈reference.storage atol=1e-10
                @test y.objective≈reference.objective atol=1e-6
            end
        end
    end
end

@testset "Campaign distinguishes false claims from unfinished search" begin
    data = FuzzWorker.FuzzCases.generate_case(1,"cascade","small")
    reference = FuzzWorker.reconstruct(data)
    options = Dict("replay_factor"=>2,"relative_gap"=>0.01)
    report() = Dict{String,Any}("findings"=>Any[])
    unfinished = SolveResult(Dict("status"=>"TIME_LIMIT"))
    r = report()
    FuzzWorker.check_result!(r,data.case,unfinished,"cold",reference,options)
    @test all(f["severity"]=="warning" for f in r["findings"])
    @test any(f["kind"]=="search_inconclusive" for f in r["findings"])
    for result in (SolveResult(Dict("status"=>"INFEASIBLE")),
        SolveResult(Dict("status"=>"TIME_LIMIT","global_bound"=>reference.objective-100.0)),
        SolveResult(Dict("status"=>"OPTIMAL","global_certificate"=>true)))
        r = report()
        FuzzWorker.check_result!(r,data.case,result,"warm",reference,options)
        @test any(f["severity"]=="error" for f in r["findings"])
    end
    # A free-commitment search can outperform the fixed-control local problem.
    # Running the stages in the opposite order must not call that a regression.
    r = report()
    r["stages"] = Dict("warm"=>Dict("objective"=>reference.objective+100,
        "checks"=>Dict("audit"=>Dict("valid"=>true))))
    local_result = SolveResult(Dict("status"=>"LOCALLY_SOLVED","accepted"=>true,
        "solution"=>reference))
    FuzzWorker.check_result!(r,data.case,local_result,"local",reference,options)
    @test !any(f["kind"]=="local_objective_regression" for f in r["findings"])
end

@testset "Campaign ledger integrates inflow changes inside a control interval" begin
    c = FuzzWorker.FuzzCases.generate_case(1,"cascade","small").case
    name = first(c.system.reservoirs).name
    input = OperationalSeries(object=name,attribute=:inflow,
        times=[0.25,0.75,1.5],values=[2.0,6.0,10.0])
    c = FuzzWorker.FuzzCases.updated(c;operations=[input])
    @test FuzzWorker.input_integral(c,name,:inflow,0.0,2.0,4.0)≈11.5
    @test FuzzWorker.input_integral(c,name,:inflow,0.5,1.0,4.0)≈2.0
end

@testset "Standalone worker preserves a full integer seed" begin
    mktempdir() do directory
        seed = "8010484400485847066"
        request = Dict("seed"=>seed,"family"=>"small_storage","size"=>"small",
            "stages"=>String[],"replay"=>nothing,"replay_factor"=>2)
        path = joinpath(directory,"request.json")
        writejson(path,request)
        result = FuzzWorker.main([path])
        @test result["seed"]==seed
        @test !any(f->f["severity"]=="error",result["findings"])
        saved = FuzzWorker.load_input(joinpath(directory,"input.json"))
        expected = FuzzWorker.FuzzCases.generate_case(parse(Int,seed),"small_storage","small")
        @test saved.q==expected.q
    end
end
