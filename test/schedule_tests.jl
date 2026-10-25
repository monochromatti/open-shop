function warm_dispatch_fixture(; rivers = true, minimum_power = 0.1)
    lake=Reservoir(name = :WarmLake, z0 = 100.0, slope = 1.0, curvature = 0.05,
        v0 = 2.0, vmin = 0.1, vmax = 4.0, inflow = 8.0, water_value = 10.0)
    generators=[Generator(name = name, plant = :WarmPlant, qmin = 1.0, qmax = 5.0,
        pmin = minimum_power, pmax = 5.0, efficiency = 0.9, qbest = 3.0,
        qcurvature = 0.0, hbest = 100.0, hcurvature = 0.0, hmin = 70.0,
        hmax = 130.0, minup = 0.0, mindown = 0.0, startup = 0.0) for name in (:A, :B)]
    reaches=rivers ? [River(name = :WarmReach, source = :WarmLake, target = :WarmSea,
        curves = RiverRouting.DelayCurve[], deterministic_delay = 0.25,
        capacity = 5.0, law = :orifice, coefficient = 0.1, water_value = 10.0)] : River[]
    system=HydroSystem(reservoirs = [lake], junctions = [Junction(name = :WarmIntake)],
        boundaries = [Boundary(name = :WarmSea, head = 0.0)],
        tunnels = [Tunnel(name = :WarmTunnel, source = :WarmLake, target = :WarmIntake,
            resistance = 0.2, capacity = 20.0)],
        plants = [Plant(name = :WarmPlant, source = :WarmIntake, target = :WarmSea,
            pmax = 10.0, tailwater_curve = TableCurve([0.0, 5.0, 10.0], [0.0, 0.5, 2.0]))],
        generators = generators, rivers = reaches)
    operations=[OperationalSeries(object = :B, attribute = :qmin, times = [0.0, 1.0],
        values = [2.5, 1.5]), OperationalSeries(object = :B, attribute = :qmax,
        times = [0.0, 1.0], values = [3.0, 2.0])]
    ScheduleCase(name = "changed commitment hydraulics", system = system, grid = [0.0, 1.0, 2.0],
        prices = [50.0, 90.0], operations = operations)
end

@testset "Preparation improves an accepted schedule and preserves its feasibility" begin
    c=warm_dispatch_fixture(;rivers=false)
    off=dispatch_from_controls(c,zeros(Int,2,2),zeros(2,2),zeros(0,2))
    @test replay_audit(c,off)["valid"]
    prepared=schedule_case(c;initial=off,time_limit=15.0,proposal_time_limit=2.0,
        nlp_time_limit=5.0,max_refinements=0,refine_grid=false)
    @test prepared["accepted"]
    if prepared["accepted"]
        @test prepared["solution"]["objective"]>off.objective+1.0
        @test validate(c,prepared["solution"])["valid"]
        @test replay_audit(c,prepared["solution"])["valid"]
    end
    damaged=OpenSHOP.schedule_dict(off)
    damaged["power"][1]=1.0
    @test_throws ArgumentError schedule_case(c;initial=damaged,time_limit=15.0)
end

@testset "Shared preparation deadline" begin
    c=warm_dispatch_fixture()
    u=zeros(Int,2,2)
    expired=solve_verified(c;u,deadline=time()-1,max_refinements=2)
    @test !expired["accepted"]
    @test expired["solution"]===nothing
    @test expired["case"].grid==c.grid

    # Auditing a supplied schedule remains necessary even when no optimization
    # time remains. A deadline must never turn stale acceptance into evidence.
    x=dispatch_from_controls(c,u,zeros(2,2),ones(1,2))
    audited=solve_verified(c;u,initial=x,deadline=time()-1,max_refinements=0)
    @test audited["accepted"]
    damaged=OpenSHOP.schedule_dict(x)
    damaged["power"][1]=1.0
    rejected=solve_verified(c;u,initial=damaged,deadline=time()-1,max_refinements=0)
    @test !rejected["accepted"]

    prepared=schedule_case(c;time_limit=1e-9)
    @test !prepared["accepted"]
    @test prepared["solution"]===nothing
    for invalid in (0.0,-1.0,NaN)
        @test_throws ArgumentError schedule_case(c;time_limit=invalid)
    end
    @test_throws ArgumentError solve_verified(c;deadline=NaN)
end

@testset "Replay margins preserve the requested decision grid" begin
    base=warm_dispatch_fixture(;rivers=false)
    lake=OpenSHOP._river_replace(only(base.system.reservoirs);inflow=80.0)
    units=[OpenSHOP._river_replace(g;pmax=2.0) for g in base.system.generators]
    system=OpenSHOP._river_replace(base.system;reservoirs=[lake],generators=units)
    c=OpenSHOP._river_replace(base;system)
    u=[1 1;0 0]
    verified=solve_verified(c;u,max_refinements=1,refine_grid=false,time_limit=10.0)
    @test verified["accepted"]
    @test verified["case"].grid==c.grid
    @test verified["solution"]["u"]==u
    @test replay_audit(c,verified["solution"])["valid"]
end

@testset "Local reconstruction corrects controls and recomputes physical values" begin
    c=warm_dispatch_fixture()
    u=zeros(Int, 2, 2)
    physical=dispatch_from_controls(c, u, zeros(2, 2), ones(1, 2))
    raw=OpenSHOP.schedule_dict(physical)
    raw["status"]="TIME_LIMIT"
    raw["generator_q"].=1e-7
    raw["gate"].+=1e-14
    raw["objective"]=123456.0
    raw["power"].=100.0
    raw["validation"]=validate(c, raw)
    @test !raw["validation"]["valid"]
    original=deepcopy(raw)
    rebuilt, correction=OpenSHOP._reconstruct_candidate(c, raw)
    @test raw==original
    @test correction.flow==1e-7
    @test correction.gate>0
    @test rebuilt["validation"]["valid"]
    report=Dict{String,Any}("status"=>raw["status"])
    raw=OpenSHOP._repair_dispatch!(c, OpenSHOP._schedule_input(raw;case=c), report)
    @test report["forward_reconstructed"]
    @test report["status"]=="TIME_LIMIT"
    @test report["raw_solver_validation"]==original["validation"]
    @test report["raw_solver_objective"]==123456.0
    @test raw["generator_q"]==zeros(2, 2)
    @test raw["gate"]==ones(1, 2)
    @test raw["power"]==zeros(2, 2)
    @test raw["V"]≈physical["V"] atol=1e-10
    @test raw["objective"]≈physical["objective"] atol=1e-8
    @test replay_audit(c, raw)["valid"]
    for (key, value) in (("generator_q", 1e-3), ("generator_q", -1e-3),
            ("gate", 1.0001), ("gate", NaN))
        invalid=deepcopy(original)
        invalid[key][1]=value
        @test_throws ArgumentError OpenSHOP._reconstruct_candidate(c, invalid)
        report=Dict{String,Any}()
        invalid=OpenSHOP._repair_dispatch!(c, OpenSHOP._schedule_input(invalid;case=c), report)
        @test !get(report, "forward_reconstructed", false)
        @test haskey(report, "forward_reconstruction_error")
        @test invalid["objective"]==original["objective"]
    end
end

@testset "Warm starts require compatible grids and finite arrays" begin
    for rivers in (false, true)
        c=warm_dispatch_fixture(; rivers)
        u=[1 1; 0 0]
        initial=dispatch_from_controls(c, u, [4.0 4.0; 0.0 0.0],
            fill(0.8, length(c.system.rivers), 2))
        compatible, reason=OpenSHOP._dispatch_start(c, u, initial)
        @test compatible.generator_q==initial.generator_q
        @test compatible.objective==initial.objective
        @test reason["used"]
        for grid in ([0.0, 0.5, 2.0], [0.0, 0.5, 1.0, 1.5, 2.0])
            wrong=OpenSHOP.schedule_dict(initial)
            wrong["grid"]=grid
            discarded, reason=OpenSHOP._dispatch_start(c, u, wrong)
            @test discarded===nothing
            @test occursin("grid", reason["ignored_reason"])
        end
        legacy=OpenSHOP.schedule_dict(initial)
        delete!(legacy, "grid")
        @test first(OpenSHOP._dispatch_start(c,u,legacy)) isa Schedule
        for key in ("generator_q", "power", "V", "gate")
            for damage in (:shape, :missing)
                wrong=deepcopy(legacy)
                damage==:shape ?
                    (wrong[key]=zeros(size(wrong[key], 1), size(wrong[key], 2)+1)) :
                    delete!(wrong, key)
                discarded, reason=OpenSHOP._dispatch_start(c, u, wrong)
                @test discarded===nothing
                @test !reason["used"]
            end
        end
        wrong=OpenSHOP.schedule_dict(initial)
        wrong["H"][1]=NaN
        @test first(OpenSHOP._dispatch_start(c, u, wrong))===nothing
        @test first(OpenSHOP._dispatch_start(c, u, nothing))===nothing
    end
end

@testset "Passive river semantics survive guarded reconstruction" begin
    c=warm_dispatch_fixture()
    weir=OpenSHOP._river_replace(only(c.system.rivers); law = :weir, coefficient = 0.001)
    c=OpenSHOP._river_replace(c; system = OpenSHOP._river_replace(c.system; rivers = [weir]))
    x=dispatch_from_controls(c, zeros(Int, 2, 2), zeros(2, 2), zeros(1, 2))
    @test !x["validation"]["valid"]
    report=Dict{String,Any}()
    repaired=OpenSHOP._repair_dispatch!(c, x, report)
    @test repaired.generator_q==x.generator_q
    @test repaired.objective==x.objective
    @test !get(report, "forward_reconstructed", false)
    @test !report["forward_reconstruction_audit"]["valid"]
    @test x["gate"]==zeros(1, 2)
end

@testset "Verified dispatch uses a fresh original audit" begin
    c=warm_dispatch_fixture()
    u=zeros(Int, 2, 2)
    physical=dispatch_from_controls(c, u, zeros(2, 2), ones(1, 2))
    @test physical["validation"]["valid"]
    for key in ("power", "V")
        bad=OpenSHOP.schedule_dict(physical)
        bad[key][end]+=1.0
        @test bad["validation"]["valid"]
        result=solve_verified(c; u, initial=bad, max_refinements=0)
        @test !result["accepted"]
        @test !result["attempts"][1]["model_valid"]
        @test !isempty(result["attempts"][1]["model_errors"])
        @test !haskey(result["audit"], "replay_grid")
    end
    for key in ("power", "gate")
        bad=OpenSHOP.schedule_dict(physical)
        bad[key]=zeros(1, 3)
        result=solve_verified(c; u, initial=bad, max_refinements=0)
        @test !result["accepted"]
        @test "missing or malformed $key" in result["attempts"][1]["model_errors"]
    end
    valid=OpenSHOP.schedule_dict(physical)
    valid["validation"]=Dict("valid"=>false, "errors"=>["stale diagnostic"])
    valid["raw_solver_validation"]=Dict("valid"=>false)
    result=solve_verified(c; u, initial=valid, max_refinements=0)
    @test result["accepted"]
    @test result["audit"]["original_audit"]["valid"]
    @test isempty(result["attempts"][1]["model_errors"])
    @test haskey(result["audit"], "replay_grid")
    wrong_grid=OpenSHOP.schedule_dict(physical)
    wrong_grid["grid"]=[0.0, 0.5, 2.0]
    @test_throws ArgumentError solve_verified(c; u, initial=wrong_grid, max_refinements=0)
end
