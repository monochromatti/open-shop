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

@testset "Local reconstruction corrects controls and recomputes physical values" begin
    c=warm_dispatch_fixture()
    u=zeros(Int, 2, 2)
    physical=dispatch_from_controls(c, u, zeros(2, 2), ones(1, 2))
    raw=deepcopy(physical)
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
    OpenSHOP._repair_dispatch!(c, raw)
    @test raw["forward_reconstructed"]
    @test raw["status"]=="TIME_LIMIT"
    @test raw["raw_solver_validation"]==original["validation"]
    @test raw["raw_solver_objective"]==123456.0
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
        OpenSHOP._repair_dispatch!(c, invalid)
        @test !get(invalid, "forward_reconstructed", false)
        @test haskey(invalid, "forward_reconstruction_error")
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
        @test compatible===initial
        @test reason["used"]
        for grid in ([0.0, 0.5, 2.0], [0.0, 0.5, 1.0, 1.5, 2.0])
            wrong=deepcopy(initial)
            wrong["grid"]=grid
            discarded, reason=OpenSHOP._dispatch_start(c, u, wrong)
            @test discarded===nothing
            @test occursin("grid", reason["ignored_reason"])
        end
        legacy=deepcopy(initial)
        delete!(legacy, "grid")
        @test first(OpenSHOP._dispatch_start(c, u, legacy))===legacy
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
        wrong=deepcopy(initial)
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
    OpenSHOP._repair_dispatch!(c, x)
    @test !get(x, "forward_reconstructed", false)
    @test !x["forward_reconstruction_audit"]["valid"]
    @test x["gate"]==zeros(1, 2)
end

@testset "Verified dispatch uses a fresh original audit" begin
    c=warm_dispatch_fixture()
    u=zeros(Int, 2, 2)
    physical=dispatch_from_controls(c, u, zeros(2, 2), ones(1, 2))
    @test physical["validation"]["valid"]
    for key in ("power", "V")
        bad=deepcopy(physical)
        bad[key][end]+=1.0
        @test bad["validation"]["valid"]
        result=solve_verified(c; u, initial=bad, max_refinements=0)
        @test !result["accepted"]
        @test !result["attempts"][1]["model_valid"]
        @test !isempty(result["attempts"][1]["model_errors"])
        @test !haskey(result["audit"], "replay_grid")
    end
    for key in ("power", "gate")
        bad=deepcopy(physical)
        bad[key]=zeros(1, 3)
        result=solve_verified(c; u, initial=bad, max_refinements=0)
        @test !result["accepted"]
        @test "missing or malformed $key" in result["attempts"][1]["model_errors"]
    end
    valid=deepcopy(physical)
    valid["validation"]=Dict("valid"=>false, "errors"=>["stale diagnostic"])
    valid["raw_solver_validation"]=Dict("valid"=>false)
    result=solve_verified(c; u, initial=valid, max_refinements=0)
    @test result["accepted"]
    @test result["audit"]["original_audit"]["valid"]
    @test isempty(result["attempts"][1]["model_errors"])
    @test haskey(result["audit"], "replay_grid")
    wrong_grid=deepcopy(physical)
    wrong_grid["grid"]=[0.0, 0.5, 2.0]
    @test_throws ArgumentError solve_verified(c; u, initial=wrong_grid, max_refinements=0)
end
