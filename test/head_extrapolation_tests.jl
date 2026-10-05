function head_extrapolation_case(policy = :error)
    table=TurbineTable(
        [50.0, 100.0],
        [0.0, 2.0, 4.0],
        [0.6 0.7; 0.7 0.8; 0.8 0.9],
        [1.0, 1.0],
        [4.0, 4.0];
        head_extrapolation = policy,
    )
    lake=Reservoir(
        name = :Lake,
        z0 = 23.0,
        slope = 1.0,
        v0 = 2.0,
        vmin = 1.0,
        vmax = 3.0,
        inflow = 2.0,
        water_value = 10.0,
    )
    unit=Generator(
        name = :Unit,
        plant = :Plant,
        qmin = 1.0,
        qmax = 4.0,
        pmin = 0.05,
        pmax = 1.0,
        hmin = 20.0,
        hmax = 30.0,
        hbest = 25.0,
        qbest = 2.0,
        qcurvature = 0.0,
        hcurvature = 0.0,
        initial_on = 1,
        minup = 0.0,
        mindown = 0.0,
        startup = 0.0,
        turbine_table = table,
    )
    system=HydroSystem(;
        reservoirs = [lake],
        junctions = Junction[],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
        tunnels = Tunnel[],
        plants = [Plant(name = :Plant, source = :Lake, target = :Sea, pmax = 1.0)],
        generators = [unit],
        rivers = River[],
    )
    ScheduleCase(;
        name = "explicit_head_extrapolation",
        system,
        grid = [0.0, 1.0],
        prices = [50.0],
    )
end

@testset "Operating turbine head extrapolation is explicit and serializable" begin
    strict=head_extrapolation_case()
    @test only(strict.system.generators).turbine_table.head_extrapolation==:error
    @test_throws ArgumentError OpenSHOP.validate_inputs(strict)
    @test_throws ArgumentError head_extrapolation_case(:clamp)

    c=head_extrapolation_case(:linear)
    @test OpenSHOP.validate_inputs(c)
    table=only(c.system.generators).turbine_table
    @test table.heads==[50.0, 100.0]
    # Raw table APIs still require an explicit extrapolation request.
    @test_throws DomainError OpenSHOP.turbine_efficiency(table, 2.0, 25.0)
    @test OpenSHOP.turbine_efficiency(table, 2.0, 25.0; extrapolation = :linear)≈0.65 atol=1e-12
    @test OpenSHOP.turbine_efficiency(table, 2.0, 125.0; extrapolation = :linear)≈0.85 atol=1e-12
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    missing=case_dict(c)
    delete!(missing["generators"][1]["turbine_table"], "head_extrapolation")
    @test_throws ArgumentError case_from_dict(missing)

    # Inflow balances discharge: actual net head is exactly25, below both
    # reference heads. Independent efficiency interpolation gives eta=.65.
    x=dispatch_from_controls(c, ones(Int, 1, 1), fill(2.0, 1, 1), zeros(0, 1))
    @test x["validation"]["valid"]
    @test x["H"][1, 1]≈25.0 atol=1e-10
    @test x["power"][1, 1]≈0.00981*2*25*0.65 atol=1e-10
    @test replay_audit(c, x)["valid"]
    @test OpenSHOP._build_dispatch(c).m isa JuMP.Model
    global_model=OpenSHOP._build_global_dispatch(c; joint = true)
    lift=OpenSHOP._lift_start!(global_model, c, x)
    @test lift["valid"]
    @test lift["objective"]≈x["objective"] atol=1e-8
end

@testset "Preparation margin tightens candidates without changing the case" begin
    c = head_extrapolation_case(:linear)
    unit = OpenSHOP._river_replace(only(c.system.generators); pmax = 0.7)
    c = OpenSHOP._river_replace(
        c;
        system = OpenSHOP._river_replace(c.system; generators = [unit]),
    )
    verified = solve_verified(
        c;
        u = ones(Int, 1, 1),
        max_refinements = 0,
        operational_margin = 0.1,
        time_limit = 10.0,
    )
    @test verified["accepted"]
    @test verified["solution"]["power"][1, 1] <= 0.60001
    @test verified["solution"]["power"][1, 1] > 0.5
    @test only(c.system.generators).pmax == 0.7
    @test validate(c, verified["solution"])["valid"]
    @test_throws ArgumentError solve_verified(c; operational_margin = -0.1)
    @test_throws ArgumentError schedule_case(c; operational_margin = Inf)
end
