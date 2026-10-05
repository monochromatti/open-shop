function outlet_test_case(receiver_head; floor = 100.0, shutdown = 0.0, periods = 1)
    source=Reservoir(
        name = :Source,
        z0 = 200.0,
        slope = 1.0,
        v0 = 2.0,
        vmin = 1.0,
        vmax = 3.0,
        inflow = 2.0,
        water_value = 10.0,
    )
    receiver=Reservoir(
        name = :Receiver,
        z0 = receiver_head-2.0,
        slope = 1.0,
        v0 = 2.0,
        vmin = 1.0,
        vmax = 3.0,
        water_value = 5.0,
    )
    plant=Plant(
        name = :Plant,
        source = :Source,
        target = :Receiver,
        pmax = 50.0,
        outlet_head_floor = floor,
        tailwater_curve = TableCurve([0.0, 5.0], [0.0, 0.5]),
    )
    unit=Generator(
        name = :Unit,
        plant = :Plant,
        qmin = 1.0,
        qmax = 5.0,
        pmin = 0.1,
        pmax = 50.0,
        efficiency = 0.9,
        qbest = 3.0,
        hbest = 100.0,
        hmin = 20.0,
        hmax = 150.0,
        qcurvature = 0.0,
        hcurvature = 0.0,
        initial_on = 1,
        initial_age = 8.0,
        minup = 0.0,
        mindown = 0.0,
        startup = 0.0,
        shutdown = shutdown,
    )
    reach=River(
        name = :Downstream,
        source = :Receiver,
        target = :Sea,
        law = :controlled,
        capacity = 4.0,
        curves = RiverRouting.DelayCurve[],
        deterministic_delay = 0.0,
        water_value = 5.0,
    )
    system=HydroSystem(;
        reservoirs = [source, receiver],
        junctions = Junction[],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
        tunnels = Tunnel[],
        plants = [plant],
        generators = [unit],
        rivers = [reach],
    )
    ScheduleCase(;
        name = "outlet_head_test",
        system,
        grid = collect(0.0:periods),
        prices = fill(50.0, periods),
    )
end

@testset "Tabulated river law needs no analytic coefficient" begin
    original = outlet_test_case(80.0)
    river = OpenSHOP._river_replace(
        only(original.system.rivers);
        law = :weir,
        coefficient = 0.0,
        discharge_curve = TableCurve([79.0, 81.0], [1.0, 3.0]),
    )
    system = OpenSHOP._river_replace(original.system; rivers = [river])
    case = OpenSHOP._river_replace(original; system)
    @test OpenSHOP.validate_inputs(case)
    solution = dispatch_from_controls(case, ones(Int, 1, 1), fill(2.0, 1, 1), ones(1, 1))
    @test solution["validation"]["valid"]
    @test solution["river_release"][1, 1]≈2.0 atol=1e-8
end

@testset "Outlet floor clamps turbine head below the floor" begin
    c=outlet_test_case(80.0)
    @test OpenSHOP.validate_inputs(c)
    x=dispatch_from_controls(c, ones(Int, 1, 1), fill(2.0, 1, 1), fill(0.5, 1, 1))
    @test x["validation"]["valid"]
    @test x["H"][2, 1]≈80.0 atol=1e-9
    @test x["power"][1, 1]≈0.00981*2*(202-100-0.2)*0.9 atol=1e-9
    @test OpenSHOP.outlet_head(only(c.system.plants), 80.0)==100.0
    local_model=OpenSHOP._build_dispatch(c)
    global_model=OpenSHOP._build_global_dispatch(c; joint = true)
    @test local_model.m isa JuMP.Model
    @test !any(
        v->startswith(JuMP.name(v), "outlet_head_"),
        JuMP.all_variables(global_model.m),
    )
    ordinary=OpenSHOP._build_global_dispatch(
        outlet_test_case(80.0; floor = nothing);
        joint = true,
    )
    @test JuMP.num_variables(global_model.m)==JuMP.num_variables(ordinary.m)
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    invalid=outlet_test_case(80.0; floor = NaN)
    @test_throws ErrorException OpenSHOP.validate_inputs(invalid)
end

@testset "Outlet floor leaves higher receiver heads unchanged" begin
    c=outlet_test_case(120.0)
    x=dispatch_from_controls(c, ones(Int, 1, 1), fill(2.0, 1, 1), fill(0.5, 1, 1))
    @test x["validation"]["valid"]
    @test x["power"][1, 1]≈0.00981*2*(202-120-0.2)*0.9 atol=1e-9
    @test OpenSHOP.outlet_head(only(c.system.plants), 120.0)==120.0
    global_model=OpenSHOP._build_global_dispatch(c; joint = true)
    @test !any(
        v->startswith(JuMP.name(v), "outlet_head_"),
        JuMP.all_variables(global_model.m),
    )
end

@testset "Outlet reference preserves receiver delivery and downstream law" begin
    floored=outlet_test_case(80.0)
    ordinary=outlet_test_case(80.0; floor = nothing)
    a=dispatch_from_controls(floored, ones(Int, 1, 1), fill(2.0, 1, 1), fill(0.5, 1, 1))
    b=dispatch_from_controls(ordinary, ones(Int, 1, 1), fill(2.0, 1, 1), fill(0.5, 1, 1))
    @test a["V"]≈b["V"] atol=1e-10
    @test a["river_release"]≈b["river_release"] atol=1e-10
    @test a["V"][2, end]≈2.0 atol=1e-10
    simulation=simulate(floored, a["generator_q"], a["gate"])
    @test only(simulation["boundary_outflow"])≈2.0 atol=1e-10
    @test a["river_release"][1, 1]≈4*0.5 atol=1e-10
    @test a["power"][1, 1]<b["power"][1, 1]
    @test replay_audit(floored, a)["valid"]
    @test only(floored.system.plants).target==:Receiver
    crossing=outlet_test_case(100.0)
    global_model=OpenSHOP._build_global_dispatch(crossing; joint = true)
    @test JuMP.variable_by_name(global_model.m, "outlet_head_1_1")!==nothing
    @test JuMP.variable_by_name(global_model.m, "outlet_head_1_1_cell[1]")!==nothing
    # Registered local operator uses the exact piecewise linear graph too.
    @test OpenSHOP._build_dispatch(crossing).m isa JuMP.Model
end

@testset "Delayed rivers discharge to boundaries with conserved inventory" begin
    c=outlet_test_case(80.0)
    reach=OpenSHOP._river_replace(only(c.system.rivers); deterministic_delay = 0.25)
    c=OpenSHOP._river_replace(
        c;
        system = OpenSHOP._river_replace(c.system; rivers = [reach]),
    )
    @test OpenSHOP.validate_inputs(c)
    @test river_order(c)==[1]
    x=dispatch_from_controls(c, ones(Int, 1, 1), fill(2.0, 1, 1), fill(0.5, 1, 1))
    @test x["validation"]["valid"]
    z=simulate(c, x["generator_q"], x["gate"])
    @test only(z["boundary_outflow"])≈1.5 atol=1e-10
    @test only(x["terminal_transit"])≈0.0036*2*0.25 atol=1e-10
    inflow=sum(r.inflow for r in c.system.reservoirs)
    @test sum(z["V"][:, end])+sum(z["transit"][:, end])-sum(z["V"][:, 1])-sum(
        z["transit"][:, 1],
    ) ≈ 0.0036*(inflow-only(z["boundary_outflow"])) atol=1e-10
    @test replay_audit(c, x)["valid"]
    b=OpenSHOP._build_global_dispatch(c; joint = true)
    @test OpenSHOP._lift_start!(b, c, x)["valid"]
    backwards=OpenSHOP._river_replace(
        c;
        system = OpenSHOP._river_replace(
            c.system;
            rivers = [OpenSHOP._river_replace(reach; source = :Sea, target = :Receiver)],
        ),
    )
    @test_throws ArgumentError river_order(backwards)
end

@testset "Low measured turbine efficiency is a declared operating choice" begin
    c=outlet_test_case(80.0)
    unit=OpenSHOP._river_replace(
        only(c.system.generators);
        efficiency = 0.65,
        min_efficiency = 0.6,
    )
    c=OpenSHOP._river_replace(
        c;
        system = OpenSHOP._river_replace(c.system; generators = [unit]),
    )
    @test OpenSHOP.validate_inputs(c)
    @test Generator(
        name = :Default,
        plant = :Plant,
        qmin = 1.0,
        qmax = 2.0,
        pmin = 0.1,
        pmax = 5.0,
        hmin = 1.0,
        hmax = 200.0,
    ).min_efficiency==0.0
    x=dispatch_from_controls(c, ones(Int, 1, 1), fill(2.0, 1, 1), fill(0.5, 1, 1))
    @test x["validation"]["valid"]
    @test x["power"][1, 1]≈0.00981*2*(202-100-0.2)*0.65 atol=1e-9
    @test replay_audit(c, x)["valid"]
    @test OpenSHOP._build_dispatch(c).m isa JuMP.Model
    b=OpenSHOP._build_global_dispatch(c; joint = true)
    @test OpenSHOP._lift_start!(b, c, x)["valid"]
    forbidden=OpenSHOP._river_replace(
        c;
        system = OpenSHOP._river_replace(
            c.system;
            generators = [OpenSHOP._river_replace(unit; min_efficiency = 0.7)],
        ),
    )
    @test !validate(forbidden, x)["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
end

@testset "Shutdown cost is charged once, including initial shutdown" begin
    c=outlet_test_case(80.0; shutdown = 37.0, periods = 2)
    c=OpenSHOP._river_replace(
        c;
        operations = [
            OperationalSeries(
                object = :Unit,
                attribute = :forced_on,
                times = [0.0],
                values = [0.0],
            ),
        ],
    )
    x=dispatch_from_controls(c, zeros(Int, 1, 2), zeros(1, 2), zeros(1, 2))
    @test x["validation"]["valid"]
    @test replay_audit(c, x)["valid"]
    income=sum(c.prices[t]*diff(c.grid)[t]*x["power"][1, t] for t in 1:2)
    water=sum(
        r.water_value*(x["V"][i, end]-r.v0) for (i, r) in enumerate(c.system.reservoirs)
    )
    @test x["objective"]≈income+water-37.0 atol=1e-10
    bad=copy(x)
    bad["objective"]+=37.0
    @test !validate(c, bad)["valid"]
    zero=OpenSHOP._river_replace(
        c;
        system = OpenSHOP._river_replace(
            c.system;
            generators = [
                OpenSHOP._river_replace(only(c.system.generators); shutdown = 0.0),
            ],
        ),
    )
    assign(v) = JuMP.name(v)=="sd[1,1]" ? 1.0 : JuMP.is_fixed(v) ? JuMP.fix_value(v) : 0.0
    for builder in (
        z->OpenSHOP._build_dispatch(z; u = zeros(Int, 1, 2)),
        z->OpenSHOP._build_global_dispatch(z; joint = true),
    )
        charged=builder(c)
        ordinary=builder(zero)
        @test JuMP.value(assign, charged.obj)-JuMP.value(assign, ordinary.obj)≈-37.0 atol=1e-10
    end
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    @test_throws ErrorException OpenSHOP.validate_inputs(
        outlet_test_case(80.0; shutdown = -1.0),
    )
end
