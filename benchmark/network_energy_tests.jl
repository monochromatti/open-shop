using Test, JuMP, OpenSHOP
include("network_energy.jl")

function energy_test_system(; reservoirs = Reservoir[], junctions = Junction[],
        boundaries = Boundary[], tunnels = Tunnel[], plants = Plant[],
        generators = Generator[], rivers = River[])
    HydroSystem(; reservoirs, junctions, boundaries, tunnels, plants, generators, rivers)
end

function energy_test_lift(c, warm; fixed_u = nothing)
    @test validate(c, warm)["valid"]
    b = OpenSHOP._build_global_dispatch(c; joint = true, u = warm["u"], warm, fixed_u)
    baseline_variables = num_variables(b.m)
    baseline_constraints = num_constraints(b.m; count_variable_in_set_constraints = false)
    profile = add_network_energy!(b, c)
    @test num_variables(b.m) == baseline_variables + profile["auxiliary_variables"]
    @test num_constraints(b.m; count_variable_in_set_constraints = false) ==
        baseline_constraints + profile["energy_rows"] + profile["support_rows"] +
        profile["mccormick_rows"] + profile["distributed_injection_rows"]
    audit = OpenSHOP._lift_start!(b, c, warm)
    @test audit["valid"]
    @test audit["assigned"] == num_variables(b.m)
    assigned = Dict(v => start_value(v) for v in all_variables(b.m))
    at(x) = x isa Number ? x : value(v -> assigned[v], x)
    diagnostic = network_energy_diagnostic(b, c, at)
    @test diagnostic["max_original_positive_violation"] <= 1e-6
    @test all(is_fixed(v) || (has_lower_bound(v) && has_upper_bound(v)) || is_binary(v)
              for v in all_variables(b.m))
    b, assigned, diagnostic
end

@testset "Fenchel hydraulic energy, both signs and zero" begin
    for resistance in (0.1, 0.25, 2.0), q in (-20.0, -10.0, -1.0, 0.0, 1.0, 10.0, 20.0)
        d = resistance * q * abs(q)
        @test _ne_flow_content(q, resistance) + _ne_head_content(d, resistance) ≈ q * d
    end
    # Unlike a fit, the content/co-content inequality is valid for every pair.
    for resistance in (0.1, 0.25, 2.0), q in range(-20.0, 20.0; length = 13),
            d in range(-150.0, 150.0; length = 17)
        @test _ne_flow_content(q, resistance) + _ne_head_content(d, resistance) - q * d >= -1e-9
    end
end

@testset "Boundary exchange, reference shifts, reverse and closed tunnels" begin
    for (ha, hb, opening, expected) in ((100.0, 0.0, 1.0, 10.0),
            (0.0, 100.0, 1.0, -10.0), (-100.0, -100.0, 1.0, 0.0),
            (100.0, -300.0, 0.0, 0.0), (100.0, 0.0, 0.25, 5.0))
        s = energy_test_system(boundaries = [Boundary(name = :a, head = ha),
            Boundary(name = :b, head = hb)], tunnels = [Tunnel(name = :pipe,
            source = :a, target = :b, resistance = 1.0, capacity = 20.0, opening = opening)])
        c = ScheduleCase(name = "signed boundary energy", system = s, grid = [0.0, 1.0], prices = [1.0])
        warm = dispatch_from_controls(c, zeros(Int, 0, 1), zeros(0, 1), zeros(0, 1))
        @test warm["tunnel_q"][1, 1] ≈ expected atol = 1e-9
        b, assigned, d = energy_test_lift(c, warm)
        @test opening == 0 ? b.m.ext[:network_energy_profile]["energy_rows"] == 0 :
            d["max_original_positive_violation"] <= 1e-8
        # A closed gate must never couple otherwise unrelated heads.
        @test b.m.ext[:network_energy_profile]["closed_tunnel_intervals_excluded"] == Int(opening == 0)
        if expected != 0
            b2 = OpenSHOP._build_global_dispatch(c; joint = true, u = warm["u"], warm)
            add_network_energy!(b2, c; head_reference = -1234.0)
            @test OpenSHOP._lift_start!(b2, c, warm)["valid"]
            at2(x) = x isa Number ? x : value(v -> start_value(v), x)
            d2 = network_energy_diagnostic(b2, c, at2)
            @test d2["largest_errors"][1]["original_energy_residual"] ≈
                d["largest_errors"][1]["original_energy_residual"] atol = 1e-8
        end
    end
end

@testset "Coupled nonlinear storage, finite-delay river and unit off" begin
    a = Reservoir(name = :a, z0 = 100.0, slope = 2.0, curvature = 0.3,
        v0 = 3.0, vmin = 1.0, vmax = 4.0, inflow = 4.0, water_value = 1.0)
    z = Reservoir(name = :z, z0 = 10.0, slope = 0.5,
        v0 = 1.0, vmin = 0.1, vmax = 4.0, water_value = 1.0)
    unit = Generator(name = :unit, plant = :plant, qmin = 1.0, qmax = 6.0,
        pmin = 0.1, pmax = 10.0, efficiency = 0.9, qbest = 4.0, hbest = 90.0,
        hmin = 50.0, hmax = 120.0, qcurvature = 0.0, hcurvature = 0.0,
        initial_on = 0, minup = 0.0, mindown = 0.0, startup = 0.0)
    river = River(name = :river, source = :a, target = :z, curves = RiverRouting.DelayCurve[],
        law = :controlled, capacity = 5.0, deterministic_delay = 0.5,
        history_grid = [-2.0, -1.0, 0.0], history_release = [2.0, 2.0], water_value = 1.0)
    s = energy_test_system(reservoirs = [a, z], junctions = [Junction(name = :intake, hmin = 0.0, hmax = 130.0)],
        boundaries = [Boundary(name = :sea, head = 0.0)],
        tunnels = [Tunnel(name = :feed, source = :a, target = :intake, resistance = 0.2, capacity = 20.0),
            Tunnel(name = :drain, source = :z, target = :sea, resistance = 0.05, capacity = 30.0)],
        plants = [Plant(name = :plant, source = :intake, target = :z, pmax = 10.0)],
        generators = [unit], rivers = [river])
    c = ScheduleCase(name = "discrete storage and delay", system = s, grid = [0.0, 1.0, 2.0], prices = [1.0, 2.0])
    warm = dispatch_from_controls(c, reshape([1, 0], 1, 2), reshape([4.0, 0.0], 1, 2), fill(0.6, 1, 2))
    @test warm["river_release"][1, :] ≈ [3.0, 3.0]
    b, assigned, diagnostic = energy_test_lift(c, warm)
    at(x) = x isa Number ? x : value(v -> assigned[v], x)
    @test at(b.arrivals[1, 1]) / 0.0036 ≈ 2.5
    @test at(b.arrivals[1, 2]) / 0.0036 ≈ 3.0
    @test warm["power"][1, 2] == 0.0
    @test b.m.ext[:network_energy_profile]["mccormick_rows"] > 0
    @test b.m.ext[:network_energy_profile]["distributed_injection_rows"] == 0
    # The observed downstream source rate must use arrivals, not this release.
    first_component = b.m.ext[:network_energy_data].plan.records[1]
    z_injection = only(i for i in first_component.injections if i.i == 2)
    rate = at(z_injection.external)
    storage_rate = (warm["V"][2, 2] - warm["V"][2, 1]) / 0.0036
    @test rate ≈ 2.5 - storage_rate atol = 1e-8
    @test abs(rate - (3.0 - storage_rate)) > 0.4
    # Quadratic head proves that silently substituting a continuous primitive is wrong.
    v0, v1 = warm["V"][1, 1], warm["V"][1, 2]
    primitive(v) = a.z0 * v + a.slope * v^2 / 2 + a.curvature * v^3 / 3
    midpoint_work = OpenSHOP.head(a, (v0 + v1) / 2) * (v1 - v0)
    @test primitive(v1) - primitive(v0) - midpoint_work ≈ a.curvature * (v1 - v0)^3 / 12 atol = 1e-12
    @test abs(primitive(v1) - primitive(v0) - midpoint_work) > 1e-9

    distributed = River(name=:river, source=:a, target=:z,
        curves=[RiverRouting.DelayCurve(0.,[.5,1.,1.5],[.5,.5]),
                RiverRouting.DelayCurve(5.,[0.,.5,1.],[.5,.5])],
        law=:controlled, capacity=5., history_grid=[-2.,-1.,0.],
        history_release=[2.,2.], water_value=1.)
    sd=energy_test_system(reservoirs=s.reservoirs,junctions=s.junctions,
        boundaries=s.boundaries,tunnels=s.tunnels,plants=s.plants,
        generators=s.generators,rivers=[distributed])
    cd=ScheduleCase(name="distributed energy arrivals",system=sd,grid=c.grid,prices=c.prices)
    wd=dispatch_from_controls(cd,warm["u"],warm["generator_q"],warm["gate"])
    bd, _, _=energy_test_lift(cd,wd)
    @test bd.m.ext[:network_energy_profile]["distributed_injection_rows"]>0
end

@testset "Off negative head/efficiency and unsupported running domains" begin
    unit = Generator(name = :unit, plant = :plant, qmin = 1.0, qmax = 6.0,
        pmin = 0.1, pmax = 10.0, efficiency = 0.9, qbest = 3.0, hbest = 100.0,
        hmin = 10.0, hmax = 110.0, qcurvature = 0.0, hcurvature = 0.4,
        initial_on = 0, minup = 0.0, mindown = 0.0)
    s = energy_test_system(boundaries = [Boundary(name = :a, head = -100.0),
        Boundary(name = :b, head = 0.0)],
        tunnels = [Tunnel(name = :pipe, source = :a, target = :b, resistance = 1.0, capacity = 20.0)],
        plants = [Plant(name = :plant, source = :a, target = :b, pmax = 10.0)], generators = [unit])
    c = ScheduleCase(name = "negative off head", system = s, grid = [0.0, 1.0], prices = [1.0])
    warm = dispatch_from_controls(c, zeros(Int, 1, 1), zeros(1, 1), zeros(0, 1))
    @test OpenSHOP.efficiency(unit, 0.0, -100.0) < 0
    b, assigned, diagnostic = energy_test_lift(c, warm; fixed_u = zeros(Int, 1, 1))
    @test b.m.ext[:network_energy_profile]["energy_rows"] == 1
    @test isempty(b.m.ext[:network_energy_profile]["skipped_component_intervals"])

    # A negative tail contribution increases turbine head beyond network drop.
    # Omitting it would make the proposed upper bound invalid, so skip that component.
    badplant = Plant(name = :plant, source = :a, target = :b, pmax = 10.0,
        tailwater_curve = TableCurve([0.0, 6.0], [-2.0, -1.0]))
    bad = ScheduleCase(name = "negative extra loss", system = energy_test_system(
        boundaries = [Boundary(name = :a, head = 100.0), Boundary(name = :b, head = 0.0)],
        plants = [badplant], generators = [unit]), grid = [0.0, 1.0], prices = [1.0])
    @test_throws ArgumentError OpenSHOP._build_global_dispatch(bad; joint = true)

    zeroeff = Generator(name=:unit, plant=:plant, qmin=1., qmax=6.,
        pmin=.1, pmax=10., efficiency=0., qbest=3., hbest=100.,
        hmin=-100., hmax=110., qcurvature=0., hcurvature=0.,
        initial_on=0, minup=0., mindown=0.)
    zero = ScheduleCase(name="zero efficiency at negative running head",
        system=energy_test_system(
            boundaries=[Boundary(name=:a,head=-100.), Boundary(name=:b,head=0.)],
            plants=[Plant(name=:plant,source=:a,target=:b,pmax=10.)],
            generators=[zeroeff]), grid=[0.,1.], prices=[1.],
        operations=[OperationalSeries(object=:unit,attribute=:pmin,times=[0.],values=[0.])])
    # Public operating minima cannot loosen the strictly positive static pmin.
    @test_throws ArgumentError OpenSHOP._build_global_dispatch(zero;joint=true)
end

@testset "Tunnel capacities focus support domains" begin
    s=energy_test_system(junctions=[Junction(name=:a,hmin=-1000.,hmax=1000.),
        Junction(name=:b,hmin=-1000.,hmax=1000.)],
        tunnels=[Tunnel(name=:pipe,source=:a,target=:b,resistance=.01,capacity=5.)])
    c=ScheduleCase(name="bounded losses",system=s,grid=[0.,1.],prices=[1.])
    b=OpenSHOP._build_global_dispatch(c;joint=true)
    e=only(only(_ne_plan(b,c).records).tunnels)
    @test e.dlo ≈ -.25 atol=1e-8
    @test e.dhi ≈ .25 atol=1e-8
end

@testset "Separates a known optimistic quadratic tunnel relaxation" begin
    s = energy_test_system(boundaries = [Boundary(name = :a, head = 100.0), Boundary(name = :b, head = 0.0)],
        tunnels = [Tunnel(name = :pipe, source = :a, target = :b, resistance = 1.0, capacity = 20.0)])
    c = ScheduleCase(name = "optimistic root point", system = s, grid = [0.0, 1.0], prices = [1.0])
    warm = dispatch_from_controls(c, zeros(Int, 0, 1), zeros(0, 1), zeros(0, 1))
    b, assigned, _ = energy_test_lift(c, warm)
    # q=5, square-lift=100 satisfies the usual graph hull q²<=lift<=20q.
    # It creates the requested 100 m drop with only half the physical flow.
    @test 5.0^2 <= 100.0 <= 20.0 * 5.0
    assigned[variable_by_name(b.m, "q[1,1]")] = 5.0 / 50
    assigned[variable_by_name(b.m, "tunnel_positive_1_1")] = 5.0
    assigned[variable_by_name(b.m, "tunnel_negative_1_1")] = 0.0
    flow_lift = only(v for v in all_variables(b.m) if startswith(name(v), "energy_flow_"))
    assigned[flow_lift] = 0.0 # best relaxed epigraph value at q=5 for these supports
    at(x) = x isa Number ? x : value(v -> assigned[v], x)
    data = b.m.ext[:network_energy_data]
    @test at(data.rows[1].lhs) - at(data.rows[1].rhs) > 160.0
    @test network_energy_diagnostic(b, c, at)["max_original_positive_violation"] > 200.0
    for ref in all_constraints(b.m; include_variable_in_set_constraints = false)
        startswith(string(ref), "energy_flow_") || continue
        object = constraint_object(ref)
        @test at(object.func) >= object.set.lower - 1e-9
    end
    @test_throws ErrorException add_network_energy!(b, c)
    @test_throws ArgumentError add_network_energy!(b, c; support_points = 1)
end
