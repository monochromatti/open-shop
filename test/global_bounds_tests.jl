function bounds_only_case(; reservoirs = Reservoir[], junctions = Junction[],
    boundaries = Boundary[], tunnels = Tunnel[], plants = Plant[],
    generators = Generator[], rivers = River[], river_junctions = RiverJunction[],
    operations = OperationalSeries[], grid = [0.0, 1.0])
    system=HydroSystem(; reservoirs, junctions, boundaries, tunnels, plants,
        generators, rivers, river_junctions)
    ScheduleCase(; name = "reachable_domains", system, grid,
        prices = fill(50.0, length(grid)-1), operations)
end

function bounds_only_domains(c; tightened = true)
    OpenSHOP.validate_inputs(c)
    exact=all(r->r.deterministic_delay!==nothing, c.system.rivers)
    nd=exact ? OpenSHOP._transport_data(c, nothing) : nothing
    rd=exact ? nothing : routing_data(c)
    OpenSHOP._global_reachable_bounds(c, nd, rd; tightened)
end

function assert_domain_enclosure(c, d, x; tolerance = 1e-7)
    @test x["validation"]["valid"]
    @test all(d.lower .<= x["V"] .+ tolerance)
    @test all(d.upper .>= x["V"] .- tolerance)
    @test all(d.tunnel_lower .<= x["tunnel_q"] .+ tolerance)
    @test all(d.tunnel_upper .>= x["tunnel_q"] .- tolerance)
    @test all(d.generator_lower .<= x["generator_q"] .+ tolerance)
    @test all(d.generator_upper .>= x["generator_q"] .- tolerance)
    @test all(d.release_lower .<= x["river_release"] .+ tolerance)
    @test all(d.release_upper .>= x["river_release"] .- tolerance)
    for (i, name) in enumerate(nodes(c.system)), t in eachindex(c.prices)
        lo, hi=d.node_head_bounds[(name, t)]
        @test lo-tolerance<=x["H"][i, t]<=hi+tolerance
    end
end

@testset "Junction continuity proves intake directions and operating flow caps" begin
    lake=Reservoir(name = :Lake, z0 = 100.0, slope = 1.0,
        v0 = 10.0, vmin = 9.0, vmax = 11.0, inflow = 3.0, water_value = 10.0)
    unit=Generator(name = :Unit, plant = :Station, qmin = 1.0, qmax = 6.0,
        pmin = 0.05, pmax = 8.0, hmin = 80.0, hmax = 120.0, hbest = 100.0,
        qbest = 3.0, qcurvature = 0.0, hcurvature = 0.0,
        minup = 0.0, mindown = 0.0, startup = 0.0)
    c=bounds_only_case(; reservoirs = [lake],
        junctions = [Junction(name = :Intake, hmin = 80.0, hmax = 120.0),
            Junction(name = :Penstock, hmin = 80.0, hmax = 120.0)],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
        tunnels = [Tunnel(name = :Main, source = :Lake, target = :Intake,
                resistance = 0.1, capacity = 20.0),
            Tunnel(name = :Branch, source = :Intake, target = :Penstock,
                resistance = 0.2, capacity = 20.0)],
        plants = [Plant(name = :Station, source = :Penstock, target = :Sea, pmax = 8.0)],
        generators = [unit], grid = [0.0, 1.0, 2.0],
        operations = [OperationalSeries(object = :Unit, attribute = :qmax,
                times = [0.0, 1.0], values = [4.0, 2.0]),
            OperationalSeries(object = :Unit, attribute = :forced_on,
                times = [0.0, 1.0], values = [1.0, 0.0])])
    d=bounds_only_domains(c)
    old=bounds_only_domains(c; tightened = false)
    @test !hasproperty(old, :tunnel_lower)
    @test d.generator_lower==[1.0 0.0]
    @test d.generator_upper==[4.0 0.0]
    @test all(d.tunnel_lower[:, 1] .>= 1.0-2e-7)
    @test all(d.tunnel_upper[:, 1] .<= 4.0+3e-7)
    @test d.tunnel_lower[:, 2]==[0.0, 0.0]
    @test d.tunnel_upper[:, 2]==[0.0, 0.0]
    @test d.node_head_bounds[(:Intake, 1)][1]>108.0
    @test d.upper[1, 2]<old.upper[1, 2]
    @test all(d.lower .>= old.lower)
    @test all(d.upper .<= old.upper)
    x=dispatch_from_controls(c, [1 0], [3.0 0.0], zeros(0, 2))
    assert_domain_enclosure(c, d, x)
end

@testset "Head-loss bounds preserve either sign and closed-head independence" begin
    for difference in (-4.0, 4.0)
        c=bounds_only_case(;
            boundaries = [Boundary(name = :A, head = 100.0+difference),
                Boundary(name = :B, head = 100.0)],
            tunnels = [Tunnel(name = :Link, source = :A, target = :B,
                resistance = 1.0, capacity = 5.0)])
        d=bounds_only_domains(c)
        flow=copysign(2.0, difference)
        @test d.tunnel_lower[1, 1]<=flow<=d.tunnel_upper[1, 1]
        @test d.tunnel_upper[1, 1]-d.tunnel_lower[1, 1]<=3e-7
        @test difference<0 ? d.tunnel_upper[1, 1]<0 : d.tunnel_lower[1, 1]>0
    end
    c=bounds_only_case(;
        boundaries = [Boundary(name = :A, head = 100.0), Boundary(name = :B, head = 104.0)],
        tunnels = [Tunnel(name = :Closed, source = :A, target = :B,
            resistance = 1.0, capacity = 5.0, opening = 0.0)])
    d=bounds_only_domains(c)
    @test d.tunnel_lower==d.tunnel_upper==zeros(1, 1)
    @test d.node_head_bounds[(:A, 1)]==(100.0, 100.0)
    @test d.node_head_bounds[(:B, 1)]==(104.0, 104.0)
end

@testset "A genuine reversible tunnel encloses a changing flow direction" begin
    reservoirs=[Reservoir(name = name, z0 = 100.0, slope = 1.0,
        v0 = 1.0, vmin = 0.5, vmax = 1.5, water_value = 10.0) for name in (:A, :B)]
    rivers=[River(name = Symbol("Drain", name), source = name, target = :Sea,
        law = :controlled, capacity = 50.0, curves = RiverRouting.DelayCurve[],
        deterministic_delay = 0.0, water_value = 0.0) for name in (:A, :B)]
    c=bounds_only_case(; reservoirs,
        boundaries = [Boundary(name = :Sea, head = 0.0)], rivers,
        tunnels = [Tunnel(name = :Reversible, source = :A, target = :B,
            resistance = 0.25, capacity = 3.0)], grid = [0.0, 1.0, 2.0])
    d=bounds_only_domains(c)
    @test all(d.tunnel_lower .< 0)
    @test all(d.tunnel_upper .> 0)
    x=dispatch_from_controls(c, zeros(Int, 0, 2), zeros(0, 2), [0.4 0.0; 0.0 1.0])
    @test x["tunnel_q"][1, 1]<0<x["tunnel_q"][1, 2]
    assert_domain_enclosure(c, d, x)
end

@testset "River-law reachability includes gates, tables, dry branches and soft minima" begin
    lake=Reservoir(name = :Lake, z0 = 10.0, slope = 0.1,
        v0 = 10.0, vmin = 9.0, vmax = 11.0, water_value = 10.0)
    for (law, curve, expected_cap) in ((:controlled, nothing, 30.0),
        (:orifice, nothing, 0.61), (:weir, nothing, 2.01),
        (:controlled, TableCurve([10.9, 11.1], [0.0, 8.0]), 1.3))
        r=River(name = :Outlet, source = :Lake, target = :Sea, law = law,
            coefficient = 2.0, crest = 10.0, capacity = 100.0,
            curves = RiverRouting.DelayCurve[], deterministic_delay = 0.0,
            water_value = 0.0, discharge_curve = curve)
        operations=[OperationalSeries(object = :Outlet, attribute = :gate_min,
            times = [0.0], values = [0.2])]
        law!=:weir && push!(operations, OperationalSeries(object = :Outlet,
            attribute = :gate_max, times = [0.0], values = [0.3]))
        c=bounds_only_case(; reservoirs = [lake], rivers = [r], operations,
            boundaries = [Boundary(name = :Sea, head = 0.0)])
        d=bounds_only_domains(c)
        @test d.release_upper[1, 1]<=expected_cap+2e-7
        @test d.release_lower[1, 1]>0
        x=dispatch_from_controls(c, zeros(Int, 0, 1), zeros(0, 1),
            fill(law==:weir ? 1.0 : 0.25, 1, 1))
        assert_domain_enclosure(c, d, x)
    end
    dry=Reservoir(name = :DryLake, z0 = -1.0, slope = 1.0,
        v0 = 0.1, vmin = 0.0, vmax = 0.2, water_value = 10.0)
    r=River(name = :DryOutlet, source = :DryLake, target = :Sea,
        law = :orifice, allow_dry = true, coefficient = 2.0,
        crest = 0.0, capacity = 100.0, curves = RiverRouting.DelayCurve[],
        deterministic_delay = 0.0, water_value = 0.0)
    c=bounds_only_case(; reservoirs = [dry], rivers = [r],
        boundaries = [Boundary(name = :Sea, head = 0.0)])
    d=bounds_only_domains(c)
    @test d.release_lower[1, 1]==0
    @test d.release_upper[1, 1]<=1e-7
    assert_domain_enclosure(c, d,
        dispatch_from_controls(c, zeros(Int, 0, 1), zeros(0, 1), ones(1, 1)))
    controlled=OpenSHOP._river_replace(r; law = :controlled, allow_dry = false)
    soft=bounds_only_case(; reservoirs = [dry], rivers = [controlled],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
        operations = [OperationalSeries(object = :DryOutlet, attribute = :min_release,
                times = [0.0], values = [2.0]),
            OperationalSeries(object = :DryOutlet, attribute = :release_penalty,
                times = [0.0], values = [10.0])])
    @test bounds_only_domains(soft).release_lower[1, 1]==0
    assert_domain_enclosure(soft, bounds_only_domains(soft),
        dispatch_from_controls(soft, zeros(Int, 0, 1), zeros(0, 1), zeros(1, 1)))
end

@testset "Delayed confluence propagates release caps and preserves cohorts" begin
    reservoirs=[Reservoir(name = name, z0 = 100.0, slope = 1.0,
        v0 = 2.0, vmin = 1.0, vmax = 3.0, water_value = 10.0) for name in (:A, :B)]
    rivers=[River(name = :First, source = :A, target = :Merge, law = :controlled,
            capacity = 3.0, curves = RiverRouting.DelayCurve[], deterministic_delay = 0.25,
            water_value = 0.0),
        River(name = :Second, source = :B, target = :Merge, law = :controlled,
            capacity = 4.0, curves = RiverRouting.DelayCurve[], deterministic_delay = 0.0,
            water_value = 0.0),
        River(name = :Downstream, source = :Merge, target = :Sea, law = :junction,
            capacity = 10.0, curves = RiverRouting.DelayCurve[], deterministic_delay = 0.5,
            water_value = 0.0)]
    c=bounds_only_case(; reservoirs, rivers,
        river_junctions = [RiverJunction(name = :Merge)],
        boundaries = [Boundary(name = :Sea, head = 0.0)], grid = [0.0, 1.0, 2.0],
        operations = [OperationalSeries(object = :First, attribute = :gate_max,
                times = [0.0], values = [0.2]),
            OperationalSeries(object = :Second, attribute = :gate_max,
                times = [0.0], values = [0.5])])
    d=bounds_only_domains(c)
    @test d.release_upper[1, 1]<=0.600001
    @test d.release_upper[2, 1]<=2.000001
    @test d.release_upper[3, 1]<=2.4501
    @test d.release_upper[3, 2]<=2.6001
    x=dispatch_from_controls(c, zeros(Int, 0, 2), zeros(0, 2),
        [0.1 0.2; 0.25 0.5; 0.0 0.0])
    assert_domain_enclosure(c, d, x)
end

@testset "Distributed delay quadratic bounds enclose independent routed schedules" begin
    lake=Reservoir(name = :Lake, z0 = 100.0, slope = 1.0,
        v0 = 2.0, vmin = 1.0, vmax = 3.0, water_value = 10.0)
    curves=[RiverRouting.DelayCurve(0.0, [0.0, 0.5, 1.0], [0.2, 0.8]),
        RiverRouting.DelayCurve(4.0, [0.0, 0.5, 1.0], [0.8, 0.2])]
    river=River(; name = :Reach, source = :Lake, target = :Sea,
        law = :controlled, capacity = 4.0, curves, water_value = 0.0)
    c=bounds_only_case(; reservoirs = [lake], rivers = [river],
        boundaries = [Boundary(name = :Sea, head = 0.0)], grid = [0.0, 1.0, 2.0],
        operations = [OperationalSeries(object = :Reach, attribute = :gate_max,
            times = [0.0], values = [0.4])])
    d=bounds_only_domains(c)
    @test maximum(d.release_upper)<=1.600001
    for gate in (0.0, 0.2, 0.4)
        x=dispatch_from_controls(c, zeros(Int, 0, 2), zeros(0, 2), fill(gate, 1, 2))
        assert_domain_enclosure(c, d, x)
        @test all(d.arrival_lower .<= x["arrival_volume"] .+ 1e-7)
        @test all(d.arrival_upper .>= x["arrival_volume"] .- 1e-7)
    end
end
