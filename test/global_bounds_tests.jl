function bounds_only_case(; reservoirs = Reservoir[], junctions = Junction[],
    boundaries = Boundary[], tunnels = Tunnel[], plants = Plant[],
    generators = Generator[], rivers = River[], river_junctions = RiverJunction[],
    operations = OperationalSeries[], grid = [0.0, 1.0])
    system=HydroSystem(; reservoirs, junctions, boundaries, tunnels, plants,
        generators, rivers, river_junctions)
    ScheduleCase(; name = "capacity_domains", system, grid,
        prices = fill(50.0, length(grid)-1), operations)
end

function bounds_only_domains(c)
    OpenSHOP.validate_inputs(c)
    exact=all(r->r.deterministic_delay!==nothing, c.system.rivers)
    nd=exact ? OpenSHOP._transport_data(c, nothing) : nothing
    rd=exact ? nothing : routing_data(c)
    OpenSHOP._global_capacity_bounds(c, nd, rd)
end

function assert_domain_enclosure(c, d, x; tolerance = 1e-7)
    @test x["validation"]["valid"]
    @test all(d.lower .<= x["V"] .+ tolerance)
    @test all(d.upper .>= x["V"] .- tolerance)
    R=length(c.system.reservoirs)
    @test all(d.hlo .<= x["H"][1:R,:] .+ tolerance)
    @test all(d.hhi .>= x["H"][1:R,:] .- tolerance)
    @test all(d.arrival_lower .<= x["arrival_volume"] .+ tolerance)
    @test all(d.arrival_upper .>= x["arrival_volume"] .- tolerance)

end

@testset "Capacity domains enclose junction dispatch and dated operating limits" begin
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
    @test d.lower[1,1]==d.upper[1,1]==lake.v0
    @test d.lower[1,2]≈lake.v0-0.0036*(20.0-3.0)-1e-7
    @test d.upper[1,2]≈lake.v0+0.0036*(20.0+3.0)+1e-7
    x=dispatch_from_controls(c, [1 0], [3.0 0.0], zeros(0, 2))
    assert_domain_enclosure(c, d, x)
end

@testset "Capacity bounds preserve either flow sign and closed-head independence" begin
    for difference in (-4.0, 4.0)
        c=bounds_only_case(;
            boundaries = [Boundary(name = :A, head = 100.0+difference),
                Boundary(name = :B, head = 100.0)],
            tunnels = [Tunnel(name = :Link, source = :A, target = :B,
                resistance = 1.0, capacity = 5.0)])
        d=bounds_only_domains(c)
        x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1))
        @test x["objective"]==0.0
        @test x["tunnel_q"][1,1]≈copysign(2.0,difference)
        assert_domain_enclosure(c,d,x)
    end
    c=bounds_only_case(;
        boundaries = [Boundary(name = :A, head = 100.0), Boundary(name = :B, head = 104.0)],
        tunnels = [Tunnel(name = :Closed, source = :A, target = :B,
            resistance = 1.0, capacity = 5.0, opening = 0.0)])
    d=bounds_only_domains(c)
    x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1))
    @test x["objective"]==0.0
    @test x["tunnel_q"]==zeros(1,1)
    @test x["H"][:,1]==[100.0,104.0]
    assert_domain_enclosure(c,d,x)
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
    x=dispatch_from_controls(c, zeros(Int, 0, 2), zeros(0, 2), [0.4 0.0; 0.0 1.0])
    @test x["tunnel_q"][1, 1]<0<x["tunnel_q"][1, 2]
    assert_domain_enclosure(c, d, x)
end

@testset "Capacity enclosures cover river laws, gates, dry branches and soft minima" begin
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
        x=dispatch_from_controls(c, zeros(Int, 0, 1), zeros(0, 1),
            fill(law==:weir ? 1.0 : 0.25, 1, 1))
        @test 0<x["river_release"][1,1]<=expected_cap
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
    assert_domain_enclosure(c, d,
        dispatch_from_controls(c, zeros(Int, 0, 1), zeros(0, 1), ones(1, 1)))
    controlled=OpenSHOP._river_replace(r; law = :controlled, allow_dry = false)
    soft=bounds_only_case(; reservoirs = [dry], rivers = [controlled],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
        operations = [OperationalSeries(object = :DryOutlet, attribute = :min_release,
                times = [0.0], values = [2.0]),
            OperationalSeries(object = :DryOutlet, attribute = :release_penalty,
                times = [0.0], values = [10.0])])
    assert_domain_enclosure(soft, bounds_only_domains(soft),
        dispatch_from_controls(soft, zeros(Int, 0, 1), zeros(0, 1), zeros(1, 1)))
end

@testset "Delayed confluence capacity bounds preserve routed cohorts" begin
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
    @test d.arrival_upper[3,1]≈0.0036*(0.25*3.0+0.5*4.0)
    @test d.arrival_upper[3,2]≈0.0036*(3.0+4.0)
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
    for gate in (0.0, 0.2, 0.4)
        x=dispatch_from_controls(c, zeros(Int, 0, 2), zeros(0, 2), fill(gate, 1, 2))
        assert_domain_enclosure(c, d, x)
    end
end

@testset "Forward reachability and future storage restrictions" begin
    for inflow in (-1000.0,1000.0)
        lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=1.0,vmin=0.0,vmax=2.0,
            inflow=inflow,water_value=0.0)
        c=bounds_only_case(;reservoirs=[lake])
        @test_throws ArgumentError bounds_only_domains(c)
    end
    lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=1.0,vmin=0.0,vmax=3.0,water_value=0.0)
    river=River(name=:Outlet,source=:Lake,target=:Sea,law=:controlled,capacity=10.0,
        curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,water_value=0.0)
    future(attribute,value)=bounds_only_case(;reservoirs=[lake],rivers=[river],
        boundaries=[Boundary(name=:Sea,head=0.0)],grid=[0.0,1.0,2.0,3.0],
        operations=[OperationalSeries(object=:Lake,attribute=attribute,times=[0.0,2.0],
            values=[attribute==:vmin ? 0.0 : 3.0,value])])
    minimum_target=bounds_only_domains(future(:vmin,0.99))
    maximum_target=bounds_only_domains(future(:vmax,0.95))
    # The backward pass carries future limits into the preceding storage vertex.
    @test minimum_target.lower[1,2]>=0.99-2e-7
    @test maximum_target.upper[1,2]<=0.95+0.0036*10.0+2e-7
    @test maximum_target.upper[1,2]<lake.v0
    @test_throws ArgumentError bounds_only_domains(future(:vmin,1.1))
    @test_throws ArgumentError bounds_only_domains(future(:vmax,0.8))
end

@testset "Polynomial and table extrema include interior points" begin
    @test OpenSHOP._global_quadratic_range(1.0,-1.0,0.0,1.0)==(0.0,0.25)
    @test OpenSHOP._global_quadratic_range(-1.0,1.0,0.0,1.0)==(-0.25,0.0)
    @test OpenSHOP._global_quadratic_range(2.0,0.0,-1.0,3.0)==(-2.0,6.0)
    @test OpenSHOP._global_quadratic_range(1.0,-1.0,1.0,2.0)==(-2.0,0.0)
    lake=Reservoir(name=:Lake,z0=100.0,slope=2.0,curvature=-1.0,v0=1.0,
        vmin=0.0,vmax=2.0,water_value=0.0)
    @test OpenSHOP._global_level_range(lake,0.0,2.0)==(100.0,101.0)
    curve=TableCurve([0.0,1.0,2.0],[100.0,104.0,101.0])
    tabulated=OpenSHOP._river_replace(lake;level_curve=curve)
    @test OpenSHOP._global_level_range(tabulated,0.5,1.5)==(102.0,104.0)
    @test OpenSHOP._global_polynomial_range((0.5,1.0,-1.0,0.0),0.0,1.0)==(0.5,0.75)
    @test OpenSHOP._global_polynomial_range((0.0,-1.0,0.0,1.0),-1.0,1.0)[2]≈2/(3sqrt(3))
end
