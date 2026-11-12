function storage_path_case(grid; fixed_outlet=false)
    reservoirs = [
        Reservoir(name=:Upper, v0=10.0, vmin=9.0, vmax=11.0, z0=100.0,
            slope=0.1, water_value=0.0),
        Reservoir(name=:Pond, v0=0.002, vmin=0.0004, vmax=0.0036, z0=50.0,
            slope=0.1, water_value=0.0),
    ]
    rivers = [
        River(name=:Inlet, source=:Upper, target=:Pond, law=:controlled,
            capacity=60.0, gate_min=1.0, water_value=0.0,
            curves=[RiverRouting.DelayCurve(0.0, [0.0, 2.0], [1.0])],
            history_grid=[-24.0, 0.0], history_release=[12.0]),
        River(name=:Outlet, source=:Pond, target=:Sea, law=:controlled,
            capacity=64.0, water_value=0.0, deterministic_delay=0.0,
            curves=RiverRouting.DelayCurve[]),
    ]
    operations = fixed_outlet ? [OperationalSeries(object=:Outlet, attribute=bound,
        times=grid[1:end-1], values=[(12 + 12 * (a + b)) / 64
            for (a,b) in zip(grid[1:end-1],grid[2:end])])
        for bound in (:gate_min, :gate_max)] : OperationalSeries[]
    ScheduleCase(;name="storage between control boundaries",
        system=HydroSystem(;reservoirs, rivers, boundaries=[Boundary(name=:Sea, head=0.0)]),
        grid, prices=zeros(length(grid)-1), operations,
        routing_grid=OpenSHOP.refined_grid(grid))
end

@testset "River volume windows agree with independent response integration" begin
    base = storage_path_case([0.0,0.4,1.0])
    inlet = OpenSHOP._river_replace(first(base.system.rivers);curves=[
        RiverRouting.DelayCurve(q,edges,[1.0]) for (q,edges) in
        ((0.0,[0.1,0.7]),(30.0,[0.05,0.4]),(60.0,[0.0,0.2]))])
    c = OpenSHOP._river_replace(base;system=OpenSHOP._river_replace(base.system;
        rivers=[inlet,last(base.system.rivers)]))
    model = Model()
    release = @variable(model,[1:2,1:2],lower_bound=0.0)
    routing = OpenSHOP.RiverTransport(model,c,release)
    for flow in (0.0,30.0,37.0,60.0)
        set_start_value.(release,[flow flow;0.0 0.0])
        for (a,b) in ((0.0,1.0),(0.2,1.0),(0.13,0.83))
            expression = OpenSHOP.routing_volume_expression(model,c,routing,1,a,b)
            expected = only(OpenSHOP.route_volumes(inlet,[a,b],c.grid,[flow,flow])) +
                only(OpenSHOP.route_volumes(inlet,[a,b],inlet.history_grid,inlet.history_release))
            @test value(start_value,expression) ≈ expected atol=1e-12
        end
    end
end

@testset "Storage paths reject an impossible hourly release" begin
    c = storage_path_case([0.0, 1.0])
    # The arriving rate is 12 + 24t m³/s. Releasing 24 returns the pond to
    # 2,000 m³ at the hour boundary but drains it below zero halfway through.
    x = dispatch_from_controls(c, zeros(Int,0,1), zeros(0,1), reshape([1.0,24/64],2,1))
    @test x.storage[2,end] ≈ 0.002 atol=1e-10
    @test !validate(c,x)["valid"]
    @test !replay_audit(c,x)["valid"]
    # Every constant release fails: the midpoint needs q ≤ 170/9, while
    # the final upper storage bound needs q ≥ 212/9 m³/s.
    proposal = OpenSHOP.propose_commitment(c;time_limit=60.0)
    @test proposal["status"] == "INFEASIBLE"
end

@testset "Small ponds retain their physical bounds between arrival knots" begin
    base=storage_path_case([0.0,1.0])
    pond=OpenSHOP._river_replace(base.system.reservoirs[2];
        v0=0.0005,vmin=0.000499,vmax=0.001499)
    inlet=OpenSHOP._river_replace(first(base.system.rivers);capacity=12.0+1/60)
    system=OpenSHOP._river_replace(base.system;
        reservoirs=[first(base.system.reservoirs),pond],rivers=[inlet,last(base.system.rivers)])
    c=OpenSHOP._river_replace(base;system,routing_grid=[0.0,1.0])
    # Arrival is 12+t/120 m³/s. Its hourly mean release restores both
    # endpoints to 500 m³, but the quadratic path reaches 496.25 m³.
    controls=reshape([1.0,(12.0+1/240)/64],2,1)
    x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),controls)
    @test x.storage[2,end]≈0.0005 atol=1e-12
    audit=validate(c,x)
    @test !audit["valid"]
    @test audit["residuals"]["storage_path_Pond"]≈2.75e-6 atol=1e-12
    @test any(t->isapprox(t,0.5;atol=1e-10),audit["storage_path_times"])
    @test !replay_audit(c,x)["valid"]
    # The same flow is safe with four extra cubic metres below the initial
    # level. Detect the excursion analytically, without adding time steps.
    safe=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(system;
        reservoirs=[first(system.reservoirs),OpenSHOP._river_replace(pond;vmin=0.000495)]))
    y=dispatch_from_controls(safe,zeros(Int,0,1),zeros(0,1),controls)
    @test validate(safe,y)["valid"]
    @test isempty(validate(safe,y)["storage_path_times"])
    changed=deepcopy(y)
    changed.storage[2,:].+=2e-6
    @test !validate(safe,changed)["valid"]
end

@testset "Storage paths retain feasible changing releases" begin
    c = storage_path_case(collect(0.0:0.25:1.0);fixed_outlet=true)
    result = solve_verified(c;u=zeros(Int,0,4),time_limit=60.0,max_refinements=0)
    @test result["accepted"]
    if result["accepted"]
        x = result["solution"]
        @test vec(x.river_release[2,:]) ≈ [15.0,21.0,27.0,33.0] atol=1e-5
        @test x.storage[2,:] ≈ fill(0.002,5) atol=1e-9
        @test replay_audit(c,x)["valid"]
    end
end

@testset "Dispatch repairs an interior storage minimum before optimizing revenue" begin
    for origin in (0.0,1000.0)
        reservoirs=[Reservoir(name=:Upper,z0=100.0,slope=1.0,v0=1.0,
            vmin=0.0,vmax=2.0,water_value=0.0),
            Reservoir(name=:Pond,z0=50.0,slope=1.0,v0=0.002,
                vmin=0.0004,vmax=0.04,water_value=0.0)]
        inlet=River(name=:Inlet,source=:Upper,target=:Pond,law=:controlled,
            capacity=36.0,gate_min=1.0,water_value=0.0,
            curves=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])],
            history_grid=[origin-2.0,origin],history_release=[12.0])
        unit=Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=40.0,
            pmin=0.0,pmax=40.0,hmin=40.0,hmax=60.0,qbest=20.0,hbest=50.0,
            efficiency=0.9,qcurvature=0.0,hcurvature=0.0,
            minup=0.0,mindown=0.0,startup=0.0)
        system=HydroSystem(;reservoirs,rivers=[inlet],generators=[unit],
            plants=[Plant(name=:Station,source=:Pond,target=:Sea,pmax=40.0)],
            boundaries=[Boundary(name=:Sea,head=0.0)])
        c=ScheduleCase(;name="interior storage minimum",system,
            grid=[origin,origin+1.0],prices=[100.0],routing_grid=[origin,origin+1.0])
        # Arrival is 12+24t m³/s. The minimum of
        # V(t)=2000+3600*((12-q)*t+12t²) m³ limits constant discharge to
        # q=12+sqrt(2*24*1600/3600), although endpoint bounds permit more.
        optimal_flow=12+sqrt(2*24*1600/3600)
        result=solve_verified(c;u=ones(Int,1,1),time_limit=20.0,max_refinements=12)
        @test result["accepted"]
        if result["accepted"]
            x=result["solution"]
            @test validate(result["problem"],x)["valid"]
            @test replay_audit(result["problem"],x;grid=OpenSHOP.refined_grid(x.grid;factor=4))["valid"]
            @test maximum(abs,x.generator_q.-optimal_flow)<=0.02
        end
    end
end
