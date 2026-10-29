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
    routing = OpenSHOP._dispatch_transport_expressions(model,c,release)
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
    @test !replay_audit(c,x)["valid"]
    # Every constant release fails: the midpoint needs q ≤ 170/9, while
    # the final upper storage bound needs q ≥ 212/9 m³/s.
    proposal = OpenSHOP.propose_commitment(c;time_limit=60.0)
    @test proposal["status"] == "INFEASIBLE"
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
