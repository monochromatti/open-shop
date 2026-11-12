@testset "Dispatch routing expressions preserve transported water" begin
    grid=[0.0, 0.2, 0.8, 1.5, 2.0]
    base=graph_case(; grid, history = true)
    rates=[3.0 7.0 2.0 6.0; 5.0 1.0 4.0 3.0; 2.0 8.0 1.0 5.0]
    for delay in (nothing, 0.35, 2.35)
        exact=delay!==nothing
        rivers=[
            if exact
                OpenSHOP._river_replace(r; curves = RiverRouting.DelayCurve[],
                    deterministic_delay = delay)
            else
                OpenSHOP._river_replace(r; curves = [
                    RiverRouting.DelayCurve(0.0, [0.1, 0.5, 3.0], [0.2, 0.8]),
                    RiverRouting.DelayCurve(30.0, [0.1, 0.5, 3.0], [0.8, 0.2]),
                ])
            end for r in base.system.rivers
        ]
        c=replace_graph(base; rivers)
        m=Model()
        @variable(m, q[1:3, 1:4])
        set_start_value.(q, rates)
        routing=OpenSHOP.RiverTransport(m,c, q)
        if !exact
            expected=OpenSHOP.routed_schedule(c,rates)
            set_start_value.(q,expected.release)
        end
        OpenSHOP.start_transport!(routing)
        evaluate(x)=value(start_value,x)
        arrivals=evaluate.(routing.arrivals)
        terminal=evaluate.(routing.terminal)
        expected=OpenSHOP.routed_schedule(c,rates)
        @test arrivals≈expected.arrival_volume atol=1e-13
        @test terminal≈expected.transit[:,end] atol=1e-13
        for (i,r) in enumerate(rivers)
            initial=OpenSHOP.remaining_volume(r,first(grid),r.history_grid,r.history_release)
            @test sum(arrivals[i,:])+terminal[i]≈
                initial+0.0036sum(diff(grid).*expected.release[i,:]) atol=1e-13
        end
        if exact
            cached=OpenSHOP.RiverTransport(m,c,q;transport=deterministic_network_data(c))
            @test evaluate.(cached.arrivals)≈arrivals atol=1e-13
            @test evaluate.(cached.terminal)≈terminal atol=1e-13
        end
    end
    empty_case=warm_dispatch_fixture(; rivers = false)
    empty_routing=OpenSHOP.RiverTransport(Model(),empty_case, zeros(0, 2))
    @test size(empty_routing.arrivals)==(0, 2)
    @test isempty(empty_routing.terminal)
end

@testset "Fixed travel references retain narrow operating intervals" begin
    duration=2.0^-12
    for origin in (0.0,nextfloat(2.0^40))
        grid=[origin,origin+duration,origin+1.0]
        river=River(name=:Reach,source=:Lake,target=:Sea,law=:controlled,
            capacity=10.0,water_value=0.0,
            curves=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0]),
                RiverRouting.DelayCurve(10.0,[1.0,2.0],[1.0])],
            history_grid=[origin-3.0,origin],history_release=[0.0])
        system=HydroSystem(reservoirs=[Reservoir(name=:Lake,v0=1.0,vmin=0.0,
            vmax=2.0,water_value=0.0)],rivers=[river],
            boundaries=[Boundary(name=:Sea,head=0.0)])
        c=ScheduleCase(;name="narrow operating interval",system,grid,prices=zeros(2))
        model=Model()
        @variable(model,0.0<=q[1:1,1:2]<=10.0)
        flows=[1.0 10.0]
        set_start_value.(q,flows)
        coordinate=(r,q,i,k,period)->OpenSHOP._proposal_delay_coordinate(r,flows[i,period])
        routing=OpenSHOP.RiverTransport(model,c,q;coordinate)
        arrival=sum(value.(start_value,routing.arrivals))
        # The first release uses 90% uniform [0,1] and 10% uniform [1,2].
        # The second release cannot arrive within the one-hour horizon.
        expected=0.9*0.0036*duration*(1-duration/2)
        @test arrival≈expected atol=1e-12
        @test arrival+only(value.(start_value,routing.terminal))≈
            0.0036sum(flows.*permutedims(diff(grid))) atol=1e-12
    end
end

@testset "Local dispatch retains time-varying outlet bounds" begin
    base=warm_dispatch_fixture()
    operations=vcat(base.operations, [
        OperationalSeries(object = :WarmReach, attribute = :capacity,
            times = [0.0, 1.0], values = [5.0, 3.0]),
        OperationalSeries(object = :WarmReach, attribute = :gate_min,
            times = [0.0, 1.0], values = [0.1, 0.2]),
        OperationalSeries(object = :WarmReach, attribute = :gate_max,
            times = [0.0, 1.0], values = [0.8, 0.6])])
    c=OpenSHOP._river_replace(base; operations)
    result=solve_case(c;u=ones(Int,2,2),time_limit=15.0)
    @test result.accepted
    if result.solution!==nothing
        @test all(result.solution.river_release .<= [5.0 3.0] .+ 1e-5)
        @test all([0.1 0.2] .- 1e-6 .<= result.solution.gate .<= [0.8 0.6] .+ 1e-6)
    end
end

@testset "A warm river network retains its transported water in dispatch" begin
    grid=[0.0,0.2,0.8,1.5,2.0]
    base=graph_case(;grid,history=true)
    curves=[RiverRouting.DelayCurve(0.0,[0.1,0.5,3.0],[0.2,0.8]),
        RiverRouting.DelayCurve(30.0,[0.1,0.5,3.0],[0.8,0.2])]
    for delay in (nothing,0.35)
        rivers=[OpenSHOP._river_replace(r;
            law=i<=2 ? :controlled : :junction,
            curves=delay===nothing ? curves : RiverRouting.DelayCurve[],
            deterministic_delay=delay) for (i,r) in enumerate(base.system.rivers)]
        c=replace_graph(base;rivers)
        u=zeros(Int,0,4)
        gates=[3.0 7.0 2.0 6.0;5.0 1.0 4.0 3.0;0.0 0.0 0.0 0.0]/30
        warm=dispatch_from_controls(c,u,zeros(0,4),gates)
        @test validate(c,warm)["valid"]
        result=solve_case(c;u,warm,replay_grid=OpenSHOP.refined_grid(grid),time_limit=15.0)
        @test result.accepted
        if result.accepted
            @test validate(c,result.solution)["valid"]
        end
        # A numerical starting point may need restoration. Initializing relay
        # equations must permit trial releases beyond the physical table domain.
        trial=deepcopy(schedule_dict(warm))
        trial["river_release"][1,1]=31.0
        @test solve_case(c;u,warm=trial,time_limit=15.0).accepted
    end
end
