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
        routing=OpenSHOP._dispatch_transport_expressions(c, q)
        evaluate(x)=value(start_value, x)
        arrivals=evaluate.(routing.arrivals)
        terminal=evaluate.(routing.terminal)
        @test routing.exact==exact
        if exact
            # Exact confluence routing carries upstream cohorts through the junction.
            expected=OpenSHOP.route_network_exact(c, rates)
            @test arrivals≈expected["arrival_volume"] atol=1e-13
            @test terminal≈expected["transit"][:, end] atol=1e-13
            cached=OpenSHOP._dispatch_transport_expressions(c, q;
                transport = routing.nd)
            @test evaluate.(cached.arrivals)≈arrivals atol=1e-13
            @test evaluate.(cached.terminal)≈terminal atol=1e-13
        else
            # Distributed routing transports each prescribed reach release independently.
            for (i, r) in enumerate(rivers)
                expected=OpenSHOP.route_volumes(r, grid, grid, rates[i, :])+
                    OpenSHOP.route_volumes(r, grid, r.history_grid, r.history_release)
                inventory=OpenSHOP.remaining_volume(r, last(grid), grid, rates[i, :])+
                    OpenSHOP.remaining_volume(r, last(grid), r.history_grid, r.history_release)
                @test arrivals[i, :]≈expected atol=1e-13
                @test terminal[i]≈inventory atol=1e-13
                initial=OpenSHOP.remaining_volume(r, first(grid), r.history_grid,
                    r.history_release)
                @test sum(arrivals[i, :])+terminal[i]≈
                    initial+0.0036sum(diff(grid).*rates[i, :]) atol=1e-13
            end
        end
    end
    empty_case=warm_dispatch_fixture(; rivers = false)
    empty_routing=OpenSHOP._dispatch_transport_expressions(empty_case, zeros(0, 2))
    @test size(empty_routing.arrivals)==(0, 2)
    @test isempty(empty_routing.terminal)
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
    b=OpenSHOP._build_dispatch(c; u = ones(Int, 2, 2))
    @test upper_bound.(b.m[:rq])≈[0.05 0.03]
    @test lower_bound.(b.a)≈[0.1 0.2]
    @test upper_bound.(b.a)≈[0.8 0.6]
end
