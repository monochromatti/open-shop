function river_reference_case(;history=[2.0,8.0])
    curves=[RiverRouting.DelayCurve(0.0,[1.0,2.0],[1.0]),
        RiverRouting.DelayCurve(10.0,[0.0,1.0],[1.0])]
    lake=Reservoir(name=:Lake,z0=100.0,v0=1.0,vmin=0.5,vmax=2.0,water_value=0.0)
    river=River(name=:Reach,source=:Lake,target=:Sea,law=:controlled,
        capacity=10.0,curves=curves,water_value=0.0,
        history_grid=[-2.0,-1.0,0.0],history_release=history)
    ScheduleCase(name="contemporaneous river reference",grid=[0.0,1.0,2.0,3.0],
        prices=zeros(3),system=HydroSystem(reservoirs=[lake],rivers=[river],
            boundaries=[Boundary(name=:Sea,head=0.0)]))
end

@testset "Each release and historical cohort selects its own delay distribution" begin
    c=river_reference_case()
    q=reshape([0.0,10.0,5.0],1,:)
    x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),q/10.0)
    # At q=2 the fast share is .2; at q=8 it is .8. Their original
    # historical pulses leave 5.6 m³/s-hours in transit at time zero.
    # Current q=10 is entirely fast; q=5 is half fast and half slow.
    arrivals=0.0036 .* [4.8,5.8,6.25]
    inventories=0.0036 .* [5.6,0.8,5.0,3.75]
    @test validate(c,x)["valid"]
    @test vec(x.river_release)≈vec(q) atol=1e-12
    @test vec(x.arrival_volume)≈arrivals atol=1e-12
    @test only(x.terminal_transit)≈last(inventories) atol=1e-12
    @test vec(x.storage)≈[1.0,1.0,0.964,0.946] atol=1e-10
    routed=route_network(c,q)
    @test vec(routed["transit"])≈inventories atol=1e-12
    @test vec(routed["arrival_volume"])≈arrivals atol=1e-12
    @test inventories[2:end]≈first(inventories).+cumsum(0.0036 .* vec(q).-arrivals) atol=1e-12
    @test [point_arrival(only(c.system.rivers),c.grid,vec(q),t) for t in (0.5,1.5,2.5)]≈
        [4.8,5.8,6.25] atol=1e-12

    # Swapping the first two releases changes arrivals, although total release
    # and the last cohort are unchanged. This is not a previous-interval lookup.
    swapped=route_network(c,reshape([10.0,0.0,5.0],1,:))
    @test vec(swapped["arrival_volume"])≈0.0036 .* [9.8,5.8,1.25] atol=1e-12
    @test vec(swapped["transit"])≈0.0036 .* [5.6,5.8,0.0,3.75] atol=1e-12
    other_history=route_network(river_reference_case(history=[8.0,2.0]),q)
    @test vec(other_history["arrival_volume"])≈0.0036 .* [1.8,5.8,6.25] atol=1e-12
    @test vec(other_history["transit"])≈0.0036 .* [2.6,0.8,5.0,3.75] atol=1e-12

    refined=with_grid(c,collect(0.0:0.5:3.0))
    fine=dispatch_from_controls(refined,zeros(Int,0,6),zeros(0,6),
        reshape([0.0,0.0,1.0,1.0,0.5,0.5],1,:))
    @test vec(sum(reshape(vec(fine.arrival_volume),2,:);dims=1))≈arrivals atol=1e-12
    @test replay_audit(c,x)["valid"]
    continued=restart_case(c,x,1.0)
    rest=dispatch_from_controls(continued,zeros(Int,0,2),zeros(0,2),reshape([1.0,0.5],1,:))
    @test validate(continued,rest)["valid"]
    @test rest.arrival_volume≈x.arrival_volume[:,2:end] atol=1e-12
    @test rest.storage≈x.storage[:,2:end] atol=1e-10
    @test rest.terminal_transit≈x.terminal_transit atol=1e-12
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
end

@testset "Delay reference endpoints belong to a declared routing domain" begin
    curves=only(river_reference_case().system.rivers).curves
    grid=[0.0,1.0,2.0,3.0]
    @test RiverRouting.route(grid,[0.0,1.0],[0.0],curves)==zeros(3)
    @test RiverRouting.route(grid,[0.0,1.0],[10.0],curves)≈[0.018,0.018,0.0] atol=1e-12
    @test RiverRouting.route(grid,[0.0,1.0],[5.0],curves)≈[0.0045,0.009,0.0045] atol=1e-12
    @test_throws ArgumentError RiverRouting.route(grid,[0.0,1.0],[10.1],curves)
    positive=[RiverRouting.DelayCurve(2.0,[1.0,2.0],[1.0]),last(curves)]
    @test_throws ArgumentError RiverRouting.route(grid,[0.0,1.0],[1.0],positive)
    # A single curve is explicitly independent of its descriptive reference.
    @test RiverRouting.route(grid,[0.0,1.0],[20.0],positive[1:1])≈[0.0,0.036,0.036] atol=1e-12
end
