@testset "Transit trajectories conserve released water on independent observation times" begin
    grid=[-3.0,-1.1,0.0,0.2,0.9,2.0]
    q=[8.0,0.0,20.0,5.0,12.0]
    times=[-4.0,-3.0,-1.7,-0.4,0.0,0.13,0.9,1.7,2.0,4.0,8.0]
    curves=[RiverRouting.DelayCurve(0.0,[0.1,0.3,2.0],[0.2,0.8]),
        RiverRouting.DelayCurve(20.0,[0.0,0.7,1.2],[0.7,0.3])]
    for delay in (nothing,0.0,0.37,3.0)
        r=River(name=:Reach,source=:Lake,target=:Sea,capacity=20.0,
            water_value=0.0,deterministic_delay=delay,
            curves=delay===nothing ? curves : RiverRouting.DelayCurve[])
        transit=OpenSHOP.remaining_volume(r,times,grid,q)
        released=[0.0036sum(q[k]*max(0.0,min(t,grid[k+1])-grid[k]) for k in eachindex(q))
            for t in times]
        arrived=cumsum(OpenSHOP.route_volumes(r,times,grid,q))
        @test transit[1]==0.0
        @test minimum(transit)>=0.0
        @test transit[2:end]+arrived≈released[2:end] atol=1e-12
        @test transit≈[OpenSHOP.remaining_volume(r,t,grid,q) for t in times] atol=1e-12
        @test OpenSHOP.remaining_volume(r,[0.13,0.13],grid,q)≈fill(transit[6],2) atol=1e-12
        @test OpenSHOP.remaining_volume(r,[0.13],grid,q)≈[transit[6]] atol=1e-12
        @test isempty(OpenSHOP.remaining_volume(r,Float64[],grid,q))
        @test_throws ArgumentError OpenSHOP.remaining_volume(r,reverse(times),grid,q)
        @test_throws ArgumentError OpenSHOP.remaining_volume(r,[NaN],grid,q)
    end
    # A one-hour release of 10^9 m³/s has just 1.8e-7 Mm³ left shortly
    # before the last arrival. Integrating the survival tail retains it.
    fixed=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])]
    t=2.0-1e-6
    expected=0.0036*1e9*(2.0-t)^2/2
    @test RiverRouting.remaining_volume(t,[0.0,1.0],[1e9],fixed)≈expected rtol=1e-12
    @test only(RiverRouting.remaining_volume([t],[0.0,1.0],[1e9],fixed))≈expected rtol=1e-12
end

@testset "Shifted deterministic releases preserve interval overlaps" begin
    grid=[-2.0,-0.2,0.0,0.1,0.8,2.0]
    q=[2.0,0.0,5.0,12.0,1.0]
    for delay in (0.0,0.37,3.0,1e16), arrivals in
        ([-5.0,-4.0],[-3.0,-0.5,0.0,0.3,1.1,4.0,8.0],[8.0,9.0])
        r=River(name=:Reach,source=:Lake,target=:Sea,capacity=20.0,
            water_value=0.0,deterministic_delay=delay,curves=RiverRouting.DelayCurve[])
        expected=[0.0036sum(q[k]*max(0.0,min(b,grid[k+1]+delay)-max(a,grid[k]+delay))
            for k in eachindex(q)) for (a,b) in zip(arrivals[1:end-1],arrivals[2:end])]
        @test OpenSHOP.route_volumes(r,arrivals,grid,q)≈expected atol=1e-12
    end
end

function transport_pulse_case(hours)
    reservoirs=[Reservoir(name=:Lake,z0=100.0,v0=100.0,vmin=0.0,vmax=200.0,water_value=0.0),
        Reservoir(name=:Receiver,z0=50.0,v0=100.0,vmin=0.0,vmax=200.0,water_value=0.0)]
    rivers=[River(name=:Source,source=:Lake,target=:Downstream,law=:controlled,
        capacity=200.0,water_value=0.0,
        curves=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])]),
        River(name=:Downstream,source=:auto,target=:Receiver,law=:junction,
            capacity=200.0,water_value=0.0,
            curves=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])])]
    ScheduleCase(name="two distributed reaches",system=HydroSystem(;reservoirs,rivers),
        grid=collect(0.0:hours),prices=zeros(hours))
end

@testset "Transport accuracy is retained when the horizon grows" begin
    for hours in (24,72,168)
        c=transport_pulse_case(hours)
        q=zeros(2,hours)
        q[1,1:3].=96.0
        controlled=route_network_controlled(c,q;absolute_tolerance=2e-4)
        @test controlled["converged"]
        @test controlled["estimated_cumulative_error_Mm3"]<=2e-4
        @test controlled["estimated_transit_error_Mm3"]<=2e-4
        # Two independent uniform [0,1] delays have the triangular density
        # f(d)=d on [0,1], 2-d on [1,2]. A three-hour rectangular release
        # therefore delivers 1/6, 5/6, 1, 5/6, 1/6 hours of flow.
        exact=0.0036*96.0.*[1/6,5/6,1.0,5/6,1/6]
        @test controlled["arrival_volume"][2,1:5]≈exact atol=2e-4
        @test sum(controlled["arrival_volume"][2,:])≈0.0036*96.0*3 atol=1e-12
    end
    c=transport_pulse_case(168)
    q=zeros(2,168);q[1,1:3].=96.0
    limited=route_network_controlled(c,q;absolute_tolerance=2e-4,max_internal_intervals=2048)
    @test !limited["converged"]
    @test_throws ArgumentError route_network_controlled(c,q;max_internal_intervals=0)
    @test_throws ArgumentError route_network_controlled(c,q;max_internal_intervals=NaN)
end
