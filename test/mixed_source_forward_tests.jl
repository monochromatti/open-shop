using Test

function mixed_forward_case(; distributed=false, fixed_reach=nothing, receiver=:Sea)
    reservoirs=[Reservoir(name=:Upper,z0=100.0,slope=0.1,v0=5.0,
        vmin=1.0,vmax=10.0,water_value=0.0),
        Reservoir(name=:Lower,z0=60.0,slope=0.1,v0=5.0,
        vmin=1.0,vmax=10.0,water_value=0.0)]
    plants=[Plant(name=:UpperPlant,source=:Upper,target=:Lower,
        discharge_river=:Main,pmax=100.0),
        Plant(name=:LowerPlant,source=:Lower,target=:Sea,
        discharge_river=:Main,pmax=100.0)]
    generators=[Generator(name=Symbol("Unit",i),plant=p.name,qmin=0.1,
        qmax=10.0,pmin=0.01,pmax=100.0,hmin=0.0,hmax=200.0,
        qbest=3.0,hbest=100.0,efficiency=0.9,qcurvature=0.0,
        hcurvature=0.0,minup=0.0,mindown=0.0) for (i,p) in enumerate(plants)]
    tunnel=Tunnel(name=:Outlet,source=:Upper,target=:Sea,
        discharge_river=:Main,resistance=20.0,capacity=10.0)
    function reach(name,source,target,capacity,delay,inflow,law)
        distributed_here=distributed && name!=fixed_reach
        curves=distributed_here ? [RiverRouting.DelayCurve(q,[delay,delay+0.5],[1.0])
            for q in (0.0,capacity)] : RiverRouting.DelayCurve[]
        River(;name,source,target,capacity,inflow,law,curves,
            deterministic_delay=distributed_here ? nothing : delay,
            water_value=0.0,history_grid=[-1.0,0.0],history_release=[0.0])
    end
    rivers=[reach(:Tributary,:Upper,:Main,20.0,0.5,1.0,:controlled),
        reach(:Main,:auto,receiver,100.0,0.25,2.0,:junction)]
    system=normalize_river_connections(HydroSystem(;reservoirs,plants,generators,
        tunnels=[tunnel],junctions=Junction[],
        boundaries=[Boundary(name=:Sea,head=0.0)],rivers))
    operations=[OperationalSeries(object=:Tributary,attribute=:inflow,
        times=[0.0,1.0,2.0],values=[1.0,2.0,0.0]),
        OperationalSeries(object=:Main,attribute=:inflow,
        times=[0.0,1.0,2.0],values=[2.0,3.0,1.0])]
    ScheduleCase(name="mixed-source forward conservation",system=system,
        grid=[0.0,1.0,2.0,3.0],prices=zeros(3),operations=operations)
end

@testset "Mixed plant, tunnel, tributary and natural-inflow conservation" begin
    for (distributed,fixed_reach) in ((false,nothing),(true,nothing),
            (true,:Tributary),(true,:Main)), receiver in (:Sea,:Lower)
        c=mixed_forward_case(;distributed,fixed_reach,receiver)
        q=[3.0 2.0 4.0;2.0 3.0 1.0]
        gate=[0.1 0.2 0.05;0.0 0.0 0.0]
        z=simulate(c,q,gate)
        @test z["converged"]
        nat_upper=[1.0,2.0,0.0]
        nat_main=[2.0,3.0,1.0]
        for t in 1:3
            controlled_release=(20.0-nat_upper[t])*gate[1,t]
            expected_upper=z["V"][1,t]-0.0036*(q[1,t]+z["tunnel_q"][1,t]+controlled_release)
            expected_lower=z["V"][2,t]-0.0036*q[2,t]+
                (receiver==:Lower ? z["arrival_volume"][2,t] : 0.0)
            @test z["V"][1,t+1]≈expected_upper atol=2e-10
            @test z["V"][2,t+1]≈expected_lower atol=2e-10
            @test z["river_release"][1,t]≈controlled_release+nat_upper[t] atol=1e-10
            @test z["river_release"][2,t]≈sum(q[:,t])+z["tunnel_q"][1,t]+
                nat_main[t]+z["arrival_volume"][1,t]/0.0036 atol=1e-9
            inventory_change=sum(z["V"][:,t+1]-z["V"][:,t])+
                sum(z["transit"][:,t+1]-z["transit"][:,t])
            @test inventory_change≈0.0036*(nat_upper[t]+nat_main[t]-
                z["boundary_outflow"][t]) atol=2e-10
            upper_head=100.0+0.1*(z["V"][1,t]+z["V"][1,t+1])/2
            lower_head=60.0+0.1*(z["V"][2,t]+z["V"][2,t+1])/2
            @test z["power"][1,t]≈0.00981*q[1,t]*(upper_head-lower_head)*0.9 atol=1e-10
        end
        if !distributed
            # The tributary's pulse traverses two delays without being rebinned
            # to its first interval average before entering the second reach.
            expected=0.0036*(0.75*(5.0+z["tunnel_q"][1,1]+2.0)+0.25*2.9)
            @test z["arrival_volume"][2,1]≈expected atol=1e-10
        end
    end
end

@testset "Natural inflow alone supplies an internal reach source" begin
    reach=River(name=:Natural,source=:auto,target=:Sea,law=:junction,
        capacity=10.0,inflow=2.0,curves=RiverRouting.DelayCurve[],
        deterministic_delay=0.25,water_value=0.0,
        history_grid=[-1.0,0.0],history_release=[0.0])
    system=normalize_river_connections(HydroSystem(reservoirs=Reservoir[],
        junctions=Junction[],boundaries=[Boundary(name=:Sea,head=0.0)],
        tunnels=Tunnel[],plants=Plant[],generators=Generator[],rivers=[reach]))
    c=ScheduleCase(name="natural inflow only",system=system,grid=[0.0,1.0,2.0],
        prices=zeros(2))
    z=simulate(c,zeros(0,2),zeros(1,2))
    @test z["converged"]
    @test vec(z["river_release"])≈[2.0,2.0] atol=1e-12
    @test vec(z["arrival_volume"])≈0.0036*[1.5,2.0] atol=1e-12
    @test vec(z["transit"])≈0.0036*[0.0,0.5,0.5] atol=1e-12
end

@testset "Routed tunnel direction is explicit and never clipped" begin
    function outlet_system(source,target)
        reach=River(name=:Reach,source=:auto,target=:Sea,law=:junction,
            capacity=20.0,curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,
            water_value=0.0,history_grid=[-1.0,0.0],history_release=[0.0])
        normalize_river_connections(HydroSystem(reservoirs=Reservoir[],
            junctions=Junction[],boundaries=[Boundary(name=:High,head=100.0),
                Boundary(name=:Low,head=0.0),Boundary(name=:Sea,head=0.0)],
            tunnels=[Tunnel(name=:Tunnel,source=source,target=target,resistance=1.0,
                capacity=20.0,discharge_river=:Reach)],plants=Plant[],
            generators=Generator[],rivers=[reach]))
    end
    forward=forward_step(outlet_system(:High,:Low),Float64[],1.0,
        Float64[],[0.0],[0.0];current_transfer=ones(1,2))
    @test forward["converged"]
    @test only(forward["tunnel_q"])≈10.0 atol=1e-12
    @test only(forward["river_release"])≈10.0 atol=1e-12
    @test forward["boundary_net_inflow"]≈[-10.0,0.0,10.0] atol=1e-12
    reverse=forward_step(outlet_system(:Low,:High),Float64[],1.0,
        Float64[],[0.0],[0.0];current_transfer=ones(1,2))
    @test !reverse["converged"]
    @test only(reverse["tunnel_q"])≈-10.0 atol=1e-12
    @test reverse["bound_violations"]["tunnel_flow"]≈10.0 atol=1e-12
end
