using Random

@testset "Batched arrival observations retain scalar rates and jump ownership" begin
    grid=[1.5,1.75,2.0]
    flows=[4.0,7.0]
    times=sort!(unique!(vcat([0.0,1.5,2.05,2.3,10.0],grid .+0.3)))
    for delay in (nothing,0.0,0.3),side in (:left,:right)
        river=River(name=:Reach,source=:Lake,target=:Sea,capacity=10.0,water_value=0.0,
            deterministic_delay=delay,curves=delay===nothing ?
                [RiverRouting.DelayCurve(0.0,[0.0,0.2,0.7],[0.4,0.6]),
                 RiverRouting.DelayCurve(10.0,[0.1,0.3,0.5],[0.7,0.3])] : RiverRouting.DelayCurve[])
        rates=OpenSHOP.cohort_arrival_rates(river,grid,flows,times;side)
        expected=[OpenSHOP.cohort_point_arrival(river,grid,flows,t;side) for t in times]
        @test rates≈expected atol=1e-12
        @test isempty(OpenSHOP.cohort_arrival_rates(river,grid,flows,Float64[];side))
    end
end

@testset "Cohort operators reproduce routed volumes and physical transit" begin
    rng=MersenneTwister(7183)
    for trial in 1:8
        release=vcat(-1.7,-1.7 .+ cumsum(0.03 .+ rand(rng,15)))
        observation=vcat(-2.0,-2.0 .+ cumsum(0.07 .+ rand(rng,21)))
        curves=[RiverRouting.DelayCurve(0.0,[0.0,0.19,0.81],[0.0,1.0]),
            RiverRouting.DelayCurve(20.0,[0.09,0.4,1.7],[0.7,0.3]),
            RiverRouting.DelayCurve(60.0,[0.0,0.08,0.43],[1.0,0.0])]
        for (delay,single) in ((nothing,false),(nothing,true),(0.0,true),(0.17,true),(20.0,true))
            r=River(name=:Reach,source=:Lake,target=:Sea,capacity=60.0,water_value=0.0,
                deterministic_delay=delay,curves=delay===nothing ?
                    (single ? curves[2:2] : curves) : RiverRouting.DelayCurve[])
            q=60 .* rand(rng,length(release)-1)
            B=OpenSHOP.RiverTransfer(r,observation,release)
            arrivals=[sum(0.0036*(release[j+1]-release[j])*
                OpenSHOP.released_fraction(r,q[j],fractions)
                for (j,fractions) in OpenSHOP.arrival_terms(B,t);init=0.0)
                for t in 1:length(observation)-1]
            expected=OpenSHOP.route_volumes(r,observation,release,q)
            @test arrivals≈expected atol=1e-12
            # Chronological propagation must give the same water as row evaluation.
            streamed=zeros(length(arrivals))
            for j in eachindex(q),(t,fractions) in OpenSHOP.release_terms(B,j)
                streamed[t]+=0.0036*(release[j+1]-release[j])*
                    OpenSHOP.released_fraction(r,q[j],fractions)
            end
            @test streamed≈expected atol=1e-12
            terminal=sum(0.0036*(release[j+1]-release[j])*
                OpenSHOP.released_fraction(r,q[j],view(B.terminal,:,j)) for j in eachindex(q))
            @test terminal≈OpenSHOP.remaining_volume(r,last(observation),release,q) atol=1e-12
            released=0.0036sum(q[j]*max(0.0,min(last(observation),release[j+1])-release[j])
                for j in eachindex(q))
            @test sum(arrivals)+terminal≈released atol=1e-12
            for (a,b) in ((-3.0,-2.0),(-0.13,0.83),(0.19,3.17),(10.0,11.0))
                row=OpenSHOP.volume_transfer_row(r,release,a,b)
                volume=sum(0.0036*(release[k+1]-release[k])*
                    OpenSHOP.released_fraction(r,q[k],view(row.fractions,:,j))
                    for (j,k) in enumerate(row.cohorts);init=0.0)
                @test volume≈only(OpenSHOP.route_volumes(r,[a,b],release,q)) atol=1e-12
            end
            lag=delay===nothing ? 0.19 : delay
            for time in (first(release),release[4],last(release),0.23,20.0,
                release[4]+lag,release[end-1]+lag),side in (:left,:right)
                row=OpenSHOP.point_transfer_row(r,release,time;side)
                arrival=sum(OpenSHOP.released_fraction(r,q[k],view(row.fractions,:,j))
                    for (j,k) in enumerate(row.cohorts);init=0.0)
                @test arrival≈OpenSHOP.cohort_point_arrival(r,release,q,time;side) atol=1e-11
            end
        end
    end
end

@testset "Atomic boundaries and small terminal inventories retain their meaning" begin
    r=River(name=:Reach,source=:Lake,target=:Sea,capacity=10.0,water_value=0.0,
        deterministic_delay=0.3,curves=RiverRouting.DelayCurve[])
    grid=[1.5,1.75,2.0]
    q=[4.0,7.0]
    for (side,expected) in ((:left,4.0),(:right,7.0))
        row=OpenSHOP.point_transfer_row(r,grid,2.05;side)
        actual=sum(OpenSHOP.released_fraction(r,q[k],view(row.fractions,:,j))
            for (j,k) in enumerate(row.cohorts);init=0.0)
        @test actual==expected
    end
    r=River(name=:Reach,source=:Lake,target=:Sea,capacity=1e9,water_value=0.0,
        curves=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])])
    time=2.0-1e-6
    B=OpenSHOP.RiverTransfer(r,[0.0,time],[0.0,1.0])
    inventory=0.0036*1e9*only(B.terminal)
    @test inventory≈0.0036*1e9*(2.0-time)^2/2 rtol=1e-12
end

@testset "Observation expressions preserve history and release-time weights" begin
    curves=[RiverRouting.DelayCurve(0.0,[0.1,0.7],[1.0]),
        RiverRouting.DelayCurve(30.0,[0.05,0.4],[1.0]),
        RiverRouting.DelayCurve(60.0,[0.0,0.2],[1.0])]
    r=River(;name=:Reach,source=:Lake,target=:Sea,law=:controlled,capacity=60.0,water_value=0.0,
        curves,history_grid=[-2.0,-0.21,0.0],history_release=[17.0,43.0])
    c=ScheduleCase(name="transfer observations",system=HydroSystem(reservoirs=[Reservoir(name=:Lake,v0=10.0,vmin=0.0,vmax=20.0,water_value=0.0)],
        rivers=[r],boundaries=[Boundary(name=:Sea,head=0.0)]),
        grid=[0.0,0.13,0.7,1.4],prices=zeros(3),routing_grid=[0.0,0.13,0.23,0.7,1.0,1.4])
    model=Model()
    release=@variable(model,[1:1,1:3],lower_bound=0.0,upper_bound=60.0)
    routing=OpenSHOP._dispatch_transport_expressions(model,c,release)
    for q in ([0.0,30.0,60.0],[17.0,43.0,3.0])
        set_start_value.(release,permutedims(q))
        for (a,b) in ((0.0,1.4),(0.13,0.83),(0.21,0.9))
            actual=value(start_value,OpenSHOP.routing_volume_expression(model,c,routing,1,a,b))
            expected=only(OpenSHOP.route_volumes(r,[a,b],c.grid,q))+
                only(OpenSHOP.route_volumes(r,[a,b],r.history_grid,r.history_release))
            @test actual≈expected atol=1e-12
        end
        for time in (0.0,0.21,0.7,1.4),side in (:left,:right)
            expression=OpenSHOP.routing_point_expression(model,c,routing,1,time;side)
            @test value(start_value,expression)≈point_arrival(r,c.grid,q,time;side) atol=1e-11
        end
        @test value(start_value,routing.terminal[1])≈OpenSHOP.remaining_volume(r,last(c.grid),c.grid,q)+
            OpenSHOP.remaining_volume(r,last(c.grid),r.history_grid,r.history_release) atol=1e-12
    end
end

@testset "Narrow travel responses retain small late-window volumes" begin
    flow=17.0
    for origin in (0.0,1000.0), hours in (6.0,168.0), width in (1e-3,1e-6)
        grid=collect(range(origin,origin+hours;length=65))
        river=River(name=:Reach,source=:Lake,target=:Sea,capacity=60.0,
            water_value=0.0,history_grid=[origin-8.0,origin],history_release=[0.0],
            curves=[RiverRouting.DelayCurve(q,[lag,lag+width],[1.0])
                for (q,lag) in zip((0.0,30.0,60.0),(0.1,0.3,0.5))])
        system=HydroSystem(reservoirs=[Reservoir(name=:Lake,v0=50.0,
            vmin=0.0,vmax=100.0,water_value=0.0)],rivers=[river],
            boundaries=[Boundary(name=:Sea,head=0.0)])
        c=ScheduleCase(;name="steady delayed flow",system,grid,prices=zeros(64))
        model=Model()
        @variable(model,0.0<=release[1:1,1:64]<=60.0)
        set_start_value.(release,flow)
        routing=OpenSHOP._dispatch_transport_expressions(model,c,release)
        for (variable,expression) in routing.starts
            set_start_value(variable,value(start_value,expression))
        end
        # A constant entry stream has reached steady flow long before these
        # observations. The expected volume needs no routing calculation.
        a=origin+0.9hours
        for duration in (0.25,1e-4,1e-8)
            b=a+duration
            expression=OpenSHOP.routing_volume_expression(model,c,routing,1,a,b)
            @test value(start_value,expression)≈0.0036flow*(b-a) atol=1e-13
        end
        for side in (:left,:right)
            expression=OpenSHOP.routing_point_expression(model,c,routing,1,a;side)
            @test value(start_value,expression)≈flow atol=1e-5
        end
        mean_delay=(13/30)*(0.1+width/2)+(17/30)*(0.3+width/2)
        @test value(start_value,only(routing.terminal))≈0.0036flow*mean_delay atol=1e-12
    end
end
