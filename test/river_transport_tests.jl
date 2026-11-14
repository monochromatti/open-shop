@testset "Transport resolution preserves control times, windows and restart releases" begin
    fixture=mixed_source_fixture(;distributed=true)
    base=fixture.case
    for factor in (1,3,7)
        mesh=OpenSHOP.refined_grid(base.grid;factor)
        c=OpenSHOP._river_replace(base;routing_grid=mesh)
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        x=dispatch_from_controls(c,fixture.u,fixture.q,fixture.gates)
        @test x.grid==base.grid
        @test x.commitment==fixture.u
        @test x.generator_q==fixture.q
        @test x.diagnostics["validation"]["valid"]
        mixed_source_conservation(c,simulate(c,fixture.q,fixture.gates))
        domains=OpenSHOP._global_capacity_bounds(c,OpenSHOP.transport_discretization(c))
        @test all(x.arrival_volume.>=domains.arrival_lower.-1e-10)
        @test all(x.arrival_volume.<=domains.arrival_upper.+1e-10)
        @test all(x.storage.>=domains.lower.-1e-10)
        @test all(x.storage.<=domains.upper.+1e-10)
        restarted=restart_case(c,x,2.0)
        rest=dispatch_from_controls(restarted,fixture.u[:,4:end],fixture.q[:,4:end],fixture.gates[:,4:end])
        @test rest.diagnostics["validation"]["valid"]
        @test rest.storage≈x.storage[:,4:end] atol=1e-9
        @test rest.arrival_volume≈x.arrival_volume[:,4:end] atol=1e-9
        @test rest.terminal_transit≈x.terminal_transit atol=1e-9
        if factor==1
            @test !validate_replay(c,x)["valid"] # Conservation alone cannot establish timing accuracy.
        end
    end
    for mesh in ([0.0,last(base.grid)],[-1.0;base.grid],[base.grid;Inf])
        invalid=OpenSHOP._river_replace(base;routing_grid=mesh)
        @test_throws ArgumentError OpenSHOP.validate_inputs(invalid)
    end
end

@testset "Prescribed transport retains its mesh under hydraulic refinement" begin
    fixture=mixed_source_small_fixture(;mixed=true)
    c=fixture.case
    for factor in (1,3)
        grid=OpenSHOP.refined_grid(c.grid;factor)
        z=simulate(c,fixture.q,fixture.gates;grid)
        @test z.converged
        physical=OpenSHOP._river_replace(c;grid,prices=zeros(length(grid)-1))
        routing_grid=OpenSHOP.river_routing_grid(c;physical_grid=grid)
        # Independently route the specified source rates on the declared mesh.
        release=zeros(length(c.system.rivers),length(grid)-1)
        for (d,r) in enumerate(c.system.rivers),t in axes(release,2)
            r.law==:controlled || continue
            k=searchsortedlast(c.grid,(grid[t]+grid[t+1])/2)
            inflow=OpenSHOP.opaverage(c,r.name,:inflow,grid[t],grid[t+1],r.inflow)
            capacity=OpenSHOP.opaverage(c,r.name,:capacity,grid[t],grid[t+1],r.capacity)
            release[d,t]=inflow+max(0,capacity-inflow)*fixture.gates[d,k]
        end
        routed=OpenSHOP._route_network(physical,release;grid=routing_grid,
            generator_q=z.generator_q)
        edges=[searchsortedfirst(routing_grid,t) for t in grid]
        arrivals=[sum(view(routed.arrival_volume,d,edges[t]:(edges[t+1]-1)))
            for d in axes(release,1),t in axes(release,2)]
        @test z.arrival_volume≈arrivals atol=1e-11
        @test z.transit≈routed.transit[:,edges] atol=1e-11
        mixed_source_conservation(physical,z)
    end
end

@testset "River junction checks retain one-sided waves and pointwise restrictions" begin
    fixture=mixed_source_small_fixture(;mixed=true)
    for offset in (0.0,1000.0)
        base=fixture.case
        rivers=[OpenSHOP._river_replace(r;history_grid=r.history_grid.+offset,
            arrival_policy=:pointwise,min_arrival=0.3r.capacity) for r in base.system.rivers]
        operations=[OpenSHOP._river_replace(z;times=z.times.+offset) for z in base.operations]
        for r in rivers
            push!(operations,OperationalSeries(object=r.name,attribute=:capacity,
                times=offset.+[0.0,0.37,1.05,2.2],values=r.capacity.*[1.0,0.15,0.8,1.0]))
        end
        c=OpenSHOP._river_replace(base;grid=base.grid.+offset,operations,
            system=OpenSHOP._river_replace(base.system;rivers))
        initial=simulate(c,fixture.q,fixture.gates)
        routed=OpenSHOP.routed_schedule(c,initial.river_release;generator_q=fixture.q)
        residuals=OpenSHOP.river_arrival_residuals(c,routed)
        for (d,r) in enumerate(rivers)
            series=routed.release_series[d]
            expected=maximum(series.release[k]-OpenSHOP.opinterval(c,r.name,:capacity,t,r.capacity)
                for k in eachindex(series.release),t in eachindex(c.prices)
                if min(series.grid[k+1],c.grid[t+1])>max(series.grid[k],c.grid[t]);init=0.0)
            upstream=findall(e->e.target==r.source,rivers)
            if r.law==:junction && !isempty(upstream)
                times=OpenSHOP.junction_arrival_times(c,[z.grid for z in routed.release_series],upstream)
                for time in times
                    sides=time==first(c.grid) ? (:right,) :
                        time==last(c.grid) ? (:left,) : (:left,:right)
                    for side in sides
                        t=OpenSHOP.river_control_interval(c,time,side)
                        rate=routed.injections[d,t]+r.allocation*sum(point_arrival(
                            rivers[i],routed.release_series[i].grid,routed.release_series[i].release,time;side)
                            for i in upstream;init=0.0)
                        expected=max(expected,rate-OpenSHOP.opvalue(c,r.name,:capacity,time,r.capacity;side))
                    end
                end
            end
            @test residuals.capacity[d]≈expected atol=1e-11
            # This must still use the target reach's releases after inspecting
            # its upstream reaches, including nonzero minimum arrival laws.
            expected_arrival=OpenSHOP.arrival_requirement_violation(c,r,series.grid,series.release)
            @test residuals.pointwise_arrival[d]≈expected_arrival atol=1e-11
        end
        @test maximum(residuals.capacity)>0
        @test maximum(residuals.pointwise_arrival)>0
    end
end

@testset "Recorded river flows give feasible cold proposals under tight storage limits" begin
    reservoirs=[Reservoir(name=:Lake,z0=100.0,slope=0.1,v0=10.0,
        vmin=9.0,vmax=11.0,water_value=0.0),
        Reservoir(name=:Receiver,z0=50.0,slope=0.1,v0=10.0,
            vmin=10.0,vmax=10.00721,water_value=10.0)]
    river=River(name=:Link,source=:Lake,target=:Receiver,law=:controlled,
        capacity=10.0,gate_min=0.2,water_value=0.0,initial_release=2.0,
        history_grid=[-4.0,0.0],history_release=[2.0],
        curves=[RiverRouting.DelayCurve(0.0,[2.0,3.0],[1.0]),
            RiverRouting.DelayCurve(10.0,[0.0,1.0],[1.0])])
    c=ScheduleCase(name="recorded river flow",system=HydroSystem(;reservoirs,rivers=[river]),
        grid=[0.0,1.0],prices=[0.0],operations=[OperationalSeries(object=:Link,
            attribute=:gate_max,times=[0.0],values=[0.2])])
    # Continuing the historical two m³/s gives exactly two m³/s of arrivals.
    # A faster nominal response would falsely overfill the receiving reservoir.
    proposal=OpenSHOP.propose_commitment(c;time_limit=10.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        result=solve_verified(c;u=proposal["u"],time_limit=15.0,
            max_refinements=0,refine_grid=false)
        @test result["accepted"]
        if result["accepted"]
            x=OpenSHOP._schedule_input(result["solution"];case=c)
            @test x.river_release[1,1]≈2.0 atol=1e-5
            @test x.arrival_volume[1,1]≈0.0072 atol=1e-8
            @test x.storage[2,end]≈10.0072 atol=1e-8
        end
    end
end

@testset "Distributed arrivals satisfy point and window laws in local and global optimization" begin
    fixture=mixed_source_small_fixture(;mixed=true)
    base=fixture.case
    initial=dispatch_from_controls(base,fixture.u,fixture.q,fixture.gates)
    routed=OpenSHOP.routed_schedule(base,initial.river_release;
        generator_q=initial.generator_q,tunnel_q=initial.tunnel_q)
    main=only(findall(r->r.name==:MainRiver,base.system.rivers))
    series=routed.release_series[main]
    minimum_rate=minimum_arrival(base.system.rivers[main],series.grid,series.release)["value"]
    @test minimum_rate>0.0
    rivers=[r.name==:MainRiver ? OpenSHOP._river_replace(r;arrival_policy=:pointwise,
        min_arrival=0.9minimum_rate) : r.name==:ReceiverOutlet ?
        OpenSHOP._river_replace(r;arrival_window_grid=[0.0,0.5,2.5],min_arrival=2.9) : r
        for r in base.system.rivers]
    c=OpenSHOP._river_replace(base;system=OpenSHOP._river_replace(base.system;rivers))
    x=dispatch_from_controls(c,fixture.u,fixture.q,fixture.gates)
    @test x.diagnostics["validation"]["valid"]
    @test validate_replay(c,x)["valid"]
    local_result=solve_verified(c;u=fixture.u,warm=x,time_limit=20.0,
        max_refinements=0,refine_grid=false)
    @test local_result["accepted"]
    global_result=OpenSHOP.solve(c;initial=x,time_limit=30.0,relative_gap=1e-4)
    @test global_result.accepted
    @test global_result.global_certificate
    @test global_result.feasible_lower_bound>=x.objective-1e-5
    @test validate(c,global_result.solution)["valid"]
    @test validate_replay(c,global_result.solution)["valid"]
end

@testset "Transport acceptance checks operating peaks between control boundaries" begin
    reservoirs=[Reservoir(name=name,z0=z,v0=10.0,vmin=9.0,vmax=11.0,
        slope=0.1,water_value=name==:Receiver ? 10.0 : 0.0) for (name,z) in ((:Lake,100.0),(:Receiver,50.0))]
    rivers=[River(name=:Source,source=:Lake,target=:Downstream,law=:controlled,
        capacity=10.0,water_value=0.0,curves=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])]),
        River(name=:Downstream,source=:auto,target=:Receiver,law=:junction,
            capacity=20.0,water_value=0.0,deterministic_delay=0.0,curves=RiverRouting.DelayCurve[])]
    c=ScheduleCase(name="junction peak",system=HydroSystem(;reservoirs,rivers),
        grid=[0.0,1.0,2.0],routing_grid=[0.0,1.0,2.0],prices=zeros(2),
        operations=[OperationalSeries(object=:Downstream,attribute=:capacity,
            times=[0.0],values=[5.5]),
            OperationalSeries(object=:Source,attribute=:gate_max,times=[1.0],values=[0.0])])
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),[1.0 0.0;0.0 0.0])
    @test !x.diagnostics["validation"]["valid"]
    @test maximum(x.river_release[2,:])<5.5
    # The source's triangular arrival has mean five on both control intervals,
    # yet reaches ten at their common boundary. Zero downstream delay makes the
    # coarse cumulative deliveries exact, hiding the capacity violation.
    validation=validate_transport(c,x)
    @test validation["converged"]
    @test validation["cumulative_schedule_discrepancy_Mm3"]<1e-10
    @test !validation["operating_valid"]
    @test !validation["accepted"]
    @test validation["violations"]["routed_release_Downstream"]>4.0
    expected=0.0036*5.5*10.0 # All water in the one-hour pulse arrives by the horizon.
    u=zeros(Int,0,2)
    local_result=solve_verified(c;u,time_limit=15.0,max_refinements=0,refine_grid=false)
    @test local_result["accepted"]
    @test local_result["solution"]["objective"]>=expected-1e-3
    result=OpenSHOP.solve(c;time_limit=30.0,relative_gap=1e-5,preparation_time_limit=0.0)
    @test result.accepted && result.global_certificate
    @test result.solution.river_release[1,1]≈5.5 atol=1e-4
    @test result.feasible_lower_bound≈expected atol=1e-5
    @test result.global_bound≈expected atol=1e-5

end
