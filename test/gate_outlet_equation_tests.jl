go_copy(x;kwargs...)=typeof(x)(;merge((; (k=>getfield(x,k) for k in fieldnames(typeof(x)))...),
    (;kwargs...))...)

function go_outlet_case(law,delta,opening;allow_dry=true)
    withdrawal=law==:weir ? 2max(delta,0)^1.5 : 2opening*sqrt(max(delta,0))
    lake=Reservoir(name=:Lake,z0=98.0+delta,slope=1.0,v0=2.0,vmin=1.5,vmax=2.5,
        inflow=withdrawal,water_value=1.0)
    river=River(name=:Outlet,source=:Lake,target=:Tail,law=law,coefficient=2.0,
        crest=100.0,allow_dry=allow_dry,capacity=10.0,inflow=0.3,
        curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,water_value=0.0)
    operations=[OperationalSeries(object=:Outlet,attribute=attribute,times=[0.0],values=[opening])
        for attribute in (:gate_min,:gate_max)]
    ScheduleCase(name="outlet thresholds",system=HydroSystem(reservoirs=[lake],rivers=[river],
        boundaries=[Boundary(name=:Tail,head=0.0)]),grid=[0.0,0.5,1.25],prices=zeros(2),operations=operations)
end

@testset "Orifice and weir thresholds retain natural river inflow" begin
    for (law,opening) in ((:orifice,0.0),(:orifice,0.5),(:orifice,1.0),(:weir,1.0)),
        delta in (-1.0,0.0,1.0)
        c=go_outlet_case(law,delta,opening)
        withdrawal=law==:weir ? 2max(delta,0)^1.5 : 2opening*sqrt(max(delta,0))
        x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),fill(opening,1,2))
        @test validate(c,x)["valid"]
        @test x.storage≈fill(2.0,1,3) atol=1e-10
        @test x.river_release≈fill(withdrawal+0.3,1,2) atol=1e-10
        @test x.arrival_volume≈reshape(0.0036.*diff(c.grid).*(withdrawal+0.3),1,:) atol=1e-10
        @test x.objective≈0.0 atol=1e-10
        @test replay_audit(c,x)["valid"]
        local_result=solve_case(c;warm=x,time_limit=30.0)
        @test local_result.status in ("LOCALLY_SOLVED","ALMOST_LOCALLY_SOLVED")
        @test local_result.accepted
        if local_result.solution!==nothing
            @test local_result.solution.river_release≈x.river_release atol=2e-5
            @test replay_audit(c,local_result.solution)["valid"]
        end
        result=solve(c;initial=x,time_limit=30.0)
        @test result.accepted
        if result.solution!==nothing
            @test result.solution.river_release≈x.river_release atol=2e-5
            @test result.solution.storage≈x.storage atol=1e-6
        end
    end
    # Without dry operation the complete static domain must stay above crest.
    @test_throws ErrorException OpenSHOP.validate_inputs(go_outlet_case(:orifice,0.0,0.5;allow_dry=false))
    @test OpenSHOP.validate_inputs(go_outlet_case(:orifice,1.0,0.5;allow_dry=false))
end

@testset "The optimizer can close a wet outlet while the reservoir fills" begin
    base=go_outlet_case(:orifice,1.0,1.0)
    lake=go_copy(only(base.system.reservoirs);water_value=10000.0)
    c=go_copy(base;operations=OperationalSeries[],system=go_copy(base.system;reservoirs=[lake]))
    initial=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),ones(1,2))
    expected=0.0036*2.0*(last(c.grid)-first(c.grid))
    for (method,result) in (("Local dispatch",solve_case(c;warm=initial,time_limit=30.0)),
        ("Global solve",solve(c;initial,time_limit=30.0,relative_gap=1e-6)))
        @testset "$method" begin
            @test result.accepted
            if result.solution!==nothing
                @test result.solution.objective≈10000.0*expected atol=1e-3
                @test maximum(abs,result.solution.gate)<2e-5
                @test result.solution.storage[1,end]≈2.0+expected atol=2e-7
                @test result.solution.river_release≈fill(0.3,1,2) atol=2e-5
                @test replay_audit(c,result.solution)["valid"]
            end
        end
    end
end

@testset "A filling reservoir crosses the dry-outlet threshold" begin
    for (law,opening) in ((:orifice,0.5),(:weir,1.0))
        base=go_outlet_case(law,-0.1,opening)
        lake=go_copy(only(base.system.reservoirs);inflow=100.0)
        c=go_copy(base;system=go_copy(base.system;reservoirs=[lake]),
            grid=[0.0,0.1,0.5,1.0],prices=zeros(3))
        volumes=[2.0];withdrawals=Float64[]
        for dt in diff(c.grid)
            previous=last(volumes);lo,hi=previous,2.5
            withdrawal(v)=law==:weir ? 2max(97.9+(previous+v)/2-100.0,0)^1.5 :
                2opening*sqrt(max(97.9+(previous+v)/2-100.0,0))
            for _ in 1:80
                midpoint=(lo+hi)/2
                residual=midpoint-previous-0.0036dt*(100.0-withdrawal(midpoint))
                residual>0 ? (hi=midpoint) : (lo=midpoint)
            end
            volume=(lo+hi)/2
            push!(volumes,volume);push!(withdrawals,withdrawal(volume))
        end
        x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),fill(opening,1,3))
        @test vec(x.storage)≈volumes atol=1e-9
        @test vec(x.river_release)≈withdrawals.+0.3 atol=1e-8
        @test x.head[1,1]<100.0<x.head[1,2]<x.head[1,3]
        @test replay_audit(c,x)["valid"]
        local_result=solve_case(c;warm=x,time_limit=30.0)
        @test local_result.status in ("LOCALLY_SOLVED","ALMOST_LOCALLY_SOLVED")
        @test local_result.accepted
        if local_result.solution!==nothing
            @test vec(local_result.solution.storage)≈volumes atol=1e-6
        end
        result=solve(c;initial=x,time_limit=30.0)
        @test result.accepted
        if result.solution!==nothing
            @test vec(result.solution.storage)≈volumes atol=1e-6
        end
    end
end

@testset "Tunnel position and resistance give signed flow in either direction" begin
    for curve in (nothing,TableCurve([0.0,1.0],[6.0,0.0])),opening in (0.0,0.5,1.0),
        delta in (-20.0,0.0,20.0)
        tunnel=Tunnel(name=:Link,source=:A,target=:B,resistance=2.0,capacity=10.0,
            opening=opening,gate_positions=[0.0,0.5,1.0],gate_loss_factor_curve=curve)
        c=ScheduleCase(name="signed gate law",system=HydroSystem(tunnels=[tunnel],
            boundaries=[Boundary(name=:A,head=100.0+delta),Boundary(name=:B,head=100.0)]),
            grid=[0.0,1.0],prices=[0.0],operations=[OperationalSeries(object=:Link,
                attribute=:opening,times=[0.0],values=[opening])])
        resistance=opening==0 ? Inf : curve===nothing ? 2.0/opening : 2.0+6.0*(1-opening)
        expected=opening==0 ? 0.0 : sign(delta)*sqrt(abs(delta)/resistance)
        x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1);
            tunnel_opening=fill(opening,1,1))
        @test validate(c,x)["valid"]
        @test only(x.tunnel_q)≈expected atol=1e-9
        @test x.head[:,1]==[100.0+delta,100.0]
        @test replay_audit(c,x)["valid"]
        result=solve(c;initial=x,time_limit=20.0)
        @test result.accepted
        if result.solution!==nothing
            @test only(result.solution.tunnel_q)≈expected atol=1e-5
        end
    end
end

@testset "Turbine outlet elevation is continuous through submergence" begin
    for receiver in (99.999,100.0,100.001)
        c=outlet_test_case(receiver)
        operations=[OperationalSeries(object=:Unit,attribute=:discharge,times=[0.0],values=[2.0]),
            OperationalSeries(object=:Downstream,attribute=:gate_min,times=[0.0],values=[0.5]),
            OperationalSeries(object=:Downstream,attribute=:gate_max,times=[0.0],values=[0.5])]
        c=go_copy(c;operations)
        expected=0.00981*2.0*(202.0-max(receiver,100.0)-0.2)*0.9
        x=dispatch_from_controls(c,ones(Int,1,1),fill(2.0,1,1),fill(0.5,1,1))
        @test only(x.power)≈expected atol=1e-10
        @test x.storage≈fill(2.0,2,2) atol=1e-10
        @test x.objective≈50.0*expected atol=1e-8
        @test replay_audit(c,x)["valid"]
        result=solve(c;initial=x,time_limit=30.0)
        @test result.accepted
        if result.solution!==nothing
            @test only(result.solution.power)≈expected atol=2e-5
            @test replay_audit(c,result.solution)["valid"]
        end
    end
end

@testset "Unsupported outlet data are rejected instead of ignored" begin
    c=go_outlet_case(:orifice,1.0,0.5)
    for (objects,attribute) in (("rivers","downstream_elevation"),("rivers","width_depth_curve"))
        data=case_dict(c);data[objects][1][attribute]=1.0
        @test_throws ArgumentError case_from_dict(data)
    end
    tunnel=Tunnel(name=:Link,source=:A,target=:B,resistance=1.0,capacity=10.0)
    tc=ScheduleCase(name="submerged tunnel",system=HydroSystem(tunnels=[tunnel],
        boundaries=[Boundary(name=:A,head=100.0),Boundary(name=:B,head=90.0)]),
        grid=[0.0,1.0],prices=[0.0])
    for attribute in ("start_height","end_height","gate_opening_curve")
        data=case_dict(tc);data["tunnels"][1][attribute]=1.0
        @test_throws ArgumentError case_from_dict(data)
    end
    uncontrolled=go_copy(go_outlet_case(:weir,1.0,1.0);operations=[
        OperationalSeries(object=:Outlet,attribute=:gate_max,times=[0.0],values=[0.5])])
    @test_throws ArgumentError OpenSHOP.validate_inputs(uncontrolled)
    # A junction cannot substitute for a reservoir-dependent outlet.
    river=go_copy(only(c.system.rivers);source=:Merge)
    bad=go_copy(c;system=go_copy(c.system;rivers=[river],river_junctions=[RiverJunction(name=:Merge)]))
    @test_throws ErrorException OpenSHOP.validate_inputs(bad)
end
