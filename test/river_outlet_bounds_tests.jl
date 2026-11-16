function outlet_bound_case(;law=:controlled,curve=nothing,minimum=0.0,
    gate_min=0.0,inflow=0.4,policy=:interval_average,operations=OperationalSeries[])
    lake=Reservoir(name=:Lake,z0=100.0,slope=5.0,v0=0.02,vmin=0.005,
        vmax=0.04,water_value=100_000.0)
    river=River(name=:Outlet,source=:Lake,target=:Sea,law=law,capacity=5.0,
        coefficient=2.0,crest=99.0,discharge_curve=curve,inflow=inflow,
        min_arrival=minimum,gate_min=gate_min,arrival_policy=policy,water_value=0.0,
        curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,history_grid=[-1.0,0.0],history_release=[inflow])
    system=HydroSystem(reservoirs=[lake],rivers=[river],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    ScheduleCase(name="reservoir outlet bounds",system=system,grid=[0.0,1.0,2.0],
        prices=zeros(2),operations=operations)
end

@testset "Outlet laws preserve nonnegative discharge when gates close" begin
    cases=[outlet_bound_case(),outlet_bound_case(law=:orifice),
        outlet_bound_case(law=:orifice,curve=TableCurve([100.0,101.0],[5.0,10.0])),
        outlet_bound_case(law=:weir,curve=TableCurve([100.0,100.2,101.0],[0.0,0.0,10.0])),
        outlet_bound_case(law=:orifice,inflow=0.0),
        outlet_bound_case(law=:weir,inflow=0.0,curve=TableCurve([100.0,100.2,101.0],[0.0,0.0,10.0]))]
    for c in cases
        result=solve_case(c;u=zeros(Int,0,2),time_limit=45.0)
        @test result.accepted
        if result.accepted
            x=result.solution
            @test x.river_release≈fill(only(c.system.rivers).inflow,1,2) atol=1e-5
            @test x.storage≈fill(0.02,1,3) atol=1e-7
            @test validate_replay(c,x;grid=OpenSHOP.refined_grid(c.grid;factor=4))["valid"]
        end
    end
end

@testset "Gate and arrival restrictions still set the minimum release" begin
    for policy in (:interval_average,:pointwise)
        c=outlet_bound_case(minimum=1.4,policy=policy)
        result=solve_case(c;u=zeros(Int,0,2),time_limit=20.0)
        @test result.accepted
        if result.accepted
            @test result.solution.river_release≈fill(1.4,1,2) atol=1e-5
            @test validate_replay(c,result.solution;grid=OpenSHOP.refined_grid(c.grid;factor=4))["valid"]
        end
    end
    c=outlet_bound_case(gate_min=0.25)
    result=solve_case(c;u=zeros(Int,0,2),time_limit=20.0)
    @test result.accepted
    if result.accepted
        @test result.solution.river_release≈fill(0.4+(5.0-0.4)*0.25,1,2) atol=1e-5
    end
    c=outlet_bound_case(policy=:pointwise)
    inset=solve_case(c;u=zeros(Int,0,2),arrival_margin=0.6,time_limit=20.0)
    @test inset.accepted
    if inset.accepted
        @test inset.solution.river_release≈fill(0.6,1,2) atol=1e-5
    end
end

@testset "A closed river capacity still permits a zero-flow schedule" begin
    operation=OperationalSeries(object=:Outlet,attribute=:capacity,times=[0.0],values=[0.0])
    c=outlet_bound_case(inflow=0.0,operations=[operation])
    result=solve_case(c;u=zeros(Int,0,2),time_limit=20.0)
    @test result.accepted
    if result.accepted
        @test result.solution.river_release≈zeros(1,2) atol=1e-8
        @test result.solution.storage≈fill(0.02,1,3) atol=1e-8
    end
end
