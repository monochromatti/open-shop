function controlled_outlet_case(;capacities=[5.0,3.0],openings=[0.4,1.0])
    reservoirs=[
        Reservoir(name=:Upper,z0=100.0,slope=1.0,v0=2.0,vmin=0.5,vmax=4.0,water_value=0.0),
        Reservoir(name=:Lower,z0=50.0,slope=1.0,v0=1.0,vmin=0.5,vmax=4.0,water_value=1000.0),
    ]
    outlet=River(name=:Outlet,source=:Upper,target=:Lower,law=:controlled,
        capacity=5.0,inflow=1.0,water_value=0.0,curves=RiverRouting.DelayCurve[],
        deterministic_delay=0.5,history_grid=[-1.0,0.0],history_release=[0.0])
    operations=[OperationalSeries(object=:Outlet,attribute=attribute,
        times=[0.0,1.0],values=values) for (attribute,values) in
        ((:capacity,capacities),(:inflow,[1.0,2.0]),(:gate_max,openings))]
    ScheduleCase(name="controlled outlet economics",grid=[0.0,1.0,2.0],prices=zeros(2),
        system=HydroSystem(;reservoirs,rivers=[outlet]),operations=operations)
end

@testset "Controlled outlets optimize release with natural inflow and dated limits" begin
    for (capacities,openings,expected) in (
        ([5.0,3.0],[0.4,1.0],[2.6,3.0]),
        ([1.0,2.0],[0.4,1.0],[1.0,2.0]),
        ([1.0,3.0],[0.4,1.0],[1.0,3.0]),
        ([5.0,3.0],[0.0,0.0],[1.0,2.0]),
    )
        c=controlled_outlet_case(;capacities,openings)
        # All water arriving in Lower earns 1000 per Mm³. The half-hour delay
        # delivers the whole first release and half the second by the horizon.
        value=3.6*(expected[1]+expected[2]/2)
        for result in (solve_case(c;time_limit=20.0),
            solve_verified(c;time_limit=20.0,max_refinements=0))
            @test result["accepted"]
            if result["accepted"]
                x=result["solution"]
                @test vec(x.river_release)≈expected atol=1e-4
                @test x.objective≈value atol=2e-4
                @test validate(c,x)["valid"]
                @test replay_audit(c,x)["valid"]
            end
        end
    end
    @test_throws ArgumentError solve_case(controlled_outlet_case(;capacities=[0.5,3.0]))
end
