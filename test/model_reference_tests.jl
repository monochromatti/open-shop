function check_reference_lift(c, seed; kwargs...)
    b=OpenSHOP._build_global_dispatch(c; kwargs...)
    audit=OpenSHOP._lift_start!(b,c,seed)
    @test audit["valid"]
    @test audit["objective"]≈seed["objective"] atol=1e-6
    b
end

@testset "Physical starts preserve gates, reverse flows and mixed equipment curves" begin
    for tail in (0.0,110.0)
        c=tg_flow_case(;tail)
        opening=reshape([0.0,0.5,1.0],1,3)
        seed=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);
            tunnel_opening=opening)
        @test seed["validation"]["valid"]
        check_reference_lift(c,seed;joint=true)
        check_reference_lift(c,seed;joint=true,tunnel_opening=opening)
    end
    for direction in (-1.0,1.0)
        c=delta_channel_case(;direction)
        family=LevelFlowCurve([100.0,110.0],[
            TableCurve([-10.0,0.0,10.0],[-20.0,0.0,20.0]),
            TableCurve([-10.0,-2.0,0.0,10.0],[-40.0,-8.0,0.0,40.0]),
        ])
        river=tg_copy(only(c.system.rivers);delta_level_ref_up_flow_curve=family)
        c=tg_copy(c;system=tg_copy(c.system;rivers=[river]))
        seed=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(1,3))
        check_reference_lift(c,seed;joint=true)
    end
    c=outlet_test_case(80.0)
    river=tg_copy(only(c.system.rivers);deterministic_delay=nothing,
        curves=multireference_curves())
    c=tg_copy(c;system=tg_copy(c.system;rivers=[river]))
    seed=dispatch_from_controls(c,ones(Int,1,1),fill(2.0,1,1),fill(0.5,1,1))
    check_reference_lift(c,seed;joint=true)

    # Shared plant head, cubic efficiency, electrical table, reservoir vertex
    # level ramps, outlet floor and a river that crosses the dry/wet boundary.
    c=outlet_test_case(100.0)
    source=tg_copy(first(c.system.reservoirs);
        level_curve=TableCurve([1.0,2.0,3.0],[201.0,202.0,203.0]),
        level_ramp_up=10.0,level_ramp_down=10.0)
    table=TurbineTable([20.0,100.0,150.0],[1.0,3.0,5.0],
        [0.82 0.85 0.86;0.88 0.92 0.93;0.86 0.89 0.90],
        fill(1.0,3),fill(5.0,3);interpolation=:pchip_discharge)
    unit=tg_copy(only(c.system.generators);turbine_table=table,
        generator_efficiency_curve=TableCurve([0.0,25.0,50.0],[0.90,0.95,0.98]))
    second=tg_copy(unit;name=:SecondUnit)
    plant=tg_copy(only(c.system.plants);minup=1.0,mindown=1.0,initial_age=8.0,
        tailwater_curve=TableCurve([0.0,10.0],[0.0,1.0]))
    dry=River(name=:DrySpill,source=:Receiver,target=:Sea,law=:weir,
        coefficient=0.1,crest=100.0,allow_dry=true,capacity=10.0,
        deterministic_delay=0.0,curves=RiverRouting.DelayCurve[],water_value=0.0)
    system=tg_copy(c.system;reservoirs=[source,last(c.system.reservoirs)],
        plants=[plant],generators=[unit,second],rivers=[only(c.system.rivers),dry])
    c=tg_copy(c;system)
    seed=dispatch_from_controls(c,reshape([1,0],2,1),reshape([2.0,0.0],2,1),
        reshape([0.5,1.0],2,1))
    @test seed["validation"]["valid"]
    for fixed in (false,true)
        check_reference_lift(c,seed;joint=true,fixed_u=fixed ? seed["u"] : nothing)
    end
end
