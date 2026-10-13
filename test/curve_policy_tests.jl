function single_head_case(; interpolation=:bilinear)
    table=TurbineTable(;heads=[100.0], discharge=[1.0,3.0], efficiency=[0.88,0.92],
        qmin=[0.5], qmax=[4.0], head_extrapolation=:constant,
        discharge_extrapolation=:linear, interpolation)
    lake=Reservoir(;name=:Lake,z0=98.0,slope=1.0,v0=2.0,vmin=1.0,vmax=3.0,
        inflow=2.0,water_value=0.0)
    unit=Generator(;name=:Unit,plant=:Plant,qmin=0.5,qmax=4.0,pmin=0.1,pmax=5.0,
        hmin=90.0,hmax=110.0,qbest=2.0,hbest=100.0,minup=0.0,mindown=0.0,
        startup=0.0,turbine_table=table,
        generator_efficiency_curve=TableCurve([0.0,5.0],[0.95,0.95]))
    system=HydroSystem(;reservoirs=[lake],generators=[unit],
        boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Plant,source=:Lake,target=:Sea,pmax=5.0)])
    ScheduleCase(;name="single_head",system,grid=[0.0,1.0,2.0],prices=[50.0,100.0])
end

@testset "Explicit single-head and discharge continuation" begin
    @test_throws ArgumentError TurbineTable([100.0],[1.0,3.0],[0.88,0.92],
        [1.0],[3.0])
    @test_throws ArgumentError TurbineTable([100.0],[1.0,3.0],[0.88,0.92],
        [1.0],[3.0];head_extrapolation=:linear)
    @test_throws ArgumentError TurbineTable([90.0,110.0],[1.0,3.0],fill(0.9,2,2),
        [1.0,1.0],[3.0,3.0];head_extrapolation=:constant)
    @test_throws ArgumentError TurbineTable([100.0],[1.0,3.0],[0.88,0.92],
        [0.5],[4.0];head_extrapolation=:constant)
    @test_throws ArgumentError TurbineTable([100.0],[1.0,3.0],[0.88,0.92],
        [-0.5],[4.0];head_extrapolation=:constant,discharge_extrapolation=:linear)
    @test_throws ArgumentError TurbineTable([100.0],[1.0,3.0],[0.88,0.92],
        [1.0],[3.0];head_extrapolation=:constant,discharge_extrapolation=:clamp)
    for interpolation in (:bilinear,:pchip_discharge)
        c=single_head_case(;interpolation)
        g=only(c.system.generators);table=g.turbine_table
        @test OpenSHOP.validate_inputs(c)
        @test table.heads==[100.0]
        @test table.discharge==[1.0,3.0]
        @test table.qmin==[0.5] && table.qmax==[4.0]
        @test g.qmin==0.5 && g.qmax==4.0
        for h in (70.0,100.0,130.0),q in (1.0,2.0,3.0)
            @test OpenSHOP.turbine_efficiency(table,q,h)≈0.86+0.02q atol=1e-12
            @test OpenSHOP.turbine_qmin(table,h)==0.5
            @test OpenSHOP.turbine_qmax(table,h)==4.0
            @test power(g,q,h)≈0.00981*q*h*(0.86+0.02q)*0.95 atol=1e-12
        end
        @test_throws DomainError OpenSHOP.turbine_efficiency(table,0.5,100.0)
        @test_throws DomainError OpenSHOP.turbine_efficiency(table,4.0,100.0)
        @test_throws DomainError OpenSHOP.turbine_efficiency(table,2.0,NaN)
        @test OpenSHOP.turbine_efficiency(table,0.5,100.0;extrapolation=:linear)≈0.87
        @test OpenSHOP.turbine_efficiency(table,4.0,100.0;extrapolation=:linear)≈0.94
        @test power(g,0.0,100.0)==0.0
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        missing=case_dict(c)
        delete!(missing["generators"][1]["turbine_table"],"head_extrapolation")
        @test_throws ArgumentError case_from_dict(missing)
        missing=case_dict(c)
        delete!(missing["generators"][1]["turbine_table"],"discharge_extrapolation")
        @test_throws ArgumentError case_from_dict(missing)

        # Independent water and energy calculation, beyond the efficiency knots.
        x=dispatch_from_controls(c,ones(Int,1,2),fill(4.0,1,2),zeros(0,2))
        expected_power=0.00981 .* 4.0 .* [99.9964,99.9892] .* 0.94 .* 0.95
        @test x["validation"]["valid"]
        @test vec(x.storage)≈[2.0,1.9928,1.9856] atol=1e-10
        @test vec(x.power)≈expected_power atol=1e-10
        @test x.objective≈sum(c.prices .* expected_power) atol=1e-8
        @test replay_audit(c,x)["valid"]
        local_result=solve_case(c;u=ones(Int,1,2),warm=x,time_limit=10.0)
        @test local_result["accepted"]
        @test local_result["objective"]≈x.objective atol=1e-3
        model=OpenSHOP._build_global_dispatch(c;joint=true)
        lift=OpenSHOP._lift_start!(model,c,x)
        @test lift["valid"]
        @test lift["objective"]≈x.objective atol=1e-8

        off=dispatch_from_controls(c,reshape([1,0],1,2),reshape([2.0,0.0],1,2),zeros(0,2))
        @test off["validation"]["valid"]
        @test off.power[1,2]==0.0
        @test OpenSHOP._lift_start!(model,c,off)["valid"]
        restarted=restart_case(c,x,1.0)
        rest=dispatch_from_controls(restarted,ones(Int,1,1),fill(4.0,1,1),zeros(0,1))
        @test rest["validation"]["valid"]
        @test rest.storage≈x.storage[:,2:end] atol=1e-10
        @test rest.power≈x.power[:,2:end] atol=1e-10
    end

    # A three-knot PCHIP has curved interior segments, but its declared extension
    # follows end secants (.15 and -.075), rather than endpoint derivatives.
    table=TurbineTable([100.0],[1.0,2.0,4.0],[0.8,0.95,0.8],[0.5],[5.0];
        head_extrapolation=:constant,discharge_extrapolation=:linear,
        interpolation=:pchip_discharge)
    @test OpenSHOP.turbine_efficiency(table,1.5,80.0)≈0.903125 atol=1e-12
    @test OpenSHOP.turbine_efficiency(table,0.0,80.0;extrapolation=:linear)≈0.65
    @test OpenSHOP.turbine_efficiency(table,5.0,120.0;extrapolation=:linear)≈0.725

    c=single_head_case()
    cubic_unit=OpenSHOP._river_replace(only(c.system.generators);turbine_table=table)
    cubic_case=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        generators=[cubic_unit]))
    cubic=dispatch_from_controls(cubic_case,ones(Int,1,2),fill(1.5,1,2),zeros(0,2))
    @test cubic["validation"]["valid"]
    @test vec(cubic.power)≈0.00981 .* 1.5 .* [100.0009,100.0027] .* 0.903125 .* 0.95 atol=1e-10
    cubic_model=OpenSHOP._build_global_dispatch(cubic_case;joint=true)
    @test !isempty(cubic_model.graphs.generator_graphs[1,1].turbine.r)
    @test OpenSHOP._lift_start!(cubic_model,cubic_case,cubic)["valid"]

    x=dispatch_from_controls(c,ones(Int,1,2),fill(4.0,1,2),zeros(0,2))
    result=OpenSHOP.solve(c;initial=x,time_limit=30.0,preparation_time_limit=0.0)
    @test result.accepted
    @test result.solution.objective≈x.objective atol=1e-3
    @test validate(c,result.solution)["valid"]
    @test replay_audit(c,result.solution)["valid"]
end
