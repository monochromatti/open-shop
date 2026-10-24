function window_rule_case(;grid=[0.0,0.5,2.0,3.0,4.0],rules=FlowRequirement[],
    level_windows=LevelWindow[],history_grid=Float64[],history_storage=Float64[],
    level_curve=nothing,inflow=0.0)
    lake=Reservoir(;name=:Lake,z0=90.0,slope=1.0,v0=10.0,vmin=8.0,vmax=12.0,
        water_value=0.0,level_windows,history_grid,history_storage,level_curve,inflow)
    unit=Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=8.0,pmin=0.0,pmax=20.0,
        hmin=60.0,hmax=160.0,hbest=100.0,qbest=4.0,efficiency=0.9,
        qcurvature=0.0,hcurvature=0.0,initial_on=1,minup=0.0,mindown=0.0,startup=0.0)
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=20.0,ramp=nothing)],generators=[unit])
    ScheduleCase(;name="window restrictions",system,grid,prices=fill(50.0,length(grid)-1),
        flow_requirements=rules)
end

window_dispatch(c,q)=dispatch_from_controls(c,ones(Int,1,length(c.prices)),
    reshape(q,1,:),zeros(0,length(c.prices)))

@testset "Window water volumes use physical overlap and survive restart" begin
    rule=FlowRequirement(name=:Gauge,generators=[:Unit],
        history_grid=[-2.0,-0.5,0.0],history_flow=[1.0,3.0],
        volume_windows=[VolumeWindow(start_time=0.25,end_time=2.5,
            min_volume=0.0288,max_volume=0.0288),
            VolumeWindow(start_time=-1.0,end_time=0.25,min_volume=0.009,max_volume=0.009)])
    c=window_rule_case(rules=[rule])
    x=window_dispatch(c,[2.0,4.0,3.0,1.0])
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    continued=restart_case(c,x,0.5)
    y=window_dispatch(continued,[4.0,3.0,1.0])
    @test validate(continued,y)["valid"]
    @test replay_audit(continued,y)["valid"]
    @test y.storage≈x.storage[:,2:end] atol=1e-10
    again=restart_case(continued,y,2.0)
    z=window_dispatch(again,[3.0,1.0])
    @test validate(again,z)["valid"]
    @test z.storage≈x.storage[:,3:end] atol=1e-10
    bad=OpenSHOP._river_replace(c;flow_requirements=[OpenSHOP._river_replace(rule;
        volume_windows=[VolumeWindow(start_time=0.25,end_time=2.5,max_volume=0.028)])])
    @test validate(bad,x)["residuals"]["volume_window_Gauge_1"]≈0.0008 atol=1e-10
end

@testset "Rolling limits catch shifted historical knots and both sides of jumps" begin
    rule=FlowRequirement(name=:Gauge,generators=[:Unit],history_grid=[-2.0,-1.25,-0.75,0.0],
        history_flow=[0.0,8.0,0.0],rolling_ramps=[RollingRamp(hours=1.0,up=1.0)])
    c=window_rule_case(grid=[0.0,2.0],rules=[rule])
    x=window_dispatch(c,[2.0])
    audit=validate(c,x)
    @test !audit["valid"]
    # Both scheduling endpoints have zero upward deviation. At t=0.25 the
    # historical pulse has just left the window: mean=0.5, current flow=2.
    @test audit["residuals"]["rolling_ramp_Gauge_1"]≈0.5 atol=1e-10
    jumps=OpenSHOP._river_replace(rule;history_grid=[-1.0,0.0],history_flow=[0.0])
    cc=window_rule_case(grid=[0.0,1.0,2.0],rules=[jumps])
    y=window_dispatch(cc,[0.0,4.0])
    @test validate(cc,y)["residuals"]["rolling_ramp_Gauge_1"]≈3.0
    down=OpenSHOP._river_replace(jumps;history_flow=[4.0],
        rolling_ramps=[RollingRamp(hours=1.0,down=1.0)])
    cdown=OpenSHOP._river_replace(cc;flow_requirements=[down])
    @test validate(cdown,y)["residuals"]["rolling_ramp_Gauge_1"]≈3.0
end

@testset "Rolling and fixed-period flow averages retain history through restarts" begin
    rule=FlowRequirement(name=:Gauge,generators=[:Unit],history_grid=[-2.0,0.0],history_flow=[2.0],
        rolling_ramps=[RollingRamp(hours=1.25,up=2.0,down=3.0)],
        period_ramps=[PeriodRamp(grid=[-1.0,0.0,2.0,4.0],up=1.5,down=1.5)])
    c=window_rule_case(rules=[rule])
    x=window_dispatch(c,[2.0,4.0,3.0,1.0])
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    strict=OpenSHOP._river_replace(c;flow_requirements=[OpenSHOP._river_replace(rule;
        period_ramps=[PeriodRamp(grid=[-1.0,0.0,2.0,4.0],up=1.0,down=1.0)])])
    # Period means are 2, 3.5 and 2; arithmetic averages would miss this.
    @test validate(strict,x)["residuals"]["period_ramp_Gauge_1"]≈0.5
    continued=restart_case(c,x,0.5)
    y=window_dispatch(continued,[4.0,3.0,1.0])
    @test validate(continued,y)["valid"]
    @test replay_audit(continued,y)["valid"]
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
end

@testset "Reservoir level ranges include history, internal extrema, and nonlinear levels" begin
    windows=[LevelWindow(start_time=-0.5,end_time=2.0,max_change=0.0302)]
    c=window_rule_case(level_windows=windows,history_grid=[-1.0,0.0],history_storage=[10.01,10.0])
    x=window_dispatch(c,[2.0,4.0,3.0,1.0])
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    continued=restart_case(c,x,0.5)
    y=window_dispatch(continued,[4.0,3.0,1.0])
    @test validate(continued,y)["valid"]
    @test replay_audit(continued,y)["valid"]
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
    strict=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        reservoirs=[OpenSHOP._river_replace(only(c.system.reservoirs);
            level_windows=[LevelWindow(start_time=-0.5,end_time=2.0,max_change=0.03)])]))
    @test validate(strict,x)["residuals"]["level_window_Lake_1"]≈0.0002 atol=1e-10
    table=TableCurve([8.0,10.0,12.0],[80.0,100.0,140.0])
    curved=window_rule_case(level_curve=table,
        level_windows=[LevelWindow(start_time=0.0,end_time=2.0,max_change=0.252)])
    z=window_dispatch(curved,[2.0,4.0,3.0,1.0])
    @test validate(curved,z)["valid"]
    # An interior rise and fall must not disappear when endpoints coincide.
    reversal=window_rule_case(grid=[0.0,1.0,2.0],inflow=4.0,
        level_windows=[LevelWindow(start_time=0.0,end_time=2.0,max_change=0.01)])
    v=window_dispatch(reversal,[0.0,8.0])
    @test v.storage[1,1]≈v.storage[1,end] atol=1e-10
    @test validate(reversal,v)["residuals"]["level_window_Lake_1"]≈0.0044 atol=1e-10
end

@testset "Window restrictions change optimized schedules in local and global solves" begin
    rule=FlowRequirement(name=:Gauge,generators=[:Unit],history_grid=[-1.0,0.0],history_flow=[1.0],
        volume_windows=[VolumeWindow(start_time=0.25,end_time=3.5,max_volume=0.036)],
        rolling_ramps=[RollingRamp(hours=1.0,up=2.0,down=2.0)],
        period_ramps=[PeriodRamp(grid=[-1.0,0.0,2.0,4.0],up=2.0,down=2.0)])
    c=window_rule_case(rules=[rule],
        level_windows=[LevelWindow(start_time=0.0,end_time=4.0,max_change=0.05)])
    warm=window_dispatch(c,ones(4))
    @test validate(c,warm)["valid"]
    local_result=solve_case(c;u=ones(Int,1,4),warm,time_limit=30.0)
    @test local_result.accepted
    if local_result.solution!==nothing
        x=local_result.solution
        @test validate(c,x)["valid"]
        @test x.objective>warm.objective+1.0
        @test replay_audit(c,x)["valid"]
        @test sum(0.0036*max(0.0,min(3.5,c.grid[t+1])-max(0.25,c.grid[t]))*x.generator_q[1,t]
            for t in eachindex(c.prices))<=0.036+1e-6
        global_result=OpenSHOP.solve(c;initial=x,time_limit=30.0,preparation_time_limit=0.0)
        @test global_result.accepted
        @test global_result.global_bound===nothing || global_result.global_bound>=x.objective-1e-3
        @test get(global_result.diagnostics,"start_audit",Dict("valid"=>false))["valid"]
    end
    # A nonlinear level table must also constrain actual optimized dispatch.
    table=TableCurve([8.0,10.0,12.0],[80.0,100.0,140.0])
    curved=window_rule_case(level_curve=table,
        level_windows=[LevelWindow(start_time=0.0,end_time=4.0,max_change=0.12)])
    idle=window_dispatch(curved,zeros(4))
    result=solve_case(curved;u=ones(Int,1,4),warm=idle,time_limit=30.0)
    @test result.accepted
    if result.solution!==nothing
        @test validate(curved,result.solution)["valid"]
        @test result.solution.storage[1,end]≈9.988 atol=2e-5
        global_result=OpenSHOP.solve(curved;initial=result.solution,time_limit=20.0,preparation_time_limit=0.0)
        @test global_result.accepted
        @test validate(curved,global_result.solution)["valid"]
    end
    constant=window_rule_case(grid=[0.0,1.0,2.0],inflow=1.0,
        level_windows=[LevelWindow(start_time=0.0,end_time=2.0,max_change=0.0)])
    reference=window_dispatch(constant,ones(2))
    result=solve_case(constant;u=ones(Int,1,2),warm=reference,time_limit=20.0)
    @test result.accepted
    @test result.solution.storage≈fill(10.0,1,3) atol=1e-8
    zero_rule=FlowRequirement(name=:Gauge,generators=[:Unit],
        history_grid=[-1.0,0.0],history_flow=[0.0],
        volume_windows=[VolumeWindow(start_time=0.0,end_time=4.0,min_volume=0.0,max_volume=0.0)],
        rolling_ramps=[RollingRamp(hours=1.0,up=0.0,down=0.0)],
        period_ramps=[PeriodRamp(grid=[-1.0,0.0,2.0,4.0],up=0.0,down=0.0)])
    stopped=window_rule_case(rules=[zero_rule])
    off=dispatch_from_controls(stopped,zeros(Int,1,4),zeros(1,4),zeros(0,4))
    result=solve_case(stopped;u=zeros(Int,1,4),warm=off,time_limit=20.0)
    @test result.accepted
    @test result.solution.generator_q==zeros(1,4)
end

@testset "Incomplete or ambiguous window input is rejected" begin
    for rule in (
        FlowRequirement(name=:Gauge,rolling_ramps=[RollingRamp(hours=1.0,up=1.0)]),
        FlowRequirement(name=:Gauge,rolling_ramps=[RollingRamp(hours=-1.0,up=1.0)]),
        FlowRequirement(name=:Gauge,volume_windows=[VolumeWindow(start_time=0.0,end_time=5.0,max_volume=1.0)]),
        FlowRequirement(name=:Gauge,volume_windows=[VolumeWindow(start_time=1.0,end_time=1.0,max_volume=1.0)]),
        FlowRequirement(name=:Gauge,period_ramps=[PeriodRamp(grid=[0.0,2.0,5.0],up=1.0)]),
    )
        @test_throws ArgumentError OpenSHOP.validate_inputs(window_rule_case(rules=[rule]))
    end
    @test_throws ArgumentError OpenSHOP.validate_inputs(window_rule_case(
        level_windows=[LevelWindow(start_time=0.25,end_time=2.0,max_change=0.1)]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(window_rule_case(
        level_windows=[LevelWindow(start_time=-1.0,end_time=2.0,max_change=0.1)]))
end

@testset "Window observations retain delayed arrivals and signed tunnel orientation" begin
    lakes=[Reservoir(name=n,z0=z,slope=0.01,v0=10.0,vmin=8.0,vmax=12.0,water_value=0.0)
        for (n,z) in ((:Upper,100.0),(:Lower,0.0))]
    river=River(name=:Reach,source=:Upper,target=:Lower,law=:controlled,capacity=4.0,
        curves=RiverRouting.DelayCurve[],deterministic_delay=0.5,water_value=0.0,
        history_grid=[-1.0,0.0],history_release=[0.0])
    c=ScheduleCase(name="arrival windows",system=HydroSystem(reservoirs=lakes,rivers=[river]),
        grid=[0.0,0.5,2.0,3.0,4.0],prices=zeros(4),flow_requirements=[FlowRequirement(
            name=:Gauge,river_arrivals=[:Reach],history_grid=[-1.0,0.0],history_flow=[0.0],
            volume_windows=[VolumeWindow(start_time=0.25,end_time=2.5,
                min_volume=0.0189,max_volume=0.0189)],
            rolling_ramps=[RollingRamp(hours=1.0,up=3.0,down=3.0)])])
    gate=reshape([1.0,0.5,0.75,0.0],1,:)
    x=dispatch_from_controls(c,zeros(Int,0,4),zeros(0,4),gate)
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    continued=restart_case(c,x,0.5)
    y=dispatch_from_controls(continued,zeros(Int,0,3),zeros(0,3),gate[:,2:end])
    @test validate(continued,y)["valid"]
    @test y.arrival_volume≈x.arrival_volume[:,2:end] atol=1e-10
    local_result=solve_case(c;u=zeros(Int,0,4),warm=x,time_limit=20.0)
    @test local_result.accepted
    # Fixed boundary heads force negative flow through the declared orientation.
    s=HydroSystem(boundaries=[Boundary(name=:A,head=0.0),Boundary(name=:B,head=4.0)],
        tunnels=[Tunnel(name=:Link,source=:A,target=:B,resistance=1.0,capacity=3.0)])
    signed=ScheduleCase(name="signed windows",system=s,grid=[0.0,1.0,2.0],prices=zeros(2),
        flow_requirements=[FlowRequirement(name=:Gauge,tunnels=[:Link],min_flow=-3.0,
            history_grid=[-1.0,0.0],history_flow=[-2.0],
            volume_windows=[VolumeWindow(start_time=-0.5,end_time=1.5,
                min_volume=-0.0144,max_volume=-0.0144)],
            rolling_ramps=[RollingRamp(hours=1.0,up=0.1,down=0.1)])])
    z=dispatch_from_controls(signed,zeros(Int,0,2),zeros(0,2),zeros(0,2))
    @test z.tunnel_q≈fill(-2.0,1,2) atol=1e-10
    @test validate(signed,z)["valid"]
    result=OpenSHOP.solve(signed;initial=z,time_limit=15.0,preparation_time_limit=0.0)
    @test result.accepted
    @test validate(signed,result.solution)["valid"]
    @test get(result.diagnostics,"start_audit",Dict("valid"=>false))["valid"]
end
