@testset "Curved storage levels have a quantified midpoint energy error" begin
    curvature=0.2
    lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,curvature=curvature,
        v0=1.0,vmin=0.5,vmax=3.5,water_value=0.0)
    generators=[Generator(name=name,plant=plant,mode=mode,qmin=0.0,qmax=600.0,
        pmin=0.0,pmax=1000.0,hmin=90.0,hmax=120.0,hbest=100.0,qbest=300.0,
        efficiency=1.0,qcurvature=0.0,hcurvature=0.0,initial_on=0,
        minup=0.0,mindown=0.0,startup=0.0) for (name,plant,mode) in
        ((:Turbine,:Generation,:generate),(:Pump,:Pumping,:pump))]
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Generation,source=:Lake,target=:Sea,pmax=1000.0,ramp=nothing),
            Plant(name=:Pumping,source=:Sea,target=:Lake,pmax=1000.0,ramp=nothing)],generators=generators)
    c=ScheduleCase(name="ideal closed water cycle",system=system,grid=[0.0,1.0,2.0,3.0],
        prices=ones(3))
    u=[0 1 1;1 0 0]
    q=[0.0 1.0 1.0;2.0 0.0 0.0]/0.0036
    coarse=dispatch_from_controls(c,u,q,zeros(0,3))
    @test validate(c,coarse)["valid"]
    @test coarse.storage≈[1.0 3.0 2.0 1.0] atol=1e-10
    # Integral of level dV returns to zero. Midpoint work differs by
    # curvature*ΔV³/12 on each step, yielding this artificial electrical gain.
    gain=(0.00981/0.0036)*curvature/2
    @test sum(net_power(c,coarse))≈gain atol=1e-9
    for factor in (2,4)
        fine=simulate(c,q,zeros(0,3);commitment=u,grid=OpenSHOP.refined_grid(c.grid;factor))
        @test fine.converged
        @test fine.storage[1,end]≈lake.v0 atol=1e-10
        @test sum(net_power(c,fine).*diff(fine.grid))≈gain/factor^2 atol=1e-9
    end
end

@testset "Passivity alone does not imply a unique midpoint step" begin
    conductances=[3.0,1.0,1/3]/0.0036
    response=LevelFlowCurve([2.25,2.5,2.75],
        [TableCurve([-2.0,0.0,2.0],[-2*k,0.0,2*k]) for k in conductances])
    system=HydroSystem(reservoirs=[
        Reservoir(name=:Upper,z0=0.0,slope=1.0,v0=3.0,vmin=0.1,vmax=4.0,water_value=0.0),
        Reservoir(name=:Lower,z0=0.0,slope=1.0,v0=1.0,vmin=0.1,vmax=4.0,water_value=0.0)],
        rivers=[River(name=:Channel,source=:Upper,target=:Lower,law=:delta,
            curves=RiverRouting.DelayCurve[],water_value=0.0,capacity=2000.0,
            deterministic_delay=0.0,delta_level_ref_up_flow_curve=response)])
    c=ScheduleCase(name="multiple midpoint roots",system=system,grid=[0.0,1.0],prices=[0.0])
    @test OpenSHOP.validate_inputs(c)
    for volume in (1.5,2.0,2.5)
        flow=(3-volume)/0.0036
        step=forward_step(system,[3.0,1.0],1.0,Float64[],[0.0],[0.0];
            initial_guess=[volume,4-volume,flow])
        @test step.converged
        @test step.storage≈[volume,4-volume] atol=1e-10
        @test step.river_release≈[flow] atol=1e-9
        @test (step.head[1]-step.head[2])*flow>0
    end
end

@testset "Replay retains a feasible pressure reference for an idle unit" begin
    system=HydroSystem(junctions=[Junction(name=:Intake,hmin=90.0,hmax=110.0)],
        boundaries=[Boundary(name=:Tail,head=100.0)],
        plants=[Plant(name=:Station,source=:Intake,target=:Tail,pmax=1.0,ramp=nothing)],
        generators=[Generator(name=:Unit,plant=:Station,pmin=0.0,qmin=0.0,pmax=1.0,qmax=1.0,
            hmin=5.0,hmax=10.0,hbest=5.0,qbest=1.0,qcurvature=0.0,hcurvature=0.0,
            initial_on=1,minup=0.0,mindown=0.0,startup=0.0)])
    c=ScheduleCase(name="idle pressure reference",system=system,grid=[0.0,1.0],prices=[0.0],
        operations=[OperationalSeries(object=:Unit,attribute=:forced_on,times=[0.0],values=[1.0])])
    reference=reshape([106.0,100.0],2,1)
    warm=dispatch_from_controls(c;commitment=ones(Int,1,1),generator_q=zeros(1,1),
        gate=zeros(0,1),head_reference=reference)
    @test validate(c,warm)["valid"]
    @test replay_audit(c,warm)["valid"]
    @test reference==reshape([106.0,100.0],2,1)
    for result in (solve_case(c;u=ones(Int,1,1),warm,time_limit=30.0),
        solve(c;initial=warm,time_limit=45.0))
        @test result.accepted
        result.solution===nothing || @test replay_audit(c,result.solution)["valid"]
    end
    invalid=dispatch_from_controls(c,ones(Int,1,1),zeros(1,1),zeros(0,1);
        head_reference=reshape([100.0,100.0],2,1))
    @test !validate(c,invalid)["valid"]
    @test_throws DimensionMismatch simulate(c,zeros(1,1),zeros(0,1);head_reference=zeros(1,1))
    @test_throws ArgumentError simulate(c,zeros(1,1),zeros(0,1);head_reference=fill(NaN,2,1))
end

@testset "A head reference cannot replace an anchored tunnel equation" begin
    system=HydroSystem(junctions=[Junction(name=:Intake,hmin=90.0,hmax=110.0)],
        boundaries=[Boundary(name=:Supply,head=100.0)],
        tunnels=[Tunnel(name=:Connection,source=:Supply,target=:Intake,resistance=1.0,capacity=1.0)])
    reference=[1000.0,100.0]
    step=forward_step(system,Float64[],1.0,Float64[],Float64[],Float64[];
        head_reference=reference)
    @test step.converged
    @test step.head≈[100.0,100.0] atol=1e-8
    @test step.tunnel_q≈[0.0] atol=1e-8
    @test reference==[1000.0,100.0]
end

@testset "A coarse midpoint step can overshoot passive level equalization" begin
    gain=1000.0
    response=LevelFlowCurve([0.0],[TableCurve([-4.0,0.0,4.0],[-4gain,0.0,4gain])])
    system=HydroSystem(reservoirs=[
        Reservoir(name=:Upper,z0=0.0,slope=1.0,v0=3.0,vmin=0.1,vmax=4.0,water_value=0.0),
        Reservoir(name=:Lower,z0=0.0,slope=1.0,v0=1.0,vmin=0.1,vmax=4.0,water_value=0.0)],
        rivers=[River(name=:Channel,source=:Upper,target=:Lower,law=:delta,
            curves=RiverRouting.DelayCurve[],water_value=0.0,capacity=4000.0,
            deterministic_delay=0.0,delta_level_ref_up_flow_curve=response)])
    coarse=forward_step(system,[3.0,1.0],1.0,Float64[],[0.0],[0.0])
    @test coarse.converged
    ratio=0.0036gain
    @test coarse.storage[1]-coarse.storage[2]≈2*(1-ratio)/(1+ratio) atol=1e-9
    @test coarse.storage[1]<coarse.storage[2]
    c=ScheduleCase(name="stiff passive equalization",system=system,grid=[0.0,1.0],prices=[0.0])
    fine=simulate(c,zeros(0,1),zeros(1,1);grid=[0.0,0.25,0.5,0.75,1.0])
    @test fine.converged
    @test fine.storage[1,end]-fine.storage[2,end]≈2*((1-ratio/4)/(1+ratio/4))^4 atol=1e-9
    @test all(fine.storage[1,:].>fine.storage[2,:])
end
