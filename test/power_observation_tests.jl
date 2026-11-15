using ForwardDiff

function drawdown_power_case(;table=nothing,pmax=0.05,initial=3.5)
    reservoir=Reservoir(name=:Lake,z0=0.0,slope=100.0,v0=initial/100,
        vmin=0.0,vmax=0.1,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Plant,qmin=0.01,qmax=10.0,pmin=0.0,pmax=pmax,
        efficiency=0.9,qcurvature=0.0,hcurvature=0.0,hbest=2.0,qbest=1.0,
        hmin=table===nothing ? 0.1 : 1.0,hmax=table===nothing ? 5.0 : 3.5,turbine_table=table,startup=0.0,minup=0.0,mindown=0.0)
    system=HydroSystem(reservoirs=[reservoir],boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Plant,source=:Lake,target=:Sea,pmax=1.0)],generators=[unit])
    ScheduleCase(name="power during drawdown",system=system,grid=[0.0,1.0],prices=[100.0])
end

@testset "Prescribed-storage pressure states advance no time" begin
    c=drawdown_power_case()
    system=OpenSHOP._river_replace(c.system;
        junctions=[Junction(name=:Intake,hmin=0.0,hmax=10.0)],
        tunnels=[Tunnel(name=:Penstock,source=:Lake,target=:Intake,resistance=0.1,capacity=10.0)],
        plants=[Plant(name=:Plant,source=:Intake,target=:Sea,pmax=1.0)])
    x=forward_step(system,[0.035],0.0,[2.0],Float64[],Float64[];commitment=[1])
    @test x.converged
    @test x.storage≈[0.035] atol=1e-12
    @test x.tunnel_q≈[2.0] atol=1e-9
    @test x.generator_head≈[3.5-0.1*2^2] atol=1e-9
    @test x.power≈[0.00981*2*(3.5-0.1*2^2)*0.9] atol=1e-10
    @test_throws ArgumentError forward_step(system,[0.035],-1.0,[2.0],Float64[],Float64[])
end

@testset "Power limits include interior efficiency maxima" begin
    # eta(h)=(4-h)/3: power at fixed flow peaks at h=2, although
    # both endpoint powers and power at the mean head 2.25 are lower.
    table=TurbineTable([1.0,3.5],[TableCurve([0.01,10.0],fill(e,2)) for e in (1.0,1/6)],
        fill(0.01,2),fill(10.0,2))
    c=drawdown_power_case(;table,pmax=0.090)
    g=only(c.system.generators)
    q=(3.5-1.0)/(100*0.0036)
    @test power(g,q,1.0)<g.pmax
    @test power(g,q,3.5)<g.pmax
    @test power(g,q,2.25)<g.pmax
    @test power(g,q,2.0)>g.pmax
    envelope=OpenSHOP.unit_head_envelope(g,q,1.0,3.5)
    @test envelope.at_upper≈2.0 atol=1e-12
    @test envelope.upper≈power(g,q,2.0)-g.pmax atol=1e-12
    @test envelope.upper>0
    x=dispatch_from_controls(c,ones(Int,1,1),reshape([q],1,1),zeros(0,1))
    @test x.storage[1,end]≈0.01 atol=1e-10
    @test !validate(c,x)["valid"]
    @test only(validate(c,x)["shortfalls"].unit_power)≈power(g,q,2.0)-g.pmax atol=1e-8
    # Algebraic extrema retain derivatives, including the movement of the
    # stationary head. Compare to finite differences of the physical envelope.
    f(q)=OpenSHOP.unit_head_envelope(g,q,1.0,3.5).upper
    @test ForwardDiff.derivative(f,q)≈(f(q+1e-4)-f(q-1e-4))/2e-4 rtol=1e-7
    @test all(isfinite,ForwardDiff.hessian(z->OpenSHOP.unit_head_envelope(g,z[1],z[2],z[3]).upper,[q,1.0,3.5]))
    result=solve_case(c;u=ones(Int,1,1),warm=x,time_limit=30.0)
    @test result.accepted
    if result.solution!==nothing
        y=result.solution
        @test y.generator_q[1,1]≈g.pmax/(0.00981*4/3) rtol=1e-3
        @test validate_replay(c,y;grid=OpenSHOP.refined_grid(y.grid;factor=16))["physically_valid"]
    end
end

@testset "Initial power, rather than average power, limits drawdown" begin
    c=drawdown_power_case(;pmax=0.05)
    g=only(c.system.generators)
    x=solve_case(c;u=ones(Int,1,1),time_limit=30.0)
    @test x.accepted
    if x.solution!==nothing
        q=x.solution.generator_q[1,1]
        @test q≈g.pmax/(0.00981*0.9*3.5) atol=2e-4
        @test x.solution.power[1,1]<g.pmax-0.001
        @test validate_replay(c,x.solution;grid=OpenSHOP.refined_grid(x.solution.grid;factor=16))["valid"]
    end
end

@testset "Pressure evaluation retains a fixed-speed pump's flow law" begin
    lower=Reservoir(name=:Lower,z0=10.0,slope=0.0,v0=1.0,vmin=0.0,vmax=2.0,water_value=0.0)
    pump=Pump(name=:Pump,plant=:Lift,qmin=1.0,qmax=20.0,pmin=0.0,pmax=20.0,
        hmin=90.0,hmax=110.0,qbest=10.0,hbest=100.0,motor_efficiency=0.95,
        characteristic=PumpCharacteristic([90.0,110.0],[10.2,9.8],[0.85,0.85]))
    system=HydroSystem(reservoirs=[lower],boundaries=[Boundary(name=:Upper,head=110.0)],
        junctions=[Junction(name=:Intake,hmin=0.0,hmax=20.0)],
        tunnels=[Tunnel(name=:Tunnel,source=:Lower,target=:Intake,resistance=0.01,capacity=20.0)],
        plants=[Plant(name=:Lift,source=:Intake,target=:Upper,pmax=20.0)],generators=[pump])
    x=forward_step(system,[1.0],0.0,[0.0],Float64[],Float64[];commitment=[1])
    q=20/(1+sqrt(1+0.008))
    @test x.converged
    @test x.storage≈[1.0] atol=1e-12
    @test only(x.generator_q)≈q atol=1e-9
    @test only(x.tunnel_q)≈q atol=1e-9
    @test only(x.generator_head)≈100+0.01q^2 atol=1e-9
end

@testset "Head-span bounds agree with physical conversion laws" begin
    generate=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=8.0,pmin=0.02,pmax=0.12,
        hmin=1.0,hmax=3.0,qbest=5.0,hbest=2.0,efficiency=0.92,
        qcurvature=0.1,hcurvature=0.2)
    electrical=OpenSHOP._river_replace(generate;
        generator_efficiency_curve=TableCurve([0.0,0.06,0.12],[0.92,0.95,0.96]))
    pump=OpenSHOP._river_replace(generate;mode=:pump,motor_efficiency=0.95)
    fixed=OpenSHOP._river_replace(pump;
        characteristic=PumpCharacteristic([1.0,2.0,3.0],[6.0,5.0,3.0],[0.8,0.9,0.82]))
    for g in (generate,electrical,pump,fixed)
        envelope=OpenSHOP.unit_head_envelope(g,5.0,1.0,3.0)
        upper=Float64[];lower=Float64[]
        for h in range(1.0,3.0;length=10001)
            q=g.characteristic===nothing ? 5.0 : OpenSHOP.pump_discharge(g.characteristic,h)
            eta=OpenSHOP.efficiency(g,q,h)
            shaft=0.00981q*h
            if g.mode==:generate
                push!(upper,shaft*eta-g.pmax/OpenSHOP.generator_efficiency(g,g.pmax))
                push!(lower,g.pmin/OpenSHOP.generator_efficiency(g,g.pmin)-shaft*eta)
            else
                push!(upper,shaft-g.pmax*g.motor_efficiency*eta)
                push!(lower,g.pmin*g.motor_efficiency*eta-shaft)
            end
        end
        @test envelope.upper≈maximum(upper) atol=1e-8
        @test envelope.lower≈maximum(lower) atol=1e-8
        f(z)=OpenSHOP.unit_head_envelope(g,z[1],z[2],z[3]).upper
        @test all(isfinite,ForwardDiff.gradient(f,[5.0,1.0,3.0]))
        @test all(isfinite,ForwardDiff.hessian(f,[5.0,1.0,3.0]))
        at_two=OpenSHOP.unit_head_envelope(g,5.0,2.0,2.0)
        eta=OpenSHOP.efficiency(g,5.0,2.0)
        expected=g.mode==:generate ? 0.00981*5*2*eta-g.pmax/OpenSHOP.generator_efficiency(g,g.pmax) :
            0.00981*5*2-g.pmax*g.motor_efficiency*eta
        @test at_two.upper≈expected atol=1e-10
    end
end
