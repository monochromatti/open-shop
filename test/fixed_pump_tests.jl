function fixed_pump_case(;network=false,characteristic=PumpCharacteristic(
    [90.0,100.0,110.0],[10.2,10.0,9.8],[0.845,0.85,0.855]))
    lower=Reservoir(name=:Lower,z0=-10.0,v0=10.0,vmin=9.0,vmax=11.0,water_value=0.0)
    upper=Reservoir(name=:Upper,z0=90.0,v0=10.0,vmin=9.0,vmax=11.0,water_value=0.0)
    pump=Pump(name=:Pump,plant=:Lift,qmin=9.0,qmax=11.0,pmin=1.0,pmax=20.0,
        characteristic=characteristic,qbest=10.0,hbest=100.0,hmin=90.0,hmax=110.0,
        motor_efficiency=0.95,min_efficiency=0.8,minup=0.5,mindown=0.75,
        startup=0.0,shutdown=0.0)
    plant=Plant(name=:Lift,source=network ? :Intake : :Lower,target=:Upper,pmax=20.0,
        ramp=nothing,tailwater_curve=TableCurve([0.0,20.0],[0.0,1.0]))
    system=HydroSystem(reservoirs=[lower,upper],plants=[plant],generators=[pump],
        junctions=network ? [Junction(name=:Intake,hmin=-20.0,hmax=20.0)] : Junction[],
        tunnels=network ? [Tunnel(name=:IntakeTunnel,source=:Lower,target=:Intake,
            resistance=0.01,capacity=12.0,gate_positions=[0.0,1.0])] : Tunnel[])
    ScheduleCase(name="fixed-speed pumping",system=system,grid=[0.0,0.5,1.25,2.0],prices=[-50.0,20.0,-50.0])
end

@testset "Fixed-speed pump characteristic follows coupled head" begin
    c=fixed_pump_case();u=[1 0 1];guess=zeros(1,3)
    x=dispatch_from_controls(c,u,guess,zeros(0,3))
    @test validate(c,x)["valid"]
    @test guess==zeros(1,3)
    step=forward_step(c.system,[10.0,10.0],0.5,[0.0],Float64[],Float64[];commitment=[1])
    @test step.converged
    @test step.generator_q≈x.generator_q[:,1] atol=1e-9
    @test x.generator_q[1,2]==x.power[1,2]==0.0
    volumes=[10.0,10.0]
    for t in 1:3
        dt=c.grid[t+1]-c.grid[t]
        h0=100.0+volumes[2]-volumes[1]
        k=0.05+0.0036dt
        q=u[t]==0 ? 0.0 : (12.0-0.02h0)/(1+0.02k)
        h=h0+k*q
        eta=0.8+0.0005h
        consumed=0.00981q*h/(eta*0.95)
        @test x.generator_q[1,t]≈q atol=1e-9
        @test x.power[1,t]≈consumed atol=1e-9
        volumes .+=0.0036dt*q .* [-1,1]
        @test x.storage[:,t+1]≈volumes atol=1e-10
    end
    @test vec(sum(x.storage;dims=1))≈fill(20.0,4) atol=1e-10
    @test net_power(c,x)≈-vec(x.power) atol=1e-12
    @test x.objective≈-sum(c.prices.*diff(c.grid).*vec(x.power)) atol=1e-9
    alternate=dispatch_from_controls(c,u,fill(10.5,1,3),zeros(0,3))
    @test alternate.generator_q≈x.generator_q atol=1e-9
    @test simulate(c,zeros(1,3),zeros(0,3);commitment=u).generator_q≈x.generator_q atol=1e-9
    @test simulate(c,zeros(Int,1,3),zeros(0,3);commitment=u).generator_q≈x.generator_q atol=1e-9
    @test simulate(c,[10.0 0.0 10.0],zeros(0,3)).generator_q≈x.generator_q atol=1e-9
    @test replay_audit(c,x)["valid"]
    refined=simulate(c,guess,zeros(0,3);commitment=u,grid=OpenSHOP.refined_grid(c.grid;factor=4))
    @test refined.converged
    @test refined.generator_q[1,1]!=refined.generator_q[1,4]
    @test maximum(abs,refined.storage[:,[1,5,9,13]]-x.storage)<1e-5
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    mktemp() do path,io
        close(io)
        writejson(path,case_dict(c))
        @test case_dict(readcase(path))==case_dict(c)
    end
    bad=schedule_dict(x);bad["generator_q"][1,1]+=0.2
    @test !validate(c,bad)["valid"]
    @test validate(c,bad)["residuals"]["pump_characteristic_Pump"]>0.19
    continued=restart_case(c,x,0.5)
    remainder=dispatch_from_controls(continued,u[:,2:3],zeros(1,2),zeros(0,2))
    @test remainder.generator_q≈x.generator_q[:,2:3] atol=1e-9
    @test remainder.storage≈x.storage[:,2:4] atol=1e-9
    @test continued.system.generators[1].characteristic!==nothing
    @test continued.system.generators[1].initial_discharge≈x.generator_q[1,1]
    @test continued.system.generators[1].initial_on==1
    @test continued.system.generators[1].initial_age==0.5

    local_result=solve_case(c;u,warm=x,time_limit=45.0)
    @test local_result.accepted
    if local_result.solution!==nothing
        @test local_result.solution.generator_q≈x.generator_q atol=1e-6
        @test local_result.solution.objective≈x.objective atol=1e-4
    end
    result=solve(c;initial=x,time_limit=90.0,relative_gap=1e-4)
    @test result.accepted
    @test result.global_certificate
    if result.solution!==nothing
        @test result.solution.commitment==u
        @test validate(c,result.solution)["valid"]
        @test replay_audit(c,result.solution)["valid"]
    end
end

@testset "Fixed-speed pumping through a selectable pressure intake" begin
    c=fixed_pump_case(;network=true);u=ones(Int,1,3)
    x=dispatch_from_controls(c,u,zeros(1,3),zeros(0,3);tunnel_opening=ones(1,3))
    @test validate(c,x)["valid"]
    @test x.tunnel_q≈x.generator_q atol=1e-9
    for t in 1:3
        dt=diff(c.grid)[t]
        h0=100.0+x.storage[2,t]-x.storage[1,t]
        b=1+0.02*(0.05+0.0036dt)
        q=2*(12-0.02h0)/(b+sqrt(b^2+4*0.0002*(12-0.02h0)))
        @test x.generator_q[1,t]≈q atol=1e-8
    end
    @test replay_audit(c,x)["valid"]
    forced=OpenSHOP._river_replace(c;operations=[OperationalSeries(object=:Pump,
        attribute=:forced_on,times=[0.0],values=[1.0])])
    proposal=OpenSHOP.propose_commitment(forced;reference=x,time_limit=20.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        @test proposal["u"]==u
        @test proposal["tunnel_opening"]==ones(1,3)
    end
    result=solve_case(c;u,warm=x,tunnel_opening=ones(1,3),time_limit=45.0)
    @test result.accepted
end

@testset "Fixed-speed pump input policies and bounds" begin
    @test_throws ArgumentError PumpCharacteristic([100.0],[10.0],[0.85])
    @test_throws ArgumentError PumpCharacteristic([100.0,90.0],[10.0,11.0],[0.8,0.85])
    @test_throws ArgumentError PumpCharacteristic([90.0,100.0],[0.0,10.0],[0.8,0.85])
    @test_throws ArgumentError PumpCharacteristic([90.0,100.0],[11.0,10.0],[0.8,1.01])
    c=fixed_pump_case()
    @test pump_discharge(c.system.generators[1].characteristic,95.0)≈10.1
    @test_throws DomainError pump_discharge(c.system.generators[1].characteristic,80.0)
    single=PumpCharacteristic([100.0],[10.0],[0.85];head_extrapolation=:constant)
    cc=fixed_pump_case(;characteristic=single)
    x=dispatch_from_controls(cc,[1 0 1],zeros(1,3),zeros(0,3))
    @test x.generator_q≈[10.0 0.0 10.0] atol=1e-9
    @test validate(cc,x)["valid"]
    @test_throws DomainError pump_discharge(single,NaN)
    linear=PumpCharacteristic([98.0,102.0],[10.04,9.96],[0.849,0.851];head_extrapolation=:linear)
    extended=fixed_pump_case(;characteristic=linear)
    xx=dispatch_from_controls(extended,[1 0 1],zeros(1,3),zeros(0,3))
    @test validate(extended,xx)["valid"]
    @test xx.generator_q[1,1]≈10/(1+0.02*(0.05+0.0036*0.5)) atol=1e-9
    profile=case_dict(c)
    profile["generators"][1]["mode"]=:generate
    @test_throws ArgumentError case_from_dict(profile)
    profile=case_dict(c)
    profile["generators"][1]["turbine_table"]=Dict("heads"=>[90.0,110.0],
        "discharge"=>[0.0,11.0],"efficiency"=>[0.85 0.85;0.85 0.85],
        "qmin"=>[9.0,9.0],"qmax"=>[11.0,11.0])
    @test_throws ArgumentError case_from_dict(profile)
    profile=case_dict(c);profile["generators"][1]["hmax"]=120.0
    @test_throws ArgumentError case_from_dict(profile)
    profile=case_dict(c)
    profile["generators"][1]["characteristic"]=Dict("heads"=>[100.0,101.0],
        "discharge"=>[1.0,0.5],"efficiency"=>[0.8,0.9],"head_extrapolation"=>"linear")
    @test_throws ArgumentError case_from_dict(profile)
    limited=OpenSHOP._river_replace(c;operations=[OperationalSeries(object=:Pump,
        attribute=:qmax,times=[0.0],values=[9.5])])
    x=dispatch_from_controls(limited,[1 0 1],zeros(1,3),zeros(0,3))
    @test !validate(limited,x)["valid"]
    @test x.generator_q[1,1]>9.5 # no clipping to make an impossible point feasible
    profile=case_dict(c);profile["generators"][1]["initial_age"]=0.25
    @test !admissible(case_from_dict(profile),[1 0 1])

    turbine=Generator(name=:Turbine,plant=:Return,qmin=9.0,qmax=11.0,pmin=1.0,pmax=20.0,
        hmin=90.0,hmax=110.0,initial_on=0,minup=0.0,mindown=0.0)
    system=OpenSHOP._river_replace(c.system;generators=[c.system.generators[1],turbine],
        plants=[c.system.plants[1],Plant(name=:Return,source=:Upper,target=:Lower,pmax=20.0,ramp=nothing)],
        commitment_exclusions=[CommitmentExclusion(name=:Reversible,groups=[[:Pump],[:Turbine]])])
    reversible=OpenSHOP._river_replace(c;system)
    @test !admissible(reversible,ones(Int,2,3))
    @test_throws ArgumentError simulate(reversible,zeros(2,3),zeros(0,3);commitment=ones(Int,2,3))
end
