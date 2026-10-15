function spinning_case(;soft=nothing,available=true,up_energy=1.0,groups=nothing)
    lake=Reservoir(name=:Lake,z0=100.0,slope=1e-8,v0=10.0,vmin=1.0,vmax=100.0,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=2.0,qmax=10.0,pmin=1.7658,pmax=8.829,
        efficiency=0.9,qcurvature=0.0,hcurvature=0.0,qbest=5.0,hbest=100.0,hmin=50.0,hmax=150.0,
        minup=0.0,mindown=0.0,initial_on=0,initial_age=8.0,startup=5.0,
        ramp_up=20.0,ramp_down=20.0,discharge_ramp_up=30.0,discharge_ramp_down=30.0,
        reserve_available=available,reserve_up_energy=up_energy,reserve_down_energy=0.8)
    station=Plant(name=:Station,source=:Lake,target=:Sea,pmax=8.829,ramp=nothing)
    system=HydroSystem(reservoirs=[lake],generators=[unit],plants=[station],boundaries=[Boundary(name=:Sea,head=0.0)])
    contract=ReserveGroup(name=:FRR,generators=[:Unit],up_requirement=3.0,down_requirement=1.0,
        up_price=100.0,down_price=25.0,up_max=4.0,down_max=4.0,shortfall_cost=soft,
        response_hours=0.25,activation_hours=0.25)
    ScheduleCase(name="spinning capacity",system=system,grid=[0.0,1.0,2.0],prices=[50.0,50.0],
        reserve_groups=groups===nothing ? [contract] : groups)
end

@testset "Spinning reserves: energy, capacity and qualified response" begin
    c=spinning_case()
    local_result=solve_case(c;u=ones(Int,1,2),time_limit=120.0)
    @test local_result.accepted
    @test !local_result.global_certificate
    x=local_result.solution
    if x!==nothing
        @test OpenSHOP.validate_reserves(c,x)["valid"]
        @test vec(x.reserves.up)≈[4.0,4.0] atol=1e-5
        @test vec(x.power)≈[4.829,4.829] atol=1e-5
        @test vec(x.reserves.down)≈fill(4.829-1.7658,2) atol=1e-5
        expected=2*(50*4.829+100*4+25*(4.829-1.7658))-5
        @test x.objective≈expected atol=1e-3
        @test replay_audit(c,x)["valid"]
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        profile=OpenSHOP.jsonready(x)
        @test validate(c,profile)["valid"]
        @test OpenSHOP._schedule_input(profile;case=c).reserves.up==x.reserves.up
        mixed=deepcopy(profile)
        mixed["reserves"]=Dict{String,Any}(mixed["reserves"])
        mixed["reserves"]["up"]=copy(x.reserves.up)
        @test OpenSHOP._schedule_input(mixed;case=c).reserves==x.reserves
        mixed=deepcopy(profile)
        mixed["reserves"]=Dict{String,Any}(mixed["reserves"])
        mixed["reserves"]["down"]=copy(x.reserves.down)
        @test OpenSHOP._schedule_input(mixed;case=c).reserves==x.reserves
        bad=schedule_dict(x)
        # Alter a capacity control without changing its physical witness.
        bad["reserves"]=deepcopy(bad["reserves"])
        bad["reserves"].up[1,1,1]+=1.0
        @test !validate(c,bad)["valid"]
        next=restart_case(c,x,1.0)
        @test next.reserve_groups==c.reserve_groups
        @test next.system.generators[1].reserve_up_energy==1.0
        refined=with_grid(c,[0.0,0.5,1.0,1.5,2.0])
        @test refined.reserve_groups==c.reserve_groups
        b=OpenSHOP._build_global_dispatch(c;joint=true)
        lift=OpenSHOP._lift_start!(b,c,x)
        @test lift["valid"]
        for (i,v) in enumerate(all_variables(b.m));set_name(v,"reserve_reference_$i");end
        @test OpenSHOP._lift_start!(b,c,x)["valid"]
        result=OpenSHOP.solve(c;initial=x,time_limit=30.0,relative_gap=1e-3)
        @test result.accepted
        @test result.solution.reserves isa ReserveAllocation
        @test result.objective>=x.objective-1e-4
    end
end

@testset "Reserve eligibility and independent shortfall accounting" begin
    c=spinning_case(;soft=100.0,available=false)
    result=solve_case(c;u=ones(Int,1,2),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        x=result.solution
        @test maximum(abs,x.reserves.up)<1e-6
        @test maximum(abs,x.reserves.down)<1e-6
        @test x.reserves.up_shortfall≈fill(3.0,1,2) atol=1e-6
        @test x.reserves.down_shortfall≈fill(1.0,1,2) atol=1e-6
        @test OpenSHOP.reserve_economics(c,x.reserves).cost≈800.0 atol=1e-5
        bad=schedule_dict(x);r=deepcopy(x.reserves)
        r.up_shortfall[1,1]=0.0;bad["reserves"]=r
        @test !validate(c,bad)["valid"]
        @test OpenSHOP.reserve_economics(c,r).cost≈800.0 atol=1e-5
        # Reserve shortfall never appears in energy or water balances.
        @test validate(c,x)["residuals"]["global_water"]<1e-6
    end
    first=only(spinning_case().reserve_groups)
    second=OpenSHOP._river_replace(first;name=:Other,up_requirement=0.0,down_requirement=0.0)
    c=spinning_case(;groups=[first,second])
    result=solve_case(c;u=ones(Int,1,2),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test sum(result.solution.reserves.up[:,:,1])<=4.0+1e-5 # shared unit MWh budget
    end
    c=spinning_case(;up_energy=0.5,soft=100.0)
    result=solve_case(c;u=ones(Int,1,2),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test vec(result.solution.reserves.up)≈[2.0,2.0] atol=1e-5
        @test vec(result.solution.reserves.up_shortfall)≈[1.0,1.0] atol=1e-5
    end
end

@testset "Reserve witnesses use shared tailwater and real efficiency laws" begin
    c=spinning_case()
    first=only(c.system.generators)
    table=TurbineTable([50.0,100.0,150.0],[0.0,5.0,10.0],
        [0.85 0.9 0.92;0.85 0.9 0.92;0.85 0.9 0.92],fill(2.0,3),fill(10.0,3);
        interpolation=:pchip_discharge)
    first=OpenSHOP._river_replace(first;turbine_table=table,
        generator_efficiency_curve=TableCurve([0.0,9.0],[0.95,0.97]),reserve_up_energy=nothing,reserve_down_energy=nothing)
    second=OpenSHOP._river_replace(first;name=:Second,reserve_available=false)
    plant=OpenSHOP._river_replace(only(c.system.plants);pmax=17.658,
        tailwater_curve=TableCurve([0.0,20.0],[0.0,4.0]))
    group=OpenSHOP._river_replace(only(c.reserve_groups);activation_hours=0.0)
    system=OpenSHOP._river_replace(c.system;generators=[first,second],plants=[plant])
    c=OpenSHOP._river_replace(c;system,reserve_groups=[group])
    result=solve_case(c;u=ones(Int,2,2),time_limit=90.0)
    @test result.accepted
    if result.solution!==nothing
        x=result.solution
        @test OpenSHOP.validate_reserves(c,x)["valid"]
        @test maximum(abs,x.reserves.up[2,:,:])<1e-6
        @test any(x.reserves.up_q[2,:] .> x.generator_q[2,:] .+ 1e-5)
        for t in 1:2
            h=x.head[1,t]-OpenSHOP.tailwater(plant,sum(x.reserves.up_q[:,t]))
            for i in 1:2
                actual=power(c.system.generators[i],x.reserves.up_q[i,t],h)
                @test actual≈x.power[i,t]+sum(x.reserves.up[i,:,t]) atol=1e-4
            end
        end
        @test replay_audit(c,x)["valid"]
        b=OpenSHOP._build_global_dispatch(c;joint=true)
        @test OpenSHOP._lift_start!(b,c,x)["valid"]
    end
end

@testset "Reserve response bounds and nonmonotonic power" begin
    c=spinning_case(;soft=100.0)
    g=OpenSHOP._river_replace(only(c.system.generators);ramp_up=4.0)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]))
    result=solve_case(c;u=ones(Int,1,2),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test vec(result.solution.reserves.up)≈[1.0,1.0] atol=1e-5
        @test vec(result.solution.reserves.up_shortfall)≈[2.0,2.0] atol=1e-5
    end
    c=spinning_case()
    g=OpenSHOP._river_replace(only(c.system.generators);qmin=1.0,pmin=0.380628,
        pmax=8.829,qcurvature=0.8,min_efficiency=0.0,ramp_up=nothing,ramp_down=nothing,
        discharge_ramp_up=nothing,discharge_ramp_down=nothing,reserve_up_energy=nothing,reserve_down_energy=nothing)
    group=OpenSHOP._river_replace(only(c.reserve_groups);activation_hours=0.0,up_requirement=1.0,
        down_requirement=0.0,up_max=6.0,down_max=0.0,up_price=100.0)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]),
        reserve_groups=[group],prices=zeros(2))
    result=solve_case(c;u=ones(Int,1,2),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        x=result.solution
        @test x.reserves.up_q[1,1]<9.0
        @test x.reserves.up[1,1,1]>power(g,g.qmax,100.0)-x.power[1,1]+1.0
        @test OpenSHOP.validate_reserves(c,x)["valid"]
    end
    missing=OpenSHOP._river_replace(g;reserve_up_energy=nothing,reserve_down_energy=nothing)
    invalid=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[missing]),
        reserve_groups=[OpenSHOP._river_replace(group;activation_hours=1.0)])
    @test_throws ArgumentError OpenSHOP.validate_inputs(invalid)
end

@testset "Reserve controls reconstruct on changing heads and refine independently" begin
    c=spinning_case()
    lake=OpenSHOP._river_replace(only(c.system.reservoirs);z0=95.0,slope=0.5)
    group=OpenSHOP._river_replace(only(c.reserve_groups);up_requirement=0.5,down_requirement=0.5)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;reservoirs=[lake]),reserve_groups=[group])
    q=fill(6.0,1,2);u=ones(Int,1,2)
    controls=ReserveAllocation(up=fill(0.5,1,1,2),down=fill(0.5,1,1,2),generator_q=q)
    x=dispatch_from_controls(c;commitment=u,generator_q=q,gate=zeros(0,2),reserves=controls)
    @test validate(c,x)["valid"]
    @test x.reserves.up_q[1,1]!=x.reserves.up_q[1,2]
    @test controls.up_q==q # reconstruction owns the inferred witnesses
    @test replay_audit(c,x)["valid"]
    @test replay_audit(c,x;grid=collect(0.0:0.25:2.0))["valid"]
    profile=OpenSHOP.jsonready(x)
    @test validate(c,profile)["valid"]
    next=restart_case(c,x,1.0)
    continuation=dispatch_from_controls(next;commitment=ones(Int,1,1),generator_q=q[:,2:2],
        gate=zeros(0,1),reserves=OpenSHOP._reserve_slice(x.reserves,[2]))
    @test continuation.diagnostics["validation"]["valid"]
    @test continuation.storage[:,end]≈x.storage[:,end] atol=1e-10
    proposal=OpenSHOP.propose_commitment(c;time_limit=30.0)
    @test proposal["u"]==ones(Int,1,2)
    qmin=fill(2.0,1,2)
    impossible=ReserveAllocation(up=zeros(1,1,2),down=fill(0.5,1,1,2),generator_q=qmin)
    @test_throws ArgumentError dispatch_from_controls(c;commitment=u,generator_q=qmin,
        gate=zeros(0,2),reserves=impossible)
end

@testset "Maintenance removes spinning eligibility without energy accounting slack" begin
    c=spinning_case(;soft=100.0)
    c=OpenSHOP._river_replace(c;operations=[OperationalSeries(object=:Unit,attribute=:maintenance,
        times=[0.0],values=[1.0])])
    result=solve_case(c;u=zeros(Int,1,2),time_limit=30.0)
    @test result.accepted
    if result.solution!==nothing
        x=result.solution
        @test x.power==zeros(1,2)
        @test x.reserves.up==zeros(1,1,2)
        @test x.reserves.down==zeros(1,1,2)
        @test x.objective≈-800.0 atol=1e-5
    end
    invalid=OpenSHOP._river_replace(c;reserve_groups=[ReserveGroup(name=:PumpReserve,generators=[:Pump])])
    @test_throws ArgumentError OpenSHOP.validate_inputs(invalid)
end

@testset "Reserve operating points are equalities and share station response" begin
    c=spinning_case()
    group=OpenSHOP._river_replace(only(c.reserve_groups);up_requirement=0.5,down_requirement=0.5)
    c=OpenSHOP._river_replace(c;reserve_groups=[group])
    q=fill(6.0,1,2)
    x=dispatch_from_controls(c;commitment=ones(Int,1,2),generator_q=q,gate=zeros(0,2),
        reserves=ReserveAllocation(up=fill(0.5,1,1,2),down=fill(0.5,1,1,2),generator_q=q))
    bad=schedule_dict(x);r=deepcopy(x.reserves)
    r.up_q.=8.0;r.down_q.=4.0;bad["reserves"]=r
    @test !validate(c,bad)["valid"]
    audit=OpenSHOP.validate_reserves(c,OpenSHOP._schedule_input(bad;case=c))
    @test audit["residuals"]["Unit_up_witness"]>1.0
    @test audit["residuals"]["Unit_down_witness"]>1.0

    g=only(c.system.generators)
    first=OpenSHOP._river_replace(g;ramp_up=100.0,discharge_ramp_up=100.0)
    second=OpenSHOP._river_replace(g;name=:Second,ramp_up=100.0,discharge_ramp_up=0.01)
    plant=OpenSHOP._river_replace(only(c.system.plants);pmax=17.658,
        tailwater_curve=TableCurve([0.0,20.0],[0.0,4.0]))
    fast=ReserveGroup(name=:Fast,generators=[:Unit],response_hours=0.01,
        up_price=100.0,up_max=1.0,down_max=0.0)
    slow=ReserveGroup(name=:Slow,generators=[:Second],response_hours=1.0,up_max=0.0,down_max=0.0)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        generators=[first,second],plants=[plant]),reserve_groups=[fast,slow])
    q=fill(4.0,2,2);up=zeros(2,2,2);up[1,1,:].=0.5
    x=dispatch_from_controls(c;commitment=ones(Int,2,2),generator_q=q,gate=zeros(0,2),
        reserves=ReserveAllocation(up=up,down=zeros(2,2,2),generator_q=q))
    @test x.reserves.up_q[2,1]-q[2,1]>second.discharge_ramp_up*fast.response_hours
    @test !validate(c,x)["valid"]
    @test OpenSHOP.validate_reserves(c,x)["residuals"]["Second_up_flow_response"]>0.001
    c=OpenSHOP._river_replace(c;operations=[OperationalSeries(object=name,
        attribute=:discharge,times=[0.0],values=[4.0]) for name in (:Unit,:Second)])
    result=solve_case(c;u=ones(Int,2,2),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test maximum(result.solution.reserves.up)<0.05
        @test OpenSHOP.validate_reserves(c,result.solution)["valid"]
    end
end

@testset "Downward capacity reconstructs at a smooth maximum-power point" begin
    c=spinning_case()
    g=OpenSHOP._river_replace(only(c.system.generators);qmin=1.0,pmin=0.380628,
        qcurvature=0.8,ramp_up=nothing,ramp_down=nothing,
        discharge_ramp_up=nothing,discharge_ramp_down=nothing,
        reserve_up_energy=nothing,reserve_down_energy=nothing)
    group=ReserveGroup(name=:Down,generators=[:Unit],up_max=0.0,
        down_requirement=0.5,down_max=0.5)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]),
        reserve_groups=[group])
    peak=(20.0+sqrt(437.5))/6.0
    q=fill(peak,1,2)
    x=dispatch_from_controls(c;commitment=ones(Int,1,2),generator_q=q,gate=zeros(0,2),
        reserves=ReserveAllocation(up=zeros(1,1,2),down=fill(0.5,1,1,2),generator_q=q))
    @test validate(c,x)["valid"]
    @test x.reserves.down_q[1,1]≈5.4991873 atol=1e-3
    @test power(g,x.reserves.down_q[1,1],x.head[1,1])≈x.power[1,1]-0.5 atol=1e-4
    @test replay_audit(c,x)["valid"]
end

@testset "Finer reserve replay retains physical discharge response limits" begin
    c=spinning_case()
    lake=OpenSHOP._river_replace(only(c.system.reservoirs);z0=0.0,slope=10.0)
    group=OpenSHOP._river_replace(only(c.reserve_groups);up_requirement=0.5,down_requirement=0.5)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;reservoirs=[lake]),
        reserve_groups=[group])
    q=fill(6.0,1,2)
    x=dispatch_from_controls(c;commitment=ones(Int,1,2),generator_q=q,gate=zeros(0,2),
        reserves=ReserveAllocation(up=fill(0.5,1,1,2),down=fill(0.5,1,1,2),generator_q=q))
    rate=(maximum(x.reserves.up_q-q)+3e-5)/group.response_hours
    g=OpenSHOP._river_replace(only(c.system.generators);discharge_ramp_up=rate)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]))
    @test validate(c,x)["valid"]
    audit=replay_audit(c,x)
    @test !audit["valid"]
    @test audit["violations"]["reserve_Unit_up_flow_response"]>1e-4
end

@testset "Nonmonotonic reserve inversion preserves an admitted supplied root" begin
    c=spinning_case()
    table=TurbineTable([50.0,100.0,150.0],[1.0,4.0,6.0,8.0,10.0],
        repeat(reshape([0.8,0.8,0.1,0.8,0.8],:,1),1,3),fill(4.0,3),fill(10.0,3))
    g=OpenSHOP._river_replace(only(c.system.generators);qmin=1.0,pmin=0.1,
        turbine_table=table,ramp_up=nothing,ramp_down=nothing,
        discharge_ramp_up=nothing,discharge_ramp_down=nothing,
        reserve_up_energy=nothing,reserve_down_energy=nothing)
    group=ReserveGroup(name=:Down,generators=[:Unit],down_max=0.5,up_max=0.0)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]),
        reserve_groups=[group])
    q=fill(7.0,1,2)
    nominal=dispatch_from_controls(OpenSHOP._river_replace(c;reserve_groups=ReserveGroup[]),
        ones(Int,1,2),q,zeros(0,2))
    down=reshape(nominal.power .- 3.0,1,1,2)
    controls=ReserveAllocation(up=zeros(1,1,2),down=down,generator_q=q)
    # This nearby physical branch is distinct from an inadmissible root below q=4.
    controls.down_q.=6.98
    x=dispatch_from_controls(c;commitment=ones(Int,1,2),generator_q=q,gate=zeros(0,2),reserves=controls)
    @test validate(c,x)["valid"]
    @test all(x.reserves.down_q .> 6.9)
    @test x.power[1,1]-x.reserves.down[1,1,1]≈3.0 atol=1e-8
    controls.down_q.=q
    reconstructed=dispatch_from_controls(c;commitment=ones(Int,1,2),generator_q=q,
        gate=zeros(0,2),reserves=controls)
    @test validate(c,reconstructed)["valid"]
    @test all(reconstructed.reserves.down_q .> 6.9)
end
@testset "Generation reserves coexist with pumping and storage penalties" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=0.01,v0=100.0,vmin=90.0,vmax=110.0,
        water_value=0.0)
    lower=Reservoir(name=:Lower,z0=-1.0,v0=1.0,vmin=0.9,vmax=2.0,water_value=0.0)
    upper=Reservoir(name=:Upper,z0=99.0,v0=1.0,vmin=0.9,vmax=2.0,water_value=0.0)
    generator=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=10.0,pmin=0.1,pmax=10.0,
        hmin=90.0,hmax=110.0,efficiency=0.9,qbest=5.5,hbest=100.0,qcurvature=0.0,hcurvature=0.0,
        initial_on=0,minup=0.0,mindown=0.0,startup=0.0,shutdown=0.0)
    pump=Pump(name=:Pump,plant=:Lift,qmin=1.0,qmax=5.0,pmin=0.1,pmax=10.0,
        hmin=90.0,hmax=110.0,efficiency=0.9,motor_efficiency=0.9,qbest=3.0,hbest=100.0,
        qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0,startup=0.0,shutdown=0.0)
    system=HydroSystem(reservoirs=[lake,lower,upper],generators=[generator,pump],
        boundaries=[Boundary(name=:Sea,head=0.0)],plants=[
            Plant(name=:Station,source=:Lake,target=:Sea,pmax=10.0,ramp=nothing),
            Plant(name=:Lift,source=:Lower,target=:Upper,pmax=10.0,ramp=nothing)])
    group=ReserveGroup(name=:Reserve,generators=[:Unit],up_requirement=1.0,
        down_requirement=0.5,up_max=1.0,down_max=0.5,up_price=20.0,down_price=10.0)
    target=SoftConstraint(name=:StorageTarget,object=:Lake,quantity=:volume,time=2.0,
        lower=100.0,max_lower_slack=0.05,lower_cost=1000.0)
    c=ScheduleCase(name="pump and generator reserve",system=system,grid=[0.0,1.0,2.0],
        prices=[50.0,50.0],reserve_groups=[group],soft_constraints=[target],operations=[
            OperationalSeries(object=:Unit,attribute=:discharge,times=[0.0],values=[5.5]),
            OperationalSeries(object=:Pump,attribute=:discharge,times=[0.0],values=[3.0])])
    q=[5.5 5.5;3.0 3.0];up=zeros(2,1,2);down=copy(up)
    up[1,1,:].=1.0;down[1,1,:].=0.5
    capacity=ReserveAllocation(;up,down,generator_q=q)
    x=dispatch_from_controls(c,ones(Int,2,2),q,zeros(0,2);reserves=capacity)
    generated=[0.00981*5.5*0.9*(101.0-0.000099*(2t-1)) for t in 1:2]
    consumed=[0.00981*3.0*(100.0+0.0108*(2t-1))/0.81 for t in 1:2]
    expected=50.0*sum(generated-consumed)+2.0*(20.0+10.0*0.5)-39.6
    @test x.power[1,:]≈generated atol=1e-9
    @test x.power[2,:]≈consumed atol=1e-9
    @test x.reserves.up_q[2,:]==x.reserves.down_q[2,:]==q[2,:]
    @test x.soft_penalty_cost≈39.6 atol=1e-8
    @test x.objective≈expected atol=1e-8
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    continued=restart_case(c,x,1.0)
    remainder=dispatch_from_controls(continued,ones(Int,2,1),q[:,2:2],zeros(0,1);
        reserves=ReserveAllocation(up=up[:,:,2:2],down=down[:,:,2:2],generator_q=q[:,2:2]))
    @test validate(continued,remainder)["valid"]
    @test remainder.reserves.up_q[2,1]==3.0
    @test remainder.soft_penalty_cost≈39.6 atol=1e-8
    for result in (solve_case(c;u=ones(Int,2,2),warm=x,time_limit=45.0),
        solve(c;initial=x,time_limit=90.0,relative_gap=1e-4))
        @test result.accepted
        if result.solution!==nothing
            @test result.solution.objective≈expected atol=2e-3
            @test validate(c,result.solution)["valid"]
            @test replay_audit(c,result.solution)["valid"]
        end
    end
    invalid=case_dict(c);invalid["reserve_groups"][1]["generators"]=["Pump"]
    @test_throws ArgumentError case_from_dict(invalid)
    bad=ScheduleCase(name=c.name,system=c.system,grid=c.grid,prices=c.prices,
        operations=c.operations,soft_constraints=c.soft_constraints,
        reserve_groups=[ReserveGroup(name=:BadReserve,generators=[:Pump])])
    @test !OpenSHOP.validate_reserves(bad,x)["valid"]
    @test !validate(bad,x)["valid"]
    @test_throws ArgumentError dispatch_from_controls(bad,ones(Int,2,2),q,zeros(0,2);reserves=capacity)
    collision=case_dict(c);collision["soft_constraints"][1]["name"]="Reserve"
    @test_throws ArgumentError case_from_dict(collision)

    # A selectable intake changes nominal head from the supplied reference.
    # Adding an independent reserve group must not pin a nonparticipating
    # pump or generating station to its reference head in the proposal.
    for mode in (:pump,:generate),with_reserves in (false,true)
        profile=case_dict(c)
        profile["junctions"]=[Dict("name"=>"Intake","hmin"=>90.0,"hmax"=>110.0)]
        profile["tunnels"]=[Dict("name"=>"IntakeTunnel","source"=>"Lake","target"=>"Intake",
            "resistance"=>0.1,"capacity"=>12.0,"opening"=>1.0,"gate_positions"=>[0.0,1.0])]
        profile["plants"][1]["source"]="Intake"
        if mode==:generate
            profile["generators"][2]["mode"]="generate"
            profile["generators"][2]["motor_efficiency"]=1.0
            profile["plants"][2]["source"]="Upper"
            profile["plants"][2]["target"]="Lower"
        end
        with_reserves || (profile["reserve_groups"]=[])
        variant=case_from_dict(profile)
        proposal=OpenSHOP.propose_commitment(variant;time_limit=20.0)
        @test haskey(proposal,"u")
        if haskey(proposal,"u")
            @test proposal["u"]==ones(Int,2,2)
        end
    end
end
