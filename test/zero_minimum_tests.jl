function zero_minimum_case(;pmin=0.0,qmin=0.0,pmax=5.0,qmax=5.0,minup=0.0,operations=OperationalSeries[])
    lake=Reservoir(name=:Lake,z0=99.0,v0=1.0,vmin=0.5,vmax=2.0,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Station,pmin=pmin,qmin=qmin,pmax=pmax,qmax=qmax,
        hmin=90.0,hmax=110.0,efficiency=0.9,qbest=3.0,hbest=100.0,
        qcurvature=0.0,hcurvature=0.0,initial_on=0,minup=minup,mindown=0.0,
        startup=3.0,shutdown=2.0)
    system=HydroSystem(reservoirs=[lake],generators=[unit],
        boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=5.0,ramp=nothing)])
    ScheduleCase(;name="zero unit minima",system=system,grid=[0.0,1.0,2.0],prices=[50.0,50.0],operations)
end

@testset "Zero unit minima preserve conversion and optimized dispatch" begin
    for (pmin,qmin) in ((0.0,1.0),(0.1,0.0),(0.0,0.0))
        c=zero_minimum_case(;pmin,qmin)
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        warm=dispatch_from_controls(c,ones(Int,1,2),fill(2.0,1,2),zeros(0,2))
        @test validate(c,warm)["valid"]
        @test warm.power≈[0.00981*2*0.9*(100-0.0036*(2t-1)) for t in 1:2]' atol=1e-9
        for result in (solve_case(c;u=ones(Int,1,2),warm,time_limit=30.0),
            solve(c;initial=warm,time_limit=45.0,relative_gap=1e-4))
            @test result.accepted
            if result.solution!==nothing
                @test result.solution.generator_q≈fill(5.0,1,2) atol=1e-4
                @test validate(c,result.solution)["valid"]
                @test replay_audit(c,result.solution)["valid"]
            end
        end
    end
end

@testset "Idle on-state is distinct from prescribed zero and maintenance" begin
    c=zero_minimum_case(minup=2.0)
    idle=dispatch_from_controls(c,ones(Int,1,2),zeros(1,2),zeros(0,2))
    off=dispatch_from_controls(c,zeros(Int,1,2),zeros(1,2),zeros(0,2))
    @test validate(c,idle)["valid"] && validate(c,off)["valid"]
    @test idle.storage==off.storage
    @test idle.objective==-3.0 && off.objective==0.0
    @test replay_audit(c,idle)["valid"]
    restarted=restart_case(c,idle,1.0)
    @test only(restarted.system.generators).initial_on==1
    @test validate(restarted,dispatch_from_controls(restarted,ones(Int,1,1),zeros(1,1),zeros(0,1)))["valid"]
    @test_throws ArgumentError dispatch_from_controls(restarted,zeros(Int,1,1),zeros(1,1),zeros(0,1))
    for attribute in (:power,:discharge,:maintenance)
        value=attribute==:maintenance ? 1.0 : 0.0
        restricted=zero_minimum_case(operations=[OperationalSeries(;object=:Unit,attribute,
            times=[0.0],values=[value])])
        initial=dispatch_from_controls(restricted,zeros(Int,1,2),zeros(1,2),zeros(0,2))
        @test_throws ArgumentError dispatch_from_controls(restricted,ones(Int,1,2),zeros(1,2),zeros(0,2))
        result=solve(restricted;initial,time_limit=30.0)
        @test result.accepted
        @test result.solution.commitment==zeros(Int,1,2)
    end
end

@testset "Unit bounds remain finite, nonnegative and ordered" begin
    for kwargs in ((;qmin=-1.0),(;pmin=-1.0),(;qmax=-1.0),(;pmax=-1.0),
        (;qmin=NaN),(;pmax=Inf))
        @test_throws ErrorException case_from_dict(case_dict(zero_minimum_case(;kwargs...)))
    end
end

@testset "Zero upper bounds do not create a nonfinite commitment proposal" begin
    c=zero_minimum_case(pmax=0.0,qmax=0.0,operations=[
        OperationalSeries(object=:Unit,attribute=:forced_on,times=[0.0],values=[1.0])])
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    proposal=OpenSHOP.propose_commitment(c;time_limit=15.0)
    @test haskey(proposal,"u")
    @test proposal["u"]==ones(Int,1,2)
    warm=dispatch_from_controls(c,ones(Int,1,2),zeros(1,2),zeros(0,2))
    for result in (solve_case(c;u=ones(Int,1,2),warm,time_limit=30.0),
        solve(c;initial=warm,time_limit=30.0))
        @test result.accepted
        @test result.solution.generator_q≈zeros(1,2) atol=1e-8
        @test result.solution.power≈zeros(1,2) atol=1e-8
        @test replay_audit(c,result.solution)["valid"]
    end
end

@testset "An idle committed unit may have a negative permitted head" begin
    for table in (nothing, TurbineTable([1.0,10.0],[0.0,1.0],fill(0.9,2,2),
        [0.0,0.0],[1.0,1.0];interpolation=:bilinear,head_extrapolation=:linear))
        system=HydroSystem(
            junctions=[Junction(name=:Intake,hmin=90.0,hmax=110.0)],
            boundaries=[Boundary(name=:Supply,head=100.0),Boundary(name=:Tail,head=105.0)],
            tunnels=[Tunnel(name=:Penstock,source=:Supply,target=:Intake,resistance=1.0,capacity=1.0)],
            plants=[Plant(name=:Station,source=:Intake,target=:Tail,pmax=1.0,ramp=nothing)],
            generators=[Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=1.0,pmin=0.0,pmax=1.0,
                hmin=-10.0,hmax=10.0,hbest=5.0,qbest=1.0,efficiency=0.9,
                qcurvature=0.0,hcurvature=0.0,initial_on=1,minup=0.0,mindown=0.0,
                startup=0.0,turbine_table=table)])
        c=ScheduleCase(name="idle negative head",system=system,grid=[0.0,1.0],prices=[10.0],
            operations=[OperationalSeries(object=:Unit,attribute=:forced_on,times=[0.0],values=[1.0])])
        warm=dispatch_from_controls(c,ones(Int,1,1),zeros(1,1),zeros(0,1))
        @test validate(c,warm)["valid"]
        for result in (solve_case(c;u=ones(Int,1,1),warm,time_limit=30.0),
            solve(c;initial=warm,time_limit=45.0))
            @test result.accepted
            if result.solution!==nothing
                @test result.solution.commitment==ones(Int,1,1)
                @test result.solution.power≈zeros(1,1) atol=1e-7
                @test result.solution.generator_q≈zeros(1,1) atol=1e-7
                @test result.solution.head[1,1]≈100.0 atol=1e-6
                @test validate(c,result.solution)["valid"]
            end
        end
    end
end
