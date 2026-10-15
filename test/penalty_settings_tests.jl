function penalty_settings_case(;soft_constraints=SoftConstraint[],reserve_groups=ReserveGroup[],operations=OperationalSeries[])
    reservoir=Reservoir(name=:Lake,z0=99.0,v0=1.0,vmin=0.95,vmax=2.0,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,pmin=0.1,pmax=5.0,
        hmin=90.0,hmax=110.0,efficiency=0.9,qcurvature=0.0,hcurvature=0.0,
        minup=0.0,mindown=0.0,initial_on=0,startup=0.0,shutdown=0.0,
        reserve_available=false)
    system=HydroSystem(reservoirs=[reservoir],generators=[unit],
        plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=5.0,ramp=nothing)],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    ScheduleCase(;name="effective penalty equations",system,grid=[0.0,0.5,2.0],
        prices=[40.0,40.0],soft_constraints,reserve_groups,operations)
end

@testset "Effective state and schedule penalties have explicit units and costs" begin
    rules=[SoftConstraint(name=:OperatingVolume,object=:Lake,quantity=:volume,time=0.5,
            lower=1.0,max_lower_slack=0.005,lower_cost=1000.0),
        SoftConstraint(name=:TerminalLevel,object=:Lake,quantity=:level,time=2.0,
            lower=100.0,max_lower_slack=0.05,lower_cost=200.0),
        SoftConstraint(name=:PowerWindow,object=:Station,quantity=:power,time=0.0,end_time=0.5,
            upper=0.0,max_upper_slack=5.0,upper_cost=10.0),
        SoftConstraint(name=:DischargeWindow,object=:Unit,quantity=:discharge,time=0.0,end_time=2.0,
            upper=2.0,max_upper_slack=1.0,upper_cost=1000.0),
        SoftConstraint(name=:HardVolume,object=:Lake,quantity=:volume,time=2.0,
            lower=0.95,lower_cost=25.0)]
    c=penalty_settings_case(;soft_constraints=rules)
    settings=penalty_settings(c)
    @test settings.defaults.soft_allowance==0.0
    @test settings.defaults.soft_coefficient==0.0
    @test settings.defaults.river_release_coefficient==0.0
    @test settings.defaults.reserve_shortfall_coefficient===nothing
    @test [s.rule for s in settings.soft]==rules
    @test all(s.active for s in settings.soft)
    @test [s.unit for s in settings.soft]==[:Mm3,:metres,:MW,:m3_per_s,:Mm3]
    @test [s.coefficient_unit for s in settings.soft]==
        [:currency_per_Mm3,:currency_per_metre,:currency_per_MWh,:currency_per_Mm3,:currency_per_Mm3]
    @test [s.cost_factor for s in settings.soft]≈[1.0,1.0,0.5,0.0072,1.0] atol=1e-15
    @test settings.soft[1].lower_origin==:explicit_penalty
    @test settings.soft[1].upper_origin==:absent
    @test settings.soft[5].lower_origin==:hard_limit
    # .5*1+1.5*3=5 m³/s-hours withdraws .018 Mm³. The first midpoint
    # head is 99.9991 m. Point costs are charged once, not per grid interval.
    x=dispatch_from_controls(c,ones(Int,1,2),reshape([1.0,3.0],1,:),zeros(0,2))
    assessment=soft_constraint_assessment(c,x)
    expected=[1.8,3.6,10.0*0.5*0.00981*0.9*99.9991,3.6,0.0]
    @test validate(c,x)["valid"]
    @test [s.lower for s in assessment.slacks]≈[0.0018,0.018,0.0,0.0,0.0] atol=1e-10
    @test [s.upper for s in assessment.slacks]≈[0.0,0.0,0.00981*0.9*99.9991,0.5,0.0] atol=1e-10
    @test [s.penalty_cost for s in assessment.slacks]≈expected atol=1e-7
    @test assessment.cost≈sum(expected) atol=1e-7
    for (setting,slack) in zip(settings.soft,assessment.slacks)
        rule=setting.rule
        @test slack.penalty_cost≈setting.cost_factor*(rule.lower_cost*slack.lower+
            rule.upper_cost*slack.upper) atol=1e-10
    end
    @test objective_components(c,x).soft_penalty≈sum(expected) atol=1e-7
    @test objective_components(c,x).total≈x.objective atol=1e-7
    refined=with_grid(c,[0.0,0.25,0.5,1.0,1.5,2.0])
    fine=dispatch_from_controls(refined,ones(Int,1,5),reshape([1.0,1.0,3.0,3.0,3.0],1,:),zeros(0,5))
    @test [s.cost_factor for s in penalty_settings(refined).soft]≈[1.0,1.0,0.5,0.0072,1.0] atol=1e-15
    @test soft_constraint_assessment(refined,fine).cost≈assessment.cost atol=1e-7
    @test replay_audit(c,x)["valid"]
    @test penalty_settings(case_from_dict(case_dict(c)))==settings
    corrupted=schedule_dict(x);corrupted["V"][1,end]+=0.001
    @test !validate(c,corrupted)["valid"]
    @test validate(c,corrupted)["residuals"]["balance_Lake"]>0.0009
    @test_throws ArgumentError restart_case(c,x,0.5) # Original output window is unfinished.
end

@testset "Terminal cost remains attached to its original absolute assessment" begin
    target=SoftConstraint(name=:Terminal,object=:Lake,quantity=:volume,time=2.0,
        lower=1.0,max_lower_slack=0.02,lower_cost=1000.0)
    c=penalty_settings_case(;soft_constraints=[target])
    x=dispatch_from_controls(c,ones(Int,1,2),reshape([1.0,3.0],1,:),zeros(0,2))
    prefix=ScheduleCase(name=c.name,system=c.system,grid=[0.0,0.5],prices=[40.0],soft_constraints=[target])
    before=dispatch_from_controls(prefix,ones(Int,1,1),ones(1,1),zeros(0,1))
    @test !only(penalty_settings(prefix).soft).active
    @test soft_constraint_assessment(prefix,before).cost==0.0
    continued=restart_case(c,x,0.5)
    after=dispatch_from_controls(continued,ones(Int,1,1),fill(3.0,1,1),zeros(0,1))
    @test only(penalty_settings(continued).soft).rule.time==2.0
    @test only(penalty_settings(continued).soft).cost_factor==1.0
    @test soft_constraint_assessment(continued,after).cost≈18.0 atol=1e-7
    @test x.objective≈before.objective+after.objective atol=1e-7
    @test penalty_settings(case_from_dict(case_dict(continued)))==penalty_settings(continued)
end

@testset "River release penalties expose zero defaults and volume conversion" begin
    lake=Reservoir(name=:Lake,z0=99.0,v0=1.0,vmin=0.95,vmax=2.0,water_value=0.0)
    river=River(name=:Reach,source=:Lake,target=:Sea,law=:controlled,capacity=10.0,
        deterministic_delay=0.0,curves=RiverRouting.DelayCurve[],water_value=0.0)
    system=HydroSystem(reservoirs=[lake],rivers=[river],boundaries=[Boundary(name=:Sea,head=0.0)])
    base=ScheduleCase(name="release policy",system=system,grid=[0.0,0.5,2.0],prices=zeros(2))
    @test all(s.coefficient==0.0 && s.coefficient_origin==:default_zero &&
        s.minimum==0.0 && s.max_shortfall==0.0 for s in penalty_settings(base).release)
    operations=[OperationalSeries(object=:Reach,attribute=:min_release,times=[0.0,0.5],values=[2.0,4.0]),
        OperationalSeries(object=:Reach,attribute=:release_penalty,times=[0.0,0.5],values=[0.0,1000.0])]
    c=ScheduleCase(name=base.name,system=system,grid=base.grid,prices=base.prices,operations=operations)
    settings=penalty_settings(c).release
    @test [s.minimum for s in settings]==[2.0,4.0]
    @test all(s.minimum_origin==:operational_series && s.coefficient_origin==:operational_series for s in settings)
    @test [s.max_shortfall for s in settings]==[0.0,4.0]
    @test all(s.unit==:m3_per_s && s.coefficient_unit==:currency_per_Mm3 for s in settings)
    @test [s.cost_factor for s in settings]≈[0.0018,0.0054] atol=1e-15
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),reshape([0.2,0.1],1,:))
    @test validate(c,x)["valid"]
    @test vec(x.shortfall_release)≈[0.0,0.0162] atol=1e-12
    @test objective_components(c,x).release_penalty≈16.2 atol=1e-10
    deficit=max(0.0,settings[2].minimum-x.river_release[1,2])
    @test settings[2].coefficient*settings[2].cost_factor*deficit≈16.2 atol=1e-10
    hard=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),zeros(1,2))
    @test !validate(c,hard)["valid"]
    @test validate(c,hard)["residuals"]["hard_release_requirements"]≈2.0 atol=1e-10
    fine=with_grid(c,[0.0,0.25,0.5,1.0,1.5,2.0])
    xf=dispatch_from_controls(fine,zeros(Int,0,5),zeros(0,5),reshape([0.2,0.2,0.1,0.1,0.1],1,:))
    @test objective_components(fine,xf).release_penalty≈16.2 atol=1e-10
    @test sum(s.cost_factor for s in penalty_settings(fine).release if s.coefficient>0)≈0.0054 atol=1e-15
    @test penalty_settings(case_from_dict(case_dict(c)))==penalty_settings(c)
end

@testset "Reserve shortfall defaults hard and explicit costs accrue per MW-hour" begin
    hard=ReserveGroup(name=:Hard,generators=[:Unit])
    soft=ReserveGroup(name=:Soft,generators=[:Unit],up_requirement=3.0,
        down_requirement=1.0,shortfall_cost=100.0)
    operation=OperationalSeries(object=:Soft,attribute=:up_requirement,times=[0.5],values=[2.0])
    c=penalty_settings_case(;reserve_groups=[hard,soft],operations=[operation])
    settings=penalty_settings(c).reserve
    hard_settings=filter(s->s.rule.name==:Hard,settings)
    soft_settings=filter(s->s.rule.name==:Soft,settings)
    @test all(s.coefficient===nothing && s.coefficient_origin==:default_hard &&
        s.max_up_shortfall==s.max_down_shortfall==0.0 for s in hard_settings)
    @test all(s.coefficient==100.0 && s.coefficient_origin==:explicit_penalty for s in soft_settings)
    @test [s.up_requirement for s in soft_settings]==[3.0,2.0]
    @test [s.up_requirement_origin for s in soft_settings]==[:group_value,:operational_series]
    @test all(s.max_up_shortfall==s.up_requirement && s.max_down_shortfall==1.0 &&
        s.unit==:MW && s.coefficient_unit==:currency_per_MWh for s in soft_settings)
    x=dispatch_from_controls(c,ones(Int,1,2),reshape([1.0,3.0],1,:),zeros(0,2))
    @test validate(c,x)["valid"]
    @test x.reserves.up_shortfall≈[0.0 0.0;3.0 2.0] atol=1e-12
    @test x.reserves.down_shortfall≈[0.0 0.0;1.0 1.0] atol=1e-12
    expected=sum(s.coefficient*s.cost_factor*(s.up_requirement+s.down_requirement) for s in soft_settings)
    @test expected≈650.0 atol=1e-12
    @test objective_components(c,x).reserve_penalty≈expected atol=1e-10
    loaded=case_from_dict(case_dict(c))
    effective=s->(s.rule.name,s.rule.generators,s.time,s.end_time,s.up_requirement,
        s.down_requirement,s.coefficient,s.coefficient_origin,s.cost_factor)
    @test map(effective,penalty_settings(loaded).reserve)==map(effective,settings)
    @test objective_components(loaded,x).reserve_penalty≈expected atol=1e-10
end
