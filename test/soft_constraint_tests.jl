include(joinpath(@__DIR__,"..","examples","soft_limits.jl"))
using .SoftLimitsExample: soft_limit_case
soft_copy(x;kwargs...)=typeof(x)(;merge((;(k=>getfield(x,k) for k in fieldnames(typeof(x)))...), (;kwargs...))...)

@testset "A terminal target has an explicit dispatch tradeoff" begin
    c=soft_limit_case();u=ones(Int,1,1);q=fill(5.0,1,1);gate=zeros(0,1)
    x=dispatch_from_controls(c,u,q,gate)
    # Exact conservation and midpoint level: V=1-.0036*5, H=100-.0018*5.
    revenue=40.0*0.00981*0.9*5.0*(100.0-0.0018*5.0)
    @test x.storage[1,end]≈0.982 atol=1e-12
    @test x.soft_slacks[1].lower≈0.018 atol=1e-12
    @test x.soft_slacks[1].unit==:Mm3
    @test x.soft_penalty_cost≈18.0 atol=1e-10
    @test x.objective≈revenue-18.0 atol=1e-9
    @test validate(c,x)["valid"]
    @test !isempty(validate(c,x)["permitted_penalties"])
    @test replay_audit(c,x)["valid"]
    local_result=solve_case(c;u,warm=x,time_limit=60.0)
    @test local_result.accepted
    if local_result.solution!==nothing
        @test local_result.solution.objective≈revenue-18.0 atol=1e-3
    end
    for (case,expected_q,expected_cost) in ((c,5.0,18.0),
        (soft_limit_case(penalty=20000.0),0.0,0.0),(soft_limit_case(allowance=0.0),0.0,0.0))
        initial=dispatch_from_controls(case,zeros(Int,1,1),zeros(1,1),gate)
        result=solve(case;initial,time_limit=60.0,relative_gap=1e-4)
        @test result.accepted
        reference_objective=expected_q==0.0 ? 0.0 : revenue-expected_cost
        if result.global_bound===nothing
            @test !result.global_certificate
        else
            @test isfinite(result.global_bound) && result.global_bound>=reference_objective-1e-3
        end
        if result.solution!==nothing
            @test result.solution.generator_q[1,1]≈expected_q atol=1e-4
            @test result.solution.soft_penalty_cost≈expected_cost atol=1e-3
            @test replay_audit(case,result.solution)["valid"]
        end
    end
    too_small=soft_limit_case(allowance=0.01)
    failed=dispatch_from_controls(too_small,u,fill(4.0,1,1),gate)
    @test !validate(too_small,failed)["valid"]
    @test validate(too_small,failed)["residuals"]["soft_allowance_TerminalStorage"]>0.004
    physical=schedule_dict(x);physical["V"][1,end]+=0.001
    @test !validate(c,physical)["valid"]
    @test validate(c,physical)["residuals"]["balance_Lake"]>0.0009
    forged=schedule_dict(x);forged["soft_penalty_cost"]=0.0
    @test !validate(c,forged)["valid"]
    missing=schedule_dict(x);delete!(missing,"soft_slacks")
    @test !validate(c,missing)["valid"]
    wrong_grid=schedule_dict(x);wrong_grid["grid"]=[0.0,2.0]
    @test !validate(c,wrong_grid)["valid"]
    # A generous tactical allowance cannot enlarge the static safety domain.
    safer=soft_copy(c;system=soft_copy(c.system;
        reservoirs=[soft_copy(only(c.system.reservoirs);vmin=0.99)]),
        soft_constraints=[soft_copy(only(c.soft_constraints);max_lower_slack=0.05)])
    initial=dispatch_from_controls(safer,u,fill(2.0,1,1),gate)
    bounded=solve_case(safer;u,warm=initial,time_limit=45.0)
    @test bounded.accepted
    if bounded.solution!==nothing
        @test bounded.solution.storage[1,end]≈0.99 atol=1e-6
        @test bounded.solution.generator_q[1,1]≈(1.0-0.99)/0.0036 atol=1e-4
    end
end

@testset "Soft output targets do not force unit states or soften hard schedules" begin
    base=soft_limit_case();gate=zeros(0,1)
    power_target=SoftConstraint(name=:PowerTarget,object=:Station,quantity=:power,time=0.0,end_time=1.0,
        upper=0.0,max_upper_slack=5.0,upper_cost=10.0)
    c=soft_copy(base;soft_constraints=[power_target])
    x=dispatch_from_controls(c,ones(Int,1,1),fill(5.0,1,1),gate)
    @test validate(c,x)["valid"]
    @test x.soft_slacks[1].unit==:MW
    @test x.soft_penalty_cost≈10.0*x.power[1,1] atol=1e-10
    result=solve(c;initial=x,time_limit=60.0,relative_gap=1e-4)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.generator_q[1,1]≈5.0 atol=1e-4
        @test result.solution.objective≈30.0*x.power[1,1] atol=1e-3
    end
    positive=SoftConstraint(name=:DesiredProduction,object=:Unit,quantity=:power,time=0.0,end_time=1.0,
        lower=2.0,max_lower_slack=2.0,lower_cost=1.0)
    costly=soft_copy(only(base.system.generators);startup=1000.0)
    offcase=soft_copy(base;system=soft_copy(base.system;generators=[costly]),soft_constraints=[positive])
    off=dispatch_from_controls(offcase,zeros(Int,1,1),zeros(1,1),gate)
    @test validate(offcase,off)["valid"]
    @test off.soft_penalty_cost≈2.0 atol=1e-12
    result=solve(offcase;initial=off,time_limit=60.0,relative_gap=1e-4)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.commitment==zeros(Int,1,1)
        @test result.solution.objective≈-2.0 atol=1e-8
    end
    # The same target must admit the off state during cold-start preparation.
    cold=solve(offcase;time_limit=60.0,preparation_time_limit=20.0,relative_gap=1e-4)
    @test cold.accepted
    if cold.solution!==nothing
        @test cold.solution.commitment==zeros(Int,1,1)
        @test cold.solution.objective≈-2.0 atol=1e-8
    end
    flow=SoftConstraint(name=:FlowTarget,object=:Unit,quantity=:discharge,time=0.0,end_time=1.0,
        upper=3.0,max_upper_slack=2.0,upper_cost=1000.0)
    flowcase=soft_copy(base;soft_constraints=[flow])
    high=dispatch_from_controls(flowcase,ones(Int,1,1),fill(5.0,1,1),gate)
    @test high.soft_slacks[1].upper≈2.0 atol=1e-12
    @test high.soft_penalty_cost≈7.2 atol=1e-10
    hard=soft_copy(flowcase;operations=[OperationalSeries(object=:Unit,attribute=:discharge,times=[0.0],values=[3.0])])
    bad=dispatch_from_controls(hard,ones(Int,1,1),fill(5.0,1,1),gate)
    @test !validate(hard,bad)["valid"]
    solved=solve_case(hard;u=ones(Int,1,1),time_limit=30.0)
    @test solved.accepted
    if solved.solution!==nothing
        @test solved.solution.generator_q[1,1]≈3.0 atol=1e-5
        @test solved.solution.soft_penalty_cost≈0.0 atol=1e-4
    end
    # Even a generous soft power band cannot waive a static generator limit.
    excessive=dispatch_from_controls(c,ones(Int,1,1),fill(7.0,1,1),gate)
    @test !validate(c,excessive)["valid"]
    @test validate(c,excessive)["residuals"]["unit_bounds_Unit"]>1.0
end

@testset "Level slacks use the state vertex and preserve exact table calibration" begin
    base=soft_limit_case();r=only(base.system.reservoirs)
    curve=TableCurve([0.95,1.0,2.0],[99.9,100.0,104.0])
    reservoir=soft_copy(r;z0=98.0,slope=2.0,level_curve=curve)
    target=SoftConstraint(name=:LevelTarget,object=:Lake,quantity=:level,time=1.0,
        lower=100.0,max_lower_slack=0.036,lower_cost=1000.0)
    c=soft_copy(base;system=soft_copy(base.system;reservoirs=[reservoir]),soft_constraints=[target])
    x=dispatch_from_controls(c,ones(Int,1,1),fill(5.0,1,1),zeros(0,1))
    @test x.soft_slacks[1].unit==:metres
    @test x.soft_slacks[1].lower≈0.036 atol=1e-10
    @test x.soft_penalty_cost≈36.0 atol=1e-8
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    local_result=solve_case(c;u=ones(Int,1,1),warm=x,time_limit=45.0)
    @test local_result.accepted
    result=solve(c;initial=x,time_limit=60.0,relative_gap=1e-4)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.objective≈x.objective atol=1e-3
    end
end

@testset "Assessment times, rate windows, restart and JSON preserve penalty costs" begin
    base=soft_limit_case()
    point=SoftConstraint(name=:OperatingStorage,object=:Lake,quantity=:volume,time=1.0,
        lower=0.99,max_lower_slack=0.008,lower_cost=500.0)
    endpoint=soft_copy(only(base.soft_constraints);time=2.0,max_lower_slack=0.036)
    c=soft_copy(base;grid=[0.0,1.0,2.0],prices=[40.0,40.0],soft_constraints=[point,endpoint])
    x=dispatch_from_controls(c,ones(Int,1,2),fill(5.0,1,2),zeros(0,2))
    @test x.soft_penalty_cost≈40.0 atol=1e-9
    @test replay_audit(c,x)["valid"]
    fine=with_grid(c,collect(0.0:0.25:2.0))
    xf=dispatch_from_controls(fine,ones(Int,1,8),fill(5.0,1,8),zeros(0,8))
    @test xf.soft_penalty_cost≈x.soft_penalty_cost atol=1e-8
    @test xf.objective≈x.objective atol=1e-8
    prefix=soft_copy(c;grid=[0.0,1.0],prices=[40.0])
    past=dispatch_from_controls(prefix,ones(Int,1,1),fill(5.0,1,1),zeros(0,1))
    continued=restart_case(c,x,1.0)
    rest=dispatch_from_controls(continued,ones(Int,1,1),fill(5.0,1,1),zeros(0,1))
    @test rest.soft_penalty_cost≈36.0 atol=1e-9
    @test x.objective≈past.objective+rest.objective atol=1e-8
    @test rest.storage≈x.storage[:,2:end] atol=1e-10
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
    flow=SoftConstraint(name=:AverageFlow,object=:Station,quantity=:discharge,time=0.0,end_time=2.0,
        lower=3.0,upper=3.0,max_lower_slack=3.0,max_upper_slack=2.0,lower_cost=1000.0,upper_cost=1000.0)
    window=soft_copy(c;soft_constraints=[flow])
    rates=reshape([1.0,5.0],1,:)
    balanced=dispatch_from_controls(window,ones(Int,1,2),rates,zeros(0,2))
    @test balanced.soft_penalty_cost≈0.0 atol=1e-10
    @test replay_audit(window,balanced)["valid"]
    @test_throws ArgumentError restart_case(window,balanced,1.0)
    first_window=soft_copy(flow;name=:FirstWindow,end_time=1.0)
    second_window=soft_copy(flow;name=:SecondWindow,time=1.0)
    split_windows=soft_copy(c;soft_constraints=[first_window,second_window])
    full=dispatch_from_controls(split_windows,ones(Int,1,2),rates,zeros(0,2))
    before=dispatch_from_controls(soft_copy(split_windows;grid=[0.0,1.0],prices=[40.0]),
        ones(Int,1,1),rates[:,1:1],zeros(0,1))
    after=dispatch_from_controls(restart_case(split_windows,full,1.0),
        ones(Int,1,1),rates[:,2:2],zeros(0,1))
    @test full.soft_penalty_cost≈14.4 atol=1e-9
    @test before.soft_penalty_cost≈7.2 atol=1e-9
    @test after.soft_penalty_cost≈7.2 atol=1e-9
    @test full.objective≈before.objective+after.objective atol=1e-8
    mktempdir() do dir
        path=joinpath(dir,"soft.json");writejson(path,case_dict(c))
        @test case_dict(readcase(path))==case_dict(c)
        output=joinpath(dir,"schedule.json");writejson(output,schedule_dict(x))
        loaded=OpenSHOP.JSON3.read(read(output,String),Dict{String,Any})
        @test validate(c,loaded)["valid"]
        @test loaded["soft_penalty_cost"]≈40.0 atol=1e-9
    end
end

@testset "Plant schedules average summed outputs over original unequal windows" begin
    base=soft_limit_case()
    second=soft_copy(only(base.system.generators);name=:Second)
    system=soft_copy(base.system;generators=[only(base.system.generators),second],
        plants=[soft_copy(only(base.system.plants);pmax=10.0)])
    volume=SoftConstraint(name=:StorageUpper,object=:Lake,quantity=:volume,time=2.0,
        upper=0.95,max_upper_slack=0.02,upper_cost=500.0)
    flow=SoftConstraint(name=:PlantFlow,object=:Station,quantity=:discharge,time=0.0,end_time=2.0,
        lower=6.0,upper=6.0,max_lower_slack=2.0,max_upper_slack=2.0,
        lower_cost=1000.0,upper_cost=1000.0)
    unit=SoftConstraint(name=:UnitPower,object=:Second,quantity=:power,time=0.0,end_time=0.5,
        upper=0.0,max_upper_slack=5.0,upper_cost=10.0)
    c=soft_copy(base;system,grid=[0.0,0.5,2.0],prices=[40.0,40.0],
        soft_constraints=[volume,flow,unit])
    rates=[1.0 3.0;2.0 4.0]
    x=dispatch_from_controls(c,ones(Int,2,2),rates,zeros(0,2))
    # Integrated aggregate discharge = .5*(1+2)+1.5*(3+4) = 12 m³/s-hours.
    assessment=soft_constraint_assessment(c,x)
    @test x.storage[1,end]≈1.0-0.0036*12.0 atol=1e-12
    @test assessment.slacks[1].upper≈0.0068 atol=1e-12
    @test assessment.slacks[1].penalty_cost≈3.4 atol=1e-9
    @test assessment.slacks[2].lower≈0.0 atol=1e-12
    @test assessment.slacks[2].upper≈0.0 atol=1e-12
    @test assessment.slacks[3].penalty_cost≈10.0*0.5*x.power[2,1] atol=1e-9
    @test validate(c,x)["valid"]
    fine=with_grid(c,[0.0,0.25,0.5,1.0,1.5,2.0])
    xf=dispatch_from_controls(fine,ones(Int,2,5),rates[:,[1,1,2,2,2]],zeros(0,5))
    @test xf.soft_penalty_cost≈x.soft_penalty_cost atol=1e-8
    @test replay_audit(c,x)["valid"]
end

@testset "Malformed or unsupported soft restrictions are rejected" begin
    c=soft_limit_case();rule=only(c.soft_constraints)
    for bad in (soft_copy(rule;object=:Missing),soft_copy(rule;quantity=:inflow),
        soft_copy(rule;max_lower_slack=Inf),soft_copy(rule;max_lower_slack=-1.0),
        soft_copy(rule;lower_cost=-1.0),soft_copy(rule;lower_cost=0.0),
        soft_copy(rule;lower=NaN),soft_copy(rule;end_time=2.0),soft_copy(rule;time=0.5),
        soft_copy(rule;upper=0.5),soft_copy(rule;max_lower_slack=1e308,lower_cost=1e308))
        @test_throws ArgumentError OpenSHOP.validate_inputs(soft_copy(c;soft_constraints=[bad]))
    end
    @test_throws ArgumentError OpenSHOP.validate_inputs(soft_copy(c;soft_constraints=[rule,rule]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(soft_copy(c;
        soft_constraints=[soft_copy(rule;name=:Unit)]))
    output=SoftConstraint(name=:Power,object=:Unit,quantity=:power,time=0.0,end_time=1.0,
        upper=1.0,max_upper_slack=1.0,upper_cost=1.0)
    for bad in (soft_copy(output;end_time=nothing),soft_copy(output;end_time=Inf),
        soft_copy(output;end_time=0.0),soft_copy(output;end_time=0.5),
        soft_copy(output;time=-1.0),soft_copy(output;lower_cost=1.0))
        @test_throws ArgumentError OpenSHOP.validate_inputs(soft_copy(c;soft_constraints=[bad]))
    end
end
