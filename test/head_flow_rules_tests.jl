function head_rule_case(;unit_kwargs=NamedTuple(),plant_kwargs=NamedTuple(),operations=OperationalSeries[])
    lake=Reservoir(name=:Lake,z0=90.0,slope=1.0,v0=10.0,vmin=8.0,vmax=12.0,water_value=0.0)
    unit=Generator(;name=:Unit,plant=:Station,qmin=1.0,qmax=4.0,pmin=0.1,pmax=10.0,
        hmin=80.0,hmax=110.0,hbest=100.0,qbest=2.0,efficiency=0.9,qcurvature=0.0,hcurvature=0.0,
        initial_on=0,minup=0.0,mindown=0.0,startup=0.0,unit_kwargs...)
    plant=Plant(;name=:Station,source=:Lake,target=:Sea,pmax=10.0,ramp=nothing,plant_kwargs...)
    ScheduleCase(;name="head restrictions",system=HydroSystem(reservoirs=[lake],
        boundaries=[Boundary(name=:Sea,head=0.0)],plants=[plant],generators=[unit]),
        grid=[0.0,1.0,2.0],prices=[100.0,100.0],operations)
end

@testset "Startup gross head uses the starting surface, not average turbine head" begin
    c=head_rule_case(unit_kwargs=(startup_gross_head_min=100.0,))
    x=dispatch_from_controls(c,ones(Int,1,2),fill(4.0,1,2),zeros(0,2))
    @test validate(c,x)["valid"]
    @test x.head[1,1]<100.0
    @test replay_audit(c,x)["valid"]
    continued=restart_case(c,x,1.0)
    y=dispatch_from_controls(continued,ones(Int,1,1),fill(4.0,1,1),zeros(0,1))
    @test validate(continued,y)["valid"] # Continuing below the startup threshold is allowed.
    @test y.storage≈x.storage[:,2:3]
    strict=head_rule_case(unit_kwargs=(startup_gross_head_min=100.01,))
    @test !validate(strict,x)["valid"]
    local_result=solve_case(c;u=ones(Int,1,2),warm=x,time_limit=20.0)
    @test local_result.accepted
    @test local_result.solution.generator_q≈fill(4.0,1,2) atol=1e-3
    off=dispatch_from_controls(strict,zeros(Int,1,2),zeros(1,2),zeros(0,2))
    global_result=OpenSHOP.solve(strict;initial=off,time_limit=30.0,preparation_time_limit=0.0)
    @test global_result.accepted
    @test global_result.global_certificate
    @test global_result.solution.commitment==zeros(Int,1,2)
    @test abs(global_result.solution.objective)<1e-4
    curved=OpenSHOP._river_replace(strict;system=OpenSHOP._river_replace(strict.system;
        reservoirs=[OpenSHOP._river_replace(only(strict.system.reservoirs);curvature=0.1)],
        generators=[OpenSHOP._river_replace(only(strict.system.generators);startup_gross_head_min=110.2)]))
    proposal=OpenSHOP.propose_commitment(curved;time_limit=10.0)
    @test haskey(proposal,"u")
    @test proposal["u"][1,1]==0
end

@testset "Gross head references and time-dependent running limits" begin
    c=head_rule_case(unit_kwargs=(gross_head_min=99.99,gross_head_max=100.0,))
    x=dispatch_from_controls(c,ones(Int,1,2),ones(1,2),zeros(0,2))
    @test validate(c,x)["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    restricted=OpenSHOP._river_replace(c;operations=[OperationalSeries(object=:Unit,
        attribute=:gross_head_min,times=[1.0],values=[100.0])])
    @test !validate(restricted,x)["valid"]
    off=dispatch_from_controls(restricted,zeros(Int,1,2),zeros(1,2),zeros(0,2))
    @test validate(restricted,off)["valid"]
    p=only(c.system.plants)
    with_junction=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        junctions=[Junction(name=:Intake,hmin=80.0,hmax=110.0)],
        tunnels=[Tunnel(name=:Penstock,source=:Lake,target=:Intake,resistance=0.5,capacity=5.0)],
        plants=[OpenSHOP._river_replace(p;source=:Intake,gross_head_source=:Lake)]))
    y=dispatch_from_controls(with_junction,ones(Int,1,2),ones(1,2),zeros(0,2))
    @test validate(with_junction,y)["valid"]
    @test y.head[2,1]≈y.head[1,1]-0.5 atol=1e-7
    @test case_dict(case_from_dict(case_dict(with_junction)))==case_dict(with_junction)
    ambiguous=OpenSHOP._river_replace(with_junction;system=OpenSHOP._river_replace(with_junction.system;
        plants=[OpenSHOP._river_replace(only(with_junction.system.plants);gross_head_source=nothing)]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(ambiguous)
end

function flow_rule_case(;river_kwargs=NamedTuple(),rules=FlowRequirement[],operations=OperationalSeries[])
    river=River(;name=:Reach,source=:Upper,target=:Lower,law=:controlled,capacity=4.0,
        curves=RiverRouting.DelayCurve[],deterministic_delay=0.5,water_value=0.0,
        history_grid=[-1.0,0.0],history_release=[0.0],river_kwargs...)
    reservoirs=[Reservoir(name=n,z0=z,slope=0.01,v0=10.0,vmin=8.0,vmax=12.0,water_value=0.0)
        for (n,z) in ((:Upper,100.0),(:Lower,0.0))]
    ScheduleCase(;name="river operating rules",system=HydroSystem(;reservoirs,rivers=[river]),
        grid=[0.0,1.0,2.0],prices=zeros(2),flow_requirements=rules,operations)
end

@testset "Flow observations distinguish delayed arrivals and upstream release" begin
    c=flow_rule_case(rules=[FlowRequirement(name=:Gauge,river_arrivals=[:Reach],max_flow=3.0,
        min_volume=0.0216,max_volume=0.0216)])
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),ones(1,2))
    @test !validate(c,x)["valid"] # Arrival rates are 2 then 4, despite the same releases.
    @test validate(c,x)["residuals"]["flow_requirement_Gauge"]≈1.0
    valid=OpenSHOP._river_replace(c;flow_requirements=[OpenSHOP._river_replace(only(c.flow_requirements);
        max_flow=1.5,min_volume=0.0108,max_volume=0.0108)])
    valid=OpenSHOP._river_replace(valid;system=OpenSHOP._river_replace(valid.system;
        reservoirs=[OpenSHOP._river_replace(r;water_value=r.name==:Upper ? 1000.0 : 0.0) for r in valid.system.reservoirs]))
    x=dispatch_from_controls(valid,zeros(Int,0,2),zeros(0,2),reshape([0.75,0.0],1,2))
    @test validate(valid,x)["valid"]
    @test replay_audit(valid,x)["valid"]
    @test x.arrival_volume≈[0.0054 0.0054] atol=1e-10
    continued=restart_case(valid,x,1.0)
    @test only(continued.flow_requirements).initial_volume≈0.0054
    y=dispatch_from_controls(continued,zeros(Int,0,1),zeros(0,1),zeros(1,1))
    @test validate(continued,y)["valid"]
    @test y.storage≈x.storage[:,2:3] atol=1e-9
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
    local_result=solve_case(valid;warm=x,time_limit=20.0)
    @test local_result.accepted
    @test sum(local_result.solution.arrival_volume)≈0.0108 atol=1e-6
    @test local_result.solution.river_release≈[3.0 0.0] atol=1e-3
    proposal=OpenSHOP.propose_commitment(valid;time_limit=10.0)
    @test haskey(proposal,"river_release")
    global_result=OpenSHOP.solve(valid;initial=x,time_limit=30.0,preparation_time_limit=0.0)
    @test global_result.accepted
    @test validate(valid,global_result.solution)["valid"]
    @test global_result.global_certificate
    @test global_result.solution.objective≈-10.8 atol=1e-3
end

@testset "River costs count rate changes and transitions without a duration multiplier" begin
    c=flow_rule_case(river_kwargs=(min_allowed_flow=1.0,flow_change_cost=2.0,
        flow_start_cost=7.0,flow_stop_cost=3.0,))
    c=OpenSHOP._river_replace(c;grid=[0.0,0.5,2.0],prices=zeros(2))
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),reshape([0.5,0.0],1,2))
    @test validate(c,x)["valid"]
    @test x.objective≈-18.0 atol=1e-9 # 2*(|2-0|+|0-2|)+7+3
    costs=objective_components(c,x)
    @test (costs.river_change_cost,costs.river_start_cost,costs.river_stop_cost)==(8.0,7.0,3.0)
    @test costs.total≈x.objective
    @test replay_audit(c,x)["valid"]
    restarted=restart_case(c,x,0.5)
    y=dispatch_from_controls(restarted,zeros(Int,0,1),zeros(0,1),zeros(1,1))
    @test validate(restarted,y)["valid"]
    @test y.objective≈-7.0
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    @test_throws ArgumentError solve_verified(c;u=zeros(Int,0,2),initial=x,
        river_on=zeros(Int,1,2),max_refinements=0)
    @test_throws ArgumentError solve_verified(c;u=zeros(Int,0,2),initial=x,
        river_on=ones(Int,1,1),max_refinements=0)
    gap=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),fill(0.1,1,2))
    @test !validate(c,gap)["valid"]
    local_result=solve_case(c;river_on=reshape([1,0],1,2),warm=x,time_limit=20.0)
    @test local_result.accepted
    @test local_result.solution.objective≈-14.0 atol=1e-3 # Flow 1, then zero.
    global_result=OpenSHOP.solve(c;initial=x,time_limit=30.0,preparation_time_limit=0.0)
    @test global_result.accepted
    @test global_result.solution.objective≈0.0 atol=1e-3
    @test global_result["start_audit"]["valid"]
end

@testset "A river startup cost changes the preferred release schedule" begin
    rule=FlowRequirement(name=:Release,rivers=[:Reach],min_volume=0.0036)
    c=flow_rule_case(river_kwargs=(min_allowed_flow=1.0,flow_start_cost=10.0,),rules=[rule])
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),reshape([0.25,0.0],1,2))
    @test validate(c,x)["valid"]
    global_result=OpenSHOP.solve(c;initial=x,time_limit=30.0,preparation_time_limit=0.0)
    @test global_result.accepted
    @test global_result.solution.objective≈-10.0 atol=1e-3
    @test sum(global_result.solution.river_release)>=1.0-1e-5
    prepared=schedule_case(c;proposal_time_limit=5.0,nlp_time_limit=15.0,
        max_refinements=0,time_limit=60.0)
    @test prepared["accepted"]
    @test prepared["solution"].objective≈-10.0 atol=1e-3
end

@testset "A storage restriction at the terminal timestamp is enforced" begin
    c=head_rule_case()
    restricted=OpenSHOP._river_replace(c;operations=[OperationalSeries(object=:Lake,
        attribute=:vmin,times=[2.0],values=[10.0])])
    x=dispatch_from_controls(c,ones(Int,1,2),ones(1,2),zeros(0,2))
    @test !validate(restricted,x)["valid"]
    result=OpenSHOP.solve(restricted;time_limit=30.0,preparation_time_limit=0.0)
    @test result.accepted
    @test result.solution.storage[1,end]>=10.0-1e-5
    @test result.solution.objective≈0.0 atol=1e-3
end

@testset "Rate-change costs favor a constant release without introducing flow commitment" begin
    c=flow_rule_case(river_kwargs=(flow_change_cost=2.0,),rules=[FlowRequirement(
        name=:Allowance,rivers=[:Reach],min_volume=0.0036,max_volume=0.0036)])
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),reshape([0.25,0.0],1,2))
    @test x.objective≈-4.0
    local_result=solve_case(c;warm=x,time_limit=20.0)
    @test local_result.accepted
    @test local_result.solution.river_release≈fill(0.5,1,2) atol=1e-3
    @test local_result.solution.objective≈-1.0 atol=1e-3
    global_result=OpenSHOP.solve(c;initial=x,time_limit=30.0,preparation_time_limit=0.0)
    @test global_result.accepted
    @test global_result.solution.objective≈-1.0 atol=1e-3
    @test global_result["native_start_audit"]["valid"]
end

@testset "Pumping startup head uses positive lift at the starting surface" begin
    c=head_rule_case(unit_kwargs=(mode=:pump,startup_gross_head_min=100.005,))
    p=only(c.system.plants)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        plants=[OpenSHOP._river_replace(p;source=:Sea,target=:Lake)]))
    x=dispatch_from_controls(c,ones(Int,1,2),fill(4.0,1,2),zeros(0,2))
    @test !validate(c,x)["valid"] # Mean lift exceeds the threshold, but initial lift does not.
    @test x.head[1,1]>100.005
    relaxed=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        generators=[OpenSHOP._river_replace(only(c.system.generators);startup_gross_head_min=100.0)]))
    @test validate(relaxed,x)["valid"]
    result=solve_case(relaxed;u=ones(Int,1,2),warm=x,time_limit=20.0)
    @test result.accepted
end

@testset "Signed tunnel observations and invalid operating inputs" begin
    c=flow_rule_case()
    link=Tunnel(name=:Return,source=:Lower,target=:Upper,resistance=100.0,capacity=2.0)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;tunnels=[link]),
        flow_requirements=[FlowRequirement(name=:SignedGauge,tunnels=[:Return],min_flow=-2.0,max_flow=-0.5)])
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),zeros(1,2))
    @test validate(c,x)["valid"]
    @test all(<(0),x.tunnel_q)
    @test abs(sum(x.storage[:,end])-sum(x.storage[:,1]))<1e-9
    strict=OpenSHOP._river_replace(c;flow_requirements=[FlowRequirement(
        name=:SignedGauge,tunnels=[:Return],min_flow=-2.0,max_flow=-1.1)])
    @test !validate(strict,x)["valid"]
    for r in (OpenSHOP._river_replace(only(c.system.rivers);flow_start_cost=1.0),
        OpenSHOP._river_replace(only(c.system.rivers);min_allowed_flow=0.0),
        OpenSHOP._river_replace(only(c.system.rivers);min_allowed_flow=1.0,initial_release=0.5))
        @test_throws ArgumentError OpenSHOP.validate_inputs(OpenSHOP._river_replace(c;
            system=OpenSHOP._river_replace(c.system;rivers=[r])))
    end
    for rule in (FlowRequirement(name=:Invalid,river_arrivals=[:Missing]),
        FlowRequirement(name=:Invalid,max_flow=-1.0),
        FlowRequirement(name=:Invalid,min_volume=2.0,max_volume=1.0))
        @test_throws ArgumentError OpenSHOP.validate_inputs(OpenSHOP._river_replace(c;flow_requirements=[rule]))
    end
end
