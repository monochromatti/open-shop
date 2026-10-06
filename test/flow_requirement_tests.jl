@testset "Aggregate flow observations" begin
    base=analytic_global_fixture(:on).case
    river=River(name=:Bypass,source=:AnalyticLake,target=:AnalyticTail,
        curves=RiverRouting.DelayCurve[],capacity=10.0,law=:controlled,
        deterministic_delay=0.0,water_value=0.0,history_grid=[-1.0,0.0],history_release=[0.0])
    sys=OpenSHOP._river_replace(base.system;rivers=[river])
    rule=FlowRequirement(name=:Reach,generators=[:AnalyticUnit],rivers=[:Bypass],min_flow=8.0)
    c=OpenSHOP._river_replace(base;system=sys,flow_requirements=[rule],operations=[
        OperationalSeries(object=:Reach,attribute=:inflow,times=[0.0],values=[1.0])])
    @test OpenSHOP.validate_inputs(c)
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    good=dispatch_from_controls(c,ones(Int,1,1),fill(5.0,1,1),fill(0.2,1,1))
    @test good["validation"]["valid"]
    @test replay_audit(c,good)["valid"]
    bad=dispatch_from_controls(c,ones(Int,1,1),fill(5.0,1,1),zeros(1,1))
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["flow_requirement_Reach"]≈2.0
    local_result=solve_case(c;u=ones(Int,1,1),time_limit=10.0)
    @test local_result["validation"]["valid"]
    @test local_result["generator_q"][1,1]+local_result["river_release"][1,1]+1>=8-1e-5
    proposal=OpenSHOP.propose_commitment(c;time_limit=5.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        @test proposal["generator_q"][1,1]+proposal["river_release"][1,1]+1>=8-1e-5
    end
    result=OpenSHOP.solve(c;initial=good,time_limit=15.0)
    @test result["accepted"]
    @test result["start_audit"]["valid"]
    # Refinement and restart retain the observation and its physical input times.
    longer=OpenSHOP._river_replace(c;grid=[0.0,1.0,2.0],prices=[100.0,100.0],operations=[
        OperationalSeries(object=:Reach,attribute=:inflow,times=[0.0,1.0],values=[1.0,2.0]),
        OperationalSeries(object=:Reach,attribute=:min_flow,times=[0.0,1.0],values=[8.0,9.0])])
    schedule=dispatch_from_controls(longer,ones(Int,1,2),fill(5.0,1,2),fill(0.2,1,2))
    @test schedule["validation"]["valid"]
    @test with_grid(longer,[0.0,0.5,1.0,1.5,2.0]).flow_requirements==[rule]
    continued=restart_case(longer,schedule,1.0)
    @test continued.flow_requirements==[rule]
    @test OpenSHOP.opinterval(continued,:Reach,:inflow,1,0.0)==2.0
    @test_throws ArgumentError OpenSHOP.validate_inputs(OpenSHOP._river_replace(c;
        flow_requirements=[FlowRequirement(name=:Reach,rivers=[:Missing])]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(OpenSHOP._river_replace(c;
        flow_requirements=[FlowRequirement(name=:Reach,generators=[:AnalyticUnit,:AnalyticUnit])]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(OpenSHOP._river_replace(c;
        flow_requirements=[FlowRequirement(name=:AnalyticUnit)]))
end

@testset "Certificates require a resolved root when LPs were solved" begin
    missing=Dict("available"=>true,"lp_iterations"=>100,"first_root_lp_upper_bound"=>nothing)
    @test OpenSHOP._unresolved_root_certificate(missing,"OPTIMAL")
    @test !OpenSHOP._unresolved_root_certificate(merge(missing,Dict("lp_iterations"=>0)),"OPTIMAL")
    @test !OpenSHOP._unresolved_root_certificate(merge(missing,Dict("first_root_lp_upper_bound"=>10.0)),"OPTIMAL")
    @test !OpenSHOP._unresolved_root_certificate(missing,"TIME_LIMIT")
end
