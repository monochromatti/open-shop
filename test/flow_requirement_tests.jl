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
    for mode in (:baseline,:tightened)
        result=OpenSHOP.solve(c;initial=good,time_limit=15.0,formulation=mode)
        @test result["accepted"]
        @test result["start_audit"]["valid"]
    end
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

@testset "Exact turbine cell ranges and affine graphs" begin
    # This polynomial has an interior maximum at discharge coordinate 1/2.
    @test OpenSHOP._global_polynomial_range((0.5,1.0,-1.0,0.0),0.0,1.0)==(0.5,0.75)
    @test OpenSHOP._global_polynomial_range((0.0,-1.0,0.0,1.0),-1.0,1.0)[2]≈2/(3sqrt(3))
    a=(0.6,0.4,-0.5,0.1); b=(0.8,-0.2,0.3,-0.05)
    lo,hi=OpenSHOP._global_turbine_cell_range(a,b,-0.5,1.5,-0.3,1.2)
    @test all(range(-0.5,1.5;length=61)) do t
        all(range(-0.3,1.2;length=31)) do u
            v=sum((a[k]+u*(b[k]-a[k]))*t^(k-1) for k in 1:4)
            lo<=v<=hi
        end
    end
    curve=TableCurve([0.0,1.0,2.0],[0.9,0.95,0.8])
    for (lo,hi) in ((0.2,0.8),(1.0,1.0),(-1.0,0.0))
        m=Model(); @variable(m,x)
        y=OpenSHOP._global_table!(m,curve,x,lo,hi)
        @test num_variables(m)==2
        @test JuMP.lower_bound(y)<=JuMP.upper_bound(y)
    end
    # Isolated off cells and on cells follow commitment, preserving head at off.
    table=TurbineTable([60.0,100.0],[2.0,4.0,8.0],[0.8 0.85;0.9 0.95;0.85 0.9],[2.0,2.0],[8.0,8.0])
    m=Model(); @variable(m,0<=q<=8);@variable(m,60<=h<=100);@variable(m,u,Bin)
    eta=OpenSHOP._global_turbine!(m,table,q,h,0.0,8.0,60.0,100.0;
        name=:eta,commitment=u,min_on_flow=3.0)
    cells=m.ext[:global_turbine_cells]["eta"]
    @test any(pair->pair[1]==(0.0,0.0),cells)
    @test all(pair->pair[1]==(0.0,0.0)||pair[1][1]>=3.0,cells)
    @test lower_bound(eta)>0.0
    @test upper_bound(eta)<1.0
end

@testset "Certificates require a resolved root when LPs were solved" begin
    missing=Dict("available"=>true,"lp_iterations"=>100,"first_root_lp_upper_bound"=>nothing)
    @test OpenSHOP._unresolved_root_certificate(missing,"OPTIMAL")
    @test !OpenSHOP._unresolved_root_certificate(merge(missing,Dict("lp_iterations"=>0)),"OPTIMAL")
    @test !OpenSHOP._unresolved_root_certificate(merge(missing,Dict("first_root_lp_upper_bound"=>10.0)),"OPTIMAL")
    @test !OpenSHOP._unresolved_root_certificate(missing,"TIME_LIMIT")
end
