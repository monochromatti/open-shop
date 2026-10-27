@testset "Objective bridge starts include constants, sense and scaling" begin
    for scale in (1.0,10000.0), nonlinear in (false,true), sense in (JuMP.MOI.MIN_SENSE,JuMP.MOI.MAX_SENSE)
        m=Model();set_silent(m)
        @variable(m,0<=x<=1)
        objective=nonlinear ? (3+x-x^2)/scale : (3+x)/scale
        set_objective_sense(m,sense)
        set_objective_function(m,objective)
        set_start_value(x,0.5)
        set_optimizer(m,optimizer_with_attributes(SCIP.Optimizer,
            "limits/time"=>10.0,"limits/nodes"=>0,"numerics/feastol"=>1e-8))
        JuMP.MOI.Utilities.attach_optimizer(backend(m))
        audit=OpenSHOP._complete_scip_start!(m)
        @test audit["valid"]
        @test isempty(audit["missing"])
        @test audit["scaled_objective"]≈(nonlinear ? 3.25 : 3.5)/scale atol=1e-14
        if nonlinear
            opt=backend(m).optimizer
            F=JuMP.MOI.get(backend(m),JuMP.MOI.ObjectiveFunctionType())
            bridge=JuMP.MOI.Bridges.bridge(opt,JuMP.MOI.ObjectiveFunction{F}())
            @test JuMP.MOI.get(JuMP.MOI.Bridges.recursive_model(opt),
                JuMP.MOI.VariablePrimalStart(),bridge.slack)≈3.25/scale atol=1e-14
        end
        optimize!(m)
        diagnostics=OpenSHOP._scip_diagnostics(m)
        @test termination_status(m)!=JuMP.MOI.TIME_LIMIT
        @test diagnostics["native_incumbent_available"]
        native_objective=diagnostics["native_incumbent_objective"]/10000
        @test sense==JuMP.MOI.MAX_SENSE ? native_objective≥audit["scaled_objective"]-1e-8 :
            native_objective≤audit["scaled_objective"]+1e-8
        if nonlinear && sense==JuMP.MOI.MAX_SENSE
            @test native_objective≈audit["scaled_objective"] atol=1e-8
        end
    end
end

@testset "A native completeness failure clears every submitted start" begin
    m=Model();set_silent(m)
    @variable(m,0<=x<=1);@variable(m,0<=unused<=1)
    @objective(m,Max,x-x^2)
    set_start_value(x,0.5)
    set_optimizer(m,SCIP.Optimizer)
    JuMP.MOI.Utilities.attach_optimizer(backend(m))
    @test !OpenSHOP._complete_scip_start!(m)["valid"]
    OpenSHOP._clear_scip_start!(m)
    @test all(JuMP.MOI.get(unsafe_backend(m),JuMP.MOI.VariablePrimalStart(),v)===nothing
        for v in JuMP.MOI.get(unsafe_backend(m),JuMP.MOI.ListOfVariableIndices()))
end

function river_value_start_case()
    c=zero_minimum_case()
    g=only(c.system.generators)
    table=OpenSHOP.synthetic_turbine_table(g;head_derating=0.0)
    river=River(name=:Spill,source=:Lake,target=:Sea,law=:controlled,capacity=1.0,
        curves=[RiverRouting.DelayCurve(0.0,[0.0,0.5,1.0],[0.5,0.5]),
            RiverRouting.DelayCurve(1.0,[0.0,1.0,2.0],[0.2,0.8])],
        water_value=50.0,history_grid=[-1.0,0.0],history_release=[0.0])
    system=OpenSHOP._river_replace(c.system;
        generators=[OpenSHOP._river_replace(g;turbine_table=table)],rivers=[river])
    OpenSHOP._river_replace(c;system)
end

@testset "Nonlinear river inventory objective retains its native start with power cuts" begin
    c=river_value_start_case()
    warm=dispatch_from_controls(c,ones(Int,1,2),fill(2.0,1,2),fill(0.3,1,2))
    @test validate(c,warm)["valid"]
    result=solve(c;initial=warm,time_limit=45.0)
    @test result.accepted
    @test result["start_audit"]["valid"]
    @test result["native_start_audit"]["valid"]
    @test result["native_start_audit"]["objective_slacks"]==1
    @test result["native_start_audit"]["scaled_objective"]≈warm.objective/10000 atol=1e-10
    @test result["scip_diagnostics"]["native_incumbent_available"]
    @test result["scip_diagnostics"]["native_incumbent_objective"]≥warm.objective-1e-5
    @test result["power_cut_statistics"]!==nothing
end

@testset "Absent and rejected starts never acquire an audited native start" begin
    c=zero_minimum_case()
    warm=dispatch_from_controls(c,ones(Int,1,2),fill(2.0,1,2),zeros(0,2))
    for result in (solve(c;preparation_time_limit=0.0,time_limit=30.0),
        solve(c;initial=warm,fixed_u=zeros(Int,1,2),time_limit=30.0))
        @test result.accepted
        @test !haskey(result,"native_start_audit")
        @test validate(c,result.solution)["valid"]
    end
end

# A controllable monotone clock makes exhaustion deterministic. Preparation is
# supplied as an ordinary callable returning the same public result type.
function exhausted_solve_clock()
    called=Ref(false)
    ()->begin
        first_call=!called[]
        called[]=true
        first_call ? 0.0 : 31.0
    end
end

@testset "An exhausted preparation budget never enters global construction" begin
    c=zero_minimum_case()
    warm=dispatch_from_controls(c,ones(Int,1,2),fill(2.0,1,2),zeros(0,2))
    for incumbent in (nothing,warm)
        elapsed=Ref(0.0)
        clock=()->elapsed[]
        prepare=(case;kwargs...)->begin
            elapsed[]=31.0
            SolveResult(Dict("status"=>"PREPARED",
                "accepted"=>incumbent!==nothing);solution=incumbent)
        end
        result=OpenSHOP._solve(c,clock,prepare;
            time_limit=30.0,preparation_time_limit=10.0)
        @test result.status=="PREPARATION_BUDGET_EXHAUSTED"
        @test result["budget_exhausted_phase"]=="preparation"
        @test result["construction_seconds"]==0.0
        @test result["solve_seconds"]==0.0
        @test !haskey(result,"model_profile")
        @test result.total_seconds==31.0
        @test result["budget_overrun_seconds"]==1.0
        @test result.global_bound===nothing && result.absolute_gap===nothing && result.relative_gap===nothing
        @test !result.global_certificate
        @test result.accepted==(incumbent!==nothing)
        if incumbent===nothing
            @test result.solution===nothing
            @test result.feasible_lower_bound===nothing
        else
            @test result.solution.objective≈warm.objective atol=1e-8
            @test validate(c,result.solution)["valid"]
            @test result["replay_audit"]["valid"]
        end
    end
    elapsed=Ref(0.0)
    prepare=(case;kwargs...)->begin
        elapsed[]=31.0
        SolveResult(Dict("status"=>"PREPARED","accepted"=>true);solution=warm)
    end
    unchecked=OpenSHOP._solve(c,()->elapsed[],prepare;time_limit=30.0,
        preparation_time_limit=10.0,replay=false)
    @test unchecked.solution!==nothing && !unchecked.accepted
    @test validate(c,unchecked.solution)["valid"]
    @test unchecked.global_bound===nothing && !unchecked.global_certificate
end

@testset "Deadline exhaustion preserves initial validation and fixed-control scope" begin
    c=zero_minimum_case()
    warm=dispatch_from_controls(c,ones(Int,1,2),fill(2.0,1,2),zeros(0,2))
    unused_prepare=(case;kwargs...)->error("preparation must be skipped for a supplied initial schedule")
    expired=OpenSHOP._solve(c,()->31.0,unused_prepare,0.0;time_limit=30.0)
    @test expired.solution===nothing && !expired.accepted
    @test expired.status=="INITIALIZATION_BUDGET_EXHAUSTED"
    @test expired["construction_seconds"]==0.0 && expired["solve_seconds"]==0.0
    @test expired.global_bound===nothing && !expired.global_certificate
    result=OpenSHOP._solve(c,exhausted_solve_clock(),unused_prepare;
        initial=warm,time_limit=30.0,fixed_u=zeros(Int,1,2))
    @test result.status=="INITIALIZATION_BUDGET_EXHAUSTED"
    @test result.solution===nothing && !result.accepted
    @test result.global_bound===nothing && !result.global_certificate
    malformed=deepcopy(warm)
    malformed.power[1,1]+=1.0
    @test_throws ArgumentError OpenSHOP._solve(c,exhausted_solve_clock(),unused_prepare;
        initial=malformed,time_limit=30.0)
    @test_throws ArgumentError OpenSHOP._solve(c,exhausted_solve_clock(),unused_prepare;
        initial=warm,time_limit=30.0,fixed_u=fill(2,1,2))
    system=OpenSHOP._river_replace(c.system;
        boundaries=[c.system.boundaries;Boundary(name=:Pool,head=0.0)],
        tunnels=[Tunnel(name=:Gate,source=:Sea,target=:Pool,resistance=1.0,capacity=1.0,
            opening=1.0,gate_positions=[0.0,1.0])])
    gated=OpenSHOP._river_replace(c;system)
    gate_warm=dispatch_from_controls(gated,ones(Int,1,2),fill(2.0,1,2),zeros(0,2);
        tunnel_opening=ones(1,2))
    restricted=OpenSHOP._solve(gated,exhausted_solve_clock(),unused_prepare;
        initial=gate_warm,time_limit=30.0,fixed_tunnel_opening=zeros(1,2))
    @test restricted.solution===nothing && !restricted.accepted
    @test restricted.global_bound===nothing && !restricted.global_certificate
end

@testset "A rejected global embedding preserves a physically accepted schedule" begin
    c=zero_minimum_case()
    u=ones(Int,1,length(c.prices))
    warm=dispatch_from_controls(c,u,fill(2.0,size(u)),zeros(0,length(c.prices)))
    unit=only(c.system.generators)
    # The physical acceptance tolerance is less strict than the native start's
    # row tolerance. This point is acceptable for use, but not for a certificate.
    unit=OpenSHOP._river_replace(unit;pmin=minimum(warm.power)+2e-5)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[unit]))
    @test replay_audit(c,warm)["valid"]
    calls=Ref(0)
    clock=()->begin
        calls[]+=1
        calls[]>=6 ? 10.0 : 0.0
    end
    prepare=(c;kwargs...)->error("supplied start needs no preparation")
    result=OpenSHOP._solve(c,clock,prepare,0.0;initial=warm,time_limit=5.0,preparation_time_limit=0.0)
    @test result.accepted
    @test result.solution!==nothing
    @test result.solution.power≈warm.power atol=1e-9
    @test replay_audit(c,result.solution)["valid"]
    @test result.feasible_lower_bound===nothing
    @test result.relative_gap===nothing && result.absolute_gap===nothing
    @test !result.global_certificate
end
