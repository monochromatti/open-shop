include("tunnel_gate_fixtures.jl")

# Independent storage integration, pressure loss, and chronological gate cost.
function tg_reference(c,opening)
    T=length(c.prices)
    V=zeros(2,T+1);V[:,1].=5.0
    power=zeros(1,T)
    prior=[1.0,0.0]
    cost=0.0
    for t in 1:T
        active=only(findall(>(0.0),opening[:,t]))
        V[:,t+1].=V[:,t]
        V[active,t+1]-=0.0036*5.0*(c.grid[t+1]-c.grid[t])
        head=(active==1 ? 99.5 : 79.5)+0.1*(V[active,t]+V[active,t+1])/2-25.0
        power[1,t]=0.00981*0.9*5.0*head
        cost+=5.0*sum(abs.(opening[:,t]-prior))
        prior=opening[:,t]
    end
    (;V,power,cost,objective=sum(c.prices .* diff(c.grid) .* vec(power))-cost)
end

@testset "Gate proposals respect the reference physics of passive pressure links" begin
    base=tunnel_gate_case()
    cross=Tunnel(name=:LakeCross,source=:HighLake,target=:LowLake,
        resistance=0.05,capacity=30.0)
    c=tg_copy(base;system=tg_copy(base.system;tunnels=[base.system.tunnels;cross]))
    # Physical flow is about 20 m³/s. Referencing this link to the station's
    # nominal 5 m³/s would require a surrogate flow beyond its 30 m³/s capacity.
    proposal=OpenSHOP.propose_commitment(c;time_limit=20.0)
    @test haskey(proposal,"u")
    @test haskey(proposal,"tunnel_opening")
    if haskey(proposal,"u")
        opening=proposal["tunnel_opening"]
        @test opening[3,:]==ones(3)
        @test opening[1:2,:]==[1.0 0.0 1.0;0.0 1.0 0.0]
        result=solve_verified(c;u=proposal["u"],tunnel_opening=opening,
            time_limit=10.0,max_refinements=0)
        @test result["accepted"]
        if result["accepted"]
            @test all(q->19.0<q<21.0,result["solution"]["tunnel_q"][3,:])
            @test replay_audit(c,result["solution"])["valid"]
        end
    end
end

function tg_flow_reference(openings,tail;curve=true)
    V=[10.0];q=Float64[]
    for opening in openings
        delta=90.0+last(V)-tail
        if opening==0.0
            flow=0.0
        else
            resistance=curve ? (opening==0.5 ? 4.0 : 1.0) : 1.0/opening
            flow=sign(delta)*(-0.0018+sqrt(0.0018^2+4resistance*abs(delta)))/(2resistance)
        end
        push!(q,flow);push!(V,last(V)-0.0036flow)
    end
    q,V
end

@testset "Discrete tunnel gates: independent enumeration and joint optimum" begin
    c=tunnel_gate_case()
    @test OpenSHOP.validate_inputs(c)
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    choices=[[i,j,k] for i in 1:2,j in 1:2,k in 1:2]
    schedules=[[Float64(choice[t]==e) for e in 1:2,t in 1:3] for choice in vec(choices)]
    expected=[tg_reference(c,a) for a in schedules]
    best_index=argmax(getproperty.(expected,:objective))
    best_opening=schedules[best_index]
    @test best_opening==[1.0 0.0 1.0;0.0 1.0 0.0]
    @test expected[best_index].cost==20.0
    u=ones(Int,1,3);q=fill(5.0,1,3);gates=zeros(0,3)
    local_objectives=Float64[]
    for (opening,reference) in zip(schedules,expected)
        x=dispatch_from_controls(c,u,q,gates;tunnel_opening=opening)
        @test x["validation"]["valid"]
        @test x["tunnel_opening"]==opening
        @test x["V"]≈reference.V atol=1e-9
        @test x["power"]≈reference.power atol=1e-8
        @test x["gate_adjustment_cost"]≈reference.cost atol=1e-10
        @test x["objective"]≈reference.objective atol=1e-8
        local_result=solve_case(c;u,warm=x,tunnel_opening=opening,time_limit=10.0)
        @test local_result["validation"]["valid"]
        @test local_result["objective"]≈reference.objective atol=1e-4
        push!(local_objectives,local_result["objective"])
    end
    @test maximum(local_objectives)≈expected[best_index].objective atol=1e-4
    initial=dispatch_from_controls(c,u,q,gates;tunnel_opening=schedules[1])
    result=OpenSHOP.solve(c;initial,time_limit=60.0,relative_gap=1e-4,absolute_gap=1e-3)
    @test result["accepted"]
    @test result["global_certificate"]
    @test !result["commitment_fixed"]
    if result["solution"]!==nothing
        solution=result["solution"]
        @test solution["tunnel_opening"]==best_opening
        @test solution["objective"]≈expected[best_index].objective atol=1e-3
        @test result["global_bound"]>=expected[best_index].objective-1e-6
        @test replay_audit(c,solution)["valid"]
        @test result["start_audit"]["valid"]
        @test result["start_audit"]["assigned"]==result["start_audit"]["variables"]
    end
    fixed_units=OpenSHOP.solve(c;initial,fixed_u=u,time_limit=20.0,relative_gap=1e-4)
    @test fixed_units["accepted"]
    @test fixed_units["commitment_fixed"]
    @test !fixed_units["tunnel_opening_fixed"]
    @test fixed_units["solution"]["tunnel_opening"]==best_opening
    fixed_gates=OpenSHOP.solve(c;initial,fixed_tunnel_opening=schedules[1],
        time_limit=20.0,relative_gap=1e-4)
    @test fixed_gates["accepted"]
    @test fixed_gates["tunnel_opening_fixed"]
    @test occursin("fixed tunnel openings",fixed_gates["certificate_scope"])
    @test fixed_gates["model_profile"]["tunnel_gate_binaries"]==0
    @test fixed_gates["solution"]["objective"]≈expected[1].objective atol=1e-4
    prepared=schedule_case(c;proposal_time_limit=5.0,nlp_time_limit=10.0,max_refinements=0)
    @test prepared["accepted"]
    @test !prepared["global_optimality_proven"]
    if prepared["accepted"]
        @test prepared["solution"]["tunnel_opening"]==best_opening
        @test prepared["solution"]["objective"]≈expected[best_index].objective atol=1e-3
        @test replay_audit(c,prepared["solution"])["valid"]
    end
    # An off-state zero-head reference must yield finite coefficients. Partial
    # reference dictionaries remain usable as in the prescribed-gate proposal.
    zerohead=tg_copy(c;system=tg_copy(c.system;
        generators=[tg_copy(only(c.system.generators);hmin=0.0)]))
    reference=Dict("H"=>copy(initial["H"]))
    reference["H"][nodeindex(zerohead.system)[:Intake],:].=0.0
    proposal=OpenSHOP.propose_commitment(zerohead;reference,time_limit=5.0)
    @test haskey(proposal,"u")
    @test haskey(proposal,"tunnel_opening")
    if haskey(proposal,"tunnel_opening")
        @test proposal["tunnel_opening"]==best_opening
    end
end

@testset "Gate movement costs, replay, restart and independent tampering" begin
    c=tunnel_gate_case();opening=[1.0 0.0 1.0;0.0 1.0 0.0]
    x=dispatch_from_controls(c,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);tunnel_opening=opening)
    @test replay_audit(c,x)["valid"]
    verified=solve_verified(c;u=ones(Int,1,3),initial=x,tunnel_opening=opening,
        max_refinements=0,time_limit=10.0)
    @test verified["accepted"]
    prefix=tg_copy(c;grid=c.grid[1:2],prices=c.prices[1:1])
    first=dispatch_from_controls(prefix,ones(Int,1,1),fill(5.0,1,1),zeros(0,1);
        tunnel_opening=opening[:,1:1])
    continued=restart_case(c,x,1.0)
    @test [e.initial_opening for e in continued.system.tunnels]==opening[:,1]
    rest=dispatch_from_controls(continued,ones(Int,1,2),fill(5.0,1,2),zeros(0,2);
        tunnel_opening=opening[:,2:3])
    @test rest["validation"]["valid"]
    @test rest["V"]≈x["V"][:,2:end] atol=1e-9
    @test rest["power"]≈x["power"][:,2:end] atol=1e-9
    @test x["objective"]≈first["objective"]+rest["objective"] atol=1e-9
    @test x["gate_adjustment_cost"]≈first["gate_adjustment_cost"]+rest["gate_adjustment_cost"] atol=1e-10
    @test replay_audit(continued,rest)["valid"]
    bad=deepcopy(x);bad["tunnel_opening"][1,2]=0.5
    @test !validate(c,bad)["valid"]
    bad=deepcopy(x);bad["tunnel_opening"][:,2].=1.0
    @test !validate(c,bad)["valid"]
    bad=deepcopy(x);bad["objective"]+=1.0
    @test !validate(c,bad)["valid"]
    bad=deepcopy(x);bad["gate_adjustment_cost"]+=1.0
    @test !validate(c,bad)["valid"]
    irregular=tg_copy(c;grid=[0.0,0.5,2.0,3.0])
    changed=dispatch_from_controls(irregular,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);
        tunnel_opening=opening)
    @test changed["validation"]["valid"]
    @test changed["gate_adjustment_cost"]≈20.0 atol=1e-10
    @test changed["objective"]≈tg_reference(irregular,opening).objective atol=1e-8
    operations=[c.operations;OperationalSeries(object=:HighIntake,attribute=:gate_adjustment_cost,
        times=[1.0,2.0],values=[7.0,9.0])]
    cc=tg_copy(c;operations)
    charged=dispatch_from_controls(cc,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);tunnel_opening=opening)
    @test charged["validation"]["valid"]
    @test charged["gate_adjustment_cost"]≈26.0 atol=1e-10
    @test charged["objective"]≈x["objective"]-6.0 atol=1e-9
end

@testset "Discrete positive positions, signed flow and closed unequal heads" begin
    for tail in (0.0,200.0), curve in (true,false)
        c=tg_flow_case(;tail,curve=curve ? TableCurve([0.0,0.5,1.0],[8.0,3.0,0.0]) : nothing)
        opening=reshape([0.0,0.5,1.0],1,:)
        expected,storage=tg_flow_reference(vec(opening),tail;curve)
        x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);tunnel_opening=opening)
        @test x["validation"]["valid"]
        @test vec(x["tunnel_q"])≈expected atol=1e-8
        @test vec(x["V"])≈storage atol=1e-9
        @test replay_audit(c,x)["valid"]
        @test abs(x["H"][1,1]-x["H"][2,1])>90.0
        b=OpenSHOP._build_global_dispatch(c;warm=x,joint=true)
        @test OpenSHOP._lift_start!(b,c,x)["valid"]
        @test OpenSHOP._lift_start!(b,c,x)["assigned"]==JuMP.num_variables(b.m)
        local_result=solve_case(c;warm=x,tunnel_opening=opening,time_limit=10.0)
        @test local_result["validation"]["valid"]
        @test vec(local_result["tunnel_q"])≈expected atol=1e-5
    end
    c=tg_flow_case(;tail=200.0,opening=0.0,initial=0.0)
    closed=zeros(1,3)
    x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);tunnel_opening=closed)
    @test x["validation"]["valid"]
    @test x["tunnel_q"]≈zeros(1,3) atol=1e-10
    @test x["V"]≈fill(10.0,1,4) atol=1e-10
    local_result=solve_case(c;warm=x,tunnel_opening=closed,time_limit=10.0)
    @test local_result["validation"]["valid"]
    # Changing positive-position resistance exercises the joint quadratic
    # gate-loss representation, with an independently known movement optimum.
    c=tg_flow_case(;tail=200.0,cost=2.0)
    initial=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=reshape([0.0,0.5,1.0],1,:))
    result=OpenSHOP.solve(c;initial,time_limit=30.0,relative_gap=1e-4)
    @test result["accepted"]
    if result["solution"]!==nothing
        @test result["solution"]["tunnel_opening"]==fill(0.5,1,3)
        @test result["objective"]≈0.0 atol=1e-5
        @test all(<(0.0),result["solution"]["tunnel_q"])
        @test replay_audit(c,result["solution"])["valid"]
    end
end

@testset "Tunnel gate input restrictions and prescribed positions" begin
    for positions in ([1.0,0.0],[0.0,0.5,0.5,1.0],[-0.1,0.5,1.0],[0.0,0.5,1.1],[0.0,NaN,1.0])
        @test_throws ArgumentError OpenSHOP.validate_inputs(tg_flow_case(;positions))
    end
    @test_throws ArgumentError OpenSHOP.validate_inputs(tg_flow_case(;opening=0.25))
    @test_throws ArgumentError OpenSHOP.validate_inputs(tg_flow_case(;initial=0.25))
    @test_throws ArgumentError OpenSHOP.validate_inputs(tg_flow_case(;cost=-1.0))
    for interlock in (TunnelInterlock(name=:Bad,tunnels=[:Unknown],max_open=1),
        TunnelInterlock(name=:Bad,tunnels=[:GateLink,:GateLink],max_open=1),
        TunnelInterlock(name=:Bad,tunnels=[:GateLink],min_open=1,max_open=0))
        @test_throws ArgumentError OpenSHOP.validate_inputs(tg_flow_case(;interlocks=[interlock]))
    end
    operation=OperationalSeries(object=:GateLink,attribute=:opening,times=[1.0],values=[0.5])
    c=tg_flow_case(;operations=[operation])
    opening=reshape([1.0,0.5,0.5],1,:)
    x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);tunnel_opening=opening)
    @test x["validation"]["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    @test_throws ArgumentError dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3))
    for bad in (reshape([1.0,1.0,0.5],1,:),fill(0.25,1,3),fill(NaN,1,3))
        @test_throws ArgumentError dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);tunnel_opening=bad)
    end
    @test_throws DimensionMismatch dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);tunnel_opening=zeros(1,2))
    tiny=tg_flow_case(;positions=[0.0,1e-12,0.5,1.0],opening=0.0,
        interlocks=[TunnelInterlock(name=:Closed,tunnels=[:GateLink],max_open=0)])
    @test OpenSHOP.validate_inputs(tiny)
    @test_throws ArgumentError dispatch_from_controls(tiny,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=fill(1e-12,1,3))
    legacy=tg_flow_case(;positions=Float64[],tail=200.0)
    @test JuMP.num_variables(OpenSHOP._build_dispatch(legacy).m)==13
    @test JuMP.num_variables(OpenSHOP._build_global_dispatch(legacy).m)==22
    singleton=tg_flow_case(;positions=[0.5],tail=200.0,cost=2.0)
    @test JuMP.num_variables(OpenSHOP._build_dispatch(singleton).m)==13
    @test JuMP.num_variables(OpenSHOP._build_global_dispatch(singleton).m)==22
    rule=TunnelInterlock(name=:ChangingCount,tunnels=[:GateLink],max_open=1)
    cc=tg_flow_case(;interlocks=[rule],operations=[OperationalSeries(object=:ChangingCount,
        attribute=:max_open,times=[1.0],values=[0.0])])
    controlled=dispatch_from_controls(cc,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=reshape([1.0,0.0,0.0],1,:))
    @test controlled["validation"]["valid"]
    @test case_dict(case_from_dict(case_dict(cc)))==case_dict(cc)
    @test_throws ArgumentError dispatch_from_controls(cc,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=ones(1,3))
end

@testset "Joint gate proposal and scheduling through a loop and delayed confluence" for distributed in (false,true)
    base=channel_network_case(;distributed)
    tunnel=tg_copy(only(base.system.tunnels);gate_positions=[0.0,0.5,1.0],
        initial_opening=0.5,gate_adjustment_cost=0.1)
    c=tg_copy(base;system=tg_copy(base.system;tunnels=[tunnel]))
    opening=reshape([0.5,0.5,0.5,0.0,0.0],1,:)
    x=dispatch_from_controls(c,ones(Int,2,5),fill(1.0,2,5),
        [zeros(3,5);fill(0.5,1,5);zeros(1,5)];tunnel_opening=opening)
    @test x["validation"]["valid"]
    @test replay_audit(c,x)["valid"]
    @test x["terminal_transit"][5]>0.0
    proposal=OpenSHOP.propose_commitment(c;reference=x,time_limit=5.0)
    @test haskey(proposal,"u")
    @test haskey(proposal,"tunnel_opening")
    if haskey(proposal,"tunnel_opening")
        @test size(proposal["tunnel_opening"])==(1,5)
        @test all(in([0.0,0.5,1.0]),proposal["tunnel_opening"])
        @test proposal["tunnel_opening"][1,4:5]==[0.0,0.0]
        @test admissible(c,proposal["u"])
    end
    prepared=schedule_case(c;initial=x,proposal_time_limit=3.0,nlp_time_limit=5.0,max_refinements=0)
    @test prepared["accepted"]
    @test !prepared["global_optimality_proven"]
    @test all(p->haskey(p,"tunnel_opening"),prepared["proposals"])
    if prepared["accepted"]
        actual=case_from_dict(prepared["case"])
        @test validate(actual,prepared["solution"])["valid"]
        @test replay_audit(actual,prepared["solution"])["valid"]
        @test haskey(prepared["solution"],"tunnel_opening")
    end
end
