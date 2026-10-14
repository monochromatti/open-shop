include(joinpath(@__DIR__,"..","examples","directional_interlocks.jl"))
using .DirectionalInterlocksExample: directional_interlock_case

di_copy(x;kwargs...)=typeof(x)(;merge((; (k=>getfield(x,k) for k in fieldnames(typeof(x)))...),
    (;kwargs...))...)

# Fixed-head boundaries isolate chronology tests from tiny-step hydraulic
# conditioning. The explicitly protected station does not depend on Flag flow.
function di_clock_case(grid)
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,pmin=0.1,pmax=5.0,
        hmin=40.0,hmax=110.0,efficiency=0.9,qbest=5.0,hbest=100.0,
        qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0)
    rule=TunnelInterlock(name=:ClosureWait,tunnels=[:Flag],forward_tunnels=[:Flag],
        plants=[:Station],forward_switch_time=1.0)
    system=HydroSystem(boundaries=[Boundary(name=:Supply,head=100.0),Boundary(name=:Tail,head=0.0)],
        tunnels=[Tunnel(name=:Flag,source=:Supply,target=:Tail,resistance=100.0,
            capacity=2.0,opening=1.0,initial_opening=1.0,gate_positions=[0.0,1.0])],
        plants=[Plant(name=:Station,source=:Supply,target=:Tail,pmax=5.0)],
        generators=[unit],tunnel_interlocks=[rule])
    ScheduleCase(name="closure clock",system=system,grid=grid,prices=fill(40.0,length(grid)-1))
end

# Independent midpoint head and storage calculation for a scheduled single
# intake. Positive prices make q=5 optimal whenever the station may run:
# throughout q∈[1,5], q*(reservoir head - q²) is strictly increasing here.
function di_reference(c,u,A;discharge=5.0)
    V=repeat([5.0,5.0],1,length(c.grid));p=zeros(size(u))
    previous=[1.0,0.0];cost=0.0
    for t in eachindex(c.prices)
        active=only(findall(>(0),A[:,t]))
        q=discharge*u[1,t]
        V[:,t+1].=V[:,t]
        V[active,t+1]-=0.0036*diff(c.grid)[t]*q
        source=(active==1 ? 99.5 : 79.5)+0.05*(V[active,t]+V[active,t+1])
        p[1,t]=0.00981*0.9*q*(source-q^2/A[active,t])
        cost+=5.0*sum(abs.(A[:,t]-previous))
        previous=A[:,t]
    end
    (;storage=V,power=p,objective=sum(c.prices.*diff(c.grid).*vec(p))-cost)
end

@testset "Endogenous gate closures trade production against waiting" begin
    base=directional_interlock_case()
    unit=di_copy(only(base.system.generators);qmax=3.0,qbest=3.0)
    operations=[OperationalSeries(object=name,attribute=:opening,times=[3.0],values=[value])
        for (name,value) in ((:HighIntake,0.0),(:LowIntake,1.0))]
    c=di_copy(base;grid=collect(0.0:6.0),prices=[20.0,1.0,1.0,100.0,60.0,50.0],
        operations,system=di_copy(base.system;generators=[unit]))
    # Each early interval may choose either intake at half or full opening.
    # Gate changes after hour 3 are prescribed. For these positive prices,
    # the physical power derivatives are positive through q=3, even at half
    # opening, and future head depletion is smaller than the current gain.
    # Thus each unblocked interval runs at q=3; enumerate gate decisions only.
    choices=((1.0,0.0),(0.5,0.0),(0.0,1.0),(0.0,0.5))
    candidates=map(Iterators.product(choices,choices,choices)) do early
        A=hcat(collect.(early)...,fill([0.0,1.0],3)...)
        u=zeros(Int,1,6);previous=(1.0,0.0);until=0.0
        for t in 1:6
            for (i,wait) in ((1,1.5),(2,0.5))
                previous[i]>0 && A[i,t]==0 && (until=max(until,c.grid[t]+wait))
            end
            u[1,t]=Int(c.grid[t]>=until)
            previous=(A[1,t],A[2,t])
        end
        (;opening=A,commitment=u,reference=di_reference(c,u,A;discharge=3.0))
    end
    best=argmax(candidate.reference.objective for candidate in candidates)
    winner=collect(candidates)[best]
    witness=dispatch_from_controls(c,winner.commitment,3.0 .*winner.commitment,zeros(0,6);
        tunnel_opening=winner.opening)
    @test validate(c,witness)["valid"]
    @test witness.objective≈winner.reference.objective atol=1e-8
    @test replay_audit(c,witness)["valid"]
    A=tunnel_opening_schedule(c);initial_u=reshape([1,1,1,0,0,1],1,:)
    initial=dispatch_from_controls(c,initial_u,1.0 .*initial_u,zeros(0,6);tunnel_opening=A)
    result=solve(c;initial,time_limit=90.0,relative_gap=1e-4)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.objective≈winner.reference.objective atol=1e-3
        @test result.solution.objective>initial.objective+100.0
        @test validate_tunnel_openings(c,result.solution.tunnel_opening,result.solution.commitment)
        @test replay_audit(c,result.solution)["valid"]
    end
end

@testset "Asymmetric plant-stop windows use physical closure times" begin
    c=directional_interlock_case();A=tunnel_opening_schedule(c)
    u=reshape([1,1,0,0,1,0,1,1],1,:);q=5.0 .*u;gate=zeros(0,8)
    @test OpenSHOP.validate_inputs(c)
    @test validate_tunnel_openings(c,A,u)
    x=dispatch_from_controls(c,u,q,gate;tunnel_opening=A)
    expected=di_reference(c,u,A)
    @test validate(c,x)["valid"]
    @test x.storage≈expected.storage atol=1e-9
    @test x.power≈expected.power atol=1e-8
    @test x.objective≈expected.objective atol=1e-8
    @test replay_audit(c,x)["valid"]
    @test simulate(c,q,gate;grid=OpenSHOP.refined_grid(c.grid;factor=5),tunnel_opening=A).converged
    # Low opens at 1h, even though production must wait until 2.5h. High opens
    # at 3h, and production may resume at 3.5h, with one gate always open.
    for t in (3,4,6)
        early=copy(u);early[1,t]=1
        @test_throws ArgumentError validate_tunnel_openings(c,A,early)
        @test_throws ArgumentError dispatch_from_controls(c,early,5.0 .*early,gate;tunnel_opening=A)
        @test_throws ArgumentError simulate(c,5.0 .*early,gate;tunnel_opening=A)
    end
    # A claimed on-state with zero discharge still violates the operational rule.
    false_state=copy(u);false_state[1,3]=1
    @test_throws ArgumentError dispatch_from_controls(c,false_state,q,gate;tunnel_opening=A)
    bad=schedule_dict(x);bad["u"][1,4]=1
    @test !validate(c,bad)["valid"]
    @test !replay_audit(c,bad)["valid"]
    @test_throws ArgumentError solve_case(c;u=false_state,tunnel_opening=A,time_limit=5.0)
    # Moving from one positive opening to another is not a closure.
    positive=[1.0 .5 1.0 1.0 1.0 1.0 1.0 1.0;zeros(8)']
    free=di_copy(c;operations=OperationalSeries[])
    @test validate_tunnel_openings(free,positive,ones(Int,1,8))
    # The clock also starts if the closing tunnel carried no water beforehand.
    no_prior_flow=zeros(Int,1,8);no_prior_flow[1,4]=1
    @test_throws ArgumentError validate_tunnel_openings(c,A,no_prior_flow)
    reopening=copy(A);reopening[:,4]=[1.0,0.0]
    @test_throws ArgumentError validate_tunnel_openings(free,reopening,no_prior_flow)
end

@testset "Optimization enforces closure waits in proposals, NLP and SCIP" begin
    c=directional_interlock_case();A=tunnel_opening_schedule(c)
    u=reshape([1,1,0,0,1,0,1,1],1,:);expected=di_reference(c,u,A)
    initial=dispatch_from_controls(c,u,1.0 .*u,zeros(0,8);tunnel_opening=A)
    local_result=solve_case(c;u,warm=initial,tunnel_opening=A,time_limit=45.0)
    @test local_result.accepted
    if local_result.solution!==nothing
        @test local_result.solution.objective≈expected.objective atol=1e-3
        @test replay_audit(c,local_result.solution)["valid"]
    end
    proposal=OpenSHOP.propose_commitment(c;reference=initial,time_limit=10.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        @test validate_tunnel_openings(c,proposal["tunnel_opening"],proposal["u"])
        @test proposal["u"]==u
    end
    result=solve(c;initial,time_limit=90.0,relative_gap=1e-4)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.commitment==u
        @test result.solution.generator_q≈5.0 .*u atol=2e-4
        @test result.solution.objective≈expected.objective atol=1e-3
        @test replay_audit(c,result.solution)["valid"]
    end
    if result.global_bound!==nothing
        @test isfinite(result.global_bound) && result.global_bound>=expected.objective-1e-3
    end
    prepared=schedule_case(c;proposal_time_limit=5.0,nlp_time_limit=15.0,
        max_refinements=0,time_limit=60.0)
    @test prepared["accepted"]
    if prepared["solution"]!==nothing
        @test validate_tunnel_openings(c,prepared["solution"].tunnel_opening,prepared["solution"].commitment)
        @test prepared["solution"].objective≈expected.objective atol=1e-3
    end
end

@testset "Restart retains both directions and a closure at the restart edge" begin
    c=directional_interlock_case();A=tunnel_opening_schedule(c)
    u=reshape([1,1,0,0,1,0,1,1],1,:);q=5.0 .*u
    x=dispatch_from_controls(c,u,q,zeros(0,8);tunnel_opening=A)
    for (time,forward,backward) in ((1.0,0.0,0.0),(2.0,0.5,0.0),(2.5,0.0,0.0),
        (3.0,0.0,0.0),(3.5,0.0,0.0))
        edge=findfirst(==(time),c.grid)
        continued=restart_case(c,x,time);rule=only(continued.system.tunnel_interlocks)
        @test rule.initial_forward_wait≈forward atol=1e-12
        @test rule.initial_backward_wait≈backward atol=1e-12
        rest=dispatch_from_controls(continued,u[:,edge:end],q[:,edge:end],zeros(0,9-edge);
            tunnel_opening=A[:,edge:end])
        @test validate(continued,rest)["valid"]
        @test replay_audit(continued,rest)["valid"]
        @test rest.storage≈x.storage[:,edge:end] atol=1e-9
        @test rest.power≈x.power[:,edge:end] atol=1e-8
        @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
        if time==2.0
            early=copy(u[:,edge:end]);early[1,1]=1
            @test_throws ArgumentError validate_tunnel_openings(continued,A[:,edge:end],early)
        end
    end
    # There is also a grid edge strictly inside the shorter backward wait.
    fine=with_grid(c,sort([c.grid;3.25]))
    indices=[searchsortedlast(c.grid,fine.grid[t]) for t in eachindex(fine.prices)]
    xf=dispatch_from_controls(fine,u[:,indices],q[:,indices],zeros(0,length(indices));
        tunnel_opening=A[:,indices])
    continued=restart_case(fine,xf,3.25)
    @test only(continued.system.tunnel_interlocks).initial_backward_wait≈0.25 atol=1e-12
    edge=findfirst(==(3.25),fine.grid)
    rest=dispatch_from_controls(continued,xf.commitment[:,edge:end],xf.generator_q[:,edge:end],
        zeros(0,length(continued.prices));tunnel_opening=xf.tunnel_opening[:,edge:end])
    @test rest.storage≈xf.storage[:,edge:end] atol=1e-9
    @test replay_audit(continued,rest)["valid"]
end

@testset "Terminal obligations and repeated closures survive continuation" begin
    base=directional_interlock_case()
    ops=[OperationalSeries(object=name,attribute=:opening,times=[0.0,1.0,3.0,4.5],values=values)
        for (name,values) in ((:HighIntake,[1.0,0.0,1.0,0.0]),(:LowIntake,[0.0,1.0,0.0,1.0]))]
    c=di_copy(base;grid=[base.grid;6.0;6.5],prices=[base.prices;40.0;40.0],operations=ops)
    A=tunnel_opening_schedule(c);u=reshape([1,1,0,0,1,0,1,0,0,1],1,:);q=5.0 .*u
    x=dispatch_from_controls(c,u,q,zeros(0,10);tunnel_opening=A)
    @test validate(c,x)["valid"]
    continued=restart_case(c,x,5.0)
    @test only(continued.system.tunnel_interlocks).initial_forward_wait≈1.0 atol=1e-12
    rest=dispatch_from_controls(continued,u[:,9:end],q[:,9:end],zeros(0,2);tunnel_opening=A[:,9:end])
    @test rest.storage≈x.storage[:,9:end] atol=1e-9
    @test replay_audit(continued,rest)["valid"]
    prefix=di_copy(c;grid=c.grid[1:9],prices=c.prices[1:8])
    past=dispatch_from_controls(prefix,u[:,1:8],q[:,1:8],zeros(0,8);tunnel_opening=A[:,1:8])
    obligations=validate(prefix,past)["terminal_tunnel_switch_obligations"]
    @test any(state.interlock==:SupplyHandover && state.forward_wait≈1.0 for state in obligations)
    @test x.objective≈past.objective+rest.objective atol=1e-8
end

@testset "A wait expired within time tolerance stays expired on restart" begin
    c=di_clock_case([0.0,1.0-6e-9,1.0-5e-9,2.0])
    A=zeros(1,3);u=reshape([0,1,1],1,:);q=5.0 .*u
    x=dispatch_from_controls(c,u,q,zeros(0,3);tunnel_opening=A)
    @test validate(c,x)["valid"]
    continued=restart_case(c,x,c.grid[3])
    @test only(continued.system.tunnel_interlocks).initial_forward_wait==0.0
    rest=dispatch_from_controls(continued,u[:,3:end],q[:,3:end],zeros(0,1);tunnel_opening=A[:,3:end])
    @test validate(continued,rest)["valid"]
    @test rest.storage≈x.storage[:,3:end] atol=1e-9
end

@testset "Control normalization preserves open and closed identity" begin
    c=di_clock_case([0.0,0.5,1.0,2.0])
    tunnel=di_copy(only(c.system.tunnels);opening=0.0,initial_opening=0.0)
    c=di_copy(c;system=di_copy(c.system;tunnels=[tunnel]))
    tiny=reshape([big"1e-500",BigFloat(0),BigFloat(0)],1,:)
    @test_throws ArgumentError dispatch_from_controls(c,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);
        tunnel_opening=tiny)
    @test_throws ArgumentError simulate(c,fill(5.0,1,3),zeros(0,3);tunnel_opening=tiny)
    @test_throws ArgumentError validate_tunnel_openings(c,tiny,ones(Int,1,3))
    above=fill(nextfloat(BigFloat(1)),1,3)
    @test_throws ArgumentError dispatch_from_controls(c,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);
        tunnel_opening=above)
    closed=zeros(BigFloat,1,3)
    x=dispatch_from_controls(c,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);tunnel_opening=closed)
    @test validate(c,x)["valid"]
end

@testset "Initial history and directional references are explicit" begin
    c=directional_interlock_case();rule=only(c.system.tunnel_interlocks)
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    mktempdir() do dir
        path=joinpath(dir,"handover.json");writejson(path,case_dict(c))
        @test case_dict(readcase(path))==case_dict(c)
    end
    for bad in (di_copy(rule;forward_switch_time=-1.0),di_copy(rule;initial_forward_wait=2.0),
        di_copy(rule;forward_tunnels=[:Unknown]),di_copy(rule;backward_tunnels=[:HighIntake]),
        di_copy(rule;forward_tunnels=[:HighIntake,:HighIntake]),di_copy(rule;plants=[:Unknown]),
        di_copy(rule;plants=Symbol[]),di_copy(rule;forward_tunnels=Symbol[]),
        di_copy(rule;initial_forward_wait=1.0))
        badcase=di_copy(c;system=di_copy(c.system;tunnel_interlocks=[bad]))
        @test_throws ArgumentError OpenSHOP.validate_inputs(badcase)
    end
    # A pre-horizon closure .5h ago leaves 1h of the forward wait. The unit
    # has already been stopped for .5h and must stay off until hour 11.
    initial=di_copy(rule;initial_forward_wait=1.0)
    off=di_copy(only(c.system.generators);initial_on=0,initial_age=0.5)
    start=di_copy(c;grid=[10.0,10.5,11.0,12.0],prices=fill(40.0,3),
        operations=OperationalSeries[],system=di_copy(c.system;generators=[off],tunnel_interlocks=[initial]))
    @test OpenSHOP.validate_inputs(start)
    A=ones(2,3);A[2,:].=0.0
    @test validate_tunnel_openings(start,A,reshape([0,0,1],1,:))
    @test_throws ArgumentError validate_tunnel_openings(start,A,ones(Int,1,3))
    recent=di_copy(start;system=di_copy(start.system;generators=[di_copy(off;initial_age=0.25)]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(recent)
    # Timer constraints apply to every member of the protected plant.
    member=di_copy(off;name=:OtherUnit)
    multi=di_copy(start;system=di_copy(start.system;generators=[off,member]))
    @test_throws ArgumentError validate_tunnel_openings(multi,A,[0 0 1;1 0 0])
end
