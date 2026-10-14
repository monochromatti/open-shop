include(joinpath(@__DIR__,"..","examples","river_splits.jl"))
using .RiverSplitsExample: river_split_case

# Independent rate integration: delay each path in continuous time, then
# integrate its piecewise constant rate at its actual change times.
function split_integral(rate,knots,a,b)
    cuts=sort!(unique([a;filter(t->a<t<b,knots);b]))
    sum(rate((cuts[k]+cuts[k+1])/2)*(cuts[k+1]-cuts[k])
        for k in 1:length(cuts)-1)
end

@testset "Several tributaries and consecutive forks" begin
    grid=[0.0,0.7,1.8,3.0]
    receiver=Reservoir(name=:Receiver,z0=0.0,slope=0.1,v0=5.0,
        vmin=1.0,vmax=10.0,water_value=0.0)
    rivers=[River(;name,source,target,allocation,inflow,deterministic_delay=delay,
        capacity=10.0,curves=RiverRouting.DelayCurve[],water_value=0.0,
        history_grid=[-3.0,0.0],history_release=[0.0])
        for (name,source,target,allocation,inflow,delay) in (
            (:A1,:auto,:X,1.0,2.0,0.2),(:A2,:auto,:X,1.0,1.0,0.4),
            (:B1,:X,:Y,0.3,0.5,0.1),(:B2,:X,:Receiver,0.7,0.1,0.3),
            (:C1,:Y,:Receiver,0.2,0.2,0.25),(:C2,:Y,:Receiver,0.8,0.0,0.6))]
    operations=[OperationalSeries(object=name,attribute=:inflow,
        times=grid[1:end-1],values=values)
        for (name,values) in ((:A1,[2.0,0.0,4.0]),(:A2,[1.0,3.0,1.0]))]
    c=ScheduleCase(name="consecutive forks",system=HydroSystem(reservoirs=[receiver],rivers=rivers),
        grid=grid,prices=zeros(3),operations=operations)
    a1(t)=t<0 ? 0.0 : [2.0,0.0,4.0][searchsortedlast(grid,t)]
    a2(t)=t<0 ? 0.0 : [1.0,3.0,1.0][searchsortedlast(grid,t)]
    b1(t)=t<0 ? 0.0 : 0.3*(a1(t-0.2)+a2(t-0.4))+0.5
    b2(t)=t<0 ? 0.0 : 0.7*(a1(t-0.2)+a2(t-0.4))+0.1
    c1(t)=t<0 ? 0.0 : 0.2*b1(t-0.1)+0.2
    c2(t)=t<0 ? 0.0 : 0.8*b1(t-0.1)
    rates=(a1,a2,b1,b2,c1,c2);delays=(0.2,0.4,0.1,0.3,0.25,0.6)
    knots=sort!(unique([t+d for t in grid for d in
        (0.0,0.1,0.2,0.3,0.4,0.5,0.25,0.35,0.55,0.75,0.6,0.7,0.9,1.1)]))
    expected=[0.0036*split_integral(rates[d],knots,grid[t]-delays[d],grid[t+1]-delays[d])
        for d in 1:6,t in 1:3]
    for junctions in (RiverJunction[],[RiverJunction(name=:X),RiverJunction(name=:Y)])
        equivalent=OpenSHOP._river_replace(c;system=HydroSystem(reservoirs=[receiver],rivers=rivers,
            river_junctions=junctions))
        x=dispatch_from_controls(equivalent,zeros(Int,0,3),zeros(0,3),zeros(6,3))
        @test x.arrival_volume≈expected atol=1e-10
        @test validate(equivalent,x)["valid"]
        @test replay_audit(equivalent,x)["valid"]
        @test case_dict(case_from_dict(case_dict(equivalent)))==case_dict(equivalent)
        @test case_dict(equivalent)["river_junctions"]==[Dict("name"=>j.name) for j in junctions]
        continued=restart_case(equivalent,x,0.7)
        remainder=dispatch_from_controls(continued,zeros(Int,0,2),zeros(0,2),zeros(6,2))
        @test remainder.storage≈x.storage[:,2:end] atol=1e-10
        @test remainder.arrival_volume≈x.arrival_volume[:,2:end] atol=1e-10
        @test remainder.terminal_transit≈x.terminal_transit atol=1e-10
    end
end

@testset "Small prescribed shares retain terminal water value" begin
    rivers=[River(;name,source,target,allocation,inflow,deterministic_delay=delay,
        capacity=1.0,curves=RiverRouting.DelayCurve[],water_value=value,
        history_grid=[-2.0,0.0],history_release=[0.0])
        for (name,source,target,allocation,inflow,delay,value) in (
            (:Incoming,:Lake,:Fork,1.0,0.0,0.0,0.0),
            (:Small,:Fork,:Sea,1e-13,0.0,2.0,1e18),
            (:Main,:Fork,:Sea,1-1e-13,0.0,0.0,0.0))]
    rivers[1]=OpenSHOP._river_replace(rivers[1];law=:controlled)
    lake=Reservoir(name=:Lake,z0=10.0,slope=1.0,v0=1.0,vmin=0.0,vmax=2.0,water_value=100.0)
    c=ScheduleCase(name="small share",system=HydroSystem(reservoirs=[lake],rivers=rivers,
        boundaries=[Boundary(name=:Sea,head=0.0)]),grid=[0.0,1.0],prices=[0.0])
    x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),reshape([0.25,0.0,0.0],3,1))
    @test x.terminal_transit[2]≈0.9e-16 rtol=1e-12
    # Releasing 1 m³/s earns 360 in transit and costs .36 in reservoir value.
    # Dropping the small share would reverse the optimal decision to no release.
    result=solve_case(c;u=zeros(Int,0,1),warm=x,time_limit=30.0)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.river_release[1,1]≈1.0 atol=1e-4
        @test result.solution.objective≈359.64 atol=1e-3
    end
end

@testset "Prescribed fork, branch-local sources and delayed reunion" begin
    c=river_split_case()
    q=reshape([4.0,0.0,3.0,2.0],1,:);u=Int.(q.>0);gates=zeros(4,4)
    x=dispatch_from_controls(c,u,q,gates)
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    # The canal supplies 1 m³/s to Left. Both local inputs enter after the fork.
    upstream(t)=t<0 ? 3.5 : q[1,searchsortedlast(c.grid,t)]
    left(t)=t<0 ? 3.4 : 0.4*upstream(t-0.25)+2.0
    right(t)=t<0 ? 2.1 : 0.6*upstream(t-0.25)
    downstream(t)=t<0 ? 6.0 : left(t-0.35)+right(t-0.8)+0.5
    rates=(upstream,left,right,downstream);delays=(0.25,0.35,0.8,0.2)
    knots=sort!(unique([t+d for t in c.grid for d in
        (0.0,0.25,0.35,0.8,0.2,0.6,1.05,0.55,1.0,0.85,1.25)]))
    arrivals=[0.0036*split_integral(rates[d],knots,
        c.grid[t]-delays[d],c.grid[t+1]-delays[d]) for d in 1:4,t in 1:4]
    transit=[0.0036*split_integral(rates[d],knots,t-delays[d],t)
        for d in 1:4,t in c.grid]
    @test x.arrival_volume≈arrivals atol=1e-10
    trajectory=simulate(c,q,gates)
    @test trajectory.transit≈transit atol=1e-10
    @test x.tunnel_q≈ones(1,4) atol=1e-10
    @test vec(x.storage[1,:])≈[5.0;5.0 .-0.0036 .*cumsum(diff(c.grid).*vec(q))] atol=1e-10
    @test vec(x.storage[2,:])≈[5.0;5.0 .+cumsum(arrivals[4,:])] atol=1e-10
    total=vec(sum(x.storage;dims=1)+sum(trajectory.transit;dims=1))
    @test total≈total[1] .+0.0036*2.5 .*c.grid atol=1e-10

    # Water before and after the fork is retained on restart, without reallocation.
    for time in (0.5,1.5)
        edge=findfirst(==(time),c.grid)
        continued=restart_case(c,x,time)
        remainder=dispatch_from_controls(continued,u[:,edge:end],q[:,edge:end],gates[:,edge:end])
        @test validate(continued,remainder)["valid"]
        @test replay_audit(continued,remainder)["valid"]
        @test remainder.storage≈x.storage[:,edge:end] atol=1e-10
        @test simulate(continued,q[:,edge:end],gates[:,edge:end]).transit≈trajectory.transit[:,edge:end] atol=1e-10
        @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)
    end
    public=case_dict(c)
    @test public["rivers"][2]["source"]==:Fork
    @test public["rivers"][4]["source"]==:Merge
    @test public["rivers"][2]["allocation"]==0.4
    @test case_dict(ScheduleCase(name=c.name,system=c.system,grid=c.grid,prices=c.prices))==public
    mktempdir() do dir
        path=joinpath(dir,"fork.json")
        writejson(path,case_dict(c))
        @test case_dict(readcase(path))==public
    end

    # A wrong split can conserve total water. Check the prescribed branch law too.
    wrong=deepcopy(x);wrong.river_release[2,2]+=0.1;wrong.river_release[3,2]-=0.1
    @test !validate(c,wrong)["valid"]
    @test validate(c,wrong)["residuals"]["river_allocation_Left"]>1e-4
end

@testset "Fork capacity limits optimization" begin
    c=river_split_case();q=ones(1,4);u=ones(Int,1,4);gates=zeros(4,4)
    initial=dispatch_from_controls(c,u,q,gates)
    # Left's 3.6 m³/s capacity = 2 local + 40% of at most 4 upstream.
    best=dispatch_from_controls(c,u,fill(4.0,1,4),gates)
    # Receiver's arrival stream is 6 until .8h, 6.2 until 1.25h, then 6.5.
    # Upper loses 4 m³/s. Their midpoint storage levels determine net head.
    receiver=[5.0;5.0 .+0.0036 .*cumsum([3.0,6.215,6.5,9.75])]
    upper=5.0 .-0.0036*4.0 .*c.grid
    nethead=50.0 .+0.1 .*(upper[1:end-1]+upper[2:end]-receiver[1:end-1]-receiver[2:end])
    expected=sum(c.prices.*diff(c.grid).*nethead)*0.00981*0.9*4.0
    @test validate(c,best)["valid"]
    @test best.objective≈expected atol=1e-9
    @test best.objective>initial.objective
    # Solver-scale roundoff is admitted without altering controls or water.
    rounded=fill(4.0+5e-7,1,4)
    near=dispatch_from_controls(c,u,rounded,gates)
    @test near.generator_q==rounded
    @test maximum(near.river_release[2,:])>3.6
    @test 0<validate(c,near)["residuals"]["river_capacity_Left"]<1e-6
    @test replay_audit(c,near)["valid"]
    @test_throws ArgumentError dispatch_from_controls(c,u,fill(5.0,1,4),gates)
    local_result=solve_case(c;u,warm=initial,time_limit=60.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    if haskey(local_result,"objective")
        @test local_result["objective"]≈expected atol=1e-3
        @test replay_audit(c,local_result)["valid"]
    end
    proposal=OpenSHOP.propose_commitment(c;reference=initial,time_limit=10.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        @test admissible(c,proposal["u"])
    end
    result=solve(c;initial,time_limit=90.0,relative_gap=1e-4)
    @test result.accepted
    if result["global_bound"]!==nothing
        @test isfinite(result["global_bound"]) && result["global_bound"]>=expected-1e-6
    end
    if result.solution!==nothing
        @test result.solution.objective≈expected atol=1e-3
        @test result.solution.generator_q≈fill(4.0,1,4) atol=2e-4
    end
end

@testset "Distributed fork responds to total branch flow" begin
    c=river_split_case();r=c.system.rivers
    curves=[RiverRouting.DelayCurve(flow,[delay,delay+0.5],[1.0])
        for (flow,delay) in ((0.0,0.0),(2.0,0.5),(3.6,1.0))]
    rivers=[OpenSHOP._river_replace(reach;history_release=[0.0],
        deterministic_delay=d==2 ? nothing : 0.0,
        curves=d==2 ? curves : RiverRouting.DelayCurve[])
        for (d,reach) in enumerate(r)]
    c=ScheduleCase(name="distributed fork",system=OpenSHOP._river_replace(c.system;rivers),
        grid=[0.0,1.0,2.0,3.0],prices=[40.0,90.0,55.0])
    q=reshape([2.0,0.0,0.0],1,:);u=Int.(q.>0);gates=zeros(4,3)
    x=dispatch_from_controls(c,u,q,gates)
    @test validate(c,x)["valid"]
    # 2 local + .4*2 upstream = 2.8. Half of its delay mass uses [.5,1]
    # and half uses [1,1.5]. Only 1/8 arrives in the first hour.
    @test x.river_release[2,:]≈[2.8,2.0,2.0] atol=1e-10
    @test x.arrival_volume[2,1]≈0.0036*2.8/8 atol=1e-10
    routed=route_network(c,zeros(4,3);generator_q=q,tunnel_q=ones(1,3))
    @test routed["arrival_volume"]≈x.arrival_volume atol=1e-10
    trajectory=simulate(c,q,gates)
    total=vec(sum(trajectory.storage;dims=1)+sum(trajectory.transit;dims=1))
    @test total≈total[1] .+0.0036*2.5 .*c.grid atol=1e-10
    continued=restart_case(c,x,1.0)
    remainder=dispatch_from_controls(continued,u[:,2:end],q[:,2:end],gates[:,2:end])
    @test remainder.storage≈x.storage[:,2:end] atol=1e-10
    @test remainder.terminal_transit≈x.terminal_transit atol=1e-10
    @test validate(continued,remainder)["valid"]
    local_result=solve_case(c;u,warm=x,feasibility_only=true,time_limit=20.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
end

@testset "Zero allocation leaves branch-local water and history intact" begin
    c=river_split_case();r=c.system.rivers
    rivers=[r[1],OpenSHOP._river_replace(r[2];allocation=0.0),
        OpenSHOP._river_replace(r[3];allocation=1.0),r[4]]
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;rivers))
    x=dispatch_from_controls(c,ones(Int,1,4),ones(1,4),zeros(4,4))
    @test validate(c,x)["valid"]
    @test x.river_release[2,:]≈fill(2.0,4) atol=1e-10
    @test x.arrival_volume[2,1]≈0.0036*(0.35*3.4+0.15*2.0) atol=1e-10
    @test replay_audit(c,x)["valid"]
end

@testset "Invalid allocations and stale transport are rejected" begin
    c=river_split_case();r=c.system.rivers
    altered(a,b)=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        rivers=[r[1],OpenSHOP._river_replace(r[2];allocation=a),
            OpenSHOP._river_replace(r[3];allocation=b),r[4]]))
    for (a,b) in ((-0.1,1.1),(NaN,0.6),(Inf,0.6),(0.4,0.5),(0.4,0.7))
        bad=altered(a,b)
        @test_throws Exception OpenSHOP.validate_inputs(bad)
        @test_throws ArgumentError simulate(bad,ones(1,4),zeros(4,4))
        @test_throws ArgumentError forward_step(bad.system,[5.0,5.0],1.0,[1.0],zeros(4),zeros(4))
    end
    changed=altered(0.5,0.5)
    data=deterministic_network_data(c)
    @test_throws ArgumentError route_network_controlled(changed,zeros(4,4);
        generator_q=ones(1,4),tunnel_q=ones(1,4),transport=data)
    # A physical outlet is controlled by its own law, not a fork allocation.
    outlet=River(name=:Outlet,source=:Upper,target=:Receiver,allocation=0.5,
        law=:controlled,capacity=10.0,curves=RiverRouting.DelayCurve[],
        deterministic_delay=0.0,water_value=0.0)
    bad=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        rivers=[r;outlet]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(bad)
    public=case_dict(c)
    public["rivers"][1]["target"]=:Left
    @test_throws ArgumentError case_from_dict(public)
    collision=case_dict(c)
    collision["rivers"][1]["target"]=Symbol("__river_node__Fork")
    for d in (2,3)
        collision["rivers"][d]["source"]=Symbol("__river_node__Fork")
    end
    @test_throws ArgumentError case_from_dict(collision)
end
