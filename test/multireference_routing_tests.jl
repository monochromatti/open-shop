function multireference_curves()
    edges=[0.1,0.6,1.3,2.4]
    [RiverRouting.DelayCurve(q,edges,w) for (q,w) in (
        (0.0,[0.0,0.2,0.8]),
        (8.0,[0.15,0.55,0.3]),
        (20.0,[0.45,0.45,0.1]),
        (50.0,[0.75,0.2,0.05]),
    )]
end

# Independent neighboring-reference interpolation, including knot values.
function multireference_expected(curves,values,q)
    length(curves)==1 && return q*only(values)
    j=findlast(c->c.reference_flow<=q,curves)
    j==length(curves) && return q*last(values)
    left,right=curves[j].reference_flow,curves[j+1].reference_flow
    q*((right-q)*values[j]+(q-left)*values[j+1])/(right-left)
end

@testset "Distributed curves accept fixed or multiple reference flows" begin
    curves=multireference_curves()
    @test isnothing(RiverRouting.check_curves(curves;capacity=40.0))
    @test isnothing(RiverRouting.check_curves(curves;capacity=50.0))
    fixed=[RiverRouting.DelayCurve(12.0,[0.0,0.8],[1.0])]
    @test isnothing(RiverRouting.check_curves(fixed;capacity=100.0))
    @test_throws ArgumentError RiverRouting.check_curves(curves;capacity=51.0)
    @test_throws ArgumentError RiverRouting.check_curves(curves[2:end];capacity=40.0)
    @test_throws ArgumentError RiverRouting.check_curves(curves;capacity=Inf)
    @test_throws ArgumentError RiverRouting.check_curves(curves;capacity=-1.0)
    @test_throws ArgumentError RiverRouting.check_curves(reverse(curves))
    @test_throws ArgumentError RiverRouting.check_curves([curves[1],curves[1]])
    @test_throws ArgumentError RiverRouting.check_curves(RiverRouting.DelayCurve[])
    @test_throws ArgumentError RiverRouting.DelayCurve(0.0,[0.0,1.0,2.0],[-0.1,1.1])
    @test_throws ArgumentError RiverRouting.DelayCurve(0.0,[0.0,1.0,2.0],[0.2,0.3])
    @test_throws ArgumentError RiverRouting.blend(fixed,NaN)
    @test_throws ArgumentError RiverRouting.blend(fixed,-1.0)
    @test_throws ArgumentError RiverRouting.blend(curves,51.0)
    @test RiverRouting.transfer_value(fixed,[0.25],40.0)==10.0
    for q in (0.0,4.0,8.0,12.0,20.0,35.0,50.0)
        weights=RiverRouting.blend(curves,q)
        @test sum(last,weights)≈1.0
        @test all(pair->pair[2]>=0,weights)
        @test length(weights)<=2
        @test length(weights)==1 || weights[2][1]==weights[1][1]+1
        for coefficients in ([0.1,0.8,0.3,0.5],[1.0,0.0,1.0,0.0])
            @test RiverRouting.transfer_value(curves,coefficients,q)≈
                multireference_expected(curves,coefficients,q) atol=1e-12
        end
    end
end

@testset "Multiple-reference cohorts conserve volume and retain restart history" begin
    curves=multireference_curves()
    grid=[-2.0,-1.7,-0.6,0.0,0.3,1.2,2.5,3.4]
    q=[0.0,4.0,8.0,12.0,20.0,35.0,50.0]
    arrivals=[-5.0,-2.0,-0.5,0.0,0.4,1.0,1.8,3.4,4.0,6.0]
    B=RiverRouting.coefficients(arrivals,grid,curves)
    expected=[sum(0.0036*(grid[k+1]-grid[k])*
        multireference_expected(curves,view(B,t,k,:),q[k]) for k in eachindex(q))
        for t in 1:(length(arrivals)-1)]
    routed=RiverRouting.route(arrivals,grid,q,curves)
    @test routed≈expected atol=1e-12
    @test minimum(routed)>=0
    @test sum(routed)≈0.0036*sum(diff(grid).*q) atol=1e-12
    for time in (-1.3,0.0,0.6,2.2,4.0,6.0)
        released=0.0036*sum(q[k]*max(0.0,min(time,grid[k+1])-grid[k])
            for k in eachindex(q))
        arrived=only(RiverRouting.route([-5.0,time],grid,q,curves))
        remaining=RiverRouting.remaining_volume(time,grid,q,curves)
        @test remaining>=0
        @test arrived+remaining≈released atol=1e-12
    end
    time=0.8
    history=release_history(grid,q,time)
    k=searchsortedlast(grid,time)
    future_grid=vcat(time,grid[(k+1):end])
    future_release=q[k:end]
    later_grid=[time,1.0,1.8,3.4,4.0,6.0]
    @test RiverRouting.route(later_grid,grid,q,curves)≈
        RiverRouting.route(later_grid,history.grid,history.release,curves)+
        RiverRouting.route(later_grid,future_grid,future_release,curves) atol=1e-12
    for t in later_grid
        @test RiverRouting.remaining_volume(t,grid,q,curves)≈
            RiverRouting.remaining_volume(t,history.grid,history.release,curves)+
            RiverRouting.remaining_volume(t,future_grid,future_release,curves) atol=1e-12
    end
    r=River(name=:distributed,source=:Lake,target=:Sea,curves=curves,
        capacity=50.0,water_value=0.0,history_grid=[-2.0,0.0],history_release=[12.0])
    @test size(OpenSHOP.point_coefficients(r,[0.0,0.7,1.4],0.9))==(2,4)
    @test point_arrival(r,[0.0,0.7,1.4],[4.0,35.0],0.9)≈
        sum(multireference_expected(curves,view(OpenSHOP.point_coefficients(r,
            [0.0,0.7,1.4],0.9),k,:),[4.0,35.0][k]) for k in 1:2)+
        multireference_expected(curves,view(OpenSHOP.point_coefficients(r,
            r.history_grid,0.9),1,:),12.0) atol=1e-12
end

@testset "Numerical routing streams only overlapping cohort arrivals" begin
    curves=[RiverRouting.DelayCurve(q,[0.0,0.05,0.2],weights) for (q,weights) in
        ((0.0,[0.1,0.9]),(20.0,[0.5,0.5]),(50.0,[0.9,0.1]))]
    releases=[-1.0,-0.4,0.0,0.01,0.07,0.6,1.2]
    flows=[0.0,8.0,20.0,31.0,50.0,12.0]
    # Include exact support edges, disjoint windows and partial overlaps.
    for arrival_grid in ([-2.0,-1.0,-0.95,-0.2,0.0,0.05,0.2,0.6,1.4,2.0],
        [-4.0,-3.0], [2.0,3.0], [0.03,0.04,0.08,0.2,0.9])
        B=RiverRouting.coefficients(arrival_grid,releases,curves)
        expected=[sum(0.0036*(releases[k+1]-releases[k])*
            multireference_expected(curves,view(B,t,k,:),flows[k])
            for k in eachindex(flows)) for t in 1:(length(arrival_grid)-1)]
        @test RiverRouting.route(arrival_grid,releases,flows,curves)≈expected atol=1e-12
    end
    # The dense optimization API remains available, but numerical replay must
    # not allocate its 1024×1024×3 (24 MiB) transfer tensor.
    grid=collect(range(0.0,24.0;length=1025))
    q=[(5.0,12.0,24.0,36.0)[mod1(k,4)] for k in 1:1024]
    RiverRouting.route(grid,grid,q,curves) # Warm this exact specialization.
    allocated=@allocated RiverRouting.route(grid,grid,q,curves)
    @test allocated<5_000_000
    volume=sum(RiverRouting.route(grid,grid,q,curves))
    @test volume+RiverRouting.remaining_volume(last(grid),grid,q,curves)≈
        0.0036*sum(diff(grid).*q) atol=1e-10
end

@testset "Piecewise-quadratic transfer bounds include stationary extrema" begin
    curves=multireference_curves()
    @test collect(RiverRouting.transfer_range(curves,[1.0,0.0,0.0,0.0],0.0,8.0))≈[0.0,2.0]
    @test collect(RiverRouting.transfer_range(curves,[-2.0,1.0,0.0,0.0],0.0,8.0))≈[-8/3,8.0]
    for coefficients in ([0.1,0.8,0.3,0.5],[-0.7,0.3,0.4,-0.2]),
        (a,b) in ((0.0,50.0),(3.0,39.0),(8.0,8.0),(12.0,20.0))
        lo,hi=RiverRouting.transfer_range(curves,coefficients,a,b)
        points=range(a,b;length=1001)
        values=[multireference_expected(curves,coefficients,x) for x in points]
        @test lo<=minimum(values)+1e-10
        @test hi>=maximum(values)-1e-10
        # Endpoint or stationary extrema are attained within the bounded cells.
        @test minimum(values)-lo<=1e-3
        @test hi-maximum(values)<=1e-3
    end
    @test_throws ArgumentError RiverRouting.transfer_range(curves,[1.0],0.0,50.0)
    @test_throws ArgumentError RiverRouting.transfer_range(curves,ones(4),-1.0,50.0)
    @test_throws ArgumentError RiverRouting.transfer_range(curves,ones(4),0.0,51.0)
end

@testset "Local and global expressions select the same reference pair" begin
    curves=multireference_curves()
    coefficients=[0.1,0.8,0.3,0.5]
    local_model=Model(OpenSHOP.Ipopt.Optimizer)
    set_silent(local_model)
    @variable(local_model,0<=q<=50)
    expression=OpenSHOP.distributed_transfer_expression(local_model,curves,q,coefficients)
    second_expression=OpenSHOP.distributed_transfer_expression(local_model,curves,q,
        copy(coefficients))
    @test length(local_model.ext[:distributed_transfer_operators])==1
    for flow in (0.0,4.0,8.0,12.0,20.0,35.0,50.0)
        actual=JuMP.value(v->flow,expression)
        @test actual≈multireference_expected(curves,coefficients,flow) atol=1e-11
        @test JuMP.value(v->flow,second_expression)≈actual atol=1e-12
    end
    fix(q,12.0;force=true)
    @variable(local_model,y)
    @constraint(local_model,y==expression)
    @objective(local_model,Min,y)
    optimize!(local_model)
    @test termination_status(local_model)==OpenSHOP.MOI.LOCALLY_SOLVED
    @test value(y)≈multireference_expected(curves,coefficients,12.0) atol=1e-8

    # An interior stationary optimum exercises derivatives of the registered
    # operator with a free discharge, rather than a fixed-variable evaluation.
    stationary_model=Model(OpenSHOP.Ipopt.Optimizer)
    set_silent(stationary_model)
    @variable(stationary_model,0<=stationary_q<=8,start=2)
    stationary_expression=OpenSHOP.distributed_transfer_expression(stationary_model,
        curves,stationary_q,[1.0,0.0,0.0,0.0])
    @objective(stationary_model,Max,stationary_expression)
    optimize!(stationary_model)
    @test termination_status(stationary_model)==OpenSHOP.MOI.LOCALLY_SOLVED
    @test value(stationary_q)≈4.0 atol=1e-6
    @test objective_value(stationary_model)≈2.0 atol=1e-8

    global_model=Model(OpenSHOP.SCIP.Optimizer)
    set_silent(global_model)
    @variable(global_model,0<=x<=50)
    coordinate=OpenSHOP.distributed_reference_coordinate!(global_model,curves,x;
        name=:routing_test)
    @test OpenSHOP.distributed_reference_coordinate!(global_model,curves,x;
        name=:routing_test)===coordinate
    @test_throws ArgumentError OpenSHOP.distributed_reference_coordinate!(
        global_model,curves,2x;name=:routing_test)
    global_expression=OpenSHOP.distributed_transfer_expression(global_model,curves,x,
        coefficients;coordinate)
    terminal_expression=OpenSHOP.distributed_transfer_expression(global_model,curves,x,
        1 .- coefficients;coordinate)
    @test length(global_model.ext[:global_tensor_coordinates])==1
    for flow in (0.0,4.0,8.0,12.0,20.0,35.0,50.0)
        weights=zeros(4)
        for (j,w) in RiverRouting.blend(curves,flow)
            weights[j]=w
        end
        assignment=Dict(x=>flow)
        @test OpenSHOP._lift_distributed_references!(global_model,assignment)===assignment
        @test [assignment[v] for v in coordinate.weights]≈weights atol=1e-12
        @test sum(weights)≈1.0
        @test sum(curves[j].reference_flow*weights[j] for j in 1:4)≈flow
        @test JuMP.value(v->assignment[v],global_expression)≈
            multireference_expected(curves,coefficients,flow) atol=1e-11
        @test JuMP.value(v->assignment[v],global_expression+terminal_expression)≈flow atol=1e-11
    end
    fix(x,12.0;force=true)
    @variable(global_model,global_arrival)
    @constraint(global_model,global_arrival==global_expression)
    @objective(global_model,Min,global_arrival)
    optimize!(global_model)
    @test termination_status(global_model)==OpenSHOP.MOI.OPTIMAL
    @test value(global_arrival)≈multireference_expected(curves,coefficients,12.0) atol=1e-8
    @test OpenSHOP._lift_distributed_references!(Model(),Dict())==Dict()

    # Two references stay a direct quadratic without auxiliary variables.
    pair=[first(curves),last(curves)]
    pair_model=Model()
    @variable(pair_model,0<=pair_q<=50)
    @test OpenSHOP.distributed_reference_coordinate!(pair_model,pair,pair_q)===nothing
    pair_expression=OpenSHOP.distributed_transfer_expression(pair_model,pair,pair_q,
        [0.1,0.5])
    @test pair_expression isa QuadExpr
    @test num_variables(pair_model)==1
    @test JuMP.value(v->12.0,pair_expression)≈
        multireference_expected(pair,[0.1,0.5],12.0)

    # A fixed curve needs no reference coordinate or nonlinear operator.
    fixed=[curves[2]]
    fixed_model=Model()
    @variable(fixed_model,0<=z<=50)
    @test OpenSHOP.distributed_reference_coordinate!(fixed_model,fixed,z)===nothing
    @test OpenSHOP.distributed_transfer_expression(fixed_model,fixed,z,[0.2])==0.2z
    @test num_variables(fixed_model)==1
    identical=[RiverRouting.DelayCurve(c.reference_flow,first(curves).edges,
        first(curves).weights) for c in curves]
    @test OpenSHOP.distributed_reference_coordinate!(fixed_model,identical,z)===nothing
    @test OpenSHOP.distributed_transfer_expression(fixed_model,identical,z,fill(0.2,4))==0.2z
    @test num_variables(fixed_model)==1
end
