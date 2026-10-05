using OpenSHOP:
    point_coefficients, point_arrival, arrival_knots, minimum_arrival, release_history

@testset "Exact instantaneous river arrivals" begin
    curve=RiverRouting.DelayCurve(0.0, [0.5, 1.5], [1.0])
    reach(; kw...) = River(
        name = :test,
        source = :source,
        target = :target,
        curves = [curve],
        capacity = 30.0,
        water_value = 1.0,
        history_grid = [-3.0, 0.0],
        history_release = [0.0];
        kw...,
    )
    r=reach()
    grid=[0.0, 1.0]
    q=[10.0]
    # Convolution of a unit-width pulse and unit-width delay bin is triangular.
    @test point_arrival(r, grid, q, 0.5)≈0.0
    @test point_arrival(r, grid, q, 1.0)≈5.0
    @test point_arrival(r, grid, q, 1.5)≈10.0
    @test point_arrival(r, grid, q, 2.0)≈5.0
    @test point_arrival(r, grid, q, 2.5)≈0.0
    @test size(point_coefficients(r, grid, 1.0))==(1, 2)
    @test point_coefficients(r, grid, 1.0; side = :left)==point_coefficients(
        r,
        grid,
        1.0;
        side = :right,
    )
    @test_throws ArgumentError point_coefficients(r, grid, 1.0; side = :invalid)
    @test_throws ArgumentError point_coefficients(r, [0.0, 0.0], 1.0)
    @test_throws ArgumentError point_arrival(r, grid, [-1.0], 1.0)
    # Independent integration over the delay density, with a nonuniform grid.
    nonuniform=[0.0, 0.7, 2.1, 3.0]
    flows=[4.0, 12.0, 2.0]
    n=20000
    for t in [0.1, 0.9, 1.3, 2.7, 3.8]
        quadrature=sum(
            begin
                d=0.5+(i-0.5)/n
                s=t-d
                k=findfirst(j->nonuniform[j]<=s<nonuniform[j + 1], eachindex(flows))
                k===nothing ? 0.0 : flows[k]
            end for i in 1:n
        )/n
        @test point_arrival(r, nonuniform, flows, t)≈quadrature atol=0.002
    end
    # Integrating point arrivals reproduces independent volume routing.
    fine=collect(0.0:0.001:3.0)
    mid=(fine[1:(end - 1)] .+ fine[2:end]) ./ 2
    @test 0.0036*0.001*sum(point_arrival(r, grid, q, t) for t in mid)≈0.036 atol=1e-10
    @test arrival_knots(r, [0.0, 3.0])==[0.0, 0.5, 1.5, 3.0]
end

@testset "Atom limits and historical releases" begin
    function atom(delay; history = [6.0])
        River(
            name = :atom,
            source = :a,
            target = :b,
            curves = RiverRouting.DelayCurve[],
            capacity = 30.0,
            water_value = 1.0,
            deterministic_delay = delay,
            history_grid = [-1.0, 0.0],
            history_release = history,
        )
    end
    g=[0.0, 0.5, 1.5, 2.0]
    q=[10.0, 20.0, 4.0]
    r=atom(0.0)
    @test point_arrival(r, g, q, 0.0; side = :left)==6.0
    @test point_arrival(r, g, q, 0.0; side = :right)==10.0
    @test point_arrival(r, g, q, 0.5; side = :left)==10.0
    @test point_arrival(r, g, q, 0.5; side = :right)==20.0
    @test point_arrival(r, g, q, 2.0; side = :left)==4.0
    @test point_arrival(r, g, q, 2.0; side = :right)==0.0
    @test minimum_arrival(r, g, q)["value"]==4.0
    rr=atom(0.25)
    @test point_arrival(rr, g, q, 0.1)==6.0
    @test point_arrival(rr, g, q, 0.25; side = :left)==6.0
    @test point_arrival(rr, g, q, 0.25; side = :right)==10.0
    @test point_arrival(rr, g, q, 0.75; side = :left)==10.0
    @test point_arrival(rr, g, q, 0.75; side = :right)==20.0
    @test arrival_knots(rr, g)==[0.0, 0.25, 0.75, 1.75, 2.0]
    @test minimum_arrival(rr, g, q)["value"]==4.0
    @test_throws ArgumentError point_coefficients(atom(-1.0), g, 0.0)
end

@testset "Flow-dependent delay and exact knot minimum" begin
    slow=RiverRouting.DelayCurve(0.0, [0.2, 0.6, 1.0], [0.0, 1.0])
    fast=RiverRouting.DelayCurve(20.0, [0.2, 0.6, 1.0], [1.0, 0.0])
    r=River(
        name = :mixed,
        source = :a,
        target = :b,
        curves = [slow, fast],
        capacity = 20.0,
        water_value = 1.0,
        history_grid = [-2.0, -0.5, 0.0],
        history_release = [10.0, 10.0],
    )
    g=[0.0, 0.7, 1.4, 2.0]
    q=[10.0, 0.0, 20.0]
    # q=10 selects a 50/50 mixture at release time, even after it ceases.
    K=point_coefficients(r, g, 0.8)
    hK=point_coefficients(r, r.history_grid, 0.8)
    expected=10*(K[1, 1]+K[1, 2])/2+20*K[3, 2]+sum(10*(hK[k, 1]+hK[k, 2])/2 for k in 1:2)
    @test point_arrival(r, g, q, 0.8)≈expected
    minimum=minimum_arrival(r, g, q)
    knots=arrival_knots(r, g)
    dense=collect(range(0.0, 2.0, length = 20001))
    @test minimum["value"]<=Base.minimum(point_arrival(r, g, q, t) for t in dense)+1e-10
    @test minimum["time"] in knots
    # Every segment is affine, so no missed interior lower value is possible.
    for (a, b) in zip(knots[1:(end - 1)], knots[2:end])
        @test point_arrival(r, g, q, (a+b)/2)≈(
            point_arrival(r, g, q, a)+point_arrival(r, g, q, b)
        )/2 atol=1e-9
    end
    @test_throws ArgumentError point_arrival(r, g, [21.0, 0.0, 0.0], 1.0)
end

@testset "Restart preserves release cohorts" begin
    g=[-2.0, -0.7, 0.0, 1.0, 3.0]
    q=[2.0, 8.0, 10.0, 4.0]
    h=release_history(g, q, 0.4)
    @test h.grid==[-2.0, -0.7, 0.0, 0.4]
    @test h.release==[2.0, 8.0, 10.0]
    @test release_history(g, q, 1.0).grid==[-2.0, -0.7, 0.0, 1.0]
    @test release_history(g, q, 4.0).grid==g
    @test isempty(release_history(g, q, -3.0).release)
    curve=RiverRouting.DelayCurve(0.0, [0.5, 1.5], [1.0])
    original=River(
        name = :r,
        source = :a,
        target = :b,
        curves = [curve],
        capacity = 20.0,
        water_value = 1.0,
        history_grid = g[1:3],
        history_release = q[1:2],
    )
    restarted=River(
        name = :r,
        source = :a,
        target = :b,
        curves = [curve],
        capacity = 20.0,
        water_value = 1.0,
        history_grid = h.grid,
        history_release = h.release,
    )
    for t in [0.4, 0.7, 1.0, 1.8, 2.8]
        @test point_arrival(original, [0.0, 1.0, 3.0], q[3:4], t)≈point_arrival(
            restarted,
            [0.4, 1.0, 3.0],
            q[3:4],
            t,
        )
    end
end
