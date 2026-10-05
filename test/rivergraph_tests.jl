using Test

function graph_case(; grid = collect(0.0:1.0:8.0), history = false)
    rs=[
        Reservoir(
            name = n,
            z0 = 10.0,
            slope = 1.0,
            v0 = 1.0,
            vmin = 0.0,
            vmax = 5.0,
            water_value = 1.0,
        ) for n in [:A, :B, :C]
    ]
    cs=[
        RiverRouting.DelayCurve(0.0, [0.25, 0.75], [1.0]),
        RiverRouting.DelayCurve(30.0, [0.25, 0.75], [1.0]),
    ]
    reaches=[
        River(
            name = n,
            source = a,
            target = b,
            curves = cs,
            capacity = 30.0,
            water_value = 1.0,
            history_grid = [-2.0, -1.0, 0.0],
            history_release = history ? [2.0, 2.0] : [0.0, 0.0],
        ) for (n, a, b) in [(:AM, :A, :M), (:BM, :B, :M), (:MC, :M, :C)]
    ]
    sys=HydroSystem(
        reservoirs = rs,
        junctions = Junction[],
        boundaries = Boundary[],
        tunnels = Tunnel[],
        plants = Plant[],
        generators = Generator[],
        river_junctions = [RiverJunction(name = :M)],
        rivers = reaches,
    )
    ScheduleCase(name = "merge", system = sys, grid = grid, prices = zeros(length(grid)-1))
end
function replace_graph(c; rivers = c.system.rivers, junctions = c.system.river_junctions)
    s=c.system
    ScheduleCase(
        name = c.name,
        grid = c.grid,
        prices = c.prices,
        system = HydroSystem(
            reservoirs = s.reservoirs,
            junctions = s.junctions,
            boundaries = s.boundaries,
            tunnels = s.tunnels,
            plants = s.plants,
            generators = s.generators,
            river_junctions = junctions,
            rivers = rivers,
        ),
    )
end

@testset "Many-to-one river graph and inventory" begin
    for grid in [collect(0.0:1.0:8.0), [0.0, 0.2, 0.8, 1.5, 2.0, 3.5, 5.0, 8.0]],
        history in (false, true)

        c=graph_case(grid = grid, history = history)
        T=length(grid)-1
        q=zeros(3, T)
        q[1, 1]=10.0
        q[2, 1]=5.0
        out=route_network(c, q)
        a=out["arrival_volume"]
        rel=out["release"]
        W=out["transit"]
        @test out["order"]==[1, 2, 3]
        @test 0.0036 .* diff(grid) .* rel[3, :]≈a[1, :]+a[2, :] atol=1e-14
        for t in 1:T
            # Internal junction flows cancel: only reservoir sources and sink
            # arrivals enter total network inventory accounting.
            entering=0.0036*sum(diff(grid)[1:t] .* (q[1, 1:t]+q[2, 1:t]))
            @test sum(W[:, t + 1])-sum(W[:, 1])≈entering-sum(a[3, 1:t]) atol=2e-14
            for d in 1:3
                @test W[d, t + 1]-W[d, t]≈0.0036*diff(grid)[t]*rel[d, t]-a[d, t] atol=2e-14
            end
        end
        data=routing_data(c)
        @test size(data["B"][1])==(T, T, 2)
        @test data["history_initial"]≈W[:, 1]
        @test all(data["history_terminal"] .== 0.0)
    end
    c=graph_case()
    q=zeros(3, 8)
    q[1, 1]=10.0
    q[2, 1]=5.0
    o=route_network(c, q)
    @test sum(o["arrival_volume"][3, :])≈0.0036*15.0 atol=1e-14
    # Exact source-cohort routing is retained under refinement.
    fine=collect(0.0:0.125:8.0)
    f=route_network(c, q; grid = fine)
    @test f["arrival_volume"][1, :]≈RiverRouting.route(
        fine,
        c.grid,
        q[1, :],
        c.system.rivers[1].curves,
    )
    @test sum(f["arrival_volume"][3, :])≈sum(o["arrival_volume"][3, :]) atol=1e-14
    @test maximum(
        abs.(
            o["arrival_volume"][3, :]-[
                sum(f["arrival_volume"][3, (8t - 7):8t]) for t in 1:8
            ],
        ),
    )>1e-5
    @test_throws ArgumentError route_network(c, q; grid = [0.0, 9.0])
    @test_throws ArgumentError route_network(c, zeros(2, 8))
    q[1, 1]=31.0
    @test_throws ArgumentError route_network(c, q)
end

@testset "Reject invalid river graphs" begin
    c=graph_case()
    r=c.system.rivers
    @test_throws ArgumentError river_order(replace_graph(c; rivers = r[1:2]))
    @test_throws ArgumentError river_order(replace_graph(c; rivers = [r[3]]))
    @test_throws ArgumentError river_order(replace_graph(c; rivers = [r; r[3]]))
    cycle=River(
        name = :cycle,
        source = :C,
        target = :A,
        curves = r[1].curves,
        capacity = 30.0,
        water_value = 1.0,
    )
    @test_throws ArgumentError river_order(replace_graph(c; rivers = [r; cycle]))
    unknown=River(
        name = :unknown,
        source = :A,
        target = :Unknown,
        curves = r[1].curves,
        capacity = 30.0,
        water_value = 1.0,
    )
    @test_throws ArgumentError river_order(replace_graph(c; rivers = [r; unknown]))
    @test_throws ArgumentError river_order(
        replace_graph(c; junctions = [RiverJunction(name = :A)]),
    )
end

@testset "Flow-dependent merge agrees with optimization transfer fractions" begin
    c=graph_case()
    cs=[
        RiverRouting.DelayCurve(0.0, [0.1, 0.5, 1.0], [0.2, 0.8]),
        RiverRouting.DelayCurve(30.0, [0.1, 0.5, 1.0], [0.8, 0.2]),
    ]
    reaches=[
        River(
            name = r.name,
            source = r.source,
            target = r.target,
            capacity = 30.0,
            curves = cs,
            history_grid = r.history_grid,
            history_release = r.history_release,
            water_value = 1.0,
        ) for r in c.system.rivers
    ]
    c=replace_graph(c; rivers = reaches)
    q=zeros(3, 8)
    q[1, 1:3]=[10.0, 20.0, 4.0]
    q[2, 1:3]=[3.0, 7.0, 0.0]
    out=route_network(c, q)
    data=routing_data(c)
    for d in 1:3
        release=out["release"][d, :]
        B=data["B"][d]
        prediction=[
            sum(
                0.0036*diff(c.grid)[k]*release[k]*(
                    B[t, k, 1]+release[k]/30*(B[t, k, 2]-B[t, k, 1])
                ) for k in 1:8
            ) for t in 1:8
        ]
        @test prediction≈out["arrival_volume"][d, :] atol=1e-14
    end
    @test sum(out["arrival_volume"][3, :])≈0.0036*44.0 atol=1e-14
end
