using Test
@testset "Independent cascade integration" begin
    sys(;
        reservoirs = Reservoir[],
        junctions = Junction[],
        boundaries = Boundary[],
        tunnels = Tunnel[],
        plants = Plant[],
        generators = Generator[],
        river_junctions = RiverJunction[],
        rivers = River[],
    ) = HydroSystem(;
        reservoirs,
        junctions,
        boundaries,
        tunnels,
        plants,
        generators,
        river_junctions,
        rivers,
    )
    for (ha, hb, q) in [(100.0, 0.0, 20.0), (0.0, 100.0, -20.0), (10.0, 10.0, 0.0)]
        s=sys(
            boundaries = [Boundary(name = :a, head = ha), Boundary(name = :b, head = hb)],
            tunnels = [
                Tunnel(
                    name = :p,
                    source = :a,
                    target = :b,
                    resistance = 0.25,
                    capacity = 50.0,
                ),
            ],
        )
        z=forward_step(s, Float64[], 1.0, Float64[], Float64[], Float64[])
        @test z["converged"]
        @test only(z["tunnel_q"])≈q atol=1e-9
        @test sum(z["boundary_net_inflow"])≈0 atol=1e-9
    end
    # Nonlinear storage-head equation coupled to an orifice; compare against
    # an independently bracketed scalar bisection, not Newton's own residual.
    r=Reservoir(
        name = :a,
        z0 = 20.0,
        slope = 3.0,
        curvature = 0.5,
        v0 = 2.0,
        vmin = 0.1,
        vmax = 4.0,
        water_value = 1.0,
        inflow = 2.0,
    )
    curve=[RiverRouting.DelayCurve(0.0, [1.0, 2.0], [1.0])]
    reach=River(
        name = :reach,
        source = :a,
        target = :out,
        curves = curve,
        capacity = 100.0,
        law = :orifice,
        coefficient = 5.0,
        crest = 10.0,
        water_value = 1.0,
    )
    s=sys(
        reservoirs = [r],
        boundaries = [Boundary(name = :out, head = 0.0)],
        rivers = [reach],
    )
    z=forward_step(s, [2.0], 1.0, Float64[], [0.7], [0.0])
    lo, hi=0.0, 4.0
    for _ in 1:100
        mid=(lo+hi)/2
        vmean=(2+mid)/2
        f=mid-2-0.0036*(2-5*0.7*sqrt(20+3vmean+0.5vmean^2-10))
        if f>0
            hi=mid
        else
            lo=mid
        end
    end
    @test z["converged"]
    @test z["Vnew"][1]≈(lo+hi)/2 atol=1e-10
    @test z["river_release"][1]≈5*0.7*sqrt(z["H"][1]-10) atol=1e-10
    shut=forward_step(s, [2.0], 1.0, Float64[], [0.0], [0.0]; initial_guess = z)
    @test shut["river_release"][1]≈0 atol=1e-10
    @test shut["Vnew"][1]≈2+0.0036*2 atol=1e-10
    # Weir is passive, independent of the gate schedule.
    wr=River(
        name = :w,
        source = :a,
        target = :out,
        curves = curve,
        capacity = 100.0,
        law = :weir,
        coefficient = 0.15,
        crest = 10.0,
        water_value = 1.0,
    )
    ws=sys(reservoirs = [r], boundaries = s.boundaries, rivers = [wr])
    a=forward_step(ws, [2.0], 1.0, Float64[], [0.0], [0.0])
    b=forward_step(ws, [2.0], 1.0, Float64[], [1.0], [0.0])
    @test a["converged"] && b["converged"]
    @test a["river_release"]≈b["river_release"] atol=1e-10
    @test a["river_release"][1]≈0.15*(a["H"][1]-10)^1.5 atol=1e-10
    # Two tributaries merge at zero-storage junction and enter another reach.
    ra=Reservoir(
        name = :a,
        z0 = 30.0,
        slope = 0.0,
        v0 = 2.0,
        vmin = 0.1,
        vmax = 4.0,
        water_value = 1.0,
    )
    rb=Reservoir(
        name = :b,
        z0 = 40.0,
        slope = 0.0,
        v0 = 2.0,
        vmin = 0.1,
        vmax = 4.0,
        water_value = 1.0,
    )
    rivers=[
        River(
            name = :one,
            source = :a,
            target = :merge,
            curves = curve,
            capacity = 100.0,
            law = :orifice,
            coefficient = 2.0,
            crest = 5.0,
            water_value = 1.0,
            history_release = fill(1.0, 4),
        ),
        River(
            name = :two,
            source = :b,
            target = :merge,
            curves = curve,
            capacity = 100.0,
            law = :orifice,
            coefficient = 2.0,
            crest = 15.0,
            water_value = 1.0,
            history_release = fill(2.0, 4),
        ),
        River(
            name = :three,
            source = :merge,
            target = :out,
            curves = curve,
            capacity = 100.0,
            water_value = 1.0,
            history_release = fill(3.0, 4),
        ),
    ]
    s=sys(
        reservoirs = [ra, rb],
        boundaries = [Boundary(name = :out, head = 0.0)],
        river_junctions = [RiverJunction(name = :merge)],
        rivers = rivers,
    )
    c=ScheduleCase(
        name = "merge",
        system = s,
        grid = collect(0.0:1.0:8.0),
        prices = zeros(8),
    )
    z=simulate(c, zeros(0, 8), ones(3, 8))
    @test z["converged"]
    @test z["river_release"][1:2, :]≈fill(10.0, 2, 8) atol=1e-10
    @test z["river_release"][3, :]≈vec(sum(z["arrival_volume"][1:2, :]; dims = 1))/0.0036 atol=1e-10
    total=vec(sum(z["V"]; dims = 1)+sum(z["transit"]; dims = 1))
    @test maximum(abs.(total[2:end]+cumsum(0.0036 .* z["boundary_outflow"]) .- total[1]))<1e-10
    for (d, reach) in enumerate(rivers)
        w=RiverRouting.remaining_volume(
            8.0,
            reach.history_grid,
            reach.history_release,
            reach.curves,
        )+RiverRouting.remaining_volume(8.0, c.grid, z["river_release"][d, :], reach.curves)
        @test z["transit"][d, end]≈w atol=1e-10
    end
    fine=simulate(c, zeros(0, 8), ones(3, 8); grid = collect(0.0:0.25:8.0))
    @test fine["converged"]
    @test sum(fine["V"][:, end])+sum(fine["transit"][:, end])+0.0036*0.25*sum(
        fine["boundary_outflow"],
    )≈total[1] atol=1e-10
    @test_throws ArgumentError simulate(
        c,
        zeros(0, 8),
        ones(3, 8);
        grid = collect(0.0:2.0:8.0),
    )
    @test_throws ArgumentError forward_step(
        s,
        [2.0, 2.0],
        0.0,
        Float64[],
        ones(3),
        zeros(3),
    )
end
