using Test
@testset "Coupled short and zero river delays" begin
    shortsys(;
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
    reservoir(name, z; v = 2.0, slope = 0.0) = Reservoir(
        name = name,
        z0 = z,
        slope = slope,
        v0 = v,
        vmin = 0.01,
        vmax = 10.0,
        water_value = 1.0,
    )
    atom(
        name,
        source,
        target;
        delay = 0.0,
        law = :orifice,
        coefficient = 2.0,
        crest = 0.0,
        capacity = 100.0,
    ) = River(
        name = name,
        source = source,
        target = target,
        curves = RiverRouting.DelayCurve[],
        deterministic_delay = delay,
        capacity = capacity,
        law = law,
        coefficient = coefficient,
        crest = crest,
        water_value = 1.0,
    )
    a=reservoir(:a, 25.0)
    b=reservoir(:b, 100.0)
    basin=reservoir(:basin, 10.0; slope = 2.0)
    downstream=reservoir(:downstream, 0.0)
    reaches=[
        atom(:tributary_a, :a, :merge),
        atom(:tributary_b, :b, :merge),
        atom(:merged, :merge, :basin; law = :junction),
        atom(:outlet, :basin, :downstream; coefficient = 3.0),
    ]
    s=shortsys(
        reservoirs = [a, b, basin, downstream],
        river_junctions = [RiverJunction(name = :merge)],
        rivers = reaches,
    )
    c=ScheduleCase(
        name = "zero-delay confluence feedback",
        system = s,
        grid = [0.0, 1.0, 2.0],
        prices = zeros(2),
    )
    z=simulate(c, zeros(0, 2), ones(4, 2))
    @test z["converged"]
    @test z["river_release"][1, :]≈fill(10.0, 2) atol=1e-10
    @test z["river_release"][2, :]≈fill(20.0, 2) atol=1e-10
    @test z["river_release"][3, :]≈fill(30.0, 2) atol=1e-10
    @test z["arrival_volume"] ≈ 0.0036 .* z["river_release"] atol=1e-10
    @test maximum(abs, z["transit"])<1e-12
    # Scalar bisection independently gives the downstream basin's first step.
    lo, hi=1.0, 3.0
    for _ in 1:100
        mid=(lo+hi)/2
        f=mid - 2.0 - 0.0036*(30.0 - 3sqrt(10.0 + 2*(2.0 + mid)/2))
        if f>0
            hi=mid
        else
            lo=mid
        end
    end
    @test z["V"][3, 2]≈(lo+hi)/2 atol=1e-10
    @test z["V"][4, 2] - 2.0 ≈ 0.0036*z["river_release"][4, 1] atol=1e-10
    @test vec(sum(z["V"]; dims = 1))≈fill(8.0, 3) atol=1e-10
    prescribed=forward_step(
        s,
        [2.0, 2.0, 2.0, 2.0],
        1.0,
        Float64[],
        ones(4),
        [1.0, 2.0, 3.0, 4.0],
    )
    @test prescribed["river_arrival"]==[1.0, 2.0, 3.0, 4.0]
    @test prescribed["river_release"][3]≈3.0 atol=1e-10
    @test prescribed["converged"]

    # Fixed quarter-hour delay: 75% of a one-hour pulse arrives immediately.
    pulse=atom(:pulse, :a, :downstream; delay = 0.25)
    p=shortsys(reservoirs = [a, downstream], rivers = [pulse])
    pc=ScheduleCase(
        name = "quarter hour pulse",
        system = p,
        grid = [0.0, 1.0, 2.0, 3.0],
        prices = zeros(3),
    )
    q=simulate(pc, zeros(0, 3), reshape([1.0, 0.0, 0.0], 1, :))
    @test q["converged"]
    @test vec(q["arrival_volume"])≈[0.027, 0.009, 0.0] atol=1e-12
    @test vec(q["transit"])≈[0.0, 0.009, 0.0, 0.0] atol=1e-12
    @test vec(sum(q["V"]; dims = 1)+sum(q["transit"]; dims = 1))≈fill(4.0, 4) atol=1e-10
    fine=simulate(
        pc,
        zeros(0, 3),
        reshape([1.0, 0.0, 0.0], 1, :);
        grid = collect(0.0:0.25:3.0),
    )
    @test fine["converged"]
    @test fine["V"][:, end]≈q["V"][:, end] atol=1e-10
    @test vec(fine["arrival_volume"])≈[0.0, fill(0.009, 4)..., zeros(7)...] atol=1e-12

    # Uniform delays over [0,.5] have the same .75 same-step fraction.
    curves=[
        RiverRouting.DelayCurve(0.0, [0.0, 0.5], [1.0]),
        RiverRouting.DelayCurve(100.0, [0.0, 0.5], [1.0]),
    ]
    spread=River(
        name = :pulse,
        source = :a,
        target = :downstream,
        curves = curves,
        capacity = 100.0,
        law = :orifice,
        coefficient = 2.0,
        crest = 0.0,
        water_value = 1.0,
    )
    sc=ScheduleCase(
        name = "short distributed pulse",
        system = shortsys(reservoirs = [a, downstream], rivers = [spread]),
        grid = pc.grid,
        prices = pc.prices,
    )
    sd=simulate(sc, zeros(0, 3), reshape([1.0, 0.0, 0.0], 1, :))
    @test sd["converged"]
    @test sd["arrival_volume"]≈q["arrival_volume"] atol=1e-12
    @test sd["transit"]≈q["transit"] atol=1e-12
    # Blended contemporary release chooses its own diagonal transfer fraction.
    flowcurves=[
        RiverRouting.DelayCurve(0.0, [0.0, 0.25], [1.0]),
        RiverRouting.DelayCurve(20.0, [0.5, 0.75], [1.0]),
    ]
    flow=River(
        name = :flow,
        source = :a,
        target = :downstream,
        curves = flowcurves,
        capacity = 20.0,
        law = :orifice,
        coefficient = 2.0,
        crest = 0.0,
        water_value = 1.0,
    )
    fc=ScheduleCase(
        name = "flow dependent short pulse",
        system = shortsys(reservoirs = [a, downstream], rivers = [flow]),
        grid = pc.grid,
        prices = pc.prices,
    )
    fd=simulate(fc, zeros(0, 3), reshape([1.0, 0.0, 0.0], 1, :))
    @test fd["converged"]
    @test vec(fd["arrival_volume"])≈[0.0225, 0.0135, 0.0] atol=1e-12
    @test vec(fd["transit"])≈[0.0, 0.0135, 0.0, 0.0] atol=1e-12
    over=forward_step(
        fc.system,
        [2.0, 2.0],
        1.0,
        Float64[],
        [1.0],
        [0.0];
        current_transfer = ones(1, 2),
    )
    @test over["converged"]
    bad=River(
        name = :bad,
        source = :a,
        target = :downstream,
        curves = flowcurves,
        capacity = 20.0,
        law = :orifice,
        coefficient = 5.0,
        crest = 0.0,
        water_value = 1.0,
    )
    invalid=forward_step(
        shortsys(reservoirs = [a, downstream], rivers = [bad]),
        [2.0, 2.0],
        1.0,
        Float64[],
        [1.0],
        [0.0];
        current_transfer = ones(1, 2),
    )
    @test !invalid["converged"]
    @test !invalid["routing_domain_valid"]
    @test invalid["river_release"][1]≈25.0 atol=1e-10
    @test invalid["bound_violations"]["river_flow"]≈5.0 atol=1e-10

    # Historical water from an atom is counted as initial transit, then drained.
    historical=River(
        name = :h,
        source = :a,
        target = :downstream,
        curves = RiverRouting.DelayCurve[],
        deterministic_delay = 0.25,
        capacity = 100.0,
        law = :orifice,
        coefficient = 2.0,
        crest = 0.0,
        water_value = 1.0,
        history_grid = [-1.0, 0.0],
        history_release = [10.0],
    )
    hc=ScheduleCase(
        name = "historical atom",
        system = shortsys(reservoirs = [a, downstream], rivers = [historical]),
        grid = pc.grid,
        prices = pc.prices,
    )
    hh=simulate(hc, zeros(0, 3), zeros(1, 3))
    @test hh["converged"]
    @test hh["transit"][1, 1]≈0.009 atol=1e-12
    @test vec(hh["arrival_volume"])≈[0.009, 0.0, 0.0] atol=1e-12
    @test sum(hh["V"][:, end])≈4.009 atol=1e-10

    # Controlled releases support prescribed schedules independently of head.
    controlled=atom(
        :controlled,
        :a,
        :downstream;
        law = :controlled,
        capacity = 40.0,
        crest = 1000.0,
    )
    cs=shortsys(reservoirs = [a, downstream], rivers = [controlled])
    cc=ScheduleCase(
        name = "controlled release",
        system = cs,
        grid = pc.grid,
        prices = pc.prices,
    )
    cz=simulate(cc, zeros(0, 3), reshape([0.5, 0.25, 0.0], 1, :))
    @test cz["converged"]
    @test vec(cz["river_release"])≈[20.0, 10.0, 0.0] atol=1e-12
    @test vec(cz["arrival_volume"]) ≈ 0.0036 .* [20.0, 10.0, 0.0] atol=1e-12
    @test cz["V"][:, end]≈[1.892, 2.108] atol=1e-10

    # Midpoint hydraulics retains negative pipe flows together with routing.
    tunnel=Tunnel(
        name = :reverse,
        source = :a,
        target = :downstream,
        resistance = 0.25,
        capacity = 100.0,
    )
    reverse=shortsys(
        reservoirs = [reservoir(:a, 0.0), reservoir(:downstream, 100.0)],
        tunnels = [tunnel],
    )
    rz=forward_step(reverse, [2.0, 2.0], 1.0, Float64[], Float64[], Float64[])
    @test rz["converged"]
    @test rz["tunnel_q"][1]≈-20.0 atol=1e-10
    @test rz["Vnew"]≈[2.072, 1.928] atol=1e-10
    @test_throws ArgumentError forward_step(
        p,
        [2.0, 2.0],
        1.0,
        Float64[],
        [1.0],
        [0.0];
        current_transfer = fill(1.1, 1, 2),
    )
end
