using ForwardDiff, SparseArrays

@testset "Midpoint hydraulics: equations and derivatives" begin
    reservoir = Reservoir(
        name = :Lake,
        z0 = 90.0,
        slope = 2.0,
        v0 = 5.0,
        vmin = 1.0,
        vmax = 10.0,
        water_value = 0.0,
    )
    tunnel =
        Tunnel(name = :Outlet, source = :Lake, target = :Sea, resistance = 1.0, capacity = 20.0)
    s = HydroSystem(
        reservoirs = [reservoir],
        tunnels = [tunnel],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
    )
    dt = 0.5
    empty = Float64[]
    data =
        OpenSHOP.HydraulicStepData(s, [5.0], dt, empty, empty, empty, [1.0], zeros(0, 2), nothing)
    work = OpenSHOP.HydraulicWorkspace(data)
    x = [4.9, 10.0]
    f = similar(x)
    OpenSHOP.hydraulic_residual!(f, x, data, work)
    expected = [(x[1]-5.0)/(0.0036dt)+x[2], 90.0+2.0*(5.0+x[1])/2-x[2]^2]
    @test f ≈ expected
    @test x == [4.9, 10.0]
    function residual(x)
        dual_work = OpenSHOP.HydraulicWorkspace(data, eltype(x))
        OpenSHOP.hydraulic_residual!(similar(x), x, data, dual_work)
    end
    @test ForwardDiff.jacobian(residual, x) ≈ [1/(0.0036dt) 1; 1 -20]
    result = forward_step(s, [5.0], dt, empty, empty, empty)
    α = 0.0036dt
    expected_q = (-α + sqrt(α^2 + 400))/2
    @test result["converged"]
    @test only(result["tunnel_q"]) ≈ expected_q atol=1e-9
    @test only(result["Vnew"]) ≈ 5.0-α*expected_q atol=1e-11
    @test only(result["boundary_net_inflow"]) ≈ expected_q
    @test_throws ArgumentError forward_step(s, [5.0], dt, empty, empty, empty; tolerance = NaN)
    @test_throws ArgumentError forward_step(s, [5.0], dt, empty, empty, empty; max_iterations = -1)
    @test_throws ArgumentError forward_step(
        s,
        [5.0],
        dt,
        empty,
        empty,
        empty;
        current_network = zeros(1, 1),
    )
    @test_throws DimensionMismatch forward_step(s, reshape([5.0], 1, 1), dt, empty, empty, empty)
    @test_throws ArgumentError forward_step(s, [5.0], big"1e-400", empty, empty, empty)
    @test_throws ArgumentError forward_step(s, [big"1e400"], dt, empty, empty, empty)
    @test_throws ArgumentError forward_step(
        s,
        [5.0],
        dt,
        empty,
        empty,
        empty;
        tunnel_opening = [big"1e-400"],
    )
    smallest_opening = [nextfloat(0.0)]
    tiny = forward_step(s, [5.0], dt, empty, empty, empty; tunnel_opening = smallest_opening)
    @test tiny["converged"]
    @test tiny["tunnel_opening"] == smallest_opening
    @test only(tiny["tunnel_q"]) > 0

end

@testset "Named controls preserve positional API behavior" begin
    f = MixedSourceNetworkExample.mixed_source_small_fixture()
    c = f.case
    q, gate, u = copy(f.q), copy(f.gates), copy(f.u)
    @test simulate(c; generator_q = q, gate) == simulate(c, q, gate)
    x = dispatch_from_controls(c; commitment = u, generator_q = q, gate)
    @test x == dispatch_from_controls(c, u, q, gate)
    @test x["validation"]["valid"]
    @test (q, gate, u) == (f.q, f.gates, f.u)
    @test replay_audit(c, x)["valid"]
    @test occursin(c.name, sprint(show, MIME"text/plain"(), c))
    a, b = HydroSystem(), HydroSystem()
    push!(a.boundaries, Boundary(name = :Sea, head = 0.0))
    @test isempty(b.boundaries)
    @test_throws DimensionMismatch simulate(c; generator_q = zeros(1, 1), gate)
end


@testset "Dense and sparse routing inputs produce the same physical trajectory" begin
    reservoirs = [
        Reservoir(
            name = :Upper,
            z0 = 90.0,
            slope = 1.0,
            v0 = 5.0,
            vmin = 1.0,
            vmax = 10.0,
            water_value = 0.0,
        ),
        Reservoir(
            name = :Lower,
            z0 = 50.0,
            slope = 1.0,
            v0 = 4.0,
            vmin = 1.0,
            vmax = 10.0,
            water_value = 0.0,
        ),
    ]
    rivers = [
        River(
            name = :Outlet,
            source = :Upper,
            target = :Merge,
            law = :controlled,
            capacity = 10.0,
            inflow = 1.0,
            water_value = 0.0,
            deterministic_delay = 0.0,
            curves = RiverRouting.DelayCurve[],
        ),
        River(
            name = :Merged,
            source = :Merge,
            target = :Lower,
            law = :junction,
            capacity = 20.0,
            inflow = 2.0,
            water_value = 0.0,
            deterministic_delay = 0.0,
            curves = RiverRouting.DelayCurve[],
        ),
    ]
    system = HydroSystem(; reservoirs, rivers, river_junctions = [RiverJunction(name = :Merge)])
    network = sparse(Float32[1 0; 1 1])
    saved = copy(network)
    gates = [0.4, 0.0]
    empty = Float64[]
    transfer = zeros(2, 1)
    data = OpenSHOP.HydraulicStepData(
        system,
        [5.0, 4.0],
        1.0,
        empty,
        gates,
        zeros(2),
        empty,
        transfer,
        network,
    )
    dense = forward_step(
        system,
        [5.0, 4.0],
        1.0,
        empty,
        gates,
        zeros(2);
        current_transfer = transfer,
        current_network = Matrix(network),
    )
    sparse_step = forward_step(
        system,
        [5.0, 4.0],
        1.0,
        empty,
        gates,
        zeros(2);
        current_transfer = transfer,
        current_network = network,
    )
    @test dense["converged"] && sparse_step["converged"]
    @test sparse_step["Vnew"] ≈ [5.0-0.0036*3.6, 4.0+0.0036*6.6] atol=1e-10
    @test sparse_step["river_release"] ≈ [4.6, 6.6]
    @test sparse_step["river_arrival"] ≈ [4.6, 6.6]
    @test sparse_step["Vnew"] ≈ dense["Vnew"] atol=1e-12
    @test network == saved
    x = [4.95, 4.1, 4.6]
    function sparse_residual(x)
        dual_work = OpenSHOP.HydraulicWorkspace(data, eltype(x))
        OpenSHOP.hydraulic_residual!(similar(x), x, data, dual_work)
    end
    @test ForwardDiff.jacobian(sparse_residual, x) ≈ [1/0.0036 0 1; 0 1/0.0036 -1; 0 0 1]
    @test_throws ArgumentError forward_step(
        system,
        [5.0, 4.0],
        1.0,
        empty,
        gates,
        zeros(2);
        current_transfer = transfer,
        current_network = fill(big"1e400", 2, 2),
    )
    @test_throws ArgumentError forward_step(
        system,
        [5.0, 4.0],
        1.0,
        empty,
        [big"1e-400", 0],
        zeros(2);
        current_transfer = transfer,
        current_network = network,
    )
    plant = Plant(name = :Station, source = :Upper, target = :Lower, pmax = 20.0)
    unit = Generator(
        name = :Unit,
        plant = :Station,
        qmin = 0.1,
        qmax = 10.0,
        pmin = 0.01,
        pmax = 20.0,
        hmin = 0.0,
        hmax = 200.0,
    )
    generating = HydroSystem(; reservoirs, plants = [plant], generators = [unit])
    @test_throws ArgumentError forward_step(
        generating,
        [5.0, 4.0],
        1.0,
        [big"1e-400"],
        empty,
        empty,
    )
end
