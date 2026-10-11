using ForwardDiff

@testset "Midpoint hydraulics: equations, derivatives and reusable storage" begin
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
    @test @inferred(OpenSHOP.hydraulic_residual!(f, x, data, work)) === f
    expected = [(x[1]-5.0)/(0.0036dt)+x[2], 90.0+2.0*(5.0+x[1])/2-x[2]^2]
    @test f ≈ expected
    saved = copy(f)
    @test (@allocated OpenSHOP.hydraulic_residual!(f, x, data, work)) == 0
    @test f == saved
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
end

@testset "Named controls and minimal model construction" begin
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
    @test isempty(HydroSystem().rivers)
    a, b = HydroSystem(), HydroSystem()
    push!(a.boundaries, Boundary(name = :Sea, head = 0.0))
    @test isempty(b.boundaries)
    @test_throws DimensionMismatch simulate(c; generator_q = zeros(1, 1), gate)
end
