@testset "Short steps respect the resolution of each reservoir's storage" begin
    for volume in (1e-5, 600.0), dt in (1.0, 1e-7)
        inflow = 2.3
        lake = Reservoir(name=:Lake, z0=100.0, v0=volume, vmin=0.0,
            vmax=volume+1.0, inflow=inflow, water_value=0.0)
        system = HydroSystem(reservoirs=[lake])
        step = forward_step(system, [volume], dt, Float64[], Float64[], Float64[])
        expected = setprecision(256) do
            Float64(BigFloat(volume) + BigFloat(0.0036)*BigFloat(dt)*BigFloat(inflow))
        end
        @test step.converged
        @test abs(only(step.storage)-expected) <= eps(expected)
        @test abs(only(step.storage)-volume-0.0036dt*inflow) <= 2eps(expected)
        @test step.residual ≈ abs((only(step.storage)-volume)/(0.0036dt)-inflow) atol=1e-12
    end

    # A large reservoir's rounding allowance must not hide a different
    # reservoir's water imbalance, even on a very short interval.
    system = HydroSystem(reservoirs=[
        Reservoir(name=:Large, z0=100.0, v0=600.0, vmin=0.0, vmax=1000.0,
            water_value=0.0),
        Reservoir(name=:Small, z0=100.0, v0=1e-5, vmin=0.0, vmax=2e-5,
            water_value=0.0),
    ])
    bad = forward_step(system, [600.0, 1e-5], 1e-7, Float64[], Float64[], Float64[];
        initial_guess=[600.0, 1.1e-5], max_iterations=0)
    @test !bad.converged
    @test bad.storage[2] > 1e-5
end

@testset "Short-step rounding does not relax a tunnel's pressure law" begin
    lake = Reservoir(name=:Lake, z0=100.0, v0=600.0, vmin=0.0, vmax=1000.0,
        water_value=0.0)
    tunnel = Tunnel(name=:Outlet, source=:Lake, target=:Sea, resistance=1.0,
        capacity=10.0)
    system = HydroSystem(reservoirs=[lake], tunnels=[tunnel],
        boundaries=[Boundary(name=:Sea, head=90.0)])
    dt = 1e-7
    bad = forward_step(system, [600.0], dt, Float64[], Float64[], Float64[];
        initial_guess=[600.0-0.0036dt, 1.0], max_iterations=0)
    @test !bad.converged
    @test abs(bad.head[1]-bad.head[2]-only(bad.tunnel_q)^2) > 1.0
end

@testset "Converged channel flows retain the coupled water balance" begin
    lakes=[Reservoir(name=name,z0=100.0-1e-5,slope=1e-6,v0=10.0,
        vmin=9.0,vmax=11.0,water_value=0.0) for name in (:A,:B,:C)]
    curve=LevelFlowCurve([100.0],[TableCurve([-1.0,0.0,1.0],[-1.0,0.0,1.0])])
    channels=[River(name=name,source=a,target=b,law=:delta,capacity=1.0,
        deterministic_delay=0.0,curves=RiverRouting.DelayCurve[],water_value=0.0,
        delta_level_ref_up_flow_curve=curve) for (name,a,b) in ((:AB,:A,:B),(:BC,:B,:C))]
    system=HydroSystem(reservoirs=lakes,rivers=channels)
    tolerance=1e-9
    flow=[0.9tolerance,-0.9tolerance]
    net=[-flow[1],flow[1]-flow[2],flow[2]]
    storage=fill(10.0,3)+0.0036net
    # Both flow-law residuals and all balances already satisfy tolerance.
    # Their errors have opposite signs at the middle reservoir.
    step=forward_step(system,fill(10.0,3),1.0,Float64[],zeros(2),zeros(2);
        initial_guess=vcat(storage,flow),max_iterations=0,tolerance)
    @test step.converged
    net=[-step.river_release[1],step.river_release[1]-step.river_release[2],step.river_release[2]]
    @test maximum(abs,(step.storage.-10.0)./0.0036-net)<=tolerance
    @test maximum(abs,step.river_release)<=tolerance
end

@testset "Open tunnels recover signed flow from zero and opposed guesses" begin
    for delta in (-4.0, 4.0), guess in (0.0, -2sign(delta), 5sign(delta))
        system = HydroSystem(boundaries=[Boundary(name=:A, head=100.0+delta),
                Boundary(name=:B, head=100.0)],
            tunnels=[Tunnel(name=:Link, source=:A, target=:B, resistance=1.0,
                capacity=10.0)])
        step = forward_step(system, Float64[], 1.0, Float64[], Float64[], Float64[];
            initial_guess=[guess])
        @test step.converged
        @test only(step.tunnel_q) ≈ sign(delta)*sqrt(abs(delta)) atol=1e-9
    end
end

@testset "A small reservoir crosses a level-table knot" begin
    curve = TableCurve([0.0, 1e-5, 2e-5], [100.0, 101.0, 110.0])
    volume, dt = 1.5e-5, 0.001
    lake = Reservoir(name=:Lake, level_curve=curve, v0=volume, vmin=0.0,
        vmax=2e-5, water_value=0.0)
    system = HydroSystem(reservoirs=[lake],
        boundaries=[Boundary(name=:Sea, head=90.0)],
        tunnels=[Tunnel(name=:Outlet, source=:Lake, target=:Sea, resistance=1.0,
            capacity=10.0)])
    # Independently integrate the two affine segments to obtain the pressure
    # law, then bracket the scalar water balance without Newton derivatives.
    level(v) = v <= 1e-5 ? 100.0+1e5v : 101.0+9e5*(v-1e-5)
    function average_level(v)
        edges = sort!(unique!([v; filter(t -> v<t<volume, curve.x); volume]))
        sum((b-a)*(level(a)+level(b))/2 for (a,b) in zip(edges[1:end-1],edges[2:end])) /
            (volume-v)
    end
    lo, hi = 0.0, volume
    for _ in 1:80
        v = (lo+hi)/2
        balance = v-volume+0.0036dt*sqrt(average_level(v)-90.0)
        balance > 0 ? (hi=v) : (lo=v)
    end
    expected = (lo+hi)/2
    step = forward_step(system, [volume], dt, Float64[], Float64[], Float64[])
    @test step.converged
    @test only(step.storage) < 1e-5 < volume
    @test only(step.storage) ≈ expected atol=1e-14
    @test only(step.tunnel_q) ≈ (volume-expected)/(0.0036dt) atol=1e-8
end

@testset "Dry and shut orifices retain their physical zero-flow branch" begin
    for opening in (0.0, 1.0), delta in (-0.001, 0.0, 0.001)
        lake = Reservoir(name=:Lake, z0=98.0+delta, slope=1.0, v0=2.0,
            vmin=1.0, vmax=3.0, inflow=100.0, water_value=0.0)
        river = River(name=:Outlet, source=:Lake, target=:Sea, law=:orifice,
            coefficient=2.0, crest=100.0, allow_dry=true, capacity=10.0,
            deterministic_delay=0.0, curves=RiverRouting.DelayCurve[], water_value=0.0)
        system = HydroSystem(reservoirs=[lake], rivers=[river],
            boundaries=[Boundary(name=:Sea, head=0.0)])
        dt = 0.01
        discharge(v) = 2opening*sqrt(max(delta+(v-2.0)/2, 0.0))
        lo, hi = 2.0, 2.0+0.0036dt*100.0
        for _ in 1:80
            v = (lo+hi)/2
            v-2.0-0.0036dt*(100.0-discharge(v)) > 0 ? (hi=v) : (lo=v)
        end
        expected = (lo+hi)/2
        step = forward_step(system, [2.0], dt, Float64[], [opening], [0.0])
        @test step.converged
        @test only(step.storage) ≈ expected atol=1e-12
        @test only(step.river_release) ≈ discharge(expected) atol=1e-8
    end
end

@testset "Dry outlets retain a nonnegative flow domain during Newton" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=2.0,vmin=1.0,vmax=3.0,
        inflow=2.0,water_value=0.0)
    rivers=[River(name=name,source=:Lake,target=:Sea,law=:weir,crest=200.0,
        coefficient=1.0,allow_dry=true,capacity=10.0,deterministic_delay=0.0,
        curves=RiverRouting.DelayCurve[],water_value=0.0) for name in (:A,:B,:C)]
    system=HydroSystem(reservoirs=[lake],rivers=rivers,
        boundaries=[Boundary(name=:Sea,head=0.0)])
    # Zero laws can acquire tiny negative values in a coupled linear solve.
    # Retain the rate domain without sacrificing the storage equation.
    for guess in ([2.0,-1e-30,1e-30,-1e-30],[2.0,0.1,0.2,0.3])
        step=forward_step(system,[2.0],1.0,Float64[],ones(3),zeros(3);initial_guess=guess)
        @test step.converged
        @test all(>=(0.0),step.river_release)
        @test maximum(abs,step.river_release)<=1e-9
        @test abs((only(step.storage)-2.0)/0.0036+sum(step.river_release)-2.0)<=1e-9
    end
end
