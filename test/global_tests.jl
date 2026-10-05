function analytic_global_fixture(kind)
    price=kind in (:off, :forced_negative) ? -20.0 : 100.0
    eta=0.9
    k=0.00981*eta
    A=80.0
    b=0.036
    water=kind==:interior ? price*k*(A-2b*8.0)/0.0036 : 1000.0
    lake=Reservoir(
        name = :AnalyticLake,
        z0 = 40.0,
        slope = 20.0,
        v0 = 3.0,
        vmin = 1.0,
        vmax = 4.0,
        water_value = water,
    )
    unit=Generator(
        name = :AnalyticUnit,
        plant = :AnalyticPlant,
        qmin = 5.0,
        qmax = 20.0,
        pmin = 1.0,
        pmax = 25.0,
        efficiency = eta,
        qbest = 10.0,
        qcurvature = 0.0,
        hbest = A,
        hcurvature = 0.0,
        hmin = 35.0,
        hmax = 105.0,
        initial_on = 0,
        initial_age = 8.0,
        minup = 0.0,
        mindown = 0.0,
        startup = 0.0,
    )
    sys=HydroSystem(
        reservoirs = [lake],
        junctions = Junction[],
        boundaries = [Boundary(name = :AnalyticTail, head = 20.0)],
        tunnels = Tunnel[],
        plants = [
            Plant(
                name = :AnalyticPlant,
                source = :AnalyticLake,
                target = :AnalyticTail,
                pmax = 25.0,
                ramp = 100.0,
            ),
        ],
        generators = [unit],
        rivers = River[],
    )
    operations=kind==:forced_negative ?
               [
        OperationalSeries(
            object = :AnalyticUnit,
            attribute = :forced_on,
            times = [0.0],
            values = [1.0],
        ),
    ] : OperationalSeries[]
    c=ScheduleCase(
        name = "analytic_$(kind)",
        system = sys,
        grid = [0.0, 1.0],
        prices = [price],
        operations = operations,
    )
    # Concave positive-price quadratic: endpoints and stationary point exhaust
    # the on domain. For negative prices the quadratic is convex,
    # so endpoints suffice; an unforced unit also has the isolated off choice.
    candidates=kind==:forced_negative ? [5.0, 20.0] : [0.0, 5.0, 20.0]
    price>0 && push!(candidates, clamp((price*k*A-0.0036*water)/(2price*k*b), 5.0, 20.0))
    f(q) = price*k*q*(A-b*q)-0.0036*water*q
    q=candidates[argmax(f.(candidates))]
    (case = c, discharge = q, objective = f(q), on = Int(q>0))
end

@testset "Native global solve against independent analytic optima" begin
    for kind in (:off, :on, :interior, :forced_negative)
        fixture=analytic_global_fixture(kind)
        c=fixture.case
        initial=kind==:on ?
                dispatch_from_controls(
            c,
            ones(Int, 1, 1),
            fill(fixture.discharge, 1, 1),
            zeros(0, 1),
        ) : nothing
        result=OpenSHOP.solve(
            c;
            time_limit = 60.0,
            relative_gap = 1e-4,
            absolute_gap = 1e-3,
            initial,
        )
        @test !result["commitment_fixed"]
        @test result["solution"]!==nothing
        if result["solution"]!==nothing
            solution=result["solution"]
            @test solution["validation"]["valid"]
            @test result["accepted"]
            @test solution["u"][1, 1]==fixture.on
            # An objective-gap certificate permits a wider flow error around a
            # flat stationary point; the objective and bounds are checked above.
            flow_tolerance = kind == :interior ? 0.2 : 1e-3
            @test solution["generator_q"][1, 1]≈fixture.discharge atol=flow_tolerance
            @test solution["objective"]≈fixture.objective atol=1e-3
            @test result["global_certificate"]
            @test result["global_bound"]!==nothing
            if result["global_bound"]!==nothing
                @test isfinite(result["global_bound"])
                @test result["global_bound"]>=fixture.objective-1e-6
                @test result["global_bound"]>=result["feasible_lower_bound"]-1e-6
                gap=max(0.0, result["global_bound"]-result["feasible_lower_bound"])
                @test result["absolute_gap"]≈gap atol=1e-12
                @test result["relative_gap"]≈gap/max(
                    1.0,
                    abs(result["feasible_lower_bound"]),
                ) atol=1e-12
                @test gap<=1e-3 || result["relative_gap"]<=1e-4
            end
            if initial!==nothing
                @test result["start_audit"]["valid"]
                @test result["start_audit"]["assigned"]==result["start_audit"]["variables"]
                @test result["start_audit"]["constraints"]>0
                @test result["start_audit"]["objective"]≈initial["objective"] atol=1e-6
            end
        end
    end
end

@testset "A fixed-commitment certificate states its conditional scope" begin
    fixture = analytic_global_fixture(:on)
    c = fixture.case
    result = solve(c; fixed_u = zeros(Int, 1, 1), time_limit = 60.0, relative_gap = 1e-4)
    @test result["accepted"]
    @test result["global_certificate"]
    @test result["commitment_fixed"]
    @test occursin("fixed commitment", result["certificate_scope"])
    @test result["objective"] ≈ 0.0 atol=1e-8
    @test fixture.objective > result["global_bound"]
end
