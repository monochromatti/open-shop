using Test, JuMP, OpenSHOP
using OpenSHOP: river_order
include("rivergraph_tests.jl")
include("forward_tests.jl")
include("short_delay_tests.jl")
include("temporal_tests.jl")
include("curve_tests.jl")
include("outlet_head_tests.jl")
include("head_extrapolation_tests.jl")

@testset "Named object round trip and bounded model" begin
    lake=Reservoir(
        name = :Lake,
        z0 = 100.0,
        slope = 1.0,
        v0 = 2.0,
        vmin = 1.0,
        vmax = 3.0,
        water_value = 10.0,
    )
    unit=Generator(
        name = :Unit,
        plant = :Plant,
        qmin = 1.0,
        qmax = 5.0,
        pmin = 0.1,
        pmax = 5.0,
        qbest = 3.0,
        hbest = 100.0,
        hmin = 90.0,
        hmax = 110.0,
        qcurvature = 0.0,
        hcurvature = 0.0,
        minup = 0.0,
        mindown = 0.0,
    )
    sys=HydroSystem(
        reservoirs = [lake],
        junctions = Junction[],
        boundaries = [Boundary(name = :Sea, head = 0.0)],
        tunnels = Tunnel[],
        plants = [Plant(name = :Plant, source = :Lake, target = :Sea, pmax = 5.0)],
        generators = [unit],
        rivers = River[],
    )
    c=ScheduleCase(name = "core", system = sys, grid = [0.0, 1.0], prices = [50.0])
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    b=OpenSHOP._build_global_dispatch(c; joint = true)
    @test all(
        v->JuMP.is_binary(v) ||
           JuMP.is_fixed(v) ||
           (JuMP.has_lower_bound(v) && JuMP.has_upper_bound(v)),
        JuMP.all_variables(b.m),
    )
    x=solve_case(c; u = ones(Int, 1, 1), time_limit = 10.0)
    @test x["validation"]["valid"]
    @test validate(c, dispatch_from_controls(c, x["u"], x["generator_q"], x["gate"]))["valid"]
end

include("global_tests.jl")
