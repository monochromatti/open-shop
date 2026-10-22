using Test, JuMP, OpenSHOP
using OpenSHOP: river_order
include("fixed_pump_tests.jl")
include("rivergraph_tests.jl")
include("forward_tests.jl")
include("hydraulic_solver_tests.jl")
include("mixed_source_forward_tests.jl")
include("short_delay_tests.jl")
include("temporal_tests.jl")
include("curve_tests.jl")
include("outlet_head_tests.jl")
include("head_extrapolation_tests.jl")
include("curve_policy_tests.jl")
include("schedule_tests.jl")
include("dispatch_data_tests.jl")
include("commitment_controls_tests.jl")
include("operating_controls_tests.jl")
include("zero_minimum_tests.jl")
include("model_consistency_tests.jl")
include("storage_integration_tests.jl")

@testset "Named object round trip and optimized dispatch" begin
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
    x=solve_case(c; u = ones(Int, 1, 1), time_limit = 10.0)
    @test x["validation"]["valid"]
    @test validate(c, dispatch_from_controls(c, x["u"], x["generator_q"], x["gate"]))["valid"]
end

include("global_tests.jl")

include("flow_requirement_tests.jl")
include("global_bounds_tests.jl")
include("diagnostics_tests.jl")

include("tensor_table_tests.jl")

include("power_bounds_tests.jl")
include("power_hull_tests.jl")
include("power_hull_start_tests.jl")
include("table_power_bounds_tests.jl")
include("table_power_coordinate_tests.jl")
include("power_cut_tests.jl")
include("solver_start_budget_tests.jl")

include("multireference_routing_tests.jl")
include("delay_policy_tests.jl")
include("mixed_source_network_tests.jl")
include("shared_tailrace_tests.jl")
include("river_split_tests.jl")

include("river_source_input_tests.jl")

include("level_flow_curves.jl")
include("tunnel_gate_losses.jl")
include("delta_channel_tests.jl")
include("tunnel_gate_controls_tests.jl")
include("tunnel_gate_regression_tests.jl")
include("directional_interlock_tests.jl")
include("gate_outlet_equation_tests.jl")
include("pumping_tests.jl")
include("terminal_value_tests.jl")

include("numerical_api_tests.jl")

include("replay_refactor_tests.jl")
include("model_reference_tests.jl")

include("result_api_tests.jl")

include("soft_constraint_tests.jl")
include("reserve_tests.jl")
include("river_reference_semantics_tests.jl")
include("penalty_settings_tests.jl")
include("unequal_turbine_grids_tests.jl")
