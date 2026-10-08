using Test, JuMP, OpenSHOP, HiGHS

# Reuse the physical fixture and keep the original certificate regressions in
# this experiment's validation run. The production test suite stays unchanged.
include(joinpath(@__DIR__, "..", "test", "table_power_bounds_tests.jl"))

const PROOF_TABLE_POLICIES = (:baseline, :head_shared, :discharge_affine, :shared)

function proof_coordinate_relaxation(nodes; affine=false, known=nothing)
    m = Model(HiGHS.Optimizer)
    set_silent(m)
    @variable(m, 0 <= u <= 1)
    known !== nothing && fix(u, known; force=true)
    weights = length(nodes) == 1 ? [1.0] :
        @variable(m, [1:length(nodes)], lower_bound=0, upper_bound=1)
    length(nodes) > 1 && @constraint(m, sum(weights) == 1)
    coordinate = (nodes=Float64.(nodes), weights=weights)
    x = sum(nodes[j] * weights[j] for j in eachindex(nodes))
    if affine
        first(nodes)==0.0 && all(>(0.0),nodes[2:end]) && @constraint(m,x<=last(nodes)*u)
        moment = x
    else
        lo, hi = extrema(nodes)
        moment = @variable(m, lower_bound=min(0.0, lo), upper_bound=max(0.0, hi))
        # The shared lift uses the existing on-state moment, whose ordinary
        # binary-product relaxation remains in the physical model.
        @constraint(m, moment >= lo * u)
        @constraint(m, moment <= hi * u)
        @constraint(m, moment >= x - hi * (1 - u))
        @constraint(m, moment <= x - lo * (1 - u))
    end
    records = OpenSHOP._TablePowerCoordinateRecord[]
    onweights = OpenSHOP._table_power_onweights!(m, coordinate, u, moment,
        "proof_coordinate", records; affine)
    (; m, u, weights, onweights, moment, records)
end

function proof_old_gate_residuals(coefficients, weights, onweights, u)
    lo, hi = extrema(coefficients)
    expression = sum(coefficients[j] * weights[j] for j in eachindex(coefficients))
    z = sum(coefficients[j] * onweights[j] for j in eachindex(coefficients))
    (lo * u - z, z - hi * u, expression - hi * (1 - u) - z,
        z - expression + lo * (1 - u))
end

@testset "Shared lifts imply every former signed scalar gate inequality" begin
    # Optimize each former facet's violation over the whole shared relaxation;
    # this checks fractional states, rather than sampling only integer starts.
    for (nodes, coefficients) in (([-3.0, -1.0, 5.0], [-2.0, 3.0, -1.0]),
        ([-3.0, -1.0, 5.0], [5.0, -4.0, 2.0]),
        ([-3.0, -1.0, 5.0], [-2.0, -2.0, -2.0]), ([-2.0], [-3.0])),
        known in (nothing, 0.0, 1.0)
        r = proof_coordinate_relaxation(nodes; known)
        for residual in proof_old_gate_residuals(coefficients, r.weights, r.onweights, r.u)
            @objective(r.m, Max, residual)
            optimize!(r.m)
            @test termination_status(r.m) == JuMP.MOI.OPTIMAL
            @test objective_value(r.m) <= 1e-9
        end
        if known === nothing && length(nodes) > 1
            @test length(r.records) == 1
        else
            @test isempty(r.records)
        end
    end

    for coefficients in ([-2.0, 3.0, -1.0], [5.0, -4.0, 2.0], [-2.0, -2.0, -2.0]),
        known in (nothing, 0.0, 1.0)
        nodes = [0.0, 2.0, 5.0]
        r = proof_coordinate_relaxation(nodes; affine=true, known)
        @test isempty(r.records)
        @test num_variables(r.m) == 4
        mass_residual = sum(r.onweights) - r.u
        moment_residual = sum(nodes[j] * r.onweights[j] for j in eachindex(nodes)) - r.moment
        for residual in (mass_residual, -mass_residual, moment_residual, -moment_residual,
            proof_old_gate_residuals(coefficients, r.weights, r.onweights, r.u)...)
            @objective(r.m, Max, residual)
            optimize!(r.m)
            @test termination_status(r.m) == JuMP.MOI.OPTIMAL
            @test objective_value(r.m) <= 1e-9
        end
    end
end

function proof_gate_witness(nodes, weights_value, u_value, moment_value, coefficients;
    shared=false, affine=false)
    m = Model(HiGHS.Optimizer)
    set_silent(m)
    @variable(m, 0 <= u <= 1)
    fix(u, u_value; force=true)
    # The commitment must stay unfixed in the helper: this is a fractional
    # relaxation witness, not a model with known binary commitment.
    unfix(u)
    @constraint(m, u == u_value)
    @variable(m, 0 <= weights[1:length(nodes)] <= 1)
    @constraint(m, sum(weights) == 1)
    for j in eachindex(nodes)
        @constraint(m, weights[j] == weights_value[j])
    end
    coordinate = (nodes=Float64.(nodes), weights=weights)
    moment = affine ? sum(nodes[j] * weights[j] for j in eachindex(nodes)) : moment_value
    records = OpenSHOP._TablePowerCoordinateRecord[]
    onweights = shared ? OpenSHOP._table_power_onweights!(m, coordinate, u, moment,
        "witness", records; affine) : nothing
    gates = OpenSHOP._TablePowerGateRecord[]
    @variable(m, -10 <= p <= 10)
    for (i, row) in enumerate(coefficients)
        intercept = shared ? sum(row[j] * onweights[j] for j in eachindex(row)) :
            OpenSHOP._table_power_intercept!(m, row, weights, u, "witness_$i", gates)
        @constraint(m, p <= intercept)
    end
    @objective(m, Max, p)
    optimize!(m)
    @test termination_status(m) == JuMP.MOI.OPTIMAL
    objective_value(m)
end

@testset "Shared weights remove strict fractional gate witnesses" begin
    # Three independent gates can reuse the same on-state mass three times.
    rows = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    baseline = proof_gate_witness([-1.0, 1.0, 3.0], fill(1 / 3, 3), .5, .5, rows)
    shared = proof_gate_witness([-1.0, 1.0, 3.0], fill(1 / 3, 3), .5, .5, rows; shared=true)
    @test baseline ≈ 1 / 3 atol=1e-9
    @test shared ≈ 1 / 6 atol=1e-9

    # The existing on-state head moment supplies information that a scalar
    # intercept gate cannot express, even with just one support row.
    rows = [[1.0, 0.0]]
    baseline = proof_gate_witness([-2.0, 3.0], [.5, .5], .5, 1.0, rows)
    shared = proof_gate_witness([-2.0, 3.0], [.5, .5], .5, 1.0, rows; shared=true)
    @test baseline ≈ .5 atol=1e-9
    @test shared ≈ .1 atol=1e-9

    # Zero discharge has a unique coordinate vector. Its on-state first
    # weight is known affinely, so a separate scalar product overestimates it.
    baseline = proof_gate_witness([0.0, 2.0], [.5, .5], .5, 1.0, rows; affine=true)
    shared = proof_gate_witness([0.0, 2.0], [.5, .5], .5, 1.0, rows;
        shared=true, affine=true)
    @test baseline ≈ .5 atol=1e-9
    @test shared ≈ 0.0 atol=1e-9
end

@testset "Affine discharge preconditions retain a general shared fallback" begin
    # A coordinate without a unique zero-flow first node cannot use the
    # affine identity. Negative nodes and positive first nodes use a lift.
    for nodes in ([1.0, 2.0], [-1.0, 0.0, 2.0])
        r = proof_coordinate_relaxation(nodes; affine=true)
        @test length(r.records) == 1
        @test num_variables(r.m) == 1 + 2 * length(nodes)
        @objective(r.m, Max, sum(r.onweights) - r.u)
        optimize!(r.m)
        @test termination_status(r.m) == JuMP.MOI.OPTIMAL
        @test objective_value(r.m) <= 1e-9
    end
end

@testset "All policies retain identical certificates and complete physical starts" begin
    for interpolation in (:bilinear, :pchip_discharge), operating in (false, true),
        on in (0, 1), fixed in (false, true)
        c = table_power_bounds_case(; interpolation, positive_electrical=operating, operating)
        seed = dispatch_from_controls(c, reshape([on, 0], 2, 1),
            reshape([8.0 * on, 0.0], 2, 1), zeros(0, 1))
        @test seed["validation"]["valid"]
        baseline_supports = nothing
        baseline_counts = nothing
        for policy in PROOF_TABLE_POLICIES
            b = OpenSHOP._build_global_dispatch(c; joint=true,
                fixed_u=fixed ? seed["u"] : nothing, table_power_policy=policy)
            profile = b.m.ext[:global_table_power_bounds_profile]
            supports = [(unit=r.unit,interval=r.interval,axis=r.axis,slopes=r.slopes,coefficients=r.coefficients)
                for r in b.m.ext[:global_table_power_supports]]
            counts = (profile["certificate_evaluations"], profile["certificate_cache_hits"])
            if policy == :baseline
                baseline_supports = supports
                baseline_counts = counts
            else
                @test isequal(supports, baseline_supports)
                @test counts == baseline_counts
            end
            audit = OpenSHOP._lift_start!(b, c, seed)
            @test audit["valid"]
            @test audit["assigned"] == audit["variables"] == num_variables(b.m)
            @test audit["objective"] ≈ seed["objective"] atol=1e-6
            for record in b.m.ext[:global_table_power_coordinates], j in eachindex(record.onweights)
                @test start_value(record.onweights[j]) ≈
                    start_value(record.u) * start_value(record.weights[j]) atol=1e-12
            end
            if !fixed && on==1
                hull = b.m.ext[:global_power_hulls]["power_hull_2_1"]
                @test start_value(hull.head) < 0
                @test start_value(hull.eta) < 0
            end
        end
    end
end

@testset "Singleton and analytic cases remain valid under every policy" begin
    table = TurbineTable([1.0, 50.0, 100.0], [5.0, 10.0, 20.0],
        [.15 .2 .25; .85 .9 .95; .75 .8 .85], fill(5.0, 3), fill(20.0, 3);
        interpolation=:pchip_discharge)
    g = Generator(name=:Unit, plant=:Plant, qmin=5.0, qmax=20.0, pmin=1.0, pmax=25.0,
        hmin=60.0, hmax=60.0, minup=0.0, mindown=0.0, turbine_table=table)
    system = HydroSystem(reservoirs=Reservoir[], junctions=Junction[],
        boundaries=[Boundary(name=:Source, head=60.0), Boundary(name=:Tail, head=0.0)],
        tunnels=Tunnel[], rivers=River[], plants=[Plant(name=:Plant, source=:Source, target=:Tail, pmax=25.0)],
        generators=[g])
    c = ScheduleCase(name="proof_singleton", system=system, grid=[0.0, 1.0], prices=[100.0],
        operations=[OperationalSeries(object=:Unit, attribute=a, times=[0.0], values=[8.0])
            for a in (:qmin, :qmax)])
    for on in (0, 1), fixed in (false, true)
        seed = dispatch_from_controls(c, fill(on, 1, 1), fill(8.0 * on, 1, 1), zeros(0, 1))
        @test seed["validation"]["valid"]
        for policy in PROOF_TABLE_POLICIES
            b = OpenSHOP._build_global_dispatch(c; joint=true,
                fixed_u=fixed ? seed["u"] : nothing, table_power_policy=policy)
            @test OpenSHOP._lift_start!(b, c, seed)["valid"]
            @test isempty(b.m.ext[:global_table_power_coordinates])
        end
    end
    c = table_power_bounds_case(; analytic=true)
    seed = dispatch_from_controls(c, reshape([1, 0], 2, 1), reshape([8.0, 0.0], 2, 1), zeros(0, 1))
    @test seed["validation"]["valid"]
    for policy in PROOF_TABLE_POLICIES
        b = OpenSHOP._build_global_dispatch(c; joint=true, table_power_policy=policy)
        @test b.m.ext[:global_table_power_bounds_profile]["skipped"]["analytic"] == 2
        @test isempty(b.m.ext[:global_table_power_gates])
        @test isempty(b.m.ext[:global_table_power_coordinates])
        @test OpenSHOP._lift_start!(b, c, seed)["valid"]
    end
    @test_throws ArgumentError OpenSHOP._build_global_dispatch(c; joint=true, table_power_policy=:unknown)
end
