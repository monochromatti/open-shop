using Test, JuMP, OpenSHOP, SCIP

if !isdefined(Main, :table_power_bounds_case)
    include("table_power_bounds_tests.jl")
end

@testset "New targeted slopes enclose original interpolation and extensions" begin
    support, _, _, _ = OpenSHOP._table_power_support_cache()
    for interpolation in (:bilinear, :pchip_discharge)
        table = TurbineTable([1.0, 2.0, 3.0], [1.0, 1.5, 2.0],
            [.7 .5 .4; .95 .8 .6; .8 .7 .4], fill(1.0, 3), fill(2.0, 3); interpolation)
        qnodes = OpenSHOP._global_tensor_nodes(table.discharge, 0.0, 3.0)
        hnodes = OpenSHOP._global_tensor_nodes(table.heads, -1.0, 5.0)
        for (qbox, hbox) in (((1.2, 1.8), (1.4, 2.6)), ((0.25, 2.5), (0.7, 3.4)),
            ((1.3, 1.3), (2.2, 2.2))), fraction in OpenSHOP._POWER_CUT_FRACTIONS
            a = .00981 * .95 * hbox[2] * fraction
            b = .00981 * .95 * qbox[2] * fraction
            head = OpenSHOP._table_power_head_coefficients(table, qbox, hbox, hnodes, .95, a, support)
            discharge = OpenSHOP._table_power_discharge_coefficients(table, qbox, hbox, qnodes, .95, b, support)
            worst_head = -Inf
            worst_discharge = -Inf
            for q in range(qbox...; length=17), h in range(hbox...; length=19)
                power = .00981 * .95 * q * h * OpenSHOP.turbine_efficiency(table, q, h; extrapolation=:linear)
                mu = OpenSHOP._global_tensor_coordinate_weights(hnodes, h)
                lambda = OpenSHOP._global_tensor_coordinate_weights(qnodes, q)
                worst_head = max(worst_head, power-a*q-sum(head .* mu))
                worst_discharge = max(worst_discharge, power-b*h-sum(discharge .* lambda))
            end
            @test worst_head <= 1e-10
            @test worst_discharge <= 1e-10
        end
    end
end

@testset "Native affine cuts preserve every complete physical start" begin
    for interpolation in (:bilinear, :pchip_discharge), fixed in (false, true)
        c = table_power_bounds_case(; interpolation)
        seed = dispatch_from_controls(c, reshape([1, 0], 2, 1), reshape([8.0, 0.0], 2, 1), zeros(0, 1))
        @test seed["validation"]["valid"]
        b = OpenSHOP._build_global_dispatch(c; joint=true, fixed_u=fixed ? seed["u"] : nothing)
        @test OpenSHOP._lift_start!(b, c, seed)["valid"]
        set_optimizer(b.m, SCIP.Optimizer)
        JuMP.MOI.Utilities.attach_optimizer(backend(b.m))
        separator = OpenSHOP._install_power_cuts!(b, c; certificate_budget=0.0)
        variables = all_variables(b.m)
        assigned = Dict(v=>start_value(v) for v in variables)
        @test all(v->assigned[v] !== nothing, variables)
        native_by_reference = Dict(SCIP.VarRef(optimizer_index(v).value)=>Float64(assigned[v]) for v in variables)
        native = [native_by_reference[reference] for reference in separator.references]
        support, _, _, _ = OpenSHOP._table_power_support_cache()
        for (index, coordinate) in enumerate(separator.coordinates), fraction in OpenSHOP._POWER_CUT_FRACTIONS
            source = coordinate.source
            slope = .00981 * source.electrical_max *
                (source.axis == :head ? source.hbox[2] : source.qbox[2]) * fraction
            coefficients = source.axis == :head ?
                OpenSHOP._table_power_head_coefficients(source.table, source.qbox, source.hbox,
                    source.nodes, source.electrical_max, slope, support) :
                OpenSHOP._table_power_discharge_coefficients(source.table, source.qbox, source.hbox,
                    source.nodes, source.electrical_max, slope, support)
            row = OpenSHOP._power_cut_row(coordinate, slope, coefficients)
            residual = 40.0 * OpenSHOP._power_cut_value(row, v->native[v])
            record = b.m.ext[:global_table_power_supports][index]
            original = record.power-slope*record.term-
                sum(coefficients[j] * record.onweights[j] for j in eachindex(coefficients))
            @test residual ≈ OpenSHOP._table_power_value(original, assigned) atol=1e-10
            @test residual <= 1e-8
            if OpenSHOP._table_power_value(record.commitment, assigned) == 0
                @test abs(residual) <= 1e-12
            end
        end
        # Off head coordinates need not lie inside the on-state certificate box.
        # The rows still vanish because their shared on-state weights vanish.
        statistics = OpenSHOP._power_cut_statistics(separator)
        @test isempty(statistics["errors"])
        @test statistics["cuts_added_to_global_pool"] == 0
    end
end

function power_cut_native_witness(; certificate_budget=1.5, max_cuts=64, corrupt=false)
    # A constant efficiency graph has an exact tangent at fraction 0.375.
    # The four pre-existing slopes leave a strict gap at this fixed q/head.
    table = TurbineTable([1.0, 3.0], [1.0, 3.0], fill(.9, 2, 2), [1.0, 1.0], [3.0, 3.0])
    qbox = (1.0, 3.0)
    hbox = (1.0, 3.0)
    head = 1.25
    discharge = 2.0
    nodes = [head]
    support, _, _, _ = OpenSHOP._table_power_support_cache()
    slopes = [.00981 * hbox[2] * f for f in OpenSHOP._TABLE_POWER_FRACTIONS]
    coefficients = [OpenSHOP._table_power_head_coefficients(table, qbox, hbox, nodes, 1.0, a, support) for a in slopes]
    old_upper = minimum(slopes[j]*discharge + coefficients[j][1] for j in eachindex(slopes))
    exact_power = .00981 * .9 * discharge * head
    flow_cost = .00981 * .9 * head
    @test old_upper > exact_power + 1e-3

    m = Model(SCIP.Optimizer)
    set_silent(m)
    @variable(m, u, Bin)
    @constraint(m, u == 1)
    @variable(m, 1 <= q <= 3)
    @variable(m, 1 <= h_on <= 3)
    @constraint(m, h_on == head)
    @variable(m, 0 <= p <= 1)
    for j in eachindex(slopes)
        @constraint(m, p <= slopes[j]*q + coefficients[j][1]*u)
    end
    # The LP prefers (dummy, auxiliary)=(.75,1.5), while the best integer
    # solution is (1,1). Ordinary bound propagation cannot fix the binary.
    # This independent witness makes SCIP call root separators.
    @variable(m, dummy, Bin)
    @variable(m, 0 <= auxiliary <= 2)
    @constraint(m, auxiliary <= 2*dummy)
    @constraint(m, auxiliary <= 3-2*dummy)
    # Pricing discharge at the exact tangent slope makes the old envelope's
    # interior kink optimal. The new row contains two unfixed LP variables.
    @objective(m, Max, p - flow_cost*q + .1*auxiliary)
    set_optimizer_attribute(m, "presolving/maxrounds", 0)
    set_optimizer_attribute(m, "misc/usesymmetry", 0)
    set_optimizer_attribute(m, "limits/time", 10.0)
    JuMP.MOI.Utilities.attach_optimizer(backend(m))
    optimizer = unsafe_backend(m)
    @test SCIP.SCIPsetHeuristics(optimizer, SCIP.SCIP_PARAMSETTING_OFF, true) == SCIP.SCIP_OKAY
    @test SCIP.SCIPsetSeparating(optimizer, SCIP.SCIP_PARAMSETTING_OFF, true) == SCIP.SCIP_OKAY
    # Separation emphasis resets the pool frequency to ten. This small witness
    # must process its global pool cut before it branches on the dummy binary.
    set_optimizer_attribute(m, "separating/minefficacyroot", 0.0)
    set_optimizer_attribute(m, "separating/poolfreq", 1)
    set_optimizer_attribute(m, "cutselection/hybrid/minorthoroot", 0.0)
    m.ext[:global_table_power_supports] = [OpenSHOP._TablePowerSupportCoordinate(
        1, 1, :head, table, qbox, hbox, 1.0, nodes, OpenSHOP._TablePowerTerm[u],
        q, p, u, q, h_on)]
    b = (m=m, GQ=reshape([q], 1, 1))
    c = (prices=[100.0], grid=[0.0, 1.0])
    separator = OpenSHOP._install_power_cuts!(b, c; certificate_budget, max_cuts)
    corrupt && (table.efficiency[1,1] = NaN)
    optimize!(m)
    @test termination_status(m) == JuMP.MOI.OPTIMAL
    (; m, separator, old_upper, exact_power, flow_cost, q, p)
end

@testset "Native root separator adds and applies a certified global cut" begin
    witness = power_cut_native_witness()
    statistics = OpenSHOP._power_cut_statistics(witness.separator)
    @test isempty(statistics["errors"])
    @test statistics["calls"] > 0
    @test 0 < statistics["cuts_added_to_global_pool"] <= statistics["max_cuts"]
    @test statistics["native_cuts_applied"] > 0
    @test statistics["rounds"] <= statistics["max_rounds"]
    @test statistics["coordinate_checks"] <= statistics["max_coordinate_checks"]
    @test statistics["certificate_seconds"] >= statistics["polynomial_certificate_seconds"]
    @test abs(objective_value(witness.m) - .1) <= 1e-7
    @test value(witness.p) ≈ witness.flow_cost*value(witness.q) atol=1e-7
    @test any(eachindex(witness.separator.certified)) do index
        witness.separator.added[index] && witness.separator.certified[index].fraction == .375
    end

    disabled = power_cut_native_witness(; certificate_budget=0.0)
    disabled_statistics = OpenSHOP._power_cut_statistics(disabled.separator)
    @test isempty(disabled_statistics["errors"])
    @test disabled_statistics["support_vectors_certified"] == 0
    @test disabled_statistics["cuts_added_to_global_pool"] == 0
    @test objective_value(disabled.m) - .1 ≈ disabled.old_upper-disabled.exact_power atol=1e-7

    capped = power_cut_native_witness(; max_cuts=0)
    capped_statistics = OpenSHOP._power_cut_statistics(capped.separator)
    @test isempty(capped_statistics["errors"])
    @test capped_statistics["coordinate_checks"] == 0
    @test capped_statistics["cuts_added_to_global_pool"] == 0
    @test objective_value(capped.m) - .1 ≈ capped.old_upper-capped.exact_power atol=1e-7
end

@testset "Callback failures reject proof while preserving the audited schedule" begin
    failed = power_cut_native_witness(; corrupt=true)
    statistics = OpenSHOP._power_cut_statistics(failed.separator)
    @test !isempty(statistics["errors"])
    @test statistics["cuts_added_to_global_pool"] == 0
    c = table_power_bounds_case()
    seed = dispatch_from_controls(c, reshape([1,0],2,1), reshape([8.0,0.0],2,1), zeros(0,1))
    @test seed["validation"]["valid"]
    # Even an apparently exact bound must be discarded after a callback error.
    result = Dict{String,Any}("global_bound"=>seed["objective"], "global_certificate"=>true,
        "power_cut_statistics"=>statistics, "accepted"=>true)
    OpenSHOP._reject_power_cut_bound!(result)
    OpenSHOP._set_incumbent!(result, seed, "initial", 1e-4, 0.0)
    @test result["global_bound"] === nothing
    @test !result["global_certificate"]
    @test result["solution"] === seed && result["validation"]["valid"]
    @test result["accepted"]
    @test result["relative_gap"] === nothing
    @test occursin("power-cut", result["bound_rejection"])
end

@testset "Public solve uses one power-cut path and skips analytic models" begin
    for analytic in (false,true)
        c = table_power_bounds_case(; analytic)
        seed = dispatch_from_controls(c, reshape([1,0],2,1), reshape([8.0,0.0],2,1), zeros(0,1))
        @test seed["validation"]["valid"]
        result = OpenSHOP.solve(c; initial=seed, time_limit=20.0, relative_gap=1e-4)
        @test result["accepted"] && result["validation"]["valid"]
        @test result["start_audit"]["valid"]
        @test result["global_bound"] === nothing || result["global_bound"] >= seed["objective"]-1e-6
        statistics = result["power_cut_statistics"]
        if analytic
            @test statistics === nothing
        else
            @test statistics !== nothing
            @test isempty(statistics["errors"]) && statistics["infeasible_flags"] == 0
            @test statistics["cuts_added_to_global_pool"] <= statistics["max_cuts"]
        end
        @test !haskey(result, "solver_error")
        @test !haskey(result, "power_cut_statistics_error")
    end
end
