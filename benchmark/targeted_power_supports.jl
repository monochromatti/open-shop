# Experimental, bounded root separation. Every row is certified over the
# original physical on-state box; LP values only select which rows to try.
using SCIP, JuMP, OpenSHOP

const TARGETED_POWER_FRACTIONS = (0.125, 0.375, 0.625, 0.875, 1.125)

struct TargetedSupportAffine
    constant::Float64
    variables::Vector{SCIP.VarRef}
    coefficients::Vector{Float64}
end

function targeted_support_affine(expression, variable_map)
    expression isa Number && return TargetedSupportAffine(Float64(expression), SCIP.VarRef[], Float64[])
    expression isa VariableRef && return TargetedSupportAffine(0.0, [variable_map[expression]], [1.0])
    expression isa GenericAffExpr || throw(ArgumentError("targeted supports require affine existing coordinates"))
    terms = Dict{SCIP.VarRef,Float64}()
    for (coefficient, variable) in linear_terms(expression)
        reference = variable_map[variable]
        terms[reference] = get(terms, reference, 0.0) + coefficient
    end
    variables = sort!([v for (v, coefficient) in terms if coefficient != 0.0]; by=v->v.val)
    TargetedSupportAffine(Float64(expression.constant), variables, [terms[v] for v in variables])
end

function targeted_support_value(expression::TargetedSupportAffine, at)
    expression.constant + sum(expression.coefficients[j] * at(expression.variables[j])
        for j in eachindex(expression.variables); init=0.0)
end

struct TargetedSupportCoordinate
    unit::Int
    interval::Int
    axis::Symbol
    table::TurbineTable
    qbox::NTuple{2,Float64}
    hbox::NTuple{2,Float64}
    electrical_max::Float64
    nodes::Vector{Float64}
    onweights::Vector{TargetedSupportAffine}
    term::TargetedSupportAffine
    power::TargetedSupportAffine
    commitment::TargetedSupportAffine
    discharge::TargetedSupportAffine
    onhead::TargetedSupportAffine
    objective_weight::Float64
end

struct TargetedCertifiedSupport
    coordinate::Int
    fraction::Float64
    slope::Float64
    coefficients::Vector{Float64}
    expression::TargetedSupportAffine
end

"""Expand signed affine on-weights, retaining their constants in the row RHS."""
function targeted_support_row(coordinate, slope, coefficients)
    terms = Dict{SCIP.VarRef,Float64}()
    constant = 0.0
    for (multiplier, expression) in ((1.0, coordinate.power), (-slope, coordinate.term))
        constant += multiplier * expression.constant
        for j in eachindex(expression.variables)
            v = expression.variables[j]
            terms[v] = get(terms, v, 0.0) + multiplier * expression.coefficients[j]
        end
    end
    for j in eachindex(coefficients)
        expression = coordinate.onweights[j]
        constant -= coefficients[j] * expression.constant
        for k in eachindex(expression.variables)
            v = expression.variables[k]
            terms[v] = get(terms, v, 0.0) - coefficients[j] * expression.coefficients[k]
        end
    end
    variables = sort!([v for (v, coefficient) in terms if coefficient != 0.0]; by=v->v.val)
    result = TargetedSupportAffine(constant / 40.0, variables, [terms[v] / 40.0 for v in variables])
    isfinite(result.constant) && all(isfinite, result.coefficients) || error("nonfinite certified cut")
    result
end

mutable struct TargetedPowerSupports{F} <: SCIP.AbstractSeparator
    optimizer::SCIP.Optimizer
    coordinates::Vector{TargetedSupportCoordinate}
    support::F
    certificate_evaluations::Base.RefValue{Int}
    certificate_hits::Base.RefValue{Int}
    polynomial_seconds::Base.RefValue{Float64}
    certified::Vector{TargetedCertifiedSupport}
    checked::Set{Int}
    added::Set{Tuple{Int,Float64}}
    certificate_budget::Float64
    max_cuts::Int
    max_rounds::Int
    max_cuts_per_round::Int
    max_coordinates_per_round::Int
    max_coordinate_checks::Int
    violation_tolerance::Float64
    calls::Int
    skipped::Int
    rounds::Int
    certificate_seconds::Float64
    callback_seconds::Float64
    maximum_violation_mw::Float64
    cuts::Vector{NamedTuple{(:unit,:interval,:axis,:fraction,:violation_mw),Tuple{Int,Int,String,Float64,Float64}}}
    infeasible_flags::Int
    errors::Vector{String}
end

"""Install the experimental separator after JuMP has attached SCIP.

Only existing shared on-state coordinates are used. There are no extra
variables, local domain certificates, or changes to the physical equations.
The certificate allowance is checked before each coefficient vector; one
already-started vector may finish beyond it. It is never a solver time limit.
"""
function install_targeted_supports(b, c; certificate_budget=1.5, max_cuts=64,
        max_rounds=4, max_cuts_per_round=16, max_coordinates_per_round=16,
        max_coordinate_checks=64, violation_tolerance=1e-5)
    isfinite(certificate_budget) && certificate_budget >= 0 || throw(ArgumentError("nonnegative certificate allowance required"))
    isfinite(violation_tolerance) && violation_tolerance >= 0 || throw(ArgumentError("nonnegative violation tolerance required"))
    all(x->x isa Integer && x >= 0,
        (max_cuts, max_rounds, max_cuts_per_round, max_coordinates_per_round, max_coordinate_checks)) ||
        throw(ArgumentError("nonnegative integer separator limits required"))
    optimizer = unsafe_backend(b.m)
    optimizer isa SCIP.Optimizer || throw(ArgumentError("attached SCIP optimizer required"))
    variable_map = Dict(v=>SCIP.VarRef(optimizer_index(v).value) for v in all_variables(b.m))
    all(reference->haskey(optimizer.inner.vars, reference), values(variable_map)) || error("SCIP variable mapping incomplete")
    coordinates = TargetedSupportCoordinate[]
    for record in get(b.m.ext, :global_table_power_supports, [])
        record.onweights === nothing && throw(ArgumentError("targeted supports require both shared coordinate lifts"))
        record.axis in (:head, :discharge) || error("unsupported power support axis")
        hull = b.m.ext[:global_power_hulls]["power_hull_$(record.unit)_$(record.interval)"]
        weight = max(0.0, c.prices[record.interval]) * (c.grid[record.interval+1] - c.grid[record.interval])
        push!(coordinates, TargetedSupportCoordinate(record.unit, record.interval, record.axis,
            record.table, Float64.(record.qbox), Float64.(record.hbox), Float64(record.em),
            Float64.(record.nodes), [targeted_support_affine(w, variable_map) for w in record.onweights],
            targeted_support_affine(record.term, variable_map), targeted_support_affine(record.power, variable_map),
            targeted_support_affine(record.u, variable_map), targeted_support_affine(b.GQ[record.unit,record.interval], variable_map),
            targeted_support_affine(hull.onhead, variable_map), weight))
    end
    support, evaluations, hits, polynomial_seconds = OpenSHOP._table_power_support_cache()
    separator = TargetedPowerSupports(optimizer, coordinates, support, evaluations, hits, polynomial_seconds,
        TargetedCertifiedSupport[], Set{Int}(), Set{Tuple{Int,Float64}}(), Float64(certificate_budget),
        Int(max_cuts), Int(max_rounds), Int(max_cuts_per_round), Int(max_coordinates_per_round),
        Int(max_coordinate_checks), Float64(violation_tolerance), 0, 0, 0, 0.0, 0.0, 0.0,
        NamedTuple{(:unit,:interval,:axis,:fraction,:violation_mw),Tuple{Int,Int,String,Float64,Float64}}[], 0, String[])
    SCIP.include_sepa(optimizer, separator; name="targeted_table_power",
        description="Bounded globally certified table-power supports", priority=100000,
        freq=0, maxbounddist=1.0, delay=false)
    separator
end

function targeted_support_add_row!(separator, expression)
    optimizer = separator.optimizer
    row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(optimizer, row, optimizer.inner.sepas[separator],
        "", -SCIP.SCIPinfinity(optimizer), -expression.constant, false, false, false)
    try
        variables = [optimizer.inner.vars[v][] for v in expression.variables]
        SCIP.@SCIP_CALL SCIP.SCIPaddVarsToRow(optimizer, row[], length(variables), variables, expression.coefficients)
        # SCIP.jl's global add_cut_sepa only adds to the pool. Also submit this
        # same global row to the current separation store so it acts at once.
        SCIP.@SCIP_CALL SCIP.SCIPaddPoolCut(optimizer, row[])
        infeasible = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPaddRow(optimizer, row[], true, infeasible)
        return Bool(infeasible[])
    finally
        SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(optimizer, row)
    end
end

function targeted_support_score(coordinate, at)
    coordinate.objective_weight <= 0 && return 0.0
    u = clamp(targeted_support_value(coordinate.commitment, at), 0.0, 1.0)
    p = targeted_support_value(coordinate.power, at)
    if u <= 1e-8
        return coordinate.objective_weight * max(0.0, p)
    end
    q = clamp(targeted_support_value(coordinate.discharge, at) / u, coordinate.qbox...)
    h = clamp(targeted_support_value(coordinate.onhead, at) / u, coordinate.hbox...)
    efficiency = OpenSHOP.turbine_efficiency(coordinate.table, q, h; extrapolation=:linear)
    physical = u * 0.00981 * coordinate.electrical_max * q * h * efficiency
    # This is a selection heuristic, never a bound. The tiny tie-breaker permits
    # checking a few positively priced units with no perspective residual.
    coordinate.objective_weight * (max(0.0, p-physical) + 1e-10 * max(0.0, p))
end

function targeted_support_certify!(separator, index)
    coordinate = separator.coordinates[index]
    push!(separator.checked, index)
    scale = 0.00981 * coordinate.electrical_max *
        (coordinate.axis == :head ? coordinate.hbox[2] : coordinate.qbox[2])
    for fraction in TARGETED_POWER_FRACTIONS
        separator.certificate_seconds >= separator.certificate_budget && break
        began_ns = time_ns()
        slope = scale * fraction
        coefficients = coordinate.axis == :head ?
            OpenSHOP._table_power_head_coefficients(coordinate.table, coordinate.qbox, coordinate.hbox,
                coordinate.nodes, coordinate.electrical_max, slope, separator.support) :
            OpenSHOP._table_power_discharge_coefficients(coordinate.table, coordinate.qbox, coordinate.hbox,
                coordinate.nodes, coordinate.electrical_max, slope, separator.support)
        all(isfinite, coefficients) && isfinite(slope) || error("nonfinite targeted support certificate")
        expression = targeted_support_row(coordinate, slope, coefficients)
        push!(separator.certified, TargetedCertifiedSupport(index, fraction, slope, coefficients, expression))
        separator.certificate_seconds += Float64(time_ns()-began_ns) / 1e9
    end
end

function SCIP.exec_lp(separator::TargetedPowerSupports)
    began_ns = time_ns()
    separator.calls += 1
    try
        optimizer = separator.optimizer
        node = SCIP.SCIPgetFocusNode(optimizer)
        if node == C_NULL || SCIP.SCIPnodeGetDepth(node) != 0 ||
                Bool(SCIP.SCIPinProbing(optimizer)) || Bool(SCIP.SCIPinDive(optimizer)) ||
                SCIP.SCIPgetLPSolstat(optimizer) != SCIP.SCIP_LPSOLSTAT_OPTIMAL ||
                !Bool(SCIP.SCIPisLPPrimalReliable(optimizer)) || !Bool(SCIP.SCIPisLPRelax(optimizer)) ||
                separator.rounds >= separator.max_rounds || length(separator.added) >= separator.max_cuts ||
                isempty(separator.coordinates) || !isempty(separator.errors)
            separator.skipped += 1
            return SCIP.SCIP_DIDNOTRUN
        end
        separator.rounds += 1
        values = Dict{SCIP.VarRef,Float64}()
        at = function(variable)
            get!(values, variable) do
                result = Float64(SCIP.SCIPgetSolVal(optimizer, C_NULL, optimizer.inner.vars[variable][]))
                isfinite(result) || error("nonfinite root LP value")
                result
            end
        end
        if separator.certificate_seconds < separator.certificate_budget &&
                length(separator.checked) < separator.max_coordinate_checks
            remaining = [j for j in eachindex(separator.coordinates) if !(j in separator.checked)]
            scores = Dict(j=>targeted_support_score(separator.coordinates[j], at) for j in remaining)
            filter!(j->scores[j] > 0.0, remaining)
            sort!(remaining; by=j->(-scores[j], j))
            count = min(separator.max_coordinates_per_round,
                separator.max_coordinate_checks-length(separator.checked), length(remaining))
            for index in first(remaining, count)
                separator.certificate_seconds >= separator.certificate_budget && break
                targeted_support_certify!(separator, index)
            end
        end
        candidates = Tuple{Float64,Float64,Int}[]
        for (index, cut) in enumerate(separator.certified)
            (cut.coordinate, cut.fraction) in separator.added && continue
            violation = 40.0 * targeted_support_value(cut.expression, at)
            violation > separator.violation_tolerance || continue
            weight = separator.coordinates[cut.coordinate].objective_weight
            push!(candidates, (weight * violation, violation, index))
        end
        sort!(candidates; by=x->(-x[1], -x[2], x[3]))
        count = min(separator.max_cuts_per_round, separator.max_cuts-length(separator.added), length(candidates))
        for (_, violation, index) in first(candidates, count)
            cut = separator.certified[index]
            coordinate = separator.coordinates[cut.coordinate]
            expression = cut.expression
            infeasible = targeted_support_add_row!(separator, expression)
            push!(separator.added, (cut.coordinate, cut.fraction))
            push!(separator.cuts, (unit=coordinate.unit, interval=coordinate.interval,
                axis=string(coordinate.axis), fraction=cut.fraction, violation_mw=violation))
            separator.maximum_violation_mw = max(separator.maximum_violation_mw, violation)
            if infeasible
                separator.infeasible_flags += 1
                error("targeted global row reported root infeasibility; experiment requires independent review")
            end
        end
        return count > 0 ? SCIP.SCIP_SEPARATED : SCIP.SCIP_DIDNOTFIND
    catch exception
        push!(separator.errors, sprint(showerror, exception))
        return SCIP.SCIP_DIDNOTFIND
    finally
        separator.callback_seconds += Float64(time_ns()-began_ns) / 1e9
    end
end

function targeted_supports_statistics(separator::TargetedPowerSupports)
    pointer = get(separator.optimizer.inner.sepas, separator, C_NULL)
    Dict("calls"=>separator.calls, "skipped"=>separator.skipped, "rounds"=>separator.rounds,
        "coordinates_available"=>length(separator.coordinates), "coordinate_checks"=>length(separator.checked),
        "support_vectors_certified"=>length(separator.certified), "cuts_added_to_global_pool"=>length(separator.added),
        "native_cuts_applied"=>pointer == C_NULL ? nothing : Int(SCIP.SCIPsepaGetNCutsApplied(pointer)),
        "certificate_evaluations"=>separator.certificate_evaluations[], "certificate_cache_hits"=>separator.certificate_hits[],
        "certificate_seconds"=>separator.certificate_seconds, "polynomial_certificate_seconds"=>separator.polynomial_seconds[],
        "certificate_budget_seconds"=>separator.certificate_budget,
        "certificate_budget_overrun_seconds"=>max(0.0, separator.certificate_seconds-separator.certificate_budget),
        "callback_seconds"=>separator.callback_seconds, "maximum_violation_mw"=>separator.maximum_violation_mw,
        "max_cuts"=>separator.max_cuts, "max_rounds"=>separator.max_rounds,
        "max_cuts_per_round"=>separator.max_cuts_per_round,
        "max_coordinates_per_round"=>separator.max_coordinates_per_round,
        "max_coordinate_checks"=>separator.max_coordinate_checks,
        "cuts"=>separator.cuts, "errors"=>copy(separator.errors),
        "infeasible_flags"=>separator.infeasible_flags,
        "cut_submission"=>"same global row added to pool and forced into current separation store; no parameter changes",
        "scope"=>"experimental root-only global cuts certified on original on-state boxes; LP values select rows only")
end
