# LP values select extra rows; polynomial certificates justify each row over
# the original on-state box. Neither probing bounds nor sampled power justify
# a cut. The exact turbine graphs and nonlinear power equations remain intact.
const _POWER_CUT_FRACTIONS = (0.125,0.375,0.625,0.875,1.125)
const _POWER_CUT_LIMITS = (certificate_seconds=1.5,cuts=64,rounds=4,
    cuts_per_round=16,coordinates_per_round=16,coordinate_checks=64,violation_mw=1e-5)

struct _PowerCutAffine
    constant::Float64
    columns::Vector{Int}
    coefficients::Vector{Float64}
end

function _power_cut_column!(variable,columns,references,pointers,optimizer)
    get!(columns,variable) do
        reference=SCIP.VarRef(optimizer_index(variable).value)
        pointer=optimizer.inner.vars[reference][]
        pointer!=C_NULL || error("missing native power-cut variable")
        push!(references,reference)
        push!(pointers,pointer)
        length(references)
    end
end

function _power_cut_affine(expression,columns,references,pointers,optimizer)
    expression isa Number && return _PowerCutAffine(Float64(expression),Int[],Float64[])
    if expression isa VariableRef
        column=_power_cut_column!(expression,columns,references,pointers,optimizer)
        return _PowerCutAffine(0.0,[column],[1.0])
    end
    expression isa AffExpr || throw(ArgumentError("power cuts require affine existing coordinates"))
    terms=Dict{Int,Float64}()
    for (coefficient,variable) in linear_terms(expression)
        column=_power_cut_column!(variable,columns,references,pointers,optimizer)
        terms[column]=get(terms,column,0.0)+coefficient
    end
    indices=sort!([j for (j,coefficient) in terms if coefficient!=0.0])
    _PowerCutAffine(Float64(expression.constant),indices,[terms[j] for j in indices])
end

_power_cut_value(expression::_PowerCutAffine,at)=expression.constant+
    sum(expression.coefficients[j]*at(expression.columns[j]) for j in eachindex(expression.columns);init=0.0)

struct _NativePowerCoordinate
    source::_TablePowerSupportCoordinate
    onweights::Vector{_PowerCutAffine}
    term::_PowerCutAffine
    power::_PowerCutAffine
    commitment::_PowerCutAffine
    discharge::_PowerCutAffine
    onhead::_PowerCutAffine
    objective_weight::Float64
end

struct _CertifiedPowerCut
    coordinate::Int
    fraction::Float64
    expression::_PowerCutAffine
end

# Discharge's first on-weight includes a constant (-1). Move every affine
# constant to the native row's RHS, including signed table intercepts.
function _power_cut_row(coordinate,slope,coefficients)
    terms=Dict{Int,Float64}()
    constant=0.0
    for (multiplier,expression) in ((1.0,coordinate.power),(-slope,coordinate.term))
        constant+=multiplier*expression.constant
        for j in eachindex(expression.columns)
            column=expression.columns[j]
            terms[column]=get(terms,column,0.0)+multiplier*expression.coefficients[j]
        end
    end
    for j in eachindex(coefficients)
        expression=coordinate.onweights[j]
        constant-=coefficients[j]*expression.constant
        for k in eachindex(expression.columns)
            column=expression.columns[k]
            terms[column]=get(terms,column,0.0)-coefficients[j]*expression.coefficients[k]
        end
    end
    columns=sort!([j for (j,coefficient) in terms if coefficient!=0.0])
    row=_PowerCutAffine(constant/40.0,columns,[terms[j]/40.0 for j in columns])
    isfinite(row.constant)&&all(isfinite,row.coefficients) || error("nonfinite certified power cut")
    row
end

mutable struct _PowerCutSeparator{F} <: SCIP.AbstractSeparator
    optimizer::SCIP.Optimizer
    coordinates::Vector{_NativePowerCoordinate}
    references::Vector{SCIP.VarRef}
    pointers::Vector{Ptr{SCIP.SCIP_VAR}}
    values::Vector{Float64}
    visits::Vector{Int}
    support::F
    certificate_evaluations::Base.RefValue{Int}
    certificate_hits::Base.RefValue{Int}
    polynomial_seconds::Base.RefValue{Float64}
    certified::Vector{_CertifiedPowerCut}
    checked::BitVector
    added::BitVector
    ranking::Vector{Tuple{Float64,Int}}
    candidates::Vector{Tuple{Float64,Float64,Int}}
    limits::typeof(_POWER_CUT_LIMITS)
    calls::Int
    skipped::Int
    rounds::Int
    coordinate_checks::Int
    cuts_added::Int
    certificate_seconds::Float64
    callback_seconds::Float64
    maximum_violation_mw::Float64
    infeasible_flags::Int
    errors::Vector{String}
end

function _install_power_cuts!(b,c;certificate_budget=_POWER_CUT_LIMITS.certificate_seconds,
        max_cuts=_POWER_CUT_LIMITS.cuts,max_rounds=_POWER_CUT_LIMITS.rounds,
        max_cuts_per_round=_POWER_CUT_LIMITS.cuts_per_round,
        max_coordinates_per_round=_POWER_CUT_LIMITS.coordinates_per_round,
        max_coordinate_checks=_POWER_CUT_LIMITS.coordinate_checks,
        violation_tolerance=_POWER_CUT_LIMITS.violation_mw)
    sources=get(b.m.ext,:global_table_power_supports,_TablePowerSupportCoordinate[])::Vector{_TablePowerSupportCoordinate}
    isempty(sources) && return nothing
    isfinite(certificate_budget)&&certificate_budget>=0 || throw(ArgumentError("invalid power-cut certificate allowance"))
    isfinite(violation_tolerance)&&violation_tolerance>=0 || throw(ArgumentError("invalid power-cut violation tolerance"))
    all(x->x isa Integer&&x>=0,(max_cuts,max_rounds,max_cuts_per_round,
        max_coordinates_per_round,max_coordinate_checks)) || throw(ArgumentError("invalid power-cut work limits"))
    optimizer=unsafe_backend(b.m)
    optimizer isa SCIP.Optimizer || throw(ArgumentError("attached SCIP optimizer required"))
    columns=Dict{VariableRef,Int}()
    references=SCIP.VarRef[]
    pointers=Ptr{SCIP.SCIP_VAR}[]
    coordinates=_NativePowerCoordinate[]
    affine(x)=_power_cut_affine(x,columns,references,pointers,optimizer)
    for source in sources
        weight=max(0.0,c.prices[source.interval])*(c.grid[source.interval+1]-c.grid[source.interval])
        push!(coordinates,_NativePowerCoordinate(source,affine.(source.onweights),affine(source.term),
            affine(source.power),affine(source.commitment),affine(source.discharge),affine(source.onhead),weight))
    end
    support,evaluations,hits,seconds=_table_power_support_cache()
    limits=(certificate_seconds=Float64(certificate_budget),cuts=Int(max_cuts),rounds=Int(max_rounds),
        cuts_per_round=Int(max_cuts_per_round),coordinates_per_round=Int(max_coordinates_per_round),
        coordinate_checks=Int(max_coordinate_checks),violation_mw=Float64(violation_tolerance))
    separator=_PowerCutSeparator(optimizer,coordinates,references,pointers,
        Vector{Float64}(undef,length(references)),zeros(Int,length(references)),support,evaluations,hits,seconds,
        _CertifiedPowerCut[],falses(length(coordinates)),BitVector(),Tuple{Float64,Int}[],
        Tuple{Float64,Float64,Int}[],limits,0,0,0,0,0,0.0,0.0,0.0,0,String[])
    SCIP.include_sepa(optimizer,separator;name="table_power",description="Certified table-power supports",
        priority=100000,freq=0,maxbounddist=1.0,delay=false)
    separator
end

function _power_cut_add_row!(separator,expression)
    optimizer=separator.optimizer
    row=Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(optimizer,row,optimizer.inner.sepas[separator],
        "",-SCIP.SCIPinfinity(optimizer),-expression.constant,false,false,false)
    try
        variables=separator.pointers[expression.columns]
        SCIP.@SCIP_CALL SCIP.SCIPaddVarsToRow(optimizer,row[],length(variables),variables,expression.coefficients)
        # The global pool alone need not process a new cut in this LP round.
        # Submit the same globally valid row to the current separation store.
        SCIP.@SCIP_CALL SCIP.SCIPaddPoolCut(optimizer,row[])
        infeasible=Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPaddRow(optimizer,row[],true,infeasible)
        Bool(infeasible[])
    finally
        SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(optimizer,row)
    end
end

function _power_cut_score(coordinate,at)
    coordinate.objective_weight<=0 && return 0.0
    source=coordinate.source
    u=clamp(_power_cut_value(coordinate.commitment,at),0.0,1.0)
    p=_power_cut_value(coordinate.power,at)
    u<=1e-8 && return coordinate.objective_weight*max(0.0,p)
    q=clamp(_power_cut_value(coordinate.discharge,at)/u,source.qbox...)
    h=clamp(_power_cut_value(coordinate.onhead,at)/u,source.hbox...)
    physical=u*0.00981*source.electrical_max*q*h*turbine_efficiency(source.table,q,h;extrapolation=:linear)
    isfinite(physical) || error("nonfinite power-cut selection value")
    # This perspective residual only ranks work. It supplies no certificate.
    coordinate.objective_weight*(max(0.0,p-physical)+1e-10*max(0.0,p))
end

function _power_cut_certify!(separator,index)
    coordinate=separator.coordinates[index]
    source=coordinate.source
    separator.checked[index]=true
    separator.coordinate_checks+=1
    scale=0.00981*source.electrical_max*(source.axis==:head ? source.hbox[2] : source.qbox[2])
    for fraction in _POWER_CUT_FRACTIONS
        separator.certificate_seconds>=separator.limits.certificate_seconds && break
        began=time_ns()
        slope=scale*fraction
        coefficients=source.axis==:head ?
            _table_power_head_coefficients(source.table,source.qbox,source.hbox,source.nodes,
                source.electrical_max,slope,separator.support) :
            _table_power_discharge_coefficients(source.table,source.qbox,source.hbox,source.nodes,
                source.electrical_max,slope,separator.support)
        all(isfinite,coefficients)&&isfinite(slope) || error("nonfinite power-cut certificate")
        push!(separator.certified,_CertifiedPowerCut(index,fraction,_power_cut_row(coordinate,slope,coefficients)))
        push!(separator.added,false)
        separator.certificate_seconds+=Float64(time_ns()-began)/1e9
    end
end

function SCIP.exec_lp(separator::_PowerCutSeparator)
    began=time_ns()
    separator.calls+=1
    try
        optimizer=separator.optimizer
        limits=separator.limits
        node=SCIP.SCIPgetFocusNode(optimizer)
        if node==C_NULL || SCIP.SCIPnodeGetDepth(node)!=0 ||
                Bool(SCIP.SCIPinProbing(optimizer)) || Bool(SCIP.SCIPinDive(optimizer)) ||
                SCIP.SCIPgetLPSolstat(optimizer)!=SCIP.SCIP_LPSOLSTAT_OPTIMAL ||
                !Bool(SCIP.SCIPisLPPrimalReliable(optimizer)) || !Bool(SCIP.SCIPisLPRelax(optimizer)) ||
                separator.rounds>=limits.rounds || separator.cuts_added>=limits.cuts || !isempty(separator.errors)
            separator.skipped+=1
            return SCIP.SCIP_DIDNOTRUN
        end
        separator.rounds+=1
        at=function(column)
            if separator.visits[column]!=separator.calls
                result=Float64(SCIP.SCIPgetSolVal(optimizer,C_NULL,separator.pointers[column]))
                isfinite(result) || error("nonfinite root LP power-cut value")
                separator.values[column]=result
                separator.visits[column]=separator.calls
            end
            separator.values[column]
        end
        if separator.certificate_seconds<limits.certificate_seconds && separator.coordinate_checks<limits.coordinate_checks
            empty!(separator.ranking)
            for j in eachindex(separator.coordinates)
                separator.checked[j] && continue
                score=_power_cut_score(separator.coordinates[j],at)
                score>0.0 && push!(separator.ranking,(score,j))
            end
            sort!(separator.ranking;by=x->(-x[1],x[2]))
            count=min(limits.coordinates_per_round,limits.coordinate_checks-separator.coordinate_checks,length(separator.ranking))
            for j in 1:count
                separator.certificate_seconds>=limits.certificate_seconds && break
                _power_cut_certify!(separator,separator.ranking[j][2])
            end
        end
        empty!(separator.candidates)
        for (index,cut) in enumerate(separator.certified)
            separator.added[index] && continue
            violation=40.0*_power_cut_value(cut.expression,at)
            violation>limits.violation_mw || continue
            weight=separator.coordinates[cut.coordinate].objective_weight
            push!(separator.candidates,(weight*violation,violation,index))
        end
        sort!(separator.candidates;by=x->(-x[1],-x[2],x[3]))
        count=min(limits.cuts_per_round,limits.cuts-separator.cuts_added,length(separator.candidates))
        for j in 1:count
            _,violation,index=separator.candidates[j]
            infeasible=_power_cut_add_row!(separator,separator.certified[index].expression)
            separator.added[index]=true
            separator.cuts_added+=1
            separator.maximum_violation_mw=max(separator.maximum_violation_mw,violation)
            if infeasible
                separator.infeasible_flags+=1
                error("certified global power cut reported root infeasibility")
            end
        end
        count>0 ? SCIP.SCIP_SEPARATED : SCIP.SCIP_DIDNOTFIND
    catch exception
        push!(separator.errors,sprint(showerror,exception))
        SCIP.SCIP_DIDNOTFIND
    finally
        separator.callback_seconds+=Float64(time_ns()-began)/1e9
    end
end

function _power_cut_statistics(separator::_PowerCutSeparator)
    pointer=get(separator.optimizer.inner.sepas,separator,C_NULL)
    limits=separator.limits
    Dict("calls"=>separator.calls,"skipped"=>separator.skipped,"rounds"=>separator.rounds,
        "coordinates_available"=>length(separator.coordinates),"coordinate_checks"=>separator.coordinate_checks,
        "support_vectors_certified"=>length(separator.certified),"cuts_added_to_global_pool"=>separator.cuts_added,
        "native_cuts_applied"=>pointer==C_NULL ? nothing : Int(SCIP.SCIPsepaGetNCutsApplied(pointer)),
        "certificate_evaluations"=>separator.certificate_evaluations[],"certificate_cache_hits"=>separator.certificate_hits[],
        "certificate_seconds"=>separator.certificate_seconds,"polynomial_certificate_seconds"=>separator.polynomial_seconds[],
        "certificate_budget_seconds"=>limits.certificate_seconds,
        "certificate_budget_overrun_seconds"=>max(0.0,separator.certificate_seconds-limits.certificate_seconds),
        "callback_seconds"=>separator.callback_seconds,"maximum_violation_mw"=>separator.maximum_violation_mw,
        "max_cuts"=>limits.cuts,"max_rounds"=>limits.rounds,"max_cuts_per_round"=>limits.cuts_per_round,
        "max_coordinates_per_round"=>limits.coordinates_per_round,"max_coordinate_checks"=>limits.coordinate_checks,
        "infeasible_flags"=>separator.infeasible_flags,"errors"=>copy(separator.errors))
end

function _reject_power_cut_bound!(result)
    statistics=get(result,"power_cut_statistics",nothing)
    if get(result,"power_cut_statistics_error",nothing)!==nothing ||
            (statistics!==nothing && (!isempty(statistics["errors"]) || statistics["infeasible_flags"]>0))
        result["rejected_global_bound"]=result["global_bound"]
        result["global_bound"]=nothing
        result["global_certificate"]=false
        result["bound_rejection"]="certified power-cut callback or statistics failed; native bound is retained only in diagnostics"
    end
    result
end
