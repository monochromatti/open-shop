using JuMP, OpenSHOP, HiGHS, SCIP

# Experimental additive envelope of the original power graph:
# P <= sum(a .* discharge_onweights) + sum(b .* head_onweights).
# Samples do not certify this inequality. Every original polynomial cell and
# linear extension is enclosed by its Bernstein coefficients. A second pass
# through the production support oracle guards the candidate LP's tolerances.
const _JOINT_POWER_LIMITS = (certificate_seconds=1.5, cuts=64, rounds=4,
    coordinates_per_round=16, coordinate_checks=64, violation_mw=1e-5)

struct _JointPowerBudgetExpired <: Exception end
_joint_power_budget(deadline)=time_ns()/1e9>=deadline ? throw(_JointPowerBudgetExpired()) : nothing

mutable struct JointPowerOracle
    table::TurbineTable
    qbox::NTuple{2,Float64}
    hbox::NTuple{2,Float64}
    qnodes::Vector{Float64}
    hnodes::Vector{Float64}
    electrical_max::Float64
    model::Model
    qcoefficients::Vector{VariableRef}
    hcoefficients::Vector{VariableRef}
    cells::Vector{NTuple{4,Float64}}
    coefficient_bound::Float64
    bernstein_rows::Int
    construction_seconds::Float64
    lp_seconds::Float64
    certificate_seconds::Float64
    solves::Int
    certificate_evaluations::Int
end

function _joint_power_nodes(nodes,box)
    values=Float64.(nodes)
    !isempty(values)&&all(isfinite,values)&&all(diff(values).>0)&&
        first(values)<=box[1]<=box[2]<=last(values) ||
        throw(ArgumentError("joint power coordinates must cover the original on-state box"))
    values
end

function _joint_power_segments(table_nodes,coordinate_nodes,box)
    nodes=sort!(unique(vcat(box[1],box[2],
        [x for x in table_nodes if box[1]<x<box[2]],
        [x for x in coordinate_nodes if box[1]<x<box[2]])))
    length(nodes)==1 ? [(only(nodes),only(nodes))] : collect(zip(nodes[1:end-1],nodes[2:end]))
end

# Same degree-(4,2) polynomial and normalized power-to-Bernstein transform as
# _power_box_upper. The added guard also matches its floating-point enclosure.
function _joint_power_bernstein(table,qlo,qhi,hlo,hhi,em)
    a=OpenSHOP._global_tensor_polynomial(table,qlo,qhi,hlo)
    z=OpenSHOP._global_tensor_polynomial(table,qlo,qhi,hhi)
    d=z.-a
    dq=qhi-qlo;dh=hhi-hlo
    C=zeros(5,3)
    for k in 1:4, (row,flow) in ((k,qlo),(k+1,dq))
        C[row,1]+=flow*hlo*a[k]
        C[row,2]+=flow*(hlo*d[k]+dh*a[k])
        C[row,3]+=flow*dh*d[k]
    end
    C .*= .00981*em
    guard=1e-9*max(1.0,sum(abs,C))
    B=zeros(5,3)
    for i in 0:4,j in 0:2
        B[i+1,j+1]=sum(C[k+1,l+1]*binomial(i,k)/binomial(4,k)*
            binomial(j,l)/binomial(2,l) for k in 0:i,l in 0:j)+guard
    end
    all(isfinite,B) || error("nonfinite joint power polynomial certificate")
    B
end

function JointPowerOracle(table,qbox,hbox,qnodes,hnodes,em;coefficient_bound=nothing,deadline=Inf)
    began=time_ns()
    qbox=Float64.((qbox[1],qbox[2]));hbox=Float64.((hbox[1],hbox[2]))
    all(isfinite,(qbox...,hbox...,em))&&qbox[1]<=qbox[2]&&hbox[1]<=hbox[2]&&em>=0 ||
        throw(ArgumentError("invalid original joint power domain"))
    qnodes=_joint_power_nodes(qnodes,qbox);hnodes=_joint_power_nodes(hnodes,hbox)
    qsegments=_joint_power_segments(table.discharge,qnodes,qbox)
    hsegments=_joint_power_segments(table.heads,hnodes,hbox)
    cells=NTuple{4,Float64}[(a,b,h0,h1) for (a,b) in qsegments for (h0,h1) in hsegments]
    polynomial_cells=Tuple{NTuple{4,Float64},Matrix{Float64}}[]
    for (qlo,qhi,hlo,hhi) in cells,k in 0:3
        _joint_power_budget(deadline)
        a=qlo+(qhi-qlo)*k/4;b=qlo+(qhi-qlo)*(k+1)/4
        push!(polynomial_cells,((a,b,hlo,hhi),_joint_power_bernstein(table,a,b,hlo,hhi,em)))
    end
    magnitude=maximum(x->maximum(abs,x[2]),polynomial_cells)
    # Bounds remove gauge/unused-node recession directions. They only restrict
    # which valid envelope we find; they never relax the certificate. In
    # particular, original off-domain coordinate nodes remain present.
    bound=coefficient_bound===nothing ? 16*max(1.0,magnitude) : Float64(coefficient_bound)
    isfinite(bound)&&bound>magnitude || throw(ArgumentError("joint coefficient bound must exceed the power enclosure"))
    m=Model(HiGHS.Optimizer);set_silent(m)
    # Do not change HiGHS' process-wide scheduler after the dispatch proposal
    # or another model has used it. The small LP itself uses serial simplex.
    set_optimizer_attribute(m,"solver","simplex")
    set_optimizer_attribute(m,"parallel","off")
    set_optimizer_attribute(m,"primal_feasibility_tolerance",1e-9)
    set_optimizer_attribute(m,"dual_feasibility_tolerance",1e-9)
    a=@variable(m,[1:length(qnodes)],lower_bound=-bound,upper_bound=bound)
    b=@variable(m,[1:length(hnodes)],lower_bound=-bound,upper_bound=bound)
    # Adding a constant to one axis and subtracting it from the other leaves
    # all physical envelopes unchanged because both on-weight sums equal u.
    gauge=clamp(searchsortedlast(qnodes,qbox[1]),1,length(qnodes))
    fix(a[gauge],0.0;force=true)
    for ((qlo,qhi,hlo,hhi),B) in polynomial_cells,i in 0:4,j in 0:2
        _joint_power_budget(deadline)
        # These are Bernstein control coordinates of an affine envelope, not
        # sampled values of power. Affine polynomials have these coefficients.
        qw=OpenSHOP._global_tensor_coordinate_weights(qnodes,qlo+(qhi-qlo)*i/4)
        hw=OpenSHOP._global_tensor_coordinate_weights(hnodes,hlo+(hhi-hlo)*j/2)
        @constraint(m,sum(qw[k]*a[k] for k in eachindex(qw))+
            sum(hw[k]*b[k] for k in eachindex(hw))>=max(0.0,B[i+1,j+1]))
    end
    JointPowerOracle(table,qbox,hbox,qnodes,hnodes,Float64(em),m,a,b,cells,bound,
        15*length(polynomial_cells),(time_ns()-began)/1e9,0.0,0.0,0,0)
end

function _joint_power_chord(nodes,coefficients,lo,hi)
    left=sum(coefficients.*OpenSHOP._global_tensor_coordinate_weights(nodes,lo))
    right=sum(coefficients.*OpenSHOP._global_tensor_coordinate_weights(nodes,hi))
    slope=lo==hi ? 0.0 : (right-left)/(hi-lo)
    slope,left-slope*lo
end

"""Certify an additive nodal envelope over the complete original on-state box.

Returns the guarded maximum residual; <=0 proves domination. The certificate
uses original PCHIP/bilinear polynomials and original linear extensions. This
is a finite-precision polynomial enclosure, with the production roundoff guard.
"""
function joint_power_certificate(oracle::JointPowerOracle,qcoefficients,hcoefficients;deadline=Inf)
    length(qcoefficients)==length(oracle.qnodes)&&length(hcoefficients)==length(oracle.hnodes)&&
        all(isfinite,qcoefficients)&&all(isfinite,hcoefficients) ||
        throw(ArgumentError("invalid joint power coefficients"))
    began=time_ns();upper=-Inf
    try
        for (qlo,qhi,hlo,hhi) in oracle.cells
            time_ns()/1e9>=deadline && return nothing
            a,qa=_joint_power_chord(oracle.qnodes,qcoefficients,qlo,qhi)
            b,hb=_joint_power_chord(oracle.hnodes,hcoefficients,hlo,hhi)
            residual=OpenSHOP._power_support(oracle.table,qlo,qhi,hlo,hhi,
                oracle.electrical_max,a,b)-qa-hb
            # For electrical efficiency e in [0,em], exact hydraulic power
            # obeys e*F <= max(0,em*F), also when an extension makes F negative.
            # A separable affine envelope attains its minimum at cell corners.
            zero_residual=-min(a*qlo,a*qhi)-min(b*hlo,b*hhi)-qa-hb
            isfinite(residual) || error("nonfinite joint power residual certificate")
            upper=max(upper,residual,zero_residual);oracle.certificate_evaluations+=1
        end
        upper
    finally
        oracle.certificate_seconds+=(time_ns()-began)/1e9
    end
end

"""Optimize one joint power envelope at existing root on-state weights.

The finite coefficient box affects strength only. HiGHS solves a sufficient
Bernstein LP, then the independent production residual oracle certifies and
guards its answer. Returns nothing on allowance expiry or no LP incumbent.
"""
function joint_power_support(oracle::JointPowerOracle,qweights,hweights;time_limit=Inf)
    length(qweights)==length(oracle.qnodes)&&length(hweights)==length(oracle.hnodes)&&
        all(isfinite,qweights)&&all(isfinite,hweights) || throw(ArgumentError("invalid joint root weights"))
    (isfinite(time_limit)&&time_limit>=0)||time_limit==Inf || throw(ArgumentError("invalid joint allowance"))
    time_limit==0 && return nothing
    deadline=time_ns()/1e9+time_limit
    m=oracle.model
    @objective(m,Min,sum(qweights[j]*oracle.qcoefficients[j] for j in eachindex(qweights))+
        sum(hweights[j]*oracle.hcoefficients[j] for j in eachindex(hweights)))
    isfinite(time_limit) ? set_time_limit_sec(m,max(1e-6,deadline-time_ns()/1e9)) : unset_time_limit_sec(m)
    began=time_ns();optimize!(m)
    oracle.lp_seconds+=(time_ns()-began)/1e9;oracle.solves+=1
    if !has_values(m)||primal_status(m)!=JuMP.MOI.FEASIBLE_POINT
        termination_status(m) in (JuMP.MOI.TIME_LIMIT,JuMP.MOI.ITERATION_LIMIT,JuMP.MOI.INTERRUPTED) && return nothing
        error("joint coefficient LP has no feasible candidate: $(termination_status(m)), $(raw_status(m))")
    end
    qa=value.(oracle.qcoefficients);hb=value.(oracle.hcoefficients)
    upper=joint_power_certificate(oracle,qa,hb;deadline)
    upper===nothing && return nothing
    # Uniform discharge correction multiplies sum(nu)=u. It is therefore
    # valid for signed intercepts and vanishes exactly when the unit is off.
    guard=max(0.0,upper)+1e-9*max(1.0,maximum(abs,qa),maximum(abs,hb))
    qa .+= guard
    (;qcoefficients=qa,hcoefficients=hb,certified=true,residual_upper=upper-guard,
        correction=guard,lp_status=string(termination_status(m)),
        value=sum(qa.*qweights)+sum(hb.*hweights),
        coefficient_bound=oracle.coefficient_bound,
        active_coefficient_bounds=count(x->abs(x)>=oracle.coefficient_bound-1e-7,vcat(qa.-guard,hb)))
end

function joint_power_support(table,qbox,hbox,qnodes,hnodes,em,qweights,hweights;kwargs...)
    joint_power_support(JointPowerOracle(table,qbox,hbox,qnodes,hnodes,em),qweights,hweights;kwargs...)
end

struct _JointNativeCoordinate
    discharge::OpenSHOP._NativePowerCoordinate
    head::OpenSHOP._NativePowerCoordinate
end

function _joint_power_row(coordinate,qa,hb)
    terms=Dict{Int,Float64}();constant=coordinate.head.power.constant
    power=coordinate.head.power
    for j in eachindex(power.columns)
        terms[power.columns[j]]=get(terms,power.columns[j],0.0)+power.coefficients[j]
    end
    for (weights,coefficients) in ((coordinate.discharge.onweights,qa),(coordinate.head.onweights,hb)),j in eachindex(weights)
        expression=weights[j];coefficient=coefficients[j]
        constant-=coefficient*expression.constant
        for k in eachindex(expression.columns)
            column=expression.columns[k]
            terms[column]=get(terms,column,0.0)-coefficient*expression.coefficients[k]
        end
    end
    columns=sort!([j for (j,v) in terms if v!=0.0])
    row=OpenSHOP._PowerCutAffine(constant/40,columns,[terms[j]/40 for j in columns])
    isfinite(row.constant)&&all(isfinite,row.coefficients) || error("nonfinite certified joint power row")
    row
end

mutable struct JointPowerSeparator <: SCIP.AbstractSeparator
    optimizer::SCIP.Optimizer
    coordinates::Vector{_JointNativeCoordinate}
    references::Vector{SCIP.VarRef}
    pointers::Vector{Ptr{SCIP.SCIP_VAR}}
    values::Vector{Float64}
    visits::Vector{Int}
    oracles::IdDict{TurbineTable,Dict{Any,JointPowerOracle}}
    limits::typeof(_JOINT_POWER_LIMITS)
    calls::Int
    skipped::Int
    rounds::Int
    coordinate_checks::Int
    cuts_added::Int
    work_seconds::Float64
    callback_seconds::Float64
    maximum_violation_mw::Float64
    infeasible_flags::Int
    timeouts::Int
    cache_hits::Int
    errors::Vector{String}
    certified::Vector{Any}
end

function install_joint_power_supports!(b,c;certificate_budget=_JOINT_POWER_LIMITS.certificate_seconds,
        max_cuts=_JOINT_POWER_LIMITS.cuts,max_rounds=_JOINT_POWER_LIMITS.rounds,
        max_coordinates_per_round=_JOINT_POWER_LIMITS.coordinates_per_round,
        max_coordinate_checks=_JOINT_POWER_LIMITS.coordinate_checks,
        violation_tolerance=_JOINT_POWER_LIMITS.violation_mw)
    sources=get(b.m.ext,:global_table_power_supports,OpenSHOP._TablePowerSupportCoordinate[])
    isempty(sources) && return nothing
    isfinite(certificate_budget)&&certificate_budget>=0&&isfinite(violation_tolerance)&&violation_tolerance>=0 ||
        throw(ArgumentError("invalid joint cut allowance"))
    all(x->x isa Integer&&x>=0,(max_cuts,max_rounds,max_coordinates_per_round,max_coordinate_checks)) ||
        throw(ArgumentError("invalid joint cut work limits"))
    optimizer=unsafe_backend(b.m)
    optimizer isa SCIP.Optimizer || throw(ArgumentError("attached SCIP optimizer required"))
    columns=Dict{VariableRef,Int}();references=SCIP.VarRef[];pointers=Ptr{SCIP.SCIP_VAR}[]
    affine(x)=OpenSHOP._power_cut_affine(x,columns,references,pointers,optimizer)
    native(source)=OpenSHOP._NativePowerCoordinate(source,affine.(source.onweights),affine(source.term),
        affine(source.power),affine(source.commitment),affine(source.discharge),affine(source.onhead),
        max(0.0,c.prices[source.interval])*(c.grid[source.interval+1]-c.grid[source.interval]))
    paired=Dict{Tuple{Int,Int},Dict{Symbol,OpenSHOP._TablePowerSupportCoordinate}}()
    for source in sources
        get!(paired,(source.unit,source.interval),Dict{Symbol,OpenSHOP._TablePowerSupportCoordinate}())[source.axis]=source
    end
    coordinates=_JointNativeCoordinate[]
    for key in sort!(collect(keys(paired)))
        pair=paired[key]
        haskey(pair,:head)&&haskey(pair,:discharge) || continue
        q=pair[:discharge];h=pair[:head]
        q.table===h.table&&q.qbox==h.qbox&&q.hbox==h.hbox&&q.electrical_max==h.electrical_max ||
            error("joint source domains differ")
        push!(coordinates,_JointNativeCoordinate(native(q),native(h)))
    end
    isempty(coordinates) && return nothing
    limits=(certificate_seconds=Float64(certificate_budget),cuts=Int(max_cuts),rounds=Int(max_rounds),
        coordinates_per_round=Int(max_coordinates_per_round),coordinate_checks=Int(max_coordinate_checks),
        violation_mw=Float64(violation_tolerance))
    separator=JointPowerSeparator(optimizer,coordinates,references,pointers,
        Vector{Float64}(undef,length(references)),zeros(Int,length(references)),
        IdDict{TurbineTable,Dict{Any,JointPowerOracle}}(),limits,
        0,0,0,0,0,0.0,0.0,0.0,0,0,0,String[],Any[])
    SCIP.include_sepa(optimizer,separator;name="joint_table_power",description="Joint certified table-power supports",
        priority=99999,freq=0,maxbounddist=1.0,delay=false)
    separator
end

function _joint_power_oracle!(separator,coordinate;deadline=Inf)
    q=coordinate.discharge.source;h=coordinate.head.source
    entries=get!(separator.oracles,q.table) do
        Dict{Any,JointPowerOracle}()
    end
    key=(q.qbox,q.hbox,Tuple(q.nodes),Tuple(h.nodes),q.electrical_max)
    if haskey(entries,key)
        separator.cache_hits+=1
        return entries[key]
    end
    oracle=JointPowerOracle(q.table,q.qbox,q.hbox,q.nodes,h.nodes,q.electrical_max;deadline)
    entries[key]=oracle
end

function SCIP.exec_lp(separator::JointPowerSeparator)
    began=time_ns();separator.calls+=1
    try
        optimizer=separator.optimizer;limits=separator.limits
        node=SCIP.SCIPgetFocusNode(optimizer)
        if node==C_NULL||SCIP.SCIPnodeGetDepth(node)!=0||Bool(SCIP.SCIPinProbing(optimizer))||
                Bool(SCIP.SCIPinDive(optimizer))||SCIP.SCIPgetLPSolstat(optimizer)!=SCIP.SCIP_LPSOLSTAT_OPTIMAL||
                !Bool(SCIP.SCIPisLPPrimalReliable(optimizer))||!Bool(SCIP.SCIPisLPRelax(optimizer))||
                separator.rounds>=limits.rounds||separator.cuts_added>=limits.cuts||
                separator.coordinate_checks>=limits.coordinate_checks||separator.work_seconds>=limits.certificate_seconds||
                !isempty(separator.errors)
            separator.skipped+=1;return SCIP.SCIP_DIDNOTRUN
        end
        separator.rounds+=1
        at=function(column)
            if separator.visits[column]!=separator.calls
                result=Float64(SCIP.SCIPgetSolVal(optimizer,C_NULL,separator.pointers[column]))
                isfinite(result) || error("nonfinite joint root value")
                separator.values[column]=result;separator.visits[column]=separator.calls
            end
            separator.values[column]
        end
        ranking=Tuple{Float64,Int}[]
        for j in eachindex(separator.coordinates)
            score=OpenSHOP._power_cut_score(separator.coordinates[j].head,at)
            score>0 && push!(ranking,(score,j))
        end
        sort!(ranking;by=x->(-x[1],x[2]))
        count=min(limits.coordinates_per_round,limits.coordinate_checks-separator.coordinate_checks,
            limits.cuts-separator.cuts_added,length(ranking))
        added=0
        for j in 1:count
            separator.work_seconds>=limits.certificate_seconds && break
            index=ranking[j][2];coordinate=separator.coordinates[index]
            separator.coordinate_checks+=1;work_began=time_ns()
            candidate=nothing
            try
                deadline=time_ns()/1e9+limits.certificate_seconds-separator.work_seconds
                oracle=_joint_power_oracle!(separator,coordinate;deadline)
                remaining=limits.certificate_seconds-separator.work_seconds-(time_ns()-work_began)/1e9
                candidate=joint_power_support(oracle,
                    [OpenSHOP._power_cut_value(x,at) for x in coordinate.discharge.onweights],
                    [OpenSHOP._power_cut_value(x,at) for x in coordinate.head.onweights];time_limit=max(0.0,remaining))
            catch exception
                exception isa _JointPowerBudgetExpired || rethrow()
            finally
                separator.work_seconds+=(time_ns()-work_began)/1e9
            end
            if candidate===nothing
                separator.timeouts+=1;continue
            end
            expression=_joint_power_row(coordinate,candidate.qcoefficients,candidate.hcoefficients)
            violation=40*OpenSHOP._power_cut_value(expression,at)
            violation>limits.violation_mw || continue
            infeasible=OpenSHOP._power_cut_add_row!(separator,expression)
            push!(separator.certified,(coordinate=index,candidate=candidate,expression=expression))
            separator.cuts_added+=1;added+=1
            separator.maximum_violation_mw=max(separator.maximum_violation_mw,violation)
            if infeasible
                separator.infeasible_flags+=1;error("certified joint cut reported root infeasibility")
            end
        end
        added>0 ? SCIP.SCIP_SEPARATED : SCIP.SCIP_DIDNOTFIND
    catch exception
        push!(separator.errors,sprint(showerror,exception));SCIP.SCIP_DIDNOTFIND
    finally
        separator.callback_seconds+=(time_ns()-began)/1e9
    end
end

function joint_power_supports_statistics(separator::JointPowerSeparator)
    pointer=get(separator.optimizer.inner.sepas,separator,C_NULL)
    oracles=[oracle for entries in values(separator.oracles) for oracle in values(entries)]
    Dict("calls"=>separator.calls,"skipped"=>separator.skipped,"rounds"=>separator.rounds,
        "coordinates_available"=>length(separator.coordinates),"coordinate_checks"=>separator.coordinate_checks,
        "cuts_added_to_global_pool"=>separator.cuts_added,"native_cuts_applied"=>pointer==C_NULL ? nothing : Int(SCIP.SCIPsepaGetNCutsApplied(pointer)),
        "oracle_models"=>length(oracles),"oracle_cache_hits"=>separator.cache_hits,
        "bernstein_lp_rows"=>sum(x->x.bernstein_rows,oracles;init=0),
        "oracle_lp_solves"=>sum(x->x.solves,oracles;init=0),
        "certificate_evaluations"=>sum(x->x.certificate_evaluations,oracles;init=0),
        "construction_seconds"=>sum(x->x.construction_seconds,oracles;init=0.0),
        "lp_seconds"=>sum(x->x.lp_seconds,oracles;init=0.0),
        "certificate_seconds"=>sum(x->x.certificate_seconds,oracles;init=0.0),
        "work_seconds"=>separator.work_seconds,"work_budget_seconds"=>separator.limits.certificate_seconds,
        "work_budget_overrun_seconds"=>max(0.0,separator.work_seconds-separator.limits.certificate_seconds),
        "callback_seconds"=>separator.callback_seconds,"maximum_violation_mw"=>separator.maximum_violation_mw,
        "max_cuts"=>separator.limits.cuts,"max_rounds"=>separator.limits.rounds,
        "max_coordinate_checks"=>separator.limits.coordinate_checks,"timeouts"=>separator.timeouts,
        "infeasible_flags"=>separator.infeasible_flags,"errors"=>copy(separator.errors))
end
