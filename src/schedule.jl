"""Schedule an input case using successive MILP proposals and nonlinear dispatch.

Each accepted nonlinear improvement supplies the next proposal's hydraulic
reference. Stop when discrete controls repeat, an accepted dispatch fails to
improve the best schedule, or the shared deadline expires. Cold dispatch restores
feasibility before optimizing revenue; a failed improvement retains that
feasible trajectory. Unsuccessful assignments are excluded from subsequent
MILP proposals, so another discrete choice can be tried. These exclusions
never restrict a later global search.

`time_limit` defaults to 60 seconds and is shared by construction, solves and
refinements. The MILP/NLP limits also bound individual solves. Budgets are
checked between phases; construction, a solver iteration or final replay can
overrun the allowance. All work is included in elapsed time. An accepted result is feasible under
the declared numerical checks, not a global optimality certificate. `operational_margin` optionally
requests a symmetric power inset during preparation, capped at half each
permitted range. Acceptance
always uses the original physical restrictions.
`initial` is an accepted schedule to retain. `warm` is an unchecked numerical
starting point after a forecast update; its controls are independently reconstructed
and checked before they can be retained. Its arrays must use this case's object order and a
hydraulic grid retaining the operating boundaries.
"""
function schedule_case(
    c::ScheduleCase;
    initial = nothing,
    warm = nothing,
    proposal_time_limit = 10.0,
    nlp_time_limit = 30.0,
    max_refinements = 48,
    operational_margin = 0.0,
    time_limit = 60.0,
    refine_grid = true,
)
    time_limit>0 && !isnan(time_limit) ||
        throw(ArgumentError("positive preparation time limit required"))
    isfinite(proposal_time_limit) && proposal_time_limit>0 ||
        throw(ArgumentError("positive finite proposal time limit required"))
    isfinite(nlp_time_limit) && nlp_time_limit>0 ||
        throw(ArgumentError("positive finite dispatch time limit required"))
    max_refinements isa Integer && max_refinements>=0 ||
        throw(ArgumentError("nonnegative integer refinement budget required"))
    isfinite(operational_margin) && operational_margin >= 0 ||
        throw(ArgumentError("nonnegative finite operational_margin required"))
    started=time()
    deadline=started+time_limit
    validate_inputs(c)
    attempts=Any[]
    best=nothing
    if initial!==nothing
        initial_schedule=_schedule_input(initial;case=c)
        mesh=TimeDiscretization(c.grid,initial_schedule.grid)
        size(initial_schedule.commitment)==(length(c.system.generators),length(mesh.dt)) ||
            throw(ArgumentError("initial schedule commitment dimensions differ from its grid"))
        u=initial_schedule.commitment[:,first.(mesh.cells)]
        v=_adaptive_dispatch(
            c;
            u,
            initial,
            tunnel_opening=initial_schedule.tunnel_opening[:,first.(mesh.cells)],
            max_refinements = 0,
            time_limit = nlp_time_limit,
            operational_margin,
            deadline,
            refine_grid,
        )
        v.accepted ||
            throw(ArgumentError("supplied initial schedule fails independent acceptance"))
        best=v
    end
    fallback=ones(Int, length(c.system.generators), length(c.prices))
    bounds=_commitment_bounds(c)
    fallback=clamp.(fallback,bounds.unit_lower,bounds.unit_upper)
    fallback_opening=tunnel_opening_schedule(c)
    fallback_gates_valid=try
        validate_tunnel_openings(c,fallback_opening,fallback)
    catch error
        error isa ArgumentError || rethrow()
        false
    end
    reference=best===nothing ? warm : best.solution
    proposals=Any[]
    control_type=Tuple{Matrix{Int},Matrix{Float64},Matrix{Int}}
    dispatched=Set{control_type}()
    excluded=control_type[]
    if best===nothing && warm!==nothing
        start=_schedule_input(warm;case=c)
        mesh=TimeDiscretization(c.grid,start.grid)
        size(start.commitment)==(size(fallback,1),length(mesh.dt)) &&
            all(x->x==0||x==1,start.commitment) ||
            throw(ArgumentError("warm commitment must be binary and unit by interval"))
        reference,reason=_dispatch_start(c,fallback,warm;grid=refine_grid ? start.grid : c.grid)
        reference!==nothing || throw(ArgumentError("incompatible warm schedule: "*
            get(reason,"ignored_reason","missing hydraulic arrays")))
    end
    while time()<deadline
        p=propose_commitment(c;excluded,time_limit=proposal_time_limit,reference,
            incumbent=best===nothing ? nothing : best.solution,deadline)
        push!(proposals,p)
        !haskey(p,"u") && best!==nothing && break
        # A missing MILP point is not evidence that the nonlinear problem is
        # infeasible. Use the ordinary admissible on-state guess in that case.
        u=get(p,"u",fallback)
        opening=get(p,"tunnel_opening",fallback_opening)
        river_on=get(p,"river_on",fixed_river_on(c,nothing,nothing))
        time()<deadline && admissible(c,u) &&
            (haskey(p,"u") || fallback_gates_valid) || break
        # Own the snapshot: later numerical starts must not change a set key.
        controls=(copy(u),copy(opening),copy(river_on))
        controls in dispatched && break
        push!(dispatched,controls)
        try
            problem=best===nothing ? c : best.case
            v=_adaptive_dispatch(problem;u,tunnel_opening=opening,river_on,
                warm=reference,max_refinements,time_limit=nlp_time_limit,
                operational_margin,deadline,refine_grid)
            push!(attempts,Dict("stage"=>"dispatch","accepted"=>v.accepted,
                "seconds"=>v.total_seconds,"steps"=>_dispatch_reports(v.attempts),"u"=>copy(u),
                "tunnel_opening"=>copy(opening),"elapsed_seconds"=>time()-started,
                "warm_start_objective"=>reference===nothing ? nothing : reference["objective"]))
            if !v.accepted
                push!(excluded,controls)
                continue
            end
            if best===nothing || v.solution.objective>best.solution.objective
                best=v
                reference=v.solution
            else
                break
            end
        catch error
            error isa Union{ArgumentError,DomainError} || rethrow()
            push!(attempts,Dict("accepted"=>false,"error"=>sprint(showerror,error)))
            break
        end
    end
    elapsed=time()-started
    SolveResult(Dict(
        "status"=>best===nothing ? "NO_ACCEPTED_SOLUTION" : "FEASIBLE",
        "accepted"=>best!==nothing,
        "problem"=>best===nothing ? c : best.case,
        "solution"=>best===nothing ? nothing : best.solution,
        "replay_validation"=>best===nothing ? nothing : _replay_report(best.replay_validation),
        "proposals"=>proposals,
        "attempts"=>attempts,
        "seconds"=>elapsed,
        "total_seconds"=>elapsed,
        "feasible_lower_bound"=>best===nothing ? nothing : best.solution.objective,
        "acceptance_scope"=>"physical equations, operating limits and finer replay",
        "global_optimality_proven"=>false,
        "requested_operational_margin_MW"=>operational_margin,
    ))
end
