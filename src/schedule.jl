"""Schedule an input case using MILP proposals and audited nonlinear dispatch.

`time_limit` caps total preparation, including construction and refinements;
the MILP/NLP limits also cap individual solves. Final replay may overrun the
allowance and is included in elapsed time. An accepted result is feasible under
the declared numerical checks, not a global optimality certificate. `operational_margin` optionally
requests a symmetric power inset during preparation, capped at half each
permitted range. Acceptance
always uses the original physical restrictions.
`initial` is an accepted schedule to retain. `warm` is an unchecked numerical
starting point after a forecast update; it is never returned without solving
and independent acceptance. Its arrays must use this case's object order/grid.
"""
function schedule_case(
    c::ScheduleCase;
    initial = nothing,
    warm = nothing,
    proposal_time_limit = 10.0,
    nlp_time_limit = 30.0,
    max_refinements = 2,
    operational_margin = 0.0,
    time_limit = Inf,
    proposal_optimizer = _PROPOSAL_OPTIMIZER,
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
        u=initial["u"]
        v=_verified_dispatch(
            c;
            u,
            initial,
            tunnel_opening=get(initial,"tunnel_opening",tunnel_opening_schedule(c)),
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
    seen=Set{Tuple{String,String,String}}()
    dispatch_options=(;started,deadline,nlp_time_limit,max_refinements,operational_margin,refine_grid)
    if best===nothing && warm!==nothing
        size(warm["u"])==size(fallback) && all(x->x==0||x==1,warm["u"]) ||
            throw(ArgumentError("warm commitment must be binary and unit by interval"))
        reference,start=_dispatch_start(c,fallback,warm)
        reference!==nothing || throw(ArgumentError("incompatible warm schedule: "*
            get(start,"ignored_reason","missing hydraulic arrays")))
        best,reference=_dispatch_candidate(c,reference.commitment,reference.tunnel_opening,
            fixed_river_on(c,nothing,reference),best,reference,attempts,seen;dispatch_options...)
    end
    radius=4
    for scale in (1.0,0.7,1.3)
        time()<deadline || break
        local_reference=best===nothing ? nothing : best.solution
        previous=best===nothing ? -Inf : best.solution.objective
        p=propose_commitment(
            c;
            water_scale = local_reference===nothing ? scale : 1.0,
            incumbent=local_reference,
            radius=local_reference===nothing ? nothing : radius,
            time_limit = proposal_time_limit,
            reference,
            deadline,
            optimizer=proposal_optimizer,
        )
        push!(proposals, p)
        # Dispatch before building the next surrogate so an accepted gain can
        # supply both its hydraulic reference and its nonlinear starting point.
        if haskey(p,"u")
            best,reference=_dispatch_candidate(c,p["u"],p["tunnel_opening"],p["river_on"],
                best,reference,attempts,seen;dispatch_options...)
        end
        # Expand after an unsuccessful neighborhood; restart locally after a
        # gain. The budget, rather than a claim of local completeness, stops us.
        radius=local_reference===nothing || (best!==nothing && best.solution.objective>previous) ? 4 : 2radius
    end
    if fallback_gates_valid && admissible(c,fallback)
        best,reference=_dispatch_candidate(c,fallback,fallback_opening,fixed_river_on(c,nothing,nothing),
            best,reference,attempts,seen;dispatch_options...)
    end
    Dict(
        "accepted"=>best!==nothing,
        "case"=>case_dict(best===nothing ? c : best.case),
        "solution"=>best===nothing ? nothing : best.solution,
        "audit"=>best===nothing ? nothing : _replay_report(best.audit),
        "proposals"=>proposals,
        "attempts"=>attempts,
        "seconds"=>time()-started,
        "global_optimality_proven"=>false,
        "requested_operational_margin_MW"=>operational_margin,
    )
end

# Fixed discrete controls define one nonlinear subproblem. The incumbent and
# hydraulic reference are explicit inputs and outputs; only the report and the
# set of previously tried controls are appended to in place.
function _dispatch_candidate(c,u,opening,river_on,best,reference,attempts,seen;
    started,deadline,nlp_time_limit,max_refinements,operational_margin,refine_grid)
    time()<deadline || return best,reference
    key=(join(vec(u)),join(vec(opening),','),join(vec(river_on)))
    key in seen && return best,reference
    push!(seen,key)
    admissible(c,u) || return best,reference
    # Find feasibility before optimizing revenue on a cold start. Retain an
    # accepted feasible trajectory if the subsequent optimization fails.
    for feasibility in (best===nothing ? (true,false) : (false,))
        time()<deadline || break
        try
            v=_verified_dispatch(c;u,tunnel_opening=opening,river_on,max_refinements,
                time_limit=nlp_time_limit,
                feasibility_only=feasibility,warm=reference,operational_margin,deadline,refine_grid)
            push!(attempts,Dict(
                "stage"=>feasibility ? "feasibility" : "dispatch",
                "accepted"=>v.accepted,"seconds"=>v.total_seconds,"steps"=>v.attempts,
                "u"=>copy(u),"tunnel_opening"=>copy(opening),"elapsed_seconds"=>time()-started,
                "warm_start_objective"=>reference===nothing ? nothing : reference["objective"],
            ))
            if v.accepted && (best===nothing || v.solution.objective>best.solution.objective)
                best=v
                reference=v.solution
            end
        catch error
            error isa Union{ArgumentError,DomainError} || rethrow()
            push!(attempts,Dict("accepted"=>false,"error"=>sprint(showerror,error)))
            break
        end
    end
    best,reference
end
