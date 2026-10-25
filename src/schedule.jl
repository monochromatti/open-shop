"""Schedule an input case using MILP proposals and audited nonlinear dispatch.

`time_limit` caps total preparation, including construction and refinements;
the MILP/NLP limits also cap individual solves. Final replay may overrun the
allowance and is included in elapsed time. An accepted result is feasible under
the declared numerical checks, not a global optimality certificate. `operational_margin` optionally
tightens candidate power/envelope constraints during preparation. Acceptance
always uses the original physical restrictions.
"""
function schedule_case(
    c::ScheduleCase;
    initial = nothing,
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
        v=solve_verified(
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
        v["accepted"] ||
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
    reference=best!==nothing && best["case"].grid==c.grid ? best["solution"] : initial
    proposals=Any[]
    seen=Set{Tuple{String,String,String}}()
    function dispatch_candidate!(u,opening,river_on)
        time()<deadline || return
        key=(join(vec(u)),join(vec(opening),','),join(vec(river_on)))
        key in seen && return
        push!(seen, key)
        admissible(c, u) || return
        # Establish a feasible trajectory before asking a cold NLP to trade
        # water against revenue. Keep it even if the subsequent polish fails.
        for feasibility in (best===nothing ? (true,false) : (false,))
            time()<deadline || break
            try
                v=solve_verified(
                    c;
                    u,
                    tunnel_opening=opening,
                    river_on,
                    max_refinements,
                    time_limit = feasibility ? nlp_time_limit/4 : nlp_time_limit,
                    feasibility_only=feasibility,
                    warm = reference,
                    operational_margin,
                    deadline,
                    refine_grid,
                )
                push!(attempts,Dict(
                    "stage"=>feasibility ? "feasibility" : "dispatch",
                    "accepted"=>v["accepted"],
                    "seconds"=>v["total_seconds"],
                    "steps"=>v["attempts"],
                    "u"=>copy(u),
                    "tunnel_opening"=>copy(opening),
                    "elapsed_seconds"=>time()-started,
                    "warm_start_objective"=>reference===nothing ? nothing : reference["objective"],
                ))
                if v["accepted"] &&
                   (best===nothing || v["solution"]["objective"]>best["solution"]["objective"])
                    best=v
                    # Refined arrays cannot anchor an original-grid proposal.
                    v["case"].grid==c.grid && (reference=v["solution"])
                end
                v["accepted"] && v["case"].grid!=c.grid && break
            catch e
                push!(attempts,Dict("accepted"=>false,"error"=>sprint(showerror,e)))
                break
            end
        end
    end
    radius=4
    for (iteration,scale) in enumerate((1.0,0.7,1.3))
        time()<deadline || break
        local_reference=iteration>1 ? reference : nothing
        previous=best===nothing ? -Inf : best["solution"]["objective"]
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
        haskey(p,"u") && dispatch_candidate!(p["u"],p["tunnel_opening"],p["river_on"])
        # Expand after an unsuccessful neighborhood; restart locally after a
        # gain. The budget, rather than a claim of local completeness, stops us.
        radius=best!==nothing && best["solution"]["objective"]>previous ? 4 : 2radius
    end
    fallback_gates_valid && admissible(c,fallback) &&
        dispatch_candidate!(fallback,fallback_opening,fixed_river_on(c,nothing,nothing))
    Dict(
        "accepted"=>best!==nothing,
        "case"=>best===nothing ? case_dict(c) : case_dict(best["case"]),
        "solution"=>best===nothing ? nothing : best["solution"],
        "audit"=>best===nothing ? nothing : best["audit"],
        "proposals"=>proposals,
        "attempts"=>attempts,
        "seconds"=>time()-started,
        "global_optimality_proven"=>false,
        "requested_operational_margin_MW"=>operational_margin,
    )
end
