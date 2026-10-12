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
        validate_tunnel_openings(c,fallback_opening)
    catch error
        error isa ArgumentError || rethrow()
        false
    end
    selectable=any(length(tunnel_allowed_positions(c,e,t))>1
        for e in c.system.tunnels,t in eachindex(c.prices))
    if best===nothing && !selectable && fallback_gates_valid && admissible(c, fallback)
        v=solve_verified(
            c;
            u = fallback,
            tunnel_opening=fallback_opening,
            max_refinements,
            time_limit = nlp_time_limit,
            feasibility_only = true,
            operational_margin,
            deadline,
            refine_grid,
        )
        push!(
            attempts,
            Dict(
                "stage"=>"feasibility",
                "accepted"=>v["accepted"],
                "seconds"=>v["total_seconds"],
                "steps"=>v["attempts"],
                "u"=>copy(fallback),
                "elapsed_seconds"=>time()-started,
            ),
        )
        v["accepted"] && (best=v)
    end
    reference=best!==nothing && best["case"].grid==c.grid ? best["solution"] : initial
    proposals=Any[]
    for scale in (1.0, 0.7, 1.3)
        time()<deadline || break
        p=propose_commitment(
            c;
            water_scale = scale,
            time_limit = proposal_time_limit,
            reference,
            deadline,
            optimizer=proposal_optimizer,
        )
        push!(proposals, p)
    end
    candidates=[(p["u"],p["tunnel_opening"]) for p in proposals if haskey(p, "u")]
    fallback_gates_valid && admissible(c,fallback) && push!(candidates,(fallback,fallback_opening))
    seen=Set{Tuple{String,String}}()
    for (u,opening) in candidates
        time()<deadline || break
        key=(join(vec(u)),join(vec(opening),','))
        key in seen && continue
        push!(seen, key)
        admissible(c, u) || continue
        try
            v=solve_verified(
                c;
                u,
                tunnel_opening=opening,
                max_refinements,
                time_limit = nlp_time_limit,
                warm = reference,
                operational_margin,
                deadline,
                refine_grid,
            )
            push!(
                attempts,
                Dict(
                    "stage"=>"dispatch",
                    "accepted"=>v["accepted"],
                    "seconds"=>v["total_seconds"],
                    "steps"=>v["attempts"],
                    "u"=>copy(u),
                    "tunnel_opening"=>copy(opening),
                    "elapsed_seconds"=>time()-started,
                    "warm_start_objective"=>reference===nothing ? nothing : reference["objective"],
                ),
            )
            if v["accepted"] &&
               (best===nothing || v["solution"]["objective"]>best["solution"]["objective"])
                best=v
            end
        catch e
            push!(attempts, Dict("accepted"=>false, "error"=>sprint(showerror, e)))
        end
    end
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
