"""Schedule an input case using MILP proposals and audited nonlinear dispatch.

Limits apply per MILP/NLP solve; construction, replay and refinements are reported
in elapsed time. An accepted result is feasible under the declared numerical
checks, not a global optimality certificate. `operational_margin` optionally
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
)
    isfinite(operational_margin) && operational_margin >= 0 ||
        throw(ArgumentError("nonnegative finite operational_margin required"))
    started=time()
    validate_inputs(c)
    attempts=Any[]
    best=nothing
    if initial!==nothing
        u=initial["u"]
        v=solve_verified(
            c;
            u,
            initial,
            max_refinements = 0,
            time_limit = nlp_time_limit,
            operational_margin,
        )
        v["accepted"] ||
            throw(ArgumentError("supplied initial schedule fails independent acceptance"))
        best=v
    end
    fallback=ones(Int, length(c.system.generators), length(c.prices))
    for (j, g) in enumerate(c.system.generators), t in eachindex(c.prices)
        forced=opinterval(c, g.name, :forced_on, t, -1.0)
        forced>=0 && (fallback[j, t]=Int(forced))
    end
    if best===nothing && admissible(c, fallback)
        v=solve_verified(
            c;
            u = fallback,
            max_refinements,
            time_limit = nlp_time_limit,
            feasibility_only = true,
            operational_margin,
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
        p=propose_commitment(
            c;
            water_scale = scale,
            time_limit = proposal_time_limit,
            reference,
        )
        push!(proposals, p)
    end
    candidates=[p["u"] for p in proposals if haskey(p, "u")]
    admissible(c, fallback) && push!(candidates, fallback)
    seen=Set{String}()
    for u in candidates
        key=join(vec(u))
        key in seen && continue
        push!(seen, key)
        admissible(c, u) || continue
        try
            v=solve_verified(
                c;
                u,
                max_refinements,
                time_limit = nlp_time_limit,
                warm = reference,
                operational_margin,
            )
            push!(
                attempts,
                Dict(
                    "stage"=>"dispatch",
                    "accepted"=>v["accepted"],
                    "seconds"=>v["total_seconds"],
                    "steps"=>v["attempts"],
                    "u"=>copy(u),
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
