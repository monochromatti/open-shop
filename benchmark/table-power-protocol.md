# Table-aware power comparison

Baseline is OpenSHOP 0.4.2 at `77ea5ff`, including its supporting power planes and tight joint power envelope. All candidates retain the original equations, unit decisions, table interpolation, numerical certificate checks and physical replay.

The first screen compares baseline, head-coordinate supports, discharge-coordinate supports, both axes, and both axes with gated intercepts. Ungated supports add no variables; gated supports use at most one continuous intercept variable per row. None adds a binary or a new table coordinate. Polynomial Bernstein enclosures certify every support over the relevant domain; sampled feasibility checks validate the implementation but do not justify the inequalities.

Each of six hosted jobs uses one frozen case and one audited initial schedule, a 120-second allowance, one repetition and corrected root-event observations. Four jobs use the normal 2/6/24-hour and seasonal 24-hour Tokke–Vinje reconstruction. Two synthetic examples guard accepted numerical certification and model overhead. Seed preparation and warm-up times are reported separately. Construction, transformation, complete start lifting, native copy and optimization consume the global allowance; extraction, audits and statistics can overrun it. Scalar summaries retain failures and hashes. Raw input, seed and LP vectors are not uploaded.

Operating comparisons use final usable upper bounds against independently reconstructed objectives. Also compare each bound to the best accepted objective within its matched job. An apparently fast termination with no usable enclosing bound is a failure. Do not infer time to optimality from improvements at a fixed allowance.

The screen can motivate a compact refinement, but does not select production. Freeze finalists before repeated confirmation at a 300-second allowance, with profile order reversed on repetition two and root callbacks disabled. Evaluate all four operating cases and the small certification guards. Seek repeatable improvement on both 24-hour cases and the unweighted mean operating gap; report any regressions, numerical failures, objective changes, model growth and small-case cost. Prefer the smaller construction when proof progress is comparable.

Publish one production policy only if the evidence supports it. Otherwise retain the baseline. Archive all experimental machinery and remove private solve/start hooks from the library before completing this phase.

## Curvature refinement screen

The six-hour job in the first screen completed before the other operating jobs. At 120 seconds, head-coordinate rows lowered the accepted gap from 10.3800% to 8.4765%; discharge-only rows gave 10.4219%, both ungated axes 8.5650%, and both gated axes 7.7476%. All used the same accepted objective. This supports examining head-coordinate coupling and avoiding an assumption that more rows help.

The next frozen 120-second screen compares baseline, head_gate, axes_gate, head_curvature, axes_curvature, and axes_hybrid across the same six case definitions. head_gate tests whether the smaller gated construction suffices. Curvature rows use guarded derivative bounds over full original coordinate cells and fixed centers; no integer variable enters a quadratic. axes_hybrid retains the linear interval-inflated rows alongside convex curvature rows because neither dominates when SOS2 adjacency is relaxed. The 652 focused checks, small native solves and complete free/fixed Tokke–Vinje lifts passed before dispatch. Confirmation selection remains subject to all operating results, not the six-hour job alone.

## Frozen confirmation

Both completed screens accepted every measurement (30 linear,36 refinement). Both gated axes won all four operating comparisons in each screen. In the refinement the mean operating gap was7.1029%, against8.9957% for baseline and7.5435% for head-only gating; curvature and the hybrid were weaker. Both gated axes also outperform head-only gating on each 24-hour case by approximately0.44–0.46 percentage points. That consistent difference warrants the additional continuous auxiliaries. Synthetic guards retain certification; the PCHIP guard has a modest measured slowdown.

Freeze baseline and axes_gate for confirmation at300 seconds, two repetitions per case, reversed order on repetition two, no root-event callbacks. Run all six cases. Do not promote on screening data alone. The exact nonlinear equalities, tables, operational limits and initial-schedule/certificate audits remain unchanged.
