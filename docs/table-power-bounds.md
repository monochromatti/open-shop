# Table-coordinate power bounds

OpenSHOP 0.4.3 strengthens SCIP's upper bound using the existing turbine-table coordinates. It retains the exact nonlinear power equality and adds no binary decisions. Repeated five-minute comparisons reduce the mean operating gap by 21.8% with unchanged delivered objectives. Quadratic curvature bounds and smaller variants remain in the experiment archive.

## Construction

A box envelope for discharge, head and efficiency permits combinations that the turbine curve cannot produce. The added rows bound the complete composed power function against the existing discharge and head weights, rather than treating efficiency as an independent factor. This follows the motivation of [composite relaxations in factorable programming](https://arxiv.org/abs/2310.07168); it is a compact construction for our existing table graph, not an implementation of that paper's full algorithm.

Let `F(q,h)=0.00981 e_max q h eta_turbine(q,h)`, with electrical-efficiency extrema covering the running power interval. For each of four fixed slopes, guarded Bernstein coefficients bound `F−a q` at each head node. The maximum residual against the node-to-node chord is then bounded over every physically usable portion of the original cell. Both endpoint coefficients are raised enough to cover that residual. An analogous construction bounds `F−b h` at discharge nodes. Original PCHIP slopes and secant extensions are preserved.

The resulting rows are `P≤a q+u Σc_j μ_j` and `P≤b h_on+u Σd_i λ_i`. Signed nodal sums use bounded continuous auxiliaries and four linear binary-product hull rows. On, the supports bound physical power. Off, discharge, power, on-head and gated intercepts vanish; the actual head and efficiency continuation remain free. Fixed commitments and constant nodal sums simplify directly. A shared support cache avoids recalculating identical polynomial bounds.

Eight power rows per applicable unit and interval supplement the existing table graph, supporting planes and joint trilinear envelope. This changes the continuous relaxation while preserving the original feasible integer schedules. Numerical guarding has the same scope as the existing power bounds; it is not an exact rounding-error proof.

## Screening and refinement

The first hosted screen compares five profiles across four Tokke–Vinje reconstructions and two synthetic certification guards. Each case uses one frozen input and one independently reconstructed audited seed, a 120-second global allowance and one repetition.

| Operating case | Baseline | Head only | Discharge only | Both axes | Both gated |
|---|---:|---:|---:|---:|---:|
| 2 h | 9.497% | 7.939% | 9.633% | 7.968% | 6.940% |
| 6 h | 10.380% | 8.477% | 10.422% | 8.565% | 7.748% |
| 24 h | 9.681% | 8.574% | 9.578% | 8.591% | 7.946% |
| Seasonal 24 h | 7.261% | 6.950% | 7.233% | 6.980% | 6.573% |
| Unweighted mean | 9.205% | 7.985% | 9.216% | 8.026% | 7.301% |

All delivered objectives are identical within each operating job. Head coupling helps consistently. Discharge rows alone do not, and adding both ungated axes slightly regresses the head-only result. Gating both axes provides the strongest bound under the common budget.

The refinement compares head-only gating, both-axis gating, head curvature, both-axis curvature, and both curvature plus the original linear rows. The convex variant bounds negative second derivatives over full original coordinate cells. Adding a fixed coordinate square makes the composed power function convex on each cell, so its endpoint chord bounds it above. Nonnegative raw intercepts and Jensen's inequality preserve off states. No binary enters the quadratic and no JuMP variable is added, although SCIP can create internal expression auxiliaries.

| Refinement profile | Mean operating gap at 120 s |
|---|---:|
| Baseline | 8.996% |
| Head gated | 7.544% |
| Both gated | 7.103% |
| Head curvature | 8.464% |
| Both curvature | 8.262% |
| Both curvature plus linear rows | 7.740% |

Both gated axes win every operating case in both screens. They outperform the smaller head-gated construction by about 0.44–0.46 percentage points on both 24-hour cases in the refinement. That repeatable difference warrants the extra continuous auxiliaries. Curvature guards are valid and can be faster on the small PCHIP case, but provide weaker operating bounds. They are not part of the production model.

## Repeated confirmation

Each profile ran twice at 300 seconds, with order reversed on the second repetition and diagnostic callbacks disabled. All reported gaps use independently accepted delivered objectives and usable native upper bounds.

| Case | Delivered objective, both profiles | Baseline gap | Selected gap | Median elapsed: baseline / selected |
|---|---:|---:|---:|---:|
| Tokke–Vinje, 2 h | 72,197.273 | 9.422% | 6.663% | 300.30 / 300.30 s |
| Tokke–Vinje, 6 h | 219,896.833 | 9.175% | 6.656% | 300.91 / 300.89 s |
| Tokke–Vinje, 24 h | 738,991.315 | 9.537% | 7.839% | 301.39 / 301.54 s |
| Seasonal 24 h | 453,194.141 | 7.261% | 6.533% | 302.44 / 302.41 s |

The unweighted mean operating gap falls from 8.849% to 6.923%, a 21.8% relative reduction. Both repetitions improve every operating case, including both 24-hour cases. Delivered objectives are identical within each paired operating job. The improvement is a tighter upper bound, not higher delivered revenue. All 90 hosted measurements across screening, refinement and confirmation are accepted, with valid starts and enclosing upper bounds. No operating case reaches the requested global gap certificate.

## Cost and limits

On the normal 24-hour case, the selected supplement adds 2,688 continuous variables and 13,440 constraints: the model grows from 24,764 variables and 31,471 rows to 27,452 and 44,911. No binary decisions or table coordinates are added. Four signed gated sums are created per axis and unit/interval, unless a fixed state or constant sum permits simplification.

On the normal 24-hour confirmation, median construction rises from 2.393 to 2.673 s while total elapsed remains around 301.5 s. Native root work dominates: the six-hour selected model spends 289.7 s in optimization-based bound tightening; the normal and seasonal 24-hour models remain at one node, spending 202.7 s and 174.7 s respectively in that propagator. Those timings can overlap LP and other native counters.

Both small guards retain numerical certification. PCHIP median elapsed rises from 0.641 to 0.797 s (+24.4%). The selected PCHIP repetitions take 0.882 s and 0.712 s; the first includes 0.170 s construction versus 0.024 s in the second. Both are retained. Native solve medians rise from 0.594 s to 0.669 s. The analytic guard adds no rows and stays essentially unchanged at 0.732 s versus 0.738 s. These costs are the tradeoff for stronger bounds on the operating cases.

A stronger upper bound at a fixed allowance is proof progress, not a demonstrated speedup to global optimality. The delivered operating schedules remain valuable lower bounds, and their gaps remain unresolved. The Tokke–Vinje reconstruction retains the importer’s documented hydraulic/generation subset and operating-rule exclusions; these experiments do not establish full SHOP parity. Certificates apply to the declared discrete model, not a continuous-time watercourse optimum.

## Validation and reproducibility

The experiment has 652 focused checks for guarded polynomial enclosures, unsafe raw nodal interpolation, full-cell curvature, clipped and extended domains, singleton coordinates, electrical efficiency maxima, off continuation and complete free/fixed starts. Small native solves and eight complete Tokke–Vinje lifts passed before refinement. An independent review found no material validity or accounting issue. The production port passes 2,076 row-equivalence checks against the frozen experimental construction, 13,961 library assertions, all seven importer checks against the pinned source, and the public example. A native Tokke–Vinje production solve is accepted; the small PCHIP case also certifies without a supplied schedule. Linux and macOS CI pass.

Each paired job freezes the same case and seed across profiles. Seed preparation and measured warmups are excluded from the global allowance. Construction, support computation, complete start lifting, native copy and optimization consume it; extraction, replay, audits and statistics can overrun it. Confirmation reverses profile order on repetition two and disables root-event callbacks. One Julia and BLAS thread is used. Cross-campaign objectives and raw timings are not paired comparisons because seeds and hosted CPUs can differ.

The [scalar measurements](../benchmark/results/table-power-bounds.json) retain all screening, refinement and confirmation records, including hashes, construction and warmup times, audit outcomes, model sizes, native bounds and plugin counters. All usable upper bounds are checked against the best independently accepted objective in their matched job. Native plugin timers can overlap and must not be summed as exclusive costs.

Root diagnostics come only from finite, reliable optimal root LP events outside probing or diving. Affine-objective events additionally check objective agreement. The last captured event is not a final dual-bound support; its discrepancies do not decompose the global gap. There is no postsolve LP reread.

The runnable alternatives and predeclared protocol are archived at [table-power-2026-10-07](https://github.com/monochromatti/open-shop/tree/table-power-2026-10-07). Hosted evidence is linked from the [linear screen](https://github.com/monochromatti/open-shop/actions/runs/37661173804), [refinement](https://github.com/monochromatti/open-shop/actions/runs/37663946086) and [confirmation](https://github.com/monochromatti/open-shop/actions/runs/37667324551). Production retains one construction and the existing [benchmark command](benchmarks.md#reproducing-the-production-benchmark).
