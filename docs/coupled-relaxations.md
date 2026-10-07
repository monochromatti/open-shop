# Coupled relaxation experiments

OpenSHOP 0.4.2 adds a joint power envelope with tighter on-state bounds. Repeated five-minute comparisons reduce the unweighted mean operating gap by 13.5%, with unchanged delivered objectives. Small cases retain certification but take longer. Every experiment retains the exact nonlinear power, hydraulic, storage and river equations. The baseline is OpenSHOP 0.4.1, including its 12 head-aware power inequalities per tabulated unit and interval.

## Constructions

**Native RLT.** Reformulation-linearization cuts multiply valid linear rows by nonnegative bound factors and relax the resulting products. Four profiles allow 10 or 20 previously unrecognized products, expose non-original rows and variables, or combine the two. Root-only separation and the 100-variable limit remain unchanged. These are bounded tests of [SCIP's RLT separation](https://arxiv.org/abs/2211.13545), rather than exhaustive parameter tuning.

**Table consistency.** Nonnegative joint discharge/head weights have row and column sums equal to the existing SOS2 axis weights; efficiency equals their weighted table values. For PCHIP tables, additional moment constraints link the cubic correction products to the same head weights. Exact points admit the product assignment of the axis weights. All new rows are linear, with no new binary decisions. The motivation is [shared structure in piecewise multilinear relaxations](https://doi.org/10.1137/22M1507486).

**Hydraulic energy coupling.** For an open passive tunnel, write its exact law as d = r q |q|. The convex functions phi(q) = r |q|³ / 3 and phi*(d) = 2 |d|^(3/2) / (3 sqrt(r)) satisfy phi(q) + phi*(d) = q d. Summing hydraulic incidence yields an energy relation between tunnel dissipation, generation and external water injection. Linear supports bound the convex terms from below; bounded McCormick products relax head times injection. This is our adaptation of [duality-based strengthening for operational water networks](https://arxiv.org/abs/2208.03551).

The injection includes the declared midpoint storage change, inflow, actual delayed river arrivals minus releases, and fixed-head boundary exchange. Closed tunnels are excluded. Components with unsupported efficiency or outlet-loss signs are skipped. The discrete head times storage change is retained: substituting a continuous potential integral would be incorrect across storage-curve knots or for nonlinear level curves.

**Joint power envelope.** A multilinear product over a box admits a vertex convex-hull formulation. At most eight corner weights represent discharge, net head and turbine efficiency together. Their mass equals the on/off decision; two auxiliary products retain full off-state head and efficiency freedom. Power is bounded using safely computed electrical-efficiency extrema. The screen separates upper-only, upper-plus-lower, and both combined with hydraulic energy. The original nonlinear equality remains. This follows the [multilinear vertex formulation](https://arxiv.org/abs/2001.00514), adapted to our unit states.

## Validation and protocol

All 2,398 focused experiment checks pass: 452 for table consistency, 791 for hydraulic energy and 1,155 for the joint product envelope, including its conditional-box refinement. They cover physical auxiliary assignments, off states, operational limits, singleton axes, polynomial extensions, PCHIP corrections, signed and closed tunnels, boundaries, finite and distributed delays, and nonlinear storage. Small LP witnesses demonstrate that each construction excludes some points allowed by independent relaxations. Complete two-hour Tokke–Vinje starts also pass the existing algebraic and physical audit.

Each hosted job freezes one input and one independently reconstructed audited seed for all profiles. The allowance includes model construction, added constraints, start lifting, native copying and optimization. Seed preparation and measured warm-up times are separate. Extraction, replay and statistics contribute to elapsed time and can overrun the allowance. One Julia and BLAS thread are used. Results compare profiles within each job; different campaigns may prepare different schedules.

The first screen has seven profiles across six cases, once each at 120 seconds. The refinement compares four product-envelope profiles under the same protocol. Repeated confirmation uses a 300-second allowance with profile order reversed on the second repetition and diagnostic callbacks disabled. Unsolved cases are judged by final usable gaps; solved cases by time to accepted numerical certification. A bound improvement at the same budget is not evidence of faster time to optimality.

## First screen

All 42 runs deliver accepted schedules and valid lifted starts. Delivered objectives are identical within each operating job. No candidate convincingly improves both 24-hour cases.

| Profile | Normal 24 h gap | Seasonal 24 h gap |
|---|---:|---:|
| Baseline | 10.076% | 7.815% |
| RLT, 10 unknown products | 10.076% | 7.815% |
| RLT, 20 unknown products | 10.122% | 7.817% |
| RLT, non-original rows/variables | 10.076% | 7.815% |
| RLT, both changes | 10.124% | 7.818% |
| Table consistency | 10.020% | 7.874% |
| Hydraulic energy | 10.035% | 7.817% |

Allowing unknown products generates useful RLT cuts on smaller cases, but no applied cuts on the normal 24-hour case and only one or two on the seasonal case. Added table consistency nearly doubles native certification time on the small PCHIP example, from 0.544 to 1.035 seconds. Neither cut counts nor stronger small LP witnesses alone establish a production benefit.

## Product-envelope screen

All 24 refinement runs pass schedule and lifted-start audits. Delivered objectives again agree within each operating job. Joint power envelopes improve all four operating gaps:

| Profile | 2 h | 6 h | Normal 24 h | Seasonal 24 h |
|---|---:|---:|---:|---:|
| Baseline | 11.467% | 11.138% | 10.820% | 7.811% |
| Upper envelope | 10.519% | 10.115% | 10.128% | 7.593% |
| Upper and lower envelopes | 10.083% | 10.125% | 10.103% | 7.592% |
| Both plus hydraulic energy | 10.616% | 10.191% | 10.096% | 7.604% |

The product profiles add 3,360 continuous variables on the normal 24-hour case, a 15.7% increase, without adding binaries. Hydraulic energy adds another 1,248 variables and about 5,251 rows, with little improvement. Upper/lower power envelopes increase native certification time on the PCHIP guard from 0.554 to 0.871 seconds in this single screen; upper-only takes 0.667 seconds. The added proof strength therefore has a measurable small-case cost.

Both 24-hour cases still exhaust the allowance at one node. All captured affine-objective snapshots pass finiteness, primal reliability, relaxation status and objective consistency checks. The last observed normal 24-hour events show maximum original power discrepancies of 20.790 MW for baseline at 33.5 seconds and 18.516 MW for the upper/lower profile at 65.9 seconds. These are observations at different times and stages, not final dual-bound supports. At the first observed event, multiplication discrepancy falls sharply while the table contribution rises; fixing one relaxed relationship lets the remaining relaxation exploit another.

The confirmation also tests a compact refinement: compute conditional on-state efficiency bounds from the exact table or analytic curve, then derive a minimum head from the positive power minimum and safe flow/efficiency maxima. This changes only the envelope box and adds no variables or rows. Full off-state bounds remain intact.

## Diagnostic correction

An older driver read the current LP after optimization ended. That read can return stale or temporary values from bound tightening, so its pointwise residuals cannot be described as the final revenue root relaxation. The previously reported 21.38 MW, 7.43 m, 1.51 m and fixed-commitment 9.32 MW residuals are withdrawn. Native global bounds, gaps, accepted schedules and callback-free selection timings remain valid.

The corrected driver captures only optimal root LP solved events outside probing and diving. It records primal reliability, relaxation status, finite values and agreement between the native original LP objective and the captured affine revenue expression. First and last mean eligible observed events, possibly across restarts; they need not be the completed root or the LP supporting the final dual bound. Power residuals are decomposed into multiplication and table contributions at these infeasible points. They are not a decomposition of the overall optimality gap.

## Repeated confirmation and decision

All 48 confirmation runs pass the complete start and schedule audits. Profile
order is reversed on the second repetition, and root callbacks are disabled.
The selected tight upper/lower envelope gives the strongest bound on each
operating case in both repetitions:

| Case | Delivered objective, all profiles | Baseline gap | Upper only | Upper/lower | Tight upper/lower | Median elapsed: baseline / selected |
|---|---:|---:|---:|---:|---:|---:|
| Tokke–Vinje, 2 h | 72,197.273 | 11.193% | 10.409% | 9.759% | 9.419% | 300.30 / 300.29 s |
| Tokke–Vinje, 6 h | 219,896.833 | 11.124% | 10.010% | 9.821% | 9.175% | 300.86 / 300.86 s |
| Tokke–Vinje, 24 h | 743,600.989 | 10.135% | 9.501% | 9.397% | 8.961% | 302.08 / 302.19 s |
| Seasonal Tokke–Vinje, 24 h | 453,194.141 | 7.815% | 7.615% | 7.599% | 7.261% | 302.29 / 303.57 s |

The unweighted mean gap falls from 10.067% to 8.704%, a 13.5% relative reduction.
These are median gaps against each profile's independently validated objective;
using the best accepted objective pooled within each job gives the same operating
comparison. Preparation differs between campaigns, so their objective levels are
not comparisons of formulation quality.

Both 24-hour cases still process one node. At 6 hours the selected envelope
also remains at the root, while baseline processes about 9,400 nodes. A tighter
bound is more valuable here than a larger node count. No operating run reaches
the requested numerical optimality certificate, and this comparison does not
measure time to optimality.

The small PCHIP case certifies in a median 0.968 seconds with the selected
envelope versus 0.632 seconds for baseline; native optimization rises from
0.556 to 0.897 seconds. The analytic case takes 0.778 versus 0.718 seconds,
with native times 0.712 versus 0.645 seconds. All small confirmation runs
certify at the requested 1e-4 relative gap. Objective differences below 0.00005
are numerical precision. The larger operating proof gain is accepted despite
this small-case overhead; production has one fixed policy.

Normal 24-hour construction rises from 3.368 to 3.425 seconds. The envelope
adds 3,360 continuous variables and 4,704 linear rows: total counts rise from
21,404 to 24,764 variables and 26,767 to 31,471 rows. It adds no binary variables.
The original nonlinear equalities, SOS2 tables, supporting power planes,
operational restrictions and certificate checks are retained.

Production includes only the selected envelope. RLT parameters remain at
SCIP's defaults; table consistency, hydraulic energy and the other product
variants are archived. Temporary solve/setup hooks and generic experimental
start callbacks are absent from the library. Concrete Julia records and reused
table lookups avoid unnecessary metadata dispatch and allocations; no speedup
is attributed to that cleanup. The production checks pass 13,744 assertions
on Linux and macOS, plus the seven importer tests and the public example.
Public benchmark runs also verify free and fixed starts on both small cases
and a freshly imported two-hour Tokke–Vinje case.

## Evidence and remaining work

All 114 measured runs in this phase are accepted, with valid complete starts,
usable enclosing bounds and available native statistics. The
[scalar measurements](../benchmark/results/coupled-relaxations.json) retain
input/seed/source hashes, CPU metadata, timing, objectives, bounds, audits,
parameter changes, model growth and selected native counters. Full logs and
native statistics remain in the hosted artifacts. The
[first screen](https://github.com/monochromatti/open-shop/actions/runs/37614790616),
[product screen](https://github.com/monochromatti/open-shop/actions/runs/37618727359)
and [confirmation](https://github.com/monochromatti/open-shop/actions/runs/37622518720)
are preserved with the runnable
[experiment tag](https://github.com/monochromatti/open-shop/tree/coupled-relaxations-2026-10-07).
The first screen's postsolve last-LP observations retain an explicit caveat.

The remaining proof cost is still predominantly native relaxation work. In the
normal 24-hour confirmation, inclusive OBBT time falls from about 253 to
183 seconds; seasonal OBBT rises from about 227 to 248 seconds. Those plugin
timers overlap LP work and cannot be summed as exclusive costs. This is not a
universal reduction in tightening cost.

The next targeted work should couple the complete turbine-table/power relation
more closely, or refine selected head/flow regions where validated event
observations show large errors. Network and storage relaxations also remain
material. Global table-consistency growth and extra solver tuning did not earn
production complexity in this phase.

Numerical certification remains limited to the declared scheduling-grid equations and solver tolerances. Floating-point guards are not directed-rounding proofs. Tokke–Vinje remains a documented hydraulic/generation reconstruction with operating-rule exclusions, rather than a complete SHOP model.
