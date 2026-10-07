# Root relaxation performance

This page records the 0.4.1 power-plane selection. The subsequent
[coupled relaxation experiments](coupled-relaxations.md) assess additional bounds.

OpenSHOP strengthens its nonlinear power equations with head-aware supporting
inequalities. These give SCIP tighter upper bounds without new variables or a
second solver mode. The full nonlinear equations, turbine interpolation and
commitment constraints remain unchanged.

The full set is the production default: 12 supporting inequalities per tabulated
unit and interval. Full and reduced sets both lower the mean operating gap by
about 5.2%. The full set gives the tighter bound on both 24-hour cases and lower
PCHIP certification cost; the reduced set is slightly stronger at 2 and 6 hours.
No formulation selector is added.

The final matched comparison uses two repetitions per profile with a 300-second
allowance, one thread, identical inputs and identical initial controls:

| Case | Delivered objective, all profiles | Baseline gap | Full planes | Reduced planes | Median elapsed: baseline / full |
|---|---:|---:|---:|---:|---:|
| Tokke–Vinje, 2 h | 72,196.856 | 11.397% | 10.974% | 10.911% | 300.18 / 300.18 s |
| Tokke–Vinje, 6 h | 217,669.103 | 12.881% | 12.284% | 12.230% | 300.55 / 300.54 s |
| Tokke–Vinje, 24 h | 743,693.455 | 10.638% | 10.139% | 10.157% | 302.14 / 302.43 s |
| Seasonal Tokke–Vinje, 24 h | 453,194.141 | 8.453% | 7.811% | 7.867% | 302.04 / 302.81 s |

Delivered operating objectives are unchanged. Both 24-hour cases still spend
the full allowance at the root node. This is stronger proof progress at the
same compute budget, not a demonstrated reduction in time to optimality.

The small PCHIP example reaches accepted numerical certification in 0.399 s
with baseline, 0.490 s with full planes and 0.598 s with reduced planes
(medians). The full-set overhead is approximately 23%. Its tiny objective
difference, about 0.000021, is numerical precision. The analytic example has no
tabulated units and receives no added planes: native optimization times remain
approximately 0.61 s. Its total-time differences do not demonstrate a formulation
gain.

Normal 24-hour model construction increases from 1.42 to 3.77 seconds. It is
included in the allowance. The variable count remains 21,404; the constraint count
increases from 22,735 to 26,767 through 4,032 additional linear inequalities.

The normal and seasonal 24-hour Tokke–Vinje cases spend nearly the entire
allowance constructing and tightening the root relaxation. Baseline Julia model
construction takes roughly one second. The slow part is SCIP's proof work.

## Diagnosis

In the first paired 120-second campaign, optimization-based bound tightening
(OBBT) consumed 98.76 seconds on the normal 24-hour case and found three domain
reductions. It consumed 108.20 seconds on the seasonal case and found eight.
Both runs processed only the root node. These are inclusive native plugin
times; LP and plugin times overlap and must not be summed as separate costs.

Fixing all unit commitments still left a 4.89% conditional gap at 60 seconds
in the local normal 24-hour case. This remains evidence that unit commitment
alone does not explain the unfinished proof. Laptop timings do not select
production changes.

A later instrument review found that the diagnostic driver reread the current
LP after optimization. That read can reflect stale or temporary bound-tightening
values. The previously reported 21.38 MW, 7.43 m, 1.51 m and fixed-commitment
9.32 MW residual claims are withdrawn. Native bounds, accepted schedules and
callback-free selection timings remain valid. The
[coupled relaxation experiments](coupled-relaxations.md) use corrected event-only
observations outside probing and diving, with reliability and objective checks.

SCIP already detects unit/network symmetries in this watercourse. Its default
dynamic handling depends on branching decisions. Static symmetry restrictions
were tested to see whether they could help before branching. Unit histories,
efficiency curves, costs and limits remain part of the symmetry detector's
model; no manual ordering treats unequal units as interchangeable.

## Interventions

The screening varied one production formulation using private experiment hooks:

- Limit root separation to 20 rounds, reduce OBBT effort, or disable OBBT.
- Add four supporting power inequalities per tabulated unit and interval.
- Invert tunnel head-loss equations to bound signed flows and propagate
  zero-storage junction balances. Rivers are skipped in that propagation.
- Use static symmetry handling, alone and with hydraulic bounds.
- Add an aggregate plant energy upper bound at the shared net head.
- Add head-aware supporting planes to the power bounds.

All physical equations and binary commitments are retained. Supporting
constraints use polynomial Bernstein enclosures with floating-point guards;
no fitted replacement curve supplies an upper bound. The plant inequality is
P <= 0.00981 H sum(a_i Q_i), where a_i bounds the admissible turbine and
electrical efficiencies. An off unit has P=Q=0. Positive output and nonnegative
efficiencies force positive head for an on unit. Unsupported negative
electrical-efficiency domains are skipped.

The head-aware planes were the only intervention to improve all four operating
bounds in the 120-second screen: 11.890% to 11.458% at 2 hours, 11.767% to
11.137% at 6 hours, 10.789% to 10.235% at normal 24 hours, and 8.453% to
7.811% at seasonal 24 hours. Delivered objectives were unchanged within each
paired job. Those screening results motivated the repeated confirmation.

The other interventions did not establish a universal gain. Disabling OBBT advanced
more nodes on larger cases but worsened the 2-hour bound. Four supporting lines
improved the seasonal gap from 8.453% to 8.228%, with almost no change at normal
24 hours and a regression at 6 hours. Limiting root separation had no effect:
the relevant runs performed fewer than 20 rounds.

Static symmetry handling and extra hydraulic domain propagation supplied no
consistent improvement. SCIP's existing automatic symmetry handling remains
enabled. The aggregate energy inequality helped at 6 hours but regressed the
seasonal case. These experiments are archived rather than exposed as production
options.

A later campaign exposed a start-lifting corner case. Tiny negative forward
integration residuals at idle tunnels selected a reverse auxiliary direction,
conflicting with a proven nonnegative flow domain. Lifting now honors an
already fixed direction and retains the complete physical and algebraic audit.
A materially inconsistent flow remains rejected.

## Protocol and evidence

Within each matrix job, profiles use identical frozen input and independently
reconstructed initial controls. The allowance includes construction, native setup
and optimization. Extraction,
replay and statistics are included in elapsed time and can overrun the allowance.
Seed preparation is measured separately. Each profile receives an excluded
10-second warm-up allowance in its first repetition; this field records the
nominal allowance, not measured elapsed warm-up. A preceding common 15-second
warm-up is also excluded and has no elapsed field. Comparisons use one thread
and alternate order when repeated.
Different jobs have different CPUs and sometimes different prepared schedules;
objectives and times are compared within paired jobs, not across campaigns.

Root snapshots and displayed progress are diagnostics. Final bounds and gap
certificates use the existing numerical certificate checks for the declared
scheduling-grid equations. Initial controls and raw LP vectors remain local;
only scalar measurements are published.

The final confirmation uses baseline, full planes (12 inequalities per tabulated
unit and interval), and a reduced set (8), twice each on six cases. Profile order
is reversed for the second repetition. Root callbacks are disabled in this
confirmation. The two small examples measure time to accepted numerical
certification; the operating cases measure the final usable gap at the same
300-second allowance.

The screenings recorded 98 runs, of which 97 were accepted. The physics campaign
is incomplete: one invalid lifted start stopped its normal 24-hour job after
two profiles. That failed record is retained. The direction-lifting fix was
then verified in the subsequent screening and final confirmation.

All 36 final confirmation runs passed schedule and lifted-start audits. Across
all four campaigns, 133 of 134 recorded measurements were accepted. The
production tests passed 13,054 assertions on both Linux and macOS, plus the
seven dataset-importer tests and the public example. Regression checks cover
polynomial extensions, electrical efficiency, tightened operating limits,
negative off-state heads and fixed tunnel directions.

[Scalar measurements](../benchmark/results/root-performance.json) include source,
case and seed hashes, CPU metadata, timings, audit results, usable bounds and
selected native statistics. Full native JSON and displayed trajectories remain
in the hosted artifacts;
the published progress rows are sampled at declared time milestones.
The runnable experimental driver is preserved at
[`root-contest-2026-10-06`](https://github.com/monochromatti/open-shop/tree/root-contest-2026-10-06).
The [confirmation workflow](https://github.com/monochromatti/open-shop/actions/runs/37524736859)
is the selection evidence. Production keeps one formulation and the original
SCIP parameter defaults.

The power inequalities are conservative polynomial enclosures in exact
arithmetic. Their implementation uses floating-point guards, consistent with the
project's numerical certificate scope; it is not a directed-rounding proof.
The substantial remaining operating gaps mean this phase establishes better
bounds, not optimal schedules or full SHOP parity. The subsequent [coupled relaxation experiments](coupled-relaxations.md) test
links between power products, table graphs and hydraulic conservation. Fixing
commitment alone does not remove the continuous-relaxation gap.
