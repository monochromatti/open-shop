# Solver workflow

OpenSHOP keeps HiGHS commitment proposals, Ipopt nonlinear dispatch and SCIP's
full nonlinear integer search. The public entry point is `solve(case)`, with one
shared time allowance. A prepared schedule supplies a complete feasible start;
its unit and gate choices remain free in SCIP.

SCIP alone is mathematically sufficient. The measured difference is whether it
finds a useful schedule within the allowance. In these experiments it found no
incumbent on either gate-network horizon or either Tokke–Vinje horizon within
60 seconds, in two repetitions each. Both prepared workflows delivered accepted
schedules. SCIP's installed build already calls Ipopt internally: the native
statistics record those solves. External dispatch adds a smaller physical model
with fixed discrete choices and direct turbine evaluation.

Replacing HiGHS with SCIP in the proposal MILP reduces some preparation times,
but does not dominate delivered objective quality. The conservative choice is
to retain HiGHS and simplify the interface. This is an empirical choice for the
current preparation algorithm, not a claim that three solver packages are
mathematically necessary.

## Matched comparison

Every arm starts from the same frozen input without a supplied schedule. The
60-second allowance includes preparation, construction, optimizer transfer and
SCIP search. Final extraction and replay can overrun it; elapsed time records
that overrun. Each case and method is warmed separately, with all warm-up
controls discarded. Measured methods run twice in rotating order, using one
Julia, BLAS and optimizer thread on the same Apple Silicon host.

The first campaign allocates at most 40% to preparation. The four synthetic
cases finish preparation within that allowance. Values are medians; higher
objective is better.

| Case | SCIP alone: objective / total seconds | HiGHS–Ipopt–SCIP: objective / total seconds | SCIP–Ipopt–SCIP: objective / total seconds |
|---|---:|---:|---:|
| Distributed rivers, 2 h | 14,860.130 / 0.443 | 14,860.130 / 0.401 | 14,860.130 / 0.413 |
| Turbine tables, 8 h | 9,851.227 / 1.613 | 9,851.174 / 1.075 | 9,851.174 / 1.066 |
| Twelve-reservoir gates, 6 h | No incumbent / 60.051 | 14,023.477 / 9.219 | 14,023.477 / 8.139 |
| Twelve-reservoir gates, 24 h | No incumbent / 60.052 | 84,634.302 / 60.169 | 84,634.302 / 60.142 |

The larger synthetic case has six units, ten tunnels, 23 river connections,
delayed confluences, selectable intake and bypass gates, interlocks, negative
prices and a signed connection between its two basins. Its SCIP proposal phase
takes 7.49 seconds versus HiGHS's 21.93 seconds. The first validated valuable
schedule appears after 3.47 versus 17.96 seconds. Both deliver the same objective
and final 1.285% gap; neither closes the requested 0.1% gap.

The 40% allowance is insufficient for a valuable 24-hour Tokke–Vinje schedule:
both methods retain a feasible fallback with objective −2,073,741.668. Giving
preparation 80% still fails when the replay correction restarts from the
optimized candidate. Retaining the original hydraulic start instead delivers
the following repeated results, with the same overall 60-second allowance:

| Tokke–Vinje | HiGHS proposals: delivered objective | SCIP proposals: delivered objective | Preparation seconds: HiGHS / SCIP |
|---|---:|---:|---:|
| 6 h | 219,711.487 | 217,625.414 | 44.12 / 34.06 |
| 24 h | 743,102.309 | 747,359.356 | 48.53 / 48.50 |

SCIP's delivered objective is 0.950% lower at six hours and 0.573% higher at
24 hours. At six hours its first MILP objective is better while its nonlinear
objective is worse: the surrogate ranks the commitments differently from the
physical problem. Solver substitution alone does not fix that approximation.

The six-hour gaps are 6.281% and 7.322%. The 24-hour gaps exceed 1,100% because
only about eleven seconds remain for global construction and search. These are
loose enclosing bounds, not estimates of actual schedule losses. This allocation
prioritizes an accepted schedule; a longer overall allowance gives SCIP more
time to improve its bound. Preparation is capped at 60 seconds even for longer
runs. The declared midpoint optimization model and its certificate scope remain
unchanged.

“First validated” can refer to a poor fallback. The retained records therefore
also report acceptance within 0.1% of the best delivered objective observed on
that input across the campaigns. That reference is empirical, not an optimum.
SCIP incumbents are audited after `optimize!` returns; native solution times do
not establish earlier physical acceptance.

## Integration

`solve(case; time_limit=60.0)` now prepares a start automatically. It reserves
at most `min(60.0, 0.8time_limit)` seconds for preparation, and gives all remaining
time to SCIP. Preparation tries the existing three water-value scales and fixes
their discrete choices during each Ipopt solve. A failed finer replay can cause
one retry with additional measured power or arrival margins. That retry retains
the decision grid and original hydraulic start. It does not alter the global
model's bounds or constraints.

The result reports `preparation_seconds`, `preparation_accepted` and
`prepared_objective`. A supplied `initial`, `fixed_u` or
`fixed_tunnel_opening` skips automatic preparation. Experts can disable the
accelerator with `preparation_time_limit=0`. The lower-level preparation and
dispatch functions remain available for explicit numerical experiments.

HiGHS's process-wide thread scheduler is left unchanged by the library. The
comparison runner configures its single-thread factories before any HiGHS
solve. Forcing a different thread count inside a library call conflicts with
previous HiGHS instances in the same process.
[HiGHS.jl thread documentation](https://github.com/jump-dev/HiGHS.jl#threads).

## Reproduction

The [runner](../benchmark/workflows.jl) and
[scalar records](../benchmark/results/solver-workflow.json) retain case hashes,
construction/search allowances, accepted objectives, enclosing bounds and
native work counters. The records contain 36 primary comparisons, four
allocation checks and eight stable-start confirmations. Raw Tokke–Vinje data
remain outside the repository. Its reconstruction uses the pinned public source
and the existing importer, including its documented operating-rule exclusions.

```sh
nix develop
python examples/tokke_vinje/fetch.py
python examples/tokke_vinje/import.py --hours 6 24
julia --pkgimages=no --project=. -e '
    include("benchmark/workflows.jl")
    benchmark_workflows("results/workflows",
        ["distributed", "tables", "gates6", "gates24", "tokke6", "tokke24"];
        budget=60.0, repeats=2, preparation_fraction=0.8)
'
```

The current source retains the original hydraulic start during replay retries.
The records distinguish the earlier retry and allocation experiments rather
than treating them as measurements of this selected policy.

The complete suite passes 19,515 assertions in 125 test sets. A fresh process
also checks the public `solve` entry point with default HiGHS settings: the
24-hour synthetic network and six-hour Tokke–Vinje use 60 seconds each, and
24-hour Tokke–Vinje uses 120 seconds. All three produce accepted schedules,
complete SCIP starts and passing fourfold replay. The latter delivers objective
742,898.077 with upper bound 797,351.914 and a 7.330% gap. Fourfold replay changes
its storage by at most 0.000650 Mm³. All 24 reconstructed delivered schedules
from the two prepared comparison methods also pass fourfold replay. These
additional checks are outside the measured comparison times.
