# Proposal solver comparison

OpenSHOP uses SCIP for commitment proposals and global nonlinear integer search,
and Ipopt for fixed-state nonlinear dispatch. HiGHS is confined to the separate
benchmark environment. The public `solve(case)` does not select among backends.
The proposal equations and preparation algorithm are unchanged.

The confirmation compares 96 preparation runs on twelve inputs, with three
repetitions per method and allowance. All 96 delivered schedules pass an
additional fourfold chronological replay. SCIP reaches a schedule within 1% of
the best observed objective in 48/48 runs; HiGHS does so in 47/48. Its one miss
is a 24-second Tokke–Vinje run that retains a feasible, poor fallback. With a
stricter 0.1% target the counts are 42/48 and 36/48 respectively. These targets
compare observed feasible schedules; they are not optimality gaps.

## Schedule quality and speed

The following are medians of three warmed repetitions with a 60-second total
preparation allowance. Positive objective changes favor SCIP. The quality time
is the first independent acceptance within 1% of the best fourfold-valid
objective observed on the same input across both methods and allowances.

| Input | SCIP objective change | Preparation seconds: HiGHS / SCIP | Time within 1%: HiGHS / SCIP |
|---|---:|---:|---:|
| distributed | +0.000% | 0.03 / 0.04 | 0.03 / 0.04 |
| gates24 | +0.000% | 20.67 / 6.78 | 17.07 / 3.17 |
| gates48 | +0.000% | 31.23 / 18.68 | 22.19 / 9.61 |
| gates6 | +0.000% | 2.38 / 1.42 | 1.86 / 0.89 |
| gates_restricted | +0.000% | 18.79 / 6.82 | 15.20 / 3.24 |
| gates_volatile | +0.018% | 16.43 / 3.64 | 15.87 / 3.09 |
| gates_wet | +0.000% | 22.00 / 7.81 | 19.07 / 4.89 |
| tables | +0.000% | 0.78 / 0.78 | 0.76 / 0.77 |
| tokke24 | +0.573% | 60.40 / 60.41 | 40.65 / 35.12 |
| tokke6 | -0.949% | 43.63 / 33.63 | 14.93 / 14.96 |
| tokke_prices | +0.128% | 60.43 / 60.39 | 38.24 / 37.93 |
| tokke_wet | +0.561% | 60.41 / 60.40 | 39.65 / 35.61 |

SCIP reaches comparable quality 5.4 times sooner on the daily gated network
and 2.3 times sooner on its 48-hour extension. Most of the difference is in the
proposal solves. On Tokke–Vinje those solves take roughly half a second at
24 hours; nonlinear dispatch dominates, so replacing the proposal engine has
much less effect on runtime. At six hours SCIP dispatches four unique
candidates instead of five, reducing total preparation by ten seconds, but
both methods find their useful schedule at about fifteen seconds.

SCIP's six-hour Tokke–Vinje objective is 2,086.073 currency units lower, or
0.949%. Its daily objective is 4,257.047 higher, or 0.573%. The physical
objective includes electricity revenue, transition/gate costs, release penalties
and changes in stored-water value. Initial stored-water value is already
subtracted; these percentages do not measure differences in generated energy.
The records retain the objective components.

Neither backend dominates every case. HiGHS saves about ten milliseconds on
the small distributed-river input, and SCIP loses some six-hour objective
quality. The substantial gated-network gains, comparable real-data behavior
and removal of a production dependency favor SCIP. A more accurate solve of
the linear surrogate need not produce a better nonlinear schedule: commitment
rankings depend on the frozen head and efficiency approximation. Both backends
use the same approximation, and every candidate faces the original physical
checks.

## Inputs and protocol

The inputs are the distributed-river and turbine-table fixtures, a connected
twelve-reservoir gated network at 6, 24 and 48 hours, three daily variations of
that network, and the Tokke–Vinje reconstruction at 6 and 24 hours with two daily
variations. These are four network families, including one reconstruction from
real data; twelve inputs do not mean twelve independent real watercourses.
The reconstruction retains its documented operating-rule exclusions and pinned
upstream revision `ba2f2fc1f95d18978a04dd2c658aeef79126981b`.

The synthetic variations alternate prices between −30 and 120, increase initial
storage by 10% of available room and inflow by 2 m³/s, or increase minimum
on/off durations to 4/3 hours and triple startup costs. The latter leaves the
best delivered commitment unchanged; it checks additional restrictions rather
than demonstrating more difficult unit switching. Tokke–Vinje variations
multiply alternate prices by 0.5/1.5, or increase initial storage by 2% of
available room and inflow by 20%. All other physical and operating data are
unchanged. The runner constructs these perturbations before measurements and
freezes each input with a content hash.

Each method uses the same three water-value scales, original hydraulic start,
fixed reference, per-solve limits, replay-driven margin correction and original
decision grid. Loading and per-network warmups are excluded; warmup schedules
are discarded. Method order rotates between repetitions. Julia, BLAS and both
proposal optimizers use one thread on the same Apple Silicon host. The installed
versions are Julia 1.13.1, JuMP 1.32.0, HiGHS.jl 1.26.0 / HiGHS 1.15.1,
SCIP.jl 0.12.8 / SCIP 10.0.3, and Ipopt.jl 1.16.0 / Ipopt 3.14.19.

Every input receives 60 seconds of preparation; the 6/24-hour real-data and
24/48-hour gated horizons also receive 24 seconds. Construction, all proposals,
dispatch and ordinary acceptance are included. Final acceptance can slightly
overrun the allowance. The additional fourfold replay is outside that timer.
The first acceptance timestamp uses the normal double-grid checks; the delivered
schedule is also checked at fourfold resolution.

The best observed objective is an empirical reference across both allowances.
This prevents a poor fallback from counting as useful simply because neither
short run improved it. The nonlinear solver can return an accepted iterate at
a time or iteration limit. Quality can therefore vary with the exact stopping
time and is not necessarily monotone with a larger allowance. Three repetitions
and synthetic perturbations support this implementation choice, not universal
solver superiority.

## Integration checks

Six additional whole-workflow checks use ordinary library settings: HiGHS keeps
its thread default, and SCIP uses one thread. The daily gated network and
six-hour Tokke–Vinje have 30-second shared allowances; daily Tokke–Vinje has
120 seconds. All six results have accepted schedules, complete valid SCIP
starts and passing fourfold replay. Their objective differences reproduce the
preparation comparison. The final daily SCIP check was repeated in isolation
following brief concurrent test startup; the discarded timing is not used.

Both methods consume their full allowance on these larger cases. Faster
preparation leaves more time for global search; it does not imply an equal
reduction in total elapsed time. The daily SCIP result has objective 747,359.356,
upper bound 798,233.448 and gap 6.807%. This does not establish improved global
proof speed: the equations, relaxation and proof strategy are unchanged.

The production environment passes 19,515 assertions in 125 test sets and all
three public examples. The constraint tests now use SCIP as well. Exhaustive
fixed-commitment checks start each independent feasibility problem with a fresh
optimizer, avoiding accumulation of MIP starts across unrelated assignments.
A fresh Julia process runs the single-reservoir example without loading HiGHS
and closes its requested gap. Including loading and example setup takes 60.61
seconds; the solve reports 43.47 seconds, including 34.15 seconds of preparation.
These first-call numbers include compilation and are separate from warmed
comparisons. The example permits 180 seconds so compilation does not exhaust
its demonstration allowance.

## Reproduction

```sh
nix develop
julia --pkgimages=no --project=benchmark -e 'using Pkg; Pkg.instantiate()'
python examples/tokke_vinje/fetch.py
python examples/tokke_vinje/import.py --hours 6 24
julia --pkgimages=no --project=benchmark benchmark/proposal-solvers.jl results/proposals
```

The [runner](../benchmark/proposal-solvers.jl) and
[scalar records](../benchmark/results/proposal-solvers.json) preserve objectives,
acceptance events, objective components, timing, input/source hashes and finer
replay results. Raw upstream inputs and control arrays remain outside the
repository. The [earlier workflow comparison](solver-workflow.md) establishes
why external Ipopt dispatch remains valuable alongside SCIP.
