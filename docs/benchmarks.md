# Solver selection

Native SCIP is the global solver used by OpenSHOP. A matched comparison with
Alpine.jl on the same bounded nonlinear equations favored SCIP for elapsed time,
bound quality and reliability under a time limit.

## Synthetic comparison

Measurements used Julia 1.13.1 on an Apple M1 Max with 32 GiB RAM. All engines
used one thread. Fresh models and solver instances received equal allowances;
compilation warmups were excluded. Construction, solution extraction and audits
are included in elapsed times. SCIP.jl 0.12.8 used SCIP 10.0; Alpine.jl 0.5.8 used
SCIP for relaxation bounds and Juniper/Ipopt for incumbents.

The first round requested a 0.1% gap with 30 seconds per run. The second round
requested 0.01% on the four smaller cases and 0.1% on the two larger cases, with
60 seconds per run. All unit commitment decisions remained free apart from a
declared unit outage.

| Case | Target | SCIP elapsed / gap | Alpine elapsed / gap |
|---|---:|---:|---:|
| Six units, two hours, distributed rivers | 0.1% | 0.255 s / 0.00947% | 22.385 s / 0.07185% |
| Six units, two hours, distributed rivers | 0.01% | 0.240 s / 0.00947% | 60.063 s / 0.07185% |
| Two stations, eight hours, PCHIP turbine tables | 0.1% | 0.726 s / 0.00701% | 5.245 s / 0.07923% |
| Two stations, eight hours, PCHIP turbine tables | 0.01% | 0.724 s / 0.00701% | 60.098 s / no returned bound |
| Six units, six hours, negative prices | 0.1% | 17.844 s / 0.08358% | 30.132 s / no validated incumbent |
| Six units, six hours, negative prices | 0.01% | 18.375 s / 0.00995% | 60.096 s / no returned bound |

Both methods' objectives agree closely when they return validated schedules:
14,860.129535 versus 14,860.129617 on the six-unit baseline, and 9,851.225564
versus 9,851.227255 on the tabulated case. The different objectives lie within
the certified tolerances. A tighter reported SCIP bound does not by itself
imply its incumbent is the best feasible schedule.

The six-hour SCIP schedules pass the discrete equations and gap checks but fail
finer replay. They are not physically accepted deliveries. The two-hour and
tabulated eight-hour cases pass both checks. The tabulated case includes four
starts and two shutdowns, so this comparison exercises actual commitment
changes rather than an all-on dispatch.

## Harder schedules

The remaining cases include a restricted 12-hour cascade, an eight-unit full-day
network with successive confluences, and a 12-hour distributed-delay network.
Neither method completed a global certificate on those cases in 30 or 60
seconds from its default start. SCIP consistently returned finite upper bounds;
several Alpine runs ended with an absent-primal result-extraction error. Such
failures and overruns are retained in the results.

A second experiment supplied the same independently validated feasible schedule
to both methods. Every algebraic start variable was populated and checked
against all constraints, bounds and binary domains before submission.

| Case | Shared seed objective | SCIP after 60 s | Alpine after 60 s |
|---|---:|---|---|
| Restricted 12-hour cascade | 9,715.775 | L=23,661.330, U=26,142.162, gap=10.48%; improved schedule fails finer replay | No returned bound; extraction failure |
| Eight-unit full-day network | 33,338.715 | L=33,338.715, U=87,718.787, gap=163.11%; seed passes replay | No returned bound; extraction failure |

Preparation cost 29.85 and 4.09 seconds respectively, charged equally to both
end-to-end workflows. The distributed 12-hour case produced no shared seed
passing finer replay, so neither seeded engine was scored on that case.

These results select SCIP as the simpler and stronger global-search foundation.
They do not establish fast full-day optimality or industrial SHOP parity. The
large remaining gaps are a concrete limit of the present formulation and search.

The complete 24 default-start and four shared-start records are in
[solver_comparison.json](../benchmark/results/solver_comparison.json). Numerical
certificates apply to the declared discrete hydraulic model. Distributed river
mixing and finer replay remain separate modeling checks.

## Reproducing the selected backend

The packaged implementation repeats the two replay-valid comparison cases
three times each. With the baseline default, median elapsed times are 0.240 seconds for the six-unit
two-hour case and 0.727 seconds for the tabulated eight-hour case. All six runs
pass the equation audit, finer replay and requested 0.1% gap. Their objectives
and bounds reproduce the original SCIP comparison. Records are in
[native_reproduction.json](../benchmark/results/native_reproduction.json).

```sh
./scripts/julia.sh benchmark/run.jl
```

The experimental tightened graph also passes all six synthetic runs. Its median
is 0.306 seconds for distributed rivers and 0.327 seconds for turbine tables;
the latter's gap is 0.000229%, versus the baseline's 0.00701%. The table case
improves, while distributed rivers become slower. All twelve scalar records are
retained in the same file, including a first-run elapsed-time outlier. The
operating Tokke–Vinje regression below prevents making this option the default.

Run the analytic example and tests with the pinned environment. The tests include
independently known off-state, on-state, interior-flow and negative-objective
optima; they check objective bounds and the reported gap.

```sh
./scripts/julia.sh examples/single_reservoir.jl
./scripts/julia.sh -e 'using Pkg; Pkg.test()'
```

The final external-data test uses the
[Tokke–Vinje reconstruction](../examples/tokke_vinje/README.md):

```sh
./scripts/julia.sh benchmark/tokke_vinje.jl
```

It freezes each operating case and one audited seed, then compares the baseline and tightened formulations with free and fixed commitment. The fourth argument sets the global allowance and the fifth sets repetitions. The records include preparation, construction, root/final bounds, search statistics, certificates and replay outcomes. To reproduce the earlier restricted input scope, import with `--profile hydraulic`.

## Operating Tokke–Vinje profile

The current importer adds nine source rules: five minimum river releases,
three seasonal minimum-storage schedules and the aggregate Vest flow minimum.
The normal cases begin at 2024-09-01 00:00 UTC. A separate 24-hour case begins
at 2024-09-29 12:00 UTC and crosses the September 30 rule changes. These inputs
have a different feasible set from the earlier hydraulic profile below.

The experiment freezes one case and one independently audited seed per horizon.
Baseline, network-domain and tightened-table formulations receive the same
controls and global allowance. The two-hour baseline and tightened pairs repeat
twice in alternating order. Other comparisons run once. Preparation, an
independent fixed-commitment dispatch probe and compilation warmups are excluded
from the scored times and recorded separately. The probe never replaces the
submitted seed; it checks that the solver bound encloses another known feasible
schedule. Fixed-commitment bounds apply only to the selected on/off decisions.
They still include tunnel-direction and table-cell choices where required; this
is a dispatch diagnostic, not a continuous convex subproblem.

The two-hour results did not justify promoting either experimental formulation.
With free commitment, the baseline gap was 52.13%, network domains 1,575.61%,
and tightened tables 2,466.94% after 60 seconds. All retained the same audited
objective, 72,191.373. Fixing commitment reduced the baseline gap to 13.95%;
tightened tables reached 23.82%. Both repetitions reproduced those bounds.

The smaller formulations performed more LP iterations but had difficulty
resolving the root relaxation. The first root LP bound was unavailable in both
experimental free-commitment runs. Variable count alone was a misleading
performance indicator: the baseline had 5,159 variables, network domains 4,363,
and tightened tables 4,229. This evidence supports retaining the baseline as
the default and keeping the two alternatives explicitly experimental.

At six hours, the baseline free/fixed gaps were 56.02% and 47.91% with a
120-second allowance. Network domains reached 55.99% with free commitment;
tightened-table gaps exceeded 1,200%. Fixing commitment therefore does not
resolve the longer case's dispatch-bound weakness.

The 24-hour baseline contains 336 unit on/off variables and 14,061 turbine-cell
selectors. Its exact piecewise data representation creates substantially more
search choices than physical unit commitment alone. A promising next experiment
must improve this representation and root-relaxation conditioning without
replacing the supplied physical curves or restricting the feasible set around
a local schedule.

The selected baseline delivered these operating schedules:

| Horizon | Preparation | SCIP + audit elapsed | Feasible objective L | Upper bound U | Free-commitment gap |
|---|---:|---:|---:|---:|---:|
| 2h | 33.71 s | 60.33 s median | 72,191.373 | 109,823.901 | 52.13% |
| 6h | 44.40 s | 121.06 s | 217,701.664 | 339,655.226 | 56.02% |
| 24h | 107.09 s | 123.12 s | 738,775.327 | 1,170,270.011 | 58.41% |

All pass the original equations and finer chronological replay. None closes the
requested 0.1% gap. The separate normal-case dispatch probe costs 18–20 seconds and is an
experimental diagnostic, not a required scheduling step. The two-hour seed is
reused from an earlier preparation; its original cost is shown above.

At 24 hours, network domains and tightened tables both reach an 866.53%
free-commitment gap. The baseline fixed-commitment run reaches 1,107.04%, versus
852.86% for tightened tables. Both lack a finite first root LP bound. These
unstable restricted runs cannot cleanly quantify the mathematical cost of
commitment; their poor bounds are retained as numerical performance evidence.

Construction costs approximately 0.2 seconds at two hours and 4 seconds at
24 hours; network-domain propagation itself takes only milliseconds. SCIP's
relaxation and search dominate elapsed time. Further Julia allocation tuning
would not address the main measured bottleneck.

The seasonal 24-hour baseline delivers an audited objective of 453,194.141,
with upper bound 743,920.105 and gap 64.15%. Seed preparation costs 92.56
seconds; the scored repeat takes 129.85 seconds with a 120-second allowance.
The dated case crosses reductions in three flow minima and the Ståvatn storage
minimum. An interrupted first pair is retained with `timing_usable=false`:
its driver and solver elapsed clocks disagree substantially. The scored repeat
inhibits idle sleep and reuses the same case and controls. No storage clipping,
conditional waiver or relaxed environmental rule is used.

Scalar records, case/control fingerprints, native root/final statistics and
start-audit summaries are in
[operating_formulations.json](../benchmark/results/operating_formulations.json).
Each comparison passed its formulation explicitly; the shipped default retains
the baseline. Upstream raw inputs and generated schedules remain external.
These measurements do not establish full SHOP compatibility or fast full-day
optimality.

An earlier numerical failure returned an apparent optimal bound below a
separately audited feasible fixed-commitment dispatch. The library now withholds
an optimality claim when LP iterations occurred without a finite first root LP
bound. The benchmark also rejects any bound below its independent probe.
These safeguards are conservative numerical checks, not rigorous arithmetic
certificates.

## Earlier hydraulic Tokke–Vinje profile

The original external-data test uses 17 reservoirs, 14 physical units, 19 tunnels
(including explicit intake/penstock loss branches), 36 retained rivers and
15 hydraulic junctions. The zero-delay mixing reach removed by the importer
is equivalent within the declared profile. All runs use hourly decision
intervals from 2024-09-01 00:00 UTC.

Local preparation requests a 0.1 MW operational margin. It addresses observed
power/envelope overshoots in finer replay, while the global SCIP model and
acceptance tolerances remain unchanged. Preparation includes a feasibility
fallback and three MILP proposals with audited Ipopt dispatch.

| Horizon | Variables | Preparation | SCIP + audit elapsed | Feasible objective L | Upper bound U | Gap |
|---|---:|---:|---:|---:|---:|---:|
| 2h | 5,159 | 39.09 s | 63.02 s | 73,512.047 | 107,110.796 | 45.71% |
| 6h | 16,043 | 19.67 s | 121.13 s | 224,519.398 | 341,577.817 | 52.14% |
| 24h | 66,251 | 71.37 s | 123.02 s | 766,840.028 | 1,188,178.992 | 54.94% |

All three delivered schedules pass the original equations, operational checks,
finer replay and conservation audit. **None closes the requested 0.1% gap.**
The global allowances are 60 seconds at two hours and 120 seconds at six and
24 hours. Construction takes 0.29, 0.65 and 2.57 seconds respectively; the
remaining elapsed time is mostly solver work. Small allowance overruns are
reported, rather than discarded.

Without the preparation margin, audited objectives were −28,915.595,
−2,289,527.146 and −2,086,600.389. Stronger local candidates were rejected for
finer-replay violations. The margin experiment improves the feasible seeds;
it does not demonstrate faster optimality proof or shrink SCIP's feasible set.

Both complete sets of scalar records, start-audit summaries and case hashes
are in [tokke_vinje.json](../benchmark/results/tokke_vinje.json). Upstream input
data and generated cases remain external. This is the restricted hydraulic
and generation reconstruction described in the example, not a full SHOP
import, a comparison with licensed SHOP, or an operational Tokke–Vinje plan.
The remaining global gaps and omitted operating rules prevent claiming a
competitive full SHOP replacement today.
