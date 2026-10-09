# Schedule generation

This study measures `schedule_case`: HiGHS proposes unit commitments, Ipopt
dispatches each unique proposal with the nonlinear equations, and independent
hydraulic replay accepts the resulting controls. Global SCIP search is outside
this experiment.

The selected change recovers usable local output and avoids a repeated original
equation audit. Existing nonlinear starting values and commitment proposals are
retained. Rebuilding hydraulic starts and refreshing the reference after every
improvement did not provide a reliable objective/time benefit.

## What changed

Local and global candidates now share the same guarded reconstruction. It permits
discharge roundoff no larger than 1e-6 m³/s and gate excursions no larger than
1e-8. Off-unit discharge becomes exactly zero. Larger control violations and
nonfinite controls remain failures. All storage, heads, tunnel flows, generation,
river arrivals and the objective are then recomputed from the corrected controls.
Raw solver diagnostics and correction magnitudes remain available.

Reconstruction is not acceptance. The reconstructed schedule must pass the
original equations and operating limits, finer chronological replay and the
independent transport check. An existing seasonal 24-hour probe is recovered
after correcting a gate excursion of 1.47e-14, with objective 453,194.141.
Two other retained probes become valid on their original grids but still fail
finer power/envelope checks and remain rejected.

`solve_verified` previously called `validate` immediately before `replay_audit`,
which called it again on the same unchanged values. Because validation includes
a complete hydraulic integration, this duplicated substantial work. The fresh
original audit inside `replay_audit` now supplies `model_valid` and `model_errors`.
Cached validation in a supplied schedule is never trusted. Finer replay and
transport checks are unchanged.

Warm values are used only when their arrays are finite and have compatible
dimensions. Explicitly mismatched grids are discarded for warm starts and
rejected for supplied initial candidates. Legacy schedules without grid metadata
remain supported. Compatible values retain the previous initialization behavior.

## Warm-start experiments

Four initial variants ran twice on each ordinary 2-, 6- and 24-hour Tokke–Vinje
case and the seasonal 24-hour case: release baseline, guarded reconstruction,
rebuilt/refreshed starts, and both changes. A second screen separated hydraulic
rebuilding from reference refresh on the 2- and 6-hour inputs.

The combined warm strategy reduced the 2-hour median from 4.792 to 1.837 seconds,
but reduced its accepted objective from 72,199.463 to 71,051.466, a 1.59% loss.
Its 6-hour result was both slower and slightly worse. Hydraulic rebuilding alone
retained objective 72,195.043 at 2 hours but was slower and worse at 6 hours.
Reference refresh alone also lost objective at 2 hours and increased 6-hour time.
These variants were removed from production.

Matching proposal and commitment hashes show that the initial comparison used
the same candidate commitments. Different starting values led Ipopt to different
local solutions. A physically consistent start is useful for auditing, but does
not guarantee a better local optimum in this nonconvex problem.

## Final operating comparison

The unchanged release and selected policy each ran twice under the same limits.
They produced the same minimum and median final objectives on every operating
case. Proposal hashes match, as do dispatched commitments; baseline dispatch
identities are inferred from proposal order, forced states and dwell checks.

| Case | Final objective, both medians | Median elapsed, release / selected | Rejected dispatches, release / selected |
|---|---:|---:|---:|
| Tokke–Vinje, 2 h | 72,199.463 | 4.795 / 4.827 s | 2 / 0 |
| Tokke–Vinje, 6 h | 217,703.587 | 14.278 / 14.110 s | 0 / 0 |
| Tokke–Vinje, 24 h | 743,663.504 | 48.874 / 49.330 s | 2 / 0 |
| Seasonal Tokke–Vinje, 24 h | 453,194.141 | 37.088 / 39.665 s | 4 / 2 |

Rejection counts sum both repetitions. Accepted output is retained even when
other attempts fail. Ordinary 24-hour objectives range from 743,645.751 to
743,681.256 in both implementations.

This is a reliability improvement with modest and mixed total-time effects.
The 6-hour total falls 1.2%; 2 hours is essentially unchanged. Ordinary and
seasonal 24-hour totals rise 0.9% and 6.9% because additional reconstructed
candidates undergo acceptance checks. The retained policy fixes false numerical
rejections, while the duplicate-audit removal limits their added cost.

At 2 hours, a recovered schedule with objective 71,025.043 is available at a
median 0.816 seconds. The release's first economic schedule is better,
72,199.463, at approximately 2.654 seconds. Attaining that same objective takes
2.730 seconds with the selected policy. Earlier usable output is therefore not
a speedup to equal quality. At 6 hours, the equal-objective time falls from
approximately 6.041 to 5.915 seconds; at ordinary 24 hours it falls from
approximately 15.286 to 14.885 seconds. Seasonal attainment rises from
approximately 25.309 to 28.011 seconds.
These events measure when the algorithm obtains an incumbent. `schedule_case`
returns after its candidate loop; this change adds no streaming result API.

Removing duplicate integration reduces 6-hour cumulative allocations from
3.986 to 3.601 GB. Recovery adds acceptance work elsewhere: 2-hour allocations
rise from 1.644 to 1.714 GB, ordinary 24 hours from 16.246 to 17.129 GB, and
seasonal 24 hours from 11.315 to 15.581 GB.

## Synthetic and supplied-schedule checks

Each additional case ran three times per implementation. The table uses the
last two repetitions: the first distributed-river call still includes about
5.4 seconds of compilation, and the first supplied-schedule call also compiles
additional code. Those first calls remain visible in the scalar evidence.

| Case | Final objective, both | Median elapsed, release / selected |
|---|---:|---:|
| Two-station PCHIP turbine tables, 8 h | 9,790.216 | 1.528 / 1.490 s |
| Six-unit distributed river network, 2 h | 14,860.129 | 0.059 / 0.057 s |
| Tokke–Vinje, 2 h, supplied schedule | 72,208.121 | 6.048 / 6.060 s |

All repetitions pass acceptance. The supplied schedule starts at objective
72,199.463; both implementations improve it. These are local preparation
objectives, not claims that the best known or global objective was reached.

The published evidence contains 74 measured scheduling calls: 32 in the initial
screen, eight in the separated warm-start screen, 16 in operating confirmation
and 18 in these additional checks. Warmup calls and the three retained raw
recovery probes are separate. The Julia suite passes 14,848 assertions,
including stale cached audits, malformed initial arrays, off-unit discharge,
bounded correction, incompatible grids and passive river laws.

## Remaining bottleneck

On ordinary 24-hour candidates, warmed construction takes roughly 0.03 seconds
and Ipopt typically consumes its 10-second allowance. At 2 and 6 hours, several
attempts reach the 1,000-iteration limit. Removing replay work cannot deliver a
large speedup while those iterations dominate.

A separate warmed 6-hour fixed-commitment profile reaches 1,000 iterations in
4.786 seconds of native Ipopt time. Its primal-dual linear-system timer accounts
for 4.103 seconds, about 86%. All function evaluations take 0.414 seconds, about
9%, and constraint Jacobian evaluation takes 0.033 seconds, less than 1%.
These timers are nested and must not be added together. No Hessian evaluations
occur: tabulated operators select limited-memory curvature. There are no
restoration iterations; mean line-search trials are 1.035, with a maximum of
seven. This profile points to linear-system work and convergence rather than
Julia derivative callback overhead as the immediate target.

The final scaled dual infeasibility is large, despite an almost-flat late
objective. This does not establish a cause, but warrants investigating
conditioning, curvature and scaling before optimizing callback allocations.
It is a solver diagnostic, not acceptance or local optimality evidence.

The next experiment should keep the same starts and equations and compare
sparse linear solvers, followed by curvature and scaling changes. The
existing tables have derivative discontinuities at knots; changing or smoothing
them would change the physical model. Analytic within-cell derivatives and
shared plant calculations remain possible experiments, with lower immediate
priority given this profile.
Ipopt.jl also provides a built-in SPRAL alternative to its default MUMPS linear
solver, with documented startup environment requirements. [Ipopt.jl linear
solvers](https://jump.dev/JuMP.jl/stable/packages/Ipopt/#Linear-Solvers).

## Measurement protocol

The operating comparison uses Julia 1.13.1 on Apple Silicon, the pinned project
manifest, one Julia thread, one BLAS thread and `OMP_NUM_THREADS=1`. Runs execute
serially. Two excluded complete 2-hour scheduling calls warm each process before
measurement. Each proposal has a 2-second allowance; each NLP has a 10-second
allowance. Refinement is disabled and the preparation margin is 0.1 MW.

These are per-solve allowances, not an end-to-end deadline. End-to-end times
include construction, optimization, reconstruction and acceptance; loading,
warmup, seed reconstruction and result serialization are excluded. Allocations
are cumulative temporary allocations, not peak resident memory. CPU/wall ratios
and compilation time are retained to expose contaminated measurements.

First economic improvement excludes the feasibility-only reference. Baseline
event times are bounded using ordered stage durations and unassigned overhead;
the selected implementation records them directly. Time to a common objective
uses the minimum final baseline objective across its repetitions. A variant
which never reaches it receives no attainment time.

Two repetitions support a descriptive comparison, not statistical significance.
Time-limited 24-hour objectives can vary slightly. These measurements establish
accepted schedules under unchanged checks; they do not establish global optima.

## Reproduction

```sh
./scripts/julia.sh benchmark/schedule.jl \
  --label selected --output results/scheduling \
  --case watercourse.json --proposal-seconds 2 --nlp-seconds 10 --repeats 2
```

Use `--warmup-case` to specify the warmup input, `--seed` immediately after its
case to benchmark an existing schedule, and `--refinements` for additional
refined-grid solves. `--help` describes the remaining measurement options.
The driver reports source, manifest, input, control and commitment hashes and
stores complete results locally. It has no production formulation selector
and never invokes global SCIP optimization.

Imported upstream inputs and schedules are not published. The scalar comparison
is retained in [scheduling.json](../benchmark/results/scheduling.json).
