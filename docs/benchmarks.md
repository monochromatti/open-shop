# Benchmarks

OpenSHOP uses one production global model: native SCIP with exact quadratic
SOS2 turbine interpolation and SOS2 graphs for every linear curve. The choice
followed two optimization rounds and a frozen, repeated seven-case comparison.

The comparison preserved the physical equations, input files and audited
initial controls. It changed their algebraic representation. Neither method
obtained a materially better delivered objective on the operating watercourses;
the deciding differences were proof progress, robustness and model size.

The architecture comparison below describes version 0.4.0. The subsequent
[root relaxation study](performance.md) adds the production power bounds and
reports their matched performance separately. Version 0.4.2 then adds the
selected joint power envelope; its [coupled-relaxation comparison](coupled-relaxations.md)
reports the proof gain and small-case overhead. The subsequent
[table-coordinate comparison](table-power-bounds.md) selects gated table supports
for version 0.4.3. Its repeated comparisons reduce the mean operating gap by
21.8% against 0.4.2, with unchanged delivered objectives.

## Final matched comparison

These are free-commitment results with a 300-second allowance and a 0.01%
requested relative gap. Each pair ran twice on the same Ubuntu host, in
alternating order. Both methods delivered schedules that passed equation checks
and finer chronological replay in all 88 confirmation runs, including the
synthetic and conditional runs.

| Case | Delivered objective, both methods | Cartesian gap | Selected SOS2 gap | Median elapsed: Cartesian / SOS2 |
|---|---:|---:|---:|---:|
| Tokke–Vinje, 2 h | 72,196.856 | 10.668% | 11.452% | 300.26 / 300.25 s |
| Tokke–Vinje, 6 h | 217,669.103 | 14.507% | 12.874% | 300.58 / 300.47 s |
| Tokke–Vinje, 24 h | 743,693.298 | 1,112.640% | 10.640% | 312.76 / 302.72 s |
| Seasonal Tokke–Vinje, 24 h | 453,194.141 | 2,451.925% | 8.443% | 302.64 / 302.23 s |
| Synthetic PCHIP watercourse, 6 h | 217,757.050 | Unusable bound | 12.877% | 129.14 / 300.36 s |

Elapsed time includes construction, optimization, extraction and audits. Those
last steps can overrun the allowance. The reported gap uses each delivered
objective and its usable final upper bound; it is not SCIP's displayed gap.
The synthetic PCHIP watercourse changes only the declared turbine interpolation
in the 6-hour reconstruction. It is not the original SHOP calibration.

The Cartesian PCHIP run terminates as `OPTIMAL`, but has no finite first root LP
bound despite performing LP iterations. Its native bound also lies below a
better independently reconstructed and fully lifted nonlinear probe. Both
repetitions therefore receive no certificate and no usable upper bound. Its
shorter termination time is not a successful proof.

The selected model contains **21,404 variables** on the normal 24-hour case,
versus Cartesian's **66,104**. It is also smaller than the previous SOS2 model's
25,159 variables. This is a structural comparison; fewer variables alone do not
guarantee faster search.

The selected method wins on larger cases and supplies usable bounds on every
confirmation case. Cartesian is slightly stronger at 2 hours and faster on the
small PCHIP example. These advantages do not justify a second production path.
The frozen selection protocol scores free-commitment runs, treats missing
usable bounds as failures and gives equal weight to cases. Fixed-commitment
proofs do not decide the scheduling architecture.

## Synthetic certification and fixed commitment

Synthetic campaigns use three repetitions, gap targets of 0.1% and 0.0001%, and
both supplied and absent initial schedules. Julia and solver paths are warmed
before measurement. An absent initial schedule retains the model's default
partial guesses; it does not mean a fresh Julia process.

| Synthetic case | Initial schedule | Cartesian median | Selected SOS2 median |
|---|---|---:|---:|
| Two-station PCHIP, 8 h | Supplied | 0.17–0.18 s | 0.36 s |
| Two-station PCHIP, 8 h | Absent | 0.99–1.01 s | 2.03–2.04 s |
| Six-unit analytic network, 2 h | Supplied | 0.63–0.69 s | 0.63–0.68 s |
| Six-unit analytic network, 2 h | Absent | 0.80–0.82 s | 0.76–0.81 s |

Every synthetic confirmation run reaches an accepted numerical gap certificate.
The PCHIP objective is approximately 9,851.227 and the analytic objective
14,860.130. Differences smaller than 0.0001 objective units are numerical
precision, not evidence of a better exact optimum.

With commitment fixed to the common seed, Cartesian reports conditional
certificates at 2 and 6 hours in approximately 13 and 58 seconds. The selected
SOS2 method instead reaches conditional gaps of 2.682% and 3.033% at the
five-minute budget. At normal and seasonal 24 hours its conditional gaps are
4.394% and 0.128%. These restrict the feasible set to one unit schedule and do
not prove optimal commitment. Earlier screening produced material Cartesian
bound reversals against checked probes; those bounds remain rejected in the
published records.

## Optimization rounds and evidence

The first round compared original cell graphs and independent axes, exact
Cartesian efficiency ranges, additional range cuts, fixed-state domain pruning
and quadratic PCHIP products. The quadratic PCHIP graph improved the original
SOS2 synthetic certification cost by about 28% in its matched campaign. Exact
Cartesian ranges improved its small-case cost by about 11%; extra cuts regressed
the free operating bound.

The second round tested both families with an effective-flow auxiliary and
SOS2 for all linear curves. Effective flow increased synthetic certification
cost by approximately 43% for Cartesian and 23% for SOS2, without a compelling
operating gain. It was removed. The all-linear SOS2 extension had essentially
unchanged small-case cost, reduced the operating model and improved conditional
bounds. It provides one coordinate primitive across the production model.

The frozen finalists were Cartesian cells with exact efficiency ranges and
fixed-state pruning, and quadratic SOS2 turbines with fixed-state pruning and
SOS2 linear curves. Confirmation added the larger horizons, a seasonal-rule
transition and the derived PCHIP watercourse. No solver-parameter tuning occurred
between freezing the finalists and confirmation.

[Scalar measurements](../benchmark/results/selection.json) preserve screening,
refinement and confirmation records, source and control hashes, CPU information,
construction/solve times, native diagnostics and probe audits. Final confirmation
used Julia 1.13.1, JuMP 1.31.2, SCIP.jl 0.12.8 and SCIP 10.0.3 with one thread.
Different matrix jobs use different CPUs; comparisons are within each paired
job, not raw times across campaigns. Local laptop timings are excluded because
machine sleep contaminated earlier measurements.

The runnable comparison and predeclared selection protocol are archived at
[`table-contest-2026-10-06`](https://github.com/monochromatti/open-shop/tree/table-contest-2026-10-06).
The [confirmation workflow](https://github.com/monochromatti/open-shop/actions/runs/37421689056)
retains the hosted evidence. Production has no formulation selector or alternate
global table implementation.

## Reproducing the production benchmark

```sh
./scripts/julia.sh benchmark/run.jl
./scripts/julia.sh benchmark/run.jl case.json results/benchmark 300 2
```

The driver freezes an input and one independently reconstructed seed, warms the
solver paths, then repeats free and fixed commitment. Seed preparation and
nonlinear probes are measured separately from the global allowance. Probes
must pass the original equations and a complete model lift before they can
invalidate an upper bound. Rounded progress logs are diagnostics only.

For the external watercourse, follow the [Tokke–Vinje example](../examples/tokke_vinje/README.md).
Upstream cases and schedules are not bundled or uploaded. The hosted Performance
workflow publishes scalar summaries only.

The remaining free-commitment gaps are substantial. This work establishes one
measured production architecture, not full SHOP parity or fast global proof on
large watercourses. Bounds are numerical certificates for the declared discrete
model; finer replay checks controls but does not certify the continuous-time
optimum. The subsequent [performance study](performance.md) identifies continuous
relaxation errors as the immediate bottleneck and measures equation-preserving
power bounds without adding another table implementation.
