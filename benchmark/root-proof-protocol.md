# Root proof experiments

The control is OpenSHOP 0.4.4: shared on-state table coordinates, the wider
certified power-cut allowance, and no extra OBBT bilinear projection LPs.
All profiles retain the exact original nonlinear equations and existing integer
and SOS2 decisions. This branch adds a private optimizer setup hook only for
matched experiments; it does not add a public solver policy selector.

## Profiles

- `current`: unchanged production configuration.
- `obbt_1`: ordinary OBBT iteration multiplier 1 instead of the native 10.
- `obbt_quarter`: multiplier 0.25 instead of 10.
- `obbt_cap30`: multiplier 1, then stop further OBBT calls at an ordinary root
  LP boundary after at least 30 cumulative native OBBT seconds.
- `joint`: current configuration plus bounded certified joint table-power cuts.
- `joint_obbt_1`: joint cuts plus multiplier 1; advanced after complementary
  early-cut and later-branching gains in the 6-hour screen.

The native 5,000-iteration minimum is retained. The multiplier determines a
**per-invocation** allowance, recomputed from root LP iterations. Zero would
mean unlimited and is never used. The cumulative stop cannot interrupt a call;
actual work and any overrun are recorded. It changes only OBBT frequency for
subsequent calls, retaining existing reductions and generalized variable bounds.

Joint cuts minimize an additive nodal envelope on the existing discharge and
head on-state weights. A small HiGHS LP imposes sufficient inequalities using
the exact power polynomial's Bernstein coefficients on every original cell and
linear extension. A separate production support-oracle pass checks the result
and adds a conservative correction. No sampled fit certifies a cut. Finite
coefficient bounds restrict which valid envelope is found; they do not relax
its validity. Construction, coefficient LPs and certification count toward the
cut work allowance and the solve's wall-clock time.

## Measurements

Each case freezes one imported input and one independently audited schedule.
Every profile and repetition receives the same schedule and starts the same
physical model. One Julia and BLAS thread, the locked Manifest, and Julia
1.13.1 are used. The second repetition reverses profile order. Warmup and common
seed preparation are excluded; model construction, copy, callback work and
native solving are included. Extraction and physical replay are included in
actual return time, which can exceed the allowance.

Only SCIP DUALBOUNDIMPROVED events and the final accepted solver bound supply
upper-bound trajectories and threshold times. They use original objective
units. A local root LP objective, converted from SCIP's transformed objective,
is recorded **only as a diagnostic**. LP snapshots require optimal primal/dual
reliable root LPs outside probing, diving and repropagation. SCIP.jl discards
the event pointer, so FIRSTLPSOLVED and LPSOLVED have separate handlers.

Root snapshots record cumulative OBBT time, completed calls, domain reductions,
LP counts/iterations, and table-power separator time/calls. OBBT's timer can
overlap probing/diving LP time; the ordinary dual-LP timer excludes probing
solves in the pinned solver. Reduction counts exclude probing changes and generalized
bounds, and cannot by themselves quantify useful tightening. Changes between
snapshots are associated with an interval of work, not a causal decomposition.

An unreached gap target is right censored, never assigned the time limit.
Callbacks, native statistics, complete-start audits and physical/replay audits
must pass. Every global upper-bound observation must enclose the best audited
schedule found anywhere in that matched job. Callback failures withhold bound
and timing results. Imported dataset files and controls remain ignored; only
scalar measurements, native statistics and diagnostic trajectories are shared.

## Stages

Screen five profiles on four operating cases and two small guards, at 120
seconds once per profile. Confirm any promising candidate against the control
at 300 seconds with two repetitions. Select using time to the same valid bound
and remaining gap at equal time, while disclosing horizon-specific tradeoffs.
A promising 120-second result alone does not change the production default.

The summary validates matched provenance and observations before computing
comparisons. It reports both final accepted bounds and the last global bound
observed within the nominal allowance. Common-bound times attained after that
allowance remain separate from times attained within it. Actual return time
also includes reconstruction and chronological replay.

```sh
python3 -m unittest discover -s benchmark -p root_proof_summary_tests.py
python3 benchmark/root_proof_summary.py results/root-proof-confirm \
  --output results/root-proof-confirm-comparison.json
```
