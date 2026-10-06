# Selecting the production table graph

Freeze this protocol before confirmation results are collected. The comparison
is between exact representations of the same physical equations, not fitted
tables or different operating domains.

1. Screen original graphs, exact Cartesian bounds/cuts, fixed-state pruning and
   quadratic SOS2 PCHIP. Test an effective-flow lift and native SOS2 linear tables
   on the strongest candidates; both refinements are available to both families.
2. Freeze one candidate from each family. Confirm on normal Tokke–Vinje cases
   at 2, 6 and 24 hours, the seasonal 24-hour case, a synthetic PCHIP variant of
   the 6-hour watercourse, and both synthetic cases. The PCHIP watercourse changes
   only the declared turbine interpolation; it is not the original SHOP case.
3. Use one CPU thread, silent solver output, warmed construction/solve paths,
   identical frozen inputs and a common audited seed. Alternate solve order.
   Operating confirmation uses two repetitions, 300-second allowances and a
   relative gap target of 0.0001. Synthetic confirmation uses three repetitions,
   gap targets of 0.001 and 0.000001, with and without a supplied incumbent.
   The latter retains the model's default partial guesses. Preparation time is
   reported separately. Each measured total includes construction and audits.
4. Require original-equation validation, finer chronological replay, complete
   supplied-start audits, and bounds consistent with known feasible schedules.
   Publish scalar diagnostics only; imported upstream data remain external.
5. Select using free-commitment results. Fixed-commitment certificates are
   conditional diagnostics, not proof of scheduling optimality.

For each case, pool the largest independently audited discrete objective L.
Compare upper bounds using g = max(0,U-L)/max(1,abs(L)); investigate any upper
bound below L by more than 0.000001. A run earns certification-time credit only
if its own delivered schedule is accepted and certified. Its loss is then
wall_time/budget. Otherwise its loss is 1+log1p(g/target), plus any fraction of
the budget exceeded. Missing accepted schedules or usable bounds are failures.
Take median losses across repetitions and an equally weighted geometric mean
across cases; warm/cold and gap-target groups receive equal shares within a
synthetic case. Report delivered objectives and bounds alongside this score.

If the difference is less than 10% or ordinary repeat variation, prefer the
simpler graph. Stop after the effective-flow round and frozen confirmation;
keep one production implementation and remove public formulation switches.
This establishes the best measured method after diminishing returns, not a
proof that no future improvement is possible. Preserve this experimental
revision and measurements as research evidence.

## Frozen finalists

`cartesian_pruned` uses exact polynomial efficiency ranges and fixed-state
pruning, without the extra cell-weighted range cuts. Exact ranges give the best
free-commitment result in the first operating screen; adding cuts regresses it.
It retains the original linear-table cell graphs.

`tensor_refined_all` uses quadratic PCHIP products, fixed-state pruning and SOS2
for every linear table. The all-table extension has essentially unchanged small
case/free operating results, reduces model size, improves the conditional bound,
and lets production use a single table-coordinate primitive.

Neither finalist includes effective flow: its measured small-case certification
cost increases and its operating improvements do not justify another variable
and three constraints per unit-period. Confirmation compares these frozen
methods without further solver-parameter tuning.

Reference probes are independently reconstructed and fully lifted into each
graph before their objectives can invalidate a bound. Raw NLP objectives are
retained as diagnostics only. Sub-0.0001 synthetic objective differences and
pooled bound reversals at this scale remain unresolved numerical precision;
they are not evidence of better exact optima. Own certificate guards remain
unchanged, and a material fully audited counterexample rejects a solver bound.
