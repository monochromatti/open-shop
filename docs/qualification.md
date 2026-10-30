# Operational qualification

`scripts/qualify.py` measures the scheduling search separately from global bound
closure. It compares cold solves at several budgets against the best accepted
schedule found with more time. It can then execute a prefix, restart with the
actual storage, river cohorts and operating history, and compare cold and warm
replanning after changed inputs. The optimizer and physical model are unchanged.

Save a campaign plan outside the checkout, for example `../work/qualification.json`:

```json
[
  {"family": "cascade", "size": "large", "seed": "101",
   "restricted": true, "budgets": [30, 90], "reference_budget": 300},
  {"family": "pressure", "size": "large", "seed": "102"},
  {"family": "integrated", "seed": "103", "hours": 72, "branches": 8,
   "pond_volume": 0.0002, "recovery": ["forecast", "outage", "restriction"]},
  {"input": "tokke_vinje_72h.json", "recovery": ["forecast", "outage"]}
]
```

Run it with the package environment instantiated:

```sh
nix develop -c python3 scripts/qualify.py ../work/qualification.json \
  --output ../work/qualification-run
```

The six [fuzz families](fuzzing.md) and their sizes are available. `restricted`
adds the witness-compatible operating restrictions used by the fuzz campaign;
it defaults to false. The integrated family combines a pressure ring, selectable
gates, tabulated turbines, delayed river splits and merges, a small forebay and
downstream generation. Seeds vary its resistances, inflows, prices and water
values. Its defaults are 48 hourly intervals, four branches and a 2,000 m³
forebay (`pond_volume` is Mm³). An `input` can be a package case JSON or a saved
fuzz input containing reference controls; paths are relative to the plan.

Defaults are two repetitions at 30 and 90 seconds, followed by 300 seconds of
reference scheduling. Short runs cap each nonlinear solve at 45 seconds; the
reference allows an individual nonlinear solve to use its whole remaining
budget. Both use five-second commitment proposals and one hydraulic refinement.
Alternate repetitions reverse budget order. Each cold run starts independently,
without a physical reference or a previous solution. Reference scheduling may
start from the best qualified short solution or generated reference controls.
`global_reference: true` adds a separate SCIP solve with the reference budget
and the best qualified schedule as its initial solution. It is optional because
this campaign primarily investigates scheduling, rather than bound closure.
`public_budgets: [90]` additionally measures cold calls to the public `solve`
entry point, using its default time allocation. These rows are identified by
method so preparation timings are not mistaken for complete solve timings.
They use the same repetition count, with no supplied initial schedule. Their
search statistics record whether the returned schedule came from preparation
or global search, together with time spent in each phase.

Each case runs in its own process. A separate small-case warmup is recorded;
remaining compilation is included in measured calls and reported separately.
Reported solve time includes the solver's acceptance checks. Additional campaign
checks are outside that time. Phase-boundary budgets can overrun; actual elapsed
time is retained. `--timeout` sets a hard process limit, including warmup and all
repetitions and recovery events (default one hour). Run plans concurrently in
separate output directories on independent machines when comparing timings.

## Reading the results

Every returned accepted schedule is checked against the original equations,
objective accounting, an independent storage-plus-transit ledger, and four-times
finer replay. `accepted` records the package's acceptance decision; `qualified`
also requires the campaign's finer replay and water ledger. Only qualified
solutions contribute to economic comparisons. Generation failures, exceptions,
timeouts and replay sensitivity remain visible; inputs are never resampled.

`value_shortfall` is the best qualified objective found on **the same input**
minus the run's objective. `relative_shortfall` divides this by
`max(1, abs(reference_objective))`. It is an empirical comparison, not a bound
on loss relative to the unknown optimum. The absolute value matters when the
reference is near zero or negative. Missing schedules have no economic score;
they must be counted separately from successful runs. Neither an unsuccessful
short search nor an unsuccessful reference search proves infeasibility.

Recovery events are sequential, with a shrinking horizon. They do not yet test
extension of a fixed-length forecast window. Each executes `advance` hours
(default six) from the previous reference schedule, limited by the remaining
grid. First the unchanged continuation must reconstruct the same storage. Then:

- `forecast` reduces future reservoir inflows by 30% and reverses the remaining
  price sequence with an added 20 currency units/MWh.
- `outage` makes an eligible generating unit unavailable after two intervals,
  for at least four hours and its minimum down time.
- `restriction` reduces an eligible unit's discharge capacity to 60%, subject
  to its minimum flow, over the same interval.

The last two events exclude units with prescribed discharge, power or commitment
series to avoid inventing contradictory requests. Cold and warm recovery use
the same new case and the largest short budget, in alternating order. Warm
recovery receives only a shifted numerical guess; it cannot retain the obsolete
schedule as a feasible fallback. A generous reference run follows. The next
event executes that reference schedule, so `short_budget_recovered` must be
examined separately: a successful long rescue does not count as short recovery.
No changed-input event is assumed feasible merely because its parent was.

Outputs follow the fuzz campaign's process isolation and checkpoint convention:
`run.json`, exact inputs, full solver results, `partial.json`, per-case
`result.json`, and an incrementally updated `summary.json`. Use `--resume` with
the same source, plan and inputs to skip completed cases. Keep external inputs
and outputs under their source's redistribution terms. This slow campaign is
deliberately separate from per-PR CI; ordinary tests cover its accounting and
state-revision behavior without running long searches.
