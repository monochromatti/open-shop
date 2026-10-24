# Long-running optimization tests

The standalone campaign generates synthetic watercourses and operating inputs,
then tries fixed-commitment nonlinear dispatch, warm-started global optimization,
and a cold solve with ordinary feasible-start preparation. It is intended for
overnight runs, before releases, and after changes to the formulation or solvers.
It is not part of the per-PR optimization suite.

From the repository root:

```sh
./scripts/julia.sh -e 'using Pkg; Pkg.instantiate()'
nix develop -c python3 scripts/fuzz.py --cases 120 --seed 20261023 \
  --output ../work/fuzz-20261023
```

An existing Julia installation with the package environment instantiated also
works. Pass its executable with `--julia /absolute/path/to/julia`. Python needs
only its standard library. Each child uses one Julia thread and one BLAS thread.

Every case runs in a fresh Julia process. A crash or hung native solver cannot
take down the campaign. `--time-limit 60` gives each of its three optimization
calls a 60-second allowance; `--timeout 900` kills the entire case after 15
minutes, including compilation, input generation and physical checks. Solver
allowances do not interrupt compilation or replay. These process-level timings
are unsuitable for comparing warm solver performance.

## Inputs

The six families are visited in order, with an independent RNG stream for each
case. After a complete family cycle, the campaign advances to the next requested
size. Small cases have 2–3 principal reservoirs and 3–5 intervals; medium cases
have 4–6 and 6–10. `--sizes small,medium,large` adds 8–12 and 12–24. The balancing
pond and reversible-machine families retain two reservoirs while increasing the
horizon; pressure and river families add their receiver reservoir. Cascade
stations can have multiple units. These are ranges, not a model-size guarantee.

| Family | Physical and numerical variation |
| --- | --- |
| `cascade` | Several stations, shared station limits, polynomial or tabulated storage levels, analytic or measured-style turbine curves on unequal grids |
| `pressure` | Multiple supply reservoirs and stations around a growing ring of pressure junctions, tunnel loops, zero/partial/full gates, gate resistance tables, selectable or prescribed positions, gate-count interlocks |
| `rivers` | Several tributaries, plant outfalls, consecutive splits and downstream merges in larger cases, unequal travel times, natural inflow, nonzero history, fixed or flow-dependent distributed responses |
| `small_storage` | A 10–20,000 m³ balancing pond between two stations, with turnover much larger than its storage and a small independent inflow |
| `channels` | Signed level-difference channels between lakes; a dry, crest-adjacent or wet orifice, an uncontrolled weir, or a tabulated outlet |
| `pumping` | Reversible commitment exclusions, adjustable or fixed-speed pumping, conversion losses and varying lift |

Prices include negative values, zero and spikes. Inflows, initial storages,
storage-to-level curves, efficiencies and capacities vary. Unit histories and
minimum on/off times accompany on/off controls. Additional restrictions include
time-dependent reservoir bounds, unit flow and power limits, gross/startup head
limits at explicit reference endpoints, river switching and flow-change costs,
arrival observations and volume allowances, ramps, maintenance,
forced states, hard discharge schedules, tunnel capacities, late gate schedules,
outlet bounds, minimum release, a flow observation and a priced terminal target.
They are constructed to admit the reference schedule. The numerical values are
illustrative physical inputs, not measured equipment or watercourse calibration.

Generation does not retry until it finds an easy case. Invalid references and
hydraulic exceptions are saved, with their phase, for inspection. They do not
establish that the optimizers excluded a feasible solution. A validated reference
does establish feasibility of the original interval model; it need not be optimal
or pass finer replay.

Routed networks start from a steady prehistory: station outfalls, natural inflow
and both sides of each split obey conservation, with station flows consistent
with their initial on/off histories. Analytic outlets in the channel family are
instantaneous, so no arbitrary pre-horizon outlet inventory is invented.

## Checks and interpretation

Each input passes through the package's JSON admission rules. Each returned
accepted schedule is checked against the original physical and
operating equations. The campaign also recomputes objective components, closes a
separate whole-basin storage-plus-transit ledger, reconstructs the controls, and
runs four-times-finer replay. The reference schedule is continued from an
intermediate grid edge to check restart histories and river cohorts.

Failures include an accepted invalid schedule, exceptions or process crashes,
false infeasibility despite a valid reference, a global upper bound below a known
feasible objective, a rejected physical warm start, changed fixed controls, and
an unsupported global-certificate claim. Each new feasible candidate is compared
with bounds from earlier stages as well as its own bound. Physical audits remain numerical checks;
the campaign is not a proof of correctness or exhaustive domain coverage.

Time limits without a schedule are `search_inconclusive` warnings. Other local
search failures and a local objective below its feasible starting value are also
recorded for review. Finer-replay failures are `replay_sensitivity` warnings: they
can reflect coarse-step accuracy or a binding operating limit, even when the
declared interval equations were solved correctly. An unfinished global gap alone
is not a failure. Warm runs record native incumbent availability separately from
an independently retained reference schedule, so a fallback does not masquerade
as successful native search.

The final exit code is 1 for failures, 0 otherwise. `--strict` also fails on
warnings. An interrupted campaign exits 130. The per-case findings and aggregate
counts distinguish these outcomes.

## Artifacts, continuation and reproduction

All output stays outside the package checkout. `run.json` records the campaign,
source fingerprint, Git revision and host; each case has:

- `request.json`: seed, family, size, solver allowances and requested stages;
- `input.json`: the exact case and controls, saved before optimization;
- `reference.json`: the reconstructed reference when available;
- `local.json`, `warm.json`, `cold.json`: complete results for the requested calls;
- `progress.json`: the current phase, including when a process times out;
- `result.json`: checks, findings, versions and elapsed times;
- `worker.log` and global solver logs/statistics for diagnosis.

`partial.json` retains completed stages if a later stage crashes. The root
`summary.json` is updated after every case. Completed cases, including failures,
are skipped on resume; an interrupted, incomplete case is rerun:

```sh
nix develop -c python3 scripts/fuzz.py --cases 240 --seed 20261023 \
  --output ../work/fuzz-20261023 --resume
```

Resume requires the same inputs, solver settings and source fingerprint. For a
new revision, replay a saved physical input into a fresh output directory:

```sh
nix develop -c python3 scripts/fuzz.py \
  --replay ../work/fuzz-20261023/000017/input.json \
  --output ../work/fuzz-reproduction --time-limit 120
```

The saved input is the reliable reproduction across Julia/RNG versions; a seed
alone is not. Solver trajectories can still vary across versions or machines.
An input is ordinary package JSON nested under `case`, with its four control
matrices under `controls`, so it can also be reduced by hand into a regression
test. No automatic reduction or alternative optimization backend is introduced.

For narrower investigations, select `--families pressure,channels`, `--sizes
small`, or `--stages local,warm`. `--generate-only` runs generation and reference
checks without optimization. The normal test suite contains a small deterministic
smoke check of each generator; it does not run a random optimization campaign.
The process supervisor has separate fast tests:

```sh
python3 -m unittest discover -s test -p 'test_fuzz_runner.py'
```
