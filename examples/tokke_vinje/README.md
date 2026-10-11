# Tokke–Vinje example

This example independently reconstructs a hydraulic and generation benchmark
from [SINTEF's public Tokke–Vinje repository](https://gitlab.sintef.no/energy/open-modelling-tools/open-datasets/tokke-vinje-watercourse/),
pinned to `ba2f2fc1f95d18978a04dd2c658aeef79126981b`. It reads public YAML, CSV
and literal notebook configuration; it does not execute upstream Python or
require a SHOP license. Upstream input files and generated cases are external
data, excluded from this package. No explicit redistribution license was found
in the reviewed upstream revision.

From the package root, run:

```sh
nix develop . -c python3 examples/tokke_vinje/fetch.py
nix develop . -c python3 examples/tokke_vinje/import.py --output examples/tokke_vinje/generated/operating
scripts/julia.sh benchmark/run.jl examples/tokke_vinje/generated/operating/tokke_vinje_24h.json results/tokke_vinje 300 2
```

The pinned Nix environment supplies Python 3 and PyYAML. Without Nix, the
importer needs Python 3 and either PyYAML or Ruby/Psych. Git is needed to
fetch and verify the pinned revision. Optional upstream Git LFS documents are
not needed. To use an existing source checkout, pass `--source /path/to/checkout`
to both scripts. The importer also accepts `--output`, `--start`, `--hours` and
`--profile operating|hydraulic`. The default operating profile imports current
hard environmental rules. The hydraulic profile reproduces the earlier
benchmark's scope without those rules.

The default uses dated historical prices, scaled historical inflows and
reconstructed initial reservoir levels. Shared intake and penstock losses are
represented explicitly. Turbine tables preserve reconstructed source knots,
using bilinear interpolation and explicit `head_extrapolation="linear"`.
Reference heads are not treated as legal operating limits. Startup and shutdown costs are
retained; missing initial unit state uses an explicit all-off policy.

The operating profile reads the current raw `qfomin`, `qmin` and `mamin` files;
the notebook's explicitly old schedules are excluded. Five river minimum
releases and three reservoir minimum-storage schedules repeat annually at
UTC midnight. Raw `mamin` thresholds are storage in Mm³, not elevations. They
combine with the existing physical minima. The source's December 31 baseline
is interpreted as January 1, and values hold until the next change.
These current source schedules are applied to historical benchmark dates;
historical rule revisions are not reconstructed.

The retained Vest reach combines Lio plant discharge, four incoming bypass/spill
arrivals and its dated natural inflow. Lio retains Bandak as its hydraulic head
reference while `discharge_river` sets its water destination. Natural inflow
enters Vest once; it is not also added to Bandak. A hard minimum on total Vest
release enforces the environmental requirement. Zero delay preserves the
earlier eliminated-reach physics. Minimum flows remain hard requirements. The
source's conditional extraordinary-inflow waivers are unsupported and never
activated implicitly. Initial historical levels are preserved; violations
cause an error rather than clipping. No current additional maximum-storage
schedule is supplied by these raw inputs.

Generate a case crossing September 30, when three release requirements and
the Ståvatn storage requirement change:

```sh
nix develop . -c python3 examples/tokke_vinje/import.py --start 2024-09-29T12:00:00Z --hours 24 --output examples/tokke_vinje/generated/seasonal24
```

Importer tests use synthetic fixtures and need no dataset. The optional
integration check uses an external pinned checkout and verifies historical
states and dated inflows for normal and seasonal cases:

```sh
nix develop . -c python3 examples/tokke_vinje/test_import.py
nix develop . -c python3 examples/tokke_vinje/test_import.py --source examples/tokke_vinje/data/tokke-vinje-watercourse
```

This is a restricted benchmark, **not a complete SHOP import or operational
schedule**. Reserve markets, river costs, conditional environmental waivers
and soft-bound penalty semantics are excluded. Reservoir and river
nominal bounds are hard restrictions. Tunnel mouths must stay submerged.
No calibrated river travel times are supplied. Terminal water values follow
an explicitly declared benchmark policy. The generated `metadata.json` records
every static source attribute, the complete source file inventory, historical
provenance, restrictions and omitted features. Unclassified static attributes,
missing history or invalid initial storage cause an error.

The benchmark requests a 0.1 MW margin during local seed preparation. This
addresses small replay overshoots at active power and turbine-envelope limits;
it does not alter the global SCIP model or the acceptance tolerances. Pass a
different `operational_margin` to `benchmark_cases` to change this preparation
margin; it is recorded in the benchmark results. The command-line runner's
third argument is the global solve allowance in seconds.

## River network validation

The 24-hour operating case was tested against main revision `61fef09` on
2026-10-10. Main's importer eliminated Vest; the current importer retains it.
The case includes 17 reservoirs, 14 generating units, 19 pressure-loss links
and 37 river reaches. Five of those links are documented physical tunnels;
the others represent shared intake and penstock losses.

Two unseeded scheduling runs used each configuration, after two excluded
two-hour warmups. Julia 1.13.1 and BLAS each used one thread, with five seconds
per MILP proposal, fifteen seconds per nonlinear solve, no additional grid
solves and a 0.1 MW preparation margin. The unchanged imported input had the
same SHA-256 in the first two configurations.

| Configuration | Accepted objective range | Median scheduling time |
| --- | ---: | ---: |
| Main, eliminated Vest | 742,568.627–742,630.850 | 73.425 s |
| New solver, same eliminated-Vest input | 742,592.638–742,605.257 | 74.415 s |
| New solver, retained Vest | 742,935.082–743,562.606 | 75.051 s |

All six schedules passed equation validation, conservation and finer replay.
The retained reach carries 8.00009–15.48258 m³/s in the first schedule,
including dated natural inflow of 2.86–3.77 m³/s. Lio operates for 22 hours,
with discharge between zero and 11.51235 m³/s. The scheduling pipeline improves
its first accepted feasibility schedule from −2,071,406.894 to 743,562.606.

The first retained-Vest schedule was reconstructed on both representations.
Objective, power and tunnel flow match exactly; maximum reservoir-storage
difference is 2.84×10⁻¹⁴ Mm³. The objective differences between unseeded runs
therefore reflect nonlinear search trajectories and time limits, rather than
a changed zero-delay feasible set.

Both models also ran native SCIP optimization with that common audited start,
free commitment and the same sixty-second allowance. Delivered objective was
743,562.606 in both; total times were 61.839 and 62.131 seconds. Neither improved
the supplied incumbent within that allowance. Their remaining optimality gaps
were 7.368% and 7.370%. The retained-reach model adds 72 variables and 96
constraints before presolving.

Both native results passed fourfold replay (96 intervals) and a restart after
twelve hours. Maximum storage discrepancy on replay was 0.000651 Mm³, maximum
water residual 8.41×10⁻¹³ Mm³ and restart storage discrepancy 4.83×10⁻¹³ Mm³;
terminal transit agreed exactly. These checks used the existing tolerances.

This validates plant outfalls, many-to-one merging, dated natural river inflow
and minimum reach release in optimization on the public reconstruction.
No river travel times or tunnel-to-river outfalls are supplied by the dataset;
those features retain their synthetic test coverage. The documented operating
exclusions still apply, and the time-limited SCIP runs do not establish global
optimality.
