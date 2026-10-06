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

The eliminated Vest reach has a hard aggregate flow observation: Lio unit
discharge plus the four incoming bypass/spill releases plus its dated local
inflow. That local inflow is already delivered to Bandak; the observation adds
no water or hydraulic state. Minimum flows remain hard requirements. The
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
third benchmark argument of `0` to reproduce preparation without that margin.
Both settings are recorded in the benchmark results.
