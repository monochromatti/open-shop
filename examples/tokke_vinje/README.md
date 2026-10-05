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
nix develop . -c python3 examples/tokke_vinje/import.py
scripts/julia.sh benchmark/tokke_vinje.jl examples/tokke_vinje/generated
```

The pinned Nix environment supplies Python 3 and PyYAML. Without Nix, the
importer needs Python 3 and either PyYAML or Ruby/Psych. Git is needed to
fetch and verify the pinned revision. Optional upstream Git LFS documents are
not needed. To use an existing source checkout, pass `--source /path/to/checkout`
to both scripts. The importer also accepts `--output`, `--start` and `--hours`.

The default uses dated historical prices, scaled historical inflows and
reconstructed initial reservoir levels. Shared intake and penstock losses are
represented explicitly. Turbine tables preserve reconstructed source knots,
using bilinear interpolation and explicit `head_extrapolation="linear"`.
Reference heads are not treated as legal operating limits. Startup and shutdown costs are
retained; missing initial unit state uses an explicit all-off policy.

This is a restricted benchmark, **not a complete SHOP import or operational
schedule**. Environmental minimum flows and storage rules, reserve markets,
river costs and soft-bound penalty semantics are excluded. Reservoir and river
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
