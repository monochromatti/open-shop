# Tokke–Vinje case

This example reconstructs a hydraulic and generation case from
[SINTEF's public dataset](https://gitlab.sintef.no/energy/open-modelling-tools/open-datasets/tokke-vinje-watercourse/),
pinned to revision `ba2f2fc1f95d18978a04dd2c658aeef79126981b`.
It reads YAML, CSV and literal notebook configuration without executing upstream
Python. Source data and generated cases remain external to the package.

From the package root:

```sh
nix develop . -c python3 examples/tokke_vinje/fetch.py
nix develop . -c python3 examples/tokke_vinje/import.py --hours 24 --output examples/tokke_vinje/generated
./scripts/julia.sh -e 'using OpenSHOP; c=readcase("examples/tokke_vinje/generated/tokke_vinje_24h.json"); r=solve(c; time_limit=300.0); writejson("schedule.json", r)'
```

The Nix shell supplies Python and PyYAML. Without Nix, use Python 3 with PyYAML
or Ruby/Psych; Git fetches and verifies the source revision. An existing checkout
can be supplied with `--source /path/to/checkout`. The importer accepts `--start`,
`--hours`, `--output` and `--profile operating|hydraulic`. The default operating
profile includes the current source's seasonal hard release and storage rules;
the hydraulic profile omits those rules. Input dates use UTC.

The reconstruction uses historical prices, scaled historical inflows and
initial reservoir levels. Shared intakes and penstocks are explicit pressure
links. Turbine knots are preserved with bilinear interpolation and declared
linear head extrapolation. Missing commitment history uses an all-off initial
policy. Terminal values use an explicit historical-price policy, recorded with
the case rather than inferred from the physical laws.

The Vest reach combines Lio discharge, incoming bypass/spill arrivals and natural
inflow. Bandak remains Lio's hydraulic reference; Vest receives the water.
Current annual minimum-release and storage tables are applied to the requested
historical dates; historical rule revisions are not reconstructed.

Conditional environmental waivers, reserve products, river costs and source
soft-bound penalties are outside this reconstruction. Nominal storage and river
capacities are hard. Tunnel mouths must remain submerged, and the dataset
supplies no calibrated river delays. Invalid initial storage fails instead of
being clipped. Generated `metadata.json` records provenance, mappings and
exclusions; consult it before interpreting a result as an operating schedule.

Synthetic importer tests need no external data:

```sh
python3 -m unittest discover -s examples/tokke_vinje -p 'test_*.py'
```

The optional integration check uses an external pinned checkout:

```sh
python3 examples/tokke_vinje/test_import.py --source examples/tokke_vinje/data/tokke-vinje-watercourse
```

No explicit redistribution license was found in the pinned source revision.
Fetch it separately and retain its original terms and provenance.
