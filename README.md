# OpenSHOP

A Julia library for nonlinear hydropower scheduling on connected watercourses.
It models reservoirs, generating and pumping units, shared pressure networks, gates,
instantaneous channels and delayed rivers. SCIP searches the mixed-integer
nonlinear model; an independent chronological simulator checks the controls.

The [paper](paper/paper.pdf) states the complete physical and operational model,
its algebraic formulation and solution method. The editable source is
[paper/paper.typ](paper/paper.typ).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/monochromatti/open-shop")
```

The tested toolchain is Julia 1.13 on macOS and Linux. For development, the
repository pins its Julia packages and provides a Nix shell:

```sh
git clone https://github.com/monochromatti/open-shop.git
cd open-shop
./scripts/julia.sh -e 'using Pkg; Pkg.instantiate()'
./scripts/julia.sh examples/single_reservoir.jl
./scripts/julia.sh -e 'using Pkg; Pkg.test()'
```

Without Nix, use `julia --project=.` instead of `./scripts/julia.sh`.

## Use

```julia
using OpenSHOP

case = readcase("watercourse.json")
result = solve(case; time_limit=60.0, relative_gap=1e-3)
if result.accepted
    schedule = result.solution
    writejson("schedule.json", schedule)
end
```

`accepted` means the equations, operating restrictions and finer replay passed.
`global_certificate` means the schedule's validated interval-model objective and
SCIP's enclosing upper bound meet the requested gap. The two checks are separate:
an interval-model certificate can coexist with a failed finer replay.
A useful schedule may be returned before that gap closes.

Construct cases with native named Julia objects or JSON. The
[API guide](docs/api.md) explains solving, replay, results and restart;
[input objects and units](docs/input.md), [operating rules](docs/operating-rules.md)
and [river networks](docs/river-networks.md) describe case data.

Examples cover a [single reservoir](examples/single_reservoir.jl),
[energy-based terminal valuation](examples/terminal_value.jl),
[station operating rules](examples/operating_rules.jl),
[many-to-one river connections](examples/mixed_source_network.jl),
[distinct penstocks with a shared tailrace](examples/shared_tailrace.jl),
[prescribed river splits](examples/river_splits.jl),
[pumped-storage arbitrage](examples/pumped_storage.jl),
[signed channels](examples/level_channels.jl) and
[discrete tunnel gates](examples/tunnel_gates.jl).
[Spinning capacity](examples/spinning_reserves.jl) co-optimizes generation and
upward/downward reservations under the [declared reserve contract](docs/reserves.md).
[Directional handovers](examples/directional_interlocks.jl) keep a station off
for supplied waiting periods after intake closures.
[Explicit penalty limits](examples/soft_limits.jl) trade generation revenue
against a bounded terminal storage shortfall while preserving physical limits.
The [penalty settings API](docs/penalties.md) exposes effective defaults,
allowances, units and cost equations before optimization.
The [Tokke–Vinje importer](examples/tokke_vinje/README.md) reconstructs a documented
case from an external public dataset; upstream data are not bundled.

## Scope

Scheduling is deterministic. Head–storage relationships, signed quadratic
pressure losses, turbine and electrical efficiencies, finite gate choices,
commitment, schedules, directional ramps, environmental requirements and river
transit are explicit parts of the model. Constitutive functions must be supplied
and calibrated. Pumps use hydraulic efficiency data and constant motor efficiency;
reversible mode groups explicitly exclude simultaneous operation. See
[pumping](docs/pumping.md) for the model and supported scope.
A spinning-capacity contract adds unit/plant operating-point witnesses,
response limits, explicit shortfalls and capacity revenue. Other reserve
products, full activation trajectories, stochastic scenarios, river geometry
and transient pressure dynamics are outside the current model.

[Curve policies](docs/curves.md) explains single-head turbine data, explicit
extrapolation and zero-flow river-delay coverage, with hand-calculated examples.

Feasible-start preparation uses SCIP MILP proposals and Ipopt nonlinear dispatch.
The full SCIP search retains the original nonlinear equations. JuMP provides
the common algebraic interface; Julia also hosts physical reconstruction and
sparse river transport. See the paper for interpolation and certificate scope.

## Paper

The repository includes the typeset article, its Typst source and bibliography.
With Typst installed, run `make -C paper`. Alternatively, use the pinned Nix shell:

```sh
nix develop .#paper -c make -C paper
```

## License

The package and paper are MIT licensed. SCIP is Apache-2.0 licensed.
External datasets retain their own terms.
