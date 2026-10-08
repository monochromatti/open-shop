# OpenSHOP

Hydropower scheduling in Julia, with nonlinear hydraulic equations and global
optimality bounds from SCIP.

OpenSHOP models reservoirs, generating units, shared hydraulic networks and
delayed rivers as named objects. It chooses unit commitment and dispatch in one
bounded mixed-integer nonlinear model. An independent chronological simulator
checks the returned controls, water balance and operating limits.

This is a generation scheduling library. It does not yet provide complete
compatibility with SINTEF SHOP. Pumping, reserve markets and stochastic scenarios
are outside the current model.

## Installation

Install from GitHub into a Julia environment:

```julia
using Pkg
Pkg.add(url="https://github.com/monochromatti/open-shop")
```

The tested toolchain is Julia 1.13 on macOS and Linux. The repository pins
its Julia packages and provides a Nix development shell for both platforms.

```sh
git clone https://github.com/monochromatti/open-shop.git
cd open-shop
./scripts/julia.sh -e 'using Pkg; Pkg.instantiate()'
./scripts/julia.sh examples/single_reservoir.jl
./scripts/julia.sh -e 'using Pkg; Pkg.test()'
```

Without Nix, use `julia --project=.` instead of `./scripts/julia.sh`.

The [single-reservoir example](examples/single_reservoir.jl) constructs a case
from Julia objects and solves it against a known optimum.

## Solving a schedule

```julia
using OpenSHOP

case = readcase("watercourse.json")
result = solve(case; time_limit=60.0, relative_gap=1e-3)

if result["accepted"]
    schedule = result["solution"]
    writejson("schedule.json", schedule)
end
```

Two checks have different meanings:

| Result | Meaning |
|---|---|
| `accepted` | The schedule passed the equations, operating constraints and finer chronological replay. |
| `global_certificate` | Its validated objective and SCIP's enclosing global upper bound meet the requested gap for the discrete optimization model. |

A time limit can leave a useful feasible schedule with an unfinished proof.
Inspect `feasible_lower_bound`, `global_bound`, `relative_gap`, `status` and
`total_seconds`. Solver termination alone is not an optimality certificate.
Construction and first-call compilation consume the time allowance; extraction
and replay can overrun it, recorded in `budget_overrun_seconds`. Passing
`fixed_u` restricts both the search and its certificate to that commitment.

The solver uses one exact table representation: SOS2 coordinates for linear
curves and shared quadratic products for discharge PCHIP interpolation. It
preserves the supplied table values and slopes. Units at the same plant share
hydraulic head and compatible table coordinates. Fixed unit states restrict
table domains to their physical branch. Shared on-state weights, a joint power
envelope and a bounded separator for certified power supports strengthen SCIP's
relaxation while retaining the exact equations.
See the [benchmark results](docs/benchmarks.md) and
[coupled relaxation experiments](docs/coupled-relaxations.md) for measured schedule,
bound and runtime tradeoffs.

For larger cases, a feasible initial schedule can help the global search:

```julia
prepared = schedule_case(case; proposal_time_limit=5.0,
    nlp_time_limit=15.0, max_refinements=0)
initial = prepared["accepted"] ? prepared["solution"] : nothing
result = solve(case; initial, time_limit=60.0)
```

Preparation uses HiGHS commitment proposals and Ipopt nonlinear dispatch. Its
objective is a feasible lower bound; SCIP supplies the global upper bound.
The proposal and NLP limits apply to each invocation: preparation tries three
MILP proposals and may dispatch several candidate commitments. Its total time
is reported separately from `solve`'s allowance. An initial schedule must
use the same grid and pass a complete algebraic feasibility check. If a better
discrete candidate fails replay, a replay-valid initial schedule is retained and
its gap is recomputed against the global bound.

Candidate preparation also accepts `operational_margin` in MW to keep
controls away from active power and turbine-envelope limits during replay.
This can make a narrow operating range infeasible; it is not a replay
guarantee. SCIP always uses the original case restrictions.

## Model and data

- Nonlinear head–storage curves, signed tunnel losses and flow conservation.
- Turbine tables, electrical efficiency curves and aggregate plant limits.
- Unit on/off states, minimum up/down times and transition costs.
- River confluences, release laws, finite travel times and environmental limits.
- Time-dependent operating restrictions, aggregate flow observations and restart state.

See [input objects and units](docs/input.md), [equations and certificate scope](docs/model.md)
and [benchmark results](docs/benchmarks.md).

The [Tokke–Vinje example](examples/tokke_vinje/README.md) fetches SINTEF's public
source dataset and reconstructs a documented generation benchmark. The default
profile includes supplied seasonal minimum-flow and reservoir storage rules. Its
mapping report identifies supported attributes, conditional waivers and exclusions.
Upstream data are not bundled with this repository.

## Why Julia and SCIP.jl?

JuMP provides the algebraic model, while SCIP.jl sends it to the native SCIP
solver. SCIP performs the nonlinear search; a bounded Julia separator supplies
certified power cuts.
Julia also hosts the typed network data, sparse transport compilation and
independent physical simulator. Keeping these in one language avoids a second
model implementation. [SCIP.jl](https://jump.dev/JuMP.jl/stable/packages/SCIP/)
exposes the solver's C API and MathOptInterface.

## License

OpenSHOP is MIT licensed. SCIP is Apache-2.0 licensed. External datasets retain
their own terms. OpenSHOP is independent of SINTEF SHOP and does not include its
source code or require its license.
