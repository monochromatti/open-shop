# Rust development

The Rust library and Python extension are being developed alongside the Julia
package. They use the same physical units except that Rust storage and river
inventory are in **m³**, rather than Julia's Mm³. The JSON schemas differ.
The Rust code does not yet implement the complete model in the paper.

## Build and run

```sh
nix develop .#rust
cargo nextest run -p openshop
cargo run -p openshop --example schedule -- examples/rust/cascade.json 10
```

The example has two headwater stations, a pressure loop, rivers that merge and
split, and a downstream station with a 1000 m³ forebay. The final argument is the
optimization time allowance in seconds. Output is JSON in physical units.

For Rust applications, deserialize a `ScheduleCase` or construct its public input
structs, then call `openshop::solve(&case, &options)`. Collection order determines
result axes. Names identify connections; input objects are never mutated.

The development shell pins Rust, rust-analyzer, Clippy, nextest, Ipopt and SCIP.
Cargo uses incremental compilation and modest optimization in development.
`nix build .#openshop` separately builds and tests the release library with crane,
which caches dependencies independently of the package source. Native solver
libraries come from Nix. No nightly Rust compiler is needed.

## Python

Build a wheel containing the extension and native solver libraries:

```sh
python3 -m venv .venv
source .venv/bin/activate
maturin build --profile dev --interpreter .venv/bin/python --auditwheel repair --out target/wheels
python -m pip install target/wheels/openshop-*.whl
python -m pip install pytest
python -m pytest python/tests
```

The wheel targets the Python version and operating system used to build it.
These are local development wheels, not published cross-platform releases.
Linux wheel repair also requires `patchelf`, included in the Nix shell.

```python
import json
from openshop import solve

with open("examples/rust/cascade.json") as source:
    case = json.load(source)
result = solve(case, time_limit=10.0)
if result["validation"]["feasible"]:
    print(result["objective"], result["unit_power"])
```

The extension copies Python inputs into owned Rust data and releases the Python
interpreter lock during optimization. Invalid input raises `ValueError`; model
construction and solver-interface failures raise `RuntimeError`. A solve that
returns an unfinished or infeasible point sets `validation.feasible` to false.

## Equations and current coverage

Reservoir storage obeys `V_next − V = 3600 Δt (inflow + incoming − outgoing)`.
The storage–elevation curve is quadratic. Its integral mean has a polynomial
expression that remains well defined when the two volumes are equal. Constant
discharge through a unit produces `P = 0.00981 Q H η(Q,H)` MW, with quadratic
efficiency and penstock losses. Unit commitment enforces either zero discharge
or the stated on-state discharge, head and power ranges. Startup costs and
terminal water values enter the revenue objective.

Tunnel flow is signed. Open tunnels satisfy
`(H_source − H_target) opening² = loss_factor Q |Q|`.
A closed tunnel has zero flow and imposes no pressure relation. Pressure
junctions conserve water without storage. Only junction heads are independent
head variables; reservoir heads and unit powers are derived expressions.

A fixed-delay river translates each release interval by its travel time.
At river junctions all arrival changes are retained, including historical
releases and plant outfalls. Outgoing reaches divide arrivals by prescribed
allocation fractions. Shared endpoint names define river junctions internally;
no additional public river-junction objects are needed. River cycles require an
intervening reservoir. River inventory includes water still in transit at the
end of the horizon, with its supplied terminal value.

Reservoir storage is piecewise linear between arrival changes. Interior states
are linear expressions of interval endpoints and incoming releases, so storage
bounds can be checked throughout an interval without adding storage variables.
Mean reservoir heads are integrated over those segments. This prevents an
arrival pulse from overflowing a small reservoir between control times.
Tunnel pressure equations still apply to interval mean heads; adaptive hydraulic
stepping and finer chronological replay remain to be ported.

Oximo builds the equations. The stable POUNCE expression evaluator supplies sparse
first and second derivatives directly in memory, including constraint curvature.
Native Ipopt solves continuous dispatch. SCIP proposes on/off schedules using
linear generation estimates and the same water balances and river transport.
Nonlinear dispatch evaluates those proposals, and search can change earlier
commitments. The proposal objective is **not a global bound**. POUNCE's optimizer
is exercised in comparison tests, rather than exposed as a user option.

`validation` evaluates water balance, pressure losses and operating limits in
physical units. It checks this interval model, not the finer replay or global
optimality provided by the complete Julia package. A feasible result alone says
nothing about optimality. `solver_status` is Ipopt's termination code;
`iterations` describes the final nonlinear solve, not the entire commitment search.

Still to be implemented are measured curve tables, flow-dependent distributed
river responses, adaptive hydraulics, pumps and reversible modes, optimized gate
positions, plant-level limits, dwell times, ramps, time-dependent operating
restrictions, restart, rolling horizons and global bounds. Full watercourses,
including Tokke–Vinje, must exercise those capabilities before any final
formulation or solver comparison. The small tests establish correctness of the
implemented equations; they cannot settle that comparison.
