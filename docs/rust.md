# Rust development

The Rust library is being developed alongside the Julia
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

Construct a `ScheduleCase` with native Rust structs and prepare a
`SchedulingProblem::new(&case)`. Preparation validates data and resolves network
incidence and transport once. Call `problem.solve(&options)` to optimize; the
prepared problem can be reused for another solve. It borrows immutable inputs,
so Rust prevents the inputs changing while its numerical structure is in use.
Collection order determines result axes. Names identify connections.

```rust,no_run
use openshop::{ScheduleCase, SchedulingProblem, SolveOptions};
fn optimize(case: &ScheduleCase) -> Result<(), openshop::Error> {
let problem = SchedulingProblem::new(case)?;
let schedule = problem.solve(&SolveOptions::default())?;
if schedule.validation.feasible {
    println!("Revenue: {}", schedule.objective);
}
    Ok(())
}
```

`crates/openshop/examples/reservoir.rs` constructs and solves a case entirely
in Rust. The `schedule` example reads JSON at the application boundary; JSON
parsing is not a library dependency. Serialization derives on data structs do
not perform I/O. The Rust API is a primary interface, not an implementation
detail behind a Python client. Python development is out of scope.

The development shell pins Rust, rust-analyzer, Clippy, nextest, Ipopt and SCIP.
Cargo uses incremental compilation and modest optimization in development.
`nix build .#openshop` separately builds and tests the release library with crane,
which caches dependencies independently of the package source. Native solver
libraries come from Nix. No nightly Rust compiler is needed.

## Unit operation

Each unit has a `UnitOperation` record. Its default has no minimum operating times
or discharge ramps, zero transition costs, and an initially stopped unit. Supply
`minimum_up_time`, `minimum_down_time` and `initial_duration` in hours. Minimum
times are measured from control edges in physical time, including an obligation
remaining from the initial state. The final state may carry an obligation beyond
the modeled horizon; this is not a requirement to finish every run before the end.

`flow_ramp.up` and `.down` bound changes in interval-average discharge per hour.
They use the distance between interval midpoints, including unequal intervals.
`initial_interval_hours` is the duration of the historical observation. An
initially stopped unit has zero initial discharge. An initially running unit can
supply `initial_flow`; without it, the first discharge ramp has no historical
reference and subsequent ramps still apply. Transition costs charge each start
and stop once. Prescribed commitments violating minimum times are rejected at
preparation. Power ramps and plant-wide operating restrictions remain to be ported.

## Numerical structure

`physics` contains scalar-generic constitutive equations used by both ordinary
evaluation and automatic derivative expressions. `model` defines physical inputs
and validates their domains. `river` prepares
conservative transport maps. `dispatch` assembles the nonlinear equations, while
`trajectory` evaluates physical states, revenue and residuals without invoking a
solver. `schedule` coordinates discrete proposals and nonlinear solves. Native
solver ownership, callbacks and mutable workspaces stay in their adapters.
File access and JSON conversion belong to the executable examples.

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
positions, plant-level limits, power ramps, time-dependent operating
restrictions, restart, rolling horizons and global bounds. Full watercourses,
including Tokke–Vinje, must exercise those capabilities before any final
formulation or solver comparison. The small tests establish correctness of the
implemented equations; they cannot settle that comparison.
