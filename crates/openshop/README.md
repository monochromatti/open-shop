# OpenSHOP

A Rust library for nonlinear hydropower scheduling on connected watercourses.
The numerical core describes water balances, hydraulic pressure laws, turbine
production and operating restrictions. Ipopt solves nonlinear dispatch; SCIP
proposes unit commitments. The library does not yet implement the complete
physical and operational model in the Julia reference.

The primary Rust interface is a prepared scheduling problem:

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

Inputs use metres, m³, m³/s, hours, MW and currency/MWh. Network incidence and
transport are prepared once; the problem borrows immutable inputs and can be
solved repeatedly. Physical validation checks equations and restrictions. It
does not establish global optimality or replace finer hydraulic replay.

From the repository's pinned Nix development shell:

```sh
nix develop .#rust
cargo run -p openshop --example reservoir
cargo run -p openshop --example schedule
cargo nextest run -p openshop
```

The first example constructs a case with native Rust structs. The second uses
the bundled synthetic watercourse, with storage and turbine tables, pressure loops,
delayed rivers, splits and merges, power ramps and a 1000 m³ forebay. It also accepts a JSON case path and a
time allowance in seconds. JSON handling belongs to the example application.

Building from source requires Rust 1.98 or newer, Ipopt, SCIP and libclang. The
Nix shell supplies their pinned versions and native header/linker paths.
The repository's `docs/rust.md` describes the supported equations and remaining
work. The package and synthetic examples are MIT licensed; solver dependencies
retain their own licenses.
