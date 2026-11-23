# Rust development

The Rust library is being developed alongside the Julia
package. They use the same physical units except that Rust storage and river
inventory are in **m³**, rather than Julia's Mm³. The JSON schemas differ.
The Rust code does not yet implement the complete model in the paper.

## Build and run

```sh
nix develop .#rust
cargo nextest run -p openshop
cargo run -p openshop --example schedule -- crates/openshop/examples/data/cascade.json 10
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
parsing stays in the example; `serde_json` is only a direct development dependency
of OpenSHOP. Serialization derives on data structs do
not perform I/O. The Rust API is a primary interface, not an implementation
detail behind a Python client. Python development is out of scope.

The development shell pins Rust, rust-analyzer, Clippy, nextest, Ipopt and SCIP.
Cargo uses incremental compilation and modest optimization in development.
`nix build .#openshop` separately builds and tests the release library with crane,
which caches dependencies independently of the package source. Native solver
libraries come from Nix. No nightly Rust compiler is needed.

## Storage curves

Use `StorageCurve::Polynomial { reference_volume, elevation, slope, curvature }`
for an analytic relationship, or `StorageCurve::Table(TableCurve::new(x, y)?)`
for measured storage/elevation pairs. Storage coordinates are m³ and elevations
are metres. Tables retain the supplied knots and linear segments: there is no
fitting or smoothing. Input data must be finite, ordered, nondecreasing in level,
and cover the complete allowed reservoir volume range.

`curve.level(volume)?`, `curve.mean_level(before, after)?` and
`curve.potential(volume)?` evaluate the physical relationship. Table lookups reject
coordinates outside the measured domain. The integral mean is computed cell by
cell using weighted slopes; it does not divide by a zero volume change or subtract
nearby cubic potentials. The potential primitive is in m⁴; multiplying its change
by water density and gravity gives the potential-energy change in joules.

The nonlinear model uses one mean-head equation per interval for each nonconstant
polynomial or measured curve. The equation integrates every segment of the reservoir's actual
piecewise-linear storage path, including changes caused by river arrivals.
Its first and second derivatives are assembled locally through the affine storage
paths. Sparse structure and chain-rule coefficients are prepared once, and numeric
callbacks reuse their buffers.
A table constant over its operating range needs no head unknown.

Trial points use end-secant continuation to let a solver recover from violated
constraints; storage bounds and the measured domain are still enforced and
validated. At knots, derivatives use the right segment except at the final knot.
These are piecewise derivatives: the table is not smoothed, and classical
second derivatives need not exist at a knot. Selected chain-rule identities are
checked in Lean; finite-difference tests also check the complete sparse Hessian.
Neither establishes a global optimality proof or verifies floating-point arithmetic.

`crates/openshop/examples/data/tabulated_cascade.json` uses synthetic storage and turbine tables with
pressure loops, delayed river splits/merges, commitment, power ramps and a 1000 m³
forebay. These numbers are synthetic, not field measurements:

```sh
cargo run -p openshop --example schedule -- crates/openshop/examples/data/tabulated_cascade.json 10
```

## Turbine curves

Use `Efficiency::constant(value)`, `Efficiency::Quadratic { ... }`, or
`Efficiency::Table(TurbineTable::new(heads, curves, flow, options)?)`.
`heads` are ordered net hydraulic heads in metres; each has its own original
`TableCurve` of discharge (m³/s) against hydraulic efficiency, and a `Bounds`
record giving its on-state discharge envelope. The curves may use different
discharge grids. The envelope is independent of the efficiency knots and adds
to the unit's static discharge limits.

```rust,no_run
use openshop::{Bounds, DischargeInterpolation, TableCurve, TurbineOptions, TurbineTable};
let table = TurbineTable::new(
    vec![80.0, 120.0],
    vec![
        TableCurve::new(vec![5.0, 12.0, 20.0], vec![0.82, 0.91, 0.87])?,
        TableCurve::new(vec![5.0, 9.0, 15.0, 20.0], vec![0.83, 0.90, 0.93, 0.89])?,
    ],
    vec![Bounds::new(5.0, 15.0), Bounds::new(5.0, 20.0)],
    TurbineOptions { interpolation: DischargeInterpolation::Pchip, ..Default::default() },
)?;
# Ok::<(), openshop::Error>(())
```

Interpolation in discharge is linear or shape-preserving cubic Hermite (PCHIP),
with weighted harmonic interior slopes following Fritsch and Butland (1984,
[DOI: 10.1137/0905021](https://doi.org/10.1137/0905021)). Each original curve is interpolated independently; the
result and the discharge envelope are linear in head. No new knots are fitted.
`table.efficiency(q, h)?` and `table.flow_limits(h)?` evaluate the declared model.
Power always uses the actual net hydraulic head, including penstock loss.

By default the data must cover the complete unit operating range. Explicit
`Extrapolation::Linear` permits end-secant continuation in head or discharge.
Discharge continuation uses the end secants even for PCHIP. A single reference
head requires `head_extrapolation: Extrapolation::Constant`, declaring efficiency
and the envelope independent of head. It does not make power independent of head.
Preparation checks efficiency throughout the allowed rectangle, including cubic
stationary points, and checks that discharge envelopes overlap at every head.
Values are neither clipped nor silently repaired.

Table units use net head and power unknowns with the physical equalities
`H = H_source − H_target − loss_factor Q²` and `P = 0.00981 Q H η(Q,H)`.
Their power limits are variable bounds and their power ramps are linear rows.
The local machine law supplies exact first and second derivatives using fixed-size
three-variable `num-dual` arithmetic; sparse chain-rule maps are prepared once. The head equality
retains the curvature of the upstream hydraulic equations. Off units have zero
flow and power and do not impose a turbine operating-head range. Trial evaluation
uses secant continuation to let a solver recover from violated constraints;
returned schedules must satisfy the declared physical domains.

## Unit operation

Each unit has a `UnitOperation` record. Its default has no minimum operating times
or ramps, zero transition costs, and an initially stopped unit. Supply
`minimum_up_time`, `minimum_down_time` and `initial_duration` in hours. Minimum
times are measured from control edges in physical time, including an obligation
remaining from the initial state. The final state may carry an obligation beyond
the modeled horizon; this is not a requirement to finish every run before the end.

`flow_ramp.up` and `.down` bound changes in interval-average discharge per hour.
They use the distance between interval midpoints, including unequal intervals.
`initial_interval_hours` is the duration of the historical observation. Supply `initial_flow` for a historical interval average; an observation for a
stopped unit must be zero. Without an observation, the first discharge ramp has
no historical reference and subsequent ramps still apply. Transition costs charge each start
and stop once. Prescribed commitments violating minimum times are rejected at
preparation.

`power_ramp.up` and `.down` similarly bound changes in generation, in MW/hour,
using the calculated nonlinear power. Starting permits an additional increase
of `power.min`; stopping permits the corresponding decrease. In symbols,
`P_t − P_prev ≤ up Δτ + P_min start` and
`P_prev − P_t ≤ down Δτ + P_min stop`, where `Δτ` is midpoint spacing.
Discharge ramps have no additional transition allowance. Supply `initial_power` for a historical average; omitting it leaves the first
power ramp without a historical reference. A supplied observation for a stopped
unit must be zero.
Historical flow and power need not lie in the new horizon's operating range;
the ramp equations determine whether the transition is possible. Plant-wide restrictions use the same operating equations on aggregate flow and
power.

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
The storage–elevation curve is quadratic or a measured piecewise-linear table.
Its exact integral mean remains well defined when the two volumes are equal. Constant
discharge through a unit produces `P = 0.00981 Q H η(Q,H)` MW, with quadratic or tabulated
efficiency and quadratic penstock losses. Unit commitment enforces either zero discharge
or the stated on-state discharge, head and power ranges. Startup costs and changes in stored and in-transit water value enter the revenue
objective. Initial inventory value is subtracted, matching the Julia convention.

Tunnel flow is signed. Open tunnels satisfy
`(H_source − H_target) opening = loss_factor Q |Q|` without a gate-loss curve.
`opening` is a head multiplier in this empirical law; it is not an area fraction.
With `gate_loss_factor_curve`, the law is
`H_source − H_target = (loss_factor + G(opening)) Q |Q|`.
The supplied table adds resistance in s²/m⁵, must cover positions [0,1], and must
have nonnegative finite values. Preparation evaluates the prescribed position
once; no gate-table lookup occurs in nonlinear callbacks.
Pressure validation reports the unscaled head error in metres, so a small opening
cannot hide a large hydraulic error. A closed tunnel has zero flow and imposes no pressure relation. Pressure
junctions conserve water without storage. Junction heads are independent
unknowns. Nonconstant reservoir mean heads and turbine net heads/powers have local
constitutive equations. All machines use the
same net-head and electrical-power unknowns with local work balances.

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
Measured storage means and turbine production supply local derivatives through
the same sparse nonlinear interface. Native Ipopt solves continuous dispatch. SCIP proposes on/off schedules using
linear generation estimates and the same water balances and river transport.
Nonlinear dispatch evaluates those proposals, and search can change earlier
commitments. The proposal objective is **not a global bound**. POUNCE's optimizer
is exercised in comparison tests, rather than exposed as a user option.

`validation` evaluates water balance, pressure losses and operating limits in
physical units. It checks this interval model, not the finer replay or global
optimality provided by the complete Julia package. A feasible result alone says
nothing about optimality. `solver_status` is Ipopt's termination code;
`iterations` describes the final nonlinear solve, not the entire commitment search.

Still required are adaptive hydraulic integration and
independent replay, tunnel interlocks,
soft constraints and reserves, restart/rolling horizons and
global bounds. Complete synthetic watercourses and Tokke–Vinje must exercise
those capabilities before the Rust implementation can replace Julia. Small
examples establish correctness of their equations; they cannot settle the final
speed and solution comparison.


## Plants and pumping

Create a `Plant::new(name, source, target)` and attach each `Unit::new(...)` by
its plant name. Units at a plant share intake, outlet level, optional outlet-head
floor and a tailwater curve of aggregate discharge. Plant flow/power bounds,
ramps and minimum operating times apply to the aggregate. A plant is on when
any of its units is on. Initial aggregate state follows the member units; supply
aggregate output history when their historical observation intervals differ.

Set `unit.mode = UnitMode::Pump` for pumping. Water endpoints always follow
water motion. Pump consumption and generator output are positive magnitudes in
`unit_power` and `plant_power`; the objective subtracts pump consumption. Pump
working head includes the elevation lift, shared tailwater and pipe losses.
`motor_efficiency` applies to pumps. `CommitmentExclusion` permits at most one
named group to operate, while allowing parallel units within the same group.

A variable-flow pump uses `Efficiency::Quadratic` or `Efficiency::Table`.
A fixed-speed pump instead supplies `PumpCharacteristic::new(heads, discharge,
efficiency, head_extrapolation)`, with one discharge and efficiency at each head.
The hydraulic solution determines its discharge; commitment decides on/off.
The single-head policy must explicitly be `Extrapolation::Constant`.

All machines satisfy local work balances. For pumping,
`0.00981 Q H = P η_hydraulic η_motor`; for generation,
`P = 0.00981 Q H η_hydraulic η_electrical(P)`.
The optional `generator_efficiency_curve` maps electrical output to efficiency;
its required shaft power must strictly increase. These equations avoid dividing
by efficiency or inverting the electrical curve during nonlinear trial steps.
Power limits and ramps remain bounds and linear rows. Small local Hessians
are differentiated with fixed-size arithmetic and assembled into the sparse
whole-system Hessian. Gate and hydraulic equations retain their own curvature.


## Time-dependent operating inputs

`ScheduleCase::operations` holds `OperationalSeries` records. An explicit
`Object::Unit`, `Object::Plant`, `Object::Reservoir`, `Object::Tunnel` or
`Object::River` reference identifies its target. `Attribute` names the physical
quantity. Times are absolute hours, values are held until the next knot, and
static data apply before the first knot. The control grid must include every
operating change inside the horizon. Series may retain data outside the horizon.

Preparation resolves these records into numerical arrays. Bounds can tighten
the static domain; they cannot extend a physical curve's domain. Unit/plant
inputs include production and discharge bounds/schedules, maintenance, forced
on/off states, costs and directional ramps. Reservoir inputs include inflows,
storage bounds and volume/level ramps. Tunnel openings and capacities, and river
inflows/capacities, can also vary in time. Continuous reservoir storage satisfies
both sides of a bound change at a scheduling vertex. Level ramps use actual
vertex levels, not interval mean levels.

A rate limit is integrated over the distance between observation midpoints.
If part of that interval is unrestricted, the comparison stays unrestricted.
The first comparison needs a supplied historical interval average. Plant flow
and power can inherit unit observations when their intervals match; explicit
aggregate observations use the plant's declared interval. Numerical evaluation
uses the prepared arrays without name lookup or input interpretation.

## Reservoir outlets and level-difference channels

`River::law` selects a physical outlet: `Controlled`, `Orifice`, `Weir`,
`Level`, or `Delta`. Ordinary outlets optimize a gate fraction held over each
control interval. `River::gate` and time-dependent `GateMin`/`GateMax` bound
that fraction. A weir has no gate control. Natural inflow enters after the
reservoir withdrawal and is included in downstream routing and capacity checks.

Orifices satisfy `Q = C a sqrt(max(H − crest, 0))`; weirs satisfy
`Q = C max(H − crest, 0)^(3/2)`. Possible drying requires `allow_dry`.
The nonlinear equations use the square root of wet depth as an unknown rather
than differentiating a square root at zero. A `Level` table supplies full-open
discharge as a nonnegative, nondecreasing function of reservoir mean head.

`LevelFlowCurve::new(reference_levels, curves)` defines a passive signed channel
response. Each curve maps upstream minus downstream level to discharge and
must be nondecreasing and contain `(0, 0)`. Curves are interpolated in either
upstream or downstream reference level; outside that reference range the nearest
curve applies. Channel equations extend each delta-level curve by its end
secants, while raw `value` lookups reject delta levels outside the supplied
domain. A channel connects two reservoirs immediately, with signed discharge,
zero delay, and no gate, natural inflow or historical inventory.


## Distributed river travel times

`RiverTransport::Fixed { delay }` translates releases exactly in hours.
`RiverTransport::Distributed { curves }` supplies `DelayCurve` records containing
a reference flow, ordered delay-bin edges and bin probabilities. Probabilities
must sum to one; they are never silently normalized. A single curve is independent
of flow. Several curves must cover zero through the reach capacity, and their
probability distributions interpolate between neighboring reference flows.

Each release retains its original reference blend. The model forms reference
flows once per release and reuses them for every arrival and inventory
observation. The last reference flow is the total minus the others, preserving
conservation algebraically. Nonlinear derivative assembly only sees the small
piecewise-quadratic reference-flow functions.

Distributed river junctions pass interval-average arrivals into downstream
reaches on a separate transport mesh. `ScheduleCase::routing_grid` can supply
that mesh; it must include all control vertices. The default has four cells per
control interval. Refinement changes arrival timing at river junctions and must
be checked for convergence; a small equation residual alone does not establish
transport accuracy. Networks containing only fixed delays retain exact shifted
release boundaries and require no mixing mesh.

`arrival_volume` and `remaining_volume` expose interval transport calculations
in m³. Arrival fractions integrate the nonnegative pieces of the trapezoidal
arrival density directly. Remaining inventory integrates survival probability
directly, avoiding cancellation when only a small tail remains in transit.


## Hydraulic operating limits

`Unit::minimum_efficiency` restricts the machine's actual hydraulic efficiency;
it does not replace low curve values with a floor. The same restriction applies
to generating, variable-speed pumping and fixed-speed characteristics.

`gross_head_min` and `gross_head_max` constrain the surface-level difference
before machine and tunnel losses. `startup_gross_head_min` applies at an
off-to-on transition and uses the surface levels at that control edge. Running
limits use the represented interval mean levels. Pumping reverses the sign so
that its gross head is a positive lift. These limits also accept time-dependent
`GrossHeadMin`, `GrossHeadMax` and `StartupGrossHeadMin` series on units.

Plant endpoints supply the default gross-head references. If an endpoint is a
pressure junction, supply `Plant::gross_head_source` and/or `gross_head_target`
explicitly: gross-head references must be reservoir surfaces or fixed boundaries.
They may differ from machine head endpoints. Startup table levels are evaluated
at the actual starting storage, without smoothing or interpolating in time.


## Flow observations and water allowances

`FlowRequirement::new(name)` defines a sum of unit discharges, river releases,
river arrivals and signed tunnel flows. Add it to `ScheduleCase::flow_requirements`.
An observation introduces no water or hydraulic node. Local `inflow` augments
the observed rate only. River arrivals retain their control-interval average
meaning when a window cuts an interval.

`min_flow`/`max_flow` bound rates in m³/s. `min_volume`/`max_volume` bound the
horizon's integrated water in m³, including `initial_volume` carried from an
earlier schedule. `VolumeWindow` bounds water over its own absolute-hour window
and counts historical water exactly once. `RollingRamp` bounds the difference
from the preceding `hours`-long average, including both sides of a jump.
`PeriodRamp` compares time-weighted averages of consecutive, possibly unequal
periods. Ramp limits are differences in m³/s, rather than rates per hour.

Windows crossing the horizon start require complete `history_grid`/`history_flow`
data. History records interval rates and must end at the horizon start. Windows
are never shortened to fit available history. Active volume windows and period
ramps must finish within the horizon. `Object::FlowRequirement(name)` accepts
time-dependent `FlowMin`, `FlowMax` and `Inflow` inputs.

Preparation integrates interval overlaps into sparse linear equations. The same
rows constrain commitment proposals and nonlinear dispatch. Results expose
`observed_flow` and `observed_volume` in requirement order; the latter includes
the previously counted water at its initial vertex.


## Reservoir level ranges

`ReservoirOperation::level_windows` contains `LevelWindow` records with absolute
start/end hours and a maximum level range in metres. Window endpoints inside
the horizon must be control vertices. A window ending after the horizon
constrains the known prefix; a window starting earlier requires complete
`StorageHistory` in `ReservoirOperation::history`. Historical storage is linearly
interpolated in time before applying the level–storage curve. Its final record
must equal the supplied initial storage.

One unknown lower level bounds all surfaces in a window. The equations include
river-arrival vertices between control times, rather than checking control
endpoints alone. Monotonicity of the level–storage curve makes those vertices
sufficient for the represented piecewise-linear storage path. Physical validation
independently calculates the largest minus smallest level, including history.
A zero range therefore enforces a constant surface without changing the curve.


## Tunnel outfalls

`Tunnel::new(name, source, target, capacity, loss_factor)` constructs an open
pressure connection. Setting `discharge_river` sends its water into that reach
while `target` remains the pressure reference. Outfall flow is nonnegative in
the declared direction; reverse flow cannot become a negative river release.
The source loses water once, and the receiving reservoir gains it through the
reach's arrival response. River inventory accounts for the intervening water.
A closed outfall tunnel still disconnects pressure. Ordinary tunnels remain
bidirectional. Flow observations can refer to the tunnel or its downstream
river arrivals explicitly.


## Discrete tunnel gates

A nonempty `Tunnel::gate_positions` lets optimization select among its ordered
positions in `[0, 1]`. Static and initial openings must be listed positions.
An `Opening` time series is a hard schedule after its first knot; preceding
intervals retain the full declared choice. Empty `gate_positions` keeps the
opening prescribed, including arbitrary continuous scheduled values.

SCIP proposes unit states and gate positions together. Pressure tangents
distinguish positions hydraulically, with inactive-row bounds derived from the
affine head and flow domains. Nonlinear dispatch solves the original pressure
law at the selected position. Search exclusions include both unit and gate
decisions, so fixing all units does not lock free gates.

`gate_adjustment_cost` charges each absolute position change, including the
change from `initial_opening` (or the static opening if omitted). Public results
include physical `tunnel_opening` and total `gate_adjustment_cost`.

Every positive position is prepared as a resistance in the equation
`H_source − H_target = resistance Q |Q|`. For the head-multiplier convention this
resistance is `loss_factor/opening`; for a gate-loss curve it is
`loss_factor + G(opening)`. This algebraic normalization keeps pressure residuals
in metres even at a small opening. Zero remains a disconnected gate.
