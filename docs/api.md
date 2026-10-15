# Julia API

Construct named objects with keyword arguments, then pass a `HydroSystem` to
`ScheduleCase`. Absent component collections default to empty arrays. Case
construction resolves direct river links without modifying the supplied system.
See [input objects and units](input.md) for fields and JSON conventions.

## Solve and inspect

```julia
case = readcase("watercourse.json")
result = solve(case; time_limit=60.0, relative_gap=1e-3)
if result.accepted
    schedule = result.solution
    println(schedule.objective)
    writejson("schedule.json", schedule)
end
```

`solve` searches unit commitment, eligible tunnel positions and nonlinear
dispatch. Feasible-start preparation uses SCIP MILP proposals and Ipopt dispatch;
SCIP searches the full nonlinear model and supplies its upper bound.
`preparation_time_limit=0` disables preparation. An accepted `initial` schedule
skips it and supplies a checked start without fixing decisions.
`fixed_u` and `fixed_tunnel_opening` restrict the search and certificate to those
supplied decisions.

Both `solve` and `solve_case` return `SolveResult`. Its stable fields are
`status`, `accepted`, `solution`, `global_certificate`, `global_bound`,
`relative_gap`, `absolute_gap` and `total_seconds`; `objective` and
`feasible_lower_bound` are derived properties. `solution` may be `nothing`.
Detailed work counters and validation reports are in `diagnostics`.
A time limit may leave an accepted schedule and an unfinished gap calculation.

For `solve`, `accepted` includes finer chronological replay.
`global_certificate` records an enclosing objective bound within the
requested gap for the declared interval model, separately from finer acceptance. For `solve_case`, `accepted`
means original-grid equations and operating restrictions passed; use
`solve_verified` for finer replay. Local solver convergence is not a global
certificate. The [paper](../paper/paper.pdf) defines the model and certificate.

## Prescribed controls and replay

```julia
local_result = solve_case(case; u=commitment, tunnel_opening=positions)
schedule = dispatch_from_controls(case; commitment=commitment,
    generator_q=discharge, gate=outlet_opening, tunnel_opening=positions)
trajectory = simulate(case; generator_q=discharge,
    gate=outlet_opening, tunnel_opening=positions)
audit = validate(case, schedule)
finer = replay_audit(case, schedule)
```

Rows follow case object order; interval columns follow `case.grid`.
`Schedule.storage` has one column per grid vertex in Mm³. `head` is metres,
`power` is MW, and `generator_q`, `tunnel_q` and `river_release` are m³/s.
Both `Generator` and `Pump` use the common `case.system.generators` row order.
Unit discharge and power are nonnegative magnitudes; `net_power(case, schedule)`
returns interval production minus consumption. See [pumping](pumping.md) for
connections, efficiency data and reversible mode exclusions.
`arrival_volume` and `terminal_transit` are Mm³. `commitment` is binary;
`gate` contains river outlet openings. Returned arrays belong to the calculation;
input arrays are read-only. `simulate` returns `HydraulicTrajectory`, and
`forward_step` returns `HydraulicStep`. Physical convergence can coexist with
an operating-limit violation.

`SoftConstraint` declares a separate bounded, priced operating or terminal
restriction; see [penalty equations](operating-rules.md#explicit-penalty-limits).
`Schedule.soft_slacks` contains `SoftConstraintSlack` records with `name`,
`lower` shortfall, `upper` excess, `unit` and positive `penalty_cost`.
`soft_penalty_cost` is their total cost, already subtracted from `objective`.
These fields are absent from the JSON profile when the case has no soft rules.
`soft_constraint_assessment(case, schedule)` recomputes canonical minimal slacks
and returns `slacks`, `cost` and per-rule `excess` beyond allowed slack.
`validate` checks this accounting independently and reports allowed deviations
in `permitted_penalties`; they do not make a schedule invalid. Bare hydraulic
simulation returns the physical trajectory without economic accounting.

An empty `Tunnel.gate_positions` leaves openings prescribed by `opening` and
operating series. A nonempty sorted list enables discrete choices. Selectable
positions must be exact listed values; a closed gate has zero flow and independent
endpoint pressures. `TunnelInterlock` counts strictly positive positions.
`gate_adjustment_cost` charges absolute position changes, including the initial
change, without a duration multiplier. Independent reconstruction of selectable
gates requires an explicit opening matrix.

Interlocks can protect named plants with asymmetric shutdown waits after gate
closures. `validate_tunnel_openings(case, positions, commitment)` checks both
gate controls and these commitment restrictions; the two-argument form checks
gate controls alone. See [gate handovers](operating-rules.md#gate-handovers-and-plant-shutdowns)
for the timing and initial-history convention.

`gate_loss_factor_curve` supplies additional resistance versus opening. The
description of positions, coefficient units, dry outlets and elevation
thresholds is in [gate and outlet equations](hydraulics.md).
`River(law=:delta)` supplies exactly one upstream- or downstream-reference
`LevelFlowCurve`; its signed flow is instantaneous. Delta segments extend
linearly outside their knots, and the nearest reference curve is held outside
its reference-level range. Cover intended operation with calibrated curves.
See [river networks](river-networks.md) and [operating rules](operating-rules.md).

## Serialization and continuation

`case_dict`, `case_from_dict`, `readcase` and `writejson` provide the JSON boundary.
Physical schedule matrices serialize as rows, retaining keys such as `V`, `H`
and `u`. Solve results contain one nested `solution` profile.
Native results expose named fields and support existing string-key reads.
They do not support dictionary mutation. `schedule_dict(schedule)` returns an
independent editable profile; validate edited schedules before reuse.

`with_grid(case, edges)` changes the decision grid while retaining physical-time
inputs. `restart_case(case, schedule, time)` requires a source schedule that passes original-grid validation
and an existing grid edge. It carries storage, executed river cohorts, unit and
plant state ages, previous outputs and gate positions. `transport_audit` checks
convergence of distributed mixing separately from hydraulic integration.
Soft state observations and output assessment windows retain their original
times. Restart removes completed assessments and rejects a restart inside an
unfinished output window, whose partial history is not represented by this
profile. A terminal target is assessed only at its original time.

## Spinning capacity

Supply `reserve_groups=[ReserveGroup(...)]` on `ScheduleCase`.
`solve` and `solve_case` retain their ordinary `SolveResult`; capacity is in
`result.solution.reserves::ReserveAllocation` when capacity contracts are present. See [reserve equations and scope](reserves.md).
