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

`solve` searches unit commitment, eligible tunnel positions, enabled river flow states and nonlinear
dispatch. Feasible-start preparation uses SCIP MILP proposals and Ipopt dispatch;
SCIP searches the full nonlinear model and supplies its upper bound.
`preparation_time_limit=0` disables preparation. An accepted `initial` schedule
skips it and supplies a checked start without fixing decisions.
`fixed_u` and `fixed_tunnel_opening` restrict the search and certificate to those
supplied decisions.

`schedule_case` runs preparation on its own. It dispatches each MILP proposal
immediately and uses the best accepted trajectory on the proposal grid to update the next proposal
and nonlinear starting point. Rejected candidates do not replace the accepted
schedule. Without an accepted trajectory, a short feasibility-only dispatch
precedes revenue optimization. Proposals with an accepted reference search joint neighborhoods of
unit states, river states and categorical gate positions around the best
accepted trajectory. A supplied `initial` must pass independent acceptance for this case;
it is not an unchecked guess for changed physical inputs. All proposal and
dispatch attempts share `time_limit`, with individual limits set by
`proposal_time_limit` and `nlp_time_limit`.
For a forecast update, `schedule_case(case; warm=previous_iterate)` instead
uses an unchecked hydraulic reference and numerical starting point. The ordinary
proposal search selects discrete states under the updated inputs.
The start must have compatible arrays and the same grid; it need not satisfy
the updated physical or operating inputs. It is never retained as a feasible
schedule without solving and independent acceptance. `initial` takes precedence
when both are supplied.
Distributed-delay proposals use the accepted trajectory's entry-cohort flows
as their reference. Without a trajectory, they use each reach's last recorded
flow, or zero when no previous flow is known. The proposal freezes only the
delay mixture; nonlinear dispatch retains the flow-dependent response.

`solve_case` fixes river operating states through `river_on[river, interval]`;
a warm schedule supplies omitted states, otherwise enabled rivers default to on.
See [river operating costs and head restrictions](operating-rules.md).

`LevelWindow`, `VolumeWindow`, `RollingRamp`, and `PeriodRamp` specify
[restrictions over time windows](operating-rules.md#restrictions-over-time-windows).
Attach level windows to a reservoir and flow windows to a `FlowRequirement`.
Their absolute times, history and physical time weights survive restart and
refinement. Rolling/period flow limits compare averages of the observation
record; reservoir windows bound the complete recorded level range.

Both `solve` and `solve_case` return `SolveResult`. Its stable fields are
`status`, `accepted`, `solution`, `global_certificate`, `global_bound`,
`relative_gap`, `absolute_gap`, `feasible_lower_bound` and `total_seconds`;
`objective` is derived from the returned schedule. `solution` may be `nothing`.
Detailed work counters and validation reports are in `diagnostics`.
A time limit may leave an accepted schedule and an unfinished gap calculation.
If preparation or initial reconstruction exhausts the allowance, `solve` skips
global construction. It retains an independently checked incumbent when one is
available, subject to the requested fixed controls and replay policy. The result
reports `budget_exhausted_phase`, `total_seconds` and `budget_overrun_seconds`;
its global bound and gaps remain unavailable when search did not supply a bound.
Time budgets are checked between phases. Model construction, a solver iteration,
physical validation and replay can overrun the allowance.

An audited initial schedule is retained independently of SCIP. After optimizer
attachment, `native_start_audit` checks that every native variable has a finite
start, including any objective-bridge slack. That slack uses the actual scaled
objective, including constants. This is a completeness check, not a new physical
feasibility proof. `scip_diagnostics.native_incumbent_available` and
`native_incumbent_objective` report SCIP's own incumbent separately; retaining a
preparation schedule does not establish that SCIP accepted it.

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

Replay reports separate arrival and discharge shortfalls (m³/s) from power
shortfalls (MW). These values describe constraint exceedances; they do not
set solver margins. After failed replay, `solve_verified` may add finer
hydraulic equations with the same unit discharges, outlet openings, commitments,
tunnel positions and reserve allocations. Fixed-speed pump flows remain hydraulic
states. Schedules and ramps use the original operating windows. The case, returned
schedule and objective retain the requested grid. `refine_grid` remains accepted
for compatibility; both values now preserve the control grid.

The retry enforces both discretizations. It can be more restrictive than either
mesh alone and is not a continuous-time error certificate. `max_refinements`
bounds this additional work. An explicitly requested `operational_margin` is a
power inset in MW, capped at half the permitted power range so a fixed-power
unit remains well-defined. It is not converted into a discharge margin.

`validate` checks supplied arrays against the equations, restrictions and river
transport directly. It does not repeat the same-grid hydraulic solve.
`replay_audit` independently integrates on a finer mesh and retains the separate
transport check. Neither accepts stored validation flags as evidence. Unexpected
programming errors in replay or reconstruction are raised to the caller.

Replay compares storage and endpoint level at common grid vertices. By default,
each storage difference must be below both 2,000 m³ and 0.1% of the reservoir's
usable volume plus one litre. Each level difference must be below 1 cm.
The small absolute allowance keeps the volume comparison meaningful for a fixed
or nearly zero storage range. Set `storage_convergence` (Mm³),
`storage_relative_convergence` (capacity fraction), and `level_convergence` (m)
on `replay_audit` or `solve_verified` to request different accuracy. An infinite
tolerance disables that comparison. The returned per-reservoir storage and level
differences follow reservoir order. Passing these comparisons estimates mesh
convergence; it does not bound the error against a continuous-time solution.

A physically accepted schedule can be returned even if its stricter global-model
start check fails. In that case it supplies neither a `feasible_lower_bound` nor
an optimality gap or certificate. Its physical acceptance and its eligibility
for a global bound are separate facts.

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

The independent hydraulic solve uses automatic derivatives of its equations.
Its `tolerance` applies directly to algebraic laws. Reservoir water balances
also allow two floating-point spacings of the endpoint storage, converted to
m³/s by dividing by `0.0036 * duration`. This accounts for rounding a small
volume change into a large stored volume. `residual` reports the raw maximum
equation error, so it can exceed `tolerance` when `converged` is true. Operating
bounds and time-discretization checks remain separate.

`ScheduleCase(routing_grid=...)` optionally supplies a river transport mesh
containing every control boundary. Distributed confluences otherwise use four
transport cells per control interval. This mesh adds no controls or reservoir
states; schedule arrays and operating windows retain `case.grid`.
`transport_audit` separately checks timing convergence and operating arrival
and capacity limits on the independently routed waves, including peaks between
control boundaries. See [river networks](river-networks.md).

`head(reservoir, volume)` returns the surface elevation. The interval equations
instead use `mean_head(reservoir, before, after)`, the volume-averaged elevation
and discrete gradient of stored potential. This method integrates polynomial
and piecewise-linear data, including crossings of table knots and equal endpoint
volumes. The same equation is used by simulation, dispatch, global optimization
and validation; see [hydraulic equations](hydraulics.md).

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

`penalty_settings(case)` returns effective `soft`, river `release` and `reserve`
penalty equations plus an explicit `defaults` tuple. Entries retain original
typed rules where present and expose physical/coefficient units, cost factors,
bounded shortfalls, activity and mathematical origins. See [effective penalties](penalties.md)
for defaults, original-window accrual and the distinction from public SHOP costs.

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
The [rolling-horizon example](../examples/rolling_horizon.jl) builds the next
forecast window from the restarted system and observation histories. It extends
the horizon, changes prices and inflows, and shifts the previous iterate into
`warm`. The executed state comes from `restart_case`, not from that guess.

## Spinning capacity

Supply `reserve_groups=[ReserveGroup(...)]` on `ScheduleCase`.
`solve` and `solve_case` retain their ordinary `SolveResult`; capacity is in
`result.solution.reserves::ReserveAllocation` when capacity contracts are present. See [reserve equations and scope](reserves.md).

## Water values and objective accounting

`energy_conversion_factor(head, efficiency)` gives MWh/Mm³ at explicit reference
conditions. `volume_water_value(energy_value, factor; downstream_value=0)` gives
currency/Mm³. `objective_components(case, schedule)` recomputes the economic
accounting. See [terminal water values](terminal-value.md) for cascade,
historical-transit and currency conventions.
