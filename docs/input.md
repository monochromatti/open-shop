# Input profile

OpenSHOP reads a JSON case with `readcase(path)` or a dictionary with `case_from_dict(data)`. `case_dict(case)` produces the public profile and `writejson(path, data)` writes it. Object names connect the watercourse; JSON arrays do not determine hydraulic connectivity.

The top-level fields are `schema_version` (use `2`), `name`, `grid`, `prices`, `reservoirs`, `junctions`, `boundaries`, `tunnels`, `plants`, `generators`, `river_junctions`, `rivers`, `operations`, and optional `flow_requirements`, `tunnel_interlocks`, `commitment_exclusions` and `soft_constraints`. Component arrays may be empty. `grid` contains strictly increasing finite time edges; `prices` has one value per interval. Matrices are stored as arrays of rows.

## Objects

The fields below are the accepted object attributes. Fields marked “required” have no constructor default. All numeric data must be finite where the profile requires physical bounds.

| Object | Required fields | Optional fields |
| --- | --- | --- |
| Reservoir | `name`, `v0`, `vmin`, `vmax`, `water_value` | `z0`, `slope`, `curvature`, `inflow`, `level_curve`, `volume_ramp_up`, `volume_ramp_down`, `level_ramp_up`, `level_ramp_down` |
| Hydraulic junction | `name` | `hmin`, `hmax` |
| Boundary | `name`, `head` | None |
| Tunnel | `name`, `source`, `target`, `resistance`, `capacity` | `opening`, `gate_positions`, `initial_opening`, `gate_adjustment_cost`, `gate_loss_factor_curve`, `discharge_river` |
| Tunnel interlock | `name`, `tunnels` | `min_open`, `max_open`, `forward_tunnels`, `backward_tunnels`, `plants`, `forward_switch_time`, `backward_switch_time`, `initial_forward_wait`, `initial_backward_wait` |
| Plant | `name`, `source`, `target`, `pmax` | `discharge_river`, `ramp`, `initial_power`, `initial_interval_hours`, `tailwater_curve`, `outlet_head_floor`, `pmin`, `qmin`, `qmax`, `ramp_up`, `ramp_down`, `discharge_ramp_up`, `discharge_ramp_down`, `initial_discharge`, `minup`, `mindown`, `initial_on`, `initial_age` |
| Unit (`Generator` or `Pump`) | `name`, `plant`, `qmin`, `qmax`, `pmin`, `pmax`, `hmin`, `hmax` | `mode`, `motor_efficiency`, `efficiency`, `min_efficiency`, `qbest`, `qcurvature`, `hbest`, `hcurvature`, `initial_on`, `initial_age`, `minup`, `mindown`, `startup`, `shutdown`, `turbine_table`, `generator_efficiency_curve`, `ramp_up`, `ramp_down`, `discharge_ramp_up`, `discharge_ramp_down`, `initial_power`, `initial_discharge`, `initial_interval_hours` |
| Commitment exclusion | `name`, `groups` | — |
| River junction | `name` | None |
| River | `name`, `target`, `curves`, `capacity`, `water_value` | `source`, `allocation`, `inflow`, `law`, `coefficient`, `crest`, `min_arrival`, `arrival_policy`, `arrival_window_grid`, `deterministic_delay`, `gate_min`, `history_grid`, `history_release`, `discharge_curve`, `delta_level_ref_up_flow_curve`, `delta_level_ref_down_flow_curve`, `allow_dry`, `ramp_up`, `ramp_down`, `initial_release`, `initial_interval_hours` |
| Operation | `object`, `attribute`, `times`, `values` | None |
| Flow requirement | `name` | `generators`, `rivers`, `inflow`, `min_flow` |
| Soft constraint | `name`, `object`, `quantity`, `time` | `end_time`, `lower`, `upper`, `max_lower_slack`, `max_upper_slack`, `lower_cost`, `upper_cost` |

Reservoir head is `z0 + slope * V + curvature * V²`, unless `level_curve` supplies the volume-to-head relationship. The relationship must increase over the storage domain. Hydraulic junctions have no storage; signed tunnel and generation flows satisfy continuity. Boundaries have fixed heads and represent external water exchange.

Positive tunnel discharge follows `source` to `target`; negative discharge reverses that direction. `opening` lies between zero and one, and zero closes the tunnel. Resistance must be positive. A plant's units share endpoints and aggregate capacity/ramp restrictions. Tailwater is additional head loss indexed by total plant discharge. Optional `outlet_head_floor` sets the outlet reference to the greater of the receiving node head and that floor. Generating head is source head minus this reference and loss; pumping lift reverses that expression, adding loss. The target receives the discharge immediately unless `discharge_river` redirects it into delayed transport. All units in a plant have the same `mode`, either `"generate"` (default) or `"pump"`; plant endpoints follow actual water motion. Unit discharge and power are positive production/consumption magnitudes. `motor_efficiency` is a constant pump efficiency, defaulting to one. See [pumping](pumping.md) for mode groups and supported scope. Explicit initial state, age and power support chronological dwell and ramp constraints. `startup` is charged once on an off-to-on transition; `shutdown` is charged once on an on-to-off transition, including the first interval relative to `initial_on`.

River targets may be reservoirs, boundaries, river junctions or other rivers. A boundary receives water leaving the modeled system; delayed arrivals, rather than upstream releases, determine that exchange.

`gate_positions` defaults to an empty array, leaving tunnel openings prescribed.
A nonempty sorted unique list enables discrete gate decisions. Positions lie in
`[0,1]`; static `opening` and explicit `initial_opening` must be listed. An
opening series prescribes positions after its first knot and leaves earlier
optimized intervals free. `TunnelInterlock` counts strictly positive member
positions in each interval, between `min_open` and `max_open`. Movement cost is
`gate_adjustment_cost * abs(opening - previous_opening)`, including the first
interval relative to `initial_opening` or the static opening. It has no duration
multiplier. See [discrete tunnel gates](api.md) for controls and solver
scope.

Optional interlock `forward_tunnels` and `backward_tunnels` identify disjoint
closure groups within `tunnels`. Closing a member keeps every unit of the listed
`plants` off for `forward_switch_time` or `backward_switch_time`, in hours.
Initial residual waits preserve pre-horizon closures; their histories must agree
with the protected units' off states and ages. See
[gate handovers](operating-rules.md#gate-handovers-and-plant-shutdowns) for the
equations, exact event semantics and continuation rules.

`min_efficiency` is a lower bound on a running unit's turbine efficiency and defaults to zero. It is independent of the electrical-efficiency curve.

River `law` is `controlled`, `orifice`, `weir`, `junction` or `delta`. The delta law uses supplied reference-level flow tables for a signed, zero-delay channel between reservoirs; see [channel and gate-loss curves](api.md). For controlled outlets, reservoir withdrawal is `(operating capacity - natural inflow) * gate`; total reach release adds natural inflow. Orifice release is `coefficient * gate * sqrt(head - crest)`. Weir release is `coefficient * (head - crest)^1.5` with an uncontrolled gate of one. `allow_dry=true` replaces negative head above crest by zero. A discharge table replaces the analytic outlet relationship. Environmental arrivals use `interval_average` or `pointwise`; supplied `arrival_window_grid` declares the windows for interval-average requirements.

## Shared streams and direct river links

Reservoir, hydraulic-junction and boundary names identify physical endpoints. Tunnel networks and instantaneous delta-level channels may contain loops. Delayed river networks must be acyclic. A zero-storage river node can have several incoming and outgoing reaches; prescribed allocation fractions determine how it distributes the arrivals. Each node needs a water supply and at least one outgoing reach. Hydraulic junctions and river nodes have different roles.

A river can name another river as its `target`. The receiving river uses `law="junction"` and omits `source` or sets it to `"auto"`. All rivers targeting that receiving river share an inferred confluence. A shared source/target label distinct from every object name also declares an internal river node; no `RiverJunction` object is required. OpenSHOP restores direct river links and node labels on serialization. Names beginning with `__river_merge__` or `__river_node__` are reserved. An explicit source on a receiving direct-link river is rejected.

Plants and tunnels can set `discharge_river` to a receiving reach while retaining
`target` as their hydraulic head reference. Routed tunnel outfalls are forward
only; other tunnel links remain reversible. Receiving reaches combine these
outfalls, upstream rivers and their own natural `inflow`. An inflow-only reach
also omits `source`. See [river networks](river-networks.md) for equations,
capacity conventions, history and mixed-delay behavior.

A shared penstock can be represented by a tunnel from a reservoir to a hydraulic junction, with multiple plants taking their source at that junction. Junction continuity sums their unit discharges, and the common tunnel head-loss law therefore depends on the combined flow. Units within one plant likewise share that plant's hydraulic endpoints and tailwater relationship.

At each river node, outgoing reaches use `law="junction"` and finite `allocation` fractions in `[0,1]` summing to one. The default is one. The sum is checked with `1e-12` absolute and relative tolerance; shares are never inferred or normalized. Fractions apply only to inherited upstream arrivals; natural inflow and plant/tunnel outfalls enter the named branch afterward. A branch's historical release already includes its allocation. Fractions are constant over the horizon, not optimization decisions or operating-series attributes. Reservoir outlets and instantaneous channels use `allocation=1`. Separate controlled reservoir outlets remain independent decisions.

## Tables and transport

A one-dimensional table contains `x` and `y` arrays, with strictly increasing coordinates. Level tables map storage to head; tailwater tables map aggregate discharge to additional head; discharge tables map source level to ungated outlet flow; generator-efficiency tables map electrical power to efficiency. Coverage and monotonicity requirements are validated according to the physical role.

A turbine table contains `heads`, `discharge`, `efficiency`, `qmin`, `qmax`, and optional `interpolation` and `head_extrapolation`. Efficiency rows correspond to discharge coordinates and columns to head coordinates. Head-indexed `qmin` and `qmax` describe the operating discharge envelope. Interpolation is `bilinear` or `pchip_discharge`; the latter is cubic in discharge and linear in head. Turbine and electrical efficiencies are separate.

`head_extrapolation` defaults to `"error"`: reference heads must cover the generator's declared operating head range. Explicit `"linear"` extends efficiency and discharge envelopes using the outermost pair of reference heads. Supplied knots remain unchanged; the generator's head, discharge, power and efficiency restrictions still apply to the extended values. Raw table-query APIs remain strict unless the caller explicitly requests linear extrapolation.

Each distributed delay curve contains `reference_flow`, `edges`, and `weights`. Edges are finite increasing nonnegative delays; weights are nonnegative bin probabilities summing to one. One distributed curve is independent of flow. Multiple curves have strictly increasing reference flows covering zero through `capacity`; their neighboring transfer distributions are blended using contemporaneous release flow. A deterministic reach instead supplies nonnegative `deterministic_delay`; `curves` may be an empty array. History edges end at the scheduling horizon's first edge and have one nonnegative release value per interval. History is physical water already in transit, not an optimization decision.

Deterministic transport preserves release cohorts through confluences. A network containing distributed reaches uses interval-average mixing at downstream junctions. Its scheduling grid is therefore part of the model; numerical refinement does not turn its discrete global bound into a continuous-time certificate.

## Operating series

An operation names an existing object and one supported attribute:

| Object | Attributes |
| --- | --- |
| Reservoir | `inflow`, `vmin`, `vmax`, `volume_ramp_up`, `volume_ramp_down`, `level_ramp_up`, `level_ramp_down` |
| Generator | `qmin`, `qmax`, `pmin`, `pmax`, `forced_on`, `maintenance`, `power`, `discharge`, `startup`, `shutdown`, `ramp_up`, `ramp_down`, `discharge_ramp_up`, `discharge_ramp_down` |
| Plant | `pmin`, `pmax`, `qmin`, `qmax`, `forced_on`, `maintenance`, `power`, `discharge`, `ramp_up`, `ramp_down`, `discharge_ramp_up`, `discharge_ramp_down` |
| Tunnel | `capacity`, `opening`, `gate_adjustment_cost` |
| Tunnel interlock | `min_open`, `max_open` |
| River | `inflow`, `capacity`, `gate_min`, `gate_max`, `min_arrival`, `min_release`, `release_penalty`, `ramp_up`, `ramp_down` |
| Flow requirement | `inflow`, `min_flow` |

`times` are strictly increasing absolute-hour knots and `values` are piecewise constant. Before the first knot the object's default applies; the last value is held afterward. Ordinary dispatch interval data use time averages. Storage edges and pointwise arrival requirements retain their corresponding endpoint/event semantics. Place knots on scheduling edges when a commitment or outage transition must occur at a particular time.

Operating minimum/maximum restrictions tighten their static object bounds. Gate and opening values lie in `[0,1]`; `forced_on` is `-1` for free commitment, `0` for off, or `1` for on. Duplicate series for one object/attribute are rejected. `min_release` is a hard requirement unless a positive `release_penalty` declares a soft shortfall cost. Soft shortfalls are measured in Mm³ and reported explicitly; they are not silently dropped.

`maintenance=1` makes a unit unavailable; a plant outage makes every member unit unavailable. Forced states, schedules and initial minimum-duration restrictions are intersected: one input cannot overwrite another. `power` and `discharge` prescribe hard interval schedules. They are inactive before their first knot; an explicit zero schedule forces the affected unit or plant off. Positive schedules require an operating unit, while positive `pmin`/`qmin` alone do not force operation.

Plant `qmax` defaults to no additional aggregate discharge restriction. Directional plant power ramps inherit `ramp` when omitted; other omitted ramp limits are disabled. See [operating rules](operating-rules.md) for timing and history conventions.

## Flow observations

A `flow_requirements` entry records a hard minimum at an operating measurement point. Its `generators` list names unit discharges; its `rivers` list names river releases. The constraint is

```math
\sum_{g\in G_o} q_{g,t}+\sum_{r\in R_o} Q^{release}_{r,t}+I_{o,t}\ge Q^{min}_{o,t}.
```

All terms are in m³/s. `inflow` and `min_flow` default to zero and accept time-dependent operations. The measurement adds no hydraulic node or water balance term. The physical network must already account for the observed water. Contributions must be unique existing names, and observation names must be globally unique. River contributions denote releases, not delayed arrivals; constrain a River’s `min_arrival` when travel time matters.

Refinement retains these rules and their physical input times. Independent equation validation and finer replay check the aggregate minimum. Flow requirements are hard constraints; conditional waivers and soft shortfall penalties are not supplied by this observation profile.

## Explicit penalty limits

`soft_constraints` contains separate, named `SoftConstraint` records. They do not
change an existing `OperationalSeries`. An empty array enables no additional
softening. A reservoir accepts `quantity="volume"` or `"level"` at an absolute
state assessment `time`; omit `end_time`. A generator or plant accepts
`quantity="power"` or `"discharge"` averaged over the original absolute window
`[time,end_time]`. Plant outputs sum every member unit before averaging.
Use `lower=upper` for a target schedule, or one or both sides for a band.

The nonnegative `max_lower_slack` and `max_upper_slack` bound the admitted
shortfall and excess; both default to zero. `lower_cost` and `upper_cost` default
to zero. A positive slack allowance requires a strictly positive cost. An absent
limit must have zero allowance and cost. Limits, allowances, costs and their
bounded monetary products must be finite. Volume, power and discharge limits
are nonnegative; level limits may use any finite datum. Names are globally
unique across physical objects, flow observations and soft constraints.

| Quantity | Slack unit | Penalty coefficient | Cost |
| --- | --- | --- | --- |
| Reservoir volume | Mm³ | currency/Mm³ | coefficient × point slack |
| Reservoir level | metres | currency/metre | coefficient × point slack |
| Unit or plant power | MW | currency/MWh | coefficient × window hours × average slack |
| Unit or plant discharge | m³/s | currency/Mm³ | coefficient × 0.0036 × window hours × average slack |

All active assessment times must be grid edges. Records wholly outside the
current horizon are inactive, allowing a forecast prefix to retain its original
future terminal target. A rate window crossing a horizon boundary is rejected.
State bands are restrictions at their declared observations; they do not impose
continuous-time bands between observations. Supply a record at each desired
state observation or output window. See [penalty equations and continuation](operating-rules.md#explicit-penalty-limits)
and the [example](../examples/soft_limits.jl).

## Units and errors

| Quantity | Unit |
| --- | --- |
| Grid, history, delay, dwell, initial age, operating times, switch times and residual waits | hours |
| Storage and transit | Mm³ |
| Discharge and inflow | m³/s |
| Head, crest and elevation | metres |
| Electrical power | MW |
| Power ramp | MW/hour |
| Discharge or release ramp | (m³/s)/hour |
| Storage ramp | Mm³/hour |
| Level ramp | metres/hour |
| Price | objective currency/MWh |
| Startup and shutdown | objective currency |
| Tunnel gate adjustment coefficient | objective currency per unit normalized position |
| Water value and release penalty | objective currency/Mm³ |
| Efficiency, gate, opening, allocation and delay probability | dimensionless |

Unknown fields at the case, object, table, delay-curve or operation level raise an error. Missing required fields, invalid bounds, unsupported laws/interpolation, unknown endpoints, duplicate names/series and invalid topology also fail input validation. Inputs use this documented profile; physical calibration must be supplied explicitly. Keep the case and returned schedule together: changing the grid changes schedule dimensions and the discrete optimization problem.
