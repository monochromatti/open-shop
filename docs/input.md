# Input profile

OpenSHOP reads a JSON case with `readcase(path)` or a dictionary with `case_from_dict(data)`. `case_dict(case)` produces the public profile and `writejson(path, data)` writes it. Object names connect the watercourse; JSON arrays do not determine hydraulic connectivity.

The top-level fields are `schema_version` (use `2`), `name`, `grid`, `prices`, `reservoirs`, `junctions`, `boundaries`, `tunnels`, `plants`, `generators`, `river_junctions`, `rivers`, `operations`, and optional `flow_requirements`. Component arrays may be empty. `grid` contains strictly increasing finite time edges; `prices` has one value per interval. Matrices are stored as arrays of rows.

## Objects

The fields below are the accepted object attributes. Fields marked “required” have no constructor default. All numeric data must be finite where the profile requires physical bounds.

| Object | Required fields | Optional fields |
| --- | --- | --- |
| Reservoir | `name`, `v0`, `vmin`, `vmax`, `water_value` | `z0`, `slope`, `curvature`, `inflow`, `level_curve` |
| Hydraulic junction | `name` | `hmin`, `hmax` |
| Boundary | `name`, `head` | None |
| Tunnel | `name`, `source`, `target`, `resistance`, `capacity` | `opening` |
| Plant | `name`, `source`, `target`, `pmax` | `ramp`, `initial_power`, `initial_interval_hours`, `tailwater_curve`, `outlet_head_floor` |
| Generator | `name`, `plant`, `qmin`, `qmax`, `pmin`, `pmax`, `hmin`, `hmax` | `efficiency`, `min_efficiency`, `qbest`, `qcurvature`, `hbest`, `hcurvature`, `initial_on`, `initial_age`, `minup`, `mindown`, `startup`, `shutdown`, `turbine_table`, `generator_efficiency_curve` |
| River junction | `name` | None |
| River | `name`, `target`, `curves`, `capacity`, `water_value` | `source`, `law`, `coefficient`, `crest`, `min_arrival`, `arrival_policy`, `arrival_window_grid`, `deterministic_delay`, `gate_min`, `history_grid`, `history_release`, `discharge_curve`, `allow_dry` |
| Operation | `object`, `attribute`, `times`, `values` | None |
| Flow requirement | `name` | `generators`, `rivers`, `inflow`, `min_flow` |

Reservoir head is `z0 + slope * V + curvature * V²`, unless `level_curve` supplies the volume-to-head relationship. The relationship must increase over the storage domain. Hydraulic junctions have no storage; signed tunnel and generation flows satisfy continuity. Boundaries have fixed heads and represent external water exchange.

Positive tunnel discharge follows `source` to `target`; negative discharge reverses that direction. `opening` lies between zero and one, and zero closes the tunnel. Resistance must be positive. A plant's units share endpoints and aggregate capacity/ramp restrictions. Tailwater is additional head loss indexed by total plant discharge. Optional `outlet_head_floor` sets the turbine outlet reference to the greater of the receiving node head and that floor; net head is source head minus this reference and tailwater loss. The receiving node still receives the original water discharge. Generators produce electrical power; pumping is not part of this profile. Explicit initial state, age and power support chronological dwell and ramp constraints. `startup` is charged once on an off-to-on transition; `shutdown` is charged once on an on-to-off transition, including the first interval relative to `initial_on`.

River targets may be reservoirs, boundaries, river junctions or other rivers. A boundary receives water leaving the modeled system; delayed arrivals, rather than upstream releases, determine that exchange.

`min_efficiency` is a lower bound on a running unit's turbine efficiency and defaults to zero. It is independent of the electrical-efficiency curve.

River `law` is `controlled`, `orifice`, `weir` or `junction`. Controlled release is operating capacity times gate opening. Orifice release is `coefficient * gate * sqrt(head - crest)`. Weir release is `coefficient * (head - crest)^1.5` with an uncontrolled gate of one. `allow_dry=true` replaces negative head above crest by zero. A discharge table replaces the analytic outlet relationship. Environmental arrivals use `interval_average` or `pointwise`; supplied `arrival_window_grid` declares the windows for interval-average requirements.

## Shared streams and direct river links

Reservoir, hydraulic-junction and boundary names identify physical endpoints. Tunnel networks may contain loops. River networks must be acyclic. Multiple incoming rivers can merge into one outgoing river through a named river junction; conservation determines the downstream release. One river junction must have incoming reaches and exactly one outgoing reach. Hydraulic junctions and river junctions have different roles.

A river can name another river as its `target`. The receiving river uses `law="junction"` and omits `source` or sets it to `"auto"`. All rivers targeting that receiving river share an inferred confluence. OpenSHOP creates its internal junction during input construction and restores the direct river links on serialization. Names beginning with `__river_merge__` are reserved for inferred junctions. An explicit source on a receiving direct-link river is rejected.

A shared penstock can be represented by a tunnel from a reservoir to a hydraulic junction, with multiple plants taking their source at that junction. Junction continuity sums their unit discharges, and the common tunnel head-loss law therefore depends on the combined flow. Units within one plant likewise share that plant's hydraulic endpoints and tailwater relationship.

River splitting and allocation among multiple downstream branches are unsupported. Declare separate controlled releases from a reservoir when that physical arrangement applies; do not assume that a zero-storage confluence allocates flow.

## Tables and transport

A one-dimensional table contains `x` and `y` arrays, with strictly increasing coordinates. Level tables map storage to head; tailwater tables map aggregate discharge to additional head; discharge tables map source level to ungated outlet flow; generator-efficiency tables map electrical power to efficiency. Coverage and monotonicity requirements are validated according to the physical role.

A turbine table contains `heads`, `discharge`, `efficiency`, `qmin`, `qmax`, and optional `interpolation` and `head_extrapolation`. Efficiency rows correspond to discharge coordinates and columns to head coordinates. Head-indexed `qmin` and `qmax` describe the operating discharge envelope. Interpolation is `bilinear` or `pchip_discharge`; the latter is cubic in discharge and linear in head. Turbine and electrical efficiencies are separate.

`head_extrapolation` defaults to `"error"`: reference heads must cover the generator's declared operating head range. Explicit `"linear"` extends efficiency and discharge envelopes using the outermost pair of reference heads. Supplied knots remain unchanged; the generator's head, discharge, power and efficiency restrictions still apply to the extended values. Raw table-query APIs remain strict unless the caller explicitly requests linear extrapolation.

Each distributed delay curve contains `reference_flow`, `edges`, and `weights`. Edges are finite increasing nonnegative delays; weights are nonnegative bin probabilities summing to one. Distributed routing requires two curves, at zero and `capacity`, and blends their transfer distributions using contemporaneous release flow. A deterministic reach instead supplies nonnegative `deterministic_delay`; `curves` may be an empty array. History edges end at the scheduling horizon's first edge and have one nonnegative release value per interval. History is physical water already in transit, not an optimization decision.

Deterministic transport preserves release cohorts through confluences. A network containing distributed reaches uses interval-average mixing at downstream junctions. Its scheduling grid is therefore part of the model; numerical refinement does not turn its discrete global bound into a continuous-time certificate.

## Operating series

An operation names an existing object and one supported attribute:

| Object | Attributes |
| --- | --- |
| Reservoir | `inflow`, `vmin`, `vmax` |
| Generator | `qmin`, `qmax`, `pmin`, `pmax`, `forced_on` |
| Plant | `pmax` |
| Tunnel | `capacity`, `opening` |
| River | `capacity`, `gate_min`, `gate_max`, `min_arrival`, `min_release`, `release_penalty` |
| Flow requirement | `inflow`, `min_flow` |

`times` are strictly increasing absolute-hour knots and `values` are piecewise constant. Before the first knot the object's default applies; the last value is held afterward. Ordinary dispatch interval data use time averages. Storage edges and pointwise arrival requirements retain their corresponding endpoint/event semantics. Place knots on scheduling edges when a commitment or outage transition must occur at a particular time.

Operating minimum/maximum restrictions tighten their static object bounds. Gate and opening values lie in `[0,1]`; `forced_on` is `-1` for free commitment, `0` for off, or `1` for on. Duplicate series for one object/attribute are rejected. `min_release` is a hard requirement unless a positive `release_penalty` declares a soft shortfall cost. Soft shortfalls are measured in Mm³ and reported explicitly; they are not silently dropped.

## Flow observations

A `flow_requirements` entry records a hard minimum at an operating measurement point. Its `generators` list names unit discharges; its `rivers` list names river releases. The constraint is

```math
\sum_{g\in G_o} q_{g,t}+\sum_{r\in R_o} Q^{release}_{r,t}+I_{o,t}\ge Q^{min}_{o,t}.
```

All terms are in m³/s. `inflow` and `min_flow` default to zero and accept time-dependent operations. The measurement adds no hydraulic node or water balance term. The physical network must already account for the observed water. Contributions must be unique existing names, and observation names must be globally unique. River contributions denote releases, not delayed arrivals; constrain a River’s `min_arrival` when travel time matters.

Refinement retains these rules and their physical input times. Independent equation validation and finer replay check the aggregate minimum. Flow requirements are hard constraints; conditional waivers and soft shortfall penalties are not supplied by this observation profile.

## Units and errors

| Quantity | Unit |
| --- | --- |
| Grid, history, delay, dwell, initial age, operating times | hours |
| Storage and transit | Mm³ |
| Discharge and inflow | m³/s |
| Head, crest and elevation | metres |
| Electrical power | MW |
| Ramp | MW/hour |
| Price | objective currency/MWh |
| Startup and shutdown | objective currency |
| Water value and release penalty | objective currency/Mm³ |
| Efficiency, gate, opening and delay probability | dimensionless |

Unknown fields at the case, object, table, delay-curve or operation level raise an error. Missing required fields, invalid bounds, unsupported laws/interpolation, unknown endpoints, duplicate names/series and invalid topology also fail input validation. OpenSHOP does not import arbitrary proprietary SHOP files or infer physical calibration from names. Keep the case and returned schedule together: changing the grid changes schedule dimensions and the discrete optimization problem.
