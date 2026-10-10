# Operating rules

Operating inputs add ordinary equality and inequality constraints to the physical
model. Local nonlinear dispatch and global optimization use the same constraint
construction. Schedule validation checks the corresponding physical values.

## Availability, commitment and schedules

`OperationalSeries` accepts unit and plant `maintenance` (0 or 1) and `forced_on`
(-1 free, 0 off, 1 on). Plant state means at least one member unit is operating.
Maintenance at a plant makes every member unavailable. A forced-on plant retains
the choice of which units run. Restrictions are intersected, including residual
initial minimum-duration obligations. Direct conflicts raise an input error;
infeasible combinations over the future horizon remain an optimization result.

`power` and `discharge` series prescribe hard interval-average output, respectively
in MW and m³/s. Plant schedules constrain member sums. A series is inactive before
its first knot, and its final value is held. Zero output requires the affected
unit or plant to be off; positive output requires operation. Unit schedules must
be zero or within the current on-state operating range. Plant schedules must be
zero or within their aggregate bounds. Both power and discharge can be supplied;
their compatibility is determined by the hydraulic equations.

```julia
operations = [
    OperationalSeries(object=:UnitB, attribute=:maintenance,
        times=[0.0, 2.0], values=[1.0, 0.0]),
    OperationalSeries(object=:Plant, attribute=:discharge,
        times=[0.0, 2.0], values=[20.0, 30.0]),
    OperationalSeries(object=:Plant, attribute=:ramp_up,
        times=[0.0, 2.0], values=[5.0, 10.0]),
]
```

Plant `pmin` and `qmin` are minimum aggregate output while operating, not forced
commitment. `pmax` always applies. Optional `qmax` adds an aggregate upper
discharge bound. Operating bound series tighten their static domains. With no
static plant `qmax`, the sum of unit capacities is the domain ceiling.

Generators retain their `minup`, `mindown`, `initial_on` and `initial_age`.
Plants accept these fields too, with zero minimum durations by default.
Plant initial state is inferred from its units; an explicit state must agree.
Positive plant minimum durations require explicit aggregate `initial_age`,
because overlapping unit histories need not identify when the plant first
started. A supplied running-plant age must cover every currently running unit's
age; an off-plant age equals the shortest member off age. Without aggregate dwell,
the known running-age lower bound suffices. Remaining unit and plant obligations
are reported at the horizon and carried by `restart_case`.

Unit `startup` and `shutdown` series specify costs at transition times. Costs
are charged once per transition, including the first interval relative to the
initial state; they are not multiplied by interval duration.

## Directional ramps

Generators and plants accept `ramp_up`/`ramp_down` for power and
`discharge_ramp_up`/`discharge_ramp_down` for discharge. Rivers use
`ramp_up`/`ramp_down` for release. All accept corresponding operating series.
An omitted directional plant power limit inherits the existing symmetric
`ramp`; other omitted limits are disabled. Zero permits no change in that
direction. Set `ramp=nothing` (JSON `null`) to disable the symmetric default.
Rate limits are nonnegative and finite where enabled.

Interval outputs are located at their midpoints, `m[t]=(grid[t]+grid[t+1])/2`.
For discharge or river release `x`, the constraints are

```math
-\int_{m_{t-1}}^{m_t} r^{down}(s)\,ds
\le x_t-x_{t-1}\le
\int_{m_{t-1}}^{m_t} r^{up}(s)\,ds.
```

Constant rates yield the familiar limit `rate*(dt[t-1]+dt[t])/2`.
Piecewise-constant rates are integrated over physical time, including historical
knots. An unrestricted part of the comparison window leaves that directional
comparison unrestricted. Power ramps additionally allow each starting unit's
current minimum power and each stopping unit's previous minimum power. These
jumps apply to both unit and aggregate plant power ramps. Discharge and river
ramps have no transition exemption.

Generators and plants can supply `initial_power`, `initial_discharge` and
`initial_interval_hours`. The historical interval ends at the horizon start;
its midpoint defines the first comparison. Missing plant history can be derived
when every member supplies the quantity over the same historical interval.
Differing member windows require explicit aggregate history. A missing quantity
omits its first comparison; later comparisons still apply. River release history
provides the previous value and duration, when present. Explicit river
`initial_release` and `initial_interval_hours` override that fallback. Restart
carries the previous control-window average separately from transport cohorts:
a delayed confluence can contain several release pulses inside that window.

Reservoirs accept `volume_ramp_up`/`volume_ramp_down` in Mm³/hour and
`level_ramp_up`/`level_ramp_down` in metres/hour. These compare storage or level
at consecutive grid vertices, with rate integrals over the full interval.
Level limits use `head(V[t+1])-head(V[t])`, including nonlinear or tabulated
head–storage curves, rather than midpoint hydraulic head. Initial storage
already supplies the first vertex. These are limits on the discrete trajectory;
they do not establish a continuous-time derivative bound.

Finer replay preserves original operating windows for schedules and ramps.
Power, discharge and release are averaged back onto those windows, while
reservoir values are sampled at their original vertices. Restart retains
previous unit/plant output, durations, state ages and river cohorts.

## SHOP correspondence and limits

The implemented rules cover the physical purposes of SHOP commitment inputs,
maintenance, production/discharge schedules and ordinary directional ramps.
See SHOP's [generator attributes](https://docs.shop.sintef.energy/objects/generator/generator.html),
[plant attributes](https://docs.shop.sintef.energy/objects/plant/plant.html) and
[ramping tutorial](https://docs.shop.sintef.energy/examples/ramping/ramping.html).

This is an explicit OpenSHOP input profile, not an arbitrary SHOP-attribute
importer. SHOP has enable flags and defaults for several inputs; OpenSHOP enables
an optional rule by supplying its value or series. Schedules and ramps here are
hard, within physical domains; SHOP also offers soft rules and penalties.
The documented SHOP hourly scaling does not establish identical behavior on
uneven grids, so midpoint timing is an explicit OpenSHOP convention.

Global optimization retains binary unit commitment throughout the horizon.
Selective fractional commitment through SHOP `mip_flag`/`mip_length` is not
implemented. The HiGHS proposal uses the same commitment, discharge schedules
and linear ramp constraints, but omits exact power schedules and nonlinear
reservoir level ramps: its approximate power conversion cannot establish their
physical feasibility. Nonlinear dispatch and validation enforce all supplied
rules. Proposal objectives and bounds remain surrogate results.

Rolling, amplitude and nonsequential ramps, soft operating penalties and complete
SHOP attribute import remain outside this feature set.
