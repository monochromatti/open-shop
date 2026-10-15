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

Unit limits accept `0 <= qmin <= qmax` and `0 <= pmin <= pmax`, with finite
bounds. A zero minimum introduces no artificial epsilon. Commitment remains a
separate operating state: a running unit may have zero flow and power when its
limits and conversion law permit that point, and still incurs transitions and
dwell obligations. An explicit zero output schedule retains its off-state
meaning. Maintenance and forced states also retain their ordinary meanings.

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

## Gate handovers and plant shutdowns

`TunnelInterlock` retains its instantaneous `min_open`/`max_open` limits.
It can also require a plant shutdown after a gate closes:

```julia
TunnelInterlock(name=:SupplyHandover, tunnels=[:HighIntake, :LowIntake],
    min_open=1, max_open=1,
    forward_tunnels=[:HighIntake], backward_tunnels=[:LowIntake],
    plants=[:Station], forward_switch_time=1.5, backward_switch_time=0.5)
```

The forward and backward lists identify two disjoint closure groups. They are
subsets of `tunnels`; their names describe connection roles, not the sign of
water flow. `plants` explicitly names the stations protected by the wait.
This avoids inferring a station from a shared or looping pressure network.
The times are finite, nonnegative hours and constant over the horizon. Zero
disables a group's wait. A positive time requires that group and protected plants.

A closure occurs at a control edge when a gate changes from a strictly positive
position to zero. It starts the corresponding wait even if the tunnel carried
no water. Positive-to-positive movement does not start a wait. Every unit of
every protected plant must stay off until expiry. Another intake may open
immediately while the plant is idle: the example keeps one gate open throughout.
Reopening a gate does not cancel its unfinished wait, and later closures start
new waits. The latest applicable deadline controls operation.

For open status `o[e,t]` and unit commitment `u[g,k]`, a closure at `τ[t]` is
enforced by the linear inequalities

```math
u_{g,k}+o_{e,t-1}-o_{e,t}\le1,
\qquad \tau_t\le\tau_k<\tau_t+L.
```

Here `L` is the group's waiting time, and `o[e,-1]` is the initial gate status.
When a gate closes, the difference is one and each protected `u[g,k]` must be
zero. In the other transitions the inequality adds no restriction. The existing
gate and unit variables suffice; no additional closure binary is needed.

These comparisons use elapsed hours, including on irregular grids, with a
`1e-8` hour boundary tolerance matching minimum dwell checks. Controls are constant
within each interval. If a wait expires inside an interval, production can resume
only at the next control edge. Include the expiry as an edge to permit an exact
restart time. Hard output schedules, minimum dwell and shutdown waits all apply;
conflicting requirements make the optimization infeasible.

`initial_forward_wait` and `initial_backward_wait` describe unexpired waits at
the horizon start; their default zero declares none. Each is bounded by its
group's switch time. For positive remaining wait `w`, the implied closure was
`L-w` hours before the horizon. All protected units must initially be off with
an off-age at least that long. Initial gate status alone cannot describe this
history, since a gate may have reopened while the plant stayed off.

`validate_tunnel_openings(case, positions, commitment)` checks the joint rule.
The two-argument form checks positions, hard schedules and open counts.
`simulate` also rejects positive protected-unit discharge during a wait;
schedule validation additionally rejects a claimed on-state with zero discharge.
Finer replay preserves gate event times. Validation reports outstanding waits
in `terminal_tunnel_switch_obligations`, and `restart_case` carries both residual
clocks along with initial openings and unit histories. Residual waits are zero
once they have expired within the same time tolerance used for admission.
Conversion to Float64 must preserve a gate's zero/nonzero state and cannot admit
controls outside their original `[0,1]` domain.
The [directional-interlock example](../examples/directional_interlocks.jl) optimizes
a two-intake station with asymmetric waits on an irregular grid.

## Explicit penalty limits

`ScheduleCase.soft_constraints` enables bounded deviations from selected
reservoir observations or unit/plant output targets. Existing hard schedules,
availability, dwell, ramps and safety limits keep their meanings. In particular,
a positive soft output target does not force a unit on, and a zero soft target
does not force it off.

For a named observation `x`, the lower and upper limits use separate slacks:

```math
x+s^-\ge L,\qquad x-s^+\le U,\qquad
0\le s^-\le\bar s^-,\qquad0\le s^+\le\bar s^+.
```

Omit an inequality when its limit is absent. A zero allowance keeps that side
hard. The cost subtracted from the objective is

```math
C=a\,(c^-s^-+c^+s^+).
```

For reservoir volume, `x=V(time)` and `a=1`; for level,
`x=head(V(time))` and `a=1`. The cost is charged once at that state observation.
A terminal target uses its original absolute terminal time. Its cost is not a
running cost per hour.

For power or discharge, `x` is the duration-weighted average over the supplied
window `[time,end_time]`. A plant observation sums its units' outputs. The
power quantity is a nonnegative magnitude: output for generators and
consumption for pumps, as in hard operating schedules. It is distinct from
signed net exchange with the grid. The
factor `a` is the window duration in hours for power, or `0.0036` times that
duration for discharge. These average targets permit one part of a window to
offset another; use separate records to restrict separate windows.

```julia
target = SoftConstraint(name=:TerminalStorage, object=:Lake,
    quantity=:volume, time=24.0, lower=5.0,
    max_lower_slack=0.1, lower_cost=1000.0)
production = SoftConstraint(name=:DesiredProduction, object=:Station,
    quantity=:power, time=8.0, end_time=9.0, lower=20.0, upper=20.0,
    max_lower_slack=5.0, max_upper_slack=2.0,
    lower_cost=200.0, upper_cost=100.0)
# Pass these separately as ScheduleCase(...; soft_constraints=[target, production]).
```

Conservation, head–storage relationships, pressure loss, turbine conversion and
all static physical domains stay hard. A reservoir's static `vmin`/`vmax`
remain safety bounds. The numerical domains used for global optimization can
intersect them with `lower-max_lower_slack` and `upper+max_upper_slack` for volume
observations. They cannot tighten to the nominal target while slack is allowed,
or extend beyond the physical domain. Level observations retain the supplied
head–storage law in nonlinear dispatch and the full global model.

`soft_constraint_assessment(case, schedule)` independently computes the minimal
actual deficit and excess, their units and costs, and any excess over the
allowances. `Schedule.soft_slacks` and `soft_penalty_cost` expose this accounting.
`validate` checks the reported accounting and objective, lists permitted
deviations in `permitted_penalties`, and rejects a deviation beyond its allowance.
An allowed penalty is compatible with a valid physical schedule.

Refinement preserves the original observation times and rate windows. It does
not add penalty points or multiply costs by the number of finer intervals.
Restart drops observations and windows already completed, retains future ones
at their original times, and therefore does not move a terminal target onto a
new horizon endpoint. Restart inside a rate window is rejected: this profile
does not carry partially assessed production history. Choose restart edges at
assessment boundaries, or declare smaller original windows.

The [penalty-cost documentation](https://docs.shop.sintef.energy/commands/penalty_cost.html)
and [individual terminal-penalty documentation](https://docs.shop.sintef.energy/commands/individual_endpenalty.html)
describe production, discharge and terminal storage penalty costs. OpenSHOP
keeps their quantity-specific units explicit and always keeps water balance hard.
The [worked example](../examples/soft_limits.jl) chooses full discharge at
1000 currency/Mm³ and saves the water at 20000 currency/Mm³.

## Optimization scope

Unit commitment remains binary throughout the horizon. The MILP proposal enforces
commitment, discharge schedules and linear ramps to find candidate states.
Nonlinear dispatch enforces physical power schedules and reservoir level ramps.
The proposal also uses soft-limit costs, with surrogate power and a secant
head–storage approximation for level observations; its candidate still needs
the exact nonlinear dispatch and independent assessment.
Proposal objectives are surrogate values, not bounds on the physical problem.

Rolling, amplitude and nonsequential ramps and conditional environmental waivers
are outside this input profile.
