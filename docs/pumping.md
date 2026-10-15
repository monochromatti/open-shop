# Pumping and reversible operating modes

`Pump` constructs the common `Generator` unit record with `mode=:pump` and
`initial_on=0`. Add generating and pumping modes to `HydroSystem.generators`
in their desired row order. Commitment, `generator_q`, and `power` share this
one unit axis. Discharge and power are nonnegative magnitudes. `net_power(case,
schedule)` returns production minus consumption in MW for each interval.
`pmin`, `pmax`, power schedules, ramps and historical power are consumption
magnitudes for pumps, with the same units and timing conventions as generation.

Each plant contains one mode, and its endpoints follow actual water motion.
A pump plant's `source` is its intake and `target` is the raised-water destination;
either endpoint can connect to the existing pressure network. A generating
plant carrying the reverse flow uses the opposite endpoints. A single reversible
machine is represented by those two mode records and an explicit exclusion:

```julia
rule = CommitmentExclusion(name=:Reversible,
    groups=[[:Turbine], [:Pump]])
system = HydroSystem(reservoirs=lakes, plants=[generating, pumping],
    generators=[turbine, pump], commitment_exclusions=[rule])
```

At most one group may run per interval. Multiple units within one group may
run together. Every cross-group pair satisfies `u[i,t] + u[j,t] <= 1`; no
additional group state is introduced. Initially running modes must also obey
the exclusion. Unit minimum up/down times, startups, shutdowns, operating
schedules and plant shutdown interlocks apply to pumps. Restart retains those
histories, initial discharge/consumption and the exclusion groups. Replay with
only discharges rejects simultaneous positive flows in excluded groups; full
schedule validation additionally checks the original binary commitments.

For head loss `L(q)` on the path and outlet reference `Hout`, the effective
head is `Hin - Hout - L(q)` for generation and `Hout - Hin + L(q)` for pumping.
Pump work and consumption satisfy

```math
C\,\eta_{hyd}(q,H)\,\eta_{motor}=0.00981\,qH,
\qquad C=\frac{0.00981\,qH}{\eta_{hyd}(q,H)\eta_{motor}}.
```

Hydraulic efficiency uses the analytic parameters or the existing `TurbineTable`,
including head-dependent discharge envelopes and declared interpolation policies.
`motor_efficiency` is a constant in `(0,1]`. Zero discharge consumes zero power.
Operating head must be positive; both flow and consumption limits hold only
while committed. The optimizer keeps the work equality in its common physical
dispatch equations. Generation power supports and cuts are omitted for pump
rows because multiplying efficiency does not bound reciprocal-efficiency
consumption. Their omission can loosen a global relaxation; it does not alter
the declared pump equations or a valid certificate.

The [synthetic arbitrage example](../examples/pumped_storage.jl) moves water
between two reservoirs across a 250 m lift at low prices and returns it through
generation at high prices. Its parameters are illustrative. Tests calculate
water balances and consumption independently and also connect pumping to a
pressure junction, use hydraulic efficiency tables, and compare restarted replay.

## Fixed-speed pumps

Supply one discharge and hydraulic efficiency at each reference head:
This follows SHOP's [binary-pump definition](https://docs.shop.sintef.energy/objects/pump/pump.html#binary-pumps).

```julia
curve = PumpCharacteristic([200.0, 250.0, 300.0],
    [22.0, 20.0, 18.0], [0.86, 0.89, 0.87])
pump = Pump(name=:Pump, plant=:Lift, characteristic=curve,
    qmin=15.0, qmax=25.0, pmin=40.0, pmax=80.0,
    hmin=200.0, hmax=300.0, motor_efficiency=0.96)
```

At each operating head, this is one point, not a variable-discharge envelope:

```math
q=uQ(h),\qquad P\eta(h)\eta_{motor}=0.00981qh.
```

Both supplied functions are piecewise linear in head. The default
`head_extrapolation=:error` requires coverage of the operating head range.
`:linear` explicitly permits end-secant continuation; its discharge and
efficiency must remain positive, with efficiency at most one, throughout that
range. A single reference head requires `head_extrapolation=:constant`, declaring
constant discharge and hydraulic efficiency while consumption still uses actual
head. Raw `pump_discharge` and `pump_efficiency` lookups reject out-of-table heads
unless `extrapolation=:linear` is passed; a declared constant characteristic has
no lookup head boundary. Heads must increase strictly. No knots are added,
sorted, fitted or clipped. A characteristic cannot accompany `turbine_table`.
Static and operational flow/power limits intersect the characteristic; an
incompatible limit makes operation infeasible.

Forward integration solves characteristic discharge together with midpoint
storage, tunnel losses and pressure. In `dispatch_from_controls`, commitment is
the fixed-speed pump control; its `generator_q` entries are numerical guesses,
and zero is a valid guess for a running pump. The returned schedule contains
solved flows. `simulate(...; commitment=u)` uses the same rule; without explicit
commitment it infers on from a positive flow entry and off from zero. Ordinary
adjustable-unit flows must be zero when their supplied commitment is off;
commitment may remain on at zero flow. Schedule validation checks stored flows
without repairing them. Finer replay
holds commitment and re-solves discharge at each finer midpoint. Restart retains
the characteristic, commitment history and preceding physical flow/consumption.

The [fixed-speed example](../examples/fixed_speed_pumping.jl) shows these flows
changing with storage and checks independently calculated consumption.
Variable-discharge pumps and fixed-flow points (`qmin == qmax`) remain supported.
Electrical motor-efficiency curves are explicitly rejected for pumps.
Mixed pumping/generation members in one plant, shared ternary-machine
short-circuit equipment characteristics, reserve services and machine-specific
direction-change delays are also outside this model. Independently connected
pumps and generators may run together when no exclusion is supplied.

The [SHOP pump definition](https://docs.shop.sintef.energy/objects/pump/pump.html)
uses reciprocal hydraulic and motor efficiencies. Its
[reversible turbine example](https://docs.shop.sintef.energy/examples/pump/reversible_turbine.html)
also imposes explicit commitment exclusion between pumping and generation;
OpenSHOP uses the grouping equation above rather than reproducing SHOP's
connection-object representation.
