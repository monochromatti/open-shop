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

Variable discharge envelopes and fixed-flow points (`qmin == qmax`) are
supported. A fixed-flow point is not a fixed rotational-speed characteristic.
SHOP's head-dependent single-point fixed-speed pump curves require a dedicated
characteristic that is outside this scope; no extra calibration points are
invented. Electrical motor-efficiency curves are explicitly rejected for pumps.
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
