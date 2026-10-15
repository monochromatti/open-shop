# Gate and outlet equations

All elevations and heads use metres in the same datum. Discharge is m³/s.
Reservoir head is evaluated at midpoint storage, `Z((V[t]+V[t+1])/2)`;
it need not equal the average of the two endpoint elevations.

## Tunnel gates

`opening`, `initial_opening` and `gate_positions` are already normalized
dimensionless positions in `[0,1]`. OpenSHOP does not infer a conversion from
actuator travel, angle or gate height. Convert those quantities before supplying
positions and evaluate a loss table on the same normalized axis.

For an open tunnel from node `s` to node `d`,

```math
H_s-H_d = K(a)q|q|,
\qquad K(a)=\begin{cases}R/a & \text{without a gate-loss table},\\
R+G(a) & \text{with a gate-loss table}.\end{cases}
```

Both `resistance=R` and the ordinates of `gate_loss_factor_curve=G` are in
s²/m⁵. A supplied gate coefficient is **added** to the tunnel coefficient;
it is not a dimensionless multiplier. The table is interpolated linearly and
must cover `[0,1]`, with nonnegative values. `R` is strictly positive.

For example, `R=2`, `G(a)=6(1-a)` and a 20 m head difference give
`q=2 m³/s` at half opening and `q=√10 m³/s` at full opening. Reversing the head
difference reverses these flows. Equal heads give zero flow.

At `a=0`, the equation is `q=0`, with independent endpoint heads. We do not divide
by zero or impose equal heads across a closed gate. Closure is exact: input
tolerance and conversion cannot turn a nonzero position into a closed gate.
Without a loss table, half opening instead gives `q=√5 m³/s` in this example.
These two parameterizations are alternatives and should not be applied together.

Signed flow is permitted between hydraulic nodes, subject to `|q|≤capacity`.
A tunnel with `discharge_river` permits only forward flow. Its hydraulic target
still supplies head; discharged water enters the specified river.

The pressure law assumes submerged tunnel connections. Tunnel entrance and exit
elevations, partial entrance submergence and inlet weir geometry are not inputs
to this model. Supplying `start_height`, `end_height` or an actuator mapping as
`gate_opening_curve` in a case file raises an input error. A dry entrance must
instead be represented by an appropriate supported outlet response; it must not
be modeled as an always-submerged tunnel.

## Reservoir outlets and spill

Natural river inflow `n` enters the reach directly. Reservoir withdrawal `d`
and total reach entry `R` satisfy `R=n+d`; only `d` leaves reservoir storage.
For head above crest `ℓ=H_source-crest`, the analytic laws are

```math
d=C a\sqrt{[\ell]_+}\quad\text{(orifice)},
\qquad d=C[\ell]_+^{3/2}\quad\text{(weir)}.
```

Here `[ℓ]+ = max(ℓ,0)`. The orifice coefficient has units m^(5/2)/s;
the weir coefficient has units m^(3/2)/s. The orifice opening is a dimensionless
fraction. An uncontrolled weir always has opening one; a restrictive gate
schedule for it is rejected.

With `allow_dry=true`, reservoir withdrawal is zero below and at the crest.
Natural inflow continues through the reach in both states. Above the crest,
withdrawal follows the selected law continuously. With `allow_dry=false`,
the complete declared storage domain must remain strictly above the crest;
a case whose domain crosses it is rejected.

Local dispatch represents the square root with a nonnegative variable `w`:

```math
w^2=\max(\ell,0).
```

For positive `ℓ`, the equation forces `w²=ℓ`; for nonpositive `ℓ`, it forces `w=0`.
Thus `d=C a w` or `d=C w³`. This exact piecewise polynomial form avoids
singular square-root and fractional-power derivatives at a dry outlet, without
adding a redundant inequality. The positive-part operator retains its corner
at the crest. Global optimization uses its
existing binary positive-part graph and a squared-root equality for the same law.

A `discharge_curve` replaces the analytic law. Its abscissae are source
elevations in metres, and its ordinates are reservoir withdrawal in m³/s.
The supplied table must include its own dry behavior, cover the storage domain,
and have nonnegative, nondecreasing ordinates. Gated withdrawal is `a D(H)`;
an uncontrolled weir uses `D(H)`. The analytic coefficient and crest do not
modify a supplied complete response curve.

These are prescribed one-level responses. The analytic laws do not calculate
downstream backwater, gate geometry or submergence corrections. Where those
effects matter, supply an externally calculated response; an instantaneous
channel whose flow depends on both reservoir levels uses `law=:delta` instead.
Unsupported geometry attributes are rejected. Reservoir-dependent outlet laws
require a reservoir source and cannot replace a zero-storage river junction.

## Turbine outlet elevation

`Plant.outlet_head_floor=z_out` changes the turbine's outlet reference:

```math
h = H_source-\max(H_target,z_out)-T(q_plant).
```

Below the outlet elevation, the target level does not affect turbine head.
At and above it, the target level supplies counterpressure. The relationship
is continuous at equality. The supplied `tailwater_curve=T` adds a
nonnegative head loss in metres as a function of total plant discharge.
Neither the elevation floor nor tailwater changes the destination of the water.

The tests compare these equations with independently calculated flows, storage,
power and objectives, including closed/half/full gates, forward/reverse flow,
dry/crest/wet outlets and receiver levels on both sides of the outlet floor.
They run nonlinear and global optimization and finer replay; they do not assert
the number or names of solver variables.

The public [tunnel reference](https://docs.shop.sintef.energy/objects/tunnel/tunnel.html)
describes additional entrance and exit geometry, and the
[river reference](https://docs.shop.sintef.energy/objects/river/river.html)
describes gate-height-dependent responses. Those attributes are not implied by
OpenSHOP's normalized positions or analytic outlet laws. The equations and
exclusions above define this package's supported hydraulic profile.
