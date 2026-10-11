# Supplied channel and gate-loss curves

OpenSHOP uses supplied hydraulic response functions. Cross-section geometry
and calibration remain outside the scheduling model.

## Tunnel gates

`Tunnel.gate_loss_factor_curve` is a `TableCurve` from the existing opening
coordinate in `[0,1]` to additional resistance in s²/m⁵. For a nonzero opening
`s`, the signed tunnel equation is

```math
H_u-H_d = \bigl(R_{tunnel}+R_{gate}(s)\bigr)q|q|.
```

The table must cover `[0,1]` and contain nonnegative coefficients. Opening zero
means closed and imposes `q=0`; it does not impose equality of endpoint heads.
An `OperationalSeries` on `opening` prescribes its chronological position.
Openings remain prescribed when `gate_positions` is empty. Nonempty listed
positions enable [discrete gate decisions and interlocks](tunnel-gates.md);
the positive-position resistance and closed-gate laws stay the same.

Without a gate-loss curve, the existing law remains
`s*(H_u-H_d)=R_tunnel*q*abs(q)`. No gate-position-to-height curve is needed.
SHOP's tunnel documentation describes the same additive resistance law, with
its gate-position mapping supplied separately. This profile uses its existing
normalized opening coordinate directly.

## Reservoir-to-reservoir channels

A `River` with `law=:delta` connects two distinct reservoirs. Its signed
upstream release is

```math
q=F(H_u-H_d,H_{ref}).
```

Set exactly one of `delta_level_ref_up_flow_curve` and
`delta_level_ref_down_flow_curve`; this selects `H_ref=H_u` or `H_ref=H_d`.
A `LevelFlowCurve` contains increasing `reference_levels` and one `TableCurve`
per reference. Each table maps level difference in metres to signed discharge
in m³/s. Delta grids can differ. Every curve must contain `(0,0)`, have
nondecreasing discharge, and preserve the sign of the level difference.
Forward and reverse responses can be asymmetric.

Discharge is linear between delta knots and between neighboring reference
levels. The physical model linearly extrapolates the outer delta segments.
SHOP also documents extrapolation outside delta ranges without specifying
the extrapolation formula. Outside the reference-level span,
OpenSHOP holds the nearest reference curve; SHOP's public documentation does
not specify this case. Cover the intended absolute-level range to avoid
relying on that extension. A single reference curve is independent of absolute
level. Public `level_flow_value` calls reject delta extrapolation unless
`extrapolation=:linear` is supplied explicitly.

```julia
response = LevelFlowCurve([100.0, 110.0], [
    TableCurve([-5.0, 0.0, 5.0], [-15.0, 0.0, 20.0]),
    TableCurve([-5.0, 0.0, 5.0], [-25.0, 0.0, 30.0]),
])
channel = River(name=:Connection, source=:Upper, target=:Lower,
    law=:delta, delta_level_ref_up_flow_curve=response,
    curves=RiverRouting.DelayCurve[], deterministic_delay=0.0,
    capacity=50.0, water_value=0.0)
```

Levels use the same midpoint storage states as other hydraulic equations.
A positive `q` withdraws from `source` and enters `target`; negative flow
reverses both transfers. Capacity limits `abs(q)` and can be tightened with an
operating series. Signed-flow ramp limits and `initial_release` retain their
usual chronological meaning.

These are instantaneous, zero-storage hydraulic channels. They require zero
delay and no delay curves, natural inflow, historical inventory, outlet gate
or environmental arrival law. Channel loops are permitted. Ordinary delayed
river networks still require an acyclic topology and nonnegative transport;
a signed channel cannot feed a confluence directly. Channel history is zero
because no water remains in transit. Restart carries reservoir states and the
previous signed flow for ramp constraints.

Local dispatch evaluates the response with explicit analytic gradients and
Hessians. The global model uses the union of delta knots and shared SOS2
coordinates to represent the same piecewise bilinear graph exactly. It adds
one discharge column per reference coordinate instead of a product variable
for every pair of weights. Commitment proposals allow signed channel flow;
the hydraulic response is enforced by the subsequent nonlinear dispatch.
Simulation and validation evaluate the supplied curves independently.

Existing `discharge_curve` inputs continue to provide upstream-level response
functions, corresponding to the single-curve `up_level_flow_curve` use case.
Existing distributed delay curves provide the `time_delay_curve` response.
This change adds no geometry calculation or gate-height model.

See the [SHOP tunnel attributes](https://docs.shop.sintef.energy/objects/tunnel/tunnel.html)
and [SHOP river flow functions](https://docs.shop.sintef.energy/objects/river/river.html).
