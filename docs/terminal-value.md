# Terminal water values

`Reservoir.water_value` and `River.water_value` are supplied constants in
currency/Mm³. They value future opportunities beyond the modeled horizon.
They are not electricity prices or marginal values calculated by the optimizer.

## Energy-equivalent conversion

At a supplied reference net head `h` metres and combined turbine/electrical
efficiency `η`, one Mm³ produces

```math
f = \frac{0.00981}{0.0036}h\eta = 2.725h\eta
\quad\text{MWh/Mm³}.
```

This follows from `P=0.00981 q h η` MW and `volume=0.0036 q Δt` Mm³.
For energy value `e` in currency/MWh, the local volume value is `e f`.

```julia
factor = energy_conversion_factor(120.0, 0.90) # 294.3 MWh/Mm³
lower_value = volume_water_value(30.0, energy_conversion_factor(80.0, 0.92))
upper_value = volume_water_value(30.0, factor; downstream_value=lower_value)
```

The lower value is `6016.8` and the upper value is `14845.8 currency/Mm³`.
The upper calculation assumes each unit of water subsequently generates at
both stations. Add only downstream opportunities not already included in
`factor`. Alternative paths and partial delivery require supplied weighted
values; neither helper traverses the network or assumes every plant will run.

The reference head and efficiency stay constant in this valuation. Dispatch
still uses actual nonlinear head and efficiency. Future losses, pumping costs,
limited generation and uncertain prices must be reflected in supplied values
where relevant; the conversion does not predict them. Efficiencies must lie
in `(0,1]`, heads must be nonnegative, and inputs/results must be finite.
Negative energy values are permitted.

## Reservoir and transit accounting

The objective values increments from the initial state:

```math
J = \text{energy revenue} - \text{costs}
  + \sum_i\lambda_i(V_{i,T}-V_{i,0})
  + \sum_r\lambda_r(W_{r,T}-W_{r,0}).
```

`W[r,0]` is historical water already in transit, not zero by convention.
`W[r,T]` includes all cohorts still in that reach at the end. A unit of water
is valued where it resides, so it is not counted in both a reservoir and transit.

Water between the stations can generate only below. The example sets that
reach's value to `lower_value`, not `upper_value`. Every reach needs an explicit
volume value; no downstream value is inferred from connectivity. Historical
transit and natural inflow may also have future value.

Absolute terminal valuation differs by the constant
`Σ λ_i V_i,0 + Σ λ_r W_r,0`. It gives the same maximizer for fixed coefficients
and a different reported objective. OpenSHOP reports the incremental form.
For the same executed controls and trajectory, with unchanged values and prices,
prefix and restarted objectives add to the
full-horizon objective. Changing valuation between horizons changes the problem.

`objective_components(case, schedule)` independently recomputes energy revenue,
reservoir/transit increments, transitions, gate costs, release penalties and
their total. It does not establish feasibility or return a dual price; use
`validate` and finer replay separately.

## A hand calculation

The [cascade example](../examples/terminal_value.jl) has a half-hour reach and
historical release `2 m³/s`, giving initial transit `0.0036 Mm³`.
Upper/lower discharges of `4/2 m³/s` in both hourly intervals give
upper storage `[2, 1.9856, 1.9712]` and lower storage `[2, 2.0036, 2.0108]` Mm³.
Terminal transit is `0.0072 Mm³`.

Reservoir-value increments are `14845.8 × (-0.0288)` and `6016.8 × 0.0108`.
The transit increment is `6016.8 × (0.0072 - 0.0036)` currency units.
Actual upper net heads are `119.991` and `119.9712 m`, and lower heads are
`80.0018` and `80.0072 m`. Revenue uses these heads and current electricity
price `100 currency/MWh`, not future value `30`. Startup charges total `8`.
Tests compare these components with reconstruction, optimization and restart.

## Currency and marginal values

OpenSHOP has no implicit NOK/EUR convention and retrieves no exchange rates.
Normalize every price, startup/shutdown charge, gate charge, penalty, reserve
payment and water value to one currency before constructing a case. A fixed
exchange rate must scale all monetary inputs, including downstream values;
it does not scale physical quantities or MWh/Mm³ factors. Tests confirm that
this scales the complete objective while preserving physical schedules.

Supplied water value is an objective coefficient. A constraint dual, when
available from a continuous solve, is a sensitivity of that solve under its
current assumptions and active constraints. Discrete commitment and
nonconvexity limit that interpretation. The helpers do not calculate such
sensitivities or claim a globally valid marginal water price.

The public [reservoir reference](https://docs.shop.sintef.energy/objects/reservoir/reservoir.html)
describes local energy-value conversion using the first downstream plant and
initial head. OpenSHOP requires reference head, combined efficiency and
downstream value explicitly. The public description does not establish every
efficiency or network aggregation choice; the equations above state this
package's independent assumptions.
