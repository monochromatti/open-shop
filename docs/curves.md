# Equipment curves and river delays

Calibration data and operating limits describe different things. A turbine's
efficiency knots describe its energy conversion. Its discharge envelope states
where it may run. OpenSHOP keeps both inputs and requires an explicit policy
when the operating domain extends beyond the calibration data.

## Turbine tables

`TurbineTable` has discharge rows and reference-head columns. Heads are metres,
discharges are m³/s, and efficiencies are fractions in `(0, 1]`. Discharge knots
must contain at least two strictly increasing points. Each reference head has
supplied `qmin` and `qmax`, independently of the efficiency knots. These limits
intersect the generator's own operating limits.

| Input | Policy | Meaning |
| --- | --- | --- |
| Multiple reference heads | `head_extrapolation=:error` (default) | Head knots must cover the operating head range. |
| Multiple reference heads | `head_extrapolation=:linear` | Continue the outer head secants. |
| One reference head | `head_extrapolation=:constant` (required) | Efficiency and the discharge envelope are independent of head. Power still uses actual head. |
| Discharge outside efficiency knots | `discharge_extrapolation=:error` (default) | Reject incomplete coverage of the table envelope or generator discharge bounds. |
| Discharge outside efficiency knots | `discharge_extrapolation=:linear` | Continue the end discharge secants without changing the operating envelope. |

No calibration points are invented. A single-head table retains its one head
and efficiency column when serialized. `:constant` is rejected for multi-head
data, where its meaning would be ambiguous. Negative or reversed discharge
limits, nonfinite points and unsupported policy names are rejected.

Inside the table, `interpolation=:bilinear` is linear in discharge and head.
`:pchip_discharge` uses shape-preserving cubic Hermite interpolation in
discharge and linear interpolation in head. Outside the discharge knots,
`:linear` always means end **secants**, including for PCHIP. It does not continue
the cubic endpoint tangent.

Raw lookup functions are strict by default. For example,
`OpenSHOP.turbine_efficiency(table, q, h; extrapolation=:linear)` requests
numerical continuation explicitly. Declaring operating extrapolation permits
that domain in a complete case but does not change the default raw lookup.
A declared single-head constant response accepts any finite lookup head.
Generator head and efficiency limits still apply while the unit runs.

Solver trial points and idle units use finite numerical extensions. These
extensions do not admit an otherwise invalid operating schedule. Optimization,
independent reconstruction and validation use the same interpolants and
operating policies. Full-model simulation rejects incomplete declared
calibration, but a converged trajectory may still violate operating limits.
Use `validate` or `dispatch_from_controls` to check those limits.

The electrical-efficiency curve is indexed by **electrical output** in MW.
Its efficiencies must be positive fractions and it must cover zero through
maximum power. Positive segment intercepts make `P/ηₑ(P)` strictly increasing,
so the implicit power equation has at most one solution. Zero discharge gives
zero power. Zero electrical efficiency, including at the zero-power knot, is
currently unsupported. Without an electrical-efficiency curve the conversion
factor is one.

### A calculation beyond the efficiency knots

The [single-head example](../examples/turbine_tables.jl) declares a 100 m
reference head, efficiency knots `(1, 0.88)` and `(3, 0.92)`, and an operating
range of `0.5–4 m³/s`. Linear continuation gives `ηₜ(4)=0.94`. Electrical
efficiency is `0.95`.

Reservoir elevation is `98 + V` metres, initial storage is `2 Mm³`, and inflow
is `2 m³/s`. Releasing `4 m³/s` in each of two one-hour intervals removes
`0.0072 Mm³` per interval. The midpoint heads are `99.9964` and `99.9892 m`.
Electrical output is `P = 0.00981 × 4 × h × 0.94 × 0.95 MW`.
At prices `50` and `100` currency/MWh, revenue is `50 P₁ + 100 P₂`.
Terminal water value and transition costs are zero in this example.
Tests compare storage and power within `1e-10` and the objective within `1e-8`
currency units. Solver comparisons use `1e-3` currency units.

## Distributed river delays

OpenSHOP uses the **strict coverage** policy. Multiple reference curves must
include a zero-flow reference and cover the reach's static capacity. Missing
zero-flow calibration is rejected even if a proposed schedule keeps every
release above the smallest reference. This applies to optimization, full-model
simulation, validation and serialized case loading. Operating capacity may
tighten the static capacity.

A single `RiverRouting.DelayCurve` declares a flow-independent distribution.
Its reference flow is descriptive: zero, smaller positive and larger flows
use that same distribution within the reach's capacity. No low-flow curve is
inferred. Raw `RiverRouting` functions can use a smaller positive reference
domain, but reject release or history flows outside it when multiple curves
are supplied.

Delay edges are **hours**. Each weight is the probability that released water
has a travel time in that bin, uniformly within it. It is not a probability
density. For edges `[0, 1, 5]` and weights `[0.25, 0.75]`, the densities are
`0.25` and `0.1875` per hour, and mean travel time is
`0.25 × 0.5 + 0.75 × 3 = 2.375 hours`.
Weights must be nonnegative and sum to one within `1e-12` absolute and relative
tolerances. OpenSHOP does not rescale them. A trailing zero weight may be
supplied when there is one weight per edge.

At an intermediate release flow, the neighboring distributions are mixed
linearly. Their bin edges may differ. The mixture is selected from total reach
entry flow, including natural inflow and upstream arrivals, and stays attached
to the released cohort. Zero flow releases no water. Historical cohorts keep
their original flows and release intervals, including water still in transit
at time zero.

For example, mix a slow uniform delay `[2, 3]` at reference flow zero and a
fast uniform delay `[0.5, 1.5]` at reference flow ten. A one-hour release at
`2.5 m³/s` uses weights `3/4` and `1/4`. By hour two, the slow fraction has
delivered nothing and the fast fraction has delivered `7/8`. Arrival volume is
`0.0036 × 2.5 × 1/4 × 7/8 = 0.00196875 Mm³`.
Of the released `0.009 Mm³`, `0.00703125 Mm³` remains in transit.
Tests check these values within `1e-12 Mm³`.

Distributed confluences conserve water but average arriving streams on the
transport grid before sending them downstream. `transport_audit` reports
cumulative-arrival and terminal-inventory differences on a finer grid
separately from equation feasibility. Conservation alone does not establish
arrival-time accuracy. Restart retains release histories, not just averaged
arrival buckets.

## Public reference and scope

The public [generator reference](https://docs.shop.sintef.energy/objects/generator/generator.html)
describes efficiency curves and head-dependent operating ranges. The public
[river reference](https://docs.shop.sintef.energy/objects/river/river.html)
describes flow-dependent wave shapes and bin volume fractions. Neither
reference establishes every endpoint policy above. The single-head assumption,
separate operating envelope, strict low-flow admission and numerical
extensions are OpenSHOP's declared semantics. They do not establish identical
behavior in another optimizer.
