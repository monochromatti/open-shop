# Equipment curves and river delays

Calibration data and operating limits describe different things. A turbine's
efficiency knots describe its energy conversion. Its discharge envelope states
where it may run. OpenSHOP keeps both inputs and requires an explicit policy
when the operating domain extends beyond the calibration data.

## Turbine tables

`TurbineTable` accepts discharge rows and reference-head columns. Heads are metres,
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

### Unequal discharge grids

Supply one `TableCurve` per reference head when the measured discharge knots
differ. Both the rectangular constructor and this headwise form use the same
numerical representation:

```julia
table = TurbineTable(heads=[100.0,200.0], curves=[
    TableCurve([1.0,3.0,5.0], [0.80,0.90,0.88]),
    TableCurve([1.0,2.0,5.0], [0.82,0.87,0.94])],
    qmin=[1.0,1.0], qmax=[5.0,5.0], interpolation=:pchip_discharge)
```

Each original curve is interpolated on its own discharge grid. At a head
between references, OpenSHOP linearly interpolates those efficiencies at the
same discharge. In the example, linear discharge interpolation gives
`η(2.5,100)=0.875` and `η(2.5,200)=0.881666…`, hence
`η(2.5,150)=0.878333…`. Original-grid PCHIP gives `η(2.5,100)=0.891875`.
Resampling that first curve at an added `q=2` knot and recalculating PCHIP
slopes would give `0.89025`; OpenSHOP preserves the original `0.891875`.

The global equations use the union of the supplied discharge knots only as a
polynomial partition. Each original cubic is re-expressed on its subintervals;
its slopes are never recalculated on the union grid. These partition nodes are
not additional calibration measurements. `table.curves` exposes the owned
original per-head `TableCurve` records; derived knots and slopes are internal.
Treat calibration tables as fixed after construction. Rectangular constructors
and JSON remain supported. Unequal data serialize as `curves` containing the
original `x`/`y` arrays in head order, without resampled columns.

Efficiency coverage is separate from the supplied `qmin`/`qmax` envelope.
A strict raw lookup at an exact reference head needs only that curve's domain.
Between heads, discharge must be covered by both participating curves.
Explicit linear continuation uses each original curve's own end secants before
interpolating in head. A complete case with `discharge_extrapolation=:error`
requires every head column participating over the unit's operating head range
to cover its declared static discharge range, including the outer pair when
head extrapolation is enabled. Otherwise declare discharge continuation
explicitly; no operating limits are clipped or inferred from efficiency knots.

The [worked example](../examples/unequal_turbine_grids.jl) checks water, power
and revenue against hand calculations and solves the case through the ordinary
API. The [SHOP generator reference](https://docs.shop.sintef.energy/objects/generator/generator.html#turb-eff-curves)
describes discharge/efficiency curves at reference heads but does not establish
this unequal-grid interpolation rule. The policy above is OpenSHOP's explicit
mathematical contract; it does not claim undocumented SHOP equivalence.

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

The reference is the flow of **that release cohort**. It is not the previous
time interval's release, the arrival flow, or a previous optimization iterate.
Writing the reference flows as `q[j]`, a cohort released at `Q` between two
references uses

```math
w=\frac{Q-q_j}{q_{j+1}-q_j},\qquad F_Q=(1-w)F_j+wF_{j+1}.
```

Here `F` is the travel-time CDF. Each historical cohort uses its own original
release and interval to select these weights. A subsequent change in flow does
not reselect the distribution of water already released. At a distributed
confluence, total entry is the allocated arriving flow plus branch-local
injections, averaged on its declared transport interval. This is the flow used
for the new downstream cohort.

Raw routing accepts the first and last reference exactly and rejects a flow
below or above a multi-curve reference range. The full case additionally needs
zero-through-capacity coverage. A single curve applies at every nonnegative
flow, subject to the reach capacity in a full case. Endpoint holding outside
a multi-curve domain is not an admitted routing policy.

Numerical trial continuations are separate: Newton holds the nearest endpoint
mixture when evaluating a trial release outside the references. Nonlinear
optimizer transfer operators instead continue the outer released-transfer
polynomial. Physical checks and optimization bounds retain the declared domain;
neither continuation admits an out-of-range release.

The nonlinear dispatch and global model both use `Q F_Q`, so the released
transfer is piecewise quadratic in `Q`. The MILP commitment proposal freezes
transfer coefficients at its reference schedule or a default capacity-based
flow as a surrogate; its objective is not a bound on this nonlinear model.
Feasibility and replay use the contemporaneous-cohort law.

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

In particular, the river reference describes SHOP optimization selecting and
interpolating the wave curves using flow from its **previous iteration**. This
is an iterative approximation, distinct from OpenSHOP's contemporaneous-flow
equations above. The public [time-delay example](https://docs.shop.sintef.energy/examples/river/advanced_time_delay.html)
illustrates historical flow and dispersed arrivals, with one flow-independent
wave curve. Matching such trajectories cannot establish equal optimization
problems, multi-reference endpoint behavior or global bounds. Synthetic
endpoint-holding simulations likewise do not establish SHOP's optimization
equations. OpenSHOP makes no optimizer-equivalence claim from these examples.
