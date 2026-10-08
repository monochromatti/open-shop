# Model

The control grid has intervals of duration Δt in hours. Storage is in Mm³,
discharge in m³/s, head in metres and electrical power in MW. Hydraulics are
quasi-steady within an interval; reservoir storage evolves with an implicit
midpoint step.

## Water and hydraulic head

For reservoir i and interval t,

```math
V_{i,t+1}=V_{i,t}+0.0036\,\Delta t_t
  (I_{i,t}+Q^{in}_{i,t}-Q^{out}_{i,t}).
```

Its hydraulic elevation is the supplied head–storage function evaluated at
`(V[i,t]+V[i,t+1])/2`. A hydraulic junction has no storage and balances signed
flows. A boundary has a prescribed head and permits external water exchange.

A tunnel connecting a to b obeys

```math
o_t(H_{a,t}-H_{b,t})=kQ_t|Q_t|.
```

Opening zero fixes its flow to zero. Capacity bounds both directions. The global
model expresses signed flow using positive and negative parts and a direction
binary, retaining the exact quadratic loss law. SCIP can infer a fixed direction
during presolve. Reversible links retain both directions.

## Generation and operation

For each unit, u is binary. Discharge and power lie between their supplied
minimum and maximum values when u=1; both are zero when u=0. Minimum up/down
times use interval edges and the supplied initial state and elapsed age.
Startup and shutdown variables record transitions.

Plant net head is upstream hydraulic head minus downstream outlet head and
the supplied flow-dependent tailwater term. An outlet-head floor gives
`max(receiver_head, outlet_head_floor)` while water still enters the original
receiving node.

```math
P=0.00981\,q\,h\,\eta_{turbine}(q,h)\,\eta_{electrical}(P).
```

Efficiency can be an analytic function or a turbine table. Turbine tables are
linear in head and linear or shape-preserving cubic in discharge. The global
model uses exact polynomial expressions within selected table cells, including
their off-state continuation; it does not substitute a linear power curve.
Head-dependent turbine flow envelopes apply to running units. A fixed off
unit has exactly zero discharge and power; its hydraulic head remains part of
the network. A fixed on unit's table domain starts at its declared minimum
flow. Efficiency bounds include endpoints, cubic stationary points and the
original table's secant extrapolation.

Aggregate plant capacity and production ramp limits apply to interval-average
power. Start/stop allowances use each unit's minimum power. If previous plant
power is supplied, the first interval is also ramp constrained. Any remaining
terminal dwell obligation is returned for the next horizon.

## Rivers

Controlled outlets, orifices, weirs and supplied discharge curves determine
source releases. Dry-capable outlets use an explicit wet/dry branch. Rivers
conserve released cohorts through their travel-time models and confluences.

For deterministic delays, transport compilation preserves original release
pulses through successive confluences. Zero and sub-interval delays are allowed;
same-interval arrivals remain coupled to hydraulic and storage equations.
Environmental minima can apply to interval averages or exact instantaneous
breakpoints. A named aggregate flow observation can additionally constrain the sum of unit discharges, river releases and exogenous inflow without adding a hydraulic node or counting the water twice.

Distributed travel-time curves blend distributions according to the release
flow. Each original cohort retains its chosen distribution. Downstream mixing
uses the declared scheduling grid. This approximation is checked with transport
and grid refinement; its optimization bound does not certify a continuous-time
mixing model.

Historical releases establish initial water in transit. Terminal water still in
transit retains the explicitly supplied river water value.

## Objective

The maximized objective is

```math
\sum_t \pi_t\Delta t_t\sum_g P_{g,t}
-\sum_{g,t}(c^{start}_g s^{up}_{g,t}+c^{stop}_g s^{down}_{g,t})
-C_{release}
+\sum_i w_i(V_{i,T}-V_{i,0})
+\sum_r w_r(W_{r,T}-W_{r,0}).
```

Release penalties use declared shortfall volumes. Constant initial inventories
are subtracted so they do not inflate the objective or relative-gap denominator.

## Feasibility and certificates

SCIP searches the bounded algebraic model and returns an upper bound U for this
maximization problem. A schedule is reconstructed from its controls by a separate
Newton-based chronological simulator, then all original equations and operating
constraints are audited. Its independently evaluated objective is the lower
bound L.

Tiny solver boundary errors are recorded and normalized only for discharge
between −1e−6 and zero, or a gate at most 1e−8 outside [0,1]. The resulting
controls receive a fresh physical reconstruction, objective and audit. No
storage, head or hydraulic state is clipped. Larger control violations fail.

The reported relative gap is `max(0,U-L)/max(1,abs(L))`, provided U encloses L
within 1e−6 objective units. A global certificate requires that interval to meet
the requested relative or absolute gap. It is a numerical solver certificate
for the discrete model, rather than an interval-arithmetic proof. When
`fixed_u` is supplied, that model restricts commitment to the supplied matrix;
its certificate is conditional on those on/off decisions. `commitment_fixed`
and `certificate_scope` identify this restriction in the returned result.

A native `OPTIMAL` result is withheld as a certificate when SCIP performed LP iterations but returned no finite first root LP bound. An observed numerical failure produced this combination and an upper bound below a separately audited feasible schedule. The native result remains in diagnostics. This conservative guard is not a general proof against floating-point solver errors.

Finer chronological replay is a separate acceptance test. An optimal discrete
schedule can fail replay when its grid is too coarse. The library reports this
failure rather than attaching the discrete certificate to different controls
or a refined model.

Supplying an initial schedule does not fix commitment. Every primary and
auxiliary start is constructed and checked against all algebraic constraints,
bounds and binary domains before SCIP receives it.

When a better discrete candidate fails replay, a supplied replay-valid initial
schedule can be retained. Its objective and gap are recomputed against the same
global bound. The rejected candidate remains in `discrete_candidate` for
diagnostics; its certificate is never transferred to the retained controls.

## Algebraic representation and diagnostics

The global model uses conservative capacity-based storage and arrival domains.
Forward and backward passes enclose every feasible storage trajectory; no local
schedule is used to restrict the feasible set. Polynomial ranges include
interior stationary points, with outward floating-point slack.

Every linear curve uses an SOS2 coordinate graph. At most two adjacent weights
are nonzero; coordinate and curve value are their weighted sums of the original
knots. Turbine tables use independent discharge weights `λ` and head weights
`μ`. Compatible head coordinates are shared by units and flow envelopes.

PCHIP discharge interpolation is represented exactly by quadratic products:

\[
w_i=\lambda_i\lambda_{i+1},\qquad
r_i=\lambda_i w_i,\qquad
 e_j=\sum_i E_{ij}\lambda_i+\sum_i\left(A_{ij}r_i+B_{ij}(w_i-r_i)\right),
\qquad \eta=\sum_j\mu_j e_j.
\]

On an active adjacent pair, `λᵢ+λᵢ₊₁=1`, so `rᵢ=λᵢ²λᵢ₊₁` and
`wᵢ−rᵢ=λᵢλᵢ₊₁²`. The coefficients reproduce the original PCHIP cubic;
clipping a domain does not recompute slopes. Bilinear tables have `A=B=0` and
need no discharge products. Products are shared across head columns. The
turbine power relation remains nonlinear. SOS2 constraints still require
integer search; they do not turn this into a continuous convex problem.

The power equation is strengthened with supporting inequalities for the original
turbine polynomial. For each on-domain, a guarded Bernstein enclosure supplies
`c >= max(P_upper(Q,H) - a Q - b H)`. The head-aware row is

\[
P\le aQ+b(H-H_*)+(c+bH_*)u,\qquad H_*=H_{\rm lo},\quad b\ge0.
\]

On, it is the supporting plane. Off, `P=Q=0` and its right-hand side is
nonnegative over the full shared-head domain, including negative heads. The
original interpolation, efficiency and nonlinear power equalities remain in
the model. These rows strengthen SCIP's relaxation; they do not replace the
physical functions. Their floating-point guards retain the model's numerical
certificate scope, rather than supplying a formal rounding-error proof.

A joint power envelope also represents discharge, on-state head and turbine
efficiency with at most eight corner weights `z`. Their mass is `u`; their
moments reproduce discharge and auxiliary on-state head and efficiency.
Linear binary-product hulls enforce `h_on = u h` and `eta_on = u eta`, using
full source bounds so an off unit retains its hydraulic head and efficiency
continuation. With corner coordinates `(q_c,h_c,eta_c)`, the added bounds are

```math
0.00981 e_{min}\sum_c q_c h_c\eta_c z_c
\le P\le
0.00981 e_{max}\sum_c q_c h_c\eta_c z_c.
```

Electrical-efficiency extrema cover the operating power interval. The on-state
box uses operational flow/head limits and exact efficiency-curve extrema,
including cubic stationary points and secant extensions. Positive minimum power
also supplies a conservative minimum head. Off-state bounds remain unchanged.
Unsupported negative electrical-efficiency domains receive no joint envelope.

Every physical on point has rank-one corner weights that reproduce its product;
an off point has zero mass and zero on-moments. The supplement therefore retains
feasible schedules and adds no binary decisions. It strengthens the relaxation
of the original nonlinear equality. Numerical guards have the same certificate
scope as the existing polynomial bounds. See the
[coupled relaxation benchmarks](coupled-relaxations.md) for the measured benefit
and small-case cost.

Table-coordinate power bounds also couple the complete power polynomial to the
existing head and discharge weights. Write
`F(q,h)=0.00981 e_max q h eta_turbine(q,h)` on the conditional running domain.
For each fixed slope `a`, head-node coefficients start as guarded upper bounds
on `max_q(F(q,h_j)−a q)`. A Bernstein bound on each original cell's residual
against the nodal chord raises both endpoint coefficients enough to cover
intermediate head values. All head supports share one vector of on-state
weights `ζ`:

```math
0\le\zeta_j\le\mu_j,\qquad
\sum_j\zeta_j=u,\qquad
\sum_j(h_j-h_1)\zeta_j=h_{on}-h_1u,
\qquad P\le a q+\sum_j c_j\zeta_j.
```

When on, equal total masses force `ζ=μ`. When off, `ζ=0` while the original
head weights remain free to represent the physical head. Fractional commitments
must use consistent on-state weights across every support and the existing
head moment. This implies the previous independent scalar gate hulls, and can
exclude fractional points those hulls allowed.

Discharge-node coefficients similarly bound `max_h(F(q_i,h)−b h)`, with certified
cell corrections. When the discharge coordinates have a unique first node at
zero and strictly positive remaining nodes, the on-state weights are affine:

```math
\nu_1=\lambda_1+u-1\ge0,\qquad \nu_i=\lambda_i\ (i>1),
\qquad P\le b h_{on}+\sum_i d_i\nu_i.
```

These weights have mass `u` and reproduce discharge. They need no new variables.
Other coordinate domains use the general shared lift; fixed states and singleton
coordinates simplify directly. The formulation preserves off-state head freedom
and adds no binary decisions. Four initial slopes on each axis give eight
supporting power rows per applicable unit and interval. Merely joining raw nodal
bounds is unsafe when power peaks between table knots; the cell correction is
part of the construction.

A bounded root separator proposes additional slopes for selected unit/interval
coordinates. Positive price, interval duration and optimistic relaxed power
guide the selection. Every added row is separately certified over its complete
original conditional running domain, including the original PCHIP polynomial
and extensions. Ranking scores and relaxed equation residuals are not proof
bounds. Rows enter both SCIP's global cut pool and its current LP. The separator
uses the existing shared weights, cached coefficients and reusable native value
buffers. It adds no partitions or binary decisions. Callback failures withhold
the reported global bound while retaining independently audited schedules.

These bounds preserve the original nonlinear power equality. See the
[proof-speed experiments](proof-speed.md) for the selected configuration and
cost, and the earlier [table-coordinate comparison](table-power-bounds.md) for
the coefficient construction.

Complete starting schedules and returned solver solutions are checked for
SOS2 adjacency and unit integrality. Tiny solver leakage at an off unit is
removed before physical reconstruction; larger violations reject the controls.
Accepted off units therefore have exactly zero flow and power.

`diagnostics_path="scip.log"` enables native SCIP progress logging and writes
`scip.log.statistics.json`. Returned `scip_statistics` includes native LP and
plugin counters, including optimization-based bound tightening. Timers can
overlap and must not be summed as exclusive costs.
`scip_diagnostics` contains native root/final bounds, node counts, LP iterations
and solution counts. Root and displayed log bounds are observations; only the
final enclosing bound and independently reconstructed objective determine the
reported gap. Native infinity sentinels are missing bounds. An apparent optimum
without a finite first root LP bound, despite LP iterations, is withheld.
SCIP's numerical termination status is not an independent proof. A checked
feasible schedule above a purported upper bound rejects that bound.
