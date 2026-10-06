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
