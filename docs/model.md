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
binary, retaining the exact quadratic loss law. Conservative head and continuity bounds can prove a direction; in that case the model uses its signed quadratic directly and omits the direction binary. Reversible links retain both directions.

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
Head-dependent turbine flow envelopes apply to running units. The experimental tightened graph uses isolated zero-discharge cells for off operation; on-cell selectors sum to unit commitment. The experimental tightened table graph bounds efficiency using endpoint and cubic stationary values, including declared extrapolation. It scales clipped cell coordinates to `[0,1]`, uses only the active branch when commitment is fixed, and represents single-segment one-dimensional tables with their exact affine graph.

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

## Formulation diagnostics

The default `solve(...; formulation=:baseline)` uses capacity-based network domains and the original exact table graphs. The experimental `:domains` and `:tightened` alternatives propagate conservative flow, head and storage intervals through conservation, river laws and signed tunnel losses before constructing the bounded model. Both prune unreachable table cells; `:tightened` also changes table graphs and strengthens turbine on/off selection. Neither uses a candidate trajectory to restrict the feasible set. Floating-point outward slack protects numerical enclosures; the certificate remains a numerical solver certificate. Matched operating Tokke–Vinje runs regressed with these alternatives, so they remain opt-in.

`formulation=:tensor` keeps baseline network domains and replaces turbine cell selection with independent discharge and head SOS2 weights. At most two adjacent weights on each axis may be nonzero. Each head column is a discharge interpolant, and efficiency is the head-weighted sum of those columns. Bilinear tables use linear discharge columns. Discharge PCHIP tables add shared adjacent-weight cubic corrections, preserving the original polynomial on clipped intervals and the original secants outside the data range. Column bounds include cubic stationary points and declared head extrapolation. Head-dependent flow envelopes reuse the same head weights; units share those weights when their physical head variable and clipped knots coincide. SOS2 constraints still require discrete search even though the graph has no explicit turbine-cell binaries. Other one-dimensional curves retain their baseline graphs.

`formulation=:domains` applies network domains and cell pruning with the original table graphs, separating those effects in an ablation. `formulation=:baseline` reproduces the earlier capacity-based domains and table graphs for matched experiments on the same physical case. `diagnostics_path="scip.log"` enables native SCIP progress logging. Returned `scip_diagnostics` include native root/final bounds, node counts, LP iterations and solution counts. Root or displayed log bounds are diagnostic observations; only the enclosing final bound and independently audited objective determine the reported gap. Native infinity sentinels are reported as missing bounds.
