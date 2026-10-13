// A self-contained article: no external Typst packages or fonts are required.
#set document(title: "Nonlinear hydropower scheduling with coupled hydraulics and delayed river networks")
#set page(paper: "a4", margin: 25mm, numbering: "1")
#set text(font: "New Computer Modern", size: 10pt, lang: "en")
#set par(justify: true, leading: 0.55em, spacing: 0.7em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)", supplement: none)
#show heading.where(level: 1): set text(size: 12pt)
#show heading.where(level: 2): set text(size: 10.5pt)
#show math.equation: set block(above: 0.65em, below: 0.65em)
#show figure.caption: set text(size: 9pt)
#show ref: it => {
  if it.element != none and it.element.func() == math.equation {
    link(it.target, numbering("(1)", ..counter(math.equation).at(it.element.location())))
  } else {
    it
  }
}

#align(center, text(size: 17pt)[
  Nonlinear hydropower scheduling with coupled hydraulics
  and delayed river networks
])
#v(0.8em)
#block(inset: (x: 8mm))[
  #align(center)[*Abstract*]
  Short-term hydropower scheduling couples discrete operation with the movement
  of water and the conversion of hydraulic energy. We give a mixed-integer
  nonlinear formulation with finite hydraulic domains for reservoirs, shared pressurized
  waterways, generating units, discrete gates, instantaneous channels and
  delayed river networks. Storage evolves by an implicit midpoint balance.
  Pressure losses and efficiency-dependent generation remain nonlinear;
  supplied response curves have explicit interpolation semantics.
  Deterministic river transport preserves release pulses through confluences,
  whereas distributed transport uses conservative interval mixing.
  The solution method combines mixed-integer linear proposals, nonlinear
  dispatch, exact algebraic interpolation graphs and spatial branch-and-bound.
  Independent chronological integration separates physical acceptance from
  the objective-gap certificate of the finite optimization model. The
  formulation is implemented in Julia using open-source optimization libraries.
]

= The scheduling problem
A generating schedule must account for water volume and hydraulic head.
Discharging a reservoir removes water, changes its level, and changes
the energy obtained from later discharges. Units also have operating
ranges: an idle unit produces nothing, while a running unit cannot
usually sustain arbitrarily small output. Head-dependent generation and
these discrete operating domains have a long history in nonlinear
hydropower scheduling~@catalao.

We consider a deterministic horizon with given prices, inflows,
equipment curves and operating restrictions. The physical model is
deliberately quasi-steady. Pressure and discharge equilibrate within
each control interval; reservoir storage evolves between intervals.
Tunnel inertia, water hammer, free-surface river momentum equations,
pumping, reserve products and stochastic decisions are outside its
scope. River geometry is represented by supplied response functions,
rather than resolved during scheduling.

The aim is a single explicit physical optimization problem. Linear
approximations may help choose initial discrete decisions, but they do
not define the delivered schedule's hydraulic feasibility. The final
equations retain head-dependent generation, signed pressure losses and
contemporaneous flow-dependent transport. Exactness below means
equivalence to these declared functions and this time discretization,
not exact reproduction of a real watercourse.

= Volumes, heads and connections
<sec:water>
Let $0 = tau_0 < dots.h.c < tau_T$ be time edges in hours,
$Delta_t = tau_(t + 1) - tau_t$, and
$m_t = ( tau_t + tau_(t + 1) ) \/ 2$, for $t = 0 , dots.h , T - 1$.
Operating inputs are piecewise constant, with changes on control edges.
Use the sets $cal(R)$ of reservoirs, $cal(J)$ of pressure junctions,
$cal(B)$ of prescribed-head boundaries, $cal(E)$ of tunnels, $cal(C)$ of
instantaneous channels, $cal(G)$ of units and $cal(D)$ of delayed
reaches. Units belong to plants $cal(G)_p$. @tab:variables fixes
the principal variables and units.

#figure(
  table(
    columns: (auto, 1fr, auto),
    align: left,
    inset: (x: 4pt, y: 3pt),
    stroke: none,
    table.hline(stroke: 0.7pt),
    table.header([Variable], [Meaning], [Unit]),
    table.hline(stroke: 0.4pt),
    [$V_(i t), W_(r t)$], [Reservoir volume, river water in transit], [$"Mm"^3$],
    [$H_(n t), h_(p t)$], [Node head, turbine net head], [m],
    [$Q_(e t), X_(c t), q_(g t)$], [Tunnel, channel, unit discharge], [$"m"^3/"s"$],
    [$R_(r t), d_(r t)$], [Reach entry flow, reservoir withdrawal], [$"m"^3/"s"$],
    [$A_(r t)$], [Volume arriving from a delayed reach], [$"Mm"^3$],
    [$P_(g t)$], [Electrical power], [MW],
    [$u_(g t), U_(p t)$], [Unit and plant operating states], [${0, 1}$],
    [$b_(e t), a_(r t)$], [Tunnel position, reservoir outlet opening], [$[0, 1]$],
    table.hline(stroke: 0.7pt),
  ),
  kind: table,
  caption: [Physical variables. Storage and river inventory use time vertices;
    other indexed quantities use control intervals.],
) <tab:variables>

Head measures mechanical energy per unit weight. At a reservoir it is
the water-surface elevation; at a flooded pressure junction it includes
pressure head. Raising a parcel by one metre and losing one metre of
head to friction have the same energy scale. Relative velocities and
fixed local losses are absorbed into the supplied equipment
relationships.

The unit conversions are

$ kappa = frac(3600, 10^6) = 0.0036, quad
  gamma = frac(rho g, 10^6) = 0.00981 \
  (rho = 1000 "kg"/"m"^3, quad g = 9.81 "m"/"s"^2). $ <eq:units>

Thus $kappa Delta q$ is a volume in $"Mm"^3$, and $gamma q h$ is
hydraulic power in MW.

Define $F_(n t)$ as the net internal rate entering node $n$, before
reservoir natural inflow. With incidence coefficients positive into a
node,

$ F_(n t) = sum_e B_(n e)^E Q_(e t)
  + sum_c B_(n c)^C X_(c t) + sum_g B_(n g)^G q_(g t)
  + sum_r D_(n r) frac(A_(r t), kappa Delta_t)
  - sum_r S_(n r) d_(r t). $ <eq:incidence>

$D_(n r) = 1$ if reach $r$ delivers into node $n$, and $S_(n r) = 1$ if
a reservoir outlet at $n$ supplies it. Other entries are zero. Tunnel
and unit incidence uses the actual water destination: a river outfall
withdraws from its source and injects the named reach, without
immediately crediting the hydraulic target.

Reservoirs and pressure junctions obey, respectively,

$ V_(i,t+1) - V_(i t) = kappa Delta_t (I_(i t) + F_(i t)), quad
  V_(i 0) = V_i^"init". $ <eq:storage>

$ F_(j t) = 0 quad (j in cal(J)), quad
  underline(V)_(i t) <= V_(i t) <= overline(V)_(i t)
  quad (t = 0, dots, T). $ <eq:junction>

A junction cannot retain water, so its equation is a balance rather than
an evolution law. Boundaries exchange water with the exterior and have
fixed heads. Reservoir hydraulic heads use midpoint storage:

$ H_(i t) = Z_i (frac(V_(i t) + V_(i,t+1), 2)), quad
  underline(H)_j <= H_(j t) <= overline(H)_j, quad
  H_(b t) = H_b^"given". $ <eq:heads>

$Z_i$ is an increasing supplied level--volume table, or
$Z_i (V) = z_i + s_i V + c_i V^2$. Midpoint storage is substituted
into $Z_i$; averaging the two endpoint levels would give a different
nonlinear model.

#figure(
  block(width: 160mm, height: 37mm)[
    #let node(x, body) = place(top + left, dx: x, dy: 22mm,
      box(width: 24mm, height: 13mm, stroke: 0.5pt, radius: 2pt,
        align(center + horizon, text(size: 8pt, body))))
    #let arrow(x) = {
      place(top + left, dx: x, dy: 28.5mm, line(length: 10mm, stroke: 0.5pt))
      place(top + left, dx: x + 8.5mm, dy: 27.7mm,
        line(end: (1.5mm, 0.8mm), stroke: 0.5pt))
      place(top + left, dx: x + 8.5mm, dy: 29.3mm,
        line(end: (1.5mm, -0.8mm), stroke: 0.5pt))
    }
    #node(0mm, [Reservoir A \ $V_A, H_A$])
    #node(34mm, [Pressure \ junction $H_J$])
    #node(68mm, [Plant \ $q_1 + q_2$])
    #node(102mm, [Reach \ inventory $W$])
    #node(136mm, [Reservoir B \ $V_B, H_B$])
    #for x in (24mm, 58mm, 92mm, 126mm) { arrow(x) }
    #place(top + left, dx: 17mm, dy: 17mm, text(size: 8pt)[tunnel $Q$])
    #place(top + left, dx: 92mm, dy: 17mm, text(size: 8pt)[release])
    #place(top + left, dx: 127mm, dy: 17mm, text(size: 8pt)[arrival])
    #place(top + left, line(start: (148mm, 22mm), end: (148mm, 8mm),
      stroke: (thickness: 0.5pt, dash: "dashed")))
    #place(top + left, line(start: (148mm, 8mm), end: (80mm, 8mm),
      stroke: (thickness: 0.5pt, dash: "dashed")))
    #place(top + left, line(start: (80mm, 8mm), end: (80mm, 22mm),
      stroke: (thickness: 0.5pt, dash: "dashed")))
    #place(top + left, line(start: (79.2mm, 20.5mm), end: (80mm, 22mm), stroke: 0.5pt))
    #place(top + left, line(start: (80.8mm, 20.5mm), end: (80mm, 22mm), stroke: 0.5pt))
    #place(top + left, dx: 86mm, dy: 2mm, text(size: 8pt)[outlet head reference])
  ],
  caption: [Water destination and hydraulic reference have different roles.
    Solid arrows move water; the dashed connection supplies the plant's outlet
    head. Water reaches B after transport, although its level already affects
    generation.],
) <fig:network>

= Pressure links and generating units
== Tunnels and gates
For an oriented tunnel $e = (v , w)$, positive $Q_e$ flows from $v$
to $w$. At positive opening,

$ H_(v t) - H_(w t) = K_e (b_(e t)) Q_(e t) abs(Q_(e t)), \
  K_e (b) = cases(
    R_e / b & "without a gate-loss curve",
    R_e + G_e (b) & "with a supplied gate-loss curve".
  ) $ <eq:tunnel>

$R_e > 0$ and $G_e >= 0$ have units $"s"^2/"m"^5$. The signed
product ensures that friction opposes flow in either direction.
Continuity and these pressure equations determine flow sharing through
branches and loops; there is no independent allocation rule.

A closed gate instead imposes

$ b_(e t) = 0 ==> Q_(e t) = 0, quad
  abs(Q_(e t)) <= C_(e t). $ <eq:closed>

It imposes no equality of endpoint heads. A directed river outfall also
has $Q_(e t) >= 0$. Tunnel positions are either prescribed or
selected from a finite declared set $cal(S)_e subset [ 0 , 1 ]$;
arbitrary continuous tunnel positions are not optimization decisions.

#block(breakable: false)[
An instantaneous channel between reservoirs has a different supplied
law:

$ X_(c t) = Phi_c (H_(v t) - H_(w t), H_("ref",t)), quad
  abs(X_(c t)) <= C_(c t). $ <eq:channel>
]

The reference is explicitly either endpoint level. Each constituent
discharge curve is nondecreasing, passes through $(0 , 0)$ and
preserves the sign of the level difference. Responses may be asymmetric.
Interpolation is linear in level difference and between reference
levels; endpoint segments continue linearly in difference, while the
nearest reference curve is held outside the reference span. Channels
have neither delay nor transit inventory.

== Energy conversion
Let $q_(p t) = sum_(g in cal(G)_p) q_(g t)$ and
$P_(p t) = sum_(g in cal(G)_p) P_(g t)$. All units at plant $p$ share

$ h_(p t) = H_("source"(p),t)
  - max{H_("target"(p),t), h_p^"floor"}
  - T_p (q_(p t)). $ <eq:nethead>

Omit the floor when absent. $T_p >= 0$ is an additional tailwater
rise, indexed by total plant discharge; it is not the complete
downstream elevation. Shared penstocks are separate tunnels, whose loss
therefore depends on the combined flow through the penstock.

Unit operation is a disjunction,

$ u_(g t) in {0, 1}, quad
  underline(q)_(g t) u_(g t) <= q_(g t) <= overline(q)_(g t) u_(g t), quad
  underline(P)_(g t) u_(g t) <= P_(g t) <= overline(P)_(g t) u_(g t). $ <eq:onoff>

$ P_(g t) = gamma q_(g t) h_(p t)
  eta_g (q_(g t), h_(p t)) e_g (P_(g t)). $ <eq:power>

On, additionally enforce

$ underline(h)_g <= h_(p t) <= overline(h)_g, quad
  eta_g^min <= eta_g (q_(g t), h_(p t)) <= 1, \
  q_g^"lo" (h_(p t)) <= q_(g t) <= q_g^"hi" (h_(p t)). $ <eq:unitdomain>

The head and efficiency bounds always apply to running units; the flow
envelopes apply when supplied by a turbine table. Off, $q = P = 0$ and
these operating inequalities are inactive; hydraulic heads remain in the
network.

Turbine efficiency is either a table or the analytic relationship

$ eta_g (q,h) = eta_g^0
  - alpha_g (frac(q - q_g^*, q_g^*))^2
  - beta_g (frac(h - h_g^*, h_g^*))^2. $ <eq:analyticeta>

Tables are linear in head and linear or shape-preserving cubic in
discharge. Electrical efficiency $e_g$ is a separate piecewise-linear
function of _electrical output_. Hence @eq:power is implicit in
$P$. The accepted electrical curves make $P \/ e_g (P)$ strictly
increasing, so there is at most one admissible scalar conversion. If no
curve is supplied, $e_g = 1$. This distinguishes turbine shaft
conversion from electrical losses.

Plant state is the logical OR of its units:

$ U_(p t) >= u_(g t) quad (g in cal(G)_p), quad
  U_(p t) <= sum_(g in cal(G)_p) u_(g t), quad U_(p t) in {0, 1}. $ <eq:plantstate>

Enforce $P_(p t) <= overline(P)_(p t)$,
$q_(p t) <= overline(q)_(p t)$ when supplied, and
$P_(p t) >= underline(P)_(p t) U_(p t)$,
$q_(p t) >= underline(q)_(p t) U_(p t)$.

= Reservoir outlets and river transport
== Entry laws and confluences
A reservoir-source reach has total entry $R_(r t) = n_(r t) + d_(r t)$,
where $n_(r t)$ is its natural inflow and $d_(r t)$ is withdrawal from
the reservoir. One chosen outlet law is

$ d_(r t) = cases(
    (C_(r t) - n_(r t)) a_(r t) & "controlled outlet",
    c_r a_(r t) sqrt(ell_(r t)) & "orifice",
    c_r ell_(r t)^(3/2) & "uncontrolled weir",
    a_(r t) D_r (H_("source"(r),t)) & "supplied discharge curve".
  ) $ <eq:outlets>

Here
$ell = H_("source") - z_r^("crest") >= 0$, or
$ell = [H_"source" - z_r^"crest"]_+$ if dry
operation is allowed. An uncontrolled weir has $a = 1$, including when a
discharge table supplies its law. Other openings satisfy their declared
bounds in $[ 0 , 1 ]$. Ordinary reaches obey
$n_(r t) <= R_(r t) <= C_(r t)$. Subtracting $d = R - n$, rather
than $R$, in @eq:incidence counts natural river inflow exactly once.

#block(breakable: false)[
For a confluence-fed reach, define local injection

$ L_(r t) = n_(r t) + sum_(g -> r) q_(g t) + sum_(e -> r) Q_(e t), quad
  R_(r t) = L_(r t) + sum_(s -> r) frac(A_(s t), kappa Delta_t). $ <eq:confluence>
]

Any number of upstream reaches and outfalls may merge. A confluence has
one outgoing reach; it does not decide how to split a stream. Delayed
reaches form a directed acyclic graph. Hydraulic loops through tunnels
and instantaneous channels remain admissible.

== Cohort response functions
Travel time separates release from arrival. The response-function
viewpoint is classical~@dooge; here its precise discretization is part
of the model. Reference flows are nonnegative, strictly ordered and
cover admissible entry flows when more than one distribution is
supplied. Delay edges are nonnegative and strictly increasing. For
reference flow $rho_(r j)$, a supplied delay distribution is uniform on
each bin $[ delta_(r j ell) , delta_(r j , ell + 1) ]$ with
probability $w_(r j ell) >= 0$, $sum_ell w_(r j ell) = 1$. Its
cumulative distribution is

$ F_(r j) (s) = sum_ell w_(r j ell)
  min{1, max{0, frac(s - delta_(r j ell), delta_(r j,ell+1) - delta_(r j ell))}}. $ <eq:cdf>

With two or more reference flows, the mixture weights satisfy

$ theta_(r k j) >= 0, quad sum_j theta_(r k j) = 1, quad
  R_(r k) = sum_j rho_(r j) theta_(r k j), quad
  theta_(r k) in op("SOS2"). $ <eq:delayweights>

Thus only neighboring distributions are blended. With a single
distribution, $theta_(r k 1) = 1$ and there is no coordinate
reproduction equation: its reference flow does not fix the release. The
selected mixture stays with that cohort. Operating capacity changes do
not rescale calibration reference flows.

For a rectangular cohort released during $[ tau_k , tau_(k + 1) )$,
the fraction arriving during interval $t$ at reference $j$ is

$ B_(r t k j) = frac(1, Delta_k) integral_(tau_k)^(tau_(k+1))
  [F_(r j) (tau_(t+1) - s) - F_(r j) (tau_t - s)] dif s. $ <eq:transfer>

The arrival volume and instantaneous arrival rate are

$ A_(r t) = A_(r t)^"hist"
  + kappa sum_k Delta_k R_(r k) sum_j theta_(r k j) B_(r t k j). $ <eq:arrival>

$ f_r (s) = f_r^"hist" (s) + sum_k R_(r k) sum_j theta_(r k j)
  [F_(r j) (s - tau_k) - F_(r j) (s - tau_(k+1))]. $ <eq:pointarrival>

The integrals are evaluated exactly for the declared uniform bins.
Interpolating the _fraction_ and then multiplying by released flow
makes the transfer piecewise quadratic in $R$; interpolating already
released volumes would define a different response.

If every delayed reach has deterministic delay $delta_r >= 0$, use
event-shaped entry streams $tilde(R)_r (s)$ and

$ f_r (s) = tilde(R)_r (s - delta_r), quad
  tilde(R)_r (s) = L_r (s) + sum_(s' -> r) f_(s') (s), \
  A_(r t) = kappa integral_(tau_t)^(tau_(t+1)) f_r (s) dif s. $ <eq:deterministic>

at confluences for $s >= 0$. For $s < 0$, $tilde(R)_r (s)$ is the
independently supplied entry history of that reach, not a backward
evaluation of the confluence recursion. On $[ tau_t , tau_(t + 1) )$,
a reservoir-source stream is $tilde(R)_r (s) = n_(r t) + d_(r t)$,
including natural inflow. Link event streams to interval variables by

$ R_(r t) = frac(1, Delta_t) integral_(tau_t)^(tau_(t+1)) tilde(R)_r (s) dif s. $ <eq:eventmean>

Original pulses are shifted through successive confluences; they are not
replaced by interval averages. Zero and sub-interval delays therefore
produce genuine same-interval coupling. Entry-capacity restrictions are
also enforced on the event stream, including both sides of its jumps.

If any reach uses a distributed response, @eq:confluence instead mixes
arrivals into new interval-average downstream cohorts. Equations
@eq:transfer and @eq:arrival remain exact for each such cohort, but
downstream timing is a grid approximation. Refining the mixing grid and
comparing cumulative arrivals and inventory estimates this error. Mass
conservation alone cannot detect an incorrectly spread pulse.

== Inventory, history and environmental flows
#block(breakable: false)[
Released water that has not arrived is still in the system:

$ W_(r,t+1) = W_(r t) + kappa Delta_t R_(r t) - A_(r t), quad
  W_(r 0) = W_r^"hist". $ <eq:inventory>
]

Historical entry cohorts supply initial inventory and future historical
arrivals: the historical terms in @eq:arrival and @eq:pointarrival use
the same formulas over supplied past intervals. For a historical cohort
$k$ on $[ s_k , e_k )$, let $F_k$ be its selected delay mixture. Then

$ W_r^"hist" = kappa sum_(k in cal(H)_r) R_(r k)
  integral_(s_k)^(e_k) [1 - F_k (-s)] dif s. $ <eq:historyinventory>

Unspecified earlier entry is zero in the declared model. Each reach's
history already includes its pre-horizon upstream contributions;
regenerating those contributions would double-count water. Restart
retains original executed cohorts and their release flows, as well as
reservoir states and operating histories. An average arrival bucket is
insufficient to reconstruct future arrival timing.

*Conservation identity.* Summing @eq:storage and @eq:inventory
cancels every internal transfer. Thus

$ sum_i (V_(i,t+1) - V_(i t)) + sum_r (W_(r,t+1) - W_(r t)) \
  = kappa Delta_t (sum_i I_(i t) + sum_r n_(r t) + E_t). $ <eq:globalwater>

where $E_t$ is net boundary exchange, positive into the modeled system.
For example, a plant outfall withdraws $q$ from its source, adds $q$ to
reach entry, and later leaves the reach as $A$; only that arrival enters
the receiving reservoir. Equation @eq:globalwater follows by this same
cancellation for each arc, including reverse instantaneous flow.

Environmental restrictions may apply to reach entry
$R_(r t) >= R_(r t)^min$, window-average arrival, or instantaneous
arrival:

$ integral_J f_r (s) dif s >= integral_J f_r^min (s) dif s, quad
  f_r (s) >= f_r^min (s). $ <eq:environment>

Only the selected arrival policy applies. For prescribed rectangular
cohorts, distributed arrivals are piecewise linear; deterministic
arrivals are piecewise constant. Checking response and requirement
breakpoints and the relevant one-sided limits enforces the pointwise
requirement for that transport model. At the horizon endpoints use the
interior limits.

A flow observation is an additional inequality, not a water injection:

$ sum_(g in cal(G)_o) q_(g t)
  + sum_(r in cal(D)_o union cal(C)_o) R_(r t) + I_(o t) >= M_(o t). $ <eq:observation>

For an instantaneous channel, its observed release is
$R_(c t) = X_(c t)$. Named contributions must refer to water already
represented in the network.

= Chronological operating constraints
Define unit starts and stops relative to the supplied initial state:

$ u_(g t) - u_(g,t-1) = s_(g t)^+ - s_(g t)^-, quad
  s_(g t)^+ + s_(g t)^- <= 1, quad s_(g t)^(plus.minus) in {0, 1}. $ <eq:transitions>

For every start at $tau_t$, enforce $u_(g k) >= s_(g t)^(+)$ whenever
$0 <= tau_k - tau_t < L_g^("up")$; for a stop, enforce
$u_(g k) <= 1 - s_(g t)^(-)$ whenever
$0 <= tau_k - tau_t < L_g^("down")$. Initial state age
fixes the remaining dwell interval. Apply the same constraints to
aggregate plant state $U$ when plant dwell is declared. Changing the
running member does not restart an already operating plant. Dwell beyond
the horizon is carried to a restart, rather than forced by an invented
terminal state.

Availability, forced states and prescribed output are hard restrictions:

$ underline(u)_(g t) <= u_(g t) <= overline(u)_(g t), quad
  P_(g t) = P_(g t)^"given", quad q_(g t) = q_(g t)^"given". $ <eq:schedules>

where supplied. Plant state also satisfies
$underline(U)_(p t) <= U_(p t) <= overline(U)_(p t)$; plant output
schedules constrain sums. Maintenance forces the relevant units off. A
zero schedule implies off; a positive schedule implies on. All
restrictions intersect, including initial dwell; no input silently
overrides another.

#block(breakable: false)[
For interval-average discharge or release $x_t$, rate limits give

$ -integral_(m_(t-1))^(m_t) r_x^- (s) dif s
  <= x_t - x_(t-1) <= integral_(m_(t-1))^(m_t) r_x^+ (s) dif s. $ <eq:ramps>
]

Unit and plant power use the same bounds with additional allowed jumps
$sum_(g in cal(G)_x) underline(P)_(g t) s_(g t)^(+)$ upward and
$sum_(g in cal(G)_x) underline(P)_(g , t - 1) s_(g t)^(-)$ downward.
These allow the minimum running output at a transition, not unrestricted
startup or shutdown. Discharge and river ramps have no jump exemption.
The first comparison uses a supplied previous-window average and
midpoint. If absent, only that first comparison is omitted.

Reservoir volume ramps use $V_(i , t + 1) - V_(i t)$ and level ramps
use $Z_i ( V_(i , t + 1) ) - Z_i ( V_(i t) )$, with rate integrals
over $[ tau_t , tau_(t + 1) ]$. They constrain vertex changes, not
continuous-time derivatives.

For a tunnel interlock $cal(E)_ell$, with
$o_(e t) = upright(bold(1)) [ b_(e t) > 0 ]$, require

$ underline(N)_(ell t) <= sum_(e in cal(E)_ell) o_(e t) <= overline(N)_(ell t). $ <eq:interlock>

A positive position counts as open even when its flow is zero.

= Economic objective
With price $pi_t$ per MWh, reservoir and reach values $v_i , v_r$ per
$"Mm"^3$, and nonnegative movement and release-shortfall costs, maximize

$ cal(V) = & sum_t pi_t Delta_t sum_g P_(g t)
  - sum_(g,t) (c_(g t)^+ s_(g t)^+ + c_(g t)^- s_(g t)^-) \
  & - sum_(e,t) c_(e t)^b abs(b_(e t) - b_(e,t-1))
  - sum_(r,t) c_(r t)^r xi_(r t) \
  & + sum_i v_i (V_(i T) - V_(i 0)) + sum_r v_r (W_(r T) - W_(r 0)). $ <eq:objective>

The first gate comparison uses its given initial position. Transition
and movement costs are charged once, without a duration multiplier. A
positive release penalty makes the entry minimum soft:

$ xi_(r t) >= 0, quad xi_(r t) >= kappa Delta_t (R_(r t)^min - R_(r t)). $ <eq:shortfall>

$xi$ is a volume in $"Mm"^3$; positive cost makes it the actual
shortfall. Without a positive penalty, the entry minimum remains hard
and $xi = 0$. Arrival minima and aggregate observations remain hard.

Terminal transit value prevents a release near the horizon from losing
all economic value merely because it has not yet arrived. These values
are explicit input assumptions. Subtracting initial inventories removes
a constant endowment from both reported value and its relative-gap
denominator.

= Finite algebraic representation
<sec:graphs>
The constraints above define a nonconvex MINLP. Discharge changes
immediate power, future head and sometimes travel response. Neither
fixed commitment nor exact interpolation makes these couplings convex.

== Exact interpolation graphs
For a piecewise-linear curve with knots $(x_i , y_i)$, use

$ lambda_i >= 0, quad sum_i lambda_i = 1, quad
  x = sum_i x_i lambda_i, quad y = sum_i y_i lambda_i, quad lambda in op("SOS2"). $ <eq:sos>

SOS2 means that at most two nonzero weights are allowed, and their
indices must be adjacent. It is a disjunctive interpolation
constraint~@huch. On one segment the coordinate equation determines
those two weights, so @eq:sos gives the original curve value exactly.

Turbine tables use independent discharge weights $lambda_i$ and head
weights $mu_j$. Let $E_(i j)$ be knot efficiencies. A shape-preserving
cubic discharge interpolant~@fritsch has the representation

$ w_i = lambda_i lambda_(i+1), quad r_i = lambda_i w_i, \
  eta_j = sum_i E_(i j) lambda_i
    + sum_i [A_(i j) r_i + B_(i j) (w_i - r_i)], quad
  eta = sum_j mu_j eta_j. $ <eq:cubic>

For segment width $delta q_i = q_(i + 1) - q_i$, secant
$d_(i j) = ( E_(i + 1 , j) - E_(i j) ) \/ delta q_i$ and original
endpoint slopes $m_(i j)$, the coefficients are

$ A_(i j) = delta q_i (m_(i j) - d_(i j)), quad
  B_(i j) = delta q_i (d_(i j) - m_(i+1,j)). $ <eq:hermitecoefficients>

On an active pair, $lambda_i + lambda_(i + 1) = 1$, so the correction
basis is $lambda_i^2 lambda_(i + 1)$ and $lambda_i lambda_(i + 1)^2$,
which spans the cubic corrections that vanish at both endpoints. Off the
active pair these products vanish. For bilinear tables $A = B = 0$.
Clipping or extending a domain does not refit the supplied slopes; outer
discharge extensions use their declared secants. Head extrapolation,
where explicitly allowed, uses the outer head secant. Off-state
efficiency extensions are finite algebraic continuations, not additional
running-domain restrictions.

Compatible units share head coordinates, and each cohort shares one flow
coordinate across arrivals, pointwise requirements and terminal
inventory. The same construction represents piecewise-bilinear channel
responses. Auxiliary equalities in @eq:cubic, the power equation and
released transport products remain nonlinear. SOS2 has not turned the
problem into a MILP. These graphs preserve the specified functions;
valid inequalities and branching improve the relaxation used to search
them.

== Quadratic tunnel and gate graphs
Signed losses can be expressed with a direction variable $z_(e t)$:

$ Q = Q^+ - Q^-, quad 0 <= Q^+ <= C z, quad 0 <= Q^- <= C(1-z), \
  S = (Q^+)^2 - (Q^-)^2, quad z in {0, 1}. $ <eq:signedloss>

Then $S = Q abs(Q)$. Directed outfalls need only $S = Q^2$; known head
ordering can fix a direction. Finite gate selectors give both $b$ and
its open state $o$ as affine expressions. Define
$d_o = o (H_v - H_w)$. For bounds $L <= H_v - H_w <= U$, its
exact binary-product hull is

$ L o <= d_o <= U o, quad
  H_v - H_w - U(1-o) <= d_o <= H_v - H_w - L(1-o). $ <eq:binaryhull>

Together with $abs(Q) <= C o$ and $d_o = K (b) S$, this retains
@eq:tunnel and @eq:closed. $K (b)$ is affine in the finite position
selectors, with any finite coefficient assigned to the closed position
since $S = 0$ there. Only quadratic product equalities are needed.
Wet/dry outlet branches similarly represent
$ell = max ( H - z^("crest") , 0 )$; an auxiliary root
$r >= 0$, $r^2 = ell$, represents the square root in @eq:outlets.

== Bounds that retain the physical equations
Spatial branch-and-bound needs enclosing relaxations. Finite equipment
and storage domains, conservation and forward/backward inventory bounds
supply domains without restricting them to a previously found schedule.
Polynomial bounds include interior extrema, not only table knots.

Let $overline(e)$ bound electrical efficiency and define the
running-domain power majorant
$p^("ub") (q , h) = gamma overline(e) thin q h eta (q , h)$,
with nonnegative $q , h , eta$. It is polynomial on each interpolation
cell, even when the electrical conversion is implicit. If
$p^("ub") (q , h) <= a q + b h + c$, with $b >= 0$
and a full head lower bound $h_(*)$, then

$ P <= a q + b(h - h_*) + (c + b h_*) u. $ <eq:support>

is valid for both states: on it bounds physical power; off $P = q = 0$
and the right side is nonnegative. Bernstein polynomial enclosures can
bound the residual over each original interpolation cell to establish
$c$. A plane fitted only at samples does not establish this inequality.

A joint box envelope couples discharge, running head and efficiency. Use
exact binary-product hulls for $h_("on") = u h$ and
$eta_("on") = u eta$. Nonnegative corner weights $z_c$ satisfy

$ sum_c z_c = u, quad sum_c q_c z_c = q, quad
  sum_c h_c z_c = h_"on", quad sum_c eta_c z_c = eta_"on", \
  gamma e_min sum_c q_c h_c eta_c z_c
  <= P <= gamma e_max sum_c q_c h_c eta_c z_c. $ <eq:box>

Up to eight corners describe a nonnegative running box, and
$e_min , e_max$ bound electrical efficiency on its power interval. A
physical on point has product corner weights reproducing $q h eta$; off
all weights vanish. The exact efficiency and power equalities remain.

Table-aware supports also share on-state head weights $zeta_j$:

$ 0 <= zeta_j <= mu_j, quad sum_j zeta_j = u, quad
  sum_j (h_j - h_1) zeta_j = h_"on" - h_1 u, \
  P <= a q + sum_j c_j zeta_j. $ <eq:coordinatecuts>

Coefficients enclose the complete polynomial power majorant, including
corrections between head nodes. On, $zeta = mu$; off, $zeta = 0$ while
the hydraulic head stays free. Analogous discharge-coordinate supports
couple the same power law along the other axis. A bounded separator may
add further valid supports. These constraints tighten relaxations; none
replaces a physical equality or adds a new operating decision.

= Solution, reconstruction and certificate
JuMP supplies the algebraic Julia interface~@jump. SCIP solves the
global nonconvex model by branch-and-cut with spatial branching,
nonlinear relaxations and bound tightening~@scip. A preliminary MILP
proposes unit states and finite gate positions using approximate
generation and hydraulics. With those discrete decisions fixed, Ipopt
solves nonlinear dispatch~@ipopt. Several proposals can be dispatched
under a shared time allowance. Only proposals passing reconstruction and
finer replay are submitted as full-model starts. The MILP objective is a
surrogate, never an enclosing bound for @eq:objective. Fixed-discrete
dispatch is a piecewise-smooth nonconvex NLP, and a stationary point is
not a global optimum.

Physical reconstruction then holds unit discharges and openings fixed
and advances the coupled equations chronologically. For each interval, a
Newton solve includes the next storage, pressure-junction heads, tunnel
flows and active outlet releases. Current arrivals remain in the solve
when travel time overlaps the same interval. Unobservable pressure
components receive one numerical reference while retaining every water
balance; admissible common pressure shifts respect junction bounds.
Electrical output is recovered from @eq:power. Reusable numerical arrays
separate this equation evaluation from input parsing and solver
reporting.

Acceptance checks the original interval equations and operating
restrictions, then repeats integration on a finer grid with the original
controls. Schedules and ramps are compared over their original operating
windows; storage vertices, total water, cumulative arrivals and terminal
transit are checked separately. This prevents an apparently converged
pressure solve from concealing an operating violation or a transport
error. For prescribed tolerances $epsilon_V , epsilon_(cal(V))$,
consistency also requires

$ max_(i,t) abs(V_i^"fine" (tau_t) - V_(i t)) <= epsilon_V, quad
  max_r abs(W_r^"fine" (tau_T) - W_(r T)) <= epsilon_V, \
  abs(cal(V)^"fine" - cal(V)) <= epsilon_(cal(V)) max(1, abs(cal(V))). $ <eq:replaycriteria>

Refinement supplies a numerical consistency test, not a rigorous
enclosure of continuous-time hydraulic error.

Let $L$ be the independently reconstructed value of a schedule feasible
for the declared interval model and $U$ a finite enclosing solver upper
bound. The reported gaps are

$ g_"abs" = max(0, U-L), quad g_"rel" = frac(max(0, U-L), max(1, abs(L))). $ <eq:gap>

A numerical certificate requires an enclosing bound and the requested
absolute or relative gap. It concerns the declared interval model,
interpolation and finite gate sets. Fixing discrete decisions makes it
conditional on those decisions. The interval-model certificate and finer
replay acceptance are separate outcomes: a schedule with a certified
objective gap can still fail finer replay. Neither certifies a
continuous-time optimum. Floating-point solver tolerances and guarded
polynomial bounds do not constitute an interval-arithmetic proof.

= Discussion
The useful separation is between the physical model, methods for finding
feasible controls, and methods for enclosing their optimum. Direct
nonlinear dispatch keeps head and efficiency changes inside the
scheduling problem; cohort accounting preserves water and its timing;
exact interpolation graphs retain the supplied response rather than
silently fitting another one. Independent reconstruction gives the
delivered controls their physical meaning, while global bounds quantify
unresolved economic improvement.

These choices do not remove nonconvexity, make coarse river mixing
exact, or guarantee a short proof on a large network. Their benefit is a
model whose approximations and certificates can be stated together. Site
calibration, unmodeled hydraulic regimes and uncertainty remain separate
questions. The formulation can be extended by adding explicit physical
laws and operating constraints while retaining the same conservation,
interpolation and acceptance principles.

#block(breakable: false)[
  #bibliography("references.bib", style: "ieee", title: "References")
]
