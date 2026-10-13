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
  A hydropower schedule chooses when to run generating units and how much water
  to release. Each release changes both the available water and the head for
  subsequent generation; in a connected watercourse it also changes downstream
  supply after a travel delay. We formulate these decisions as a mixed-integer
  nonlinear optimization problem with reservoirs, shared tunnels, finite gate
  choices and river networks. Reservoir balances use midpoint storage to
  evaluate head. Pressure losses and efficiency-dependent power remain
  nonlinear. River transport accounts for both arrivals and water still in
  transit, with an explicit distinction between fixed travel times and
  distributed delays. Supplied curves are represented by their specified
  interpolants. Feasible schedules are found using discrete proposals and
  nonlinear dispatch; global search provides an upper bound on their possible
  improvement. Independent integration checks the returned controls on the
  optimization grid and a finer grid. The model and solution method are
  implemented in Julia with open-source solvers.
]

= The scheduling problem

A reservoir's water can produce electricity now or be retained for later
use. Releasing it changes the reservoir level, and therefore changes the
energy available from the next release. In a cascade the same water may
pass through several plants, with pressure losses along tunnels and travel
delays along rivers. These effects connect decisions across the network and
over time. Unit operation adds a discrete choice: an idle unit produces
zero, whereas a running unit must remain within its operating range.
Models combining these choices with head-dependent generation have a long
history~@catalao.

The given data are the network, initial volumes and operating histories,
prices, inflows, equipment curves and operating limits. The decisions are
unit states, unit discharges, outlet openings and tunnel positions chosen
from finite sets. Reservoir volumes, tunnel flows, hydraulic heads, river
arrivals and electrical output must then satisfy the equations below. The
objective is revenue from generation, less operating penalties, plus the
value assigned to water remaining at the end of the horizon.

We assume deterministic inputs and quasi-steady hydraulics: tunnel pressures
and flows are in equilibrium within each scheduling interval, while reservoir
storage changes over time. This excludes tunnel inertia and water hammer.
River geometry enters through supplied discharge and travel-time response
curves; the optimization does not solve river momentum equations. Pumping,
reserve products and stochastic decisions are outside the present model.

The final formulation retains pressure losses, changing head and efficiency
in its nonlinear equations. Linear approximations help propose initial
discrete decisions. We first state the water balances, hydraulic laws and
operating constraints, then explain how they are represented and solved.
An exact representation of a supplied curve means the same interpolated
function; the curve itself remains an approximation of the watercourse.

= Water balance and hydraulic head
<sec:water>

Let $0=tau_0<dots.h.c<tau_T$ be time edges in hours,
$Delta_t=tau_(t+1)-tau_t$ the interval duration, and
$m_t=(tau_t+tau_(t+1)) \/ 2$ its midpoint, for $t=0,dots,T-1$.
Interval inputs are piecewise constant, with changes at these edges.
Volumes are recorded at the edges; the hydraulic equations describe each
interval using its midpoint reservoir volumes.

Hydraulic head is mechanical energy per unit weight, expressed in metres.
At a reservoir it is the water-surface elevation. At a flooded junction it
also includes pressure divided by $rho g$, where $rho$ is water density and
$g$ is gravitational acceleration. A head difference drives flow: some of
that energy is dissipated by friction and some may be converted to electricity.
Velocity corrections and local losses are included in the supplied equipment
relationships rather than modeled as separate states.

#block(breakable: false)[
We measure volume in
$"Mm"^3$, here denoting one million cubic metres ($10^6 "m"^3$), discharge
in $"m"^3/"s"$ and time in hours. The conversion factors are

$ kappa = frac(3600, 10^6) = 0.0036, quad
  gamma = frac(rho g, 10^6) = 0.00981 \
  (rho = 1000 "kg"/"m"^3, quad g = 9.81 "m"/"s"^2). $ <eq:units>
]

Thus $kappa Delta q$ converts a discharge over an interval into a volume,
and $gamma q h$ is hydraulic power at discharge $q$ and net head $h$.
For example, a release of $10 "m"^3/"s"$ for one hour removes
$0.036 "Mm"^3$ of water. A reservoir's volume change is this conversion
applied to its net inflow.

To write the same accounting for a network, let $cal(R)$, $cal(J)$ and
$cal(B)$ index reservoirs, pressure junctions and prescribed-head boundaries.
A node $n$ is any of these three. Tunnels, instantaneous channels, generating
units and delayed river reaches have index sets $cal(E)$, $cal(C)$, $cal(G)$
and $cal(D)$, respectively. The units at plant $p$ form $cal(G)_p$.
@tab:variables lists the variables. Underlines and overlines denote supplied
lower and upper bounds.

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

Let $F_(n t)$ be net flow into node $n$ from the modeled connections,
excluding natural reservoir inflow $I_(i t)$. The matrices $B^E$, $B^C$ and
$B^G$ record the tunnel, channel and unit connections. For a connection
between nodes, its column is $-1$ at the source, $+1$ at the destination,
and zero elsewhere. Negative discharge reverses these contributions.
Arrival volumes from rivers and withdrawals into them give the remaining terms:

$ F_(n t) = sum_e B_(n e)^E Q_(e t)
  + sum_c B_(n c)^C X_(c t) + sum_g B_(n g)^G q_(g t)
  + sum_r D_(n r) frac(A_(r t), kappa Delta_t)
  - sum_r S_(n r) d_(r t). $ <eq:incidence>

$D_(n r)=1$ when reach $r$ delivers into node $n$, and $S_(n r)=1$ when an
outlet at $n$ withdraws water into that reach; other entries are zero. Arrival
volume $A_(r t)$ is divided by $kappa Delta_t$ to give a flow rate.
A tunnel or plant discharging into a river has its withdrawal in $F$, while
its injection enters the river balance below. Its hydraulic outlet reference
need not be its immediate water destination: a downstream reservoir can
determine the plant's outlet head before the released water arrives there,
as @fig:network illustrates.

Water stored at the next time edge equals the present volume plus net inflow
over the interval. A pressure junction stores no water, so its net inflow
must be zero:

$ V_(i,t+1) - V_(i t) = kappa Delta_t (I_(i t) + F_(i t)), quad
  V_(i 0) = V_i^"init". $ <eq:storage>

$ F_(j t) = 0 quad (j in cal(J)), quad
  underline(V)_(i t) <= V_(i t) <= overline(V)_(i t)
  quad (t = 0, dots, T). $ <eq:junction>

Initial volumes and volume limits are given data. Junction balances determine
how flows meet and divide. Boundary nodes have prescribed heads and may
exchange water with the exterior. Reservoir head depends on midpoint volume:

$ H_(i t) = Z_i (frac(V_(i t) + V_(i,t+1), 2)), quad
  underline(H)_j <= H_(j t) <= overline(H)_j, quad
  H_(b t) = H_b^"given". $ <eq:heads>

$Z_i$ converts volume to surface elevation. It is an increasing supplied
table, or $Z_i (V)=z_i+s_i V+c_i V^2$ with given coefficients.
Equations @eq:storage and @eq:heads form an implicit midpoint step: the ending
volume helps determine the head that determines the interval flows.
For nonlinear $Z_i$, applying the function to the average volume differs
from averaging the two endpoint elevations.

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

= Tunnel hydraulics and generation

== Tunnels and gates

For an oriented tunnel $e=(v,w)$, positive $Q_e$ flows from $v$ to $w$.
With opening $b>0$, the head difference equals friction loss:

$ H_(v t) - H_(w t) = K_e (b_(e t)) Q_(e t) abs(Q_(e t)), \
  K_e (b) = cases(
    R_e / b & "without a gate-loss curve",
    R_e + G_e (b) & "with a supplied gate-loss curve".
  ) $ <eq:tunnel>

$R_e>0$ is tunnel resistance and $G_e>=0$ is any additional gate resistance;
both have units $"s"^2/"m"^5$. For positive $Q$, upstream head exceeds
downstream head; for negative $Q$, the difference reverses. The loss magnitude
grows as $Q^2$, and energy is dissipated in the direction of flow. These head
equations and the junction balances jointly determine flow sharing through
branches and loops.

A closed gate permits a head difference while stopping flow:

$ b_(e t) = 0 ==> Q_(e t) = 0, quad
  abs(Q_(e t)) <= C_(e t). $ <eq:closed>

Equation @eq:tunnel therefore applies only to positive openings. Capacity
$C_(e t)$ bounds flow in either direction; an outfall permitted to discharge
only into a river also has $Q_(e t)>=0$. Tunnel positions are prescribed
or chosen from a finite supplied set $cal(S)_e subset [0,1]$. A prescribed
position may be fractional; arbitrary continuous tunnel positions are not
optimization decisions in this model.

#block(breakable: false)[
An instantaneous channel between reservoirs follows a supplied discharge law:

$ X_(c t) = Phi_c (H_(v t) - H_(w t), H_("ref",t)), quad
  abs(X_(c t)) <= C_(c t). $ <eq:channel>
]

$Phi_c$ is a family of curves indexed by a specified upstream or downstream
reservoir level. Each curve gives discharge as a function of level difference,
is nondecreasing and passes through $(0,0)$. Positive and negative responses
need not be symmetric. We interpolate linearly along each curve and between
reference levels. Outside the tabulated differences we extend the end segments;
outside the reference levels we use the nearest curve. The channel transfers
water in the same interval and has no separate volume in transit.

== Energy conversion

Write $q_(p t)=sum_(g in cal(G)_p) q_(g t)$ and
$P_(p t)=sum_(g in cal(G)_p) P_(g t)$ for plant totals. All units at a plant
share the same net head:

$ h_(p t) = H_("source"(p),t)
  - max{H_("target"(p),t), h_p^"floor"}
  - T_p (q_(p t)). $ <eq:nethead>

The target head is the outlet reference in @fig:network. A supplied floor
prevents the effective outlet level from falling below $h_p^"floor"$; the
maximum is omitted when no floor is given. $T_p (q_(p t))>=0$ is the additional
tailwater rise caused by total plant discharge. Shared penstocks are separate
tunnels, so their friction depends on the combined flow of the units they supply.

The binary state $u_(g t)$ selects an idle unit or its running range:

$ u_(g t) in {0, 1}, quad
  underline(q)_(g t) u_(g t) <= q_(g t) <= overline(q)_(g t) u_(g t), quad
  underline(P)_(g t) u_(g t) <= P_(g t) <= overline(P)_(g t) u_(g t). $ <eq:onoff>

The electrical output is the available hydraulic power multiplied by turbine
and electrical efficiencies:

$ P_(g t) = gamma q_(g t) h_(p t)
  eta_g (q_(g t), h_(p t)) e_g (P_(g t)). $ <eq:power>

#block(breakable: false)[
When $u_(g t)=1$, the unit must also satisfy head and efficiency limits and,
when supplied by a turbine table, head-dependent discharge limits:

$ underline(h)_g <= h_(p t) <= overline(h)_g, quad
  eta_g^min <= eta_g (q_(g t), h_(p t)) <= 1, \
  q_g^"lo" (h_(p t)) <= q_(g t) <= q_g^"hi" (h_(p t)). $ <eq:unitdomain>
]

When $u_(g t)=0$, @eq:onoff gives $q=P=0$ and the operating inequalities in
@eq:unitdomain are inactive. The plant head is still determined by the
hydraulic network, even when none of its units runs.

Turbine efficiency is a supplied table or the analytic relationship

$ eta_g (q,h) = eta_g^0
  - alpha_g (frac(q - q_g^*, q_g^*))^2
  - beta_g (frac(h - h_g^*, h_g^*))^2. $ <eq:analyticeta>

Here $q_g^*$ and $h_g^*$ are reference operating conditions, $eta_g^0$ is
their turbine efficiency, and nonnegative $alpha_g,beta_g$ describe its
decrease away from them. Tables interpolate linearly in head and either
linearly or with shape-preserving cubics in discharge.

Turbine efficiency converts hydraulic power to shaft power. Electrical
efficiency $e_g$ accounts for the remaining conversion to electrical output.
It is a separate piecewise-linear function of _electrical output_ $P$, so
@eq:power must be solved for $P$ even when $q,h$ are known. Accepted curves
make $P \/ e_g (P)$ strictly increasing, giving at most one admissible output
for a given shaft power. With no electrical-efficiency curve, we use $e_g=1$.

Plant state $U_(p t)$ is one whenever at least one unit runs:

$ U_(p t) >= u_(g t) quad (g in cal(G)_p), quad
  U_(p t) <= sum_(g in cal(G)_p) u_(g t), quad U_(p t) in {0, 1}. $ <eq:plantstate>

Plant limits apply to aggregate output and discharge:
$P_(p t)<=overline(P)_(p t)$ and $q_(p t)<=overline(q)_(p t)$ when supplied,
and $P_(p t)>=underline(P)_(p t) U_(p t)$,
$q_(p t)>=underline(q)_(p t) U_(p t)$.

= Outlets and river transport

== Entry laws and confluences

A river _reach_ is a link with a defined entry and downstream destination.
For a reach fed by a reservoir outlet, total entry is
$R_(r t)=n_(r t)+d_(r t)$: natural reach inflow $n$ plus reservoir withdrawal
$d$. With outlet opening $a_(r t)$, one chosen law determines withdrawal:

$ d_(r t) = cases(
    (C_(r t) - n_(r t)) a_(r t) & "controlled outlet",
    c_r a_(r t) sqrt(ell_(r t)) & "orifice",
    c_r ell_(r t)^(3/2) & "uncontrolled weir",
    a_(r t) D_r (H_("source"(r),t)) & "supplied discharge curve".
  ) $ <eq:outlets>

$ell$ is the head above the outlet crest, $H_"source"-z_r^"crest"$.
It must be nonnegative unless dry operation is allowed, in which case we use
$[H_"source"-z_r^"crest"]_+$, with $[x]_+=max(x,0)$.
The coefficient $c_r$ or curve $D_r$ is supplied. An uncontrolled weir has
$a=1$, also when its law is given as a table. Other openings satisfy their
bounds in $[0,1]$. Ordinary reaches obey $n_(r t)<=R_(r t)<=C_(r t)$.
Only $d=R-n$ is subtracted from reservoir storage: the natural reach inflow
enters the river directly and has never been in the reservoir.

#block(breakable: false)[
At a confluence, water arriving from upstream reaches joins natural inflow
and any plant or tunnel outfalls. Let $L_(r t)$ collect these local injections:

$ L_(r t) = n_(r t) + sum_(g -> r) q_(g t) + sum_(e -> r) Q_(e t), quad
  R_(r t) = L_(r t) + sum_(s -> r) frac(A_(s t), kappa Delta_t). $ <eq:confluence>
]

The notation $g->r$, $e->r$ and $s->r$ identifies connections into reach $r$.
Any number of upstream reaches may merge into one outgoing reach. No allocation
decision is needed at such a junction. The delayed reach network is directed
and acyclic, so upstream responses can be propagated downstream in order.
The separate tunnel and instantaneous-channel network may contain loops.

== Release cohorts and travel-time responses

A _cohort_ is the water entering a reach during one release interval.
With a fixed travel time, that release simply appears downstream later.
With a distributed delay, different fractions of the cohort arrive at
different times. The distribution describes how water is spread in time;
it does not introduce uncertainty into the scheduling problem. Supplied
response functions take the place of a resolved river-flow model~@dooge.

For reference entry flow $rho_(r j)$, let fraction $w_(r j ell)$ have a
travel time uniformly distributed over
$[delta_(r j ell),delta_(r j,ell+1)]$. The delay edges are nonnegative and
strictly increasing; $w_(r j ell)>=0$ and $sum_ell w_(r j ell)=1$.
The cumulative fraction that has arrived after elapsed time $s$ is

$ F_(r j) (s) = sum_ell w_(r j ell)
  min{1, max{0, frac(s - delta_(r j ell), delta_(r j,ell+1) - delta_(r j ell))}}. $ <eq:cdf>

Entry flow may change the travel response. For two or more nonnegative, strictly increasing reference
flows covering the admissible entry range, interpolate between the neighboring
distributions. In the following equations, SOS2 is a _special ordered set of
type 2_: at most two weights may be nonzero, and their indices must be adjacent.

$ theta_(r k j) >= 0, quad sum_j theta_(r k j) = 1, quad
  R_(r k) = sum_j rho_(r j) theta_(r k j), quad
  theta_(r k) in op("SOS2"). $ <eq:delayweights>

The weights reproduce the cohort's total entry flow, including natural inflow
and upstream arrivals. For example, halfway between two reference flows,
the two distributions receive equal weights. With only one distribution,
$theta_(r k 1)=1$ and the flow-reproduction equation is omitted; the release
can still vary. The selected mixture remains attached to the cohort after
release. Changing an operating capacity does not change these reference flows.

Consider water entering at time $s$ in interval $k$. The fraction of it arriving
during interval $t$ is the cumulative fraction at age $tau_(t+1)-s$ minus
that at age $tau_t-s$. Averaging over the release interval gives

$ B_(r t k j) = frac(1, Delta_k) integral_(tau_k)^(tau_(k+1))
  [F_(r j) (tau_(t+1) - s) - F_(r j) (tau_t - s)] dif s. $ <eq:transfer>

Multiplying by released volume and summing over cohorts gives arrival volume
$A_(r t)$:

$ A_(r t) = A_(r t)^"hist"
  + kappa sum_k Delta_k R_(r k) sum_j theta_(r k j) B_(r t k j). $ <eq:arrival>

Water arriving at clock time $s$ from cohort $k$ must have a travel time between
$s-tau_(k+1)$ and $s-tau_k$. Integrating the delay density over this range
gives the difference of cumulative fractions below. Multiplying by the
cohort's release rate gives its contribution to instantaneous arrival rate:

$ f_r (s) = f_r^"hist" (s) + sum_k R_(r k) sum_j theta_(r k j)
  [F_(r j) (s - tau_k) - F_(r j) (s - tau_(k+1))]. $ <eq:pointarrival>

The historical terms account for releases before $tau_0$. The integrals are
exact for the supplied uniform bins. The interpolated quantity is the
_fraction_ that arrives. Multiplication by the flow-dependent released volume
then makes the transfer piecewise quadratic in $R$. Interpolating released
volumes directly would change this travel response.

When every reach has a fixed travel time $delta_r>=0$, we retain the full
piecewise-constant entry stream $tilde(R)_r (s)$ instead of averaging it at
each confluence. Shifting a stream by its delay and adding upstream arrivals gives

$ f_r (s) = tilde(R)_r (s - delta_r), quad
  tilde(R)_r (s) = L_r (s) + sum_(s' -> r) f_(s') (s), \
  A_(r t) = kappa integral_(tau_t)^(tau_(t+1)) f_r (s) dif s. $ <eq:deterministic>

The confluence relation applies for $s>=0$. At a reservoir source,
$tilde(R)_r (s)=n_(r t)+d_(r t)$ on $[tau_t,tau_(t+1))$.
For $s<0$, it is the independently supplied entry history of that reach.
Interval entry variables are the averages of these streams:

$ R_(r t) = frac(1, Delta_t) integral_(tau_t)^(tau_(t+1)) tilde(R)_r (s) dif s. $ <eq:eventmean>

Keeping the stream preserves the time of every release change. For example,
a one-hour release passing through two reaches with half-hour delays must
arrive as a one-hour release shifted by one hour. On an hourly mixing grid,
averaging after the first reach would spread some arrival before and after
that final interval.
Zero delays and delays shorter than a control interval allow newly released
water to affect downstream flow in that same interval. Entry capacities are
therefore also checked on the stream, including both sides of each jump.

If any reach has a distributed delay, @eq:confluence instead uses the mean
upstream arrivals as the next reach's interval release. Equations @eq:transfer
and @eq:arrival are exact for those rectangular cohorts, but this downstream
averaging approximates the continuous timing. A finer transport grid tests
its effect by comparing cumulative arrivals and water in transit. Conservation
alone does not reveal a pulse that has arrived at the wrong time.

== Water in transit, history and environmental flows

#block(breakable: false)[
Water released but not yet delivered belongs to reach inventory $W$.
Its change equals entry volume minus arrival volume:

$ W_(r,t+1) = W_(r t) + kappa Delta_t R_(r t) - A_(r t), quad
  W_(r 0) = W_r^"hist". $ <eq:inventory>
]

Pre-horizon entry records determine both initial inventory and future
historical arrivals. For a past cohort $k$ released on $[s_k,e_k)$,
let $F_k$ be its selected delay mixture. Water released at $s<0$ has age
$-s$ at the horizon start, so the fraction still in transit is $1-F_k (-s)$:

$ W_r^"hist" = kappa sum_(k in cal(H)_r) R_(r k)
  integral_(s_k)^(e_k) [1 - F_k (-s)] dif s. $ <eq:historyinventory>

The set $cal(H)_r$ contains the supplied historical cohorts of reach $r$;
unspecified earlier entry is zero. Each reach's history already includes its
past upstream contributions. Propagating those contributions again would
count them twice. Restarting from an executed schedule therefore retains the
original entry cohorts and their flows, along with reservoir states and
operating histories. A past average arrival alone does not identify the
release times needed to predict remaining arrivals.

Summing reservoir and reach balances cancels every internal transfer:

$ sum_i (V_(i,t+1) - V_(i t)) + sum_r (W_(r,t+1) - W_(r t)) \
  = kappa Delta_t (sum_i I_(i t) + sum_r n_(r t) + E_t). $ <eq:globalwater>

Here $E_t$ is net boundary exchange, positive into the modeled system.
For a plant outfall, the source reservoir loses $q$, the reach receives $q$,
and the receiving reservoir gains the later arrival $A$. These internal
movements cancel; only natural inflows and external exchanges change total
stored and in-transit water. The same identity includes reversed flow on
instantaneous connections.

Environmental requirements can constrain reach entry
$R_(r t)>=R_(r t)^min$, average arrival over a specified time window $J$,
or arrival at every instant:

$ integral_J f_r (s) dif s >= integral_J f_r^min (s) dif s, quad
  f_r (s) >= f_r^min (s). $ <eq:environment>

The supplied policy chooses the window or pointwise arrival constraint.
For rectangular cohorts, distributed arrival rates are piecewise linear;
fixed-delay streams are piecewise constant. Between combined response and
requirement breakpoints, the difference $f_r-f_r^min$ is affine, so its
minimum occurs at an endpoint. Endpoint checks thus enforce the pointwise
constraint; jumps require both one-sided limits. At the horizon endpoints,
use the limit from inside the horizon.

An aggregate flow requirement can also combine named unit discharges, river
entries, channel flows and a given inflow contribution $I_(o t)$:

$ sum_(g in cal(G)_o) q_(g t)
  + sum_(r in cal(D)_o union cal(C)_o) R_(r t) + I_(o t) >= M_(o t). $ <eq:observation>

The sets $cal(G)_o$, $cal(D)_o$ and $cal(C)_o$ select the contributions, and
$M_(o t)$ is the required minimum. For a channel, $R_(c t)=X_(c t)$ retains
its sign. This constraint observes flows already accounted for in the network;
it adds no water to any balance.

= Operating constraints over time

A unit's initial operating state and its duration in that state are given.
Binary starts $s^+$ and stops $s^-$ record changes between consecutive intervals:

$ u_(g t) - u_(g,t-1) = s_(g t)^+ - s_(g t)^-, quad
  s_(g t)^+ + s_(g t)^- <= 1, quad s_(g t)^(plus.minus) in {0, 1}. $ <eq:transitions>

#block(breakable: false)[
Minimum running time $L_g^"up"$ means that a start at $tau_t$ requires
$u_(g k)>=s_(g t)^+$ for each later interval with
#box[$0<=tau_k-tau_t<L_g^"up"$]. Minimum stopped time $L_g^"down"$ similarly
requires $u_(g k)<=1-s_(g t)^-$ whenever
#box[$0<=tau_k-tau_t<L_g^"down"$]. These tests use elapsed hours, including on
an irregular grid. Initial state age determines any remaining mandatory
running or stopped time.
]

The same constraints apply to $U$ when minimum plant operating times are
specified. Replacing one running unit with another does not restart an
already operating plant. An unfinished minimum duration at the horizon end
is carried into a restart rather than imposed through an artificial terminal state.

#block(breakable: false)[
Availability, forced states and prescribed output restrict the same decisions:

$ underline(u)_(g t) <= u_(g t) <= overline(u)_(g t), quad
  P_(g t) = P_(g t)^"given", quad q_(g t) = q_(g t)^"given". $ <eq:schedules>
]

Only supplied schedules are imposed. Plant state also satisfies
$underline(U)_(p t)<=U_(p t)<=overline(U)_(p t)$, and plant schedules constrain
the corresponding sums. Maintenance forces units off; zero and positive
output schedules require the appropriate off and on states. All restrictions,
including initial minimum durations, apply together. Inconsistent restrictions
make the model infeasible.

#block(breakable: false)[
For interval-average discharge or release $x_t$, supplied upward and downward
rates $r_x^+,r_x^-$ limit changes between interval midpoints:

$ -integral_(m_(t-1))^(m_t) r_x^- (s) dif s
  <= x_t - x_(t-1) <= integral_(m_(t-1))^(m_t) r_x^+ (s) dif s. $ <eq:ramps>
]

A running unit cannot enter its positive operating range from zero when
the allowed inter-interval increase is smaller than its minimum output.
Power ramps therefore
allow the minimum-output jump at startup and shutdown: add
$sum_(g in cal(G)_x) underline(P)_(g t) s_(g t)^+$ to the upward allowance
and $sum_(g in cal(G)_x) underline(P)_(g,t-1) s_(g t)^-$ to the downward
allowance. Here $cal(G)_x$ contains that unit or all units at the constrained
plant. Discharge and river-release
ramps have no jump allowance. The first interval uses a supplied previous
window's average and midpoint; if absent, only that comparison is omitted.

Reservoir volume and level ramps instead constrain vertex changes
$V_(i,t+1)-V_(i t)$ and $Z_i (V_(i,t+1))-Z_i (V_(i t))$, with rate integrals
over $[tau_t,tau_(t+1)]$. These bounds limit changes between recorded states,
not derivatives at every instant.

A tunnel interlock limits how many gates in a group $cal(E)_ell$ may be
open at once. With $o_(e t)=bold(1)[b_(e t)>0]$,

$ underline(N)_(ell t) <= sum_(e in cal(E)_ell) o_(e t) <= overline(N)_(ell t). $ <eq:interlock>

A positive position counts as open even when its flow is zero.

= Economic objective

Generation at price $pi_t$ per MWh earns $pi_t Delta_t P_(g t)$.
Let $c^+,c^-$ be unit transition costs, $c^b$ nonnegative gate-movement costs,
$c^r$ nonnegative release-shortfall costs, and $v_i,v_r$ the values assigned
to water in reservoirs and reaches, per $"Mm"^3$. We maximize

$ cal(V) = & sum_t pi_t Delta_t sum_g P_(g t)
  - sum_(g,t) (c_(g t)^+ s_(g t)^+ + c_(g t)^- s_(g t)^-) \
  & - sum_(e,t) c_(e t)^b abs(b_(e t) - b_(e,t-1))
  - sum_(r,t) c_(r t)^r xi_(r t) \
  & + sum_i v_i (V_(i T) - V_(i 0)) + sum_r v_r (W_(r T) - W_(r 0)). $ <eq:objective>

Starts, stops and gate movements are charged once per event, without a duration
multiplier. The first gate comparison uses its supplied initial position.
The last two terms value changes in reservoir storage and reach inventory.
Valuing water in transit allows a late release to retain downstream economic
value even when it has not yet arrived. These terminal values are input
assumptions about use beyond the horizon.

A positive release penalty permits a deficit below the entry minimum.
The deficit variable is measured as a volume:

$ xi_(r t) >= 0, quad xi_(r t) >= kappa Delta_t (R_(r t)^min - R_(r t)). $ <eq:shortfall>

A positive cost makes $xi$ equal the actual shortfall at optimum. Without
that penalty, the entry minimum is hard and $xi=0$. Arrival minima and
aggregate flow requirements remain hard. Subtracting initial water inventories
in @eq:objective removes a constant endowment from the reported value while
leaving the maximizing decisions unchanged.

= Representing curves and strengthening bounds
<sec:graphs>

Together, the preceding equations and inequalities define a mixed-integer
nonlinear program (MINLP). The head, efficiency and travel-response couplings
are nonconvex even when unit states and gate choices are fixed. We distinguish
two tasks: representing the supplied functions without changing them, and
constructing relaxations that bound the objective during global search.

== Representing supplied curves exactly

Direct evaluation of a curve is convenient for nonlinear dispatch. Global
search also needs its inputs and outputs exposed as algebraic variables.
The _graph_ of a function is simply the set of input-output pairs that obey
it. For a piecewise-linear curve with knots $(x_i,y_i)$, its graph is

$ lambda_i >= 0, quad sum_i lambda_i = 1, quad
  x = sum_i x_i lambda_i, quad y = sum_i y_i lambda_i, quad lambda in op("SOS2"). $ <eq:sos>

The SOS2 condition selects one segment~@huch. Its two neighboring weights
sum to one, and their weighted input coordinate fixes the interpolation
fraction. The weighted output therefore gives exactly the supplied linear
interpolant. Allowing nonadjacent weights would also allow points between
distant knots and would change the curve.

For a turbine table, discharge weights $lambda_i$ and head weights $mu_j$
each satisfy the corresponding SOS2 and coordinate constraints. Let $E_(i j)$
be efficiencies at the knots. Linear interpolation in head combines the
values $eta_j$ on the neighboring head curves. A shape-preserving cubic in
discharge~@fritsch is represented by adding corrections to the linear curve:

$ w_i = lambda_i lambda_(i+1), quad r_i = lambda_i w_i, \
  eta_j = sum_i E_(i j) lambda_i
    + sum_i [A_(i j) r_i + B_(i j) (w_i - r_i)], quad
  eta = sum_j mu_j eta_j. $ <eq:cubic>

For segment width $delta q_i=q_(i+1)-q_i$, secant slope
$d_(i j)=(E_(i+1,j)-E_(i j)) \/ delta q_i$ and original endpoint derivatives
$m_(i j)$, use

$ A_(i j) = delta q_i (m_(i j) - d_(i j)), quad
  B_(i j) = delta q_i (d_(i j) - m_(i+1,j)). $ <eq:hermitecoefficients>

On the active segment, the two correction terms are multiples of
$lambda_i^2 lambda_(i+1)$ and $lambda_i lambda_(i+1)^2$. Both vanish at the
endpoints; their coefficients make the endpoint derivatives equal the supplied
ones. The auxiliary products $w_i,r_i$ implement these terms using quadratic
equalities. Setting $A=B=0$ gives bilinear interpolation. Restricting a domain
does not refit the original cubic slopes. Outside the discharge table we use
the declared end secants; head extrapolation, when allowed, uses its outer
secant. An idle unit's network head may lie outside the head range of its
running table. A finite extension keeps the algebraic relation defined there;
$q=P=0$, and the running head and efficiency limits remain inactive.

Units sharing the same net-head variable and compatible table grids can reuse
head weights. A release cohort likewise uses one set of flow weights for its arrivals, environmental
checks and terminal inventory, ensuring one travel response throughout.
Channel responses use the analogous piecewise-bilinear graph. The product
equalities in @eq:cubic, @eq:power and @eq:arrival remain nonlinear: these
representations preserve the curves rather than turning the model into a MILP.

== Quadratic equations for tunnels and gates

A binary direction variable separates positive and negative tunnel flow:

$ Q = Q^+ - Q^-, quad 0 <= Q^+ <= C z, quad 0 <= Q^- <= C(1-z), \
  S = (Q^+)^2 - (Q^-)^2, quad z in {0, 1}. $ <eq:signedloss>

If $z=1$, only $Q^+$ can be positive; if $z=0$, only $Q^-$ can be positive.
Thus $S=Q abs(Q)$ in both cases. A directed outfall needs only $S=Q^2$;
known head ordering can also fix a direction.

For each finite gate position, a binary selector identifies whether that
position is chosen. Exactly one selector is active, so $b$, the open state
$o$ and the selected resistance $K(b)$ are weighted sums of known values.
To impose the loss equation only when the gate is open, define
$d_o=o(H_v-H_w)$. If $L<=H_v-H_w<=U$, this binary product is represented by

$ L o <= d_o <= U o, quad
  H_v - H_w - U(1-o) <= d_o <= H_v - H_w - L(1-o). $ <eq:binaryhull>

At $o=0$ these inequalities force $d_o=0$; at $o=1$ they force
$d_o=H_v-H_w$. Together with $abs(Q)<=C o$ and $d_o=K(b)S$, they recover
@eq:tunnel and @eq:closed. Any finite resistance may be assigned to the
closed selector because $S=0$ there. This construction needs only quadratic
product equalities. Similarly, a wet/dry binary branch represents
$ell=max(H-z^"crest",0)$, and a nonnegative auxiliary root $r$ with $r^2=ell$
represents the square root in the outlet law @eq:outlets.

== Upper bounds for global optimization

#block(breakable: false)[
A _relaxation_ allows additional combinations of variables while retaining
every feasible schedule. Its maximum objective is therefore an upper bound
on what the physical model can earn. Relaxed flow, head and efficiency may
combine to predict physically impossible power. A _valid inequality_ excludes
some of these combinations without excluding a feasible schedule; it tightens
the bound while leaving the physical equations unchanged.
]

Finite storage and equipment limits provide variable domains. Conservation
and inventory propagation can narrow them using the initial and terminal
restrictions. These bounds cover all admissible schedules, rather than a
neighborhood of one previously found solution.

For clarity, omit unit and time indices in the following bounds.
For a running unit with nonnegative $q,h,eta$, let $overline(e)$ bound its
electrical efficiency. Then $p^"ub" (q,h)=gamma overline(e) q h eta(q,h)$
is an upper bound on physical electrical power. It is polynomial within
each interpolation cell even when @eq:power is implicit in $P$.
Suppose a plane $a q+b h+c$ bounds this polynomial throughout the running
domain, with $b>=0$. For a head lower bound $h_*$ valid in both operating
states, the following inequality includes the idle state:

$ P <= a q + b(h - h_*) + (c + b h_*) u. $ <eq:support>

When the unit runs, it reduces to $P<=a q+b h+c$. When it is idle,
$P=q=0$ and the right side is nonnegative. To establish $c$, we bound the
residual $p^"ub"-a q-b h$ over every original interpolation cell.
In a Bernstein representation of a polynomial on a box, the basis weights
are nonnegative and sum to one, so its value cannot exceed its largest
coefficient. This supplies a bound over the entire cell, including interior
extrema that a plane fitted at a few samples might miss.

Another relaxation retains joint information about flow, head and efficiency.
Consider a box containing their admissible running values, with up to eight
corners $(q_c,h_c,eta_c)$. Exact binary-product inequalities, as in
@eq:binaryhull, define $h_"on"=u h$ and $eta_"on"=u eta$.
Nonnegative weights $z_c$ on the corners satisfy

$ sum_c z_c = u, quad sum_c q_c z_c = q, quad
  sum_c h_c z_c = h_"on", quad sum_c eta_c z_c = eta_"on", \
  gamma e_min sum_c q_c h_c eta_c z_c
  <= P <= gamma e_max sum_c q_c h_c eta_c z_c. $ <eq:box>

Here $e_min,e_max$ bound electrical efficiency over the power range.
For any point in the box, multilinear interpolation weights reproduce its
coordinates and the product $q h eta$. Hence the inequalities include every
physical running point. When $u=0$, all corner weights vanish. The original
efficiency and power equalities still determine the physical output.

For a tabulated turbine with head knots $h_j$, head weights $mu_j$ already
available in its graph can
also strengthen power bounds. We introduce $zeta_j$ to equal the head
weights while running and vanish while idle:

$ 0 <= zeta_j <= mu_j, quad sum_j zeta_j = u, quad
  sum_j (h_j - h_1) zeta_j = h_"on" - h_1 u, \
  P <= a q + sum_j c_j zeta_j. $ <eq:coordinatecuts>

The coefficients $c_j$ are chosen so that the weighted expression bounds the
polynomial upper bound $p^"ub"$, including its variation between head knots.
The constraints give $zeta=mu$ when $u=1$, and $zeta=0$ when $u=0$ without
fixing the idle hydraulic head. Discharge weights support analogous bounds
along the other table axis. Additional planes may be generated during the
search, provided each is checked over its complete domain. All these
inequalities strengthen the relaxation while retaining the original equations.

= Solving the model and checking a schedule

== Finding feasible controls

JuMP provides the algebraic interface in Julia~@jump. SCIP is used for both
linear mixed-integer proposals and global nonlinear search~@scip; Ipopt
solves continuous nonlinear dispatch with discrete choices fixed~@ipopt.
The preparatory step uses a MILP approximation to propose unit states and
finite gate positions. Nonlinear dispatch then adjusts continuous decisions
under the original physical laws for those choices. Several proposals can
be tried within a shared time allowance. Candidates passing the independent
checks below are submitted as starts for the full SCIP model, where the
discrete decisions are free again.

This preparation seeks feasible schedules, not bounds on the best possible
objective. Its linear objective is a surrogate for @eq:objective.
Dispatch with these choices fixed remains a piecewise-smooth, nonconvex nonlinear
program, so convergence to a stationary point does not establish global
optimality.

== Reconstructing the water trajectory

An independent integration holds the selected unit discharges and openings
fixed and recovers the resulting water trajectory. At each interval, Newton's
method solves jointly for ending storage, junction heads, tunnel flows and
active outlet releases. With a delay shorter than the interval, some current
release arrives before the interval ends; that arrival must participate in
the same solve rather than being computed only from past releases. Electrical
output is then recovered from @eq:power.

Some isolated pressure components determine only head differences. Adding
the same constant to all their heads leaves the flow equations unchanged.
One numerical reference removes this indeterminacy during solution, while
every water balance is retained. The recovered heads can be shifted together
to respect the junction bounds when such a shift is admissible.

Acceptance first requires the reconstructed schedule to satisfy the original
interval equations and operating restrictions. Integration is then repeated
on a finer grid with the same controls. Schedules and ramps are still checked
over their original operating windows; storage at the original vertices,
total water, cumulative arrivals and terminal inventory are compared separately.
For supplied volume and relative objective tolerances
$epsilon_V,epsilon_(cal(V))$, require also

$ max_(i,t) abs(V_i^"fine" (tau_t) - V_(i t)) <= epsilon_V, quad
  max_r abs(W_r^"fine" (tau_T) - W_(r T)) <= epsilon_V, \
  abs(cal(V)^"fine" - cal(V)) <= epsilon_(cal(V)) max(1, abs(cal(V))). $ <eq:replaycriteria>

This refinement tests numerical consistency with the original controls.
It estimates the effect of the time grid rather than rigorously bounding
continuous-time hydraulic error.

== Bounding the remaining objective improvement

A feasible schedule gives an achievable objective value $L$. Global search
provides an upper bound $U$ that no feasible schedule in the declared model
can exceed. SCIP obtains such bounds from relaxations, tightens them and
branches on discrete choices or continuous variable ranges. A region can
be discarded when its bound cannot improve the current feasible value.
A good schedule therefore helps the search even before it is proved optimal.

Here $L$ is the independently reconstructed value of a schedule feasible for
the interval model, and $U$ must be a finite enclosing bound. The absolute
and relative limits on possible improvement are

$ g_"abs" = max(0, U-L), quad g_"rel" = frac(max(0, U-L), max(1, abs(L))). $ <eq:gap>

For example, $L=100$ and $U=103$ in the same monetary units limit improvement
to at most $3$, or $3\%$ under this relative-gap definition. They do not imply
that a better schedule exists. Subtracting the initial water endowment in
@eq:objective leaves maximizing decisions unchanged but affects a relative
gap's denominator; both $L$ and $U$ use that same objective convention.

A certificate records that the requested absolute or relative gap has been
met for the chosen interval model, interpolation rules and finite gate sets.
If discrete choices are fixed, it is conditional on those choices.
The finer-grid check answers a separate question about the physical trajectory:
a schedule can have a certified interval-model gap and still fail finer
integration. Neither result certifies a continuous-time optimum.
These are floating-point numerical certificates, interpreted at solver
tolerances. They do not include a rigorous enclosure of every rounding error
in the solver and polynomial bound calculations.

#block(breakable: false)[
= Discussion

The formulation connects storage, pressure loss, generation and travel time
within one scheduling problem. Head and efficiency change as decisions
change, rather than being fixed beforehand in a production approximation.
River inventory gives released water an explicit place in the balance until
it arrives. Supplied curves remain the functions chosen by the modeler,
including their stated interpolation and extension rules.

The solution method uses this same physical model in distinct ways. Nonlinear
dispatch finds feasible controls; independent integration checks what those
controls do; global bounds limit the economic improvement left unresolved.
Stronger relaxations help the global search without replacing hydraulic laws.
These choices preserve the model's couplings, but they also preserve its
nonconvexity. A large network may still require substantial time to close
the optimality gap.

The numerical model has a stated time resolution, river response approximation
and finite set of tunnel positions. Calibration, transient hydraulic regimes
and uncertain future inputs remain separate modeling questions. Extensions
can add explicit physical laws and operating constraints while using the
same water accounting and distinguishing feasible operation from objective
bounds and time-grid consistency.
]

#block(breakable: false)[
  #bibliography("references.bib", style: "ieee", title: "References")
]
