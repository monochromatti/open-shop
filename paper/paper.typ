#import "style.typ": article, title-block, abstract-block, accent, muted, rule
#show: article

#title-block[
  Nonlinear hydropower scheduling \
  with coupled hydraulics \
  and delayed river networks
]
#abstract-block[
  Short-term hydropower scheduling must decide when to generate or pump and how much
  water to release. The water available to downstream plants depends on
  earlier releases and river travel times. Reservoir levels and shared tunnel
  losses also affect the power obtained from each unit. We formulate the
  schedule as a mixed-integer nonlinear program with unit commitment, finite
  gate choices and delayed river transport. Reservoir balances use an implicit
  discrete-gradient step that respects stored potential energy. Pressure losses
  and power conversion are nonlinear, and
  tabulated equipment curves keep their specified interpolants. For fixed
  river travel times, releases are shifted through the network without
  averaging at confluences. Distributed delays use interval cohorts and track
  water still in transit. SCIP searches the full model using feasible starts
  prepared by mixed-integer proposals and nonlinear dispatch. A separate
  integration checks the controls on the optimization grid and a finer grid.
  Optional spinning capacity uses feasible unit operating points at scheduled
  hydraulic node heads, with shared tailwater and qualified response limits.
  We state the limits of the resulting objective-gap certificate. The
  implementation uses Julia and open-source solvers.
]

= The scheduling problem

A reservoir can generate now or save water for later. A release lowers its
level and changes the energy available from the next release. In a cascade,
the water may pass through several plants before leaving the system. The
upstream schedule then determines when water becomes available downstream.
Shared tunnels impose another dependency, since a plant's discharge affects
pressure losses for other plants on the waterway. A unit may be stopped,
but when running it must stay within its operating range. These discrete
choices and the dependence of generation on head have long been studied
together~@catalao.

We take the network, initial volumes, operating histories, prices, inflows,
equipment curves and operating limits as given. The optimizer chooses unit
states, continuous unit discharges, outlet openings and optional spinning
reserve allocations. It selects each
adjustable tunnel position from a finite set. The associated volumes, flows,
heads, river arrivals and electrical outputs must satisfy the equations
below. We maximize net electrical revenue
less operating penalties, with a terminal value for the water left in the
system.

We assume deterministic inputs and quasi-steady tunnel hydraulics: pressures
and flows are in equilibrium within each scheduling interval while reservoir
storage evolves over time. Tunnel inertia and water hammer are excluded.
Discharge and travel-time response curves describe the rivers, without
solving river momentum equations during scheduling. We do
not include stochastic decisions. Pumps use supplied hydraulic
efficiency functions and constant motor efficiency; generating and pumping
modes retain separate positive flow and power magnitudes.
The reserve capacity model
below does not simulate hydraulic activation across the watercourse.

The optimization uses nonlinear pressure-loss and power equations. Linear
approximations provide initial discrete choices for the solution procedure
in @sec:solution.

= Water balance and hydraulic head
<sec:water>

Let $0=tau_0<dots.h.c<tau_T$ be time edges in hours,
$Delta_t=tau_(t+1)-tau_t$ the interval duration, and
$m_t=(tau_t+tau_(t+1)) \/ 2$ its midpoint, for $t=0,dots,T-1$.
Interval inputs are piecewise constant, with changes at these edges.
We record volumes at the edges and evaluate interval hydraulics using
volume-averaged reservoir heads.

Hydraulic head is mechanical energy per unit weight, expressed in metres.
At a reservoir it is the water-surface elevation. At a flooded junction it
includes elevation and pressure divided by $rho g$, where $rho$ is water
density and $g$ is gravitational acceleration. Water flows down a head
gradient. Friction dissipates energy along the waterway, and a turbine can
extract power from the flow. The equipment relationships account for
velocity corrections and local losses.

#block(breakable: false)[
Volume is measured in $"Mm"^3$ (one million cubic metres, $10^6 "m"^3$),
discharge in $"m"^3/"s"$ and time in hours. The conversion factors are

$ kappa = frac(3600, 10^6) = 0.0036, quad
  gamma = frac(rho g, 10^6) = 0.00981 \
  (rho = 1000 "kg"/"m"^3, quad g = 9.81 "m"/"s"^2). $ <eq:units>
]

Thus $kappa Delta q$ converts a discharge over an interval into a volume,
and $gamma q h$ is hydraulic power at discharge $q$ and net head $h$.
For example, a release of $10 "m"^3/"s"$ for one hour removes
$0.036 "Mm"^3$ of water. We use the same conversion for net inflow in the reservoir balance.

Let $cal(R)$, $cal(J)$ and
$cal(B)$ index reservoirs, pressure junctions and prescribed-head boundaries.
A node $n$ is any of these three. Tunnels, instantaneous channels, generating and pumping
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
    [$H_(n t), h_(p t)$], [Node head, effective generating head or pumping lift], [m],
    [$Q_(e t), X_(c t), q_(g t)$], [Tunnel, channel, unit discharge], [$"m"^3/"s"$],
    [$R_(r t), d_(r t)$], [Reach entry flow, reservoir withdrawal], [$"m"^3/"s"$],
    [$A_(r t)$], [Volume arriving from a delayed reach], [$"Mm"^3$],
    [$P_(g t)$], [Positive production or consumption magnitude], [MW],
    [$r_(g k t)^(plus.minus), sigma_(k t)^(plus.minus)$], [Spinning allocation, group shortfall], [MW],
    [$u_(g t), U_(p t)$], [Unit and plant operating states], [${0, 1}$],
    [$b_(e t), a_(r t)$], [Tunnel position, reservoir outlet opening], [$[0, 1]$],
    table.hline(stroke: 0.7pt),
  ),
  kind: table,
  caption: [Variables and units. Storage and river inventory are recorded at time
    vertices. The other indexed quantities refer to control intervals.],
) <tab:variables>

Let $F_(n t)$ be net flow into node $n$ from the modeled connections,
excluding natural reservoir inflow $I_(i t)$. The matrices $B^E$, $B^C$ and
$B^G$ record the tunnel, channel and unit connections. For a connection
between nodes, its column is $-1$ at the source, $+1$ at the destination,
and zero elsewhere. Negative tunnel/channel discharge reverses these contributions.
Unit endpoints always follow their actual water direction, so a pump withdrawing
from the lower reservoir and delivering to the upper one uses nonnegative $q$.
Arrival volumes from rivers and withdrawals into them give the remaining terms:

$ F_(n t) = sum_e B_(n e)^E Q_(e t)
  + sum_c B_(n c)^C X_(c t) + sum_g B_(n g)^G q_(g t)
  + sum_r D_(n r) frac(A_(r t), kappa Delta_t)
  - sum_r S_(n r) d_(r t). $ <eq:incidence>

$D_(n r)=1$ if reach $r$ delivers into node $n$, and $S_(n r)=1$ if an outlet
at $n$ withdraws water into that reach. The other entries are zero.
Dividing $A_(r t)$ by $kappa Delta_t$ converts arrival volume to flow rate.
For a tunnel or plant outfall, $F$ includes the withdrawal from the source
node. The outfall enters the river balance below. Its outlet head can depend
on a downstream reservoir even though the water reaches that reservoir
later (@fig:network).

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
exchange water with the exterior. Let $Z_(i)(V)$ be reservoir surface elevation
and $Psi_(i)$ an antiderivative of $Z_(i)$. The interval head is
the volume average of that elevation:

$ overline(Z)_(i) (a,b)=cases(
  (Psi_(i) (b)-Psi_(i) (a))/(b-a) & "if" a != b,
  Z_(i) (a) & "if" a=b), quad
  H_(i t)=overline(Z)_(i) (V_(i t),V_(i,t+1)). $ <eq:heads>

$ underline(H)_j <= H_(j t) <= overline(H)_j, quad
  H_(b t) = H_b^"given". $ <eq:boundaryheads>

$Z_(i)$ converts volume to surface elevation. It is an increasing supplied
table, or $Z_(i) (V)=z_i+s_i V+c_i V^2$ with given coefficients.
The ending volume is unknown in @eq:storage and also enters @eq:heads.
We solve these coupled equations implicitly. The head is a discrete gradient
of stored potential~@celledoni: multiplying it by the interval volume change
gives exactly $Psi_(i) (b)-Psi_(i) (a)$. In MWh, the gravitational energy change is
$(0.00981/kappa)(Psi_(i) (b)-Psi_(i) (a))$. A common elevation datum cancels from
the balance of a closed water network.

For a quadratic level curve, the formula requires no division:

$ overline(Z)_(i) (a,b)=z_i+s_i (a+b)/2
  +c_i ((a+b)^2/4+(b-a)^2/12). $ <eq:polynomialhead>

For a tabulated curve, we integrate its original linear segments, splitting
only at crossed knots. An endpoint trapezoid would be wrong if an intervening
knot changes the slope. The implementation also avoids subtracting nearby
potentials when $a$ and $b$ almost coincide. Affine curves retain the usual
midpoint head. With curvature, this head differs from $Z_(i) ((a+b)/2)$ by
$c_(i) (b-a)^2/12$.

#figure(
  block(width: 144mm, height: 40mm)[
    #let node(x, body) = place(top + left, dx: x, dy: 23mm,
      box(width: 24mm, height: 14mm, stroke: 0.5pt + rule,
        fill: rgb("#f3f6f7"), radius: 1pt,
        align(center + horizon, text(size: 9pt, body))))
    #let arrow(x) = {
      place(top + left, dx: x, dy: 30mm, line(length: 6mm, stroke: 0.6pt + accent))
      place(top + left, dx: x + 4.5mm, dy: 29.2mm,
        line(end: (1.5mm, 0.8mm), stroke: 0.6pt + accent))
      place(top + left, dx: x + 4.5mm, dy: 30.8mm,
        line(end: (1.5mm, -0.8mm), stroke: 0.6pt + accent))
    }
    #node(0mm, [Reservoir A \ $V_A, H_A$])
    #node(30mm, [Pressure \ junction $H_J$])
    #node(60mm, [Plant \ $q_1 + q_2$])
    #node(90mm, [Reach \ inventory $W$])
    #node(120mm, [Reservoir B \ $V_B, H_B$])
    #for x in (24mm, 54mm, 84mm, 114mm) { arrow(x) }
    #place(top + left, dx: 17mm, dy: 17mm, text(size: 8pt, fill: muted)[tunnel $Q$])
    #place(top + left, dx: 81mm, dy: 17mm, text(size: 8pt, fill: muted)[release])
    #place(top + left, dx: 112mm, dy: 17mm, text(size: 8pt, fill: muted)[arrival])
    #place(top + left, line(start: (132mm, 23mm), end: (132mm, 8mm),
      stroke: (thickness: 0.5pt, paint: muted, dash: "dashed")))
    #place(top + left, line(start: (132mm, 8mm), end: (72mm, 8mm),
      stroke: (thickness: 0.5pt, paint: muted, dash: "dashed")))
    #place(top + left, line(start: (72mm, 8mm), end: (72mm, 23mm),
      stroke: (thickness: 0.5pt, paint: muted, dash: "dashed")))
    #place(top + left, line(start: (71.2mm, 21.5mm), end: (72mm, 23mm), stroke: 0.5pt + muted))
    #place(top + left, line(start: (72.8mm, 21.5mm), end: (72mm, 23mm), stroke: 0.5pt + muted))
    #place(top + left, dx: 78mm, dy: 2mm, text(size: 8pt, fill: muted)[outlet head reference])
  ],
  caption: [Solid arrows show where the water goes. The dashed connection specifies
    the plant's outlet head. Water reaches reservoir B after traveling through
    the river, but B's present level already affects generation.],
) <fig:network>

= Tunnel hydraulics and energy conversion

== Tunnels and gates

For an oriented tunnel $e=(v,w)$, positive $Q_e$ flows from $v$ to $w$.
With opening $b>0$, the head difference equals friction loss:

$ H_(v t) - H_(w t) = K_e (b_(e t)) Q_(e t) abs(Q_(e t)), \
  K_e (b) = cases(
    R_e / b & "without a gate-loss curve",
    R_e + G_e (b) & "with a supplied gate-loss curve".
  ) $ <eq:tunnel>

The tunnel resistance $R_e>0$ and additional gate resistance $G_e>=0$ have
units $"s"^2/"m"^5$. Positive $Q$ requires $H_v>H_w$, and negative $Q$
requires $H_v<H_w$. In either direction, friction dissipates energy with a
loss proportional to $Q^2$. Flows through branches and loops follow from
these head equations and the junction balances.

A closed gate permits a head difference while stopping flow:

$ b_(e t) = 0 ==> Q_(e t) = 0, quad
  abs(Q_(e t)) <= C_(e t). $ <eq:closed>

Equation @eq:tunnel applies when the gate is open. Capacity $C_(e t)$ bounds
flow in either direction. River outfalls that allow only discharge have
$Q_(e t)>=0$. A tunnel position may be prescribed or chosen from a finite
set $cal(S)_e subset [0,1]$. Prescribed positions may take any value in
$[0,1]$, but the optimizer searches only the listed positions.

#block(breakable: false)[
An instantaneous channel between reservoirs follows a supplied discharge law:

$ X_(c t) = Phi_c (H_(v t) - H_(w t), H_("ref",t)), quad
  abs(X_(c t)) <= C_(c t). $ <eq:channel>
]

$Phi_c$ is a family of curves indexed by a specified upstream or downstream
reservoir level. Each curve gives discharge as a function of level difference,
is nondecreasing and passes through $(0,0)$. Positive and negative responses
need not be symmetric. We interpolate linearly along each curve and between
reference levels. We extend the end segments beyond the tabulated level differences. Outside
the reference-level range, we use the nearest curve. Water passes through
the channel within the interval, without a separate transit inventory.

== Energy conversion

Write $q_(p t)=sum_(g in cal(G)_p) q_(g t)$ and
$P_(p t)=sum_(g in cal(G)_p) P_(g t)$ for plant totals. All units at a plant
share the same effective head. Each plant has one operating mode, with
$sigma_p=1$ for generation and $sigma_p=-1$ for pumping:

$ h_(p t) = sigma_p (H_("source"(p),t)
  - max{H_("target"(p),t), h_p^"floor"}
  - T_p (q_(p t))). $ <eq:nethead>

The target head is the outlet reference in @fig:network. If an outlet-level
floor is specified, the maximum in @eq:nethead enforces it. Otherwise we use
the target head alone. Total plant discharge raises the tailwater by
$T_p (q_(p t))>=0$. We model shared penstocks as tunnels, whose friction
loss depends on the combined discharge of the units they supply. For pumping,
@eq:nethead is destination head minus intake head plus the path loss.

Distinct penstocks can feed units that share a pressurized tailrace. Give each
penstock its own intake node and connect the plant outlets to one junction.
Conservation makes the tailrace tunnel carry their combined discharge.
Its quadratic loss raises the common outlet head, reducing every unit's head
once. These hydraulic supply groups remain separate plants for operating
rules; a common outlet does not combine their power caps or commitment rules.

#block(breakable: false)[
The binary state $u_(g t)$ distinguishes a stopped unit from a committed unit:

$ u_(g t) in {0, 1}, quad
  underline(q)_(g t) u_(g t) <= q_(g t) <= overline(q)_(g t) u_(g t), quad
  underline(P)_(g t) u_(g t) <= P_(g t) <= overline(P)_(g t) u_(g t). $ <eq:onoff>
]

Lower flow and power limits may be zero. Commitment still records startup
and dwell obligations when a committed unit produces no power.

The electrical output is the available hydraulic power multiplied by turbine
and electrical efficiencies:

$ P_(g t) = gamma q_(g t) h_(p t)
  eta_g (q_(g t), h_(p t)) e_g (P_(g t)). $ <eq:power>

For a pump, the hydraulic work is supplied by consumed electrical power.
With constant motor efficiency $e_g^"motor" in (0,1]$, the equation is

$ P_(g t) eta_g(q_(g t),h_(p t)) e_g^"motor"
  = gamma q_(g t) h_(p t). $ <eq:pump>

Lower efficiency therefore increases consumption. The same on/off and
operating-domain constraints apply, with positive pumping lift required.
Hydraulic efficiency tables retain their declared interpolation. Electrical
motor-efficiency curves are outside this formulation. For a fixed-speed pump,
one supplied discharge and hydraulic efficiency at each reference head define
piecewise-linear functions $Q_g(h)$ and $eta_g(h)$. Its equations are

$ q_(g t)=u_(g t) Q_g(h_(p t)), quad
  P_(g t) eta_g(h_(p t)) e_g^"motor"=gamma q_(g t) h_(p t). $

No variable-discharge efficiency surface is inferred. Head coverage is required
unless linear end-secant extrapolation is explicitly supplied. A single
reference head requires an explicit constant characteristic. Static and
operative unit bounds intersect this point, and an incompatible bound makes
operation infeasible. Independent integration retains commitment and solves
characteristic flow jointly with storage and pressure, including on finer
replay grids. Schedule validation checks the original flow against the
characteristic without replacing it. Equal discharge bounds on an ordinary
variable-discharge unit instead prescribe a fixed-flow operating point.

#block(breakable: false)[
When $u_(g t)=1$, the unit must also satisfy head and efficiency limits and,
when supplied by a turbine table, head-dependent discharge limits:

$ underline(h)_g <= h_(p t) <= overline(h)_g, quad
  eta_g^min <= eta_g (q_(g t), h_(p t)) <= 1, \
  q_g^"lo" (h_(p t)) <= q_(g t) <= q_g^"hi" (h_(p t)). $ <eq:unitdomain>
]

When $u_(g t)=0$, @eq:onoff gives $q=P=0$ and the operating inequalities in
@eq:unitdomain are inactive. Plant head remains tied to the hydraulic node
heads when all units are off. Isolated pressure components may retain a common
reference freedom, as discussed below.

Turbine efficiency is a supplied table or the analytic relationship

$ eta_g (q,h) = eta_g^0
  - alpha_g (frac(q - q_g^*, q_g^*))^2
  - beta_g (frac(h - h_g^*, h_g^*))^2. $ <eq:analyticeta>

Here $q_g^*$ and $h_g^*$ are reference operating conditions, $eta_g^0$ is
their turbine efficiency. Nonnegative $alpha_g,beta_g$ describe a decrease
away from them; signed coefficients are also admitted, subject to the same
operating efficiency bounds. At reference head $h_j$, a supplied efficiency curve
$f_j(q)$ has its own discharge knots. These grids may differ between heads.
Interpolate each original curve linearly or with its original shape-preserving
cubic slopes, then interpolate its value at the same discharge between heads:

$ eta_g(q,h) = (1-b) f_j(q) + b f_(j+1)(q), quad
  b = (h-h_j)/(h_(j+1)-h_j). $ <eq:headwiseefficiency>

A strict lookup between heads requires both contributing curves to cover $q$.
At an exact reference head, only that curve contributes. For complete-case
strict coverage, every column participating over the declared operating head
range must cover the unit's static discharge range. When discharge continuation
is explicitly allowed, each curve uses its own outer secant before head
interpolation. Efficiency coverage does not define or modify the operating
discharge envelope.
A single reference head requires an explicit assumption that efficiency and
the discharge envelope are independent of head. The actual hydraulic head
still enters @eq:power. Operating limits are supplied separately from
efficiency knots, with explicit linear continuation required outside the data.

Turbine efficiency converts hydraulic power to shaft power. Electrical
efficiency $e_g$ accounts for the remaining conversion to electrical output.
It is a separate piecewise-linear function of _electrical output_ $P$, so
@eq:power must be solved for $P$ even when $q,h$ are known. We require
$P \/ e_g (P)$ to be strictly increasing over the admissible output range.
The scalar power equation then has at most one admissible solution.
If no electrical-efficiency curve is given, $e_g=1$.

Plant state $U_(p t)$ is one whenever at least one unit runs:

$ U_(p t) >= u_(g t) quad (g in cal(G)_p), quad
  U_(p t) <= sum_(g in cal(G)_p) u_(g t), quad U_(p t) in {0, 1}. $ <eq:plantstate>

Plant limits apply to aggregate output and discharge:
$P_(p t)<=overline(P)_(p t)$ and $q_(p t)<=overline(q)_(p t)$ when supplied,
and $P_(p t)>=underline(P)_(p t) U_(p t)$,
$q_(p t)>=underline(q)_(p t) U_(p t)$.

A single reversible machine has generating and pumping mode records on
oppositely connected plants. An exclusion consists of supplied groups of unit
names. For any two different groups $K_a,K_b$, we require

$ u_(i t)+u_(j t)<=1 quad (i in K_a, j in K_b, a!=b). $ <eq:modeexclusion>

Parallel units within one group may run together. Initially running states
must satisfy the same exclusion. Existing unit dwell histories, power and
discharge histories carry into a restarted horizon without a separate pump
state. Independently connected pumps and turbines may run simultaneously when
no exclusion is supplied; shared ternary-machine characteristics are not modeled.

= Outlets and river transport

== Entry laws, confluences and splits

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
The reservoir loses $d=R-n$. Natural reach inflow enters the river directly
and must not be subtracted from reservoir storage.

#block(breakable: false)[
The square-root derivative is singular at a dry outlet. Local dispatch instead
uses $w>=0$ with $w^2=[ell]_+$. This piecewise polynomial equality gives
exactly $w=sqrt([ell]_+)$, so withdrawal is $c a w$ or $c w^3$.
The global model uses a binary positive-part graph followed by $w^2=[ell]_+$.
Both formulations retain the crest boundary without smoothing the physical law.
]
The analytic outlet laws are one-level responses; they do not calculate
downstream backwater or gate geometry. A supplied one-level curve is valid only
where the omitted conditions are fixed or accounted for in its calibration;
it cannot describe independently varying downstream counterpressure. Directed
river connections must be physically applicable over the operating range.
Transport conservation does not enforce downhill flow or a gravitational
energy balance along a routed reach. A plant or tunnel's hydraulic outlet
reference must also represent its actual outfall head. Tunnel connections
assume submerged entrances and exits.

#block(breakable: false)[
At a zero-storage river node, incoming water is distributed among outgoing
reaches with prescribed shares $alpha_r$. Local natural inflow and plant or
tunnel outfalls enter each branch after this allocation. Let $L_(r t)$ collect
these local injections:

$ L_(r t) = n_(r t) + sum_(g -> r) q_(g t) + sum_(e -> r) Q_(e t), quad
  R_(r t) = L_(r t) + alpha_r sum_(s -> r) frac(A_(s t), kappa Delta_t). $ <eq:confluence>
]

Here $g->r$ and $e->r$ identify local outfalls, while $s->r$ identifies
reaches arriving at the source node of $r$. At every zero-storage node,
$0<=alpha_r<=1$ and the outgoing shares sum to one. These supplied constants
hold throughout the horizon. A confluence with one outgoing reach has
$alpha_r=1$. A 40% branch receiving 5 m³/s upstream and 1 m³/s locally
carries 3 m³/s. Its capacity must accommodate both sources; spare capacity
on another branch cannot change the prescribed split. Delayed reaches form
a directed acyclic network, so we propagate arrivals in upstream order.
Tunnels and instantaneous channels may form loops.

== Release cohorts and travel-time responses

A _cohort_ is the water entering a reach during one release interval.
With a fixed travel time, that release simply appears downstream later.
With a distributed delay, different fractions of the cohort arrive at
different times. The response functions are deterministic approximations
of river transport~@dooge.

For reference entry flow $rho_(r j)$, let fraction $w_(r j ell)$ have a
travel time uniformly distributed over
$[delta_(r j ell),delta_(r j,ell+1)]$. Delay edges must be nonnegative and
strictly increasing. The fractions obey
$w_(r j ell)>=0$ and $sum_ell w_(r j ell)=1$.
The cumulative fraction that has arrived after elapsed time $s$ is

$ F_(r j) (s) = sum_ell w_(r j ell)
  min{1, max{0, frac(s - delta_(r j ell), delta_(r j,ell+1) - delta_(r j ell))}}. $ <eq:cdf>

Entry flow may change the travel response. For two or more nonnegative, strictly increasing reference
flows covering zero through reach capacity, interpolate between the neighboring
distributions. A missing zero-flow reference is rejected. Zero flow releases
no water. In the following equations, SOS2 is a _special ordered set of
type 2_: at most two weights may be nonzero, and their indices must be adjacent.

$ theta_(r k j) >= 0, quad sum_j theta_(r k j) = 1, quad
  R_(r k) = sum_j rho_(r j) theta_(r k j), quad
  theta_(r k) in op("SOS2"). $ <eq:delayweights>

The weights reproduce the cohort's total entry flow, including natural inflow
and upstream arrivals. For example, halfway between two reference flows,
the two distributions receive equal weights. A single distribution needs
only $theta_(r k 1)=1$. We then omit the flow-reproduction equation, so the
reference flow does not fix the release.
Once released, each cohort keeps its chosen mixture. Reference flows are
calibration data and are unaffected by changes in operating capacity.

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

The historical terms account for releases before $tau_0$. We evaluate the
integrals exactly for the uniform bins. Since $theta$ depends linearly on
$R$ within a reference-flow interval, multiplying the arrival fraction by
released volume gives a quadratic dependence on $R$. Interpolating arrival
volumes directly would give a different travel response.

Each cohort selects $theta$ from its own contemporaneous total release,
including historical cohorts at their original flows. The nonlinear dispatch,
global model and physical reconstruction use this same law. A frozen-reference
MILP proposes commitments but provides no bound for this nonlinear problem.
Physical releases stay within the supplied reference domain;
endpoint holding in Newton trials and outer-polynomial extensions in nonlinear
operators are numerical continuations only.

#block(breakable: false)[
When all reaches have fixed travel times $delta_r>=0$, we propagate the
piecewise-constant entry stream $tilde(R)_r (s)$ through the network. Each
reach shifts its stream by its delay. At a river node, each outgoing stream
receives its share of the shifted upstream streams:

$ f_r (s) = tilde(R)_r (s - delta_r), quad
  tilde(R)_r (s) = L_r (s) + alpha_r sum_(s' -> r) f_(s') (s), \
  A_(r t) = kappa integral_(tau_t)^(tau_(t+1)) f_r (s) dif s. $ <eq:deterministic>
]

The allocation relation applies for $s>=0$. At a reservoir source,
$tilde(R)_r (s)=n_(r t)+d_(r t)$ on $[tau_t,tau_(t+1))$.
For $s<0$, it is the independently supplied entry history of that reach,
already allocated to its branch. Upstream historical water arriving at a node
after zero is allocated once, on entry to the outgoing reaches.
Interval entry variables are the averages of these streams:

$ R_(r t) = frac(1, Delta_t) integral_(tau_t)^(tau_(t+1)) tilde(R)_r (s) dif s. $ <eq:eventmean>

For example,
a one-hour release passing through two reaches with half-hour delays must
arrive as a one-hour release shifted by one hour. On an hourly mixing grid,
averaging after the first reach would spread some arrival before and after
that final interval.
A release can reach downstream within its control interval when the delay
is zero or shorter than that interval. We check entry capacities on the
stream, including both sides of each jump.

If any reach has a distributed delay, @eq:confluence instead uses the mean
upstream arrivals as the next reach's interval release. Equations @eq:transfer
and @eq:arrival are exact for those rectangular cohorts, but this downstream
averaging approximates the continuous timing. We assess this error on a
finer transport grid by comparing cumulative
arrivals and water in transit. A volume balance cannot detect water that
arrives at the wrong time.

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

The historical cohorts of reach $r$ form $cal(H)_r$. We assume zero entry
before the supplied history. A reach's records already include past upstream
contributions, so propagating them again would count the water twice. To
restart after executing part of a schedule, we save its entry cohorts and
flows as well as reservoir states and operating histories. Average arrivals
alone cannot determine when the water still in transit will arrive.

Summing reservoir and reach balances cancels every internal transfer:

$ sum_i (V_(i,t+1) - V_(i t)) + sum_r (W_(r,t+1) - W_(r t)) \
  = kappa Delta_t (sum_i I_(i t) + sum_r n_(r t) + E_t). $ <eq:globalwater>

Here $E_t$ is net boundary exchange, positive into the modeled system.
For a plant outfall, the source reservoir loses $q$ and the reach receives
it. Arrival $A$ later leaves the reach and enters the receiving reservoir.
Both transfers cancel in the sum. Total stored and in-transit water changes
only through natural inflows and external exchanges, including when an
instantaneous connection carries reversed flow.

#block(breakable: false)[
Environmental requirements can constrain reach entry
$R_(r t)>=R_(r t)^min$, average arrival over a specified time window $J$,
or arrival at every instant:

$ integral_J f_r (s) dif s >= integral_J f_r^min (s) dif s, quad
  f_r (s) >= f_r^min (s). $ <eq:environment>
]

The supplied policy chooses the window or pointwise arrival constraint.
Distributed arrival rates from rectangular cohorts are piecewise linear.
Fixed-delay streams are piecewise constant. Between combined response and
requirement breakpoints, the difference $f_r-f_r^min$ is affine, so its
minimum occurs at an endpoint. It is enough to check these endpoints, taking both one-sided limits at a
jump. At the horizon endpoints, use the limit from inside the horizon.

An aggregate flow requirement can also combine named unit discharges, river
entries, channel flows and a given inflow contribution $I_(o t)$:

$ sum_(g in cal(G)_o) q_(g t)
  + sum_(r in cal(D)_o union cal(C)_o) R_(r t) + I_(o t) >= M_(o t). $ <eq:observation>

The sets $cal(G)_o$, $cal(D)_o$ and $cal(C)_o$ select the contributions, and
$M_(o t)$ is the required minimum. Channel flow $R_(c t)=X_(c t)$ is signed.
The inequality constrains existing
flows and adds no water to the balances.

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

Apply these minimum-duration constraints to $U$ if they are also specified
for the plant. A handover between running units does not restart the plant.
Any minimum duration unfinished at the horizon end becomes a requirement
on the next schedule.

#block(breakable: false)[
Availability, forced states and prescribed output restrict the same decisions:

$ underline(u)_(g t) <= u_(g t) <= overline(u)_(g t), quad
  P_(g t) = P_(g t)^"given", quad q_(g t) = q_(g t)^"given". $ <eq:schedules>
]

We impose an output schedule only where one is supplied. Plant schedules
constrain the sums of unit outputs, and plant state obeys
$underline(U)_(p t)<=U_(p t)<=overline(U)_(p t)$. Maintenance requires units
to be off. Zero scheduled output requires the off state, and positive output
requires the on state. These constraints also apply during mandatory initial
running or stopped periods. A conflict makes the model infeasible.

Optional gross-head limits use named reservoir or fixed-level references,
before waterway losses. With the same mode sign $sigma_g$ as the power
balance, let $h_(g t)^G=sigma_g (overline(Z)_(s t)-overline(Z)_(d t))$.
Running and starting restrictions have different assessment times:

$ u_(g t)=1 & ==> underline(h)_(g t)^G <= h_(g t)^G <= overline(h)_(g t)^G, \
  s_(g t)^+=1 & ==> sigma_g (Z_s (V_(s t))-Z_d (V_(d t))) >= h_(g t)^"start". $

A fixed boundary contributes its prescribed level in either expression.
The startup condition uses the surface at the control edge, not the subsequent
interval average. It does not stop a running unit when the level later falls
below the startup threshold. References remain fixed when tunnel gates move;
we require explicit references where a plant endpoint is a pressure junction.

#block(breakable: false)[
For interval-average discharge or release $x_t$, supplied upward and downward
rates $r_x^+,r_x^-$ limit changes between interval midpoints:

$ -integral_(m_(t-1))^(m_t) r_x^- (s) dif s
  <= x_t - x_(t-1) <= integral_(m_(t-1))^(m_t) r_x^+ (s) dif s. $ <eq:ramps>
]

A running unit cannot enter its positive operating range from zero when
the allowed inter-interval increase is smaller than its minimum output.
We therefore allow a minimum-output jump at startup and shutdown in the
power-ramp bounds. Add
$sum_(g in cal(G)_x) underline(P)_(g t) s_(g t)^+$ to the upward allowance
and $sum_(g in cal(G)_x) underline(P)_(g,t-1) s_(g t)^-$ to the downward
allowance. Here $cal(G)_x$ contains that unit or all units at the constrained
plant. Discharge and river-release ramps have no jump allowance.
For the first interval, use the average and midpoint of the supplied previous
window. Omit that comparison if no previous window is supplied.

Reservoir volume and level ramps instead constrain vertex changes
$V_(i,t+1)-V_(i t)$ and $Z_(i) (V_(i,t+1))-Z_(i) (V_(i t))$, with rate integrals
over $[tau_t,tau_(t+1)]$. The ramp limits apply between recorded states.

#block(breakable: false)[
A tunnel interlock limits how many gates in a group $cal(E)_ell$ may be
open at once. With $o_(e t)=bold(1)[b_(e t)>0]$,

$ underline(N)_(ell t) <= sum_(e in cal(E)_ell) o_(e t) <= overline(N)_(ell t). $ <eq:interlock>
]

A positive position counts as open even when its flow is zero.

Some handovers also require a plant shutdown after a gate closes. For each
interlock, specify protected plants $cal(P)_ell$ and disjoint closure groups
$cal(E)_ell^d subset.eq cal(E)_ell$, with $d$ denoting forward or backward.
Each group has a supplied waiting time $L_ell^d>=0$ in hours. Closing its gate
at $tau_t$ requires every protected unit to be off during the wait:

$ u_(g k) + o_(e,t-1) - o_(e t) <= 1, quad
  tau_t <= tau_k < tau_t + L_ell^d, \
  e in cal(E)_ell^d, quad g in union_(p in cal(P)_ell) cal(G)_p. $ <eq:switchwait>

At closure, the gate-status difference is one, so the inequality forces
$u_(g k)=0$. For every other transition it adds no restriction. It uses the
existing variables, with initial gate status supplied before the first interval.
Another intake may open immediately while the plant remains idle. Movement
between positive positions is not a closure, and reopening does not cancel a
wait. Controls change only at grid edges; a wait expiring inside an interval
permits generation at the next edge.

A pre-horizon closure can leave a residual wait $w_(ell d)^0$:

$ 0 <= w_(ell d)^0 <= L_ell^d, quad
  u_(g k)=0 quad "when" quad tau_k-tau_0 < w_(ell d)^0. $ <eq:initialswitchwait>

For a positive residual, the protected units must initially be off with off-age
at least $L_ell^d-w_(ell d)^0$. This makes their histories consistent with the
implied closure time. Zero declares no pending wait. On restart we retain
the latest unexpired closure deadline in each group, together with the gate
positions and unit histories. Closing a gate at the restart edge starts a new
wait in the continuation.

== Flow observations and river operating states

An operating observation can sum unit discharges, tunnel flows, upstream reach
entries and downstream arrivals. Its interval-average rate is

$ F_(o t) = sum_(g in cal(G)_o) q_(g t)
  + sum_(e in cal(E)_o) Q_(e t) + sum_(r in cal(R)_o^"in") R_(r t)
  + sum_(r in cal(R)_o^"out") frac(A_(r t), kappa Delta_t) + I_(o t). $

Tunnel contributions retain their declared direction. Each listed contribution
must describe water already present in the physical network. The observation
adds no water balance. Bounds may constrain this rate and its horizon volume:

$ underline(F)_(o t) <= F_(o t) <= overline(F)_(o t), quad
  underline(B)_o <= B_o^0 + kappa sum_t Delta_t F_(o t) <= overline(B)_o. $

Only supplied upper or volume bounds apply. A continuation adds executed
observed volume to $B_o^0$ and retains the terminal allowance. Downstream
observations use delayed arrivals even when no release schedule is supplied.

Some restrictions apply to a day or a preceding time window rather than to
one interval. Extend each observation as the constant $F_(o t)$ over its
interval and retain the supplied pre-horizon history. Let
$C_o(s)=kappa integral_(tau_1)^s F_o(v) dif v$. Its vertex values satisfy

$ C_(o,t+1)=C_(o t)+kappa Delta_t F_(o t), quad C_(o 1)=0. $

Linear interpolation gives $C_o$ inside an interval. Historical values are
known constants. Bounds on a window's volume constrain $C_o(b)-C_o(a)$.
For a rolling duration $w>0$, define

$ M_o(s;w)=frac(C_o(s)-C_o(s-w), kappa w), quad
  -d_o <= F_o(s)-M_o(s;w) <= u_o. $

These limits are flow differences, not rates of change per hour. The
deviation is affine between changes of $F_o$ and those changes shifted by
$w$. We therefore test both sides of each such knot. A pulse leaving the
history window can make a constraint bind between scheduling boundaries.
This construction adds no time steps. Cumulative variables keep the
integral equations sparse. Fixed-period ramps instead bound differences
between $(C_o(b)-C_o(a))/(kappa (b-a))$ for consecutive declared periods.
All averages use elapsed time, including when periods have unequal lengths.

The constant record is an operating observation. For downstream arrivals it
retains their original interval means, not their instantaneous waveform.
Cutting such an interval apportions its recorded volume in proportion to
time. Finer replay checks flow windows on these original intervals.

For a reservoir level range $D$ over a declared window, introduce one lower
level $ell$ and impose

$ ell <= Z_i(V_(i t)) <= ell+D $

at every vertex in the window, together with its historical extrema. This
has the same feasible projection as bounding every pairwise level difference
by $D$, with a linear rather than quadratic number of rows. Window endpoints
inside the horizon must be vertices. The increasing level law and the
piecewise-linear storage path place all extrema at these vertices; finer
replay also checks its additional vertices. A range may constrain the known
prefix of a window extending past the horizon. A volume total or fixed-period
mean requires its full period within the known past and scheduling horizon.
Restart retains the executed history and the original window times.

A routed river may have a minimum flowing rate $R_r^"min">0$. Binary state
$z_(r t)$ then imposes

$ R_r^"min" z_(r t) <= R_(r t) <= C_(r t) z_(r t). $

Starts $a_(r t)^+$ and stops $a_(r t)^-$ obey the transition equations above
with $z$ in place of $u$, and no minimum-duration requirement. Costs depend on
these transitions and on changes in upstream rate:

$ C^"river" = sum_(r,t) (c_(r t)^"change" abs(R_(r t)-R_(r,t-1))
  + c_(r t)^"start" a_(r t)^+ + c_(r t)^"stop" a_(r t)^-). $

The rate-change term also applies without a river operating state. Each cost
is charged once, without a duration multiplier. The first comparison uses
supplied previous flow; an unspecified previous flow omits that rate-change
comparison. Transition costs require a known initial state and a positive
minimum rate. Signed instantaneous channels have no such state. Finer replay
assesses these costs and flow observations over the original control windows.

== Explicitly penalized operating limits

A separate named rule may permit bounded deviations from a reservoir
observation or a unit or plant output target. Existing operating schedules
remain hard unless replaced explicitly by such a rule. For an observed
quantity $x_ell$, supplied limits $L_ell,U_ell$ use nonnegative slacks:

$ x_ell + sigma_ell^- >= L_ell, quad
  x_ell - sigma_ell^+ <= U_ell, quad
  0 <= sigma_ell^(plus.minus) <= overline(sigma)_ell^(plus.minus). $ <eq:softlimits>

Omit a side when its limit is absent. A zero allowance makes that side hard.
Positive allowances have strictly positive penalty coefficients. A soft output
target adds no commitment implication; generation can remain off and incur
the permitted shortfall.

Reservoir volume observations use $x_ell=V_i (theta_ell)$ and level
observations use $x_ell=Z_(i) (V_i (theta_ell))$ at a declared absolute grid
vertex $theta_ell$. Their slacks have units $"Mm"^3$ and metres. A terminal
target is an observation at the original terminal time, which is retained
when the case is continued. These rules constrain the recorded states and
make no assertion about a continuous-time band between observations.

For power and discharge, observe the average over an original window
$[a_ell,b_ell]$ whose endpoints are grid edges:

$ x_ell = 1/(b_ell-a_ell)
  sum_(t: a_ell <= tau_t < tau_(t+1) <= b_ell) Delta_t y_t. $ <eq:softaverage>

Here $y_t$ is one unit's output or the sum of the plant's member outputs.
The slacks have units MW or $"m"^3$/s. Separate windows are required when
one part of a window must not offset another.

The nonnegative cost is

$ C_ell = alpha_ell
  (c_ell^- sigma_ell^- + c_ell^+ sigma_ell^+), quad
  alpha_ell = cases(
    1 & "reservoir volume or level observation",
    b_ell-a_ell & "power window",
    kappa (b_ell-a_ell) & "discharge window",
  ). $ <eq:softcost>

Costs are therefore supplied per $"Mm"^3$, per metre, per MWh or per
$"Mm"^3$, respectively. State observation costs are charged once at the
observation, without a duration multiplier. Finer replay keeps these same
observations and windows and does not create additional penalty terms.
Restart drops completed assessments, retains future times and rejects a
restart inside an output window because partial assessment history is not
part of this input profile.

Water balance and every constitutive relation remain hard, together with
the static domains needed for valid equipment calibration and global bounds.
A reservoir volume rule can tighten its domain only to the expanded limits
$L_ell-overline(sigma)_ell^-$ and $U_ell+overline(sigma)_ell^+$, intersected
with its hard safety domain. The nominal limits alone cannot enclose all
admitted solutions. Independent assessment reports the actual minimal
slacks $max(0,L_ell-x_ell)$ and $max(0,x_ell-U_ell)$ and their costs. A
permitted positive slack is compatible with an accepted schedule; exceeding
an allowance or violating the physics is not.

All permitted deviations and their costs are explicit case data.
Zero allowances retain hard limits; no implicit endpoint penalty is added.
= Spinning reserve capacity

A reserve group $k in cal(K)$ has eligible generating units $cal(G)_k$,
upward and downward requirements $D_(k t)^(plus.minus)$, capacity ceilings
$overline(R)_(k t)^(plus.minus)$ and response time $theta_k>0$ in hours.
Allocation $r_(g k t)^(plus.minus)$ is MW. It is zero outside $cal(G)_k$ and
obeys reserve availability $a_(g t)^"res" in {0,1}$:

$ 0 <= r_(g k t)^(plus.minus) <= a_(g t)^"res" overline(P)_g u_(g t), quad
  sum_g r_(g k t)^(plus.minus) <= overline(R)_(k t)^(plus.minus). $ <eq:reserveallocation>

Here $overline(P)_g$ is nameplate power. Omitted group ceilings use the sum
of eligible nameplate power limits. A group
requirement is hard unless a positive shortfall cost is explicitly supplied:

$ sum_g r_(g k t)^(plus.minus) + sigma_(k t)^(plus.minus) >= D_(k t)^(plus.minus), quad
  0 <= sigma_(k t)^(plus.minus) <= D_(k t)^(plus.minus). $ <eq:reserveobligation>

For a hard requirement $sigma=0$. Positive shortfall costs make it equal
$[D_(k t)^(plus.minus)-sum_g r_(g k t)^(plus.minus)]_+$ at optimum. Independent
assessment recomputes that deficit from allocations and checks the stored value.

The following reserve operating-point equations apply to generating stations.
Pumps have zero reserve allocation and retain their nominal flow and
consumption, $q_(g t)^(plus.minus)=q_(g t)$ and $P_(g t)^(plus.minus)=P_(g t)$.
All groups share each generating unit's physical headroom:

$ R_(g t)^(plus.minus) = sum_k r_(g k t)^(plus.minus), quad
  P_(g t)^+ = P_(g t) + R_(g t)^+, quad
  P_(g t)^- = P_(g t) - R_(g t)^-, \
  q_(g t)^+ >= q_(g t), quad q_(g t)^- <= q_(g t). $ <eq:reservepoints>

Both operating points retain the nominal commitment $u$. They satisfy
@eq:onoff, @eq:power and @eq:unitdomain with their respective $q^+,P^+,h^+$
or $q^-,P^-,h^-$. The plant totals obey the same operative power/discharge
bounds and plant state $U$ as nominal generation. Downward capacity therefore
cannot stop a unit or cross its minimum operating output. No monotonic
relationship between power and discharge is assumed.

The two points hold the scheduled hydraulic node heads fixed, while shared
station tailwater is re-evaluated:

$ h_(p t)^(plus.minus) = H_("source"(p),t)
  - max{H_("target"(p),t), h_p^"floor"}
  - T_p (sum_(g in cal(G)_p) q_(g t)^(plus.minus)). $ <eq:reservehead>

Every operating unit in the station participates in these power equations.
A unit with no allocation may change its witness flow to hold nominal power
when another unit changes shared tailwater. Summing allocations in
@eq:reservepoints prevents selling the same headroom in multiple groups.

Let $theta_p$ be the minimum $theta_k$ among groups with at least one eligible
unit in the station. For a unit or plant output
$x in {P,q}$, response rates use the same supplied directional rates as
@eq:ramps, with no startup/shutdown jump:

$ x_t^+ - x_t <= integral_(m_t)^(m_t+theta_p) r_x^+ (s) dif s, quad
  x_t - x_t^- <= integral_(m_t)^(m_t+theta_p) r_x^- (s) dif s. $ <eq:reserveresponse>

The station-wide response is conservative for products with different
response times and also limits units compensating another group's tailwater.
Finer replay retains these physical response limits at its operating points.

Optional sustained activation duration $delta_k^"act">=0$ is distinct from
response time. A positive duration requires supplied unit qualifications
$E_g^+,E_g^-$ in MWh:

$ sum_k delta_k^"act" r_(g k t)^+ <= E_g^+, quad
  sum_k delta_k^"act" r_(g k t)^- <= E_g^-. $ <eq:reservequalification>

These qualifications apply independently to each interval. They are not a
consumable reservoir inventory and remain unchanged at restart. Capacity
operating points add no flow or energy activation to nominal water balances.
They do not re-solve tunnel losses, storage evolution, river transport or
environmental flows for an activated trajectory. Cascade activation requires
a separate treatment, such as the two-stage stochastic and robust models
studied by Naversen et al.~@naversen. A certificate here concerns only the
declared capacity model at scheduled node heads. Pumping reserves, nonspinning
starts/stops, frequency-response droop, forbidden-zone transitions and
reserve market purchases are outside its scope.

#block(breakable: false)[
= Economic objective

Generation at price $pi_t$ per MWh earns $pi_t Delta_t P_(g t)$; pumping
consumption incurs that amount. Let $sigma_g$ be the sign of the unit's plant
mode. Net electrical exchange is $sum_g sigma_g P_(g t)$.
Let $c^+,c^-$ be unit transition costs, $c^b$ nonnegative gate-movement costs,
$c^r$ nonnegative release-shortfall costs, and $v_i,v_r$ the values assigned
to water in reservoirs and reaches, per $"Mm"^3$. We maximize

$ cal(V) = & sum_t pi_t Delta_t sum_g sigma_g P_(g t)
  - sum_(g,t) (c_(g t)^+ s_(g t)^+ + c_(g t)^- s_(g t)^-) \
  & - sum_(e,t) c_(e t)^b abs(b_(e t) - b_(e,t-1))
  - sum_(r,t) c_(r t)^r xi_(r t) - sum_ell C_ell - C^"river" \
  & + sum_i v_i (V_(i T) - V_(i 0)) + sum_r v_r (W_(r T) - W_(r 0)) \
  & + sum_(k,t) Delta_t (pi_(k t)^+ sum_g r_(g k t)^+
      + pi_(k t)^- sum_g r_(g k t)^-
      - c_k^"res" (sigma_(k t)^+ + sigma_(k t)^-)). $ <eq:objective>
]

Starts, stops and gate movements are charged once per event, without a duration
multiplier. The first gate comparison uses its supplied initial position.
The reservoir-storage and reach-inventory terms value changes in stored water.
A late release can still have value downstream after the horizon ends.
The terminal values $v_i,v_r$ are inputs that account for this later use.
Reserve prices $pi_(k t)^(plus.minus)$ and positive optional shortfall costs
$c_k^"res"$ are currency per MW-hour. Allocated MW earns capacity price
times interval hours, without implying actual activation energy. Hard reserve
requirements have zero shortfall and no penalty term.

An energy value $e_i$ in currency per MWh can be converted at a supplied
reference net head and combined efficiency:

$ f_i=2.725 h_i^"ref" eta_i^"ref", quad
  v_i=e_i f_i+v_i^"downstream". $ <eq:energyvalue>

$f_i$ is MWh per $"Mm"^3$. The downstream term includes only opportunities
not already represented by $f_i$. Water that has passed one station retains
only its remaining opportunities. Heads, efficiencies, paths and exchange
rates are not inferred. All monetary inputs use one currency.
These constants describe future use; they are not marginal values calculated
by optimization. Historical transit is included in $W_(r 0)$, and subtracting
initial inventory values preserves additive accounting across restart.

A positive release penalty permits a deficit below the entry minimum.
The deficit variable is measured as a volume:

$ xi_(r t) >= 0, quad xi_(r t) >= kappa Delta_t (R_(r t)^min - R_(r t)). $ <eq:shortfall>

A positive cost makes $xi$ equal the actual shortfall at optimum. Without
that penalty, the entry minimum is hard and $xi=0$. Arrival minima and
aggregate flow requirements remain hard. Subtracting the values of the
initial inventories removes a constant from the reported objective and
leaves the optimal decisions unchanged.

= Solving the model and checking a schedule
<sec:solution>

== Finding feasible controls

We build the models in JuMP~@jump. SCIP first solves a MILP approximation
to propose unit states and gate positions. We fix those choices and use
Ipopt~@ipopt to solve for the continuous decisions under the nonlinear
hydraulic laws. We dispatch each proposal before constructing the next one.
Without an accepted trajectory, a short solve with constant objective first
seeks a feasible point. Revenue optimization follows, whether or not that
attempt succeeds. An accepted feasibility point is retained if optimization fails.
The best accepted trajectory on the proposal grid supplies the reference for the next
approximation and the starting point for its nonlinear dispatch. All attempts
share a time budget. Candidates that pass the checks below become starting solutions
for the full SCIP search~@scip, which can change the discrete choices.

The proposal also selects enabled river operating states; nonlinear dispatch
holds those states fixed. The MILP objective approximates @eq:objective and gives no enclosing bound.
When the proposal omits tunnel pressure equations, the permitted head ranges
and gate positions supply necessary flow bounds. We do not intersect these
exact-law bounds with a secant pressure approximation, whose approximate
flows need not lie in the same interval. Later proposals limit the total
number of changed unit states, river states and gate positions around the
accepted trajectory. A different gate position counts as one change. The
neighborhood starts at four changes, doubles after an unsuccessful attempt,
and resets after a gain.

Its level observations use first-order affine approximations at the reference
storage, or at initial storage when no reference is available. Each approximation
equals the physical level at its reference point. Its power observations use the
proposal's power surrogate. The
nonlinear dispatch and global model retain the original observation laws.
Dispatch with the discrete choices fixed is a piecewise-smooth, nonconvex
problem. A stationary solution need not be globally optimal.

== Reconstructing the water trajectory

An independent integration holds the selected unit discharges and openings
fixed and recovers the resulting water trajectory. At each interval, Newton's
method solves jointly for ending storage, junction heads, tunnel flows and
active outlet releases. Column-pivoted QR handles rank deficiency in the
linearized equations, including cycles with zero tunnel flow.
Nonnegative outlet and pump flow domains are imposed on Newton trial iterates;
convergence still requires all the original equations. The jointly converged
flows are retained, since reevaluating separate outlet laws can disturb a
water balance by accumulating their rounding errors. With a delay
shorter than the interval, some current release arrives before the interval
ends and belongs in the same hydraulic
solve. Equation @eq:power then gives electrical output.

In some isolated pressure components, the equations determine only head
differences. A common shift in all heads leaves the flows unchanged.
We choose one reference head for the numerical solve and keep every water
balance. We can then shift the recovered heads to satisfy the junction
bounds if the bounds admit such a shift.
Replay retains the candidate's pressure reference where the equations permit
this freedom. An idle committed unit can still have an operating head bound
relative to another component, so choosing a different reference can wrongly
make its operating point appear inadmissible. The reference does not replace
any hydraulic equation or absolute head restriction.

The nonlinear step need not have a unique solution. For example, a channel
can be passive at each reference level while its conductance decreases so
strongly with that level that the interval equations have several roots.
Newton convergence alone does not select the continuous physical trajectory.
Finer time steps are needed to assess such cases; a failed replay does not
establish mathematical infeasibility.

Even a constant channel conductance $K$ can cause coarse-grid overshoot.
With equal affine level slopes $s$, the difference $d$ in reservoir levels
obeys $d_"new"=(1-theta)/(1+theta) d_"old"$, where
$theta=kappa Delta K s$. For $theta>1$, the levels cross during a discrete
step; continuous passive equalization approaches equilibrium without crossing.

We accept a reconstructed schedule only if it satisfies the interval equations
and operating restrictions. We repeat the integration on a finer grid using
the same controls. Schedules and ramps are checked over their original
operating windows. We compare storage at the original vertices, total water,
cumulative arrivals and terminal inventory. For the specified volume
tolerance $epsilon_V$ and relative objective tolerance $epsilon_(cal(V))$,
the reconstruction must satisfy

$ max_(i,t) abs(V_i^"fine" (tau_t) - V_(i t)) <= epsilon_V, quad
  max_r abs(W_r^"fine" (tau_T) - W_(r T)) <= epsilon_V, \
  abs(cal(V)^"fine" - cal(V)) <= epsilon_(cal(V)) max(1, abs(cal(V))). $ <eq:replaycriteria>

The comparison estimates the effect of time discretization for the chosen
controls. It does not give a rigorous bound on continuous-time hydraulic error.

The discrete-gradient chain rule removes artificial storage-cycle work.
For direct internal transfers, multiply each node balance by its own interval
head and sum: the internal edge terms become minus transferred volume times
source-minus-target head. Passive tunnels and channels dissipate this work.
Ideal pumping and generating units exchange it without loss; efficiencies
below one add conversion losses. Summing potential changes around a closed
storage cycle gives zero, even when the interval lengths differ.

This identity concerns storage and direct transfers. A delayed river carries
no mechanical-energy state in this model. Its release and arrival energy
ports do not cancel as an instantaneous transfer, and the storage chain rule
does not establish an energy balance for the omitted river dynamics.
Nor does energy consistency remove coarse-step overshoot or nonunique roots.
The interval-model objective and its bound still concern an approximation to
continuous physical time.

== Bounding the remaining objective improvement

A feasible schedule gives an achievable objective value $L$. Global search
provides an upper bound $U$ that no feasible schedule in the declared model
can exceed. SCIP obtains such bounds from relaxations, tightens them and
branches on discrete choices or continuous variable ranges. A region can
be discarded when its bound cannot improve the current feasible value.

#block(breakable: false)[
A _relaxation_ enlarges the feasible set. Its maximum is an upper bound on
the scheduling objective because every feasible schedule is still in that set.
The enlarged set may allow flow, head and efficiency values whose power
exceeds the physical conversion law. A _valid inequality_ holds for every
feasible schedule and can exclude some of these excessive values.
]

The algebraic curve representations and bound constructions are given in
@sec:graphs.

For numerical bounding, each scalar equality $f(x)=b$ is enclosed by
$b-delta <= f(x) <= b+delta$, with $delta=10^(-9)$ in the equation's model
coordinates. The solver feasibility tolerance is $10^(-10)$. This explicit
band prevents small inconsistencies between floating-point curve graphs from
making the bounding problem empty. It enlarges the feasible set, so its
upper bound also bounds exact equality solutions. The band is used only in
global search; reconstruction and acceptance use the original equations.
It does not enclose every error in coefficient arithmetic or in the solver.

Here $L$ is the independently reconstructed value of a schedule feasible for
the interval model, and $U$ must be a finite enclosing bound. The absolute
and relative limits on possible improvement are

$ g_"abs" = max(0, U-L), quad g_"rel" = frac(max(0, U-L), max(1, abs(L))). $ <eq:gap>

If $L=100$ and $U=103$, improvement is limited to $3$ monetary units, or
$3\%$ under this definition. There may be no better schedule at all.
The constant subtraction of initial water values in @eq:objective affects
the relative gap's denominator, so $L$ and $U$ must use that same objective.

A certificate states that the requested absolute or relative gap has been
met for the interval model, interpolation rules and finite gate sets.
Fixing discrete choices makes the certificate conditional on those choices.
A schedule may meet the interval-model gap and fail the finer-grid checks.
Neither test certifies a continuous-time optimum. The numerical certificate
uses solver tolerances and floating-point bounds. It does not rigorously
enclose every rounding error in the solver and polynomial calculations.

#block(breakable: false)[
= Discussion

The hydraulic equations let changes in reservoir levels and shared tunnel
flows affect dispatch during optimization. River inventory accounts for
released water until it reaches the next node. Using these nonlinear laws
can be computationally expensive. Even with commitment fixed, the continuous
problem is nonconvex, and a large network may take substantial time to close
the objective gap.

Model accuracy also depends on the input curves. Finer integration can expose
an inadequate time grid, but it cannot correct a miscalibrated efficiency or
travel-time response. Tunnel positions absent from the declared finite set
are excluded from the search. The objective certificate applies to this
finite model, with its specified curves and control grid. Changing those
assumptions requires solving a new problem.
]

#block(above: 1.6em, below: 0pt, sticky: true)[
  #line(length: 100%, stroke: 0.5pt + rule)
  #v(6pt)
  #text(size: 9pt, fill: muted, smallcaps[Appendix])
]
#counter(heading).update(0)
#set heading(numbering: "A.1", supplement: [Appendix])

= Representing curves and strengthening bounds
<sec:graphs>

The scheduling problem is a mixed-integer nonlinear program (MINLP).
Fixing unit states and gate positions still leaves a nonconvex problem in
flow, head, efficiency and travel response. SCIP needs algebraic equations
for these functions and relaxations from which to compute upper objective
bounds.

== Representing supplied curves exactly

For a curve $y=f(x)$, the _graph_ is the set of pairs $(x,y)$ satisfying
that equation. We reproduce the specified interpolant in the optimization
model. Exactness here refers to the interpolant, not to the watercourse it
approximates. A piecewise-linear curve with knots $(x_i,y_i)$ has the graph

$ lambda_i >= 0, quad sum_i lambda_i = 1, quad
  x = sum_i x_i lambda_i, quad y = sum_i y_i lambda_i, quad lambda in op("SOS2"). $ <eq:sos>

SOS2 restricts the weights to one segment~@huch. The two neighboring weights
sum to one, and the equation for $x$ fixes the interpolation fraction.
The equation for $y$ then gives the interpolated output. Without adjacency,
weights on distant knots could produce points away from the curve.

For a reservoir level table, integrate each linear cell to obtain a quadratic
potential. The endpoint level and potential share the same coordinate weights.
Writing $delta v_k=v_(k+1)-v_k$ and
$s_k=(Z_(k+1)-Z_k)/delta v_k$, the potential graph is

$ Psi(v)=sum_k Psi(v_k) lambda_k
  - sum_k (s_k delta v_k^2)/2 lambda_k lambda_(k+1). $ <eq:potentialgraph>

The correction removes the quadratic's chord error. For endpoint volumes
$a,b$, impose the work equation and a bounded divided-difference row:

$ (b-a) H=Psi(b)-Psi(a), quad H=Z(a)+D(b-a), quad
  (min_k s_k)/2 <= D <= (max_k s_k)/2. $ <eq:storagegraph>

The slope bounds range over cells intersecting the permitted storage domain.
They follow from
$D=integral_0^1 (1-theta) Z'(a+theta(b-a)) dif theta$.
For $a != b$, the first row fixes the mean head. At $a=b$, it is only $0=0$;
the second row fixes $H=Z(a)$ without a direction binary or a volume-change
threshold. Local dispatch and replay evaluate the segment integral directly.
The global model subtracts one common elevation reference from $Z$ and
integrates that shifted level, which leaves the equations unchanged while
reducing cancellation. Polynomial storage uses @eq:polynomialhead directly.

For a turbine table, let $q_i$ be the sorted union of the original discharge
knots. These coordinates partition the interpolants; they are not newly
measured calibration points. Discharge weights $lambda_i$ and head weights
$mu_j$ each satisfy the corresponding SOS2 and coordinate constraints.
Let $E_(i j)=f_j(q_i)$ be the original curve's value at a partition point.
Linear interpolation in head combines the
values $eta_j$ on the neighboring head curves. For head-independent data,
each head vertex uses the same efficiency column. A shape-preserving cubic in
discharge~@fritsch is represented by adding corrections to the linear curve:

$ w_i = lambda_i lambda_(i+1), quad r_i = lambda_i w_i, \
  eta_j = sum_i E_(i j) lambda_i
    + sum_i [A_(i j) r_i + B_(i j) (w_i - r_i)], quad
  eta = sum_j mu_j eta_j. $ <eq:cubic>

For segment width $delta q_i=q_(i+1)-q_i$, secant slope
$d_(i j)=(E_(i+1,j)-E_(i j)) \/ delta q_i$ and original polynomial derivatives
$m_(i j)^"left",m_(i j)^"right"$ at the cell's two endpoints, use

$ A_(i j) = delta q_i (m_(i j)^"left" - d_(i j)), quad
  B_(i j) = delta q_i (d_(i j) - m_(i j)^"right"). $ <eq:hermitecoefficients>

On the active segment, the two correction terms are multiples of
$lambda_i^2 lambda_(i+1)$ and $lambda_i lambda_(i+1)^2$. They vanish at the
endpoints and adjust the derivatives to match the
original slopes. The variables $w_i,r_i$ express the corrections through
quadratic equalities. Setting $A=B=0$ gives bilinear interpolation.
Each original cubic is re-expressed on its union-grid subcells, retaining the
derivatives from that polynomial. Recalculating PCHIP slopes from the
partition-point values would generally change the interpolant and is not
performed. Endpoint derivatives are taken from within each selected cell,
including the linear continuation outside its original discharge data.
Restricting the curve's domain leaves its cubic slopes unchanged. We use the
declared end secants for discharge outside the table and for head extrapolation
where allowed.

An off unit can have a network head outside its running table. A finite
extension lets us evaluate the graph there. Its running head and efficiency
limits are inactive, and $q=P=0$.

Units with the same net-head variable and compatible table grids can share
head weights. Each release cohort uses the same flow weights for arrivals,
environmental checks and terminal inventory. Channel responses use a
piecewise-bilinear graph. The product equalities in @eq:cubic, @eq:power and
@eq:arrival remain nonlinear.

== Quadratic equations for tunnels and gates

A binary direction variable separates positive and negative tunnel flow:

$ Q = Q^+ - Q^-, quad 0 <= Q^+ <= C z, quad 0 <= Q^- <= C(1-z), \
  S = (Q^+)^2 - (Q^-)^2, quad z in {0, 1}. $ <eq:signedloss>

At $z=1$, flow is nonnegative and $Q^-=0$. At $z=0$, flow is nonpositive
and $Q^+=0$. This gives $S=Q abs(Q)$ in either case. Directed outfalls use
$S=Q^2$. Known head ordering can also fix the direction.

For each finite gate position, a binary selector identifies whether that
position is chosen. Exactly one selector is active, so $b$, the open state
$o$ and the selected resistance $K(b)$ are weighted sums of known values.
To impose the loss equation only when the gate is open, define
$d_o=o(H_v-H_w)$. If $L<=H_v-H_w<=U$, this binary product is represented by

$ L o <= d_o <= U o, quad
  H_v - H_w - U(1-o) <= d_o <= H_v - H_w - L(1-o). $ <eq:binaryhull>

When $o=0$, @eq:binaryhull forces $d_o=0$. When $o=1$, it forces
$d_o=H_v-H_w$. Adding $abs(Q)<=C o$ and $d_o=K(b)S$ gives the open and
closed gate laws in @eq:tunnel and @eq:closed. The resistance assigned to
the closed selector may be any finite value, since $S=0$ there. All the
product equalities are quadratic.

The outlet law also has a binary wet/dry branch for
$ell=max(H-z^"crest",0)$. Its square root in @eq:outlets can be written
using $r>=0$ and $r^2=ell$.

== Upper bounds for global optimization

We bound variables using storage and equipment limits. Water balances
propagate initial and terminal restrictions through the horizon to narrow
these ranges. The bounds must hold for every admissible schedule.

#block(breakable: false)[
We omit unit and time indices in the following bounds.
These strengthening rows apply only to generating modes. Pump rows retain
the original bounded work equality @eq:pump, without applying generation
supports to consumption.
For a running generating unit with nonnegative $q h eta$, let $overline(e)$ bound its
electrical efficiency. Then $p^"ub" (q,h)=gamma overline(e) q h eta(q,h)$
is an upper bound on physical electrical power. It is polynomial within
each interpolation cell even when @eq:power is implicit in $P$.
Suppose a plane $a q+b h+c$ bounds this polynomial throughout the running
domain, with #box[$b>=0$]. For a head lower bound $h_*$ valid in both operating
states, the following inequality includes the off state:

$ P <= a q + b(h - h_*) + (c + b h_*) u. $ <eq:support>
]

When the unit runs, it reduces to $P<=a q+b h+c$. When it is off,
$P=q=0$ and the right side is nonnegative. To establish $c$, we bound the
residual $p^"ub"-a q-b h$ over every original interpolation cell.
In a Bernstein representation of a polynomial on a box, the basis weights
are nonnegative and sum to one, so its value cannot exceed its largest
coefficient. The coefficient bound covers the whole cell. Fitting a plane at a few
samples alone could miss a larger residual between them.

A box containing the admissible running values of $q,h,eta$ has up to eight
corners $(q_c,h_c,eta_c)$. We couple these variables with weights on the
corners. Exact binary-product inequalities, as in
@eq:binaryhull, define $h_"on"=u h$ and $eta_"on"=u eta$.
Nonnegative weights $z_c$ on the corners satisfy

$ sum_c z_c = u, quad sum_c q_c z_c = q, quad
  sum_c h_c z_c = h_"on", quad sum_c eta_c z_c = eta_"on", \
  gamma e_min sum_c q_c h_c eta_c z_c
  <= P <= gamma e_max sum_c q_c h_c eta_c z_c. $ <eq:box>

Here $e_min,e_max$ bound electrical efficiency over the power range.
Their use assumes $q h eta>=0$, which follows from @eq:power for nonnegative
power and positive electrical efficiency. Commitment alone does not imply
$h>=0$: an idle unit with zero minimum power may have any head in its declared
running domain. The corner box retains that domain, and bounds on $h_"on"$
include both zero and any negative running heads.
For any point in the box, multilinear interpolation weights reproduce its
coordinates and the product $q h eta$. Hence the inequalities include every
physical running point. When $u=0$, all corner weights vanish. The original
efficiency and power equalities still determine the physical output.

For a tabulated turbine, we can use its head knots $h_j$ and weights $mu_j$
to bound power. Define $zeta_j=mu_j$ for a running unit and $zeta_j=0$ for
an off unit through the following constraints:

$ 0 <= zeta_j <= mu_j, quad sum_j zeta_j = u, quad
  sum_j (h_j - h_1) zeta_j = h_"on" - h_1 u, \
  P <= a q + sum_j c_j zeta_j. $ <eq:coordinatecuts>

The coefficients $c_j$ are chosen so that the weighted expression bounds the
polynomial upper bound $p^"ub"$, including its variation between head knots.
The constraints give $zeta=mu$ when $u=1$, and $zeta=0$ when $u=0$ without
fixing the head of an off unit. Discharge weights support analogous bounds
along the other table axis. We may add further planes during the search after verifying each over its
complete domain.

#block(breakable: false)[
  #bibliography("references.bib", style: "ieee", title: "References")
]
