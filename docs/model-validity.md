# Mathematical scope and validity

OpenSHOP solves a deterministic, finite-horizon, discrete-time model. Its
nonlinear equations describe quasi-steady hydraulics at midpoint storage.
They are not an exact time integration of the watercourse. This distinction
matters even when a solver closes the objective gap.

The [paper](../paper/paper.pdf) gives the equations. The
[Lean project](../test/formal/README.md) checks selected identities over the
real numbers. Julia regression tests connect those statements to numerical
examples. The Julia implementation, its translation through JuMP, and the
solvers have not been formally verified.

Input validation does not promise a nonempty feasible set. Physical limits,
water availability and operating requirements can conflict in a particular
case. That is different from a redundant formulation inequality incorrectly
removing a feasible point.

## Conservation and passive hydraulic laws

An internal connection withdraws and credits the same signed water flow.
A plant or tunnel outfall credits its receiving river instead of its hydraulic
head-reference node. A river's natural inflow enters that river directly;
it is not withdrawn from reservoir storage. Zero-storage confluences allocate
upstream arrivals with prescribed shares summing to one. Local branch injections
are added after this allocation.

Adding all reservoir and river-inventory balances cancels internal transfers.
Only external inflows and boundary exchanges change total water inventory.
This identity does not require a tree network, a preferred tunnel direction,
or equal branch parameters. Historical river inventory belongs in the initial
state; omitting it changes the problem.

For an open passive tunnel, `ΔH = K q |q|` with `K>0`, so
`ΔH q = K |q|³ ≥ 0`. A closed gate imposes `q=0` without imposing `ΔH=0`.
The supplied delta-level curves are nondecreasing and pass through `(0,0)`.
Their convex interpolation preserves `ΔH q ≥ 0`. These laws cannot create
hydraulic work at their evaluation points. Tunnel laws assume submerged
connections and omit inertia, water hammer and transient pressure storage.

At the same head and transferred volume, efficiencies in `(0,1]` ensure that
generation cannot exceed ideal hydraulic work and pumping consumes at least
that work. Pump consumption divides by hydraulic and motor efficiency. These
pointwise inequalities do **not** establish an energy balance over a
nonlinear, time-discretized storage cycle.

## Midpoint energy error

Let `Z(V)=z+sV+cV²` and let `Ψ(V)=zV+sV²/2+cV³/3` be its potential integral.
For an interval with endpoint volumes `a,b`, real arithmetic gives

```math
\Psi(b)-\Psi(a)
=(b-a)Z\!\left(\frac{a+b}{2}\right)+\frac{c(b-a)^3}{12}.
```

Thus midpoint work is exact for affine storage-level curves, but generally
not for curved ones. On the cycle `0 → 2 → 1 → 0`, midpoint pumping work minus
generating work is `-c/2`, despite returning to the initial storage. Adding a
constant to every volume leaves this cycle error unchanged.

The regression test uses an ideal pump and turbine, a curved reservoir and a
closed water cycle `1 → 3 → 2 → 1` Mm³. Its coarse-grid artificial electrical
gain is 0.2725 MWh. Halving each interval reduces that gain by a factor of four;
quartering reduces it by sixteen. Losses can mask this error, so a plausible
round-trip efficiency is not a discretization test.

Finer replay estimates grid sensitivity for the chosen controls. One comparison
is not an error enclosure or a convergence theorem. For sensitive cases, compare
successively finer grids and judge absolute energy as well as relative objective
change. An energy-preserving discrete storage head would require a formulation
change; it is not supplied by the present midpoint equations.

## Existence, uniqueness and pressure references

A nonlinear solver can find a root without establishing existence or uniqueness
throughout the permitted domain. Passivity alone does not imply that an implicit
midpoint step has a unique root. Reference-level-dependent channel curves are
one concrete example.

Take two reservoirs with `Z(V)=V`, initial volumes `(3,1)`, and total volume four.
Let `x` be the first reservoir's ending volume. In volume-per-interval units,
use a channel with `q=K(reference) ΔH` and linearly interpolate the positive
conductances `3,1,1/3` at reference levels `9/4,5/2,11/4`. Its midpoint equation is

```math
x=3-K\!\left(\frac{3+x}{2}\right)(x-1).
```

The distinct values `x=1.5,2,2.5` all solve it. Every reference curve is passive
and increasing in delta head. Lean verifies this counterexample; Julia evaluates
all three roots with zero hydraulic residual. Their admissibility as discrete
roots does not make them the same continuous trajectory.
This is adversarial synthetic input, not a calibrated channel: the strong decrease
in conductance is admitted by the input contract, which checks passivity in delta
head but does not impose joint monotonicity in both reservoir levels.

Even a constant conductance can produce coarse-grid overshoot. For equal affine
level slopes `s`, write `r=0.0036 Δt K s`. The difference in reservoir levels
satisfies `d_new=(1-r)/(1+r) d_old`. If `r>1`, the discrete levels cross, whereas
continuous passive equalization approaches equilibrium without crossing it.
The tests verify this effect and remove the reversal by refining the step.

Consequently, a replay failure can mean a different root or a Newton convergence
failure, rather than mathematical infeasibility. Acceptance checks fail closed;
they can reject feasible discrete points. Strong level dependence, long intervals
and stiff outlets call for grid refinement. No general uniqueness or solver
convergence guarantee is claimed.

An isolated junction component can also have an arbitrary common pressure shift.
That shift preserves tunnel head differences, but need not preserve an idle
unit's operating head relative to another component. `head_reference` supplies
Newton guesses and selects this pressure reference during reconstruction.
Validation and finer replay retain the candidate's reference. Every balance,
hydraulic equation and absolute operating bound remains checked; supplied heads
cannot override an anchored hydraulic solution.

## Commitment and exact algebraic representations

Commitment and physical output are separate. With zero flow and power minima,
a unit may be on and idle. A positive minimum power, nonnegative flow and
nonnegative efficiency imply positive generating head. Commitment by itself
does not. In particular, an idle unit can have negative head if its declared
operating domain permits it.

The global power envelope therefore retains the full running head interval.
Its on-head moment uses bounds including zero and both ends of that interval.
The four binary-product inequalities remain valid for signed heads. Rank-one
corner weights reproduce `q h η` for signed head coordinates as well.

Cubic discharge interpolation is retained exactly in real arithmetic through
SOS2 coordinates and product corrections; it is not replaced by a chord.
The Lean identities check the Hermite coefficients, affine cell changes,
cubic correction and auxiliary bounds. They do not verify the SOS2 solver
implementation or the floating-point calculation of knots and coefficients.

Piecewise interpolation is not globally smooth. A derivative at a knot uses
the selected segment. The dry-outlet equation `w²=max(ΔH,0)` avoids evaluating
a square-root derivative at zero, but can still have a singular constraint
Jacobian at `w=0`. Ipopt's usual smoothness and constraint-qualification
assumptions need not hold there. A successful local termination is not a proof
of feasibility or optimality; independent residual checks remain necessary.

## Transport, operating rules and objective bounds

Delayed reaches are directed transport responses, not pressure-flow laws. Their
directions and one-level outlet responses must be physically applicable over the
declared operating range. Input validation does not check downstream backwater
or gravitational elevation compatibility. A supplied uphill river connection can
conserve water and pass finer replay while creating unphysical potential energy.
Unlike tunnels and delta-level channels, ordinary routed reaches have no
head-difference equation that would prevent this. A one-level discharge curve
cannot represent independently varying downstream counterpressure. Likewise,
a plant/tunnel hydraulic target must supply the actual outfall head, even when
its water enters a separately named delayed reach. These are calibration and
connection assumptions, not consequences of the conservation proofs.

Each release cohort keeps the travel-response mixture selected by its own
release rate. A convex mixture of causal probability responses remains causal
and volume preserving. Deterministic delays preserve source pulses through
confluences. Distributed routing at confluences instead rebins arrivals uniformly
on each interval: it conserves volume but approximates timing.

The mathematical identities assume response weights and branch shares sum
exactly to one. Input checks admit tiny floating-point discrepancies in those
sums, and compiled transport uses floating-point arithmetic. Numerical water
conservation is therefore checked with tolerances, not asserted as a bitwise or
interval-arithmetic guarantee.

Dwell times, ramps and observation windows use physical hours. Terminal dwell
obligations are reported for continuation; they do not prove that a future
continuation exists. Terminal water values are supplied economic assumptions.
Reserve capacity witnesses hold scheduled hydraulic node heads fixed and do not
certify a physically activated tunnel/river trajectory.

For maximization, a genuinely feasible value `L` and an enclosing real upper
bound `U` imply that no feasible improvement exceeds `U-L`. MILP proposals and
local stationary points supply no such upper bound. Tightening a relaxation is
valid only if it retains every feasible point, including idle and off states.

SCIP bounds, polynomial coefficient bounds and schedule checks use floating-point
arithmetic and feasibility tolerances. Guard margins reduce rounding exposure;
they are not directed-rounding proofs. `global_certificate` is a numerical gap
certificate for the declared interval model and any stated fixed decisions.
It is neither a Lean certificate of a solver result nor a bound on the optimum
of the continuous physical system. Original-grid validation, finer replay and
the global gap answer different questions.
