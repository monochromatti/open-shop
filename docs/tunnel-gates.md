# Discrete tunnel gate decisions

An empty `Tunnel.gate_positions` keeps the gate prescribed by `opening` and its
operating series. A nonempty vector lets the scheduler choose one listed
position in each interval. The positions must be finite, strictly increasing,
unique and within `[0,1]`. Both `opening` and any explicit `initial_opening`
must be listed. This is discrete gate optimization; intermediate unlisted
positions are not admissible controls.

```julia
high = Tunnel(name=:HighIntake, source=:HighLake, target=:Intake,
    resistance=1.0, capacity=6.0, opening=1.0,
    gate_positions=[0.0, 1.0], initial_opening=1.0,
    gate_adjustment_cost=5.0)
low = Tunnel(name=:LowIntake, source=:LowLake, target=:Intake,
    resistance=1.0, capacity=6.0, opening=0.0,
    gate_positions=[0.0, 1.0], initial_opening=0.0,
    gate_adjustment_cost=5.0)
interlock = TunnelInterlock(name=:IntakeSelection,
    tunnels=[:HighIntake, :LowIntake], min_open=1, max_open=1)
```

Pass these tunnels and `tunnel_interlocks=[interlock]` to `HydroSystem`.
An interlock counts listed positions strictly greater than zero: a gate at
`0.5` is one open gate. `min_open` defaults to zero; `max_open` defaults to the
number of members.
The members must be distinct existing tunnel names. Each interval must satisfy
the lower and upper counts. Interlocks add operating restrictions; hydraulic
junction continuity still determines how water flows through open paths.
Their `min_open` and `max_open` counts also accept integer operating series.

## Hydraulics and movement cost

Opening zero enforces zero tunnel discharge and leaves the endpoint heads
independent. It never imposes equal pressure across a closed gate. With no
gate-loss table, positive position `s` satisfies

```math
s(H_u-H_d)=R q|q|.
```

With `gate_loss_factor_curve`, the positive-position law is

```math
H_u-H_d=\bigl(R+R_{gate}(s)\bigr)q|q|.
```

The supplied curve maps normalized position to additional resistance; see
[channel and gate-loss curves](channel-curves.md). Unrouted tunnels can reverse
flow. A `discharge_river` outfall retains its forward-only discharge domain.

The global model represents `K` positions with `K−1` binary selectors; the
first position is implicit. It retains one physical tunnel flow per interval.
Prescribed positions introduce no gate decision variables, and local nonlinear
dispatch fixes the selected positions before building its equations.

For a selectable closed/open gate, an auxiliary head difference equals
`open * (H_u-H_d)`. The four inequalities of the bounded binary-product hull
enforce this identity with linear constraints. Flow capacity is also multiplied
by the open state. The pressure-loss equation then remains quadratic and
imposes no pressure relation when the gate is closed.

When positive positions have different resistances, the selected resistance is
an affine expression in the selectors. A single signed-square variable
`w = q|q| / 2500` gives the quadratic equation
`head_difference_when_open = 2500 * selected_resistance * w`. The factor 2500
scales flow by 50 m³/s; it does not change the physical law. This avoids cubic
products and separate flows for each position. Reversible flows retain the
existing positive/negative split; a head domain that fixes the flow direction
needs no additional direction binary.

Movement cost is charged once per position change:

```math
C_{gate}=\sum_{e,t} c_{e,t}|s_{e,t}-s_{e,t-1}|.
```

The first comparison uses `initial_opening`, or `opening` when it is omitted.
The coefficient is `gate_adjustment_cost`, optionally changed by an
`OperationalSeries` of the same name, evaluated at the interval start.
Cost is in objective currency per unit
normalized position; it is not multiplied by interval duration. Moving `0 → 0.5`
costs half as much as `0 → 1` under the same coefficient. The solution reports
the total in `gate_adjustment_cost`, which is subtracted from its objective.

## Prescribing and solving gates

An `opening` operating series prescribes positions from its first knot onward.
On an optimized tunnel, the intervals before that knot remain free. On a tunnel
with empty `gate_positions`, the ordinary default opening applies beforehand.
Use scheduling edges for opening knots, and supply listed positions for an
optimized gate. Prescribed openings and interlocks must be mutually feasible.

`solve(case)` optimizes gate positions and unit commitment together. Fixing
`fixed_u` restricts unit commitment but still lets eligible tunnel gates vary.
Passing `fixed_tunnel_opening` fixes gates and makes the corresponding
certificate conditional on those openings.
A complete accepted solution contains `tunnel_opening`, with one row per tunnel
and one column per interval, including prescribed tunnels.
Discrete positions are extracted from integral selectors, so arbitrarily close
listed positions retain their identity. Supplied fixed and physical control
schedules must use exact declared positions; tolerance never turns a closed
gate into an open one.

Local nonlinear dispatch fixes both discrete control families. Pass
`tunnel_opening` to `solve_case` or `solve_verified`; without it they use the
prescribed/default opening schedule. `propose_commitment` returns both `u`
and `tunnel_opening`, and `schedule_case` dispatches the paired controls.
Preparation provides a feasible schedule when accepted, not a global bound.
For selectable gates, the MILP proposal uses affine reservoir heads, secant
tunnel pressure losses and a first-order power approximation with frozen
efficiency. This lets it value alternative water supplies. The subsequent
nonlinear dispatch and physical replay use the original hydraulic equations.
Passive links between known reservoir/boundary heads derive their reference
flow from the pressure-loss law rather than from a station's nominal discharge.
Cases with prescribed gates retain their existing proposal model.

Closing gates can leave junction components without a reservoir or boundary
pressure reference. Independent integration fixes one numerical reference per
unobservable pressure component and retains every water-balance equation in a
rectangular Newton solve. It then selects a common pressure shift compatible
with the junction bounds. Replay compares such heads modulo that common shift;
head differences, operational bounds, flows and generation are still audited.
Components whose pressure is observable through active generation are not
assigned an independent reference.

```julia
prepared = schedule_case(case; proposal_time_limit=5.0,
    nlp_time_limit=15.0, max_refinements=0)
initial = prepared["accepted"] ? prepared["solution"] : nothing
result = solve(case; initial, time_limit=60.0, relative_gap=1e-3)
```

For independent reconstruction, pass the same matrix to
`dispatch_from_controls(case, u, generator_q, river_gate; tunnel_opening)` or
`simulate(case, generator_q, river_gate; tunnel_opening)`. These interfaces
check positions, prescriptions and interlocks rather than clipping controls.
The lower-level `forward_step` accepts one interval's opening vector and
integrates its physical equations; full chronological restrictions require the
case-level interfaces. Candidate validation and finer
replay check the actual gate decisions. `restart_case` stores the last actual
opening as the next initial position, so movement costs remain chronological.

`accepted` means the equations, operating rules and finer replay passed.
`global_certificate` additionally requires an enclosing SCIP bound within the
requested objective gap. Gate optimization does not change that distinction.
The certificate covers the declared interval model and listed gate positions,
not arbitrary continuous gate movement. A withheld bound remains unavailable
even when the reconstructed schedule matches a separately calculated optimum.

Run the [reservoir supply example](../examples/tunnel_gates.jl):

```sh
./scripts/julia.sh examples/tunnel_gates.jl
```

Two reservoirs supply one station at a required discharge. Positive prices
favor the higher supply; a negative-price interval favors the lower supply.
The optimized High–Low–High switch pays both gates' movement costs. Enumerating
the eight admissible supply schedules gives an independent check of its
objective, hydraulic heads and storage depletion.

## Network validation

The [network benchmark](../benchmark/tunnel-gates.jl) joins two six-reservoir
basins with a signed level-flow channel. It has six units, ten tunnels, two
cross-linked supply networks and twenty-three river reaches. Each basin selects
one of two intakes at different elevations. A minimum high-station release
requires generation even during negative-price hours. A three-position bypass
gate and delayed plant/tributary confluences complete the network.

The runs below used Julia 1.13.1 on an Apple M1, with one Julia and one BLAS
thread, five seconds per MILP proposal and twenty seconds per nonlinear solve.
A two-hour warmup excluded compilation. These are single-run scheduling
measurements, without global optimization; they establish feasibility and
model integration, not an optimality gap or a speed comparison with SHOP.

| Horizon | Scheduling | Gate changes | Objective | Global model variables | Global construction |
| --- | ---: | ---: | ---: | ---: | ---: |
| 48 hours | 33.185 s | 18 | 169,096.105 | 16,476 | 0.491 s |
| 168 hours | 76.526 s | 58 | 590,761.726 | 57,680 | 3.424 s |

All three joint proposals were feasible in both runs. Both schedules passed
fourfold chronological replay. Maximum storage difference was below
`9.7e-7` Mm³, and maximum water-balance residual was below `7.2e-13` Mm³.
Every global model variable received an independently checked feasible start;
the week-long model has 1,344 gate-position binaries before presolve.
Gate movement costs were 180 and 580, including initial position changes.
The [scalar records](../benchmark/results/tunnel-gates.json) retain the limits,
model sizes, costs and checks. Run the same cases with:

```sh
./scripts/julia.sh benchmark/tunnel-gates.jl 48 168
```

The unchanged retained-Vest Tokke–Vinje 24-hour input also passed native SCIP
optimization with its existing audited initial schedule and a thirty-second
allowance. Delivered objective remained `743,562.606`; reconstructed objective,
power, tunnel flow and storage matched the previous schedule exactly. Its
26,213 variables and 36,640 constraints are unchanged, with zero gate-decision
binaries. Fourfold replay and a twelve-hour restart passed. The remaining
7.927% objective gap is not a certificate; this run checks compatibility with
prescribed gates. See the [regression record](../benchmark/results/tunnel-gates-tokke.json)
and [dataset scope](../examples/tokke_vinje/README.md).
