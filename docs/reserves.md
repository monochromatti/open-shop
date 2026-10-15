# Spinning reserve capacity

OpenSHOP co-optimizes generation and upward/downward spinning reserve capacity.
`ReserveGroup` specifies eligible generators, MW requirements, capacity prices,
market depth and response time. Units remain committed in both reserve operating
points; downward reserve cannot stop a unit or bypass its minimum output.

```julia
reserve = ReserveGroup(name=:Spinning, generators=[:Unit],
    up_requirement=3.0, down_requirement=1.0,
    up_price=100.0, down_price=25.0, up_max=4.0, down_max=4.0,
    response_hours=0.25)
case = ScheduleCase(name="day", system=system, grid=grid, prices=prices,
    reserve_groups=[reserve])
result = solve(case)
capacity = result.solution.reserves
```

Prices and `shortfall_cost` are currency/(MW·hour). Capacity revenue equals
reserved MW × the scheduling interval's hours × its capacity price; it is not
payment for generated activation energy. `shortfall_cost=nothing` makes both
requirements hard. An explicit positive cost permits bounded MW shortfalls up
to the corresponding requirement. Independent validation checks their values
against the actual allocation. They never supply energy or water.

`up_max` and `down_max` bound each group's total offered capacity; omitted values
use the sum of eligible units' nameplate limits. Capacity prices apply to the
whole allocated amount, including capacity covering a requirement. There is no
separate reserve purchase, sale or energy-settlement model. Requirements,
prices and group maxima can use `OperationalSeries` at absolute hour knots.

## Physical operating points

For unit `g`, allocations across all groups share the same physical headroom:

\[
 R^+_{gt}=\sum_k r^+_{gkt},\qquad R^-_{gt}=\sum_k r^-_{gkt},
 \qquad p^+_{gt}=p_{gt}+R^+_{gt},\quad p^-_{gt}=p_{gt}-R^-_{gt}.
\]

Both operating points obey the supplied turbine and electrical power equations,
head-dependent turbine flow limits, efficiency limits and operative unit/plant
power and flow bounds. Their flows satisfy `q⁺ ≥ q` and `q⁻ ≤ q`. The plant's
net head at each point is

\[
 h^{\pm}_{pt}=H_{\mathrm{source},t}
 -\max(H_{\mathrm{target},t},h_{\mathrm{outlet\ floor}})
 -\mathrm{tailwater}_p\!\left(\sum_{g\in p}q^{\pm}_{gt}\right).
\]

The outlet floor is omitted when not supplied. Every unit in the station shares
this head. A unit that provides no reserve can adjust its operating-point flow
to hold its scheduled power when another unit changes plant tailwater.

The formulation uses feasible power/flow witnesses instead of assuming maximum
power occurs at maximum discharge. Analytic efficiencies, bilinear tables,
discharge PCHIP and electrical efficiency tables retain their existing exact
representations. `reserve_available=false` (or its operational series) removes
a unit from reserve allocation while still permitting ordinary generation.
Maintenance and an off commitment remove spinning capacity.

## Response and sustained activation

`response_hours` is time to reach the reserved operating point. Existing unit
and plant directional power/discharge rates limit that change, integrated from
the interval midpoint over the response window. For overlapping groups, their
shared unit/plant operating point must respond within the fastest applicable
response time. This is conservative when several products have different
response requirements. No startup or shutdown jump allowance is credited.

`activation_hours` is a separate sustained duration, zero by default. A positive
value requires explicit `Generator.reserve_up_energy` and
`reserve_down_energy` qualifications in MWh. The model enforces

\[
 \sum_k \tau_k r^+_{gkt}\le E^+_g,\qquad
 \sum_k \tau_k r^-_{gkt}\le E^-_g.
\]

These are supplied per-interval qualification budgets, not an internally
calculated reservoir energy inventory. They remain unchanged by restart because
this model contains no executed activation. Do not substitute `response_hours`
for sustained duration or treat a qualification budget as a hydraulic proof.

## Results, reconstruction and scope

`ReserveAllocation.up/down[unit,group,interval]` are MW. Its `up_q/down_q`
matrices are physical m³/s witnesses and `up_shortfall/down_shortfall` matrices
are MW by group and interval. JSON stores capacity as one unit-by-interval
matrix per group. Named group order follows `case.reserve_groups`.

Prescribed capacity can be reconstructed with ordinary physical controls:

```julia
capacity = ReserveAllocation(up=up, down=down, generator_q=q)
schedule = dispatch_from_controls(case; commitment=u, generator_q=q,
    gate=river_gates, reserves=capacity)
```

Forward reconstruction infers witness flows through the coupled plant power
laws; independent validation checks the stored witnesses without repairing
them. Finer replay repeats capacity controls and reconstructs witnesses at the
replayed node heads. Reconstruction retains supplied physical branches and uses
bounded Newton steps, with scalar bisection when a sign bracket is available.
It solves power residuals to 1e-8 MW where possible; boundary points use the same
1e-4 MW acceptance tolerance as independent power validation. Finer replay
retains reserve response limits even when ordinary scheduling ramps are checked
on their original intervals. A failed physical witness is a failed replay. Restart and
grid refinement retain reserve contract data and absolute operating times.

**This is a capacity model at scheduled hydraulic node heads.** Nominal tunnel,
river and storage balances contain only scheduled generation. The reserve
operating points re-evaluate plant tailwater but do not re-solve an activated
watercourse's tunnel losses, storage evolution, river routing or environmental
response. A global certificate concerns this declared capacity model, not
whole-watercourse activation deliverability.

Pumping reserves, nonspinning startup/stop delivery, FCR droop/frequency response,
forbidden-zone transitions, reserve minimum-lot binaries, market purchases and
stochastic or robust activation trajectories are outside this contract.

The SHOP [reserve tutorial](https://docs.shop.sintef.energy/examples/reserve_capacity/reserve_capacity.html),
[reserve-group reference](https://docs.shop.sintef.energy/objects/reserve_group/reserve_group.html),
[generator capabilities](https://docs.shop.sintef.energy/objects/gen_reserve_capability/gen_reserve_capability.html)
and [market reference](https://docs.shop.sintef.energy/objects/market/market.html)
describe SHOP's broader product-specific behavior. SINTEF's
[activation study](https://www.sintef.no/en/publications/publication/1793080/)
shows why hydraulic activation across a cascade requires a separate stochastic
or robust model. These product names do not confer equivalent compliance on
OpenSHOP's simpler spinning-capacity contract.
