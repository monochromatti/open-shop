# Effective penalty settings

`penalty_settings(case)` exposes the effective operating and terminal penalty
equations before a solve. It returns `defaults`, `soft`, `release` and `reserve`.
It does not change the case. All costs use the same currency as `case.prices`;
no currency conversion or hidden large penalty is applied.

```julia
settings = penalty_settings(case)
println(settings.defaults)
for entry in settings.soft
    println(entry.rule.name, ": ", entry.unit, ", ", entry.coefficient_unit,
        ", factor=", entry.cost_factor, ", active=", entry.active)
end
assessment = soft_constraint_assessment(case, schedule)
println(assessment.slacks) # physical lower/upper deviations and each cost
println(objective_components(case, schedule)) # separate economic components
```

## Defaults and origins

| Input | OpenSHOP default | Consequence |
| --- | --- | --- |
| `soft_constraints` | empty | No additional targets or bands |
| Each soft allowance | `0.0` in the quantity's unit | A supplied limit remains hard |
| Each soft coefficient | `0.0` | Positive allowance requires an explicit positive cost |
| River `min_release` | `0.0 m³/s` | No additional entry minimum |
| River `release_penalty` | `0.0 currency/Mm³` | A supplied entry minimum is hard |
| Reserve `shortfall_cost` | `nothing` | A supplied reserve requirement is hard |

The returned `defaults` tuple records the numeric defaults. `soft` entries
retain the original typed `rule`, including its limits, bounded allowances,
coefficients and absolute assessment time/window. The side origins are
`absent`, `hard_limit` or `explicit_penalty`. A positive allowance can only have
an explicit positive coefficient. A zero allowance has no effective penalty,
even if its unused coefficient is positive. The normalized `SoftConstraint`
record cannot distinguish an omitted unused zero coefficient from an explicitly
supplied zero; those spellings describe the same equation.

Each `release` entry identifies a river and control interval. `minimum` and
`coefficient` are the effective interval values; their origins distinguish
`default_zero` from `operational_series`, including an explicitly supplied zero.
A series that starts inside an interval has `default_and_series` origin and
uses the same interval average as the equations. Valid scheduling grids include
interior operating changes. Signed delta-level channels have no entry-shortfall
policy and are absent from this collection.

Each `reserve` entry retains its original group `rule` and gives effective
interval up/down requirements. Requirement origins are `group_value` or
`operational_series` (or `default_and_series`). The coefficient origin is
`default_hard` for `nothing`, otherwise `explicit_penalty`; a non-null coefficient
must be strictly positive. `group_value` describes the normalized group value,
without claiming whether a default-valued keyword was supplied.

## Quantities, allowances and accrual

For either lower or upper side of a soft rule, the nonnegative cost is

```math
C=c\,\alpha\,\sigma,
\qquad \sigma^-=\max(0,L-x),\quad \sigma^+=\max(0,x-U).
```

The two sides add. `cost_factor` is `α` and `coefficient_unit` gives the unit of
`c`. The actual minimal physical deviation `σ` is distinct from its permitted
maximum and its monetary cost.

| Policy | Physical deviation (`unit`) | Coefficient unit | `cost_factor` | Maximum allowed |
| --- | --- | --- | --- | --- |
| Reservoir volume | Mm³ at one state vertex | currency/Mm³ | `1` | Original side allowance |
| Reservoir level | metres at one state vertex | currency/metre | `1` | Original side allowance |
| Power target/band | MW of original-window average | currency/MWh | window hours | Original side allowance |
| Discharge target/band | m³/s of original-window average | currency/Mm³ | `0.0036 × window hours` | Original side allowance |
| River entry minimum | m³/s below interval minimum | currency/Mm³ | `0.0036 × interval hours` | Whole minimum if coefficient is positive; otherwise zero |
| Reserve requirement | MW of up/down capacity deficit | currency/MWh | interval hours | Whole directional requirement if coefficient is positive; otherwise zero |

State costs accrue once at their declared absolute time. A terminal target is
a state observation at the original terminal time, not every restarted endpoint.
An inactive future observation costs zero. Output windows retain their original
endpoints and use the time-weighted average of the unit output, or the sum of
plant member outputs. Refinement creates no extra assessments. Separate windows
are necessary if deviations in different periods must not offset each other.
Restart drops completed assessments and rejects a restart inside an unfinished
output window. [Operating rules](operating-rules.md#explicit-penalty-limits)
give the shared optimization and independent assessment equations.

For example, an average discharge excess of `2 m³/s` over `1.5 hours` at
`1000 currency/Mm³` costs `1000 × 0.0036 × 1.5 × 2 = 10.8 currency`.
The allowed average rate and integrated `0.0108 Mm³` are different quantities.
`Schedule.shortfall_release` reports river deficits as integrated Mm³;
`ReserveAllocation.up_shortfall` and `down_shortfall` report MW.
`SoftConstraintSlack` reports the original physical unit and separate cost.
`objective_components` recomputes `release_penalty`, `soft_penalty` and
`reserve_penalty` from physical results, independently of reported cost fields.
Use `validate` separately to establish feasibility, including hard requirements
and deviations beyond allowances.

Conservation, constitutive equations, static safety/calibration domains,
ordinary `OperationalSeries` schedules, arrival minima and aggregate flow
requirements remain hard. There is no default reservoir or terminal slack,
and setting a soft target cannot override a hard static reservoir bound.

## Public SHOP reference and compatibility

SHOP's [global penalty command](https://docs.shop.sintef.energy/commands/penalty_cost.html)
documents defaults of 10,000,000 currency/Mm³ for reservoir/endpoint penalties,
50,000 currency/MWh for power limits, 100 currency/MWh for soft power,
100,000 currency/Mm³ for soft discharge and 1,000 currency/MWh for reserves.
Its currency is the configured base currency despite the table's NOK labels.
These defaults are not installed by OpenSHOP.

The [individual endpoint command](https://docs.shop.sintef.energy/commands/individual_endpenalty.html)
is deprecated and describes separate reservoir overrides and endpoint up/down
costs. It states currency/Mm³ for endpoint penalties but currency/MWh for ordinary
reservoir penalties, unlike the global command's reservoir unit. These pages
do not establish identical slack admission, time accrual or settings for every
licensed formulation. OpenSHOP therefore exposes its own explicit equations
and defaults rather than claiming automatic penalty compatibility. Matching
a coefficient alone does not establish equal feasible sets or objectives.
