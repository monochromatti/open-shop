# Numerical design

OpenSHOP keeps named watercourse objects at the public interface and ordinary
arrays inside numerical calculations. Object order determines array rows; time
intervals determine columns. Storage has one additional column for the initial
state. Units are explicit: hours, Mm³, m³/s, metres and MW.

## Equations and algorithms

`HydraulicTopology` resolves names to integer connectivity. It contains no
operating policies, solver objects or cached decisions. Tunnel and generator
incidence distinguishes the water destination from the hydraulic head reference:
a plant can discharge into a river while using a reservoir as its outlet head.
Both optimization formulations use `constrain_water_balance!` for the same
storage and junction equations. Their nonlinear representations remain distinct:
Ipopt uses direct smooth expressions; SCIP uses the bounded polynomial model.

The chronological simulator is an independent equation evaluation.
`HydraulicStepData` describes one midpoint step with prescribed controls.
`hydraulic_residual!` writes four contiguous equation groups:

1. Reservoir storage change divided by interval duration, minus net inflow.
2. Zero-storage junction flow balance.
3. Tunnel head difference minus signed quadratic pressure loss.
4. Outlet release minus its level-dependent discharge law.

`HydraulicWorkspace` holds the head, net-flow and river arrays needed to evaluate
those equations. The residual does not construct dictionaries or solver objects.
`solve_hydraulics!` owns central differences, linear solves and backtracking.
Its buffers are allocated once per step, not once per residual evaluation.
Floating pressure components fix one head coordinate while retaining every
physical balance row. Incompatible water injections still fail convergence.

Public callers pass physical objects and control arrays. These internal types
hold the normalized data for a numerical step. Normalization checks the
Float64 values and preserves exact zero/nonzero control states; overflow and
control underflow raise an input error. Supplied sparse routing matrices retain
their sparse storage. Replay passes interval inflows and capacities as numerical
vectors while reusing fixed connectivity. It does not reconstruct physical
objects each time step. A `HydraulicStep` owns its arrays; a trajectory stores a
concrete vector of these steps. Later steps cannot overwrite saved results.

## Extending the model

Keep each physical law next to its numerical evaluation and interpolation rules.
Add the direct nonlinear and bounded global representations where needed, and
test against independent physical values or conservation identities. A new
water-balance term belongs in the common optimization balance and independent
replay, rather than a generic constraint registry.

Use concrete storage for repeated calculations and generic array arguments where
the operation permits them. The hydraulic residual can use an alternate scalar
type through its workspace; its derivative test uses automatic differentiation.
The production solver normalizes physical inputs to Float64; alternate scalar
types are supported at the residual interface.

Ordinary Julia loops are suitable numerical code. Broadcast when it makes an
array equation clearer or fuses elementwise work; use views or reusable buffers
when measurement justifies them. Avoid unconditional `@inbounds`, speculative
caches and whole-network temporary arrays in inner iterations. Preserve useful
physical variable names rather than abbreviating every intermediate or turning
every expression into a helper.

## Public use

Omit object collections that are absent from the system:

```julia
system = HydroSystem(reservoirs=[lake], plants=[station], generators=[unit],
    boundaries=[tailwater])
case = ScheduleCase(name="day", system=system, grid=collect(0.0:24.0), prices=prices)
```

Use `solve(case)` for the mixed-integer problem, `solve_case(case; u=commitment)`
for continuous dispatch with prescribed discrete states, and `simulate` for
physical replay. Named controls make their meaning visible:

```julia
trajectory = simulate(case; generator_q=q, gate=river_gates,
    tunnel_opening=tunnel_gates)
schedule = dispatch_from_controls(case; commitment=u, generator_q=q,
    gate=river_gates, tunnel_opening=tunnel_gates)
valid = schedule.diagnostics["validation"]["valid"]
```

The existing positional forms remain supported. The named forms use the same
implementation and validation. Inputs are read-only; selectable tunnel gates
still require an explicit schedule. Simulation convergence alone does not assert
that operating restrictions have been satisfied.

## Numerical results

`simulate` returns a `HydraulicTrajectory`; `forward_step` returns a
`HydraulicStep`. `dispatch_from_controls` returns a `Schedule`. Both `solve`
and fixed-commitment `solve_case` return a `SolveResult`, whose `solution` is a
schedule or `nothing`. Arrays retain case object order and the units above.

```julia
result = solve(case; time_limit=60.0)
if result.accepted
    schedule = result.solution
    storage = schedule.storage
    commitment = schedule.commitment
    writejson("schedule.json", schedule)
end
```

Numerical calculations use fields directly. Solver work counters and detailed
validation reports remain in `diagnostics` dictionaries outside the kernels.
Existing string-key reads, such as `schedule["V"]`, remain available. Records
have no dictionary mutation interface. Use `schedule_dict(schedule)` to obtain
an independent editable profile, then validate it before using it as input.

Case construction resolves direct river connections for native Julia objects
and JSON inputs. The supplied `HydroSystem` is unchanged; confluence nodes are
internal to the constructed case. A caller need not perform a dictionary round
trip or explicitly normalize these connections.

Global equation builders retain references to their primary variables and
auxiliary interpolation, tunnel-loss and outlet graphs. Warm starts populate
those references and check every variable and constraint. Variable names are
only diagnostic labels. The optimization builders share one economic objective;
physical reconstruction and validation still evaluate economics independently.

### Migration from 0.7

Case and physical schedule JSON retain their attributes and array orientation.
`writejson(path, result)` now uses one nested `solution` profile for both local
and global solve results. Readers of the former flattened local solve report
should read `report["solution"]`, or write `result.solution` to obtain the same
physical schedule profile as before. Omitted selectable gate decisions remain
an error; omitted optional accounting and malformed controls retain their
validation meaning. Native commitment values are checked without Float64
rounding, and nonzero controls may not underflow to zero during normalization.

## Design references

The [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/)
explains generic interfaces, mutation naming and consistent library conventions.
The [performance manual](https://docs.julialang.org/en/v1/manual/performance-tips/)
describes concrete storage, specialization barriers and selective preallocation.
[JuMP's larger-model guide](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/design_patterns_for_larger_models/)
shows how to keep model-building functions explicit while separating input data
from expression construction. Its [performance guidance](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/performance_tips/)
explains expression construction and avoiding unnecessary temporary expressions.

The resulting design favors a few substantial mathematical operations with
explicit inputs, concrete data and reusable numerical storage.

## Refactor measurements

The [recorded comparison](../benchmark/results/hydraulics-refactor.json) uses
nine warm replay calls on the six-reservoir mixed-source fixture, with irregular
control intervals, tunnel loops, plant/tunnel river outfalls, natural inflow and
confluences. The two variants use fixed and flow-dependent distributed delays.
The baseline is commit `2ffe0f7`; both runs use Julia 1.13.1 on an Apple M1 with
one Julia and one BLAS thread. Package loading and compilation are excluded.

| Replay | Before | After | Allocated before | Allocated after |
|---|---:|---:|---:|---:|
| Fixed delay | 9.66 ms | 1.71 ms | 10.76 MB | 2.23 MB |
| Distributed delay | 13.13 ms | 1.59 ms | 17.51 MB | 1.87 MB |

Storage, heads and power match exactly in these cases. The largest flow
difference is 3.6e-15 m³/s. Local and global variable/constraint counts are
unchanged. These short replay timings are not end-to-end optimization speedups;
the allocation reduction is the more reproducible result. The residual's
analytic-equation test also checks type inference and zero warm allocation.

The [numerical-data comparison](../benchmark/results/numerical-data-refactor.json)
uses the same network and nine warmed calls against commit `d93fa7a`, after the
proposal solver consolidation. Passing interval coefficients instead of copying
physical objects, reusing connectivity and storing typed steps gives:

| Replay | Before | After | Allocated before | Allocated after |
|---|---:|---:|---:|---:|
| Fixed delay | 1.79 ms | 0.81 ms | 2.23 MB | 0.87 MB |
| Distributed delay | 1.61 ms | 0.79 ms | 1.87 MB | 0.55 MB |

All compared physical arrays and objectives match exactly. Both model sizes
are unchanged, and complete SCIP starts satisfy every constraint within 2.6e-12.
These are replay measurements on one machine, with compilation excluded;
solver search remains a separate cost. Run `benchmark/hydraulics.jl` under each
revision's project to measure replay, allocations and model sizes.

The [Tokke–Vinje comparison](../benchmark/results/numerical-data-tokke.json)
uses identical imported cases, commitment and warm controls under each revision.
Each case has one excluded warmup and three measured `solve_case` calls, with
a 120-second allowance and the unchanged 1,000-iteration cap. Input fingerprints
and individual measurements are recorded; upstream data remain outside the
repository.

| Fixed-commitment dispatch | Objective, both revisions | Before | After |
|---|---:|---:|---:|
| 6 h | 219,719.686 | 4.876 s | 4.854 s |
| 24 h | 743,435.449 | 18.056 s | 18.011 s |

These medians measure optimization, excluding construction and validation.
Every measured call passes original-grid validation. Objectives and iteration
counts match exactly; both runs reach the iteration cap. The small timing
differences support no optimization speedup claim.

The global models retain 6,439 and 26,213 variables. Their complete starts pass
every constraint within 6e-12. Global solves with 30- and 60-second allowances
retain the supplied schedules, pass fourfold chronological replay and report
gaps of 6.2% and 7.1%. Those calls exercise the complete result and acceptance
path; they do not compare global search performance against the baseline.
