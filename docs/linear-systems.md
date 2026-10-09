# Local linear systems

This study targets the nonlinear dispatch inside `schedule_case`. HiGHS still
proposes commitments, and independent hydraulic and river replay still accepts
controls. SCIP's global search is outside these measurements.

## Where the time goes

Ipopt repeatedly solves a sparse symmetric indefinite Newton system containing
the constraint Jacobian, Lagrangian curvature and barrier terms. The matrix's
sparsity, numerical conditioning and repeated right-hand sides matter more than
dense matrix multiplication alone.

The initial six-hour profile spent 4.103 of 4.786 seconds in the primal-dual
system solver. Constraint Jacobians took 0.033 seconds. The run reached the
1,000-iteration cap; this is a bottleneck diagnosis, not evidence of a local or
global optimum.

The pinned Ipopt 3.14.19 MUMPS adapter loops over multiple right-hand sides,
calling MUMPS separately for each column. Its numeric-factorization timer is
missing, so a reported factorization time of zero does not mean no work occurred.
SPRAL's adapter batches the right-hand sides and times numeric factorization.
Their factorization subtimers therefore cannot be compared directly.
[MUMPS adapter](https://github.com/coin-or/Ipopt/blob/releases/3.14.19/src/Algorithm/LinearSolvers/IpMumpsSolverInterface.cpp),
[SPRAL adapter](https://github.com/coin-or/Ipopt/blob/releases/3.14.19/src/Algorithm/LinearSolvers/IpSpralSolverInterface.cpp).

The turbine efficiency operator currently supplies values for two inputs.
Ipopt uses limited-memory curvature when the nonlinear evaluator cannot supply
a Hessian. Its default Sherman–Morrison implementation adds back-solves for the
low-rank updates. The alternative `extended` implementation places those updates
in a larger bordered system. It can avoid back-solves but increase fill.
[Low-rank solver](https://github.com/coin-or/Ipopt/blob/releases/3.14.19/src/Algorithm/IpLowRankAugSystemSolver.cpp),
[extended solver](https://github.com/coin-or/Ipopt/blob/releases/3.14.19/src/Algorithm/IpLowRankSSAugSystemSolver.cpp).

## Matched kernel experiments

The initial screen used the same first commitment proposal and frozen starting
values at two and six hours, two measured repetitions per configuration, and
30-second dispatch allowances. All 52 dispatches passed original checks and
finer replay. Values below are medians; they measure `optimize!`, excluding
construction and acceptance. Higher accepted objective is better.

| Configuration | 2 h seconds | 2 h objective | 6 h seconds | 6 h objective |
|---|---:|---:|---:|---:|
| Bundled MUMPS, default curvature | 0.447 | 71,025.043 | 4.868 | 217,703.587 |
| Extended curvature system | 1.077 | 71,068.915 | 2.872 | 217,705.844 |
| Extended, history 3 | 0.800 | 71,055.754 | 2.382 | 217,701.245 |
| Extended, history 12 | 0.543 | 71,062.061 | 4.252 | 217,693.598 |
| SPRAL | 1.218 | 71,061.682 | 9.100 | 217,702.985 |
| SPRAL, extended | 1.609 | 71,055.395 | 10.449 | 217,694.486 |
| SPRAL, METIS and Ruiz scaling | 0.706 | 71,055.608 | 4.921 | 217,696.594 |
| MUMPS, METIS ordering | 1.672 | 71,055.841 | 5.007 | 217,702.564 |
| MUMPS, pivot threshold 0.0001 | 0.725 | 71,067.112 | 4.840 | 217,702.552 |
| Exact table derivatives | 0.797 | 71,055.924 | 3.695 | 217,585.907 |
| Apple Accelerate, native 32-bit BLAS | 1.345 | 71,063.046 | 3.962 | 217,700.271 |
| Accelerate, extended | 1.100 | 71,062.312 | 2.959 | 217,700.636 |
| Extended, four BLAS threads | 1.084 | 71,068.915 | 2.885 | 217,705.844 |

The six-hour default spent 4.162 seconds in the primal-dual solver, including
3.001 seconds of back-solves. Extended curvature reduced these to 2.169 and
0.667 seconds. History 3 reduced total kernel time further, but a faster kernel
alone does not establish a faster or better complete schedule.

The exact-derivative prototype supplied the gradient and lower-triangular
Hessian of the selected table cell, preserving the original value function and
its existing directional derivative convention at the final discharge knot.
It passed 5,288 derivative and composed-power checks. It reduced back-solves and
enabled sparse Hessian evaluation, but lowered the six-hour objective and failed
to improve complete scheduling. The tables retain their original piecewise
smoothness; no smoothing was introduced.

## Complete scheduling

The full pipeline retains feasibility preparation, all three commitment
proposals, nonlinear dispatch, guarded reconstruction and independent acceptance.
Each configuration used two measured repetitions after two excluded two-hour
warmups, a two-second allowance per proposal, ten seconds per NLP, no additional
refinement and a 0.1 MW operating margin. These allowances are per invocation.

| Configuration | 2 h total seconds | 2 h objective | 6 h total seconds | 6 h objective |
|---|---:|---:|---:|---:|
| Production control | 4.848 | 72,199.463 | 14.124 | 217,703.587 |
| Extended | 4.405 | 72,154.397 | 8.840 | 217,701.593 |
| Extended, history 3 | 1.687 | 72,194.465 | 8.442 | 217,688.962 |
| Extended, history 12 | 4.139 | 72,145.984 | 13.906 | 217,694.775 |
| Exact derivatives | 4.723 | 72,193.252 | 15.804 | 217,674.954 |

History 3 saves about 65% and 40% on the short cases, with objectives about
0.007% lower. Extended history 6 loses 0.062% at two hours. The first synthetic
calls also include compilation; their timing is retained separately and cannot
be presented as warmed solver performance. Longer-case confirmation rejects history 3 as a universal default:

| Case | Production / history 3 seconds | Production objective | History 3 objective |
|---|---:|---:|---:|
| Ordinary 24 h | 49.505 / 54.493 | 743,677.334 | 743,301.565 |
| Seasonal 24 h | 39.792 / 36.101 | 453,194.141–465,269.764 | 453,194.141 |
| Synthetic turbine tables, warmed repeat | 1.473 / 1.791 | 9,790.216 | 9,790.216 |
| Synthetic distributed rivers, warmed repeat | 0.058 / 0.057 | 14,860.129 | 14,860.129 |

The ordinary 24-hour case is about 10% slower and 0.051% lower in objective.
The seasonal control found a materially better schedule in one of its two
repetitions; its range is shown instead of hiding that variation in a median.
The synthetic timings use the second, compilation-free repetition only and
support no statistical significance claim. All eight confirmation calls passed
physical and replay acceptance.

## Native right-hand-side batching

A separate experiment rebuilt Ipopt 3.14.19 twice with the same Apple Clang
options and existing MUMPS 5.8.1, METIS, ASL and libblastrampoline libraries.
Both builds add the missing factorization timer. Only the second changes
MUMPS's dense right-hand-side handling. The rebuilds disable the same optional
backends, use identical isolated Julia projects and a shared frozen OpenSHOP
source, and leave installed artifacts and global preferences untouched.

The [batch patch](../benchmark/patches/ipopt-3.14.19-batch-rhs.patch) replaces
one solve per column with one `NRHS`/`LRHS` solve of the existing column-major
block. The separate [timer patch](../benchmark/patches/ipopt-3.14.19-mumps-timer.patch)
is required in both controls for comparable native counters. These are patches
against Ipopt's EPL-licensed source, not part of OpenSHOP's production solver.

A small indefinite-matrix check reused one factorization for 1, 3 and 12 right-hand
sides, compared against repeated singleton solves, and checked a singleton after
a batch. All 12 checks passed: maximum scaled residual was 6.8e-17, and block
and singleton results agreed exactly. This checks the native block interface;
the hydropower comparisons exercise the patched Ipopt adapter.

Three measured dispatch repetitions followed the same excluded warmups:

| Case | Rebuilt control / batch seconds | Control / batch iterations | Control / batch objective |
|---|---:|---:|---:|
| 2 h | 1.634 / 0.296 | 1,000 / 237 | 71,056.326 / 71,051.832 |
| 6 h | 4.854 / 3.364 | 1,000 / 1,000 | 217,700.328 / 217,689.834 |

The six-hour kernel is 30.7% faster, with back-solve time falling from 2.989 to
1.523 seconds. The two-hour 81.9% reduction also includes earlier convergence;
it must not be attributed entirely to faster triangular solves. All 12
hydropower dispatches passed physical and replay acceptance.

The rebuilt control differs numerically from the bundled binary: its two-hour
run reaches 1,000 iterations rather than 277. Comparing only the patched rebuild
against the bundled binary would confound batching with compilation. Even the
matched rebuild changes its trajectory slightly under batching, so delivered
objectives remain part of the comparison.

The first complete batching campaign passed two-hour scheduling, then crashed
while unwinding a C++ exception during the six-hour case. The crash report showed
Julia's `libunwind` alongside the macOS C++ exception runtime. The custom linker's
log confirms it accidentally picked up that unwinder from the Julia library
search path. Both rebuilds were changed to the system unwinder, leaving source,
compiler settings and numerical libraries unchanged. The repeated complete
comparison then passed all eight calls:

| Case | Rebuilt control / batch total seconds | Control / batch objective |
|---|---:|---:|
| 2 h | 5.701 / 3.706 | 72,196.243 / 72,190.213 |
| 6 h | 18.304 / 13.140 | 217,702.364 / 217,701.635 |

Batching saves 35.0% and 28.2% against the matched rebuilds, with objective losses
of 0.0084% and 0.00034%. These are not reductions against the shipped binary:
that control's pipeline took 4.848 and 14.124 seconds, with objectives
72,199.463 and 217,703.587. Against it, the patched rebuild is approximately 24%
and 7% faster, with slightly lower revenue. The initial crash and the two sets
of binary fingerprints remain in the evidence.

## Production decision

OpenSHOP retains bundled Ipopt/MUMPS and its default curvature policy. The
portable configuration changes did not produce a consistent gain in speed and
revenue across the longer operating and synthetic cases. No alternate production
solver mode, native override or new dependency is introduced.

Batching is the strongest targeted linear-algebra result. It merits further
work on a normally distributed Ipopt build and longer-horizon, cross-platform
confirmation before adoption. Maintaining a private native solver build for
the measured short-case gains would add a substantial packaging obligation.
The patches remain separate benchmark artifacts.

The remaining local cost combines repeated sparse back-solves with stalled
Newton iterations. Faster linear algebra reduces each iteration's cost, while
conditioning and curvature determine how many iterations are needed and which
schedule is delivered. Replacing Julia's BLAS or raising thread counts did not
resolve that second issue in this study.

[Scalar evidence](../benchmark/results/linear-systems.json) records 64 accepted
fixed-commitment dispatches, 58 accepted completed scheduling calls, the failed
native campaign, derivative checks and matrix checks. It preserves input/start
hashes, compiler and library fingerprints, iterations, native counters and
physical/replay acceptance without distributing upstream cases or schedules.
The toolchain was Julia 1.13.1, JuMP 1.32.0, Ipopt.jl 1.16.0 and native Ipopt
3.14.19 on Apple Silicon, with one thread except the declared BLAS experiment.

## Research basis

Tasseff, Coffrin, Wächter and Laird compare Ipopt's sparse backends across NLP
classes. Their results show that ordering, solver choice and useful thread counts
depend on matrix structure and problem size; parallelism alone is insufficient.
[Exploring Benefits of Linear Solver Parallelism](https://arxiv.org/abs/1909.08104).

Hogg and Scott examine scaling, stability and fill. Barrier systems can become
ill-conditioned near a solution; scaling may improve robustness and reduce work,
but computing it also costs time. Their experiments motivate matched tests rather
than unconditional scaling or smaller pivot tolerances.
[On the effects of scaling on the performance of Ipopt](https://arxiv.org/abs/1301.7283).

Pacaud and colleagues reformulate KKT systems for sparse Cholesky factorization
on GPUs. Their 2026 revision studies conditioning as well as performance and
reports less robust edge cases on CUTEst. This is a possible route for much larger
models on NVIDIA hardware. It requires a solver/modeling-stack change and is not
a drop-in acceleration for this Apple Silicon experiment.
[Condensed Interior-Point Methods](https://arxiv.org/html/2405.14236v3).

## Reproducing a native comparison

Freeze a case, commitment and warm schedule once. For example:

```julia
using OpenSHOP, TOML
c = readcase("benchmark/cases/turbine-tables.json")
solve_verified(c; feasibility_only=true, time_limit=45, max_refinements=0)
reference = solve_verified(c; feasibility_only=true, time_limit=45,
                           max_refinements=0)
reference["accepted"] || error("No accepted reference")
warm = reference["solution"]
writejson("results/linear-fixture.json",
          Dict("case"=>case_dict(c), "u"=>warm["u"], "warm"=>warm))
open(io -> TOML.print(io, Dict()), "results/defaults.toml", "w")
open(io -> TOML.print(io, Dict("limited_memory_aug_solver"=>"extended")),
     "results/extended.toml", "w")
```

Run each configuration serially in a fresh process and output directory:

```sh
OMP_NUM_THREADS=1 ./scripts/julia.sh benchmark/linear-systems.jl \
  results/linear-fixture.json results/defaults.toml results/linear-default 3
OMP_NUM_THREADS=1 ./scripts/julia.sh benchmark/linear-systems.jl \
  results/linear-fixture.json results/extended.toml results/linear-extended 3
```

The runner excludes two warmups, records actual start values before Ipopt's
interior push, and saves native logs, iteration callbacks, original validation
and finer replay. Compare hashes, accepted objectives and total elapsed time.
`optimize_seconds` measures the dispatch kernel; acceptance and construction are
separate. Julia allocations exclude native solver storage.

SPRAL additionally requires `OMP_CANCELLATION=TRUE` and `OMP_PROC_BIND=TRUE`
before Julia starts. Ipopt's BLAS uses 32-bit integer indices; changing Julia's
64-bit BLAS alone does not change that backend. Thread experiments must account
for both wall and process CPU limits.
[Ipopt.jl linear solvers and BLAS](https://jump.dev/JuMP.jl/stable/packages/Ipopt/#Linear-Solvers).
