# Real-arithmetic model proofs

These Lean 4 proofs establish selected identities and counterexamples for the
mathematical model. They do not verify the Julia program, floating-point
arithmetic, JuMP, SCIP or Ipopt. Assumptions are explicit arguments to each
theorem. In particular, efficiency bounds, nonnegative volumes, causal response
functions and exact probability normalization are hypotheses, not inferred
properties of arbitrary input files.

Install [Lean through Elan](https://lean-lang.org/install/manual/), then run:

```sh
cd test/formal
lake exe cache get
lake build --wfail
lake env lean CheckAxioms.lean
```

`lean-toolchain` pins Lean 4.24.0. The Lake manifest pins mathlib and its
transitive dependencies. Neither is needed to install or use the Julia package.
CI builds these proofs and checks every declaration in the `OpenSHOP` namespace,
including generated auxiliaries. There are 46 explicitly stated theorems.
Only Lean's standard logical axioms (`propext`, `Classical.choice`, `Quot.sound`)
are allowed. Proof placeholders and additional axioms fail that check.

| Proof file | Statements checked | Julia equations |
| --- | --- | --- |
| `OpenSHOP/Hydraulics.lean` | Internal incidence cancellation; storage plus transit conservation; confluence allocation; signed tunnel dissipation; passive channel interpolation; pointwise generation/pumping work bounds | `hydraulic_topology.jl`, `river_sources.jl`, `tunnel_physics.jl`, `channel_physics.jl`, `curve_physics.jl` |
| `OpenSHOP/Hydraulics.lean` | Midpoint potential remainder and artificial cycle work; passive equalization overshoot; multiple roots with passive reference-dependent channels; pressure-shift invariance and bound translation | `forward.jl`, `hydraulic_residual.jl`, `pressure_references.jl` |
| `OpenSHOP/Interpolation.lean` | Hermite coefficients; SOS2 cubic correction and product bounds; affine cell reparameterization; signed tensor-product moment; corner partition; coefficient enclosure implication; implicit electrical-power inversion | `curves.jl`, `polynomial_tables.jl`, `tensor_tables.jl`, `power_hull.jl`, `power_bounds.jl` |
| `OpenSHOP/Commitment.lean` | Exact binary products for signed variables; positive-power head implication and its idle counterexample; binary transition rows; storage interval propagation; signed on-head bounds; conditional objective-gap implication | `algebra.jl`, `joint_states.jl`, `global_bounds.jl`, `power_hull.jl`, `solve.jl` |
| `OpenSHOP/Routing.lean` | Convex transfer fractions; quadratic released transfers; cohort inventory partition; causal rectangular-cohort arrivals; window additivity; constant-release refinement | `routing.jl`, `temporal.jl`, `distributed_routing.jl` |

This mapping is reviewed by humans; there is no automatic equivalence proof
between Lean definitions and Julia source. The interpolation proofs concern a
selected cell and its algebra, not the solver's SOS2 selection machinery. The
coefficient enclosure theorem assumes valid coefficient bounds; it does not
certify their floating-point computation. The routing proofs assume the supplied
CDF is causal and monotone; they do not verify the complete routing integrator.
The gap theorem assumes a genuine enclosing real bound and does not certify
SCIP output.

The [model validity guide](../../docs/model-validity.md) explains the physical
assumptions and numerical limits. `../model_consistency_tests.jl` contains the
midpoint energy and nonuniqueness examples and the pressure-reference regression.
`../zero_minimum_tests.jl` checks that both local and global optimization retain
an idle committed unit at negative permitted head.
