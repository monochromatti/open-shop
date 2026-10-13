# Contributing

Run the test suite before submitting a change:

```sh
./scripts/julia.sh -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Physical changes need an independent conservation or equation check, not only
a test of the JuMP model that implements them. Solver changes should preserve
the distinction between feasibility, finer replay and global bound quality.
Record the case, control grid, solver version, limits and actual elapsed time
when reporting a performance change.

Do not add external watercourse data without checking its redistribution terms.

## Julia and numerical code

Keep the equations visible, with documented units and array axes. Separate input
conversion and result serialization from repeated numerical evaluation. Use
concrete array storage and ordinary functions; introduce a type when it represents
physical data, an operator or reusable numerical workspace.

Keep related equations together. Extract a shared mathematical operation when it
removes duplicated knowledge, rather than splitting calculations into preparation
and processing layers. Document public methods with signatures, dimensions,
mutation behavior and the distinction between convergence and feasibility.

Refactors should preserve physical results and model semantics. Test conservation,
known solutions, invalid inputs and difficult domains through public operations.
Use kernel inference/allocation checks only where they protect a measured numerical
property. Report compilation separately from warm timings, and do not describe a
replay improvement as an end-to-end solver speedup.

The [paper](paper/paper.typ) states the model and solution method. Update it
when changing physical equations, interpolation or certificate scope. Keep
experiment outputs outside the repository.
