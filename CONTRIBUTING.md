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
