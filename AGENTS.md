# Working on OpenSHOP

OpenSHOP is a scientific library for hydropower scheduling. Its code should
read like the physical equations and numerical methods it implements.

## Language and names

- Use plain language in explanations, documentation, comments and pull requests.
  State the physical or numerical meaning before introducing specialist terms.
  Explain unfamiliar terms and abbreviations on first use.
- Prefer established terms from hydraulic engineering, physics, numerical
  analysis and optimization. Name the quantity or calculation, not a software
  workflow around it. For example, use `validate_replay`, `constraint_residuals`,
  `release_grid`, `trial_solution` and `reserve_operating_points`.
- Describe water released over a time interval as an interval release. A time
  series of those releases is a release series. Avoid calling either a
  "cohort". Where rivers meet, say "river junction" or explain "confluence".
- Use validation or residual evaluation for checking equations and restrictions.
  Avoid auditing, certification or qualification language for ordinary numerical
  checks. Reserve a mathematical certificate for an actual stated bound or proof.
- Replace vague bookkeeping names such as `candidate` and `witness` with the
  specific mathematical object: solver solution, trial point, commitment
  proposal, refined-grid state or reserve operating point.
- Keep conventional mathematical names when they make the equations clearer.
  `V`, `Q`, `head`, `residual`, `Jacobian`, `incumbent` and solver diagnostics can
  be appropriate. Explain them to users rather than inventing unusual synonyms.
  Document units, array axes and time-grid conventions where they matter.
- Keep names consistent across numerical routines, public APIs, serialized
  results, examples, tests and the paper. Document the current API directly.
- The project has no external users yet. Make API changes directly and update
  their uses. Do not add migration guides, deprecated aliases or compatibility
  layers unless explicitly requested. Prefer one clear current interface.

## Numerical code

- Keep related equations together in substantive modules. Extract a function
  when it names a mathematical operation, a reusable numerical method or a
  public boundary. Do not fragment a derivation into many preparation helpers.
- Let physical dependencies and the numerical method determine control flow.
  Do not introduce managers, registries, candidate-ranking systems or decisions
  driven by diagnostic labels around the equations.
- Separate input/output and solver statistics from physical calculations.
  Validation evaluates the model; it must not silently change the model or
  relax its physical restrictions.
- Use idiomatic Julia arrays, views, in-place operations and structured linear
  algebra where they improve clarity and measured performance. Preserve the
  distinction between control intervals, hydraulic steps and transport cells.
- A terminology refactor must preserve equations, tolerances, solver settings,
  physical units and numerical behavior. Avoid combining it with algorithm changes.

## Verification and documentation

- Test observable behavior: conservation, physical laws, operating restrictions,
  returned solutions, numerical convergence and supported input/output contracts.
  Avoid assertions about private names, variable counts or exact code structure.
- Run the relevant Julia tests and the complete suite for changes spanning
  numerical modules. The normal package command is
  `./scripts/julia.sh -e 'using Pkg; Pkg.test()'`.
- If formal statements change, build `test/formal` with `lake build --wfail`
  and run `lake env lean CheckAxioms.lean`. State proof assumptions and limits.
- Update examples and documentation with public API changes. Keep the paper's
  terminology aligned with the implementation, while preserving its equations
  and the meaning of cited research.
- Keep long experiment outputs and performance reports outside the repository.
