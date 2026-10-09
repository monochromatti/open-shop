# Ipopt MUMPS benchmark patches

These patches target Ipopt 3.14.19's EPL-licensed source. See the
[upstream license](https://github.com/coin-or/Ipopt/blob/releases/3.14.19/LICENSE)
and the [linear-system study](../../docs/linear-systems.md).
They are not applied by OpenSHOP or its Nix development shell.

Apply the timer patch to both source trees; apply batching only to the second:

```sh
patch -p1 < ipopt-3.14.19-mumps-timer.patch
patch -p1 < ipopt-3.14.19-batch-rhs.patch
```

Build both with the same compiler, options and native dependencies. Select the
rebuilt libraries through isolated Ipopt_jll preferences, including the AMPL
interface; verify that loading SCIP does not load a second Ipopt library.
The macOS experiment also required consistent C++ exception/unwind libraries;
the study records the initial mixed-runtime crash and the corrected comparison.
Use the same frozen fixture with `benchmark/linear-systems.jl`. Compare
accepted objectives and complete scheduling as well as kernel time. The study's
scalar evidence records the source, patch, dependency and binary checksums.
