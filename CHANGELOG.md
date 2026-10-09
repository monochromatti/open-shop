# Changelog

## Unreleased

- Added a fixed-commitment linear-system profiler with frozen starts, native timings, library fingerprints and independent replay acceptance.
- Benchmarked sparse backends, curvature history, exact derivatives, BLAS and MUMPS right-hand-side batching; retained production defaults after longer-case regressions.
- Shared guarded control reconstruction between local and global candidates, retaining corrected controls, recomputed objectives and raw solver diagnostics.
- Removed a duplicate hydraulic audit in fixed-commitment verification while preserving fresh original validation, finer replay and transport acceptance.
- Added compatible-grid checks for starts, iteration and stage timings, and a separate schedule-generation benchmark.
- Tested hydraulic start reconstruction and reference refresh, retaining the existing initialization after objective/time regressions.

## 0.4.5

- Reduced SCIP's ordinary OBBT iteration multiplier from 10 to 1, retaining its filtering, minimum allowance and generalized bounds.
- Selected the change after 78 accepted matched measurements. Time to the control's common final bound falls from 299.51 to 21.10 seconds at 2 hours and from 292.61 to 126.03 seconds at 6 hours, with unchanged objectives.
- Published the nearly unchanged 24-hour results and small ordinary-horizon regression. Operating global proofs remain unfinished.
- Tested certified joint power-cut LPs and a cumulative root budget, retaining the simpler production path and archiving the alternatives.

## 0.4.4

- Replaced independent power gates with shared on-state turbine coordinates, preserving binary branches and the exact nonlinear equations.
- Added a bounded root separator for globally certified power supports, with cached coefficients and reusable native-value buffers.
- Retained ordinary OBBT while disabling its extra bilinear projection LPs. Selected one wider cut budget after paired, repeated six-case comparisons.
- Reduced the mean five-minute operating gap from 7.183% to 6.494% with unchanged delivered objectives. The 2-hour time to 7% falls from 99.79 to 1.46 seconds; larger operating cases remain uncertified.
- Published bound trajectories, native statistics and the documented smaller-budget tradeoff. Added regression tests for shared lifts, native cut application and callback failure containment.

## 0.4.3

- Added gated power supports tied to existing head and discharge table coordinates, preserving the exact nonlinear equations and adding no binary decisions.
- Selected one construction after 90 accepted matched measurements. Repeated five-minute operating comparisons reduced the mean gap by 21.8% with unchanged objectives; small PCHIP certification takes about 0.16 seconds longer.
- Added typed gate records, complete start lifting and regression checks for polynomial extensions, operational limits, electrical efficiency maxima and singleton domains.
- Published the table-coordinate comparison and archived quadratic curvature and smaller alternatives without public solver modes.

## 0.4.2

- Strengthened exact nonlinear power with a joint discharge/head/efficiency envelope and safe conditional on-state bounds, without new binary decisions.
- Selected one fixed envelope after 114 accepted matched benchmark runs. Repeated five-minute operating comparisons reduced the mean gap by 13.5% with unchanged objectives; small-case certification takes longer.
- Added complete corner-weight start lifting and regression checks for analytic curves, PCHIP tables, operational limits and off-state continuation.
- Published the coupled-relaxation measurements and archived alternatives without adding public solver modes.
- Corrected earlier pointwise root-relaxation claims from a postsolve LP read. Native global bounds, accepted schedules and callback-free comparisons remain valid.

## 0.4.1

- Strengthened the original nonlinear generation equations with head-aware supporting power bounds.
- Added optional native SCIP JSON statistics for LP and bound-tightening costs.
- Preserved fixed tunnel direction choices when lifting numerically zero-flow starts.
- Added regression checks for operating limits, off-state heads, electrical efficiency and polynomial extensions.
- Published 134 performance measurements, including repeated Tokke–Vinje confirmation and the retained failed-start record.

## 0.4.0

- Settled on exact quadratic SOS2 tables for the production SCIP model after repeated seven-case confirmation.
- Removed the formulation keyword, Cartesian graphs and experimental network-bound and effective-flow implementations.
- Shared SOS2 coordinates across compatible linear curves and preserved fixed-state turbine domains.
- Added independently reconstructed probe audits, exact zero-flow delivery for off units, randomized table tests and empty-reservoir accounting.
- Published the frozen screening, refinement and confirmation measurements and one production benchmark command.

## 0.3.0

- Made the exact turbine-table graph the default with independent SOS2 axes, shared head weights and discharge PCHIP corrections.
- Retained the baseline graph for small or fixed-commitment workloads where it can be faster.
- Published matched hosted Tokke–Vinje benchmarks at 2, 6 and 24 hours, including the seasonal operating case.
- Added complete start lifting and adjacency audits for SOS2 constraints, plus paired table-graph benchmarks.

## 0.2.0

- Added aggregate minimum-flow observations to dispatch, commitment proposals, replay, serialization and restart.
- Imported the current Tokke–Vinje seasonal minimum-flow and storage rules, with dated normal and seasonal benchmark inputs.
- Added experimental conservative network domains, proven tunnel directions, table-cell pruning and exact turbine cell bounds.
- Added an experimental table graph with commitment coupling, scaled cell coordinates and direct affine segments.
- Added native SCIP diagnostics, conservative certificate checks and paired benchmarks with frozen cases and shared audited starts.

## 0.1.0

- Added named-object hydropower cases with reservoirs, shared plants, signed tunnels and delayed river networks.
- Added operational series, outages, commitment history, startup and shutdown costs and release-shortfall accounting.
- Added analytic and tabulated hydraulic and efficiency relationships, including discharge PCHIP interpolation and explicit head extrapolation.
- Added native SCIP global scheduling and numerical bounds, Ipopt fixed-commitment dispatch and independent physical replay.
- Added JSON case serialization, restart support, core tests and a pinned Nix development environment.
