# Changelog

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
