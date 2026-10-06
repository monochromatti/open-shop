# Changelog

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
