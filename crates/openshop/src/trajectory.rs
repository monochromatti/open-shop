//! Physical trajectories, revenue and residuals. No solver or I/O operations.

use crate::{ScheduleCase, model::Topology, river::Routing};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Schedule {
    /// Reservoir by grid vertex, m³.
    pub storage: Vec<Vec<f64>>,
    /// Unit by control interval, m³/s.
    pub unit_flow: Vec<Vec<f64>>,
    /// Unit by control interval, MW.
    pub unit_power: Vec<Vec<f64>>,
    pub tunnel_flow: Vec<Vec<f64>>,
    /// River by control interval. Withdrawal excludes local river inflow/outfalls.
    pub river_withdrawal: Vec<Vec<f64>>,
    pub river_release: Vec<Vec<f64>>,
    /// Downstream arrival volume, m³, by river and control interval.
    pub river_arrival_volume: Vec<Vec<f64>>,
    /// Water in transit, m³, by river and grid vertex.
    pub river_storage: Vec<Vec<f64>>,
    /// Node by control interval, mean head in metres.
    pub head: Vec<Vec<f64>>,
    pub commitment: Vec<Vec<bool>>,
    pub objective: f64,
    pub solver_status: i32,
    pub iterations: usize,
    pub validation: Validation,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Validation {
    pub feasible: bool,
    /// Maximum unadjusted water-balance error, m³/s.
    pub water_balance: f64,
    /// Maximum pressure-law error, metres.
    pub pressure_balance: f64,
    /// Maximum bound violation divided by the quantity's scale.
    pub operating_limits: f64,
}

pub(crate) fn evaluate(
    case: &ScheduleCase,
    topology: &Topology,
    routing: &Routing,
    mut schedule: Schedule,
    junction_head: Vec<Vec<f64>>,
    tolerance: f64,
) -> Schedule {
    let s = &case.system;
    let columns: Vec<_> = schedule
        .unit_flow
        .iter()
        .flatten()
        .copied()
        .chain(
            routing
                .outlets
                .iter()
                .flat_map(|&r| schedule.river_withdrawal[r].iter().copied()),
        )
        .collect();
    schedule.river_release = routing
        .releases
        .iter()
        .map(|series| {
            case.grid
                .windows(2)
                .map(|w| series.integral(w[0], w[1]).value(&columns) / (w[1] - w[0]))
                .collect()
        })
        .collect();
    schedule.river_arrival_volume = routing
        .arrivals
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(t, rate)| 3600. * (case.grid[t + 1] - case.grid[t]) * rate.value(&columns))
                .collect()
        })
        .collect();
    schedule.river_storage = routing
        .inventory
        .iter()
        .map(|row| row.iter().map(|v| v.value(&columns)).collect())
        .collect();
    schedule.head = s
        .reservoirs
        .iter()
        .enumerate()
        .map(|(r, reservoir)| {
            schedule.storage[r]
                .windows(2)
                .enumerate()
                .map(|(t, v)| {
                    let vertices: Vec<_> = routing.storage_shape[r][t]
                        .iter()
                        .map(|(f, c)| (*f, (1. - f) * v[0] + f * v[1] + c.value(&columns)))
                        .collect();
                    vertices
                        .windows(2)
                        .map(|w| {
                            (w[1].0 - w[0].0) * reservoir.level.trial_mean_level(w[0].1, w[1].1)
                        })
                        .sum()
                })
                .collect()
        })
        .chain(junction_head)
        .chain(s.boundaries.iter().map(|b| vec![b.head; case.prices.len()]))
        .collect();
    schedule.unit_power = s
        .units
        .iter()
        .enumerate()
        .map(|(i, u)| {
            let (source, target) = topology.unit_nodes[i];
            schedule.unit_flow[i]
                .iter()
                .enumerate()
                .map(|(t, &q)| {
                    u.trial_generation(
                        q,
                        schedule.head[source][t] - schedule.head[target][t] - u.loss_factor * q * q,
                    )
                })
                .collect()
        })
        .collect();
    schedule.objective = objective(case, &schedule);
    schedule.validation = validate(case, topology, routing, &schedule, tolerance);
    schedule
}

fn objective(case: &ScheduleCase, schedule: &Schedule) -> f64 {
    let revenue = case
        .system
        .units
        .iter()
        .enumerate()
        .map(|(i, u)| {
            case.prices
                .iter()
                .enumerate()
                .map(|(t, &price)| {
                    let (start, stop) = u.operation.transitions(&schedule.commitment[i], t);
                    price * schedule.unit_power[i][t] * (case.grid[t + 1] - case.grid[t])
                        - if start { u.operation.startup_cost } else { 0. }
                        - if stop { u.operation.shutdown_cost } else { 0. }
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    revenue
        + case
            .system
            .rivers
            .iter()
            .enumerate()
            .map(|(r, river)| river.water_value * schedule.river_storage[r][case.prices.len()])
            .sum::<f64>()
        + case
            .system
            .reservoirs
            .iter()
            .enumerate()
            .map(|(r, res)| res.water_value * schedule.storage[r][case.prices.len()])
            .sum::<f64>()
}

fn validate(
    case: &ScheduleCase,
    topology: &Topology,
    routing: &Routing,
    schedule: &Schedule,
    tolerance: f64,
) -> Validation {
    let mut result = Validation {
        feasible: true,
        ..Default::default()
    };
    let violation = |value: f64, lo: f64, hi: f64, scale: f64| {
        if value.is_finite() {
            (lo - value).max(value - hi).max(0.) / scale
        } else {
            f64::INFINITY
        }
    };
    let columns: Vec<_> = schedule
        .unit_flow
        .iter()
        .flatten()
        .copied()
        .chain(
            routing
                .outlets
                .iter()
                .flat_map(|&r| schedule.river_withdrawal[r].iter().copied()),
        )
        .collect();
    for (r, series) in routing.releases.iter().enumerate() {
        for rate in &series.rates {
            result.operating_limits = result.operating_limits.max(violation(
                rate.value(&columns),
                0.,
                case.system.rivers[r].capacity,
                case.system.rivers[r].capacity.max(1.),
            ));
        }
        for t in 0..case.prices.len() {
            let seconds = 3600. * (case.grid[t + 1] - case.grid[t]);
            let balance = (schedule.river_storage[r][t + 1] - schedule.river_storage[r][t])
                / seconds
                - schedule.river_release[r][t]
                + schedule.river_arrival_volume[r][t] / seconds;
            result.water_balance = result.water_balance.max(if balance.is_finite() {
                balance.abs()
            } else {
                f64::INFINITY
            });
        }
    }
    for (r, reservoir) in case.system.reservoirs.iter().enumerate() {
        let scale = (reservoir.storage.max - reservoir.storage.min).max(1.);
        result.operating_limits = result
            .operating_limits
            .max((schedule.storage[r][0] - reservoir.initial_storage).abs() / scale);
        for &volume in &schedule.storage[r] {
            result.operating_limits = result.operating_limits.max(violation(
                volume,
                reservoir.storage.min,
                reservoir.storage.max,
                scale,
            ));
        }
        for (t, vertices) in routing.storage_shape[r].iter().enumerate() {
            for (f, c) in vertices {
                let volume = (1. - f) * schedule.storage[r][t]
                    + f * schedule.storage[r][t + 1]
                    + c.value(&columns);
                result.operating_limits = result.operating_limits.max(violation(
                    volume,
                    reservoir.storage.min,
                    reservoir.storage.max,
                    scale,
                ));
            }
        }
        for t in 0..case.prices.len() {
            let net = case.inflow(r, t)
                + routing.reservoir_net[r][t].value(&columns)
                + topology.incidence[r]
                    .iter()
                    .map(|&(u, sign)| sign * schedule.unit_flow[u][t])
                    .sum::<f64>()
                + topology.tunnel_incidence[r]
                    .iter()
                    .map(|&(e, sign)| sign * schedule.tunnel_flow[e][t])
                    .sum::<f64>();
            let rate = (schedule.storage[r][t + 1] - schedule.storage[r][t])
                / (3600. * (case.grid[t + 1] - case.grid[t]));
            let error = (rate - net).abs();
            result.water_balance = result.water_balance.max(if error.is_finite() {
                error
            } else {
                f64::INFINITY
            });
        }
    }
    let nr = case.system.reservoirs.len();
    for (j, junction) in case.system.junctions.iter().enumerate() {
        for t in 0..case.prices.len() {
            let node = nr + j;
            let net = topology.incidence[node]
                .iter()
                .map(|&(u, sign)| sign * schedule.unit_flow[u][t])
                .sum::<f64>()
                + topology.tunnel_incidence[node]
                    .iter()
                    .map(|&(e, sign)| sign * schedule.tunnel_flow[e][t])
                    .sum::<f64>();
            result.water_balance = result.water_balance.max(net.abs());
            result.operating_limits = result.operating_limits.max(violation(
                schedule.head[node][t],
                junction.head.min,
                junction.head.max,
                junction.head.max.abs().max(1.),
            ));
        }
    }
    for (e, tunnel) in case.system.tunnels.iter().enumerate() {
        let law = &topology.tunnel_laws[e];
        for t in 0..case.prices.len() {
            let (a, b) = topology.tunnel_nodes[e];
            let q = schedule.tunnel_flow[e][t];
            let error = if tunnel.opening == 0. {
                // Flow is checked by its zero-capacity bound below; closed gates
                // impose no pressure law.
                0.
            } else {
                ((schedule.head[a][t] - schedule.head[b][t]) * law.head_factor
                    - law.resistance * q * q.abs())
                .abs()
                    / law.head_factor
            };
            result.pressure_balance = result.pressure_balance.max(error);
            let capacity = if tunnel.opening == 0. {
                0.
            } else {
                tunnel.capacity
            };
            result.operating_limits = result.operating_limits.max(violation(
                q,
                -capacity,
                capacity,
                tunnel.capacity.max(1.),
            ));
        }
    }
    for (i, u) in case.system.units.iter().enumerate() {
        for row in &topology.unit_commitment[i] {
            result.operating_limits = result
                .operating_limits
                .max(row.residual(&schedule.commitment[i]).max(0.));
        }
        for t in 0..case.prices.len() {
            let previous = if t == 0 {
                u.operation.previous_flow()
            } else {
                Some(schedule.unit_flow[i][t - 1])
            };
            if let Some(previous) = previous {
                let change = schedule.unit_flow[i][t] - previous;
                let dt = u.operation.ramp_hours(&case.grid, t);
                let lo = u
                    .operation
                    .flow_ramp
                    .down
                    .map_or(f64::NEG_INFINITY, |v| -v * dt);
                let hi = u.operation.flow_ramp.up.map_or(f64::INFINITY, |v| v * dt);
                result.operating_limits =
                    result
                        .operating_limits
                        .max(violation(change, lo, hi, u.flow.max.max(1.)));
            }
        }
        for t in 0..case.prices.len() {
            let previous = if t == 0 {
                u.operation.previous_power()
            } else {
                Some(schedule.unit_power[i][t - 1])
            };
            if let Some(previous) = previous {
                let (start, stop) = u.operation.transitions(&schedule.commitment[i], t);
                let dt = u.operation.ramp_hours(&case.grid, t);
                let lo = u.operation.power_ramp.down.map_or(f64::NEG_INFINITY, |v| {
                    -v * dt - if stop { u.power.min } else { 0. }
                });
                let hi = u.operation.power_ramp.up.map_or(f64::INFINITY, |v| {
                    v * dt + if start { u.power.min } else { 0. }
                });
                let change = schedule.unit_power[i][t] - previous;
                result.operating_limits =
                    result
                        .operating_limits
                        .max(violation(change, lo, hi, u.power.max.max(1.)));
            }
        }
        let (a, b) = topology.unit_nodes[i];
        for t in 0..case.prices.len() {
            let on = schedule.commitment[i][t];
            let q = schedule.unit_flow[i][t];
            let p = schedule.unit_power[i][t];
            let h = schedule.head[a][t] - schedule.head[b][t] - u.loss_factor * q * q;
            result.operating_limits = result.operating_limits.max(violation(
                q,
                if on { u.flow.min } else { 0. },
                if on { u.flow.max } else { 0. },
                u.flow.max.max(1.),
            ));
            result.operating_limits = result.operating_limits.max(violation(
                p,
                if on { u.power.min } else { 0. },
                if on { u.power.max } else { 0. },
                u.power.max.max(1.),
            ));
            if on && let crate::Efficiency::Table(table) = &u.efficiency {
                let minimum = table.envelope_data(h, false).value;
                let maximum = table.envelope_data(h, true).value;
                result.operating_limits =
                    result
                        .operating_limits
                        .max(violation(q, minimum, maximum, u.flow.max.max(1.)));
            }
            if on {
                result.operating_limits = result.operating_limits.max(violation(
                    h,
                    u.head.min,
                    u.head.max,
                    u.head.max.max(1.),
                ));
            }
        }
    }
    result.feasible = result.water_balance <= tolerance
        && result.pressure_balance <= 100. * tolerance
        && result.operating_limits <= tolerance
        && schedule.objective.is_finite();
    result
}
