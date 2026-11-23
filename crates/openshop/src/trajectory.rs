//! Physical trajectories, revenue and residuals. No solver or I/O operations.

use crate::ScheduleCase;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Schedule {
    /// Reservoir by grid vertex, m³.
    pub storage: Vec<Vec<f64>>,
    /// Unit by control interval, m³/s.
    pub unit_flow: Vec<Vec<f64>>,
    /// Unit by control interval, MW.
    pub unit_power: Vec<Vec<f64>>,
    /// Net working head or pumping lift, unit by control interval, m.
    /// An idle generating unit can have negative head if its declared range permits it.
    pub unit_head: Vec<Vec<f64>>,
    pub plant_flow: Vec<Vec<f64>>,
    pub plant_power: Vec<Vec<f64>>,
    pub tunnel_flow: Vec<Vec<f64>>,
    pub tunnel_opening: Vec<Vec<f64>>,
    pub gate_adjustment_cost: f64,
    /// River by control interval. Withdrawal excludes local river inflow/outfalls.
    pub river_withdrawal: Vec<Vec<f64>>,
    /// Reservoir outlet gate fraction, river by control interval; zero for signed channels and river junctions.
    pub river_opening: Vec<Vec<f64>>,
    pub river_release: Vec<Vec<f64>>,
    /// Downstream arrival volume, m³, by river and control interval.
    pub river_arrival_volume: Vec<Vec<f64>>,
    /// Water in transit, m³, by river and grid vertex.
    pub river_storage: Vec<Vec<f64>>,
    /// Node by control interval, mean head in metres.
    pub head: Vec<Vec<f64>>,
    pub commitment: Vec<Vec<bool>>,
    /// Named flow requirement by control interval, m³/s.
    pub observed_flow: Vec<Vec<f64>>,
    /// Previously counted water plus integrated observations, m³, at vertices.
    pub observed_volume: Vec<Vec<f64>>,
    pub objective: f64,
    pub solver_status: i32,
    pub iterations: usize,
    pub validation: Validation,
}

impl Schedule {
    /// Production minus pumping consumption, MW, in control-interval order.
    pub fn net_power(&self, case: &ScheduleCase) -> Result<Vec<f64>, crate::Error> {
        if self.unit_power.len() != case.system.units.len()
            || self
                .unit_power
                .iter()
                .any(|row| row.len() != case.prices.len())
        {
            return Err(crate::Error::Input(
                "unit power axes do not match the case".into(),
            ));
        }
        Ok((0..case.prices.len())
            .map(|t| {
                case.system
                    .units
                    .iter()
                    .enumerate()
                    .map(|(g, u)| u.mode.power_sign() * self.unit_power[g][t])
                    .sum()
            })
            .collect())
    }
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

/// Interval heads integrated over the represented piecewise-linear storage paths.
pub(crate) fn hydraulic_heads(
    case: &ScheduleCase,
    routing: &crate::river::Routing,
    storage: &[Vec<f64>],
    columns: &[f64],
    junction_head: Vec<Vec<f64>>,
) -> Vec<Vec<f64>> {
    case.system
        .reservoirs
        .iter()
        .enumerate()
        .map(|(r, reservoir)| {
            storage[r]
                .windows(2)
                .enumerate()
                .map(|(t, v)| {
                    let vertices: Vec<_> = routing.storage_shape[r][t]
                        .iter()
                        .map(|(f, c)| (*f, (1. - f) * v[0] + f * v[1] + c.value(columns)))
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
        .chain(
            case.system
                .boundaries
                .iter()
                .map(|b| vec![b.head; case.prices.len()]),
        )
        .collect()
}

pub(crate) fn evaluate(
    problem: &crate::SchedulingProblem<'_>,
    mut schedule: Schedule,
    junction_head: Vec<Vec<f64>>,
    tolerance: f64,
) -> Schedule {
    let (case, topology, routing, operating) = (
        problem.case,
        &problem.topology,
        &problem.routing,
        &problem.operating,
    );
    let s = &case.system;
    let mut columns: Vec<_> = schedule
        .unit_flow
        .iter()
        .flatten()
        .copied()
        .chain(schedule.tunnel_flow.iter().flatten().copied())
        .chain(
            routing
                .outlets
                .iter()
                .flat_map(|&r| schedule.river_withdrawal[r].iter().copied()),
        )
        .collect();
    routing.complete_columns(&mut columns);
    schedule.observed_flow = problem
        .observations
        .rates
        .iter()
        .map(|row| row.iter().map(|f| f.value(&columns)).collect())
        .collect();
    schedule.observed_volume = case
        .flow_requirements
        .iter()
        .zip(&schedule.observed_flow)
        .map(|(r, rates)| {
            let mut total = r.initial_volume;
            let mut volumes = vec![total];
            for (q, w) in rates.iter().zip(case.grid.windows(2)) {
                total += q * 3600. * (w[1] - w[0]);
                volumes.push(total);
            }
            volumes
        })
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
    schedule.head = hydraulic_heads(case, routing, &schedule.storage, &columns, junction_head);
    schedule.plant_flow = topology
        .plant_units
        .iter()
        .map(|members| {
            (0..case.prices.len())
                .map(|t| members.iter().map(|&g| schedule.unit_flow[g][t]).sum())
                .collect()
        })
        .collect();
    schedule.unit_head = s
        .units
        .iter()
        .enumerate()
        .map(|(g, u)| {
            let (source, target) = topology.unit_nodes[g];
            let p = topology.unit_plants[g];
            let plant = &s.plants[p];
            (0..case.prices.len())
                .map(|t| {
                    let outlet = schedule.head[target][t]
                        .max(plant.outlet_head_floor.unwrap_or(f64::NEG_INFINITY));
                    let tailwater = plant
                        .tailwater_curve
                        .as_ref()
                        .map_or(0., |c| c.extended_value(schedule.plant_flow[p][t]));
                    (schedule.head[source][t]
                        - outlet
                        - tailwater
                        - u.loss_factor * schedule.unit_flow[g][t].powi(2))
                        * u.mode.power_sign()
                })
                .collect()
        })
        .collect();
    schedule.unit_power = s
        .units
        .iter()
        .enumerate()
        .map(|(g, u)| {
            schedule.unit_flow[g]
                .iter()
                .zip(&schedule.unit_head[g])
                .map(|(&q, &h)| u.trial_power(q, h))
                .collect()
        })
        .collect();
    schedule.plant_power = topology
        .plant_units
        .iter()
        .map(|members| {
            (0..case.prices.len())
                .map(|t| members.iter().map(|&g| schedule.unit_power[g][t]).sum())
                .collect()
        })
        .collect();
    schedule.gate_adjustment_cost = gate_adjustment_cost(&s.tunnels, &schedule.tunnel_opening);
    schedule.objective = objective(case, &schedule, operating);
    schedule.validation = validate(problem, &schedule, tolerance);
    schedule
}

fn objective(
    case: &ScheduleCase,
    schedule: &Schedule,
    operating: &crate::operating_data::OperatingData,
) -> f64 {
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
                    price
                        * u.mode.power_sign()
                        * schedule.unit_power[i][t]
                        * (case.grid[t + 1] - case.grid[t])
                        - if start {
                            operating.units[i][t].startup_cost
                        } else {
                            0.
                        }
                        - if stop {
                            operating.units[i][t].shutdown_cost
                        } else {
                            0.
                        }
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    revenue - schedule.gate_adjustment_cost
        + case
            .system
            .rivers
            .iter()
            .enumerate()
            .map(|(r, river)| {
                river.water_value
                    * (schedule.river_storage[r][case.prices.len()] - schedule.river_storage[r][0])
            })
            .sum::<f64>()
        + case
            .system
            .reservoirs
            .iter()
            .enumerate()
            .map(|(r, res)| {
                res.water_value * (schedule.storage[r][case.prices.len()] - res.initial_storage)
            })
            .sum::<f64>()
}

fn validate(
    problem: &crate::SchedulingProblem<'_>,
    schedule: &Schedule,
    tolerance: f64,
) -> Validation {
    let (case, topology, routing, operating) = (
        problem.case,
        &problem.topology,
        &problem.routing,
        &problem.operating,
    );
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
    let mut columns: Vec<_> = schedule
        .unit_flow
        .iter()
        .flatten()
        .copied()
        .chain(schedule.tunnel_flow.iter().flatten().copied())
        .chain(
            routing
                .outlets
                .iter()
                .flat_map(|&r| schedule.river_withdrawal[r].iter().copied()),
        )
        .collect();
    routing.complete_columns(&mut columns);
    result.operating_limits = result
        .operating_limits
        .max(problem.observations.residual(&columns))
        .max(
            problem
                .interlocks
                .residual(&schedule.commitment, &schedule.tunnel_opening),
        )
        .max(
            problem
                .level_windows
                .residual(case, &schedule.storage, &columns),
        );
    for (r, series) in routing.releases.iter().enumerate() {
        for (k, rate) in series.rates.iter().enumerate() {
            let period = case
                .grid
                .partition_point(|time| *time <= series.grid[k])
                .saturating_sub(1)
                .min(case.prices.len() - 1);
            result.operating_limits = result.operating_limits.max(violation(
                rate.value(&columns),
                if case.system.rivers[r].law.is_delta() {
                    -operating.river_capacity[r][period]
                } else {
                    0.
                },
                operating.river_capacity[r][period],
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
    for &e in &routing.outlets {
        let river = &case.system.rivers[e];
        let source = case
            .system
            .reservoirs
            .iter()
            .position(|r| r.name == river.source)
            .unwrap();
        for t in 0..case.prices.len() {
            let h = schedule.head[source][t];
            let gate = schedule.river_opening[e][t];
            let expected = match &river.law {
                crate::RiverLaw::Delta { curve, reference } => {
                    let target = case
                        .system
                        .reservoirs
                        .iter()
                        .position(|r| r.name == river.target)
                        .unwrap();
                    let other = schedule.head[target][t];
                    curve.extended_value(
                        h - other,
                        match reference {
                            crate::ReferenceLevel::Upstream => h,
                            crate::ReferenceLevel::Downstream => other,
                        },
                    )
                }
                _ => river.law.trial_outlet_flow(
                    h,
                    gate,
                    operating.river_capacity[e][t],
                    operating.river_inflow[e][t],
                ),
            };
            result.water_balance = result
                .water_balance
                .max((schedule.river_withdrawal[e][t] - expected).abs());
            if !river.law.is_delta() {
                let limits = operating.river_gate[e][t];
                result.operating_limits = result
                    .operating_limits
                    .max(violation(gate, limits.min, limits.max, 1.));
            }
        }
    }
    for (r, reservoir) in case.system.reservoirs.iter().enumerate() {
        let scale = (reservoir.storage.max - reservoir.storage.min).max(1.);
        result.operating_limits = result
            .operating_limits
            .max((schedule.storage[r][0] - reservoir.initial_storage).abs() / scale);
        for (t, &volume) in schedule.storage[r].iter().enumerate() {
            result.operating_limits = result.operating_limits.max(violation(
                volume,
                operating.reservoirs[r].storage[t].min,
                operating.reservoirs[r].storage[t].max,
                scale,
            ));
        }
        for t in 0..case.prices.len() {
            let volume = operating.reservoirs[r].volume_change[t];
            let level = operating.reservoirs[r].level_change[t];
            let dv = schedule.storage[r][t + 1] - schedule.storage[r][t];
            let dz = reservoir.level.trial_level(schedule.storage[r][t + 1])
                - reservoir.level.trial_level(schedule.storage[r][t]);
            result.operating_limits = result.operating_limits.max(violation(
                dv,
                volume.down.map_or(f64::NEG_INFINITY, |x| -x),
                volume.up.unwrap_or(f64::INFINITY),
                scale,
            ));
            result.operating_limits = result.operating_limits.max(violation(
                dz,
                level.down.map_or(f64::NEG_INFINITY, |x| -x),
                level.up.unwrap_or(f64::INFINITY),
                100.,
            ));
        }
        for (t, vertices) in routing.storage_shape[r].iter().enumerate() {
            for (f, c) in vertices {
                let volume = (1. - f) * schedule.storage[r][t]
                    + f * schedule.storage[r][t + 1]
                    + c.value(&columns);
                result.operating_limits = result.operating_limits.max(violation(
                    volume,
                    operating.reservoirs[r].interval_storage[t].min,
                    operating.reservoirs[r].interval_storage[t].max,
                    scale,
                ));
            }
        }
        for t in 0..case.prices.len() {
            let net = operating.reservoirs[r].inflow[t]
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
        for t in 0..case.prices.len() {
            let interval = &operating.tunnels[e][t];
            let Some(law) = interval
                .positions
                .iter()
                .find(|p| p.opening == schedule.tunnel_opening[e][t])
            else {
                result.operating_limits = f64::INFINITY;
                continue;
            };
            let (a, b) = topology.tunnel_nodes[e];
            let q = schedule.tunnel_flow[e][t];
            let error = if law.opening == 0. {
                // Flow is checked by its zero-capacity bound below; closed gates
                // impose no pressure law.
                0.
            } else {
                ((schedule.head[a][t] - schedule.head[b][t]) - law.resistance * q * q.abs()).abs()
            };
            result.pressure_balance = result.pressure_balance.max(error);
            let capacity = if law.opening == 0. {
                0.
            } else {
                interval.capacity
            };
            result.operating_limits = result.operating_limits.max(violation(
                q,
                if topology.tunnel_outfalls[e].is_some() {
                    0.
                } else {
                    -capacity
                },
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
            let limits = &operating.units[i][t];
            let previous_q = if t == 0 {
                u.operation.previous_flow()
            } else {
                Some(schedule.unit_flow[i][t - 1])
            };
            if let Some(previous) = previous_q {
                result.operating_limits = result.operating_limits.max(violation(
                    schedule.unit_flow[i][t] - previous,
                    limits.flow_change.down.map_or(f64::NEG_INFINITY, |x| -x),
                    limits.flow_change.up.unwrap_or(f64::INFINITY),
                    u.flow.max.max(1.),
                ));
            }
            let previous_p = if t == 0 {
                u.operation.previous_power()
            } else {
                Some(schedule.unit_power[i][t - 1])
            };
            if let Some(previous) = previous_p {
                let (start, stop) = u.operation.transitions(&schedule.commitment[i], t);
                result.operating_limits = result.operating_limits.max(violation(
                    schedule.unit_power[i][t] - previous,
                    limits.power_change.down.map_or(f64::NEG_INFINITY, |x| {
                        -x - if stop { limits.stopping_power_min } else { 0. }
                    }),
                    limits.power_change.up.map_or(f64::INFINITY, |x| {
                        x + if start { limits.power.min } else { 0. }
                    }),
                    u.power.max.max(1.),
                ));
            }
            if limits
                .forced
                .is_some_and(|state| state != schedule.commitment[i][t])
            {
                result.operating_limits = f64::INFINITY;
            }
            if let Some(p) = limits.fixed_power {
                result.operating_limits = result
                    .operating_limits
                    .max((schedule.unit_power[i][t] - p).abs() / u.power.max.max(1.));
            }
            if let Some(q) = limits.fixed_flow {
                result.operating_limits = result
                    .operating_limits
                    .max((schedule.unit_flow[i][t] - q).abs() / u.flow.max.max(1.));
            }
        }
        for t in 0..case.prices.len() {
            let limits = &operating.units[i][t];
            let on = schedule.commitment[i][t];
            let q = schedule.unit_flow[i][t];
            let p = schedule.unit_power[i][t];
            let h = schedule.unit_head[i][t];
            result.operating_limits = result.operating_limits.max(violation(
                q,
                if on { limits.flow.min } else { 0. },
                if on { limits.flow.max } else { 0. },
                u.flow.max.max(1.),
            ));
            result.operating_limits = result.operating_limits.max(violation(
                p,
                if on { limits.power.min } else { 0. },
                if on { limits.power.max } else { 0. },
                u.power.max.max(1.),
            ));
            if on && let Some(characteristic) = &u.pump_characteristic {
                result.operating_limits = result
                    .operating_limits
                    .max((q - characteristic.trial_discharge(h)).abs() / u.flow.max.max(1.));
            }
            if on
                && u.pump_characteristic.is_none()
                && let crate::Efficiency::Table(table) = &u.efficiency
            {
                let minimum = table.envelope_data(h, false).value;
                let maximum = table.envelope_data(h, true).value;
                result.operating_limits =
                    result
                        .operating_limits
                        .max(violation(q, minimum, maximum, u.flow.max.max(1.)));
            }
            if on {
                result.operating_limits = result
                    .operating_limits
                    .max((-q * h).max(0.) / (u.flow.max.max(1.) * 100.));
                let eta = crate::machine::hydraulic_efficiency(u, q, h);
                result.operating_limits =
                    result
                        .operating_limits
                        .max(violation(eta, u.minimum_efficiency, 1., 1.));
                if let Some(reference) = &topology.gross_heads[i] {
                    let gross = (schedule.head[reference.source][t]
                        - schedule.head[reference.target][t])
                        * u.mode.power_sign();
                    result.operating_limits = result.operating_limits.max(violation(
                        gross,
                        limits.gross_head_min.unwrap_or(f64::NEG_INFINITY),
                        limits.gross_head_max.unwrap_or(f64::INFINITY),
                        100.,
                    ));
                    if u.operation.transitions(&schedule.commitment[i], t).0
                        && let Some(lo) = limits.startup_gross_head_min
                    {
                        let at_vertex = |node: usize| {
                            if node < case.system.reservoirs.len() {
                                case.system.reservoirs[node]
                                    .level
                                    .trial_level(schedule.storage[node][t])
                            } else {
                                case.system.boundaries[node
                                    - case.system.reservoirs.len()
                                    - case.system.junctions.len()]
                                .head
                            }
                        };
                        let initial = (at_vertex(reference.source) - at_vertex(reference.target))
                            * u.mode.power_sign();
                        result.operating_limits = result.operating_limits.max(violation(
                            initial,
                            lo,
                            f64::INFINITY,
                            100.,
                        ));
                    }
                }
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
    for i in 0..case.system.plants.len() {
        let members = &topology.plant_units[i];
        let on: Vec<_> = (0..case.prices.len())
            .map(|t| members.iter().any(|&g| schedule.commitment[g][t]))
            .collect();
        for row in &topology.plant_commitment[i] {
            result.operating_limits = result.operating_limits.max(row.residual(&on).max(0.));
        }
        let q_scale = members
            .iter()
            .map(|&g| case.system.units[g].flow.max)
            .sum::<f64>()
            .max(1.);
        let p_scale = members
            .iter()
            .map(|&g| case.system.units[g].power.max)
            .sum::<f64>()
            .max(1.);
        for (t, &running) in on.iter().enumerate() {
            let limits = &operating.plants[i][t];
            if limits.forced.is_some_and(|state| state != running) {
                result.operating_limits = f64::INFINITY;
            }
            if let Some(p) = limits.fixed_power {
                result.operating_limits = result
                    .operating_limits
                    .max((schedule.plant_power[i][t] - p).abs() / p_scale);
            }
            if let Some(q) = limits.fixed_flow {
                result.operating_limits = result
                    .operating_limits
                    .max((schedule.plant_flow[i][t] - q).abs() / q_scale);
            }
            result.operating_limits = result.operating_limits.max(violation(
                schedule.plant_flow[i][t],
                if running { limits.flow.min } else { 0. },
                limits.flow.max,
                q_scale,
            ));
            result.operating_limits = result.operating_limits.max(violation(
                schedule.plant_power[i][t],
                if running { limits.power.min } else { 0. },
                limits.power.max,
                p_scale,
            ));
            let prior_q = if t == 0 {
                operating.plant_initial_flow[i]
            } else {
                Some(schedule.plant_flow[i][t - 1])
            };
            if let Some(prior) = prior_q {
                let change = schedule.plant_flow[i][t] - prior;
                result.operating_limits = result.operating_limits.max(violation(
                    change,
                    limits.flow_change.down.map_or(f64::NEG_INFINITY, |x| -x),
                    limits.flow_change.up.unwrap_or(f64::INFINITY),
                    q_scale,
                ));
            }
            let prior_p = if t == 0 {
                operating.plant_initial_power[i]
            } else {
                Some(schedule.plant_power[i][t - 1])
            };
            if let Some(prior) = prior_p {
                let jump_up: f64 = members
                    .iter()
                    .filter(|&&g| {
                        case.system.units[g]
                            .operation
                            .transitions(&schedule.commitment[g], t)
                            .0
                    })
                    .map(|&g| operating.units[g][t].power.min)
                    .sum();
                let jump_down: f64 = members
                    .iter()
                    .filter(|&&g| {
                        case.system.units[g]
                            .operation
                            .transitions(&schedule.commitment[g], t)
                            .1
                    })
                    .map(|&g| operating.units[g][t].stopping_power_min)
                    .sum();
                result.operating_limits = result.operating_limits.max(violation(
                    schedule.plant_power[i][t] - prior,
                    limits
                        .power_change
                        .down
                        .map_or(f64::NEG_INFINITY, |x| -x - jump_down),
                    limits
                        .power_change
                        .up
                        .map_or(f64::INFINITY, |x| x + jump_up),
                    p_scale,
                ));
            }
        }
    }
    for rule in &topology.exclusions {
        for t in 0..case.prices.len() {
            if rule
                .iter()
                .filter(|group| group.iter().any(|&g| schedule.commitment[g][t]))
                .count()
                > 1
            {
                result.operating_limits = f64::INFINITY;
            }
        }
    }
    let finite = schedule
        .storage
        .iter()
        .chain(&schedule.unit_flow)
        .chain(&schedule.unit_power)
        .chain(&schedule.unit_head)
        .chain(&schedule.head)
        .chain(&schedule.tunnel_flow)
        .chain(&schedule.tunnel_opening)
        .chain(&schedule.river_storage)
        .chain(&schedule.observed_flow)
        .chain(&schedule.observed_volume)
        .chain(&schedule.river_withdrawal)
        .chain(&schedule.river_arrival_volume)
        .chain(&schedule.river_release)
        .chain(&schedule.river_opening)
        .all(|row| row.iter().all(|v| v.is_finite()));
    result.feasible = finite
        && result.water_balance <= tolerance
        && result.pressure_balance <= 100. * tolerance
        && result.operating_limits <= tolerance
        && schedule.objective.is_finite();
    result
}

/// Total variation of physical gate positions, with the declared initial position.
pub(crate) fn gate_adjustment_cost(tunnels: &[crate::Tunnel], openings: &[Vec<f64>]) -> f64 {
    tunnels
        .iter()
        .zip(openings)
        .map(|(e, row)| {
            let mut previous = e.initial_opening.unwrap_or(e.opening);
            row.iter()
                .map(|&a| {
                    let cost = (a - previous).abs() * e.gate_adjustment_cost;
                    previous = a;
                    cost
                })
                .sum::<f64>()
        })
        .sum()
}
