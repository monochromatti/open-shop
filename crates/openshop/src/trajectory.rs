//! Physical trajectories, revenue and residuals. No solver or I/O operations.

use crate::ScheduleCase;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Schedule {
    /// Hydraulic time edges, hours; state and power arrays use these cells.
    pub grid: Vec<f64>,
    /// River mixing time edges, hours; empty for exact fixed-delay routing.
    pub routing_grid: Vec<f64>,
    /// Original decision times, hours; controls remain held on these intervals.
    pub operating_grid: Vec<f64>,
    /// Reservoir by grid vertex, m³.
    pub storage: Vec<Vec<f64>>,
    /// Unit by hydraulic cell, m³/s.
    pub unit_flow: Vec<Vec<f64>>,
    /// Unit by hydraulic cell, MW.
    pub unit_power: Vec<Vec<f64>>,
    /// Net working head or pumping lift, unit by hydraulic cell, m.
    /// An idle generating unit can have negative head if its declared range permits it.
    pub unit_head: Vec<Vec<f64>>,
    pub plant_flow: Vec<Vec<f64>>,
    pub plant_power: Vec<Vec<f64>>,
    pub tunnel_flow: Vec<Vec<f64>>,
    pub tunnel_opening: Vec<Vec<f64>>,
    pub gate_adjustment_cost: f64,
    /// River by hydraulic cell. Withdrawal excludes local river inflow/outfalls.
    pub river_withdrawal: Vec<Vec<f64>>,
    /// Reservoir outlet gate fraction, river by hydraulic cell; zero for signed channels and river junctions.
    pub river_opening: Vec<Vec<f64>>,
    pub river_release: Vec<Vec<f64>>,
    /// Downstream arrival volume, m³, by river and hydraulic cell.
    pub river_arrival_volume: Vec<Vec<f64>>,
    /// Water in transit, m³, by river and grid vertex.
    pub river_storage: Vec<Vec<f64>>,
    /// Node by hydraulic cell, mean head in metres.
    pub head: Vec<Vec<f64>>,
    pub commitment: Vec<Vec<bool>>,
    /// Named flow requirement by control interval, m³/s.
    pub observed_flow: Vec<Vec<f64>>,
    /// Previously counted water plus integrated observations, m³, at vertices.
    pub observed_volume: Vec<Vec<f64>>,
    /// Minimal deviations from the declared operating targets, in physical units.
    pub soft_slacks: Vec<crate::SoftSlack>,
    /// Total target-deviation cost, already subtracted from the objective.
    pub soft_penalty_cost: f64,
    /// Missing release volume, m³, by river and original operating interval.
    /// Only priced minimum releases are included here.
    pub release_shortfall_volume: Vec<Vec<f64>>,
    /// Shortfall cost subtracted from the objective.
    pub release_penalty_cost: f64,
    pub objective: f64,
    /// Ipopt termination code for an optimized result; None for pure reconstruction.
    pub solver_status: Option<i32>,
    pub iterations: usize,
    pub validation: Validation,
    /// Independent finer-grid differences for an optimized schedule. None means
    /// that time convergence has not been checked for these controls and states.
    pub time_error: Option<crate::TimeIntegrationError>,
}

impl Schedule {
    /// Production minus pumping consumption, MW, in control-interval order.
    pub fn net_power(&self, case: &ScheduleCase) -> Result<Vec<f64>, crate::Error> {
        let mesh = crate::time_discretization::TimeDiscretization::new(&case.grid, &self.grid)?;
        if self.operating_grid != case.grid {
            return Err(crate::Error::Input(
                "schedule operating grid does not match the case".into(),
            ));
        }
        if self.unit_power.len() != case.system.units.len()
            || self.unit_power.iter().any(|row| row.len() != mesh.dt.len())
        {
            return Err(crate::Error::Input(
                "unit power axes do not match the case".into(),
            ));
        }
        let power = mesh.means(&self.unit_power, 0.)?;
        Ok((0..case.prices.len())
            .map(|t| {
                case.system
                    .units
                    .iter()
                    .enumerate()
                    .map(|(g, u)| u.mode.power_sign() * power[g][t])
                    .sum()
            })
            .collect())
    }
}

/// Original-window quantities used by operating equations. Nonlinear physical
/// laws are evaluated on hydraulic cells, before these averages are taken.
struct OperatingValues {
    storage: Vec<Vec<f64>>,
    unit_flow: Vec<Vec<f64>>,
    unit_power: Vec<Vec<f64>>,
    plant_flow: Vec<Vec<f64>>,
    plant_power: Vec<Vec<f64>>,
    river_release: Vec<Vec<f64>>,
    head: Vec<Vec<f64>>,
    controls: crate::Controls,
}
impl OperatingValues {
    fn new(
        mesh: &crate::time_discretization::TimeDiscretization,
        s: &Schedule,
    ) -> Result<Self, crate::Error> {
        Ok(Self {
            storage: mesh.states(&s.storage)?,
            unit_flow: mesh.means(&s.unit_flow, 0.)?,
            unit_power: mesh.means(&s.unit_power, 0.)?,
            plant_flow: mesh.means(&s.plant_flow, 0.)?,
            plant_power: mesh.means(&s.plant_power, 0.)?,
            river_release: mesh.means(&s.river_release, 0.)?,
            head: mesh.means(&s.head, 0.)?,
            controls: s.controls()?,
        })
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
    mesh: &crate::time_discretization::TimeDiscretization,
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
                .map(|b| vec![b.head; mesh.dt.len()]),
        )
        .collect()
}

pub(crate) fn evaluate(
    problem: &crate::SchedulingProblem<'_>,
    mut schedule: Schedule,
    junction_head: Vec<Vec<f64>>,
    tolerance: f64,
) -> Result<Schedule, crate::Error> {
    schedule.time_error = None;
    let mesh = &problem.mesh;
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
            mesh.grid
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
                .map(|(t, rate)| 3600. * mesh.dt[t] * rate.value(&columns))
                .collect()
        })
        .collect();
    schedule.river_storage = routing
        .inventory
        .iter()
        .map(|row| row.iter().map(|v| v.value(&columns)).collect())
        .collect();
    schedule.head = hydraulic_heads(
        case,
        routing,
        mesh,
        &schedule.storage,
        &columns,
        junction_head,
    );
    schedule.plant_flow = topology
        .plant_units
        .iter()
        .map(|members| {
            (0..mesh.dt.len())
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
            (0..mesh.dt.len())
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
            (0..mesh.dt.len())
                .map(|t| members.iter().map(|&g| schedule.unit_power[g][t]).sum())
                .collect()
        })
        .collect();
    schedule.gate_adjustment_cost = gate_adjustment_cost(&s.tunnels, &schedule.tunnel_opening);
    let values = OperatingValues::new(mesh, &schedule)?;
    let (slacks, penalty_cost, soft_excess) =
        problem
            .soft_limits
            .assess(case, &values.storage, &values.unit_power, &values.unit_flow);
    schedule.soft_slacks = slacks;
    schedule.soft_penalty_cost = penalty_cost;
    schedule.release_shortfall_volume = operating
        .river_release_penalty
        .iter()
        .enumerate()
        .map(|(r, prices)| {
            prices
                .iter()
                .enumerate()
                .map(|(t, &price)| {
                    if price == 0. {
                        0.
                    } else {
                        (operating.river_min_release[r][t] - values.river_release[r][t]).max(0.)
                            * (3600. * mesh.duration[t])
                    }
                })
                .collect()
        })
        .collect();
    schedule.release_penalty_cost = schedule
        .release_shortfall_volume
        .iter()
        .zip(&operating.river_release_penalty)
        .flat_map(|(volumes, prices)| volumes.iter().zip(prices).map(|(v, p)| v * p))
        .sum();
    schedule.objective = objective(case, &schedule, &values, operating);
    schedule.validation = validate(problem, &schedule, &values, tolerance);
    schedule.validation.operating_limits = schedule.validation.operating_limits.max(soft_excess);
    schedule.validation.feasible &= soft_excess <= tolerance;
    Ok(schedule)
}

fn objective(
    case: &ScheduleCase,
    schedule: &Schedule,
    values: &OperatingValues,
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
                    let (start, stop) = u.operation.transitions(&values.controls.commitment[i], t);
                    price
                        * u.mode.power_sign()
                        * values.unit_power[i][t]
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
    let river_cost: f64 = case
        .system
        .rivers
        .iter()
        .enumerate()
        .map(|(e, r)| {
            (0..case.prices.len())
                .map(|t| {
                    let prior = if t == 0 {
                        r.previous_release().map(|(q, _)| q)
                    } else {
                        Some(values.river_release[e][t - 1])
                    };
                    prior.map_or(0., |q| {
                        let release = values.river_release[e][t];
                        let before = r.flow_on(q);
                        let after = r.flow_on(release);
                        operating.river_change_cost[e][t] * (release - q).abs()
                            + if after && !before {
                                operating.river_startup_cost[e][t]
                            } else {
                                0.
                            }
                            + if before && !after {
                                operating.river_shutdown_cost[e][t]
                            } else {
                                0.
                            }
                    })
                })
                .sum::<f64>()
        })
        .sum();
    revenue
        - schedule.gate_adjustment_cost
        - river_cost
        - schedule.soft_penalty_cost
        - schedule.release_penalty_cost
        + case
            .system
            .rivers
            .iter()
            .enumerate()
            .map(|(r, river)| {
                river.water_value
                    * (schedule.river_storage[r][schedule.grid.len() - 1]
                        - schedule.river_storage[r][0])
            })
            .sum::<f64>()
        + case
            .system
            .reservoirs
            .iter()
            .enumerate()
            .map(|(r, res)| {
                res.water_value * (values.storage[r][case.prices.len()] - res.initial_storage)
            })
            .sum::<f64>()
}

/// Express the change in stored-water value through conserved water fluxes.
/// On the balance equations this equals terminal value minus initial value.
/// It removes the large, physically redundant storage-range gradient from the
/// numerical objective; internal transfers carry differences in water value.
pub(crate) fn water_value_change(
    case: &ScheduleCase,
    topology: &crate::model::Topology,
    routing: &crate::river::Routing,
    operating: &crate::operating_data::OperatingData,
    mesh: &crate::time_discretization::TimeDiscretization,
) -> crate::river::LinearForm {
    let s = &case.system;
    let nt = mesh.dt.len();
    let mut value = crate::river::LinearForm::default();
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        for t in 0..nt {
            let mut net = routing.reservoir_net[r][t].clone();
            net.constant += operating.reservoirs[r].inflow[mesh.interval[t]];
            for &(g, sign) in &topology.incidence[r] {
                *net.coefficients.entry(g * nt + t).or_default() += sign;
            }
            for &(e, sign) in &topology.tunnel_incidence[r] {
                *net.coefficients
                    .entry((s.units.len() + e) * nt + t)
                    .or_default() += sign;
            }
            value.add_scaled(&net, 3600. * mesh.dt[t] * reservoir.water_value);
        }
    }
    for (r, river) in s.rivers.iter().enumerate() {
        value.add_scaled(&routing.inventory[r][nt], river.water_value);
        value.add_scaled(&routing.inventory[r][0], -river.water_value);
    }
    value
}

fn validate(
    problem: &crate::SchedulingProblem<'_>,
    schedule: &Schedule,
    values: &OperatingValues,
    tolerance: f64,
) -> Validation {
    let (case, topology, routing, operating) = (
        problem.case,
        &problem.topology,
        &problem.routing,
        &problem.operating,
    );
    let mesh = &problem.mesh;
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
                .residual(&values.controls.commitment, &values.controls.tunnel_opening),
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
        for t in 0..mesh.dt.len() {
            let seconds = 3600. * mesh.dt[t];
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
        for t in 0..mesh.dt.len() {
            let j = mesh.interval[t];
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
                    operating.river_capacity[e][j],
                    operating.river_inflow[e][j],
                ),
            };
            result.water_balance = result
                .water_balance
                .max((schedule.river_withdrawal[e][t] - expected).abs());
            if !river.law.is_delta() {
                let limits = operating.river_gate[e][j];
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
            let bounds = match mesh.vertices.binary_search(&t) {
                Ok(j) => operating.reservoirs[r].storage[j],
                Err(_) => operating.reservoirs[r].interval_storage[mesh.interval[t]],
            };
            result.operating_limits = result
                .operating_limits
                .max(violation(volume, bounds.min, bounds.max, scale));
        }
        for t in 0..case.prices.len() {
            let volume = operating.reservoirs[r].volume_change[t];
            let level = operating.reservoirs[r].level_change[t];
            let dv = values.storage[r][t + 1] - values.storage[r][t];
            let dz = reservoir.level.trial_level(values.storage[r][t + 1])
                - reservoir.level.trial_level(values.storage[r][t]);
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
                    operating.reservoirs[r].interval_storage[mesh.interval[t]].min,
                    operating.reservoirs[r].interval_storage[mesh.interval[t]].max,
                    scale,
                ));
            }
        }
        for t in 0..mesh.dt.len() {
            let j = mesh.interval[t];
            let net = operating.reservoirs[r].inflow[j]
                + routing.reservoir_net[r][t].value(&columns)
                + topology.incidence[r]
                    .iter()
                    .map(|&(u, sign)| sign * schedule.unit_flow[u][t])
                    .sum::<f64>()
                + topology.tunnel_incidence[r]
                    .iter()
                    .map(|&(e, sign)| sign * schedule.tunnel_flow[e][t])
                    .sum::<f64>();
            let rate = (schedule.storage[r][t + 1] - schedule.storage[r][t]) / (3600. * mesh.dt[t]);
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
        for t in 0..mesh.dt.len() {
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
        for t in 0..mesh.dt.len() {
            let j = mesh.interval[t];
            let interval = &operating.tunnels[e][j];
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
                .max(row.residual(&values.controls.commitment[i]).max(0.));
        }
        for t in 0..case.prices.len() {
            let limits = &operating.units[i][t];
            let previous_q = if t == 0 {
                u.operation.previous_flow()
            } else {
                Some(values.unit_flow[i][t - 1])
            };
            if let Some(previous) = previous_q {
                result.operating_limits = result.operating_limits.max(violation(
                    values.unit_flow[i][t] - previous,
                    limits.flow_change.down.map_or(f64::NEG_INFINITY, |x| -x),
                    limits.flow_change.up.unwrap_or(f64::INFINITY),
                    u.flow.max.max(1.),
                ));
            }
            let previous_p = if t == 0 {
                u.operation.previous_power()
            } else {
                Some(values.unit_power[i][t - 1])
            };
            if let Some(previous) = previous_p {
                let (start, stop) = u.operation.transitions(&values.controls.commitment[i], t);
                result.operating_limits = result.operating_limits.max(violation(
                    values.unit_power[i][t] - previous,
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
                .is_some_and(|state| state != values.controls.commitment[i][t])
            {
                result.operating_limits = f64::INFINITY;
            }
            if let Some(p) = limits.fixed_power {
                result.operating_limits = result
                    .operating_limits
                    .max((values.unit_power[i][t] - p).abs() / u.power.max.max(1.));
            }
            if let Some(q) = limits.fixed_flow {
                result.operating_limits = result
                    .operating_limits
                    .max((values.unit_flow[i][t] - q).abs() / u.flow.max.max(1.));
            }
        }
        for t in 0..mesh.dt.len() {
            let j = mesh.interval[t];
            let limits = &operating.units[i][j];
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
                    let gross = (values.head[reference.source][j]
                        - values.head[reference.target][j])
                        * u.mode.power_sign();
                    result.operating_limits = result.operating_limits.max(violation(
                        gross,
                        limits.gross_head_min.unwrap_or(f64::NEG_INFINITY),
                        limits.gross_head_max.unwrap_or(f64::INFINITY),
                        100.,
                    ));
                    if t == mesh.cells[j].start
                        && u.operation.transitions(&values.controls.commitment[i], j).0
                        && let Some(lo) = limits.startup_gross_head_min
                    {
                        let at_vertex = |node: usize| {
                            if node < case.system.reservoirs.len() {
                                case.system.reservoirs[node]
                                    .level
                                    .trial_level(values.storage[node][j])
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
            .map(|t| members.iter().any(|&g| values.controls.commitment[g][t]))
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
                    .max((values.plant_power[i][t] - p).abs() / p_scale);
            }
            if let Some(q) = limits.fixed_flow {
                result.operating_limits = result
                    .operating_limits
                    .max((values.plant_flow[i][t] - q).abs() / q_scale);
            }
            result.operating_limits = result.operating_limits.max(violation(
                values.plant_flow[i][t],
                if running { limits.flow.min } else { 0. },
                limits.flow.max,
                q_scale,
            ));
            result.operating_limits = result.operating_limits.max(violation(
                values.plant_power[i][t],
                if running { limits.power.min } else { 0. },
                limits.power.max,
                p_scale,
            ));
            let prior_q = if t == 0 {
                operating.plant_initial_flow[i]
            } else {
                Some(values.plant_flow[i][t - 1])
            };
            if let Some(prior) = prior_q {
                let change = values.plant_flow[i][t] - prior;
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
                Some(values.plant_power[i][t - 1])
            };
            if let Some(prior) = prior_p {
                let jump_up: f64 = members
                    .iter()
                    .filter(|&&g| {
                        case.system.units[g]
                            .operation
                            .transitions(&values.controls.commitment[g], t)
                            .0
                    })
                    .map(|&g| operating.units[g][t].power.min)
                    .sum();
                let jump_down: f64 = members
                    .iter()
                    .filter(|&&g| {
                        case.system.units[g]
                            .operation
                            .transitions(&values.controls.commitment[g], t)
                            .1
                    })
                    .map(|&g| operating.units[g][t].stopping_power_min)
                    .sum();
                result.operating_limits = result.operating_limits.max(violation(
                    values.plant_power[i][t] - prior,
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
                .filter(|group| group.iter().any(|&g| values.controls.commitment[g][t]))
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

#[cfg(test)]
mod tests {
    use crate::*;
    use std::collections::BTreeMap;

    #[test]
    fn conserved_flux_value_matches_terminal_storage_and_historical_river_inventory() {
        let case = ScheduleCase {
            system: HydroSystem {
                reservoirs: vec![Reservoir {
                    name: "Lake".into(),
                    storage: Bounds::new(0., 1e6),
                    initial_storage: 5e5,
                    level: StorageCurve::Polynomial {
                        reference_volume: 0.,
                        elevation: 100.,
                        slope: 0.,
                        curvature: 0.,
                    },
                    water_value: 0.002,
                    operation: Default::default(),
                }],
                boundaries: vec![Boundary {
                    name: "Sea".into(),
                    head: 20.,
                }],
                plants: vec![Plant::new("Station", "Lake", "Sea")],
                units: vec![Unit::new(
                    "Unit",
                    "Station",
                    Bounds::new(0., 30.),
                    Bounds::new(0., 100.),
                    Bounds::new(0., 200.),
                    Efficiency::constant(0.9),
                )],
                rivers: vec![River {
                    operation: Default::default(),
                    name: "Reach".into(),
                    source: "Lake".into(),
                    target: "Sea".into(),
                    capacity: 20.,
                    transport: RiverTransport::Fixed { delay: 0.5 },
                    allocation: 1.,
                    inflow: 2.,
                    water_value: 0.003,
                    history_grid: vec![-1., 0.],
                    history_release: vec![6.],
                    law: RiverLaw::Controlled,
                    gate: Bounds::new(0., 1.),
                }],
                ..Default::default()
            },
            grid: vec![0., 1., 1.5],
            prices: vec![100.; 2],
            inflows: BTreeMap::from([("Lake".into(), vec![3.; 2])]),
            commitment: Default::default(),
            operations: vec![],
            soft_constraints: vec![],
            flow_requirements: vec![],
            routing_grid: vec![],
        };
        let problem = SchedulingProblem::new(&case).unwrap();
        let physical = problem
            .reconstruct(&Controls {
                head_guess: None,
                commitment: vec![vec![true; 2]],
                unit_flow: vec![vec![10.; 2]],
                river_opening: vec![vec![0.25, 0.75]],
                tunnel_opening: vec![],
            })
            .unwrap();
        assert!(physical.validation.feasible, "{:?}", physical.validation);
        let mut columns: Vec<_> = physical
            .unit_flow
            .iter()
            .flatten()
            .copied()
            .chain(physical.river_withdrawal.iter().flatten().copied())
            .collect();
        problem.routing.complete_columns(&mut columns);
        let terminal = 0.002 * (physical.storage[0][2] - physical.storage[0][0])
            + 0.003 * (physical.river_storage[0][2] - physical.river_storage[0][0]);
        assert!((terminal + 105.3).abs() < 1e-8);
        assert!((problem.water_value.value(&columns) - terminal).abs() < 1e-8);
    }
}
