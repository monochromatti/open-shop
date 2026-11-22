//! Nonlinear dispatch and independent evaluation of its physical equations.

use crate::{
    Error, ScheduleCase, SolveOptions,
    model::Topology,
    nlp::{self, Backend},
    river::Routing,
};
use oximo::{
    expr::{Expr, VarId},
    prelude::*,
};
use serde::{Deserialize, Serialize};
use std::time::Instant;

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

struct Dispatch {
    problem: pounce_nl::nl_reader::NlTnlp,
    storage: Vec<Vec<VarId>>,
    flow: Vec<Vec<VarId>>,
    tunnel_flow: Vec<Vec<VarId>>,
    river_withdrawal: Vec<Vec<VarId>>,
    junction_head: Vec<Vec<VarId>>,
    flow_scale: f64,
}

/// Optimize on/off schedules with SCIP proposals and exact nonlinear dispatch.
/// Prescribed commitments remain fixed; other units can change commitment
/// between nonlinear solves. This search does not establish a global bound.
/// An infeasible or unfinished solve is reported through `validation.feasible`;
/// a solver's success flag alone does not establish physical feasibility.
pub fn solve(case: &ScheduleCase, options: &SolveOptions) -> Result<Schedule, Error> {
    solve_with_backend(case, options, Backend::Ipopt)
}

fn solve_with_backend(
    case: &ScheduleCase,
    options: &SolveOptions,
    backend: Backend,
) -> Result<Schedule, Error> {
    options.validate()?;
    let topology = case.topology()?;
    let routing = Routing::new(case, &topology)?;
    if case.commitment.len() == case.system.units.len() {
        let commitment = case
            .system
            .units
            .iter()
            .map(|u| case.commitment[&u.name].clone())
            .collect();
        return dispatch(
            case, &topology, &routing, commitment, None, options, backend,
        );
    }
    let started = Instant::now();
    let mut excluded = Vec::new();
    let mut best: Option<Schedule> = None;
    let mut last = None;
    while options.time_limit - started.elapsed().as_secs_f64() > 0.01 {
        let remaining = options.time_limit - started.elapsed().as_secs_f64();
        let Some(proposal) = crate::commitment::propose(
            case,
            &topology,
            &routing,
            best.as_ref(),
            &excluded,
            remaining.min(options.time_limit / 4.),
        )?
        else {
            break;
        };
        excluded.push(proposal.on.clone());
        let mut budget = options.clone();
        budget.time_limit =
            (options.time_limit - started.elapsed().as_secs_f64()).min(options.time_limit / 4.);
        if budget.time_limit <= 0. {
            break;
        }
        budget.max_iterations = options.max_iterations.min(200);
        let result = dispatch(
            case,
            &topology,
            &routing,
            proposal.on.clone(),
            Some(&proposal),
            &budget,
            backend,
        )?;
        if result.validation.feasible
            && best.as_ref().is_none_or(|p| result.objective > p.objective)
        {
            best = Some(result.clone());
        }
        last = Some(result);
    }
    let result = best.or(last).ok_or_else(|| {
        Error::Solver("no commitment produced a nonlinear trial point within the time limit".into())
    })?;
    let remaining = options.time_limit - started.elapsed().as_secs_f64();
    if result.validation.feasible && remaining > 0.01 {
        let start = crate::commitment::Proposal {
            on: result.commitment.clone(),
            storage: result.storage.clone(),
            flow: result.unit_flow.clone(),
            tunnel_flow: result.tunnel_flow.clone(),
            river_withdrawal: result.river_withdrawal.clone(),
        };
        let mut budget = options.clone();
        budget.time_limit = remaining;
        let refined = dispatch(
            case,
            &topology,
            &routing,
            result.commitment.clone(),
            Some(&start),
            &budget,
            backend,
        )?;
        if refined.validation.feasible && refined.objective >= result.objective {
            return Ok(refined);
        }
    }
    Ok(result)
}

fn dispatch(
    case: &ScheduleCase,
    topology: &Topology,
    routing: &Routing,
    commitment: Vec<Vec<bool>>,
    start: Option<&crate::commitment::Proposal>,
    options: &SolveOptions,
    backend: Backend,
) -> Result<Schedule, Error> {
    let dispatch = build_dispatch(case, topology, routing, &commitment, start)?;
    let point = nlp::solve(dispatch.problem, options, backend)?;
    let storage = dispatch
        .storage
        .iter()
        .enumerate()
        .map(|(r, row)| {
            let bounds = case.system.reservoirs[r].storage;
            row.iter()
                .map(|i| bounds.min + (bounds.max - bounds.min) * point.x[i.index()])
                .collect()
        })
        .collect();
    let unit_flow = dispatch
        .flow
        .iter()
        .map(|row| {
            row.iter()
                .map(|i| dispatch.flow_scale * point.x[i.index()])
                .collect()
        })
        .collect();
    let tunnel_flow = dispatch
        .tunnel_flow
        .iter()
        .map(|row| {
            row.iter()
                .map(|i| dispatch.flow_scale * point.x[i.index()])
                .collect()
        })
        .collect();
    let mut schedule = Schedule {
        storage,
        unit_flow,
        tunnel_flow,
        river_withdrawal: vec![vec![0.; case.prices.len()]; case.system.rivers.len()],
        river_release: Vec::new(),
        river_arrival_volume: Vec::new(),
        river_storage: Vec::new(),
        unit_power: Vec::new(),
        head: Vec::new(),
        commitment,
        objective: 0.,
        solver_status: point.status,
        iterations: point.iterations,
        validation: Validation::default(),
    };
    let s = &case.system;
    for (k, &river) in routing.outlets.iter().enumerate() {
        schedule.river_withdrawal[river] = dispatch.river_withdrawal[k]
            .iter()
            .map(|i| dispatch.flow_scale * point.x[i.index()])
            .collect();
    }
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
                        .map(|w| (w[1].0 - w[0].0) * reservoir.level.mean_level(w[0].1, w[1].1))
                        .sum()
                })
                .collect()
        })
        .chain(
            dispatch
                .junction_head
                .iter()
                .map(|row| row.iter().map(|i| 100. * point.x[i.index()]).collect()),
        )
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
                    u.generation(
                        q,
                        schedule.head[source][t] - schedule.head[target][t] - u.loss_factor * q * q,
                    )
                })
                .collect()
        })
        .collect();
    schedule.objective = objective(case, &schedule);
    schedule.validation = validate(
        case,
        topology,
        routing,
        &schedule,
        options.tolerance.max(1e-7),
    );
    Ok(schedule)
}

fn build_dispatch(
    case: &ScheduleCase,
    topology: &Topology,
    routing: &Routing,
    commitment: &[Vec<bool>],
    start: Option<&crate::commitment::Proposal>,
) -> Result<Dispatch, Error> {
    let s = &case.system;
    let (nr, nu, nt) = (s.reservoirs.len(), s.units.len(), case.prices.len());
    let m = Model::new("hydropower dispatch");
    variable!(m,0.0 <= v[r in 0..nr,t in 0..nt+1] <= 1.0);
    variable!(m,q[u in 0..nu,t in 0..nt] >= 0.0);
    variable!(m,qt[e in 0..s.tunnels.len(),t in 0..nt]);
    variable!(m,hj[j in 0..s.junctions.len(),t in 0..nt]);
    variable!(m,qr[e in 0..routing.outlets.len(),t in 0..nt]>=0.);
    let flow_scale = s
        .units
        .iter()
        .map(|u| u.flow.max)
        .chain(s.tunnels.iter().map(|e| e.capacity))
        .chain(s.rivers.iter().map(|r| r.capacity))
        .fold(1.0, f64::max);
    let zero = || -> Expr<'_> { Expr::constant(v[[0, 0]].arena(), 0.0).into() };
    let mut river_columns = Vec::new();
    for u in 0..nu {
        for t in 0..nt {
            river_columns.push(Expr::from(q[[u, t]] * flow_scale));
        }
    }
    for (k, &e) in routing.outlets.iter().enumerate() {
        for t in 0..nt {
            m.set_bounds(qr[[k, t]], 0., s.rivers[e].capacity / flow_scale)
                .map_err(|e| Error::Model(e.to_string()))?;
            m.set_initial(
                qr[[k, t]],
                start
                    .map_or(0., |p| p.river_withdrawal[e][t] / flow_scale)
                    .clamp(0., s.rivers[e].capacity / flow_scale),
            )
            .map_err(|e| Error::Model(e.to_string()))?;
            river_columns.push(Expr::from(qr[[k, t]] * flow_scale));
        }
    }
    for (i, series) in routing.releases.iter().enumerate() {
        for rate in &series.rates {
            let release = rate.expression(&river_columns, zero());
            constraint!(m, release / flow_scale <= s.rivers[i].capacity / flow_scale);
        }
    }
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        let range = reservoir.storage.max - reservoir.storage.min;
        if range <= 0. {
            return Err(Error::Input(format!(
                "{}: positive storage range required",
                reservoir.name
            )));
        }
        let initial = (reservoir.initial_storage - reservoir.storage.min) / range;
        m.fix(v[[r, 0]], initial)
            .map_err(|e| Error::Model(e.to_string()))?;
        for t in 0..=nt {
            let guess = start
                .map_or(initial, |p| {
                    (p.storage[r][t] - reservoir.storage.min) / range
                })
                .clamp(0., 1.);
            m.set_initial(v[[r, t]], guess)
                .map_err(|e| Error::Model(e.to_string()))?;
        }
        for t in 0..nt {
            let mut net = zero() + case.inflow(r, t);
            for &(u, sign) in &topology.incidence[r] {
                net = net + q[[u, t]] * (sign * flow_scale);
            }
            for &(e, sign) in &topology.tunnel_incidence[r] {
                net = net + qt[[e, t]] * (sign * flow_scale);
            }
            net = net + routing.reservoir_net[r][t].expression(&river_columns, zero());
            let rate =
                (v[[r, t + 1]] - v[[r, t]]) * (range / (3600. * (case.grid[t + 1] - case.grid[t])));
            constraint!(
                m,
                name = format!("water_balance[{r},{t}]"),
                (rate - net) / flow_scale == 0.0
            );
        }
    }
    for (j, junction) in s.junctions.iter().enumerate() {
        let node = nr + j;
        for t in 0..nt {
            m.set_bounds(
                hj[[j, t]],
                junction.head.min / 100.,
                junction.head.max / 100.,
            )
            .map_err(|e| Error::Model(e.to_string()))?;
            m.set_initial(hj[[j, t]], (junction.head.min + junction.head.max) / 200.)
                .map_err(|e| Error::Model(e.to_string()))?;
            let mut net = zero();
            for &(u, sign) in &topology.incidence[node] {
                net = net + q[[u, t]] * sign;
            }
            for &(e, sign) in &topology.tunnel_incidence[node] {
                net = net + qt[[e, t]] * sign;
            }
            constraint!(m, name = format!("junction_balance[{j},{t}]"), net == 0.);
        }
    }
    let mut head = Vec::new();
    for t in 0..nt {
        let mut values: Vec<Expr<'_>> = Vec::new();
        for (r, reservoir) in s.reservoirs.iter().enumerate() {
            let range = reservoir.storage.max - reservoir.storage.min;
            let c = reservoir.level;
            let mut vertices = Vec::new();
            for (fraction, correction) in &routing.storage_shape[r][t] {
                let volume = v[[r, t]] * (1. - fraction)
                    + v[[r, t + 1]] * *fraction
                    + correction.expression(&river_columns, zero()) / range;
                if *fraction > 0. && *fraction < 1. {
                    constraint!(m, volume >= 0.);
                    constraint!(m, volume <= 1.);
                }
                vertices.push((
                    *fraction,
                    volume * range + (reservoir.storage.min - c.reference_volume),
                ));
            }
            let mut level = zero();
            for w in vertices.windows(2) {
                let (a, b) = (w[0].1, w[1].1);
                level = level
                    + ((a + b) * (c.slope / 2.)
                        + (a.square() + a * b + b.square()) * (c.curvature / 3.)
                        + c.elevation)
                        * (w[1].0 - w[0].0);
            }
            values.push(level);
        }
        values.extend(
            s.junctions
                .iter()
                .enumerate()
                .map(|(j, _)| Expr::from(hj[[j, t]] * 100.)),
        );
        values.extend(s.boundaries.iter().map(|b| zero() + b.head));
        head.push(values);
    }
    let initial_head: Vec<_> = s
        .reservoirs
        .iter()
        .map(|r| r.level.level(r.initial_storage))
        .chain(s.junctions.iter().map(|j| (j.head.min + j.head.max) / 2.))
        .chain(s.boundaries.iter().map(|b| b.head))
        .collect();
    for (e, tunnel) in s.tunnels.iter().enumerate() {
        let (source, target) = topology.tunnel_nodes[e];
        let dh = initial_head[source] - initial_head[target];
        let guess = dh.signum() * (dh.abs() / tunnel.loss_factor).sqrt() * tunnel.opening;
        for t in 0..nt {
            let capacity = if tunnel.opening == 0. {
                0.
            } else {
                tunnel.capacity / flow_scale
            };
            m.set_bounds(qt[[e, t]], -capacity, capacity)
                .map_err(|e| Error::Model(e.to_string()))?;
            m.set_initial(
                qt[[e, t]],
                start
                    .map_or(guess / flow_scale, |p| p.tunnel_flow[e][t] / flow_scale)
                    .clamp(-capacity, capacity),
            )
            .map_err(|e| Error::Model(e.to_string()))?;
            if tunnel.opening > 0. {
                let flow = qt[[e, t]] * flow_scale;
                let pressure = (head[t][source] - head[t][target]) * tunnel.opening.powi(2)
                    - flow * flow.abs() * tunnel.loss_factor;
                constraint!(
                    m,
                    name = format!("pressure_loss[{e},{t}]"),
                    pressure / 100. == 0.
                );
            }
        }
    }
    let mut revenue = zero();
    for (i, u) in s.units.iter().enumerate() {
        let (source, target) = topology.unit_nodes[i];
        for t in 0..nt {
            let on = commitment[i][t];
            let (lo, hi) = if on {
                (u.flow.min / flow_scale, u.flow.max / flow_scale)
            } else {
                (0., 0.)
            };
            m.set_bounds(q[[i, t]], lo, hi)
                .map_err(|e| Error::Model(e.to_string()))?;
            m.set_initial(
                q[[i, t]],
                start
                    .map_or((lo + hi) / 2., |p| p.flow[i][t] / flow_scale)
                    .clamp(lo, hi),
            )
            .map_err(|e| Error::Model(e.to_string()))?;
            let flow = q[[i, t]] * flow_scale;
            let net_head = head[t][source] - head[t][target] - flow.square() * u.loss_factor;
            let eta = -(flow - u.efficiency.best_flow).square() * u.efficiency.flow_curvature
                - (net_head - u.efficiency.best_head).square() * u.efficiency.head_curvature
                + u.efficiency.peak;
            let power = flow * net_head * eta * 0.00981;
            if on {
                constraint!(
                    m,
                    name = format!("head_min[{i},{t}]"),
                    net_head / 100. >= u.head.min / 100.
                );
                constraint!(
                    m,
                    name = format!("head_max[{i},{t}]"),
                    net_head / 100. <= u.head.max / 100.
                );
                constraint!(
                    m,
                    name = format!("power_min[{i},{t}]"),
                    power / u.power.max.max(1.) >= u.power.min / u.power.max.max(1.)
                );
                constraint!(
                    m,
                    name = format!("power_max[{i},{t}]"),
                    power / u.power.max.max(1.) <= u.power.max / u.power.max.max(1.)
                );
            }
            revenue = revenue + power * (case.prices[t] * (case.grid[t + 1] - case.grid[t]));
            if on
                && !(if t == 0 {
                    u.initial_on
                } else {
                    commitment[i][t - 1]
                })
            {
                revenue = revenue - u.startup_cost;
            }
        }
    }
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        revenue = revenue
            + (v[[r, nt]] * (reservoir.storage.max - reservoir.storage.min)
                + reservoir.storage.min)
                * reservoir.water_value;
    }
    for (r, river) in s.rivers.iter().enumerate() {
        revenue = revenue
            + routing.inventory[r][nt].expression(&river_columns, zero()) * river.water_value;
    }
    objective!(m, Max, revenue / 1000.);
    let storage = (0..nr)
        .map(|r| (0..=nt).map(|t| v[[r, t]].var_id().unwrap()).collect())
        .collect();
    let flow = (0..nu)
        .map(|u| (0..nt).map(|t| q[[u, t]].var_id().unwrap()).collect())
        .collect();
    let tunnel_flow = (0..s.tunnels.len())
        .map(|e| (0..nt).map(|t| qt[[e, t]].var_id().unwrap()).collect())
        .collect();
    let junction_head = (0..s.junctions.len())
        .map(|j| (0..nt).map(|t| hj[[j, t]].var_id().unwrap()).collect())
        .collect();
    let river_withdrawal = (0..routing.outlets.len())
        .map(|r| (0..nt).map(|t| qr[[r, t]].var_id().unwrap()).collect())
        .collect();
    Ok(Dispatch {
        problem: nlp::prepare(&m)?,
        storage,
        flow,
        tunnel_flow,
        junction_head,
        river_withdrawal,
        flow_scale,
    })
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
                    let start = schedule.commitment[i][t]
                        && !(if t == 0 {
                            u.initial_on
                        } else {
                            schedule.commitment[i][t - 1]
                        });
                    price * schedule.unit_power[i][t] * (case.grid[t + 1] - case.grid[t])
                        - if start { u.startup_cost } else { 0. }
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
        for t in 0..case.prices.len() {
            let (a, b) = topology.tunnel_nodes[e];
            let q = schedule.tunnel_flow[e][t];
            let error = if tunnel.opening == 0. {
                q.abs()
            } else {
                ((schedule.head[a][t] - schedule.head[b][t]) * tunnel.opening.powi(2)
                    - tunnel.loss_factor * q * q.abs())
                .abs()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Boundary, Bounds, Efficiency, HydroSystem, Reservoir, StorageCurve, Unit};
    use std::collections::BTreeMap;

    fn reservoir_case() -> ScheduleCase {
        ScheduleCase {
            system: HydroSystem {
                reservoirs: vec![Reservoir {
                    name: "Lake".into(),
                    storage: Bounds::new(1e6, 4e6),
                    initial_storage: 3e6,
                    level: StorageCurve {
                        reference_volume: 0.,
                        elevation: 40.,
                        slope: 20e-6,
                        curvature: 0.,
                    },
                    water_value: 0.001,
                }],
                boundaries: vec![Boundary {
                    name: "Sea".into(),
                    head: 20.,
                }],
                units: vec![Unit {
                    name: "Unit".into(),
                    source: "Lake".into(),
                    target: "Sea".into(),
                    discharge_river: None,
                    flow: Bounds::new(5., 20.),
                    power: Bounds::new(1., 25.),
                    head: Bounds::new(35., 105.),
                    efficiency: Efficiency::constant(0.9),
                    loss_factor: 0.,
                    startup_cost: 0.,
                    initial_on: false,
                }],
                ..Default::default()
            },
            grid: vec![0., 1.],
            prices: vec![100.],
            inflows: BTreeMap::new(),
            commitment: BTreeMap::from([("Unit".into(), vec![true])]),
        }
    }

    #[test]
    fn complete_reservoir_optimization() {
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let case = reservoir_case();
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!(
                (result.unit_flow[0][0] - 20.).abs() < 1e-5,
                "{backend:?}: {:?}",
                result.unit_flow
            );
            assert!((result.storage[0][1] - 2_928_000.).abs() < 0.05);
            let expected_power =
                0.00981 * 20. * 0.9 * (40. + 20e-6 * (3e6 + 2_928_000.) / 2. - 20.);
            let expected = 100. * expected_power + 0.001 * 2_928_000.;
            assert!(
                (result.objective - expected).abs() < 0.01,
                "{backend:?}: {} vs {expected}",
                result.objective
            );
        }
    }

    #[test]
    fn off_units_have_exactly_zero_discharge_and_generation() {
        let mut case = reservoir_case();
        case.commitment.get_mut("Unit").unwrap()[0] = false;
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.unit_flow[0][0], 0.);
        assert_eq!(result.unit_power[0][0], 0.);
        assert!((result.storage[0][1] - 3e6).abs() < 1e-6);
    }

    #[test]
    fn storage_head_is_energy_consistent_even_for_equal_volumes() {
        let c = StorageCurve {
            reference_volume: 1e6,
            elevation: 75.,
            slope: 2e-5,
            curvature: 1e-11,
        };
        for (a, b) in [(1e6, 2e6), (2e6, 1e6), (1e6, 1e6), (1e6, 1e6 + 1e-5)] {
            let mean = c.mean_level(a, b);
            assert!(mean.is_finite());
            if a == b {
                assert_eq!(mean, c.level(a));
            } else {
                assert!(
                    ((b - a) * mean - (c.potential(b) - c.potential(a))).abs()
                        < 1e-7 * (1. + ((b - a) * mean).abs())
                );
            }
        }
    }

    #[test]
    fn commitment_can_turn_off_during_negative_prices() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 1.5];
        case.prices = vec![-100., 150.];
        case.commitment.clear();
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(
                &case,
                &SolveOptions {
                    time_limit: 5.,
                    ..Default::default()
                },
                backend,
            )
            .unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert_eq!(result.commitment[0], vec![false, true]);
            assert_eq!(result.unit_flow[0][0], 0.);
            assert!((result.unit_flow[0][1] - 20.).abs() < 1e-5);
        }
    }

    #[test]
    fn pressure_loop_conserves_water_and_handles_reversed_tunnel_orientation() {
        use crate::{Junction, Tunnel};
        let mut case = reservoir_case();
        case.system.junctions.push(Junction {
            name: "Intake".into(),
            head: Bounds::new(25., 120.),
        });
        case.system.units[0].source = "Intake".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        case.system.tunnels = vec![
            Tunnel {
                name: "Main".into(),
                source: "Lake".into(),
                target: "Intake".into(),
                capacity: 20.,
                loss_factor: 0.1,
                opening: 1.,
            },
            Tunnel {
                name: "Parallel".into(),
                source: "Intake".into(),
                target: "Lake".into(),
                capacity: 20.,
                loss_factor: 0.4,
                opening: 1.,
            },
        ];
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.unit_flow[0][0] - 20.).abs() < 1e-4);
            let (a, b) = (result.tunnel_flow[0][0], result.tunnel_flow[1][0]);
            assert!(a > 0. && b < 0.);
            assert!((a / (-b) - 2.).abs() < 1e-6);
            assert!((a - b - result.unit_flow[0][0]).abs() < 1e-6);
        }
    }

    #[test]
    fn small_reservoir_couples_two_plants_without_small_time_steps() {
        let mut case = reservoir_case();
        case.system.reservoirs.push(Reservoir {
            name: "Forebay".into(),
            storage: Bounds::new(0., 1000.),
            initial_storage: 500.,
            level: StorageCurve {
                reference_volume: 0.,
                elevation: 45.,
                slope: 0.01,
                curvature: 0.,
            },
            water_value: 0.,
        });
        case.system.units[0].target = "Forebay".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        let mut downstream = case.system.units[0].clone();
        downstream.name = "Downstream".into();
        downstream.source = "Forebay".into();
        downstream.target = "Sea".into();
        case.system.units.push(downstream);
        case.grid = vec![0., 1., 3.];
        case.prices = vec![30., 150.];
        case.commitment = BTreeMap::from([
            ("Unit".into(), vec![true, true]),
            ("Downstream".into(), vec![true, true]),
        ]);
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            for (t, dt) in [1., 2.].iter().enumerate() {
                let change = 3600. * dt * (result.unit_flow[0][t] - result.unit_flow[1][t]);
                assert!((change - (result.storage[1][t + 1] - result.storage[1][t])).abs() < 0.005);
            }
        }
    }

    fn delayed_forebay_case() -> ScheduleCase {
        let mut case = reservoir_case();
        case.system.reservoirs.push(Reservoir {
            name: "Forebay".into(),
            storage: Bounds::new(0., 1000.),
            initial_storage: 500.,
            level: StorageCurve {
                reference_volume: 0.,
                elevation: 45.,
                slope: 0.01,
                curvature: 0.,
            },
            water_value: 0.,
        });
        case.system.boundaries.push(Boundary {
            name: "Tailwater".into(),
            head: 60.,
        });
        case.system.units[0].target = "Tailwater".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        case.system.units[0].discharge_river = Some("River".into());
        let mut lower = case.system.units[0].clone();
        lower.name = "Lower".into();
        lower.source = "Forebay".into();
        lower.target = "Sea".into();
        lower.discharge_river = None;
        case.system.units.push(lower);
        let mut river = crate::River {
            name: "River".into(),
            source: "Outfall".into(),
            target: "Forebay".into(),
            delay: 0.5,
            capacity: 30.,
            allocation: 1.,
            inflow: 0.,
            water_value: 0.,
            history_grid: Vec::new(),
            history_release: Vec::new(),
        };
        river.history_grid = vec![-1., 0.];
        river.history_release = vec![5.];
        case.system.rivers.push(river);
        case.grid = vec![0., 0.5, 1.5];
        case.prices = vec![-100., 150.];
        case.commitment = BTreeMap::from([
            ("Unit".into(), vec![true, true]),
            ("Lower".into(), vec![true, true]),
        ]);
        case
    }

    #[test]
    fn delayed_outfall_optimization_respects_small_reservoir_between_control_times() {
        let case = delayed_forebay_case();
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.river_arrival_volume[0][0] - 9000.).abs() < 0.01);
            // Within the second interval the arrival changes at hour 1.0.
            let middle =
                result.storage[1][1] + 1800. * (result.unit_flow[0][0] - result.unit_flow[1][1]);
            assert!(
                (-0.01..=1000.01).contains(&middle),
                "{backend:?}: middle storage {middle}"
            );
            let mean_volume = (result.storage[1][1] + 2. * middle + result.storage[1][2]) / 4.;
            assert!((result.head[1][1] - (45. + 0.01 * mean_volume)).abs() < 1e-8);
            for t in 0..case.prices.len() {
                let seconds = 3600. * (case.grid[t + 1] - case.grid[t]);
                let change = result.river_arrival_volume[0][t] - seconds * result.unit_flow[1][t];
                assert!((result.storage[1][t + 1] - result.storage[1][t] - change).abs() < 0.01);
            }
        }
    }

    #[test]
    fn endpoint_balance_cannot_hide_overflow_between_control_times() {
        let mut case = delayed_forebay_case();
        case.grid = vec![0., 1.];
        case.prices = vec![100.];
        case.system.rivers[0].history_release = vec![20.];
        case.commitment =
            BTreeMap::from([("Unit".into(), vec![false]), ("Lower".into(), vec![true])]);
        let result = solve(
            &case,
            &SolveOptions {
                max_iterations: 100,
                ..Default::default()
            },
        )
        .unwrap();
        // Half an hour of 20 m³/s followed by zero cannot be buffered in 1000 m³
        // with a single constant downstream flow, although endpoint balance can.
        assert!(!result.validation.feasible);
    }

    #[test]
    fn commitment_search_includes_delayed_water_and_river_inventory() {
        let mut case = delayed_forebay_case();
        case.commitment.clear();
        let result = solve(
            &case,
            &SolveOptions {
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        // In-flight water forces the downstream unit on in the first interval.
        assert!(result.commitment[1][0]);
        let releases: f64 = result.unit_flow[0]
            .iter()
            .enumerate()
            .map(|(t, q)| q * 3600. * (case.grid[t + 1] - case.grid[t]))
            .sum();
        let arrivals: f64 = result.river_arrival_volume[0].iter().sum();
        assert!(
            (result.river_storage[0][0] + releases - arrivals - result.river_storage[0][2]).abs()
                < 0.01
        );
    }

    #[test]
    fn nonlinear_efficiency_can_make_partial_load_optimal() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level = StorageCurve {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        case.system.reservoirs[0].water_value = 0.;
        let unit = &mut case.system.units[0];
        unit.flow = Bounds::new(5., 15.);
        unit.head = Bounds::new(90., 110.);
        unit.efficiency = Efficiency {
            peak: 0.9,
            best_flow: 10.,
            flow_curvature: 0.01,
            best_head: 100.,
            head_curvature: 0.,
        };
        let expected = (0.4 + 0.148_f64.sqrt()) / 0.06;
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(result.validation.feasible);
            assert!(
                (result.unit_flow[0][0] - expected).abs() < 1e-5,
                "{backend:?}: {:?}",
                result.unit_flow
            );
        }
    }

    #[test]
    fn a_closed_tunnel_disconnects_pressure_without_equating_heads() {
        let mut case = reservoir_case();
        case.system.tunnels.push(crate::Tunnel {
            name: "Closed".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 50.,
            loss_factor: 0.01,
            opening: 0.,
        });
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.tunnel_flow[0][0], 0.);
        assert!(result.head[0][0] - result.head[1][0] > 70.);
    }

    #[test]
    fn an_interrupted_solve_does_not_claim_physical_feasibility() {
        let result = solve(
            &reservoir_case(),
            &SolveOptions {
                max_iterations: 0,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(!result.validation.feasible);
        assert!(result.validation.water_balance > 1.);
    }
}
