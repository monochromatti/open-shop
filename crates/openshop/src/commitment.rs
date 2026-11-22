//! Unit commitment proposals with linearized generation and exact water balance.
//! These schedules initialize nonlinear dispatch; their objectives are not bounds
//! on the nonlinear problem. No-good cuts let search change earlier commitments.

use crate::{Error, Schedule, ScheduleCase, model::Topology, river::Routing};
use oximo::{ScipOptions, prelude::*, solvers::Scip};
use std::time::Duration;

pub(crate) struct Proposal {
    pub on: Vec<Vec<bool>>,
    pub storage: Vec<Vec<f64>>,
    pub flow: Vec<Vec<f64>>,
    pub tunnel_flow: Vec<Vec<f64>>,
    pub river_withdrawal: Vec<Vec<f64>>,
}

pub(crate) fn propose(
    case: &ScheduleCase,
    topology: &Topology,
    routing: &Routing,
    previous: Option<&Schedule>,
    excluded: &[Vec<Vec<bool>>],
    seconds: f64,
) -> Result<Option<Proposal>, Error> {
    let s = &case.system;
    let (nr, nu, nt) = (s.reservoirs.len(), s.units.len(), case.prices.len());
    let m = Model::new("unit commitment");
    variable!(m,0.0<=v[r in 0..nr,t in 0..nt+1]<=1.0);
    variable!(m,q[g in 0..nu,t in 0..nt]>=0.0);
    variable!(m,qt[e in 0..s.tunnels.len(),t in 0..nt]);
    variable!(m,qr[e in 0..routing.outlets.len(),t in 0..nt]>=0.);
    variable!(m,on[g in 0..nu,t in 0..nt],Bin);
    variable!(m,start[g in 0..nu,t in 0..nt],Bin);
    let zero = || -> Expr<'_> { Expr::constant(v[[0, 0]].arena(), 0.).into() };
    let flow_scale = s
        .units
        .iter()
        .map(|u| u.flow.max)
        .chain(s.tunnels.iter().map(|e| e.capacity))
        .chain(s.rivers.iter().map(|r| r.capacity))
        .fold(1., f64::max);
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
            river_columns.push(Expr::from(qr[[k, t]] * flow_scale));
        }
    }
    for (i, series) in routing.releases.iter().enumerate() {
        for rate in &series.rates {
            let release = crate::nlp::linear_expression(rate, &river_columns, zero());
            constraint!(m, release / flow_scale <= s.rivers[i].capacity / flow_scale);
        }
    }
    for (e, tunnel) in s.tunnels.iter().enumerate() {
        for t in 0..nt {
            let capacity = if tunnel.opening == 0. {
                0.
            } else {
                tunnel.capacity / flow_scale
            };
            m.set_bounds(qt[[e, t]], -capacity, capacity)
                .map_err(|e| Error::Model(e.to_string()))?;
        }
    }
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        let range = reservoir.storage.max - reservoir.storage.min;
        m.fix(
            v[[r, 0]],
            (reservoir.initial_storage - reservoir.storage.min) / range,
        )
        .map_err(|e| Error::Model(e.to_string()))?;
        for t in 0..nt {
            let mut net = zero() + case.inflow(r, t);
            for &(g, sign) in &topology.incidence[r] {
                net = net + q[[g, t]] * (sign * flow_scale);
            }
            for &(e, sign) in &topology.tunnel_incidence[r] {
                net = net + qt[[e, t]] * (sign * flow_scale);
            }
            net = net
                + crate::nlp::linear_expression(
                    &routing.reservoir_net[r][t],
                    &river_columns,
                    zero(),
                );
            constraint!(
                m,
                (v[[r, t + 1]] - v[[r, t]]) * (range / (3600. * (case.grid[t + 1] - case.grid[t])))
                    / flow_scale
                    - net / flow_scale
                    == 0.
            );
            for (fraction, correction) in &routing.storage_shape[r][t] {
                if *fraction > 0. && *fraction < 1. {
                    let volume = v[[r, t]] * (1. - fraction)
                        + v[[r, t + 1]] * *fraction
                        + crate::nlp::linear_expression(correction, &river_columns, zero()) / range;
                    constraint!(m, volume >= 0.);
                    constraint!(m, volume <= 1.);
                }
            }
        }
    }
    for j in 0..s.junctions.len() {
        for t in 0..nt {
            let node = nr + j;
            let mut net = zero();
            for &(g, sign) in &topology.incidence[node] {
                net = net + q[[g, t]] * sign;
            }
            for &(e, sign) in &topology.tunnel_incidence[node] {
                net = net + qt[[e, t]] * sign;
            }
            constraint!(m, net == 0.);
        }
    }
    let initial_head: Vec<_> = s
        .reservoirs
        .iter()
        .map(|r| r.level.trial_level(r.initial_storage))
        .chain(s.junctions.iter().map(|j| (j.head.min + j.head.max) / 2.))
        .chain(s.boundaries.iter().map(|b| b.head))
        .collect();
    let mut revenue = zero();
    for (g, unit) in s.units.iter().enumerate() {
        for row in &topology.unit_commitment[g] {
            let lhs = row
                .coefficients
                .iter()
                .fold(zero(), |sum, &(t, a)| sum + on[[g, t]] * a);
            constraint!(m, lhs >= row.rhs);
        }
        let (a, b) = topology.unit_nodes[g];
        for t in 0..nt {
            if let Some(fixed) = case.commitment.get(&unit.name) {
                m.fix(on[[g, t]], if fixed[t] { 1. } else { 0. })
                    .map_err(|e| Error::Model(e.to_string()))?;
            }
            let head = previous
                .map_or(initial_head[a] - initial_head[b], |p| {
                    p.head[a][t] - p.head[b][t]
                })
                .clamp(unit.head.min, unit.head.max);
            let flow = previous
                .map_or((unit.flow.min + unit.flow.max) / 2., |p| p.unit_flow[g][t])
                .clamp(unit.flow.min, unit.flow.max);
            let head = (head - unit.loss_factor * flow * flow).clamp(unit.head.min, unit.head.max);
            let coefficient = 0.00981 * head * unit.efficiency.value(flow, head);
            constraint!(m, q[[g, t]] >= on[[g, t]] * (unit.flow.min / flow_scale));
            constraint!(m, q[[g, t]] <= on[[g, t]] * (unit.flow.max / flow_scale));
            constraint!(
                m,
                q[[g, t]] * (flow_scale * coefficient) >= on[[g, t]] * unit.power.min
            );
            constraint!(
                m,
                q[[g, t]] * (flow_scale * coefficient) <= on[[g, t]] * unit.power.max
            );
            if t == 0 {
                constraint!(
                    m,
                    start[[g, t]] >= on[[g, t]] - if unit.operation.initial_on { 1. } else { 0. }
                );
            } else {
                constraint!(m, start[[g, t]] >= on[[g, t]] - on[[g, t - 1]]);
            }
            let previous_on = if t == 0 {
                zero() + if unit.operation.initial_on { 1. } else { 0. }
            } else {
                on[[g, t - 1]].into()
            };
            constraint!(m, start[[g, t]] <= on[[g, t]]);
            constraint!(m, start[[g, t]] <= 1. - previous_on);
            let stop = start[[g, t]] - on[[g, t]] + previous_on;
            let previous_flow = if t == 0 {
                unit.operation
                    .previous_flow()
                    .map(|v| zero() + v / flow_scale)
            } else {
                Some(q[[g, t - 1]].into())
            };
            if let Some(previous) = previous_flow {
                let change = q[[g, t]] - previous;
                let dt = unit.operation.ramp_hours(&case.grid, t);
                if let Some(up) = unit.operation.flow_ramp.up {
                    constraint!(m, change <= up * dt / flow_scale);
                }
                if let Some(down) = unit.operation.flow_ramp.down {
                    constraint!(m, -change <= down * dt / flow_scale);
                }
            }
            revenue = revenue
                + q[[g, t]]
                    * (flow_scale
                        * coefficient
                        * case.prices[t]
                        * (case.grid[t + 1] - case.grid[t]))
                - start[[g, t]] * unit.operation.startup_cost
                - stop * unit.operation.shutdown_cost;
        }
    }
    for commitment in excluded {
        let mut difference = zero();
        for g in 0..nu {
            for t in 0..nt {
                difference = if commitment[g][t] {
                    difference + 1. - on[[g, t]]
                } else {
                    difference + on[[g, t]]
                };
            }
        }
        constraint!(m, difference >= 1.);
    }
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        revenue = revenue
            + (v[[r, nt]] * (reservoir.storage.max - reservoir.storage.min)
                + reservoir.storage.min)
                * reservoir.water_value;
    }
    for (r, river) in s.rivers.iter().enumerate() {
        revenue = revenue
            + crate::nlp::linear_expression(&routing.inventory[r][nt], &river_columns, zero())
                * river.water_value;
    }
    objective!(m, Max, revenue / 1000.);
    let options = ScipOptions::default()
        .time_limit(Duration::from_secs_f64(seconds))
        .threads(1)
        .verbose(false);
    let result = Scip::new()
        .solve(&m, &options)
        .map_err(|e| Error::Solver(e.to_string()))?;
    if !result.has_solution() {
        return Ok(None);
    }
    let value = |variable| {
        result
            .value_of(variable)
            .map_err(|e| Error::Solver(e.to_string()))?
            .ok_or_else(|| Error::Solver("SCIP returned an incomplete solution".into()))
    };
    let mut commitment = vec![vec![false; nt]; nu];
    let mut storage = vec![vec![0.; nt + 1]; nr];
    let mut flow = vec![vec![0.; nt]; nu];
    let mut tunnel_flow = vec![vec![0.; nt]; s.tunnels.len()];
    let mut river_withdrawal = vec![vec![0.; nt]; s.rivers.len()];
    for r in 0..nr {
        for t in 0..=nt {
            let volume = value(v[[r, t]])?;
            storage[r][t] = s.reservoirs[r].storage.min
                + volume * (s.reservoirs[r].storage.max - s.reservoirs[r].storage.min);
        }
    }
    for g in 0..nu {
        for t in 0..nt {
            commitment[g][t] = value(on[[g, t]])? > 0.5;
            flow[g][t] = flow_scale * value(q[[g, t]])?;
        }
    }
    for e in 0..s.tunnels.len() {
        for t in 0..nt {
            tunnel_flow[e][t] = flow_scale * value(qt[[e, t]])?;
        }
    }
    for (k, &r) in routing.outlets.iter().enumerate() {
        for t in 0..nt {
            river_withdrawal[r][t] = flow_scale * value(qr[[k, t]])?;
        }
    }
    Ok(Some(Proposal {
        on: commitment,
        storage,
        flow,
        tunnel_flow,
        river_withdrawal,
    }))
}
