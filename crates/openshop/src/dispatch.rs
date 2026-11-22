//! Simultaneous nonlinear water balances, pressure laws and operating equations.

use crate::nlp;
use crate::{Error, ScheduleCase, model::Topology, river::Routing};
use oximo::{
    expr::{Expr, VarId},
    prelude::*,
};

pub(crate) struct Dispatch {
    pub problem: crate::constitutive::Nlp,
    pub storage: Vec<Vec<VarId>>,
    pub flow: Vec<Vec<VarId>>,
    pub tunnel_flow: Vec<Vec<VarId>>,
    pub river_withdrawal: Vec<Vec<VarId>>,
    pub junction_head: Vec<Vec<VarId>>,
    pub flow_scale: f64,
}

pub(crate) fn build_dispatch(
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
    let mut curves = Vec::new();
    let table_index: Vec<_> = s
        .reservoirs
        .iter()
        .map(|r| {
            if let crate::StorageCurve::Table(curve) = &r.level
                && r.level.trial_level(r.storage.min) != r.level.trial_level(r.storage.max)
            {
                let index = curves.len();
                curves.push(std::sync::Arc::new(curve.clone()));
                Some(index)
            } else {
                None
            }
        })
        .collect();
    variable!(m,hm[k in 0..curves.len(),t in 0..nt]);
    let mut turbine_curves = Vec::new();
    let turbine_index: Vec<_> = s
        .units
        .iter()
        .map(|u| match &u.efficiency {
            crate::Efficiency::Table(c) => {
                let k = turbine_curves.len();
                turbine_curves.push(std::sync::Arc::new(c.clone()));
                Some(k)
            }
            _ => None,
        })
        .collect();
    variable!(m,ht[k in 0..turbine_curves.len(),t in 0..nt]);
    variable!(m,pt[k in 0..turbine_curves.len(),t in 0..nt]);
    let mut laws = Vec::new();
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
            let release = crate::nlp::linear_expression(rate, &river_columns, zero());
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
            net = net
                + crate::nlp::linear_expression(
                    &routing.reservoir_net[r][t],
                    &river_columns,
                    zero(),
                );
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
            let c = &reservoir.level;
            let mut vertices = Vec::new();
            for (fraction, correction) in &routing.storage_shape[r][t] {
                let volume = v[[r, t]] * (1. - fraction)
                    + v[[r, t + 1]] * *fraction
                    + crate::nlp::linear_expression(correction, &river_columns, zero()) / range;
                if *fraction > 0. && *fraction < 1. {
                    constraint!(m, volume >= 0.);
                    constraint!(m, volume <= 1.);
                }
                vertices.push((*fraction, volume * range + reservoir.storage.min));
            }
            let level = match c {
                crate::StorageCurve::Polynomial {
                    reference_volume,
                    elevation,
                    slope,
                    curvature,
                } => {
                    let mut level = zero();
                    for w in vertices.windows(2) {
                        level = level
                            + crate::physics::polynomial_mean_level(
                                *reference_volume,
                                *elevation,
                                *slope,
                                *curvature,
                                w[0].1,
                                w[1].1,
                            ) * (w[1].0 - w[0].0);
                    }
                    level
                }
                crate::StorageCurve::Table(_) => {
                    let lo = c.trial_level(reservoir.storage.min);
                    let hi = c.trial_level(reservoir.storage.max);
                    if let Some(k) = table_index[r] {
                        m.set_bounds(hm[[k, t]], lo / 100., hi / 100.)
                            .map_err(|e| Error::Model(e.to_string()))?;
                        let (before, after) = start.map_or(
                            (reservoir.initial_storage, reservoir.initial_storage),
                            |p| (p.storage[r][t], p.storage[r][t + 1]),
                        );
                        m.set_initial(
                            hm[[k, t]],
                            c.trial_mean_level(before, after).clamp(lo, hi) / 100.,
                        )
                        .map_err(|e| Error::Model(e.to_string()))?;
                        let path = {
                            let arena = m.arena();
                            vertices
                                .iter()
                                .map(|(fraction, volume)| {
                                    let linear = oximo::expr::extract_linear(&arena, volume.id())
                                        .ok_or_else(|| {
                                        Error::Model("storage path must be affine".into())
                                    })?;
                                    Ok((
                                        *fraction,
                                        crate::river::LinearForm {
                                            constant: linear.constant,
                                            coefficients: linear
                                                .coeffs
                                                .iter()
                                                .map(|(v, c)| (v.index(), *c))
                                                .collect(),
                                        },
                                    ))
                                })
                                .collect::<Result<Vec<_>, Error>>()?
                        };
                        laws.push(
                            crate::constitutive::StorageHeadEquation {
                                curve: curves[k].clone(),
                                vertices: path,
                                head: hm[[k, t]].var_id().unwrap().index(),
                                head_scale: 100.,
                            }
                            .into(),
                        );
                        Expr::from(hm[[k, t]] * 100.)
                    } else {
                        // A nondecreasing curve with equal endpoint levels is constant
                        // throughout the allowed storage range; no unknown head is needed.
                        zero() + lo
                    }
                }
            };
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
        .map(|r| r.level.trial_level(r.initial_storage))
        .chain(s.junctions.iter().map(|j| (j.head.min + j.head.max) / 2.))
        .chain(s.boundaries.iter().map(|b| b.head))
        .collect();
    for (e, tunnel) in s.tunnels.iter().enumerate() {
        let (source, target) = topology.tunnel_nodes[e];
        let dh = initial_head[source] - initial_head[target];
        let law = &topology.tunnel_laws[e];
        let guess = dh.signum() * (law.head_factor * dh.abs() / law.resistance).sqrt();
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
                let pressure = (head[t][source] - head[t][target]) * law.head_factor
                    - flow * flow.abs() * law.resistance;
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
        let mut powers: Vec<Expr<'_>> = Vec::with_capacity(nt);
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
            let previous = if t == 0 {
                u.operation.previous_flow().map(|v| zero() + v / flow_scale)
            } else {
                Some(q[[i, t - 1]].into())
            };
            if let Some(previous) = previous {
                let change = q[[i, t]] - previous;
                let dt = u.operation.ramp_hours(&case.grid, t);
                if let Some(up) = u.operation.flow_ramp.up {
                    constraint!(m, change <= up * dt / flow_scale);
                }
                if let Some(down) = u.operation.flow_ramp.down {
                    constraint!(m, -change <= down * dt / flow_scale);
                }
            }
            let flow: Expr<'_> = (q[[i, t]] * flow_scale).into();
            let net_head = head[t][source] - head[t][target] - flow.square() * u.loss_factor;
            let power: Expr<'_> = match &u.efficiency {
                crate::Efficiency::Quadratic {
                    peak,
                    best_flow,
                    flow_curvature,
                    best_head,
                    head_curvature,
                } => {
                    let eta = crate::physics::quadratic_efficiency(
                        *peak,
                        *best_flow,
                        *flow_curvature,
                        *best_head,
                        *head_curvature,
                        flow,
                        net_head,
                    );
                    crate::physics::generation(flow, net_head, eta)
                }
                crate::Efficiency::Table(_) => {
                    let k = turbine_index[i].unwrap();
                    let scale = u.power.max.max(1.);
                    if on {
                        m.set_bounds(ht[[k, t]], u.head.min / 100., u.head.max / 100.)
                            .map_err(|e| Error::Model(e.to_string()))?;
                        m.set_bounds(pt[[k, t]], u.power.min / scale, u.power.max / scale)
                            .map_err(|e| Error::Model(e.to_string()))?;
                        let h0 = (u.head.min + u.head.max) / 2.;
                        let q0 = start
                            .map_or((u.flow.min + u.flow.max) / 2., |p| p.flow[i][t])
                            .clamp(u.flow.min, u.flow.max);
                        m.set_initial(ht[[k, t]], h0 / 100.)
                            .map_err(|e| Error::Model(e.to_string()))?;
                        m.set_initial(
                            pt[[k, t]],
                            u.trial_generation(q0, h0).clamp(u.power.min, u.power.max) / scale,
                        )
                        .map_err(|e| Error::Model(e.to_string()))?;
                        constraint!(m, ht[[k, t]] - net_head / 100. == 0.);
                        laws.extend(crate::constitutive::Equation::turbine(
                            turbine_curves[k].clone(),
                            crate::river::LinearForm::variable(
                                q[[i, t]].var_id().unwrap().index(),
                                flow_scale,
                            ),
                            crate::river::LinearForm::variable(
                                ht[[k, t]].var_id().unwrap().index(),
                                100.,
                            ),
                            crate::river::LinearForm::variable(
                                pt[[k, t]].var_id().unwrap().index(),
                                scale,
                            ),
                            flow_scale,
                            scale,
                        ));
                    } else {
                        m.fix(ht[[k, t]], 0.)
                            .map_err(|e| Error::Model(e.to_string()))?;
                        m.fix(pt[[k, t]], 0.)
                            .map_err(|e| Error::Model(e.to_string()))?;
                    }
                    (pt[[k, t]] * scale).into()
                }
            };
            let (starts, stops) = u.operation.transitions(&commitment[i], t);
            let previous = if t == 0 {
                u.operation.previous_power().map(|p| zero() + p)
            } else {
                Some(powers[t - 1])
            };
            if let Some(previous) = previous {
                let dt = u.operation.ramp_hours(&case.grid, t);
                let scale = u.power.max.max(1.);
                let change = (power - previous) / scale;
                if let Some(up) = u.operation.power_ramp.up {
                    constraint!(
                        m,
                        change <= (up * dt + if starts { u.power.min } else { 0. }) / scale
                    );
                }
                if let Some(down) = u.operation.power_ramp.down {
                    constraint!(
                        m,
                        -change <= (down * dt + if stops { u.power.min } else { 0. }) / scale
                    );
                }
            }
            powers.push(power);
            if on && turbine_index[i].is_none() {
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
            let (start, stop) = u.operation.transitions(&commitment[i], t);
            revenue = revenue
                - if start { u.operation.startup_cost } else { 0. }
                - if stop { u.operation.shutdown_cost } else { 0. };
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
            + crate::nlp::linear_expression(&routing.inventory[r][nt], &river_columns, zero())
                * river.water_value;
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
        problem: nlp::prepare(&m, laws)?,
        storage,
        flow,
        tunnel_flow,
        junction_head,
        river_withdrawal,
        flow_scale,
    })
}
