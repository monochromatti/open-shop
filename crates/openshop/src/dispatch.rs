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
    variable!(m,ht[g in 0..nu,t in 0..nt]);
    variable!(m,pt[g in 0..nu,t in 0..nt]);
    let machines: Vec<_> = s.units.iter().cloned().map(std::sync::Arc::new).collect();
    let curved_plants: Vec<_> = s
        .plants
        .iter()
        .enumerate()
        .filter_map(|(i, p)| p.tailwater_curve.as_ref().map(|_| i))
        .collect();
    variable!(m,hw[k in 0..curved_plants.len(),t in 0..nt]);
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
                rate - net == 0.0
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
    let mut tailwater = vec![vec![zero(); nt]; s.plants.len()];
    for (k, &p) in curved_plants.iter().enumerate() {
        let curve = std::sync::Arc::new(s.plants[p].tailwater_curve.as_ref().unwrap().clone());
        for t in 0..nt {
            let mut total = crate::river::LinearForm::default();
            for &g in &topology.plant_units[p] {
                total
                    .coefficients
                    .insert(q[[g, t]].var_id().unwrap().index(), flow_scale);
            }
            m.set_initial(hw[[k, t]], curve.values()[0] / 100.)
                .map_err(|e| Error::Model(e.to_string()))?;
            laws.push(crate::constitutive::Equation {
                linear: crate::river::LinearForm::variable(
                    hw[[k, t]].var_id().unwrap().index(),
                    -1.,
                ),
                terms: vec![crate::constitutive::Term {
                    function: crate::constitutive::PhysicalFunction::Table(curve.clone()),
                    inputs: [total, Default::default(), Default::default()],
                    weight: 1. / 100.,
                }],
                bounds: crate::Bounds::new(0., 0.),
            });
            tailwater[p][t] = (hw[[k, t]] * 100.).into();
        }
    }
    let mut all_power = Vec::new();
    let mut all_flow = Vec::new();
    let mut all_on = Vec::new();
    let mut all_starts = Vec::new();
    let mut all_stops = Vec::new();
    let mut revenue = zero();
    for (i, u) in s.units.iter().enumerate() {
        let (source, target) = topology.unit_nodes[i];
        let mut powers = Vec::new();
        let mut flows = Vec::new();
        let mut states = Vec::new();
        let mut starts = Vec::new();
        let mut stops = Vec::new();
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
            let flow: Expr<'_> = (q[[i, t]] * flow_scale).into();
            let plant = topology.unit_plants[i];
            let p = &s.plants[plant];
            let outlet = p.outlet_head_floor.map_or(head[t][target], |floor| {
                (head[t][target] + floor + (head[t][target] - floor).abs()) * 0.5
            });
            let net_head =
                (head[t][source] - outlet - tailwater[plant][t] - flow.square() * u.loss_factor)
                    * u.mode.power_sign();
            let scale = u.power.max.max(1.);
            if on {
                m.set_bounds(ht[[i, t]], u.head.min / 100., u.head.max / 100.)
                    .map_err(|e| Error::Model(e.to_string()))?;
                m.set_bounds(pt[[i, t]], u.power.min / scale, u.power.max / scale)
                    .map_err(|e| Error::Model(e.to_string()))?;
                let h0 = ((initial_head[source] - initial_head[target]) * u.mode.power_sign())
                    .clamp(u.head.min, u.head.max);
                let q0 = start
                    .map_or((u.flow.min + u.flow.max) / 2., |p| p.flow[i][t])
                    .clamp(u.flow.min, u.flow.max);
                m.set_initial(ht[[i, t]], h0 / 100.)
                    .map_err(|e| Error::Model(e.to_string()))?;
                let guess = u.trial_power(q0, h0);
                let p0 = if guess.is_finite() {
                    guess.clamp(u.power.min, u.power.max)
                } else {
                    (u.power.min + u.power.max) / 2.
                };
                m.set_initial(pt[[i, t]], p0 / scale)
                    .map_err(|e| Error::Model(e.to_string()))?;
                constraint!(m, ht[[i, t]] - net_head / 100. == 0.);
                laws.extend(crate::constitutive::Equation::production(
                    machines[i].clone(),
                    crate::river::LinearForm::variable(
                        q[[i, t]].var_id().unwrap().index(),
                        flow_scale,
                    ),
                    crate::river::LinearForm::variable(ht[[i, t]].var_id().unwrap().index(), 100.),
                    crate::river::LinearForm::variable(pt[[i, t]].var_id().unwrap().index(), scale),
                    flow_scale,
                    scale,
                ));
            } else {
                m.fix(ht[[i, t]], 0.)
                    .map_err(|e| Error::Model(e.to_string()))?;
                m.fix(pt[[i, t]], 0.)
                    .map_err(|e| Error::Model(e.to_string()))?;
            }
            let power: Expr<'_> = (pt[[i, t]] * scale).into();
            let (starts_now, stops_now) = u.operation.transitions(&commitment[i], t);
            powers.push(power);
            flows.push(flow);
            states.push(zero() + if on { 1. } else { 0. });
            starts.push(zero() + if starts_now { 1. } else { 0. });
            stops.push(zero() + if stops_now { 1. } else { 0. });
            revenue = revenue
                + power
                    * (u.mode.power_sign() * case.prices[t] * (case.grid[t + 1] - case.grid[t]));
            let (start, stop) = u.operation.transitions(&commitment[i], t);
            revenue = revenue
                - if start { u.operation.startup_cost } else { 0. }
                - if stop { u.operation.shutdown_cost } else { 0. };
        }
        all_power.push(powers);
        all_flow.push(flows);
        all_on.push(states);
        all_starts.push(starts);
        all_stops.push(stops);
    }
    let plant_on: Vec<Vec<_>> = topology
        .plant_units
        .iter()
        .map(|members| {
            (0..nt)
                .map(|t| {
                    zero()
                        + if members.iter().any(|&g| commitment[g][t]) {
                            1.
                        } else {
                            0.
                        }
                })
                .collect()
        })
        .collect();
    crate::operating_equations::constrain(
        &m,
        case,
        topology,
        &crate::operating_equations::OperatingExpressions {
            power: all_power,
            flow: all_flow,
            on: all_on,
            starts: all_starts,
            stops: all_stops,
            plant_on,
        },
    );
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
