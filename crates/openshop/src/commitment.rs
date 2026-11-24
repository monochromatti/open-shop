//! Unit and gate proposals with linearized hydraulics and exact water balance.
//! These schedules initialize nonlinear dispatch; their objectives are not bounds
//! on the nonlinear problem. No-good cuts let search change earlier unit and gate decisions.

use crate::{Error, Schedule};
use oximo::{ScipOptions, prelude::*, solvers::Scip};
use std::time::Duration;

#[derive(Clone, PartialEq, Eq)]
pub(crate) struct DiscreteControls {
    pub unit_on: Vec<Vec<bool>>,
    pub gate_position: Vec<Vec<usize>>,
    /// Only reaches with a positive on-state minimum, in prepared domain order.
    pub river_on: Vec<Vec<bool>>,
}
pub(crate) struct Proposal {
    pub controls: DiscreteControls,
    pub storage: Vec<Vec<f64>>,
    pub flow: Vec<Vec<f64>>,
    pub tunnel_flow: Vec<Vec<f64>>,
    pub river_withdrawal: Vec<Vec<f64>>,
    pub river_opening: Vec<Vec<f64>>,
    pub head: Vec<Vec<f64>>,
}

impl Proposal {
    /// Bounded physical control guesses for reconstruction. This initializes a
    /// nonlinear solve; acceptance still checks the original physical limits.
    pub fn physical_controls(&self, p: &crate::SchedulingProblem<'_>) -> crate::Controls {
        crate::Controls {
            commitment: self.controls.unit_on.clone(),
            unit_flow: self
                .flow
                .iter()
                .enumerate()
                .map(|(g, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(t, &q)| {
                            if self.controls.unit_on[g][t] {
                                let limits = p.operating.units[g][t].flow;
                                q.clamp(limits.min, limits.max)
                            } else {
                                0.
                            }
                        })
                        .collect()
                })
                .collect(),
            river_opening: self
                .river_opening
                .iter()
                .enumerate()
                .map(|(e, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(t, &a)| {
                            let limits = p.operating.river_gate[e][t];
                            a.clamp(limits.min, limits.max)
                        })
                        .collect()
                })
                .collect(),
            tunnel_opening: p
                .operating
                .tunnels
                .iter()
                .enumerate()
                .map(|(e, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(t, interval)| {
                            interval.positions[self.controls.gate_position[e][t]].opening
                        })
                        .collect()
                })
                .collect(),
        }
    }
    pub fn from_schedule(
        problem: &crate::SchedulingProblem<'_>,
        result: &Schedule,
    ) -> Result<Self, Error> {
        let gate_position = problem
            .operating
            .tunnels
            .iter()
            .enumerate()
            .map(|(e, row)| {
                row.iter()
                    .enumerate()
                    .map(|(t, v)| {
                        v.positions
                            .iter()
                            .position(|p| p.opening == result.tunnel_opening[e][t])
                    })
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| Error::Model("result contains an undeclared tunnel opening".into()))?;
        Ok(Self {
            controls: DiscreteControls {
                unit_on: result.commitment.clone(),
                gate_position,
                river_on: problem
                    .observations
                    .domains
                    .iter()
                    .map(|d| {
                        result.river_release[d.river]
                            .iter()
                            .map(|&q| problem.case.system.rivers[d.river].flow_on(q))
                            .collect()
                    })
                    .collect(),
            },
            storage: result.storage.clone(),
            flow: result.unit_flow.clone(),
            tunnel_flow: result.tunnel_flow.clone(),
            river_withdrawal: result.river_withdrawal.clone(),
            river_opening: result.river_opening.clone(),
            head: result.head.clone(),
        })
    }
}

pub(crate) fn propose(
    problem: &crate::SchedulingProblem<'_>,
    previous: Option<&Schedule>,
    excluded: &[DiscreteControls],
    seconds: f64,
) -> Result<Option<Proposal>, Error> {
    let (case, topology, routing, operating) = (
        problem.case,
        &problem.topology,
        &problem.routing,
        &problem.operating,
    );
    let s = &case.system;
    let (nr, nu, nt) = (s.reservoirs.len(), s.units.len(), case.prices.len());
    let m = Model::new("unit commitment");
    variable!(m,0.0<=v[r in 0..nr,t in 0..nt+1]<=1.0);
    variable!(m,q[g in 0..nu,t in 0..nt]>=0.0);
    variable!(m,qt[e in 0..s.tunnels.len(),t in 0..nt]);
    variable!(m,hj[j in 0..s.junctions.len(),t in 0..nt]);
    variable!(m,qr[e in 0..routing.outlets.len(),t in 0..nt]);
    variable!(m,0.<=ga[e in 0..routing.outlets.len(),t in 0..nt]<=1.);
    variable!(m,on[g in 0..nu,t in 0..nt],Bin);
    variable!(m,0.0<=start[g in 0..nu,t in 0..nt]<=1.0);
    variable!(m,0.0<=plant_on[p in 0..s.plants.len(),t in 0..nt]<=1.0);
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
    for e in 0..s.tunnels.len() {
        for t in 0..nt {
            river_columns.push((qt[[e, t]] * flow_scale).into());
        }
    }
    for (k, &e) in routing.outlets.iter().enumerate() {
        for t in 0..nt {
            let lower = if s.rivers[e].law.is_delta() {
                -operating.river_capacity[e][t] / flow_scale
            } else {
                0.
            };
            m.set_bounds(
                qr[[k, t]],
                lower,
                operating.river_capacity[e][t] / flow_scale,
            )?;
            river_columns.push(Expr::from(qr[[k, t]] * flow_scale));
        }
    }
    variable!(m,rf[k in 0..routing.transfers.len()]);
    let mut reference_columns: Vec<_> = previous.map_or(vec![0.; routing.primary_columns], |p| {
        p.unit_flow
            .iter()
            .flatten()
            .copied()
            .chain(p.tunnel_flow.iter().flatten().copied())
            .chain(
                routing
                    .outlets
                    .iter()
                    .flat_map(|&e| p.river_withdrawal[e].iter().copied()),
            )
            .collect()
    });
    for (k, e) in routing.transfers.iter().enumerate() {
        let q0 = e.flow.value(&reference_columns);
        let d = e.fraction.data(q0);
        reference_columns.push(d.value);
        let input = crate::nlp::linear_expression(&e.flow, &river_columns, zero());
        constraint!(
            m,
            rf[[k]] - (input * d.gradient[0] + d.value - q0 * d.gradient[0]) / flow_scale == 0.
        );
        river_columns.push((rf[[k]] * flow_scale).into());
    }
    problem.observations.constrain(&m, &river_columns, zero());
    let nd = problem.observations.domains.len();
    variable!(m, river_state[k in 0..nd,t in 0..nt], Bin);
    variable!(m, 0. <= river_start[k in 0..nd,t in 0..nt] <= 1.);
    let river_on = (0..nd)
        .map(|k| (0..nt).map(|t| Expr::from(river_state[[k, t]])).collect())
        .collect::<Vec<_>>();
    let river_starts = (0..nd)
        .map(|k| (0..nt).map(|t| Expr::from(river_start[[k, t]])).collect())
        .collect::<Vec<_>>();
    let river_transition_cost =
        problem
            .observations
            .constrain_states(&m, &river_columns, &river_on, &river_starts, zero());
    for (i, series) in routing.releases.iter().enumerate() {
        for (k, rate) in series.rates.iter().enumerate() {
            let period = case
                .grid
                .partition_point(|time| *time <= series.grid[k])
                .saturating_sub(1)
                .min(nt - 1);
            let release = crate::nlp::linear_expression(rate, &river_columns, zero());
            constraint!(
                m,
                release / flow_scale <= operating.river_capacity[i][period] / flow_scale
            );
        }
    }
    let choices: usize = operating
        .tunnels
        .iter()
        .flatten()
        .filter(|p| p.positions.len() > 1)
        .map(|p| p.positions.len())
        .sum();
    variable!(m,position[k in 0..choices],Bin);
    let mut gate_weights = Vec::new();
    let mut gate_opening = Vec::new();
    let mut gate_on = Vec::new();
    let mut index = 0;
    for (e, intervals) in operating.tunnels.iter().enumerate() {
        let mut weights_row = Vec::new();
        let mut opening_row = Vec::new();
        let mut on_row = Vec::new();
        for (t, interval) in intervals.iter().enumerate() {
            let weights: Vec<Expr<'_>> = if interval.positions.len() == 1 {
                vec![zero() + 1.]
            } else {
                let row = (index..index + interval.positions.len())
                    .map(|j| position[[j]].into())
                    .collect();
                index += interval.positions.len();
                row
            };
            let sum = weights.iter().fold(zero(), |sum, &w| sum + w);
            constraint!(m, sum == 1.);
            let open = interval
                .positions
                .iter()
                .zip(&weights)
                .filter(|(p, _)| p.opening > 0.)
                .fold(zero(), |sum, (_, w)| sum + *w);
            let opening = interval
                .positions
                .iter()
                .zip(&weights)
                .fold(zero(), |sum, (p, w)| sum + *w * p.opening);
            let cap = interval.capacity / flow_scale;
            m.set_bounds(
                qt[[e, t]],
                if topology.tunnel_outfalls[e].is_some() {
                    0.
                } else {
                    -cap
                },
                cap,
            )?;
            constraint!(m, qt[[e, t]] <= open * cap);
            if topology.tunnel_outfalls[e].is_none() {
                constraint!(m, qt[[e, t]] >= -open * cap);
            }
            on_row.push(open);
            weights_row.push(weights);
            opening_row.push(opening);
        }
        gate_on.push(on_row);
        gate_weights.push(weights_row);
        gate_opening.push(opening_row);
    }
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        let range = reservoir.storage.max - reservoir.storage.min;
        for (t, bounds) in operating.reservoirs[r].storage.iter().enumerate() {
            m.set_bounds(
                v[[r, t]],
                (bounds.min - reservoir.storage.min) / range,
                (bounds.max - reservoir.storage.min) / range,
            )?;
        }
        m.fix(
            v[[r, 0]],
            (reservoir.initial_storage - reservoir.storage.min) / range,
        )?;
        for t in 0..nt {
            let change = (v[[r, t + 1]] - v[[r, t]]) * range;
            if let Some(up) = operating.reservoirs[r].volume_change[t].up {
                constraint!(m, change <= up);
            }
            if let Some(down) = operating.reservoirs[r].volume_change[t].down {
                constraint!(m, -change <= down);
            }
            let mut net = zero() + operating.reservoirs[r].inflow[t];
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
                    constraint!(
                        m,
                        volume
                            >= (operating.reservoirs[r].interval_storage[t].min
                                - reservoir.storage.min)
                                / range
                    );
                    constraint!(
                        m,
                        volume
                            <= (operating.reservoirs[r].interval_storage[t].max
                                - reservoir.storage.min)
                                / range
                    );
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
    for (j, node) in s.junctions.iter().enumerate() {
        for t in 0..nt {
            m.set_bounds(hj[[j, t]], node.head.min / 100., node.head.max / 100.)?;
        }
    }
    // Tangent hydraulic laws guide integer proposals. Nonlinear dispatch
    // subsequently enforces the original equations, not these tangents.
    let mut mean_head = vec![vec![zero(); nr]; nt];
    let mut reference_head = vec![vec![0.; nr]; nt];
    let mut head_bounds = vec![vec![crate::Bounds::new(0., 0.); nr]; nt];
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        let range = reservoir.storage.max - reservoir.storage.min;
        for t in 0..nt {
            let (a, b) = previous.map_or(
                (reservoir.initial_storage, reservoir.initial_storage),
                |p| (p.storage[r][t], p.storage[r][t + 1]),
            );
            let d = reservoir.level.mean_data(a, b);
            let before = operating.reservoirs[r].storage[t];
            let after = operating.reservoirs[r].storage[t + 1];
            let da = [
                d.gradient[0] * (before.min - a),
                d.gradient[0] * (before.max - a),
            ];
            let db = [
                d.gradient[1] * (after.min - b),
                d.gradient[1] * (after.max - b),
            ];
            head_bounds[t][r] = crate::Bounds::new(
                d.value + da[0].min(da[1]) + db[0].min(db[1]),
                d.value + da[0].max(da[1]) + db[0].max(db[1]),
            );
            reference_head[t][r] = d.value;
            mean_head[t][r] = zero()
                + d.value
                + (v[[r, t]] * range + reservoir.storage.min - a) * d.gradient[0]
                + (v[[r, t + 1]] * range + reservoir.storage.min - b) * d.gradient[1];
        }
    }
    // Conditional pressure tangents give each positive gate position its own
    // hydraulic meaning in the MILP. Inactive row bounds follow the actual
    // affine head and discharge domains, without an arbitrary big-M constant.
    for (e, intervals) in operating.tunnels.iter().enumerate() {
        let (source, target) = topology.tunnel_nodes[e];
        for (t, interval) in intervals.iter().enumerate() {
            let level = |node: usize| {
                if node < nr {
                    mean_head[t][node]
                } else if node < nr + s.junctions.len() {
                    (hj[[node - nr, t]] * 100.).into()
                } else {
                    zero() + s.boundaries[node - nr - s.junctions.len()].head
                }
            };
            let bounds = |node: usize| {
                if node < nr {
                    head_bounds[t][node]
                } else if node < nr + s.junctions.len() {
                    s.junctions[node - nr].head
                } else {
                    let h = s.boundaries[node - nr - s.junctions.len()].head;
                    crate::Bounds::new(h, h)
                }
            };
            let (a, b) = (bounds(source), bounds(target));
            let delta = level(source) - level(target);
            let dh = previous.map_or(initial_head[source] - initial_head[target], |p| {
                p.head[source][t] - p.head[target][t]
            });
            for (j, law) in interval.positions.iter().enumerate() {
                if law.opening == 0. {
                    continue;
                }
                let q0 = previous
                    .filter(|p| p.tunnel_opening[e][t] == law.opening)
                    .map_or_else(
                        || dh.signum() * (dh.abs() / law.resistance).sqrt(),
                        |p| p.tunnel_flow[e][t],
                    )
                    .clamp(
                        if topology.tunnel_outfalls[e].is_some() {
                            0.
                        } else {
                            -interval.capacity
                        },
                        interval.capacity,
                    );
                let derivative = 2. * law.resistance * q0.abs();
                let constant = law.resistance * q0 * q0.abs();
                let residual = delta - qt[[e, t]] * (flow_scale * derivative) + constant;
                let lo = (a.min - b.max) - derivative * interval.capacity + constant;
                let hi = (a.max - b.min) + derivative * interval.capacity + constant;
                let selected = gate_weights[e][t][j];
                constraint!(m, residual / 100. >= lo.min(0.) / 100. * (1. - selected));
                constraint!(m, residual / 100. <= hi.max(0.) / 100. * (1. - selected));
            }
        }
    }
    for (k, &e) in routing.outlets.iter().enumerate() {
        let river = &s.rivers[e];
        let source = s
            .reservoirs
            .iter()
            .position(|r| r.name == river.source)
            .unwrap();
        for t in 0..nt {
            let gate = operating.river_gate[e][t];
            m.set_bounds(ga[[k, t]], gate.min, gate.max)?;
            let h0 = reference_head[t][source];
            let h = mean_head[t][source];
            let flow: Expr<'_> = match &river.law {
                crate::RiverLaw::Delta { curve, reference } => {
                    let target = s
                        .reservoirs
                        .iter()
                        .position(|r| r.name == river.target)
                        .unwrap();
                    let other0 = reference_head[t][target];
                    let other = mean_head[t][target];
                    let (ref0, href) = match reference {
                        crate::ReferenceLevel::Upstream => (h0, h),
                        crate::ReferenceLevel::Downstream => (other0, other),
                    };
                    let d = curve.data(h0 - other0, ref0);
                    m.fix(ga[[k, t]], 0.)?;
                    zero()
                        + d.value
                        + (h - other - (h0 - other0)) * d.gradient[0]
                        + (href - ref0) * d.gradient[1]
                }
                _ => {
                    let coefficient = river.law.outlet_flow(
                        h0,
                        1.,
                        operating.river_capacity[e][t],
                        operating.river_inflow[e][t],
                    )?;
                    let slope = match &river.law {
                        crate::RiverLaw::Orifice {
                            coefficient, crest, ..
                        } if h0 > *crest => coefficient / (2. * (h0 - crest).sqrt()),
                        crate::RiverLaw::Weir {
                            coefficient, crest, ..
                        } if h0 > *crest => 1.5 * coefficient * (h0 - crest).sqrt(),
                        crate::RiverLaw::Level { curve } => {
                            curve.mean_data(h0, h0).gradient.iter().sum()
                        }
                        _ => 0.,
                    };
                    let g0 = previous.map_or((gate.min + gate.max) / 2., |p| p.river_opening[e][t]);
                    ga[[k, t]] * coefficient + (h - h0) * (slope * g0)
                }
            };
            constraint!(m, qr[[k, t]] - flow / flow_scale == 0.);
        }
    }
    variable!(m,band_floor[k in 0..problem.level_windows.bands.len()]);
    for (k, band) in problem.level_windows.bands.iter().enumerate() {
        let i = band.reservoir;
        let r = &s.reservoirs[i];
        let range = r.storage.max - r.storage.min;
        m.set_bounds(
            band_floor[[k]],
            band.floor.min / 100.,
            band.floor.max / 100.,
        )?;
        for point in &band.points {
            let (t, fraction) = (point.interval, point.fraction);
            let volume = ((v[[i, t]] * (1. - fraction) + v[[i, t + 1]] * fraction) * range
                + r.storage.min)
                + crate::nlp::linear_expression(&point.correction, &river_columns, zero());
            let v0 = previous.map_or(r.initial_storage, |p| {
                p.storage[i][t] * (1. - fraction) + p.storage[i][t + 1] * fraction
            }) + point.correction.value(&reference_columns);
            let d = r.level.mean_data(v0, v0);
            let level = (volume - v0) * (d.gradient[0] + d.gradient[1]) + d.value;
            constraint!(m, level / 100. - band_floor[[k]] >= 0.);
            constraint!(m, level / 100. - band_floor[[k]] <= band.max_change / 100.);
        }
    }
    let mut all_power = Vec::new();
    let mut all_flow = Vec::new();
    let mut all_on = Vec::new();
    let mut all_starts = Vec::new();
    let mut all_stops = Vec::new();
    let mut revenue = zero();
    for (g, unit) in s.units.iter().enumerate() {
        let (a, b) = topology.unit_nodes[g];
        let mut powers = Vec::new();
        let mut flows = Vec::new();
        let mut states = Vec::new();
        let mut starts = Vec::new();
        let mut stops = Vec::new();
        for t in 0..nt {
            let limits = &operating.units[g][t];
            if let Some(fixed) = case.commitment.get(&unit.name) {
                m.fix(on[[g, t]], if fixed[t] { 1. } else { 0. })?;
            }
            let plant = &s.plants[topology.unit_plants[g]];
            let total_flow = topology.plant_units[topology.unit_plants[g]]
                .iter()
                .map(|&j| {
                    previous.map_or((s.units[j].flow.min + s.units[j].flow.max) / 2., |p| {
                        p.unit_flow[j][t]
                    })
                })
                .sum::<f64>();
            let outlet = previous
                .map_or(initial_head[b], |p| p.head[b][t])
                .max(plant.outlet_head_floor.unwrap_or(f64::NEG_INFINITY));
            let tailwater = plant
                .tailwater_curve
                .as_ref()
                .map_or(0., |c| c.extended_value(total_flow));
            let flow = previous
                .map_or((unit.flow.min + unit.flow.max) / 2., |p| p.unit_flow[g][t])
                .clamp(unit.flow.min, unit.flow.max);
            let head = ((previous.map_or(initial_head[a], |p| p.head[a][t])
                - outlet
                - tailwater
                - unit.loss_factor * flow * flow)
                * unit.mode.power_sign())
            .clamp(unit.head.min, unit.head.max);
            let eta = crate::machine::hydraulic_efficiency(unit, flow, head);
            let coefficient = match unit.mode {
                crate::UnitMode::Pump => 0.00981 * head / (eta * unit.motor_efficiency),
                crate::UnitMode::Generate => {
                    let p_ref = previous
                        .map_or((unit.power.min + unit.power.max) / 2., |p| {
                            p.unit_power[g][t]
                        })
                        .clamp(unit.power.min, unit.power.max);
                    let electrical = unit
                        .generator_efficiency_curve
                        .as_ref()
                        .map_or(1., |c| c.extended_value(p_ref));
                    0.00981 * head * eta * electrical
                }
            };
            if let Some(reference) = &topology.gross_heads[g] {
                let interval_level = |node: usize| {
                    if node < nr {
                        mean_head[t][node]
                    } else {
                        zero() + s.boundaries[node - nr - s.junctions.len()].head
                    }
                };
                let gross = (interval_level(reference.source) - interval_level(reference.target))
                    * unit.mode.power_sign();
                if let Some(lo) = limits.gross_head_min {
                    constraint!(
                        m,
                        gross / 100.
                            >= (lo - (lo - reference.range.min).max(0.) * (1. - on[[g, t]])) / 100.
                    );
                }
                if let Some(hi) = limits.gross_head_max {
                    constraint!(
                        m,
                        gross / 100.
                            <= (hi + (reference.range.max - hi).max(0.) * (1. - on[[g, t]])) / 100.
                    );
                }
                if let Some(lo) = limits.startup_gross_head_min {
                    let vertex_level = |node: usize| {
                        if node < nr {
                            let r = &s.reservoirs[node];
                            let volume = previous.map_or(r.initial_storage, |p| p.storage[node][t]);
                            let d = r.level.mean_data(volume, volume);
                            zero()
                                + d.value
                                + (v[[node, t]] * (r.storage.max - r.storage.min) + r.storage.min
                                    - volume)
                                    * (d.gradient[0] + d.gradient[1])
                        } else {
                            zero() + s.boundaries[node - nr - s.junctions.len()].head
                        }
                    };
                    let initial = (vertex_level(reference.source) - vertex_level(reference.target))
                        * unit.mode.power_sign();
                    constraint!(
                        m,
                        initial / 100.
                            >= (lo - (lo - reference.range.min).max(0.) * (1. - start[[g, t]]))
                                / 100.
                    );
                }
            }
            if let crate::Efficiency::Table(table) = &unit.efficiency {
                let minimum = table.envelope_data(head, false).value.max(unit.flow.min);
                let maximum = table.envelope_data(head, true).value.min(unit.flow.max);
                constraint!(m, q[[g, t]] >= on[[g, t]] * (minimum / flow_scale));
                constraint!(m, q[[g, t]] <= on[[g, t]] * (maximum / flow_scale));
            }
            constraint!(m, q[[g, t]] >= on[[g, t]] * (limits.flow.min / flow_scale));
            constraint!(m, q[[g, t]] <= on[[g, t]] * (limits.flow.max / flow_scale));
            constraint!(
                m,
                q[[g, t]] * (flow_scale * coefficient) >= on[[g, t]] * limits.power.min
            );
            constraint!(
                m,
                q[[g, t]] * (flow_scale * coefficient) <= on[[g, t]] * limits.power.max
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
            let power: Expr<'_> = (q[[g, t]] * (flow_scale * coefficient)).into();
            powers.push(power);
            flows.push((q[[g, t]] * flow_scale).into());
            states.push(on[[g, t]].into());
            starts.push(start[[g, t]].into());
            stops.push(stop);
            revenue = revenue
                + q[[g, t]]
                    * (flow_scale
                        * coefficient
                        * unit.mode.power_sign()
                        * case.prices[t]
                        * (case.grid[t + 1] - case.grid[t]))
                - start[[g, t]] * limits.startup_cost
                - stop * limits.shutdown_cost;
        }
        all_power.push(powers);
        all_flow.push(flows);
        all_on.push(states);
        all_starts.push(starts);
        all_stops.push(stops);
    }
    let mut plant_states = Vec::new();
    for (p, members) in topology.plant_units.iter().enumerate() {
        let mut states = Vec::new();
        for t in 0..nt {
            let sum = members.iter().fold(zero(), |sum, &g| sum + on[[g, t]]);
            for &g in members {
                constraint!(m, plant_on[[p, t]] >= on[[g, t]]);
            }
            constraint!(m, plant_on[[p, t]] <= sum);
            states.push(plant_on[[p, t]].into());
        }
        plant_states.push(states);
    }
    problem
        .interlocks
        .constrain(&m, &all_on, &gate_on, zero(), nt);
    crate::operating_equations::constrain(
        &m,
        case,
        topology,
        operating,
        &crate::operating_equations::OperatingExpressions {
            power: all_power,
            flow: all_flow,
            on: all_on,
            starts: all_starts,
            stops: all_stops,
            plant_on: plant_states,
        },
    );
    for commitment in excluded {
        let mut difference = zero();
        for g in 0..nu {
            for t in 0..nt {
                difference = if commitment.unit_on[g][t] {
                    difference + 1. - on[[g, t]]
                } else {
                    difference + on[[g, t]]
                };
            }
        }
        for (e, row) in gate_weights.iter().enumerate() {
            for (t, weights) in row.iter().enumerate() {
                difference = difference + 1. - weights[commitment.gate_position[e][t]];
            }
        }
        for (k, row) in commitment.river_on.iter().enumerate() {
            for (t, &state) in row.iter().enumerate() {
                difference = difference
                    + if state {
                        1. - river_state[[k, t]]
                    } else {
                        river_state[[k, t]]
                    };
            }
        }
        constraint!(m, difference >= 1.);
    }
    let moving: Vec<_> = s
        .tunnels
        .iter()
        .enumerate()
        .filter_map(|(e, t)| (t.gate_adjustment_cost > 0.).then_some(e))
        .collect();
    variable!(m,0.<=movement[k in 0..moving.len(),t in 0..nt]<=1.);
    for (k, &e) in moving.iter().enumerate() {
        for t in 0..nt {
            let prior = if t == 0 {
                zero() + s.tunnels[e].initial_opening.unwrap_or(s.tunnels[e].opening)
            } else {
                gate_opening[e][t - 1]
            };
            let change = gate_opening[e][t] - prior;
            constraint!(m, movement[[k, t]] >= change);
            constraint!(m, movement[[k, t]] >= -change);
            revenue = revenue - movement[[k, t]] * s.tunnels[e].gate_adjustment_cost;
        }
    }
    revenue = revenue + crate::nlp::linear_expression(&problem.water_value, &river_columns, zero());
    revenue = revenue - problem.observations.cost(&m, &river_columns, zero());
    revenue = revenue - river_transition_cost;
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
    let mut river_opening = vec![vec![0.; nt]; s.rivers.len()];
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
            river_opening[r][t] = value(ga[[k, t]])?;
        }
    }
    let mut index = 0;
    let gate_position = operating
        .tunnels
        .iter()
        .map(|row| {
            row.iter()
                .map(|interval| {
                    if interval.positions.len() == 1 {
                        return Ok(0);
                    }
                    let values = (index..index + interval.positions.len())
                        .map(|j| value(position[[j]]))
                        .collect::<Result<Vec<_>, _>>()?;
                    index += interval.positions.len();
                    Ok(values
                        .iter()
                        .enumerate()
                        .max_by(|(_, a), (_, b)| a.total_cmp(b))
                        .unwrap()
                        .0)
                })
                .collect::<Result<Vec<_>, Error>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let junction_head = (0..s.junctions.len())
        .map(|j| {
            (0..nt)
                .map(|t| value(hj[[j, t]]).map(|h| 100. * h))
                .collect::<Result<Vec<_>, Error>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut columns: Vec<_> = flow
        .iter()
        .flatten()
        .copied()
        .chain(tunnel_flow.iter().flatten().copied())
        .chain(
            routing
                .outlets
                .iter()
                .flat_map(|&r| river_withdrawal[r].iter().copied()),
        )
        .collect();
    routing.complete_columns(&mut columns);
    let head = crate::trajectory::hydraulic_heads(case, routing, &storage, &columns, junction_head);
    Ok(Some(Proposal {
        controls: DiscreteControls {
            unit_on: commitment,
            gate_position,
            river_on: (0..nd)
                .map(|k| {
                    (0..nt)
                        .map(|t| value(river_state[[k, t]]).map(|v| v > 0.5))
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Vec<_>, _>>()?,
        },
        storage,
        flow,
        tunnel_flow,
        river_withdrawal,
        river_opening,
        head,
    }))
}
