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
    operating: &crate::operating_data::OperatingData,
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
            )
            .map_err(|e| Error::Model(e.to_string()))?;
            river_columns.push(Expr::from(qr[[k, t]] * flow_scale));
        }
    }
    variable!(m,rf[k in 0..routing.transfers.len()]);
    let mut reference_columns: Vec<_> = previous.map_or(vec![0.; routing.primary_columns], |p| {
        p.unit_flow
            .iter()
            .flatten()
            .copied()
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
    for e in 0..s.tunnels.len() {
        for t in 0..nt {
            let capacity = if operating.tunnels[e][t].opening == 0. {
                0.
            } else {
                operating.tunnels[e][t].capacity / flow_scale
            };
            m.set_bounds(qt[[e, t]], -capacity, capacity)
                .map_err(|e| Error::Model(e.to_string()))?;
        }
    }
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        let range = reservoir.storage.max - reservoir.storage.min;
        for (t, bounds) in operating.reservoirs[r].storage.iter().enumerate() {
            m.set_bounds(
                v[[r, t]],
                (bounds.min - reservoir.storage.min) / range,
                (bounds.max - reservoir.storage.min) / range,
            )
            .map_err(|e| Error::Model(e.to_string()))?;
        }
        m.fix(
            v[[r, 0]],
            (reservoir.initial_storage - reservoir.storage.min) / range,
        )
        .map_err(|e| Error::Model(e.to_string()))?;
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
    // Tangent hydraulic laws guide integer proposals. Nonlinear dispatch
    // subsequently enforces the original equations, not these tangents.
    let mut mean_head = vec![vec![zero(); nr]; nt];
    let mut reference_head = vec![vec![0.; nr]; nt];
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        let range = reservoir.storage.max - reservoir.storage.min;
        for t in 0..nt {
            let (a, b) = previous.map_or(
                (reservoir.initial_storage, reservoir.initial_storage),
                |p| (p.storage[r][t], p.storage[r][t + 1]),
            );
            let d = reservoir.level.mean_data(a, b);
            reference_head[t][r] = d.value;
            mean_head[t][r] = zero()
                + d.value
                + (v[[r, t]] * range + reservoir.storage.min - a) * d.gradient[0]
                + (v[[r, t + 1]] * range + reservoir.storage.min - b) * d.gradient[1];
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
            m.set_bounds(ga[[k, t]], gate.min, gate.max)
                .map_err(|e| Error::Model(e.to_string()))?;
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
                    m.fix(ga[[k, t]], 0.)
                        .map_err(|e| Error::Model(e.to_string()))?;
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
    let mut all_power = Vec::new();
    let mut all_flow = Vec::new();
    let mut all_on = Vec::new();
    let mut all_starts = Vec::new();
    let mut all_stops = Vec::new();
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
        let mut powers = Vec::new();
        let mut flows = Vec::new();
        let mut states = Vec::new();
        let mut starts = Vec::new();
        let mut stops = Vec::new();
        for t in 0..nt {
            let limits = &operating.units[g][t];
            if let Some(fixed) = case.commitment.get(&unit.name) {
                m.fix(on[[g, t]], if fixed[t] { 1. } else { 0. })
                    .map_err(|e| Error::Model(e.to_string()))?;
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
    let initial_value: f64 = s
        .reservoirs
        .iter()
        .map(|r| r.initial_storage * r.water_value)
        .sum::<f64>()
        + s.rivers
            .iter()
            .enumerate()
            .map(|(i, r)| routing.inventory[i][0].constant * r.water_value)
            .sum::<f64>();
    revenue = revenue - initial_value;
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
