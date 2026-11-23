//! Equipment limits at hydraulic vertices, with controls held over each
//! operating interval. Instantaneous pressure equations do not advance storage.

use crate::{
    Bounds, Error, Schedule, SchedulingProblem,
    constitutive::{Equation, PhysicalFunction, Term},
    nlp,
    river::{LinearForm, Routing},
};
use oximo::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
pub(crate) struct OperatingPoint {
    pub interval: usize,
    pub fraction: f64,
    correction: Vec<LinearForm>,
}
#[derive(Clone)]
pub(crate) struct OperatingPoints {
    points: Vec<OperatingPoint>,
    required: Vec<[bool; 4]>,
}

impl OperatingPoint {
    fn new(routing: &Routing, interval: usize, fraction: f64) -> Self {
        let correction = routing
            .storage_shape
            .iter()
            .map(|r| {
                let path = &r[interval];
                let j = path
                    .partition_point(|v| v.0 <= fraction)
                    .saturating_sub(1)
                    .min(path.len() - 2);
                let w = (fraction - path[j].0) / (path[j + 1].0 - path[j].0);
                let mut value = path[j].1.clone().scaled(1. - w);
                value.add_scaled(&path[j + 1].1, w);
                value
            })
            .collect();
        Self {
            interval,
            fraction,
            correction,
        }
    }
}

impl OperatingPoints {
    pub fn new(case: &crate::ScheduleCase, routing: &Routing) -> Self {
        let mut points = Vec::new();
        for t in 0..case.prices.len() {
            let mut fractions: Vec<_> = routing
                .storage_shape
                .iter()
                .flat_map(|r| r[t].iter().map(|p| p.0))
                .chain([0., 1.])
                .collect();
            fractions.sort_by(f64::total_cmp);
            fractions.dedup();
            for fraction in fractions {
                points.push(OperatingPoint::new(routing, t, fraction));
            }
        }
        Self {
            points,
            required: case
                .system
                .units
                .iter()
                .map(crate::head_limits::required_constraints)
                .collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Add times at which storage paths cross measured level knots or an
    /// outlet-head floor. These pressure observations do not advance storage.
    pub fn with_level_crossings(
        &self,
        p: &SchedulingProblem<'_>,
        storage: &[Vec<f64>],
        columns: &[f64],
    ) -> Self {
        let mut fractions = vec![Vec::new(); p.case.prices.len()];
        for point in &self.points {
            fractions[point.interval].push(point.fraction);
        }
        for (r, reservoir) in p.case.system.reservoirs.iter().enumerate() {
            let mut cuts = match &reservoir.level {
                crate::StorageCurve::Table(curve) => curve.coordinates().to_vec(),
                _ => Vec::new(),
            };
            for plant in &p.case.system.plants {
                if plant.target != reservoir.name {
                    continue;
                }
                let Some(floor) = plant.outlet_head_floor else {
                    continue;
                };
                let (mut lo, mut hi) = (reservoir.storage.min, reservoir.storage.max);
                if !(reservoir.level.trial_level(lo) < floor
                    && floor < reservoir.level.trial_level(hi))
                {
                    continue;
                }
                // Monotonicity is validated at preparation. Bisection stops
                // when no representable storage lies between its endpoints.
                loop {
                    let middle = lo / 2. + hi / 2.;
                    if !(lo < middle && middle < hi) {
                        break;
                    }
                    if reservoir.level.trial_level(middle) < floor {
                        lo = middle;
                    } else {
                        hi = middle;
                    }
                }
                cuts.push(lo / 2. + hi / 2.);
            }
            for (t, path) in p.routing.storage_shape[r].iter().enumerate() {
                let value = |v: &(f64, LinearForm)| {
                    storage[r][t] * (1. - v.0) + storage[r][t + 1] * v.0 + v.1.value(columns)
                };
                for segment in path.windows(2) {
                    let a = value(&segment[0]);
                    let b = value(&segment[1]);
                    for &knot in &cuts {
                        if a.min(b) < knot && knot < a.max(b) {
                            fractions[t].push(
                                segment[0].0
                                    + (segment[1].0 - segment[0].0) * ((knot - a) / (b - a)),
                            );
                        }
                    }
                }
            }
        }
        let mut points = Vec::new();
        for (t, fractions) in fractions.iter_mut().enumerate() {
            fractions.sort_by(f64::total_cmp);
            fractions.dedup_by(|a, b| (*a - *b).abs() <= 32. * f64::EPSILON);
            points.extend(
                fractions
                    .iter()
                    .map(|&f| OperatingPoint::new(&p.routing, t, f)),
            );
        }
        Self {
            points,
            required: self.required.clone(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn constrain(
        &self,
        p: &SchedulingProblem<'_>,
        m: &Model,
        storage: &[Vec<Expr<'_>>],
        flow: &[Vec<Expr<'_>>],
        columns: &[Expr<'_>],
        controls: &crate::commitment::DiscreteControls,
        start: Option<&crate::commitment::Proposal>,
        laws: &mut Vec<Equation>,
    ) -> Result<(), Error> {
        let s = &p.case.system;
        let (nr, nj, nu, np) = (
            s.reservoirs.len(),
            s.junctions.len(),
            s.units.len(),
            self.points.len(),
        );
        let mzero = storage[0][0] * 0.;
        variable!(m,hp[k in 0..np,n in 0..nr+nj]);
        variable!(m,tq[k in 0..np,e in 0..s.tunnels.len()]);
        variable!(m,uq[k in 0..np,g in 0..nu]);
        variable!(m,uh[k in 0..np,g in 0..nu]);
        variable!(m,up[k in 0..np,g in 0..nu]);
        variable!(m,tw[k in 0..np,j in 0..s.plants.len()]);
        let initial_controls = start.map(|v| crate::Controls {
            commitment: controls.unit_on.clone(),
            unit_flow: v
                .flow
                .iter()
                .enumerate()
                .map(|(g, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(t, &q)| if controls.unit_on[g][t] { q } else { 0. })
                        .collect()
                })
                .collect(),
            river_opening: v.river_opening.clone(),
            tunnel_opening: p
                .operating
                .tunnels
                .iter()
                .enumerate()
                .map(|(e, row)| {
                    row.iter()
                        .enumerate()
                        .map(|(t, v)| v.positions[controls.gate_position[e][t]].opening)
                        .collect()
                })
                .collect(),
        });
        let initial_columns = start.map(|v| {
            let mut columns: Vec<_> = v
                .flow
                .iter()
                .flatten()
                .copied()
                .chain(v.tunnel_flow.iter().flatten().copied())
                .chain(
                    p.routing
                        .outlets
                        .iter()
                        .flat_map(|&e| v.river_withdrawal[e].iter().copied()),
                )
                .collect();
            p.routing.complete_columns(&mut columns);
            columns
        });
        for (k, point) in self.points.iter().enumerate() {
            let (t, f) = (point.interval, point.fraction);
            let initial = if let (Some(start), Some(controls), Some(columns)) =
                (start, initial_controls.as_ref(), initial_columns.as_ref())
            {
                let volume: Vec<_> = start
                    .storage
                    .iter()
                    .enumerate()
                    .map(|(r, v)| {
                        v[t] * (1. - f) + v[t + 1] * f + point.correction[r].value(columns)
                    })
                    .collect();
                let reference = crate::hydraulics::HydraulicState {
                    storage: volume.clone(),
                    head: start.head.iter().map(|r| r[t]).collect(),
                    tunnel_flow: start.tunnel_flow.iter().map(|r| r[t]).collect(),
                    river_withdrawal: p
                        .routing
                        .outlets
                        .iter()
                        .map(|&e| start.river_withdrawal[e][t])
                        .collect(),
                    unit_flow: start.flow.iter().map(|r| r[t]).collect(),
                    iterations: 0,
                };
                match crate::hydraulics::hydraulic_step(
                    p,
                    controls,
                    t,
                    &volume,
                    0.,
                    &columns[..p.routing.primary_columns],
                    Some(&reference),
                ) {
                    Ok(state) => Some(state),
                    Err(Error::Solver(_)) => None,
                    Err(error) => return Err(error),
                }
            } else {
                None
            };
            let mut head = Vec::new();
            for (r, reservoir) in s.reservoirs.iter().enumerate() {
                let lo = reservoir.level.trial_level(reservoir.storage.min);
                let hi = reservoir.level.trial_level(reservoir.storage.max);
                m.set_bounds(hp[[k, r]], lo / 100., hi / 100.)
                    .map_err(|e| Error::Model(e.to_string()))?;
                let v0 = start.map_or(reservoir.initial_storage, |v| {
                    v.storage[r][t] * (1. - f) + v.storage[r][t + 1] * f
                });
                m.set_initial(
                    hp[[k, r]],
                    initial
                        .as_ref()
                        .map_or_else(|| reservoir.level.trial_level(v0), |v| v.head[r])
                        .clamp(lo, hi)
                        / 100.,
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                if lo < hi {
                    let volume = storage[r][t] * (1. - f)
                        + storage[r][t + 1] * f
                        + nlp::linear_expression(&point.correction[r], columns, mzero);
                    let v = nlp::linear_form(m, volume)?;
                    laws.push(Equation {
                        linear: LinearForm::variable(hp[[k, r]].var_id().unwrap().index(), -1.),
                        terms: vec![Term {
                            function: PhysicalFunction::StorageMean(Arc::new(
                                reservoir.level.clone(),
                            )),
                            inputs: [v.clone(), v, Default::default()],
                            weight: 1. / 100.,
                        }],
                        bounds: Bounds::new(0., 0.),
                    });
                }
                head.push(Expr::from(hp[[k, r]] * 100.));
            }
            for (j, junction) in s.junctions.iter().enumerate() {
                m.set_bounds(
                    hp[[k, nr + j]],
                    junction.head.min / 100.,
                    junction.head.max / 100.,
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                m.set_initial(
                    hp[[k, nr + j]],
                    initial
                        .as_ref()
                        .map_or_else(
                            || {
                                start.map_or((junction.head.min + junction.head.max) / 2., |v| {
                                    v.head[nr + j][t]
                                })
                            },
                            |v| v.head[nr + j],
                        )
                        .clamp(junction.head.min, junction.head.max)
                        / 100.,
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                head.push((hp[[k, nr + j]] * 100.).into());
            }
            head.extend(s.boundaries.iter().map(|b| mzero + b.head));
            let mut unit_flow = Vec::new();
            for (g, u) in s.units.iter().enumerate() {
                let on = controls.unit_on[g][t];
                let limits = &p.operating.units[g][t];
                let scale = u.flow.max.max(1.);
                m.set_bounds(
                    uq[[k, g]],
                    if on { limits.flow.min / scale } else { 0. },
                    if on { limits.flow.max / scale } else { 0. },
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                m.set_initial(
                    uq[[k, g]],
                    if on {
                        initial
                            .as_ref()
                            .map_or_else(
                                || {
                                    start.map_or((limits.flow.min + limits.flow.max) / 2., |v| {
                                        v.flow[g][t]
                                    })
                                },
                                |v| v.unit_flow[g],
                            )
                            .clamp(limits.flow.min, limits.flow.max)
                            / scale
                    } else {
                        0.
                    },
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                let q: Expr<'_> = (uq[[k, g]] * scale).into();
                if u.pump_characteristic.is_none() {
                    constraint!(m, (q - flow[g][t]) / scale == 0.);
                }
                unit_flow.push(q);
            }
            for (e, &(a, b)) in p.topology.tunnel_nodes.iter().enumerate() {
                let interval = &p.operating.tunnels[e][t];
                let law = &interval.positions[controls.gate_position[e][t]];
                let scale = s.tunnels[e].capacity.max(1.);
                let cap = if law.opening == 0. {
                    0.
                } else {
                    interval.capacity / scale
                };
                m.set_bounds(
                    tq[[k, e]],
                    if p.topology.tunnel_outfalls[e].is_some() {
                        0.
                    } else {
                        -cap
                    },
                    cap,
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                m.set_initial(
                    tq[[k, e]],
                    initial
                        .as_ref()
                        .map_or_else(
                            || start.map_or(0., |v| v.tunnel_flow[e][t]),
                            |v| v.tunnel_flow[e],
                        )
                        .clamp(
                            if p.topology.tunnel_outfalls[e].is_some() {
                                0.
                            } else {
                                -cap * scale
                            },
                            cap * scale,
                        )
                        / scale,
                )
                .map_err(|e| Error::Model(e.to_string()))?;
                if law.opening > 0. {
                    let q = tq[[k, e]] * scale;
                    constraint!(
                        m,
                        (head[a] - head[b] - q * q.abs() * law.resistance) / 100. == 0.
                    );
                }
            }
            for node in nr..nr + nj {
                let mut net = mzero;
                let mut scale: f64 = 1.;
                for &(g, sign) in &p.topology.incidence[node] {
                    net = net + unit_flow[g] * sign;
                    scale = scale.max(s.units[g].flow.max);
                }
                for &(e, sign) in &p.topology.tunnel_incidence[node] {
                    let qscale = s.tunnels[e].capacity.max(1.);
                    net = net + tq[[k, e]] * (sign * qscale);
                    scale = scale.max(qscale);
                }
                constraint!(m, net / scale == 0.);
            }
            let mut power: Vec<Expr<'_>> = Vec::new();
            for (j, plant) in s.plants.iter().enumerate() {
                if let Some(curve) = &plant.tailwater_curve {
                    let q = p.topology.plant_units[j]
                        .iter()
                        .fold(mzero, |a, &g| a + unit_flow[g]);
                    laws.push(Equation {
                        linear: LinearForm::variable(tw[[k, j]].var_id().unwrap().index(), -1.),
                        terms: vec![Term {
                            function: PhysicalFunction::Table(Arc::new(curve.clone())),
                            inputs: [
                                nlp::linear_form(m, q)?,
                                Default::default(),
                                Default::default(),
                            ],
                            weight: 1. / 100.,
                        }],
                        bounds: Bounds::new(0., 0.),
                    });
                    let q0 = p.topology.plant_units[j]
                        .iter()
                        .map(|&g| {
                            start.map_or((s.units[g].flow.min + s.units[g].flow.max) / 2., |v| {
                                v.flow[g][t]
                            })
                        })
                        .sum();
                    m.set_initial(tw[[k, j]], curve.extended_value(q0) / 100.)
                        .map_err(|e| Error::Model(e.to_string()))?;
                } else {
                    m.fix(tw[[k, j]], 0.)
                        .map_err(|e| Error::Model(e.to_string()))?;
                }
            }
            for (g, u) in s.units.iter().enumerate() {
                let limits = &p.operating.units[g][t];
                let scale = u.power.max.max(1.);
                if controls.unit_on[g][t] {
                    let plant = p.topology.unit_plants[g];
                    let spec = &s.plants[plant];
                    let (a, b) = p.topology.unit_nodes[g];
                    let outlet = spec.outlet_head_floor.map_or(head[b], |floor| {
                        (head[b] + floor + (head[b] - floor).abs()) * 0.5
                    });
                    let q = unit_flow[g];
                    let h = (head[a] - outlet - tw[[k, plant]] * 100. - q.square() * u.loss_factor)
                        * u.mode.power_sign();
                    m.set_bounds(uh[[k, g]], u.head.min / 100., u.head.max / 100.)
                        .map_err(|e| Error::Model(e.to_string()))?;
                    m.set_bounds(
                        up[[k, g]],
                        limits.power.min / scale,
                        limits.power.max / scale,
                    )
                    .map_err(|e| Error::Model(e.to_string()))?;
                    let h0 = initial
                        .as_ref()
                        .map_or((u.head.min + u.head.max) / 2., |v| {
                            crate::hydraulics::machine_head(p, g, &v.head, &v.unit_flow)
                        })
                        .clamp(u.head.min, u.head.max);
                    let q0 = initial.as_ref().map_or_else(
                        || start.map_or((u.flow.min + u.flow.max) / 2., |v| v.flow[g][t]),
                        |v| v.unit_flow[g],
                    );
                    m.set_initial(uh[[k, g]], h0 / 100.)
                        .map_err(|e| Error::Model(e.to_string()))?;
                    let guess = u.trial_power(q0, h0);
                    let p0 = if guess.is_finite() {
                        guess.clamp(limits.power.min, limits.power.max)
                    } else {
                        (limits.power.min + limits.power.max) / 2.
                    };
                    m.set_initial(up[[k, g]], p0 / scale)
                        .map_err(|e| Error::Model(e.to_string()))?;
                    constraint!(m, uh[[k, g]] - h / 100. == 0.);
                    laws.extend(Equation::production(
                        Arc::new(u.clone()),
                        nlp::linear_form(m, q)?,
                        LinearForm::variable(uh[[k, g]].var_id().unwrap().index(), 100.),
                        LinearForm::variable(up[[k, g]].var_id().unwrap().index(), scale),
                        u.flow.max.max(1.),
                        scale,
                    ));
                } else {
                    m.fix(uh[[k, g]], 0.)
                        .map_err(|e| Error::Model(e.to_string()))?;
                    m.fix(up[[k, g]], 0.)
                        .map_err(|e| Error::Model(e.to_string()))?;
                }
                power.push(Expr::from(up[[k, g]] * scale));
            }
            for (j, members) in p.topology.plant_units.iter().enumerate() {
                let limits = &p.operating.plants[j][t];
                let on = members.iter().any(|&g| controls.unit_on[g][t]);
                let total_q = members.iter().fold(mzero, |v, &g| v + unit_flow[g]);
                let total_p = members.iter().fold(mzero, |v, &g| v + power[g]);
                let qs = members
                    .iter()
                    .map(|&g| s.units[g].flow.max)
                    .sum::<f64>()
                    .max(1.);
                let ps = members
                    .iter()
                    .map(|&g| s.units[g].power.max)
                    .sum::<f64>()
                    .max(1.);
                constraint!(
                    m,
                    total_q / qs >= if on { limits.flow.min / qs } else { 0. }
                );
                constraint!(m, total_q / qs <= limits.flow.max / qs);
                constraint!(
                    m,
                    total_p / ps >= if on { limits.power.min / ps } else { 0. }
                );
                constraint!(m, total_p / ps <= limits.power.max / ps);
            }
        }
        for k in 0..self.points.len().saturating_sub(1) {
            let t = self.points[k].interval;
            if self.points[k + 1].interval != t {
                continue;
            }
            for (g, u) in s.units.iter().enumerate() {
                if !controls.unit_on[g][t] {
                    continue;
                }
                let limits = &p.operating.units[g][t];
                let envelope = Arc::new(crate::head_limits::HeadEnvelope {
                    unit: Arc::new(u.clone()),
                    power: limits.power,
                    flow: limits.flow,
                });
                use crate::head_limits::HeadLimit::*;
                for (index, (component, scale)) in [
                    (UpperPower, u.power.max.max(1.)),
                    (LowerPower, u.power.max.max(1.)),
                    (Efficiency, 1.),
                    (Flow, u.flow.max.max(1.)),
                ]
                .into_iter()
                .enumerate()
                {
                    if !self.required[g][index] {
                        continue;
                    }
                    laws.push(Equation {
                        linear: Default::default(),
                        terms: vec![Term {
                            function: PhysicalFunction::HeadSpan(envelope.clone(), component),
                            inputs: [
                                nlp::linear_form(m, flow[g][t])?,
                                LinearForm::variable(uh[[k, g]].var_id().unwrap().index(), 100.),
                                LinearForm::variable(
                                    uh[[k + 1, g]].var_id().unwrap().index(),
                                    100.,
                                ),
                            ],
                            weight: 1. / scale,
                        }],
                        bounds: Bounds::new(f64::NEG_INFINITY, 0.),
                    });
                }
            }
        }
        Ok(())
    }

    pub fn residual(&self, p: &SchedulingProblem<'_>, schedule: &Schedule) -> Result<f64, Error> {
        let s = &p.case.system;
        let controls = schedule.controls();
        let mut columns: Vec<_> = schedule
            .unit_flow
            .iter()
            .flatten()
            .copied()
            .chain(schedule.tunnel_flow.iter().flatten().copied())
            .chain(
                p.routing
                    .outlets
                    .iter()
                    .flat_map(|&e| schedule.river_withdrawal[e].iter().copied()),
            )
            .collect();
        p.routing.complete_columns(&mut columns);
        let mut excess: f64 = 0.;
        let violation = |x: f64, lo: f64, hi: f64, scale: f64| {
            if x.is_finite() {
                ((lo - x).max(x - hi)).max(0.) / scale
            } else {
                f64::INFINITY
            }
        };
        let mut observed_heads = vec![vec![0.; s.units.len()]; self.points.len()];
        for (k, point) in self.points.iter().enumerate() {
            let (t, f) = (point.interval, point.fraction);
            let volume: Vec<_> = schedule
                .storage
                .iter()
                .enumerate()
                .map(|(r, v)| v[t] * (1. - f) + v[t + 1] * f + point.correction[r].value(&columns))
                .collect();
            let state = crate::hydraulics::hydraulic_step(
                p,
                &controls,
                t,
                &volume,
                0.,
                &columns[..p.routing.primary_columns],
                None,
            )?;
            let power: Vec<_> = s
                .units
                .iter()
                .enumerate()
                .map(|(g, u)| {
                    u.trial_power(
                        state.unit_flow[g],
                        crate::hydraulics::machine_head(p, g, &state.head, &state.unit_flow),
                    )
                })
                .collect();
            for (g, u) in s.units.iter().enumerate() {
                let on = controls.commitment[g][t];
                let limits = &p.operating.units[g][t];
                let q = state.unit_flow[g];
                observed_heads[k][g] =
                    crate::hydraulics::machine_head(p, g, &state.head, &state.unit_flow);

                excess = excess
                    .max(violation(
                        q,
                        if on { limits.flow.min } else { 0. },
                        if on { limits.flow.max } else { 0. },
                        u.flow.max.max(1.),
                    ))
                    .max(violation(
                        power[g],
                        if on { limits.power.min } else { 0. },
                        if on { limits.power.max } else { 0. },
                        u.power.max.max(1.),
                    ));
                if on {
                    let h = crate::hydraulics::machine_head(p, g, &state.head, &state.unit_flow);
                    excess = excess
                        .max(violation(h, u.head.min, u.head.max, 100.))
                        .max(violation(
                            crate::machine::hydraulic_efficiency(u, q, h),
                            u.minimum_efficiency,
                            1.,
                            1.,
                        ));
                    if u.pump_characteristic.is_none()
                        && let crate::Efficiency::Table(table) = &u.efficiency
                    {
                        excess = excess.max(violation(
                            q,
                            table.envelope_data(h, false).value,
                            table.envelope_data(h, true).value,
                            u.flow.max.max(1.),
                        ));
                    }
                }
            }
            for (j, members) in p.topology.plant_units.iter().enumerate() {
                let limits = &p.operating.plants[j][t];
                let on = members.iter().any(|&g| controls.commitment[g][t]);
                let q = members.iter().map(|&g| state.unit_flow[g]).sum();
                let power_sum = members.iter().map(|&g| power[g]).sum();
                let qs = members
                    .iter()
                    .map(|&g| s.units[g].flow.max)
                    .sum::<f64>()
                    .max(1.);
                let ps = members
                    .iter()
                    .map(|&g| s.units[g].power.max)
                    .sum::<f64>()
                    .max(1.);
                excess = excess
                    .max(violation(
                        q,
                        if on { limits.flow.min } else { 0. },
                        limits.flow.max,
                        qs,
                    ))
                    .max(violation(
                        power_sum,
                        if on { limits.power.min } else { 0. },
                        limits.power.max,
                        ps,
                    ));
            }
            for (e, &q) in state.tunnel_flow.iter().enumerate() {
                let cap = if controls.tunnel_opening[e][t] == 0. {
                    0.
                } else {
                    p.operating.tunnels[e][t].capacity
                };
                excess = excess.max(violation(
                    q,
                    if p.topology.tunnel_outfalls[e].is_some() {
                        0.
                    } else {
                        -cap
                    },
                    cap,
                    s.tunnels[e].capacity.max(1.),
                ));
            }
            for (j, node) in s.junctions.iter().enumerate() {
                excess = excess.max(violation(
                    state.head[s.reservoirs.len() + j],
                    node.head.min,
                    node.head.max,
                    100.,
                ));
            }
        }
        for k in 0..self.points.len().saturating_sub(1) {
            let t = self.points[k].interval;
            if self.points[k + 1].interval != t {
                continue;
            }
            for (g, u) in s.units.iter().enumerate() {
                if !schedule.commitment[g][t] {
                    continue;
                }
                let limits = &p.operating.units[g][t];
                let envelope = crate::head_limits::HeadEnvelope {
                    unit: Arc::new(u.clone()),
                    power: limits.power,
                    flow: limits.flow,
                };
                let span = envelope.range(
                    schedule.unit_flow[g][t],
                    observed_heads[k][g],
                    observed_heads[k + 1][g],
                );
                excess = excess
                    .max(span.efficiency.max(0.))
                    .max(span.flow.max(0.) / u.flow.max.max(1.));
                for h in [span.at_upper, span.at_lower] {
                    let q = u
                        .pump_characteristic
                        .as_ref()
                        .map_or(schedule.unit_flow[g][t], |c| c.trial_discharge(h));
                    let power = u.trial_power(q, h);
                    excess = excess.max(violation(
                        power,
                        limits.power.min,
                        limits.power.max,
                        u.power.max.max(1.),
                    ));
                }
            }
        }
        Ok(excess)
    }
}
