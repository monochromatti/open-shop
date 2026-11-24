//! Simultaneous nonlinear water balances, pressure laws and operating equations.

use crate::Error;
use crate::nlp;
use oximo::{
    expr::{Expr, VarId},
    prelude::*,
};

pub(crate) struct Dispatch {
    pub problem: crate::constitutive::Nlp,
    pub variables: DispatchVariables,
}

pub(crate) struct DispatchVariables {
    pub storage: Vec<Vec<VarId>>,
    pub flow: Vec<Vec<VarId>>,
    pub tunnel_flow: Vec<Vec<VarId>>,
    pub river_withdrawal: Vec<Vec<VarId>>,
    pub river_opening: Vec<Vec<VarId>>,
    pub junction_head: Vec<Vec<VarId>>,
    pub unit_flow_scale: Vec<f64>,
    pub tunnel_flow_scale: Vec<f64>,
    pub river_flow_scale: Vec<f64>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum DispatchObjective {
    Revenue,
    FlowDistance,
}

pub(crate) fn build_dispatch(
    problem: &crate::SchedulingProblem<'_>,
    controls: &crate::commitment::DiscreteControls,
    start: Option<&crate::commitment::Proposal>,
    points: &crate::operating_points::OperatingPoints,
    objective: DispatchObjective,
) -> Result<Dispatch, Error> {
    let (case, topology, routing, operating) = (
        problem.case,
        &problem.topology,
        &problem.routing,
        &problem.operating,
    );
    let mesh = &problem.mesh;
    let nc = case.prices.len();
    let initial = start.map(|p| p.on_grid(&mesh.grid)).transpose()?;
    let start = initial.as_deref();
    let commitment = &controls.unit_on;
    let s = &case.system;
    let (nr, nu, nt) = (s.reservoirs.len(), s.units.len(), mesh.dt.len());
    let m = Model::new("hydropower dispatch");
    variable!(m,0.0 <= v[r in 0..nr,t in 0..nt+1] <= 1.0);
    let adjustable: Vec<_> = s
        .units
        .iter()
        .enumerate()
        .filter_map(|(g, u)| u.pump_characteristic.is_none().then_some(g))
        .collect();
    let pumps: Vec<_> = s
        .units
        .iter()
        .enumerate()
        .filter_map(|(g, u)| u.pump_characteristic.is_some().then_some(g))
        .collect();
    variable!(m,qc[g in 0..adjustable.len(),j in 0..nc] >= 0.0);
    variable!(m,qp[g in 0..pumps.len(),t in 0..nt] >= 0.0);
    // Adjustable discharges are controls; repeated hydraulic cells refer to
    // the same variable. Fixed-speed pump flows follow the pressure equation.
    let q: Vec<Vec<_>> = (0..nu)
        .map(|g| {
            if let Some(i) = adjustable.iter().position(|&i| i == g) {
                mesh.interval.iter().map(|&j| qc[[i, j]]).collect()
            } else {
                let i = pumps.iter().position(|&i| i == g).unwrap();
                (0..nt).map(|t| qp[[i, t]]).collect()
            }
        })
        .collect();
    variable!(m,qt[e in 0..s.tunnels.len(),t in 0..nt]);
    variable!(m,hj[j in 0..s.junctions.len(),t in 0..nt]);
    variable!(m,qr[e in 0..routing.outlets.len(),t in 0..nt]);
    variable!(m,0.<=gc[e in 0..routing.outlets.len(),j in 0..nc]<=1.);
    let ga: Vec<Vec<_>> = (0..routing.outlets.len())
        .map(|e| mesh.interval.iter().map(|&j| gc[[e, j]]).collect())
        .collect();
    let mut curves = Vec::new();
    let curve_index: Vec<_> = s
        .reservoirs
        .iter()
        .map(|r| {
            if r.level.trial_level(r.storage.min) != r.level.trial_level(r.storage.max) {
                let index = curves.len();
                curves.push(std::sync::Arc::new(r.level.clone()));
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
    let level_reservoirs: Vec<_> = s
        .reservoirs
        .iter()
        .enumerate()
        .filter_map(|(r, reservoir)| {
            if matches!(reservoir.level, crate::StorageCurve::Table(_))
                && (problem.soft_limits.needs_level(r)
                    || operating.reservoirs[r]
                        .level_change
                        .iter()
                        .any(|c| c.up.is_some() || c.down.is_some()))
            {
                Some(r)
            } else {
                None
            }
        })
        .collect();
    variable!(m,hv[k in 0..level_reservoirs.len(),t in 0..nc+1]);
    let mut laws = Vec::new();
    let flow_scale = s
        .units
        .iter()
        .map(|u| u.flow.max)
        .chain(s.tunnels.iter().map(|e| e.capacity))
        .chain(s.rivers.iter().map(|r| r.capacity))
        .fold(1.0, f64::max);
    let unit_scale: Vec<_> = s.units.iter().map(|u| u.flow.max.max(1.)).collect();
    let tunnel_scale: Vec<_> = s.tunnels.iter().map(|e| e.capacity.max(1.)).collect();
    let river_scale: Vec<_> = s.rivers.iter().map(|e| e.capacity.max(1.)).collect();
    let zero = || -> Expr<'_> { Expr::constant(v[[0, 0]].arena(), 0.0).into() };
    let mut river_columns = Vec::new();
    for (u, row) in q.iter().enumerate() {
        for &value in row {
            river_columns.push(Expr::from(value * unit_scale[u]));
        }
    }
    for e in 0..s.tunnels.len() {
        for t in 0..nt {
            river_columns.push((qt[[e, t]] * tunnel_scale[e]).into());
        }
    }
    for (k, &e) in routing.outlets.iter().enumerate() {
        for t in 0..nt {
            let lower = if s.rivers[e].law.is_delta() {
                -operating.river_capacity[e][mesh.interval[t]] / river_scale[e]
            } else {
                0.
            };
            if s.rivers[e].law.is_delta() {
                m.set_bounds(
                    qr[[k, t]],
                    lower,
                    operating.river_capacity[e][mesh.interval[t]] / river_scale[e],
                )?;
            }
            // A nonnegative gate and rating already give nonnegative outlet
            // withdrawal. Total-release bounds below impose capacity. Repeating
            // these as withdrawal bounds makes closed outlets degenerate.
            m.set_initial(
                qr[[k, t]],
                start
                    .map_or(0., |p| p.river_withdrawal[e][t] / river_scale[e])
                    .clamp(
                        lower,
                        operating.river_capacity[e][mesh.interval[t]] / river_scale[e],
                    ),
            )?;
            river_columns.push(Expr::from(qr[[k, t]] * river_scale[e]));
        }
    }
    variable!(m,rf[k in 0..routing.transfers.len()]);
    let mut reference_columns: Vec<_> = start.map_or(vec![0.; routing.primary_columns], |p| {
        p.flow
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
        let value = e.fraction.data(q0).value;
        reference_columns.push(value);
        m.set_initial(rf[[k]], value / flow_scale)?;
        let input = nlp::linear_expression(&e.flow, &river_columns, zero());
        laws.push(crate::constitutive::Equation {
            linear: crate::river::LinearForm::variable(
                rf[[k]].var_id().unwrap().index(),
                -flow_scale,
            ),
            terms: vec![crate::constitutive::Term {
                function: crate::constitutive::PhysicalFunction::Transfer(e.fraction.clone()),
                inputs: [
                    nlp::linear_form(&m, input)?,
                    Default::default(),
                    Default::default(),
                ],
                weight: 1.,
            }],
            bounds: crate::Bounds::new(0., 0.),
        });
        river_columns.push((rf[[k]] * flow_scale).into());
    }
    problem.observations.constrain(&m, &river_columns, zero());
    let river_on: Vec<Vec<_>> = controls
        .river_on
        .iter()
        .map(|row| {
            row.iter()
                .map(|&on| zero() + if on { 1. } else { 0. })
                .collect()
        })
        .collect();
    let river_starts: Vec<Vec<_>> = controls
        .river_on
        .iter()
        .enumerate()
        .map(|(k, row)| {
            row.iter()
                .enumerate()
                .map(|(t, &on)| {
                    let prior = if t == 0 {
                        problem.observations.domains[k].initial_on
                    } else {
                        row[t - 1]
                    };
                    zero() + if on && !prior { 1. } else { 0. }
                })
                .collect()
        })
        .collect();
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
                .min(nc - 1);
            let release = crate::nlp::linear_expression(rate, &river_columns, zero());
            constraint!(
                m,
                release / flow_scale <= operating.river_capacity[i][period] / flow_scale
            );
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
        for t in 0..=nt {
            let bounds = match mesh.vertices.binary_search(&t) {
                Ok(j) => operating.reservoirs[r].storage[j],
                Err(_) => operating.reservoirs[r].interval_storage[mesh.interval[t]],
            };
            m.set_bounds(
                v[[r, t]],
                (bounds.min - reservoir.storage.min) / range,
                (bounds.max - reservoir.storage.min) / range,
            )?;
        }
        let initial = (reservoir.initial_storage - reservoir.storage.min) / range;
        m.fix(v[[r, 0]], initial)?;
        for t in 0..=nt {
            let guess = start
                .map_or(initial, |p| {
                    (p.storage[r][t] - reservoir.storage.min) / range
                })
                .clamp(0., 1.);
            m.set_initial(v[[r, t]], guess)?;
        }
        for t in 0..nc {
            let change = (v[[r, mesh.vertices[t + 1]]] - v[[r, mesh.vertices[t]]]) * range;
            if let Some(up) = operating.reservoirs[r].volume_change[t].up {
                constraint!(m, change <= up);
            }
            if let Some(down) = operating.reservoirs[r].volume_change[t].down {
                constraint!(m, -change <= down);
            }
        }
        for t in 0..nt {
            let mut net = zero() + operating.reservoirs[r].inflow[mesh.interval[t]];
            for &(u, sign) in &topology.incidence[r] {
                net = net + q[u][t] * (sign * unit_scale[u]);
            }
            for &(e, sign) in &topology.tunnel_incidence[r] {
                net = net + qt[[e, t]] * (sign * tunnel_scale[e]);
            }
            net = net
                + crate::nlp::linear_expression(
                    &routing.reservoir_net[r][t],
                    &river_columns,
                    zero(),
                );
            let rate = (v[[r, t + 1]] - v[[r, t]]) * (range / (3600. * mesh.dt[t]));
            constraint!(
                m,
                name = format!("water_balance[{r},{t}]"),
                rate - net == 0.0
            );
        }
    }
    let mut vertex_levels = vec![vec![zero(); nc + 1]; nr];
    for (r, reservoir) in s.reservoirs.iter().enumerate() {
        if !operating.reservoirs[r]
            .level_change
            .iter()
            .any(|c| c.up.is_some() || c.down.is_some())
            && !problem.soft_limits.needs_level(r)
        {
            continue;
        }
        let range = reservoir.storage.max - reservoir.storage.min;
        let mut levels: Vec<Expr<'_>> = Vec::new();
        for t in 0..=nc {
            let volume: Expr<'_> =
                (v[[r, mesh.vertices[t]]] * range + reservoir.storage.min).into();
            let level = match &reservoir.level {
                crate::StorageCurve::Polynomial {
                    reference_volume,
                    elevation,
                    slope,
                    curvature,
                } => crate::physics::polynomial_mean_level(
                    *reference_volume,
                    *elevation,
                    *slope,
                    *curvature,
                    volume,
                    volume,
                ),
                crate::StorageCurve::Table(c) => {
                    let k = level_reservoirs.iter().position(|&i| i == r).unwrap();
                    m.set_initial(
                        hv[[k, t]],
                        reservoir.level.trial_level(reservoir.initial_storage) / 100.,
                    )?;
                    laws.push(crate::constitutive::Equation {
                        linear: crate::river::LinearForm::variable(
                            hv[[k, t]].var_id().unwrap().index(),
                            -1.,
                        ),
                        terms: vec![crate::constitutive::Term {
                            function: crate::constitutive::PhysicalFunction::Table(
                                std::sync::Arc::new(c.clone()),
                            ),
                            inputs: [
                                crate::river::LinearForm {
                                    constant: reservoir.storage.min,
                                    coefficients: std::collections::BTreeMap::from([(
                                        v[[r, mesh.vertices[t]]].var_id().unwrap().index(),
                                        range,
                                    )]),
                                },
                                Default::default(),
                                Default::default(),
                            ],
                            weight: 1. / 100.,
                        }],
                        bounds: crate::Bounds::new(0., 0.),
                    });
                    (hv[[k, t]] * 100.).into()
                }
            };
            levels.push(level);
        }
        vertex_levels[r].clone_from(&levels);
        for t in 0..nc {
            let change = (levels[t + 1] - levels[t]) / 100.;
            if let Some(up) = operating.reservoirs[r].level_change[t].up {
                constraint!(m, change <= up / 100.);
            }
            if let Some(down) = operating.reservoirs[r].level_change[t].down {
                constraint!(m, -change <= down / 100.);
            }
        }
    }
    for (j, junction) in s.junctions.iter().enumerate() {
        let node = nr + j;
        for t in 0..nt {
            m.set_bounds(
                hj[[j, t]],
                junction.head.min / 100.,
                junction.head.max / 100.,
            )?;
            m.set_initial(
                hj[[j, t]],
                start
                    .map_or((junction.head.min + junction.head.max) / 2., |p| {
                        p.head[node][t]
                    })
                    .clamp(junction.head.min, junction.head.max)
                    / 100.,
            )?;
            let mut net = zero();
            for &(u, sign) in &topology.incidence[node] {
                net = net + q[u][t] * (sign * unit_scale[u]);
            }
            for &(e, sign) in &topology.tunnel_incidence[node] {
                net = net + qt[[e, t]] * (sign * tunnel_scale[e]);
            }
            constraint!(
                m,
                name = format!("junction_balance[{j},{t}]"),
                net / flow_scale == 0.
            );
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
                    constraint!(
                        m,
                        volume
                            >= (operating.reservoirs[r].interval_storage[mesh.interval[t]].min
                                - reservoir.storage.min)
                                / range
                    );
                    constraint!(
                        m,
                        volume
                            <= (operating.reservoirs[r].interval_storage[mesh.interval[t]].max
                                - reservoir.storage.min)
                                / range
                    );
                }
                vertices.push((*fraction, volume * range + reservoir.storage.min));
            }
            let level = {
                let lo = c.trial_level(reservoir.storage.min);
                let hi = c.trial_level(reservoir.storage.max);
                if let Some(k) = curve_index[r] {
                    m.set_bounds(hm[[k, t]], lo / 100., hi / 100.)?;
                    let (before, after) = start.map_or(
                        (reservoir.initial_storage, reservoir.initial_storage),
                        |p| (p.storage[r][t], p.storage[r][t + 1]),
                    );
                    m.set_initial(
                        hm[[k, t]],
                        start
                            .map_or(c.trial_mean_level(before, after), |p| p.head[r][t])
                            .clamp(lo, hi)
                            / 100.,
                    )?;
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
    let operating_head: Vec<Vec<_>> = mesh
        .cells
        .iter()
        .enumerate()
        .map(|(j, cells)| {
            (0..(nr + s.junctions.len() + s.boundaries.len()))
                .map(|n| {
                    cells.clone().fold(zero(), |sum, k| {
                        sum + head[k][n] * (mesh.dt[k] / mesh.duration[j])
                    })
                })
                .collect()
        })
        .collect();
    // Reservoir outlet laws share the same interval mean heads as machines
    // and tunnels. Square-root lifts keep wet outlet derivatives finite.
    variable!(m,root[k in 0..routing.outlets.len(),t in 0..nt]>=0.);
    for (k, &e) in routing.outlets.iter().enumerate() {
        let river = &s.rivers[e];
        let source = s
            .reservoirs
            .iter()
            .position(|r| r.name == river.source)
            .unwrap();
        for t in 0..nt {
            let gate = operating.river_gate[e][mesh.interval[t]];
            m.set_bounds(ga[k][t], gate.min, gate.max)?;
            m.set_initial(
                ga[k][t],
                start
                    .map_or((gate.min + gate.max) / 2., |p| p.river_opening[e][t])
                    .clamp(gate.min, gate.max),
            )?;
            let h = head[t][source];
            let withdrawal = qr[[k, t]] * river_scale[e];
            match &river.law {
                crate::RiverLaw::Controlled => {
                    constraint!(
                        m,
                        withdrawal
                            - ga[k][t]
                                * (operating.river_capacity[e][mesh.interval[t]]
                                    - operating.river_inflow[e][mesh.interval[t]])
                            == 0.
                    );
                    m.fix(root[[k, t]], 0.)?;
                }
                crate::RiverLaw::Delta { curve, reference } => {
                    let target = s
                        .reservoirs
                        .iter()
                        .position(|r| r.name == river.target)
                        .unwrap();
                    let href = match reference {
                        crate::ReferenceLevel::Upstream => h,
                        crate::ReferenceLevel::Downstream => head[t][target],
                    };
                    laws.push(crate::constitutive::Equation {
                        linear: crate::river::LinearForm::variable(
                            qr[[k, t]].var_id().unwrap().index(),
                            -river_scale[e],
                        ),
                        terms: vec![crate::constitutive::Term {
                            function: crate::constitutive::PhysicalFunction::LevelFlow(
                                std::sync::Arc::new(curve.clone()),
                            ),
                            inputs: [
                                nlp::linear_form(&m, h - head[t][target])?,
                                nlp::linear_form(&m, href)?,
                                Default::default(),
                            ],
                            weight: 1.,
                        }],
                        bounds: crate::Bounds::new(0., 0.),
                    });
                    m.fix(ga[k][t], 0.)?;
                    m.fix(root[[k, t]], 0.)?;
                }
                crate::RiverLaw::Level { curve } => {
                    laws.push(crate::constitutive::Equation {
                        linear: crate::river::LinearForm::variable(
                            qr[[k, t]].var_id().unwrap().index(),
                            -river_scale[e],
                        ),
                        terms: vec![crate::constitutive::Term {
                            function: crate::constitutive::PhysicalFunction::OutletTable(
                                std::sync::Arc::new(curve.clone()),
                            ),
                            inputs: [
                                nlp::linear_form(&m, h)?,
                                nlp::linear_form(&m, ga[k][t].into())?,
                                Default::default(),
                            ],
                            weight: 1.,
                        }],
                        bounds: crate::Bounds::new(0., 0.),
                    });
                    m.fix(root[[k, t]], 0.)?;
                }
                crate::RiverLaw::Orifice {
                    coefficient,
                    crest,
                    allow_dry,
                }
                | crate::RiverLaw::Weir {
                    coefficient,
                    crest,
                    allow_dry,
                } => {
                    let reservoir = &s.reservoirs[source];
                    let lo = (reservoir.level.trial_level(reservoir.storage.min) - crest)
                        .max(0.)
                        .sqrt();
                    let hi = (reservoir.level.trial_level(reservoir.storage.max) - crest)
                        .max(0.)
                        .sqrt();
                    m.set_bounds(root[[k, t]], lo, hi)?;
                    m.set_initial(
                        root[[k, t]],
                        start
                            .map_or((lo + hi) / 2., |p| {
                                (p.head[source][t] - crest).max(0.).sqrt()
                            })
                            .clamp(lo, hi),
                    )?;
                    let depth = h - *crest;
                    let wet_depth = if *allow_dry {
                        (depth + depth.abs()) * 0.5
                    } else {
                        depth
                    };
                    if lo < hi {
                        constraint!(m, (root[[k, t]].square() - wet_depth) / 100. == 0.);
                    }
                    let discharge: Expr<'_> = if matches!(river.law, crate::RiverLaw::Weir { .. }) {
                        ((root[[k, t]] * root[[k, t]] * root[[k, t]]) * *coefficient).into()
                    } else {
                        (root[[k, t]] * ga[k][t] * *coefficient).into()
                    };
                    constraint!(m, withdrawal - discharge == 0.);
                }
            }
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
        m.set_initial(
            band_floor[[k]],
            r.level
                .trial_level(r.initial_storage)
                .clamp(band.floor.min, band.floor.max)
                / 100.,
        )?;
        let curve = std::sync::Arc::new(r.level.clone());
        for point in &band.points {
            let (t, fraction) = (point.interval, point.fraction);
            let volume = ((v[[i, t]] * (1. - fraction) + v[[i, t + 1]] * fraction) * range
                + r.storage.min)
                + nlp::linear_expression(&point.correction, &river_columns, zero());
            let input = nlp::linear_form(&m, volume)?;
            laws.push(crate::constitutive::Equation {
                linear: crate::river::LinearForm::variable(
                    band_floor[[k]].var_id().unwrap().index(),
                    -1.,
                ),
                terms: vec![crate::constitutive::Term {
                    function: crate::constitutive::PhysicalFunction::StorageMean(curve.clone()),
                    inputs: [input.clone(), input, Default::default()],
                    weight: 1. / 100.,
                }],
                bounds: crate::Bounds::new(0., band.max_change / 100.),
            });
        }
    }
    let initial_head: Vec<_> = s
        .reservoirs
        .iter()
        .map(|r| r.level.trial_level(r.initial_storage))
        .chain(s.junctions.iter().map(|j| (j.head.min + j.head.max) / 2.))
        .chain(s.boundaries.iter().map(|b| b.head))
        .collect();
    for e in 0..s.tunnels.len() {
        let (source, target) = topology.tunnel_nodes[e];
        let dh = initial_head[source] - initial_head[target];
        for t in 0..nt {
            let interval = &operating.tunnels[e][mesh.interval[t]];
            let law = &interval.positions[controls.gate_position[e][mesh.interval[t]]];
            let guess = dh.signum() * (dh.abs() / law.resistance).sqrt();
            let capacity = if law.opening == 0. {
                0.
            } else {
                interval.capacity / tunnel_scale[e]
            };
            m.set_bounds(
                qt[[e, t]],
                if topology.tunnel_outfalls[e].is_some() {
                    0.
                } else {
                    -capacity
                },
                capacity,
            )?;
            m.set_initial(
                qt[[e, t]],
                start
                    .map_or(guess / tunnel_scale[e], |p| {
                        p.tunnel_flow[e][t] / tunnel_scale[e]
                    })
                    .clamp(
                        if topology.tunnel_outfalls[e].is_some() {
                            0.
                        } else {
                            -capacity
                        },
                        capacity,
                    ),
            )?;
            if law.opening > 0. {
                let flow = qt[[e, t]] * tunnel_scale[e];
                let pressure =
                    (head[t][source] - head[t][target]) - flow * flow.abs() * law.resistance;
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
                    .insert(q[g][t].var_id().unwrap().index(), unit_scale[g]);
            }
            m.set_initial(
                hw[[k, t]],
                start.map_or(curve.values()[0], |start| {
                    curve.extended_value(
                        topology.plant_units[p]
                            .iter()
                            .map(|&g| start.flow[g][t])
                            .sum(),
                    )
                }) / 100.,
            )?;
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
            let limits = &operating.units[i][mesh.interval[t]];
            let on = commitment[i][mesh.interval[t]];
            let (lo, hi) = if on {
                (
                    limits.flow.min / unit_scale[i],
                    limits.flow.max / unit_scale[i],
                )
            } else {
                (0., 0.)
            };
            m.set_bounds(q[i][t], lo, hi)?;
            m.set_initial(
                q[i][t],
                start
                    .map_or((lo + hi) / 2., |p| p.flow[i][t] / unit_scale[i])
                    .clamp(lo, hi),
            )?;
            let flow: Expr<'_> = (q[i][t] * unit_scale[i]).into();
            let plant = topology.unit_plants[i];
            let p = &s.plants[plant];
            let outlet = p.outlet_head_floor.map_or(head[t][target], |floor| {
                (head[t][target] + floor + (head[t][target] - floor).abs()) * 0.5
            });
            let net_head =
                (head[t][source] - outlet - tailwater[plant][t] - flow.square() * u.loss_factor)
                    * u.mode.power_sign();
            if t == mesh.cells[mesh.interval[t]].start
                && let Some(reference) = &topology.gross_heads[i]
            {
                let gross = (operating_head[mesh.interval[t]][reference.source]
                    - operating_head[mesh.interval[t]][reference.target])
                    * u.mode.power_sign();
                if on {
                    if let Some(lo) = limits.gross_head_min {
                        constraint!(m, gross / 100. >= lo / 100.);
                    }
                    if let Some(hi) = limits.gross_head_max {
                        constraint!(m, gross / 100. <= hi / 100.);
                    }
                }
                if u.operation.transitions(&commitment[i], mesh.interval[t]).0
                    && let Some(lo) = limits.startup_gross_head_min
                {
                    let mut equation = crate::constitutive::Equation {
                        linear: Default::default(),
                        terms: Vec::new(),
                        bounds: crate::Bounds::new(lo / 100., f64::INFINITY),
                    };
                    for (node, sign) in [
                        (reference.source, u.mode.power_sign()),
                        (reference.target, -u.mode.power_sign()),
                    ] {
                        if node < nr {
                            let r = &s.reservoirs[node];
                            let range = r.storage.max - r.storage.min;
                            let volume = nlp::linear_form(
                                &m,
                                (v[[node, mesh.vertices[mesh.interval[t]]]] * range
                                    + r.storage.min)
                                    .into(),
                            )?;
                            equation.terms.push(crate::constitutive::Term {
                                function: crate::constitutive::PhysicalFunction::StorageMean(
                                    std::sync::Arc::new(r.level.clone()),
                                ),
                                inputs: [volume.clone(), volume, Default::default()],
                                weight: sign / 100.,
                            });
                        } else {
                            equation.linear.constant +=
                                sign * s.boundaries[node - nr - s.junctions.len()].head / 100.;
                        }
                    }
                    laws.push(equation);
                }
            }
            let scale = u.power.max.max(1.);
            if on {
                m.set_bounds(ht[[i, t]], u.head.min / 100., u.head.max / 100.)?;
                m.set_bounds(
                    pt[[i, t]],
                    limits.power.min / scale,
                    limits.power.max / scale,
                )?;
                let upstream = start.map_or(initial_head[source], |p| p.head[source][t]);
                let downstream = start
                    .map_or(initial_head[target], |p| p.head[target][t])
                    .max(p.outlet_head_floor.unwrap_or(f64::NEG_INFINITY));
                let total_flow = topology.plant_units[plant]
                    .iter()
                    .map(|&g| {
                        start.map_or((s.units[g].flow.min + s.units[g].flow.max) / 2., |p| {
                            p.flow[g][t]
                        })
                    })
                    .sum();
                let tail = p
                    .tailwater_curve
                    .as_ref()
                    .map_or(0., |c| c.extended_value(total_flow));
                let q0 = start
                    .map_or((u.flow.min + u.flow.max) / 2., |p| p.flow[i][t])
                    .clamp(u.flow.min, u.flow.max);
                let h0 = ((upstream - downstream - tail - u.loss_factor * q0 * q0)
                    * u.mode.power_sign())
                .clamp(u.head.min, u.head.max);
                m.set_initial(ht[[i, t]], h0 / 100.)?;
                let guess = u.trial_power(q0, h0);
                let p0 = if guess.is_finite() {
                    guess.clamp(limits.power.min, limits.power.max)
                } else {
                    (u.power.min + u.power.max) / 2.
                };
                m.set_initial(pt[[i, t]], p0 / scale)?;
                constraint!(m, ht[[i, t]] - net_head / 100. == 0.);
                if u.head.min < 0. {
                    constraint!(m, flow * ht[[i, t]] / unit_scale[i] >= 0.);
                }
                laws.extend(crate::constitutive::Equation::production(
                    machines[i].clone(),
                    crate::river::LinearForm::variable(
                        q[i][t].var_id().unwrap().index(),
                        unit_scale[i],
                    ),
                    crate::river::LinearForm::variable(ht[[i, t]].var_id().unwrap().index(), 100.),
                    crate::river::LinearForm::variable(pt[[i, t]].var_id().unwrap().index(), scale),
                    unit_scale[i],
                    scale,
                ));
            } else {
                m.fix(ht[[i, t]], 0.)?;
                m.fix(pt[[i, t]], 0.)?;
            }
            let power: Expr<'_> = (pt[[i, t]] * scale).into();
            powers.push(power);
            flows.push(flow);
            revenue = revenue
                + power * (u.mode.power_sign() * case.prices[mesh.interval[t]] * mesh.dt[t]);
            let (start, stop) = u.operation.transitions(&commitment[i], mesh.interval[t]);
            if t == mesh.cells[mesh.interval[t]].start {
                revenue = revenue
                    - if start { limits.startup_cost } else { 0. }
                    - if stop { limits.shutdown_cost } else { 0. };
            }
        }
        for j in 0..nc {
            let (start, stop) = u.operation.transitions(&commitment[i], j);
            states.push(zero() + if commitment[i][j] { 1. } else { 0. });
            starts.push(zero() + if start { 1. } else { 0. });
            stops.push(zero() + if stop { 1. } else { 0. });
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
            (0..nc)
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
    let gate_on: Vec<Vec<_>> = operating
        .tunnels
        .iter()
        .enumerate()
        .map(|(e, row)| {
            row.iter()
                .enumerate()
                .map(|(t, interval)| {
                    zero()
                        + if interval.positions[controls.gate_position[e][t]].opening > 0. {
                            1.
                        } else {
                            0.
                        }
                })
                .collect()
        })
        .collect();
    problem
        .interlocks
        .constrain(&m, &all_on, &gate_on, zero(), nc);
    let volumes: Vec<Vec<Expr<'_>>> = s
        .reservoirs
        .iter()
        .enumerate()
        .map(|(r, reservoir)| {
            (0..=nt)
                .map(|t| {
                    (v[[r, t]] * (reservoir.storage.max - reservoir.storage.min)
                        + reservoir.storage.min)
                        .into()
                })
                .collect()
        })
        .collect();
    points.constrain(
        problem,
        &m,
        &volumes,
        &all_flow,
        &river_columns,
        controls,
        start,
        &mut laws,
    )?;
    revenue = revenue
        - problem.soft_limits.constrain(
            &m,
            &mesh.states(&volumes)?,
            &vertex_levels,
            &mesh.means(&all_power, zero())?,
            &mesh.means(&all_flow, zero())?,
            zero(),
        )?;
    crate::operating_equations::constrain(
        &m,
        case,
        topology,
        operating,
        &crate::operating_equations::OperatingExpressions {
            power: mesh.means(&all_power, zero())?,
            flow: mesh.means(&all_flow, zero())?,
            on: all_on,
            starts: all_starts,
            stops: all_stops,
            plant_on,
        },
    );
    revenue = revenue + nlp::linear_expression(&problem.water_value, &river_columns, zero());
    revenue = revenue - problem.observations.cost(&m, &river_columns, zero());
    revenue = revenue - river_transition_cost;
    let openings: Vec<Vec<f64>> = operating
        .tunnels
        .iter()
        .enumerate()
        .map(|(e, row)| {
            row.iter()
                .enumerate()
                .map(|(t, p)| p.positions[controls.gate_position[e][t]].opening)
                .collect()
        })
        .collect();
    revenue = revenue - crate::trajectory::gate_adjustment_cost(&s.tunnels, &openings);
    match objective {
        DispatchObjective::Revenue => {
            objective!(m, Max, revenue / 1000.);
        }
        DispatchObjective::FlowDistance => {
            // The closest water plan is measured in one common discharge unit.
            // A spillway's large capacity must not make spilling inexpensive.
            let mut distance = zero();
            for t in 0..nt {
                let duration = mesh.dt[t];
                for (g, u) in s.units.iter().enumerate() {
                    let target = start.map_or_else(
                        || {
                            if controls.unit_on[g][mesh.interval[t]] {
                                (u.flow.min + u.flow.max) / 2.
                            } else {
                                0.
                            }
                        },
                        |p| p.flow[g][t],
                    );
                    let change = (river_columns[g * nt + t] - target) / 50.;
                    distance = distance + change.square() * duration;
                }
                for series in &routing.releases {
                    let rate = series
                        .integral(mesh.grid[t], mesh.grid[t + 1])
                        .scaled(1. / duration);
                    let target = rate.value(&reference_columns);
                    let change =
                        (nlp::linear_expression(&rate, &river_columns, zero()) - target) / 50.;
                    distance = distance + change.square() * duration;
                }
            }
            objective!(m, Min, distance);
        }
    }
    let storage = (0..nr)
        .map(|r| (0..=nt).map(|t| v[[r, t]].var_id().unwrap()).collect())
        .collect();
    let flow = (0..nu)
        .map(|u| (0..nt).map(|t| q[u][t].var_id().unwrap()).collect())
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
        variables: DispatchVariables {
            storage,
            flow,
            tunnel_flow,
            junction_head,
            river_withdrawal,
            river_opening: (0..routing.outlets.len())
                .map(|k| (0..nt).map(|t| ga[k][t].var_id().unwrap()).collect())
                .collect(),
            unit_flow_scale: unit_scale,
            tunnel_flow_scale: tunnel_scale,
            river_flow_scale: river_scale,
        },
    })
}
