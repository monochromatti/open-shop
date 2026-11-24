//! Flow observations and exact time-window integrals of interval averages.
//!
//! An observation constrains existing flows; it adds no water or hydraulic node.
use crate::{
    Bounds, Error, ScheduleCase,
    river::{LinearForm, Routing},
};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VolumeWindow {
    pub start_time: f64,
    pub end_time: f64,
    #[serde(default)]
    pub min_volume: Option<f64>,
    #[serde(default)]
    pub max_volume: Option<f64>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RollingRamp {
    pub hours: f64,
    #[serde(default)]
    pub up: Option<f64>,
    #[serde(default)]
    pub down: Option<f64>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PeriodRamp {
    pub grid: Vec<f64>,
    #[serde(default)]
    pub up: Option<f64>,
    #[serde(default)]
    pub down: Option<f64>,
}
/// Named sum of existing rates, m³/s; volume bounds and counted water use m³.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FlowRequirement {
    pub name: String,
    #[serde(default)]
    pub units: Vec<String>,
    #[serde(default)]
    pub rivers: Vec<String>,
    #[serde(default)]
    pub river_arrivals: Vec<String>,
    #[serde(default)]
    pub tunnels: Vec<String>,
    #[serde(default)]
    pub inflow: f64,
    #[serde(default)]
    pub min_flow: f64,
    #[serde(default)]
    pub max_flow: Option<f64>,
    #[serde(default)]
    pub min_volume: Option<f64>,
    #[serde(default)]
    pub max_volume: Option<f64>,
    #[serde(default)]
    pub initial_volume: f64,
    #[serde(default)]
    pub volume_windows: Vec<VolumeWindow>,
    #[serde(default)]
    pub rolling_ramps: Vec<RollingRamp>,
    #[serde(default)]
    pub period_ramps: Vec<PeriodRamp>,
    #[serde(default)]
    pub history_grid: Vec<f64>,
    #[serde(default)]
    pub history_flow: Vec<f64>,
}
impl FlowRequirement {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            units: vec![],
            rivers: vec![],
            river_arrivals: vec![],
            tunnels: vec![],
            inflow: 0.,
            min_flow: 0.,
            max_flow: None,
            min_volume: None,
            max_volume: None,
            initial_volume: 0.,
            volume_windows: vec![],
            rolling_ramps: vec![],
            period_ramps: vec![],
            history_grid: vec![],
            history_flow: vec![],
        }
    }
}
pub(crate) struct FlowInterval {
    pub inflow: f64,
    pub bounds: Bounds,
}
#[derive(Clone)]
pub(crate) struct FlowRow {
    pub expression: LinearForm,
    pub bounds: Bounds,
    pub scale: f64,
}
pub(crate) struct FlowDomain {
    pub river: usize,
    pub minimum: f64,
    pub initial_on: bool,
    rates: Vec<LinearForm>,
    capacity: Vec<f64>,
    startup_cost: Vec<f64>,
    shutdown_cost: Vec<f64>,
}
#[derive(Default)]
pub(crate) struct FlowEquations {
    pub rates: Vec<Vec<LinearForm>>,
    pub rows: Vec<FlowRow>,
    changes: Vec<(LinearForm, f64, f64)>,
    pub domains: Vec<FlowDomain>,
}
fn ordered(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && b > a && (b - a).is_finite()
}
fn bounds_valid(lo: Option<f64>, hi: Option<f64>) -> bool {
    lo.into_iter().chain(hi).all(|v| v.is_finite()) && !matches!((lo,hi),(Some(a),Some(b)) if a>b)
}
fn ramp_valid(up: Option<f64>, down: Option<f64>) -> bool {
    (up.is_some() || down.is_some()) && up.into_iter().chain(down).all(|v| v.is_finite() && v >= 0.)
}

impl FlowEquations {
    pub fn new(
        case: &ScheduleCase,
        routing: &Routing,
        operating: &crate::operating_data::OperatingData,
        mesh: &crate::time_discretization::TimeDiscretization,
    ) -> Result<Self, Error> {
        let mut data = Self {
            rows: routing.arrival_limits.clone(),
            ..Default::default()
        };
        let s = &case.system;
        let nt = case.prices.len();
        let (start, end) = (case.grid[0], case.grid[nt]);
        let nf = mesh.dt.len();
        let tunnel_offset = s.units.len() * nf;
        let mut names = HashSet::new();
        // Environmental release minima constrain the control-interval average,
        // including natural inflow and routed arrivals at a river junction.
        for (e, series) in routing.releases.iter().enumerate() {
            let river = &s.rivers[e];
            let rates: Vec<_> = case
                .grid
                .windows(2)
                .map(|w| series.integral(w[0], w[1]).scaled(1. / (w[1] - w[0])))
                .collect();
            if let Some(minimum) = river.operation.minimum_on_flow {
                data.domains.push(FlowDomain {
                    river: e,
                    minimum,
                    initial_on: river.flow_on(river.previous_release().unwrap().0),
                    rates: rates.clone(),
                    capacity: operating.river_capacity[e].clone(),
                    startup_cost: operating.river_startup_cost[e].clone(),
                    shutdown_cost: operating.river_shutdown_cost[e].clone(),
                });
            }
            for t in 0..nt {
                let lower = operating.river_min_release[e][t];
                if lower > 0. {
                    data.rows.push(FlowRow {
                        expression: rates[t].clone(),
                        bounds: Bounds::new(lower, f64::INFINITY),
                        scale: 100.,
                    });
                }
                let prior = if t == 0 {
                    river.previous_release().map(|(q, _)| LinearForm {
                        constant: q,
                        ..Default::default()
                    })
                } else {
                    Some(rates[t - 1].clone())
                };
                if let Some(prior) = prior {
                    let mut difference = rates[t].clone();
                    difference.add_scaled(&prior, -1.);
                    let ramp = operating.river_flow_change[e][t];
                    let scale = river.capacity.max(1.);
                    if ramp.up.is_some() || ramp.down.is_some() {
                        data.rows.push(FlowRow {
                            expression: difference.clone(),
                            bounds: Bounds::new(
                                ramp.down.map_or(f64::NEG_INFINITY, |q| -q),
                                ramp.up.unwrap_or(f64::INFINITY),
                            ),
                            scale,
                        });
                    }
                    let cost = operating.river_change_cost[e][t];
                    if cost > 0. {
                        data.changes.push((difference, cost, scale));
                    }
                }
            }
        }
        for (i, rule) in case.flow_requirements.iter().enumerate() {
            let invalid = |message: &str| Error::Input(format!("{}: {message}", rule.name));
            if rule.name.is_empty()
                || !names.insert(&rule.name)
                || !rule.inflow.is_finite()
                || rule.inflow < 0.
                || !rule.min_flow.is_finite()
                || !rule.initial_volume.is_finite()
                || !bounds_valid(Some(rule.min_flow), rule.max_flow)
                || !bounds_valid(rule.min_volume, rule.max_volume)
            {
                return Err(invalid("invalid flow observation bounds or name"));
            }
            let indices =
                |requested: &[String], available: Vec<&str>| -> Result<Vec<usize>, Error> {
                    let mut seen = HashSet::new();
                    requested
                        .iter()
                        .map(|name| {
                            if seen.insert(name) {
                                available
                                    .iter()
                                    .position(|n| *n == name)
                                    .ok_or_else(|| invalid("unknown flow observation contribution"))
                            } else {
                                Err(invalid("duplicate flow observation contribution"))
                            }
                        })
                        .collect()
                };
            let units = indices(
                &rule.units,
                s.units.iter().map(|u| u.name.as_str()).collect(),
            )?;
            let rivers = indices(
                &rule.rivers,
                s.rivers.iter().map(|r| r.name.as_str()).collect(),
            )?;
            let arrivals = indices(
                &rule.river_arrivals,
                s.rivers.iter().map(|r| r.name.as_str()).collect(),
            )?;
            let tunnels = indices(
                &rule.tunnels,
                s.tunnels.iter().map(|t| t.name.as_str()).collect(),
            )?;
            let empty_history = rule.history_grid.is_empty() && rule.history_flow.is_empty();
            if !empty_history
                && (rule.history_grid.len() != rule.history_flow.len() + 1
                    || rule.history_grid.len() < 2
                    || rule.history_grid.last() != Some(&start)
                    || !rule.history_grid.windows(2).all(|w| ordered(w[0], w[1]))
                    || !rule.history_flow.iter().all(|q| q.is_finite()))
            {
                return Err(invalid(
                    "flow history must be finite, ordered and end at the horizon start",
                ));
            }
            let history = |a: f64| -> Result<(), Error> {
                if a < start && (rule.history_grid.is_empty() || rule.history_grid[0] > a) {
                    Err(invalid(
                        "observation window requires complete earlier flow history",
                    ))
                } else {
                    Ok(())
                }
            };
            let rates: Vec<_> = (0..nt)
                .map(|t| {
                    let mut rate = LinearForm {
                        constant: operating.observed_flow[i][t].inflow,
                        ..Default::default()
                    };
                    for &u in &units {
                        for k in mesh.cells[t].clone() {
                            rate.add_scaled(
                                &LinearForm::variable(u * nf + k, 1.),
                                mesh.dt[k] / mesh.duration[t],
                            );
                        }
                    }
                    for &e in &rivers {
                        rate.add_scaled(
                            &routing.releases[e].integral(case.grid[t], case.grid[t + 1]),
                            1. / (case.grid[t + 1] - case.grid[t]),
                        );
                    }
                    for &e in &arrivals {
                        for k in mesh.cells[t].clone() {
                            rate.add_scaled(&routing.arrivals[e][k], mesh.dt[k] / mesh.duration[t]);
                        }
                    }
                    for &e in &tunnels {
                        for k in mesh.cells[t].clone() {
                            rate.add_scaled(
                                &LinearForm::variable(tunnel_offset + e * nf + k, 1.),
                                mesh.dt[k] / mesh.duration[t],
                            );
                        }
                    }
                    rate
                })
                .collect();
            for (t, rate) in rates.iter().enumerate() {
                data.rows.push(FlowRow {
                    expression: rate.clone(),
                    bounds: operating.observed_flow[i][t].bounds,
                    scale: 100.,
                });
            }
            let integral = |a: f64, b: f64| {
                let mut volume = LinearForm::default();
                for (t, w) in case.grid.windows(2).enumerate() {
                    volume.add_scaled(&rates[t], 3600. * (b.min(w[1]) - a.max(w[0])).max(0.));
                }
                for (t, w) in rule.history_grid.windows(2).enumerate() {
                    volume.constant +=
                        3600. * (b.min(w[1]) - a.max(w[0])).max(0.) * rule.history_flow[t];
                }
                volume
            };
            if rule.min_volume.is_some() || rule.max_volume.is_some() {
                let mut total = integral(start, end);
                total.constant += rule.initial_volume;
                data.rows.push(FlowRow {
                    expression: total,
                    bounds: Bounds::new(
                        rule.min_volume.unwrap_or(f64::NEG_INFINITY),
                        rule.max_volume.unwrap_or(f64::INFINITY),
                    ),
                    scale: 1e6,
                });
            }
            for w in &rule.volume_windows {
                if !ordered(w.start_time, w.end_time)
                    || !bounds_valid(w.min_volume, w.max_volume)
                    || (w.min_volume.is_none() && w.max_volume.is_none())
                {
                    return Err(invalid(
                        "volume window needs ordered finite times and at least one consistent bound",
                    ));
                }
                if w.start_time >= end || w.end_time <= start {
                    continue;
                }
                if w.end_time > end {
                    return Err(invalid(
                        "an active volume window must finish within the horizon",
                    ));
                }
                history(w.start_time)?;
                data.rows.push(FlowRow {
                    expression: integral(w.start_time, w.end_time),
                    bounds: Bounds::new(
                        w.min_volume.unwrap_or(f64::NEG_INFINITY),
                        w.max_volume.unwrap_or(f64::INFINITY),
                    ),
                    scale: 1e6,
                });
            }
            for r in &rule.rolling_ramps {
                if !ramp_valid(r.up, r.down)
                    || !r.hours.is_finite()
                    || r.hours <= 0.
                    || !(3600. * r.hours).is_finite()
                    || !ordered(start - r.hours, start)
                    || !ordered(end, end + r.hours)
                    || !(start - r.hours).is_finite()
                    || !(end + r.hours).is_finite()
                {
                    return Err(invalid(
                        "rolling ramp needs positive finite duration and a nonnegative difference limit",
                    ));
                }
                history(start - r.hours)?;
                let mut knots = case.grid.clone();
                knots.extend(
                    rule.history_grid
                        .iter()
                        .chain(&case.grid)
                        .map(|t| t + r.hours)
                        .filter(|t| start < *t && *t < end),
                );
                knots.sort_by(f64::total_cmp);
                knots.dedup();
                for time in knots {
                    let average = integral(time - r.hours, time).scaled(1. / (3600. * r.hours));
                    let k = case
                        .grid
                        .partition_point(|v| *v <= time)
                        .saturating_sub(1)
                        .min(nt - 1);
                    let mut intervals = vec![k];
                    if time > start && time < end && case.grid[k] == time {
                        intervals.push(k - 1);
                    }
                    for k in intervals {
                        let mut change = rates[k].clone();
                        change.add_scaled(&average, -1.);
                        data.rows.push(FlowRow {
                            expression: change,
                            bounds: Bounds::new(
                                r.down.map_or(f64::NEG_INFINITY, |v| -v),
                                r.up.unwrap_or(f64::INFINITY),
                            ),
                            scale: 100.,
                        });
                    }
                }
            }
            for r in &rule.period_ramps {
                if !ramp_valid(r.up, r.down)
                    || r.grid.len() < 3
                    || !r.grid.windows(2).all(|w| ordered(w[0], w[1]))
                {
                    return Err(invalid(
                        "period ramp needs two ordered finite periods and a nonnegative difference limit",
                    ));
                }
                for w in r.grid.windows(3) {
                    let (a, b, c) = (w[0], w[1], w[2]);
                    if b >= end || c <= start {
                        continue;
                    }
                    if c > end {
                        return Err(invalid(
                            "an active ramp period must finish within the horizon",
                        ));
                    }
                    history(a)?;
                    let mut change = integral(b, c).scaled(1. / (3600. * (c - b)));
                    change.add_scaled(&integral(a, b), -1. / (3600. * (b - a)));
                    data.rows.push(FlowRow {
                        expression: change,
                        bounds: Bounds::new(
                            r.down.map_or(f64::NEG_INFINITY, |v| -v),
                            r.up.unwrap_or(f64::INFINITY),
                        ),
                        scale: 100.,
                    });
                }
            }
            data.rates.push(rates);
        }
        let finite = |expression: &LinearForm| {
            expression.constant.is_finite()
                && expression.coefficients.values().all(|c| c.is_finite())
        };
        if data.rows.iter().any(|r| !finite(&r.expression))
            || data
                .changes
                .iter()
                .any(|(expression, cost, scale)| !finite(expression) || !(cost * scale).is_finite())
        {
            return Err(Error::Input(
                "flow observation integral exceeds finite arithmetic".into(),
            ));
        }
        Ok(data)
    }
    pub fn constrain<'a>(
        &self,
        model: &oximo::Model,
        columns: &[oximo::expr::Expr<'a>],
        zero: oximo::expr::Expr<'a>,
    ) {
        use oximo::prelude::*;
        for row in &self.rows {
            let expression =
                crate::nlp::linear_expression(&row.expression, columns, zero) / row.scale;
            if row.bounds.min == row.bounds.max {
                constraint!(model, expression == row.bounds.min / row.scale);
            } else {
                if row.bounds.min.is_finite() {
                    constraint!(model, expression >= row.bounds.min / row.scale);
                }
                if row.bounds.max.is_finite() {
                    constraint!(model, expression <= row.bounds.max / row.scale);
                }
            }
        }
    }
    pub fn residual(&self, columns: &[f64]) -> f64 {
        let mut excess = self
            .rows
            .iter()
            .map(|r| {
                let value = r.expression.value(columns);
                if !value.is_finite() {
                    f64::INFINITY
                } else {
                    (r.bounds.min - value).max(value - r.bounds.max).max(0.) / r.scale
                }
            })
            .fold(0., f64::max);
        for domain in &self.domains {
            for rate in &domain.rates {
                let q = rate.value(columns);
                excess = excess.max(if q.is_finite() {
                    q.abs().min((domain.minimum - q).max(0.)) / domain.minimum.max(1.)
                } else {
                    f64::INFINITY
                });
            }
        }
        excess
    }

    /// On-state flow bounds and exact linear transition equations. The same
    /// rows accept binary proposal states or prescribed nonlinear states.
    pub fn constrain_states<'a>(
        &self,
        model: &oximo::Model,
        columns: &[oximo::expr::Expr<'a>],
        on: &[Vec<oximo::expr::Expr<'a>>],
        starts: &[Vec<oximo::expr::Expr<'a>>],
        zero: oximo::expr::Expr<'a>,
    ) -> oximo::expr::Expr<'a> {
        use oximo::prelude::*;
        let mut cost = zero;
        for (k, domain) in self.domains.iter().enumerate() {
            for t in 0..domain.rates.len() {
                let scale = domain.capacity[t].max(domain.minimum).max(1.);
                let q = crate::nlp::linear_expression(&domain.rates[t], columns, zero) / scale;
                let state = on[k][t];
                let start = starts[k][t];
                let prior = if t == 0 {
                    zero + if domain.initial_on { 1. } else { 0. }
                } else {
                    on[k][t - 1]
                };
                constraint!(model, q >= state * (domain.minimum / scale));
                constraint!(model, q <= state * (domain.capacity[t] / scale));
                constraint!(model, start >= state - prior);
                constraint!(model, start <= state);
                constraint!(model, start <= 1. - prior);
                let stop = start - state + prior;
                cost = cost + start * domain.startup_cost[t] + stop * domain.shutdown_cost[t];
            }
        }
        cost
    }

    /// Convex epigraph of total variation. A change is charged once, without
    /// multiplying by an interval duration. Zero-cost changes add no unknowns.
    pub fn cost<'a>(
        &self,
        model: &'a oximo::Model,
        columns: &[oximo::expr::Expr<'a>],
        zero: oximo::expr::Expr<'a>,
    ) -> oximo::expr::Expr<'a> {
        use oximo::prelude::*;
        variable!(model, change[k in 0..self.changes.len()] >= 0.);
        let mut cost = zero;
        for (k, (difference, coefficient, scale)) in self.changes.iter().enumerate() {
            let delta = crate::nlp::linear_expression(difference, columns, zero) / *scale;
            constraint!(model, change[[k]] >= delta);
            constraint!(model, change[[k]] >= -delta);
            cost = cost + change[[k]] * (coefficient * scale);
        }
        cost
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Boundary, Efficiency, HydroSystem, Plant, Reservoir, StorageCurve, Tunnel, Unit};
    #[test]
    fn flow_observations_keep_their_original_windows_on_finer_hydraulic_cells() {
        let system = HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(0., 1e6),
                initial_storage: 5e5,
                water_value: 0.,
                operation: Default::default(),
                level: StorageCurve::Polynomial {
                    reference_volume: 0.,
                    elevation: 100.,
                    slope: 0.,
                    curvature: 0.,
                },
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            plants: vec![Plant::new("Station", "Lake", "Sea")],
            units: vec![Unit::new(
                "Unit",
                "Station",
                Bounds::new(0., 30.),
                Bounds::new(0., 100.),
                Bounds::new(1., 100.),
                Efficiency::constant(0.9),
            )],
            tunnels: vec![Tunnel::new("Link", "Lake", "Sea", 100., 1.)],
            ..Default::default()
        };
        let mut c = ScheduleCase::new(system, vec![0., 1., 3.], vec![100.; 2]);
        let mut rule = FlowRequirement::new("Combined");
        rule.units.push("Unit".into());
        rule.tunnels.push("Link".into());
        c.flow_requirements.push(rule);
        let topology = c.topology().unwrap();
        let operating = crate::operating_data::OperatingData::new(&c, &topology).unwrap();
        let mesh =
            crate::time_discretization::TimeDiscretization::new(&c.grid, &[0., 0.25, 1., 2., 3.])
                .unwrap();
        let routing = Routing::new(&c, &topology, &operating, &mesh, &c.routing_grid).unwrap();
        let equations = FlowEquations::new(&c, &routing, &operating, &mesh).unwrap();
        let values = [4., 8., 12., 20., 1., 1., 2., 2.];
        let rates: Vec<_> = equations.rates[0]
            .iter()
            .map(|q| q.value(&values))
            .collect();
        assert!((rates[0] - 8.).abs() < 1e-12 && (rates[1] - 18.).abs() < 1e-12);
        let water: f64 = rates
            .iter()
            .zip(&mesh.duration)
            .map(|(q, dt)| 3600. * q * dt)
            .sum();
        assert!((water - 3600. * 44.).abs() < 1e-8);
    }
}
