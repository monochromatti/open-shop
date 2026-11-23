//! Environmental arrival requirements on the represented release waves.
//!
//! Uniform delay bins produce piecewise-linear arrival rates. Their minima
//! occur at shifted release edges and delay-bin edges; fixed delays also
//! require both sides of a jump. The equations reuse the transport's reference
//! fluxes, rather than introducing new water or independent mixing weights.

use crate::{
    Bounds, Error, ScheduleCase,
    flow_requirements::FlowRow,
    operating_data::OperatingData,
    river::{LinearForm, ReleaseSeries, transport_series},
    river_transport::{Observation, TimeSide},
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArrivalPolicy {
    /// Minimum duration-weighted arrival rate on each observation window.
    #[default]
    IntervalAverage,
    /// Minimum instantaneous rate throughout the horizon, including jumps.
    Pointwise,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RiverOperation {
    /// Minimum downstream arrival, m³/s; includes historical releases.
    pub minimum_arrival: f64,
    pub arrival_policy: ArrivalPolicy,
    /// Observation-window edges, hours. Empty uses the control grid.
    /// A supplied grid must span the horizon; it need not match control edges.
    pub arrival_window_grid: Vec<f64>,
    /// Changes in upstream interval-average discharge, (m³/s)/hour.
    pub flow_ramp: crate::RampRates,
    /// Previous upstream interval-average release, m³/s. If omitted, use the
    /// last historical release and its duration; absent history skips this ramp.
    pub initial_release: Option<f64>,
    pub initial_interval_hours: f64,
    /// Currency per m³/s of absolute change, charged once per control change.
    pub flow_change_cost: f64,
}

impl Default for RiverOperation {
    fn default() -> Self {
        Self {
            minimum_arrival: 0.,
            arrival_policy: Default::default(),
            arrival_window_grid: vec![],
            flow_ramp: Default::default(),
            initial_release: None,
            initial_interval_hours: 1.,
            flow_change_cost: 0.,
        }
    }
}

impl crate::River {
    pub(crate) fn previous_release(&self) -> Option<(f64, f64)> {
        self.operation
            .initial_release
            .map(|q| (q, self.operation.initial_interval_hours))
            .or_else(|| {
                self.history_release
                    .last()
                    .zip(self.history_grid.windows(2).last())
                    .map(|(&q, w)| (q, w[1] - w[0]))
            })
    }
}

impl RiverOperation {
    pub(crate) fn validate(&self, grid: &[f64]) -> Result<(), Error> {
        if !self.minimum_arrival.is_finite()
            || self.minimum_arrival < 0.
            || !self.initial_interval_hours.is_finite()
            || self.initial_interval_hours <= 0.
            || self.initial_release.is_some_and(|q| !q.is_finite())
            || !self.flow_change_cost.is_finite()
            || self.flow_change_cost < 0.
            || self
                .flow_ramp
                .up
                .into_iter()
                .chain(self.flow_ramp.down)
                .any(|q| !q.is_finite() || q < 0.)
            || (!self.arrival_window_grid.is_empty()
                && (self.arrival_window_grid.len() < 2
                    || self.arrival_window_grid.first() != grid.first()
                    || self.arrival_window_grid.last() != grid.last()
                    || !self.arrival_window_grid.iter().all(|t| t.is_finite())
                    || !self
                        .arrival_window_grid
                        .windows(2)
                        .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())))
        {
            return Err(Error::Input(
                "river limits, ramp rates and costs must be finite and nonnegative; observation windows must span the horizon with finite ordered edges".into(),
            ));
        }
        Ok(())
    }
}

/// Compile linear observations of the already lifted reference fluxes.
/// The constraints are linear in these columns even when the release-dependent
/// mixture defining each reference flux is nonlinear.
pub(crate) fn arrival_limits(
    case: &ScheduleCase,
    operating: &OperatingData,
    river: usize,
    release: (&ReleaseSeries, &[Vec<LinearForm>]),
    history: (&ReleaseSeries, &[Vec<LinearForm>]),
) -> Vec<FlowRow> {
    let r = &case.system.rivers[river];
    let minimum = &operating.river_min_arrival[river];
    if minimum.iter().all(|q| *q == 0.) {
        return vec![];
    }
    let observed = |observation| {
        let mut value = transport_series(&r.transport, release.0, release.1, observation);
        value.add_scaled(
            &transport_series(&r.transport, history.0, history.1, observation),
            1.,
        );
        value
    };
    let nt = case.prices.len();
    let (start, end) = (case.grid[0], case.grid[nt]);
    let mut rows = vec![];
    let mut add = |expression: LinearForm, lower: f64| {
        if lower > 0. {
            rows.push(FlowRow {
                expression,
                bounds: Bounds::new(lower, f64::INFINITY),
                scale: r.capacity.max(1.),
            });
        }
    };
    match r.operation.arrival_policy {
        ArrivalPolicy::IntervalAverage => {
            let grid = if r.operation.arrival_window_grid.is_empty() {
                &case.grid
            } else {
                &r.operation.arrival_window_grid
            };
            for w in grid.windows(2) {
                let required: f64 = case
                    .grid
                    .windows(2)
                    .zip(minimum)
                    .map(|(cell, q)| {
                        q * (w[1].min(cell[1]) - w[0].max(cell[0])).max(0.) / (w[1] - w[0])
                    })
                    .sum();
                add(
                    observed(Observation::Arrival(w[0], w[1])).scaled(1. / (w[1] - w[0])),
                    required,
                );
            }
        }
        ArrivalPolicy::Pointwise => {
            let delays = match &r.transport {
                crate::RiverTransport::Fixed { delay } => vec![*delay],
                crate::RiverTransport::Distributed { curves } => curves
                    .iter()
                    .flat_map(|c| c.edges().iter().copied())
                    .collect(),
            };
            // Include control edges even when no arrival changes there: an
            // environmental minimum itself can change at that time.
            let mut knots = case.grid.clone();
            for &t in release.0.grid.iter().chain(&history.0.grid) {
                for &delay in &delays {
                    let time = t + delay;
                    if start <= time && time <= end {
                        knots.push(time);
                    }
                }
            }
            knots.sort_by(f64::total_cmp);
            knots.dedup();
            for time in knots {
                for side in [TimeSide::Left, TimeSide::Right] {
                    if (time == start && matches!(side, TimeSide::Left))
                        || (time == end && matches!(side, TimeSide::Right))
                    {
                        continue;
                    }
                    let t = case
                        .grid
                        .partition_point(|v| match side {
                            TimeSide::Left => *v < time,
                            TimeSide::Right => *v <= time,
                        })
                        .saturating_sub(1)
                        .min(nt - 1);
                    add(observed(Observation::Point(time, side)), minimum[t]);
                }
            }
        }
    }
    rows
}
