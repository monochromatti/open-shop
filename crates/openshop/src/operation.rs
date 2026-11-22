//! Unit commitment and interval-average discharge restrictions in physical time.

use crate::Error;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RampRates {
    /// Rate increase per hour; None leaves that direction unrestricted.
    pub up: Option<f64>,
    pub down: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct UnitOperation {
    pub initial_on: bool,
    /// Time already spent in the initial state, hours.
    pub initial_duration: f64,
    pub minimum_up_time: f64,
    pub minimum_down_time: f64,
    pub startup_cost: f64,
    pub shutdown_cost: f64,
    /// Previous interval-average discharge, m³/s. An initially off unit has zero flow.
    pub initial_flow: Option<f64>,
    pub initial_interval_hours: f64,
    /// Restrictions on changes between interval-average discharges, (m³/s)/hour.
    pub flow_ramp: RampRates,
}
impl Default for UnitOperation {
    fn default() -> Self {
        Self {
            initial_on: false,
            initial_duration: 0.,
            minimum_up_time: 0.,
            minimum_down_time: 0.,
            startup_cost: 0.,
            shutdown_cost: 0.,
            initial_flow: None,
            initial_interval_hours: 1.,
            flow_ramp: RampRates::default(),
        }
    }
}

/// Sum of the indexed on/off states times their coefficients must be >= rhs.
pub(crate) struct CommitmentInequality {
    pub coefficients: Vec<(usize, f64)>,
    pub rhs: f64,
}
impl CommitmentInequality {
    pub fn residual(&self, on: &[bool]) -> f64 {
        self.rhs
            - self
                .coefficients
                .iter()
                .map(|&(t, a)| a * if on[t] { 1. } else { 0. })
                .sum::<f64>()
    }
}
impl UnitOperation {
    pub(crate) fn validate(&self, name: &str) -> Result<(), Error> {
        if ![
            self.initial_duration,
            self.minimum_up_time,
            self.minimum_down_time,
            self.startup_cost,
            self.shutdown_cost,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.)
            || !self.initial_interval_hours.is_finite()
            || self.initial_interval_hours <= 0.
            || !self
                .initial_flow
                .iter()
                .chain(self.flow_ramp.up.iter())
                .chain(self.flow_ramp.down.iter())
                .all(|v| v.is_finite() && *v >= 0.)
            || (!self.initial_on && self.initial_flow.is_some_and(|q| q != 0.))
        {
            return Err(Error::Input(format!(
                "{name}: invalid unit operating restrictions or initial state"
            )));
        }
        Ok(())
    }
    pub(crate) fn commitment_inequalities(&self, grid: &[f64]) -> Vec<CommitmentInequality> {
        let nt = grid.len() - 1;
        let mut rows = Vec::new();
        let minimum = if self.initial_on {
            self.minimum_up_time
        } else {
            self.minimum_down_time
        };
        for t in 0..nt {
            if grid[t] - grid[0] + self.initial_duration < minimum {
                rows.push(CommitmentInequality {
                    coefficients: vec![(t, if self.initial_on { 1. } else { -1. })],
                    rhs: if self.initial_on { 1. } else { 0. },
                });
            }
            for k in t + 1..nt {
                if grid[k] - grid[t] >= self.minimum_up_time.max(self.minimum_down_time) {
                    break;
                }
                if grid[k] - grid[t] < self.minimum_up_time {
                    let mut coefficients = vec![(k, 1.), (t, -1.)];
                    let rhs = if t == 0 {
                        if self.initial_on { -1. } else { 0. }
                    } else {
                        coefficients.push((t - 1, 1.));
                        0.
                    };
                    rows.push(CommitmentInequality { coefficients, rhs });
                }
                if grid[k] - grid[t] < self.minimum_down_time {
                    let mut coefficients = vec![(k, -1.), (t, 1.)];
                    let rhs = if t == 0 {
                        if self.initial_on { 0. } else { -1. }
                    } else {
                        coefficients.push((t - 1, -1.));
                        -1.
                    };
                    rows.push(CommitmentInequality { coefficients, rhs });
                }
            }
        }
        rows
    }
    pub(crate) fn previous_flow(&self) -> Option<f64> {
        if self.initial_on {
            self.initial_flow
        } else {
            Some(0.)
        }
    }
    pub(crate) fn ramp_hours(&self, grid: &[f64], t: usize) -> f64 {
        let previous = if t == 0 {
            self.initial_interval_hours
        } else {
            grid[t] - grid[t - 1]
        };
        (previous + (grid[t + 1] - grid[t])) / 2.
    }
    pub(crate) fn transitions(&self, on: &[bool], t: usize) -> (bool, bool) {
        let previous = if t == 0 { self.initial_on } else { on[t - 1] };
        (on[t] && !previous, previous && !on[t])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn minimum_time_inequalities_match_elapsed_time_for_every_schedule() {
        let grid = [0., 0.25, 0.75, 1.5, 2.5, 4., 5.];
        for initial_on in [false, true] {
            for initial_duration in [0., 0.5, 3.] {
                let op = UnitOperation {
                    initial_on,
                    initial_duration,
                    minimum_up_time: 1.5,
                    minimum_down_time: 1.,
                    ..Default::default()
                };
                let rows = op.commitment_inequalities(&grid);
                for bits in 0..64 {
                    let on: Vec<_> = (0..6).map(|t| bits & (1 << t) != 0).collect();
                    let mut state = initial_on;
                    let mut changed = -initial_duration;
                    let mut physically_allowed = true;
                    for t in 0..6 {
                        if state != on[t] {
                            let required = if state {
                                op.minimum_up_time
                            } else {
                                op.minimum_down_time
                            };
                            physically_allowed &= grid[t] - changed >= required;
                            state = on[t];
                            changed = grid[t];
                        }
                    }
                    assert_eq!(
                        rows.iter().all(|r| r.residual(&on) <= 0.),
                        physically_allowed,
                        "{initial_on} {initial_duration} {on:?}"
                    );
                }
            }
        }
    }
}
