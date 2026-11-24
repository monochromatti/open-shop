//! Independent step doubling with unchanged operating decisions.
//!
//! The differences are numerical convergence estimates, not rigorous bounds on
//! continuous-time error. Physical feasibility is checked on the finer replay.

use crate::{Error, Schedule, SchedulingProblem};
use serde::{Deserialize, Serialize};

/// Accuracy of hydraulic time integration. Volumes are m³ and levels metres.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TimeAccuracy {
    /// Absolute storage, river inventory and cumulative-arrival tolerance, m³.
    pub absolute_storage: f64,
    /// Fraction of usable reservoir storage. The smaller of this and the
    /// absolute tolerance is used, with a one-litre rounding allowance.
    pub relative_storage: f64,
    /// Maximum difference in reservoir level, m.
    pub level: f64,
    /// Objective tolerance relative to max(1, abs(original objective)).
    pub relative_objective: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;

    #[test]
    fn exact_stationary_replay_and_invalid_numerical_inputs() {
        let system = HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(0., 1000.),
                initial_storage: 500.,
                level: StorageCurve::Polynomial {
                    reference_volume: 0.,
                    elevation: 100.,
                    slope: 0.,
                    curvature: 0.,
                },
                water_value: 0.,
                operation: Default::default(),
            }],
            ..Default::default()
        };
        let case = ScheduleCase::new(system, vec![0., 0.5, 2.], vec![50., 100.]);
        let problem = SchedulingProblem::new(&case).unwrap();
        let controls = Controls {
            head_guess: None,
            commitment: vec![],
            unit_flow: vec![],
            river_opening: vec![],
            tunnel_opening: vec![],
        };
        let schedule = problem.reconstruct(&controls).unwrap();
        let (fine, error) = problem
            .validate_replay(&schedule, &TimeAccuracy::default())
            .unwrap();
        assert!(fine.validation.feasible && error.converged());
        assert_eq!(error.storage_difference, vec![0.]);
        assert_eq!(error.level_difference, vec![0.]);
        assert_eq!(error.objective_difference, 0.);
        for value in [0., -1., f64::NAN, f64::INFINITY] {
            let accuracy = TimeAccuracy {
                level: value,
                ..Default::default()
            };
            assert!(matches!(
                problem.validate_replay(&schedule, &accuracy),
                Err(Error::Input(_))
            ));
        }
        let mut invalid = schedule.clone();
        invalid.storage[0].pop();
        assert!(matches!(
            problem.validate_replay(&invalid, &TimeAccuracy::default()),
            Err(Error::Input(_))
        ));
        let mut invalid = schedule;
        invalid.storage[0][1] = f64::NAN;
        assert!(matches!(
            problem.validate_replay(&invalid, &TimeAccuracy::default()),
            Err(Error::Input(_))
        ));
    }
}
impl Default for TimeAccuracy {
    fn default() -> Self {
        Self {
            absolute_storage: 2000.,
            relative_storage: 1e-3,
            level: 1e-2,
            relative_objective: 1e-3,
        }
    }
}
impl TimeAccuracy {
    fn validate(&self) -> Result<(), Error> {
        if [
            self.absolute_storage,
            self.relative_storage,
            self.level,
            self.relative_objective,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
        {
            return Err(Error::Input(
                "time-integration tolerances must be positive and finite".into(),
            ));
        }
        Ok(())
    }
}

/// Differences between a schedule and independently reconstructed finer steps.
/// Storage, level and river inventory differences are maxima at common time
/// vertices. Arrival differences compare cumulative volumes. Volumes are m³,
/// levels metres and objective difference is in currency units.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TimeIntegrationError {
    pub storage_difference: Vec<f64>,
    pub level_difference: Vec<f64>,
    pub river_storage_difference: Vec<f64>,
    pub arrival_volume_difference: Vec<f64>,
    pub objective_difference: f64,
    /// Largest difference divided by its declared tolerance; at most one passes.
    pub normalized_max: f64,
}
impl TimeIntegrationError {
    pub fn converged(&self) -> bool {
        self.normalized_max.is_finite() && self.normalized_max <= 1.
    }
}

pub(crate) fn bisect(grid: &[f64]) -> Vec<f64> {
    let mut fine = Vec::with_capacity(2 * grid.len() - 1);
    for w in grid.windows(2) {
        fine.push(w[0]);
        let middle = w[0] + (w[1] - w[0]) / 2.;
        if w[0] < middle && middle < w[1] {
            fine.push(middle);
        }
    }
    fine.push(*grid.last().unwrap());
    fine
}

impl SchedulingProblem<'_> {
    /// Reconstruct a schedule with twice as many hydraulic and river-mixing
    /// cells, holding the original commitments, discharges and gate positions.
    /// The returned schedule has independently checked physical residuals.
    /// Its feasibility and the returned error's convergence are separate checks.
    /// Neither check establishes global optimality or a continuous-time bound.
    pub fn validate_replay(
        &self,
        schedule: &Schedule,
        accuracy: &TimeAccuracy,
    ) -> Result<(Schedule, TimeIntegrationError), Error> {
        accuracy.validate()?;
        let mesh =
            crate::time_discretization::TimeDiscretization::new(&self.case.grid, &schedule.grid)?;
        let s = &self.case.system;
        if schedule.operating_grid != self.case.grid
            || schedule.storage.len() != s.reservoirs.len()
            || schedule.river_storage.len() != s.rivers.len()
            || schedule.river_arrival_volume.len() != s.rivers.len()
            || schedule
                .river_arrival_volume
                .iter()
                .any(|row| row.len() != mesh.dt.len() || row.iter().any(|v| !v.is_finite()))
            || (schedule.routing_grid.is_empty()
                && s.rivers.iter().any(|r| r.transport.fixed_delay().is_none()))
            || schedule
                .storage
                .iter()
                .chain(&schedule.river_storage)
                .any(|row| row.len() != mesh.grid.len() || row.iter().any(|x| !x.is_finite()))
            || !schedule.objective.is_finite()
        {
            return Err(Error::Input("replay requires finite storage, river inventory and objective on the schedule grid".into()));
        }
        let grid = bisect(&schedule.grid);
        if grid == schedule.grid {
            return Err(Error::Solver(
                "hydraulic time grid cannot be refined further in floating-point arithmetic".into(),
            ));
        }
        let mixing = if schedule.routing_grid.is_empty() {
            vec![]
        } else {
            crate::time_discretization::TimeDiscretization::new(
                &self.case.grid,
                &schedule.routing_grid,
            )?;
            bisect(&schedule.routing_grid)
        };
        let problem = self.on_mesh(&grid, &mixing)?;
        let replay = problem.reconstruct(&schedule.controls()?)?;
        let common = crate::time_discretization::TimeDiscretization::new(&schedule.grid, &grid)?;
        let mut error = TimeIntegrationError {
            storage_difference: vec![0.; s.reservoirs.len()],
            level_difference: vec![0.; s.reservoirs.len()],
            river_storage_difference: vec![0.; s.rivers.len()],
            arrival_volume_difference: vec![0.; s.rivers.len()],
            objective_difference: (replay.objective - schedule.objective).abs(),
            normalized_max: 0.,
        };
        for (r, reservoir) in s.reservoirs.iter().enumerate() {
            for (j, &k) in common.vertices.iter().enumerate() {
                let (a, b) = (schedule.storage[r][j], replay.storage[r][k]);
                error.storage_difference[r] = error.storage_difference[r].max((a - b).abs());
                let difference =
                    (reservoir.level.trial_level(a) - reservoir.level.trial_level(b)).abs();
                error.level_difference[r] =
                    error.level_difference[r].max(if difference.is_finite() {
                        difference
                    } else {
                        f64::INFINITY
                    });
            }
            let storage_limit = accuracy.absolute_storage.min(
                1e-3 + accuracy.relative_storage * (reservoir.storage.max - reservoir.storage.min),
            );
            error.normalized_max = error
                .normalized_max
                .max(error.storage_difference[r] / storage_limit)
                .max(error.level_difference[r] / accuracy.level);
        }
        for r in 0..s.rivers.len() {
            let (mut coarse_arrival, mut fine_arrival) = (0., 0.);
            for (j, cells) in common.cells.iter().enumerate() {
                coarse_arrival += schedule.river_arrival_volume[r][j];
                fine_arrival += cells
                    .clone()
                    .map(|k| replay.river_arrival_volume[r][k])
                    .sum::<f64>();
                error.arrival_volume_difference[r] =
                    error.arrival_volume_difference[r].max((coarse_arrival - fine_arrival).abs());
                let k = common.vertices[j + 1];
                let difference =
                    (replay.river_storage[r][k] - schedule.river_storage[r][j + 1]).abs();
                error.river_storage_difference[r] =
                    error.river_storage_difference[r].max(difference);
            }
            error.normalized_max = error
                .normalized_max
                .max(error.river_storage_difference[r] / accuracy.absolute_storage)
                .max(error.arrival_volume_difference[r] / accuracy.absolute_storage);
        }

        error.normalized_max = error.normalized_max.max(
            error.objective_difference
                / (accuracy.relative_objective * schedule.objective.abs().max(1.)),
        );
        Ok((replay, error))
    }
}
