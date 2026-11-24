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
    pub(crate) fn validate(&self) -> Result<(), Error> {
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
    /// Largest normalized error estimate; at most one passes. Hydraulic step
    /// doubling uses the second-order Richardson factor 4/3. A distributed
    /// river mesh uses factor two, allowing first-order mixing convergence.
    pub normalized_max: f64,
}
impl TimeIntegrationError {
    pub(crate) fn between(
        case: &crate::ScheduleCase,
        coarse: &Schedule,
        fine: &Schedule,
        accuracy: &TimeAccuracy,
    ) -> Result<Self, Error> {
        let common = crate::time_discretization::TimeDiscretization::new(&coarse.grid, &fine.grid)?;
        let s = &case.system;
        let mut error = TimeIntegrationError {
            storage_difference: vec![0.; s.reservoirs.len()],
            level_difference: vec![0.; s.reservoirs.len()],
            river_storage_difference: vec![0.; s.rivers.len()],
            arrival_volume_difference: vec![0.; s.rivers.len()],
            objective_difference: (fine.objective - coarse.objective).abs(),
            normalized_max: 0.,
        };
        for (r, reservoir) in s.reservoirs.iter().enumerate() {
            for (j, &k) in common.vertices.iter().enumerate() {
                let (a, b) = (coarse.storage[r][j], fine.storage[r][k]);
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
                coarse_arrival += coarse.river_arrival_volume[r][j];
                fine_arrival += cells
                    .clone()
                    .map(|k| fine.river_arrival_volume[r][k])
                    .sum::<f64>();
                error.arrival_volume_difference[r] =
                    error.arrival_volume_difference[r].max((coarse_arrival - fine_arrival).abs());
                let k = common.vertices[j + 1];
                let difference = (fine.river_storage[r][k] - coarse.river_storage[r][j + 1]).abs();
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
                / (accuracy.relative_objective * coarse.objective.abs().max(1.)),
        );
        error.normalized_max *= if coarse.routing_grid.is_empty() {
            4. / 3.
        } else {
            2.
        };
        Ok(error)
    }
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

/// Dörfler marking: bisect cells carrying half the squared local error.
/// Scaling before squaring avoids overflow and leaves the marked set unchanged.
pub(crate) fn bulk_refined_grid(grid: &[f64], defects: &[f64]) -> Vec<f64> {
    let largest = defects.iter().copied().fold(0., f64::max);
    if largest == 0. {
        return grid.to_vec();
    }
    let weights: Vec<_> = defects
        .iter()
        .map(|&d| {
            if largest.is_infinite() {
                if d.is_infinite() { 1. } else { 0. }
            } else {
                (d / largest).powi(2)
            }
        })
        .collect();
    let total: f64 = weights.iter().sum();
    let mut order: Vec<_> = (0..defects.len()).collect();
    order.sort_by(|&a, &b| weights[b].total_cmp(&weights[a]).then(a.cmp(&b)));
    let mut marked = vec![false; defects.len()];
    let mut covered = 0.;
    for k in order {
        if weights[k] == 0. {
            break;
        }
        marked[k] = true;
        covered += weights[k];
        if largest.is_finite() && covered >= total / 2. {
            break;
        }
    }
    let mut fine = Vec::with_capacity(grid.len() + marked.iter().filter(|&&m| m).count());
    for (k, w) in grid.windows(2).enumerate() {
        fine.push(w[0]);
        let middle = w[0] + (w[1] - w[0]) / 2.;
        if marked[k] && w[0] < middle && middle < w[1] {
            fine.push(middle);
        }
    }
    fine.push(*grid.last().unwrap());
    fine
}

impl SchedulingProblem<'_> {
    /// Two half steps start from each represented storage state. This locates
    /// local integration error without propagating earlier state errors through
    /// every later cell. Whole-horizon replay determines numerical acceptance.
    pub(crate) fn hydraulic_defects(
        &self,
        schedule: &Schedule,
        accuracy: &TimeAccuracy,
    ) -> Result<Vec<f64>, Error> {
        let grid = bisect(&schedule.grid);
        let fine = self.on_mesh(&grid, &schedule.routing_grid)?;
        let controls = schedule.controls()?;
        let point = crate::commitment::Proposal::from_schedule(self, schedule)?;
        let initial = point.on_grid(&grid)?;
        let s = &self.case.system;
        let nt = grid.len() - 1;
        let base: Vec<_> = initial
            .flow
            .iter()
            .flatten()
            .copied()
            .chain(initial.tunnel_flow.iter().flatten().copied())
            .chain(
                fine.routing
                    .outlets
                    .iter()
                    .flat_map(|&e| initial.river_withdrawal[e].iter().copied()),
            )
            .collect();
        let mut columns = base.clone();
        let cells = crate::time_discretization::TimeDiscretization::new(&schedule.grid, &grid)?;
        let mut defects = vec![0.; schedule.grid.len() - 1];
        for (t, children) in cells.cells.iter().enumerate() {
            columns.copy_from_slice(&base);
            let mut previous: Vec<_> = schedule.storage.iter().map(|row| row[t]).collect();
            let mut energy = 0.;
            for k in children.clone() {
                let state = match crate::hydraulics::hydraulic_step(
                    &fine,
                    &controls,
                    k,
                    &previous,
                    3600. * fine.mesh.dt[k],
                    &columns,
                    None,
                ) {
                    Ok(state) => state,
                    Err(Error::Solver(_)) => {
                        defects[t] = f64::INFINITY;
                        break;
                    }
                    Err(error) => return Err(error),
                };
                for (g, u) in s.units.iter().enumerate() {
                    let head =
                        crate::hydraulics::machine_head(&fine, g, &state.head, &state.unit_flow);
                    energy += u.mode.power_sign()
                        * u.trial_power(state.unit_flow[g], head)
                        * fine.mesh.dt[k];
                    columns[g * nt + k] = state.unit_flow[g];
                }
                for e in 0..s.tunnels.len() {
                    columns[(s.units.len() + e) * nt + k] = state.tunnel_flow[e];
                }
                for j in 0..fine.routing.outlets.len() {
                    columns[(s.units.len() + s.tunnels.len() + j) * nt + k] =
                        state.river_withdrawal[j];
                }
                // A changing emitted flux can affect distant arrivals even
                // when its local storage defect is small.
                for i in (k..fine.routing.primary_columns).step_by(nt) {
                    let volume = 3600. * fine.mesh.dt[k] * (columns[i] - base[i]).abs();
                    defects[t] = defects[t].max(volume / (accuracy.absolute_storage / 10.));
                }
                previous = state.storage;
            }
            if defects[t].is_infinite() {
                continue;
            }
            for (r, reservoir) in s.reservoirs.iter().enumerate() {
                let difference = (previous[r] - schedule.storage[r][t + 1]).abs();
                let level = (reservoir.level.trial_level(previous[r])
                    - reservoir.level.trial_level(schedule.storage[r][t + 1]))
                .abs();
                let limit = accuracy.absolute_storage.min(
                    1e-3 + accuracy.relative_storage
                        * (reservoir.storage.max - reservoir.storage.min),
                );
                defects[t] = defects[t]
                    .max(difference / limit)
                    .max(if level.is_finite() {
                        level / accuracy.level
                    } else {
                        f64::INFINITY
                    });
            }
            let original: f64 = s
                .units
                .iter()
                .enumerate()
                .map(|(g, u)| u.mode.power_sign() * schedule.unit_power[g][t])
                .sum();
            let objective = ((energy - original * cells.duration[t])
                * self.case.prices[self.mesh.interval[t]])
                .abs();
            defects[t] = defects[t].max(if objective.is_finite() {
                objective / (accuracy.relative_objective * schedule.objective.abs().max(1.))
            } else {
                f64::INFINITY
            });
        }
        Ok(defects)
    }

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
        let controls = self.held_controls(schedule, &mesh)?;
        let problem = self.on_mesh(&grid, &mixing)?;
        let replay = problem.reconstruct(&controls)?;
        let mut error = TimeIntegrationError::between(self.case, schedule, &replay, accuracy)?;
        if !mixing.is_empty() {
            // Isolate river-mixing error so it cannot cancel a hydraulic error.
            let transport = self
                .on_mesh(&schedule.grid, &mixing)?
                .reconstruct(&controls)?;
            let transport_error =
                TimeIntegrationError::between(self.case, schedule, &transport, accuracy)?;
            for (a, b) in error
                .storage_difference
                .iter_mut()
                .zip(transport_error.storage_difference)
            {
                *a = a.max(b);
            }
            for (a, b) in error
                .level_difference
                .iter_mut()
                .zip(transport_error.level_difference)
            {
                *a = a.max(b);
            }
            for (a, b) in error
                .river_storage_difference
                .iter_mut()
                .zip(transport_error.river_storage_difference)
            {
                *a = a.max(b);
            }
            for (a, b) in error
                .arrival_volume_difference
                .iter_mut()
                .zip(transport_error.arrival_volume_difference)
            {
                *a = a.max(b);
            }
            error.objective_difference = error
                .objective_difference
                .max(transport_error.objective_difference);
            error.normalized_max = error.normalized_max.max(transport_error.normalized_max);
        }
        Ok((replay, error))
    }
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
