//! Nonlinear dispatch and independent evaluation of its physical equations.

use crate::dispatch::{DispatchObjective, build_dispatch};
use crate::trajectory::{Schedule, Validation};
use crate::{
    Error, ScheduleCase, SolveOptions,
    model::Topology,
    nlp::{self, Backend},
    river::Routing,
};
use std::{sync::Arc, time::Instant};

/// Optimize unit states and tunnel gates with SCIP proposals and nonlinear dispatch.
/// Prescribed commitments remain fixed; other units can change commitment
/// between nonlinear solves. This search does not establish a global bound.
/// An infeasible or unfinished solve is reported through `validation.feasible`;
/// a solver's success flag alone does not establish physical feasibility.
pub struct SchedulingProblem<'a> {
    pub(crate) case: &'a ScheduleCase,
    pub(crate) mesh: crate::time_discretization::TimeDiscretization,
    pub(crate) topology: Arc<Topology>,
    pub(crate) routing: Routing,
    pub(crate) operating: Arc<crate::operating_data::OperatingData>,
    pub(crate) observations: crate::flow_requirements::FlowEquations,
    pub(crate) level_windows: crate::level_windows::LevelEquations,
    pub(crate) interlocks: Arc<crate::gate_interlocks::GateInterlocks>,
    pub(crate) operating_points: crate::operating_points::OperatingPoints,
    pub(crate) soft_limits: Arc<crate::soft_constraints::SoftLimits>,
    pub(crate) water_value: crate::river::LinearForm,
}

impl<'a> SchedulingProblem<'a> {
    /// Validate inputs and prepare network incidence and transport once.
    pub fn new(case: &'a ScheduleCase) -> Result<Self, Error> {
        let topology = Arc::new(case.topology()?);
        let mesh = crate::time_discretization::TimeDiscretization::new(&case.grid, &case.grid)?;
        let soft_limits = Arc::new(crate::soft_constraints::SoftLimits::new(case, &topology)?);
        let operating = Arc::new(crate::operating_data::OperatingData::new(case, &topology)?);
        let routing = Routing::new(case, &topology, &operating, &mesh, &case.routing_grid)?;
        let observations =
            crate::flow_requirements::FlowEquations::new(case, &routing, &operating, &mesh)?;
        let level_windows = crate::level_windows::LevelEquations::new(case, &routing, &mesh)?;
        let interlocks = Arc::new(crate::gate_interlocks::GateInterlocks::new(
            case, &topology,
        )?);
        let operating_points = crate::operating_points::OperatingPoints::new(case, &routing, &mesh);
        let water_value =
            crate::trajectory::water_value_change(case, &topology, &routing, &operating, &mesh);
        Ok(Self {
            mesh,
            soft_limits,
            water_value,
            operating_points,
            interlocks,
            level_windows,
            observations,
            case,
            topology,
            routing,
            operating,
        })
    }
    pub(crate) fn on_mesh(&self, grid: &[f64], mixing_grid: &[f64]) -> Result<Self, Error> {
        let case = self.case;
        let mesh = crate::time_discretization::TimeDiscretization::new(&case.grid, grid)?;
        let routing = Routing::new(case, &self.topology, &self.operating, &mesh, mixing_grid)?;
        let observations =
            crate::flow_requirements::FlowEquations::new(case, &routing, &self.operating, &mesh)?;
        let level_windows = crate::level_windows::LevelEquations::new(case, &routing, &mesh)?;
        let operating_points = crate::operating_points::OperatingPoints::new(case, &routing, &mesh);
        let water_value = crate::trajectory::water_value_change(
            case,
            &self.topology,
            &routing,
            &self.operating,
            &mesh,
        );
        Ok(Self {
            case,
            mesh,
            routing,
            observations,
            level_windows,
            operating_points,
            water_value,
            topology: Arc::clone(&self.topology),
            operating: Arc::clone(&self.operating),
            interlocks: Arc::clone(&self.interlocks),
            soft_limits: Arc::clone(&self.soft_limits),
        })
    }
    pub fn case(&self) -> &ScheduleCase {
        self.case
    }
    pub fn solve(&self, options: &SolveOptions) -> Result<Schedule, Error> {
        self.solve_with_backend(options, Backend::Ipopt)
    }
    /// Improve continuous dispatch from physical controls while holding unit
    /// commitments, river flow states and selectable tunnel positions fixed. A feasible supplied
    /// schedule is retained if refinement fails or reduces its revenue.
    /// This local refinement does not establish global optimality.
    pub fn refine(
        &self,
        controls: &crate::Controls,
        options: &SolveOptions,
    ) -> Result<Schedule, Error> {
        options.validate()?;
        let started = Instant::now();
        let mut physical = self.reconstruct(controls)?;
        let start = crate::commitment::Proposal::from_schedule(self, &physical)?;
        let mut budget = options.clone();
        budget.time_limit -= started.elapsed().as_secs_f64();
        if budget.time_limit <= 0.01 {
            physical.validation.feasible = false;
            return Ok(physical);
        }
        let refined = self.dispatch(
            start.controls.clone(),
            Some(&start),
            &budget,
            Backend::Ipopt,
        )?;
        Ok(refined)
    }

    fn solve_with_backend(
        &self,
        options: &SolveOptions,
        backend: Backend,
    ) -> Result<Schedule, Error> {
        options.validate()?;
        let started = Instant::now();
        let mut excluded = Vec::new();
        let mut dispatched = Vec::new();
        let mut best: Option<Schedule> = None;
        let mut last = None;
        while options.time_limit - started.elapsed().as_secs_f64() > 0.01 {
            let remaining = options.time_limit - started.elapsed().as_secs_f64();
            let Some(proposal) = crate::commitment::propose(
                self,
                best.as_ref(),
                &excluded,
                remaining.min(options.time_limit / 4.),
            )?
            else {
                break;
            };
            // Re-solving an unchanged integer assignment adds no new search
            // direction. Use the remaining allowance to refine its dispatch.
            if dispatched.contains(&proposal.controls) {
                break;
            }
            dispatched.push(proposal.controls.clone());
            let mut budget = options.clone();
            budget.time_limit =
                (options.time_limit - started.elapsed().as_secs_f64()).min(options.time_limit / 4.);
            if budget.time_limit <= 0. {
                break;
            }
            budget.max_iterations = options.max_iterations.min(200);
            let result =
                self.dispatch(proposal.controls.clone(), Some(&proposal), &budget, backend)?;
            if result.validation.feasible {
                if best.as_ref().is_none_or(|p| result.objective > p.objective) {
                    best = Some(result.clone());
                } else {
                    break;
                }
            } else {
                // This excludes a failed local trial only from this proposal
                // search; it is not a proof that the assignment is infeasible.
                excluded.push(proposal.controls.clone());
            }
            last = Some(result);
        }
        // An infeasible linearization does not prove nonlinear infeasibility.
        // If every discrete control is prescribed, solve those equations directly.
        if last.is_none() {
            let unit_on = self
                .case
                .system
                .units
                .iter()
                .enumerate()
                .map(|(g, u)| {
                    self.operating.units[g]
                        .iter()
                        .enumerate()
                        .map(|(t, interval)| {
                            self.case
                                .commitment
                                .get(&u.name)
                                .map(|on| on[t])
                                .or(interval.forced)
                        })
                        .collect::<Option<Vec<_>>>()
                })
                .collect::<Option<Vec<_>>>();
            let gate_position = self
                .operating
                .tunnels
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|interval| (interval.positions.len() == 1).then_some(0))
                        .collect::<Option<Vec<_>>>()
                })
                .collect::<Option<Vec<_>>>();
            if let (Some(unit_on), Some(gate_position)) = (unit_on, gate_position)
                && self.observations.domains.is_empty()
            {
                let mut budget = options.clone();
                budget.time_limit -= started.elapsed().as_secs_f64();
                if budget.time_limit > 0. {
                    return self.dispatch(
                        crate::commitment::DiscreteControls {
                            unit_on,
                            gate_position,
                            river_on: vec![],
                        },
                        None,
                        &budget,
                        backend,
                    );
                }
            }
        }
        let result = best.or(last).ok_or_else(|| {
            Error::Solver(
                "no unit and gate schedule produced a nonlinear trial point within the time limit"
                    .into(),
            )
        })?;
        let remaining = options.time_limit - started.elapsed().as_secs_f64();
        if remaining > 0.01 {
            let start = crate::commitment::Proposal::from_schedule(self, &result)?;
            let mut budget = options.clone();
            budget.time_limit = remaining;
            let refined = self.dispatch(start.controls.clone(), Some(&start), &budget, backend)?;
            if refined.validation.feasible
                && (!result.validation.feasible || refined.objective >= result.objective)
            {
                return Ok(refined);
            }
        }
        Ok(result)
    }
}

impl SchedulingProblem<'_> {
    fn dispatch(
        &self,
        controls: crate::commitment::DiscreteControls,
        start: Option<&crate::commitment::Proposal>,
        options: &SolveOptions,
        backend: Backend,
    ) -> Result<Schedule, Error> {
        let started = Instant::now();
        let mut mesh = start
            .filter(|p| p.grid != self.mesh.grid || p.routing_grid != self.routing.mixing_grid)
            .map(|p| self.on_mesh(&p.grid, &p.routing_grid))
            .transpose()?;
        let mut warm = None;
        let mut current: Option<Schedule> = None;
        let mut best: Option<Schedule> = None;
        let mut last = None;
        let mut polished = false;
        loop {
            let problem = mesh.as_ref().unwrap_or(self);
            let mut budget = options.clone();
            budget.time_limit -= started.elapsed().as_secs_f64();
            let restoring = current.is_none();
            if restoring {
                if budget.time_limit <= 0.01 {
                    return best.or(last).ok_or_else(|| Error::Solver("dispatch time allowance expired before a physical trajectory was available".into()));
                }
                let initial = warm.as_ref().or(start);
                let physical = if let Some(point) = initial {
                    match problem.reconstruct(&point.physical_controls(problem)?) {
                        Ok(s) => Some(s),
                        Err(Error::Input(_) | Error::Solver(_)) => None,
                        Err(error) => return Err(error),
                    }
                } else {
                    None
                };
                let projected = physical
                    .as_ref()
                    .map(|s| crate::commitment::Proposal::from_schedule(problem, s))
                    .transpose()?;
                let feasible = physical.as_ref().is_some_and(|s| s.validation.feasible)
                    && projected.as_ref().is_some_and(|p| p.controls == controls);
                budget.time_limit = options.time_limit - started.elapsed().as_secs_f64();
                current = if feasible {
                    physical
                } else if budget.time_limit > 0.01 {
                    // Restore the water plan before optimizing revenue or
                    // interpreting integration differences as mesh error.
                    Some(problem.dispatch_with_objective(
                        controls.clone(),
                        projected.as_ref().or(initial),
                        &budget,
                        backend,
                        DispatchObjective::FlowDistance,
                    )?)
                } else {
                    physical
                };
            }
            let Some(mut result) = current.take() else {
                return best.or(last).ok_or_else(|| {
                    Error::Solver("dispatch time allowance expired during reconstruction".into())
                });
            };
            // Step differences estimate discretization error only after the
            // algebraic water and pressure equations have been solved. A
            // violated operating bound alone does not prevent that estimate.
            let tolerance = options.tolerance.max(1e-7);
            if result.validation.water_balance > tolerance
                || result.validation.pressure_balance > 100. * tolerance
            {
                return Ok(best.unwrap_or(result));
            }
            if !result.validation.feasible
                && !restoring
                && options.time_limit - started.elapsed().as_secs_f64() > 0.01
            {
                // Restore the improved water plan with its current crossing
                // observations before accepting or discarding its revenue.
                warm = Some(crate::commitment::Proposal::from_schedule(
                    problem, &result,
                )?);
                last = Some(result);
                continue;
            }
            let replay = match problem.validate_replay(&result, &options.time_accuracy) {
                Ok((fine, error)) => {
                    let converged =
                        result.validation.feasible && fine.validation.feasible && error.converged();
                    if !result.validation.feasible && !fine.validation.feasible && error.converged()
                    {
                        return Ok(best.unwrap_or(result));
                    }
                    result.time_error = Some(error);
                    Some((fine, converged))
                }
                Err(Error::Solver(_)) => None,
                Err(error) => return Err(error),
            };
            let converged = replay.as_ref().is_some_and(|(_, converged)| *converged);
            if converged {
                let improved = best.as_ref().is_none_or(|s| {
                    result.objective > s.objective + options.tolerance * s.objective.abs().max(1.)
                });
                if best
                    .as_ref()
                    .is_none_or(|s| result.objective >= s.objective)
                {
                    best = Some(result.clone());
                }
                budget.time_limit = options.time_limit - started.elapsed().as_secs_f64();
                if budget.time_limit <= 0.01
                    || (polished && (!improved || matches!(result.solver_status, Some(0 | 1))))
                {
                    return Ok(best.unwrap());
                }
                let point = crate::commitment::Proposal::from_schedule(problem, &result)?;
                // Revenue solves yield checkpoints while time remains to
                // reconstruct their controls and restore finer-grid feasibility.
                budget.time_limit /= 2.;
                current = match problem.dispatch_with_objective(
                    controls.clone(),
                    Some(&point),
                    &budget,
                    backend,
                    DispatchObjective::Revenue,
                ) {
                    Ok(s) => Some(s),
                    Err(Error::Solver(_)) => return Ok(best.unwrap()),
                    Err(error) => return Err(error),
                };
                polished = true;
                continue;
            }
            // A finite-grid solution is not an accepted time trajectory until
            // the independently finer physics and numerical comparison pass.
            result.validation.feasible = false;
            if options.time_limit - started.elapsed().as_secs_f64() <= 0.01 {
                return Ok(best.unwrap_or(result));
            }
            let mut grid = if result.time_error.as_ref().is_some_and(|e| !e.converged()) {
                let defects = problem.hydraulic_defects(&result, &options.time_accuracy)?;
                crate::time_integration::bulk_refined_grid(&result.grid, &defects)
            } else {
                crate::time_integration::bisect(&result.grid)
            };
            let mut mixing = result.routing_grid.clone();
            if !mixing.is_empty() && result.time_error.as_ref().is_some_and(|e| !e.converged()) {
                let finer_mixing = crate::time_integration::bisect(&mixing);
                let transport = self
                    .on_mesh(&result.grid, &finer_mixing)?
                    .reconstruct(&result.controls()?);
                match transport {
                    Ok(transport) => {
                        let error = crate::TimeIntegrationError::between(
                            self.case,
                            &result,
                            &transport,
                            &options.time_accuracy,
                        )?;
                        if !error.converged() {
                            // Resolve river averaging before adding hydraulic
                            // states that cannot correct a transport error.
                            grid.clone_from(&result.grid);
                            mixing = finer_mixing;
                        }
                    }
                    Err(Error::Solver(_)) => mixing = finer_mixing,
                    Err(error) => return Err(error),
                }
            }
            if grid == result.grid && mixing == result.routing_grid {
                grid = crate::time_integration::bisect(&result.grid);
            }
            if grid == result.grid && mixing == result.routing_grid {
                return Ok(best.unwrap_or(result));
            }
            let initial = replay.as_ref().map_or(&result, |(fine, _)| fine);
            warm = Some(crate::commitment::Proposal::from_schedule(
                problem, initial,
            )?);
            mesh = Some(self.on_mesh(&grid, &mixing)?);
            last = Some(result);
            polished = false;
        }
    }

    fn dispatch_with_objective(
        &self,
        controls: crate::commitment::DiscreteControls,
        start: Option<&crate::commitment::Proposal>,
        options: &SolveOptions,
        backend: Backend,
        objective: DispatchObjective,
    ) -> Result<Schedule, Error> {
        let started = Instant::now();
        let initial = start.map(|s| s.on_grid(&self.mesh.grid)).transpose()?;
        let start = initial.as_deref();
        let mut points = self.operating_points.clone();
        if let Some(start) = start {
            let mut columns: Vec<_> = start
                .flow
                .iter()
                .flatten()
                .copied()
                .chain(start.tunnel_flow.iter().flatten().copied())
                .chain(
                    self.routing
                        .outlets
                        .iter()
                        .flat_map(|&e| start.river_withdrawal[e].iter().copied()),
                )
                .collect();
            self.routing.complete_columns(&mut columns);
            points = points.with_level_crossings(self, &start.storage, &columns);
        }
        let mut warm = None;
        let mut last = None;
        let mut best_trial: Option<Schedule> = None;
        loop {
            let mut budget = options.clone();
            budget.time_limit -= started.elapsed().as_secs_f64();
            if budget.time_limit <= 0. {
                return best_trial.or(last).ok_or_else(|| Error::Solver(
                    "nonlinear dispatch exceeded its time allowance before producing a trial point".into()));
            }
            let result = self.dispatch_on_points(
                controls.clone(),
                warm.as_ref().or(start),
                &points,
                &budget,
                backend,
                objective,
            )?;
            let tolerance = options.tolerance.max(1e-7);
            if objective == DispatchObjective::Revenue
                && result.validation.water_balance <= tolerance
                && result.validation.pressure_balance <= 100. * tolerance
                && result.objective.is_finite()
                && best_trial
                    .as_ref()
                    .is_none_or(|s| result.objective > s.objective)
            {
                // Keep a useful water plan for restoration if a later solve
                // with added crossing observations is interrupted. Its physical
                // validation remains unchanged; this does not accept the trial.
                best_trial = Some(result.clone());
            }
            if result.validation.feasible {
                return Ok(result);
            }
            if options.time_limit - started.elapsed().as_secs_f64() <= 0.01 {
                return Ok(best_trial.unwrap_or(result));
            }
            let mut columns: Vec<_> = result
                .unit_flow
                .iter()
                .flatten()
                .copied()
                .chain(result.tunnel_flow.iter().flatten().copied())
                .chain(
                    self.routing
                        .outlets
                        .iter()
                        .flat_map(|&e| result.river_withdrawal[e].iter().copied()),
                )
                .collect();
            self.routing.complete_columns(&mut columns);
            let refined = points.with_level_crossings(self, &result.storage, &columns);
            if refined.len() <= points.len() {
                return Ok(best_trial.unwrap_or(result));
            }
            points = refined;
            warm = Some(crate::commitment::Proposal::from_schedule(self, &result)?);
            last = Some(result);
        }
    }

    fn dispatch_on_points(
        &self,
        controls: crate::commitment::DiscreteControls,
        start: Option<&crate::commitment::Proposal>,
        points: &crate::operating_points::OperatingPoints,
        options: &SolveOptions,
        backend: Backend,
        objective: DispatchObjective,
    ) -> Result<Schedule, Error> {
        let (case, routing) = (self.case, &self.routing);
        let dispatch = build_dispatch(self, &controls, start, points, objective)?;
        let termination = match objective {
            DispatchObjective::Revenue => nlp::Termination::Optimality,
            DispatchObjective::FlowDistance => nlp::Termination::Feasibility,
        };
        let point = nlp::solve(dispatch.problem, options, backend, termination)?;
        let storage = dispatch
            .storage
            .iter()
            .enumerate()
            .map(|(r, row)| {
                let bounds = case.system.reservoirs[r].storage;
                row.iter()
                    .map(|i| bounds.min + (bounds.max - bounds.min) * point.x[i.index()])
                    .collect()
            })
            .collect();
        // Project only solver-sized bound roundoff in the physical controls.
        // Independent reconstruction then checks all equations and restrictions.
        let roundoff_to_bounds = |x: f64, bounds: crate::Bounds, scale: f64| {
            if x >= bounds.min - options.tolerance * scale
                && x <= bounds.max + options.tolerance * scale
            {
                x.clamp(bounds.min, bounds.max)
            } else {
                x
            }
        };
        let unit_flow = dispatch
            .flow
            .iter()
            .enumerate()
            .map(|(g, row)| {
                row.iter()
                    .enumerate()
                    .map(|(t, i)| {
                        if controls.unit_on[g][self.mesh.interval[t]] {
                            roundoff_to_bounds(
                                dispatch.unit_flow_scale[g] * point.x[i.index()],
                                self.operating.units[g][self.mesh.interval[t]].flow,
                                dispatch.unit_flow_scale[g],
                            )
                        } else {
                            0.
                        }
                    })
                    .collect()
            })
            .collect();
        let tunnel_flow = dispatch
            .tunnel_flow
            .iter()
            .enumerate()
            .map(|(e, row)| {
                row.iter()
                    .map(|i| dispatch.tunnel_flow_scale[e] * point.x[i.index()])
                    .collect()
            })
            .collect();
        let mut schedule = Schedule {
            grid: self.mesh.grid.clone(),
            routing_grid: self.routing.mixing_grid.clone(),
            operating_grid: case.grid.clone(),
            storage,
            unit_flow,
            tunnel_flow,
            tunnel_opening: self.mesh.hold(
                &self
                    .operating
                    .tunnels
                    .iter()
                    .enumerate()
                    .map(|(e, row)| {
                        row.iter()
                            .enumerate()
                            .map(|(t, p)| p.positions[controls.gate_position[e][t]].opening)
                            .collect()
                    })
                    .collect::<Vec<_>>(),
            )?,
            gate_adjustment_cost: 0.,
            river_withdrawal: vec![vec![0.; self.mesh.dt.len()]; case.system.rivers.len()],
            river_opening: vec![vec![0.; self.mesh.dt.len()]; case.system.rivers.len()],
            observed_flow: Vec::new(),
            observed_volume: Vec::new(),
            river_release: Vec::new(),
            river_arrival_volume: Vec::new(),
            river_storage: Vec::new(),
            unit_power: Vec::new(),
            unit_head: Vec::new(),
            plant_flow: Vec::new(),
            plant_power: Vec::new(),
            head: Vec::new(),
            commitment: self.mesh.hold(&controls.unit_on)?,
            soft_slacks: vec![],
            soft_penalty_cost: 0.,
            release_shortfall_volume: vec![],
            release_penalty_cost: 0.,
            objective: 0.,
            solver_status: Some(point.status),
            iterations: point.iterations,
            validation: Validation::default(),
            time_error: None,
        };
        for (k, &river) in routing.outlets.iter().enumerate() {
            schedule.river_opening[river] = dispatch.river_opening[k]
                .iter()
                .enumerate()
                .map(|(t, i)| {
                    roundoff_to_bounds(
                        point.x[i.index()],
                        self.operating.river_gate[river][self.mesh.interval[t]],
                        1.,
                    )
                })
                .collect();
            schedule.river_withdrawal[river] = dispatch.river_withdrawal[k]
                .iter()
                .map(|i| dispatch.river_flow_scale[river] * point.x[i.index()])
                .collect();
        }
        let junction_head = dispatch
            .junction_head
            .iter()
            .map(|row| row.iter().map(|i| 100. * point.x[i.index()]).collect())
            .collect();
        let mut schedule = crate::trajectory::evaluate(
            self,
            schedule,
            junction_head,
            options.tolerance.max(1e-7),
        )?;
        match self.reconstruct(&schedule.controls()?) {
            Ok(mut physical) => {
                physical.solver_status = schedule.solver_status;
                physical.iterations = schedule.iterations;
                Ok(physical)
            }
            Err(Error::Solver(_) | Error::Input(_)) => {
                schedule.validation.feasible = false;
                Ok(schedule)
            }
            Err(error) => Err(error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Boundary, Bounds, Efficiency, HydroSystem, Plant, Reservoir, StorageCurve, Unit};
    use std::collections::BTreeMap;

    fn solve(case: &ScheduleCase, options: &SolveOptions) -> Result<Schedule, Error> {
        SchedulingProblem::new(case)?.solve(options)
    }
    fn solve_with_backend(
        case: &ScheduleCase,
        options: &SolveOptions,
        backend: Backend,
    ) -> Result<Schedule, Error> {
        SchedulingProblem::new(case)?.solve_with_backend(options, backend)
    }

    fn reservoir_case() -> ScheduleCase {
        ScheduleCase {
            soft_constraints: vec![],
            flow_requirements: vec![],
            routing_grid: vec![],
            operations: Default::default(),
            system: HydroSystem {
                reservoirs: vec![Reservoir {
                    operation: Default::default(),
                    name: "Lake".into(),
                    storage: Bounds::new(1e6, 4e6),
                    initial_storage: 3e6,
                    level: StorageCurve::Polynomial {
                        reference_volume: 0.,
                        elevation: 40.,
                        slope: 20e-6,
                        curvature: 0.,
                    },
                    water_value: 0.001,
                }],
                boundaries: vec![Boundary {
                    name: "Sea".into(),
                    head: 20.,
                }],
                plants: vec![Plant::new("Station", "Lake", "Sea")],
                units: vec![Unit::new(
                    "Unit",
                    "Station",
                    Bounds::new(5., 20.),
                    Bounds::new(1., 25.),
                    Bounds::new(35., 105.),
                    Efficiency::constant(0.9),
                )],
                ..Default::default()
            },
            grid: vec![0., 1.],
            prices: vec![100.],
            inflows: BTreeMap::new(),
            commitment: BTreeMap::from([("Unit".into(), vec![true])]),
        }
    }

    #[test]
    fn refined_hydraulics_preserve_original_setpoints_ramps_targets_and_energy() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 2.];
        case.prices = vec![80., 160.];
        case.commitment.insert("Unit".into(), vec![true; 2]);
        case.system.units[0].operation.initial_on = true;
        case.system.units[0].operation.initial_flow = Some(5.);
        case.system.units[0].operation.flow_ramp.up = Some(10.);
        case.system.units[0].operation.startup_cost = 100.;
        case.operations.push(crate::OperationalSeries {
            object: crate::Object::Unit("Unit".into()),
            attribute: crate::Attribute::Discharge,
            times: vec![0., 0.5],
            values: vec![5., 15.],
        });
        case.soft_constraints.push(crate::SoftConstraint {
            name: "Final storage".into(),
            target: crate::SoftTarget::Storage {
                reservoir: "Lake".into(),
                time: 2.,
            },
            lower: Some(crate::PenaltyLimit {
                value: 3e6,
                allowance: 1e5,
                cost: 0.001,
            }),
            upper: None,
        });
        case.soft_constraints.push(crate::SoftConstraint {
            name: "Mean power".into(),
            target: crate::SoftTarget::UnitPower {
                unit: "Unit".into(),
                start_time: 0.,
                end_time: 2.,
            },
            lower: Some(crate::PenaltyLimit {
                value: 10.,
                allowance: 10.,
                cost: 2.,
            }),
            upper: None,
        });
        let controls = crate::Controls {
            head_guess: None,
            commitment: vec![vec![true; 2]],
            unit_flow: vec![vec![5., 15.]],
            river_opening: vec![],
            tunnel_opening: vec![],
        };
        let p = SchedulingProblem::new(&case).unwrap();
        let coarse = p.reconstruct(&controls).unwrap();
        let fine = p
            .on_mesh(&[0., 0.125, 0.5, 0.75, 1.5, 2.], &p.routing.mixing_grid)
            .unwrap()
            .reconstruct(&controls)
            .unwrap();
        assert!(coarse.validation.feasible, "{:?}", coarse.validation);
        assert!(fine.validation.feasible, "{:?}", fine.validation);
        assert!((coarse.objective - fine.objective).abs() < 1e-7);
        assert!((fine.storage[0].last().unwrap() - (3e6 - 90000.)).abs() < 1e-7);
        for (a, b) in coarse
            .net_power(&case)
            .unwrap()
            .iter()
            .zip(fine.net_power(&case).unwrap())
        {
            assert!((a - b).abs() < 1e-9);
        }
        for (a, b) in coarse.soft_slacks.iter().zip(&fine.soft_slacks) {
            assert!(
                (a.lower - b.lower).abs() < 1e-8 && (a.penalty_cost - b.penalty_cost).abs() < 1e-8
            );
        }
        assert_eq!(fine.controls().unwrap().unit_flow, controls.unit_flow);
        case.system.units[0].operation.flow_ramp.up = Some(9.);
        let p = SchedulingProblem::new(&case).unwrap();
        let invalid = p
            .on_mesh(&fine.grid, &p.routing.mixing_grid)
            .unwrap()
            .reconstruct(&controls)
            .unwrap();
        assert!(!invalid.validation.feasible);
    }

    #[test]
    fn refined_dispatch_holds_discharge_and_charges_startup_once() {
        let mut case = reservoir_case();
        case.grid = vec![0., 1., 2.];
        case.prices = vec![50., 100.];
        case.commitment.insert("Unit".into(), vec![false, true]);
        case.system.units[0].startup_gross_head_min = Some(79.9);
        case.system.units[0].operation.startup_cost = 123.;
        // Running flow is ten m³/s. Head decreases from 80 to 79.28 m,
        // so the original interval's requested power is its integrated mean.
        let expected_power = 0.00981 * 10. * 0.9 * (80. - 0.72 / 2.);
        case.operations.push(crate::OperationalSeries {
            object: crate::Object::Unit("Unit".into()),
            attribute: crate::Attribute::Power,
            times: vec![0., 1.],
            values: vec![0., expected_power],
        });
        let controls = crate::Controls {
            head_guess: None,
            commitment: vec![vec![false, true]],
            unit_flow: vec![vec![0., 10.]],
            river_opening: vec![],
            tunnel_opening: vec![],
        };
        let p = SchedulingProblem::new(&case).unwrap();
        let fine = p
            .on_mesh(&[0., 0.3, 1., 1.2, 1.5, 2.], &p.routing.mixing_grid)
            .unwrap();
        let result = fine
            .refine(
                &controls,
                &SolveOptions {
                    time_limit: 3.,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        for &q in &result.unit_flow[0][2..] {
            assert!((q - 10.).abs() < 1e-5);
        }
        let power = result.net_power(&case).unwrap();
        assert!(power[0].abs() < 1e-8 && (power[1] - expected_power).abs() < 1e-6);
        assert!((result.unit_power[0][2] - result.unit_power[0][4]).abs() > 0.01);
        let expected = expected_power * 100. - 123. - 36.;
        assert!(
            (result.objective - expected).abs() < 1e-4,
            "{} vs {expected}",
            result.objective
        );
        let replay = fine.reconstruct(&result.controls().unwrap()).unwrap();
        assert!(replay.validation.feasible, "{:?}", replay.validation);
        assert!((replay.objective - result.objective).abs() < 1e-7);
    }

    #[test]
    fn fixed_speed_pump_flows_follow_each_hydraulic_step_without_new_controls() {
        let mut case = reservoir_case();
        case.system.plants[0].source = "Sea".into();
        case.system.plants[0].target = "Lake".into();
        let pump = &mut case.system.units[0];
        pump.mode = crate::UnitMode::Pump;
        pump.head = Bounds::new(60., 100.);
        pump.flow = Bounds::new(8., 12.);
        pump.pump_characteristic = Some(
            crate::PumpCharacteristic::new(
                vec![60., 80., 100.],
                vec![11., 10., 9.],
                vec![0.85; 3],
                crate::Extrapolation::Error,
            )
            .unwrap(),
        );
        let controls = crate::Controls {
            head_guess: None,
            commitment: vec![vec![true]],
            unit_flow: vec![vec![10.]],
            river_opening: vec![],
            tunnel_opening: vec![],
        };
        let p = SchedulingProblem::new(&case).unwrap();
        let fine = p
            .on_mesh(&[0., 0.1, 0.4, 1.], &p.routing.mixing_grid)
            .unwrap();
        let result = fine
            .refine(
                &controls,
                &SolveOptions {
                    time_limit: 3.,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        let mut storage = 3e6;
        for (t, w) in result.grid.windows(2).enumerate() {
            let h0 = 20. + 20e-6 * storage;
            let q = (14. - 0.05 * h0) / (1. + 0.05 * 20e-6 * 3600. * (w[1] - w[0]) / 2.);
            assert!((result.unit_flow[0][t] - q).abs() < 1e-6);
            storage += q * 3600. * (w[1] - w[0]);
            assert!((result.storage[0][t + 1] - storage).abs() < 1e-5);
        }
        assert!((result.unit_flow[0][0] - result.unit_flow[0][2]).abs() > 0.01);
        assert!(result.objective < 0.);
    }

    #[test]
    fn refined_dispatch_restores_feasibility_in_a_small_reservoir_with_delayed_arrivals() {
        let mut case = delayed_forebay_case();
        case.system.reservoirs[1].level = StorageCurve::Table(
            crate::TableCurve::new(vec![0., 400., 600., 1000.], vec![45., 48., 51., 55.]).unwrap(),
        );
        let controls = crate::Controls {
            head_guess: None,
            commitment: vec![vec![true; 2]; 2],
            unit_flow: vec![vec![20.; 2], vec![10.; 2]],
            river_opening: vec![vec![0.; 2]],
            tunnel_opening: vec![],
        };
        let p = SchedulingProblem::new(&case).unwrap();
        let fine = p
            .on_mesh(&[0., 0.1, 0.5, 0.75, 1., 1.25, 1.5], &p.routing.mixing_grid)
            .unwrap();
        let initial = fine.reconstruct(&controls).unwrap();
        assert!(!initial.validation.feasible);
        let result = fine
            .refine(
                &controls,
                &SolveOptions {
                    time_limit: 4.,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        assert!(
            result.storage[1]
                .iter()
                .all(|v| (-1e-4..=1000.0001).contains(v))
        );
        let held = fine
            .mesh
            .hold(&result.controls().unwrap().unit_flow)
            .unwrap();
        for (row, expected) in result.unit_flow.iter().zip(held) {
            for (q, value) in row.iter().zip(expected) {
                assert!((q - value).abs() < 1e-9);
            }
        }
        let replay = fine.reconstruct(&result.controls().unwrap()).unwrap();
        assert!(replay.validation.feasible, "{:?}", replay.validation);
        assert!((replay.objective - result.objective).abs() < 1e-6);
    }

    #[test]
    fn moving_storage_knot_constraints_limit_actual_peak_generation() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.05];
        case.prices = vec![100.];
        let upper = &mut case.system.reservoirs[0];
        upper.storage = Bounds::new(0., 1000.);
        upper.initial_storage = 1000.;
        upper.water_value = 0.;
        upper.level = StorageCurve::Table(
            crate::TableCurve::new(vec![0., 600., 1000.], vec![30., 90., 92.]).unwrap(),
        );
        case.system.reservoirs.push(Reservoir {
            name: "Lower lake".into(),
            storage: Bounds::new(0., 1000.),
            initial_storage: 1000.,
            water_value: 0.,
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 20.,
                slope: 0.03,
                curvature: 0.,
            },
            operation: Default::default(),
        });
        case.system.plants[0].target = "Lower lake".into();
        case.system
            .plants
            .push(Plant::new("Lower station", "Lower lake", "Sea"));
        let unit = &mut case.system.units[0];
        unit.flow = Bounds::new(2., 5.);
        unit.power = Bounds::new(0., 1.6);
        unit.head = Bounds::new(1., 100.);
        case.system.units.push(Unit::new(
            "Lower unit",
            "Lower station",
            Bounds::new(6., 6.),
            Bounds::new(0., 10.),
            Bounds::new(1., 100.),
            Efficiency::constant(0.9),
        ));
        case.commitment.insert("Lower unit".into(), vec![true]);
        let mut controls = crate::Controls {
            head_guess: None,
            commitment: vec![vec![true]; 2],
            unit_flow: vec![vec![4.2], vec![6.]],
            river_opening: vec![],
            tunnel_opening: vec![],
        };
        let p = SchedulingProblem::new(&case).unwrap();
        // Both interval endpoints obey the upper machine's power limit, but
        // its working head peaks when the upper storage crosses 600 m³.
        let physical = p.reconstruct(&controls).unwrap();
        for k in [0, physical.grid.len() - 1] {
            let head = case.system.reservoirs[0]
                .level
                .level(physical.storage[0][k])
                .unwrap()
                - case.system.reservoirs[1]
                    .level
                    .level(physical.storage[1][k])
                    .unwrap();
            assert!(
                case.system.units[0]
                    .power(controls.unit_flow[0][0], head)
                    .unwrap()
                    < case.system.units[0].power.max
            );
        }
        assert!(!physical.validation.feasible);
        controls.unit_flow[0][0] = 3.;
        let result = p
            .refine(
                &controls,
                &SolveOptions {
                    time_limit: 3.,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        let q = result.controls().unwrap().unit_flow[0][0];
        // At the moving crossing, t=400/Q seconds and Q*H=28Q+72.
        let optimum = (1.6 / (0.00981 * 0.9) - 72.) / 28.;
        let peak_power = 0.00981 * 0.9 * (28. * q + 72.);
        assert!(peak_power <= 1.6 + 1e-6, "{peak_power}");
        assert!((q - optimum).abs() < 1e-3, "Q={q}, expected {optimum}");
        let (fine, error) = p
            .validate_replay(&result, &crate::TimeAccuracy::default())
            .unwrap();
        assert!(fine.validation.feasible && error.converged());
    }

    #[test]
    fn small_reservoir_outlet_converges_to_analytic_decay_with_held_gate() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5];
        case.prices = vec![0.];
        case.commitment.clear();
        case.system.units.clear();
        case.system.plants.clear();
        let r = &mut case.system.reservoirs[0];
        r.storage = Bounds::new(0., 1000.);
        r.initial_storage = 500.;
        r.water_value = 0.;
        r.level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 10.,
            slope: 0.01,
            curvature: 0.,
        };
        case.system.rivers.push(crate::River {
            name: "Outlet".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 5.,
            law: crate::RiverLaw::Level {
                curve: crate::TableCurve::new(vec![10., 20.], vec![0., 1.]).unwrap(),
            },
            gate: Bounds::new(1., 1.),
            operation: Default::default(),
            transport: Default::default(),
            allocation: 1.,
            inflow: 0.,
            water_value: 0.,
            history_grid: vec![],
            history_release: vec![],
        });
        let controls = crate::Controls {
            head_guess: None,
            commitment: vec![],
            unit_flow: vec![],
            river_opening: vec![vec![1.]],
            tunnel_opening: vec![],
        };
        let p = SchedulingProblem::new(&case).unwrap();
        let exact = 500. * (-1.8f64).exp();
        let mut previous_error = f64::INFINITY;
        for steps in [1, 2, 4, 8, 16] {
            let grid: Vec<_> = (0..=steps).map(|k| 0.5 * k as f64 / steps as f64).collect();
            let result = p
                .on_mesh(&grid, &p.routing.mixing_grid)
                .unwrap()
                .reconstruct(&controls)
                .unwrap();
            assert!(result.validation.feasible, "{:?}", result.validation);
            let (replay, estimate) = p
                .validate_replay(&result, &crate::TimeAccuracy::default())
                .unwrap();
            assert!(replay.validation.feasible, "{:?}", replay.validation);
            assert_eq!(
                replay.controls().unwrap().river_opening,
                controls.river_opening
            );
            if steps == 1 {
                assert!(!estimate.converged());
            }
            if steps == 16 {
                assert!(estimate.converged(), "{estimate:?}");
            }
            let end = *result.storage[0].last().unwrap();
            let error = (end - exact).abs();
            assert!(error < previous_error);
            if steps >= 4 {
                assert!(previous_error / error > 3.5);
            }
            previous_error = error;
            let released: f64 = result.river_release[0]
                .iter()
                .zip(grid.windows(2))
                .map(|(q, w)| 3600. * q * (w[1] - w[0]))
                .sum();
            assert!((released + end - 500.).abs() < 1e-6);
            assert_eq!(
                result.controls().unwrap().river_opening,
                controls.river_opening
            );
        }
        assert!(previous_error < 0.2);
        let optimized = p
            .refine(
                &controls,
                &SolveOptions {
                    time_limit: 2.,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(optimized.validation.feasible, "{:?}", optimized.validation);
        assert!(optimized.time_error.as_ref().unwrap().converged());
        assert!(
            (optimized.storage[0].last().unwrap() - exact).abs() < 1.,
            "end {:?}, exact {exact}, error {:?}",
            optimized.storage[0].last(),
            optimized.time_error
        );
        let (checked, error) = p
            .validate_replay(&optimized, &crate::TimeAccuracy::default())
            .unwrap();
        assert!(checked.validation.feasible && error.converged());
        // The one-hour midpoint equation has negative final storage although
        // the continuous outlet law preserves positivity. Grid error must be
        // resolved before declaring the physical controls infeasible.
        case.grid[1] = 1.;
        let p = SchedulingProblem::new(&case).unwrap();
        assert!(!p.reconstruct(&controls).unwrap().validation.feasible);
        let result = p
            .refine(
                &controls,
                &SolveOptions {
                    time_limit: 3.,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        assert!(result.time_error.as_ref().unwrap().converged());
        let exact = 500. * (-3.6f64).exp();
        assert!((result.storage[0].last().unwrap() - exact).abs() < 1.);
        assert!(result.storage[0].iter().all(|&v| v > 0.));
    }

    #[test]
    fn prepared_problem_can_be_solved_repeatedly_without_changing_inputs() {
        let case = reservoir_case();
        let problem = SchedulingProblem::new(&case).unwrap();
        let first = problem.solve(&SolveOptions::default()).unwrap();
        let second = problem.solve(&SolveOptions::default()).unwrap();
        assert!(first.validation.feasible && second.validation.feasible);
        assert!((first.objective - second.objective).abs() < 1e-6);
        assert_eq!(problem.case().system.reservoirs[0].initial_storage, 3e6);
    }

    #[test]
    fn preparation_rejects_invalid_physical_data_before_solving() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].initial_storage = -1.;
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
    }

    #[test]
    fn elapsed_time_and_discharge_ramps_are_enforced_on_unequal_intervals() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 1.5, 3.];
        case.prices = vec![-100., 200., 200.];
        case.commitment.clear();
        case.system.units[0].operation = crate::UnitOperation {
            initial_on: true,
            initial_duration: 0.,
            minimum_up_time: 1.,
            initial_flow: Some(5.),
            initial_interval_hours: 0.5,
            flow_ramp: crate::RampRates {
                up: Some(2.),
                down: Some(3.),
            },
            ..Default::default()
        };
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(
                &case,
                &SolveOptions {
                    time_limit: 5.,
                    ..Default::default()
                },
                backend,
            )
            .unwrap();
            assert!(result.validation.feasible, "{:?}", result.validation);
            assert_eq!(result.commitment[0], vec![true, true, true]);
            for (&q, expected) in result.unit_flow[0].iter().zip([6., 7.5, 10.]) {
                assert!(
                    (q - expected).abs() < 1e-4,
                    "{backend:?}: {:?}",
                    result.unit_flow
                );
            }
        }
    }

    #[test]
    fn power_ramps_use_nonlinear_generation_and_unequal_interval_midpoints() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 1.5, 3.];
        case.prices = vec![200.; 3];
        case.commitment.clear();
        case.system.units[0].operation = crate::UnitOperation {
            initial_on: true,
            minimum_up_time: 3.,
            initial_power: Some(5.),
            initial_interval_hours: 0.5,
            power_ramp: crate::RampRates {
                up: Some(2.),
                down: Some(3.),
            },
            ..Default::default()
        };
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(
                &case,
                &SolveOptions {
                    time_limit: 5.,
                    ..Default::default()
                },
                backend,
            )
            .unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert_eq!(result.commitment[0], vec![true; 3]);
            for (&power, expected) in result.unit_power[0].iter().zip([6., 7.5, 10.]) {
                assert!(
                    (power - expected).abs() < 1e-4,
                    "{backend:?}: {:?}",
                    result.unit_power
                );
            }
        }
    }

    #[test]
    fn minimum_power_is_allowed_at_startup_and_shutdown_with_zero_normal_ramp() {
        let mut case = reservoir_case();
        case.grid = vec![0., 1., 2.];
        case.prices = vec![200., -100.];
        case.commitment.insert("Unit".into(), vec![true, false]);
        case.system.reservoirs[0].level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        let unit = &mut case.system.units[0];
        unit.flow.min = 1.;
        unit.power.min = 0.00981 * 100. * 0.9;
        unit.operation.power_ramp = crate::RampRates {
            up: Some(0.),
            down: Some(0.),
        };
        let minimum_power = unit.power.min;
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.unit_power[0][0] - minimum_power).abs() < 1e-6);
            assert!((result.unit_flow[0][0] - 1.).abs() < 1e-6);
            assert_eq!(result.unit_power[0][1], 0.);
        }
    }

    #[test]
    fn historical_output_can_exceed_new_limits_but_transition_must_obey_ramps() {
        let mut case = reservoir_case();
        let unit = &mut case.system.units[0];
        unit.power.max = 8.;
        unit.operation = crate::UnitOperation {
            initial_on: true,
            initial_flow: Some(25.),
            initial_power: Some(15.),
            flow_ramp: crate::RampRates {
                down: Some(20.),
                ..Default::default()
            },
            power_ramp: crate::RampRates {
                // The instantaneous maximum is eight MW; declining head makes
                // the feasible interval average slightly smaller than eight.
                down: Some(7.1),
                ..Default::default()
            },
            ..Default::default()
        };
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        assert!((7.9..8.).contains(&result.unit_power[0][0]));
        let initial_power = case.system.units[0]
            .power(result.unit_flow[0][0], 80.)
            .unwrap();
        assert!((initial_power - 8.).abs() < 1e-5);
        case.system.units[0].operation.power_ramp.down = Some(6.);
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(!result.validation.feasible);
    }

    #[test]
    fn a_stopped_unit_cannot_supply_nonzero_historical_generation() {
        let mut case = reservoir_case();
        case.system.units[0].operation.initial_power = Some(1.);
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
    }

    #[test]
    fn a_start_requires_remaining_on_through_the_minimum_up_time() {
        let mut case = reservoir_case();
        case.grid = vec![0., 1., 2., 3., 4.];
        case.prices = vec![-100., 200., -100., 200.];
        case.commitment.clear();
        case.system.units[0].operation.minimum_up_time = 2.;
        let result = solve(
            &case,
            &SolveOptions {
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.commitment[0], vec![false, true, true, true]);
        assert!((result.unit_flow[0][2] - 5.).abs() < 1e-4);
    }

    #[test]
    fn initial_down_time_prevents_an_early_start() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 1.5, 3.];
        case.prices = vec![200.; 3];
        case.commitment.clear();
        case.system.units[0].operation.minimum_down_time = 1.5;
        case.system.units[0].operation.initial_duration = 0.5;
        let result = solve(
            &case,
            &SolveOptions {
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.commitment[0], vec![false, false, true]);
    }

    #[test]
    fn shutdown_is_charged_once_and_prescribed_dwell_violations_are_rejected() {
        let mut case = reservoir_case();
        case.commitment.get_mut("Unit").unwrap()[0] = false;
        case.system.units[0].operation.initial_on = true;
        case.system.units[0].operation.shutdown_cost = 137.;
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert!((result.objective - (-137.)).abs() < 1e-6);
        case.system.units[0].operation.minimum_up_time = 2.;
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
    }

    #[test]
    fn measured_storage_curve_optimizes_across_a_knot_with_exact_energy() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level = StorageCurve::Table(
            crate::TableCurve::new(vec![1e6, 2_950_000., 3e6, 4e6], vec![60., 99., 102., 110.])
                .unwrap(),
        );
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.unit_flow[0][0] - 20.).abs() < 1e-5);
            let mean = (22000. * (98.56 + 99.) / 2. + 50000. * (99. + 102.) / 2.) / 72000.;
            let expected = 100. * 0.00981 * 20. * 0.9 * (mean - 20.) - 72.;
            assert!((result.head[0][0] - mean).abs() < 1e-6);
            assert!(
                (result.objective - expected).abs() < 0.01,
                "{backend:?}: {} vs {expected}",
                result.objective
            );
            let after = result.storage[0][1];
            let curve = &case.system.reservoirs[0].level;
            let potential_change = curve.potential(3e6).unwrap() - curve.potential(after).unwrap();
            assert!((potential_change - (3e6 - after) * result.head[0][0]).abs() < 1e-5);
        }
    }

    #[test]
    fn measured_storage_curve_handles_no_volume_change_at_a_knot() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level = StorageCurve::Table(
            crate::TableCurve::new(vec![1e6, 3e6, 4e6], vec![60., 100., 105.]).unwrap(),
        );
        case.commitment.get_mut("Unit").unwrap()[0] = false;
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.storage[0][1] - 3e6).abs() < 0.001);
            assert!((result.head[0][0] - 100.).abs() < 1e-8);
        }
    }

    #[test]
    fn measured_forebay_curve_integrates_the_complete_delayed_storage_path() {
        let mut case = delayed_forebay_case();
        case.system.reservoirs[1].level = StorageCurve::Table(
            crate::TableCurve::new(vec![0., 400., 600., 1000.], vec![45., 48., 51., 55.]).unwrap(),
        );
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            let middle =
                result.storage[1][1] + 1800. * (result.unit_flow[0][0] - result.unit_flow[1][1]);
            assert!((-0.01..=1000.01).contains(&middle));
            let curve = &case.system.reservoirs[1].level;
            let mean = |a: f64, b: f64| {
                if a == b {
                    return curve.trial_level(a);
                }
                let (lo, hi) = (a.min(b), a.max(b));
                let mut knots = vec![lo];
                knots.extend([400., 600.].into_iter().filter(|&v| lo < v && v < hi));
                knots.push(hi);
                knots
                    .windows(2)
                    .map(|w| {
                        (w[1] - w[0]) * (curve.trial_level(w[0]) + curve.trial_level(w[1])) / 2.
                    })
                    .sum::<f64>()
                    / (hi - lo)
            };
            let expected =
                0.5 * (mean(result.storage[1][1], middle) + mean(middle, result.storage[1][2]));
            assert!((result.head[1][1] - expected).abs() < 1e-7);
        }
    }

    #[test]
    fn storage_tables_must_be_monotone_and_cover_the_operating_domain() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level =
            StorageCurve::Table(crate::TableCurve::new(vec![2e6, 4e6], vec![60., 100.]).unwrap());
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
        case.system.reservoirs[0].level =
            StorageCurve::Table(crate::TableCurve::new(vec![1e6, 4e6], vec![100., 60.]).unwrap());
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
    }

    #[test]
    fn measured_watercourse_optimizes_with_pressure_loops_river_splits_and_commitment() {
        let case: ScheduleCase =
            serde_json::from_str(include_str!("../examples/data/tabulated_cascade.json")).unwrap();
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(
                &case,
                &SolveOptions {
                    time_limit: 5.,
                    ..Default::default()
                },
                backend,
            )
            .unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            // Count every parcel once: in a reservoir or in transit through a reach.
            for t in 0..case.prices.len() {
                let seconds = 3600. * (case.grid[t + 1] - case.grid[t]);
                let stored: f64 = result.storage.iter().map(|row| row[t + 1] - row[t]).sum();
                let in_transit: f64 = result
                    .river_storage
                    .iter()
                    .map(|row| row[t + 1] - row[t])
                    .sum();
                let natural: f64 = case.inflows.values().map(|row| row[t]).sum();
                let exported = result.unit_flow[2][t];
                assert!((stored + in_transit - seconds * (natural - exported)).abs() < 0.01);
            }
        }
    }

    #[test]
    fn a_table_flat_over_the_operating_range_has_the_correct_constant_head() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level = StorageCurve::Table(
            crate::TableCurve::new(vec![0., 1e6, 4e6, 5e6], vec![50., 100., 100., 110.]).unwrap(),
        );
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(result.validation.feasible);
            assert_eq!(result.head[0][0], 100.);
            assert!((result.unit_flow[0][0] - 20.).abs() < 1e-5);
            assert!((result.objective - (100. * 0.00981 * 20. * 80. * 0.9 - 72.)).abs() < 0.01);
        }
    }

    #[test]
    fn complete_reservoir_optimization() {
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let case = reservoir_case();
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!(
                (result.unit_flow[0][0] - 20.).abs() < 1e-5,
                "{backend:?}: {:?}",
                result.unit_flow
            );
            assert!((result.storage[0][1] - 2_928_000.).abs() < 0.05);
            let expected_power =
                0.00981 * 20. * 0.9 * (40. + 20e-6 * (3e6 + 2_928_000.) / 2. - 20.);
            let expected = 100. * expected_power + 0.001 * (2_928_000. - 3e6);
            assert!(
                (result.objective - expected).abs() < 0.01,
                "{backend:?}: {} vs {expected}",
                result.objective
            );
        }
    }

    #[test]
    fn off_units_have_exactly_zero_discharge_and_generation() {
        let mut case = reservoir_case();
        case.commitment.get_mut("Unit").unwrap()[0] = false;
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.unit_flow[0][0], 0.);
        assert_eq!(result.unit_power[0][0], 0.);
        assert!((result.storage[0][1] - 3e6).abs() < 1e-6);
    }

    #[test]
    fn storage_head_is_energy_consistent_even_for_equal_volumes() {
        let c = StorageCurve::Polynomial {
            reference_volume: 1e6,
            elevation: 75.,
            slope: 2e-5,
            curvature: 1e-11,
        };
        for (a, b) in [(1e6, 2e6), (2e6, 1e6), (1e6, 1e6), (1e6, 1e6 + 1e-5)] {
            let mean = c.mean_level(a, b).unwrap();
            assert!(mean.is_finite());
            if a == b {
                assert_eq!(mean, c.level(a).unwrap());
            } else {
                assert!(
                    ((b - a) * mean - (c.potential(b).unwrap() - c.potential(a).unwrap())).abs()
                        < 1e-7 * (1. + ((b - a) * mean).abs())
                );
            }
        }
    }

    #[test]
    fn commitment_can_turn_off_during_negative_prices() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 1.5];
        case.prices = vec![-100., 150.];
        case.commitment.clear();
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(
                &case,
                &SolveOptions {
                    time_limit: 5.,
                    ..Default::default()
                },
                backend,
            )
            .unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert_eq!(result.commitment[0], vec![false, true]);
            assert_eq!(result.unit_flow[0][0], 0.);
            assert!((result.unit_flow[0][1] - 20.).abs() < 1e-5);
        }
    }

    #[test]
    fn pressure_loop_conserves_water_and_handles_reversed_tunnel_orientation() {
        use crate::{Junction, Tunnel};
        let mut case = reservoir_case();
        case.system.junctions.push(Junction {
            name: "Intake".into(),
            head: Bounds::new(25., 120.),
        });
        case.system.plants[0].source = "Intake".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        case.system.tunnels = vec![
            Tunnel {
                name: "Main".into(),
                source: "Lake".into(),
                target: "Intake".into(),
                capacity: 20.,
                loss_factor: 0.1,
                opening: 1.,
                discharge_river: None,
                gate_positions: vec![],
                initial_opening: None,
                gate_adjustment_cost: 0.,
                gate_loss_factor_curve: None,
            },
            Tunnel {
                name: "Parallel".into(),
                source: "Intake".into(),
                target: "Lake".into(),
                capacity: 20.,
                loss_factor: 0.4,
                opening: 1.,
                discharge_river: None,
                gate_positions: vec![],
                initial_opening: None,
                gate_adjustment_cost: 0.,
                gate_loss_factor_curve: None,
            },
        ];
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.unit_flow[0][0] - 20.).abs() < 1e-4);
            let (a, b) = (result.tunnel_flow[0][0], result.tunnel_flow[1][0]);
            assert!(a > 0. && b < 0.);
            assert!((a / (-b) - 2.).abs() < 1e-6);
            assert!((a - b - result.unit_flow[0][0]).abs() < 1e-6);
        }
    }

    #[test]
    fn small_reservoir_couples_two_plants_without_small_time_steps() {
        let mut case = reservoir_case();
        case.system.reservoirs.push(Reservoir {
            operation: Default::default(),
            name: "Forebay".into(),
            storage: Bounds::new(0., 1000.),
            initial_storage: 500.,
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 45.,
                slope: 0.01,
                curvature: 0.,
            },
            water_value: 0.,
        });
        case.system.plants[0].target = "Forebay".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        let mut downstream = case.system.units[0].clone();
        downstream.name = "Downstream".into();
        downstream.plant = "Downstream".into();
        case.system
            .plants
            .push(Plant::new("Downstream", "Forebay", "Sea"));
        case.system.units.push(downstream);
        case.grid = vec![0., 1., 3.];
        case.prices = vec![30., 150.];
        case.commitment = BTreeMap::from([
            ("Unit".into(), vec![true, true]),
            ("Downstream".into(), vec![true, true]),
        ]);
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            for (t, dt) in [1., 2.].iter().enumerate() {
                let change = 3600. * dt * (result.unit_flow[0][t] - result.unit_flow[1][t]);
                assert!((change - (result.storage[1][t + 1] - result.storage[1][t])).abs() < 0.005);
            }
        }
    }

    fn delayed_forebay_case() -> ScheduleCase {
        let mut case = reservoir_case();
        case.system.reservoirs.push(Reservoir {
            operation: Default::default(),
            name: "Forebay".into(),
            storage: Bounds::new(0., 1000.),
            initial_storage: 500.,
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 45.,
                slope: 0.01,
                curvature: 0.,
            },
            water_value: 0.,
        });
        case.system.boundaries.push(Boundary {
            name: "Tailwater".into(),
            head: 60.,
        });
        case.system.plants[0].target = "Tailwater".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        case.system.plants[0].discharge_river = Some("River".into());
        let mut lower = case.system.units[0].clone();
        lower.name = "Lower".into();
        lower.plant = "Downstream".into();
        case.system
            .plants
            .push(Plant::new("Downstream", "Forebay", "Sea"));
        case.system.units.push(lower);
        let mut river = crate::River {
            operation: Default::default(),
            law: Default::default(),
            gate: Bounds::new(0., 1.),
            name: "River".into(),
            source: "Outfall".into(),
            target: "Forebay".into(),
            transport: crate::RiverTransport::Fixed { delay: 0.5 },
            capacity: 30.,
            allocation: 1.,
            inflow: 0.,
            water_value: 0.,
            history_grid: Vec::new(),
            history_release: Vec::new(),
        };
        river.history_grid = vec![-1., 0.];
        river.history_release = vec![5.];
        case.system.rivers.push(river);
        case.grid = vec![0., 0.5, 1.5];
        case.prices = vec![-100., 150.];
        case.commitment = BTreeMap::from([
            ("Unit".into(), vec![true, true]),
            ("Lower".into(), vec![true, true]),
        ]);
        case
    }

    #[test]
    fn delayed_outfall_optimization_respects_small_reservoir_between_control_times() {
        let case = delayed_forebay_case();
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.river_arrival_volume[0][0] - 9000.).abs() < 0.01);
            // Within the second interval the arrival changes at hour 1.0.
            let middle =
                result.storage[1][1] + 1800. * (result.unit_flow[0][0] - result.unit_flow[1][1]);
            assert!(
                (-0.01..=1000.01).contains(&middle),
                "{backend:?}: middle storage {middle}"
            );
            let mean_volume = (result.storage[1][1] + 2. * middle + result.storage[1][2]) / 4.;
            assert!((result.head[1][1] - (45. + 0.01 * mean_volume)).abs() < 1e-8);
            for t in 0..case.prices.len() {
                let seconds = 3600. * (case.grid[t + 1] - case.grid[t]);
                let change = result.river_arrival_volume[0][t] - seconds * result.unit_flow[1][t];
                assert!((result.storage[1][t + 1] - result.storage[1][t] - change).abs() < 0.01);
            }
        }
    }

    #[test]
    fn endpoint_balance_cannot_hide_overflow_between_control_times() {
        let mut case = delayed_forebay_case();
        case.grid = vec![0., 1.];
        case.prices = vec![100.];
        case.system.rivers[0].history_release = vec![20.];
        case.commitment =
            BTreeMap::from([("Unit".into(), vec![false]), ("Lower".into(), vec![true])]);
        let result = solve(
            &case,
            &SolveOptions {
                max_iterations: 100,
                ..Default::default()
            },
        )
        .unwrap();
        // Half an hour of 20 m³/s followed by zero cannot be buffered in 1000 m³
        // with a single constant downstream flow, although endpoint balance can.
        assert!(!result.validation.feasible);
    }

    #[test]
    fn commitment_search_includes_delayed_water_and_river_inventory() {
        let mut case = delayed_forebay_case();
        case.commitment.clear();
        let result = solve(
            &case,
            &SolveOptions {
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        // In-flight water forces the downstream unit on in the first interval.
        assert!(result.commitment[1][0]);
        let releases: f64 = result.unit_flow[0]
            .iter()
            .enumerate()
            .map(|(t, q)| q * 3600. * (case.grid[t + 1] - case.grid[t]))
            .sum();
        let arrivals: f64 = result.river_arrival_volume[0].iter().sum();
        assert!(
            (result.river_storage[0][0] + releases - arrivals - result.river_storage[0][2]).abs()
                < 0.01
        );
    }

    #[test]
    fn constant_turbine_tables_match_analytic_generation_with_changing_head_and_losses() {
        let mut case = reservoir_case();
        case.grid = vec![0., 0.5, 1.5, 3.];
        case.prices = vec![-100., 200., 100.];
        case.commitment.insert("Unit".into(), vec![true; 3]);
        case.system.units[0].loss_factor = 0.005;
        for backend in [Backend::Ipopt, Backend::Pounce] {
            case.system.units[0].efficiency = Efficiency::constant(0.9);
            let analytic = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            let unit = &mut case.system.units[0];
            unit.efficiency = Efficiency::Table(
                crate::TurbineTable::new(
                    vec![unit.head.min, unit.head.max],
                    vec![
                        crate::TableCurve::new(vec![5., 12., 20.], vec![0.9; 3]).unwrap(),
                        crate::TableCurve::new(vec![5., 9., 15., 20.], vec![0.9; 4]).unwrap(),
                    ],
                    vec![unit.flow; 2],
                    crate::TurbineOptions {
                        interpolation: crate::DischargeInterpolation::Pchip,
                        ..Default::default()
                    },
                )
                .unwrap(),
            );
            let tabulated = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                analytic.validation.feasible && tabulated.validation.feasible,
                "{backend:?}: {:?}",
                tabulated.validation
            );
            assert!(
                (tabulated.objective - analytic.objective).abs()
                    < 1e-8 * analytic.objective.abs().max(1.),
                "{backend:?}: analytic {}, table {}, flows {:?} {:?}, iterations {} {}, status {:?} {:?}",
                analytic.objective,
                tabulated.objective,
                analytic.unit_flow,
                tabulated.unit_flow,
                analytic.iterations,
                tabulated.iterations,
                analytic.solver_status,
                tabulated.solver_status
            );
            for (&a, &b) in analytic.unit_flow[0].iter().zip(&tabulated.unit_flow[0]) {
                assert!((a - b).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn measured_turbine_envelope_and_nonlinear_power_ramp_limit_generation() {
        let mut case = reservoir_case();
        case.grid = vec![0., 1., 2.];
        case.prices = vec![100.; 2];
        case.commitment.insert("Unit".into(), vec![true; 2]);
        case.system.reservoirs[0].level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        let unit = &mut case.system.units[0];
        unit.efficiency = Efficiency::Table(
            crate::TurbineTable::new(
                vec![80., 120.],
                vec![
                    crate::TableCurve::new(vec![5., 12., 20.], vec![0.9; 3]).unwrap(),
                    crate::TableCurve::new(vec![5., 9., 15., 20.], vec![0.9; 4]).unwrap(),
                ],
                vec![Bounds::new(5., 10.), Bounds::new(5., 20.)],
                crate::TurbineOptions::default(),
            )
            .unwrap(),
        );
        unit.head = Bounds::new(80., 120.);
        unit.operation.initial_on = true;
        unit.operation.initial_power = Some(5.);
        unit.operation.power_ramp.up = Some(3.);
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!(
                (result.unit_power[0][0] - 8.).abs() < 1e-5,
                "{backend:?}: {:?}, status {:?}, iterations {}",
                result.unit_power,
                result.solver_status,
                result.iterations
            );
            assert!((result.unit_power[0][1] - 11.).abs() < 1e-5);
        }
        case.system.units[0].operation.power_ramp.up = None;
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        // At 100 m the linearly interpolated upper discharge is 15 m³/s.
        assert!((result.unit_flow[0][0] - 15.).abs() < 1e-5);
        assert!((result.unit_power[0][0] - 0.00981 * 15. * 100. * 0.9).abs() < 1e-5);
    }

    #[test]
    fn electrical_efficiency_balance_recovers_when_a_trial_shaft_power_has_no_root() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        case.system.units[0].power = Bounds::new(1., 10.);
        case.system.units[0].generator_efficiency_curve =
            Some(crate::TableCurve::new(vec![0., 10.], vec![0.2, 0.9]).unwrap());
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(
                result.validation.feasible,
                "{backend:?}: {:?}",
                result.validation
            );
            assert!((result.unit_power[0][0] - 10.).abs() < 1e-5);
            assert!((result.unit_flow[0][0] - 10. / (0.00981 * 100. * 0.9 * 0.9)).abs() < 1e-5);
        }
    }

    #[test]
    fn nonlinear_efficiency_can_make_partial_load_optimal() {
        let mut case = reservoir_case();
        case.system.reservoirs[0].level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        case.system.reservoirs[0].water_value = 0.;
        let unit = &mut case.system.units[0];
        unit.flow = Bounds::new(5., 15.);
        unit.head = Bounds::new(90., 110.);
        unit.efficiency = Efficiency::Quadratic {
            peak: 0.9,
            best_flow: 10.,
            flow_curvature: 0.01,
            best_head: 100.,
            head_curvature: 0.,
        };
        let expected = (0.4 + 0.148_f64.sqrt()) / 0.06;
        for backend in [Backend::Ipopt, Backend::Pounce] {
            let result = solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
            assert!(result.validation.feasible);
            assert!(
                (result.unit_flow[0][0] - expected).abs() < 1e-5,
                "{backend:?}: {:?}",
                result.unit_flow
            );
        }
    }

    #[test]
    fn partial_tunnel_openings_and_gate_resistance_match_hand_calculated_flows() {
        let mut case = reservoir_case();
        case.commitment.insert("Unit".into(), vec![false]);
        case.system.reservoirs[0].level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        case.system.tunnels = vec![crate::Tunnel {
            name: "Gate".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 50.,
            loss_factor: 1.,
            opening: 0.25,
            discharge_river: None,
            gate_positions: vec![],
            initial_opening: None,
            gate_adjustment_cost: 0.,
            gate_loss_factor_curve: None,
        }];
        for backend in [Backend::Ipopt, Backend::Pounce] {
            for tabulated in [false, true] {
                case.system.tunnels[0].gate_loss_factor_curve = if tabulated {
                    Some(crate::TableCurve::new(vec![0., 0.25, 1.], vec![30., 15., 0.]).unwrap())
                } else {
                    None
                };
                for reverse in [false, true] {
                    let tunnel = &mut case.system.tunnels[0];
                    tunnel.source = if reverse { "Sea" } else { "Lake" }.into();
                    tunnel.target = if reverse { "Lake" } else { "Sea" }.into();
                    let result =
                        solve_with_backend(&case, &SolveOptions::default(), backend).unwrap();
                    assert!(
                        result.validation.feasible,
                        "{backend:?}: {:?}",
                        result.validation
                    );
                    // No table: .25 * 100 = Q². Table: 100 = (1 + 15) Q².
                    let magnitude = if tabulated { 2.5 } else { 5. };
                    let signed = if reverse { -magnitude } else { magnitude };
                    assert!((result.tunnel_flow[0][0] - signed).abs() < 1e-6);
                    assert!((result.storage[0][1] - (3e6 - 3600. * magnitude)).abs() < 0.01);
                }
            }
        }
    }

    #[test]
    fn a_small_opening_does_not_hide_pressure_errors_in_validation() {
        let mut case = reservoir_case();
        case.commitment.insert("Unit".into(), vec![false]);
        case.system.reservoirs[0].level = StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: 120.,
            slope: 0.,
            curvature: 0.,
        };
        case.system.tunnels = vec![crate::Tunnel {
            name: "Gate".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 50.,
            loss_factor: 1.,
            opening: 1e-4,
            discharge_river: None,
            gate_positions: vec![],
            initial_opening: None,
            gate_adjustment_cost: 0.,
            gate_loss_factor_curve: None,
        }];
        let problem = SchedulingProblem::new(&case).unwrap();
        let mut result = problem.solve(&SolveOptions::default()).unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        result.tunnel_flow[0][0] += 1e-5;
        result.storage[0][1] -= 3600. * 1e-5;
        let checked = crate::trajectory::evaluate(&problem, result, vec![], 1e-7).unwrap();
        assert!(!checked.validation.feasible);
        assert!(checked.validation.water_balance < 1e-8);
        assert!(checked.validation.pressure_balance > 0.01);
    }

    #[test]
    fn closed_gates_ignore_the_pressure_law_even_with_zero_added_resistance() {
        let mut case = reservoir_case();
        case.system.tunnels = vec![crate::Tunnel {
            name: "Closed".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 0.,
            loss_factor: 1.,
            opening: 0.,
            discharge_river: None,
            gate_positions: vec![],
            initial_opening: None,
            gate_adjustment_cost: 0.,
            gate_loss_factor_curve: Some(
                crate::TableCurve::new(vec![0., 1.], vec![0., 0.]).unwrap(),
            ),
        }];
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.tunnel_flow[0][0], 0.);
        assert!(result.head[0][0] - result.head[1][0] > 70.);
        case.system.tunnels[0].gate_loss_factor_curve =
            Some(crate::TableCurve::new(vec![0., 1.], vec![-1., 0.]).unwrap());
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
        case.system.tunnels[0].gate_loss_factor_curve =
            Some(crate::TableCurve::new(vec![0.1, 1.], vec![2., 0.]).unwrap());
        assert!(matches!(
            SchedulingProblem::new(&case),
            Err(Error::Input(_))
        ));
    }

    #[test]
    fn a_closed_tunnel_disconnects_pressure_without_equating_heads() {
        let mut case = reservoir_case();
        case.system.tunnels.push(crate::Tunnel {
            name: "Closed".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 50.,
            loss_factor: 0.01,
            opening: 0.,
            discharge_river: None,
            gate_positions: vec![],
            initial_opening: None,
            gate_adjustment_cost: 0.,
            gate_loss_factor_curve: None,
        });
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.tunnel_flow[0][0], 0.);
        assert!(result.head[0][0] - result.head[1][0] > 70.);
    }

    #[test]
    fn an_interrupted_solve_reports_recomputed_water_balance() {
        let case = reservoir_case();
        let result = solve(
            &case,
            &SolveOptions {
                max_iterations: 0,
                ..Default::default()
            },
        )
        .unwrap();
        let balance = (result.storage[0][1] - case.system.reservoirs[0].initial_storage) / 3600.
            + result.unit_flow[0][0];
        assert!((result.validation.water_balance - balance.abs()).abs() < 1e-9);
        if balance.abs() > 1e-7 {
            assert!(!result.validation.feasible);
        }
    }
}
