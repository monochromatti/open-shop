//! Nonlinear dispatch and independent evaluation of its physical equations.

use crate::dispatch::build_dispatch;
use crate::trajectory::{Schedule, Validation};
use crate::{
    Error, ScheduleCase, SolveOptions,
    model::Topology,
    nlp::{self, Backend},
    river::Routing,
};
use std::time::Instant;

/// Optimize on/off schedules with SCIP proposals and exact nonlinear dispatch.
/// Prescribed commitments remain fixed; other units can change commitment
/// between nonlinear solves. This search does not establish a global bound.
/// An infeasible or unfinished solve is reported through `validation.feasible`;
/// a solver's success flag alone does not establish physical feasibility.
pub struct SchedulingProblem<'a> {
    case: &'a ScheduleCase,
    topology: Topology,
    routing: Routing,
}

impl<'a> SchedulingProblem<'a> {
    /// Validate inputs and prepare network incidence and transport once.
    pub fn new(case: &'a ScheduleCase) -> Result<Self, Error> {
        let topology = case.topology()?;
        let routing = Routing::new(case, &topology)?;
        Ok(Self {
            case,
            topology,
            routing,
        })
    }
    pub fn case(&self) -> &ScheduleCase {
        self.case
    }
    pub fn solve(&self, options: &SolveOptions) -> Result<Schedule, Error> {
        self.solve_with_backend(options, Backend::Ipopt)
    }
    fn solve_with_backend(
        &self,
        options: &SolveOptions,
        backend: Backend,
    ) -> Result<Schedule, Error> {
        options.validate()?;
        let (case, topology, routing) = (self.case, &self.topology, &self.routing);
        if case.commitment.len() == case.system.units.len() {
            let commitment = case
                .system
                .units
                .iter()
                .map(|u| case.commitment[&u.name].clone())
                .collect();
            return dispatch(case, topology, routing, commitment, None, options, backend);
        }
        let started = Instant::now();
        let mut excluded = Vec::new();
        let mut best: Option<Schedule> = None;
        let mut last = None;
        while options.time_limit - started.elapsed().as_secs_f64() > 0.01 {
            let remaining = options.time_limit - started.elapsed().as_secs_f64();
            let Some(proposal) = crate::commitment::propose(
                case,
                topology,
                routing,
                best.as_ref(),
                &excluded,
                remaining.min(options.time_limit / 4.),
            )?
            else {
                break;
            };
            excluded.push(proposal.on.clone());
            let mut budget = options.clone();
            budget.time_limit =
                (options.time_limit - started.elapsed().as_secs_f64()).min(options.time_limit / 4.);
            if budget.time_limit <= 0. {
                break;
            }
            budget.max_iterations = options.max_iterations.min(200);
            let result = dispatch(
                case,
                topology,
                routing,
                proposal.on.clone(),
                Some(&proposal),
                &budget,
                backend,
            )?;
            if result.validation.feasible
                && best.as_ref().is_none_or(|p| result.objective > p.objective)
            {
                best = Some(result.clone());
            }
            last = Some(result);
        }
        let result = best.or(last).ok_or_else(|| {
            Error::Solver(
                "no commitment produced a nonlinear trial point within the time limit".into(),
            )
        })?;
        let remaining = options.time_limit - started.elapsed().as_secs_f64();
        if result.validation.feasible && remaining > 0.01 {
            let start = crate::commitment::Proposal {
                on: result.commitment.clone(),
                storage: result.storage.clone(),
                flow: result.unit_flow.clone(),
                tunnel_flow: result.tunnel_flow.clone(),
                river_withdrawal: result.river_withdrawal.clone(),
            };
            let mut budget = options.clone();
            budget.time_limit = remaining;
            let refined = dispatch(
                case,
                topology,
                routing,
                result.commitment.clone(),
                Some(&start),
                &budget,
                backend,
            )?;
            if refined.validation.feasible && refined.objective >= result.objective {
                return Ok(refined);
            }
        }
        Ok(result)
    }
}

fn dispatch(
    case: &ScheduleCase,
    topology: &Topology,
    routing: &Routing,
    commitment: Vec<Vec<bool>>,
    start: Option<&crate::commitment::Proposal>,
    options: &SolveOptions,
    backend: Backend,
) -> Result<Schedule, Error> {
    let dispatch = build_dispatch(case, topology, routing, &commitment, start)?;
    let point = nlp::solve(dispatch.problem, options, backend)?;
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
    let unit_flow = dispatch
        .flow
        .iter()
        .map(|row| {
            row.iter()
                .map(|i| dispatch.flow_scale * point.x[i.index()])
                .collect()
        })
        .collect();
    let tunnel_flow = dispatch
        .tunnel_flow
        .iter()
        .map(|row| {
            row.iter()
                .map(|i| dispatch.flow_scale * point.x[i.index()])
                .collect()
        })
        .collect();
    let mut schedule = Schedule {
        storage,
        unit_flow,
        tunnel_flow,
        river_withdrawal: vec![vec![0.; case.prices.len()]; case.system.rivers.len()],
        river_release: Vec::new(),
        river_arrival_volume: Vec::new(),
        river_storage: Vec::new(),
        unit_power: Vec::new(),
        head: Vec::new(),
        commitment,
        objective: 0.,
        solver_status: point.status,
        iterations: point.iterations,
        validation: Validation::default(),
    };
    for (k, &river) in routing.outlets.iter().enumerate() {
        schedule.river_withdrawal[river] = dispatch.river_withdrawal[k]
            .iter()
            .map(|i| dispatch.flow_scale * point.x[i.index()])
            .collect();
    }
    let junction_head = dispatch
        .junction_head
        .iter()
        .map(|row| row.iter().map(|i| 100. * point.x[i.index()]).collect())
        .collect();
    let schedule = crate::trajectory::evaluate(
        case,
        topology,
        routing,
        schedule,
        junction_head,
        options.tolerance.max(1e-7),
    );
    Ok(schedule)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Boundary, Bounds, Efficiency, HydroSystem, Reservoir, StorageCurve, Unit};
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
            system: HydroSystem {
                reservoirs: vec![Reservoir {
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
                units: vec![Unit {
                    name: "Unit".into(),
                    source: "Lake".into(),
                    target: "Sea".into(),
                    discharge_river: None,
                    flow: Bounds::new(5., 20.),
                    power: Bounds::new(1., 25.),
                    head: Bounds::new(35., 105.),
                    efficiency: Efficiency::constant(0.9),
                    loss_factor: 0.,
                    operation: crate::UnitOperation::default(),
                }],
                ..Default::default()
            },
            grid: vec![0., 1.],
            prices: vec![100.],
            inflows: BTreeMap::new(),
            commitment: BTreeMap::from([("Unit".into(), vec![true])]),
        }
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
        assert!((result.objective - (3000. - 137.)).abs() < 1e-6);
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
            let expected = 100. * 0.00981 * 20. * 0.9 * (mean - 20.) + 2928.;
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
            serde_json::from_str(include_str!("../examples/data/measured_cascade.json")).unwrap();
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
            assert!((result.objective - (100. * 0.00981 * 20. * 80. * 0.9 + 2928.)).abs() < 0.01);
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
            let expected = 100. * expected_power + 0.001 * 2_928_000.;
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
        case.system.units[0].source = "Intake".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        case.system.tunnels = vec![
            Tunnel {
                name: "Main".into(),
                source: "Lake".into(),
                target: "Intake".into(),
                capacity: 20.,
                loss_factor: 0.1,
                opening: 1.,
            },
            Tunnel {
                name: "Parallel".into(),
                source: "Intake".into(),
                target: "Lake".into(),
                capacity: 20.,
                loss_factor: 0.4,
                opening: 1.,
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
        case.system.units[0].target = "Forebay".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        let mut downstream = case.system.units[0].clone();
        downstream.name = "Downstream".into();
        downstream.source = "Forebay".into();
        downstream.target = "Sea".into();
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
        case.system.units[0].target = "Tailwater".into();
        case.system.units[0].head = Bounds::new(1., 120.);
        case.system.units[0].discharge_river = Some("River".into());
        let mut lower = case.system.units[0].clone();
        lower.name = "Lower".into();
        lower.source = "Forebay".into();
        lower.target = "Sea".into();
        lower.discharge_river = None;
        case.system.units.push(lower);
        let mut river = crate::River {
            name: "River".into(),
            source: "Outfall".into(),
            target: "Forebay".into(),
            delay: 0.5,
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
        unit.efficiency = Efficiency {
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
    fn a_closed_tunnel_disconnects_pressure_without_equating_heads() {
        let mut case = reservoir_case();
        case.system.tunnels.push(crate::Tunnel {
            name: "Closed".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            capacity: 50.,
            loss_factor: 0.01,
            opening: 0.,
        });
        let result = solve(&case, &SolveOptions::default()).unwrap();
        assert!(result.validation.feasible);
        assert_eq!(result.tunnel_flow[0][0], 0.);
        assert!(result.head[0][0] - result.head[1][0] > 70.);
    }

    #[test]
    fn an_interrupted_solve_does_not_claim_physical_feasibility() {
        let result = solve(
            &reservoir_case(),
            &SolveOptions {
                max_iterations: 0,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(!result.validation.feasible);
        assert!(result.validation.water_balance > 1.);
    }
}
