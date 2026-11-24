use openshop::*;
use std::collections::BTreeMap;

fn pressure_loop() -> ScheduleCase {
    ScheduleCase {
        system: HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(0., 1e6),
                initial_storage: 5e5,
                water_value: 0.,
                operation: Default::default(),
                level: StorageCurve::Polynomial {
                    reference_volume: 0.,
                    elevation: 100.,
                    slope: 1e-4,
                    curvature: 0.,
                },
            }],
            junctions: vec![Junction {
                name: "Intake".into(),
                head: Bounds::new(-1000., 1000.),
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            plants: vec![Plant::new("Station", "Intake", "Sea")],
            units: vec![Unit::new(
                "Unit",
                "Station",
                Bounds::new(1., 40.),
                Bounds::new(0., 100.),
                Bounds::new(1., 300.),
                Efficiency::constant(0.9),
            )],
            tunnels: vec![
                Tunnel::new("Left", "Lake", "Intake", 100., 0.01),
                Tunnel::new("Right", "Intake", "Lake", 100., 0.04),
            ],
            ..Default::default()
        },
        grid: vec![0., 1., 1.5],
        prices: vec![100., 20.],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
        operations: vec![],
        soft_constraints: vec![],
        flow_requirements: vec![],
        routing_grid: vec![],
    }
}

#[test]
fn prescribed_controls_reconstruct_pressure_loops_energy_and_reversed_flows() {
    let case = pressure_loop();
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![30.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let schedule = SchedulingProblem::new(&case)
        .unwrap()
        .reconstruct(&controls)
        .unwrap();
    assert!(schedule.validation.feasible, "{:?}", schedule.validation);
    for t in 0..2 {
        assert!((schedule.tunnel_flow[0][t] - 20.).abs() < 1e-8);
        assert!((schedule.tunnel_flow[1][t] + 10.).abs() < 1e-8);
        let lake = &case.system.reservoirs[0];
        let mean = lake
            .level
            .mean_level(schedule.storage[0][t], schedule.storage[0][t + 1])
            .unwrap();
        assert!((schedule.head[1][t] - (mean - 4.)).abs() < 1e-8);
        assert!((schedule.unit_power[0][t] - 0.00981 * 30. * (mean - 4.) * 0.9).abs() < 1e-8);
    }
    assert!((schedule.storage[0][2] - 338000.).abs() < 1e-6);
    assert_eq!(controls.unit_flow, vec![vec![30.; 2]]);
}

#[test]
fn reconstruction_reports_overdraw_instead_of_clipping_storage() {
    let mut case = pressure_loop();
    case.system.reservoirs[0].storage = Bounds::new(0., 1000.);
    case.system.reservoirs[0].initial_storage = 500.;
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![30.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let schedule = SchedulingProblem::new(&case)
        .unwrap()
        .reconstruct(&controls)
        .unwrap();
    assert!(!schedule.validation.feasible);
    assert!((schedule.storage[0][2] + 161500.).abs() < 1e-6);
    assert!(schedule.validation.water_balance < 1e-8);
}

#[test]
fn inconsistent_injection_at_a_disconnected_junction_cannot_be_hidden_by_a_pressure_reference() {
    let mut case = pressure_loop();
    case.system.tunnels.clear();
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![1.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    assert!(matches!(
        SchedulingProblem::new(&case)
            .unwrap()
            .reconstruct(&controls),
        Err(Error::Solver(_))
    ));
}

#[test]
fn optimized_storage_and_transport_agree_with_independent_reconstruction() {
    let case: ScheduleCase =
        serde_json::from_str(include_str!("../examples/data/tabulated_cascade.json")).unwrap();
    let problem = SchedulingProblem::new(&case).unwrap();
    let solution = problem
        .solve(&SolveOptions {
            time_limit: 5.,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    let reconstructed = problem.reconstruct(&solution.controls().unwrap()).unwrap();
    assert!(
        reconstructed.validation.feasible,
        "{:?}",
        reconstructed.validation
    );
    assert!((reconstructed.objective - solution.objective).abs() < 0.02);
    for (a, b) in reconstructed
        .storage
        .iter()
        .flatten()
        .zip(solution.storage.iter().flatten())
    {
        assert!((a - b).abs() < 0.1, "storage differs by {} m³", a - b);
    }
    for (a, b) in reconstructed
        .river_arrival_volume
        .iter()
        .flatten()
        .zip(solution.river_arrival_volume.iter().flatten())
    {
        assert!((a - b).abs() < 0.1);
    }
}

#[test]
fn invalid_axes_and_off_unit_flows_are_rejected_before_reconstruction() {
    let case = pressure_loop();
    let problem = SchedulingProblem::new(&case).unwrap();
    let mut controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![30.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    controls.commitment[0][0] = false;
    assert!(matches!(
        problem.reconstruct(&controls),
        Err(Error::Input(_))
    ));
    controls.unit_flow[0][0] = 0.;
    controls.tunnel_opening[0].pop();
    assert!(matches!(
        problem.reconstruct(&controls),
        Err(Error::Input(_))
    ));
}

#[test]
fn average_power_below_a_limit_does_not_hide_an_initial_power_peak() {
    let mut case = pressure_loop();
    case.system.units[0].power.max = 13.1;
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![10.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let physical = SchedulingProblem::new(&case)
        .unwrap()
        .reconstruct(&controls)
        .unwrap();
    assert!(physical.unit_power[0][0] < 13.1);
    assert!(!physical.validation.feasible);
    assert!(physical.validation.operating_limits > 0.001);
}

#[test]
fn optimization_respects_instantaneous_power_as_well_as_interval_average_power() {
    let mut case = pressure_loop();
    case.system.units[0].power.max = 13.1;
    case.commitment.insert("Unit".into(), vec![true; 2]);
    let problem = SchedulingProblem::new(&case).unwrap();
    let solution = problem
        .solve(&SolveOptions {
            time_limit: 5.,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    let q = solution.unit_flow[0][0];
    let first_head = 150. - 0.01 * (q * 2. / 3.).powi(2);
    let initial_power = 0.00981 * q * first_head * 0.9;
    assert!((initial_power - 13.1).abs() < 1e-5);
    assert!(solution.unit_power[0][0] < initial_power - 0.1);
    assert!(
        problem
            .reconstruct(&solution.controls().unwrap())
            .unwrap()
            .validation
            .feasible
    );
}

#[test]
fn turbine_table_power_peaks_between_observed_heads_without_extra_time_steps() {
    let mut case = pressure_loop();
    case.system.junctions.clear();
    case.system.tunnels.clear();
    case.system.plants[0].source = "Lake".into();
    case.grid = vec![0., 1.];
    case.prices = vec![100.];
    let reservoir = &mut case.system.reservoirs[0];
    reservoir.initial_storage = 36000.;
    reservoir.level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 130.,
        slope: 100. / 36000.,
        curvature: 0.,
    };
    let unit = &mut case.system.units[0];
    unit.power.max = 11.;
    unit.head = Bounds::new(100., 300.);
    unit.efficiency = Efficiency::Table(
        TurbineTable::new(
            vec![100., 200.],
            vec![
                TableCurve::new(vec![1., 40.], vec![1.; 2]).unwrap(),
                TableCurve::new(vec![1., 40.], vec![0.5; 2]).unwrap(),
            ],
            vec![unit.flow; 2],
            TurbineOptions {
                head_extrapolation: Extrapolation::Linear,
                ..Default::default()
            },
        )
        .unwrap(),
    );
    case.commitment.insert("Unit".into(), vec![true]);
    let problem = SchedulingProblem::new(&case).unwrap();
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true]],
        unit_flow: vec![vec![10.]],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    let physical = problem.reconstruct(&controls).unwrap();
    // Both endpoint powers and the mean-step power lie below eleven MW.
    for head in [130., 180., 230.] {
        assert!(case.system.units[0].power(10., head).unwrap() < 11.);
    }
    assert!(physical.unit_power[0][0] < 11.);
    assert!(!physical.validation.feasible);
    assert!(physical.validation.operating_limits > 0.003);
    let solution = problem
        .solve(&SolveOptions {
            time_limit: 5.,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    let peak = case.system.units[0]
        .power(solution.unit_flow[0][0], 150.)
        .unwrap();
    assert!((peak - 11.).abs() < 1e-5);
    assert!(solution.unit_power[0][0] < 10.7);
}

#[test]
fn fixed_power_limits_allow_a_constant_head_operating_schedule() {
    let mut case = pressure_loop();
    case.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 150.,
        slope: 0.,
        curvature: 0.,
    };
    case.system.units[0].power = Bounds::new(10., 10.);
    case.commitment.insert("Unit".into(), vec![true; 2]);
    let solution = SchedulingProblem::new(&case)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 5.,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    assert!(
        solution.unit_power[0]
            .iter()
            .all(|p| (p - 10.).abs() < 1e-6)
    );
}

#[test]
fn a_measured_level_knot_cannot_hide_a_head_peak_between_control_edges() {
    let mut case = pressure_loop();
    case.system.junctions.clear();
    case.system.tunnels.clear();
    case.system.reservoirs[0].storage = Bounds::new(0., 1e6);
    case.system.reservoirs[0].initial_storage = 72000.;
    case.system.reservoirs[0].level = StorageCurve::Table(
        TableCurve::new(
            vec![0., 36000., 54000., 72000., 1e6],
            vec![200., 200., 250., 250., 250.],
        )
        .unwrap(),
    );
    let mut lower = case.system.reservoirs[0].clone();
    lower.name = "Lower".into();
    lower.level = StorageCurve::Table(
        TableCurve::new(
            vec![0., 36000., 54000., 72000., 1e6],
            vec![100., 100., 100., 150., 150.],
        )
        .unwrap(),
    );
    case.system.reservoirs.push(lower);
    case.system.plants[0] = Plant::new("Station", "Lake", "Lower");
    case.system
        .plants
        .push(Plant::new("Downstream", "Lower", "Sea"));
    case.system.units[0].power.max = 12.;
    case.system.units.push(Unit::new(
        "Downstream unit",
        "Downstream",
        Bounds::new(20., 20.),
        Bounds::new(0., 100.),
        Bounds::new(1., 200.),
        Efficiency::constant(0.9),
    ));
    case.grid = vec![0., 1.];
    case.prices = vec![100.];
    case.commitment.insert("Unit".into(), vec![true]);
    case.commitment.insert("Downstream unit".into(), vec![true]);
    let problem = SchedulingProblem::new(&case).unwrap();
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true]; 2],
        unit_flow: vec![vec![10.], vec![20.]],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    let physical = problem.reconstruct(&controls).unwrap();
    assert!(physical.unit_power[0][0] < 12.);
    assert!(!physical.validation.feasible);
    assert!(physical.validation.operating_limits > 0.1);
    let solution = problem
        .solve(&SolveOptions {
            time_limit: 10.,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    assert!((solution.unit_flow[0][0] - 12. / (0.00981 * 150. * 0.9)).abs() < 1e-5);
    assert!(
        problem
            .reconstruct(&solution.controls().unwrap())
            .unwrap()
            .validation
            .feasible
    );
}

#[test]
fn floating_pressure_reference_respects_absolute_bounds_without_changing_flow() {
    let mut case = pressure_loop();
    case.system.units.clear();
    case.system.plants.clear();
    case.system.junctions = vec![
        Junction {
            name: "A".into(),
            head: Bounds::new(400., 500.),
        },
        Junction {
            name: "B".into(),
            head: Bounds::new(450., 550.),
        },
    ];
    case.system.tunnels = vec![Tunnel::new("Link", "A", "B", 20., 0.1)];
    let controls = Controls {
        head_guess: None,
        commitment: vec![],
        unit_flow: vec![],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]],
    };
    let physical = SchedulingProblem::new(&case)
        .unwrap()
        .reconstruct(&controls)
        .unwrap();
    assert!(physical.validation.feasible, "{:?}", physical.validation);
    for t in 0..2 {
        assert!(physical.head[1][t] >= 450. && physical.head[1][t] <= 500.);
        assert!((physical.head[1][t] - physical.head[2][t]).abs() < 1e-9);
        assert!(physical.tunnel_flow[0][t].abs() < 1e-9);
    }
}

#[test]
fn crossing_an_outlet_head_floor_cannot_hide_an_interior_power_peak() {
    let mut case = pressure_loop();
    case.system.junctions.clear();
    case.system.tunnels.clear();
    let lake = &mut case.system.reservoirs[0];
    lake.initial_storage = 72000.;
    lake.level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 100.,
        slope: 50. / 36000.,
        curvature: 0.,
    };
    let mut lower = lake.clone();
    lower.name = "Lower".into();
    lower.level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: -80.,
        slope: 100. / 36000.,
        curvature: 0.,
    };
    case.system.reservoirs.push(lower);
    case.system.plants[0] = Plant::new("Station", "Lake", "Lower");
    case.system.plants[0].outlet_head_floor = Some(100.);
    case.system
        .plants
        .push(Plant::new("Downstream", "Lower", "Sea"));
    case.system.units[0].power.max = 7.2;
    case.system.units.push(Unit::new(
        "Downstream unit",
        "Downstream",
        Bounds::new(20., 20.),
        Bounds::new(0., 100.),
        Bounds::new(1., 200.),
        Efficiency::constant(0.9),
    ));
    case.grid = vec![0., 1.];
    case.prices = vec![100.];
    case.commitment.insert("Unit".into(), vec![true]);
    case.commitment.insert("Downstream unit".into(), vec![true]);
    let problem = SchedulingProblem::new(&case).unwrap();
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true]; 2],
        unit_flow: vec![vec![10.], vec![20.]],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    let physical = problem.reconstruct(&controls).unwrap();
    for head in [50., 75., 80.] {
        assert!(case.system.units[0].power(10., head).unwrap() < 7.2);
    }
    assert!(physical.unit_power[0][0] < 7.2);
    assert!(!physical.validation.feasible);
    assert!(physical.validation.operating_limits > 0.1);
    let solution = problem
        .solve(&SolveOptions {
            time_limit: 10.,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    let q = solution.unit_flow[0][0];
    let fraction = 2. / (20. - q);
    let peak = 0.00981 * q * (100. - 5. * q * fraction) * 0.9;
    assert!(
        (peak - 7.2).abs() < 1e-5,
        "q {q}, peak {peak}, objective {}, iterations {}, status {:?}",
        solution.objective,
        solution.iterations,
        solution.solver_status
    );
}

#[test]
fn reconstruction_checks_prescribed_commitment_in_addition_to_water_and_pressure() {
    let mut case = pressure_loop();
    case.commitment.insert("Unit".into(), vec![true; 2]);
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![false; 2]],
        unit_flow: vec![vec![0.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let physical = SchedulingProblem::new(&case)
        .unwrap()
        .reconstruct(&controls)
        .unwrap();
    assert!(!physical.validation.feasible);
    assert!(physical.validation.water_balance < 1e-8);
    assert!(physical.validation.pressure_balance < 1e-8);
    assert!(physical.validation.operating_limits > 0.);
}

#[test]
fn a_committed_idle_unit_can_have_negative_head_without_transferring_water_uphill() {
    let mut case = pressure_loop();
    case.system.junctions.clear();
    case.system.tunnels.clear();
    case.system.plants[0].source = "Lake".into();
    case.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: -10.,
        slope: 0.,
        curvature: 0.,
    };
    let u = &mut case.system.units[0];
    u.flow = Bounds::new(0., 20.);
    u.power = Bounds::new(0., 20.);
    u.head = Bounds::new(-20., 20.);
    case.commitment.insert("Unit".into(), vec![true; 2]);
    let problem = SchedulingProblem::new(&case).unwrap();
    let solution = problem
        .solve(&SolveOptions {
            time_limit: 5.,
            tolerance: 1e-11,
            ..Default::default()
        })
        .unwrap();
    assert!(solution.validation.feasible, "{:?}", solution.validation);
    assert!(solution.commitment[0].iter().all(|on| *on));
    assert!(solution.unit_flow[0].iter().all(|q| q.abs() < 1e-8));
    assert!(solution.unit_head[0].iter().all(|h| (h + 10.).abs() < 1e-8));
    assert_eq!(case.system.units[0].power(0., -10.).unwrap(), 0.);
}

#[test]
fn zero_efficiency_cannot_hide_uphill_transfer_through_a_generating_unit() {
    let mut case = pressure_loop();
    case.system.junctions.clear();
    case.system.tunnels.clear();
    case.system.plants[0].source = "Lake".into();
    case.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: -10.,
        slope: 0.,
        curvature: 0.,
    };
    let u = &mut case.system.units[0];
    u.flow = Bounds::new(0., 20.);
    u.power = Bounds::new(0., 20.);
    u.head = Bounds::new(-20., 20.);
    u.efficiency = Efficiency::Quadratic {
        peak: 0.9,
        best_flow: 10.,
        flow_curvature: 0.009,
        best_head: 100.,
        head_curvature: 0.,
    };
    case.commitment.insert("Unit".into(), vec![true; 2]);
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![20.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    let problem = SchedulingProblem::new(&case).unwrap();
    let physical = problem.reconstruct(&controls).unwrap();
    assert!(!physical.validation.feasible);
    assert!(physical.unit_power[0].iter().all(|p| p.abs() < 1e-10));
    assert!(physical.validation.operating_limits >= 0.1);
    assert!(case.system.units[0].power(20., -10.).is_err());
}

#[test]
fn continuous_refinement_improves_revenue_without_changing_discrete_controls() {
    let case = pressure_loop();
    let problem = SchedulingProblem::new(&case).unwrap();
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![2.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let initial = problem.reconstruct(&controls).unwrap();
    assert!(initial.validation.feasible);
    let result = problem
        .refine(
            &controls,
            &SolveOptions {
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    assert!(result.objective > initial.objective + 100.);
    assert_eq!(result.commitment, controls.commitment);
    assert_eq!(result.tunnel_opening, controls.tunnel_opening);
    let interrupted = problem
        .refine(
            &result.controls().unwrap(),
            &SolveOptions {
                max_iterations: 0,
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(interrupted.validation.feasible);
    assert!(interrupted.objective >= result.objective - 1e-6);
}

#[test]
fn continuous_refinement_cannot_remove_a_commitment_time_restriction() {
    let mut case = pressure_loop();
    case.system.units[0].operation.minimum_up_time = 2.;
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true, false]],
        unit_flow: vec![vec![5., 0.]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let result = SchedulingProblem::new(&case)
        .unwrap()
        .refine(
            &controls,
            &SolveOptions {
                time_limit: 3.,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(!result.validation.feasible);
    assert_eq!(result.commitment, controls.commitment);
}

#[test]
fn continuous_refinement_restores_infeasible_flows_before_improving_revenue() {
    let mut case = pressure_loop();
    case.system.units[0].power = Bounds::new(10., 20.);
    let controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![40.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let problem = SchedulingProblem::new(&case).unwrap();
    assert!(!problem.reconstruct(&controls).unwrap().validation.feasible);
    let result = problem
        .refine(
            &controls,
            &SolveOptions {
                time_limit: 5.,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    assert!(
        result.unit_power[0]
            .iter()
            .all(|p| *p >= 10. - 1e-5 && *p <= 20. + 1e-5)
    );
    assert!(result.unit_flow[0].iter().all(|q| *q < 40.));
    assert_eq!(result.commitment, controls.commitment);
    assert_eq!(result.tunnel_opening, controls.tunnel_opening);
}

fn serial_turbines() -> ScheduleCase {
    let mut c = pressure_loop();
    c.grid = vec![0., 1.];
    c.prices = vec![100.];
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 100.,
        slope: 0.,
        curvature: 0.,
    };
    c.system.junctions = vec![Junction {
        name: "Intermediate".into(),
        head: Bounds::new(21., 99.),
    }];
    c.system.boundaries[0].head = 20.;
    c.system.tunnels.clear();
    c.system.plants = vec![
        Plant::new("Upper", "Lake", "Intermediate"),
        Plant::new("Lower", "Intermediate", "Sea"),
    ];
    c.system.units = [0.9, 0.8]
        .into_iter()
        .zip(["Upper", "Lower"])
        .map(|(eta, name)| {
            Unit::new(
                name,
                name,
                Bounds::new(20., 20.),
                Bounds::new(0., 100.),
                Bounds::new(1., 79.),
                Efficiency::constant(eta),
            )
        })
        .collect();
    c.commitment = BTreeMap::from([("Upper".into(), vec![true]), ("Lower".into(), vec![true])]);
    c
}

#[test]
fn optimized_serial_turbines_retain_their_pressure_allocation_during_reconstruction() {
    let c = serial_turbines();
    let p = SchedulingProblem::new(&c).unwrap();
    let s = p
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap();
    assert!(s.validation.feasible, "{:?}", s.validation);
    // A metre of head produces more power in the upper, more efficient machine.
    let expected = 100. * 0.00981 * 20. * (79. * 0.9 + 0.8);
    assert!(
        (s.objective - expected).abs() < 1e-5,
        "revenue {} instead of {expected}",
        s.objective
    );
    let replay = p.reconstruct(&s.controls().unwrap()).unwrap();
    assert!(replay.validation.feasible);
    assert!((replay.objective - s.objective).abs() < 1e-7);
}

#[test]
fn specified_intermediate_pressure_is_retained_when_series_discharge_does_not_determine_it() {
    for pressure in [21., 99.] {
        let mut c = serial_turbines();
        c.system.junctions[0].head = Bounds::new(pressure, pressure);
        let input = Controls {
            head_guess: None,
            commitment: vec![vec![true]; 2],
            unit_flow: vec![vec![20.]; 2],
            river_opening: vec![],
            tunnel_opening: vec![],
        };
        let s = SchedulingProblem::new(&c)
            .unwrap()
            .reconstruct(&input)
            .unwrap();
        assert!(s.validation.feasible, "{:?}", s.validation);
        assert!((s.head[1][0] - pressure).abs() < 1e-9);
        let expected = 100. * 0.00981 * 20. * ((100. - pressure) * 0.9 + (pressure - 20.) * 0.8);
        assert!((s.objective - expected).abs() < 1e-8);
    }
}

#[test]
fn pressure_guesses_select_a_free_allocation_without_overriding_fixed_heads() {
    let c = serial_turbines();
    let input = Controls {
        head_guess: Some(HeadGuess {
            grid: c.grid.clone(),
            head: vec![vec![999.], vec![40.], vec![-999.]],
        }),
        commitment: vec![vec![true]; 2],
        unit_flow: vec![vec![20.]; 2],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    let s = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&input)
        .unwrap();
    assert!(s.validation.feasible);
    assert!((s.head[0][0] - 100.).abs() < 1e-9 && (s.head[2][0] - 20.).abs() < 1e-9);
    assert!((s.objective - 100. * 0.00981 * 20. * (60. * 0.9 + 20. * 0.8)).abs() < 1e-8);
}

#[test]
fn an_incorrect_pressure_guess_cannot_replace_an_anchored_tunnel_pressure_equation() {
    let c = pressure_loop();
    let input = Controls {
        head_guess: Some(HeadGuess {
            grid: c.grid.clone(),
            head: vec![vec![1000.; 2]; 3],
        }),
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![30.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let s = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&input)
        .unwrap();
    assert!(s.validation.feasible, "{:?}", s.validation);
    let r = &c.system.reservoirs[0];
    for t in 0..2 {
        let mean = r
            .level
            .mean_level(s.storage[0][t], s.storage[0][t + 1])
            .unwrap();
        assert!((s.head[1][t] - (mean - 4.)).abs() < 1e-8);
    }
    let mut incompatible = c.clone();
    incompatible.system.junctions[0].head = Bounds::new(0., 10.);
    let outside = SchedulingProblem::new(&incompatible)
        .unwrap()
        .reconstruct(&input)
        .unwrap();
    assert!(!outside.validation.feasible);
    assert!(
        outside.head[1][0] > 10.,
        "the solved state must not be clipped to the initial guess bounds"
    );
    let mut invalid = input;
    invalid.head_guess.as_mut().unwrap().head[0][0] = f64::NAN;
    assert!(matches!(
        SchedulingProblem::new(&c).unwrap().reconstruct(&invalid),
        Err(Error::Input(_))
    ));
}

#[test]
fn finer_state_arrays_recover_held_controls_and_original_energy_windows() {
    let c = pressure_loop();
    let coarse_controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![30.; 2]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 2]; 2],
    };
    let coarse = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&coarse_controls)
        .unwrap();
    let mut fine_case = c.clone();
    fine_case.grid = vec![0., 0.25, 1., 1.25, 1.5];
    fine_case.prices = vec![100., 100., 20., 20.];
    let fine_controls = Controls {
        head_guess: None,
        commitment: vec![vec![true; 4]],
        unit_flow: vec![vec![30.; 4]],
        river_opening: vec![],
        tunnel_opening: vec![vec![1.; 4]; 2],
    };
    let mut fine = SchedulingProblem::new(&fine_case)
        .unwrap()
        .reconstruct(&fine_controls)
        .unwrap();
    assert!(coarse.validation.feasible && fine.validation.feasible);
    fine.operating_grid = c.grid.clone();
    let controls = fine.controls().unwrap();
    let replay = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&controls)
        .unwrap();
    assert!(replay.validation.feasible);
    assert!((replay.objective - coarse.objective).abs() < 1e-7);
    let mean_power = fine.net_power(&c).unwrap();
    for (mean, expected) in mean_power.iter().zip(&coarse.unit_power[0]) {
        assert!((mean - expected).abs() < 1e-9);
    }
    let energy_revenue: f64 = mean_power
        .iter()
        .zip(&c.prices)
        .zip(c.grid.windows(2))
        .map(|((p, price), w)| p * price * (w[1] - w[0]))
        .sum();
    assert!((energy_revenue - fine.objective).abs() < 1e-7);
    fine.commitment[0][1] = false;
    assert!(matches!(fine.controls(), Err(Error::Input(_))));
}

#[test]
fn head_starting_profiles_use_time_window_means_on_different_partitions() {
    let mut c = serial_turbines();
    c.grid = vec![0., 0.5, 1.];
    c.prices = vec![100.; 2];
    c.commitment = BTreeMap::from([
        ("Upper".into(), vec![true; 2]),
        ("Lower".into(), vec![true; 2]),
    ]);
    let mut input = Controls {
        head_guess: Some(HeadGuess {
            grid: vec![0., 0.25, 1.],
            head: vec![vec![100.; 2], vec![21., 73.], vec![20.; 2]],
        }),
        commitment: vec![vec![true; 2]; 2],
        unit_flow: vec![vec![20.; 2]; 2],
        river_opening: vec![],
        tunnel_opening: vec![],
    };
    let p = SchedulingProblem::new(&c).unwrap();
    let s = p.reconstruct(&input).unwrap();
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.head[1][0] - 47.).abs() < 1e-9 && (s.head[1][1] - 73.).abs() < 1e-9);
    assert!((s.objective - 1334.16).abs() < 1e-8);
    input.head_guess.as_mut().unwrap().grid[1] = 1.;
    assert!(matches!(p.reconstruct(&input), Err(Error::Input(_))));
}
