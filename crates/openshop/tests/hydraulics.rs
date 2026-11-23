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
        flow_requirements: vec![],
        routing_grid: vec![],
    }
}

#[test]
fn prescribed_controls_reconstruct_pressure_loops_energy_and_reversed_flows() {
    let case = pressure_loop();
    let controls = Controls {
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
    let reconstructed = problem.reconstruct(&solution.controls()).unwrap();
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
            .reconstruct(&solution.controls())
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
            .reconstruct(&solution.controls())
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
        "q {q}, peak {peak}, objective {}, iterations {}, status {}",
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
