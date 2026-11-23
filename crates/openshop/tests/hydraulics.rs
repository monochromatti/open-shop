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
