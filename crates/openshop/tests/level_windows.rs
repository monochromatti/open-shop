use openshop::*;
use std::collections::BTreeMap;
fn case() -> ScheduleCase {
    ScheduleCase {
        flow_requirements: vec![],
        routing_grid: vec![],
        operations: vec![],
        system: HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(0., 1e6),
                initial_storage: 5e5,
                level: StorageCurve::Polynomial {
                    reference_volume: 0.,
                    elevation: 50.,
                    slope: 1e-4,
                    curvature: 0.,
                },
                water_value: 0.,
                operation: Default::default(),
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            plants: vec![Plant::new("Station", "Lake", "Sea")],
            units: vec![Unit::new(
                "Unit",
                "Station",
                Bounds::new(1., 20.),
                Bounds::new(0., 30.),
                Bounds::new(1., 200.),
                Efficiency::constant(0.9),
            )],
            ..Default::default()
        },
        grid: vec![0., 1., 3.],
        prices: vec![100., 10.],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::from([("Unit".into(), vec![true; 2])]),
    }
}
fn solve(c: &ScheduleCase) -> Schedule {
    SchedulingProblem::new(c)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap()
}
#[test]
fn level_ranges_constrain_polynomial_and_tabulated_surfaces() {
    for table in [false, true] {
        let mut c = case();
        if table {
            c.system.reservoirs[0].level = StorageCurve::Table(
                TableCurve::new(vec![0., 5e5, 1e6], vec![50., 100., 150.]).unwrap(),
            );
        }
        c.system.reservoirs[0]
            .operation
            .level_windows
            .push(LevelWindow {
                start_time: 0.,
                end_time: 3.,
                max_change: 3.,
            });
        let s = solve(&c);
        assert!(s.validation.feasible, "{:?}", s.validation);
        assert!((s.storage[0][2] - 470000.).abs() < 0.1);
        assert!((s.unit_flow[0][0] - 19. / 3.).abs() < 1e-4);
    }
}
#[test]
fn a_zero_level_range_keeps_inflow_and_discharge_balanced() {
    let mut c = case();
    c.inflows.insert("Lake".into(), vec![4., 0.]);
    c.system.reservoirs[0]
        .operation
        .level_windows
        .push(LevelWindow {
            start_time: 0.,
            end_time: 1.,
            max_change: 0.,
        });
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_flow[0][0] - 4.).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 20.).abs() < 1e-5);
}
#[test]
fn level_window_history_interpolates_storage_at_the_window_start() {
    let mut c = case();
    let op = &mut c.system.reservoirs[0].operation;
    op.history = Some(StorageHistory {
        times: vec![-1., 0.],
        storage: vec![520000., 500000.],
    });
    op.level_windows.push(LevelWindow {
        start_time: -0.5,
        end_time: 1.,
        max_change: 2.,
    });
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_flow[0][0] - 25. / 9.).abs() < 1e-4);
    c.system.reservoirs[0].operation.level_windows[0].max_change = 0.5;
    assert!(SchedulingProblem::new(&c).is_err());
}
#[test]
fn a_level_window_extending_past_the_horizon_constrains_its_known_prefix() {
    let mut c = case();
    c.system.reservoirs[0]
        .operation
        .level_windows
        .push(LevelWindow {
            start_time: 0.,
            end_time: 10.,
            max_change: 2.,
        });
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.storage[0][2] - 480000.).abs() < 0.1);
}
#[test]
fn identical_endpoint_levels_cannot_hide_an_interior_river_pulse() {
    let mut c = case();
    c.grid = vec![0., 1.];
    c.prices = vec![10.];
    c.commitment.insert("Unit".into(), vec![true]);
    let r = &mut c.system.reservoirs[0];
    r.storage = Bounds::new(0., 2000.);
    r.initial_storage = 500.;
    r.level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 95.,
        slope: 0.01,
        curvature: 0.,
    };
    r.operation.level_windows.push(LevelWindow {
        start_time: 0.,
        end_time: 1.,
        max_change: 7.,
    });
    c.system.units[0].flow = Bounds::new(0.1, 1.);
    c.operations.push(OperationalSeries {
        object: Object::Unit("Unit".into()),
        attribute: Attribute::Discharge,
        times: vec![0.],
        values: vec![0.3],
    });
    c.system.rivers.push(River {
        name: "Pulse".into(),
        source: "Upstream".into(),
        target: "Lake".into(),
        capacity: 2.,
        transport: RiverTransport::Fixed { delay: 0.5 },
        allocation: 1.,
        inflow: 0.1,
        water_value: 0.,
        law: RiverLaw::Controlled,
        gate: Bounds::new(0., 1.),
        history_grid: vec![-0.5, -0.25, 0.],
        history_release: vec![1., 0.],
    });
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.storage[0][0] - s.storage[0][1]).abs() < 1e-4);
    c.system.reservoirs[0].operation.level_windows[0].max_change = 6.;
    assert!(!solve(&c).validation.feasible);
}
#[test]
fn level_windows_require_complete_history_and_represented_endpoints() {
    let mut c = case();
    c.system.reservoirs[0]
        .operation
        .level_windows
        .push(LevelWindow {
            start_time: -1.,
            end_time: 1.,
            max_change: 3.,
        });
    assert!(SchedulingProblem::new(&c).is_err());
    c.system.reservoirs[0].operation.level_windows[0].start_time = 0.25;
    assert!(SchedulingProblem::new(&c).is_err());
}
