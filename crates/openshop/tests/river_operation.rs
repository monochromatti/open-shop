use openshop::*;
use std::collections::BTreeMap;

fn case(transport: RiverTransport) -> ScheduleCase {
    ScheduleCase {
        system: HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(0., 1e6),
                initial_storage: 5e5,
                level: StorageCurve::Polynomial {
                    reference_volume: 0.,
                    elevation: 100.,
                    slope: 0.,
                    curvature: 0.,
                },
                water_value: 1.,
                operation: Default::default(),
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            rivers: vec![River {
                name: "Reach".into(),
                source: "Lake".into(),
                target: "Sea".into(),
                capacity: 20.,
                transport,
                allocation: 1.,
                inflow: 0.,
                water_value: 0.,
                history_grid: vec![-2., 0.],
                history_release: vec![0.],
                law: RiverLaw::Controlled,
                gate: Bounds::new(0., 1.),
                operation: RiverOperation {
                    minimum_arrival: 3.,
                    ..Default::default()
                },
            }],
            ..Default::default()
        },
        grid: vec![0., 1., 2.],
        prices: vec![0.; 2],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
        operations: vec![],
        flow_requirements: vec![],
        routing_grid: vec![],
    }
}

fn controls(opening: Vec<f64>) -> Controls {
    Controls {
        commitment: vec![],
        unit_flow: vec![],
        river_opening: vec![opening],
        tunnel_opening: vec![],
    }
}

#[test]
fn instantaneous_arrivals_preserve_fixed_delay_jumps_and_distributed_waves() {
    let fixed = RiverTransport::Fixed { delay: 0.25 };
    for (time, side, expected) in [
        (0.25, TimeSide::Left, 0.),
        (0.25, TimeSide::Right, 10.),
        (1.25, TimeSide::Left, 10.),
        (1.25, TimeSide::Right, 0.),
    ] {
        assert_eq!(
            fixed.arrival_rate((0., 1.), 10., time, side).unwrap(),
            expected
        );
    }
    let distributed = RiverTransport::Distributed {
        curves: vec![DelayCurve::new(0., vec![0., 1.], vec![1.]).unwrap()],
    };
    for (time, expected) in [(0., 0.), (0.5, 5.), (1., 10.), (1.5, 5.), (2., 0.)] {
        for side in [TimeSide::Left, TimeSide::Right] {
            assert!(
                (distributed.arrival_rate((0., 1.), 10., time, side).unwrap() - expected).abs()
                    < 1e-12
            );
        }
    }
    assert!(
        distributed
            .arrival_rate((0., 1.), -1., 0.5, TimeSide::Right)
            .is_err()
    );
    assert!(
        distributed
            .arrival_rate((0., 1.), 1., f64::NAN, TimeSide::Right)
            .is_err()
    );
}

#[test]
fn interval_minima_do_not_imply_an_instantaneous_minimum() {
    let mut c = case(RiverTransport::Fixed { delay: 0.5 });
    let u = controls(vec![0.4, 0.]);
    let average = SchedulingProblem::new(&c).unwrap().reconstruct(&u).unwrap();
    assert!(average.validation.feasible, "{:?}", average.validation);
    for volume in &average.river_arrival_volume[0] {
        assert!((volume / 3600. - 4.).abs() < 1e-8);
    }
    c.system.rivers[0].operation.arrival_policy = ArrivalPolicy::Pointwise;
    assert!(
        !SchedulingProblem::new(&c)
            .unwrap()
            .reconstruct(&u)
            .unwrap()
            .validation
            .feasible
    );
}

#[test]
fn arrival_windows_integrate_the_wave_instead_of_rebinning_control_means() {
    let mut c = case(RiverTransport::Fixed { delay: 0.5 });
    c.system.rivers[0].operation.arrival_window_grid = vec![0., 0.75, 2.];
    let u = controls(vec![0.4, 0.]);
    // The first window contains 0.25 h at 8 m³/s: its mean is 8/3,
    // although the first full control interval has a mean of 4 m³/s.
    assert!(
        !SchedulingProblem::new(&c)
            .unwrap()
            .reconstruct(&u)
            .unwrap()
            .validation
            .feasible
    );
    c.system.rivers[0].operation.minimum_arrival = 2.5;
    assert!(
        SchedulingProblem::new(&c)
            .unwrap()
            .reconstruct(&u)
            .unwrap()
            .validation
            .feasible
    );
}

#[test]
fn both_sides_of_a_requirement_change_use_their_own_minimum() {
    let mut c = case(RiverTransport::Fixed { delay: 0. });
    c.system.rivers[0].operation.arrival_policy = ArrivalPolicy::Pointwise;
    c.operations.push(OperationalSeries {
        object: Object::River("Reach".into()),
        attribute: Attribute::MinArrival,
        times: vec![0., 1.],
        values: vec![3., 8.],
    });
    let schedule = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&controls(vec![0.15, 0.4]))
        .unwrap();
    assert!(schedule.validation.feasible, "{:?}", schedule.validation);
    assert!(
        !SchedulingProblem::new(&c)
            .unwrap()
            .reconstruct(&controls(vec![0.15, 0.3]))
            .unwrap()
            .validation
            .feasible
    );
}

#[test]
fn optimized_pointwise_minima_include_history_and_flow_dependent_delay_blends() {
    for transport in [
        RiverTransport::Fixed { delay: 0.5 },
        RiverTransport::Distributed {
            curves: vec![
                DelayCurve::new(0., vec![0., 1.], vec![1.]).unwrap(),
                DelayCurve::new(20., vec![0., 0.5], vec![1.]).unwrap(),
            ],
        },
    ] {
        let mut c = case(transport);
        c.system.rivers[0].history_release = vec![3.];
        c.system.rivers[0].operation.arrival_policy = ArrivalPolicy::Pointwise;
        let p = SchedulingProblem::new(&c).unwrap();
        let result = p
            .solve(&SolveOptions {
                time_limit: 5.,
                ..Default::default()
            })
            .unwrap();
        assert!(result.validation.feasible, "{:?}", result.validation);
        let r = &c.system.rivers[0];
        for j in 0..=200 {
            let time = 2. * j as f64 / 200.;
            let side = if j == 200 {
                TimeSide::Left
            } else {
                TimeSide::Right
            };
            let mut rate = r.transport.arrival_rate((-2., 0.), 3., time, side).unwrap();
            for t in 0..2 {
                rate += r
                    .transport
                    .arrival_rate(
                        (c.grid[t], c.grid[t + 1]),
                        result.river_release[0][t],
                        time,
                        side,
                    )
                    .unwrap();
            }
            assert!(rate >= 3. - 2e-6, "arrival {rate} at {time}");
        }
    }
}

#[test]
fn invalid_arrival_inputs_are_rejected_before_optimization() {
    for grid in [
        vec![0.1, 2.],
        vec![0., 2.1],
        vec![0., 1., 0.5, 2.],
        vec![0.],
    ] {
        let mut c = case(RiverTransport::default());
        c.system.rivers[0].operation.arrival_window_grid = grid;
        assert!(matches!(SchedulingProblem::new(&c), Err(Error::Input(_))));
    }
    let mut c = case(RiverTransport::default());
    c.system.rivers[0].operation.minimum_arrival = -1.;
    assert!(matches!(SchedulingProblem::new(&c), Err(Error::Input(_))));
}
