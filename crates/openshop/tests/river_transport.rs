use openshop::*;
use std::collections::BTreeMap;
fn response() -> RiverTransport {
    RiverTransport::Distributed {
        curves: vec![
            DelayCurve::new(0., vec![0., 0.4, 1.2], vec![0.2, 0.8]).unwrap(),
            DelayCurve::new(40., vec![0., 0.2, 0.8], vec![0.6, 0.4]).unwrap(),
        ],
    }
}
#[test]
fn release_volume_partitions_into_arrivals_and_transit_on_unequal_grids() {
    for transport in [RiverTransport::Fixed { delay: 0.3 }, response()] {
        for q in [0., 3., 17., 40.] {
            for release in [(-0.8, -0.1), (0., 0.01), (0.7, 1.9)] {
                let end = 1.1;
                let grid = [-2., -0.65, 0., 0.07, 0.83, end];
                let arrived: f64 = grid
                    .windows(2)
                    .map(|w| transport.arrival_volume(release, q, (w[0], w[1])).unwrap())
                    .sum();
                let transit = transport.remaining_volume(release, q, end).unwrap();
                let released = (release.1.min(end) - release.0).max(0.) * q * 3600.;
                assert!(
                    (arrived + transit - released).abs() < 1e-8,
                    "{arrived} + {transit} != {released}"
                );
            }
        }
    }
}
#[test]
fn narrow_releases_and_small_late_tails_do_not_lose_volume_to_cancellation() {
    let transport = RiverTransport::Distributed {
        curves: vec![DelayCurve::new(0., vec![0., 1.], vec![1.]).unwrap()],
    };
    let tiny = transport
        .arrival_volume((0., 1e-12), 1., (0., 1e-6))
        .unwrap();
    assert!((tiny / (3600. * 1e-12 * (1e-6 - 0.5e-12)) - 1.).abs() < 1e-9);
    let tail = transport.remaining_volume((0., 1.), 1., 2. - 1e-7).unwrap();
    assert!((tail / (1800. * 1e-14) - 1.).abs() < 2e-8, "tail {tail}");
    assert!(
        transport
            .arrival_volume((0., 1.), 1., (100., 101.))
            .unwrap()
            == 0.
    );
}
fn reach(name: &str, source: &str, target: &str, transport: RiverTransport) -> River {
    River {
        operation: Default::default(),
        name: name.into(),
        source: source.into(),
        target: target.into(),
        capacity: 40.,
        transport,
        allocation: 1.,
        law: RiverLaw::Controlled,
        gate: Bounds::new(0., 1.),
        inflow: 0.,
        water_value: 0.,
        history_grid: vec![],
        history_release: vec![],
    }
}
fn network(cells: usize) -> ScheduleCase {
    let mut left = reach("Left", "Join", "Merge", response());
    left.allocation = 0.3;
    let mut right = reach(
        "Right",
        "Join",
        "Merge",
        RiverTransport::Fixed { delay: 0.17 },
    );
    right.allocation = 0.7;
    let rivers = vec![
        reach("Inlet", "Lake", "Join", response()),
        left,
        right,
        reach("Outlet", "Merge", "Sea", response()),
    ];
    let mut operations = vec![];
    for attribute in [Attribute::GateMin, Attribute::GateMax] {
        operations.push(OperationalSeries {
            object: Object::River("Inlet".into()),
            attribute,
            times: vec![0., 0.5, 1.5],
            values: vec![0.5, 0., 0.25],
        });
    }
    let grid = vec![0., 0.5, 1.5, 3.];
    let mut routing_grid = vec![];
    for w in grid.windows(2) {
        for k in 0..cells {
            routing_grid.push(w[0] + (w[1] - w[0]) * k as f64 / cells as f64);
        }
    }
    routing_grid.push(3.);
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
                water_value: 0.,
                operation: Default::default(),
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            rivers,
            ..Default::default()
        },
        grid,
        prices: vec![0.; 3],
        soft_constraints: vec![],
        flow_requirements: vec![],
        routing_grid,
        operations,
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
    }
}
fn solve(c: &ScheduleCase) -> Schedule {
    let result = SchedulingProblem::new(c)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 10.,
            ..Default::default()
        })
        .unwrap();
    assert!(
        result.validation.feasible,
        "{:?}, status {:?}",
        result.validation, result.solver_status
    );
    result
}

#[test]
fn minimum_arrival_can_require_a_faster_response_than_at_zero_release() {
    let mut c = network(4);
    c.grid = vec![0., 1.];
    c.prices = vec![100.];
    c.operations.clear();
    c.routing_grid.clear();
    c.system.rivers = vec![reach(
        "Tailrace",
        "Tailrace junction",
        "Sea",
        RiverTransport::Distributed {
            curves: vec![
                DelayCurve::new(0., vec![2., 3.], vec![1.]).unwrap(),
                DelayCurve::new(40., vec![0., 0.25], vec![1.]).unwrap(),
            ],
        },
    )];
    let mut plant = Plant::new("Station", "Lake", "Sea");
    plant.discharge_river = Some("Tailrace".into());
    c.system.plants.push(plant);
    c.system.units.push(Unit::new(
        "Unit",
        "Station",
        Bounds::new(5., 20.),
        Bounds::new(1., 30.),
        Bounds::new(90., 110.),
        Efficiency::constant(0.9),
    ));
    c.system.rivers[0].operation.minimum_arrival = 2.;
    let result = solve(&c);
    assert!(result.controls().unwrap().commitment[0][0]);
    let arrival: f64 = result.river_arrival_volume[0].iter().sum();
    assert!(arrival >= 7200. - 1e-5, "{arrival}");
    let release: f64 = result.unit_flow[0]
        .iter()
        .zip(result.grid.windows(2))
        .map(|(q, interval)| q * (interval[1] - interval[0]) * 3600.)
        .sum();
    assert!((arrival + result.river_storage[0].last().unwrap() - release).abs() < 1e-6);
}

#[test]
fn flow_dependent_split_merge_transport_conserves_every_reservoir_withdrawal() {
    let c = network(4);
    let s = solve(&c);
    for t in 0..3 {
        let seconds = 3600. * (c.grid[t + 1] - c.grid[t]);
        assert!(
            (s.river_release[1][t] - 0.3 * s.river_arrival_volume[0][t] / seconds).abs() < 1e-7
        );
        assert!(
            (s.river_release[2][t] - 0.7 * s.river_arrival_volume[0][t] / seconds).abs() < 1e-7
        );
        assert!(
            (s.river_release[3][t]
                - (s.river_arrival_volume[1][t] + s.river_arrival_volume[2][t]) / seconds)
                .abs()
                < 1e-7
        );
        let water = s.storage[0][t + 1]
            + s.river_storage.iter().map(|r| r[t + 1]).sum::<f64>()
            + s.river_arrival_volume[3][..=t].iter().sum::<f64>();
        assert!((water - 5e5).abs() < 1e-5, "{water}");
    }
}
#[test]
fn transport_cell_refinement_converges_for_changing_flows_at_multiple_river_junctions() {
    let totals: Vec<_> = [2, 4, 8, 16]
        .iter()
        .map(|&cells| {
            solve(&network(cells)).river_arrival_volume[3]
                .iter()
                .sum::<f64>()
        })
        .collect();
    assert!(
        (totals[2] - totals[3]).abs() < (totals[0] - totals[3]).abs() * 0.4,
        "{totals:?}"
    );
}

#[test]
fn finer_replay_refines_river_mixing_without_new_gate_decisions() {
    let c = network(1);
    let problem = SchedulingProblem::new(&c).unwrap();
    let controls = Controls {
        head_guess: None,
        commitment: vec![],
        unit_flow: vec![],
        tunnel_opening: vec![],
        river_opening: vec![vec![0.5, 0., 0.25], vec![0.; 3], vec![0.; 3], vec![0.; 3]],
    };
    let mut schedule = problem.reconstruct(&controls).unwrap();
    let mut errors = Vec::new();
    for _ in 0..4 {
        let (fine, error) = problem
            .validate_replay(&schedule, &TimeAccuracy::default())
            .unwrap();
        assert!(fine.validation.feasible, "{:?}", fine.validation);
        assert_eq!(
            fine.controls().unwrap().river_opening,
            controls.river_opening
        );
        let end = fine.grid.len() - 1;
        let accounted = fine.storage[0][end]
            + fine.river_storage.iter().map(|r| r[end]).sum::<f64>()
            + fine.river_arrival_volume[3].iter().sum::<f64>();
        assert!((accounted - 5e5).abs() < 1e-5);
        errors.push(error.arrival_volume_difference[3]);
        schedule = fine;
    }
    assert!(errors[0] > 1. && errors[3] < errors[0] * 0.1, "{errors:?}");
}
#[test]
fn invalid_delay_domains_probabilities_and_meshes_are_rejected() {
    assert!(DelayCurve::new(0., vec![0., 1.], vec![0.9]).is_err());
    assert!(DelayCurve::new(0., vec![0., 1.], vec![1., 0.1]).is_err());
    let mut c = network(4);
    c.routing_grid = vec![0., 3.];
    assert!(SchedulingProblem::new(&c).is_err());
    c.routing_grid.clear();
    c.system.rivers[0].transport = RiverTransport::Distributed {
        curves: vec![
            DelayCurve::new(1., vec![0., 1.], vec![1.]).unwrap(),
            DelayCurve::new(40., vec![0., 1.], vec![1.]).unwrap(),
        ],
    };
    assert!(SchedulingProblem::new(&c).is_err());
}

#[test]
fn commitment_optimization_with_outfalls_uses_the_nonlinear_river_responses() {
    let mut c = network(4);
    c.prices = vec![60., 90., 40.];
    let mut plant = Plant::new("Station", "Lake", "Sea");
    plant.discharge_river = Some("Inlet".into());
    c.system.plants.push(plant);
    c.system.units.push(Unit::new(
        "Unit",
        "Station",
        Bounds::new(1., 5.),
        Bounds::new(0.5, 10.),
        Bounds::new(90., 110.),
        Efficiency::constant(0.9),
    ));
    let s = solve(&c);
    assert!(s.commitment[0].iter().all(|&on| on));
    assert!(s.unit_flow[0].iter().all(|q| (*q - 5.).abs() < 1e-5));
    for t in 0..3 {
        let water = s.storage[0][t + 1]
            + s.river_storage.iter().map(|r| r[t + 1]).sum::<f64>()
            + s.river_arrival_volume[3][..=t].iter().sum::<f64>();
        assert!((water - 5e5).abs() < 1e-4);
    }
}
