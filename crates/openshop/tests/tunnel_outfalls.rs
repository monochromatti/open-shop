use openshop::*;
use std::collections::BTreeMap;
fn lake(name: &str, head: f64, volume: f64) -> Reservoir {
    Reservoir {
        name: name.into(),
        storage: Bounds::new(0., 1e6),
        initial_storage: volume,
        level: StorageCurve::Polynomial {
            reference_volume: 0.,
            elevation: head,
            slope: 0.,
            curvature: 0.,
        },
        water_value: 0.,
        operation: Default::default(),
    }
}
fn case(transport: RiverTransport) -> ScheduleCase {
    let mut tunnel = Tunnel::new("Outlet", "Upper", "Tail", 10., 5.);
    tunnel.discharge_river = Some("Reach".into());
    ScheduleCase {
        soft_constraints: vec![],
        flow_requirements: vec![],
        routing_grid: vec![],
        operations: vec![],
        system: HydroSystem {
            reservoirs: vec![lake("Upper", 100., 5e5), lake("Lower", 50., 1e5)],
            boundaries: vec![Boundary {
                name: "Tail".into(),
                head: 20.,
            }],
            tunnels: vec![tunnel],
            rivers: vec![River {
                operation: Default::default(),
                name: "Reach".into(),
                source: "Outfall".into(),
                target: "Lower".into(),
                capacity: 10.,
                transport,
                allocation: 1.,
                inflow: 0.,
                water_value: 0.,
                law: RiverLaw::Controlled,
                gate: Bounds::new(0., 1.),
                history_grid: vec![],
                history_release: vec![],
            }],
            ..Default::default()
        },
        grid: vec![0., 1., 3.],
        prices: vec![100.; 2],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
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
fn tunnel_outfall_keeps_pressure_reference_and_routes_water_only_once() {
    for transport in [
        RiverTransport::Fixed { delay: 0.5 },
        RiverTransport::Distributed {
            curves: vec![
                DelayCurve::new(0., vec![0.25, 0.75], vec![1.]).unwrap(),
                DelayCurve::new(10., vec![0., 0.5], vec![1.]).unwrap(),
            ],
        },
    ] {
        let c = case(transport);
        let s = solve(&c);
        assert!(s.validation.feasible, "{:?}", s.validation);
        for t in 0..2 {
            assert!((s.tunnel_flow[0][t] - 4.).abs() < 1e-6);
            assert!((s.river_release[0][t] - 4.).abs() < 1e-6);
            let expected: f64 = c
                .grid
                .windows(2)
                .map(|w| {
                    c.system.rivers[0]
                        .transport
                        .arrival_volume((w[0], w[1]), 4., (c.grid[t], c.grid[t + 1]))
                        .unwrap()
                })
                .sum();
            assert!((s.river_arrival_volume[0][t] - expected).abs() < 0.01);
            let total = s.storage.iter().map(|v| v[t + 1]).sum::<f64>() + s.river_storage[0][t + 1];
            assert!((total - 6e5).abs() < 0.01);
        }
    }
}
#[test]
fn a_pressure_junction_balances_a_tunnel_outfall() {
    let mut c = case(RiverTransport::Fixed { delay: 0.5 });
    c.system.junctions.push(Junction {
        name: "Hub".into(),
        head: Bounds::new(20., 100.),
    });
    c.system.tunnels[0].source = "Hub".into();
    c.system.tunnels[0].loss_factor = 3.;
    c.system
        .tunnels
        .insert(0, Tunnel::new("Feeder", "Upper", "Hub", 10., 2.));
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!(
        s.tunnel_flow
            .iter()
            .flatten()
            .all(|q| (*q - 4.).abs() < 1e-6)
    );
    assert!(s.head[2].iter().all(|h| (*h - 68.).abs() < 1e-6));
}
#[test]
fn a_directed_outfall_cannot_create_reverse_river_release() {
    let mut c = case(RiverTransport::Fixed { delay: 0.5 });
    c.system.tunnels[0].source = "Tail".into();
    c.system.tunnels[0].target = "Upper".into();
    assert!(!solve(&c).validation.feasible);
    c.system.tunnels[0].opening = 0.;
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!(s.tunnel_flow[0].iter().all(|q| q.abs() < 1e-8));
}
#[test]
fn commitment_and_flow_observations_include_tunnel_outfalls() {
    let mut c = case(RiverTransport::Fixed { delay: 0.5 });
    c.system.plants.push(Plant::new("Station", "Lower", "Tail"));
    c.system.units.push(Unit::new(
        "Unit",
        "Station",
        Bounds::new(1., 5.),
        Bounds::new(0., 5.),
        Bounds::new(20., 40.),
        Efficiency::constant(0.9),
    ));
    let mut rule = FlowRequirement::new("Gauge");
    rule.tunnels.push("Outlet".into());
    rule.river_arrivals.push("Reach".into());
    rule.max_flow = Some(8.);
    c.flow_requirements.push(rule);
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!(s.commitment[0].iter().all(|on| *on));
    assert!(s.unit_flow[0].iter().all(|q| (*q - 5.).abs() < 1e-5));
    assert!((s.observed_flow[0][0] - 6.).abs() < 1e-5);
    assert!((s.observed_flow[0][1] - 8.).abs() < 1e-5);
}

#[test]
fn restarting_preserves_delayed_arrivals_and_counted_observation_water() {
    for transport in [
        RiverTransport::Fixed { delay: 0.5 },
        RiverTransport::Distributed {
            curves: vec![
                DelayCurve::new(0., vec![0.25, 1.75], vec![1.]).unwrap(),
                DelayCurve::new(10., vec![0.25, 1.25], vec![1.]).unwrap(),
            ],
        },
    ] {
        let mut original = case(transport);
        original.system.rivers[0].history_grid = vec![-1., -0.5];
        original.system.rivers[0].history_release = vec![2.];
        original.system.reservoirs[1]
            .operation
            .level_windows
            .push(LevelWindow {
                start_time: 0.,
                end_time: 3.,
                max_change: 1.,
            });
        let mut observed = FlowRequirement::new("Downstream gauge");
        observed.river_arrivals.push("Reach".into());
        observed.initial_volume = 3000.;
        observed.min_volume = Some(3000.);
        observed.volume_windows.push(VolumeWindow {
            start_time: 0.5,
            end_time: 3.,
            min_volume: Some(0.),
            max_volume: Some(60_000.),
        });
        original.flow_requirements.push(observed);
        let problem = SchedulingProblem::new(&original).unwrap();
        let full = solve(&original);
        let next = problem.restart_case(&full, 1.).unwrap();
        let restarted = SchedulingProblem::new(&next).unwrap();
        let suffix = restarted
            .reconstruct(&Controls {
                head_guess: None,
                commitment: vec![],
                unit_flow: vec![],
                river_opening: vec![vec![0.]],
                tunnel_opening: vec![vec![1.]],
            })
            .unwrap();
        assert!(suffix.validation.feasible, "{:?}", suffix.validation);
        for (expected, actual) in full.storage.iter().zip(&suffix.storage) {
            assert!((expected.last().unwrap() - actual.last().unwrap()).abs() < 1e-6);
        }
        assert!(
            (full.river_storage[0].last().unwrap() - suffix.river_storage[0].last().unwrap()).abs()
                < 1e-6
        );
        assert!(
            (full.observed_volume[0].last().unwrap() - suffix.observed_volume[0].last().unwrap())
                .abs()
                < 1e-6
        );
        assert!(problem.restart_case(&full, 0.5).is_err());
        assert!(problem.restart_case(&full, 3.).is_err());
    }
}
