use openshop::*;
use std::collections::BTreeMap;
fn case() -> ScheduleCase {
    ScheduleCase {
        soft_constraints: vec![],
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
            plants: vec![Plant::new("Station", "Lake", "Sea")],
            units: vec![Unit::new(
                "Unit",
                "Station",
                Bounds::new(1., 20.),
                Bounds::new(0., 30.),
                Bounds::new(90., 110.),
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
fn gauge() -> FlowRequirement {
    let mut r = FlowRequirement::new("Gauge");
    r.units.push("Unit".into());
    r
}
fn solve(c: &ScheduleCase) -> Schedule {
    let s = SchedulingProblem::new(c)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap();
    assert!(s.validation.feasible, "{:?}", s.validation);
    s
}
#[test]
fn rate_and_total_volume_limits_preserve_previously_counted_water() {
    let mut c = case();
    let mut r = gauge();
    r.max_flow = Some(8.);
    r.initial_volume = 5000.;
    r.max_volume = Some(86000.);
    c.flow_requirements.push(r);
    let s = solve(&c);
    assert!((s.unit_flow[0][0] - 8.).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 7.25).abs() < 1e-5);
    assert!((s.observed_volume[0][0] - 5000.).abs() < 1e-8);
    assert!((s.observed_volume[0][2] - 86000.).abs() < 0.01);
}
#[test]
fn volume_window_integrates_partial_intervals_without_counting_initial_volume_twice() {
    let mut c = case();
    let mut r = gauge();
    r.initial_volume = 5000.;
    r.volume_windows.push(VolumeWindow {
        start_time: 0.25,
        end_time: 2.5,
        min_volume: Some(40000.),
        max_volume: Some(40000.),
    });
    c.flow_requirements.push(r);
    let s = solve(&c);
    assert!((s.unit_flow[0][0] - 34600. / 2700.).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 1.).abs() < 1e-5);
}
#[test]
fn rolling_average_ramps_check_both_sides_of_a_control_change() {
    let mut c = case();
    let mut r = gauge();
    r.history_grid = vec![-2., -1., 0.];
    r.history_flow = vec![2., 4.];
    r.rolling_ramps.push(RollingRamp {
        hours: 1.,
        up: Some(1.),
        down: None,
    });
    c.flow_requirements.push(r);
    let s = solve(&c);
    assert!((s.unit_flow[0][0] - 5.).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 6.).abs() < 1e-5);
}
#[test]
fn period_ramps_use_time_weighted_historical_averages() {
    let mut c = case();
    let mut r = gauge();
    r.history_grid = vec![-1., -0.25, 0.];
    r.history_flow = vec![2., 10.];
    r.period_ramps.push(PeriodRamp {
        grid: vec![-1., 0., 1., 3.],
        up: Some(2.),
        down: None,
    });
    c.flow_requirements.push(r);
    let s = solve(&c);
    assert!((s.unit_flow[0][0] - 6.).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 8.).abs() < 1e-5);
}
#[test]
fn a_window_crossing_the_start_uses_the_complete_flow_history() {
    let mut c = case();
    let mut r = gauge();
    r.history_grid = vec![-1., -0.25, 0.];
    r.history_flow = vec![2., 10.];
    r.volume_windows.push(VolumeWindow {
        start_time: -0.5,
        end_time: 1.5,
        min_volume: None,
        max_volume: Some(21600.),
    });
    c.flow_requirements.push(r);
    let s = solve(&c);
    assert!((s.unit_flow[0][0] - 2.5).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 1.).abs() < 1e-5);
}
#[test]
fn observations_retain_signed_tunnel_flow() {
    let mut c = case();
    c.system.units.clear();
    c.system.plants.clear();
    c.commitment.clear();
    c.system.tunnels.push(Tunnel {
        name: "Link".into(),
        source: "Sea".into(),
        target: "Lake".into(),
        capacity: 10.,
        loss_factor: 6.25,
        opening: 1.,
        discharge_river: None,
        gate_positions: vec![],
        initial_opening: None,
        gate_adjustment_cost: 0.,
        gate_loss_factor_curve: None,
    });
    let mut r = FlowRequirement::new("Gauge");
    r.tunnels.push("Link".into());
    r.min_flow = -5.;
    r.max_flow = Some(-3.);
    c.flow_requirements.push(r);
    let s = solve(&c);
    assert!(s.observed_flow[0].iter().all(|q| (*q + 4.).abs() < 1e-6));
    assert!((s.observed_volume[0][2] + 43200.).abs() < 0.01);
}
#[test]
fn observed_inflow_and_limits_can_change_in_time() {
    let mut c = case();
    c.flow_requirements.push(gauge());
    for (attribute, values) in [
        (Attribute::Inflow, vec![3., 5.]),
        (Attribute::FlowMax, vec![10., 10.]),
    ] {
        c.operations.push(OperationalSeries {
            object: Object::FlowRequirement("Gauge".into()),
            attribute,
            times: vec![0., 1.],
            values,
        });
    }
    let s = solve(&c);
    assert!((s.unit_flow[0][0] - 7.).abs() < 1e-5);
    assert!((s.unit_flow[0][1] - 5.).abs() < 1e-5);
    assert!(s.observed_flow[0].iter().all(|q| (*q - 10.).abs() < 1e-6));
}
#[test]
fn missing_history_duplicate_flows_and_unfinished_windows_are_rejected() {
    let mut c = case();
    let mut r = gauge();
    r.rolling_ramps.push(RollingRamp {
        hours: 1.,
        up: Some(1.),
        down: None,
    });
    c.flow_requirements.push(r);
    assert!(SchedulingProblem::new(&c).is_err());
    c.flow_requirements = vec![gauge()];
    c.flow_requirements[0].units.push("Unit".into());
    assert!(SchedulingProblem::new(&c).is_err());
    c.flow_requirements = vec![gauge()];
    c.flow_requirements[0].volume_windows.push(VolumeWindow {
        start_time: 0.,
        end_time: 4.,
        min_volume: None,
        max_volume: Some(1e5),
    });
    assert!(SchedulingProblem::new(&c).is_err());
}
