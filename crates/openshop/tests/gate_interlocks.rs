use openshop::*;
use std::collections::BTreeMap;
fn lake() -> Reservoir {
    Reservoir {
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
    }
}
fn case() -> ScheduleCase {
    let e = Tunnel::new("Link", "Lake", "Sea", 10., 5.);
    ScheduleCase {
        soft_constraints: vec![],
        flow_requirements: vec![],
        routing_grid: vec![],
        operations: vec![],
        system: HydroSystem {
            reservoirs: vec![lake()],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            plants: vec![Plant::new("Station", "Lake", "Sea")],
            units: vec![Unit::new(
                "Unit",
                "Station",
                Bounds::new(1., 5.),
                Bounds::new(0., 10.),
                Bounds::new(90., 110.),
                Efficiency::constant(0.9),
            )],
            tunnels: vec![e],
            ..Default::default()
        },
        grid: vec![0., 0.5, 1., 1.5, 2.],
        prices: vec![100.; 4],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
    }
}
fn series(c: &mut ScheduleCase, name: &str, times: Vec<f64>, values: Vec<f64>) {
    c.operations.push(OperationalSeries {
        object: Object::Tunnel(name.into()),
        attribute: Attribute::Opening,
        times,
        values,
    });
}
fn wait(c: &mut ScheduleCase, hours: f64, remaining: f64) {
    let mut rule = TunnelInterlock::new("Switch", vec!["Link".into()]);
    rule.closures.push(ClosureWait {
        tunnels: vec!["Link".into()],
        plants: vec!["Station".into()],
        hours,
        initial_remaining: remaining,
    });
    c.system.tunnel_interlocks.push(rule);
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
fn closing_a_gate_stops_protected_units_and_reopening_does_not_cancel_the_wait() {
    let mut c = case();
    series(&mut c, "Link", vec![0.5, 1.], vec![0., 1.]);
    wait(&mut c, 1., 0.);
    let s = solve(&c);
    assert_eq!(s.commitment[0], vec![true, false, false, true]);
}
#[test]
fn closure_waits_interact_with_minimum_operating_times() {
    let mut c = case();
    c.system.units[0].operation.minimum_up_time = 0.75;
    series(&mut c, "Link", vec![0.5, 1.], vec![0., 1.]);
    wait(&mut c, 1., 0.);
    let s = solve(&c);
    assert_eq!(s.commitment[0], vec![false, false, false, true]);
}
#[test]
fn a_pending_wait_uses_history_and_survives_reopening() {
    let mut c = case();
    c.grid = vec![0., 0.5, 1., 1.5];
    c.prices = vec![100.; 3];
    c.system.tunnels[0].opening = 0.;
    c.system.units[0].operation.initial_duration = 0.25;
    series(&mut c, "Link", vec![0.5], vec![1.]);
    wait(&mut c, 1., 0.75);
    let s = solve(&c);
    assert_eq!(s.commitment[0], vec![false, false, true]);
    c.system.units[0].operation.initial_duration = 0.1;
    assert!(SchedulingProblem::new(&c).is_err());
}
#[test]
fn changing_between_positive_positions_is_not_a_closure() {
    let mut c = case();
    series(&mut c, "Link", vec![0.5], vec![0.5]);
    wait(&mut c, 1., 0.);
    let s = solve(&c);
    assert!(s.commitment[0].iter().all(|on| *on));
}
#[test]
fn positive_positions_count_as_open_at_zero_flow() {
    let mut c = case();
    c.system.units.clear();
    c.system.plants.clear();
    c.grid = vec![0., 1.];
    c.prices = vec![0.];
    c.system.junctions.push(Junction {
        name: "Hub".into(),
        head: Bounds::new(0., 200.),
    });
    c.system.tunnels[0].target = "Hub".into();
    c.system.tunnels[0].opening = 0.;
    c.system.tunnels[0].gate_positions = vec![0., 1.];
    let mut other = c.system.tunnels[0].clone();
    other.name = "Other".into();
    c.system.tunnels.push(other);
    let mut r = TunnelInterlock::new("OneOpen", vec!["Link".into(), "Other".into()]);
    r.min_open = 1;
    r.max_open = Some(1);
    c.system.tunnel_interlocks.push(r);
    let s = solve(&c);
    assert_eq!(s.tunnel_opening.iter().filter(|row| row[0] > 0.).count(), 1);
    assert!(s.tunnel_flow.iter().flatten().all(|q| q.abs() < 1e-7));
}
#[test]
fn time_dependent_open_counts_choose_one_then_two_parallel_tunnels() {
    let mut c = case();
    c.system.units.clear();
    c.system.plants.clear();
    c.grid = vec![0., 1., 2.];
    c.prices = vec![0.; 2];
    c.system.tunnels[0].opening = 0.;
    c.system.tunnels[0].gate_positions = vec![0., 1.];
    let mut other = c.system.tunnels[0].clone();
    other.name = "Other".into();
    c.system.tunnels.push(other);
    c.system.tunnel_interlocks.push(TunnelInterlock::new(
        "Count",
        vec!["Link".into(), "Other".into()],
    ));
    for attribute in [Attribute::MinOpen, Attribute::MaxOpen] {
        c.operations.push(OperationalSeries {
            object: Object::TunnelInterlock("Count".into()),
            attribute,
            times: vec![0., 1.],
            values: vec![1., 2.],
        });
    }
    let s = solve(&c);
    assert_eq!(s.tunnel_opening.iter().filter(|row| row[0] > 0.).count(), 1);
    assert_eq!(s.tunnel_opening.iter().filter(|row| row[1] > 0.).count(), 2);
}
#[test]
fn closing_a_zero_flow_gate_still_triggers_a_wait() {
    let mut c = case();
    c.system.junctions.push(Junction {
        name: "Hub".into(),
        head: Bounds::new(0., 200.),
    });
    c.system.tunnels[0].target = "Hub".into();
    let mut other = Tunnel::new("Other", "Lake", "Hub", 10., 5.);
    other.opening = 0.;
    c.system.tunnels.push(other);
    series(&mut c, "Link", vec![0.5], vec![0.]);
    series(&mut c, "Other", vec![0.5], vec![1.]);
    wait(&mut c, 1., 0.);
    let s = solve(&c);
    assert_eq!(s.commitment[0], vec![true, false, false, true]);
    assert!(s.tunnel_flow.iter().flatten().all(|q| q.abs() < 1e-7));
}
#[test]
fn inconsistent_counts_invalid_groups_and_fractional_count_schedules_are_rejected() {
    let mut c = case();
    let mut r = TunnelInterlock::new("Count", vec!["Link".into()]);
    r.min_open = 2;
    c.system.tunnel_interlocks.push(r);
    assert!(SchedulingProblem::new(&c).is_err());
    c.system.tunnel_interlocks[0].min_open = 0;
    c.operations.push(OperationalSeries {
        object: Object::TunnelInterlock("Count".into()),
        attribute: Attribute::MinOpen,
        times: vec![0.],
        values: vec![0.5],
    });
    assert!(SchedulingProblem::new(&c).is_err());
    c.operations.clear();
    c.system.tunnel_interlocks[0].closures.push(ClosureWait {
        tunnels: vec!["Unknown".into()],
        plants: vec!["Station".into()],
        hours: 1.,
        initial_remaining: 0.,
    });
    assert!(SchedulingProblem::new(&c).is_err());
}

#[test]
fn prescribed_units_cannot_bypass_a_gate_closure_wait() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![true; 4]);
    series(&mut c, "Link", vec![0.5], vec![0.]);
    wait(&mut c, 1., 0.);
    let s = SchedulingProblem::new(&c)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap();
    assert!(!s.validation.feasible);
}
