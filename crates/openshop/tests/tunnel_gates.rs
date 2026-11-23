use openshop::*;
use std::collections::BTreeMap;
fn lake(name: &str, head: f64, volume: f64, value: f64) -> Reservoir {
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
        water_value: value,
        operation: Default::default(),
    }
}
fn case() -> ScheduleCase {
    let mut e = Tunnel::new("Link", "Upper", "Lower", 10., 5.);
    e.opening = 0.;
    e.gate_positions = vec![0., 0.25, 1.];
    e.gate_adjustment_cost = 10.;
    ScheduleCase {
        flow_requirements: vec![],
        routing_grid: vec![],
        operations: vec![],
        system: HydroSystem {
            reservoirs: vec![lake("Upper", 100., 5e5, 0.), lake("Lower", 20., 1e5, 0.01)],
            tunnels: vec![e],
            ..Default::default()
        },
        grid: vec![0., 1., 2.],
        prices: vec![0.; 2],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
    }
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
fn discrete_gates_are_optimized_even_when_no_unit_commitments_are_free() {
    let c = case();
    let s = solve(&c);
    assert_eq!(s.tunnel_opening[0], vec![1., 1.]);
    assert!(s.tunnel_flow[0].iter().all(|q| (*q - 4.).abs() < 1e-6));
    assert!((s.gate_adjustment_cost - 10.).abs() < 1e-8);
    assert!((s.objective - 278.).abs() < 1e-4);
}
#[test]
fn adjustment_cost_can_make_remaining_closed_optimal() {
    let mut c = case();
    c.system.tunnels[0].gate_adjustment_cost = 1000.;
    let s = solve(&c);
    assert_eq!(s.tunnel_opening[0], vec![0., 0.]);
    assert!(s.objective.abs() < 1e-4);
}
#[test]
fn opening_schedule_is_hard_after_its_first_knot_and_free_before_it() {
    let mut c = case();
    c.grid = vec![0., 1., 2., 3.];
    c.prices = vec![0.; 3];
    c.operations.push(OperationalSeries {
        object: Object::Tunnel("Link".into()),
        attribute: Attribute::Opening,
        times: vec![1., 2.],
        values: vec![1., 0.],
    });
    let s = solve(&c);
    assert_eq!(s.tunnel_opening[0], vec![1., 1., 0.]);
    assert!((s.gate_adjustment_cost - 20.).abs() < 1e-8);
}
#[test]
fn prescribed_positions_charge_initial_and_later_absolute_changes() {
    let mut c = case();
    let e = &mut c.system.tunnels[0];
    e.gate_positions.clear();
    e.opening = 0.5;
    e.gate_adjustment_cost = 10.;
    c.operations.push(OperationalSeries {
        object: Object::Tunnel("Link".into()),
        attribute: Attribute::Opening,
        times: vec![0., 1.],
        values: vec![1., 0.25],
    });
    let s = solve(&c);
    assert_eq!(s.tunnel_opening[0], vec![1., 0.25]);
    assert!((s.gate_adjustment_cost - 12.5).abs() < 1e-8);
    assert!((s.tunnel_flow[0][0] - 4.).abs() < 1e-6);
    assert!((s.tunnel_flow[0][1] - 2.).abs() < 1e-6);
    assert!((s.objective - 203.5).abs() < 1e-4);
}
#[test]
fn gate_resistance_curves_remain_physical_during_position_optimization() {
    let mut c = case();
    let e = &mut c.system.tunnels[0];
    e.gate_positions = vec![0., 0.5, 1.];
    e.gate_loss_factor_curve = Some(TableCurve::new(vec![0., 0.5, 1.], vec![15., 5., 0.]).unwrap());
    let s = solve(&c);
    assert_eq!(s.tunnel_opening[0], vec![1., 1.]);
    assert!((s.objective - 278.).abs() < 1e-4);
}
#[test]
fn undeclared_schedules_and_invalid_initial_gate_positions_are_rejected() {
    let mut c = case();
    c.operations.push(OperationalSeries {
        object: Object::Tunnel("Link".into()),
        attribute: Attribute::Opening,
        times: vec![10.],
        values: vec![0.5],
    });
    assert!(SchedulingProblem::new(&c).is_err());
    c.operations.clear();
    c.system.tunnels[0].initial_opening = Some(0.5);
    assert!(SchedulingProblem::new(&c).is_err());
}

#[test]
fn a_small_opening_retains_pressure_accuracy_while_storage_changes() {
    let mut c = case();
    c.grid = vec![0., 1.];
    c.prices = vec![0.];
    let e = &mut c.system.tunnels[0];
    e.gate_positions.clear();
    e.opening = 1e-8;
    e.gate_adjustment_cost = 0.;
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 5e5,
        elevation: 100.,
        slope: 1e-4,
        curvature: 0.,
    };
    let s = solve(&c);
    let expected = 160. / (0.18 + (0.18_f64.powi(2) + 4. * 5e8 * 80.).sqrt());
    assert!((s.tunnel_flow[0][0] - expected).abs() < 1e-11);
    assert!(s.validation.pressure_balance < 1e-6);
}
