use openshop::*;
use std::collections::BTreeMap;

fn case() -> ScheduleCase {
    ScheduleCase {
        operations: vec![],
        system: HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(1e6, 4e6),
                initial_storage: 3e6,
                operation: Default::default(),
                water_value: 0.001,
                level: StorageCurve::Polynomial {
                    reference_volume: 0.,
                    elevation: 40.,
                    slope: 20e-6,
                    curvature: 0.,
                },
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 20.,
            }],
            plants: vec![Plant::new("Station", "Lake", "Sea")],
            units: vec![Unit::new(
                "Unit",
                "Station",
                Bounds::new(5., 20.),
                Bounds::new(1., 25.),
                Bounds::new(35., 105.),
                Efficiency::constant(0.9),
            )],
            ..Default::default()
        },
        grid: vec![0., 1., 2.],
        prices: vec![100., 100.],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
    }
}
fn series(
    object: Object,
    attribute: Attribute,
    times: Vec<f64>,
    values: Vec<f64>,
) -> OperationalSeries {
    OperationalSeries {
        object,
        attribute,
        times,
        values,
    }
}
fn solve(c: &ScheduleCase) -> Schedule {
    let result = SchedulingProblem::new(c)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    result
}
#[test]
fn requested_production_and_maintenance_control_commitment() {
    let mut c = case();
    c.grid = vec![0., 1., 2., 3.];
    c.prices = vec![100.; 3];
    c.operations = vec![
        series(
            Object::Unit("Unit".into()),
            Attribute::Power,
            vec![0., 1., 2.],
            vec![6., 0., 8.],
        ),
        series(
            Object::Unit("Unit".into()),
            Attribute::Maintenance,
            vec![1., 2.],
            vec![1., 0.],
        ),
    ];
    let result = solve(&c);
    assert_eq!(result.commitment[0], vec![true, false, true]);
    for (&p, wanted) in result.unit_power[0].iter().zip([6., 0., 8.]) {
        assert!((p - wanted).abs() < 1e-5);
    }
}
#[test]
fn changing_power_limits_are_enforced_in_nonlinear_dispatch() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![true; 2]);
    c.operations = vec![series(
        Object::Unit("Unit".into()),
        Attribute::PowerMax,
        vec![0., 1.],
        vec![8., 10.],
    )];
    let result = solve(&c);
    assert!((result.unit_power[0][0] - 8.).abs() < 1e-5);
    assert!((result.unit_power[0][1] - 10.).abs() < 1e-5);
}
#[test]
fn ramp_rates_are_integrated_between_observation_midpoints() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![true; 2]);
    c.system.units[0].operation.initial_on = true;
    c.system.units[0].operation.initial_flow = Some(5.);
    c.operations = vec![series(
        Object::Unit("Unit".into()),
        Attribute::FlowRampUp,
        vec![-0.5, 0., 1.],
        vec![2., 4., 6.],
    )];
    let result = solve(&c);
    assert!((result.unit_flow[0][0] - 8.).abs() < 1e-5);
    assert!((result.unit_flow[0][1] - 13.).abs() < 1e-5);
}
#[test]
fn omitted_historical_output_does_not_invent_a_previous_average() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![true; 2]);
    c.system.units[0].operation.flow_ramp.up = Some(0.);
    let result = solve(&c);
    assert!((result.unit_flow[0][0] - 20.).abs() < 1e-5);
}
#[test]
fn changing_storage_bounds_and_inflow_change_the_schedule() {
    let mut c = case();
    c.operations = vec![
        series(
            Object::Reservoir("Lake".into()),
            Attribute::StorageMin,
            vec![1.],
            vec![2.99e6],
        ),
        series(
            Object::Reservoir("Lake".into()),
            Attribute::Inflow,
            vec![0., 1.],
            vec![0., 20.],
        ),
    ];
    let result = solve(&c);
    assert_eq!(result.commitment[0], vec![false, true]);
    assert!((result.unit_flow[0][1] - 20.).abs() < 1e-5);
    assert!(result.storage[0].iter().all(|v| *v >= 2.99e6 - 0.001));
}
#[test]
fn tabulated_vertex_levels_enforce_level_ramps() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![true; 2]);
    c.system.reservoirs[0].level = StorageCurve::Table(
        TableCurve::new(vec![1e6, 2.8e6, 3.2e6, 4e6], vec![40., 96., 104., 120.]).unwrap(),
    );
    c.system.reservoirs[0].operation.level_ramp.down = Some(0.5);
    let result = solve(&c);
    for &q in &result.unit_flow[0] {
        assert!((q - 0.5 / (3600. * 20e-6)).abs() < 1e-4);
    }
}
#[test]
fn scheduled_gate_closure_disconnects_heads() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![false; 2]);
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 120.,
        slope: 0.,
        curvature: 0.,
    };
    c.system.tunnels.push(Tunnel {
        name: "Gate".into(),
        source: "Lake".into(),
        target: "Sea".into(),
        capacity: 50.,
        loss_factor: 1.,
        opening: 1.,
        gate_loss_factor_curve: None,
    });
    c.operations = vec![series(
        Object::Tunnel("Gate".into()),
        Attribute::Opening,
        vec![0., 1.],
        vec![0.25, 0.],
    )];
    let result = solve(&c);
    assert!((result.tunnel_flow[0][0] - 5.).abs() < 1e-6);
    assert_eq!(result.tunnel_flow[0][1], 0.);
    assert_eq!(result.head[0][1] - result.head[1][1], 100.);
}
#[test]
fn river_inflows_follow_their_original_release_intervals() {
    let mut c = case();
    c.system.units.clear();
    c.system.plants.clear();
    c.system.rivers.push(River {
        law: Default::default(),
        gate: Bounds::new(0., 1.),
        name: "River".into(),
        source: "Natural".into(),
        target: "Sea".into(),
        capacity: 10.,
        delay: 0.5,
        allocation: 1.,
        inflow: 0.,
        water_value: 0.,
        history_grid: vec![],
        history_release: vec![],
    });
    c.operations = vec![
        series(
            Object::River("River".into()),
            Attribute::Inflow,
            vec![0., 1.],
            vec![3., 7.],
        ),
        series(
            Object::River("River".into()),
            Attribute::Capacity,
            vec![0., 1.],
            vec![5., 10.],
        ),
    ];
    let result = solve(&c);
    assert_eq!(result.river_release[0], vec![3., 7.]);
    assert!((result.river_arrival_volume[0][0] - 5400.).abs() < 1e-7);
    assert!((result.river_storage[0][2] - 12600.).abs() < 1e-7);
}
#[test]
fn inconsistent_states_and_unrepresented_operating_changes_are_rejected() {
    let mut c = case();
    c.operations = vec![series(
        Object::Unit("Unit".into()),
        Attribute::PowerMax,
        vec![0.5],
        vec![8.],
    )];
    assert!(SchedulingProblem::new(&c).is_err());
    c.operations = vec![
        series(
            Object::Unit("Unit".into()),
            Attribute::Maintenance,
            vec![0.],
            vec![1.],
        ),
        series(
            Object::Unit("Unit".into()),
            Attribute::Power,
            vec![0.],
            vec![8.],
        ),
    ];
    assert!(SchedulingProblem::new(&c).is_err());
}

#[test]
fn plant_power_and_flow_can_use_different_historical_intervals() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![true; 2]);
    let unit = &mut c.system.units[0];
    unit.operation.initial_on = true;
    unit.operation.initial_flow = Some(5.);
    unit.operation.initial_power = Some(5.);
    unit.operation.initial_interval_hours = 0.5;
    let plant = &mut c.system.plants[0];
    plant.operation.initial_power = Some(5.);
    plant.operation.initial_interval_hours = 2.;
    plant.operation.power_ramp.up = Some(2.);
    plant.operation.flow_ramp.up = Some(2.);
    let result = solve(&c);
    assert!((result.unit_flow[0][0] - 6.5).abs() < 1e-5);
    c.system.plants[0].operation.flow_ramp.up = None;
    let result = solve(&c);
    assert!((result.unit_power[0][0] - 8.).abs() < 1e-5);
}
#[test]
fn shutdown_uses_the_previous_minimum_power_at_a_bound_change() {
    let mut c = case();
    c.commitment.insert("Unit".into(), vec![false; 2]);
    let unit = &mut c.system.units[0];
    unit.operation.initial_on = true;
    unit.operation.initial_power = Some(5.);
    unit.operation.power_ramp.down = Some(0.);
    c.operations = vec![series(
        Object::Unit("Unit".into()),
        Attribute::PowerMin,
        vec![0.],
        vec![10.],
    )];
    let result = SchedulingProblem::new(&c)
        .unwrap()
        .solve(&SolveOptions::default())
        .unwrap();
    assert!(!result.validation.feasible);
}
#[test]
fn future_operating_bounds_cannot_extend_a_static_physical_domain() {
    let mut c = case();
    c.operations = vec![series(
        Object::Unit("Unit".into()),
        Attribute::PowerMax,
        vec![100.],
        vec![100.],
    )];
    assert!(SchedulingProblem::new(&c).is_err());
}
