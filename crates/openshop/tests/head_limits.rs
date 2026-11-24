use openshop::*;
use std::collections::BTreeMap;
fn case() -> ScheduleCase {
    let mut unit = Unit::new(
        "Unit",
        "Station",
        Bounds::new(1., 20.),
        Bounds::new(0., 30.),
        Bounds::new(1., 150.),
        Efficiency::constant(0.9),
    );
    unit.loss_factor = 0.01;
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
                    elevation: 90.,
                    slope: 1e-4,
                    curvature: 0.,
                },
                water_value: 0.,
                operation: Default::default(),
            }],
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 20.,
            }],
            plants: vec![Plant::new("Station", "Lake", "Sea")],
            units: vec![unit],
            ..Default::default()
        },
        grid: vec![0., 1.],
        prices: vec![100.],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::from([("Unit".into(), vec![true])]),
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
fn running_gross_head_uses_surface_levels_before_machine_losses() {
    let mut c = case();
    c.system.units[0].gross_head_min = Some(118.);
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_flow[0][0] - 100. / 9.).abs() < 1e-4);
    assert!((s.head[0][0] - s.head[1][0] - 118.).abs() < 1e-5);
    assert!(s.unit_head[0][0] < 117.);
}
#[test]
fn startup_head_uses_the_actual_starting_surface_instead_of_the_interval_mean() {
    for table in [false, true] {
        let mut c = case();
        if table {
            c.system.reservoirs[0].level = StorageCurve::Table(
                TableCurve::new(vec![0., 5e5, 1e6], vec![90., 140., 190.]).unwrap(),
            );
        }
        c.system.units[0].startup_gross_head_min = Some(120.);
        let s = solve(&c);
        assert!(s.validation.feasible, "{:?}", s.validation);
        assert!(s.unit_flow[0][0] > 19.99);
        assert!(s.head[0][0] - s.head[1][0] < 120.);
        c.system.units[0].startup_gross_head_min = Some(120.1);
        assert!(!solve(&c).validation.feasible);
    }
}
#[test]
fn explicit_gross_head_references_can_differ_from_machine_endpoints() {
    let mut c = case();
    c.system.boundaries.push(Boundary {
        name: "UpperGauge".into(),
        head: 150.,
    });
    c.system.plants[0].gross_head_source = Some("UpperGauge".into());
    c.system.units[0].gross_head_min = Some(129.);
    c.system.units[0].gross_head_max = Some(131.);
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!(s.unit_flow[0][0] > 19.99);
    c.system.plants[0].gross_head_source = Some("Unknown".into());
    assert!(SchedulingProblem::new(&c).is_err());
}
#[test]
fn minimum_efficiency_restricts_dispatch_instead_of_clipping_the_curve() {
    let mut c = case();
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 120.,
        slope: 0.,
        curvature: 0.,
    };
    let u = &mut c.system.units[0];
    u.loss_factor = 0.;
    u.minimum_efficiency = 0.88;
    u.efficiency = Efficiency::Quadratic {
        peak: 0.9,
        best_flow: 10.,
        flow_curvature: 0.0008,
        best_head: 100.,
        head_curvature: 0.,
    };
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_flow[0][0] - 15.).abs() < 1e-4);
    assert!(
        (c.system.units[0]
            .efficiency
            .value(s.unit_flow[0][0], s.unit_head[0][0])
            .unwrap()
            - 0.88)
            .abs()
            < 1e-7
    );
}
#[test]
fn changing_gross_head_limits_control_commitment() {
    let mut c = case();
    c.commitment.clear();
    c.grid = vec![0., 1., 2.];
    c.prices = vec![100.; 2];
    c.operations.push(OperationalSeries {
        object: Object::Unit("Unit".into()),
        attribute: Attribute::GrossHeadMin,
        times: vec![0., 1.],
        values: vec![121., 100.],
    });
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert_eq!(s.commitment[0], vec![false, true]);
}

#[test]
fn head_dependent_flow_limits_can_exclude_part_of_a_units_declared_head_range() {
    let mut c = case();
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 140.,
        slope: 0.,
        curvature: 0.,
    };
    let u = &mut c.system.units[0];
    u.flow = Bounds::new(5., 20.);
    u.efficiency = Efficiency::Table(
        TurbineTable::new(
            vec![100., 120.],
            vec![TableCurve::new(vec![5., 20.], vec![0.8, 0.9]).unwrap(); 2],
            vec![Bounds::new(20., 20.), Bounds::new(5., 20.)],
            TurbineOptions {
                head_extrapolation: Extrapolation::Linear,
                ..Default::default()
            },
        )
        .unwrap(),
    );
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_flow[0][0] - 20.).abs() < 1e-4);
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 100.,
        slope: 0.,
        curvature: 0.,
    };
    assert!(!solve(&c).validation.feasible);
}
#[test]
fn efficiency_bounds_exclude_nonphysical_extrapolated_points_without_rejecting_the_table() {
    let mut c = case();
    c.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 420.,
        slope: 0.,
        curvature: 0.,
    };
    let u = &mut c.system.units[0];
    u.flow = Bounds::new(1., 10.);
    u.head = Bounds::new(100., 400.);
    u.loss_factor = 0.;
    u.efficiency = Efficiency::Table(
        TurbineTable::new(
            vec![100., 200.],
            vec![
                TableCurve::new(vec![0., 5., 10.], vec![0.7, 0.9, 0.7]).unwrap(),
                TableCurve::new(vec![0., 5., 10.], vec![0.6, 0.95, 0.6]).unwrap(),
            ],
            vec![Bounds::new(0., 10.); 2],
            TurbineOptions {
                interpolation: DischargeInterpolation::Pchip,
                head_extrapolation: Extrapolation::Linear,
                ..Default::default()
            },
        )
        .unwrap(),
    );
    c.operations.push(OperationalSeries {
        object: Object::Unit("Unit".into()),
        attribute: Attribute::Discharge,
        times: vec![0.],
        values: vec![5.],
    });
    assert!(!solve(&c).validation.feasible);
    c.operations[0].values[0] = 10.;
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_power[0][0] - 15.696).abs() < 1e-4);
}
