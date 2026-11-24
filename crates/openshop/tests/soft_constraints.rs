use openshop::*;
use std::collections::BTreeMap;

fn case() -> ScheduleCase {
    ScheduleCase {
        system: HydroSystem {
            reservoirs: vec![Reservoir {
                name: "Lake".into(),
                storage: Bounds::new(1e6, 4e6),
                initial_storage: 3e6,
                water_value: 0.,
                operation: Default::default(),
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
        grid: vec![0., 0.5, 2.],
        prices: vec![100.; 2],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::from([("Unit".into(), vec![true; 2])]),
        operations: vec![],
        routing_grid: vec![],
        flow_requirements: vec![],
        soft_constraints: vec![],
    }
}
fn lower(target: SoftTarget, value: f64, allowance: f64, cost: f64) -> SoftConstraint {
    SoftConstraint {
        name: "Target".into(),
        target,
        lower: Some(PenaltyLimit {
            value,
            allowance,
            cost,
        }),
        upper: None,
    }
}
fn controls() -> Controls {
    Controls {
        head_guess: vec![],
        commitment: vec![vec![true; 2]],
        unit_flow: vec![vec![5., 15.]],
        river_opening: vec![],
        tunnel_opening: vec![],
    }
}
fn solve(case: &ScheduleCase) -> Schedule {
    SchedulingProblem::new(case)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap()
}

#[test]
fn storage_penalty_trades_revenue_against_water_without_relaxing_safety_bounds() {
    let mut c = case();
    c.soft_constraints.push(lower(
        SoftTarget::Storage {
            reservoir: "Lake".into(),
            time: 2.,
        },
        3e6,
        2e5,
        0.1,
    ));
    let expensive = solve(&c);
    assert!(expensive.validation.feasible, "{:?}", expensive.validation);
    for &q in &expensive.unit_flow[0] {
        assert!((q - 5.).abs() < 1e-4);
    }
    assert!((expensive.soft_slacks[0].lower - 36000.).abs() < 1.);
    assert!((expensive.soft_penalty_cost - 3600.).abs() < 0.1);
    c.soft_constraints[0].lower.as_mut().unwrap().cost = 0.001;
    let cheap = solve(&c);
    assert!(cheap.validation.feasible, "{:?}", cheap.validation);
    for &q in &cheap.unit_flow[0] {
        assert!((q - 20.).abs() < 1e-4);
    }
    assert!((cheap.soft_slacks[0].lower - 144000.).abs() < 1.);
    assert!((cheap.soft_penalty_cost - 144.).abs() < 0.01);
    assert!(cheap.objective > expensive.objective);
    assert!(cheap.storage[0].iter().all(|v| (1e6..=4e6).contains(v)));
}

#[test]
fn penalty_costs_affect_integer_proposals_and_nonlinear_revenue_consistently() {
    let mut c = case();
    c.commitment.clear();
    c.soft_constraints.push(lower(
        SoftTarget::Storage {
            reservoir: "Lake".into(),
            time: 2.,
        },
        3e6,
        2e5,
        0.1,
    ));
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!(s.commitment[0].iter().all(|on| !on));
    assert!(s.objective.abs() < 1e-5);
    assert!(s.soft_penalty_cost.abs() < 1e-5);
}

#[test]
fn discharge_and_power_penalties_use_duration_weighted_assessment_windows() {
    let mut c = case();
    c.soft_constraints.push(lower(
        SoftTarget::UnitDischarge {
            unit: "Unit".into(),
            start_time: 0.,
            end_time: 2.,
        },
        14.,
        2.,
        0.01,
    ));
    c.soft_constraints.push(SoftConstraint {
        name: "Plant power".into(),
        target: SoftTarget::PlantPower {
            plant: "Station".into(),
            start_time: 0.,
            end_time: 2.,
        },
        lower: Some(PenaltyLimit {
            value: 10.,
            allowance: 10.,
            cost: 2.,
        }),
        upper: None,
    });
    let s = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&controls())
        .unwrap();
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.soft_slacks[0].lower - 1.5).abs() < 1e-9);
    assert!((s.soft_slacks[0].penalty_cost - 108.).abs() < 1e-8);
    let mean_power = (0.5 * s.unit_power[0][0] + 1.5 * s.unit_power[0][1]) / 2.;
    let power_penalty = (10. - mean_power).max(0.) * 2. * 2.;
    assert!((s.soft_slacks[1].penalty_cost - power_penalty).abs() < 1e-8);
    assert!((s.soft_penalty_cost - 108. - power_penalty).abs() < 1e-8);
    let gross = 100. * (0.5 * s.unit_power[0][0] + 1.5 * s.unit_power[0][1]);
    assert!((s.objective - (gross - s.soft_penalty_cost)).abs() < 1e-8);
}

#[test]
fn measured_level_targets_use_actual_storage_curve_and_bounded_deviations() {
    let mut c = case();
    c.system.reservoirs[0].level = StorageCurve::Table(
        TableCurve::new(vec![1e6, 2e6, 3e6, 4e6], vec![60., 80., 100., 120.]).unwrap(),
    );
    c.soft_constraints.push(lower(
        SoftTarget::Level {
            reservoir: "Lake".into(),
            time: 2.,
        },
        100.,
        5.,
        10000.,
    ));
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.unit_flow[0][0] - 5.).abs() < 1e-4 && (s.unit_flow[0][1] - 5.).abs() < 1e-4);
    assert!((s.soft_slacks[0].lower - 0.72).abs() < 1e-5);
    assert!((s.soft_penalty_cost - 7200.).abs() < 0.1);
    c.soft_constraints[0].lower.as_mut().unwrap().allowance = 0.5;
    let failed = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&s.controls())
        .unwrap();
    assert!(!failed.validation.feasible);
    assert!(failed.validation.operating_limits > 0.);
}

#[test]
fn upper_targets_and_inactive_assessments_report_only_minimal_physical_slacks() {
    let mut c = case();
    c.soft_constraints.push(SoftConstraint {
        name: "Upper storage target".into(),
        target: SoftTarget::Storage {
            reservoir: "Lake".into(),
            time: 0.,
        },
        lower: None,
        upper: Some(PenaltyLimit {
            value: 2.9e6,
            allowance: 1e5,
            cost: 0.01,
        }),
    });
    c.soft_constraints.push(lower(
        SoftTarget::Storage {
            reservoir: "Lake".into(),
            time: -1.,
        },
        4e6,
        0.,
        0.,
    ));
    let s = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&controls())
        .unwrap();
    assert!(s.validation.feasible);
    assert!((s.soft_slacks[0].upper - 1e5).abs() < 1e-8);
    assert!((s.soft_penalty_cost - 1000.).abs() < 1e-8);
    assert_eq!(s.soft_slacks[1].lower, 0.);
    assert_eq!(s.soft_slacks[1].penalty_cost, 0.);
    c.soft_constraints[0].upper.as_mut().unwrap().allowance = 0.;
    assert!(
        !SchedulingProblem::new(&c)
            .unwrap()
            .reconstruct(&controls())
            .unwrap()
            .validation
            .feasible
    );
}

#[test]
fn invalid_targets_and_unpriced_deviations_are_rejected_before_optimization() {
    let base = lower(
        SoftTarget::Storage {
            reservoir: "Lake".into(),
            time: 2.,
        },
        3e6,
        2e5,
        0.1,
    );
    for change in 0..7 {
        let mut c = case();
        let mut rule = base.clone();
        match change {
            0 => rule.lower.as_mut().unwrap().cost = 0.,
            1 => rule.lower.as_mut().unwrap().allowance = -1.,
            2 => {
                rule.target = SoftTarget::Storage {
                    reservoir: "missing".into(),
                    time: 2.,
                }
            }
            3 => {
                rule.target = SoftTarget::Storage {
                    reservoir: "Lake".into(),
                    time: 1.,
                }
            }
            4 => {
                rule.target = SoftTarget::UnitPower {
                    unit: "Unit".into(),
                    start_time: -1.,
                    end_time: 0.5,
                }
            }
            5 => {
                rule.upper = Some(PenaltyLimit {
                    value: 1.,
                    allowance: 0.,
                    cost: 0.,
                })
            }
            _ => rule.lower.as_mut().unwrap().cost = f64::INFINITY,
        }
        c.soft_constraints.push(rule);
        assert!(matches!(SchedulingProblem::new(&c), Err(Error::Input(_))));
    }
}

#[test]
fn aggregate_discharge_targets_sum_all_units_before_time_averaging() {
    let mut c = case();
    let mut second = c.system.units[0].clone();
    second.name = "Second".into();
    c.system.units.push(second);
    c.commitment.insert("Second".into(), vec![true; 2]);
    c.soft_constraints.push(lower(
        SoftTarget::PlantDischarge {
            plant: "Station".into(),
            start_time: 0.,
            end_time: 2.,
        },
        24.,
        2.,
        0.01,
    ));
    let mut input = controls();
    input.commitment.push(vec![true; 2]);
    input.unit_flow.push(vec![10.; 2]);
    let s = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&input)
        .unwrap();
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.soft_slacks[0].lower - 1.5).abs() < 1e-9);
    assert!((s.soft_penalty_cost - 108.).abs() < 1e-8);
}

#[test]
fn zero_allowance_is_a_hard_target_and_large_allowance_cannot_relax_reservoir_safety() {
    let mut c = case();
    c.soft_constraints.push(lower(
        SoftTarget::Storage {
            reservoir: "Lake".into(),
            time: 2.,
        },
        2.94e6,
        0.,
        0.,
    ));
    let s = solve(&c);
    assert!(s.validation.feasible, "{:?}", s.validation);
    assert!((s.storage[0][2] - 2.94e6).abs() < 0.1);
    assert!(s.soft_slacks[0].lower < 0.1 && s.soft_penalty_cost == 0.);
    c.system.reservoirs[0].storage.min = 2.98e6;
    c.soft_constraints[0].lower = Some(PenaltyLimit {
        value: 3e6,
        allowance: 2e5,
        cost: 0.001,
    });
    let failed = SchedulingProblem::new(&c)
        .unwrap()
        .reconstruct(&s.controls())
        .unwrap();
    assert!(!failed.validation.feasible);
    assert!(failed.storage[0][2] < c.system.reservoirs[0].storage.min);
}
