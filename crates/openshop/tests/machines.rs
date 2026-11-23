use openshop::*;
use std::collections::BTreeMap;

fn pumped_case() -> ScheduleCase {
    let reservoirs = vec![
        Reservoir {
            operation: Default::default(),
            name: "Upper".into(),
            storage: Bounds::new(1e6, 20e6),
            initial_storage: 5e6,
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 299.5,
                slope: 0.1e-6,
                curvature: 0.,
            },
            water_value: 0.018,
        },
        Reservoir {
            operation: Default::default(),
            name: "Lower".into(),
            storage: Bounds::new(1e6, 20e6),
            initial_storage: 15e6,
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 48.5,
                slope: 0.1e-6,
                curvature: 0.,
            },
            water_value: 0.001,
        },
    ];
    let mut turbine = Unit::new(
        "Turbine",
        "Generating",
        Bounds::new(20., 20.),
        Bounds::new(40., 60.),
        Bounds::new(200., 300.),
        Efficiency::constant(0.91),
    );
    turbine.operation.startup_cost = 10.;
    let mut pump = Unit::new(
        "Pump",
        "Pumping",
        Bounds::new(20., 20.),
        Bounds::new(50., 65.),
        Bounds::new(200., 300.),
        Efficiency::constant(0.88),
    );
    pump.mode = UnitMode::Pump;
    pump.motor_efficiency = 0.96;
    pump.operation.startup_cost = 10.;
    ScheduleCase {
        flow_requirements: vec![],
        routing_grid: vec![],
        operations: Default::default(),
        system: HydroSystem {
            reservoirs,
            plants: vec![
                Plant::new("Generating", "Upper", "Lower"),
                Plant::new("Pumping", "Lower", "Upper"),
            ],
            units: vec![turbine, pump],
            commitment_exclusions: vec![CommitmentExclusion {
                name: "Reversible".into(),
                groups: vec![vec!["Turbine".into()], vec!["Pump".into()]],
            }],
            ..Default::default()
        },
        grid: vec![0., 1., 2., 3., 4.],
        prices: vec![10., 15., 110., 120.],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::from([
            ("Turbine".into(), vec![false, false, true, true]),
            ("Pump".into(), vec![true, true, false, false]),
        ]),
    }
}
#[test]
fn pumped_cycle_conserves_water_and_pays_for_electrical_input() {
    let case = pumped_case();
    let result = SchedulingProblem::new(&case)
        .unwrap()
        .solve(&SolveOptions::default())
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    let expected = [5e6, 5.072e6, 5.144e6, 5.072e6, 5e6];
    for (t, &volume) in expected.iter().enumerate() {
        assert!(
            (result.storage[0][t] - volume).abs() < 1e-4,
            "{:?}",
            result.storage
        );
        assert!((result.storage[0][t] + result.storage[1][t] - 20e6).abs() < 1e-5);
    }
    let mut value = -20.;
    let mut produced = 0.;
    let mut consumed = 0.;
    for t in 0..4 {
        let head = 251.
            + 0.1e-6
                * (result.storage[0][t] + result.storage[0][t + 1]
                    - result.storage[1][t]
                    - result.storage[1][t + 1])
                / 2.;
        if t < 2 {
            let p = 0.00981 * 20. * head / (0.88 * 0.96);
            assert!((result.unit_power[1][t] - p).abs() < 1e-7);
            value -= case.prices[t] * p;
            consumed += p;
        } else {
            let p = 0.00981 * 20. * head * 0.91;
            assert!((result.unit_power[0][t] - p).abs() < 1e-7);
            value += case.prices[t] * p;
            produced += p;
        }
    }
    assert!(produced < consumed);
    let net = result.net_power(&case).unwrap();
    assert!(net[0] < 0. && net[2] > 0.);
    assert!((result.objective - value).abs() < 1e-5);
}
#[test]
fn commitment_search_can_choose_pumping_then_generation() {
    let mut case = pumped_case();
    case.commitment.clear();
    let result = SchedulingProblem::new(&case)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 4.,
            ..Default::default()
        })
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    assert_eq!(result.commitment[0], vec![false, false, true, true]);
    assert_eq!(result.commitment[1], vec![true, true, false, false]);
}
#[test]
fn fixed_speed_pump_flow_follows_storage_and_shared_tailwater() {
    let mut case = pumped_case();
    case.system.units.remove(0);
    case.system.plants.remove(0);
    case.system.commitment_exclusions.clear();
    case.system.reservoirs[0].initial_storage = 10e6;
    case.system.reservoirs[1].initial_storage = 10e6;
    case.system.reservoirs[0].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: 90.,
        slope: 1e-6,
        curvature: 0.,
    };
    case.system.reservoirs[1].level = StorageCurve::Polynomial {
        reference_volume: 0.,
        elevation: -10.,
        slope: 1e-6,
        curvature: 0.,
    };
    case.system.plants[0].tailwater_curve =
        Some(TableCurve::new(vec![0., 20.], vec![0., 1.]).unwrap());
    let pump = &mut case.system.units[0];
    pump.flow = Bounds::new(9., 11.);
    pump.power = Bounds::new(1., 20.);
    pump.head = Bounds::new(90., 110.);
    pump.motor_efficiency = 0.95;
    pump.pump_characteristic = Some(
        PumpCharacteristic::new(
            vec![90., 100., 110.],
            vec![10.2, 10., 9.8],
            vec![0.845, 0.85, 0.855],
            Extrapolation::Error,
        )
        .unwrap(),
    );
    case.grid = vec![0., 0.5, 1.25, 2.];
    case.prices = vec![-50., 20., -50.];
    case.commitment = BTreeMap::from([("Pump".into(), vec![true, false, true])]);
    let result = SchedulingProblem::new(&case)
        .unwrap()
        .solve(&SolveOptions::default())
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    let mut upper = 10e6;
    let mut lower = 10e6;
    for t in 0..3 {
        let dt = case.grid[t + 1] - case.grid[t];
        let h0 = 100. + (upper - lower) * 1e-6;
        let k = 0.05 + 0.0036 * dt;
        let q = if t == 1 {
            0.
        } else {
            (12. - 0.02 * h0) / (1. + 0.02 * k)
        };
        let h = h0 + k * q;
        let eta = 0.8 + 0.0005 * h;
        let p = 0.00981 * q * h / (eta * 0.95);
        assert!((result.unit_flow[0][t] - q).abs() < 1e-6);
        assert!((result.unit_power[0][t] - p).abs() < 1e-6);
        upper += 3600. * dt * q;
        lower -= 3600. * dt * q;
        assert!((result.storage[0][t + 1] - upper).abs() < 0.01);
        assert!((result.storage[1][t + 1] - lower).abs() < 0.01);
    }
}
#[test]
fn plant_caps_apply_to_parallel_units_and_share_the_outlet_level() {
    let mut case = pumped_case();
    case.system.units.remove(1);
    case.system.plants.remove(1);
    case.system.commitment_exclusions.clear();
    case.grid = vec![0., 1.];
    case.prices = vec![100.];
    case.commitment.clear();
    case.system.plants[0].flow_max = Some(25.);
    case.system.plants[0].power_max = Some(55.);
    case.system.plants[0].tailwater_curve =
        Some(TableCurve::new(vec![0., 40.], vec![0., 4.]).unwrap());
    let unit = &mut case.system.units[0];
    unit.flow = Bounds::new(5., 20.);
    unit.power = Bounds::new(1., 60.);
    let mut parallel = unit.clone();
    parallel.name = "Parallel".into();
    case.system.units.push(parallel);
    let result = SchedulingProblem::new(&case)
        .unwrap()
        .solve(&SolveOptions {
            time_limit: 3.,
            ..Default::default()
        })
        .unwrap();
    assert!(result.validation.feasible, "{:?}", result.validation);
    assert!((result.plant_power[0][0] - 55.).abs() < 1e-4);
    assert!(result.plant_flow[0][0] <= 25. + 1e-6);
    assert!((result.unit_head[0][0] - result.unit_head[1][0]).abs() < 1e-8);
}
#[test]
fn invalid_modes_and_conflicting_prescribed_commitments_are_rejected() {
    let mut case = pumped_case();
    case.commitment.insert("Turbine".into(), vec![true; 4]);
    assert!(SchedulingProblem::new(&case).is_err());
    let mut case = pumped_case();
    case.system.units[1].generator_efficiency_curve =
        Some(TableCurve::new(vec![0., 65.], vec![0.9, 0.95]).unwrap());
    assert!(SchedulingProblem::new(&case).is_err());
}
