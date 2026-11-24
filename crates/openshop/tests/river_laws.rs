use openshop::*;
use std::collections::BTreeMap;
fn reservoir(name: &str, head: f64, volume: f64, slope: f64) -> Reservoir {
    Reservoir {
        name: name.into(),
        storage: Bounds::new(0., 2e6),
        initial_storage: volume,
        level: StorageCurve::Polynomial {
            reference_volume: volume,
            elevation: head,
            slope,
            curvature: 0.,
        },
        water_value: 0.,
        operation: Default::default(),
    }
}
fn river(source: &str, target: &str, law: RiverLaw) -> River {
    River {
        operation: Default::default(),
        name: "Reach".into(),
        source: source.into(),
        target: target.into(),
        law,
        gate: Bounds::new(0., 1.),
        capacity: 100.,
        transport: RiverTransport::Fixed { delay: 0. },
        allocation: 1.,
        inflow: 0.,
        water_value: 0.,
        history_grid: vec![],
        history_release: vec![],
    }
}
fn case(reservoirs: Vec<Reservoir>, r: River) -> ScheduleCase {
    ScheduleCase {
        system: HydroSystem {
            reservoirs,
            boundaries: vec![Boundary {
                name: "Sea".into(),
                head: 0.,
            }],
            rivers: vec![r],
            ..Default::default()
        },
        grid: vec![0., 0.5, 1.5],
        prices: vec![0.; 2],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
        soft_constraints: vec![],
        flow_requirements: vec![],
        routing_grid: vec![],
        operations: vec![],
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
    assert!(
        s.validation.feasible,
        "{:?}, status {:?}",
        s.validation, s.solver_status
    );
    s
}
fn gate(c: &mut ScheduleCase, values: Vec<f64>) {
    for attribute in [Attribute::GateMin, Attribute::GateMax] {
        c.operations.push(OperationalSeries {
            object: Object::River("Reach".into()),
            attribute,
            times: vec![0., 0.5],
            values: values.clone(),
        });
    }
}
#[test]
fn outlet_laws_match_heads_and_gate_schedules_without_losing_natural_inflow() {
    let laws = [
        RiverLaw::Controlled,
        RiverLaw::Orifice {
            coefficient: 2.,
            crest: 90.,
            allow_dry: false,
        },
        RiverLaw::Weir {
            coefficient: 0.1,
            crest: 90.,
            allow_dry: false,
        },
        RiverLaw::Level {
            curve: TableCurve::new(vec![80., 100., 120.], vec![1., 8., 12.]).unwrap(),
        },
    ];
    for law in laws {
        let weir = matches!(law, RiverLaw::Weir { .. });
        let mut c = case(
            vec![reservoir("Lake", 105., 1e6, 1e-5)],
            river("Lake", "Sea", law),
        );
        c.system.rivers[0].inflow = 3.;
        if !weir {
            gate(&mut c, vec![0., 0.6]);
        }
        let s = solve(&c);
        for t in 0..2 {
            let g = if weir { 1. } else { [0., 0.6][t] };
            let q = c.system.rivers[0]
                .law
                .outlet_flow(s.head[0][t], g, 100., 3.)
                .unwrap();
            assert!((s.river_withdrawal[0][t] - q).abs() < 1e-6);
            assert!((s.river_release[0][t] - q - 3.).abs() < 1e-6);
            assert!(
                (s.river_arrival_volume[0][t]
                    - 3600. * (c.grid[t + 1] - c.grid[t]) * s.river_release[0][t])
                    .abs()
                    < 1e-5
            );
            assert!(
                (s.storage[0][t + 1] - s.storage[0][t] + 3600. * (c.grid[t + 1] - c.grid[t]) * q)
                    .abs()
                    < 1e-4
            );
        }
    }
}
#[test]
fn dry_outlets_have_zero_discharge_without_singular_root_evaluation() {
    for law in [
        RiverLaw::Orifice {
            coefficient: 2.,
            crest: 120.,
            allow_dry: true,
        },
        RiverLaw::Weir {
            coefficient: 0.1,
            crest: 120.,
            allow_dry: true,
        },
    ] {
        let c = case(
            vec![reservoir("Lake", 100., 1e6, 1e-5)],
            river("Lake", "Sea", law),
        );
        let s = solve(&c);
        assert!(s.river_withdrawal[0].iter().all(|q| q.abs() < 1e-7));
    }
}
#[test]
fn channels_follow_signed_level_differences_and_conserve_water_in_both_directions() {
    for direction in [-1., 0., 1.] {
        for reference in [ReferenceLevel::Upstream, ReferenceLevel::Downstream] {
            let response = LevelFlowCurve::new(
                vec![80., 120.],
                vec![
                    TableCurve::new(vec![-20., 0., 20.], vec![-8., 0., 8.]).unwrap(),
                    TableCurve::new(vec![-20., 0., 20.], vec![-16., 0., 16.]).unwrap(),
                ],
            )
            .unwrap();
            let c = case(
                vec![
                    reservoir("A", 100. + 10. * direction, 1e6, 1e-5),
                    reservoir("B", 100., 1e6, 2e-5),
                ],
                river(
                    "A",
                    "B",
                    RiverLaw::Delta {
                        curve: response.clone(),
                        reference,
                    },
                ),
            );
            let s = solve(&c);
            for t in 0..2 {
                let (a, b) = (s.head[0][t], s.head[1][t]);
                let h = match reference {
                    ReferenceLevel::Upstream => a,
                    ReferenceLevel::Downstream => b,
                };
                assert!(
                    (s.river_withdrawal[0][t] - response.value(a - b, h).unwrap()).abs() < 1e-6
                );
                assert!((s.storage[0][t + 1] + s.storage[1][t + 1] - 2e6).abs() < 1e-4);
                assert!(s.river_storage[0][t].abs() < 1e-6);
                assert!(s.river_withdrawal[0][t] * direction >= -1e-7);
            }
        }
    }
}
#[test]
fn passive_channel_data_and_incompatible_operating_inputs_are_rejected() {
    assert!(
        LevelFlowCurve::new(
            vec![100.],
            vec![TableCurve::new(vec![-1., 0., 1.], vec![1., 0., 2.]).unwrap()]
        )
        .is_err()
    );
    let response = LevelFlowCurve::new(
        vec![100.],
        vec![TableCurve::new(vec![-1., 0., 1.], vec![-1., 0., 1.]).unwrap()],
    )
    .unwrap();
    let mut c = case(
        vec![
            reservoir("A", 100., 1e6, 1e-5),
            reservoir("B", 90., 1e6, 1e-5),
        ],
        river(
            "A",
            "B",
            RiverLaw::Delta {
                curve: response,
                reference: ReferenceLevel::Upstream,
            },
        ),
    );
    c.system.rivers[0].transport = RiverTransport::Fixed { delay: 0.1 };
    assert!(SchedulingProblem::new(&c).is_err());
    c.system.rivers[0].transport = RiverTransport::default();
    c.operations.push(OperationalSeries {
        object: Object::River("Reach".into()),
        attribute: Attribute::Inflow,
        times: vec![0.],
        values: vec![0.],
    });
    assert!(SchedulingProblem::new(&c).is_err());
}

#[test]
fn minimum_release_includes_natural_water_before_reservoir_withdrawal() {
    let mut c = case(
        vec![reservoir("Lake", 100., 1e6, 0.)],
        river("Lake", "Sea", RiverLaw::Controlled),
    );
    c.system.reservoirs[0].water_value = 0.001;
    c.system.rivers[0].inflow = 1.;
    c.operations.push(OperationalSeries {
        object: Object::River("Reach".into()),
        attribute: Attribute::MinRelease,
        times: vec![0., 0.5],
        values: vec![2., 5.],
    });
    let s = solve(&c);
    assert!((s.river_release[0][0] - 2.).abs() < 1e-5);
    assert!((s.river_release[0][1] - 5.).abs() < 1e-5);
    assert!((s.river_withdrawal[0][0] - 1.).abs() < 1e-5);
    assert!((s.river_withdrawal[0][1] - 4.).abs() < 1e-5);
}
#[test]
fn a_river_junction_release_minimum_is_an_interval_average() {
    let mut c = case(
        vec![reservoir("Lake", 100., 1e6, 0.)],
        river("Lake", "Join", RiverLaw::Controlled),
    );
    c.grid = vec![0., 1., 2.];
    c.prices = vec![0.; 2];
    c.system.reservoirs[0].water_value = 0.001;
    c.system.rivers[0].name = "Upper".into();
    c.system.rivers[0].capacity = 10.;
    c.system.rivers[0].transport = RiverTransport::Fixed { delay: 0.5 };
    let mut lower = river("Join", "Sea", RiverLaw::Controlled);
    lower.name = "Lower".into();
    lower.capacity = 10.;
    c.system.rivers.push(lower);
    c.operations.push(OperationalSeries {
        object: Object::River("Lower".into()),
        attribute: Attribute::MinRelease,
        times: vec![0., 1.],
        values: vec![2., 0.],
    });
    let s = solve(&c);
    assert!((s.river_withdrawal[0][0] - 4.).abs() < 1e-5);
    assert!(s.river_withdrawal[0][1].abs() < 1e-5);
    assert!((s.river_release[1][0] - 2.).abs() < 1e-5);
    c.operations[0].values[0] = 11.;
    assert!(SchedulingProblem::new(&c).is_err());
}
