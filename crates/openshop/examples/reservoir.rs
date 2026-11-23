//! Construct and optimize a watercourse using the native Rust API.

use openshop::{
    Boundary, Bounds, Efficiency, HydroSystem, Plant, Reservoir, ScheduleCase, SchedulingProblem,
    SolveOptions, StorageCurve, Unit,
};
use std::collections::BTreeMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let system = HydroSystem {
        reservoirs: vec![Reservoir {
            operation: Default::default(),
            name: "Lake".into(),
            storage: Bounds::new(1e6, 4e6),
            initial_storage: 3e6,
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 40.,
                slope: 20e-6,
                curvature: 0.,
            },
            water_value: 0.001,
        }],
        boundaries: vec![Boundary {
            name: "Sea".into(),
            head: 20.,
        }],
        plants: vec![Plant::new("Station", "Lake", "Sea")],
        units: vec![Unit {
            name: "Station".into(),
            plant: "Station".into(),
            mode: openshop::UnitMode::Generate,
            pump_characteristic: None,
            motor_efficiency: 1.,
            generator_efficiency_curve: None,
            flow: Bounds::new(5., 20.),
            power: Bounds::new(1., 25.),
            head: Bounds::new(35., 105.),
            efficiency: Efficiency::constant(0.9),
            loss_factor: 0.,
            operation: openshop::UnitOperation {
                startup_cost: 100.,
                ..Default::default()
            },
        }],
        ..Default::default()
    };
    let case = ScheduleCase {
        routing_grid: vec![],
        operations: Vec::new(),
        system,
        grid: vec![0., 1., 2.],
        prices: vec![-100., 150.],
        inflows: BTreeMap::new(),
        commitment: BTreeMap::new(),
    };
    let problem = SchedulingProblem::new(&case)?;
    let schedule = problem.solve(&SolveOptions {
        time_limit: 5.,
        ..Default::default()
    })?;
    if !schedule.validation.feasible {
        return Err("no feasible schedule found".into());
    }
    println!("Revenue: {:.2}", schedule.objective);
    for (unit, flows) in problem.case().system.units.iter().zip(&schedule.unit_flow) {
        println!("{} discharge (m³/s): {flows:?}", unit.name);
    }
    Ok(())
}
