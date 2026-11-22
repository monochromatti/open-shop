//! Construct and optimize a watercourse using the native Rust API.

use openshop::{
    Boundary, Bounds, Efficiency, HydroSystem, Reservoir, ScheduleCase, SchedulingProblem,
    SolveOptions, StorageCurve, Unit,
};
use std::collections::BTreeMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let system = HydroSystem {
        reservoirs: vec![Reservoir {
            name: "Lake".into(),
            storage: Bounds::new(1e6, 4e6),
            initial_storage: 3e6,
            level: StorageCurve {
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
        units: vec![Unit {
            name: "Station".into(),
            source: "Lake".into(),
            target: "Sea".into(),
            discharge_river: None,
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
