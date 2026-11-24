//! Optimize a JSON watercourse and write its schedule to standard output.

use openshop::{ScheduleCase, SchedulingProblem, SolveOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = std::env::args().skip(1);
    let path = arguments.next();
    let seconds = arguments
        .next()
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(10.);
    let case: ScheduleCase = if let Some(path) = path {
        serde_json::from_reader(std::fs::File::open(path)?)?
    } else {
        serde_json::from_str(include_str!("data/tabulated_cascade.json"))?
    };
    let schedule = SchedulingProblem::new(&case)?.solve(&SolveOptions {
        time_limit: seconds,
        ..Default::default()
    })?;
    println!("{}", serde_json::to_string_pretty(&schedule)?);
    if !schedule.validation.feasible {
        return Err("no physically feasible schedule was found".into());
    }
    Ok(())
}
