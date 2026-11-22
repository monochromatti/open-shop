//! Nonlinear hydropower scheduling.
//!
//! Physical inputs use metres, cubic metres, cubic metres per second, hours,
//! megawatts and currency per megawatt-hour. Optimization variables are scaled
//! internally; returned schedules use the input units.

mod commitment;
mod dispatch;
mod model;
mod native_ipopt;
mod nlp;
mod physics;
mod river;
mod schedule;
mod trajectory;

pub use model::{
    Boundary, Bounds, Efficiency, HydroSystem, Junction, Reservoir, River, ScheduleCase,
    StorageCurve, Tunnel, Unit,
};
pub use schedule::SchedulingProblem;
pub use trajectory::{Schedule, Validation};

use serde::{Deserialize, Serialize};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid input: {0}")]
    Input(String),
    #[error("nonlinear model: {0}")]
    Model(String),
    #[error("solver: {0}")]
    Solver(String),
}

/// Shared time and accuracy limits. Time is wall-clock seconds.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct SolveOptions {
    pub time_limit: f64,
    pub max_iterations: i32,
    pub tolerance: f64,
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            time_limit: 60.0,
            max_iterations: 1000,
            tolerance: 1e-8,
        }
    }
}

impl SolveOptions {
    fn validate(&self) -> Result<(), Error> {
        if !self.time_limit.is_finite()
            || self.time_limit <= 0.0
            || std::time::Duration::try_from_secs_f64(self.time_limit).is_err()
            || !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.max_iterations < 0
        {
            return Err(Error::Input(
                "positive finite time and tolerance, and a nonnegative iteration limit, required"
                    .into(),
            ));
        }
        Ok(())
    }
}
