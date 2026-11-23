//! Nonlinear hydropower scheduling.
//!
//! Physical inputs use metres, cubic metres, cubic metres per second, hours,
//! megawatts and currency per megawatt-hour. Optimization variables are scaled
//! internally; returned schedules use the input units.

mod commitment;
mod constitutive;
mod curve;
mod dispatch;
mod flow_requirements;
mod level_windows;
mod machine;
mod model;
mod native_ipopt;
mod nlp;
mod operating_data;
mod operating_equations;
mod operation;
mod physics;
mod river;
mod river_laws;
mod river_transport;
mod schedule;
mod trajectory;
mod turbine;

pub use curve::TableCurve;
pub use flow_requirements::{FlowRequirement, PeriodRamp, RollingRamp, VolumeWindow};
pub use level_windows::{LevelWindow, StorageHistory};
pub use machine::PumpCharacteristic;
pub use model::{
    Boundary, Bounds, CommitmentExclusion, Efficiency, HydroSystem, Junction, Plant, Reservoir,
    River, ScheduleCase, StorageCurve, Tunnel, Unit, UnitMode,
};
pub use operating_data::{Attribute, Object, OperationalSeries, ReservoirOperation};
pub use operation::{PlantOperation, RampRates, UnitOperation};
pub use river_laws::{LevelFlowCurve, ReferenceLevel, RiverLaw};
pub use river_transport::{DelayCurve, RiverTransport};
pub use schedule::SchedulingProblem;
pub use trajectory::{Schedule, Validation};
pub use turbine::{DischargeInterpolation, Extrapolation, TurbineOptions, TurbineTable};

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
