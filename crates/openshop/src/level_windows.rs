//! Reservoir surface ranges over absolute-hour windows, including prior storage.
use crate::{
    Bounds, Error, ScheduleCase,
    river::{LinearForm, Routing},
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LevelWindow {
    pub start_time: f64,
    pub end_time: f64,
    pub max_change: f64,
}
/// Observed storage vertices, m³, connected by linear interpolation in time.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StorageHistory {
    pub times: Vec<f64>,
    pub storage: Vec<f64>,
}
pub(crate) struct StoragePoint {
    pub interval: usize,
    pub fraction: f64,
    pub correction: LinearForm,
}
pub(crate) struct LevelBand {
    pub reservoir: usize,
    pub max_change: f64,
    pub floor: Bounds,
    pub points: Vec<StoragePoint>,
    pub known_levels: Vec<f64>,
}
#[derive(Default)]
pub(crate) struct LevelEquations {
    pub bands: Vec<LevelBand>,
}
impl LevelEquations {
    pub fn new(
        case: &ScheduleCase,
        routing: &Routing,
        mesh: &crate::time_discretization::TimeDiscretization,
    ) -> Result<Self, Error> {
        let mut bands = Vec::new();
        let (start, end) = (case.grid[0], *case.grid.last().unwrap());
        for (i, r) in case.system.reservoirs.iter().enumerate() {
            let invalid = |message: &str| Error::Input(format!("{}: {message}", r.name));
            if let Some(h) = &r.operation.history {
                if h.times.is_empty()
                    || h.times.len() != h.storage.len()
                    || h.times.last() != Some(&start)
                    || h.storage.last() != Some(&r.initial_storage)
                    || !h.times.iter().all(|v| v.is_finite())
                    || !h
                        .times
                        .windows(2)
                        .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())
                    || !h.storage.iter().all(|v| v.is_finite() && *v >= 0.)
                {
                    return Err(invalid(
                        "storage history must be finite, ordered and end at the initial storage and time",
                    ));
                }
                let lo = h.storage.iter().copied().fold(r.storage.min, f64::min);
                let hi = h.storage.iter().copied().fold(r.storage.max, f64::max);
                r.level.validate(Bounds::new(lo, hi), &r.name)?;
            }
            for w in &r.operation.level_windows {
                if !w.start_time.is_finite()
                    || !w.end_time.is_finite()
                    || w.end_time <= w.start_time
                    || !(w.end_time - w.start_time).is_finite()
                    || !w.max_change.is_finite()
                    || w.max_change < 0.
                {
                    return Err(invalid(
                        "level window requires ordered finite times and nonnegative maximum range",
                    ));
                }
                if w.start_time >= end || w.end_time <= start {
                    continue;
                }
                for time in [w.start_time, w.end_time] {
                    if start <= time && time <= end && !case.grid.contains(&time) {
                        return Err(invalid(
                            "control grid must include level-window endpoints inside the horizon",
                        ));
                    }
                }
                let mut known_levels = Vec::new();
                if w.start_time < start {
                    let h = r
                        .operation
                        .history
                        .as_ref()
                        .filter(|h| h.times[0] <= w.start_time)
                        .ok_or_else(|| {
                            invalid("level window requires complete earlier storage history")
                        })?;
                    let a = w.start_time;
                    let b = w.end_time.min(start);
                    let mut times = vec![a];
                    times.extend(h.times.iter().copied().filter(|t| a < *t && *t < b));
                    times.push(b);
                    for time in times {
                        let k = h
                            .times
                            .partition_point(|t| *t <= time)
                            .saturating_sub(1)
                            .min(h.times.len() - 2);
                        let fraction = (time - h.times[k]) / (h.times[k + 1] - h.times[k]);
                        known_levels.push(
                            r.level.level(
                                h.storage[k] * (1. - fraction) + h.storage[k + 1] * fraction,
                            )?,
                        );
                    }
                }
                if w.start_time <= start {
                    known_levels.push(r.level.level(r.initial_storage)?);
                }
                let lower = r.level.level(r.storage.min)?;
                let upper = r.level.level(r.storage.max)?;
                let mut floor = Bounds::new(lower - w.max_change, upper);
                for &level in &known_levels {
                    floor.min = floor.min.max(level - w.max_change);
                    floor.max = floor.max.min(level);
                }
                floor.validate(&r.name).map_err(|_| {
                    invalid("known history violates the level range or its numeric domain")
                })?;
                if lower == upper {
                    continue;
                }
                let mut points = Vec::new();
                let mut previous = None;
                for (t, path) in routing.storage_shape[i].iter().enumerate() {
                    for (fraction, correction) in path {
                        let time = if *fraction == 0. {
                            mesh.grid[t]
                        } else if *fraction == 1. {
                            mesh.grid[t + 1]
                        } else {
                            mesh.grid[t] + fraction * (mesh.grid[t + 1] - mesh.grid[t])
                        };
                        if time < start.max(w.start_time)
                            || time > end.min(w.end_time)
                            || time == start
                            || previous == Some(time)
                        {
                            continue;
                        }
                        points.push(StoragePoint {
                            interval: t,
                            fraction: *fraction,
                            correction: correction.clone(),
                        });
                        previous = Some(time);
                    }
                }
                bands.push(LevelBand {
                    reservoir: i,
                    max_change: w.max_change,
                    floor,
                    points,
                    known_levels,
                });
            }
        }
        Ok(Self { bands })
    }
    pub fn residual(&self, case: &ScheduleCase, storage: &[Vec<f64>], columns: &[f64]) -> f64 {
        self.bands
            .iter()
            .map(|band| {
                let r = &case.system.reservoirs[band.reservoir];
                let mut low = f64::INFINITY;
                let mut high = f64::NEG_INFINITY;
                for level in band
                    .known_levels
                    .iter()
                    .copied()
                    .chain(band.points.iter().map(|p| {
                        r.level
                            .trial_level(p.volume(&storage[band.reservoir], columns))
                    }))
                {
                    if !level.is_finite() {
                        return f64::INFINITY;
                    }
                    low = low.min(level);
                    high = high.max(level);
                }
                (high - low - band.max_change).max(0.) / 100.
            })
            .fold(0., f64::max)
    }
}
impl StoragePoint {
    pub fn volume(&self, storage: &[f64], columns: &[f64]) -> f64 {
        let t = self.interval;
        storage[t] * (1. - self.fraction)
            + storage[t + 1] * self.fraction
            + self.correction.value(columns)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{HydroSystem, Reservoir, StorageCurve};
    #[test]
    fn level_windows_check_intermediate_storage_after_hydraulic_refinement() {
        let mut reservoir = Reservoir {
            name: "Pool".into(),
            storage: Bounds::new(0., 1e6),
            initial_storage: 5e5,
            water_value: 0.,
            operation: Default::default(),
            level: StorageCurve::Polynomial {
                reference_volume: 0.,
                elevation: 0.,
                slope: 1e-4,
                curvature: 0.,
            },
        };
        reservoir.operation.level_windows.push(LevelWindow {
            start_time: 0.,
            end_time: 1.,
            max_change: 5.,
        });
        let c = ScheduleCase::new(
            HydroSystem {
                reservoirs: vec![reservoir],
                ..Default::default()
            },
            vec![0., 1.],
            vec![0.],
        );
        let topology = c.topology().unwrap();
        let operating = crate::operating_data::OperatingData::new(&c, &topology).unwrap();
        let mesh =
            crate::time_discretization::TimeDiscretization::new(&c.grid, &[0., 0.5, 1.]).unwrap();
        let routing = Routing::new(&c, &topology, &operating, &mesh, &c.routing_grid).unwrap();
        let equations = LevelEquations::new(&c, &routing, &mesh).unwrap();
        assert!(equations.residual(&c, &[vec![5e5, 6e5, 5e5]], &[]) > 0.);
        assert_eq!(equations.residual(&c, &[vec![5e5, 5e5, 5e5]], &[]), 0.);
    }
}
