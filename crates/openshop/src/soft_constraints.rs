//! Bounded, explicitly priced deviations from operating targets.
//!
//! Safety bounds and physical equations remain hard. These limits add only
//! slack variables to their own named target inequalities.
use crate::{Error, ScheduleCase};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashSet,
    ops::{Add, Mul},
};

/// A target uses one state assessment or an average over a declared time window.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "quantity", rename_all = "snake_case", deny_unknown_fields)]
pub enum SoftTarget {
    Storage {
        reservoir: String,
        time: f64,
    },
    Level {
        reservoir: String,
        time: f64,
    },
    UnitPower {
        unit: String,
        start_time: f64,
        end_time: f64,
    },
    PlantPower {
        plant: String,
        start_time: f64,
        end_time: f64,
    },
    UnitDischarge {
        unit: String,
        start_time: f64,
        end_time: f64,
    },
    PlantDischarge {
        plant: String,
        start_time: f64,
        end_time: f64,
    },
}

/// One nominal bound and the largest permitted deviation, in physical units.
/// Costs are currency/m³ for storage or discharge, currency/m for level, and
/// currency/MWh for power. Discharge penalties integrate the missing volume.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PenaltyLimit {
    pub value: f64,
    #[serde(default)]
    pub allowance: f64,
    #[serde(default)]
    pub cost: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SoftConstraint {
    pub name: String,
    pub target: SoftTarget,
    #[serde(default)]
    pub lower: Option<PenaltyLimit>,
    #[serde(default)]
    pub upper: Option<PenaltyLimit>,
}

/// Minimal deviations from the nominal targets, in the target's physical units.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SoftSlack {
    pub name: String,
    pub lower: f64,
    pub upper: f64,
    pub penalty_cost: f64,
}

pub(crate) enum Target {
    Storage(usize, usize),
    Level(usize, usize),
    Output {
        units: Vec<usize>,
        weights: Vec<(usize, f64)>,
        power: bool,
    },
}
impl Target {
    fn value<D: Copy + Add<Output = D> + Mul<f64, Output = D>>(
        &self,
        storage: &[Vec<D>],
        level: &[Vec<D>],
        power: &[Vec<D>],
        flow: &[Vec<D>],
        zero: D,
    ) -> D {
        match self {
            Self::Storage(r, k) => storage[*r][*k],
            Self::Level(r, k) => level[*r][*k],
            Self::Output {
                units,
                weights,
                power: is_power,
            } => {
                let rates = if *is_power { power } else { flow };
                units
                    .iter()
                    .flat_map(|&g| weights.iter().map(move |&(t, w)| rates[g][t] * w))
                    .fold(zero, |sum, value| sum + value)
            }
        }
    }
}
struct Limit {
    name: String,
    target: Option<Target>,
    lower: Option<PenaltyLimit>,
    upper: Option<PenaltyLimit>,
    scale: f64,
    factor: f64,
}
#[derive(Default)]
pub(crate) struct SoftLimits {
    limits: Vec<Limit>,
}
impl SoftLimits {
    pub fn new(case: &ScheduleCase, topology: &crate::model::Topology) -> Result<Self, Error> {
        let system = &case.system;
        let mut names: HashSet<&str> = system
            .reservoirs
            .iter()
            .map(|o| o.name.as_str())
            .chain(system.junctions.iter().map(|o| o.name.as_str()))
            .chain(system.boundaries.iter().map(|o| o.name.as_str()))
            .chain(system.plants.iter().map(|o| o.name.as_str()))
            .chain(system.units.iter().map(|o| o.name.as_str()))
            .chain(system.tunnels.iter().map(|o| o.name.as_str()))
            .chain(system.rivers.iter().map(|o| o.name.as_str()))
            .chain(system.tunnel_interlocks.iter().map(|o| o.name.as_str()))
            .chain(case.flow_requirements.iter().map(|o| o.name.as_str()))
            .collect();
        let mut limits = Vec::new();
        let (start, end) = (case.grid[0], *case.grid.last().unwrap());
        for rule in &case.soft_constraints {
            let invalid = |reason: &str| Error::Input(format!("{}: {reason}", rule.name));
            if rule.name.is_empty()
                || !names.insert(&rule.name)
                || (rule.lower.is_none() && rule.upper.is_none())
            {
                return Err(invalid(
                    "soft constraints need a unique nonempty name and at least one bound",
                ));
            }
            let is_level = matches!(rule.target, SoftTarget::Level { .. });
            for bound in rule.lower.into_iter().chain(rule.upper) {
                if !bound.value.is_finite()
                    || (!is_level && bound.value < 0.)
                    || !bound.allowance.is_finite()
                    || bound.allowance < 0.
                    || !bound.cost.is_finite()
                    || bound.cost < 0.
                    || (bound.allowance > 0. && bound.cost <= 0.)
                {
                    return Err(invalid(
                        "bounds and allowances must be physical and finite; permitted deviations need positive costs",
                    ));
                }
            }
            if matches!((rule.lower, rule.upper), (Some(a), Some(b)) if a.value > b.value) {
                return Err(invalid("lower target exceeds upper target"));
            }
            let (target, scale, factor) = match &rule.target {
                SoftTarget::Storage { reservoir, time } | SoftTarget::Level { reservoir, time } => {
                    let r = system
                        .reservoirs
                        .iter()
                        .position(|r| r.name == *reservoir)
                        .ok_or_else(|| invalid("unknown reservoir"))?;
                    if !time.is_finite() {
                        return Err(invalid("assessment time must be finite"));
                    }
                    let target = if start <= *time && *time <= end {
                        let k = case.grid.iter().position(|t| t == time).ok_or_else(|| {
                            invalid("control grid must include the state assessment time")
                        })?;
                        Some(if is_level {
                            Target::Level(r, k)
                        } else {
                            Target::Storage(r, k)
                        })
                    } else {
                        None
                    };
                    let scale = if is_level {
                        100.
                    } else {
                        (system.reservoirs[r].storage.max - system.reservoirs[r].storage.min)
                            .max(1.)
                    };
                    (target, scale, 1.)
                }
                SoftTarget::UnitPower {
                    unit: object,
                    start_time: a,
                    end_time: b,
                }
                | SoftTarget::UnitDischarge {
                    unit: object,
                    start_time: a,
                    end_time: b,
                }
                | SoftTarget::PlantPower {
                    plant: object,
                    start_time: a,
                    end_time: b,
                }
                | SoftTarget::PlantDischarge {
                    plant: object,
                    start_time: a,
                    end_time: b,
                } => {
                    let power = matches!(
                        rule.target,
                        SoftTarget::UnitPower { .. } | SoftTarget::PlantPower { .. }
                    );
                    let units = if matches!(
                        rule.target,
                        SoftTarget::UnitPower { .. } | SoftTarget::UnitDischarge { .. }
                    ) {
                        vec![
                            system
                                .units
                                .iter()
                                .position(|u| u.name == *object)
                                .ok_or_else(|| invalid("unknown unit"))?,
                        ]
                    } else {
                        let p = system
                            .plants
                            .iter()
                            .position(|p| p.name == *object)
                            .ok_or_else(|| invalid("unknown plant"))?;
                        topology.plant_units[p].clone()
                    };
                    let duration = b - a;
                    if !a.is_finite() || !b.is_finite() || !duration.is_finite() || duration <= 0. {
                        return Err(invalid(
                            "output assessment requires a finite positive time window",
                        ));
                    }
                    let target = if *b <= start || *a >= end {
                        None
                    } else {
                        if *a < start || *b > end {
                            return Err(invalid(
                                "assessment window crosses the scheduling horizon",
                            ));
                        }
                        let lo = case.grid.iter().position(|t| t == a).ok_or_else(|| {
                            invalid("control grid must include assessment window edges")
                        })?;
                        let hi = case.grid.iter().position(|t| t == b).ok_or_else(|| {
                            invalid("control grid must include assessment window edges")
                        })?;
                        Some(Target::Output {
                            units,
                            weights: (lo..hi)
                                .map(|t| (t, (case.grid[t + 1] - case.grid[t]) / duration))
                                .collect(),
                            power,
                        })
                    };
                    (target, 100., duration * if power { 1. } else { 3600. })
                }
            };
            if !factor.is_finite()
                || rule.lower.into_iter().chain(rule.upper).any(|b| {
                    !(scale * factor * b.cost).is_finite()
                        || !(factor * b.cost * b.allowance).is_finite()
                })
            {
                return Err(invalid(
                    "penalty coefficients and bounded costs must be finite",
                ));
            }
            limits.push(Limit {
                name: rule.name.clone(),
                target,
                lower: rule.lower,
                upper: rule.upper,
                scale,
                factor,
            });
        }
        Ok(Self { limits })
    }
    pub fn needs_level(&self, reservoir: usize) -> bool {
        self.limits
            .iter()
            .any(|l| matches!(l.target, Some(Target::Level(r, _)) if r == reservoir))
    }
    pub fn constrain<'a>(
        &self,
        model: &'a oximo::Model,
        storage: &[Vec<oximo::expr::Expr<'a>>],
        level: &[Vec<oximo::expr::Expr<'a>>],
        power: &[Vec<oximo::expr::Expr<'a>>],
        flow: &[Vec<oximo::expr::Expr<'a>>],
        zero: oximo::expr::Expr<'a>,
    ) -> Result<oximo::expr::Expr<'a>, Error> {
        use oximo::prelude::*;
        if self.limits.is_empty() {
            return Ok(zero);
        }
        variable!(model, slack[k in 0..self.limits.len(), side in 0..2] >= 0.);
        let mut cost = zero;
        for (k, limit) in self.limits.iter().enumerate() {
            for (side, bound) in [limit.lower, limit.upper].into_iter().enumerate() {
                let allowance = if limit.target.is_some() {
                    bound.map_or(0., |b| b.allowance)
                } else {
                    0.
                };
                model.set_bounds(slack[[k, side]], 0., allowance / limit.scale)?;
                if let (Some(target), Some(bound)) = (&limit.target, bound) {
                    let value = target.value(storage, level, power, flow, zero) / limit.scale;
                    if side == 0 {
                        constraint!(model, value + slack[[k, side]] >= bound.value / limit.scale);
                    } else {
                        constraint!(model, value - slack[[k, side]] <= bound.value / limit.scale);
                    }
                    cost = cost + slack[[k, side]] * (limit.scale * limit.factor * bound.cost);
                }
            }
        }
        Ok(cost)
    }
    pub fn assess(
        &self,
        case: &ScheduleCase,
        storage: &[Vec<f64>],
        power: &[Vec<f64>],
        flow: &[Vec<f64>],
    ) -> (Vec<SoftSlack>, f64, f64) {
        if self.limits.is_empty() {
            return (vec![], 0., 0.);
        }
        let levels: Vec<Vec<_>> = case
            .system
            .reservoirs
            .iter()
            .zip(storage)
            .map(|(r, row)| row.iter().map(|&v| r.level.trial_level(v)).collect())
            .collect();
        let mut cost = 0.;
        let mut excess: f64 = 0.;
        let slacks = self
            .limits
            .iter()
            .map(|l| {
                let value = l
                    .target
                    .as_ref()
                    .map(|t| t.value(storage, &levels, power, flow, 0.));
                let lower = value
                    .zip(l.lower)
                    .map_or(0., |(v, b)| (b.value - v).max(0.));
                let upper = value
                    .zip(l.upper)
                    .map_or(0., |(v, b)| (v - b.value).max(0.));
                if value.is_some_and(|v| !v.is_finite()) {
                    excess = f64::INFINITY;
                }
                for (deviation, bound) in [(lower, l.lower), (upper, l.upper)] {
                    if let Some(b) = bound {
                        excess = excess.max((deviation - b.allowance).max(0.) / l.scale);
                    }
                }
                let penalty_cost = l.factor
                    * (lower * l.lower.map_or(0., |b| b.cost)
                        + upper * l.upper.map_or(0., |b| b.cost));
                cost += penalty_cost;
                SoftSlack {
                    name: l.name.clone(),
                    lower,
                    upper,
                    penalty_cost,
                }
            })
            .collect();
        (slacks, cost, excess)
    }
}
