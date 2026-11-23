//! Open-gate counts and plant shutdown obligations after a gate closes.
use crate::{Error, ScheduleCase};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeSet, HashSet},
    ops::Range,
};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClosureWait {
    pub tunnels: Vec<String>,
    pub plants: Vec<String>,
    pub hours: f64,
    #[serde(default)]
    pub initial_remaining: f64,
}
/// Count positive gate positions, regardless of whether a tunnel carries water.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TunnelInterlock {
    pub name: String,
    pub tunnels: Vec<String>,
    #[serde(default)]
    pub min_open: usize,
    #[serde(default)]
    pub max_open: Option<usize>,
    #[serde(default)]
    pub closures: Vec<ClosureWait>,
}
impl TunnelInterlock {
    pub fn new(name: impl Into<String>, tunnels: Vec<String>) -> Self {
        Self {
            name: name.into(),
            tunnels,
            min_open: 0,
            max_open: None,
            closures: vec![],
        }
    }
}
struct OpenCount {
    tunnels: Vec<usize>,
    bounds: Vec<(usize, usize)>,
}
struct ClosurePeriods {
    tunnels: Vec<usize>,
    units: Vec<usize>,
    initial: Range<usize>,
    after_closure: Vec<Range<usize>>,
}
#[derive(Default)]
pub(crate) struct GateInterlocks {
    counts: Vec<OpenCount>,
    waits: Vec<ClosurePeriods>,
    closing_tunnels: Vec<usize>,
    closing_index: Vec<Option<usize>>,
    initial_open: Vec<bool>,
}
impl GateInterlocks {
    pub fn new(case: &ScheduleCase, topology: &crate::model::Topology) -> Result<Self, Error> {
        let s = &case.system;
        let nt = case.prices.len();
        let start = case.grid[0];
        let mut data = Self::default();
        let mut names = HashSet::new();
        let mut closing = BTreeSet::new();
        data.initial_open = s
            .tunnels
            .iter()
            .map(|e| e.initial_opening.unwrap_or(e.opening) > 0.)
            .collect();
        for rule in &s.tunnel_interlocks {
            let invalid = |message: &str| Error::Input(format!("{}: {message}", rule.name));
            if rule.name.is_empty() || !names.insert(&rule.name) {
                return Err(invalid("interlock needs a unique nonempty name"));
            }
            let mut seen = HashSet::new();
            let tunnels = rule
                .tunnels
                .iter()
                .map(|name| {
                    if !seen.insert(name) {
                        return Err(invalid("repeated interlock tunnel"));
                    }
                    s.tunnels
                        .iter()
                        .position(|e| &e.name == name)
                        .ok_or_else(|| invalid("unknown interlock tunnel"))
                })
                .collect::<Result<Vec<_>, _>>()?;
            let object = crate::Object::TunnelInterlock(rule.name.clone());
            let lo_series = case
                .operations
                .iter()
                .find(|op| op.object == object && op.attribute == crate::Attribute::MinOpen);
            let hi_series = case
                .operations
                .iter()
                .find(|op| op.object == object && op.attribute == crate::Attribute::MaxOpen);
            let maximum = rule.max_open.unwrap_or(tunnels.len());
            if rule.min_open > maximum || maximum > tunnels.len() {
                return Err(invalid("open-gate count bounds exclude every count"));
            }
            let mut bounds = Vec::new();
            for &time in &case.grid[..nt] {
                let lo = lo_series
                    .and_then(|op| op.at(time, None, false))
                    .unwrap_or(rule.min_open as f64);
                let hi = hi_series
                    .and_then(|op| op.at(time, None, false))
                    .unwrap_or(maximum as f64);
                if lo > hi || hi > tunnels.len() as f64 {
                    return Err(invalid("inconsistent scheduled open-gate counts"));
                }
                bounds.push((lo as usize, hi as usize));
            }
            data.counts.push(OpenCount {
                tunnels: tunnels.clone(),
                bounds,
            });
            for wait in &rule.closures {
                if !wait.hours.is_finite()
                    || wait.hours < 0.
                    || !wait.initial_remaining.is_finite()
                    || wait.initial_remaining < 0.
                    || wait.initial_remaining > wait.hours
                    || !(*case.grid.last().unwrap() + wait.hours).is_finite()
                {
                    return Err(invalid(
                        "closure waits must be finite and nonnegative, with remaining wait no longer than the full wait",
                    ));
                }
                let mut seen = HashSet::new();
                let members = wait
                    .tunnels
                    .iter()
                    .map(|name| {
                        if !seen.insert(name) {
                            return Err(invalid("repeated tunnel in a closure group"));
                        }
                        let e = s
                            .tunnels
                            .iter()
                            .position(|e| &e.name == name)
                            .ok_or_else(|| invalid("unknown closure tunnel"))?;
                        if !tunnels.contains(&e) {
                            return Err(invalid(
                                "closure tunnel is not a member of this interlock",
                            ));
                        }
                        Ok(e)
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let mut seen = HashSet::new();
                let mut units = Vec::new();
                for name in &wait.plants {
                    if !seen.insert(name) {
                        return Err(invalid("repeated protected plant"));
                    }
                    let p = s
                        .plants
                        .iter()
                        .position(|p| &p.name == name)
                        .ok_or_else(|| invalid("unknown protected plant"))?;
                    units.extend_from_slice(&topology.plant_units[p]);
                }
                if wait.hours > 0. && (members.is_empty() || units.is_empty()) {
                    return Err(invalid(
                        "a positive closure wait needs tunnels and protected plants",
                    ));
                }
                if wait.initial_remaining > 0.
                    && units.iter().any(|&g| {
                        s.units[g].operation.initial_on
                            || s.units[g].operation.initial_duration
                                < wait.hours - wait.initial_remaining
                    })
                {
                    return Err(invalid(
                        "initial protected-unit history does not cover the pending closure wait",
                    ));
                }
                let initial = 0..case.grid[..nt]
                    .partition_point(|time| *time - start < wait.initial_remaining);
                let after_closure = (0..nt)
                    .map(|t| {
                        t..case.grid[..nt].partition_point(|time| *time - case.grid[t] < wait.hours)
                    })
                    .collect();
                closing.extend(members.iter().copied());
                data.waits.push(ClosurePeriods {
                    tunnels: members,
                    units,
                    initial,
                    after_closure,
                });
            }
        }
        data.closing_tunnels = closing.into_iter().collect();
        data.closing_index = vec![None; s.tunnels.len()];
        for (k, &e) in data.closing_tunnels.iter().enumerate() {
            data.closing_index[e] = Some(k);
        }
        Ok(data)
    }
    pub fn constrain<'a>(
        &self,
        model: &oximo::Model,
        on: &[Vec<oximo::expr::Expr<'a>>],
        open: &[Vec<oximo::expr::Expr<'a>>],
        zero: oximo::expr::Expr<'a>,
        intervals: usize,
    ) {
        use oximo::prelude::*;
        for rule in &self.counts {
            for (t, &(lo, hi)) in rule.bounds.iter().enumerate() {
                let count = rule.tunnels.iter().fold(zero, |sum, &e| sum + open[e][t]);
                constraint!(model, count >= lo as f64);
                constraint!(model, count <= hi as f64);
            }
        }
        variable!(model,0.<=closed[k in 0..self.closing_tunnels.len(),t in 0..intervals]<=1.);
        for (k, &e) in self.closing_tunnels.iter().enumerate() {
            for t in 0..intervals {
                let prior = if t == 0 {
                    zero + if self.initial_open[e] { 1. } else { 0. }
                } else {
                    open[e][t - 1]
                };
                constraint!(model, closed[[k, t]] >= prior - open[e][t]);
                constraint!(model, closed[[k, t]] <= prior);
                constraint!(model, closed[[k, t]] <= 1. - open[e][t]);
            }
        }
        for wait in &self.waits {
            for &g in &wait.units {
                for t in wait.initial.clone() {
                    constraint!(model, on[g][t] == 0.);
                }
                for &e in &wait.tunnels {
                    let k = self.closing_index[e].unwrap();
                    for (t, periods) in wait.after_closure.iter().enumerate() {
                        for p in periods.clone() {
                            constraint!(model, on[g][p] + closed[[k, t]] <= 1.);
                        }
                    }
                }
            }
        }
    }
    pub fn residual(&self, on: &[Vec<bool>], opening: &[Vec<f64>]) -> f64 {
        let mut residual: f64 = 0.;
        for rule in &self.counts {
            for (t, &(lo, hi)) in rule.bounds.iter().enumerate() {
                let count = rule.tunnels.iter().filter(|&&e| opening[e][t] > 0.).count();
                residual = residual.max((lo as f64 - count as f64).max(count as f64 - hi as f64));
            }
        }
        for wait in &self.waits {
            for &g in &wait.units {
                for t in wait.initial.clone() {
                    if on[g][t] {
                        residual = residual.max(1.);
                    }
                }
                for &e in &wait.tunnels {
                    let mut prior = self.initial_open[e];
                    for (t, periods) in wait.after_closure.iter().enumerate() {
                        let current = opening[e][t] > 0.;
                        if prior && !current && periods.clone().any(|p| on[g][p]) {
                            residual = residual.max(1.);
                        }
                        prior = current;
                    }
                }
            }
        }
        residual.max(0.)
    }
}
