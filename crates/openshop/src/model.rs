//! Physical input data. Names are resolved once before numerical assembly.

use crate::Error;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// Finite lower and upper bounds, in the quantity's documented units.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bounds {
    pub min: f64,
    pub max: f64,
}

impl Bounds {
    pub fn new(min: f64, max: f64) -> Self {
        Self { min, max }
    }
    pub fn contains(self, value: f64) -> bool {
        value.is_finite() && self.min <= value && value <= self.max
    }
    fn validate(self, name: &str) -> Result<(), Error> {
        if self.min.is_finite() && self.max.is_finite() && self.min <= self.max {
            Ok(())
        } else {
            Err(Error::Input(format!(
                "{name}: finite ordered bounds required"
            )))
        }
    }
}

/// z(V) = elevation + slope (V-reference_volume) + curvature (V-reference_volume)².
/// Volume is m³ and elevation is metres.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StorageCurve {
    pub reference_volume: f64,
    pub elevation: f64,
    pub slope: f64,
    pub curvature: f64,
}

impl StorageCurve {
    pub fn level(self, volume: f64) -> f64 {
        let v = volume - self.reference_volume;
        self.elevation + self.slope * v + self.curvature * v * v
    }
    /// Exact integral mean, including the equal-volume limit. This head makes
    /// the storage step conserve potential energy without subtracting cubics.
    pub fn mean_level(self, before: f64, after: f64) -> f64 {
        crate::physics::mean_level(&self, before, after)
    }
    /// Integral of level from reference_volume, in m⁴.
    pub fn potential(self, volume: f64) -> f64 {
        let v = volume - self.reference_volume;
        self.elevation * v + self.slope * v * v / 2.0 + self.curvature * v * v * v / 3.0
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reservoir {
    pub name: String,
    pub storage: Bounds,
    pub initial_storage: f64,
    pub level: StorageCurve,
    /// Terminal marginal value, in currency/m³.
    #[serde(default)]
    pub water_value: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Boundary {
    pub name: String,
    pub head: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Junction {
    pub name: String,
    pub head: Bounds,
}

/// Signed flow through a passive pressure connection. A closed opening
/// disconnects pressure and fixes flow to zero; it does not equate the heads.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Tunnel {
    pub name: String,
    pub source: String,
    pub target: String,
    pub capacity: f64,
    pub loss_factor: f64,
    pub opening: f64,
}

/// Turbine hydraulic efficiency, before electrical conversion.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Efficiency {
    pub peak: f64,
    #[serde(default)]
    pub best_flow: f64,
    #[serde(default)]
    pub flow_curvature: f64,
    #[serde(default)]
    pub best_head: f64,
    #[serde(default)]
    pub head_curvature: f64,
}

impl Efficiency {
    pub fn constant(value: f64) -> Self {
        Self {
            peak: value,
            best_flow: 0.,
            flow_curvature: 0.,
            best_head: 0.,
            head_curvature: 0.,
        }
    }
    pub fn value(self, flow: f64, head: f64) -> f64 {
        crate::physics::efficiency(&self, flow, head)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Unit {
    pub name: String,
    pub source: String,
    pub target: String,
    /// Water destination when discharge enters a river instead of the head node.
    #[serde(default)]
    pub discharge_river: Option<String>,
    /// On-state discharge, m³/s. Off-state discharge is exactly zero.
    pub flow: Bounds,
    pub power: Bounds,
    pub head: Bounds,
    pub efficiency: Efficiency,
    /// Head loss coefficient, in metres/(m³/s)².
    #[serde(default)]
    pub loss_factor: f64,
    #[serde(default)]
    pub operation: crate::UnitOperation,
}

impl Unit {
    /// MW, with water density 1000 kg/m³ and gravity 9.81 m/s².
    pub fn generation(&self, flow: f64, head: f64) -> f64 {
        crate::physics::generation(flow, head, self.efficiency.value(flow, head))
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HydroSystem {
    pub reservoirs: Vec<Reservoir>,
    #[serde(default)]
    pub junctions: Vec<Junction>,
    #[serde(default)]
    pub boundaries: Vec<Boundary>,
    pub units: Vec<Unit>,
    #[serde(default)]
    pub tunnels: Vec<Tunnel>,
    #[serde(default)]
    pub rivers: Vec<River>,
}

/// A river reach with a fixed travel time. Reservoir sources have optimized
/// releases. Other source/target names denote internal river junctions; all
/// arriving water is divided according to the outgoing allocation fractions.
/// Discharge is m³/s, delay and history times are hours, inventory is m³.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct River {
    pub name: String,
    pub source: String,
    pub target: String,
    pub capacity: f64,
    pub delay: f64,
    #[serde(default = "one")]
    pub allocation: f64,
    /// Natural inflow entering the top of this reach, after allocation.
    #[serde(default)]
    pub inflow: f64,
    /// Terminal inventory value, currency/m³.
    #[serde(default)]
    pub water_value: f64,
    /// Historical release intervals; omitted history means zero past release.
    #[serde(default)]
    pub history_grid: Vec<f64>,
    #[serde(default)]
    pub history_release: Vec<f64>,
}

fn one() -> f64 {
    1.
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScheduleCase {
    pub system: HydroSystem,
    /// Control-interval boundaries, hours. Durations may differ.
    pub grid: Vec<f64>,
    /// Currency/MWh, one value per control interval.
    pub prices: Vec<f64>,
    /// Reservoir inflows in m³/s; omitted reservoirs have zero inflow.
    #[serde(default)]
    pub inflows: BTreeMap<String, Vec<f64>>,
    /// Prescribed on/off states. Unspecified units are optimized.
    #[serde(default)]
    pub commitment: BTreeMap<String, Vec<bool>>,
}

pub(crate) struct Topology {
    pub unit_nodes: Vec<(usize, usize)>,
    pub incidence: Vec<Vec<(usize, f64)>>,
    pub tunnel_nodes: Vec<(usize, usize)>,
    pub tunnel_incidence: Vec<Vec<(usize, f64)>>,
    pub unit_outfalls: Vec<Option<usize>>,
    pub unit_commitment: Vec<Vec<crate::operation::CommitmentInequality>>,
}

impl ScheduleCase {
    pub(crate) fn topology(&self) -> Result<Topology, Error> {
        let s = &self.system;
        let periods = self.prices.len();
        if periods == 0
            || self.grid.len() != periods + 1
            || !self.grid.iter().all(|v| v.is_finite())
            || !self
                .grid
                .windows(2)
                .all(|w| w[1] > w[0] && (3600. * (w[1] - w[0])).is_finite())
            || !self.prices.iter().all(|v| v.is_finite())
        {
            return Err(Error::Input(
                "finite increasing grid and one price per interval required".into(),
            ));
        }
        if s.reservoirs.is_empty() {
            return Err(Error::Input("at least one reservoir required".into()));
        }
        let mut nodes = HashMap::new();
        for name in s
            .reservoirs
            .iter()
            .map(|r| &r.name)
            .chain(s.junctions.iter().map(|j| &j.name))
            .chain(s.boundaries.iter().map(|b| &b.name))
        {
            let i = nodes.len();
            if name.is_empty() || nodes.insert(name.as_str(), i).is_some() {
                return Err(Error::Input(
                    "node names must be nonempty and unique".into(),
                ));
            }
        }
        for r in &s.reservoirs {
            r.storage.validate(&r.name)?;
            let c = r.level;
            if r.storage.min < 0.
                || r.storage.max <= r.storage.min
                || !r.storage.contains(r.initial_storage)
                || ![
                    c.reference_volume,
                    c.elevation,
                    c.slope,
                    c.curvature,
                    r.water_value,
                ]
                .iter()
                .all(|v| v.is_finite())
                || [r.storage.min, r.storage.max]
                    .iter()
                    .any(|v| c.slope + 2. * c.curvature * (v - c.reference_volume) < 0.)
                || ![c.level(r.storage.min), c.level(r.storage.max)]
                    .iter()
                    .all(|v| v.is_finite())
            {
                return Err(Error::Input(format!(
                    "{}: invalid storage, level curve or initial state",
                    r.name
                )));
            }
        }
        if s.boundaries.iter().any(|b| !b.head.is_finite()) {
            return Err(Error::Input("finite boundary heads required".into()));
        }
        for j in &s.junctions {
            j.head.validate(&j.name)?;
        }
        let mut incidence = vec![Vec::new(); s.reservoirs.len() + s.junctions.len()];
        let mut unit_nodes = Vec::new();
        let mut unit_names = HashMap::new();
        let mut unit_outfalls = Vec::new();
        for (i, u) in s.units.iter().enumerate() {
            if u.name.is_empty() || unit_names.insert(u.name.as_str(), i).is_some() {
                return Err(Error::Input(
                    "unit names must be nonempty and unique".into(),
                ));
            }
            for (name, bounds) in [("flow", u.flow), ("power", u.power), ("head", u.head)] {
                bounds.validate(&format!("{} {name}", u.name))?;
            }
            u.operation.validate(&u.name)?;
            let e = u.efficiency;
            if u.flow.min < 0.
                || u.flow.max <= 0.
                || u.power.min < 0.
                || u.head.min < 0.
                || ![
                    e.peak,
                    e.best_flow,
                    e.flow_curvature,
                    e.best_head,
                    e.head_curvature,
                    u.loss_factor,
                ]
                .iter()
                .all(|v| v.is_finite())
                || e.peak <= 0.
                || e.peak > 1.
                || e.flow_curvature < 0.
                || e.head_curvature < 0.
                || u.loss_factor < 0.
                || [u.flow.min, u.flow.max].iter().any(|&q| {
                    [u.head.min, u.head.max]
                        .iter()
                        .any(|&h| e.value(q, h) <= 0.)
                })
            {
                return Err(Error::Input(format!(
                    "{}: invalid unit characteristic",
                    u.name
                )));
            }
            let source = *nodes
                .get(u.source.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown source {}", u.name, u.source)))?;
            let target = *nodes
                .get(u.target.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown target {}", u.name, u.target)))?;
            if source == target {
                return Err(Error::Input(format!(
                    "{}: identical source and target",
                    u.name
                )));
            }
            let outfall = u
                .discharge_river
                .as_ref()
                .map(|name| {
                    s.rivers
                        .iter()
                        .position(|r| &r.name == name)
                        .ok_or_else(|| {
                            Error::Input(format!("{}: unknown discharge river {name}", u.name))
                        })
                })
                .transpose()?;
            for (node, sign) in [(source, -1.), (target, 1.)] {
                if sign > 0. && outfall.is_some() {
                    continue;
                }
                if node < incidence.len() {
                    incidence[node].push((i, sign));
                }
            }
            unit_nodes.push((source, target));
            unit_outfalls.push(outfall);
        }
        for (name, values) in &self.inflows {
            if !s.reservoirs.iter().any(|r| &r.name == name)
                || values.len() != periods
                || !values.iter().all(|v| v.is_finite() && *v >= 0.)
            {
                return Err(Error::Input(format!(
                    "{name}: reservoir inflow must contain one finite nonnegative rate per interval"
                )));
            }
        }
        for (name, values) in &self.commitment {
            if !unit_names.contains_key(name.as_str()) || values.len() != periods {
                return Err(Error::Input(format!(
                    "{name}: commitment must contain one state per interval"
                )));
            }
        }
        let mut tunnel_nodes = Vec::new();
        let mut tunnel_incidence = vec![Vec::new(); incidence.len()];
        let mut tunnel_names = HashMap::new();
        for (i, tunnel) in s.tunnels.iter().enumerate() {
            if tunnel.name.is_empty()
                || tunnel_names.insert(tunnel.name.as_str(), i).is_some()
                || !tunnel.capacity.is_finite()
                || tunnel.capacity <= 0.
                || !tunnel.loss_factor.is_finite()
                || tunnel.loss_factor <= 0.
                || !tunnel.opening.is_finite()
                || !(0. ..=1.).contains(&tunnel.opening)
            {
                return Err(Error::Input(format!(
                    "{}: invalid tunnel capacity, resistance or opening",
                    tunnel.name
                )));
            }
            let source = *nodes
                .get(tunnel.source.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown source", tunnel.name)))?;
            let target = *nodes
                .get(tunnel.target.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown target", tunnel.name)))?;
            if source == target {
                return Err(Error::Input(format!(
                    "{}: identical tunnel endpoints",
                    tunnel.name
                )));
            }
            for (node, sign) in [(source, -1.), (target, 1.)] {
                if node < tunnel_incidence.len() {
                    tunnel_incidence[node].push((i, sign));
                }
            }
            tunnel_nodes.push((source, target));
        }
        let unit_commitment: Vec<_> = s
            .units
            .iter()
            .map(|u| u.operation.commitment_inequalities(&self.grid))
            .collect();
        for (i, u) in s.units.iter().enumerate() {
            if let Some(on) = self.commitment.get(&u.name)
                && unit_commitment[i].iter().any(|row| row.residual(on) > 0.)
            {
                return Err(Error::Input(format!(
                    "{}: prescribed commitment violates minimum operating times",
                    u.name
                )));
            }
            if let Some(q) = u.operation.initial_flow
                && u.operation.initial_on
                && !u.flow.contains(q)
            {
                return Err(Error::Input(format!(
                    "{}: initial discharge outside on-state range",
                    u.name
                )));
            }
        }
        Ok(Topology {
            unit_nodes,
            incidence,
            tunnel_nodes,
            tunnel_incidence,
            unit_outfalls,
            unit_commitment,
        })
    }
    pub(crate) fn inflow(&self, r: usize, t: usize) -> f64 {
        self.inflows
            .get(&self.system.reservoirs[r].name)
            .map_or(0., |v| v[t])
    }
}
