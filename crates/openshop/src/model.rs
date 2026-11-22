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
    pub(crate) fn validate(self, name: &str) -> Result<(), Error> {
        if self.min.is_finite() && self.max.is_finite() && self.min <= self.max {
            Ok(())
        } else {
            Err(Error::Input(format!(
                "{name}: finite ordered bounds required"
            )))
        }
    }
}

/// Reservoir elevation in metres as a function of storage in m³.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum StorageCurve {
    Polynomial {
        reference_volume: f64,
        elevation: f64,
        slope: f64,
        curvature: f64,
    },
    Table(crate::TableCurve),
}
impl StorageCurve {
    pub fn level(&self, volume: f64) -> Result<f64, Error> {
        if let Self::Table(c) = self {
            return c.value(volume);
        }
        let value = self.trial_level(volume);
        if !volume.is_finite() || !value.is_finite() {
            return Err(Error::Input("nonfinite level coordinate or value".into()));
        }
        Ok(value)
    }
    pub(crate) fn trial_level(&self, volume: f64) -> f64 {
        match self {
            Self::Polynomial {
                reference_volume,
                elevation,
                slope,
                curvature,
            } => {
                let v = volume - reference_volume;
                elevation + slope * v + curvature * v * v
            }
            Self::Table(c) => c.extended_value(volume),
        }
    }
    /// Exact integral mean, symmetric and well-defined at equal volumes.
    pub fn mean_level(&self, before: f64, after: f64) -> Result<f64, Error> {
        if let Self::Table(c) = self {
            return c.mean_value(before, after);
        }
        let value = self.trial_mean_level(before, after);
        if !before.is_finite() || !after.is_finite() || !value.is_finite() {
            return Err(Error::Input(
                "nonfinite mean-level coordinate or value".into(),
            ));
        }
        Ok(value)
    }
    pub(crate) fn trial_mean_level(&self, before: f64, after: f64) -> f64 {
        match self {
            Self::Polynomial {
                reference_volume,
                elevation,
                slope,
                curvature,
            } => crate::physics::polynomial_mean_level(
                *reference_volume,
                *elevation,
                *slope,
                *curvature,
                before,
                after,
            ),
            Self::Table(c) => c.mean_data(before, after).value,
        }
    }
    pub(crate) fn mean_data(&self, before: f64, after: f64) -> crate::curve::Derivatives2 {
        match self {
            Self::Table(c) => c.mean_data(before, after),
            Self::Polynomial {
                reference_volume,
                slope,
                curvature,
                ..
            } => {
                let a = before - reference_volume;
                let b = after - reference_volume;
                crate::curve::Derivatives2 {
                    value: self.trial_mean_level(before, after),
                    gradient: [
                        slope / 2. + curvature * (2. * a + b) / 3.,
                        slope / 2. + curvature * (a + 2. * b) / 3.,
                    ],
                    hessian: [
                        [2. * curvature / 3., curvature / 3.],
                        [curvature / 3., 2. * curvature / 3.],
                    ],
                }
            }
        }
    }
    /// Integral of level from the polynomial reference or first table coordinate, m⁴.
    pub fn potential(&self, volume: f64) -> Result<f64, Error> {
        self.level(volume)?;
        Ok(match self {
            Self::Polynomial {
                reference_volume,
                elevation,
                slope,
                curvature,
            } => {
                let v = volume - reference_volume;
                elevation * v + slope * v * v / 2. + curvature * v * v * v / 3.
            }
            Self::Table(c) => {
                c.mean_value(c.coordinates()[0], volume)? * (volume - c.coordinates()[0])
            }
        })
    }
    fn validate(&self, bounds: Bounds, name: &str) -> Result<(), Error> {
        let valid = match self {
            Self::Polynomial {
                reference_volume,
                elevation,
                slope,
                curvature,
            } => {
                [reference_volume, elevation, slope, curvature]
                    .iter()
                    .all(|v| v.is_finite())
                    && [bounds.min, bounds.max].iter().all(|v| {
                        let derivative = slope + 2. * curvature * (v - reference_volume);
                        derivative.is_finite() && derivative >= 0.
                    })
            }
            Self::Table(c) => {
                c.coordinates()[0] <= bounds.min
                    && bounds.max <= c.coordinates()[c.coordinates().len() - 1]
                    && c.values().windows(2).all(|w| w[1] >= w[0])
            }
        };
        if !valid
            || ![bounds.min, bounds.max]
                .iter()
                .all(|&v| self.level(v).is_ok())
        {
            return Err(Error::Input(format!(
                "{name}: storage level curve must cover storage bounds and be finite and nondecreasing"
            )));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reservoir {
    pub name: String,
    pub storage: Bounds,
    pub initial_storage: f64,
    #[serde(default)]
    pub operation: crate::ReservoirOperation,
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
    /// Head multiplier in [0,1] when no gate-loss curve is supplied.
    pub opening: f64,
    /// Additional signed-flow resistance as a function of gate position.
    #[serde(default)]
    pub gate_loss_factor_curve: Option<crate::TableCurve>,
}

/// Turbine hydraulic efficiency, before electrical conversion.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Efficiency {
    Quadratic {
        peak: f64,
        #[serde(default)]
        best_flow: f64,
        #[serde(default)]
        flow_curvature: f64,
        #[serde(default)]
        best_head: f64,
        #[serde(default)]
        head_curvature: f64,
    },
    Table(crate::TurbineTable),
}
impl Efficiency {
    pub fn constant(value: f64) -> Self {
        Self::Quadratic {
            peak: value,
            best_flow: 0.,
            flow_curvature: 0.,
            best_head: 0.,
            head_curvature: 0.,
        }
    }
    pub fn value(&self, flow: f64, head: f64) -> Result<f64, Error> {
        if !flow.is_finite() || !head.is_finite() {
            return Err(Error::Input("nonfinite efficiency coordinate".into()));
        }
        match self {
            Self::Quadratic { .. } => Ok(self.trial_value(flow, head)),
            Self::Table(c) => c.efficiency(flow, head),
        }
    }
    pub(crate) fn trial_value(&self, flow: f64, head: f64) -> f64 {
        match self {
            Self::Quadratic {
                peak,
                best_flow,
                flow_curvature,
                best_head,
                head_curvature,
            } => crate::physics::quadratic_efficiency(
                *peak,
                *best_flow,
                *flow_curvature,
                *best_head,
                *head_curvature,
                flow,
                head,
            ),
            Self::Table(c) => c.trial_efficiency(flow, head),
        }
    }
    fn validate_domain(&self, flow: Bounds, head: Bounds) -> Result<(), Error> {
        match self {
            Self::Quadratic {
                peak,
                best_flow,
                flow_curvature,
                best_head,
                head_curvature,
            } => {
                if ![peak, best_flow, flow_curvature, best_head, head_curvature]
                    .iter()
                    .all(|v| v.is_finite())
                    || *peak <= 0.
                    || *peak > 1.
                    || *flow_curvature < 0.
                    || *head_curvature < 0.
                    || [flow.min, flow.max].iter().any(|&q| {
                        [head.min, head.max]
                            .iter()
                            .any(|&h| self.trial_value(q, h) <= 0.)
                    })
                {
                    return Err(Error::Input("invalid quadratic turbine efficiency".into()));
                }
                Ok(())
            }
            Self::Table(c) => c.validate_domain(flow, head),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UnitMode {
    #[default]
    Generate,
    Pump,
}
impl UnitMode {
    pub fn power_sign(self) -> f64 {
        match self {
            Self::Generate => 1.,
            Self::Pump => -1.,
        }
    }
}

/// Common intake, outlet and operating restrictions for units at one station.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Plant {
    pub name: String,
    pub source: String,
    pub target: String,
    #[serde(default)]
    pub discharge_river: Option<String>,
    #[serde(default)]
    pub flow_min: f64,
    #[serde(default)]
    pub flow_max: Option<f64>,
    #[serde(default)]
    pub power_min: f64,
    #[serde(default)]
    pub power_max: Option<f64>,
    #[serde(default)]
    pub tailwater_curve: Option<crate::TableCurve>,
    #[serde(default)]
    pub outlet_head_floor: Option<f64>,
    #[serde(default)]
    pub operation: crate::PlantOperation,
}
impl Plant {
    pub fn new(
        name: impl Into<String>,
        source: impl Into<String>,
        target: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            source: source.into(),
            target: target.into(),
            discharge_river: None,
            flow_min: 0.,
            flow_max: None,
            power_min: 0.,
            power_max: None,
            tailwater_curve: None,
            outlet_head_floor: None,
            operation: Default::default(),
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CommitmentExclusion {
    pub name: String,
    pub groups: Vec<Vec<String>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Unit {
    pub name: String,
    pub plant: String,
    #[serde(default)]
    pub mode: UnitMode,
    /// On-state magnitudes; pumping consumption is positive here and in results.
    pub flow: Bounds,
    pub power: Bounds,
    pub head: Bounds,
    pub efficiency: Efficiency,
    #[serde(default)]
    pub pump_characteristic: Option<crate::PumpCharacteristic>,
    #[serde(default = "one")]
    pub motor_efficiency: f64,
    #[serde(default)]
    pub generator_efficiency_curve: Option<crate::TableCurve>,
    #[serde(default)]
    pub loss_factor: f64,
    #[serde(default)]
    pub operation: crate::UnitOperation,
}
impl Unit {
    pub fn new(
        name: impl Into<String>,
        plant: impl Into<String>,
        flow: Bounds,
        power: Bounds,
        head: Bounds,
        efficiency: Efficiency,
    ) -> Self {
        Self {
            name: name.into(),
            plant: plant.into(),
            mode: UnitMode::Generate,
            flow,
            power,
            head,
            efficiency,
            pump_characteristic: None,
            motor_efficiency: 1.,
            generator_efficiency_curve: None,
            loss_factor: 0.,
            operation: Default::default(),
        }
    }
    /// Electrical output or consumption, MW, at discharge and net working head.
    pub fn power(&self, flow: f64, head: f64) -> Result<f64, Error> {
        if !flow.is_finite() || !head.is_finite() {
            return Err(Error::Input("nonfinite power coordinates".into()));
        }
        if flow == 0. {
            return Ok(0.);
        }
        let eta = if let Some(c) = &self.pump_characteristic {
            c.efficiency(head)?
        } else {
            self.efficiency.value(flow, head)?
        };
        let hydraulic = 0.00981 * flow * head;
        match self.mode {
            UnitMode::Pump => Ok(hydraulic / (eta * self.motor_efficiency)),
            UnitMode::Generate => match &self.generator_efficiency_curve {
                Some(c) => crate::machine::electrical_power(c, hydraulic * eta),
                None => Ok(hydraulic * eta),
            },
        }
    }
    pub(crate) fn trial_power(&self, flow: f64, head: f64) -> f64 {
        if flow == 0. {
            return 0.;
        }
        let eta = if let Some(c) = &self.pump_characteristic {
            c.trial_efficiency(head)
        } else {
            self.efficiency.trial_value(flow, head)
        };
        let hydraulic = 0.00981 * flow * head;
        match self.mode {
            UnitMode::Pump => hydraulic / (eta * self.motor_efficiency),
            UnitMode::Generate => self
                .generator_efficiency_curve
                .as_ref()
                .map_or(hydraulic * eta, |c| {
                    crate::machine::electrical_power(c, hydraulic * eta).unwrap_or(f64::NAN)
                }),
        }
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HydroSystem {
    #[serde(default)]
    pub plants: Vec<Plant>,
    #[serde(default)]
    pub commitment_exclusions: Vec<CommitmentExclusion>,
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
    #[serde(default)]
    pub law: crate::RiverLaw,
    #[serde(default = "gate_range")]
    pub gate: Bounds,
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

fn gate_range() -> Bounds {
    Bounds::new(0., 1.)
}

fn one() -> f64 {
    1.
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScheduleCase {
    #[serde(default)]
    pub operations: Vec<crate::OperationalSeries>,
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
    pub unit_plants: Vec<usize>,
    pub plant_units: Vec<Vec<usize>>,
    pub plant_operation: Vec<crate::UnitOperation>,
    pub plant_commitment: Vec<Vec<crate::operation::CommitmentInequality>>,
    pub exclusions: Vec<Vec<Vec<usize>>>,
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
            r.level.validate(r.storage, &r.name)?;
            if r.storage.min < 0.
                || r.storage.max <= r.storage.min
                || !r.storage.contains(r.initial_storage)
                || !r.water_value.is_finite()
            {
                return Err(Error::Input(format!(
                    "{}: invalid storage or initial state",
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
        let mut plant_names = HashMap::new();
        let mut plant_nodes = Vec::new();
        for (i, p) in s.plants.iter().enumerate() {
            if p.name.is_empty() || plant_names.insert(p.name.as_str(), i).is_some() {
                return Err(Error::Input(
                    "plant names must be nonempty and unique".into(),
                ));
            }
            let source = *nodes
                .get(p.source.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown intake", p.name)))?;
            let target = *nodes
                .get(p.target.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown outlet", p.name)))?;
            if source == target
                || !p.flow_min.is_finite()
                || p.flow_min < 0.
                || !p.power_min.is_finite()
                || p.power_min < 0.
                || p.flow_max.is_some_and(|x| !x.is_finite() || x < p.flow_min)
                || p.power_max
                    .is_some_and(|x| !x.is_finite() || x < p.power_min)
                || p.outlet_head_floor.is_some_and(|x| !x.is_finite())
            {
                return Err(Error::Input(format!("{}: invalid plant limits", p.name)));
            }
            plant_nodes.push((source, target));
        }
        let mut unit_plants = Vec::new();
        let mut plant_units = vec![Vec::new(); s.plants.len()];
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
            if u.flow.min < 0.
                || u.flow.max <= 0.
                || u.power.min < 0.
                || u.head.min < 0.
                || !u.loss_factor.is_finite()
                || u.loss_factor < 0.
            {
                return Err(Error::Input(format!(
                    "{}: invalid unit characteristic",
                    u.name
                )));
            }
            if u.pump_characteristic.is_none() {
                u.efficiency
                    .validate_domain(u.flow, u.head)
                    .map_err(|e| Error::Input(format!("{}: {e}", u.name)))?;
            }
            crate::machine::validate_unit(u)?;
            let plant = *plant_names
                .get(u.plant.as_str())
                .ok_or_else(|| Error::Input(format!("{}: unknown plant {}", u.name, u.plant)))?;
            unit_plants.push(plant);
            plant_units[plant].push(i);
            let (source, target) = plant_nodes[plant];
            if source == target {
                return Err(Error::Input(format!(
                    "{}: identical source and target",
                    u.name
                )));
            }
            let outfall = s.plants[plant]
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
        let mut plant_operation = Vec::new();
        let mut plant_commitment = Vec::new();
        for (i, p) in s.plants.iter().enumerate() {
            if plant_units[i].is_empty() {
                return Err(Error::Input(format!("{}: plant has no units", p.name)));
            }
            let members: Vec<_> = plant_units[i].iter().map(|&g| &s.units[g]).collect();
            if members.iter().any(|g| g.mode != members[0].mode) {
                return Err(Error::Input(format!(
                    "{}: generating and pumping units need separate, oppositely connected plants",
                    p.name
                )));
            }
            if let Some(c) = &p.tailwater_curve {
                let maximum: f64 = members.iter().map(|g| g.flow.max).sum();
                if c.coordinates()[0] != 0.
                    || *c.coordinates().last().unwrap() < maximum
                    || c.values().iter().any(|x| *x < 0.)
                    || c.values().windows(2).any(|w| w[1] < w[0])
                {
                    return Err(Error::Input(format!(
                        "{}: tailwater must be nonnegative, nondecreasing and cover aggregate flow",
                        p.name
                    )));
                }
            }
            let operation = p.operation.resolve(&members, &p.name)?;
            plant_commitment.push(operation.commitment_inequalities(&self.grid));
            plant_operation.push(operation);
        }
        let mut exclusions = Vec::new();
        let mut exclusion_names = std::collections::HashSet::new();
        for rule in &s.commitment_exclusions {
            if rule.name.is_empty()
                || !exclusion_names.insert(&rule.name)
                || rule.groups.len() < 2
                || rule.groups.iter().any(Vec::is_empty)
            {
                return Err(Error::Input(
                    "mode exclusion needs a unique name and at least two nonempty groups".into(),
                ));
            }
            let mut seen = std::collections::HashSet::new();
            let mut groups = Vec::new();
            for group in &rule.groups {
                let mut units = Vec::new();
                for name in group {
                    let &g = unit_names.get(name.as_str()).ok_or_else(|| {
                        Error::Input(format!("{}: unknown exclusion unit {name}", rule.name))
                    })?;
                    if !seen.insert(g) {
                        return Err(Error::Input(format!(
                            "{}: repeated exclusion unit",
                            rule.name
                        )));
                    }
                    units.push(g);
                }
                groups.push(units);
            }
            if groups
                .iter()
                .filter(|group| group.iter().any(|&g| s.units[g].operation.initial_on))
                .count()
                > 1
            {
                return Err(Error::Input(format!(
                    "{}: initial units violate mode exclusion",
                    rule.name
                )));
            }
            for t in 0..periods {
                if groups
                    .iter()
                    .filter(|group| {
                        group.iter().any(|&g| {
                            self.commitment
                                .get(&s.units[g].name)
                                .is_some_and(|on| on[t])
                        })
                    })
                    .count()
                    > 1
                {
                    return Err(Error::Input(format!(
                        "{}: prescribed units violate mode exclusion",
                        rule.name
                    )));
                }
            }
            exclusions.push(groups);
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
                || tunnel.capacity < 0.
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
            if let Some(curve) = &tunnel.gate_loss_factor_curve
                && (curve.coordinates()[0] > 0.
                    || *curve.coordinates().last().unwrap() < 1.
                    || curve.values().iter().any(|v| *v < 0.))
            {
                return Err(Error::Input(format!(
                    "{}: gate resistance table must cover [0,1] and be nonnegative",
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
        }
        for (p, members) in plant_units.iter().enumerate() {
            if members
                .iter()
                .all(|&g| self.commitment.contains_key(&s.units[g].name))
            {
                let on: Vec<_> = (0..periods)
                    .map(|t| {
                        members
                            .iter()
                            .any(|&g| self.commitment[&s.units[g].name][t])
                    })
                    .collect();
                if plant_commitment[p].iter().any(|row| row.residual(&on) > 0.) {
                    return Err(Error::Input(format!(
                        "{}: prescribed commitment violates plant operating times",
                        s.plants[p].name
                    )));
                }
            }
        }
        Ok(Topology {
            unit_nodes,
            unit_plants,
            plant_units,
            plant_operation,
            plant_commitment,
            exclusions,
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
