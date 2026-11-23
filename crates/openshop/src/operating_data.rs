//! Resolve named time series into immutable interval and vertex data.
//! Name lookup and input interpretation stop at this preparation boundary.

use crate::{Bounds, Error, RampRates, ScheduleCase, model::Topology};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Object {
    Reservoir(String),
    Unit(String),
    Plant(String),
    Tunnel(String),
    River(String),
    FlowRequirement(String),
    TunnelInterlock(String),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Attribute {
    StorageMin,
    StorageMax,
    FlowMin,
    FlowMax,
    PowerMin,
    PowerMax,
    ForcedOn,
    Maintenance,
    Power,
    Discharge,
    StartupCost,
    ShutdownCost,
    PowerRampUp,
    PowerRampDown,
    FlowRampUp,
    FlowRampDown,
    VolumeRampUp,
    VolumeRampDown,
    LevelRampUp,
    LevelRampDown,
    Opening,
    Capacity,
    Inflow,
    GateMin,
    GateMax,
    GrossHeadMin,
    GrossHeadMax,
    StartupGrossHeadMin,
    MinOpen,
    MinRelease,
    MinArrival,
    FlowChangeCost,
    MaxOpen,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OperationalSeries {
    pub object: Object,
    pub attribute: Attribute,
    pub times: Vec<f64>,
    pub values: Vec<f64>,
}
impl OperationalSeries {
    pub(crate) fn at(&self, time: f64, default: Option<f64>, left: bool) -> Option<f64> {
        let k = self
            .times
            .partition_point(|&t| if left { t < time } else { t <= time });
        if k == 0 {
            default
        } else {
            Some(self.values[k - 1])
        }
    }
    fn integral(&self, a: f64, b: f64, default: Option<f64>) -> Option<f64> {
        let mut k = self.times.partition_point(|&t| t <= a);
        let mut rate = self.at(a, default, false)?;
        let mut previous = a;
        let mut integral = 0.;
        while k < self.times.len() && self.times[k] < b {
            integral += rate * (self.times[k] - previous);
            previous = self.times[k];
            rate = self.values[k];
            k += 1;
        }
        Some(integral + rate * (b - previous))
    }
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ReservoirOperation {
    pub volume_ramp: RampRates,
    pub level_ramp: RampRates,
    pub level_windows: Vec<crate::LevelWindow>,
    pub history: Option<crate::StorageHistory>,
}
#[derive(Clone, Copy, Default)]
pub(crate) struct ChangeLimit {
    pub up: Option<f64>,
    pub down: Option<f64>,
}
pub(crate) struct MachineInterval {
    pub flow: Bounds,
    pub power: Bounds,
    pub forced: Option<bool>,
    pub fixed_power: Option<f64>,
    pub fixed_flow: Option<f64>,
    pub startup_cost: f64,
    pub shutdown_cost: f64,
    pub stopping_power_min: f64,
    pub gross_head_min: Option<f64>,
    pub gross_head_max: Option<f64>,
    pub startup_gross_head_min: Option<f64>,
    pub flow_change: ChangeLimit,
    pub power_change: ChangeLimit,
}
pub(crate) struct PlantInterval {
    pub flow: Bounds,
    pub power: Bounds,
    pub forced: Option<bool>,
    pub fixed_power: Option<f64>,
    pub fixed_flow: Option<f64>,
    pub flow_change: ChangeLimit,
    pub power_change: ChangeLimit,
}
pub(crate) struct ReservoirData {
    pub storage: Vec<Bounds>,
    pub interval_storage: Vec<Bounds>,
    pub inflow: Vec<f64>,
    pub volume_change: Vec<ChangeLimit>,
    pub level_change: Vec<ChangeLimit>,
}
#[derive(Clone, Copy)]
pub(crate) struct TunnelPosition {
    pub opening: f64,
    pub resistance: f64,
}
pub(crate) struct TunnelInterval {
    pub capacity: f64,
    pub positions: Vec<TunnelPosition>,
}
pub(crate) struct OperatingData {
    pub units: Vec<Vec<MachineInterval>>,
    pub plants: Vec<Vec<PlantInterval>>,
    pub plant_initial_flow: Vec<Option<f64>>,
    pub plant_initial_power: Vec<Option<f64>>,
    pub reservoirs: Vec<ReservoirData>,
    pub tunnels: Vec<Vec<TunnelInterval>>,
    pub river_capacity: Vec<Vec<f64>>,
    pub river_inflow: Vec<Vec<f64>>,
    pub river_min_release: Vec<Vec<f64>>,
    pub river_min_arrival: Vec<Vec<f64>>,
    pub river_flow_change: Vec<Vec<ChangeLimit>>,
    pub river_change_cost: Vec<Vec<f64>>,
    pub river_gate: Vec<Vec<Bounds>>,
    pub observed_flow: Vec<Vec<crate::flow_requirements::FlowInterval>>,
}
type SeriesMap<'a> = HashMap<Attribute, &'a OperationalSeries>;
fn value(
    map: Option<&SeriesMap<'_>>,
    attribute: Attribute,
    t: f64,
    default: Option<f64>,
    left: bool,
) -> Option<f64> {
    map.and_then(|m| m.get(&attribute))
        .map_or(default, |s| s.at(t, default, left))
}
fn change(
    map: Option<&SeriesMap<'_>>,
    attribute: Attribute,
    a: f64,
    b: f64,
    rate: Option<f64>,
) -> Option<f64> {
    map.and_then(|m| m.get(&attribute))
        .map_or_else(|| rate.map(|r| r * (b - a)), |s| s.integral(a, b, rate))
}
fn forced(map: Option<&SeriesMap<'_>>, time: f64) -> Result<Option<bool>, Error> {
    let on = value(map, Attribute::ForcedOn, time, None, false)
        .and_then(|v| if v < 0. { None } else { Some(v != 0.) });
    let maintenance = value(map, Attribute::Maintenance, time, Some(0.), false).unwrap() != 0.;
    if maintenance && on == Some(true) {
        return Err(Error::Input(
            "forced operation conflicts with maintenance".into(),
        ));
    }
    Ok(if maintenance { Some(false) } else { on })
}
impl OperatingData {
    pub fn new(case: &ScheduleCase, topology: &Topology) -> Result<Self, Error> {
        use Attribute::*;
        let system = &case.system;
        let nt = case.prices.len();
        let mut series: HashMap<&Object, SeriesMap<'_>> = HashMap::new();
        let mut seen = HashSet::new();
        for s in &case.operations {
            if !seen.insert((&s.object, s.attribute))
                || s.times.is_empty()
                || s.times.len() != s.values.len()
                || !s.times.iter().chain(&s.values).all(|v| v.is_finite())
                || !s.times.windows(2).all(|w| w[1] > w[0])
            {
                return Err(Error::Input(
                    "operating series need unique attributes and paired, ordered finite data"
                        .into(),
                ));
            }
            if s.times
                .iter()
                .any(|&t| case.grid[0] < t && t < case.grid[nt] && !case.grid.contains(&t))
            {
                return Err(Error::Input(
                    "control grid must include every operating change inside the horizon".into(),
                ));
            }
            let allowed = match &s.object {
                Object::Unit(name) => {
                    system.units.iter().any(|x| &x.name == name)
                        && matches!(
                            s.attribute,
                            GrossHeadMin
                                | GrossHeadMax
                                | StartupGrossHeadMin
                                | FlowMin
                                | FlowMax
                                | PowerMin
                                | PowerMax
                                | ForcedOn
                                | Maintenance
                                | Power
                                | Discharge
                                | StartupCost
                                | ShutdownCost
                                | PowerRampUp
                                | PowerRampDown
                                | FlowRampUp
                                | FlowRampDown
                        )
                }
                Object::Plant(name) => {
                    system.plants.iter().any(|x| &x.name == name)
                        && matches!(
                            s.attribute,
                            FlowMin
                                | FlowMax
                                | PowerMin
                                | PowerMax
                                | ForcedOn
                                | Maintenance
                                | Power
                                | Discharge
                                | PowerRampUp
                                | PowerRampDown
                                | FlowRampUp
                                | FlowRampDown
                        )
                }
                Object::Reservoir(name) => {
                    system.reservoirs.iter().any(|x| &x.name == name)
                        && matches!(
                            s.attribute,
                            StorageMin
                                | StorageMax
                                | Inflow
                                | VolumeRampUp
                                | VolumeRampDown
                                | LevelRampUp
                                | LevelRampDown
                        )
                }
                Object::Tunnel(name) => {
                    system.tunnels.iter().any(|x| &x.name == name)
                        && matches!(s.attribute, Capacity | Opening)
                }
                Object::TunnelInterlock(name) => {
                    system.tunnel_interlocks.iter().any(|r| &r.name == name)
                        && matches!(s.attribute, MinOpen | MaxOpen)
                }
                Object::FlowRequirement(name) => {
                    case.flow_requirements.iter().any(|r| &r.name == name)
                        && matches!(s.attribute, FlowMin | FlowMax | Inflow)
                }
                Object::River(name) => {
                    system.rivers.iter().any(|x| {
                        &x.name == name
                            && (!x.law.is_delta()
                                || matches!(
                                    s.attribute,
                                    Capacity | FlowRampUp | FlowRampDown | FlowChangeCost
                                ))
                            && (system.reservoirs.iter().any(|r| r.name == x.source)
                                || !matches!(s.attribute, GateMin | GateMax))
                    }) && matches!(
                        s.attribute,
                        Capacity
                            | Inflow
                            | GateMin
                            | GateMax
                            | MinRelease
                            | MinArrival
                            | FlowRampUp
                            | FlowRampDown
                            | FlowChangeCost
                    )
                }
            };
            if !allowed {
                return Err(Error::Input(format!(
                    "unsupported operating attribute {:?} on {:?}",
                    s.attribute, s.object
                )));
            }
            if let Object::Tunnel(name) = &s.object
                && s.attribute == Opening
            {
                let e = system.tunnels.iter().find(|e| &e.name == name).unwrap();
                if !e.gate_positions.is_empty()
                    && s.values.iter().any(|v| !e.gate_positions.contains(v))
                {
                    return Err(Error::Input(format!(
                        "{name}: opening schedule must select declared gate positions"
                    )));
                }
            }
            let valid = if matches!(s.object, Object::FlowRequirement(_))
                && matches!(s.attribute, FlowMin | FlowMax)
            {
                true
            } else {
                match s.attribute {
                    MinOpen | MaxOpen => s.values.iter().all(|v| *v >= 0. && v.fract() == 0.),
                    ForcedOn => s.values.iter().all(|v| [-1., 0., 1.].contains(v)),
                    Maintenance => s.values.iter().all(|v| [0., 1.].contains(v)),
                    Opening | GateMin | GateMax => s.values.iter().all(|v| (0. ..=1.).contains(v)),
                    _ => s.values.iter().all(|v| *v >= 0.),
                }
            };
            if !valid {
                return Err(Error::Input(
                    "invalid operating bound, rate or state".into(),
                ));
            }
            let static_range = match &s.object {
                Object::Unit(name) => {
                    let u = system.units.iter().find(|u| &u.name == name).unwrap();
                    match s.attribute {
                        FlowMin => Some((u.flow.min, f64::INFINITY)),
                        FlowMax => Some((0., u.flow.max)),
                        PowerMin => Some((u.power.min, f64::INFINITY)),
                        PowerMax => Some((0., u.power.max)),
                        _ => None,
                    }
                }
                Object::Reservoir(name) => {
                    let r = system.reservoirs.iter().find(|r| &r.name == name).unwrap();
                    match s.attribute {
                        StorageMin => Some((r.storage.min, f64::INFINITY)),
                        StorageMax => Some((0., r.storage.max)),
                        _ => None,
                    }
                }
                Object::Tunnel(name) if s.attribute == Capacity => Some((
                    0.,
                    system
                        .tunnels
                        .iter()
                        .find(|r| &r.name == name)
                        .unwrap()
                        .capacity,
                )),
                Object::TunnelInterlock(name) => Some((
                    0.,
                    system
                        .tunnel_interlocks
                        .iter()
                        .find(|r| &r.name == name)
                        .unwrap()
                        .tunnels
                        .len() as f64,
                )),
                Object::FlowRequirement(_) => None,
                Object::River(name) => {
                    let r = system.rivers.iter().find(|r| &r.name == name).unwrap();
                    match s.attribute {
                        MinRelease => Some((0., r.capacity)),
                        MinArrival => Some((0., f64::INFINITY)),
                        Capacity => Some((0., r.capacity)),
                        GateMin => Some((r.gate.min, 1.)),
                        GateMax => Some((0., r.gate.max)),
                        _ => None,
                    }
                }
                Object::Plant(name) => {
                    let p = system.plants.iter().find(|p| &p.name == name).unwrap();
                    match s.attribute {
                        FlowMin => Some((p.flow_min, f64::INFINITY)),
                        PowerMin => Some((p.power_min, f64::INFINITY)),
                        FlowMax => Some((
                            0.,
                            p.flow_max.unwrap_or_else(|| {
                                system
                                    .units
                                    .iter()
                                    .filter(|u| u.plant == *name)
                                    .map(|u| u.flow.max)
                                    .sum()
                            }),
                        )),
                        PowerMax => Some((
                            0.,
                            p.power_max.unwrap_or_else(|| {
                                system
                                    .units
                                    .iter()
                                    .filter(|u| u.plant == *name)
                                    .map(|u| u.power.max)
                                    .sum()
                            }),
                        )),
                        _ => None,
                    }
                }
                _ => None,
            };
            if static_range.is_some_and(|(lo, hi)| s.values.iter().any(|v| *v < lo || *v > hi)) {
                return Err(Error::Input(
                    "operating bounds may only tighten the static domain, including future data"
                        .into(),
                ));
            }
            series.entry(&s.object).or_default().insert(s.attribute, s);
        }
        let mut units = Vec::new();
        for u in &system.units {
            let object = Object::Unit(u.name.clone());
            let m = series.get(&object);
            let mut intervals: Vec<MachineInterval> = Vec::new();
            for t in 0..nt {
                let time = case.grid[t];
                let b = time + (case.grid[t + 1] - time) / 2.;
                let a = if t == 0 {
                    time - u.operation.initial_interval_hours / 2.
                } else {
                    (case.grid[t - 1] + time) / 2.
                };
                let flow = Bounds::new(
                    value(m, FlowMin, time, Some(u.flow.min), false).unwrap(),
                    value(m, FlowMax, time, Some(u.flow.max), false).unwrap(),
                );
                let power = Bounds::new(
                    value(m, PowerMin, time, Some(u.power.min), false).unwrap(),
                    value(m, PowerMax, time, Some(u.power.max), false).unwrap(),
                );
                if flow.min < u.flow.min
                    || flow.max > u.flow.max
                    || flow.min > flow.max
                    || power.min < u.power.min
                    || power.max > u.power.max
                    || power.min > power.max
                {
                    return Err(Error::Input(format!(
                        "{}: operating bounds may only tighten the static domain",
                        u.name
                    )));
                }
                let fixed_power = value(m, Power, time, None, false);
                let fixed_flow = value(m, Discharge, time, None, false);
                let mut state = forced(m, time)?;
                if fixed_power.is_some_and(|p| p > 0.) || fixed_flow.is_some_and(|q| q > 0.) {
                    if state == Some(false) {
                        return Err(Error::Input(format!(
                            "{}: positive output requested while stopped",
                            u.name
                        )));
                    }
                    state = Some(true);
                }
                if let Some(fixed) = case.commitment.get(&u.name) {
                    if state.is_some_and(|on| on != fixed[t]) {
                        return Err(Error::Input(format!(
                            "{}: prescribed state conflicts with operating data",
                            u.name
                        )));
                    }
                    state = Some(fixed[t]);
                }
                if fixed_power.is_some_and(|p| p != 0. && !power.contains(p))
                    || fixed_flow.is_some_and(|q| q != 0. && !flow.contains(q))
                {
                    return Err(Error::Input(format!(
                        "{}: requested output outside operating range",
                        u.name
                    )));
                }
                let stopping_power_min = if t == 0 {
                    value(m, PowerMin, time, Some(u.power.min), true).unwrap()
                } else {
                    intervals[t - 1].power.min
                };
                let gross_head_min = value(m, GrossHeadMin, time, u.gross_head_min, false);
                let gross_head_max = value(m, GrossHeadMax, time, u.gross_head_max, false);
                if matches!((gross_head_min,gross_head_max),(Some(a),Some(b)) if a>b) {
                    return Err(Error::Input(format!(
                        "{}: inconsistent gross head limits",
                        u.name
                    )));
                }
                intervals.push(MachineInterval {
                    gross_head_min,
                    gross_head_max,
                    startup_gross_head_min: value(
                        m,
                        StartupGrossHeadMin,
                        time,
                        u.startup_gross_head_min,
                        false,
                    ),
                    flow,
                    power,
                    forced: state,
                    fixed_power,
                    fixed_flow,
                    startup_cost: value(
                        m,
                        StartupCost,
                        time,
                        Some(u.operation.startup_cost),
                        false,
                    )
                    .unwrap(),
                    shutdown_cost: value(
                        m,
                        ShutdownCost,
                        time,
                        Some(u.operation.shutdown_cost),
                        false,
                    )
                    .unwrap(),
                    stopping_power_min,
                    flow_change: ChangeLimit {
                        up: change(m, FlowRampUp, a, b, u.operation.flow_ramp.up),
                        down: change(m, FlowRampDown, a, b, u.operation.flow_ramp.down),
                    },
                    power_change: ChangeLimit {
                        up: change(m, PowerRampUp, a, b, u.operation.power_ramp.up),
                        down: change(m, PowerRampDown, a, b, u.operation.power_ramp.down),
                    },
                });
            }
            units.push(intervals);
        }
        let mut plants = Vec::new();
        let mut plant_initial_flow = Vec::new();
        let mut plant_initial_power = Vec::new();
        for (i, p) in system.plants.iter().enumerate() {
            let object = Object::Plant(p.name.clone());
            let m = series.get(&object);
            let mut intervals = Vec::new();
            let qmax = p.flow_max.unwrap_or_else(|| {
                topology.plant_units[i]
                    .iter()
                    .map(|&g| system.units[g].flow.max)
                    .sum()
            });
            let pmax = p.power_max.unwrap_or_else(|| {
                topology.plant_units[i]
                    .iter()
                    .map(|&g| system.units[g].power.max)
                    .sum()
            });
            let op = &topology.plant_operation[i];
            let history = |power: bool| -> Result<(Option<f64>, f64), Error> {
                let (explicit, up, down, up_attribute, down_attribute) = if power {
                    (
                        p.operation.initial_power,
                        op.power_ramp.up,
                        op.power_ramp.down,
                        PowerRampUp,
                        PowerRampDown,
                    )
                } else {
                    (
                        p.operation.initial_flow,
                        op.flow_ramp.up,
                        op.flow_ramp.down,
                        FlowRampUp,
                        FlowRampDown,
                    )
                };
                if explicit.is_some() {
                    return Ok((explicit, p.operation.initial_interval_hours));
                }
                let enabled = up.is_some()
                    || down.is_some()
                    || m.is_some_and(|m| {
                        m.contains_key(&up_attribute) || m.contains_key(&down_attribute)
                    });
                if !enabled {
                    return Ok((None, p.operation.initial_interval_hours));
                }
                let values: Option<Vec<_>> = topology.plant_units[i]
                    .iter()
                    .map(|&g| {
                        if power {
                            system.units[g].operation.initial_power
                        } else {
                            system.units[g].operation.initial_flow
                        }
                    })
                    .collect();
                let Some(values) = values else {
                    return Ok((None, p.operation.initial_interval_hours));
                };
                let hours = system.units[topology.plant_units[i][0]]
                    .operation
                    .initial_interval_hours;
                if topology.plant_units[i]
                    .iter()
                    .any(|&g| system.units[g].operation.initial_interval_hours != hours)
                {
                    return Err(Error::Input(format!(
                        "{}: supply aggregate history when observation durations differ",
                        p.name
                    )));
                }
                Ok((Some(values.iter().sum()), hours))
            };
            let (initial_q, q_hours) = history(false)?;
            let (initial_p, p_hours) = history(true)?;
            plant_initial_flow.push(initial_q);
            plant_initial_power.push(initial_p);
            for t in 0..nt {
                let time = case.grid[t];
                let b = time + (case.grid[t + 1] - time) / 2.;
                let a = if t == 0 {
                    time - q_hours / 2.
                } else {
                    (case.grid[t - 1] + time) / 2.
                };
                let flow = Bounds::new(
                    value(m, FlowMin, time, Some(p.flow_min), false).unwrap(),
                    value(m, FlowMax, time, Some(qmax), false).unwrap(),
                );
                let power = Bounds::new(
                    value(m, PowerMin, time, Some(p.power_min), false).unwrap(),
                    value(m, PowerMax, time, Some(pmax), false).unwrap(),
                );
                if flow.min < p.flow_min
                    || flow.max > qmax
                    || flow.min > flow.max
                    || power.min < p.power_min
                    || power.max > pmax
                    || power.min > power.max
                {
                    return Err(Error::Input(format!(
                        "{}: invalid plant operating bounds",
                        p.name
                    )));
                }
                let power_a = if t == 0 { time - p_hours / 2. } else { a };
                let fixed_power = value(m, Power, time, None, false);
                let fixed_flow = value(m, Discharge, time, None, false);
                let mut state = forced(m, time)?;
                if fixed_power.is_some_and(|v| v > 0.) || fixed_flow.is_some_and(|v| v > 0.) {
                    if state == Some(false) {
                        return Err(Error::Input(format!(
                            "{}: positive output requested while stopped",
                            p.name
                        )));
                    }
                    state = Some(true);
                }
                if fixed_power.is_some_and(|v| v != 0. && !power.contains(v))
                    || fixed_flow.is_some_and(|v| v != 0. && !flow.contains(v))
                {
                    return Err(Error::Input(format!(
                        "{}: requested output outside plant limits",
                        p.name
                    )));
                }
                let members = &topology.plant_units[i];
                if members
                    .iter()
                    .all(|&g| case.commitment.contains_key(&system.units[g].name))
                {
                    let prescribed = members
                        .iter()
                        .any(|&g| case.commitment[&system.units[g].name][t]);
                    if state.is_some_and(|state| state != prescribed) {
                        return Err(Error::Input(format!(
                            "{}: prescribed states conflict with plant operation",
                            p.name
                        )));
                    }
                }
                intervals.push(PlantInterval {
                    flow,
                    power,
                    forced: state,
                    fixed_power,
                    fixed_flow,
                    flow_change: ChangeLimit {
                        up: change(m, FlowRampUp, a, b, op.flow_ramp.up),
                        down: change(m, FlowRampDown, a, b, op.flow_ramp.down),
                    },
                    power_change: ChangeLimit {
                        up: change(m, PowerRampUp, power_a, b, op.power_ramp.up),
                        down: change(m, PowerRampDown, power_a, b, op.power_ramp.down),
                    },
                });
            }
            plants.push(intervals);
        }
        let mut reservoirs = Vec::new();
        for (r, reservoir) in system.reservoirs.iter().enumerate() {
            let object = Object::Reservoir(reservoir.name.clone());
            let m = series.get(&object);
            let mut storage = Vec::new();
            for &time in &case.grid {
                let min = value(m, StorageMin, time, Some(reservoir.storage.min), false).unwrap();
                let max = value(m, StorageMax, time, Some(reservoir.storage.max), false).unwrap();
                let bound = if time == case.grid[0] {
                    Bounds::new(min, max)
                } else {
                    Bounds::new(
                        min.max(
                            value(m, StorageMin, time, Some(reservoir.storage.min), true).unwrap(),
                        ),
                        max.min(
                            value(m, StorageMax, time, Some(reservoir.storage.max), true).unwrap(),
                        ),
                    )
                };
                if bound.min < reservoir.storage.min
                    || bound.max > reservoir.storage.max
                    || bound.min > bound.max
                {
                    return Err(Error::Input(format!(
                        "{}: invalid storage limits",
                        reservoir.name
                    )));
                }
                storage.push(bound);
            }
            if !storage[0].contains(reservoir.initial_storage) {
                return Err(Error::Input(format!(
                    "{}: initial volume outside operating limits",
                    reservoir.name
                )));
            }
            let mut interval_storage = Vec::new();
            let mut inflow = Vec::new();
            let mut volume_change = Vec::new();
            let mut level_change = Vec::new();
            for t in 0..nt {
                let (a, b) = (case.grid[t], case.grid[t + 1]);
                interval_storage.push(Bounds::new(
                    value(m, StorageMin, a, Some(reservoir.storage.min), false).unwrap(),
                    value(m, StorageMax, a, Some(reservoir.storage.max), false).unwrap(),
                ));
                inflow.push(value(m, Inflow, a, Some(case.inflow(r, t)), false).unwrap());
                volume_change.push(ChangeLimit {
                    up: change(m, VolumeRampUp, a, b, reservoir.operation.volume_ramp.up),
                    down: change(
                        m,
                        VolumeRampDown,
                        a,
                        b,
                        reservoir.operation.volume_ramp.down,
                    ),
                });
                level_change.push(ChangeLimit {
                    up: change(m, LevelRampUp, a, b, reservoir.operation.level_ramp.up),
                    down: change(m, LevelRampDown, a, b, reservoir.operation.level_ramp.down),
                });
            }
            reservoirs.push(ReservoirData {
                storage,
                interval_storage,
                inflow,
                volume_change,
                level_change,
            });
        }
        let mut tunnels = Vec::new();
        for tunnel in &system.tunnels {
            let object = Object::Tunnel(tunnel.name.clone());
            let m = series.get(&object);
            let mut intervals = Vec::new();
            for &time in &case.grid[..nt] {
                let capacity = value(m, Capacity, time, Some(tunnel.capacity), false).unwrap();
                if capacity > tunnel.capacity {
                    return Err(Error::Input(format!(
                        "{}: operating capacity exceeds static domain",
                        tunnel.name
                    )));
                }
                let permitted = if tunnel.gate_positions.is_empty() {
                    vec![value(m, Opening, time, Some(tunnel.opening), false).unwrap()]
                } else if let Some(opening) = value(m, Opening, time, None, false) {
                    vec![opening]
                } else {
                    tunnel.gate_positions.clone()
                };
                let mut positions = Vec::new();
                for opening in permitted {
                    let resistance = if opening == 0. {
                        tunnel.loss_factor
                    } else if let Some(c) = &tunnel.gate_loss_factor_curve {
                        tunnel.loss_factor + c.value(opening)?
                    } else {
                        tunnel.loss_factor / opening
                    };
                    if !resistance.is_finite() {
                        return Err(Error::Input(format!(
                            "{}: nonfinite tunnel resistance",
                            tunnel.name
                        )));
                    }
                    positions.push(TunnelPosition {
                        opening,
                        resistance,
                    });
                }
                intervals.push(TunnelInterval {
                    capacity,
                    positions,
                });
            }
            tunnels.push(intervals);
        }
        let mut river_gate = Vec::new();
        let mut river_capacity = Vec::new();
        let mut river_inflow = Vec::new();
        let mut river_min_release = Vec::new();
        let mut river_min_arrival = Vec::new();
        let mut river_flow_change = Vec::new();
        let mut river_change_cost = Vec::new();
        for river in &system.rivers {
            river.operation.validate(&case.grid)?;
            let object = Object::River(river.name.clone());
            let m = series.get(&object);
            let mut gates = Vec::new();
            let mut cap = Vec::new();
            let mut inflow = Vec::new();
            let mut minimum = Vec::new();
            let mut arrival = Vec::new();
            let mut ramps = Vec::new();
            let mut costs = Vec::new();
            let prior = river.previous_release();
            if prior.is_some_and(|(q, _)| q < 0.) && !river.law.is_delta() {
                return Err(Error::Input(
                    "directed river initial releases must be nonnegative".into(),
                ));
            }
            for (t, &time) in case.grid[..nt].iter().enumerate() {
                let c = value(m, Capacity, time, Some(river.capacity), false).unwrap();
                let q = value(m, Inflow, time, Some(river.inflow), false).unwrap();
                let minimum_release = value(m, MinRelease, time, Some(0.), false).unwrap();
                if c > river.capacity || q > c || minimum_release > c {
                    return Err(Error::Input(format!(
                        "{}: invalid river capacity/inflow",
                        river.name
                    )));
                }
                let lo = value(m, GateMin, time, Some(river.gate.min), false).unwrap();
                let hi = value(m, GateMax, time, Some(river.gate.max), false).unwrap();
                if lo < river.gate.min
                    || hi > river.gate.max
                    || lo > hi
                    || (matches!(river.law, crate::RiverLaw::Weir { .. }) && hi < 1.)
                {
                    return Err(Error::Input(format!(
                        "{}: invalid river gate limits",
                        river.name
                    )));
                }
                gates.push(if matches!(river.law, crate::RiverLaw::Weir { .. }) {
                    Bounds::new(1., 1.)
                } else {
                    Bounds::new(lo, hi)
                });
                minimum.push(minimum_release);
                arrival.push(
                    value(
                        m,
                        MinArrival,
                        time,
                        Some(river.operation.minimum_arrival),
                        false,
                    )
                    .unwrap(),
                );
                let a = if t == 0 {
                    time - prior.map_or(river.operation.initial_interval_hours, |(_, dt)| dt) / 2.
                } else {
                    case.grid[t - 1] + (time - case.grid[t - 1]) / 2.
                };
                let b = time + (case.grid[t + 1] - time) / 2.;
                let ramp = ChangeLimit {
                    up: change(m, FlowRampUp, a, b, river.operation.flow_ramp.up),
                    down: change(m, FlowRampDown, a, b, river.operation.flow_ramp.down),
                };
                if !a.is_finite()
                    || !b.is_finite()
                    || b <= a
                    || ramp.up.into_iter().chain(ramp.down).any(|v| !v.is_finite())
                {
                    return Err(Error::Input(
                        "river ramp observation times and allowances must be finite".into(),
                    ));
                }
                ramps.push(ramp);
                costs.push(
                    value(
                        m,
                        FlowChangeCost,
                        time,
                        Some(river.operation.flow_change_cost),
                        false,
                    )
                    .unwrap(),
                );
                cap.push(c);
                inflow.push(q);
            }
            river_gate.push(gates);
            river_capacity.push(cap);
            river_inflow.push(inflow);
            river_min_release.push(minimum);
            river_min_arrival.push(arrival);
            river_flow_change.push(ramps);
            river_change_cost.push(costs);
        }
        let mut observed_flow = Vec::new();
        for r in &case.flow_requirements {
            let object = Object::FlowRequirement(r.name.clone());
            let m = series.get(&object);
            let mut intervals = Vec::new();
            for &time in &case.grid[..nt] {
                let bounds = Bounds::new(
                    value(m, FlowMin, time, Some(r.min_flow), false).unwrap(),
                    value(m, FlowMax, time, r.max_flow, false).unwrap_or(f64::INFINITY),
                );
                if bounds.min > bounds.max {
                    return Err(Error::Input(format!(
                        "{}: inconsistent observed flow bounds",
                        r.name
                    )));
                }
                intervals.push(crate::flow_requirements::FlowInterval {
                    inflow: value(m, Inflow, time, Some(r.inflow), false).unwrap(),
                    bounds,
                });
            }
            observed_flow.push(intervals);
        }
        Ok(Self {
            observed_flow,
            units,
            plants,
            plant_initial_flow,
            plant_initial_power,
            reservoirs,
            tunnels,
            river_capacity,
            river_inflow,
            river_min_release,
            river_min_arrival,
            river_flow_change,
            river_change_cost,
            river_gate,
        })
    }
}
