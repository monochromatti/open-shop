//! Conservative transport through river networks.
//!
//! A release interval is translated in time, not spread over the next control
//! interval. At a river junction we retain every arriving interval boundary.
//! The resulting sparse linear maps are shared by commitment and dispatch.

use crate::{Error, ScheduleCase, model::Topology};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug, Default)]
pub(crate) struct LinearForm {
    pub constant: f64,
    pub coefficients: BTreeMap<usize, f64>,
}

impl LinearForm {
    pub fn variable(index: usize, coefficient: f64) -> Self {
        Self {
            constant: 0.,
            coefficients: BTreeMap::from([(index, coefficient)]),
        }
    }
    pub fn scaled(mut self, scale: f64) -> Self {
        self.constant *= scale;
        for a in self.coefficients.values_mut() {
            *a *= scale;
        }
        self
    }

    pub(crate) fn add_scaled(&mut self, other: &Self, scale: f64) {
        if scale == 0. {
            return;
        }
        self.constant += scale * other.constant;
        for (&i, &a) in &other.coefficients {
            let value = self.coefficients.get(&i).copied().unwrap_or(0.) + scale * a;
            if value == 0. {
                self.coefficients.remove(&i);
            } else {
                self.coefficients.insert(i, value);
            }
        }
    }
    pub fn value(&self, x: &[f64]) -> f64 {
        self.constant
            + self
                .coefficients
                .iter()
                .map(|(&i, &a)| a * x[i])
                .sum::<f64>()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ReleaseSeries {
    pub grid: Vec<f64>,
    pub rates: Vec<LinearForm>,
}

impl ReleaseSeries {
    fn at(&self, time: f64) -> LinearForm {
        let k = self.grid.partition_point(|&t| t <= time);
        if k == 0 || k >= self.grid.len() {
            LinearForm::default()
        } else {
            self.rates[k - 1].clone()
        }
    }
    /// Integral in (m³/s)·hours. Conversion to m³ is explicit at call sites.
    pub fn integral(&self, start: f64, end: f64) -> LinearForm {
        let mut volume = LinearForm::default();
        for (k, pair) in self.grid.windows(2).enumerate() {
            let overlap = end.min(pair[1]) - start.max(pair[0]);
            if overlap > 0. {
                volume.add_scaled(&self.rates[k], overlap);
            }
        }
        volume
    }
}

pub(crate) struct TransferEquation {
    pub flow: LinearForm,
    pub fraction: std::sync::Arc<crate::river_transport::ReleaseFraction>,
}
pub(crate) struct Routing {
    pub primary_columns: usize,
    /// Local released-transfer equations, in upstream-to-downstream order.
    pub transfers: Vec<TransferEquation>,
    /// Reservoir-source reach indices, in decision-variable order.
    pub outlets: Vec<usize>,
    /// Columns are unit flows, tunnel flows, then reservoir river releases;
    /// each block is interval-major. Shared reference fluxes follow these inputs.
    pub releases: Vec<ReleaseSeries>,
    /// Mean downstream arrival rate, reach by control interval, m³/s.
    pub arrivals: Vec<Vec<LinearForm>>,
    pub arrival_limits: Vec<crate::flow_requirements::FlowRow>,
    /// River inventory at each control vertex, m³.
    pub inventory: Vec<Vec<LinearForm>>,
    /// River contributions to reservoir balances, m³/s.
    pub reservoir_net: Vec<Vec<LinearForm>>,
    /// Interior storage relative to straight interpolation of interval endpoints.
    /// The correction is a volume (m³), and vanishes at both endpoints.
    pub storage_shape: Vec<Vec<Vec<(f64, LinearForm)>>>,
}

impl Routing {
    pub fn new(
        case: &ScheduleCase,
        topology: &Topology,
        operating: &crate::operating_data::OperatingData,
    ) -> Result<Self, Error> {
        let s = &case.system;
        let nt = case.prices.len();
        let (start, end) = (case.grid[0], case.grid[nt]);
        let reservoirs: HashMap<_, _> = s
            .reservoirs
            .iter()
            .enumerate()
            .map(|(i, r)| (r.name.as_str(), i))
            .collect();
        let mut names = HashMap::new();
        let mut outlets = Vec::new();
        let mut outlet_index = vec![None; s.rivers.len()];
        let mut upstream = vec![Vec::new(); s.rivers.len()];
        let mut outgoing: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
        let mut history = Vec::new();
        for (i, r) in s.rivers.iter().enumerate() {
            if r.name.is_empty()
                || names.insert(r.name.as_str(), i).is_some()
                || r.source.is_empty()
                || r.target.is_empty()
                || r.source == r.target
                || s.reservoirs.iter().any(|n| n.name == r.name)
                || s.boundaries.iter().any(|n| n.name == r.name)
                || s.junctions
                    .iter()
                    .any(|n| n.name == r.source || n.name == r.target || n.name == r.name)
                || s.boundaries.iter().any(|n| n.name == r.source)
                || !r.capacity.is_finite()
                || r.capacity <= 0.
                || !r.allocation.is_finite()
                || !(0. ..=1.).contains(&r.allocation)
                || !r.inflow.is_finite()
                || r.inflow < 0.
                || !r.water_value.is_finite()
            {
                return Err(Error::Input(format!(
                    "{}: invalid river endpoints, capacity, delay or allocation",
                    r.name
                )));
            }
            r.law.validate(r, s)?;
            r.transport.validate(r.capacity, end)?;
            let empty = r.history_grid.is_empty() && r.history_release.is_empty();
            if !empty
                && (r.history_grid.len() != r.history_release.len() + 1
                    || r.history_grid.len() < 2
                    || !r.history_grid.iter().all(|t| t.is_finite() && *t <= start)
                    || !r.history_grid.windows(2).all(|w| w[1] > w[0])
                    || !r.history_release.iter().all(|q| q.is_finite() && *q >= 0.))
            {
                return Err(Error::Input(format!(
                    "{}: invalid historical river releases",
                    r.name
                )));
            }
            if let crate::RiverTransport::Distributed { curves } = &r.transport
                && curves.len() > 1
                && r.history_release
                    .iter()
                    .any(|q| *q > curves.last().unwrap().reference_flow())
            {
                return Err(Error::Input(format!(
                    "{}: historical release exceeds delay reference domain",
                    r.name
                )));
            }
            history.push(ReleaseSeries {
                grid: r.history_grid.clone(),
                rates: r
                    .history_release
                    .iter()
                    .map(|&q| LinearForm {
                        constant: q,
                        ..Default::default()
                    })
                    .collect(),
            });
            if reservoirs.contains_key(r.source.as_str()) {
                if r.allocation != 1. {
                    return Err(Error::Input(format!(
                        "{}: allocation applies only at river junctions",
                        r.name
                    )));
                }
                outlet_index[i] = Some(outlets.len());
                outlets.push(i);
            } else {
                outgoing.entry(&r.source).or_default().push(i);
                upstream[i] = s
                    .rivers
                    .iter()
                    .enumerate()
                    .filter_map(|(j, u)| (u.target == r.source).then_some(j))
                    .collect();
            }
        }
        for (name, edges) in &outgoing {
            if (edges.iter().map(|&i| s.rivers[i].allocation).sum::<f64>() - 1.).abs() > 1e-12 {
                return Err(Error::Input(format!(
                    "{name}: outgoing river allocations must sum to one"
                )));
            }
        }
        for r in &s.rivers {
            if !reservoirs.contains_key(r.target.as_str())
                && !s.boundaries.iter().any(|b| b.name == r.target)
                && !outgoing.contains_key(r.target.as_str())
            {
                return Err(Error::Input(format!(
                    "{}: river target {} has no outgoing reach",
                    r.name, r.target
                )));
            }
        }
        for (i, r) in s.rivers.iter().enumerate() {
            if outlet_index[i].is_none()
                && upstream[i].is_empty()
                && operating.river_inflow[i].iter().all(|q| *q == 0.)
                && !topology.unit_outfalls.contains(&Some(i))
                && !topology.tunnel_outfalls.contains(&Some(i))
            {
                return Err(Error::Input(format!(
                    "{}: river source {} has no incoming reach, outfall or inflow",
                    r.name, r.source
                )));
            }
        }
        let mut order = Vec::new();
        let mut visited = vec![false; s.rivers.len()];
        while order.len() < s.rivers.len() {
            let before = order.len();
            for i in 0..s.rivers.len() {
                if !visited[i] && upstream[i].iter().all(|&j| visited[j]) {
                    order.push(i);
                    visited[i] = true;
                }
            }
            if order.len() == before {
                return Err(Error::Input(
                    "river cycles require an intervening reservoir".into(),
                ));
            }
        }
        if s.rivers.iter().any(|r| r.transport.fixed_delay().is_none()) {
            return Self::distributed(
                case,
                topology,
                operating,
                outlets,
                outlet_index,
                upstream,
                order,
                history,
            );
        }
        let mut releases: Vec<_> = (0..s.rivers.len())
            .map(|_| ReleaseSeries {
                grid: Vec::new(),
                rates: Vec::new(),
            })
            .collect();
        for i in order {
            let river = &s.rivers[i];
            let mut grid = case.grid.clone();
            for &j in &upstream[i] {
                for &t in releases[j].grid.iter().chain(&history[j].grid) {
                    let arrival = t + s.rivers[j].transport.fixed_delay().unwrap();
                    if arrival > start && arrival < end {
                        grid.push(arrival);
                    }
                }
            }
            grid.sort_by(f64::total_cmp);
            grid.dedup();
            let rates = grid
                .windows(2)
                .map(|w| {
                    let time = w[0] + (w[1] - w[0]) / 2.;
                    let t = case.grid.partition_point(|&v| v <= time) - 1;
                    let mut rate = LinearForm {
                        constant: operating.river_inflow[i][t],
                        ..Default::default()
                    };
                    if let Some(k) = outlet_index[i] {
                        rate.coefficients
                            .insert((s.units.len() + s.tunnels.len()) * nt + k * nt + t, 1.);
                    }
                    for &j in &upstream[i] {
                        let released = time - s.rivers[j].transport.fixed_delay().unwrap();
                        let previous = if released < start {
                            &history[j]
                        } else {
                            &releases[j]
                        };
                        rate.add_scaled(&previous.at(released), river.allocation);
                    }
                    for (u, outfall) in topology.unit_outfalls.iter().enumerate() {
                        if *outfall == Some(i) {
                            rate.coefficients.insert(u * nt + t, 1.);
                        }
                    }
                    for (e, outfall) in topology.tunnel_outfalls.iter().enumerate() {
                        if *outfall == Some(i) {
                            rate.coefficients.insert((s.units.len() + e) * nt + t, 1.);
                        }
                    }
                    rate
                })
                .collect();
            releases[i] = ReleaseSeries { grid, rates };
        }
        let mut arrivals = Vec::new();
        let mut arrival_limits = Vec::new();
        let mut inventory = Vec::new();
        let mut reservoir_net = vec![vec![LinearForm::default(); nt]; s.reservoirs.len()];
        for (i, r) in s.rivers.iter().enumerate() {
            // Fixed-delay reference partitions have one flux per release and
            // add no nonlinear transfer variables.
            let mut transfers = Vec::new();
            let flux = reference_fluxes(&r.transport, &releases[i], 0, &mut transfers);
            let past_flux = reference_fluxes(&r.transport, &history[i], 0, &mut transfers);
            arrival_limits.extend(crate::river_operation::arrival_limits(
                case,
                operating,
                i,
                (&releases[i], &flux),
                (&history[i], &past_flux),
            ));
            let arrived: Vec<_> = case
                .grid
                .windows(2)
                .map(|w| {
                    let mut volume = history[i].integral(
                        w[0] - r.transport.fixed_delay().unwrap(),
                        w[1] - r.transport.fixed_delay().unwrap(),
                    );
                    volume.add_scaled(
                        &releases[i].integral(
                            w[0] - r.transport.fixed_delay().unwrap(),
                            w[1] - r.transport.fixed_delay().unwrap(),
                        ),
                        1.,
                    );
                    let mut mean = LinearForm::default();
                    mean.add_scaled(&volume, 1. / (w[1] - w[0]));
                    mean
                })
                .collect();
            let transit = case
                .grid
                .iter()
                .map(|&time| {
                    let mut volume = LinearForm::default();
                    volume.add_scaled(
                        &history[i].integral(time - r.transport.fixed_delay().unwrap(), time),
                        3600.,
                    );
                    volume.add_scaled(
                        &releases[i].integral(time - r.transport.fixed_delay().unwrap(), time),
                        3600.,
                    );
                    volume
                })
                .collect();
            if let Some(&ridx) = reservoirs.get(r.target.as_str()) {
                for (t, row) in arrived.iter().enumerate() {
                    reservoir_net[ridx][t].add_scaled(row, 1.);
                }
            }
            if let Some(&ridx) = reservoirs.get(r.source.as_str()) {
                let k = outlet_index[i].unwrap();
                for (t, row) in reservoir_net[ridx].iter_mut().enumerate() {
                    *row.coefficients
                        .entry((s.units.len() + s.tunnels.len()) * nt + k * nt + t)
                        .or_default() -= 1.;
                }
            }
            arrivals.push(arrived);
            inventory.push(transit);
        }
        let mut storage_shape = vec![Vec::new(); s.reservoirs.len()];
        for (r, reservoir) in s.reservoirs.iter().enumerate() {
            let incoming: Vec<_> = s
                .rivers
                .iter()
                .enumerate()
                .filter_map(|(i, e)| (e.target == reservoir.name).then_some(i))
                .collect();
            for (t, w) in case.grid.windows(2).enumerate() {
                let mut times = vec![w[0], w[1]];
                for &i in &incoming {
                    for &release in releases[i].grid.iter().chain(&history[i].grid) {
                        let arrival = release + s.rivers[i].transport.fixed_delay().unwrap();
                        if arrival > w[0] && arrival < w[1] {
                            times.push(arrival);
                        }
                    }
                }
                times.sort_by(f64::total_cmp);
                times.dedup();
                let vertices = times
                    .iter()
                    .map(|&time| {
                        let fraction = (time - w[0]) / (w[1] - w[0]);
                        let mut correction = LinearForm::default();
                        if time > w[0] && time < w[1] {
                            for &i in &incoming {
                                let delay = s.rivers[i].transport.fixed_delay().unwrap();
                                correction.add_scaled(
                                    &history[i].integral(w[0] - delay, time - delay),
                                    3600.,
                                );
                                correction.add_scaled(
                                    &releases[i].integral(w[0] - delay, time - delay),
                                    3600.,
                                );
                                correction.add_scaled(&arrivals[i][t], -3600. * (time - w[0]));
                            }
                        }
                        (fraction, correction)
                    })
                    .collect();
                storage_shape[r].push(vertices);
            }
        }
        Ok(Self {
            primary_columns: (s.units.len() + s.tunnels.len() + outlets.len()) * nt,
            transfers: Vec::new(),
            outlets,
            releases,
            arrivals,
            arrival_limits,
            inventory,
            reservoir_net,
            storage_shape,
        })
    }
    /// Complete physical routing values in dependency order. This is also used
    /// to validate a solver result independently of its lifted transfer variables.
    pub fn complete_columns(&self, columns: &mut Vec<f64>) {
        debug_assert_eq!(columns.len(), self.primary_columns);
        for e in &self.transfers {
            columns.push(e.fraction.data(e.flow.value(columns)).value);
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn distributed(
        case: &ScheduleCase,
        topology: &Topology,
        operating: &crate::operating_data::OperatingData,
        outlets: Vec<usize>,
        outlet_index: Vec<Option<usize>>,
        upstream: Vec<Vec<usize>>,
        order: Vec<usize>,
        history: Vec<ReleaseSeries>,
    ) -> Result<Self, Error> {
        let s = &case.system;
        let nt = case.prices.len();
        let nr = s.rivers.len();
        let primary_columns = (s.units.len() + s.tunnels.len() + outlets.len()) * nt;
        let mut transfers = Vec::new();
        let mesh = if case.routing_grid.is_empty() {
            let mut grid = Vec::new();
            for w in case.grid.windows(2) {
                for k in 0..4 {
                    grid.push(w[0] + (w[1] - w[0]) * k as f64 / 4.);
                }
            }
            grid.push(*case.grid.last().unwrap());
            grid
        } else {
            case.routing_grid.clone()
        };
        let mut releases = vec![
            ReleaseSeries {
                grid: vec![],
                rates: vec![]
            };
            nr
        ];
        let mut mesh_arrival = vec![vec![]; nr];
        let mut arrivals = vec![vec![]; nr];
        let mut arrival_limits = Vec::new();
        let mut inventory = vec![vec![]; nr];
        for i in order {
            let r = &s.rivers[i];
            let grid = if upstream[i].is_empty() {
                case.grid.clone()
            } else {
                mesh.clone()
            };
            let mut rates = Vec::new();
            for (k, w) in grid.windows(2).enumerate() {
                let t = case
                    .grid
                    .partition_point(|v| *v <= w[0])
                    .saturating_sub(1)
                    .min(nt - 1);
                let mut q = LinearForm {
                    constant: operating.river_inflow[i][t],
                    ..Default::default()
                };
                if let Some(j) = outlet_index[i] {
                    q.coefficients
                        .insert((s.units.len() + s.tunnels.len()) * nt + j * nt + t, 1.);
                }
                for &j in &upstream[i] {
                    q.add_scaled(&mesh_arrival[j][k], r.allocation);
                }
                for (u, outfall) in topology.unit_outfalls.iter().enumerate() {
                    if *outfall == Some(i) {
                        *q.coefficients.entry(u * nt + t).or_default() += 1.;
                    }
                }
                for (e, outfall) in topology.tunnel_outfalls.iter().enumerate() {
                    if *outfall == Some(i) {
                        *q.coefficients
                            .entry((s.units.len() + e) * nt + t)
                            .or_default() += 1.;
                    }
                }
                rates.push(q);
            }
            releases[i] = ReleaseSeries { grid, rates };
            let flux =
                reference_fluxes(&r.transport, &releases[i], primary_columns, &mut transfers);
            let past_flux =
                reference_fluxes(&r.transport, &history[i], primary_columns, &mut transfers);
            arrival_limits.extend(crate::river_operation::arrival_limits(
                case,
                operating,
                i,
                (&releases[i], &flux),
                (&history[i], &past_flux),
            ));
            mesh_arrival[i] = mesh
                .windows(2)
                .map(|w| {
                    let observation = crate::river_transport::Observation::Arrival(w[0], w[1]);
                    let mut mean =
                        transport_series(&r.transport, &history[i], &past_flux, observation);
                    mean.add_scaled(
                        &transport_series(&r.transport, &releases[i], &flux, observation),
                        1.,
                    );
                    mean.scaled(1. / (w[1] - w[0]))
                })
                .collect();
            arrivals[i] = case
                .grid
                .windows(2)
                .map(|w| {
                    let mut mean = LinearForm::default();
                    for (k, cell) in mesh.windows(2).enumerate() {
                        if w[0] <= cell[0] && cell[1] <= w[1] {
                            mean.add_scaled(
                                &mesh_arrival[i][k],
                                (cell[1] - cell[0]) / (w[1] - w[0]),
                            );
                        }
                    }
                    mean
                })
                .collect();
            inventory[i] = case
                .grid
                .iter()
                .map(|&time| {
                    let observation = crate::river_transport::Observation::Transit(time);
                    let mut volume =
                        transport_series(&r.transport, &history[i], &past_flux, observation);
                    volume.add_scaled(
                        &transport_series(&r.transport, &releases[i], &flux, observation),
                        1.,
                    );
                    volume.scaled(3600.)
                })
                .collect();
        }
        let mut reservoir_net = vec![vec![LinearForm::default(); nt]; s.reservoirs.len()];
        let mut storage_shape = vec![vec![]; s.reservoirs.len()];
        for (r, reservoir) in s.reservoirs.iter().enumerate() {
            let incoming: Vec<_> = s
                .rivers
                .iter()
                .enumerate()
                .filter_map(|(i, e)| (e.target == reservoir.name).then_some(i))
                .collect();
            for &i in &incoming {
                for t in 0..nt {
                    reservoir_net[r][t].add_scaled(&arrivals[i][t], 1.);
                }
            }
            for (i, e) in s.rivers.iter().enumerate() {
                if e.source == reservoir.name {
                    let j = outlet_index[i].unwrap();
                    for (t, row) in reservoir_net[r].iter_mut().enumerate() {
                        *row.coefficients
                            .entry((s.units.len() + s.tunnels.len()) * nt + j * nt + t)
                            .or_default() -= 1.;
                    }
                }
            }
            for (t, w) in case.grid.windows(2).enumerate() {
                let mut vertices = vec![(0., LinearForm::default())];
                let mut correction = LinearForm::default();
                for (k, cell) in mesh.windows(2).enumerate() {
                    if cell[0] < w[0] || cell[1] > w[1] {
                        continue;
                    }
                    for &i in &incoming {
                        correction.add_scaled(&mesh_arrival[i][k], 3600. * (cell[1] - cell[0]));
                        correction.add_scaled(&arrivals[i][t], -3600. * (cell[1] - cell[0]));
                    }
                    let fraction = (cell[1] - w[0]) / (w[1] - w[0]);
                    vertices.push((
                        fraction,
                        if cell[1] == w[1] {
                            LinearForm::default()
                        } else {
                            correction.clone()
                        },
                    ));
                }
                storage_shape[r].push(vertices);
            }
        }
        Ok(Self {
            primary_columns,
            transfers,
            outlets,
            releases,
            arrivals,
            arrival_limits,
            inventory,
            reservoir_net,
            storage_shape,
        })
    }
}

/// Partition each release into reference-distribution fluxes. Only J−1
/// nonlinear fluxes are needed: the last equals total release minus the others.
/// Every arrival and inventory observation reuses this partition, so its size
/// grows with releases and reference curves rather than arrival/release pairs.
fn reference_fluxes(
    transport: &crate::RiverTransport,
    series: &ReleaseSeries,
    primary_columns: usize,
    transfers: &mut Vec<TransferEquation>,
) -> Vec<Vec<LinearForm>> {
    let references = match transport {
        crate::RiverTransport::Fixed { .. } => vec![0.],
        crate::RiverTransport::Distributed { curves } => curves
            .iter()
            .map(crate::DelayCurve::reference_flow)
            .collect(),
    };
    series
        .rates
        .iter()
        .map(|q| {
            let mut columns = Vec::new();
            let mut remainder = q.clone();
            for j in 0..references.len() - 1 {
                let mut coefficients = vec![0.; references.len()];
                coefficients[j] = 1.;
                let fraction = crate::river_transport::ReleaseFraction {
                    references: references.clone(),
                    coefficients,
                };
                let column = if q.coefficients.is_empty() {
                    LinearForm {
                        constant: fraction.data(q.constant).value,
                        ..Default::default()
                    }
                } else {
                    let index = primary_columns + transfers.len();
                    transfers.push(TransferEquation {
                        flow: q.clone(),
                        fraction: std::sync::Arc::new(fraction),
                    });
                    LinearForm::variable(index, 1.)
                };
                remainder.add_scaled(&column, -1.);
                columns.push(column);
            }
            columns.push(remainder);
            columns
        })
        .collect()
}

/// Integral in (m³/s)·hours, retaining each release's original reference blend.
pub(crate) fn transport_series(
    transport: &crate::RiverTransport,
    series: &ReleaseSeries,
    flux: &[Vec<LinearForm>],
    observation: crate::river_transport::Observation,
) -> LinearForm {
    let mut total = LinearForm::default();
    for (j, w) in series.grid.windows(2).enumerate() {
        let fraction = transport.fraction(observation, w[0], w[1]);
        for (column, &coefficient) in flux[j].iter().zip(&fraction.coefficients) {
            let duration = if matches!(
                observation,
                crate::river_transport::Observation::Point(_, _)
            ) {
                1.
            } else {
                w[1] - w[0]
            };
            total.add_scaled(column, coefficient * duration);
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Boundary, Bounds, HydroSystem, Reservoir, River, StorageCurve};

    fn reach(name: &str, source: &str, target: &str, delay: f64) -> River {
        River {
            operation: Default::default(),
            law: Default::default(),
            gate: Bounds::new(0., 1.),
            name: name.into(),
            source: source.into(),
            target: target.into(),
            transport: crate::RiverTransport::Fixed { delay },
            capacity: 100.,
            allocation: 1.,
            inflow: 0.,
            water_value: 0.,
            history_grid: Vec::new(),
            history_release: Vec::new(),
        }
    }
    fn case(rivers: Vec<River>) -> ScheduleCase {
        ScheduleCase {
            soft_constraints: vec![],
            flow_requirements: vec![],
            routing_grid: vec![],
            operations: Default::default(),
            system: HydroSystem {
                reservoirs: vec![Reservoir {
                    operation: Default::default(),
                    name: "Lake".into(),
                    storage: Bounds::new(0., 1e6),
                    initial_storage: 5e5,
                    level: StorageCurve::Polynomial {
                        reference_volume: 0.,
                        elevation: 100.,
                        slope: 0.,
                        curvature: 0.,
                    },
                    water_value: 0.,
                }],
                boundaries: vec![Boundary {
                    name: "Sea".into(),
                    head: 0.,
                }],
                rivers,
                ..Default::default()
            },
            grid: vec![0., 0.5, 1.5, 2.5],
            prices: vec![0.; 3],
            inflows: BTreeMap::new(),
            commitment: BTreeMap::new(),
        }
    }

    #[test]
    fn delayed_pulse_and_history_conserve_water_on_unequal_intervals() {
        let mut river = reach("River", "Lake", "Sea", 0.75);
        river.history_grid = vec![-1., 0.];
        river.history_release = vec![2.];
        let case = case(vec![river]);
        let routing = Routing::new(
            &case,
            &case.topology().unwrap(),
            &crate::operating_data::OperatingData::new(&case, &case.topology().unwrap()).unwrap(),
        )
        .unwrap();
        let flow = [4., 8., 0.];
        let arrivals: Vec<_> = routing.arrivals[0]
            .iter()
            .enumerate()
            .map(|(t, a)| a.value(&flow) * 3600. * (case.grid[t + 1] - case.grid[t]))
            .collect();
        for (&actual, expected) in arrivals.iter().zip([3600., 16200., 21600.]) {
            assert!((actual - expected).abs() < 1e-9);
        }
        let initial = routing.inventory[0][0].value(&flow);
        let terminal = routing.inventory[0][3].value(&flow);
        assert!((initial - 5400.).abs() < 1e-9);
        assert!((initial + 36000. - arrivals.iter().sum::<f64>() - terminal).abs() < 1e-9);
    }

    #[test]
    fn splitting_and_merging_preserves_pulses_and_total_volume() {
        let mut left = reach("Left", "Fork", "Sea", 0.25);
        left.allocation = 0.4;
        let mut right = reach("Right", "Fork", "Sea", 0.5);
        right.allocation = 0.6;
        // Deliberately provide reaches out of topological order.
        let case = case(vec![
            left,
            reach("A", "Lake", "Merge", 0.25),
            right,
            reach("Combined", "Merge", "Fork", 0.),
            reach("B", "Lake", "Merge", 0.75),
        ]);
        let routing = Routing::new(
            &case,
            &case.topology().unwrap(),
            &crate::operating_data::OperatingData::new(&case, &case.topology().unwrap()).unwrap(),
        )
        .unwrap();
        let flow = [4., 0., 8., 0., 6., 0.];
        let merged = &routing.releases[3];
        for (time, expected) in [(0.6, 4.), (0.9, 0.), (1.4, 6.), (2., 14.)] {
            assert!((merged.at(time).value(&flow) - expected).abs() < 1e-12);
        }
        let sea_arrivals: [f64; 2] = [0, 2].map(|r| {
            routing.arrivals[r]
                .iter()
                .enumerate()
                .map(|(t, a)| a.value(&flow) * 3600. * (case.grid[t + 1] - case.grid[t]))
                .sum()
        });
        let remaining: f64 = routing.inventory.iter().map(|r| r[3].value(&flow)).sum();
        assert!(
            (3600. * (4. * 0.5 + 8. + 6.) - sea_arrivals.iter().sum::<f64>() - remaining).abs()
                < 1e-8
        );
    }

    #[test]
    fn invalid_splits_and_river_cycles_are_rejected() {
        let mut a = reach("A", "Fork", "Sea", 1.);
        a.allocation = 0.4;
        let case = case(vec![a]);
        assert!(matches!(
            Routing::new(
                &case,
                &case.topology().unwrap(),
                &crate::operating_data::OperatingData::new(&case, &case.topology().unwrap())
                    .unwrap()
            ),
            Err(Error::Input(_))
        ));
        let case = super::tests::case(vec![
            reach("A", "One", "Two", 1.),
            reach("B", "Two", "One", 1.),
        ]);
        assert!(matches!(
            Routing::new(
                &case,
                &case.topology().unwrap(),
                &crate::operating_data::OperatingData::new(&case, &case.topology().unwrap())
                    .unwrap()
            ),
            Err(Error::Input(_))
        ));
    }

    #[test]
    fn a_returning_reach_does_not_erase_other_reservoir_fluxes() {
        let case = case(vec![
            reach("Return", "Bend", "Lake", 0.),
            reach("Outlet", "Lake", "Bend", 0.),
        ]);
        let routing = Routing::new(
            &case,
            &case.topology().unwrap(),
            &crate::operating_data::OperatingData::new(&case, &case.topology().unwrap()).unwrap(),
        )
        .unwrap();
        for row in &routing.reservoir_net[0] {
            assert!(row.value(&[4., 7., 2.]).abs() < 1e-12);
        }
    }
}
