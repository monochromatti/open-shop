//! Independent hydraulic reconstruction by damped Newton iteration.
//!
//! Adjustable machine flows and gate openings are controls. Storage, pressure
//! junction heads, tunnel flows, reservoir outlets and fixed-speed pump flows
//! follow from conservation and the constitutive laws. No operating bound is
//! imposed by clipping a physical state.

use crate::{Error, Schedule, SchedulingProblem};
use nalgebra::{Const, DMatrix, DVector, U1};
use num_dual::{DualNum, DualSVec64};
use serde::{Deserialize, Serialize};

/// Physical controls and numerical starting values, with objects as rows and
/// operating intervals as columns.
/// Adjustable unit discharge is in m³/s. For a running fixed-speed pump its
/// discharge entry is only an initial guess; the characteristic determines flow.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Controls {
    /// Optional node-by-interval head guesses, metres; empty uses the hydraulic
    /// initialization. Only junction entries are used. These choose a pressure
    /// solution when the equations leave it undetermined, and never replace a
    /// pressure or conservation equation.
    #[serde(default)]
    pub head_guess: Vec<Vec<f64>>,
    pub commitment: Vec<Vec<bool>>,
    pub unit_flow: Vec<Vec<f64>>,
    pub river_opening: Vec<Vec<f64>>,
    pub tunnel_opening: Vec<Vec<f64>>,
}

impl Schedule {
    /// Copy physical controls and head guesses for independent reconstruction.
    /// The guesses preserve an optimized pressure allocation between machines
    /// when discharges alone do not determine the intermediate pressure.
    pub fn controls(&self) -> Controls {
        Controls {
            head_guess: self.head.clone(),
            commitment: self.commitment.clone(),
            unit_flow: self.unit_flow.clone(),
            river_opening: self.river_opening.clone(),
            tunnel_opening: self.tunnel_opening.clone(),
        }
    }
}

impl Controls {
    fn validate(&self, p: &SchedulingProblem<'_>) -> Result<(), Error> {
        let s = &p.case.system;
        let nt = p.case.prices.len();
        if self.commitment.len() != s.units.len()
            || self.commitment.iter().any(|r| r.len() != nt)
            || [
                (&self.unit_flow, s.units.len()),
                (&self.river_opening, s.rivers.len()),
                (&self.tunnel_opening, s.tunnels.len()),
            ]
            .iter()
            .any(|(a, n)| {
                a.len() != *n
                    || a.iter()
                        .any(|r| r.len() != nt || r.iter().any(|x| !x.is_finite()))
            })
        {
            return Err(Error::Input(
                "control arrays must have the case's object and interval axes and finite values"
                    .into(),
            ));
        }
        let nodes = s.reservoirs.len() + s.junctions.len() + s.boundaries.len();
        if !self.head_guess.is_empty()
            && (self.head_guess.len() != nodes
                || self
                    .head_guess
                    .iter()
                    .any(|row| row.len() != nt || row.iter().any(|h| !h.is_finite())))
        {
            return Err(Error::Input(
                "head guesses must be finite node-by-interval arrays, or empty".into(),
            ));
        }
        for (g, row) in self.unit_flow.iter().enumerate() {
            for (t, &q) in row.iter().enumerate() {
                if q < 0. || (!self.commitment[g][t] && q != 0.) {
                    return Err(Error::Input(
                        "unit discharges must be nonnegative and zero when stopped".into(),
                    ));
                }
            }
        }
        if self
            .river_opening
            .iter()
            .flatten()
            .any(|a| !(0. ..=1.).contains(a))
        {
            return Err(Error::Input("river gate openings must be in [0,1]".into()));
        }
        for (e, row) in self.tunnel_opening.iter().enumerate() {
            for (t, a) in row.iter().enumerate() {
                if !p.operating.tunnels[e][t]
                    .positions
                    .iter()
                    .any(|v| v.opening == *a)
                {
                    return Err(Error::Input(
                        "tunnel opening is not an allowed position in this interval".into(),
                    ));
                }
            }
        }
        Ok(())
    }
}

pub(crate) struct HydraulicState {
    pub storage: Vec<f64>,
    pub head: Vec<f64>,
    pub tunnel_flow: Vec<f64>,
    pub river_withdrawal: Vec<f64>,
    pub unit_flow: Vec<f64>,
    pub iterations: usize,
}

struct Workspace<D> {
    head: Vec<D>,
    flow: Vec<D>,
    columns: Vec<D>,
}
impl<D: DualNum<Primitive = f64> + Copy> Workspace<D> {
    fn new(p: &SchedulingProblem<'_>) -> Self {
        Self {
            head: vec![D::from(0.); p.topology.incidence.len() + p.case.system.boundaries.len()],
            flow: vec![D::from(0.); p.case.system.units.len()],
            columns: Vec::with_capacity(p.routing.primary_columns + p.routing.transfers.len()),
        }
    }
}

/// Data and coordinate scaling for one hydraulic step. A zero duration solves
/// the stationary pressure problem at prescribed storage, without a water step.
struct HydraulicStep<'a, 'b> {
    p: &'a SchedulingProblem<'b>,
    controls: &'a Controls,
    interval: usize,
    previous: &'a [f64],
    duration: f64,
    columns: &'a [f64],
    pumps: Vec<usize>,
    resistance: Vec<Option<f64>>,
    scales: Vec<f64>,
    rows: Vec<f64>,
    nonnegative: Vec<usize>,
}

fn linear<D: DualNum<Primitive = f64> + Copy>(f: &crate::river::LinearForm, x: &[D]) -> D {
    f.coefficients
        .iter()
        .fold(D::from(f.constant), |v, (&i, &a)| v + x[i] * a)
}

impl HydraulicStep<'_, '_> {
    fn storage_count(&self) -> usize {
        if self.duration == 0. {
            0
        } else {
            self.previous.len()
        }
    }
    fn residual<D: DualNum<Primitive = f64> + Copy>(
        &self,
        x: &[D],
        w: &mut Workspace<D>,
        f: &mut [D],
    ) {
        let p = self.p;
        let (s, topology, routing) = (&p.case.system, &p.topology, &p.routing);
        let (nr, nj, nt) = (s.reservoirs.len(), s.junctions.len(), p.case.prices.len());
        let nv = self.storage_count();
        let tunnel = nv + nj;
        let outlet = tunnel + s.tunnels.len();
        let pump = outlet + routing.outlets.len();
        let t = self.interval;
        w.columns.clear();
        w.columns.extend(self.columns.iter().map(|&v| D::from(v)));
        for (g, q) in w.flow.iter_mut().enumerate() {
            *q = D::from(self.controls.unit_flow[g][t]);
        }
        for (k, &g) in self.pumps.iter().enumerate() {
            w.flow[g] = x[pump + k] * self.scales[pump + k];
            w.columns[g * nt + t] = w.flow[g];
        }
        for e in 0..s.tunnels.len() {
            w.columns[(s.units.len() + e) * nt + t] = x[tunnel + e] * self.scales[tunnel + e];
        }
        for k in 0..routing.outlets.len() {
            w.columns[(s.units.len() + s.tunnels.len() + k) * nt + t] =
                x[outlet + k] * self.scales[outlet + k];
        }
        // First-order lifting of the scalar transfer supplies its exact chain
        // rule to Newton's Jacobian, including shared reference fluxes.
        for transfer in &routing.transfers {
            let q = linear(&transfer.flow, &w.columns);
            let d = transfer.fraction.data(q.re());
            w.columns.push((q - q.re()) * d.gradient[0] + d.value);
        }
        for (r, reservoir) in s.reservoirs.iter().enumerate() {
            if nv == 0 {
                // The coincident-endpoint mean has the level's full derivative.
                w.head[r] = D::from(reservoir.level.trial_level(self.previous[r]));
            } else {
                let after = x[r] * self.scales[r] + reservoir.storage.min;
                let path = &routing.storage_shape[r][t];
                let mut mean = D::from(0.);
                for vertices in path.windows(2) {
                    let volume = |v: &(f64, crate::river::LinearForm)| {
                        after * v.0 + self.previous[r] * (1. - v.0) + linear(&v.1, &w.columns)
                    };
                    let a = volume(&vertices[0]);
                    let b = volume(&vertices[1]);
                    let d = reservoir.level.mean_data(a.re(), b.re());
                    let h = (a - a.re()) * d.gradient[0] + (b - b.re()) * d.gradient[1] + d.value;
                    mean += h * (vertices[1].0 - vertices[0].0);
                }
                w.head[r] = mean;
            }
        }
        for j in 0..nj {
            w.head[nr + j] = x[nv + j] * self.scales[nv + j];
        }
        for (j, b) in s.boundaries.iter().enumerate() {
            w.head[nr + nj + j] = D::from(b.head);
        }
        for node in 0..nr + nj {
            if nv == 0 && node < nr {
                continue;
            }
            let mut net = D::from(if node < nr {
                p.operating.reservoirs[node].inflow[t]
            } else {
                0.
            });
            for &(g, sign) in &topology.incidence[node] {
                net += w.flow[g] * sign;
            }
            for &(e, sign) in &topology.tunnel_incidence[node] {
                net += x[tunnel + e] * (sign * self.scales[tunnel + e]);
            }
            let row = if node < nr { node } else { nv + node - nr };
            if node < nr {
                net += linear(&routing.reservoir_net[node][t], &w.columns);
                let change = x[node] * self.scales[node] + s.reservoirs[node].storage.min
                    - self.previous[node];
                f[row] = (change / self.duration - net) / self.rows[row];
            } else {
                f[row] = -net / self.rows[row];
            }
        }
        for (e, &(a, b)) in topology.tunnel_nodes.iter().enumerate() {
            let q = x[tunnel + e] * self.scales[tunnel + e];
            f[tunnel + e] = if let Some(r) = self.resistance[e] {
                (w.head[a] - w.head[b] - q * q.abs() * r) / self.rows[tunnel + e]
            } else {
                q / self.rows[tunnel + e]
            };
        }
        for (k, &e) in routing.outlets.iter().enumerate() {
            let river = &s.rivers[e];
            let source = s
                .reservoirs
                .iter()
                .position(|r| r.name == river.source)
                .unwrap();
            let h = w.head[source];
            let g = self.controls.river_opening[e][t];
            let rate = match &river.law {
                crate::RiverLaw::Controlled => {
                    D::from(g * (p.operating.river_capacity[e][t] - p.operating.river_inflow[e][t]))
                }
                crate::RiverLaw::Level { curve } => crate::machine::table_value(curve, h) * g,
                crate::RiverLaw::Delta { curve, reference } => {
                    let target = s
                        .reservoirs
                        .iter()
                        .position(|r| r.name == river.target)
                        .unwrap();
                    let href = match reference {
                        crate::ReferenceLevel::Upstream => h,
                        crate::ReferenceLevel::Downstream => w.head[target],
                    };
                    curve.extended_value(h - w.head[target], href)
                }
                crate::RiverLaw::Orifice {
                    coefficient, crest, ..
                } => {
                    if h.re() <= *crest {
                        D::from(0.)
                    } else {
                        (h - *crest).sqrt() * (coefficient * g)
                    }
                }
                crate::RiverLaw::Weir {
                    coefficient, crest, ..
                } => {
                    if h.re() <= *crest {
                        D::from(0.)
                    } else {
                        let d = h - *crest;
                        d * d.sqrt() * *coefficient
                    }
                }
            };
            f[outlet + k] =
                (x[outlet + k] * self.scales[outlet + k] - rate) / self.rows[outlet + k];
        }
        for (k, &g) in self.pumps.iter().enumerate() {
            let unit = &s.units[g];
            let h = machine_head(p, g, &w.head, &w.flow);
            let d = unit
                .pump_characteristic
                .as_ref()
                .unwrap()
                .discharge_data(h.re());
            let expected = (h - h.re()) * d.gradient[1] + d.value;
            f[pump + k] = (w.flow[g] - expected) / self.rows[pump + k];
        }
    }

    fn converged(&self, x: &[f64], f: &[f64], tolerance: f64) -> bool {
        f.iter().enumerate().all(|(i, r)| {
            let rounding = if i < self.storage_count() {
                let v = x[i] * self.scales[i] + self.p.case.system.reservoirs[i].storage.min;
                let v = v.abs().max(self.previous[i].abs());
                2. * (v.next_up() - v) / self.duration
            } else {
                0.
            };
            r.is_finite() && (r * self.rows[i]).abs() <= tolerance + rounding
        })
    }
}

pub(crate) fn machine_head<D: DualNum<Primitive = f64> + Copy>(
    p: &SchedulingProblem<'_>,
    g: usize,
    head: &[D],
    flow: &[D],
) -> D {
    let topology = &p.topology;
    let s = &p.case.system;
    let unit = &s.units[g];
    let plant = topology.unit_plants[g];
    let spec = &s.plants[plant];
    let (a, b) = topology.unit_nodes[g];
    let target = if spec
        .outlet_head_floor
        .is_some_and(|floor| head[b].re() < floor)
    {
        D::from(spec.outlet_head_floor.unwrap())
    } else {
        head[b]
    };
    let total = topology.plant_units[plant]
        .iter()
        .fold(D::from(0.), |v, &g| v + flow[g]);
    let tail = spec
        .tailwater_curve
        .as_ref()
        .map_or(D::from(0.), |c| crate::machine::table_value(c, total));
    (head[a] - target - tail - flow[g] * flow[g] * unit.loss_factor) * unit.mode.power_sign()
}

fn newton(step: &HydraulicStep<'_, '_>, x: &mut [f64]) -> Result<usize, Error> {
    let n = x.len();
    if n == 0 {
        return Ok(0);
    }
    let mut w = Workspace::new(step.p);
    let mut dw = Workspace::<DualSVec64<8>>::new(step.p);
    let mut dual = vec![DualSVec64::<8>::from_re(0.); n];
    let mut df = dual.clone();
    let mut f = vec![0.; n];
    let mut trial_f = f.clone();
    let mut trial = f.clone();
    let mut jacobian = DMatrix::zeros(n, n);
    for iteration in 0..=100 {
        step.residual(x, &mut w, &mut f);
        if step.converged(x, &f, 1e-8) {
            return Ok(iteration);
        }
        for start in (0..n).step_by(8) {
            for (i, v) in dual.iter_mut().enumerate() {
                *v = DualSVec64::from_re(x[i]);
                if i >= start && i < start + 8 {
                    *v = v.derivative(i - start);
                }
            }
            step.residual(&dual, &mut dw, &mut df);
            for (i, v) in df.iter().enumerate() {
                let d = v.eps.unwrap_generic(Const::<8>, U1);
                for j in start..(start + 8).min(n) {
                    jacobian[(i, j)] = d[j - start];
                }
            }
        }
        let rhs = DVector::from_iterator(n, f.iter().map(|r| -r));
        let direction = jacobian
            .clone()
            .full_piv_lu()
            .solve(&rhs)
            .or_else(|| {
                let svd = jacobian.clone().svd(true, true);
                let threshold = f64::EPSILON * n as f64 * svd.singular_values.amax();
                svd.solve(&rhs, threshold).ok()
            })
            .ok_or_else(|| Error::Solver("singular hydraulic Jacobian".into()))?;
        let merit: f64 = f.iter().map(|r| r * r).sum();
        let mut alpha = 1.;
        let mut accepted = false;
        for _ in 0..45 {
            for i in 0..n {
                trial[i] = x[i] + alpha * direction[i];
            }
            for &i in &step.nonnegative {
                trial[i] = trial[i].max(0.);
            }
            step.residual(&trial, &mut w, &mut trial_f);
            let next: f64 = trial_f.iter().map(|r| r * r).sum();
            if next.is_finite() && (next < merit || step.converged(&trial, &trial_f, 1e-8)) {
                x.copy_from_slice(&trial);
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted || iteration == 100 {
            break;
        }
    }
    step.residual(x, &mut w, &mut f);
    let residual = f
        .iter()
        .zip(&step.rows)
        .map(|(f, scale)| (f * scale).abs())
        .fold(0., f64::max);
    Err(Error::Solver(format!(
        "hydraulic reconstruction did not converge in interval {} (maximum equation residual {residual})",
        step.interval
    )))
}

/// A passive pressure component without a fixed-head node admits one constant
/// head shift. Choose that reference within the intersection of absolute head
/// bounds. No head difference, water equation or machine work is changed.
fn pressure_references(
    p: &SchedulingProblem<'_>,
    controls: &Controls,
    interval: usize,
    head: &mut [f64],
    reference: Option<&HydraulicState>,
) {
    let s = &p.case.system;
    let nr = s.reservoirs.len();
    let nj = s.junctions.len();
    let mut component: Vec<_> = (0..head.len()).collect();
    for (e, &(a, b)) in p.topology.tunnel_nodes.iter().enumerate() {
        if controls.tunnel_opening[e][interval] == 0. {
            continue;
        }
        let (a, b) = (component[a], component[b]);
        for root in &mut component {
            if *root == b {
                *root = a;
            }
        }
    }
    let mut anchored = vec![false; head.len()];
    for node in (0..nr).chain(nr + nj..head.len()) {
        anchored[component[node]] = true;
    }
    for (g, &(a, b)) in p.topology.unit_nodes.iter().enumerate() {
        if !controls.commitment[g][interval] {
            continue;
        }
        let plant = &s.plants[p.topology.unit_plants[g]];
        let (a, b) = (component[a], component[b]);
        if a != b || plant.outlet_head_floor.is_some() {
            anchored[a] = true;
            anchored[b] = true;
        }
    }
    for (root, &is_anchored) in anchored.iter().enumerate() {
        if is_anchored {
            continue;
        }
        let nodes: Vec<_> = (nr..nr + nj).filter(|&j| component[j] == root).collect();
        let Some(&first) = nodes.first() else {
            continue;
        };
        let lo = nodes
            .iter()
            .map(|&j| s.junctions[j - nr].head.min - head[j])
            .fold(f64::NEG_INFINITY, f64::max);
        let hi = nodes
            .iter()
            .map(|&j| s.junctions[j - nr].head.max - head[j])
            .fold(f64::INFINITY, f64::min);
        if lo > hi {
            continue;
        }
        let desired = controls
            .head_guess
            .get(first)
            .map(|row| row[interval])
            .or_else(|| reference.map(|v| v.head[first]));
        let shift = desired.map_or(lo / 2. + hi / 2., |h| (h - head[first]).clamp(lo, hi));
        for j in nodes {
            head[j] += shift;
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn hydraulic_step(
    p: &SchedulingProblem<'_>,
    controls: &Controls,
    interval: usize,
    previous: &[f64],
    duration: f64,
    columns: &[f64],
    reference: Option<&HydraulicState>,
) -> Result<HydraulicState, Error> {
    let s = &p.case.system;
    let topology = &p.topology;
    let nr = if duration == 0. {
        0
    } else {
        s.reservoirs.len()
    };
    let nj = s.junctions.len();
    let ne = s.tunnels.len();
    let pumps: Vec<_> = s
        .units
        .iter()
        .enumerate()
        .filter_map(|(g, u)| {
            (u.pump_characteristic.is_some() && controls.commitment[g][interval]).then_some(g)
        })
        .collect();
    let known_head: Vec<_> = s
        .reservoirs
        .iter()
        .enumerate()
        .map(|(r, v)| v.level.trial_level(previous[r]))
        .chain(s.boundaries.iter().map(|b| b.head))
        .collect();
    let mean = known_head.iter().sum::<f64>() / known_head.len().max(1) as f64;
    let head: Vec<_> = known_head[..s.reservoirs.len()]
        .iter()
        .copied()
        .chain(s.junctions.iter().enumerate().map(|(j, junction)| {
            // Bounds initialize the numerical guess; Newton still solves
            // the original equations and can reveal a bound violation.
            controls
                .head_guess
                .get(s.reservoirs.len() + j)
                .map(|row| row[interval])
                .or_else(|| reference.map(|v| v.head[s.reservoirs.len() + j]))
                .unwrap_or(mean)
                .clamp(junction.head.min, junction.head.max)
        }))
        .chain(s.boundaries.iter().map(|b| b.head))
        .collect();
    let resistance: Vec<_> = p
        .operating
        .tunnels
        .iter()
        .enumerate()
        .map(|(e, row)| {
            let a = controls.tunnel_opening[e][interval];
            (a > 0.).then(|| {
                row[interval]
                    .positions
                    .iter()
                    .find(|v| v.opening == a)
                    .unwrap()
                    .resistance
            })
        })
        .collect();
    let mut x = Vec::new();
    let mut scales = Vec::new();
    let mut rows = Vec::new();
    for (r, v) in s.reservoirs.iter().enumerate().take(nr) {
        let scale = v.storage.max - v.storage.min;
        scales.push(scale);
        x.push((previous[r] - v.storage.min) / scale);
    }
    for j in 0..nj {
        scales.push(100.);
        x.push(head[s.reservoirs.len() + j] / 100.);
    }
    for node in (0..nr).chain(s.reservoirs.len()..s.reservoirs.len() + nj) {
        let scale = topology.incidence[node]
            .iter()
            .map(|&(g, _)| s.units[g].flow.max)
            .chain(
                topology.tunnel_incidence[node]
                    .iter()
                    .map(|&(e, _)| s.tunnels[e].capacity),
            )
            .fold(1., f64::max);
        rows.push(scale);
    }
    let mut nonnegative = Vec::new();
    for (e, &(a, b)) in topology.tunnel_nodes.iter().enumerate() {
        let scale = s.tunnels[e].capacity.max(1.);
        let q = resistance[e].map_or(0., |r| {
            let dh = head[a] - head[b];
            let physical = dh.signum() * (dh.abs() / r).sqrt();
            reference.map_or(physical, |v| {
                let q = v.tunnel_flow[e];
                if q.signum() == dh.signum() && q != 0. {
                    q
                } else {
                    physical
                }
            })
        });
        x.push(q / scale);
        scales.push(scale);
        rows.push(if resistance[e].is_some() { 100. } else { scale });
        if topology.tunnel_outfalls[e].is_some() {
            nonnegative.push(nr + nj + e);
        }
    }
    for (k, &e) in p.routing.outlets.iter().enumerate() {
        let river = &s.rivers[e];
        let source = s
            .reservoirs
            .iter()
            .position(|r| r.name == river.source)
            .unwrap();
        let scale = river.capacity.max(1.);
        let q = if let crate::RiverLaw::Delta { curve, reference } = &river.law {
            let target = s
                .reservoirs
                .iter()
                .position(|r| r.name == river.target)
                .unwrap();
            let href = match reference {
                crate::ReferenceLevel::Upstream => head[source],
                crate::ReferenceLevel::Downstream => head[target],
            };
            curve.extended_value(head[source] - head[target], href)
        } else {
            river.law.trial_outlet_flow(
                head[source],
                controls.river_opening[e][interval],
                p.operating.river_capacity[e][interval],
                p.operating.river_inflow[e][interval],
            )
        };
        x.push(q / scale);
        scales.push(scale);
        rows.push(scale);
        if !river.law.is_delta() {
            nonnegative.push(nr + nj + ne + k);
        }
    }
    for (k, &g) in pumps.iter().enumerate() {
        let scale = s.units[g].flow.max.max(1.);
        let q = reference.map_or(controls.unit_flow[g][interval], |v| v.unit_flow[g]);
        x.push(q.max(0.) / scale);
        scales.push(scale);
        rows.push(scale);
        nonnegative.push(nr + nj + ne + p.routing.outlets.len() + k);
    }
    let step = HydraulicStep {
        p,
        controls,
        interval,
        previous,
        duration,
        columns,
        pumps,
        resistance,
        scales,
        rows,
        nonnegative,
    };
    let iterations = newton(&step, &mut x)?;
    let mut w = Workspace::new(p);
    step.residual(&x, &mut w, &mut vec![0.; x.len()]);
    pressure_references(p, controls, interval, &mut w.head, reference);
    Ok(HydraulicState {
        storage: if nr == 0 {
            previous.to_vec()
        } else {
            (0..nr)
                .map(|r| s.reservoirs[r].storage.min + x[r] * step.scales[r])
                .collect()
        },
        head: w.head,
        tunnel_flow: (0..ne)
            .map(|e| x[nr + nj + e] * step.scales[nr + nj + e])
            .collect(),
        river_withdrawal: (0..p.routing.outlets.len())
            .map(|k| x[nr + nj + ne + k] * step.scales[nr + nj + ne + k])
            .collect(),
        unit_flow: w.flow,
        iterations,
    })
}

impl SchedulingProblem<'_> {
    /// Reconstruct physical states from controls without an optimization solver.
    /// The returned validation reports operating violations; reconstruction does
    /// not clip storage, heads, capacities or controls to make a schedule feasible.
    pub fn reconstruct(&self, controls: &Controls) -> Result<Schedule, Error> {
        controls.validate(self)?;
        let s = &self.case.system;
        let nt = self.case.prices.len();
        let mut schedule = Schedule {
            storage: s
                .reservoirs
                .iter()
                .map(|r| vec![r.initial_storage; nt + 1])
                .collect(),
            unit_flow: controls.unit_flow.clone(),
            tunnel_opening: controls.tunnel_opening.clone(),
            river_opening: controls.river_opening.clone(),
            commitment: controls.commitment.clone(),
            tunnel_flow: vec![vec![0.; nt]; s.tunnels.len()],
            river_withdrawal: vec![vec![0.; nt]; s.rivers.len()],
            gate_adjustment_cost: 0.,
            observed_flow: vec![],
            observed_volume: vec![],
            river_release: vec![],
            river_arrival_volume: vec![],
            river_storage: vec![],
            unit_power: vec![],
            unit_head: vec![],
            plant_flow: vec![],
            plant_power: vec![],
            head: vec![],
            soft_slacks: vec![],
            soft_penalty_cost: 0.,
            objective: 0.,
            solver_status: None,
            iterations: 0,
            validation: Default::default(),
        };
        let mut columns: Vec<_> = controls.unit_flow.iter().flatten().copied().collect();
        columns.resize(self.routing.primary_columns, 0.);
        let mut previous: Vec<_> = s.reservoirs.iter().map(|r| r.initial_storage).collect();
        let mut reference = None;
        let mut junctions = vec![vec![0.; nt]; s.junctions.len()];
        for t in 0..nt {
            let state = hydraulic_step(
                self,
                controls,
                t,
                &previous,
                3600. * (self.case.grid[t + 1] - self.case.grid[t]),
                &columns,
                reference.as_ref(),
            )?;
            for (r, row) in schedule.storage.iter_mut().enumerate() {
                row[t + 1] = state.storage[r];
            }
            for (g, row) in schedule.unit_flow.iter_mut().enumerate() {
                row[t] = state.unit_flow[g];
                columns[g * nt + t] = row[t];
            }
            for (e, row) in schedule.tunnel_flow.iter_mut().enumerate() {
                row[t] = state.tunnel_flow[e];
                columns[(s.units.len() + e) * nt + t] = row[t];
            }
            for (k, &e) in self.routing.outlets.iter().enumerate() {
                schedule.river_withdrawal[e][t] = state.river_withdrawal[k];
                columns[(s.units.len() + s.tunnels.len() + k) * nt + t] = state.river_withdrawal[k];
            }
            for (j, row) in junctions.iter_mut().enumerate() {
                row[t] = state.head[s.reservoirs.len() + j];
            }
            schedule.iterations += state.iterations;
            previous.clone_from(&state.storage);
            reference = Some(state);
        }
        let mut physical = crate::trajectory::evaluate(self, schedule, junctions, 1e-7);
        self.routing.complete_columns(&mut columns);
        let points = self
            .operating_points
            .with_level_crossings(self, &physical.storage, &columns);
        let excess = points.residual(self, &physical)?;
        physical.validation.operating_limits = physical.validation.operating_limits.max(excess);
        physical.validation.feasible &= excess <= 1e-7;
        Ok(physical)
    }
}
