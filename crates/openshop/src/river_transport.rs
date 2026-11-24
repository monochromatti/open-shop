//! Conservative transport of interval releases through delay distributions.
//!
//! A delay bin is uniform in travel time, with the supplied probability mass.
//! Flow selects a linear blend of adjacent response distributions at release.
use crate::{Error, curve::Derivatives2};
use serde::{Deserialize, Serialize};

/// One-sided value at a change in arrival rate. Distributed arrivals are
/// continuous; a fixed travel time preserves jumps in the release series.
#[derive(Clone, Copy, Debug)]
pub enum TimeSide {
    Left,
    Right,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "DelayData", into = "DelayData")]
pub struct DelayCurve {
    data: DelayData,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DelayData {
    reference_flow: f64,
    edges: Vec<f64>,
    weights: Vec<f64>,
}
impl From<DelayCurve> for DelayData {
    fn from(c: DelayCurve) -> Self {
        c.data
    }
}
impl TryFrom<DelayData> for DelayCurve {
    type Error = Error;
    fn try_from(d: DelayData) -> Result<Self, Error> {
        Self::new(d.reference_flow, d.edges, d.weights)
    }
}
impl DelayCurve {
    pub fn new(reference_flow: f64, edges: Vec<f64>, mut weights: Vec<f64>) -> Result<Self, Error> {
        if weights.len() == edges.len() && weights.last() == Some(&0.) {
            weights.pop();
        }
        if !reference_flow.is_finite()
            || reference_flow < 0.
            || edges.len() < 2
            || !edges.iter().all(|e| e.is_finite() && *e >= 0.)
            || !edges
                .windows(2)
                .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())
            || weights.len() + 1 != edges.len()
            || !weights.iter().all(|w| w.is_finite() && *w >= 0.)
            || (weights.iter().sum::<f64>() - 1.).abs() > 2e-12
        {
            return Err(Error::Input("delay distributions need ordered nonnegative finite edges and nonnegative bin probabilities summing to one".into()));
        }
        Ok(Self {
            data: DelayData {
                reference_flow,
                edges,
                weights,
            },
        })
    }
    pub fn reference_flow(&self) -> f64 {
        self.data.reference_flow
    }
    pub fn edges(&self) -> &[f64] {
        &self.data.edges
    }
    pub fn weights(&self) -> &[f64] {
        &self.data.weights
    }
    // Sum of two uniform intervals has a trapezoidal density. Integrate its
    // nonnegative pieces directly; differences of nearly equal CDFs lose tails.
    fn arrival_fraction(&self, a: f64, b: f64, s: f64, e: f64) -> f64 {
        self.data
            .edges
            .windows(2)
            .zip(&self.data.weights)
            .map(|(edges, &weight)| {
                let width = edges[1] - edges[0];
                let short = (e - s).min(width);
                let long = (e - s).max(width);
                let ramp = |lo: f64, hi: f64| {
                    let l = lo.max(0.);
                    let h = hi.min(short);
                    if h > l {
                        ((h - l) / long) * (l / short + h / short) / 2.
                    } else {
                        0.
                    }
                };
                let lo = a - s - edges[0];
                let hi = b - s - edges[0];
                weight
                    * (ramp(lo, hi)
                        + (hi.min(long) - lo.max(short)).max(0.) / long
                        + ramp(e - b + edges[1], e - a + edges[1]))
            })
            .sum::<f64>()
            .clamp(0., 1.)
    }
    // Integrate survival directly; a small remaining inventory must not be
    // recovered by subtracting two much larger water volumes.
    fn transit_fraction(&self, time: f64, s: f64, e: f64) -> f64 {
        let released = time.min(e);
        if released <= s {
            return 0.;
        }
        let lo = time - released;
        let hi = time - s;
        self.data
            .edges
            .windows(2)
            .zip(&self.data.weights)
            .map(|(edges, &weight)| {
                let prefix = (hi.min(edges[0]) - lo).max(0.);
                let l = lo.max(edges[0]);
                let h = hi.min(edges[1]);
                let ramp = if h > l {
                    (h - l) * ((edges[1] - l) + (edges[1] - h)) / (2. * (edges[1] - edges[0]))
                } else {
                    0.
                };
                weight * (prefix + ramp) / (e - s)
            })
            .sum()
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RiverTransport {
    Fixed { delay: f64 },
    Distributed { curves: Vec<DelayCurve> },
}
impl Default for RiverTransport {
    fn default() -> Self {
        Self::Fixed { delay: 0. }
    }
}
impl RiverTransport {
    pub(crate) fn fixed_delay(&self) -> Option<f64> {
        match self {
            Self::Fixed { delay } => Some(*delay),
            _ => None,
        }
    }
    pub(crate) fn validate(&self, capacity: f64, end: f64) -> Result<(), Error> {
        let valid = match self {
            Self::Fixed { delay } => delay.is_finite() && *delay >= 0. && (end + delay).is_finite(),
            Self::Distributed { curves } => {
                !curves.is_empty()
                    && curves.windows(2).all(|c| {
                        c[1].reference_flow() > c[0].reference_flow()
                            && (c[1].reference_flow() - c[0].reference_flow()).is_finite()
                    })
                    && (curves.len() == 1
                        || (curves[0].reference_flow() == 0.
                            && curves.last().unwrap().reference_flow() >= capacity))
                    && curves
                        .iter()
                        .all(|c| (end + c.edges().last().unwrap()).is_finite())
            }
        };
        if valid {
            Ok(())
        } else {
            Err(Error::Input("river delay references must cover zero through capacity; fixed travel times must be finite and nonnegative".into()))
        }
    }
    pub(crate) fn fraction(&self, window: Observation, s: f64, e: f64) -> ReleaseFraction {
        let at = |c: &DelayCurve| match window {
            Observation::Arrival(a, b) => c.arrival_fraction(a, b, s, e),
            Observation::Transit(time) => c.transit_fraction(time, s, e),
            Observation::Point(time, _) => c
                .edges()
                .windows(2)
                .zip(c.weights())
                .map(|(w, &weight)| {
                    weight * (w[1].min(time - s) - w[0].max(time - e)).max(0.) / (w[1] - w[0])
                })
                .sum(),
        };
        match self {
            Self::Fixed { delay } => ReleaseFraction {
                references: vec![0.],
                coefficients: vec![match window {
                    Observation::Arrival(a, b) => {
                        (b.min(e + delay) - a.max(s + delay)).max(0.) / (e - s)
                    }
                    Observation::Transit(time) => {
                        (time.min(e) - s.max(time - delay)).max(0.) / (e - s)
                    }
                    Observation::Point(time, side) => {
                        let inside = match side {
                            TimeSide::Left => s + delay < time && time <= e + delay,
                            TimeSide::Right => s + delay <= time && time < e + delay,
                        };
                        if inside { 1. } else { 0. }
                    }
                }],
            },
            Self::Distributed { curves } => ReleaseFraction {
                references: curves.iter().map(DelayCurve::reference_flow).collect(),
                coefficients: curves.iter().map(at).collect(),
            },
        }
    }
    /// Instantaneous downstream rate, m³/s, from one constant release.
    /// At a fixed-delay jump, `side` selects the value before or after it.
    pub fn arrival_rate(
        &self,
        release: (f64, f64),
        flow: f64,
        time: f64,
        side: TimeSide,
    ) -> Result<f64, Error> {
        self.check_observation(release, flow)?;
        if !time.is_finite() {
            return Err(Error::Input(
                "finite arrival observation time required".into(),
            ));
        }
        Ok(self
            .fraction(Observation::Point(time, side), release.0, release.1)
            .data(flow)
            .value)
    }
    /// Volume (m³) from one constant release arriving within an observation window.
    pub fn arrival_volume(
        &self,
        release: (f64, f64),
        flow: f64,
        window: (f64, f64),
    ) -> Result<f64, Error> {
        self.check_observation(release, flow)?;
        if !window.0.is_finite() || !window.1.is_finite() || window.1 < window.0 {
            return Err(Error::Input(
                "ordered finite arrival window required".into(),
            ));
        }
        Ok(self
            .fraction(
                Observation::Arrival(window.0, window.1),
                release.0,
                release.1,
            )
            .data(flow)
            .value
            * (release.1 - release.0)
            * 3600.)
    }
    /// Already released water (m³) still in transit at `time`.
    pub fn remaining_volume(
        &self,
        release: (f64, f64),
        flow: f64,
        time: f64,
    ) -> Result<f64, Error> {
        self.check_observation(release, flow)?;
        if !time.is_finite() {
            return Err(Error::Input("finite inventory time required".into()));
        }
        Ok(self
            .fraction(Observation::Transit(time), release.0, release.1)
            .data(flow)
            .value
            * (release.1 - release.0)
            * 3600.)
    }
    fn check_observation(&self, release: (f64, f64), flow: f64) -> Result<(), Error> {
        self.validate(flow.max(0.), release.1)?;
        if !release.0.is_finite()
            || !release.1.is_finite()
            || release.1 <= release.0
            || !(release.1 - release.0).is_finite()
            || !flow.is_finite()
            || flow < 0.
        {
            return Err(Error::Input(
                "finite nonnegative release and increasing finite release interval required".into(),
            ));
        }
        Ok(())
    }
}
#[derive(Clone, Copy)]
pub(crate) enum Observation {
    Arrival(f64, f64),
    Transit(f64),
    Point(f64, TimeSide),
}
#[derive(Clone)]
pub(crate) struct ReleaseFraction {
    pub references: Vec<f64>,
    pub coefficients: Vec<f64>,
}
impl ReleaseFraction {
    pub(crate) fn linear_coefficient(&self) -> Option<f64> {
        self.coefficients
            .iter()
            .all(|c| *c == self.coefficients[0])
            .then_some(self.coefficients[0])
    }
    pub(crate) fn data(&self, q: f64) -> Derivatives2 {
        if let Some(c) = self.linear_coefficient() {
            return Derivatives2 {
                value: c * q,
                gradient: [c, 0.],
                hessian: [[0.; 2]; 2],
            };
        }
        let j = self
            .references
            .partition_point(|v| *v <= q)
            .saturating_sub(1)
            .min(self.references.len() - 2);
        let slope = (self.coefficients[j + 1] - self.coefficients[j])
            / (self.references[j + 1] - self.references[j]);
        let fraction = self.coefficients[j] + slope * (q - self.references[j]);
        Derivatives2 {
            value: q * fraction,
            gradient: [fraction + q * slope, 0.],
            hessian: [[2. * slope, 0.], [0.; 2]],
        }
    }

    /// Range of the transported flux minus slope × release over the reference
    /// domain. Each interval is quadratic, with extrema at its endpoints
    /// or the stationary point.
    pub(crate) fn affine_residual_bounds(&self, slope: f64) -> crate::Bounds {
        let mut lower = f64::INFINITY;
        let mut upper = f64::NEG_INFINITY;
        for (q, weights) in self.references.windows(2).zip(self.coefficients.windows(2)) {
            let a = (weights[1] - weights[0]) / (q[1] - q[0]);
            let stationary = if a == 0. {
                q[0]
            } else {
                (0.5 * (q[0] + (slope - weights[0]) / a)).clamp(q[0], q[1])
            };
            for release in [q[0], q[1], stationary] {
                let residual = release * (weights[0] + a * (release - q[0]) - slope);
                lower = lower.min(residual);
                upper = upper.max(residual);
            }
        }
        crate::Bounds::new(lower, upper)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn released_transfer_derivatives_include_reference_blend_curvature() {
        let f = ReleaseFraction {
            references: vec![0., 10., 40.],
            coefficients: vec![0.8, 0.6, 0.2],
        };
        for q in [-2., 3., 17., 45.] {
            let d = f.data(q);
            let h = 1e-4;
            let left = f.data(q - h);
            let right = f.data(q + h);
            assert!(((right.value - left.value) / (2. * h) - d.gradient[0]).abs() < 1e-8);
            assert!(
                ((right.gradient[0] - left.gradient[0]) / (2. * h) - d.hessian[0][0]).abs() < 1e-8
            );
        }
    }
}
