//! Measured turbine efficiency and head-dependent operating envelopes.
//! Original discharge grids are retained; interpolation is linear in head.

use crate::{Bounds, Error, TableCurve, curve::Derivatives2};
use num_dual::DualNum;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DischargeInterpolation {
    #[default]
    Linear,
    /// Shape-preserving piecewise cubic Hermite interpolation (PCHIP).
    Pchip,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Extrapolation {
    #[default]
    Error,
    Linear,
    /// Only for a single reference head; discharge and efficiency ignore head.
    Constant,
}
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TurbineOptions {
    pub interpolation: DischargeInterpolation,
    pub head_extrapolation: Extrapolation,
    pub discharge_extrapolation: Extrapolation,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "TurbineData", into = "TurbineData")]
pub struct TurbineTable {
    data: TurbineData,
    // Coefficients in the local normalized coordinate (q - x_i)/(x_{i+1}-x_i).
    coefficients: Vec<Vec<[f64; 4]>>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TurbineData {
    heads: Vec<f64>,
    curves: Vec<TableCurve>,
    flow: Vec<Bounds>,
    #[serde(default)]
    options: TurbineOptions,
}
impl From<TurbineTable> for TurbineData {
    fn from(table: TurbineTable) -> Self {
        table.data
    }
}
impl TryFrom<TurbineData> for TurbineTable {
    type Error = Error;
    fn try_from(data: TurbineData) -> Result<Self, Error> {
        Self::new(data.heads, data.curves, data.flow, data.options)
    }
}
impl TurbineTable {
    /// Heads in metres; one efficiency curve and discharge range (m³/s) per head.
    pub fn new(
        heads: Vec<f64>,
        curves: Vec<TableCurve>,
        flow: Vec<Bounds>,
        options: TurbineOptions,
    ) -> Result<Self, Error> {
        if heads.is_empty()
            || curves.len() != heads.len()
            || flow.len() != heads.len()
            || !heads.iter().all(|h| h.is_finite() && *h > 0.)
            || !heads
                .windows(2)
                .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())
            || !curves
                .iter()
                .all(|c| c.coordinates()[0] >= 0. && c.values().iter().all(|e| *e > 0. && *e <= 1.))
            || !flow
                .iter()
                .all(|b| b.min.is_finite() && b.max.is_finite() && 0. <= b.min && b.min <= b.max)
            || (heads.len() == 1 && options.head_extrapolation != Extrapolation::Constant)
            || (heads.len() > 1 && options.head_extrapolation == Extrapolation::Constant)
            || options.discharge_extrapolation == Extrapolation::Constant
        {
            return Err(Error::Input("turbine table needs ordered positive heads, paired curves/envelopes and valid extension policies".into()));
        }
        if options.discharge_extrapolation == Extrapolation::Error
            && curves.iter().zip(&flow).any(|(c, b)| {
                c.coordinates()[0] > b.min || *c.coordinates().last().unwrap() < b.max
            })
        {
            return Err(Error::Input("turbine efficiency curves must cover discharge envelopes, or declare linear discharge extension".into()));
        }
        let mut coefficients = Vec::with_capacity(curves.len());
        for curve in &curves {
            let x = curve.coordinates();
            let y = curve.values();
            let widths: Vec<_> = x.windows(2).map(|w| w[1] - w[0]).collect();
            let secants: Vec<_> = y
                .windows(2)
                .zip(&widths)
                .map(|(w, h)| (w[1] - w[0]) / h)
                .collect();
            let mut slopes = vec![0.; x.len()];
            if x.len() == 2 {
                slopes.fill(secants[0]);
            } else {
                for i in 1..x.len() - 1 {
                    let (a, b) = (secants[i - 1], secants[i]);
                    if a.signum() == b.signum() && a != 0. && b != 0. {
                        let w1 = 2. * widths[i] + widths[i - 1];
                        let w2 = widths[i] + 2. * widths[i - 1];
                        slopes[i] = (w1 + w2) / (w1 / a + w2 / b);
                    }
                }
                let endpoint = |h1: f64, h2: f64, d1: f64, d2: f64| {
                    let value = ((2. * h1 + h2) * d1 - h1 * d2) / (h1 + h2);
                    if value.signum() != d1.signum() {
                        0.
                    } else if d1.signum() != d2.signum() && value.abs() > 3. * d1.abs() {
                        3. * d1
                    } else {
                        value
                    }
                };
                slopes[0] = endpoint(widths[0], widths[1], secants[0], secants[1]);
                let i = widths.len() - 1;
                slopes[i + 1] = endpoint(widths[i], widths[i - 1], secants[i], secants[i - 1]);
            }
            let cells: Vec<_> = (0..widths.len())
                .map(|i| match options.interpolation {
                    DischargeInterpolation::Linear => [y[i], y[i + 1] - y[i], 0., 0.],
                    DischargeInterpolation::Pchip => {
                        let a = widths[i] * slopes[i];
                        let b = widths[i] * slopes[i + 1];
                        let delta = y[i + 1] - y[i];
                        [y[i], a, 3. * delta - 2. * a - b, -2. * delta + a + b]
                    }
                })
                .collect();
            if !cells.iter().flatten().all(|v| v.is_finite()) {
                return Err(Error::Input(
                    "nonfinite turbine interpolation coefficients".into(),
                ));
            }
            coefficients.push(cells);
        }
        Ok(Self {
            data: TurbineData {
                heads,
                curves,
                flow,
                options,
            },
            coefficients,
        })
    }
    pub fn heads(&self) -> &[f64] {
        &self.data.heads
    }
    pub fn curves(&self) -> &[TableCurve] {
        &self.data.curves
    }
    pub fn flow_envelopes(&self) -> &[Bounds] {
        &self.data.flow
    }
    pub fn options(&self) -> TurbineOptions {
        self.data.options
    }
    fn head_weights(&self, h: f64) -> (usize, usize, f64) {
        if self.data.heads.len() == 1 {
            return (0, 0, 0.);
        }
        let j = self
            .data
            .heads
            .partition_point(|v| *v <= h)
            .saturating_sub(1)
            .min(self.data.heads.len() - 2);
        (
            j,
            j + 1,
            (h - self.data.heads[j]) / (self.data.heads[j + 1] - self.data.heads[j]),
        )
    }
    fn check_head(&self, h: f64) -> Result<(), Error> {
        if !h.is_finite()
            || (self.data.options.head_extrapolation == Extrapolation::Error
                && (h < self.data.heads[0] || h > *self.data.heads.last().unwrap()))
        {
            return Err(Error::Input("head outside turbine table domain".into()));
        }
        Ok(())
    }
    /// Efficiency within the declared data/extension domain, without clipping.
    pub fn efficiency(&self, q: f64, h: f64) -> Result<f64, Error> {
        self.check_head(h)?;
        if !q.is_finite() {
            return Err(Error::Input("nonfinite turbine discharge".into()));
        }
        let (j, k, b) = self.head_weights(h);
        if self.data.options.discharge_extrapolation == Extrapolation::Error {
            for col in [if b == 1. { k } else { j }, if b == 0. { j } else { k }] {
                let x = self.data.curves[col].coordinates();
                if q < x[0] || q > *x.last().unwrap() {
                    return Err(Error::Input(
                        "discharge outside turbine efficiency data".into(),
                    ));
                }
            }
        }
        Ok(self.trial_efficiency(q, h))
    }
    /// Additional on-state discharge range at the actual net head.
    pub fn flow_limits(&self, h: f64) -> Result<Bounds, Error> {
        self.check_head(h)?;
        let (j, k, b) = self.head_weights(h);
        Ok(Bounds::new(
            self.data.flow[j].min + b * (self.data.flow[k].min - self.data.flow[j].min),
            self.data.flow[j].max + b * (self.data.flow[k].max - self.data.flow[j].max),
        ))
    }
    // Return a cell polynomial in its own coordinate, using end secants outside data.
    fn column_polynomial(&self, col: usize, q: f64) -> (f64, f64, [f64; 4]) {
        let c = &self.data.curves[col];
        let x = c.coordinates();
        let y = c.values();
        let i = x
            .partition_point(|v| *v <= q)
            .saturating_sub(1)
            .min(x.len() - 2);
        let coefficients = if q < x[0] || q > *x.last().unwrap() {
            [y[i], y[i + 1] - y[i], 0., 0.]
        } else {
            self.coefficients[col][i]
        };
        (x[i], x[i + 1] - x[i], coefficients)
    }
    fn column<D: DualNum<Primitive = f64> + Copy>(&self, col: usize, q: D) -> D {
        let (origin, width, c) = self.column_polynomial(col, q.re());
        let t = (q - origin) / width;
        ((t * c[3] + c[2]) * t + c[1]) * t + c[0]
    }
    pub(crate) fn trial_efficiency<D: DualNum<Primitive = f64> + Copy>(&self, q: D, h: D) -> D {
        let (j, k, _) = self.head_weights(h.re());
        let left = self.column(j, q);
        if j == k {
            return left;
        }
        let b = (h - self.data.heads[j]) / (self.data.heads[k] - self.data.heads[j]);
        left + b * (self.column(k, q) - left)
    }
    pub(crate) fn envelope_data(&self, h: f64, maximum: bool) -> Derivatives2 {
        let (j, k, b) = self.head_weights(h);
        let value = |i: usize| {
            if maximum {
                self.data.flow[i].max
            } else {
                self.data.flow[i].min
            }
        };
        let slope = if j == k {
            0.
        } else {
            (value(k) - value(j)) / (self.data.heads[k] - self.data.heads[j])
        };
        Derivatives2 {
            value: value(j) + b * (value(k) - value(j)),
            gradient: [0., slope],
            hessian: [[0.; 2]; 2],
        }
    }
    pub(crate) fn validate_domain(&self, flow: Bounds, head: Bounds) -> Result<(), Error> {
        self.check_head(head.min)?;
        self.check_head(head.max)?;
        let mut hh: Vec<_> = self
            .data
            .heads
            .iter()
            .copied()
            .filter(|h| head.min < *h && *h < head.max)
            .collect();
        hh.extend([head.min, head.max]);
        let mut qq: Vec<_> = self
            .data
            .curves
            .iter()
            .flat_map(|c| c.coordinates().iter().copied())
            .filter(|q| flow.min < *q && *q < flow.max)
            .collect();
        qq.extend([flow.min, flow.max]);
        qq.sort_by(f64::total_cmp);
        qq.dedup();
        for h in hh {
            let limits = self.flow_limits(h)?;
            if !limits.min.is_finite()
                || !limits.max.is_finite()
                || limits.min < 0.
                || limits.min > limits.max
                || limits.min > flow.max
                || limits.max < flow.min
            {
                return Err(Error::Input(
                    "turbine and unit discharge envelopes do not overlap throughout operating head"
                        .into(),
                ));
            }
            for q in [flow.min, flow.max] {
                let efficiency = self.efficiency(q, h)?;
                if !efficiency.is_finite() || efficiency <= 0. || efficiency > 1. + 1e-14 {
                    return Err(Error::Input(
                        "extended turbine efficiency outside (0,1]".into(),
                    ));
                }
            }
            // Each discharge cell is cubic and each head cell is affine. Its extrema
            // occur at head edges and discharge endpoints or stationary points.
            for w in qq.windows(2) {
                let (j, k, b) = self.head_weights(h);
                let transformed = |col| {
                    let (origin, width, c) = self.column_polynomial(col, (w[0] + w[1]) / 2.);
                    let a = (w[0] - origin) / width;
                    let r = (w[1] - w[0]) / width;
                    [
                        ((c[3] * a + c[2]) * a + c[1]) * a + c[0],
                        r * ((3. * c[3] * a + 2. * c[2]) * a + c[1]),
                        r * r * (c[2] + 3. * c[3] * a),
                        r * r * r * c[3],
                    ]
                };
                let left = transformed(j);
                let right = transformed(k);
                let c: [f64; 4] = std::array::from_fn(|i| left[i] + b * (right[i] - left[i]));
                let valid = |t: f64| {
                    let v = ((c[3] * t + c[2]) * t + c[1]) * t + c[0];
                    v.is_finite() && v > 0. && v <= 1. + 1e-14
                };
                if !valid(0.) || !valid(1.) {
                    return Err(Error::Input(
                        "extended turbine efficiency outside (0,1]".into(),
                    ));
                }
                let (a, b, d) = (3. * c[3], 2. * c[2], c[1]);
                let roots = if a == 0. {
                    [-d / b, f64::NAN]
                } else {
                    let discriminant = b * b - 4. * a * d;
                    let q = -0.5 * (b + discriminant.sqrt().copysign(b));
                    [q / a, d / q]
                };
                if roots.iter().any(|t| *t > 0. && *t < 1. && !valid(*t)) {
                    return Err(Error::Input(
                        "extended turbine efficiency outside (0,1] between knots".into(),
                    ));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_dual::{second_derivative, second_partial_derivative};
    fn generation_data(t: &TurbineTable, q: f64, h: f64) -> Derivatives2 {
        let (value, dq, dh, cross) =
            second_partial_derivative(|(q, h)| q * h * t.trial_efficiency(q, h) * 0.00981, (q, h));
        let (_, _, qq) =
            second_derivative(|q| q * h * t.trial_efficiency(q, h.into()) * 0.00981, q);
        let (_, _, hh) =
            second_derivative(|h| h * q * t.trial_efficiency(q.into(), h) * 0.00981, h);
        Derivatives2 {
            value,
            gradient: [dq, dh],
            hessian: [[qq, cross], [cross, hh]],
        }
    }
    fn table(interpolation: DischargeInterpolation) -> TurbineTable {
        TurbineTable::new(
            vec![100., 200.],
            vec![
                TableCurve::new(vec![1., 4., 10.], vec![0.75, 0.91, 0.8]).unwrap(),
                TableCurve::new(vec![0., 6., 14., 20.], vec![0.7, 0.9, 0.93, 0.82]).unwrap(),
            ],
            vec![Bounds::new(1., 10.), Bounds::new(2., 18.)],
            TurbineOptions {
                interpolation,
                ..Default::default()
            },
        )
        .unwrap()
    }
    #[test]
    fn unequal_discharge_grids_interpolate_independently_in_head() {
        let t = table(DischargeInterpolation::Linear);
        let a = 0.91 + (5. - 4.) * (0.8 - 0.91) / 6.;
        let b = 0.7 + 5. * (0.9 - 0.7) / 6.;
        assert!((t.efficiency(5., 150.).unwrap() - (a + b) / 2.).abs() < 1e-14);
        let limits = t.flow_limits(150.).unwrap();
        assert_eq!(limits.min, 1.5);
        assert_eq!(limits.max, 14.);
        // At an exact supplied head only that head's original discharge curve is needed.
        assert!(t.efficiency(15., 200.).is_ok());
        assert!(t.efficiency(15., 150.).is_err());
        assert!(t.efficiency(5., 50.).is_err());
    }
    #[test]
    fn cubic_interpolation_preserves_knots_and_does_not_overshoot() {
        let t = table(DischargeInterpolation::Pchip);
        for (&head, c) in t.heads().iter().zip(t.curves()) {
            for (&q, &e) in c.coordinates().iter().zip(c.values()) {
                assert!((t.efficiency(q, head).unwrap() - e).abs() < 1e-14);
            }
            for i in 0..c.coordinates().len() - 1 {
                for k in 0..101 {
                    let q = c.coordinates()[i]
                        + (c.coordinates()[i + 1] - c.coordinates()[i]) * k as f64 / 100.;
                    let e = t.efficiency(q, head).unwrap();
                    assert!(
                        e >= c.values()[i].min(c.values()[i + 1]) - 1e-14
                            && e <= c.values()[i].max(c.values()[i + 1]) + 1e-14
                    );
                }
            }
        }
    }
    #[test]
    fn automatic_generation_derivatives_match_finite_differences() {
        for interpolation in [
            DischargeInterpolation::Linear,
            DischargeInterpolation::Pchip,
        ] {
            let t = table(interpolation);
            for point in [[2.3, 135.], [7.2, 180.], [11., 220.], [-0.5, 80.]] {
                let data = generation_data(&t, point[0], point[1]);
                for axis in 0..2 {
                    let step = 1e-4;
                    let mut p = point;
                    let mut n = point;
                    p[axis] += step;
                    n[axis] -= step;
                    let plus = generation_data(&t, p[0], p[1]);
                    let minus = generation_data(&t, n[0], n[1]);
                    assert!(
                        ((plus.value - minus.value) / (2. * step) - data.gradient[axis]).abs()
                            < 1e-8
                    );
                    for row in 0..2 {
                        assert!(
                            ((plus.gradient[row] - minus.gradient[row]) / (2. * step)
                                - data.hessian[row][axis])
                                .abs()
                                < 1e-8
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn single_head_and_extension_policies_are_explicit() {
        let curve = TableCurve::new(vec![1., 5., 10.], vec![0.8, 0.9, 0.85]).unwrap();
        assert!(
            TurbineTable::new(
                vec![100.],
                vec![curve.clone()],
                vec![Bounds::new(1., 10.)],
                TurbineOptions::default()
            )
            .is_err()
        );
        let t = TurbineTable::new(
            vec![100.],
            vec![curve],
            vec![Bounds::new(1., 10.)],
            TurbineOptions {
                head_extrapolation: Extrapolation::Constant,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            t.efficiency(5., 50.).unwrap(),
            t.efficiency(5., 200.).unwrap()
        );
        assert!(
            (generation_data(&t, 5., 200.).value - 4. * generation_data(&t, 5., 50.).value).abs()
                < 1e-14
        );
        let serialized = serde_json::to_string(&t).unwrap();
        let restored: TurbineTable = serde_json::from_str(&serialized).unwrap();
        assert_eq!(
            t.efficiency(7., 30.).unwrap(),
            restored.efficiency(7., 30.).unwrap()
        );
        assert!(serde_json::from_str::<TurbineTable>(&serialized.replace("0.9", "1.9")).is_err());
    }
    #[test]
    fn linear_discharge_extension_requires_an_explicit_valid_physical_domain() {
        let curve = TableCurve::new(vec![5., 10.], vec![0.8, 0.82]).unwrap();
        let options = TurbineOptions {
            head_extrapolation: Extrapolation::Constant,
            ..Default::default()
        };
        assert!(
            TurbineTable::new(
                vec![100.],
                vec![curve.clone()],
                vec![Bounds::new(0., 20.)],
                options
            )
            .is_err()
        );
        let table = TurbineTable::new(
            vec![100.],
            vec![curve],
            vec![Bounds::new(0., 20.)],
            TurbineOptions {
                discharge_extrapolation: Extrapolation::Linear,
                ..options
            },
        )
        .unwrap();
        assert!((table.efficiency(20., 80.).unwrap() - 0.86).abs() < 1e-14);
        assert!(
            table
                .validate_domain(Bounds::new(0., 20.), Bounds::new(80., 120.))
                .is_ok()
        );
        // A fixed discharge can still have an invalid extended efficiency.
        assert!(
            table
                .validate_domain(Bounds::new(100., 100.), Bounds::new(80., 120.))
                .is_err()
        );
    }

    #[test]
    fn extended_efficiency_is_checked_between_knots_not_only_at_corners() {
        let curves = vec![
            TableCurve::new(vec![0., 5., 10.], vec![0.7, 0.9, 0.7]).unwrap(),
            TableCurve::new(vec![0., 5., 10.], vec![0.6, 0.95, 0.6]).unwrap(),
        ];
        let t = TurbineTable::new(
            vec![100., 200.],
            curves,
            vec![Bounds::new(0., 10.); 2],
            TurbineOptions {
                interpolation: DischargeInterpolation::Pchip,
                head_extrapolation: Extrapolation::Linear,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            t.validate_domain(Bounds::new(0., 10.), Bounds::new(100., 200.))
                .is_ok()
        );
        // At 400 m the endpoint efficiencies are 0.4, but the interior exceeds one.
        assert!(
            t.validate_domain(Bounds::new(0., 10.), Bounds::new(100., 400.))
                .is_err()
        );
    }
}
