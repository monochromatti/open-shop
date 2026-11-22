//! Measured piecewise-linear curves and their exact integral means.

use crate::Error;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "CurveData")]
pub struct TableCurve {
    x: Vec<f64>,
    y: Vec<f64>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CurveData {
    x: Vec<f64>,
    y: Vec<f64>,
}
impl TryFrom<CurveData> for TableCurve {
    type Error = Error;
    fn try_from(data: CurveData) -> Result<Self, Error> {
        Self::new(data.x, data.y)
    }
}

/// A scalar value and its first and second derivatives in two physical inputs.
pub(crate) struct Derivatives2 {
    pub value: f64,
    pub gradient: [f64; 2],
    pub hessian: [[f64; 2]; 2],
}
impl TableCurve {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Result<Self, Error> {
        if x.len() < 2
            || x.len() != y.len()
            || !x.iter().chain(&y).all(|v| v.is_finite())
            || !x
                .windows(2)
                .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())
            || !x
                .windows(2)
                .zip(y.windows(2))
                .all(|(a, b)| ((b[1] - b[0]) / (a[1] - a[0])).is_finite())
        {
            return Err(Error::Input(
                "curve requires finite paired data and strictly increasing coordinates".into(),
            ));
        }
        Ok(Self { x, y })
    }
    pub fn coordinates(&self) -> &[f64] {
        &self.x
    }
    pub fn values(&self) -> &[f64] {
        &self.y
    }
    fn segment(&self, x: f64) -> usize {
        self.x
            .partition_point(|&v| v <= x)
            .saturating_sub(1)
            .min(self.x.len() - 2)
    }
    fn slope(&self, i: usize) -> f64 {
        (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
    }
    pub(crate) fn extended_value(&self, x: f64) -> f64 {
        let i = self.segment(x);
        self.y[i] + (x - self.x[i]) * self.slope(i)
    }
    fn contains(&self, x: f64) -> bool {
        x.is_finite() && self.x[0] <= x && x <= self.x[self.x.len() - 1]
    }
    /// Evaluate within the supplied domain. Knot derivatives use the right segment.
    pub fn value(&self, x: f64) -> Result<f64, Error> {
        if !self.contains(x) {
            return Err(Error::Input("coordinate outside curve domain".into()));
        }
        Ok(self.extended_value(x))
    }
    /// Exact integral average over the volume interval, including equal endpoints.
    /// Coordinates must be within the measured domain; reversing them is allowed.
    pub fn mean_value(&self, before: f64, after: f64) -> Result<f64, Error> {
        if !self.contains(before) || !self.contains(after) {
            return Err(Error::Input("mean endpoints outside curve domain".into()));
        }
        Ok(self.mean_data(before, after).value)
    }
    /// End-secant continuation is used only for nonlinear solver trial points.
    /// Physical bounds on the table domain are enforced by the reservoir model.
    pub(crate) fn mean_data(&self, a: f64, b: f64) -> Derivatives2 {
        if a > b {
            let v = self.mean_data(b, a);
            return Derivatives2 {
                value: v.value,
                gradient: [v.gradient[1], v.gradient[0]],
                hessian: [
                    [v.hessian[1][1], v.hessian[1][0]],
                    [v.hessian[0][1], v.hessian[0][0]],
                ],
            };
        }
        let i = self.segment(a);
        let j = self.segment(b);
        if i == j || a == b {
            let slope = self.slope(i);
            return Derivatives2 {
                value: self.y[i] + slope * ((a - self.x[i]) + (b - self.x[i])) / 2.,
                gradient: [slope / 2.; 2],
                hessian: [[0.; 2]; 2],
            };
        }
        let width = b - a;
        let mut gradient = [0.; 2];
        let mut left = a;
        for k in i..=j {
            let right = if k == j { b } else { self.x[k + 1] };
            let weight = (right - left) / width;
            let centre = ((left - a) / width + (right - a) / width) / 2.;
            let slope = self.slope(k);
            gradient[0] += slope * weight * (1. - centre);
            gradient[1] += slope * weight * centre;
            left = right;
        }
        let cross = (gradient[1] - gradient[0]) / width;
        Derivatives2 {
            value: self.extended_value(a) + width * gradient[0],
            gradient,
            hessian: [
                [(2. * gradient[0] - self.slope(i)) / width, cross],
                [cross, (self.slope(j) - 2. * gradient[1]) / width],
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn measured_curve_integral_is_exact_symmetric_and_stable_at_equal_volumes() {
        let curve = TableCurve::new(vec![0., 1., 3.], vec![10., 12., 20.]).unwrap();
        // Areas are 11 on [0,1] and 32 on [1,3].
        assert!((curve.mean_value(0., 3.).unwrap() - 43. / 3.).abs() < 1e-14);
        assert_eq!(
            curve.mean_value(0., 3.).unwrap(),
            curve.mean_value(3., 0.).unwrap()
        );
        for x in [0., 0.5, 1., 2., 3.] {
            assert_eq!(curve.mean_value(x, x).unwrap(), curve.value(x).unwrap());
        }
        let a = 2.0_f64;
        let b = f64::from_bits(a.to_bits() + 1);
        assert!((curve.mean_value(a, b).unwrap() - 16.).abs() < 1e-14);
    }
    #[test]
    fn integral_derivatives_include_constraint_curvature_and_reverse_endpoints() {
        let curve = TableCurve::new(vec![0., 1., 3., 4.], vec![10., 12., 20., 21.]).unwrap();
        for (a, b) in [(0.2, 2.4), (2.4, 0.2), (1.2, 2.4), (0.4, 3.8)] {
            let result = curve.mean_data(a, b);
            let step = 1e-5;
            for axis in 0..2 {
                let mut plus = [a, b];
                plus[axis] += step;
                let mut minus = [a, b];
                minus[axis] -= step;
                let p = curve.mean_data(plus[0], plus[1]);
                let m = curve.mean_data(minus[0], minus[1]);
                assert!(((p.value - m.value) / (2. * step) - result.gradient[axis]).abs() < 1e-8);
                for row in 0..2 {
                    assert!(
                        ((p.gradient[row] - m.gradient[row]) / (2. * step)
                            - result.hessian[row][axis])
                            .abs()
                            < 1e-8
                    );
                }
            }
        }
    }
    #[test]
    fn invalid_curve_data_and_extrapolation_are_rejected() {
        assert!(TableCurve::new(vec![1., 1.], vec![1., 2.]).is_err());
        let curve = TableCurve::new(vec![1., 2.], vec![3., 4.]).unwrap();
        assert!(curve.value(0.).is_err());
        assert!(curve.mean_value(1., 3.).is_err());
        assert!(serde_json::from_str::<TableCurve>(r#"{"x":[0,0],"y":[1,2]}"#).is_err());
    }
}
