//! Hydraulic work balances, pump characteristics and electrical conversion.

use crate::{Error, Extrapolation, TableCurve, Unit, UnitMode, curve::Derivatives2};
use num_dual::DualNum;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "PumpData", into = "PumpData")]
pub struct PumpCharacteristic {
    data: PumpData,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PumpData {
    heads: Vec<f64>,
    discharge: Vec<f64>,
    efficiency: Vec<f64>,
    #[serde(default)]
    head_extrapolation: Extrapolation,
}
impl From<PumpCharacteristic> for PumpData {
    fn from(c: PumpCharacteristic) -> Self {
        c.data
    }
}
impl TryFrom<PumpData> for PumpCharacteristic {
    type Error = Error;
    fn try_from(d: PumpData) -> Result<Self, Error> {
        Self::new(d.heads, d.discharge, d.efficiency, d.head_extrapolation)
    }
}
impl PumpCharacteristic {
    pub fn new(
        heads: Vec<f64>,
        discharge: Vec<f64>,
        efficiency: Vec<f64>,
        head_extrapolation: Extrapolation,
    ) -> Result<Self, Error> {
        if heads.is_empty()
            || heads.len() != discharge.len()
            || heads.len() != efficiency.len()
            || !heads
                .iter()
                .chain(&discharge)
                .all(|v| v.is_finite() && *v > 0.)
            || !efficiency
                .iter()
                .all(|e| e.is_finite() && *e > 0. && *e <= 1.)
            || !heads
                .windows(2)
                .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())
            || (heads.len() == 1 && head_extrapolation != Extrapolation::Constant)
            || (heads.len() > 1 && head_extrapolation == Extrapolation::Constant)
        {
            return Err(Error::Input("pump characteristic needs ordered positive heads, positive flow/efficiency and an explicit single-head policy".into()));
        }
        if heads.len() > 1
            && [&discharge, &efficiency].iter().any(|v| {
                v.windows(2)
                    .zip(heads.windows(2))
                    .any(|(y, x)| !((y[1] - y[0]) / (x[1] - x[0])).is_finite())
            })
        {
            return Err(Error::Input("nonfinite pump characteristic slope".into()));
        }
        Ok(Self {
            data: PumpData {
                heads,
                discharge,
                efficiency,
                head_extrapolation,
            },
        })
    }
    pub fn heads(&self) -> &[f64] {
        &self.data.heads
    }
    fn check(&self, h: f64) -> Result<(), Error> {
        if !h.is_finite()
            || (self.data.head_extrapolation == Extrapolation::Error
                && (h < self.data.heads[0] || h > *self.data.heads.last().unwrap()))
        {
            return Err(Error::Input("head outside pump characteristic".into()));
        }
        Ok(())
    }
    fn coefficients(&self, h: f64, efficiency: bool) -> (f64, f64, f64) {
        let y = if efficiency {
            &self.data.efficiency
        } else {
            &self.data.discharge
        };
        let x = &self.data.heads;
        if x.len() == 1 {
            return (x[0], y[0], 0.);
        }
        let i = x
            .partition_point(|v| *v <= h)
            .saturating_sub(1)
            .min(x.len() - 2);
        (x[i], y[i], (y[i + 1] - y[i]) / (x[i + 1] - x[i]))
    }
    pub fn discharge(&self, h: f64) -> Result<f64, Error> {
        self.check(h)?;
        Ok(self.trial_discharge(h))
    }
    pub fn efficiency(&self, h: f64) -> Result<f64, Error> {
        self.check(h)?;
        Ok(self.trial_efficiency(h))
    }
    pub(crate) fn trial_discharge(&self, h: f64) -> f64 {
        let (x, y, s) = self.coefficients(h, false);
        y + (h - x) * s
    }
    pub(crate) fn discharge_data(&self, h: f64) -> Derivatives2 {
        let (x, y, s) = self.coefficients(h, false);
        Derivatives2 {
            value: y + (h - x) * s,
            gradient: [0., s],
            hessian: [[0.; 2]; 2],
        }
    }
    pub(crate) fn trial_efficiency<D: DualNum<Primitive = f64> + Copy>(&self, h: D) -> D {
        let (x, y, s) = self.coefficients(h.re(), true);
        (h - x) * s + y
    }
}

pub(crate) fn table_value<D: DualNum<Primitive = f64> + Copy>(c: &TableCurve, x: D) -> D {
    let coordinates = c.coordinates();
    let values = c.values();
    let i = coordinates
        .partition_point(|v| *v <= x.re())
        .saturating_sub(1)
        .min(coordinates.len() - 2);
    (x - coordinates[i]) * ((values[i + 1] - values[i]) / (coordinates[i + 1] - coordinates[i]))
        + values[i]
}

pub(crate) fn hydraulic_efficiency<D: DualNum<Primitive = f64> + Copy>(u: &Unit, q: D, h: D) -> D {
    if let Some(c) = &u.pump_characteristic {
        return c.trial_efficiency(h);
    }
    match &u.efficiency {
        crate::Efficiency::Quadratic {
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
            q,
            h,
        ),
        crate::Efficiency::Table(c) => c.trial_efficiency(q, h),
    }
}

/// Power and flow are positive magnitudes. Pumping adds work to the water;
/// generation extracts work. No reciprocal efficiency occurs in this residual.
pub(crate) fn power_balance<D: DualNum<Primitive = f64> + Copy>(u: &Unit, q: D, h: D, p: D) -> D {
    let hydraulic = q * h * 0.00981;
    let eta = hydraulic_efficiency(u, q, h);
    match u.mode {
        UnitMode::Pump => hydraulic - p * eta * u.motor_efficiency,
        UnitMode::Generate => {
            let electrical = u
                .generator_efficiency_curve
                .as_ref()
                .map_or(D::from(1.), |c| table_value(c, p));
            hydraulic * eta * electrical - p
        }
    }
}

pub(crate) fn electrical_power(c: &TableCurve, shaft: f64) -> Result<f64, Error> {
    if shaft == 0. {
        return Ok(0.);
    }
    let x = c.coordinates();
    let y = c.values();
    let i = x
        .iter()
        .zip(y)
        .position(|(&p, &e)| p / e > shaft)
        .unwrap_or(x.len() - 1)
        .saturating_sub(1)
        .min(x.len() - 2);
    let slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    let intercept = y[i] - slope * x[i];
    let denominator = 1. - shaft * slope;
    if !shaft.is_finite() || denominator <= 0. {
        return Err(Error::Input(
            "generator efficiency extension has no finite positive root".into(),
        ));
    }
    Ok(shaft * intercept / denominator)
}

pub(crate) fn validate_unit(u: &Unit) -> Result<(), Error> {
    let invalid = |message: &str| Error::Input(format!("{}: {message}", u.name));
    if !u.motor_efficiency.is_finite() || !(0. < u.motor_efficiency && u.motor_efficiency <= 1.) {
        return Err(invalid("motor efficiency must be in (0,1]"));
    }
    match u.mode {
        UnitMode::Pump => {
            if u.head.min <= 0. || u.generator_efficiency_curve.is_some() {
                return Err(invalid(
                    "pumping needs positive head and constant motor efficiency",
                ));
            }
        }
        UnitMode::Generate => {
            if u.motor_efficiency != 1. || u.pump_characteristic.is_some() {
                return Err(invalid("pump parameters cannot apply to a generating unit"));
            }
        }
    }
    if let Some(c) = &u.pump_characteristic {
        if matches!(u.efficiency, crate::Efficiency::Table(_)) {
            return Err(invalid(
                "supply a fixed-speed characteristic or a variable-flow efficiency table",
            ));
        }
        c.check(u.head.min)?;
        c.check(u.head.max)?;
        for h in c
            .heads()
            .iter()
            .copied()
            .filter(|h| u.head.min < *h && *h < u.head.max)
            .chain([u.head.min, u.head.max])
        {
            let q = c.trial_discharge(h);
            let e = c.trial_efficiency(h);
            if !q.is_finite() || q <= 0. || !e.is_finite() || e <= 0. || e > 1. {
                return Err(invalid(
                    "extended pump characteristic has invalid flow or efficiency",
                ));
            }
        }
    }
    if let Some(c) = &u.generator_efficiency_curve {
        let x = c.coordinates();
        let y = c.values();
        if x[0] != 0.
            || *x.last().unwrap() < u.power.max
            || y.iter().any(|e| *e <= 0. || *e > 1.)
            || (0..x.len() - 1).any(|i| y[i] - x[i] * (y[i + 1] - y[i]) / (x[i + 1] - x[i]) <= 0.)
        {
            return Err(invalid(
                "electrical efficiency must cover [0,p_max] and give increasing required shaft power",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Bounds, Efficiency};
    use nalgebra::SVector;
    use num_dual::hessian;
    #[test]
    fn work_balance_derivatives_include_power_and_efficiency_couplings() {
        let mut unit = Unit::new(
            "Machine",
            "Station",
            Bounds::new(1., 20.),
            Bounds::new(0., 30.),
            Bounds::new(80., 120.),
            Efficiency::Quadratic {
                peak: 0.9,
                best_flow: 12.,
                flow_curvature: 0.001,
                best_head: 100.,
                head_curvature: 0.00001,
            },
        );
        for mode in [UnitMode::Generate, UnitMode::Pump] {
            unit.mode = mode;
            unit.motor_efficiency = if mode == UnitMode::Pump { 0.95 } else { 1. };
            unit.generator_efficiency_curve = if mode == UnitMode::Generate {
                Some(TableCurve::new(vec![0., 15., 30.], vec![0.85, 0.94, 0.97]).unwrap())
            } else {
                None
            };
            for point in [[7.0_f64, 93., 8.], [15., 107., 22.]] {
                let (value, gradient, h) = hessian(
                    |x| power_balance(&unit, x[0], x[1], x[2]),
                    &SVector::from(point),
                );
                assert!((value - power_balance(&unit, point[0], point[1], point[2])).abs() < 1e-12);
                for axis in 0..3 {
                    let mut plus = point;
                    let mut minus = point;
                    plus[axis] += 1e-4;
                    minus[axis] -= 1e-4;
                    let (p, gp, _) = hessian(
                        |x| power_balance(&unit, x[0], x[1], x[2]),
                        &SVector::from(plus),
                    );
                    let (n, gn, _) = hessian(
                        |x| power_balance(&unit, x[0], x[1], x[2]),
                        &SVector::from(minus),
                    );
                    assert!(((p - n) / 2e-4 - gradient[axis]).abs() < 1e-8);
                    for row in 0..3 {
                        assert!(((gp[row] - gn[row]) / 2e-4 - h[(row, axis)]).abs() < 1e-8);
                    }
                }
            }
        }
    }
}
