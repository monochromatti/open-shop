//! Passive level-difference channels and reservoir outlet laws.

use crate::{Error, TableCurve, curve::Derivatives2};
use num_dual::DualNum;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ReferenceLevel {
    #[default]
    Upstream,
    Downstream,
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RiverLaw {
    #[default]
    Controlled,
    Orifice {
        coefficient: f64,
        crest: f64,
        #[serde(default)]
        allow_dry: bool,
    },
    Weir {
        coefficient: f64,
        crest: f64,
        #[serde(default)]
        allow_dry: bool,
    },
    Level {
        curve: TableCurve,
    },
    Delta {
        curve: LevelFlowCurve,
        #[serde(default)]
        reference: ReferenceLevel,
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(try_from = "LevelFlowData", into = "LevelFlowData")]
pub struct LevelFlowCurve {
    data: LevelFlowData,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LevelFlowData {
    reference_levels: Vec<f64>,
    curves: Vec<TableCurve>,
}
impl From<LevelFlowCurve> for LevelFlowData {
    fn from(c: LevelFlowCurve) -> Self {
        c.data
    }
}
impl TryFrom<LevelFlowData> for LevelFlowCurve {
    type Error = Error;
    fn try_from(d: LevelFlowData) -> Result<Self, Error> {
        Self::new(d.reference_levels, d.curves)
    }
}
impl LevelFlowCurve {
    pub fn new(reference_levels: Vec<f64>, curves: Vec<TableCurve>) -> Result<Self, Error> {
        if reference_levels.is_empty()
            || reference_levels.len() != curves.len()
            || !reference_levels.iter().all(|v| v.is_finite())
            || !reference_levels
                .windows(2)
                .all(|w| w[1] > w[0] && (w[1] - w[0]).is_finite())
        {
            return Err(Error::Input(
                "level-flow curves need paired ordered finite reference levels".into(),
            ));
        }
        for c in &curves {
            if c.values().windows(2).any(|w| w[1] < w[0])
                || !c
                    .coordinates()
                    .iter()
                    .zip(c.values())
                    .any(|(&x, &y)| x == 0. && y == 0.)
                || c.coordinates()
                    .iter()
                    .zip(c.values())
                    .any(|(&x, &y)| (x < 0. && y > 0.) || (x > 0. && y < 0.))
            {
                return Err(Error::Input(
                    "level-difference curves must be nondecreasing, passive, and contain (0,0)"
                        .into(),
                ));
            }
        }
        Ok(Self {
            data: LevelFlowData {
                reference_levels,
                curves,
            },
        })
    }
    pub fn reference_levels(&self) -> &[f64] {
        &self.data.reference_levels
    }
    pub fn curves(&self) -> &[TableCurve] {
        &self.data.curves
    }
    fn reference(&self, h: f64) -> (usize, usize, f64, f64) {
        let refs = &self.data.reference_levels;
        let n = refs.len();
        if n == 1 || h < refs[0] {
            return (0, 0, 0., 0.);
        }
        if h > refs[n - 1] {
            return (n - 1, n - 1, 0., 0.);
        }
        let j = refs
            .partition_point(|v| *v <= h)
            .saturating_sub(1)
            .min(n - 2);
        let width = refs[j + 1] - refs[j];
        (j, j + 1, (h - refs[j]) / width, 1. / width)
    }
    /// Strict delta-domain lookup; reference levels outside data hold the nearest curve.
    pub fn value(&self, delta: f64, reference: f64) -> Result<f64, Error> {
        if !delta.is_finite() || !reference.is_finite() {
            return Err(Error::Input("nonfinite level-flow coordinate".into()));
        }
        let (j, k, w, _) = self.reference(reference);
        if j == k || w == 0. {
            return self.data.curves[j].value(delta);
        }
        if w == 1. {
            return self.data.curves[k].value(delta);
        }
        Ok(self.data.curves[j].value(delta)? * (1. - w) + self.data.curves[k].value(delta)? * w)
    }
    pub(crate) fn extended_value<D: DualNum<Primitive = f64> + Copy>(
        &self,
        delta: D,
        reference: D,
    ) -> D {
        let (j, k, _, _) = self.reference(reference.re());
        let left = crate::machine::table_value(&self.data.curves[j], delta);
        if j == k {
            return left;
        }
        let refs = &self.data.reference_levels;
        let weight = (reference - refs[j]) / (refs[k] - refs[j]);
        left + weight * (crate::machine::table_value(&self.data.curves[k], delta) - left)
    }
    pub(crate) fn data(&self, delta: f64, reference: f64) -> Derivatives2 {
        let (j, k, w, scale) = self.reference(reference);
        let column = |index: usize| {
            let c = &self.data.curves[index];
            let x = c.coordinates();
            let i = x
                .partition_point(|v| *v <= delta)
                .saturating_sub(1)
                .min(x.len() - 2);
            (
                c.extended_value(delta),
                (c.values()[i + 1] - c.values()[i]) / (x[i + 1] - x[i]),
            )
        };
        let (left, dl) = column(j);
        let (right, dr) = column(k);
        let cross = (dr - dl) * scale;
        Derivatives2 {
            value: left + w * (right - left),
            gradient: [dl + w * (dr - dl), (right - left) * scale],
            hessian: [[0., cross], [cross, 0.]],
        }
    }
}
impl RiverLaw {
    /// Outlet withdrawal in m³/s; delta channels use `LevelFlowCurve` instead.
    pub fn outlet_flow(
        &self,
        head: f64,
        opening: f64,
        capacity: f64,
        inflow: f64,
    ) -> Result<f64, Error> {
        if !head.is_finite()
            || !opening.is_finite()
            || !(0. ..=1.).contains(&opening)
            || !capacity.is_finite()
            || !inflow.is_finite()
            || inflow < 0.
            || capacity < inflow
        {
            return Err(Error::Input("invalid river outlet coordinate".into()));
        }
        Ok(match self {
            Self::Controlled => (capacity - inflow) * opening,
            Self::Orifice {
                coefficient, crest, ..
            } => coefficient * opening * (head - crest).max(0.).sqrt(),
            Self::Weir {
                coefficient, crest, ..
            } => coefficient * (head - crest).max(0.).powf(1.5),
            Self::Level { curve } => opening * curve.value(head)?,
            Self::Delta { .. } => {
                return Err(Error::Input(
                    "a signed channel has no gate-controlled outlet".into(),
                ));
            }
        })
    }
    // Numerical evaluation never clips a returned control. Bounds are checked
    // separately, so solver roundoff remains visible in the residuals.
    pub(crate) fn trial_outlet_flow(
        &self,
        head: f64,
        opening: f64,
        capacity: f64,
        inflow: f64,
    ) -> f64 {
        match self {
            Self::Controlled => (capacity - inflow) * opening,
            Self::Orifice {
                coefficient, crest, ..
            } => coefficient * opening * (head - crest).max(0.).sqrt(),
            Self::Weir {
                coefficient, crest, ..
            } => coefficient * (head - crest).max(0.).powf(1.5),
            Self::Level { curve } => opening * curve.extended_value(head),
            Self::Delta { .. } => f64::NAN,
        }
    }
    pub(crate) fn validate(
        &self,
        river: &crate::River,
        system: &crate::HydroSystem,
    ) -> Result<(), Error> {
        river.gate.validate(&river.name)?;
        if river.gate.min < 0. || river.gate.max > 1. {
            return Err(Error::Input("river gate bounds must lie in [0,1]".into()));
        }
        let source = system.reservoirs.iter().find(|r| r.name == river.source);
        let bad = |message: &str| Error::Input(format!("{}: {message}", river.name));
        match self {
            Self::Controlled => {}
            Self::Delta { .. } => {
                if source.is_none()
                    || !system.reservoirs.iter().any(|r| r.name == river.target)
                    || river.transport.fixed_delay() != Some(0.)
                    || river.inflow != 0.
                    || river.history_release.iter().any(|&q| q != 0.)
                    || river.gate.min != 0.
                    || river.gate.max != 1.
                    || system
                        .tunnels
                        .iter()
                        .any(|t| t.discharge_river.as_deref() == Some(&river.name))
                    || system
                        .plants
                        .iter()
                        .any(|p| p.discharge_river.as_deref() == Some(&river.name))
                {
                    return Err(bad(
                        "a signed channel requires two reservoirs, zero delay and no inflow, outfall, inventory or gate restrictions",
                    ));
                }
            }
            Self::Orifice {
                coefficient,
                crest,
                allow_dry,
            }
            | Self::Weir {
                coefficient,
                crest,
                allow_dry,
            } => {
                let r = source.ok_or_else(|| bad("an outlet law requires a reservoir source"))?;
                if !coefficient.is_finite()
                    || *coefficient <= 0.
                    || !crest.is_finite()
                    || (!allow_dry && r.level.trial_level(r.storage.min) <= *crest)
                {
                    return Err(bad(
                        "invalid outlet coefficient/crest, or possible drying without allow_dry",
                    ));
                }
                if matches!(self, Self::Weir { .. }) && river.gate.max < 1. {
                    return Err(bad("a weir cannot be closed by a gate"));
                }
            }
            Self::Level { curve } => {
                let r = source.ok_or_else(|| bad("an outlet curve requires a reservoir source"))?;
                if curve.values().iter().any(|v| *v < 0.)
                    || curve.values().windows(2).any(|w| w[1] < w[0])
                    || curve.value(r.level.trial_level(r.storage.min)).is_err()
                    || curve.value(r.level.trial_level(r.storage.max)).is_err()
                {
                    return Err(bad(
                        "outlet curve must be nonnegative, nondecreasing and cover reservoir levels",
                    ));
                }
            }
        }
        if source.is_none() && (river.gate.min != 0. || river.gate.max != 1.) {
            return Err(bad("river junctions have no independent gates"));
        }
        Ok(())
    }

    pub(crate) fn is_delta(&self) -> bool {
        matches!(self, Self::Delta { .. })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn interpolated_channel_derivatives_match_both_physical_inputs() {
        let c = LevelFlowCurve::new(
            vec![80., 120.],
            vec![
                TableCurve::new(vec![-5., 0., 5.], vec![-4., 0., 10.]).unwrap(),
                TableCurve::new(vec![-5., 0., 5.], vec![-8., 0., 15.]).unwrap(),
            ],
        )
        .unwrap();
        for point in [[-2., 95.], [2., 105.], [8., 130.], [-8., 70.]] {
            let d = c.data(point[0], point[1]);
            for j in 0..2 {
                let h = 1e-4;
                let mut left = point;
                let mut right = point;
                left[j] -= h;
                right[j] += h;
                let fd = (c.extended_value(right[0], right[1])
                    - c.extended_value(left[0], left[1]))
                    / (2. * h);
                assert!((fd - d.gradient[j]).abs() < 1e-8);
                let dl = c.data(left[0], left[1]);
                let dr = c.data(right[0], right[1]);
                for i in 0..2 {
                    assert!(
                        ((dr.gradient[i] - dl.gradient[i]) / (2. * h) - d.hessian[i][j]).abs()
                            < 1e-8
                    );
                }
            }
        }
        assert!(c.value(8., 100.).is_err());
        assert!((c.value(2., 70.).unwrap() - 4.).abs() < 1e-12);
    }
}
