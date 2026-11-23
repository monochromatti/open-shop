//! Exact machine operating limits over a continuous span of net heads.
//!
//! On each calibration segment, efficiency is quadratic at most in head.
//! Hydraulic work and the transformed power inequalities are cubic at most.
//! Their extrema occur at endpoints, data knots, or derivative roots. This
//! requires no extra hydraulic time steps and retains derivatives of the
//! stationary head when the held discharge changes.

use crate::{Bounds, Efficiency, Unit, UnitMode};
use num_dual::DualNum;
use std::sync::Arc;

#[derive(Clone, Copy)]
pub(crate) enum HeadLimit {
    UpperPower,
    LowerPower,
    Efficiency,
    Flow,
}
pub(crate) struct HeadEnvelope {
    pub unit: Arc<Unit>,
    pub power: Bounds,
    pub flow: Bounds,
}

/// Retain only head-span inequalities not already implied by the endpoint
/// equations and bounds. Polynomial ranges cover the complete declared box,
/// including original calibration cells and explicit extrapolation segments.
pub(crate) fn required_constraints(u: &Unit) -> [bool; 4] {
    let (mut emin, mut emax, mut increasing) = (f64::INFINITY, f64::NEG_INFINITY, true);
    let monotone = |values: Vec<f64>| {
        values.windows(2).all(|v| v[0] <= v[1]) || values.windows(2).all(|v| v[0] >= v[1])
    };
    let flow_monotone;
    if let Some(c) = &u.pump_characteristic {
        increasing = false;
        for h in [u.head.min, u.head.max]
            .into_iter()
            .chain(c.heads().iter().copied().filter(|h| u.head.contains(*h)))
        {
            let e = c.trial_efficiency(h);
            emin = emin.min(e);
            emax = emax.max(e);
        }
        flow_monotone = c.monotone_flow();
    } else if let Efficiency::Table(table) = &u.efficiency {
        let mut qs: Vec<_> = table
            .curves()
            .iter()
            .flat_map(|c| c.coordinates().iter().copied())
            .filter(|q| u.flow.min < *q && *q < u.flow.max)
            .chain([u.flow.min, u.flow.max])
            .collect();
        let mut hs: Vec<_> = table
            .heads()
            .iter()
            .copied()
            .filter(|h| u.head.min < *h && *h < u.head.max)
            .chain([u.head.min, u.head.max])
            .collect();
        for values in [&mut qs, &mut hs] {
            values.sort_by(f64::total_cmp);
            values.dedup();
            if values.len() == 1 {
                values.push(values[0]);
            }
        }
        for h in hs.windows(2) {
            let (j, k, _) = table.head_weights(h[0]);
            for q in qs.windows(2) {
                let polynomial = |col| {
                    let (origin, width, c) = table.column_polynomial(col, q[0] / 2. + q[1] / 2.);
                    let a = (q[0] - origin) / width;
                    let r = (q[1] - q[0]) / width;
                    [
                        ((c[3] * a + c[2]) * a + c[1]) * a + c[0],
                        r * ((3. * c[3] * a + 2. * c[2]) * a + c[1]),
                        r * r * (c[2] + 3. * c[3] * a),
                        r * r * r * c[3],
                    ]
                };
                let left = polynomial(j);
                let right = polynomial(k);
                let slope: [f64; 4] = std::array::from_fn(|i| {
                    if j == k {
                        0.
                    } else {
                        (right[i] - left[i]) / (table.heads()[k] - table.heads()[j])
                    }
                });
                for &head in h {
                    let eta =
                        std::array::from_fn(|i| left[i] + (head - table.heads()[j]) * slope[i]);
                    let (a, b, _, _) = cubic_range(eta, 0., 1.);
                    emin = emin.min(a);
                    emax = emax.max(b);
                    let factor = match u.mode {
                        UnitMode::Generate => 2. * head - table.heads()[j],
                        UnitMode::Pump => -table.heads()[j],
                    };
                    let derivative = std::array::from_fn(|i| left[i] + factor * slope[i]);
                    let (lower, _, _, _) = cubic_range(derivative, 0., 1.);
                    increasing &= lower
                        > 128.
                            * f64::EPSILON
                            * derivative.iter().map(|v| v.abs()).sum::<f64>().max(1.);
                }
            }
        }
        flow_monotone = monotone(table.flow_envelopes().iter().map(|b| b.min).collect())
            && monotone(table.flow_envelopes().iter().map(|b| b.max).collect());
    } else {
        let Efficiency::Quadratic {
            peak,
            best_flow,
            flow_curvature,
            best_head,
            head_curvature,
        } = u.efficiency
        else {
            unreachable!()
        };
        for q in [u.flow.min, u.flow.max] {
            for h in [u.head.min, u.head.max] {
                let e = u.efficiency.trial_value(q, h);
                emin = emin.min(e);
            }
        }
        emax = u.efficiency.trial_value(
            best_flow.clamp(u.flow.min, u.flow.max),
            best_head.clamp(u.head.min, u.head.max),
        );
        for q in [u.flow.min, u.flow.max] {
            let e0 = peak
                - flow_curvature * (q - best_flow).powi(2)
                - head_curvature * best_head * best_head;
            let e1 = 2. * head_curvature * best_head;
            let e2 = -head_curvature;
            let derivative = match u.mode {
                UnitMode::Generate => [e0, 2. * e1, 3. * e2, 0.],
                UnitMode::Pump => [e0, 0., -e2, 0.],
            };
            let (lower, _, _, _) = cubic_range(derivative, u.head.min, u.head.max);
            increasing &= lower
                > 128. * f64::EPSILON * derivative.iter().map(|v| v.abs()).sum::<f64>().max(1.);
        }
        flow_monotone = true;
    }
    let guard = 128. * f64::EPSILON * emin.abs().max(emax.abs()).max(1.);
    let efficiency_bounded = emin > u.minimum_efficiency + guard && emax < 1. - guard;
    [
        !increasing,
        !increasing,
        !efficiency_bounded,
        !flow_monotone,
    ]
}
pub(crate) struct HeadRange<D> {
    pub upper: D,
    pub lower: D,
    pub efficiency: D,
    pub flow: D,
    pub at_upper: D,
    pub at_lower: D,
}
impl<D: Copy> HeadRange<D> {
    pub fn component(&self, component: HeadLimit) -> D {
        match component {
            HeadLimit::UpperPower => self.upper,
            HeadLimit::LowerPower => self.lower,
            HeadLimit::Efficiency => self.efficiency,
            HeadLimit::Flow => self.flow,
        }
    }
}

fn max<D: DualNum<Primitive = f64> + Copy>(a: D, b: D) -> D {
    if a.re() >= b.re() { a } else { b }
}
fn quadratic_roots<D: DualNum<Primitive = f64> + Copy>(a: D, b: D, c: D) -> [Option<D>; 2] {
    if a.re() == 0. {
        return [if b.re() == 0. { None } else { Some(-c / b) }, None];
    }
    let discriminant = b * b - a * c * 4.;
    // A double derivative root has no strict extremum. Avoid differentiating
    // sqrt at zero; the endpoints already give the range in that case.
    if discriminant.re() <= 0. {
        return [None, None];
    }
    let q = (b + discriminant.sqrt() * if b.re() < 0. { -1. } else { 1. }) * -0.5;
    [Some(q / a), Some(c / q)]
}

fn cubic_range<D: DualNum<Primitive = f64> + Copy>(c: [D; 4], a: D, b: D) -> (D, D, D, D) {
    let value = |h: D| ((c[3] * h + c[2]) * h + c[1]) * h + c[0];
    let mut lo = value(a);
    let mut hi = lo;
    let (mut at_lo, mut at_hi) = (a, a);
    let roots = quadratic_roots(c[3] * 3., c[2] * 2., c[1]);
    for h in [Some(b), roots[0], roots[1]].into_iter().flatten() {
        if !h.re().is_finite() || h.re() < a.re() || h.re() > b.re() {
            continue;
        }
        let v = value(h);
        if v.re() < lo.re() {
            lo = v;
            at_lo = h;
        }
        if v.re() > hi.re() {
            hi = v;
            at_hi = h;
        }
    }
    (lo, hi, at_lo, at_hi)
}

impl HeadEnvelope {
    pub fn range<D: DualNum<Primitive = f64> + Copy>(&self, q: D, a: D, b: D) -> HeadRange<D> {
        let u = &self.unit;
        let (lo, hi) = if a.re() <= b.re() { (a, b) } else { (b, a) };
        let data = if let Some(c) = &u.pump_characteristic {
            c.heads()
        } else if let Efficiency::Table(c) = &u.efficiency {
            c.heads()
        } else {
            &[]
        };
        let mut result = HeadRange {
            upper: D::from(f64::NEG_INFINITY),
            lower: D::from(f64::NEG_INFINITY),
            efficiency: D::from(f64::NEG_INFINITY),
            flow: D::from(f64::NEG_INFINITY),
            at_upper: lo,
            at_lower: lo,
        };
        let mut left = lo;
        for right in data
            .iter()
            .copied()
            .filter(|h| lo.re() < *h && *h < hi.re())
            .map(D::from)
            .chain(std::iter::once(hi))
        {
            let (q0, q1, e) = if let Some(c) = &u.pump_characteristic {
                let (h, y, s) = c.coefficients(left.re(), false);
                let (he, ye, se) = c.coefficients(left.re(), true);
                (
                    D::from(y - s * h),
                    D::from(s),
                    [D::from(ye - se * he), D::from(se), D::from(0.)],
                )
            } else if let Efficiency::Table(c) = &u.efficiency {
                for h in [left, right] {
                    let minimum = c.envelope_data(h.re(), false);
                    let maximum = c.envelope_data(h.re(), true);
                    let low = (h - h.re()) * minimum.gradient[1] + minimum.value;
                    let high = (h - h.re()) * maximum.gradient[1] + maximum.value;
                    result.flow = max(result.flow, max(low - q, q - high));
                }
                (q, D::from(0.), c.head_polynomial(q, left.re()))
            } else {
                let Efficiency::Quadratic {
                    peak,
                    best_flow,
                    flow_curvature,
                    best_head,
                    head_curvature,
                } = u.efficiency
                else {
                    unreachable!()
                };
                let dq = q - best_flow;
                (
                    q,
                    D::from(0.),
                    [
                        dq * dq * (-flow_curvature) + peak - head_curvature * best_head * best_head,
                        D::from(2. * head_curvature * best_head),
                        D::from(-head_curvature),
                    ],
                )
            };
            let (emin, emax, _, _) = cubic_range([e[0], e[1], e[2], D::from(0.)], left, right);
            result.efficiency = max(
                result.efficiency,
                max(D::from(u.minimum_efficiency) - emin, emax - 1.),
            );
            for h in [left, right] {
                let flow = q0 + q1 * h;
                result.flow = max(
                    result.flow,
                    max(D::from(self.flow.min) - flow, flow - self.flow.max),
                );
            }
            let (upper, lower, at_upper, at_lower) = match u.mode {
                UnitMode::Generate => {
                    let (minimum, maximum, at_min, at_max) = cubic_range(
                        [
                            D::from(0.),
                            q * e[0] * 0.00981,
                            q * e[1] * 0.00981,
                            q * e[2] * 0.00981,
                        ],
                        left,
                        right,
                    );
                    let electrical = |p| {
                        u.generator_efficiency_curve
                            .as_ref()
                            .map_or(1., |c| c.extended_value(p))
                    };
                    (
                        maximum - self.power.max / electrical(self.power.max),
                        D::from(self.power.min / electrical(self.power.min)) - minimum,
                        at_max,
                        at_min,
                    )
                }
                UnitMode::Pump => {
                    let up = self.power.max * u.motor_efficiency;
                    let down = self.power.min * u.motor_efficiency;
                    let (_, maximum, _, at_max) = cubic_range(
                        [
                            -e[0] * up,
                            q0 * 0.00981 - e[1] * up,
                            q1 * 0.00981 - e[2] * up,
                            D::from(0.),
                        ],
                        left,
                        right,
                    );
                    let (_, minimum, _, at_min) = cubic_range(
                        [
                            e[0] * down,
                            e[1] * down - q0 * 0.00981,
                            e[2] * down - q1 * 0.00981,
                            D::from(0.),
                        ],
                        left,
                        right,
                    );
                    (maximum, minimum, at_max, at_min)
                }
            };
            if upper.re() > result.upper.re() {
                result.upper = upper;
                result.at_upper = at_upper;
            }
            if lower.re() > result.lower.re() {
                result.lower = lower;
                result.at_lower = at_lower;
            }
            left = right;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::SVector;
    use num_dual::hessian;

    #[test]
    fn interior_head_peak_and_its_derivatives_include_the_moving_stationary_point() {
        let unit = Unit::new(
            "Unit",
            "Plant",
            Bounds::new(1., 20.),
            Bounds::new(0., 20.),
            Bounds::new(1., 200.),
            Efficiency::Quadratic {
                peak: 0.9,
                best_flow: 10.,
                flow_curvature: 0.001,
                best_head: 100.,
                head_curvature: 0.0001,
            },
        );
        let envelope = HeadEnvelope {
            unit: Arc::new(unit),
            power: Bounds::new(0., 20.),
            flow: Bounds::new(1., 20.),
        };
        let x = SVector::from([8., 80., 160.]);
        let f = |x: SVector<f64, 3>| envelope.range(x[0], x[1], x[2]).upper;
        let (value, gradient, h) = hessian(|x| envelope.range(x[0], x[1], x[2]).upper, &x);
        let range = envelope.range(x[0], x[1], x[2]);
        assert!(range.at_upper > 100. && range.at_upper < 160.);
        assert!((value - f(x)).abs() < 1e-12);
        let step = 0.001;
        let mut a = x;
        let mut b = x;
        a[0] += step;
        b[0] -= step;
        assert!((gradient[0] - (f(a) - f(b)) / (2. * step)).abs() < 1e-8);
        assert!((h[(0, 0)] - (f(a) - 2. * f(x) + f(b)) / (step * step)).abs() < 1e-6);
        assert!(gradient[1].abs() + gradient[2].abs() < 1e-12);
    }
}
