//! Constitutive equations shared by numerical evaluation and equation assembly.
//! Scalar arithmetic permits ordinary values and automatic derivative expressions.
//! There is no solver state, graph preparation or I/O in these functions.

use crate::{Efficiency, StorageCurve};
use std::ops::{Add, Mul, Sub};

pub(crate) trait Scalar:
    Copy
    + Add<Output = Self>
    + Mul<Output = Self>
    + Sub<f64, Output = Self>
    + Add<f64, Output = Self>
    + Mul<f64, Output = Self>
{
}
impl<T> Scalar for T where
    T: Copy
        + Add<Output = T>
        + Mul<Output = T>
        + Sub<f64, Output = T>
        + Add<f64, Output = T>
        + Mul<f64, Output = T>
{
}

pub(crate) fn mean_level<T: Scalar>(curve: &StorageCurve, before: T, after: T) -> T {
    let a = before - curve.reference_volume;
    let b = after - curve.reference_volume;
    (a + b) * (curve.slope / 2.)
        + (a * a + a * b + b * b) * (curve.curvature / 3.)
        + curve.elevation
}

pub(crate) fn efficiency<T: Scalar>(curve: &Efficiency, flow: T, head: T) -> T {
    let dq = flow - curve.best_flow;
    let dh = head - curve.best_head;
    (dq * dq) * (-curve.flow_curvature) + (dh * dh) * (-curve.head_curvature) + curve.peak
}

pub(crate) fn generation<T: Scalar>(flow: T, head: T, efficiency: T) -> T {
    flow * head * efficiency * 0.00981
}
