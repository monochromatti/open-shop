//! Constitutive equations shared by numerical evaluation and equation assembly.
//! Scalar arithmetic permits ordinary values and automatic derivative expressions.
//! There is no solver state, graph preparation or I/O in these functions.

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

pub(crate) fn polynomial_mean_level<T: Scalar>(
    reference_volume: f64,
    elevation: f64,
    slope: f64,
    curvature: f64,
    before: T,
    after: T,
) -> T {
    let a = before - reference_volume;
    let b = after - reference_volume;
    (a + b) * (slope / 2.) + (a * a + a * b + b * b) * (curvature / 3.) + elevation
}

pub(crate) fn quadratic_efficiency<T: Scalar>(
    peak: f64,
    best_flow: f64,
    flow_curvature: f64,
    best_head: f64,
    head_curvature: f64,
    flow: T,
    head: T,
) -> T {
    let dq = flow - best_flow;
    let dh = head - best_head;
    (dq * dq) * (-flow_curvature) + (dh * dh) * (-head_curvature) + peak
}
