//! Sparse assembly of local physical functions composed with affine inputs.
//! Structure and chain-rule coefficients are prepared once. Evaluation only
//! accumulates values and exact derivatives into existing numeric buffers.

use crate::{
    Bounds, Error, TableCurve, TurbineTable, Unit, curve::Derivatives2, river::LinearForm,
};
use nalgebra::SVector;
use num_dual::hessian;
use pounce_nl::nl_reader::NlTnlp;
use pounce_nlp::tnlp::{
    BoundsInfo, IndexStyle, IpoptCq, IpoptData, NlpInfo, Solution, SparsityRequest, StartingPoint,
    TNLP,
};
use std::{collections::BTreeMap, sync::Arc};

pub(crate) struct StorageHeadEquation {
    pub curve: Arc<crate::StorageCurve>,
    /// Fractional time and affine physical volume (m³), including both endpoints.
    pub vertices: Vec<(f64, LinearForm)>,
    pub head: usize,
    pub head_scale: f64,
}

pub(crate) enum PhysicalFunction {
    StorageMean(Arc<crate::StorageCurve>),
    Production(Arc<Unit>),
    Transfer(Arc<crate::river_transport::ReleaseFraction>),
    Table(Arc<TableCurve>),
    LevelFlow(Arc<crate::LevelFlowCurve>),
    OutletTable(Arc<TableCurve>),
    Characteristic(Arc<crate::PumpCharacteristic>),
    MinimumFlow(Arc<TurbineTable>),
    MaximumFlow(Arc<TurbineTable>),
}
struct Derivatives {
    gradient: [f64; 3],
    hessian: [[f64; 3]; 3],
}
impl From<Derivatives2> for Derivatives {
    fn from(d: Derivatives2) -> Self {
        Self {
            gradient: [d.gradient[0], d.gradient[1], 0.],
            hessian: [
                [d.hessian[0][0], d.hessian[0][1], 0.],
                [d.hessian[1][0], d.hessian[1][1], 0.],
                [0.; 3],
            ],
        }
    }
}
impl PhysicalFunction {
    fn is_linear(&self) -> bool {
        matches!(self,Self::StorageMean(c) if matches!(&**c,crate::StorageCurve::Polynomial{curvature,..} if *curvature==0.))
    }
    fn value(&self, x: [f64; 3]) -> f64 {
        match self {
            Self::Transfer(c) => c.data(x[0]).value,
            Self::Production(u) => crate::machine::power_balance(u, x[0], x[1], x[2]),
            Self::StorageMean(c) => c.mean_data(x[0], x[1]).value,
            Self::Table(c) => c.extended_value(x[0]),
            Self::LevelFlow(c) => c.extended_value(x[0], x[1]),
            Self::OutletTable(c) => c.extended_value(x[0]) * x[1],
            Self::Characteristic(c) => c.trial_discharge(x[1]),
            Self::MinimumFlow(c) => c.envelope_data(x[1], false).value,
            Self::MaximumFlow(c) => c.envelope_data(x[1], true).value,
        }
    }
    fn data(&self, x: [f64; 3]) -> Derivatives {
        match self {
            Self::Production(u) => {
                let (_, g, h) = hessian(
                    |x| crate::machine::power_balance(u, x[0], x[1], x[2]),
                    &SVector::from(x),
                );
                Derivatives {
                    gradient: std::array::from_fn(|i| g[i]),
                    hessian: std::array::from_fn(|i| std::array::from_fn(|j| h[(i, j)])),
                }
            }
            Self::Transfer(c) => c.data(x[0]).into(),
            Self::StorageMean(c) => c.mean_data(x[0], x[1]).into(),
            Self::Characteristic(c) => c.discharge_data(x[1]).into(),
            Self::MinimumFlow(c) => c.envelope_data(x[1], false).into(),
            Self::MaximumFlow(c) => c.envelope_data(x[1], true).into(),
            Self::LevelFlow(c) => c.data(x[0], x[1]).into(),
            Self::OutletTable(c) => {
                let (_, g, h) = hessian(
                    |x| crate::machine::table_value(c, x[0]) * x[1],
                    &SVector::from(x),
                );
                Derivatives {
                    gradient: std::array::from_fn(|i| g[i]),
                    hessian: std::array::from_fn(|i| std::array::from_fn(|j| h[(i, j)])),
                }
            }
            Self::Table(c) => {
                let coordinates = c.coordinates();
                let i = coordinates
                    .partition_point(|v| *v <= x[0])
                    .saturating_sub(1)
                    .min(coordinates.len() - 2);
                let slope =
                    (c.values()[i + 1] - c.values()[i]) / (coordinates[i + 1] - coordinates[i]);
                Derivatives {
                    gradient: [slope, 0., 0.],
                    hessian: [[0.; 3]; 3],
                }
            }
        }
    }
}
pub(crate) struct Term {
    pub function: PhysicalFunction,
    pub inputs: [LinearForm; 3],
    pub weight: f64,
}
/// lower <= linear(x) + Σ weight f(a(x),b(x)) <= upper.
pub(crate) struct Equation {
    pub linear: LinearForm,
    pub terms: Vec<Term>,
    pub bounds: Bounds,
}
impl Equation {
    /// Electrical work balance and the machine's additional flow restrictions.
    pub fn production(
        unit: Arc<Unit>,
        flow: LinearForm,
        head: LinearForm,
        power: LinearForm,
        flow_scale: f64,
        power_scale: f64,
    ) -> Vec<Self> {
        let mut equations = vec![Self {
            linear: LinearForm::default(),
            terms: vec![Term {
                function: PhysicalFunction::Production(unit.clone()),
                inputs: [flow.clone(), head.clone(), power],
                weight: 1. / power_scale,
            }],
            bounds: Bounds::new(0., 0.),
        }];
        let mut flow_equation = |function, bounds| {
            equations.push(Self {
                linear: flow.clone().scaled(-1. / flow_scale),
                terms: vec![Term {
                    function,
                    inputs: [flow.clone(), head.clone(), LinearForm::default()],
                    weight: 1. / flow_scale,
                }],
                bounds,
            })
        };
        if let Some(c) = &unit.pump_characteristic {
            flow_equation(
                PhysicalFunction::Characteristic(Arc::new(c.clone())),
                Bounds::new(0., 0.),
            );
        } else if let crate::Efficiency::Table(c) = &unit.efficiency {
            let c = Arc::new(c.clone());
            flow_equation(
                PhysicalFunction::MinimumFlow(c.clone()),
                Bounds::new(f64::NEG_INFINITY, 0.),
            );
            flow_equation(
                PhysicalFunction::MaximumFlow(c),
                Bounds::new(0., f64::INFINITY),
            );
        }
        equations
    }
}

impl From<StorageHeadEquation> for Equation {
    fn from(law: StorageHeadEquation) -> Self {
        Self {
            linear: LinearForm {
                constant: 0.,
                coefficients: BTreeMap::from([(law.head, -1.)]),
            },
            terms: law
                .vertices
                .windows(2)
                .map(|w| Term {
                    function: PhysicalFunction::StorageMean(law.curve.clone()),
                    inputs: [w[0].1.clone(), w[1].1.clone(), LinearForm::default()],
                    weight: (w[1].0 - w[0].0) / law.head_scale,
                })
                .collect(),
            bounds: Bounds::new(0., 0.),
        }
    }
}
// Contiguous sparse coefficients for repeated local evaluations. LinearForm's
// ordered maps are used only while assembling/coalescing the equations.
struct LinearFunction {
    constant: f64,
    coefficients: Vec<(usize, f64)>,
}
impl From<LinearForm> for LinearFunction {
    fn from(form: LinearForm) -> Self {
        Self {
            constant: form.constant,
            coefficients: form.coefficients.into_iter().collect(),
        }
    }
}
impl LinearFunction {
    fn value(&self, x: &[f64]) -> f64 {
        self.constant
            + self
                .coefficients
                .iter()
                .map(|&(i, a)| a * x[i])
                .sum::<f64>()
    }
}
struct PreparedTerm {
    function: PhysicalFunction,
    inputs: [LinearFunction; 3],
    weight: f64,
    jacobian: [Vec<(usize, f64)>; 3],
    /// Sparse position and coefficients of H_aa, H_ab and H_bb.
    hessian: Vec<(usize, [f64; 6])>,
}
struct PreparedEquation {
    linear: LinearFunction,
    bounds: Bounds,
    jacobian_columns: Vec<usize>,
    linear_jacobian: Vec<(usize, f64)>,
    terms: Vec<PreparedTerm>,
}
pub(crate) struct Nlp {
    base: NlTnlp,
    base_info: NlpInfo,
    equations: Vec<PreparedEquation>,
    jacobian_size: usize,
    hessian_coordinates: Vec<(i32, i32)>,
    base_hessian_positions: Vec<usize>,
    base_hessian_values: Vec<f64>,
    #[cfg(test)]
    final_x: Option<Vec<f64>>,
}
impl Nlp {
    pub fn new(mut base: NlTnlp, equations: Vec<Equation>) -> Result<Self, Error> {
        let info = base
            .get_nlp_info()
            .ok_or_else(|| Error::Model("invalid nonlinear dimensions".into()))?;
        if info.index_style != IndexStyle::C {
            return Err(Error::Model("zero-based sparse indices required".into()));
        }
        let nh = info.nnz_h_lag as usize;
        let mut rows = vec![0; nh];
        let mut cols = vec![0; nh];
        if nh > 0
            && !base.eval_h(
                None,
                false,
                1.,
                None,
                false,
                SparsityRequest::Structure {
                    irow: &mut rows,
                    jcol: &mut cols,
                },
            )
        {
            return Err(Error::Model("invalid Hessian structure".into()));
        }
        let mut coordinates = BTreeMap::new();
        for (&r, &c) in rows.iter().zip(&cols) {
            coordinates.insert((r.max(c), r.min(c)), 0);
        }
        for equation in &equations {
            if equation.bounds.min.is_nan()
                || equation.bounds.max.is_nan()
                || equation.bounds.min > equation.bounds.max
            {
                return Err(Error::Model("invalid constitutive bounds".into()));
            }
            let forms = std::iter::once(&equation.linear)
                .chain(equation.terms.iter().flat_map(|t| t.inputs.iter()));
            for form in forms {
                if !form.constant.is_finite()
                    || form
                        .coefficients
                        .iter()
                        .any(|(&i, &a)| i >= info.n as usize || !a.is_finite())
                {
                    return Err(Error::Model("invalid constitutive affine input".into()));
                }
            }
            for term in &equation.terms {
                if !term.weight.is_finite() {
                    return Err(Error::Model("invalid constitutive weight".into()));
                }
                if term.function.is_linear() {
                    continue;
                }
                let variables: BTreeMap<_, _> = term
                    .inputs
                    .iter()
                    .flat_map(|v| v.coefficients.keys())
                    .map(|&i| (i, ()))
                    .collect();
                let vars: Vec<_> = variables.keys().copied().collect();
                for (j, &a) in vars.iter().enumerate() {
                    for &b in &vars[..=j] {
                        coordinates.insert((a as i32, b as i32), 0);
                    }
                }
            }
        }
        for (i, position) in coordinates.values_mut().enumerate() {
            *position = i;
        }
        let base_hessian_positions = rows
            .iter()
            .zip(&cols)
            .map(|(&r, &c)| coordinates[&(r.max(c), r.min(c))])
            .collect();
        let mut compiled = Vec::new();
        let mut jacobian_size = info.nnz_jac_g as usize;
        for equation in equations {
            let mut columns: BTreeMap<_, usize> = std::iter::once(&equation.linear)
                .chain(equation.terms.iter().flat_map(|t| t.inputs.iter()))
                .flat_map(|v| v.coefficients.keys())
                .map(|&i| (i, 0))
                .collect();
            for (k, pos) in columns.values_mut().enumerate() {
                *pos = k;
            }
            let mut terms = Vec::new();
            for term in equation.terms {
                let variables: BTreeMap<_, _> = term
                    .inputs
                    .iter()
                    .flat_map(|v| v.coefficients.keys())
                    .map(|&i| (i, ()))
                    .collect();
                let vars: Vec<_> = variables.keys().copied().collect();
                let mut hessian = Vec::new();
                if !term.function.is_linear() {
                    for (j, &a) in vars.iter().enumerate() {
                        for &b in &vars[..=j] {
                            let ja = term
                                .inputs
                                .each_ref()
                                .map(|v| v.coefficients.get(&a).copied().unwrap_or(0.));
                            let jb = term
                                .inputs
                                .each_ref()
                                .map(|v| v.coefficients.get(&b).copied().unwrap_or(0.));
                            let mut weights = [0.; 6];
                            let mut k = 0;
                            for u in 0..3 {
                                for v in 0..=u {
                                    weights[k] = term.weight
                                        * if u == v {
                                            ja[u] * jb[u]
                                        } else {
                                            ja[u] * jb[v] + ja[v] * jb[u]
                                        };
                                    k += 1;
                                }
                            }
                            hessian.push((coordinates[&(a as i32, b as i32)], weights));
                        }
                    }
                }
                terms.push(PreparedTerm {
                    jacobian: term.inputs.each_ref().map(|v| {
                        v.coefficients
                            .iter()
                            .map(|(&i, &a)| (columns[&i], term.weight * a))
                            .collect()
                    }),
                    function: term.function,
                    inputs: term.inputs.map(LinearFunction::from),
                    weight: term.weight,
                    hessian,
                });
            }
            jacobian_size += columns.len();
            compiled.push(PreparedEquation {
                linear_jacobian: equation
                    .linear
                    .coefficients
                    .iter()
                    .map(|(&i, &a)| (columns[&i], a))
                    .collect(),
                linear: equation.linear.into(),
                bounds: equation.bounds,
                jacobian_columns: columns.keys().copied().collect(),
                terms,
            });
        }
        Ok(Self {
            base,
            base_info: info,
            equations: compiled,
            jacobian_size,
            hessian_coordinates: coordinates.keys().copied().collect(),
            base_hessian_positions,
            base_hessian_values: vec![0.; nh],
            #[cfg(test)]
            final_x: None,
        })
    }
    #[cfg(test)]
    pub fn final_x(&self) -> Option<&[f64]> {
        self.final_x.as_deref()
    }
}
impl TNLP for Nlp {
    fn get_nlp_info(&mut self) -> Option<NlpInfo> {
        Some(NlpInfo {
            m: self.base_info.m + self.equations.len() as i32,
            nnz_jac_g: self.jacobian_size as i32,
            nnz_h_lag: self.hessian_coordinates.len() as i32,
            ..self.base_info
        })
    }
    fn get_bounds_info(&mut self, b: BoundsInfo<'_>) -> bool {
        let m = self.base_info.m as usize;
        for (k, e) in self.equations.iter().enumerate() {
            b.g_l[m + k] = e.bounds.min;
            b.g_u[m + k] = e.bounds.max;
        }
        self.base.get_bounds_info(BoundsInfo {
            x_l: b.x_l,
            x_u: b.x_u,
            g_l: &mut b.g_l[..m],
            g_u: &mut b.g_u[..m],
        })
    }
    fn get_starting_point(&mut self, s: StartingPoint<'_>) -> bool {
        let m = self.base_info.m as usize;
        if s.init_lambda {
            s.lambda[m..].fill(0.);
        }
        self.base.get_starting_point(StartingPoint {
            lambda: &mut s.lambda[..m],
            ..s
        })
    }
    fn eval_f(&mut self, x: &[f64], new_x: bool) -> Option<f64> {
        self.base.eval_f(x, new_x)
    }
    fn eval_grad_f(&mut self, x: &[f64], new_x: bool, gradient: &mut [f64]) -> bool {
        self.base.eval_grad_f(x, new_x, gradient)
    }
    fn eval_g(&mut self, x: &[f64], new_x: bool, g: &mut [f64]) -> bool {
        let m = self.base_info.m as usize;
        if !self.base.eval_g(x, new_x, &mut g[..m]) {
            return false;
        }
        for (k, e) in self.equations.iter().enumerate() {
            g[m + k] = e.linear.value(x)
                + e.terms
                    .iter()
                    .map(|t| t.weight * t.function.value(t.inputs.each_ref().map(|v| v.value(x))))
                    .sum::<f64>();
        }
        g.iter().all(|v| v.is_finite())
    }
    fn eval_jac_g(&mut self, x: Option<&[f64]>, new_x: bool, mode: SparsityRequest<'_>) -> bool {
        let mut nz = self.base_info.nnz_jac_g as usize;
        match mode {
            SparsityRequest::Structure { irow, jcol } => {
                if !self.base.eval_jac_g(
                    None,
                    new_x,
                    SparsityRequest::Structure {
                        irow: &mut irow[..nz],
                        jcol: &mut jcol[..nz],
                    },
                ) {
                    return false;
                }
                for (k, e) in self.equations.iter().enumerate() {
                    for &col in &e.jacobian_columns {
                        irow[nz] = self.base_info.m + k as i32;
                        jcol[nz] = col as i32;
                        nz += 1;
                    }
                }
            }
            SparsityRequest::Values { values } => {
                if !self.base.eval_jac_g(
                    x,
                    new_x,
                    SparsityRequest::Values {
                        values: &mut values[..nz],
                    },
                ) {
                    return false;
                }
                let Some(x) = x else {
                    return false;
                };
                for e in &self.equations {
                    let out = &mut values[nz..nz + e.jacobian_columns.len()];
                    out.fill(0.);
                    for &(i, a) in &e.linear_jacobian {
                        out[i] += a;
                    }
                    for t in &e.terms {
                        let g = t
                            .function
                            .data(t.inputs.each_ref().map(|v| v.value(x)))
                            .gradient;
                        for (input, &derivative) in t.jacobian.iter().zip(&g) {
                            for &(i, a) in input {
                                out[i] += a * derivative;
                            }
                        }
                    }
                    nz += out.len();
                }
                return values.iter().all(|v| v.is_finite());
            }
        }
        true
    }
    fn eval_h(
        &mut self,
        x: Option<&[f64]>,
        new_x: bool,
        factor: f64,
        lambda: Option<&[f64]>,
        new_lambda: bool,
        mode: SparsityRequest<'_>,
    ) -> bool {
        match mode {
            SparsityRequest::Structure { irow, jcol } => {
                for (k, &(r, c)) in self.hessian_coordinates.iter().enumerate() {
                    irow[k] = r;
                    jcol[k] = c;
                }
            }
            SparsityRequest::Values { values } => {
                let (Some(x), Some(lambda)) = (x, lambda) else {
                    return false;
                };
                if !self.base_hessian_values.is_empty()
                    && !self.base.eval_h(
                        Some(x),
                        new_x,
                        factor,
                        Some(&lambda[..self.base_info.m as usize]),
                        new_lambda,
                        SparsityRequest::Values {
                            values: &mut self.base_hessian_values,
                        },
                    )
                {
                    return false;
                }
                values.fill(0.);
                for (&k, &v) in self
                    .base_hessian_positions
                    .iter()
                    .zip(&self.base_hessian_values)
                {
                    values[k] += v;
                }
                for (k, e) in self.equations.iter().enumerate() {
                    for t in &e.terms {
                        if t.hessian.is_empty() {
                            continue;
                        }
                        let h = t
                            .function
                            .data(t.inputs.each_ref().map(|v| v.value(x)))
                            .hessian;
                        let multiplier = lambda[self.base_info.m as usize + k];
                        for &(position, a) in &t.hessian {
                            values[position] += multiplier
                                * (a[0] * h[0][0]
                                    + a[1] * h[1][0]
                                    + a[2] * h[1][1]
                                    + a[3] * h[2][0]
                                    + a[4] * h[2][1]
                                    + a[5] * h[2][2]);
                        }
                    }
                }
                return values.iter().all(|v| v.is_finite());
            }
        }
        true
    }
    fn get_constraints_linearity(&mut self, types: &mut [pounce_nlp::tnlp::Linearity]) -> bool {
        let m = self.base_info.m as usize;
        if !self.base.get_constraints_linearity(&mut types[..m]) {
            return false;
        }
        for (kind, equation) in types[m..].iter_mut().zip(&self.equations) {
            *kind = if equation.terms.iter().all(|t| t.function.is_linear()) {
                pounce_nlp::tnlp::Linearity::Linear
            } else {
                pounce_nlp::tnlp::Linearity::NonLinear
            };
        }
        true
    }
    fn get_objective_variables_linearity(
        &mut self,
        types: &mut [pounce_nlp::tnlp::Linearity],
    ) -> bool {
        self.base.get_objective_variables_linearity(types)
    }
    fn derivative_proofs(&mut self) -> pounce_nlp::constant_derivatives::DerivativeProofs {
        if self.equations.is_empty() {
            self.base.derivative_proofs()
        } else {
            Default::default()
        }
    }
    fn finalize_solution(&mut self, s: Solution<'_>, data: &IpoptData, cq: &IpoptCq) {
        #[cfg(test)]
        {
            self.final_x = Some(s.x.to_vec());
        }
        let m = self.base_info.m as usize;
        self.base.finalize_solution(
            Solution {
                g: &s.g[..m],
                lambda: &s.lambda[..m],
                ..s
            },
            data,
            cq,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oximo::prelude::*;

    fn lagrangian_gradient(
        problem: &mut dyn TNLP,
        x: &[f64],
        factor: f64,
        lambda: &[f64],
    ) -> Vec<f64> {
        let info = problem.get_nlp_info().unwrap();
        let mut gradient = vec![0.; x.len()];
        assert!(problem.eval_grad_f(x, true, &mut gradient));
        for g in &mut gradient {
            *g *= factor;
        }
        let mut rows = vec![0; info.nnz_jac_g as usize];
        let mut cols = rows.clone();
        let mut values = vec![0.; rows.len()];
        assert!(problem.eval_jac_g(
            None,
            false,
            SparsityRequest::Structure {
                irow: &mut rows,
                jcol: &mut cols
            }
        ));
        assert!(problem.eval_jac_g(
            Some(x),
            true,
            SparsityRequest::Values {
                values: &mut values
            }
        ));
        for k in 0..values.len() {
            gradient[cols[k] as usize] += lambda[rows[k] as usize] * values[k];
        }
        gradient
    }
    #[test]
    fn affine_paths_and_algebraic_terms_share_correct_sparse_second_derivatives() {
        let m = Model::new("constitutive chain rule");
        variable!(m,0.0<=a<=4.0,initial=0.2);
        variable!(m,0.0<=b<=4.0,initial=2.4);
        variable!(m,0.0<=head<=30.0,initial=14.0);
        constraint!(m, a.square() + a * b <= 5.0);
        objective!(m, Min, a.square() + b * head + head.square());
        let ai = a.var_id().unwrap().index();
        let bi = b.var_id().unwrap().index();
        let law = StorageHeadEquation {
            curve: Arc::new(crate::StorageCurve::Table(
                TableCurve::new(vec![0., 1., 3., 4.], vec![10., 12., 20., 21.]).unwrap(),
            )),
            vertices: vec![
                (
                    0.,
                    LinearForm {
                        constant: 0.,
                        coefficients: BTreeMap::from([(ai, 1.)]),
                    },
                ),
                (
                    0.4,
                    LinearForm {
                        constant: 0.2,
                        coefficients: BTreeMap::from([(ai, 0.3), (bi, 0.7)]),
                    },
                ),
                (
                    1.,
                    LinearForm {
                        constant: 0.,
                        coefficients: BTreeMap::from([(bi, 1.)]),
                    },
                ),
            ],
            head: head.var_id().unwrap().index(),
            head_scale: 1.,
        };
        let mut problem = crate::nlp::prepare(&m, vec![law.into()]).unwrap();
        check_derivatives(
            &mut problem,
            &[[0.2, 2.4, 14.], [2.6, 0.4, 14.], [1.2, 2.5, 14.]],
            &[0.6, 1.3],
        );
    }
    #[test]
    fn turbine_rows_preserve_sparse_second_derivatives_with_coupled_inputs() {
        let m = Model::new("turbine chain rule");
        variable!(m,0.0<=a<=4.0,initial=0.2);
        variable!(m,0.0<=b<=4.0,initial=2.4);
        variable!(m,0.0<=power<=30.0,initial=14.0);
        constraint!(m, a.square() + a * b <= 5.0);
        objective!(m, Min, a.square() + b * power + power.square());
        let curve = Arc::new(
            TurbineTable::new(
                vec![80., 120., 160.],
                vec![
                    TableCurve::new(vec![0., 6., 11., 20.], vec![0.6, 0.85, 0.9, 0.82]).unwrap(),
                    TableCurve::new(vec![0., 7., 13., 20.], vec![0.65, 0.88, 0.92, 0.85]).unwrap(),
                    TableCurve::new(vec![0., 4., 10., 20.], vec![0.6, 0.86, 0.94, 0.87]).unwrap(),
                ],
                vec![
                    Bounds::new(2., 12.),
                    Bounds::new(3., 15.),
                    Bounds::new(4., 20.),
                ],
                crate::TurbineOptions {
                    interpolation: crate::DischargeInterpolation::Pchip,
                    ..Default::default()
                },
            )
            .unwrap(),
        );
        let flow = LinearForm {
            constant: 2.,
            coefficients: BTreeMap::from([
                (a.var_id().unwrap().index(), 3.),
                (b.var_id().unwrap().index(), 1.),
            ]),
        };
        let head = LinearForm {
            constant: 80.,
            coefficients: BTreeMap::from([
                (a.var_id().unwrap().index(), 10.),
                (b.var_id().unwrap().index(), 3.),
            ]),
        };
        let mut unit = crate::Unit::new(
            "Turbine",
            "Station",
            Bounds::new(0., 20.),
            Bounds::new(0., 30.),
            Bounds::new(80., 160.),
            crate::Efficiency::Table((*curve).clone()),
        );
        unit.generator_efficiency_curve =
            Some(TableCurve::new(vec![0., 15., 30.], vec![0.88, 0.95, 0.97]).unwrap());
        let equations = Equation::production(
            Arc::new(unit),
            flow,
            head,
            LinearForm::variable(power.var_id().unwrap().index(), 10.),
            20.,
            10.,
        );
        let mut problem = crate::nlp::prepare(&m, equations).unwrap();
        check_derivatives(
            &mut problem,
            &[[0.2, 2.4, 14.], [2.6, 0.4, 14.], [1.2, 2.5, 14.]],
            &[0.6, 1.3, -0.2, 0.4],
        );
    }
    #[test]
    fn river_law_jacobians_and_lagrangian_hessians_include_shared_inputs() {
        let m = Model::new("river chain rule");
        variable!(m,0.<=a<=4.,initial=0.2);
        variable!(m,0.<=b<=4.,initial=2.4);
        variable!(m,0.<=z<=30.,initial=14.);
        constraint!(m, a.square() + a * b <= 5.);
        objective!(m, Min, a.square() + b * z + z.square());
        let flow = crate::nlp::linear_form(&m, (a * 3. + b + 2.).into()).unwrap();
        let reference = crate::nlp::linear_form(&m, (a * 10. + b * 3. + 80.).into()).unwrap();
        let curves = vec![
            TableCurve::new(vec![-20., 0., 20.], vec![-4., 0., 8.]).unwrap(),
            TableCurve::new(vec![-20., 0., 20.], vec![-8., 0., 16.]).unwrap(),
        ];
        let functions = [
            PhysicalFunction::Transfer(Arc::new(crate::river_transport::ReleaseFraction {
                references: vec![0., 10., 40.],
                coefficients: vec![0.8, 0.6, 0.2],
            })),
            PhysicalFunction::OutletTable(Arc::new(
                TableCurve::new(vec![0., 10., 30.], vec![1., 4., 8.]).unwrap(),
            )),
            PhysicalFunction::LevelFlow(Arc::new(
                crate::LevelFlowCurve::new(vec![80., 160.], curves).unwrap(),
            )),
        ];
        let equations = functions
            .into_iter()
            .map(|function| Equation {
                linear: LinearForm::variable(z.var_id().unwrap().index(), -2.),
                terms: vec![Term {
                    function,
                    inputs: [flow.clone(), reference.clone(), Default::default()],
                    weight: 0.3,
                }],
                bounds: Bounds::new(0., 0.),
            })
            .collect();
        let mut problem = crate::nlp::prepare(&m, equations).unwrap();
        check_derivatives(
            &mut problem,
            &[[0.2, 2.4, 14.], [2.6, 0.4, 14.], [1.2, 2.5, 14.]],
            &[0.6, 1.3, -0.2, 0.4],
        );
    }
    fn check_derivatives(problem: &mut dyn TNLP, points: &[[f64; 3]], lambda: &[f64]) {
        let info = problem.get_nlp_info().unwrap();
        let mut hr = vec![0; info.nnz_h_lag as usize];
        let mut hc = hr.clone();
        assert!(problem.eval_h(
            None,
            false,
            0.7,
            None,
            false,
            SparsityRequest::Structure {
                irow: &mut hr,
                jcol: &mut hc
            }
        ));
        for &point in points {
            let mut values = vec![0.; hr.len()];
            assert!(problem.eval_h(
                Some(&point),
                true,
                0.7,
                Some(lambda),
                true,
                SparsityRequest::Values {
                    values: &mut values
                }
            ));
            let mut dense = [[0.; 3]; 3];
            for (k, &v) in values.iter().enumerate() {
                dense[hr[k] as usize][hc[k] as usize] = v;
                dense[hc[k] as usize][hr[k] as usize] = v;
            }
            for column in 0..3 {
                let mut plus = point;
                plus[column] += 1e-5;
                let mut minus = point;
                minus[column] -= 1e-5;
                let mut gp = vec![0.; info.m as usize];
                let mut gn = gp.clone();
                assert!(problem.eval_g(&plus, true, &mut gp));
                assert!(problem.eval_g(&minus, true, &mut gn));
                let mut jr = vec![0; info.nnz_jac_g as usize];
                let mut jc = jr.clone();
                let mut jv = vec![0.; jr.len()];
                assert!(problem.eval_jac_g(
                    None,
                    false,
                    SparsityRequest::Structure {
                        irow: &mut jr,
                        jcol: &mut jc
                    }
                ));
                assert!(problem.eval_jac_g(
                    Some(&point),
                    true,
                    SparsityRequest::Values { values: &mut jv }
                ));
                for row in 0..info.m as usize {
                    let expected: f64 = (0..jv.len())
                        .filter(|&k| jr[k] as usize == row && jc[k] as usize == column)
                        .map(|k| jv[k])
                        .sum();
                    assert!(((gp[row] - gn[row]) / 2e-5 - expected).abs() < 1e-7);
                }
                let p = lagrangian_gradient(problem, &plus, 0.7, lambda);
                let n = lagrangian_gradient(problem, &minus, 0.7, lambda);
                for row in 0..3 {
                    assert!(
                        ((p[row] - n[row]) / 2e-5 - dense[row][column]).abs() < 1e-7,
                        "{point:?}, ({row},{column})"
                    );
                }
            }
        }
    }
}
