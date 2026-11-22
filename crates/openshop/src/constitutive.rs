//! Sparse constitutive equations appended to an algebraic nonlinear model.
//! Affine storage paths are composed with local integral derivatives. Sparse
//! matrix coordinates and chain-rule coefficients are prepared once; callbacks
//! only evaluate curves and accumulate into existing numeric buffers.

use crate::{Error, TableCurve, river::LinearForm};
use pounce_nl::nl_reader::NlTnlp;
use pounce_nlp::tnlp::{
    BoundsInfo, IndexStyle, IpoptCq, IpoptData, NlpInfo, Solution, SparsityRequest, StartingPoint,
    TNLP,
};
use std::{collections::BTreeMap, sync::Arc};

pub(crate) struct StorageHeadEquation {
    pub curve: Arc<TableCurve>,
    /// Fractional time and affine physical volume (m³), including both endpoints.
    pub vertices: Vec<(f64, LinearForm)>,
    pub head: usize,
    pub head_scale: f64,
}
struct Segment {
    before: LinearForm,
    after: LinearForm,
    weight: f64,
    before_jacobian: Vec<(usize, f64)>,
    after_jacobian: Vec<(usize, f64)>,
    /// Sparse position and coefficients of H_aa, H_ab and H_bb.
    hessian: Vec<(usize, [f64; 3])>,
}
struct StorageHead {
    curve: Arc<TableCurve>,
    head: usize,
    jacobian_columns: Vec<usize>,
    head_position: usize,
    segments: Vec<Segment>,
}

pub(crate) struct Nlp {
    base: NlTnlp,
    base_info: NlpInfo,
    laws: Vec<StorageHead>,
    jacobian_size: usize,
    hessian_coordinates: Vec<(i32, i32)>,
    base_hessian_positions: Vec<usize>,
    base_hessian_values: Vec<f64>,
    #[cfg(test)]
    final_x: Option<Vec<f64>>,
}
impl Nlp {
    pub fn new(mut base: NlTnlp, laws: Vec<StorageHeadEquation>) -> Result<Self, Error> {
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
        for law in &laws {
            if law.head >= info.n as usize
                || law.vertices.len() < 2
                || law.head_scale <= 0.
                || !law.head_scale.is_finite()
                || law.vertices.first().unwrap().0 != 0.
                || law.vertices.last().unwrap().0 != 1.
                || !law.vertices.windows(2).all(|w| w[0].0 < w[1].0)
            {
                return Err(Error::Model("invalid storage-head integration path".into()));
            }
            for pair in law.vertices.windows(2) {
                let mut variables: BTreeMap<_, _> = pair[0]
                    .1
                    .coefficients
                    .keys()
                    .chain(pair[1].1.coefficients.keys())
                    .map(|&i| (i, ()))
                    .collect();
                variables.remove(&law.head);
                let vars: Vec<_> = variables.keys().copied().collect();
                for (j, &a) in vars.iter().enumerate() {
                    if a >= info.n as usize {
                        return Err(Error::Model("invalid storage path variable".into()));
                    }
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
        for law in laws {
            let mut columns: BTreeMap<_, usize> = law
                .vertices
                .iter()
                .flat_map(|(_, v)| v.coefficients.keys().copied())
                .chain([law.head])
                .map(|i| (i, 0))
                .collect();
            for (k, pos) in columns.values_mut().enumerate() {
                *pos = k;
            }
            let mut segments = Vec::new();
            for pair in law.vertices.windows(2) {
                let before = pair[0].1.clone();
                let after = pair[1].1.clone();
                if before.coefficients.contains_key(&law.head)
                    || after.coefficients.contains_key(&law.head)
                {
                    return Err(Error::Model(
                        "storage path cannot depend on its own mean head".into(),
                    ));
                }
                let weight = (pair[1].0 - pair[0].0) / law.head_scale;
                let mut variables: BTreeMap<_, _> = before
                    .coefficients
                    .keys()
                    .chain(after.coefficients.keys())
                    .map(|&i| (i, ()))
                    .collect();
                variables.remove(&law.head);
                let vars: Vec<_> = variables.keys().copied().collect();
                let mut hessian = Vec::new();
                for (j, &a) in vars.iter().enumerate() {
                    for &b in &vars[..=j] {
                        let aj = before.coefficients.get(&a).copied().unwrap_or(0.);
                        let ak = before.coefficients.get(&b).copied().unwrap_or(0.);
                        let bj = after.coefficients.get(&a).copied().unwrap_or(0.);
                        let bk = after.coefficients.get(&b).copied().unwrap_or(0.);
                        hessian.push((
                            coordinates[&(a as i32, b as i32)],
                            [
                                weight * aj * ak,
                                weight * (aj * bk + bj * ak),
                                weight * bj * bk,
                            ],
                        ));
                    }
                }
                segments.push(Segment {
                    before_jacobian: before
                        .coefficients
                        .iter()
                        .map(|(&i, &v)| (columns[&i], weight * v))
                        .collect(),
                    after_jacobian: after
                        .coefficients
                        .iter()
                        .map(|(&i, &v)| (columns[&i], weight * v))
                        .collect(),
                    before,
                    after,
                    weight,
                    hessian,
                });
            }
            jacobian_size += columns.len();
            compiled.push(StorageHead {
                curve: law.curve,
                head: law.head,
                head_position: columns[&law.head],
                jacobian_columns: columns.keys().copied().collect(),
                segments,
            });
        }
        Ok(Self {
            base,
            base_info: info,
            laws: compiled,
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
            m: self.base_info.m + self.laws.len() as i32,
            nnz_jac_g: self.jacobian_size as i32,
            nnz_h_lag: self.hessian_coordinates.len() as i32,
            ..self.base_info
        })
    }
    fn get_bounds_info(&mut self, b: BoundsInfo<'_>) -> bool {
        let m = self.base_info.m as usize;
        b.g_l[m..].fill(0.);
        b.g_u[m..].fill(0.);
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
        for (k, law) in self.laws.iter().enumerate() {
            g[m + k] = law
                .segments
                .iter()
                .map(|s| {
                    s.weight
                        * law
                            .curve
                            .mean_data(s.before.value(x), s.after.value(x))
                            .value
                })
                .sum::<f64>()
                - x[law.head];
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
                for (k, law) in self.laws.iter().enumerate() {
                    for &col in &law.jacobian_columns {
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
                for law in &self.laws {
                    let out = &mut values[nz..nz + law.jacobian_columns.len()];
                    out.fill(0.);
                    out[law.head_position] = -1.;
                    for s in &law.segments {
                        let g = law
                            .curve
                            .mean_data(s.before.value(x), s.after.value(x))
                            .gradient;
                        for &(i, a) in &s.before_jacobian {
                            out[i] += a * g[0];
                        }
                        for &(i, a) in &s.after_jacobian {
                            out[i] += a * g[1];
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
                for (k, law) in self.laws.iter().enumerate() {
                    for s in &law.segments {
                        let h = law
                            .curve
                            .mean_data(s.before.value(x), s.after.value(x))
                            .hessian;
                        let multiplier = lambda[self.base_info.m as usize + k];
                        for &(position, a) in &s.hessian {
                            values[position] +=
                                multiplier * (a[0] * h[0][0] + a[1] * h[1][0] + a[2] * h[1][1]);
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
        types[m..].fill(pounce_nlp::tnlp::Linearity::NonLinear);
        true
    }
    fn get_objective_variables_linearity(
        &mut self,
        types: &mut [pounce_nlp::tnlp::Linearity],
    ) -> bool {
        self.base.get_objective_variables_linearity(types)
    }
    fn derivative_proofs(&mut self) -> pounce_nlp::constant_derivatives::DerivativeProofs {
        if self.laws.is_empty() {
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
            curve: Arc::new(
                TableCurve::new(vec![0., 1., 3., 4.], vec![10., 12., 20., 21.]).unwrap(),
            ),
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
        let mut problem = crate::nlp::prepare(&m, vec![law]).unwrap();
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
        for point in [[0.2, 2.4, 14.], [2.6, 0.4, 14.], [1.2, 2.5, 14.]] {
            let lambda = [0.6, 1.3];
            let mut values = vec![0.; hr.len()];
            assert!(problem.eval_h(
                Some(&point),
                true,
                0.7,
                Some(&lambda),
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
                let p = lagrangian_gradient(&mut problem, &plus, 0.7, &lambda);
                let n = lagrangian_gradient(&mut problem, &minus, 0.7, &lambda);
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
