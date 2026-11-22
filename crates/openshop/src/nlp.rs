//! Translate the algebraic model once, then reuse exact sparse derivatives.
//!
//! The public expression interfaces avoid a text/file round trip and preserve
//! variable ordering. Linear conservation rows bypass automatic differentiation.

use std::sync::Arc;
#[cfg(test)]
use std::{cell::RefCell, rc::Rc};

use crate::constitutive::{Nlp, StorageHeadEquation};
use oximo::{
    Model,
    core::{Domain, ObjectiveSense},
    expr::{ExprArena, ExprId, ExprNode, UnaryOp as OximoUnary, extract_linear},
};
use pounce_nl::nl_reader::{BinOp, Expr, NlBody, NlProblem, NlProblemParts, NlTnlp, UnaryOp};
#[cfg(test)]
use pounce_nlp::tnlp::TNLP;

use crate::{Error, SolveOptions};

#[derive(Clone, Copy, Debug)]
pub(crate) enum Backend {
    Ipopt,
    #[cfg(test)]
    Pounce,
}

pub(crate) struct NlpSolution {
    pub x: Vec<f64>,
    pub status: i32,
    pub iterations: usize,
}

/// Convert a prepared sparse linear form at the equation-assembly boundary.
pub(crate) fn linear_expression<'a>(
    form: &crate::river::LinearForm,
    columns: &[oximo::expr::Expr<'a>],
    zero: oximo::expr::Expr<'a>,
) -> oximo::expr::Expr<'a> {
    form.coefficients
        .iter()
        .fold(zero + form.constant, |sum, (&i, &a)| sum + columns[i] * a)
}

pub(crate) fn prepare(model: &Model, laws: Vec<StorageHeadEquation>) -> Result<Nlp, Error> {
    let variables = model.resolved_variables();
    if variables.iter().any(|v| v.domain != Domain::Real) {
        return Err(Error::Model(
            "nonlinear dispatch requires fixed integer decisions".into(),
        ));
    }
    let arena = model.arena();
    let mut expressions = vec![None; arena.len()];
    let objective = model
        .try_objective()
        .map_err(|e| Error::Model(e.to_string()))?;
    let rows = model.constraints();
    let rows: Vec<_> = rows.algebraic().iter().filter(|row| row.active).collect();
    let mut problem = NlProblem::from_expressions(NlProblemParts {
        minimize: objective.sense == ObjectiveSense::Minimize,
        objective: lower(&arena, objective.expr, &mut expressions)?,
        obj_constant: 0.0,
        constraints: rows
            .iter()
            .map(|row| lower(&arena, row.lhs, &mut expressions))
            .collect::<Result<_, _>>()?,
        x_l: variables.iter().map(|v| v.lb).collect(),
        x_u: variables.iter().map(|v| v.ub).collect(),
        x0: variables
            .iter()
            .map(|v| {
                v.initial.unwrap_or_else(|| {
                    if v.lb.is_finite() && v.ub.is_finite() {
                        0.5 * v.lb + 0.5 * v.ub
                    } else {
                        0.0_f64.max(v.lb).min(v.ub)
                    }
                })
            })
            .collect(),
        g_l: rows.iter().map(|row| row.lower).collect(),
        g_u: rows.iter().map(|row| row.upper).collect(),
        var_names: variables.iter().map(|v| v.name.to_string()).collect(),
        con_names: rows.iter().map(|row| row.name.to_string()).collect(),
    })
    .map_err(Error::Model)?;
    for (i, row) in rows.iter().enumerate() {
        if let Some(linear) = extract_linear(&arena, row.lhs) {
            problem.con_linear[i] = linear.coeffs.iter().map(|(v, c)| (v.index(), *c)).collect();
            problem.con_nonlinear[i] = NlBody::Tree(Expr::Const(0.0));
            problem.g_l[i] -= linear.constant;
            problem.g_u[i] -= linear.constant;
        }
    }
    if let Some(linear) = extract_linear(&arena, objective.expr) {
        problem.obj_linear = linear.coeffs.iter().map(|(v, c)| (v.index(), *c)).collect();
        problem.obj_nonlinear = NlBody::Tree(Expr::Const(0.0));
        problem.obj_constant = linear.constant;
    }
    Nlp::new(NlTnlp::try_new(problem).map_err(Error::Model)?, laws)
}

fn lower(arena: &ExprArena, id: ExprId, memo: &mut [Option<Arc<Expr>>]) -> Result<Expr, Error> {
    if let Some(expr) = &memo[id.index()] {
        return Ok(Expr::Cse(expr.clone()));
    }
    let binary = |op, a, b| Expr::Binary(op, Box::new(a), Box::new(b));
    let result = match arena.get(id) {
        ExprNode::Const(c) => Expr::Const(*c),
        ExprNode::Var(v) => Expr::Var(v.index()),
        ExprNode::Param(p) => Expr::Const(
            arena
                .try_param_value(*p)
                .ok_or_else(|| Error::Model("missing parameter".into()))?,
        ),
        ExprNode::Add(children) => Expr::Sum(
            children
                .iter()
                .map(|&c| lower(arena, c, memo))
                .collect::<Result<_, _>>()?,
        ),
        ExprNode::Mul(children) => {
            let mut terms = children.iter().map(|&c| lower(arena, c, memo));
            let mut product = terms.next().transpose()?.unwrap_or(Expr::Const(1.0));
            for term in terms {
                product = binary(BinOp::Mul, product, term?);
            }
            product
        }
        ExprNode::Div(a, b) => binary(BinOp::Div, lower(arena, *a, memo)?, lower(arena, *b, memo)?),
        ExprNode::Pow(a, b) => binary(BinOp::Pow, lower(arena, *a, memo)?, lower(arena, *b, memo)?),
        ExprNode::Atan2(a, b) => binary(
            BinOp::Atan2,
            lower(arena, *a, memo)?,
            lower(arena, *b, memo)?,
        ),
        ExprNode::Min(children) => Expr::MinList(
            children
                .iter()
                .map(|&c| lower(arena, c, memo))
                .collect::<Result<_, _>>()?,
        ),
        ExprNode::Max(children) => Expr::MaxList(
            children
                .iter()
                .map(|&c| lower(arena, c, memo))
                .collect::<Result<_, _>>()?,
        ),
        ExprNode::Linear { coeffs, constant } => Expr::Sum(
            std::iter::once(Expr::Const(*constant))
                .chain(
                    coeffs
                        .iter()
                        .map(|(v, c)| binary(BinOp::Mul, Expr::Const(*c), Expr::Var(v.index()))),
                )
                .collect(),
        ),
        ExprNode::Unary(op, x) => {
            let op = match op {
                OximoUnary::Neg => UnaryOp::Neg,
                OximoUnary::Abs => UnaryOp::Abs,
                OximoUnary::Sqrt => UnaryOp::Sqrt,
                OximoUnary::Exp => UnaryOp::Exp,
                OximoUnary::Log => UnaryOp::Log,
                OximoUnary::Log10 => UnaryOp::Log10,
                OximoUnary::Sin => UnaryOp::Sin,
                OximoUnary::Cos => UnaryOp::Cos,
                OximoUnary::Tan => UnaryOp::Tan,
                OximoUnary::Asin => UnaryOp::Asin,
                OximoUnary::Acos => UnaryOp::Acos,
                OximoUnary::Atan => UnaryOp::Atan,
                OximoUnary::Sinh => UnaryOp::Sinh,
                OximoUnary::Cosh => UnaryOp::Cosh,
                OximoUnary::Tanh => UnaryOp::Tanh,
                OximoUnary::Asinh => UnaryOp::Asinh,
                OximoUnary::Acosh => UnaryOp::Acosh,
                OximoUnary::Atanh => UnaryOp::Atanh,
                other => {
                    return Err(Error::Model(format!(
                        "unsupported derivative operation: {other:?}"
                    )));
                }
            };
            Expr::Unary(op, Box::new(lower(arena, *x, memo)?))
        }
    };
    let result = Arc::new(result);
    memo[id.index()] = Some(result.clone());
    Ok(Expr::Cse(result))
}

pub(crate) fn solve(
    mut problem: Nlp,
    options: &SolveOptions,
    backend: Backend,
) -> Result<NlpSolution, Error> {
    options.validate()?;
    match backend {
        Backend::Ipopt => crate::native_ipopt::solve(&mut problem, options),
        #[cfg(test)]
        Backend::Pounce => {
            let mut app = pounce_rs::IpoptApplication::new();
            let opts = app.options_mut();
            for (name, value) in [
                ("tol", options.tolerance),
                ("constr_viol_tol", options.tolerance),
                ("bound_relax_factor", 0.0),
                ("max_wall_time", options.time_limit),
            ] {
                opts.set_numeric_value(name, value, true, true)
                    .map_err(|e| Error::Solver(e.to_string()))?;
            }
            for (name, value) in [("print_level", 0), ("max_iter", options.max_iterations)] {
                opts.set_integer_value(name, value, true, true)
                    .map_err(|e| Error::Solver(e.to_string()))?;
            }
            // Match native Ipopt's algorithm and avoid implicit NLP reformulations.
            opts.set_string_value("mu_strategy", "monotone", true, true)
                .map_err(|e| Error::Solver(e.to_string()))?;
            app.initialize().map_err(|e| Error::Solver(e.to_string()))?;
            let problem = Rc::new(RefCell::new(problem));
            let status =
                app.optimize_tnlp_without_presolve(Rc::clone(&problem) as Rc<RefCell<dyn TNLP>>);
            let x = problem
                .borrow()
                .final_x()
                .ok_or_else(|| Error::Solver(format!("POUNCE returned no point: {status:?}")))?
                .to_vec();
            Ok(NlpSolution {
                x,
                status: status as i32,
                iterations: app.statistics().iteration_count.max(0) as usize,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oximo::prelude::*;
    use pounce_nlp::tnlp::SparsityRequest;

    #[test]
    fn exact_derivatives_include_curvature_of_constraint_rows() {
        let m = Model::new("pressure and energy derivatives");
        variable!(m,-2.0<=q<=2.0,initial=0.0);
        variable!(m,0.5<=h<=3.0,initial=1.0);
        constraint!(m, pressure, q * q.abs() + h.square() <= 5.0);
        constraint!(m, conservation, 2.0 * q - h + 3.0 == 0.0);
        objective!(m, Min, q.square() * h + q.exp());
        let mut evaluator = prepare(&m, Vec::new()).unwrap();
        let info = evaluator.get_nlp_info().unwrap();
        let mut rows = vec![0; info.nnz_jac_g as usize];
        let mut cols = rows.clone();
        assert!(evaluator.eval_jac_g(
            None,
            false,
            SparsityRequest::Structure {
                irow: &mut rows,
                jcol: &mut cols
            }
        ));
        let mut hr = vec![0; info.nnz_h_lag as usize];
        let mut hc = hr.clone();
        assert!(evaluator.eval_h(
            None,
            false,
            1.,
            None,
            false,
            SparsityRequest::Structure {
                irow: &mut hr,
                jcol: &mut hc
            }
        ));
        // Structure was prepared at zero flow; both flow directions must work.
        for flow in [-0.7, 0.6] {
            let point = [flow, 1.4];
            let mut gradient = [0.; 2];
            assert!(evaluator.eval_grad_f(&point, true, &mut gradient));
            assert!((gradient[0] - (2. * flow * 1.4 + flow.exp())).abs() < 1e-10);
            assert!((gradient[1] - flow * flow).abs() < 1e-10);
            let mut jac = vec![0.; rows.len()];
            assert!(evaluator.eval_jac_g(
                Some(&point),
                true,
                SparsityRequest::Values { values: &mut jac }
            ));
            let mut dense = vec![vec![0.; 2]; 2];
            for (k, &value) in jac.iter().enumerate() {
                dense[rows[k] as usize][cols[k] as usize] = value;
            }
            assert!((dense[0][0] - 2. * flow.abs()).abs() < 1e-10);
            assert!((dense[0][1] - 2.8).abs() < 1e-10);
            assert!((dense[1][0] - 2.).abs() < 1e-10);
            assert!((dense[1][1] + 1.).abs() < 1e-10);
            let mut hessian = vec![0.; hr.len()];
            let multipliers = [0.7, -0.4];
            assert!(evaluator.eval_h(
                Some(&point),
                true,
                1.,
                Some(&multipliers),
                true,
                SparsityRequest::Values {
                    values: &mut hessian
                }
            ));
            let mut dense = vec![vec![0.; 2]; 2];
            for (k, &value) in hessian.iter().enumerate() {
                let (i, j) = (hr[k] as usize, hc[k] as usize);
                dense[i][j] = value;
                dense[j][i] = value;
            }
            assert!((dense[0][0] - (2.8 + flow.exp() + 1.4 * flow.signum())).abs() < 1e-10);
            assert!((dense[0][1] - 2. * flow).abs() < 1e-10);
            assert!((dense[1][1] - 1.4).abs() < 1e-10);
        }
    }

    #[test]
    fn curve_segments_are_evaluated_again_when_the_operating_point_moves() {
        let m = Model::new("piecewise quadratic curve");
        variable!(m,0.0<=x<=3.0,initial=0.5);
        objective!(
            m,
            Min,
            (x - 1.)
                .max(oximo::expr::Expr::constant(x.arena(), 0.))
                .square()
        );
        let mut evaluator = prepare(&m, Vec::new()).unwrap();
        for (x, expected_value, expected_gradient) in [(0.5, 0., 0.), (2., 1., 2.), (0.25, 0., 0.)]
        {
            assert!((evaluator.eval_f(&[x], true).unwrap() - expected_value).abs() < 1e-12);
            let mut gradient = [0.];
            assert!(evaluator.eval_grad_f(&[x], true, &mut gradient));
            assert!((gradient[0] - expected_gradient).abs() < 1e-12);
        }
    }
}
