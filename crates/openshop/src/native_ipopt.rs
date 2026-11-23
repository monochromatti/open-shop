//! Ownership and callback boundary for the system Ipopt C interface.
//!
//! Sparse derivatives come from the same TNLP evaluator used by POUNCE. Rust
//! panics cannot cross the C boundary. No native problem or callback outlives
//! the borrowed evaluator, and the problem is freed on every exit path.

use crate::{Error, SolveOptions, nlp::NlpSolution};
use ipopt_bindgen as ffi;
use pounce_nlp::tnlp::{BoundsInfo, SparsityRequest, StartingPoint, TNLP};
use std::{
    ffi::{CStr, c_void},
    panic::{AssertUnwindSafe, catch_unwind},
    slice,
};

struct Problem(ffi::IpoptProblem);
impl Drop for Problem {
    fn drop(&mut self) {
        // SAFETY: CreateIpoptProblem returned this owned non-null handle.
        unsafe {
            ffi::FreeIpoptProblem(self.0);
        }
    }
}

struct Callbacks<'a> {
    evaluator: &'a mut dyn TNLP,
    problem: ffi::IpoptProblem,
    panicked: bool,
    iterations: usize,
    current_x: Vec<f64>,
    best_x: Vec<f64>,
    best_objective: f64,
    lower_violation: Vec<f64>,
    upper_violation: Vec<f64>,
    constraint_violation: Vec<f64>,
    tolerance: f64,
}

pub(crate) fn solve(
    evaluator: &mut dyn TNLP,
    options: &SolveOptions,
) -> Result<NlpSolution, Error> {
    let info = evaluator
        .get_nlp_info()
        .ok_or_else(|| Error::Model("invalid nonlinear dimensions".into()))?;
    let (n, m) = (info.n as usize, info.m as usize);
    let (mut xl, mut xu, mut gl, mut gu) = (vec![0.; n], vec![0.; n], vec![0.; m], vec![0.; m]);
    if !evaluator.get_bounds_info(BoundsInfo {
        x_l: &mut xl,
        x_u: &mut xu,
        g_l: &mut gl,
        g_u: &mut gu,
    }) {
        return Err(Error::Model("invalid nonlinear bounds".into()));
    }
    let (mut x, mut zl, mut zu, mut lambda) = (vec![0.; n], vec![0.; n], vec![0.; n], vec![0.; m]);
    if !evaluator.get_starting_point(StartingPoint {
        init_x: true,
        x: &mut x,
        init_z: false,
        z_l: &mut zl,
        z_u: &mut zu,
        init_lambda: false,
        lambda: &mut lambda,
    }) {
        return Err(Error::Model("invalid nonlinear start".into()));
    }
    // SAFETY: The declared lengths match all bound buffers, and every callback
    // follows Ipopt's C signatures generated from this installation's headers.
    let raw = unsafe {
        ffi::CreateIpoptProblem(
            info.n,
            xl.as_mut_ptr(),
            xu.as_mut_ptr(),
            info.m,
            gl.as_mut_ptr(),
            gu.as_mut_ptr(),
            info.nnz_jac_g,
            info.nnz_h_lag,
            info.index_style as i32,
            Some(eval_f),
            Some(eval_g),
            Some(eval_grad),
            Some(eval_jac),
            Some(eval_h),
        )
    };
    if raw.is_null() {
        return Err(Error::Solver("Ipopt could not create the problem".into()));
    }
    let problem = Problem(raw);
    for (name, value) in [
        (c"tol", options.tolerance),
        (c"constr_viol_tol", options.tolerance),
        (c"bound_relax_factor", 0.0),
        (c"max_wall_time", options.time_limit),
    ] {
        // SAFETY: Ipopt copies the option and does not mutate its name.
        if !unsafe { ffi::AddIpoptNumOption(problem.0, name.as_ptr().cast_mut(), value) } {
            return Err(Error::Solver(format!(
                "Ipopt rejected {}",
                name.to_string_lossy()
            )));
        }
    }
    for (name, value) in [(c"print_level", 0), (c"max_iter", options.max_iterations)] {
        // SAFETY: The handle is live and the option name is NUL terminated.
        if !unsafe { ffi::AddIpoptIntOption(problem.0, name.as_ptr().cast_mut(), value) } {
            return Err(Error::Solver(format!(
                "Ipopt rejected {}",
                name.to_string_lossy()
            )));
        }
    }
    // A library solve must not depend on an unrelated ipopt.opt in its caller's
    // working directory. All numerical options are set at this boundary.
    set_string(&problem, c"option_file_name", c"")?;
    set_string(&problem, c"sb", c"yes")?;
    set_string(&problem, c"mu_strategy", c"monotone")?;
    set_string(&problem, c"hessian_approximation", c"limited-memory")?;
    // SAFETY: The callback remains valid for the entire synchronous solve.
    unsafe {
        ffi::SetIntermediateCallback(problem.0, Some(intermediate));
    }
    let mut g = vec![0.; m];
    let feasible = |x: &[f64], g: &[f64]| {
        x.iter().zip(&xl).zip(&xu).all(|((&v, &lo), &hi)| {
            v.is_finite() && v >= lo - options.tolerance && v <= hi + options.tolerance
        }) && g.iter().zip(&gl).zip(&gu).all(|((&v, &lo), &hi)| {
            v.is_finite() && v >= lo - options.tolerance && v <= hi + options.tolerance
        })
    };
    let initial_objective = evaluator.eval_f(&x, true).filter(|v| v.is_finite());
    let initial_feasible =
        initial_objective.is_some() && evaluator.eval_g(&x, false, &mut g) && feasible(&x, &g);
    let mut callbacks = Callbacks {
        evaluator,
        problem: problem.0,
        panicked: false,
        iterations: 0,
        current_x: vec![0.; n],
        best_x: x.clone(),
        best_objective: if initial_feasible {
            initial_objective.unwrap()
        } else {
            f64::INFINITY
        },
        lower_violation: vec![0.; n],
        upper_violation: vec![0.; n],
        constraint_violation: vec![0.; m],
        tolerance: options.tolerance,
    };
    let mut objective = 0.0;
    // SAFETY: All buffers match CreateIpoptProblem dimensions. The context and
    // evaluator stay exclusively borrowed until IpoptSolve returns.
    let status = unsafe {
        ffi::IpoptSolve(
            problem.0,
            x.as_mut_ptr(),
            g.as_mut_ptr(),
            &mut objective,
            lambda.as_mut_ptr(),
            zl.as_mut_ptr(),
            zu.as_mut_ptr(),
            (&mut callbacks as *mut Callbacks<'_>).cast(),
        )
    };
    if callbacks.panicked {
        return Err(Error::Solver("a nonlinear callback panicked".into()));
    }
    let final_objective = callbacks.evaluator.eval_f(&x, true);
    let final_feasible = callbacks.evaluator.eval_g(&x, false, &mut g) && feasible(&x, &g);
    if callbacks.best_objective.is_finite()
        && (!final_feasible
            || final_objective.is_none_or(|v| !v.is_finite() || callbacks.best_objective < v))
    {
        x.copy_from_slice(&callbacks.best_x);
    }
    Ok(NlpSolution {
        x,
        status,
        iterations: callbacks.iterations,
    })
}

fn set_string(problem: &Problem, name: &CStr, value: &CStr) -> Result<(), Error> {
    // SAFETY: The handle is live, strings are terminated and copied by Ipopt.
    if unsafe {
        ffi::AddIpoptStrOption(
            problem.0,
            name.as_ptr().cast_mut(),
            value.as_ptr().cast_mut(),
        )
    } {
        Ok(())
    } else {
        Err(Error::Solver(format!(
            "Ipopt rejected {}",
            name.to_string_lossy()
        )))
    }
}

// SAFETY: Called only by Ipopt with the live exclusive context passed to
// IpoptSolve; f cannot retain references to that context or to C buffers.
unsafe fn callback(data: *mut c_void, f: impl FnOnce(&mut Callbacks<'_>) -> bool) -> bool {
    // SAFETY: The caller upholds the context ownership contract above.
    let context = unsafe { &mut *data.cast::<Callbacks<'_>>() };
    if context.panicked {
        return false;
    }
    match catch_unwind(AssertUnwindSafe(|| f(context))) {
        Ok(ok) => ok,
        Err(_) => {
            context.panicked = true;
            false
        }
    }
}

unsafe fn input<'a, T>(p: *const T, n: i32) -> &'a [T] {
    if n == 0 {
        &[]
    } else {
        // SAFETY: Ipopt supplies n initialized values for each input buffer.
        unsafe { slice::from_raw_parts(p, n as usize) }
    }
}
unsafe fn output<'a, T>(p: *mut T, n: i32) -> &'a mut [T] {
    if n == 0 {
        &mut []
    } else {
        // SAFETY: Ipopt supplies n exclusive writable values for each output.
        unsafe { slice::from_raw_parts_mut(p, n as usize) }
    }
}

unsafe extern "C" fn eval_f(
    n: i32,
    x: *mut f64,
    new_x: bool,
    value: *mut f64,
    data: *mut c_void,
) -> bool {
    // SAFETY: Ipopt provides valid buffers and the context for this callback.
    unsafe {
        callback(data, |c| match c.evaluator.eval_f(input(x, n), new_x) {
            Some(v) if v.is_finite() => {
                *value = v;
                true
            }
            _ => false,
        })
    }
}
unsafe extern "C" fn eval_grad(
    n: i32,
    x: *mut f64,
    new_x: bool,
    grad: *mut f64,
    data: *mut c_void,
) -> bool {
    // SAFETY: The input and output arrays have n entries and do not alias.
    unsafe {
        callback(data, |c| {
            let out = output(grad, n);
            c.evaluator.eval_grad_f(input(x, n), new_x, out) && out.iter().all(|v| v.is_finite())
        })
    }
}
unsafe extern "C" fn eval_g(
    n: i32,
    x: *mut f64,
    new_x: bool,
    m: i32,
    g: *mut f64,
    data: *mut c_void,
) -> bool {
    // SAFETY: Ipopt supplies n inputs and m writable constraint values.
    unsafe {
        callback(data, |c| {
            let out = output(g, m);
            c.evaluator.eval_g(input(x, n), new_x, out) && out.iter().all(|v| v.is_finite())
        })
    }
}
unsafe extern "C" fn eval_jac(
    n: i32,
    x: *mut f64,
    new_x: bool,
    _m: i32,
    nz: i32,
    rows: *mut i32,
    cols: *mut i32,
    values: *mut f64,
    data: *mut c_void,
) -> bool {
    // SAFETY: Ipopt supplies either two structure buffers or one value buffer,
    // each with nz entries; x may be null only on the structure call.
    unsafe {
        callback(data, |c| {
            if values.is_null() {
                c.evaluator.eval_jac_g(
                    None,
                    new_x,
                    SparsityRequest::Structure {
                        irow: output(rows, nz),
                        jcol: output(cols, nz),
                    },
                )
            } else {
                let out = output(values, nz);
                c.evaluator.eval_jac_g(
                    Some(input(x, n)),
                    new_x,
                    SparsityRequest::Values { values: out },
                ) && out.iter().all(|v| v.is_finite())
            }
        })
    }
}
unsafe extern "C" fn eval_h(
    n: i32,
    x: *mut f64,
    new_x: bool,
    factor: f64,
    m: i32,
    lambda: *mut f64,
    new_lambda: bool,
    nz: i32,
    rows: *mut i32,
    cols: *mut i32,
    values: *mut f64,
    data: *mut c_void,
) -> bool {
    // SAFETY: Same structure/value convention as eval_jac; lambda has m entries.
    unsafe {
        callback(data, |c| {
            if values.is_null() {
                c.evaluator.eval_h(
                    None,
                    new_x,
                    factor,
                    None,
                    new_lambda,
                    SparsityRequest::Structure {
                        irow: output(rows, nz),
                        jcol: output(cols, nz),
                    },
                )
            } else {
                let out = output(values, nz);
                c.evaluator.eval_h(
                    Some(input(x, n)),
                    new_x,
                    factor,
                    Some(input(lambda, m)),
                    new_lambda,
                    SparsityRequest::Values { values: out },
                ) && out.iter().all(|v| v.is_finite())
            }
        })
    }
}
unsafe extern "C" fn intermediate(
    mode: i32,
    iter: i32,
    objective: f64,
    _pr: f64,
    _du: f64,
    _mu: f64,
    _dn: f64,
    _reg: f64,
    _ad: f64,
    _ap: f64,
    _ls: i32,
    data: *mut c_void,
) -> bool {
    // SAFETY: This synchronous callback borrows only its live context.
    unsafe {
        callback(data, |c| {
            c.iterations = iter.max(0) as usize;
            if mode != 0 || !objective.is_finite() || objective >= c.best_objective {
                return true;
            }
            let n = c.current_x.len() as i32;
            let m = c.constraint_violation.len() as i32;
            // SAFETY: The live problem owns the current iterate. All requested
            // buffers have its original dimensions; null outputs are optional.
            let checked = ffi::GetIpoptCurrentViolations(
                c.problem,
                false,
                n,
                c.lower_violation.as_mut_ptr(),
                c.upper_violation.as_mut_ptr(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
                m,
                c.constraint_violation.as_mut_ptr(),
                std::ptr::null_mut(),
            );
            if !checked
                || c.lower_violation
                    .iter()
                    .chain(&c.upper_violation)
                    .chain(&c.constraint_violation)
                    .any(|v| !v.is_finite() || *v > c.tolerance)
            {
                return true;
            }
            // SAFETY: This is the same live problem and n-sized buffer. Reading
            // the iterate does not change the solver or retain its C storage.
            let read = ffi::GetIpoptCurrentIterate(
                c.problem,
                false,
                n,
                c.current_x.as_mut_ptr(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
                m,
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            );
            if read && c.current_x.iter().all(|v| v.is_finite()) {
                c.best_x.copy_from_slice(&c.current_x);
                c.best_objective = objective;
            }
            true
        })
    }
}
