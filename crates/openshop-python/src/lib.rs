//! Python conversion and ownership boundary. Numerical calculations stay in Rust.

use openshop::{Error, ScheduleCase, SolveOptions};
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};

/// Optimize a watercourse supplied as a Python dictionary.
///
/// All inputs are copied into an owned Rust case before releasing Python.
/// Results contain physical quantities in the same units as the Rust API.
#[pyfunction]
#[pyo3(signature=(case, *, time_limit=60.0, max_iterations=1000, tolerance=1e-8))]
fn solve(
    py: Python<'_>,
    case: &Bound<'_, PyAny>,
    time_limit: f64,
    max_iterations: i32,
    tolerance: f64,
) -> PyResult<Py<PyAny>> {
    let case: ScheduleCase =
        pythonize::depythonize(case).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let options = SolveOptions {
        time_limit,
        max_iterations,
        tolerance,
    };
    let schedule = py
        .detach(move || openshop::solve(&case, &options))
        .map_err(|e| match e {
            Error::Input(_) => PyValueError::new_err(e.to_string()),
            _ => PyRuntimeError::new_err(e.to_string()),
        })?;
    Ok(pythonize::pythonize(py, &schedule)?.unbind())
}

#[pymodule]
fn _openshop(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(solve, module)?)?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
