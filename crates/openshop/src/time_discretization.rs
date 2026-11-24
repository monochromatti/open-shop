//! Hydraulic cells inside unchanged operating intervals.
//!
//! Rates are projected by duration-weighted means and storage states by their
//! original time vertices. Holding a control
//! and projecting it back is the identity.
use crate::Error;
use std::ops::{Add, Mul, Range};

pub(crate) fn check_grid(grid: &[f64]) -> Result<(), Error> {
    if grid.len() < 2
        || !grid.iter().all(|t| t.is_finite())
        || !grid
            .windows(2)
            .all(|w| w[1] > w[0] && (3600. * (w[1] - w[0])).is_finite())
    {
        Err(Error::Input(
            "time grids must have finite ordered edges and durations".into(),
        ))
    } else {
        Ok(())
    }
}

pub(crate) struct TimeDiscretization {
    pub grid: Vec<f64>,
    pub interval: Vec<usize>,
    pub vertices: Vec<usize>,
    pub cells: Vec<Range<usize>>,
    pub dt: Vec<f64>,
    pub duration: Vec<f64>,
}
impl TimeDiscretization {
    pub fn new(operating: &[f64], hydraulic: &[f64]) -> Result<Self, Error> {
        check_grid(operating)?;
        check_grid(hydraulic)?;
        if operating.first() != hydraulic.first() || operating.last() != hydraulic.last() {
            return Err(Error::Input(
                "hydraulic and operating grids must span the same horizon".into(),
            ));
        }
        let vertices = operating
            .iter()
            .map(|&time| {
                let k = hydraulic.partition_point(|&t| t < time);
                if hydraulic.get(k) != Some(&time) {
                    Err(Error::Input(
                        "hydraulic grid must retain every operating boundary".into(),
                    ))
                } else {
                    Ok(k)
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let cells: Vec<_> = vertices.windows(2).map(|v| v[0]..v[1]).collect();
        let mut interval = vec![0; hydraulic.len() - 1];
        for (j, cells) in cells.iter().enumerate() {
            interval[cells.clone()].fill(j);
        }
        Ok(Self {
            grid: hydraulic.to_vec(),
            interval,
            vertices,
            cells,
            dt: hydraulic.windows(2).map(|w| w[1] - w[0]).collect(),
            duration: operating.windows(2).map(|w| w[1] - w[0]).collect(),
        })
    }
    pub fn means<D: Copy + Add<Output = D> + Mul<f64, Output = D>>(
        &self,
        rates: &[Vec<D>],
        zero: D,
    ) -> Result<Vec<Vec<D>>, Error> {
        if rates.iter().any(|row| row.len() != self.dt.len()) {
            return Err(Error::Input(
                "rate columns must match hydraulic cells".into(),
            ));
        }
        Ok(rates
            .iter()
            .map(|row| {
                self.cells
                    .iter()
                    .enumerate()
                    .map(|(j, cells)| {
                        cells.clone().fold(zero, |sum, k| {
                            sum + row[k] * (self.dt[k] / self.duration[j])
                        })
                    })
                    .collect()
            })
            .collect())
    }
    pub fn states<D: Clone>(&self, states: &[Vec<D>]) -> Result<Vec<Vec<D>>, Error> {
        if states.iter().any(|row| row.len() != self.grid.len()) {
            return Err(Error::Input(
                "state columns must match hydraulic vertices".into(),
            ));
        }
        Ok(states
            .iter()
            .map(|row| self.vertices.iter().map(|&k| row[k].clone()).collect())
            .collect())
    }
    pub fn hold<D: Clone>(&self, controls: &[Vec<D>]) -> Result<Vec<Vec<D>>, Error> {
        if controls.iter().any(|row| row.len() != self.duration.len()) {
            return Err(Error::Input(
                "control columns must match operating intervals".into(),
            ));
        }
        Ok(controls
            .iter()
            .map(|row| self.interval.iter().map(|&j| row[j].clone()).collect())
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn held_rates_preserve_original_discharge_means_and_water_volume() {
        let m = TimeDiscretization::new(&[0., 0.75, 2.], &[0., 0.25, 0.5, 0.75, 1.25, 2.]).unwrap();
        let controls = vec![vec![12., 30.], vec![7., 9.]];
        let fine = m.hold(&controls).unwrap();
        let means = m.means(&fine, 0.).unwrap();
        for (row, original) in means.iter().zip(&controls) {
            for (q, expected) in row.iter().zip(original) {
                assert!((q - expected).abs() < 1e-12);
            }
        }
        for (row, original) in fine.iter().zip(&controls) {
            let fine_volume: f64 = row.iter().zip(&m.dt).map(|(q, dt)| 3600. * q * dt).sum();
            let original_volume: f64 = original
                .iter()
                .zip(&m.duration)
                .map(|(q, dt)| 3600. * q * dt)
                .sum();
            assert!((fine_volume - original_volume).abs() < 1e-8);
        }
    }
    #[test]
    fn time_averages_preserve_energy_and_state_projection_retains_operating_edges() {
        let m = TimeDiscretization::new(&[0., 1., 3.], &[0., 0.25, 1., 2., 3.]).unwrap();
        let power = vec![vec![4., 8., 12., 20.]];
        let average = m.means(&power, 0.).unwrap();
        assert!((average[0][0] - 7.).abs() < 1e-12 && (average[0][1] - 16.).abs() < 1e-12);
        let storage = vec![m.grid.iter().map(|t| 1e6 - 3600. * 10. * t).collect()];
        let states = m.states(&storage).unwrap();
        assert!((states[0][0] - states[0][1] - 36000.).abs() < 1e-8);
        assert!((states[0][1] - states[0][2] - 72000.).abs() < 1e-8);
    }
    #[test]
    fn invalid_partitions_and_incompatible_quantity_axes_are_rejected() {
        for grid in [
            vec![0., 2.],
            vec![0., 1., 1., 2.],
            vec![0., f64::NAN, 2.],
            vec![0., 1., 3.],
        ] {
            assert!(TimeDiscretization::new(&[0., 1., 2.], &grid).is_err());
        }
        let m = TimeDiscretization::new(&[0., 2.], &[0., 1., 2.]).unwrap();
        assert!(m.means(&[vec![1.]], 0.).is_err());
        assert!(m.states(&[vec![1., 2.]]).is_err());
        assert!(m.hold(&[vec![1., 2.]]).is_err());
    }
}
