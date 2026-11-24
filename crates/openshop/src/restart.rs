//! Continue a scheduling horizon from independently reconstructed history.

use crate::{Error, Schedule, ScheduleCase, SchedulingProblem, SoftTarget, StorageHistory};

/// Execute the prefix of an on/off schedule and return its state and elapsed time.
fn state_at(
    grid: &[f64],
    mut on: bool,
    age: f64,
    states: impl Iterator<Item = bool>,
    time: f64,
) -> (bool, f64) {
    let mut since = grid[0] - age;
    for (&edge, state) in grid.iter().zip(states) {
        if state != on {
            since = edge;
            on = state;
        }
    }
    (on, time - since)
}

impl SchedulingProblem<'_> {
    /// Continue from an interior operating boundary, using the supplied controls
    /// on their hydraulic and river grids. The executed prefix becomes storage,
    /// release, flow, unit and gate history. Times remain absolute.
    ///
    /// The source controls must reconstruct a feasible schedule. An arrival
    /// averaging window or priced target cannot straddle the restart boundary.
    /// The new objective measures changes in water value from the restart state.
    pub fn restart_case(&self, schedule: &Schedule, time: f64) -> Result<ScheduleCase, Error> {
        let original = self.case;
        let edge = original
            .grid
            .iter()
            .position(|&x| x == time)
            .filter(|&i| i > 0 && i < original.prices.len())
            .ok_or_else(|| {
                Error::Input("restart time must be an interior operating boundary".into())
            })?;
        if schedule.operating_grid != original.grid {
            return Err(Error::Input(
                "restart schedule has a different operating grid".into(),
            ));
        }
        let problem = self.on_mesh(&schedule.grid, &schedule.routing_grid)?;
        let controls = problem.held_controls(schedule, &problem.mesh)?;
        let physical = problem.reconstruct(&controls)?;
        if !physical.validation.feasible {
            return Err(Error::Input(
                "restart controls violate physical equations or operating limits".into(),
            ));
        }
        let mesh = &problem.mesh;
        let vertex = mesh.vertices[edge];
        let previous = edge - 1;
        let duration = mesh.duration[previous];
        let mut next = original.clone();
        next.grid.drain(..edge);
        next.prices.drain(..edge);
        for row in next.inflows.values_mut() {
            row.drain(..edge);
        }
        for row in next.commitment.values_mut() {
            row.drain(..edge);
        }
        next.routing_grid = physical
            .routing_grid
            .iter()
            .copied()
            .filter(|t| *t >= time)
            .collect();

        // Exact history needs the routed release series, rather than hydraulic
        // cell averages, because a delayed inflow can change inside a cell.
        let mut columns: Vec<_> = physical
            .unit_flow
            .iter()
            .flatten()
            .copied()
            .chain(physical.tunnel_flow.iter().flatten().copied())
            .chain(
                problem
                    .routing
                    .outlets
                    .iter()
                    .flat_map(|&e| physical.river_withdrawal[e].iter().copied()),
            )
            .collect();
        problem.routing.complete_columns(&mut columns);
        for (r, reservoir) in next.system.reservoirs.iter_mut().enumerate() {
            reservoir.initial_storage = physical.storage[r][vertex];
            let windows = &mut reservoir.operation.level_windows;
            windows.retain(|w| w.end_time > time);
            if windows.iter().any(|w| w.start_time < time) || reservoir.operation.history.is_some()
            {
                let history = reservoir
                    .operation
                    .history
                    .get_or_insert_with(|| StorageHistory {
                        times: vec![mesh.grid[0]],
                        storage: vec![physical.storage[r][0]],
                    });
                for k in 0..vertex {
                    for (fraction, correction) in &problem.routing.storage_shape[r][k] {
                        if *fraction == 0. {
                            continue;
                        }
                        let end = if *fraction == 1. {
                            mesh.grid[k + 1]
                        } else {
                            mesh.grid[k] + fraction * mesh.dt[k]
                        };
                        if end <= *history.times.last().unwrap() {
                            continue;
                        }
                        let volume = if *fraction == 1. {
                            physical.storage[r][k + 1]
                        } else {
                            (1. - fraction) * physical.storage[r][k]
                                + fraction * physical.storage[r][k + 1]
                                + correction.value(&columns)
                        };
                        history.times.push(end);
                        history.storage.push(volume);
                    }
                }
            }
        }

        let unit_flow = mesh.means(&physical.unit_flow, 0.)?;
        let unit_power = mesh.means(&physical.unit_power, 0.)?;
        for (g, unit) in next.system.units.iter_mut().enumerate() {
            let op = &mut unit.operation;
            (op.initial_on, op.initial_duration) = state_at(
                &original.grid,
                op.initial_on,
                op.initial_duration,
                controls.commitment[g][..edge].iter().copied(),
                time,
            );
            op.initial_flow = Some(unit_flow[g][previous]);
            op.initial_power = Some(unit_power[g][previous]);
            op.initial_interval_hours = duration;
        }
        for (p, plant) in next.system.plants.iter_mut().enumerate() {
            let members = &self.topology.plant_units[p];
            let prior = plant.operation.resolve(
                &members
                    .iter()
                    .map(|&g| &original.system.units[g])
                    .collect::<Vec<_>>(),
                &plant.name,
            )?;
            let states = (0..edge).map(|t| members.iter().any(|&g| controls.commitment[g][t]));
            let (_, age) = state_at(
                &original.grid,
                prior.initial_on,
                prior.initial_duration,
                states,
                time,
            );
            plant.operation.initial_duration = Some(age);
            plant.operation.initial_flow =
                Some(members.iter().map(|&g| unit_flow[g][previous]).sum());
            plant.operation.initial_power =
                Some(members.iter().map(|&g| unit_power[g][previous]).sum());
            plant.operation.initial_interval_hours = duration;
        }

        for (e, tunnel) in next.system.tunnels.iter_mut().enumerate() {
            tunnel.initial_opening = Some(controls.tunnel_opening[e][previous]);
        }
        for rule in &mut next.system.tunnel_interlocks {
            for wait in &mut rule.closures {
                let mut until = original.grid[0] + wait.initial_remaining;
                for name in &wait.tunnels {
                    let e = original
                        .system
                        .tunnels
                        .iter()
                        .position(|t| &t.name == name)
                        .unwrap();
                    let tunnel = &original.system.tunnels[e];
                    let mut open = tunnel.initial_opening.unwrap_or(tunnel.opening) > 0.;
                    for (&edge, &opening) in original
                        .grid
                        .iter()
                        .zip(&controls.tunnel_opening[e][..previous + 1])
                    {
                        if open && opening == 0. {
                            until = until.max(edge + wait.hours);
                        }
                        open = opening > 0.;
                    }
                }
                wait.initial_remaining = (until - time).clamp(0., wait.hours);
            }
        }

        let releases = mesh.means(&physical.river_release, 0.)?;
        for (r, river) in next.system.rivers.iter_mut().enumerate() {
            let op = &mut river.operation;
            if !op.arrival_window_grid.is_empty() {
                let edge = op
                    .arrival_window_grid
                    .iter()
                    .position(|&x| x == time)
                    .ok_or_else(|| {
                        Error::Input(format!(
                            "restart crosses an arrival window for {}",
                            river.name
                        ))
                    })?;
                op.arrival_window_grid.drain(..edge);
            }
            op.initial_release = Some(releases[r][previous]);
            op.initial_interval_hours = duration;
            if river.law.is_delta() {
                continue;
            }
            let series = &problem.routing.releases[r];
            let start = series.grid[0];
            if river.history_grid.is_empty() {
                river.history_grid.push(start);
            } else if river.history_grid.last() != Some(&start) {
                river.history_grid.push(start);
                river.history_release.push(0.);
            }
            for (t, interval) in series.grid.windows(2).enumerate() {
                if interval[0] >= time {
                    break;
                }
                river.history_grid.push(interval[1].min(time));
                river
                    .history_release
                    .push(series.rates[t].value(&columns).max(0.));
            }
        }

        for (i, rule) in next.flow_requirements.iter_mut().enumerate() {
            if rule.min_volume.is_some() || rule.max_volume.is_some() {
                rule.initial_volume = physical.observed_volume[i][edge];
            }
            rule.volume_windows.retain(|w| w.end_time > time);
            if !rule.history_grid.is_empty()
                || !rule.volume_windows.is_empty()
                || !rule.rolling_ramps.is_empty()
                || !rule.period_ramps.is_empty()
            {
                let first = usize::from(!rule.history_grid.is_empty());
                rule.history_grid
                    .extend_from_slice(&original.grid[first..=edge]);
                rule.history_flow
                    .extend_from_slice(&physical.observed_flow[i][..edge]);
            }
        }
        for rule in &next.soft_constraints {
            let (a, b) = assessment_window(&rule.target);
            if a < time && time < b {
                return Err(Error::Input(format!(
                    "restart crosses priced target window {}",
                    rule.name
                )));
            }
        }
        next.soft_constraints
            .retain(|rule| assessment_window(&rule.target).1 > time);
        SchedulingProblem::new(&next)?;
        Ok(next)
    }
}

fn assessment_window(target: &SoftTarget) -> (f64, f64) {
    match target {
        SoftTarget::Storage { time, .. } | SoftTarget::Level { time, .. } => (*time, *time),
        SoftTarget::UnitPower {
            start_time,
            end_time,
            ..
        }
        | SoftTarget::PlantPower {
            start_time,
            end_time,
            ..
        }
        | SoftTarget::UnitDischarge {
            start_time,
            end_time,
            ..
        }
        | SoftTarget::PlantDischarge {
            start_time,
            end_time,
            ..
        } => (*start_time, *end_time),
    }
}
