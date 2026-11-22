//! The same unit and aggregate operating equations for commitment and dispatch.

use crate::{ScheduleCase, model::Topology};
use oximo::prelude::*;

pub(crate) struct OperatingExpressions<'a> {
    pub power: Vec<Vec<Expr<'a>>>,
    pub flow: Vec<Vec<Expr<'a>>>,
    pub on: Vec<Vec<Expr<'a>>>,
    pub starts: Vec<Vec<Expr<'a>>>,
    pub stops: Vec<Vec<Expr<'a>>>,
    pub plant_on: Vec<Vec<Expr<'a>>>,
}
pub(crate) fn constrain(
    m: &Model,
    case: &ScheduleCase,
    topology: &Topology,
    variables: &OperatingExpressions<'_>,
) {
    let OperatingExpressions {
        power,
        flow,
        on,
        starts,
        stops,
        plant_on,
    } = variables;
    let nt = case.prices.len();
    for (g, u) in case.system.units.iter().enumerate() {
        let scale_q = u.flow.max.max(1.);
        let scale_p = u.power.max.max(1.);
        let zero = power[g][0] * 0.;
        for t in 0..nt {
            let dt = u.operation.ramp_hours(&case.grid, t);
            let prior_q = if t == 0 {
                u.operation.previous_flow().map(|q| zero + q)
            } else {
                Some(flow[g][t - 1])
            };
            if let Some(prior) = prior_q {
                let change = (flow[g][t] - prior) / scale_q;
                if let Some(up) = u.operation.flow_ramp.up {
                    constraint!(m, change <= up * dt / scale_q);
                }
                if let Some(down) = u.operation.flow_ramp.down {
                    constraint!(m, -change <= down * dt / scale_q);
                }
            }
            let prior_p = if t == 0 {
                u.operation.previous_power().map(|p| zero + p)
            } else {
                Some(power[g][t - 1])
            };
            if let Some(prior) = prior_p {
                let change = (power[g][t] - prior) / scale_p;
                if let Some(up) = u.operation.power_ramp.up {
                    constraint!(
                        m,
                        change <= (up * dt + starts[g][t] * u.power.min) / scale_p
                    );
                }
                if let Some(down) = u.operation.power_ramp.down {
                    constraint!(
                        m,
                        -change <= (down * dt + stops[g][t] * u.power.min) / scale_p
                    );
                }
            }
        }
    }
    for (i, p) in case.system.plants.iter().enumerate() {
        let members = &topology.plant_units[i];
        let zero = power[members[0]][0] * 0.;
        let total_q: Vec<_> = (0..nt)
            .map(|t| members.iter().fold(zero, |sum, &g| sum + flow[g][t]))
            .collect();
        let total_p: Vec<_> = (0..nt)
            .map(|t| members.iter().fold(zero, |sum, &g| sum + power[g][t]))
            .collect();
        let operation = &topology.plant_operation[i];
        let scale_q = p
            .flow_max
            .unwrap_or_else(|| members.iter().map(|&g| case.system.units[g].flow.max).sum())
            .max(1.);
        let scale_p = p
            .power_max
            .unwrap_or_else(|| {
                members
                    .iter()
                    .map(|&g| case.system.units[g].power.max)
                    .sum()
            })
            .max(1.);
        for row in &topology.plant_commitment[i] {
            let lhs = row
                .coefficients
                .iter()
                .fold(zero, |sum, &(t, a)| sum + plant_on[i][t] * a);
            constraint!(m, lhs >= row.rhs);
        }
        for t in 0..nt {
            constraint!(
                m,
                total_q[t] / scale_q >= plant_on[i][t] * (p.flow_min / scale_q)
            );
            constraint!(
                m,
                total_p[t] / scale_p >= plant_on[i][t] * (p.power_min / scale_p)
            );
            if let Some(hi) = p.flow_max {
                constraint!(m, total_q[t] / scale_q <= hi / scale_q);
            }
            if let Some(hi) = p.power_max {
                constraint!(m, total_p[t] / scale_p <= hi / scale_p);
            }
            let dt = operation.ramp_hours(&case.grid, t);
            let prior_q = if t == 0 {
                operation.previous_flow().map(|q| zero + q)
            } else {
                Some(total_q[t - 1])
            };
            if let Some(prior) = prior_q {
                let change = (total_q[t] - prior) / scale_q;
                if let Some(up) = operation.flow_ramp.up {
                    constraint!(m, change <= up * dt / scale_q);
                }
                if let Some(down) = operation.flow_ramp.down {
                    constraint!(m, -change <= down * dt / scale_q);
                }
            }
            let prior_p = if t == 0 {
                operation.previous_power().map(|v| zero + v)
            } else {
                Some(total_p[t - 1])
            };
            if let Some(prior) = prior_p {
                let change = (total_p[t] - prior) / scale_p;
                let jump_up = members.iter().fold(zero, |sum, &g| {
                    sum + starts[g][t] * case.system.units[g].power.min
                });
                let jump_down = members.iter().fold(zero, |sum, &g| {
                    sum + stops[g][t] * case.system.units[g].power.min
                });
                if let Some(up) = operation.power_ramp.up {
                    constraint!(m, change <= (up * dt + jump_up) / scale_p);
                }
                if let Some(down) = operation.power_ramp.down {
                    constraint!(m, -change <= (down * dt + jump_down) / scale_p);
                }
            }
        }
    }
    // Each group's maximum on-state plus every competing group's maximum is <=1.
    // Pairwise rows avoid adding an integer variable for a derived mode indicator.
    for rule in &topology.exclusions {
        for (i, a) in rule.iter().enumerate() {
            for b in &rule[..i] {
                for &g in a {
                    for &h in b {
                        for (&a, &b) in on[g].iter().zip(&on[h]) {
                            constraint!(m, a + b <= 1.);
                        }
                    }
                }
            }
        }
    }
}
