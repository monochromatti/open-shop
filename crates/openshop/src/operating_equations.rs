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
    operating: &crate::operating_data::OperatingData,
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
        for row in &topology.unit_commitment[g] {
            let lhs = row
                .coefficients
                .iter()
                .fold(zero, |sum, &(t, a)| sum + on[g][t] * a);
            constraint!(m, lhs >= row.rhs);
        }
        for t in 0..nt {
            let limits = &operating.units[g][t];
            if let Some(fixed) = limits.forced {
                constraint!(m, on[g][t] == if fixed { 1. } else { 0. });
            }
            if let Some(p) = limits.fixed_power {
                constraint!(m, power[g][t] / scale_p == p / scale_p);
            }
            if let Some(q) = limits.fixed_flow {
                constraint!(m, flow[g][t] / scale_q == q / scale_q);
            }
            let prior_q = if t == 0 {
                u.operation.previous_flow().map(|q| zero + q)
            } else {
                Some(flow[g][t - 1])
            };
            if let Some(prior) = prior_q {
                let change = (flow[g][t] - prior) / scale_q;
                if let Some(up) = limits.flow_change.up {
                    constraint!(m, change <= up / scale_q);
                }
                if let Some(down) = limits.flow_change.down {
                    constraint!(m, -change <= down / scale_q);
                }
            }
            let prior_p = if t == 0 {
                u.operation.previous_power().map(|p| zero + p)
            } else {
                Some(power[g][t - 1])
            };
            if let Some(prior) = prior_p {
                let change = (power[g][t] - prior) / scale_p;
                if let Some(up) = limits.power_change.up {
                    constraint!(
                        m,
                        change <= (up + starts[g][t] * limits.power.min) / scale_p
                    );
                }
                if let Some(down) = limits.power_change.down {
                    constraint!(
                        m,
                        -change <= (down + stops[g][t] * limits.stopping_power_min) / scale_p
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
            let limits = &operating.plants[i][t];
            if let Some(fixed) = limits.forced {
                constraint!(m, plant_on[i][t] == if fixed { 1. } else { 0. });
            }
            if let Some(p) = limits.fixed_power {
                constraint!(m, total_p[t] / scale_p == p / scale_p);
            }
            if let Some(q) = limits.fixed_flow {
                constraint!(m, total_q[t] / scale_q == q / scale_q);
            }
            constraint!(
                m,
                total_q[t] / scale_q >= plant_on[i][t] * (limits.flow.min / scale_q)
            );
            constraint!(
                m,
                total_p[t] / scale_p >= plant_on[i][t] * (limits.power.min / scale_p)
            );
            if let Some(hi) = Some(limits.flow.max) {
                constraint!(m, total_q[t] / scale_q <= hi / scale_q);
            }
            if let Some(hi) = Some(limits.power.max) {
                constraint!(m, total_p[t] / scale_p <= hi / scale_p);
            }

            let prior_q = if t == 0 {
                operating.plant_initial_flow[i].map(|q| zero + q)
            } else {
                Some(total_q[t - 1])
            };
            if let Some(prior) = prior_q {
                let change = (total_q[t] - prior) / scale_q;
                if let Some(up) = limits.flow_change.up {
                    constraint!(m, change <= up / scale_q);
                }
                if let Some(down) = limits.flow_change.down {
                    constraint!(m, -change <= down / scale_q);
                }
            }
            let prior_p = if t == 0 {
                operating.plant_initial_power[i].map(|v| zero + v)
            } else {
                Some(total_p[t - 1])
            };
            if let Some(prior) = prior_p {
                let change = (total_p[t] - prior) / scale_p;
                let jump_up = members.iter().fold(zero, |sum, &g| {
                    sum + starts[g][t] * operating.units[g][t].power.min
                });
                let jump_down = members.iter().fold(zero, |sum, &g| {
                    sum + stops[g][t] * operating.units[g][t].stopping_power_min
                });
                if let Some(up) = limits.power_change.up {
                    constraint!(m, change <= (up + jump_up) / scale_p);
                }
                if let Some(down) = limits.power_change.down {
                    constraint!(m, -change <= (down + jump_down) / scale_p);
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
