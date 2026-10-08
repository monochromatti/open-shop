#!/usr/bin/env python3
"""Plot matched global-bound trajectories (Matplotlib 3.9.4)."""
import argparse
import pathlib
import statistics
import matplotlib.pyplot as plt
from root_proof_summary import load_rows, summarize

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('folder')
parser.add_argument('--output', required=True)
parser.add_argument('--profiles', default='current,obbt_1,joint')
args = parser.parse_args()
rows = load_rows(args.folder)
summarize(rows)  # Require valid, matched evidence before plotting.
profiles = args.profiles.split(',')
labels = dict(current='Current', obbt_1='OBBT ×1', obbt_quarter='OBBT ×0.25',
              obbt_cap30='OBBT boundary budget', joint='Joint power cuts',
              joint_obbt_1='Joint cuts + OBBT ×1')
colors = ['#48535d', '#087b9b', '#c46220', '#537b35', '#8566ab']
fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), sharex=True, sharey=True)
for ax, case in zip(axes.flat, ('operating-2h', 'operating-6h', 'operating-24h', 'seasonal-24h')):
    selected = [r for r in rows if r['benchmark_label'] == case]
    if not selected:
        ax.set_visible(False)
        continue
    allowance = max(r['allowance_seconds'] for r in selected)
    selected = [r for r in selected if r['allowance_seconds'] == allowance]
    lower = selected[0]['threshold_audited_lower_bound']
    for profile, color in zip(profiles, colors):
        runs = [r for r in selected if r['profile'] == profile]
        if not runs:
            continue
        times = sorted({0., float(allowance)} | {p['seconds'] for r in runs
            for p in r['global_bound_trajectory'] if p['seconds'] <= allowance})
        def at(run, t):
            points = [p['upper'] for p in run['global_bound_trajectory'] if p['seconds'] <= t]
            return 100 * (points[-1] - lower) / max(1, abs(lower)) if points else float('inf')
        values = [statistics.median(at(r,t) for r in runs) for t in times]
        ax.step(times, values, where='post', color=color, linewidth=1.6, label=labels[profile])
    ax.set_title(case.replace('operating-', 'Operating ').replace('seasonal-', 'Seasonal '), loc='left')
    ax.grid(axis='y', color='#e5e8eb', linewidth=.7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlim(0, allowance)
    ax.set_ylim(0, 11)
for ax in axes[:,0]:
    ax.set_ylabel('Global bound gap (%)')
for ax in axes[-1,:]:
    ax.set_xlabel('Elapsed proof time (seconds)')
handles, legend_labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc='lower center', ncol=len(handles), frameon=False, bbox_to_anchor=(.5,.025))
fig.suptitle('SCIP proof progress on matched operating cases', x=.065, ha='left', fontsize=15)
fig.text(.065, .005, 'Median trajectories against one frozen audited objective per case. Only globally valid bound observations.', fontsize=9, color='#48535d')
fig.tight_layout(rect=(0,.07,1,.95))
path = pathlib.Path(args.output)
path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(path, bbox_inches='tight')
