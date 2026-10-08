"""Plot native globally valid bound trajectories for paired operating cases."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

folder, profile, output = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
cases = [('operating-2h', 'Tokke–Vinje, 2 h'), ('operating-6h', 'Tokke–Vinje, 6 h'),
         ('operating-24h', 'Tokke–Vinje, 24 h'), ('seasonal-24h', 'Seasonal Tokke–Vinje, 24 h')]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'svg.fonttype': 'none', 'axes.spines.top': False,
                     'axes.spines.right': False})
figure, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
for axis, (case, title) in zip(axes.flat, cases):
    rows = json.loads((folder / ('proof-speed-' + case) / 'summary.json').read_text())
    lower = rows[0]['threshold_audited_lower_bound']
    displayed = [r for r in rows if r['profile'] in ('baseline', profile)]
    bottom = min(100 * max(0, r['global_bound'] - lower) / max(1, abs(lower)) for r in displayed) - .2
    for row in displayed:
        baseline = row['profile'] == 'baseline'
        points = row['global_bound_trajectory']
        seconds = [p['seconds'] for p in points]
        gaps = [100 * max(0, p['upper'] - lower) / max(1, abs(lower)) for p in points]
        axis.step(seconds, gaps, where='post', color='#696969' if baseline else '#087f8c',
                  linewidth=1.6, alpha=.85, linestyle='-' if row['repeat'] == 1 else '--',
                  label=('v0.4.3 baseline' if baseline else 'Selected configuration') if row['repeat'] == 1 else None)
    axis.set(title=title, xlim=(0, 300), ylim=(bottom, 10),
             xlabel='Elapsed time (s)', ylabel='Gap above audited objective (%)')
    axis.set_xticks([0, 60, 120, 180, 240, 300])
    axis.grid(True, color='#dddddd', linewidth=.5)
handles, labels = axes[0, 0].get_legend_handles_labels()
figure.legend(handles, labels, loc='outside upper center', ncol=2, frameon=False,
              title='Two repetitions per configuration: solid and dashed')
output.parent.mkdir(parents=True, exist_ok=True)
figure.savefig(output, metadata={'Date': None})
figure.savefig(output.with_suffix('.png'), dpi=160)
