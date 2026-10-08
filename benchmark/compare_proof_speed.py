"""Compare paired proof-speed records; never substitute limits for censored times."""
import collections
import json
import statistics
import sys
from pathlib import Path


def gap(upper, lower):
    return max(0.0, upper - lower) / max(1.0, abs(lower))


def bound_at(row, seconds):
    observed = [p['upper'] for p in row['global_bound_trajectory'] if p['seconds'] <= seconds]
    return min(observed) if observed else None


def mean_gap(row, lower, start, end):
    upper = bound_at(row, start)
    if upper is None or end <= start:
        return None
    previous = start
    integral = 0.0
    for point in row['global_bound_trajectory']:
        seconds = point['seconds']
        if seconds <= start:
            continue
        if seconds >= end:
            break
        integral += (seconds - previous) * gap(upper, lower)
        upper = min(upper, point['upper'])
        previous = seconds
    integral += (end - previous) * gap(upper, lower)
    return integral / (end - start)


output = []
for path in sorted(Path(sys.argv[1]).rglob('summary.json')):
    rows = json.loads(path.read_text())
    lower = rows[0]['threshold_audited_lower_bound']
    # All profiles are judged over the same interval after each has a finite
    # globally valid bound. No value is invented before the first observation.
    start = max(r['global_bound_trajectory'][0]['seconds'] for r in rows)
    end = min(r['proof_observation_end_seconds'] for r in rows)
    ten_percent_times = [next((p['seconds'] for p in r['global_bound_trajectory']
        if gap(p['upper'], lower) <= 0.10), None) for r in rows]
    useful_start = max(ten_percent_times) if all(t is not None for t in ten_percent_times) else None
    groups = collections.defaultdict(list)
    for row in rows:
        groups[row['profile']].append(row)
    profiles = {}
    for profile, group in sorted(groups.items()):
        thresholds = []
        for index, threshold in enumerate(group[0]['time_to_gap']):
            samples = [r['time_to_gap'][index] for r in group]
            times = [s['first_observed_seconds'] for s in samples]
            thresholds.append(dict(gap_percent=threshold['gap_percent'], times=times,
                all_attained=all(s['attained'] for s in samples),
                median_seconds=statistics.median(times) if all(t is not None for t in times) else None))
        profiles[profile] = dict(
            objectives=[r['feasible_lower_bound'] for r in group],
            final_gap_percent=[100 * gap(r['global_bound'], lower) for r in group],
            median_final_gap_percent=statistics.median(100 * gap(r['global_bound'], lower) for r in group),
            certificates=[r['global_certificate'] for r in group],
            median_total_seconds=statistics.median(r['total_seconds'] for r in group),
            median_proof_seconds=statistics.median(r['proof_observation_end_seconds'] for r in group),
            median_mean_gap_percent=statistics.median(100 * mean_gap(r, lower, start, end) for r in group),
            median_mean_gap_percent_after_all_reach_10_percent=(
                statistics.median(100 * mean_gap(r, lower, useful_start, end) for r in group)
                if useful_start is not None and useful_start < end else None),
            thresholds=thresholds)
    output.append(dict(artifact=path.parent.name, common_audited_objective=lower,
        common_finite_bound_interval_seconds=[start, end],
        common_10_percent_bound_interval_seconds=[useful_start, end], profiles=profiles))
print(json.dumps(output, indent=2))
