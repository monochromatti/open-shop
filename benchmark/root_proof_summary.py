#!/usr/bin/env python3
"""Validate and summarize matched root-proof scalar measurements."""
import argparse
import json
import math
import pathlib
import statistics
from collections import defaultdict


def first_upper(points, target):
    return next((p['seconds'] for p in points if p['upper'] <= target), None)


def load_rows(folder):
    rows = []
    for path in sorted(pathlib.Path(folder).rglob('summary.json')):
        data = json.loads(path.read_text())
        if isinstance(data, list):
            for row in data:
                row["benchmark_label"] = row.get("benchmark_label", path.parent.name.removeprefix("root-proof-"))
            rows.extend(data)
    if not rows:
        raise ValueError('no benchmark records')
    return rows


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        if row.get('experiment_error'):
            raise ValueError(row['experiment_error'])
        if not row['accepted'] or not row['start_audit']['valid']:
            raise ValueError('missing accepted schedule or complete start')
        if not math.isfinite(row['global_bound']):
            raise ValueError('missing finite global bound')
        if row.get('bound_event_errors') or row.get('bound_event_monotonic') is False:
            raise ValueError('invalid global-bound observations')
        if row['root_progress']['errors']:
            raise ValueError('root observer failed')
        if row['global_bound'] < row['matched_job_best_accepted_objective'] - 1e-6:
            raise ValueError('bound excludes an audited matched schedule')
        points = row['global_bound_trajectory']
        if any(p['upper'] < row['matched_job_best_accepted_objective'] - 1e-6 for p in points):
            raise ValueError('trajectory excludes an audited matched schedule')
        for key in ('power_cut_statistics', 'joint_supports'):
            stat = row.get(key)
            if stat and (stat['errors'] or stat.get('infeasible_flags', 0)):
                raise ValueError('separator failed')
        groups[(row.get('benchmark_label',row['case']), row['case_sha256'], row['allowance_seconds'])].append(row)
    results = []
    for (case, case_hash, allowance), group in sorted(groups.items()):
        for key in ('case_sha256', 'seed_controls_sha256', 'manifest_sha256',
                    'source_sha256', 'experiment_sha256', 'threshold_audited_lower_bound',
                    'julia_version', 'threads', 'blas_threads'):
            if len({r[key] for r in group}) != 1:
                raise ValueError(f'matched records differ in {key}: {case}')
        seen = [(r['profile'], r['repeat']) for r in group]
        if len(set(seen)) != len(seen):
            raise ValueError('duplicate profile/repetition')
        controls = [r for r in group if r['profile'] == 'current']
        if not controls:
            raise ValueError('missing production control')
        lower = group[0]['threshold_audited_lower_bound']
        # Use one bound attainable by every control repetition; report the
        # exact common target, never compare policy-specific gap denominators.
        common_upper = max(r['global_bound'] for r in controls)
        by_profile = defaultdict(list)
        for row in group:
            by_profile[row['profile']].append(row)
        for profile, records in sorted(by_profile.items()):
            if {r['repeat'] for r in records} != {r['repeat'] for r in controls}:
                raise ValueError('unmatched repetitions')
            times = [first_upper(r['global_bound_trajectory'], common_upper) for r in records]
            threshold_rows = []
            for target in (.10, .08, .07, .06, .05, .01, .0001):
                ts = [first_upper(r['global_bound_trajectory'], lower + target * max(1, abs(lower)))
                      for r in records]
                threshold_rows.append(dict(gap_percent=100*target,
                    attained=sum(t is not None for t in ts), repetitions=len(ts),
                    median_seconds=statistics.median(ts) if all(t is not None for t in ts) else None))
            first_branches = [r['root_progress']['first_root_branch'] for r in records]
            result = dict(case=case, case_name=group[0]["case"], case_sha256=case_hash, allowance_seconds=allowance, profile=profile,
                repetitions=len(records), audited_common_lower=lower,
                median_upper=statistics.median(r['global_bound'] for r in records),
                median_frozen_gap_percent=statistics.median(100*r['gap_against_frozen_audited_objective'] for r in records),
                median_delivered_objective=statistics.median(r['feasible_lower_bound'] for r in records),
                median_return_seconds=statistics.median(r['total_seconds'] for r in records),
                common_control_upper_target=common_upper,
                common_bound_attained=sum(t is not None for t in times),
                median_seconds_to_common_bound=statistics.median(times) if all(t is not None for t in times) else None,
                first_branch_attained=sum(b is not None for b in first_branches),
                median_seconds_to_first_branch=statistics.median(b['seconds'] for b in first_branches) if all(b is not None for b in first_branches) else None,
                median_nodes=statistics.median(r['scip_diagnostics']['nodes'] for r in records),
                median_root_observer_seconds=statistics.median(r['root_progress']['capture_seconds'] for r in records),
                time_to_gap=threshold_rows)
            def obbt(r):
                return r['scip_statistics']['propagator']['plugins']['obbt']
            result['median_obbt_seconds'] = statistics.median(obbt(r)['propagation_time'] for r in records)
            result['median_obbt_calls'] = statistics.median(obbt(r)['calls'] for r in records)
            result['obbt_boundary_stops'] = sum(r['root_progress']['obbt_disabled_at_boundary'] for r in records)
            if profile.startswith('joint'):
                result['median_joint_cuts'] = statistics.median((r.get('joint_supports') or {}).get('cuts_added_to_global_pool', 0) for r in records)
                result['median_joint_work_seconds'] = statistics.median((r.get('joint_supports') or {}).get('work_seconds', 0) for r in records)
            results.append(result)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    rows = load_rows(args.folder)
    comparison = summarize(rows)
    pathlib.Path(args.output).write_text(json.dumps(comparison, indent=2) + '\n')
    for r in comparison:
        print(f"{r['case']:48s} {r['profile']:16s} gap={r['median_frozen_gap_percent']:.4f}% "
              f"obbt={r['median_obbt_seconds']:.2f}s nodes={r['median_nodes']} "
              f"same-bound={r['median_seconds_to_common_bound']}")
