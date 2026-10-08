"""Checks that proof reports distinguish valid evidence from censored results.

Run with ``python3 -m unittest discover -s benchmark -p root_proof_summary_tests.py``.
Synthetic records use the same evidence schema as the matched benchmark jobs.
"""
import json
import pathlib
import tempfile
import unittest

from root_proof_summary import load_rows, summarize


def record(profile="current", repeat=1, upper=110.0, *, lower=100.0,
           delivered=None, allowance=100.0, end=100.0, observations=None,
           label="operating-24h", case_hash="ordinary-case", certificate=False):
    delivered = lower if delivered is None else delivered
    if observations is None:
        observations = [(1.0, 130.0), (end, upper)]
    return {
        "benchmark_label": label,
        # The importer gives both ordinary and seasonal cases this same name.
        "case": "Tokke–Vinje, 24 hours",
        "case_sha256": case_hash,
        "allowance_seconds": allowance,
        "seed_controls_sha256": "frozen-controls",
        "manifest_sha256": "locked-manifest",
        "source_sha256": "physical-source",
        "experiment_sha256": "experiment-source",
        "threshold_audited_lower_bound": lower,
        "julia_version": "1.13.1",
        "threads": 1,
        "blas_threads": 1,
        "commitment": "free",
        "scope": "whole declared discrete model",
        "native_scip_version": "10.0.3",
        "native_lp_solver": "SoPlex 8.0.3",
        "variable_count": 100,
        "constraint_count": 200,
        "profile": profile,
        "repeat": repeat,
        "accepted": True,
        "start_audit": {"valid": True},
        "global_bound": upper,
        "feasible_lower_bound": delivered,
        "matched_job_best_accepted_objective": delivered,
        "gap_against_frozen_audited_objective": (upper - lower) / max(1.0, abs(lower)),
        "relative_gap": (upper - delivered) / max(1.0, abs(delivered)),
        "global_certificate": certificate,
        "total_seconds": end + 0.5,
        "proof_observation_end_seconds": end,
        "bound_event_errors": [],
        "bound_event_monotonic": True,
        "global_bound_trajectory": [
            {"seconds": seconds, "upper": value, "native_seconds": seconds - 0.1, "run": 1}
            for seconds, value in observations
        ],
        "root_progress": {
            "errors": [], "first_root_branch": None,
            "capture_seconds": 0.001, "obbt_disabled_at_boundary": False,
        },
        "scip_diagnostics": {"nodes": 1},
        "scip_statistics": {"propagator": {"plugins": {"obbt": {
            "propagation_time": 10.0, "calls": 2,
        }}}},
        "power_cut_statistics": {"errors": [], "infeasible_flags": 0},
        "joint_supports": None,
    }


def profile_result(rows, profile="current"):
    return next(result for result in summarize(rows) if result["profile"] == profile)


class MatchedEvidenceTests(unittest.TestCase):
    def test_imported_names_do_not_merge_ordinary_and_seasonal_cases(self):
        ordinary = record()
        seasonal = record(case_hash="seasonal-case", label="seasonal-24h", lower=80.0, upper=90.0)
        del ordinary["benchmark_label"]
        del seasonal["benchmark_label"]
        with tempfile.TemporaryDirectory() as temporary:
            folder = pathlib.Path(temporary)
            for name, row in (("operating-24h", ordinary), ("seasonal-24h", seasonal)):
                artifact = folder / f"root-proof-{name}"
                artifact.mkdir()
                (artifact / "summary.json").write_text(json.dumps([row]))
            results = summarize(load_rows(folder))
        self.assertEqual(len(results), 2)
        self.assertEqual({r["case"] for r in results}, {"operating-24h", "seasonal-24h"})
        self.assertEqual(len({r["case_name"] for r in results}), 1)
        self.assertEqual({r["audited_common_lower"] for r in results}, {80.0, 100.0})

    def test_case_hash_separates_same_label_and_name(self):
        results = summarize([record(), record(case_hash="other-case", lower=80.0, upper=90.0)])
        self.assertEqual(len(results), 2)
        self.assertEqual({r["case_sha256"] for r in results}, {"ordinary-case", "other-case"})

    def test_mismatched_provenance_is_rejected(self):
        changes = {
            "seed_controls_sha256": "other-controls", "manifest_sha256": "other-lock",
            "source_sha256": "other-model", "experiment_sha256": "other-experiment",
            "threshold_audited_lower_bound": 99.0, "julia_version": "1.12.0",
            "threads": 2, "blas_threads": 2, "commitment": "fixed",
            "scope": "frozen commitment", "native_scip_version": "10.1.0",
            "native_lp_solver": "other-LP", "variable_count": 101, "constraint_count": 201,
        }
        for field, value in changes.items():
            with self.subTest(field=field):
                candidate = record("obbt_1", upper=108.0)
                candidate[field] = value
                with self.assertRaises(ValueError):
                    summarize([record(), candidate])

    def test_missing_control_is_rejected(self):
        with self.assertRaises(ValueError):
            summarize([record("obbt_1")])

    def test_duplicate_profile_repetition_is_rejected(self):
        with self.assertRaises(ValueError):
            summarize([record(), record()])

    def test_unmatched_repetitions_are_rejected(self):
        with self.assertRaises(ValueError):
            summarize([record(repeat=1), record(repeat=2), record("obbt_1", repeat=1)])

    def test_unaccepted_schedule_or_incomplete_start_is_rejected(self):
        for accepted, start_valid in ((False, True), (True, False)):
            with self.subTest(accepted=accepted, start_valid=start_valid):
                row = record()
                row["accepted"] = accepted
                row["start_audit"]["valid"] = start_valid
                with self.assertRaises(ValueError):
                    summarize([row])


class BoundValidityTests(unittest.TestCase):
    def test_final_bound_must_enclose_best_matched_audited_schedule(self):
        row = record(upper=110.0)
        row["matched_job_best_accepted_objective"] = 111.0
        with self.assertRaises(ValueError):
            summarize([row])

    def test_every_event_must_enclose_matched_schedule(self):
        row = record(observations=[(1.0, 130.0), (50.0, 99.0)])
        with self.assertRaises(ValueError):
            summarize([row])

    def test_nonfinite_bounds_and_timestamps_are_rejected(self):
        for field in ("upper", "seconds"):
            for value in (float("nan"), float("inf"), -float("inf")):
                with self.subTest(field=field, value=value):
                    row = record()
                    row["global_bound_trajectory"][0][field] = value
                    with self.assertRaises(ValueError):
                        summarize([row])
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.subTest(final_bound=value):
                row = record()
                row["global_bound"] = value
                with self.assertRaises(ValueError):
                    summarize([row])

    def test_appended_final_point_cannot_increase_bound(self):
        # The native event flag can be true while the separately appended final
        # observation is inconsistent. Check the actual sequence as well.
        row = record(upper=112.0, observations=[(1.0, 130.0), (90.0, 110.0), (100.0, 112.0)])
        self.assertTrue(row["bound_event_monotonic"])
        with self.assertRaises(ValueError):
            summarize([row])

    def test_out_of_order_event_time_is_rejected(self):
        with self.assertRaises(ValueError):
            summarize([record(observations=[(2.0, 130.0), (1.0, 110.0)])])

    def test_failed_callbacks_and_native_infeasibility_are_rejected(self):
        bad_rows = []
        row = record()
        row["experiment_error"] = "callback failed"
        bad_rows.append(row)
        row = record()
        row["bound_event_errors"] = ["observer failed"]
        bad_rows.append(row)
        row = record()
        row["bound_event_monotonic"] = False
        bad_rows.append(row)
        row = record()
        row["root_progress"]["errors"] = ["budget mutation failed"]
        bad_rows.append(row)
        row = record()
        row["power_cut_statistics"]["errors"] = ["certificate failed"]
        bad_rows.append(row)
        row = record()
        row["power_cut_statistics"]["infeasible_flags"] = 1
        bad_rows.append(row)
        row = record()
        row["joint_supports"] = {"errors": ["oracle failed"], "infeasible_flags": 0}
        bad_rows.append(row)
        row = record()
        row["joint_supports"] = {"errors": [], "infeasible_flags": 1}
        bad_rows.append(row)
        for number, row in enumerate(bad_rows):
            with self.subTest(failure=number):
                with self.assertRaises(ValueError):
                    summarize([row])


class TimingAndCertificateTests(unittest.TestCase):
    def test_unreached_threshold_is_censored_not_assigned_allowance(self):
        controls = [record(repeat=1), record(repeat=2)]
        candidates = [record("obbt_1", repeat=1, upper=106.0, end=30.0),
                      record("obbt_1", repeat=2, upper=111.0, end=80.0)]
        result = profile_result(controls + candidates, "obbt_1")
        self.assertEqual(result["common_bound_attained"], 1)
        self.assertIsNone(result["median_seconds_to_common_bound"])
        threshold = next(t for t in result["time_to_gap"] if abs(t["gap_percent"] - 7.0) < 1e-12)
        self.assertEqual(threshold["attained"], 1)
        self.assertEqual(threshold["repetitions"], 2)
        self.assertIsNone(threshold["median_seconds"])

    def test_attained_times_use_one_common_control_bound(self):
        rows = [record(repeat=1, upper=108.0), record(repeat=2, upper=110.0),
                record("obbt_1", repeat=1, upper=109.0, end=20.0),
                record("obbt_1", repeat=2, upper=109.0, end=40.0)]
        result = profile_result(rows, "obbt_1")
        self.assertEqual(result["common_control_upper_target"], 110.0)
        self.assertEqual(result["common_bound_attained"], 2)
        self.assertEqual(result["median_seconds_to_common_bound"], 30.0)

    def test_root_certificate_is_success_despite_frozen_seed_gap(self):
        # A guard finds a better schedule than its common seed. Its remaining
        # physical objective gap is tiny even though the frozen-seed gap is not.
        row = record(upper=100.620001, delivered=100.62, certificate=True)
        result = profile_result([row])
        self.assertAlmostEqual(result["median_frozen_gap_percent"], 0.620001)
        self.assertLess(result["median_actual_gap_percent"], 0.00001)
        self.assertEqual(result["global_certificates"], 1)
        self.assertEqual(result["first_branch_attained"], 0)
        self.assertIsNone(result["median_seconds_to_first_branch"])

    def test_nominal_allowance_excludes_later_native_completion(self):
        control = record(upper=105.0, end=101.0,
                         observations=[(1.0, 130.0), (99.0, 108.0), (101.0, 105.0)])
        candidate = record("obbt_1", upper=104.0, end=105.0,
                           observations=[(1.0, 130.0), (95.0, 110.0), (105.0, 104.0)])
        candidate["total_seconds"] = 109.0
        result = profile_result([control, candidate], "obbt_1")
        self.assertEqual(result["median_upper"], 104.0)
        self.assertEqual(result["median_upper_at_nominal_allowance"], 110.0)
        self.assertEqual(result["median_frozen_gap_at_nominal_allowance_percent"], 10.0)
        self.assertEqual(result["common_bound_attained"], 1)
        self.assertEqual(result["median_seconds_to_common_bound"], 105.0)
        self.assertEqual(result["common_bound_attained_within_allowance"], 0)
        self.assertIsNone(result["median_seconds_to_common_bound_within_allowance"])
        self.assertEqual(result["median_proof_seconds"], 105.0)
        self.assertEqual(result["median_return_seconds"], 109.0)

    def test_no_bound_before_allowance_is_unavailable(self):
        row = record(upper=105.0, end=101.0, observations=[(101.0, 105.0)])
        result = profile_result([row])
        self.assertIsNone(result["median_upper_at_nominal_allowance"])
        self.assertIsNone(result["median_frozen_gap_at_nominal_allowance_percent"])
        self.assertEqual(result["common_bound_attained_within_allowance"], 0)


if __name__ == "__main__":
    unittest.main()
