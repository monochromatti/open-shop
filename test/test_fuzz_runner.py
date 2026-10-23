"""Behavior of the slow-run supervisor; fake executables never invoke a solver."""

import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("fuzz", Path(__file__).parents[1] / "scripts" / "fuzz.py")
fuzz = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fuzz)


class FuzzRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.julia = self.directory / "fake-julia"
        self.julia.write_text(f"#!{sys.executable}\n" + '''
import json, os, pathlib, sys, time
request_path = pathlib.Path(sys.argv[-1])
request = json.loads(request_path.read_text())
mode = os.environ.get("FUZZ_TEST_MODE", "pass")
if mode == "timeout":
    time.sleep(10)
if mode == "crash":
    (request_path.parent / "partial.json").write_text('{"stages": {"local": "saved"}}')
    sys.exit(17)
if mode == "missing":
    sys.exit(0)
result = {"outcome": mode, "seed": request["seed"], "findings": []}
if mode == "review":
    result["findings"] = [{"kind": "search_inconclusive", "severity": "warning"}]
(request_path.parent / "result.json").write_text(json.dumps(result))
''')
        self.julia.chmod(0o755)

    def tearDown(self):
        self.temporary.cleanup()

    def run_campaign(self, *extra):
        args = ["--julia", str(self.julia), "--output", str(self.directory / "run")]
        with contextlib.redirect_stdout(io.StringIO()):
            return fuzz.main(args + list(extra))

    def test_resume_preserves_completed_inputs_and_extends_campaign(self):
        self.assertEqual(self.run_campaign("--cases", "2"), 0)
        request = (self.directory / "run" / "000001" / "request.json").read_bytes()
        result = (self.directory / "run" / "000001" / "result.json").read_bytes()
        self.assertEqual(self.run_campaign("--cases", "3", "--resume"), 0)
        self.assertEqual((self.directory / "run" / "000001" / "request.json").read_bytes(), request)
        self.assertEqual((self.directory / "run" / "000001" / "result.json").read_bytes(), result)
        summary = json.loads((self.directory / "run" / "summary.json").read_text())
        self.assertEqual(summary["completed"], 3)
        with self.assertRaises(SystemExit):
            self.run_campaign("--cases", "3", "--resume", "--seed", "12")

    def test_timeout_is_recorded_and_next_case_runs(self):
        with patch.dict(os.environ, FUZZ_TEST_MODE="timeout"):
            self.assertEqual(self.run_campaign("--cases", "2", "--timeout", "0.15"), 1)
        for i in (1, 2):
            result = json.loads((self.directory / "run" / f"{i:06d}" / "result.json").read_text())
            self.assertEqual(result["findings"][0]["kind"], "hard_timeout")

    def test_crash_retains_completed_stage_evidence(self):
        with patch.dict(os.environ, FUZZ_TEST_MODE="crash"):
            self.assertEqual(self.run_campaign("--cases", "1"), 1)
        result = json.loads((self.directory / "run" / "000001" / "result.json").read_text())
        self.assertEqual(result["findings"][0]["kind"], "process_exit")
        self.assertEqual(result["partial"]["stages"]["local"], "saved")

    def test_successful_exit_without_result_is_a_failure(self):
        with patch.dict(os.environ, FUZZ_TEST_MODE="missing"):
            self.assertEqual(self.run_campaign("--cases", "1"), 1)
        result = json.loads((self.directory / "run" / "000001" / "result.json").read_text())
        self.assertEqual(result["findings"][0]["kind"], "missing_result")

    def test_inconclusive_search_is_distinct_from_a_failure(self):
        with patch.dict(os.environ, FUZZ_TEST_MODE="review"):
            self.assertEqual(self.run_campaign("--cases", "1"), 0)
            self.assertEqual(self.run_campaign("--cases", "1", "--resume", "--strict"), 1)

    def test_rejects_destructive_output_reuse(self):
        self.assertEqual(self.run_campaign("--cases", "1"), 0)
        with self.assertRaises(SystemExit):
            self.run_campaign("--cases", "1")
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            fuzz.arguments(["--julia", str(self.julia), "--output", str(fuzz.ROOT / "results")])


if __name__ == "__main__":
    unittest.main()
