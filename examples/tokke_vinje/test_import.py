#!/usr/bin/env python3
"""Repeatable importer checks; synthetic rule fixtures need no external data.

Pass --source to also check the pinned historical operating and seasonal cases.
The integration check reads external data without running Julia or SHOP.
"""
import argparse
import copy
import importlib.util
import pathlib
import sys
import tempfile
import unittest

HERE = pathlib.Path(__file__).resolve().parent
sys.dont_write_bytecode = True
SPEC = importlib.util.spec_from_file_location("tokke_importer", HERE / "import.py")
IMPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPORTER)
SOURCE = None


class RuleTests(unittest.TestCase):
    def rule(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "rule.csv"
            path.write_text(text)
            return IMPORTER.read_seasonal_rule(path)

    def test_annual_baseline_and_inclusive_utc_transition(self):
        rule = self.rule("tid,verdi\n1990-12-31,3\n1991-06-17,7\n1991-11-18,3\n")
        for year in (2023, 2024, 2025):
            self.assertEqual(
                IMPORTER.seasonal_value(
                    rule, IMPORTER.timestamp(f"{year}-01-01T00:00Z")
                ),
                3,
            )
            self.assertEqual(
                IMPORTER.seasonal_value(
                    rule, IMPORTER.timestamp(f"{year}-06-16T23:00Z")
                ),
                3,
            )
            self.assertEqual(
                IMPORTER.seasonal_value(
                    rule, IMPORTER.timestamp(f"{year}-06-17T00:00Z")
                ),
                7,
            )
        # A timezone offset changes the effective UTC date.
        self.assertEqual(
            IMPORTER.seasonal_value(rule, IMPORTER.timestamp("2024-06-17T01:00+02:00")),
            3,
        )

    def test_interior_steps_and_terminal_boundary(self):
        rule = self.rule("tid,verdi\n1990-12-31,3\n1991-09-16,7\n")
        at = IMPORTER.timestamp("2024-09-15T12:00Z")
        self.assertEqual(IMPORTER.seasonal_series(rule, at, 48), ([0.0, 12.0], [3, 7]))
        # A change at the final instant belongs to the next horizon.
        self.assertEqual(IMPORTER.seasonal_series(rule, at, 12), ([0.0], [3]))
        at = IMPORTER.timestamp("2024-12-31T12:00Z")
        self.assertEqual(IMPORTER.seasonal_series(rule, at, 24), ([0.0, 12.0], [7, 3]))

    def test_reject_invalid_source_rules(self):
        for text in (
            "date,value\n1990-12-31,3\n",
            "tid,verdi\n1991-06-17,3\n",
            "tid,verdi\n1990-12-31,3\n1991-06-17,7\n1992-06-17,9\n",
            "tid,verdi\n1990-12-31,nan\n",
            "tid,verdi\n1990-12-31,-1\n",
        ):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.rule(text)

    def fixture(self, source):
        folder = source / "data/input_data/SEnDHub"
        folder.mkdir(parents=True)
        for filename in IMPORTER.RIVER_MINIMUM_FILES.values():
            (folder / filename).write_text("tid,verdi\n1990-12-31,3\n1991-09-16,7\n")
        for filename in IMPORTER.RESERVOIR_MINIMUM_FILES.values():
            (folder / filename).write_text(
                "tid,verdi\n1990-12-31,0\n1991-07-01,10\n1991-09-16,0\n"
            )
        (folder / "qmin_7891.csv").write_text("tid,verdi\n1990-12-31,6\n1991-09-16,9\n")
        return self.case_fixture()

    def case_fixture(self):
        rivers = set(IMPORTER.RIVER_MINIMUM_FILES) | set(IMPORTER.VEST_RIVERS)
        return {
            "name": "fixture",
            "operations": [],
            "rivers": [{"name": n} for n in sorted(rivers)],
            "reservoirs": [
                {"name": n, "v0": 20.0, "vmin": 2.0, "vmax": 50.0}
                for n in IMPORTER.RESERVOIR_MINIMUM_FILES
            ],
            "generators": [{"name": "Lio_G1", "plant": "Lio"}],
        }

    def test_storage_units_and_observation_without_second_water_flow(self):
        with tempfile.TemporaryDirectory() as directory:
            source = pathlib.Path(directory)
            case = self.fixture(source)
            before = copy.deepcopy(case)
            inflows = {IMPORTER.VEST_NAME: [1.25] * 48}
            coverage = IMPORTER.apply_operating_rules(
                case, source, IMPORTER.timestamp("2024-09-15T12:00Z"), 48, inflows
            )
            self.assertEqual(len(coverage), 9)
            self.assertEqual(case["reservoirs"], before["reservoirs"])
            self.assertEqual(case["rivers"], before["rivers"])
            self.assertEqual(case["generators"], before["generators"])
            storage = next(x for x in case["operations"] if x["object"] == "Totak")
            self.assertEqual(storage["values"], [10.0, 2.0])
            observation = case["flow_requirements"][0]
            self.assertEqual(observation["generators"], ["Lio_G1"])
            self.assertEqual(observation["rivers"], list(IMPORTER.VEST_RIVERS))
            self.assertEqual(observation["inflow"], 1.25)
            self.assertEqual(observation["min_flow"], 6)
            local = next(
                x
                for x in case["operations"]
                if x["object"] == IMPORTER.VEST_NAME and x["attribute"] == "inflow"
            )
            self.assertEqual(local["values"], inflows[IMPORTER.VEST_NAME])

    def test_invalid_initial_state_is_rejected_without_clipping(self):
        with tempfile.TemporaryDirectory() as directory:
            source = pathlib.Path(directory)
            case = self.fixture(source)
            case["reservoirs"][0]["v0"] = 9.0
            with self.assertRaisesRegex(ValueError, "no waiver applied"):
                IMPORTER.apply_operating_rules(
                    case,
                    source,
                    IMPORTER.timestamp("2024-07-02T00:00Z"),
                    2,
                    {IMPORTER.VEST_NAME: [1.0, 1.0]},
                )
            self.assertEqual(case["reservoirs"][0]["v0"], 9.0)

    def test_observation_rejects_missing_inflow_or_named_terms(self):
        with tempfile.TemporaryDirectory() as directory:
            source = pathlib.Path(directory)
            case = self.fixture(source)
            with self.assertRaisesRegex(ValueError, "dated Vest local inflow"):
                IMPORTER.apply_operating_rules(
                    case,
                    source,
                    IMPORTER.timestamp("2024-09-01T00:00Z"),
                    2,
                    {IMPORTER.VEST_NAME: [1.0]},
                )
            case = self.case_fixture()
            case["generators"] = []
            with self.assertRaisesRegex(ValueError, "observation terms"):
                IMPORTER.apply_operating_rules(
                    case,
                    source,
                    IMPORTER.timestamp("2024-09-01T00:00Z"),
                    2,
                    {IMPORTER.VEST_NAME: [1.0, 1.0]},
                )


class ExternalInputTests(unittest.TestCase):
    def test_real_initial_states_and_seasonal_crossing(self):
        if SOURCE is None:
            self.skipTest("pass --source for the external-data integration check")
        data, params, _ = IMPORTER.topology(SOURCE)
        for start, hours in (("2024-09-01T00:00Z", 24), ("2024-09-29T12:00Z", 48)):
            with self.subTest(start=start):
                at = IMPORTER.timestamp(start)
                levels, _ = IMPORTER.level_samples(SOURCE, at)
                prices, inflows, _, terminal = IMPORTER.historical_inputs(
                    SOURCE, at, hours
                )
                case, _ = IMPORTER.compile_case(
                    data, params, levels, prices, inflows, terminal, hours, SOURCE
                )
                before = copy.deepcopy(case)
                coverage = IMPORTER.apply_operating_rules(
                    case, SOURCE, at, hours, inflows
                )
                self.assertEqual(case["reservoirs"], before["reservoirs"])
                self.assertEqual(len(coverage), 9)
                self.assertEqual(len(case["flow_requirements"]), 1)
                self.assertEqual(len(case["generators"]), 14)
                self.assertEqual(len(case["reservoirs"]), 17)
                if hours == 48:
                    rule = next(
                        x
                        for x in case["operations"]
                        if x["object"] == "b_Bandak" and x["attribute"] == "min_release"
                    )
                    self.assertEqual(rule["times"], [0.0, 12.0])
                    self.assertEqual(rule["values"], [5.0, 4.0])
                    rule = next(
                        x
                        for x in case["operations"]
                        if x["object"] == "Staavatn" and x["attribute"] == "vmin"
                    )
                    self.assertEqual(rule["values"], [46.32, 0.0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--source", type=pathlib.Path)
    args, remaining = parser.parse_known_args()
    SOURCE = args.source
    unittest.main(argv=[sys.argv[0], *remaining])
