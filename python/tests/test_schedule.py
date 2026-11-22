from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import pytest

from openshop import solve


def watercourse():
    return {
        "system": {
            "reservoirs": [{
                "name": "Lake",
                "storage": {"min": 1_000_000.0, "max": 4_000_000.0},
                "initial_storage": 3_000_000.0,
                "level": {"reference_volume": 0.0, "elevation": 40.0,
                          "slope": 0.000020, "curvature": 0.0},
                "water_value": 0.001,
            }],
            "boundaries": [{"name": "Sea", "head": 20.0}],
            "units": [{
                "name": "Station", "source": "Lake", "target": "Sea",
                "flow": {"min": 5.0, "max": 20.0},
                "power": {"min": 1.0, "max": 25.0},
                "head": {"min": 35.0, "max": 105.0},
                "efficiency": {"peak": 0.9},
            }],
        },
        "grid": [0.0, 1.0, 2.0],
        "prices": [-100.0, 150.0],
    }


def test_python_optimizes_commitment_and_conserves_water():
    case = watercourse()
    result = solve(case, time_limit=5.0)
    assert result["validation"]["feasible"]
    assert result["commitment"] == [[False, True]]
    assert result["unit_flow"][0] == pytest.approx([0.0, 20.0], abs=1e-5)
    for t, flow in enumerate(result["unit_flow"][0]):
        assert result["storage"][0][t + 1] == pytest.approx(
            result["storage"][0][t] - 3600.0 * flow, abs=0.01
        )
    assert case == watercourse()


def test_python_rejects_invalid_input_without_starting_a_solver():
    case = watercourse()
    case["system"]["reservoirs"][0]["initial_storage"] = -1.0
    with pytest.raises(ValueError, match="initial state"):
        solve(case)
    with pytest.raises(ValueError, match="unknown field"):
        solve({**watercourse(), "misspelled_prices": []})


def test_independent_python_calls_own_their_input_and_solver_state():
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: solve(watercourse(), time_limit=5.0), range(2)))
    assert all(result["validation"]["feasible"] for result in results)
    assert results[0]["objective"] == pytest.approx(results[1]["objective"], abs=0.01)


def test_connected_watercourse_conserves_reservoir_and_river_water():
    path = Path(__file__).resolve().parents[2] / "examples/rust/cascade.json"
    case = json.loads(path.read_text())
    result = solve(case, time_limit=10.0)
    assert result["validation"]["feasible"]
    for t, (start, end) in enumerate(zip(case["grid"], case["grid"][1:])):
        seconds = 3600.0 * (end - start)
        stored = sum(row[t + 1] - row[t] for row in result["storage"])
        in_transit = sum(row[t + 1] - row[t] for row in result["river_storage"])
        inflow = sum(row[t] for row in case["inflows"].values())
        discharge_to_sea = result["unit_flow"][2][t]
        assert stored + in_transit == pytest.approx(
            seconds * (inflow - discharge_to_sea), abs=0.01
        )


@pytest.mark.parametrize("seconds", [0.0, float("nan"), float("inf"), 1e100])
def test_invalid_time_allowance_is_an_input_error(seconds):
    with pytest.raises(ValueError, match="time"):
        solve(watercourse(), time_limit=seconds)
