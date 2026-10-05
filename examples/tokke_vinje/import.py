#!/usr/bin/env python3
"""Build an explicitly scoped OpenSHOP benchmark from SINTEF's public data.

This independent importer reads facts and literal notebook configuration; it
never imports or executes upstream Python/pySHOP. Upstream files stay external.
The resulting generation/hydraulic profile is NOT a complete SHOP case import.
See the generated metadata for every restriction, policy and excluded attribute.
"""
from __future__ import annotations

import argparse
import ast
import bisect
import collections
import copy
import csv
import datetime as dt
import hashlib
import json
import math
import pathlib
import re
import subprocess

SOURCE_URL = "https://gitlab.sintef.no/energy/open-modelling-tools/open-datasets/tokke-vinje-watercourse"
PINNED_COMMIT = "ba2f2fc1f95d18978a04dd2c658aeef79126981b"
UTC = dt.timezone.utc


def yaml_data(path):
    try:
        import yaml
    except ImportError:
        # Ruby/Psych is part of the host's standard runtime; no YAML parser is
        # reimplemented and no scientific Python environment is required.
        script = "puts JSON.generate(YAML.safe_load(File.read(ARGV[0]), permitted_classes: [], aliases: false))"
        return json.loads(
            subprocess.check_output(
                ["ruby", "-ryaml", "-rjson", "-e", script, str(path)]
            )
        )
    return yaml.safe_load(path.read_text())


def timestamp(value):
    parsed = dt.datetime.fromisoformat(value.strip().replace("Z", "+00:00"))
    return (
        parsed.replace(tzinfo=UTC) if parsed.tzinfo is None else parsed.astimezone(UTC)
    )


def linear(xs, ys, value):
    i = max(0, min(len(xs) - 2, bisect.bisect_right(xs, value) - 1))
    return ys[i] + (value - xs[i]) * (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])


def cubic_four_points(xs, ys, value):
    # Four points with not-a-knot end conditions determine a single cubic.
    # This independently reproduces the documented four-point spline samples.
    return sum(
        y * math.prod((value - xj) / (xi - xj) for j, xj in enumerate(xs) if j != i)
        for i, (xi, y) in enumerate(zip(xs, ys))
    )


def literal(node, env):
    """Accept data literals/arithmetic only, never arbitrary notebook code."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return copy.deepcopy(env[node.id])
    if isinstance(node, (ast.List, ast.Tuple)):
        value = [literal(x, env) for x in node.elts]
        return tuple(value) if isinstance(node, ast.Tuple) else value
    if isinstance(node, ast.Dict):
        return {
            literal(k, env): literal(v, env) for k, v in zip(node.keys, node.values)
        }
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = literal(node.operand, env)
        return -value if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.BinOp):
        a, b = literal(node.left, env), literal(node.right, env)
        if isinstance(node.op, ast.Add):
            return a + b
        if isinstance(node.op, ast.Sub):
            return a - b
        if isinstance(node.op, ast.Mult):
            return a * b
        if isinstance(node.op, ast.Div):
            return a / b
    raise ValueError("not literal configuration")


def notebook_nodes(path):
    for block in re.findall(r"```python\s*\n(.*?)```", path.read_text(), flags=re.S):
        for node in ast.parse(block).body:
            yield node


def deep_update(target, updates):
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            deep_update(target[key], value)
        elif isinstance(value, dict) and isinstance(target.get(key), list):
            deep_update(target[key][0], value)
        else:
            target[key] = copy.deepcopy(value)


def topology(source):
    data = yaml_data(source / "data/input_data/SEnDHub/draft_topology.yaml")
    original_counts = {k: len(v) for k, v in data["model"].items()}
    env = {}
    generator_parameters = {}
    for node in notebook_nodes(source / "Tokke_topology_manual_modification.md"):
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
        ):
            try:
                env[node.targets[0].id] = literal(node.value, env)
            except (ValueError, KeyError, TypeError):
                # The only whitelisted calculation call produces numerical
                # efficiency samples from literal public source measurements.
                call = node.value
                if (
                    isinstance(call, ast.Call)
                    and isinstance(call.func, ast.Attribute)
                    and call.func.attr == "create_base_model_turb_efficiency"
                ):
                    xs, ys = [literal(x, env) for x in call.args[:2]]
                    kw = {x.arg: literal(x.value, env) for x in call.keywords}
                    xx = [x / xs[-1] * 100 for x in xs]
                    percentages = list(
                        range(
                            kw["start_percent"], kw["end_percent"] + 1, kw["interval"]
                        )
                    )
                    env[node.targets[0].id] = {
                        "x": percentages,
                        "y": [cubic_four_points(xx, ys, x) for x in percentages],
                    }
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        if not isinstance(call.func, ast.Attribute):
            continue
        operation = call.func.attr
        if (
            isinstance(call.func.value, ast.Name)
            and call.func.value.id == "create_final_topology"
        ):
            kw = {}
            for key in call.keywords:
                try:
                    kw[key.arg] = literal(key.value, env)
                except (ValueError, KeyError, TypeError):
                    pass
            if operation == "remove_objects_from_yaml":
                removed = kw["objects_to_remove"]
                for objects in data["model"].values():
                    for name in list(objects):
                        if any(rem in name for rem in removed):
                            del objects[name]
                data["connections"] = [
                    c
                    for c in data["connections"]
                    if not any(rem in c["from"] or rem in c["to"] for rem in removed)
                ]
            elif operation == "split_reservoirs_from_yaml":
                for instruction in kw["split_instructions"]:
                    name = instruction["reservoir_name"]
                    old = data["model"]["reservoir"].pop(name)
                    for new, fraction, lo, hi in zip(
                        instruction["new_names"],
                        instruction["percentages"],
                        instruction["lrls"],
                        instruction["hrls"],
                    ):
                        obj = copy.deepcopy(old)
                        obj.update(
                            max_vol=round(old["max_vol"] * fraction, 1), lrl=lo, hrl=hi
                        )
                        curve = obj.setdefault("vol_head", {"ref": 0})
                        if (
                            old.get("vol_head")
                            and lo == old.get("lrl")
                            and hi == old.get("hrl")
                        ):
                            curve["x"] = [round(x * fraction, 1) for x in curve["x"]]
                        else:
                            original_x = old.get("vol_head", {}).get(
                                "x", [0.0, old["max_vol"]]
                            )
                            curve.update(
                                x=[
                                    round(original_x[0] * fraction, 1),
                                    round(original_x[-1] * fraction, 1),
                                ],
                                y=[lo, hi],
                            )
                        data["model"]["reservoir"][new] = obj
                    data["connections"] = [
                        c
                        for c in data["connections"]
                        if not (
                            c["from_type"] == "reservoir"
                            and c["from"] == name
                            or c["to_type"] == "reservoir"
                            and c["to"] == name
                        )
                    ]
            elif operation == "replace_names_in_yaml":
                mappings = {
                    "reservoir": kw["reservoir_mapping"],
                    "plant": kw["plant_mapping"],
                    "river": kw["river_mapping"],
                }

                def renamed(kind, name):
                    if kind == "generator":
                        for old, new in mappings["plant"].items():
                            name = name.replace(old, new)
                    elif kind == "river":
                        name = mappings[kind].get(name, name)
                        for old, new in mappings["reservoir"].items():
                            name = name.replace(old, new)
                    else:
                        name = mappings.get(kind, {}).get(name, name)
                    return name

                for kind, objects in data["model"].items():
                    data["model"][kind] = {
                        renamed(kind, name): value for name, value in objects.items()
                    }
                for connection in data["connections"]:
                    for side in ("from", "to"):
                        connection[side] = renamed(
                            connection[side + "_type"], connection[side]
                        )
            elif operation == "update_attributes_in_yaml":
                for update in kw["updates"]:
                    obj = (
                        data["model"]
                        .setdefault(update["type"], {})
                        .setdefault(update["name"], {})
                    )
                    deep_update(obj, update["attributes"])
            elif operation == "apply_reservoir_heuristics_in_yaml":
                for obj in data["model"]["reservoir"].values():
                    curve = obj["vol_head"]
                    x0, x1 = curve["x"][-2:]
                    y0, y1 = curve["y"][-2:]
                    extra = 0.1 * x1
                    curve["x"].append(round(x1 + extra, 2 if extra > 0.01 else 4))
                    curve["y"].append(
                        round(y1 + 0.99 * (y1 - y0) / (x1 - x0) * extra, 2)
                    )
            elif operation == "apply_river_heuristics_in_yaml":
                base = kw["base_model_flow_cost"]
                x0, x1 = base["x"][-2:]
                y0, y1 = base["y"][-2:]
                curvature = (y1 - y0) / (x1 - x0) ** 2
                edges = kw["q_intervals"]
                marginal = [1.0] + [
                    round(curvature * (a + b - 2 * x0))
                    for a, b in zip(edges[1:-1], edges[2:])
                ]
                for name in kw["river_names"]:
                    obj = data["model"]["river"][name]
                    obj.pop("max_flow_const", None)
                    obj.pop("flow_cost_const", None)
                    obj["flow_cost_curve"] = {
                        "ref": 0,
                        "x": edges,
                        "y": marginal + [marginal[-1]],
                    }
            elif operation == "apply_spillway_heuristics_in_yaml":
                for name, crest in kw["crest_lengths"].items():
                    curve = data["model"]["reservoir"][name]["vol_head"]
                    x0, x1 = curve["x"][-2:]
                    y0, y1 = curve["y"][-2:]
                    lengths = crest if isinstance(crest, list) else [crest]
                    suffixes = (
                        ["_" + x for x in kw["downstream_mapping"][name]]
                        if isinstance(crest, list)
                        else [""]
                    )
                    for width, suffix in zip(lengths, suffixes):
                        discharge = round(3 * width * (y1 - y0) ** 1.5, 2)
                        discharge = max(discharge, round((x1 - x0) / 0.0864, 2))
                        data["model"]["river"]["s_" + name + suffix][
                            "up_head_flow_curve"
                        ] = [{"ref": 0, "x": [y0, y1], "y": [0.0, discharge]}]
            elif operation == "apply_generator_heuristics_in_yaml":
                generator_parameters = {
                    x["generator_name"]: x for x in kw["generators_attributes"]
                }
                for name, params in generator_parameters.items():
                    obj = data["model"]["generator"][name]
                    nominal = obj["turb_eff_curves"][0]["ref"]
                    qmax = obj["turb_eff_curves"][0]["x"][-1]
                    base = kw["efficiency_types"][params["efficiency_type"]]
                    age = kw["current_year"] - params["install_year"]
                    degradation = params["recovery_factor"] ** (
                        age // params["refurbishment_interval"]
                    ) * math.exp(-0.001 * (age % params["refurbishment_interval"]))
                    flows = [qmax * x / 100 for x in base["x"]]
                    efficiency = [y * degradation for y in base["y"]]
                    curves = [
                        {
                            "ref": nominal,
                            "x": [round(x, 2) for x in flows],
                            "y": [round(y, 2) for y in efficiency],
                        }
                    ]
                    lo, hi = params["min_head"], params["max_head"]
                    if lo != nominal:
                        start = 1 if params["efficiency_type"] == 1 else 0
                        curves.append(
                            {
                                "ref": lo,
                                "x": [round(x, 2) for x in flows[start:]],
                                "y": [round(0.95 * y, 2) for y in efficiency[start:]],
                            }
                        )
                    if hi != nominal:
                        scale = (
                            1.0
                            if nominal == lo
                            else 1 - (hi - nominal) / (nominal - lo) * 0.05
                        )
                        curves.append(
                            {
                                "ref": hi,
                                "x": [round(x, 2) for x in flows],
                                "y": [round(scale * y, 2) for y in efficiency],
                            }
                        )
                    obj["turb_eff_curves"] = sorted(curves, key=lambda c: c["ref"])
                    obj["p_min"] = round(
                        obj["p_max"]
                        * base["y"][0]
                        * flows[0]
                        / (base["y"][-1] * flows[-1]),
                        1,
                    )
                    obj["gen_eff_curve"]["x"][0] = obj["p_min"]
        elif operation == "connect_to":

            def endpoint(expression):
                if not isinstance(expression, ast.Subscript) or not isinstance(
                    expression.value, ast.Attribute
                ):
                    raise ValueError("unsupported literal connection")
                return expression.value.attr, literal(expression.slice, env)

            a_kind, a_name = endpoint(call.func.value)
            b_kind, b_name = endpoint(call.args[0])
            data["connections"].append(
                {
                    "from": a_name,
                    "from_type": a_kind,
                    "to": b_name,
                    "to_type": b_kind,
                    "connection_type": "connection_standard",
                    "order": 0,
                }
            )
    # Repeated notebook connections denote one graph edge; no duplicate flow.
    unique = {}
    for connection in data["connections"]:
        unique[tuple(connection[k] for k in ("from_type", "from", "to_type", "to"))] = (
            connection
        )
    data["connections"] = list(unique.values())
    expected = {"reservoir": 17, "plant": 9, "generator": 14, "river": 37, "tunnel": 5}
    reconstructed = {k: len(v) for k, v in data["model"].items()}
    if reconstructed != expected:
        raise ValueError(f"reviewed topology inventory changed: {reconstructed}")
    return data, generator_parameters, original_counts


def level_samples(source, at):
    result = {}
    evidence = {}
    for path in (source / "data/input_data/Sildre/reservoir_level").glob("*.csv"):
        with path.open(encoding="utf-8-sig") as stream:
            name = stream.readline().lstrip("#").split(",")[0].split(" ")[0]
            name = name.translate(
                str.maketrans(
                    {"æ": "ae", "ø": "oe", "å": "aa", "Æ": "Ae", "Ø": "Oe", "Å": "Aa"}
                )
            )
            if name == "Kjela":
                name = "Kjelavatn"
            before = after = None
            for row in csv.DictReader(stream, delimiter=";"):
                value = row.get("Vannstand (m)", "").replace(",", ".")
                if not value:
                    continue
                time = timestamp(row["Tidspunkt"])
                sample = (time, float(value))
                if time <= at and (before is None or time > before[0]):
                    before = sample
                if time >= at and (after is None or time < after[0]):
                    after = sample
            if (
                before is None
                or after is None
                or (after[0] - before[0]).total_seconds() > 31 * 86400
            ):
                continue
            value = (
                before[1]
                if before[0] == after[0]
                else before[1]
                + (after[1] - before[1])
                * (at - before[0]).total_seconds()
                / (after[0] - before[0]).total_seconds()
            )
            result[name] = round(value, 2)
            evidence[name] = {
                "file": str(path.relative_to(source)),
                "bracket": [x[0].isoformat() for x in (before, after)],
                "interpolated_level_m": result[name],
            }
    result["Hyljelihyl"] = 704.2
    result["Vatjern"] = 836.5
    result["Vaamarvatn"] = round(
        677.0 + (result["Totak"] - 680.0) / (687.3 - 680.0) * (687.3 - 677.0), 2
    )
    evidence.update(
        {
            "Hyljelihyl": {"policy": "upstream documented constant 704.2 m"},
            "Vatjern": {"policy": "upstream documented constant 836.5 m"},
            "Vaamarvatn": {
                "policy": "upstream documented proportional Totak extension"
            },
        }
    )
    return result, evidence


def historical_inputs(source, at, hours):
    env = {}
    for node in notebook_nodes(source / "Tokke_time_series_manual_modification.md"):
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "reservoir_mapping"
        ):
            env = literal(node.value, {})
            break
    if len(env) != 18:
        raise ValueError("unexpected public reservoir inflow mapping")
    end = at + dt.timedelta(hours=max(hours, 168))
    prices = {}
    with (
        source / "data/input_data/SEnDHub/Historical_day_ahead_price.csv"
    ).open() as stream:
        for row in csv.DictReader(stream):
            time = timestamp(row["date"])
            if at <= time < end:
                prices.setdefault(time, float(row["Day-ahead price (EUR/MWh)"]))
    inflows = {}
    factors = {}
    for station in sorted({s for _, s in env.values()}):
        sums = collections.defaultdict(float)
        hourly = {}
        path = source / f"data/input_data/SEnDHub/{station}_combined_DETD_Sildre.csv"
        with path.open() as stream:
            for row in csv.DictReader(stream):
                time = timestamp(row["date"])
                if 1991 <= time.year <= 2020:
                    sums[time.year] += 0.0036 * float(row["DETD inflow (m3/s)"] or 0.0)
                if (
                    at <= time < at + dt.timedelta(hours=hours)
                    and row["Combined inflow (m3/s)"]
                ):
                    hourly[time] = float(row["Combined inflow (m3/s)"])
        mean = sum(sums.values()) / len(sums)
        factors[station] = {
            "mean_annual_DETD_Mm3_1991_2020": mean,
            "source_file": str(path.relative_to(source)),
        }
        for name, (annual, mapped_station) in env.items():
            if mapped_station == station:
                inflows[name] = [
                    round(annual / mean * hourly[at + dt.timedelta(hours=t)], 2)
                    for t in range(hours)
                ]
    selected_prices = [prices[at + dt.timedelta(hours=t)] for t in range(hours)]
    next_week = [prices[at + dt.timedelta(hours=t)] for t in range(168)]
    if any(x < 0 for values in inflows.values() for x in values):
        raise ValueError("negative historical inflow")
    return selected_prices, inflows, factors, sum(next_week) / 168


def compile_case(data, params, levels, prices, inflows, terminal_price, hours, source):
    model = data["model"]
    edges = data["connections"]
    case = {
        "schema_version": 2,
        "name": f"tokke_vinje_hydraulic_generation_{hours}h",
        "grid": list(map(float, range(hours + 1))),
        "prices": prices[:hours],
        "reservoirs": [],
        "junctions": [],
        "boundaries": [],
        "tunnels": [],
        "plants": [],
        "generators": [],
        "river_junctions": [],
        "rivers": [],
        "operations": [],
    }
    origins = {
        "Bandak": "HOGGA",
        "Vinjevatn": "TOKKE",
        "Byrtevatn": "LIO",
        "Botnedalsvatn": "BYRTE",
        "Totak": "VINJE",
        "Vaamarvatn": "VINJE",
        "Langeidvatn": "HAUKELI",
        "Vatjern": "HAUKELI",
        "Songavatn": "SONGA+BITDAL",
        "Bitdalsvatn": "SONGA+BITDAL",
        "Bordalsvatn": "KJELA",
        "Foersvatn": "KJELA",
        "Kjelavatn": "VESLE KJELA",
        "Staavatn": "STAAVATN",
        "Venemo": "VENEMO",
        "Langesae": "LANGESAE",
        "Hyljelihyl": "HYLJELIHYL",
    }
    with (source / "data/input_data/SEnDHub/moduledata.csv").open() as stream:
        energy = {
            r["modnavn"]: float(r["enekv_global"]) * 1000
            for r in csv.DictReader(stream)
        }
    # A restricted pressurized-tunnel profile avoids unsupported emergence.
    minimum_wet_head = collections.defaultdict(float)
    for name, obj in model["tunnel"].items():
        for side, height in (("from", obj["start_height"]), ("to", obj["end_height"])):
            matching = [
                e
                for e in edges
                if e["to" if side == "from" else "from"] == name
                and e["to_type" if side == "from" else "from_type"] == "tunnel"
                and e[side + "_type"] == "reservoir"
            ]
            for edge in matching:
                minimum_wet_head[edge[side]] = max(minimum_wet_head[edge[side]], height)
    restrictions = []
    for name, obj in model["reservoir"].items():
        curve = obj["vol_head"]
        initial = linear(curve["y"], curve["x"], levels[name])
        lower_head = max(obj["lrl"], minimum_wet_head[name])
        vmin = linear(curve["y"], curve["x"], lower_head)
        vmax = min(obj["max_vol"], linear(curve["y"], curve["x"], obj["hrl"]))
        if not vmin <= initial <= vmax:
            raise ValueError(
                f"{name}: historical storage {initial} outside restricted [{vmin},{vmax}]"
            )
        if lower_head > obj["lrl"]:
            restrictions.append(
                {
                    "object": name,
                    "attribute": "vmin",
                    "source_lrl_m": obj["lrl"],
                    "restricted_min_head_m": lower_head,
                    "reason": "keep documented tunnel mouths submerged; emergence is not modeled",
                }
            )
        rates = inflows[name][:hours]
        if name == "Bandak":
            rates = [a + b for a, b in zip(rates, inflows["r_Vest_Vassdraget"][:hours])]
        case["reservoirs"].append(
            {
                "name": name,
                "v0": initial,
                "vmin": vmin,
                "vmax": vmax,
                "inflow": rates[0],
                "water_value": terminal_price * energy[origins[name]],
                "level_curve": {"x": curve["x"], "y": curve["y"]},
            }
        )
        case["operations"].append(
            {
                "object": name,
                "attribute": "inflow",
                "times": list(map(float, range(hours))),
                "values": rates,
            }
        )
    reservoirs = {r["name"]: r for r in case["reservoirs"]}
    node_head_bounds = {
        name: tuple(
            linear(r["level_curve"]["x"], r["level_curve"]["y"], v)
            for v in (r["vmin"], r["vmax"])
        )
        for name, r in reservoirs.items()
    }
    intake_bound_provenance = []

    def junction(name, lo=0.0, hi=1100.0):
        if not any(j["name"] == name for j in case["junctions"]):
            case["junctions"].append({"name": name, "hmin": lo, "hmax": hi})
        node_head_bounds[name] = (lo, hi)
        return name

    def loss_intake(name, source_node, resistance, cap):
        source_lo, source_hi = node_head_bounds[source_node]
        lo = source_lo - resistance * cap**2
        junction(name, lo, source_hi)
        intake_bound_provenance.append(
            {
                "node": name,
                "source_node": source_node,
                "source_head_bounds_m": [source_lo, source_hi],
                "aggregate_discharge_upper_m3_s": cap,
                "loss_factor": resistance,
                "head_bounds_m": [lo, source_hi],
                "formula": "H_intake = H_source - k*Q^2; 0 <= Q <= sum(unit qmax)",
            }
        )
        return name

    def tunnel(name, a, b, loss, cap):
        case["tunnels"].append(
            {
                "name": name,
                "source": a,
                "target": b,
                "resistance": loss,
                "capacity": cap,
                "opening": 1.0,
            }
        )

    group_map = {}
    for name, plant in model["plant"].items():
        units = [
            n
            for n in model["generator"]
            if any(
                e["from"] == n and e["to"] == name and e["from_type"] == "generator"
                for e in edges
            )
        ]
        inputs = [
            e
            for e in edges
            if e["to"] == name
            and e["to_type"] == "plant"
            and e["from_type"] in ("reservoir", "tunnel")
        ]
        if len(inputs) == 1 and inputs[0]["from_type"] == "reservoir":
            upstream = inputs[0]["from"]
        elif inputs and all(e["from_type"] == "tunnel" for e in inputs):
            feeding = [
                e["from"]
                for e in edges
                if e["from_type"] == "reservoir"
                and e["to_type"] == "tunnel"
                and e["to"] in {x["from"] for x in inputs}
            ]
            # The feeding tunnels' outlet mouths must stay submerged in this
            # declared profile. A flowing merge can sit below its source heads.
            floor = max(model["tunnel"][x["from"]]["end_height"] for x in inputs)
            ceiling = max(node_head_bounds[n][1] for n in feeding)
            upstream = junction("__" + name + "_tunnel_merge", floor, ceiling)
            intake_bound_provenance.append(
                {
                    "node": upstream,
                    "source_reservoirs": feeding,
                    "head_bounds_m": [floor, ceiling],
                    "formula": "lower = maximum feeding outlet-mouth elevation (pressurized profile); upper = maximum supplying reservoir head; nonnegative aggregate generation outflow",
                }
            )
        else:
            raise ValueError(f"ambiguous plant input: {name}")
        supply_reservoirs = (
            [upstream]
            if upstream in reservoirs
            else [
                e["from"]
                for e in edges
                if e["from_type"] == "reservoir"
                and e["to_type"] == "tunnel"
                and e["to"] in {x["from"] for x in inputs}
            ]
        )
        if not supply_reservoirs:
            raise ValueError(f"missing reservoir head bound: {name}")
        source_head_max = max(
            linear(
                reservoirs[n]["level_curve"]["x"],
                reservoirs[n]["level_curve"]["y"],
                reservoirs[n]["vmax"],
            )
            for n in supply_reservoirs
        )
        outputs = [
            e
            for e in edges
            if e["from"] == name
            and e["from_type"] == "plant"
            and e["to_type"] in ("reservoir", "river")
        ]
        if len(outputs) == 1:
            downstream = (
                "Bandak"
                if outputs[0]["to"] == "r_Vest_Vassdraget"
                else outputs[0]["to"]
            )
        elif name == "Hogga" and not outputs:
            downstream = "__Hogga_tail_boundary"
            case["boundaries"].append(
                {"name": downstream, "head": plant["outlet_line"]}
            )
        else:
            raise ValueError(f"ambiguous plant output: {name}")
        if downstream in reservoirs:
            receiver = reservoirs[downstream]
            tail_head_min = linear(
                receiver["level_curve"]["x"],
                receiver["level_curve"]["y"],
                receiver["vmin"],
            )
        else:
            tail_head_min = plant["outlet_line"]
        # Positive generation flow only loses head through shared waterways.
        # The upstream HRL and receiver minimum therefore give a conservative
        # physical upper bound, independently of turbine reference heads.
        operating_head_max = source_head_max - max(tail_head_min, plant["outlet_line"])
        if operating_head_max <= 0:
            raise ValueError(f"nonpositive physical head domain: {name}")
        total = sum(
            max(c["x"][-1] for c in model["generator"][u]["turb_eff_curves"])
            for u in units
        )
        loss = plant["main_loss"]
        if len(loss) != 1 or loss[0] <= 0:
            raise ValueError(f"unsupported shared main loss: {name}")
        intake = loss_intake("__" + name + "_main_intake", upstream, loss[0], total)
        tunnel("__" + name + "_main_loss", upstream, intake, loss[0], total)
        groups = collections.defaultdict(list)
        for unit in units:
            groups[model["generator"][unit]["penstock"]].append(unit)
        group_map[name] = []
        for number, members in sorted(groups.items()):
            coefficient = (
                plant["penstock_loss"][number - 1]
                if number <= len(plant["penstock_loss"])
                else 0.0
            )
            group_name = name if len(groups) == 1 else name + f"_penstock_{number}"
            group_map[name].append(group_name)
            node = intake
            if coefficient > 0:
                cap = sum(
                    max(c["x"][-1] for c in model["generator"][u]["turb_eff_curves"])
                    for u in members
                )
                node = loss_intake(
                    "__" + group_name + "_intake", intake, coefficient, cap
                )
                tunnel("__" + group_name + "_loss", intake, node, coefficient, cap)
            case["plants"].append(
                {
                    "name": group_name,
                    "source": node,
                    "target": downstream,
                    "pmax": sum(model["generator"][u]["p_max"] for u in members),
                    "ramp": 1e12,
                    "outlet_head_floor": plant["outlet_line"],
                }
            )
            for unit in members:
                obj = model["generator"][unit]
                curves = obj["turb_eff_curves"]
                discharges = sorted({x for curve in curves for x in curve["x"]})
                efficiencies = [
                    [linear(c["x"], c["y"], q) / 100 for c in curves]
                    for q in discharges
                ]
                if not all(0 < eta <= 1 for row in efficiencies for eta in row):
                    raise ValueError(f"invalid reconstructed efficiency {unit}")
                electrical = obj["gen_eff_curve"]
                ex = list(electrical["x"])
                ey = [y / 100 for y in electrical["y"]]
                if ex[0] > 0:
                    ey.insert(0, linear(ex, ey, 0.0))
                    ex.insert(0, 0.0)
                case["generators"].append(
                    {
                        "name": unit,
                        "plant": group_name,
                        "qmin": discharges[0],
                        "qmax": discharges[-1],
                        "pmin": obj["p_min"],
                        "pmax": obj["p_max"],
                        "efficiency": max(max(row) for row in efficiencies),
                        "min_efficiency": 0.0,
                        "qbest": 0.85 * discharges[-1],
                        "qcurvature": 0.0,
                        "hbest": model["generator"][unit]["turb_eff_curves"][
                            1 if len(curves) > 2 else 0
                        ]["ref"],
                        "hcurvature": 0.0,
                        "hmin": 1e-3,
                        "hmax": operating_head_max,
                        "initial_on": 0,
                        "initial_age": 1000.0,
                        "minup": 0.0,
                        "mindown": 0.0,
                        "startup": obj["startcost_const"],
                        "shutdown": obj["stopcost_const"],
                        "turbine_table": {
                            "heads": [c["ref"] for c in curves],
                            "discharge": discharges,
                            "efficiency": efficiencies,
                            "qmin": [c["x"][0] for c in curves],
                            "qmax": [c["x"][-1] for c in curves],
                            "interpolation": "bilinear",
                            "head_extrapolation": "linear",
                        },
                        "generator_efficiency_curve": {"x": ex, "y": ey},
                    }
                )
    for name, obj in model["tunnel"].items():
        incoming = [e for e in edges if e["to"] == name and e["to_type"] == "tunnel"]
        outgoing = [
            e for e in edges if e["from"] == name and e["from_type"] == "tunnel"
        ]
        if len(incoming) != 1 or len(outgoing) != 1:
            raise ValueError(f"ambiguous tunnel: {name}")
        a = incoming[0]["from"]
        b = (
            outgoing[0]["to"]
            if outgoing[0]["to_type"] == "reservoir"
            else "__" + outgoing[0]["to"] + "_tunnel_merge"
        )

        def head_range(node):
            if node in reservoirs:
                r = reservoirs[node]
                c = r["level_curve"]
                return linear(c["x"], c["y"], r["vmin"]), linear(
                    c["x"], c["y"], r["vmax"]
                )
            j = next(j for j in case["junctions"] if j["name"] == node)
            return j["hmin"], j["hmax"]

        alo, ahi = head_range(a)
        blo, bhi = head_range(b)
        cap = math.sqrt(max(abs(ahi - blo), abs(bhi - alo)) / obj["loss_factor"]) + 1e-6
        tunnel(name, a, b, obj["loss_factor"], cap)
    sink_name = "__downstream_river_boundary"
    river_meta = []
    for name, obj in model["river"].items():
        if name == "r_Vest_Vassdraget":
            continue
        incoming = [e for e in edges if e["to"] == name and e["to_type"] == "river"]
        outgoing = [e for e in edges if e["from"] == name and e["from_type"] == "river"]
        if len(incoming) != 1 or incoming[0]["from_type"] != "reservoir":
            raise ValueError(f"unsupported river source: {name}")
        src = incoming[0]["from"]
        target = outgoing[0]["to"] if len(outgoing) == 1 else sink_name
        if target == "r_Vest_Vassdraget":
            target = "Bandak"
        if target not in reservoirs and target != sink_name:
            raise ValueError(f"unsupported river target: {name}")
        discharge_curve = obj.get("up_head_flow_curve")
        if discharge_curve:
            dc = discharge_curve[0]
            r = reservoirs[src]
            c = r["level_curve"]
            low = linear(c["x"], c["y"], r["vmin"])
            xx = list(dc["x"])
            yy = list(dc["y"])
            if low < xx[0]:
                xx.insert(0, low)
                yy.insert(0, 0.0)
            cap = max(1.0, max(yy))
            law = "weir"
            curve = {"x": xx, "y": yy}
        else:
            cap = max(1.0, obj.get("max_flow_const", 1000.0))
            law = "controlled"
            curve = None
        river = {
            "name": name,
            "source": src,
            "target": target,
            "curves": [],
            "capacity": cap,
            "law": law,
            "deterministic_delay": 0.0,
            "water_value": 0.0,
            "history_grid": [-1.0, 0.0],
            "history_release": [0.0],
            "discharge_curve": curve,
        }
        case["rivers"].append(river)
        if obj.get("max_flow_const") == 0:
            case["operations"].append(
                {
                    "object": name,
                    "attribute": "gate_max",
                    "times": [0.0],
                    "values": [0.0],
                }
            )
        river_meta.append(
            {
                "name": name,
                "source": src,
                "target": target,
                "source_attributes": list(obj),
                "flow_costs_excluded": any(
                    k in obj for k in ("flow_cost_const", "flow_cost_curve")
                ),
                "source_soft_max_treated_as_hard": "max_flow_const" in obj,
            }
        )
    if any(r["target"] == sink_name for r in case["rivers"]):
        # River arrival is delivered out of the modeled watercourse. This head
        # is unused by the supplied upstream-only river discharge laws.
        case["boundaries"].append(
            {"name": sink_name, "head": model["plant"]["Hogga"]["outlet_line"]}
        )
    return case, {
        "pressurized_domain_restrictions": restrictions,
        "hydraulic_head_bound_provenance": intake_bound_provenance,
        "original_to_compiled_plants": group_map,
        "river_coverage": river_meta,
        "downstream_river_boundary": sink_name,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = pathlib.Path(__file__).resolve().parent
    parser.add_argument(
        "--source", type=pathlib.Path, default=root / "data/tokke-vinje-watercourse"
    )
    parser.add_argument("--output", type=pathlib.Path, default=root / "generated")
    parser.add_argument("--start", default="2024-09-01T00:00:00+00:00")
    parser.add_argument("--hours", type=int, nargs="+", default=[2, 6, 24])
    args = parser.parse_args()
    at = timestamp(args.start)
    commit = subprocess.check_output(
        ["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True
    ).strip()
    if commit != PINNED_COMMIT:
        raise ValueError(
            f"expected reviewed source commit {PINNED_COMMIT}, got {commit}"
        )
    unchanged = subprocess.run(
        [
            "git",
            "-C",
            str(args.source),
            "-c",
            "filter.lfs.smudge=",
            "-c",
            "filter.lfs.process=",
            "-c",
            "filter.lfs.required=false",
            "diff",
            "--quiet",
            "HEAD",
            "--",
        ]
    ).returncode
    if unchanged:
        raise ValueError("reviewed source checkout has modified tracked files")
    data, params, original = topology(args.source)
    levels, level_provenance = level_samples(args.source, at)
    prices, inflows, stations, terminal_price = historical_inputs(
        args.source, at, max(args.hours)
    )
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = {
        "source_url": SOURCE_URL,
        "source_commit": commit,
        "profile": "restricted hydraulic-generation benchmark; not full SHOP equivalence",
        "source_license": "no explicit LICENSE/COPYING/NOTICE found; external data not bundled with published library",
        "start_utc": at.isoformat(),
        "source_draft_inventory": original,
        "documented_final_inventory": {k: len(v) for k, v in data["model"].items()},
        "historical_level_provenance": level_provenance,
        "inflow_scaling": stations,
        "terminal_policy": {
            "method": "next 168-hour mean historical price × NVE DETD enekv_global",
            "forward_mean_price_EUR_MWh": terminal_price,
            "unit_conversion": "kWh/m3 × 1000 = MWh/Mm3",
            "equivalence": "explicit benchmark policy; not SHOP's internal energy_value_input conversion or the upstream filling-percentile policy",
        },
        "declarations": [
            "All 17 documented reservoirs, 14 physical generating units and 5 documented tunnels are retained.",
            "Turbine efficiency curves are reconstructed from the upstream documented heuristic; they are not measured plant efficiency data.",
            "Supplied turbine knots are preserved on a union discharge grid using piecewise-linear discharge/head interpolation; SHOP spline semantics are not claimed.",
            "All source main losses and shared penstocks are compiled into exact quadratic tunnel/junction blocks; per-penstock plant groups preserve the sum of original unit power limits.",
            "Source outlet_line is represented by the exact outlet_head_floor field of OpenSHOP.",
            "Startup and shutdown costs are preserved; initial commitment is an explicit all-off benchmark policy with no dwell restrictions because source history is absent.",
            "No physical travel-time attribute is supplied; zero delay is used, not an invented calibrated delay.",
            "The zero-delay r_Vest_Vassdraget mixing reach is eliminated exactly for the declared profile; Lio output and four bypass/spill arrivals go to Bandak, and local inflow is added to Bandak.",
            "Environmental constraints, reserves, river flow costs and source soft-bound penalty semantics are excluded from this generation/hydraulic profile.",
            "Nominal reservoir volume bounds and source river maximum flows are enforced as hard benchmark restrictions; spillway tables remain but flooding beyond HRL is not allowed.",
            "Storage lower bounds are tightened where needed to keep documented inter-reservoir tunnel mouths submerged; unsupported emergence is excluded, not approximated.",
            "No artificial minimum efficiency is imposed; all supplied physical low-efficiency points are retained.",
            "Historical water levels use documented <=31-day interpolation plus Hyljelihyl/Vatjern constants and the Totak-derived Vaamarvatn extension.",
            "The dataset's final output_data is gitignored; topology is independently reconstructed from literal public notebook configuration without licensed SHOP.",
            "Downstream rivers terminate at a fixed external Boundary; its head is unused by the supplied upstream-only river discharge laws.",
        ],
        "attribute_coverage": {
            "reservoir": {
                "max_vol": "hard nominal storage bound, restricted versus SHOP soft/flood behavior",
                "lrl": "minimum storage via inverse vol_head, tightened for submerged mouths where declared",
                "hrl": "hard nominal maximum via inverse vol_head",
                "vol_head": "preserved source knots plus documented extra flood point; historic-max extension not needed for selected start",
            },
            "plant": {
                "outlet_line": "exact outlet_head_floor",
                "main_loss": "exact additional shared quadratic tunnel",
                "penstock_loss": "exact additional shared-penstock tunnel and grouped physical units",
            },
            "generator": {
                "penstock": "exact unit membership in shared loss groups",
                "p_min": "preserved documented heuristic result",
                "p_max": "preserved",
                "p_nom": "equal to p_max in this dataset; metadata only",
                "startcost_const": "exact startup",
                "stopcost_const": "exact shutdown",
                "gen_eff_curve": "on-domain points preserved; linear numerical extension added at zero MW",
                "turb_eff_curves": "preserved reconstructed knots; bilinear interpolation policy explicitly selected",
            },
            "tunnel": {
                "length": "unused by steady-state pressure equation, metadata only",
                "start_height": "pressurized-domain restriction",
                "end_height": "pressurized-domain restriction",
                "loss_factor": "exact q|q| loss equation",
            },
            "river": {
                "upstream_elevation": "bypass-bed domain implicit in source reservoir LRL; spill curve defines crest",
                "max_flow_const": "treated as hard; source soft-bound semantics excluded",
                "flow_cost_const": "excluded",
                "flow_cost_curve": "excluded",
                "main_river": "metadata only; routing is zero delay",
                "up_head_flow_curve": "exact supplied piecewise-linear nonnegative discharge curve with zero below crest",
            },
        },
        "cases": [],
    }
    metadata["declarations"].append(
        "Turbine reference heads are not operating limits. Head bounds use positive head and a conservative upstream-HRL minus minimum outlet bound; values outside reference heads use explicit linear extrapolation, whose equivalence to SHOP interpolation is not claimed."
    )
    metadata["declarations"].append(
        "Shared main and penstock intake bounds propagate H_intake = H_source - k*Q^2 with 0 <= Q <= sum(unit qmax): lower = source_lower - k*Qmax^2, upper = source_upper. Exact equations and feasible physical domain are preserved; every interval and formula is recorded."
    )
    metadata["attribute_coverage"]["generator"][
        "turb_eff_curves"
    ] += "; linear extrapolation outside reference heads, positive physical operating-head domain"
    # Classification is exhaustive for the reconstructed static attributes.
    # A newly introduced upstream field must be reviewed before proceeding.
    metadata["source_object_attribute_inventory"] = {}
    for kind, objects in data["model"].items():
        metadata["source_object_attribute_inventory"][kind] = {
            name: sorted(obj) for name, obj in objects.items()
        }
        actual = {key for obj in objects.values() for key in obj}
        unknown = actual - metadata["attribute_coverage"][kind].keys()
        if unknown:
            raise ValueError(
                f"unclassified source {kind} attributes: {sorted(unknown)}"
            )
    metadata["time_series_coverage"] = {
        "Historical_day_ahead_price.csv": "exact dated hourly price samples",
        "combined_DETD_Sildre.csv": "documented station assignment and annual inflow scaling",
        "reservoir_level/*.csv": "dated level samples and explicit initial-state reconstruction",
        "environmental_constraints": "excluded: time-varying river minimum flow and minimum/maximum storage",
        "initial_unit_state": "absent from public source; explicit all-off benchmark policy",
        "terminal_water_value": "explicit next-week historical mean × public global energy-equivalent benchmark policy",
        "reserves": "excluded from this profile",
    }
    tracked = (
        subprocess.check_output(
            ["git", "-C", str(args.source), "ls-files", "-z"], text=True
        )
        .rstrip("\0")
        .split("\0")
    )
    metadata["source_file_inventory"] = [
        {"path": name, "checkout_bytes": (args.source / name).stat().st_size}
        for name in tracked
    ]
    for hours in args.hours:
        case, coverage = compile_case(
            data, params, levels, prices, inflows, terminal_price, hours, args.source
        )
        path = args.output / f"tokke_vinje_{hours}h.json"
        text = json.dumps(case, indent=2, allow_nan=False) + "\n"
        path.write_text(text)
        metadata["cases"].append(
            {
                "hours": hours,
                "file": str(path),
                "sha256": hashlib.sha256(text.encode()).hexdigest(),
                "compiled_inventory": {
                    k: len(case[k])
                    for k in (
                        "reservoirs",
                        "junctions",
                        "boundaries",
                        "tunnels",
                        "plants",
                        "generators",
                        "rivers",
                    )
                },
                **coverage,
            }
        )
    (args.output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {
                "start": metadata["start_utc"],
                "source_commit": commit,
                "cases": [
                    {
                        "hours": c["hours"],
                        "file": c["file"],
                        "inventory": c["compiled_inventory"],
                    }
                    for c in metadata["cases"]
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
