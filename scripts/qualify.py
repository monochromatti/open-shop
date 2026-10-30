#!/usr/bin/env python3
"""Compare scheduling budgets and operational recovery on a saved campaign plan."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import sys
from datetime import datetime, timezone

import fuzz


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path, help="JSON list of case requests; see docs/qualification.md")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--timeout", type=float, default=3600.0, help="hard seconds per case, including compilation")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    if not args.julia or not Path(args.julia).is_file():
        parser.error("use nix develop, or pass --julia /path/to/julia")
    if not 0 < args.timeout < float("inf"):
        parser.error("--timeout must be finite and positive")
    args.output = args.output.resolve()
    if args.output == fuzz.ROOT or fuzz.ROOT in args.output.parents:
        parser.error("keep run artifacts outside the package checkout")
    args.julia = str(Path(args.julia).resolve())
    requests = json.loads(args.plan.read_text())
    if not isinstance(requests, list) or not requests or any(not isinstance(r, dict) for r in requests):
        parser.error("plan must be a nonempty JSON list of case requests")
    for request in requests:
        if request.get("input"):
            request["input"] = str((args.plan.resolve().parent / request["input"]).resolve())
        request.setdefault("seed", "1")
        request["seed"] = str(request["seed"])
        request.setdefault("family", "external" if request.get("input") else "integrated")
        request.setdefault("size", "medium")
    source = fuzz.source_identity()
    digest = hashlib.sha256(source["source_sha256"].encode())
    for path in (Path(__file__), fuzz.ROOT / "scripts" / "qualify.jl",
                 fuzz.ROOT / "examples" / "rolling_horizon.jl", fuzz.ROOT / "Project.toml"):
        digest.update(path.read_bytes())
    source["source_sha256"] = digest.hexdigest()
    inputs = {r["input"]: hashlib.sha256(Path(r["input"]).read_bytes()).hexdigest()
              for r in requests if r.get("input")}
    configuration = {"requests": requests, "inputs": inputs, "source": source,
                     "julia": args.julia, "timeout": args.timeout}
    metadata = args.output / "run.json"
    if args.resume:
        if not metadata.is_file() or json.loads(metadata.read_text())["configuration"] != configuration:
            parser.error("resume needs the same plan, saved inputs, source and runtime settings")
    elif args.output.exists() and any(args.output.iterdir()):
        parser.error("output directory is not empty; use --resume or another directory")
    else:
        args.output.mkdir(parents=True, exist_ok=True)
        fuzz.atomic_json(metadata, {"configuration": configuration,
            "created_utc": datetime.now(timezone.utc).isoformat(), "platform": platform.platform()})
    results = []
    try:
        for index, request in enumerate(requests, 1):
            directory = args.output / f"{index:06d}"
            print(f"[{index}/{len(requests)}] {request['family']} {request['size']} seed={request['seed']}", flush=True)
            if args.resume and (directory / "result.json").is_file():
                result = json.loads((directory / "result.json").read_text())
            else:
                result = fuzz.run_case(args, directory, request, worker="qualify.jl")
            results.append(result)
            fuzz.atomic_json(args.output / "summary.json", results)
            print(f"  {result['outcome']}", flush=True)
    except KeyboardInterrupt:
        return 130
    return int(any(r["outcome"] == "fail" for r in results))


if __name__ == "__main__":
    sys.exit(main())
