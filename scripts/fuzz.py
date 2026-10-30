#!/usr/bin/env python3
"""Run isolated, reproducible optimization cases. No third-party Python packages."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ("cascade", "pressure", "rivers", "small_storage", "channels", "pumping")
SIZES = ("small", "medium", "large")
STAGES = ("local", "warm", "cold")


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def selection(value, allowed):
    result = value.split(",")
    if not result or len(set(result)) != len(result) or any(v not in allowed for v in result):
        raise argparse.ArgumentTypeError("choose comma-separated values from " + ",".join(allowed))
    return result


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=100, help="total cases, including resumed cases")
    parser.add_argument("--seed", type=int, default=20261023, help="campaign seed")
    parser.add_argument("--families", type=lambda s: selection(s, FAMILIES), default=list(FAMILIES))
    parser.add_argument("--sizes", type=lambda s: selection(s, SIZES), default=["small", "medium"])
    parser.add_argument("--stages", type=lambda s: selection(s, STAGES), default=list(STAGES))
    parser.add_argument("--generate-only", action="store_true", help="check references without optimization")
    parser.add_argument("--time-limit", type=float, default=60.0, help="seconds per optimization")
    parser.add_argument("--timeout", type=float, default=900.0, help="hard seconds per isolated case")
    parser.add_argument("--relative-gap", type=float, default=0.01)
    parser.add_argument("--replay-factor", type=int, default=4)
    parser.add_argument("--julia", default=shutil.which("julia"), help="Julia executable; run inside nix develop")
    parser.add_argument("--output", type=Path, help="artifact directory, outside the package checkout")
    parser.add_argument("--resume", action="store_true", help="skip completed cases in --output")
    parser.add_argument("--replay", type=Path, help="re-run an exact saved input.json")
    parser.add_argument("--strict", action="store_true", help="also exit nonzero for warnings")
    args = parser.parse_args(argv)
    if args.cases < 1 or not 0 <= args.seed < 2**63:
        parser.error("--cases must be positive and --seed must fit a nonnegative signed 64-bit integer")
    if not 0 < args.time_limit < float("inf") or not 0 < args.timeout < float("inf"):
        parser.error("time limits must be finite and positive")
    if not 0 <= args.relative_gap < 1 or args.replay_factor < 2:
        parser.error("--relative-gap must be in [0,1) and --replay-factor must be at least 2")
    if not args.julia:
        parser.error("Julia was not found; use nix develop -c python3 scripts/fuzz.py")
    args.julia = str(Path(args.julia).resolve())
    if not Path(args.julia).is_file():
        parser.error("--julia must name an executable file")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    args.output = (args.output or ROOT.parent / "work" / ("fuzz-" + stamp)).resolve()
    if args.output == ROOT or ROOT in args.output.parents:
        parser.error("keep run artifacts outside the package checkout")
    if args.resume and not (args.output / "run.json").is_file():
        parser.error("--resume needs an existing --output with run.json")
    if args.replay:
        args.replay = args.replay.resolve()
        if not args.replay.is_file():
            parser.error("--replay must name a saved input.json")
        args.cases = 1
    if args.generate_only:
        args.stages = []
    return args


def source_identity():
    paths = sorted((ROOT / "src").glob("*.jl")) + [ROOT / "Manifest.toml"]
    paths += [ROOT / "scripts" / name for name in ("fuzz.py", "fuzz_cases.jl", "fuzz_worker.jl")]
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(ROOT)).encode())
        digest.update(path.read_bytes())
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                              text=True, capture_output=True, check=True).stdout.strip()
    return {"git_revision": revision, "source_sha256": digest.hexdigest()}


def case_request(args, index):
    # Every case has its own RNG stream. Resuming or reordering a campaign cannot
    # change subsequent inputs. Exact saved inputs also survive RNG version changes.
    digest = hashlib.sha256(f"{args.seed}:{index}".encode()).digest()
    seed = int.from_bytes(digest[:8], "little") & (2**63 - 1)
    request = {"seed": str(seed), "family": args.families[(index - 1) % len(args.families)],
               "size": args.sizes[((index - 1) // len(args.families)) % len(args.sizes)],
               "stages": args.stages, "time_limit": args.time_limit,
               "relative_gap": args.relative_gap, "replay_factor": args.replay_factor,
               "replay": str(args.replay) if args.replay else None}
    if args.replay and (args.replay.parent / "request.json").is_file():
        original = json.loads((args.replay.parent / "request.json").read_text())
        request.update({key: original[key] for key in ("seed", "family", "size")})
        request["seed"] = str(request["seed"])
    return request


def run_case(args, directory, request, worker="fuzz_worker.jl"):
    directory.mkdir(parents=True, exist_ok=True)
    atomic_json(directory / "request.json", request)
    command = [args.julia, "--startup-file=no", "--pkgimages=no", f"--project={ROOT}",
               str(ROOT / "scripts" / worker), str(directory / "request.json")]
    environment = {**os.environ, "JULIA_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    environment.setdefault("JULIA_DEPOT_PATH", str(ROOT / ".depot"))
    started = time.monotonic()
    failure = None
    with (directory / "worker.log").open("w") as log:
        try:
            process = subprocess.run(command, cwd=ROOT, env=environment, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
            if process.returncode:
                failure = ("process_exit", f"worker exited with code {process.returncode}")
        except subprocess.TimeoutExpired:
            failure = ("hard_timeout", f"worker exceeded {args.timeout:g} seconds")
        except OSError as error:
            failure = ("launch_error", str(error))
    result_path = directory / "result.json"
    if failure is None:
        try:
            result = json.loads(result_path.read_text())
            if not isinstance(result, dict) or result.get("outcome") not in ("pass", "review", "fail"):
                raise ValueError("missing or invalid worker outcome")
            if not isinstance(result.get("findings"), list) or any(
                not isinstance(f, dict) or not isinstance(f.get("kind"), str)
                or f.get("severity") not in ("warning", "error") for f in result["findings"]
            ):
                raise ValueError("missing or malformed worker findings")
        except (OSError, ValueError) as error:
            failure = ("missing_result", str(error))
    if failure is not None:
        result = {"seed": request["seed"], "family": request["family"], "size": request["size"],
                  "outcome": "fail", "findings": [{"severity": "error", "kind": failure[0],
                                                       "stage": "process", "detail": failure[1]}]}
        # Retain completed stages if a later stage crashed or timed out.
        partial = directory / "partial.json"
        if partial.is_file():
            result["partial"] = json.loads(partial.read_text())
    result["wall_seconds"] = time.monotonic() - started
    atomic_json(result_path, result)
    return result


def main(argv=None):
    args = arguments(argv)
    configuration = {**case_request(args, 1), "seed": args.seed, "families": args.families,
                     "sizes": args.sizes, "timeout": args.timeout, "julia": args.julia,
                     "source": source_identity()}
    if args.replay:
        configuration["replay_sha256"] = hashlib.sha256(args.replay.read_bytes()).hexdigest()
    metadata_path = args.output / "run.json"
    if args.resume:
        old = json.loads(metadata_path.read_text())
        if old["configuration"] != configuration:
            raise SystemExit("Resume configuration or source changed. Start a new run or use --replay.")
    elif args.output.exists() and any(args.output.iterdir()):
        raise SystemExit("Output directory is not empty. Use --resume or choose another directory.")
    else:
        args.output.mkdir(parents=True, exist_ok=True)
        atomic_json(metadata_path, {"configuration": configuration,
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                    "python": platform.python_version(), "platform": platform.platform()})
    counts = {"pass": 0, "review": 0, "fail": 0}
    findings = {}
    print(f"Artifacts: {args.output}", flush=True)
    try:
        for index in range(1, args.cases + 1):
            directory = args.output / f"{index:06d}"
            if args.resume and (directory / "result.json").is_file():
                result = json.loads((directory / "result.json").read_text())
            else:
                request = case_request(args, index)
                print(f"[{index}/{args.cases}] {request['family']} {request['size']} "
                      f"seed={request['seed']}", flush=True)
                result = run_case(args, directory, request)
            counts[result["outcome"]] += 1
            for finding in result["findings"]:
                kind = finding["kind"]
                findings[kind] = findings.get(kind, 0) + 1
            atomic_json(args.output / "summary.json", {"completed": sum(counts.values()),
                        "outcomes": counts, "findings": findings})
            print(f"  {result['outcome']}: " + ", ".join(
                sorted({f['kind'] for f in result['findings']})), flush=True)
    except KeyboardInterrupt:
        print("\nInterrupted. Completed cases are saved; use --resume to continue.", file=sys.stderr)
        return 130
    print(json.dumps({"outcomes": counts, "findings": findings}, indent=2))
    return 1 if counts["fail"] or (args.strict and counts["review"]) else 0


if __name__ == "__main__":
    sys.exit(main())
