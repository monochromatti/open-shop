#!/usr/bin/env python3
"""Fetch the reviewed SINTEF checkout as external, ignored example data.

This script grants no rights to redistribute upstream material. It fetches a
pinned public revision and bypasses optional Git LFS assets, which this example
does not need. Existing directories are verified, never reset or overwritten.
"""
import argparse
import pathlib
import subprocess

SOURCE_URL = "https://gitlab.sintef.no/energy/open-modelling-tools/open-datasets/tokke-vinje-watercourse.git"
PINNED_COMMIT = "ba2f2fc1f95d18978a04dd2c658aeef79126981b"


def git(*args):
    return subprocess.check_output(["git", *map(str, args)], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent
        / "data/tokke-vinje-watercourse",
    )
    args = parser.parse_args()
    if args.source.exists():
        commit = git("-C", args.source, "rev-parse", "HEAD")
        if commit != PINNED_COMMIT:
            raise ValueError(
                f"existing checkout is {commit}; expected {PINNED_COMMIT}; use another --source directory"
            )
    else:
        args.source.parent.mkdir(parents=True, exist_ok=True)
        git("clone", "--no-checkout", "--depth", "1", SOURCE_URL, args.source)
        # The default branch may have advanced since this example was reviewed.
        git("-C", args.source, "fetch", "--depth", "1", "origin", PINNED_COMMIT)
        git(
            "-C",
            args.source,
            "-c",
            "filter.lfs.smudge=",
            "-c",
            "filter.lfs.process=",
            "-c",
            "filter.lfs.required=false",
            "checkout",
            "--detach",
            PINNED_COMMIT,
        )
    # The imported CSV/YAML/Markdown inputs must be regular checked-out files.
    for name in (
        "README.md",
        "Tokke_topology_manual_modification.md",
        "Tokke_time_series_manual_modification.md",
        "data/input_data/SEnDHub/draft_topology.yaml",
    ):
        if not (args.source / name).is_file():
            raise ValueError(f"missing checked-out input: {name}")
    print(f"External source: {args.source.resolve()}\nReviewed commit: {PINNED_COMMIT}")


if __name__ == "__main__":
    main()
