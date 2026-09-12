#!/usr/bin/env python3
"""CTest adapter: regenerate every coupled-oracle dependency in fresh scratch."""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
from manifest import generate_manifest


def run(command, cwd):
    result = subprocess.run([str(value) for value in command], cwd=cwd)
    if result.returncode:
        raise RuntimeError(f"current command failed rc={result.returncode}: {command[0]}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--evolution-executable", type=Path, required=True)
    parser.add_argument("--trackr-executable", type=Path, required=True)
    parser.add_argument("--structural-executable", type=Path, required=True)
    parser.add_argument("--phase5b-executable", type=Path, required=True)
    parser.add_argument("--production-executable", type=Path, required=True)
    args = parser.parse_args()
    source = args.source_root.resolve()
    scratch_parent = args.scratch_root.resolve()
    scratch_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="coupled-oracles-", dir=scratch_parent) as name:
        scratch = Path(name)
        entry = scratch / "entry-hashes.json"
        generate_manifest(source, entry)
        qualification = scratch / "qualification"
        run([
            sys.executable, source / "tests/rotochemical/produce_qualification.py",
            "--source-root", source,
            "--trackr-executable", args.trackr_executable,
            "--structural-executable", args.structural_executable,
            "--phase5b-executable", args.phase5b_executable,
            "--production-executable", args.production_executable,
            "--output", qualification,
        ], source)
        run([
            sys.executable, source / "tests/rotochemical/qualified_suite.py",
            args.evolution_executable, "oracles", "--source-root", source,
            "--build-dir", args.build_dir, "--entry-manifest", entry,
            "--qualification", qualification, "--output", scratch / "oracles",
        ], source)
    print("PASS coupled oracles are fresh-context reproducible")


if __name__ == "__main__":
    main()
