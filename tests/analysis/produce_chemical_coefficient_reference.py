#!/usr/bin/env python3
"""Fresh governed Phase-5C corrected-coefficient artifact producer.

Every invocation constructs and characterizes the Structure-1 fixture, then
executes the compiled production G/Q/Z/I/W path. The reviewed candidate and
governed baseline are never inputs. Compiler provenance comes from the
compiled production executable and is retained without normalization.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


GOVERNED_CLASSIFICATION = (
    "GOVERNED REGRESSION BASELINE; human-ratified generic/free-gas "
    "Phase-5C coefficient scope"
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_logged(command, log_path, *, cwd, env):
    with log_path.open("w") as log:
        return subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=cwd,
            env=env,
        ).returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--characterization-executable", required=True, type=Path)
    parser.add_argument("--production-executable", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    analysis = source_root / "tests" / "analysis"

    characterization_root = output / "characterization"
    characterization_rc = run_logged(
        [
            sys.executable,
            str(analysis / "chemical_trackr_budget.py"),
            str(args.characterization_executable.resolve()),
            str(characterization_root),
        ],
        output / "characterization.log",
        cwd=source_root,
        env=environment,
    )
    if characterization_rc:
        raise SystemExit(
            f"Phase-5C characterization failed; see {output / 'characterization.log'}"
        )

    characterizations = list(characterization_root.glob("run-*/result.json"))
    if len(characterizations) != 1:
        raise SystemExit("Expected exactly one fresh characterization result")

    computed_root = output / "computed"
    production_rc = run_logged(
        [
            sys.executable,
            str(analysis / "chemical_production_evidence.py"),
            str(args.production_executable.resolve()),
            str(characterizations[0]),
            str(computed_root),
        ],
        output / "production.log",
        cwd=source_root,
        env=environment,
    )
    if production_rc:
        raise SystemExit(
            f"Phase-5C production failed; see {output / 'production.log'}"
        )

    computed_candidate = computed_root / "candidate.json"
    artifact = json.loads(computed_candidate.read_text())
    if artifact.get("candidate_only") is not True:
        raise SystemExit("Fresh computed payload lacks candidate-only source classification")
    if "governed_baseline" in artifact or "classification" in artifact:
        raise SystemExit("Unexpected pre-existing governed classification fields")
    compiler = artifact.get("provenance", {}).get("toolchain", {}).get("compiler")
    if not isinstance(compiler, str) or not compiler.strip():
        raise SystemExit("Compiled production path emitted no compiler provenance")

    computed_candidate_sha = digest(computed_candidate)
    artifact["candidate_only"] = False
    artifact["governed_baseline"] = True
    artifact["classification"] = GOVERNED_CLASSIFICATION
    target = output / "phase5c_chemical_coefficients.json"
    target.write_text(
        json.dumps(artifact, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )

    report = {
        "characterization_raw_rc": characterization_rc,
        "compiler": compiler,
        "compiler_provenance_source": "production executable __VERSION__",
        "computed_candidate_sha256": computed_candidate_sha,
        "governed_artifact_sha256": digest(target),
        "production_raw_rc": production_rc,
    }
    (output / "producer-report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(f"compiler={compiler}")
    print(f"computed_candidate_sha256={computed_candidate_sha}")
    print(f"governed_artifact_sha256={digest(target)}")
    print(target)


if __name__ == "__main__":
    main()
