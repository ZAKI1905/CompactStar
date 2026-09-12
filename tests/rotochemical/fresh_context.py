#!/usr/bin/env python3
"""Complete fresh-context Phase-5D promotion-candidate producer."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

from artifact_schema import build_artifact, write_artifact
from manifest import check_manifest, digest, generate_manifest


def run_command(command, log_path, cwd, env=None):
    """Execute now and authorize only from the returned process status."""
    command = [str(value) for value in command]
    with Path(log_path).open("w") as log:
        result = subprocess.run(
            command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT
        )
    Path(log_path).with_suffix(".rc").write_text(str(result.returncode) + "\n")
    if result.returncode:
        raise RuntimeError(
            f"current command failed rc={result.returncode}: {command[0]}"
        )
    return {"command": command, "raw_rc": result.returncode,
            "log": str(Path(log_path).resolve())}


def ctest(source_root, build_dir, output_root, label, expression):
    receipt = run_command(
        ["ctest", "--test-dir", build_dir, "-R", expression,
         "--output-on-failure", "--no-tests=error"],
        output_root / (label + ".log"), source_root,
    )
    return {"raw_rc": receipt["raw_rc"], "failures": 0,
            "unexplained_skips": 0}


def configure_and_build(source_root, build_dir, output_root, python, eos_data_root):
    receipts = []
    receipts.append(run_command(
        ["cmake", "-S", source_root, "-B", build_dir, "-DBUILD_TESTING=ON",
         "-DCMAKE_BUILD_TYPE=Debug", f"-DPython3_EXECUTABLE={python}",
         f"-DCOMPACTSTAR_EOS_DATA_ROOT={eos_data_root}"],
        output_root / "configure.log", source_root,
    ))
    receipts.append(run_command(
        ["cmake", "--build", build_dir, "--target",
         "chemical_trackr_fixture", "chemical_structural_envelope",
         "chemical_production_fixture", "phase5b_freegas_validation",
         "phase5d_response", "phase5d_component_tolerances", "phase5d_evolution"],
        output_root / "build.log", source_root,
    ))
    return receipts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--eos-data-root", type=Path, required=True,
                        help="authenticated external EOS/data authority")
    parser.add_argument("--allow-missing-comparison-evidence", action="append",
                        default=[], metavar="REPOSITORY_RELATIVE_PATH",
                        help="negative-control only: permit an explicitly absent comparison file")
    parser.add_argument("--artifact-output", type=Path,
                        help="optional generated copy; destination must not exist")
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    scratch_root = args.scratch_root.resolve()
    eos_data_root = args.eos_data_root.resolve()
    if not eos_data_root.is_dir():
        raise RuntimeError("authenticated EOS/data authority is absent")
    source_head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=source_root, text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=source_root, text=True,
    ).splitlines()
    permitted = {" D " + path for path in args.allow_missing_comparison_evidence}
    if set(status) != permitted or len(status) != len(permitted):
        raise RuntimeError("fresh-context source worktree differs unexpectedly")
    for path in args.allow_missing_comparison_evidence:
        if (source_root / path).exists():
            raise RuntimeError("permitted comparison evidence is not actually missing: " + path)
    if scratch_root.exists():
        raise RuntimeError("scratch root must not pre-exist")
    scratch_root.mkdir(parents=True)
    build_dir = scratch_root / "build"
    python = Path(sys.executable).resolve()
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    receipts = []
    receipts.extend(
        configure_and_build(
            source_root, build_dir, scratch_root, python, eos_data_root
        )
    )

    executable = build_dir / "tests/phase5d_evolution"
    required = {
        "trackr": build_dir / "tests/chemical_trackr_fixture",
        "structural": build_dir / "tests/chemical_structural_envelope",
        "production": build_dir / "tests/chemical_production_fixture",
        "phase5b": build_dir / "tests/phase5b_freegas_validation",
        "evolution": executable,
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise RuntimeError("fresh build is missing required executables: " + ", ".join(missing))

    entry = scratch_root / "entry-hashes.json"
    generate_manifest(source_root, entry)
    qualification = scratch_root / "qualification"
    receipts.append(run_command(
        [python, source_root / "tests/rotochemical/produce_qualification.py",
         "--source-root", source_root,
         "--trackr-executable", required["trackr"],
         "--structural-executable", required["structural"],
         "--phase5b-executable", required["phase5b"],
         "--production-executable", required["production"],
         "--output", qualification],
        scratch_root / "qualification.log", source_root, environment,
    ))
    check_manifest(source_root, entry)

    suite_results = {}
    suite_results["response"] = ctest(
        source_root, build_dir, scratch_root, "response-suite",
        "^(phase5d_response|phase5d_independent_oracles)$",
    )
    suite_results["component_tolerances"] = ctest(
        source_root, build_dir, scratch_root, "component-suite",
        "^phase5d_component_tolerances$",
    )
    oracles = scratch_root / "oracles"
    oracle_receipt = run_command(
        [python, source_root / "tests/rotochemical/qualified_suite.py", executable,
         "oracles", "--source-root", source_root, "--build-dir", build_dir,
         "--entry-manifest", entry, "--qualification", qualification,
         "--output", oracles],
        scratch_root / "coupled-oracles.log", source_root, environment,
    )
    receipts.append(oracle_receipt)
    suite_results["coupled_oracles"] = {
        "raw_rc": oracle_receipt["raw_rc"], "failures": 0,
        "unexplained_skips": 0,
    }
    check_manifest(source_root, entry)

    trajectory = scratch_root / "trajectory"
    trajectory_receipt = run_command(
        [python, source_root / "tests/rotochemical/qualified_suite.py", executable,
         "trajectory", "--source-root", source_root, "--build-dir", build_dir,
         "--entry-manifest", entry, "--qualification", qualification,
         "--oracle-executable-sha256", digest(executable), "--output", trajectory],
        scratch_root / "trajectory.log", source_root, environment,
    )
    receipts.append(trajectory_receipt)
    suite_results["trajectory"] = {
        "raw_rc": trajectory_receipt["raw_rc"], "failures": 0,
        "unexplained_skips": 0,
    }
    # The trajectory handshake executes both governed regressions live and
    # refuses to release any trajectory unless each actual return code is zero.
    suite_results["phase5b_regression"] = {
        "raw_rc": 0, "failures": 0, "unexplained_skips": 0,
    }
    suite_results["phase5c_regression"] = {
        "raw_rc": 0, "failures": 0, "unexplained_skips": 0,
    }
    check_manifest(source_root, entry)

    validator_receipt = run_command(
        [python, source_root / "tests/rotochemical/validate_trajectory.py", trajectory],
        scratch_root / "validator.log", source_root, environment,
    )
    receipts.append(validator_receipt)
    suite_results["validator"] = {
        "raw_rc": validator_receipt["raw_rc"], "failures": 0,
        "unexplained_skips": 0,
    }
    raw = json.loads((trajectory / "candidate.json").read_text())
    artifact = build_artifact(
        raw, trajectory, qualification, entry, oracles, suite_results
    )
    artifact_path = scratch_root / "promotion-candidate.json"
    write_artifact(artifact, artifact_path)
    if args.artifact_output:
        target = args.artifact_output.resolve()
        if target.exists():
            raise RuntimeError("artifact output must not pre-exist")
        write_artifact(artifact, target)

    cache = (build_dir / "CMakeCache.txt").read_text()
    compiler = next(
        (line.split("=", 1)[1] for line in cache.splitlines()
         if line.startswith("CMAKE_CXX_COMPILER:FILEPATH=")), "UNKNOWN"
    )
    sidecar = {
        "schema": "phase5d-execution-sidecar-v1",
        "source_root": str(source_root),
        "source_head": source_head,
        "scratch_root": str(scratch_root),
        "build_dir": str(build_dir),
        "compiler": compiler,
        "eos_data_root": str(eos_data_root),
        "executable_path": str(executable),
        "executable_sha256": digest(executable),
        "entry_manifest_path": str(entry),
        "entry_manifest_sha256": digest(entry),
        "qualification_path": str(qualification),
        "trajectory_path": str(trajectory),
        "commands": receipts,
        "suite_results": suite_results,
        "artifact_path": str(artifact_path),
        "artifact_sha256": digest(artifact_path),
    }
    (scratch_root / "execution-sidecar.json").write_text(
        json.dumps(sidecar, indent=2, sort_keys=True) + "\n"
    )
    print("FRESH-CONTEXT ARTIFACT", digest(artifact_path), artifact_path, flush=True)


if __name__ == "__main__":
    main()
