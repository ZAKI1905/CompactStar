#!/usr/bin/env python3
"""Fresh governed Phase-5D controlled-evolution regression and controls."""

import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True

from compare_artifacts import (  # noqa: E402
    compare_baseline,
    compare_governed,
    validate_governed,
)


class ComparisonFailure(RuntimeError):
    pass


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_comparison(generated_path, baseline_path):
    generated_path = Path(generated_path).resolve()
    baseline_path = Path(baseline_path).resolve()
    if generated_path == baseline_path:
        raise ComparisonFailure("refusing generated/baseline self-comparison")
    generated = json.loads(generated_path.read_text())
    validate_governed(generated)
    if not baseline_path.is_file():
        raise ComparisonFailure("expected governed baseline is absent after generation")
    baseline = json.loads(baseline_path.read_text())
    return compare_governed(generated, baseline)


def validate_sidecar(sidecar_path, generated_path, scratch_root):
    sidecar_path = Path(sidecar_path).resolve()
    generated_path = Path(generated_path).resolve()
    scratch_root = Path(scratch_root).resolve()
    if not sidecar_path.is_file():
        raise ComparisonFailure("execution-provenance sidecar is absent")
    sidecar = json.loads(sidecar_path.read_text())
    required = {
        "schema", "source_root", "source_head", "scratch_root", "build_dir",
        "compiler", "eos_data_root", "executable_path", "executable_sha256",
        "entry_manifest_path", "entry_manifest_sha256", "qualification_path",
        "trajectory_path", "commands", "suite_results", "artifact_path",
        "artifact_sha256", "artifact_state",
    }
    if set(sidecar) != required:
        raise ComparisonFailure("execution-provenance sidecar schema differs")
    if sidecar["schema"] != "phase5d-execution-sidecar-v1":
        raise ComparisonFailure("execution-provenance sidecar identity differs")
    if sidecar["artifact_state"] != "governed":
        raise ComparisonFailure("execution sidecar does not describe governed production")
    if (not isinstance(sidecar["compiler"], str)
            or not sidecar["compiler"].strip()
            or sidecar["compiler"] == "UNKNOWN"):
        raise ComparisonFailure("execution sidecar compiler provenance is incomplete")
    if Path(sidecar["scratch_root"]).resolve() != scratch_root:
        raise ComparisonFailure("execution sidecar scratch provenance differs")
    if Path(sidecar["artifact_path"]).resolve() != generated_path:
        raise ComparisonFailure("execution sidecar artifact path differs")
    if sidecar["artifact_sha256"] != digest(generated_path):
        raise ComparisonFailure("execution sidecar artifact hash differs")
    if not isinstance(sidecar["commands"], list) or not sidecar["commands"]:
        raise ComparisonFailure("execution sidecar command provenance is incomplete")
    if not isinstance(sidecar["suite_results"], dict):
        raise ComparisonFailure("execution sidecar suite provenance is incomplete")
    return sidecar


def must_refuse(label, operation):
    try:
        operation()
    except (ComparisonFailure, RuntimeError, FileNotFoundError):
        print("negative_" + label + "=PASS")
        return "PASS"
    raise ComparisonFailure("negative control escaped: " + label)


def mutation_controls(baseline):
    controls = {}

    def reject(label, path, value):
        mutated = copy.deepcopy(baseline)
        target = mutated
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        controls[label] = must_refuse(
            label, lambda: compare_governed(mutated, baseline)
        )

    checkpoint = baseline["trajectory_checkpoints"][1]
    numeric_key = next(key for key, value in checkpoint.items()
                       if isinstance(value, float))
    reject(
        "numeric_trajectory",
        ["trajectory_checkpoints", 1, numeric_key],
        math.nextafter(checkpoint[numeric_key], math.inf),
    )
    w_value = baseline["coefficient_authority"]["W_I"]["W"][0]
    reject(
        "Z_or_W", ["coefficient_authority", "W_I", "W", 0],
        math.nextafter(w_value, math.inf),
    )
    ltilde_process = sorted(baseline["coefficient_authority"]["Ltilde"])[0]
    ltilde_value = baseline["coefficient_authority"]["Ltilde"][ltilde_process]["value"]
    reject(
        "Ltilde", ["coefficient_authority", "Ltilde", ltilde_process, "value"],
        math.nextafter(ltilde_value, math.inf),
    )
    reject(
        "fixture_metadata", ["fixture", "rho_c_g_cm3"],
        math.nextafter(baseline["fixture"]["rho_c_g_cm3"], math.inf),
    )
    reject(
        "solver_tolerance", ["solver", "configuration", "rtol"], 2.0e-7,
    )
    source_path = sorted(
        baseline["source_provenance"]["scientific_production_source_hashes"]
    )[0]
    reject(
        "source_provenance",
        ["source_provenance", "scientific_production_source_hashes", source_path],
        "0" * 64,
    )
    reject(
        "envelope_label", ["thermal_authority", "envelope_provenance"],
        "iron Potekhin1997",
    )
    reject(
        "physical_spin_interpretation",
        ["classification", "physical_spin_interpretation"], True,
    )

    promotion = copy.deepcopy(baseline)
    promotion["classification"].update({
        "classification": "promotion_candidate",
        "candidate_only": True,
        "governed_baseline": False,
    })
    differences = compare_baseline(promotion, baseline)
    expected = {
        ("classification", "classification"),
        ("classification", "candidate_only"),
        ("classification", "governed_baseline"),
    }
    if {tuple(item["path"]) for item in differences} != expected:
        raise ComparisonFailure("classification transition inventory differs")
    controls["classification_transition"] = "PASS"
    print("positive_classification_transition=PASS")
    return controls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--eos-data-root", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    run_root = Path(tempfile.mkdtemp(prefix="regression-", dir=output_root))
    scratch_root = run_root / "fresh-governed"
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    command = [
        sys.executable, str(args.producer.resolve()),
        "--source-root", str(args.source_root.resolve()),
        "--scratch-root", str(scratch_root),
        "--eos-data-root", str(args.eos_data_root.resolve()),
        "--artifact-state", "governed",
    ]
    with (run_root / "producer.log").open("w") as log:
        result = subprocess.run(
            command, cwd=args.source_root.resolve(), env=environment,
            stdout=log, stderr=subprocess.STDOUT,
        )
    (run_root / "producer.rc").write_text(str(result.returncode) + "\n")
    if result.returncode:
        raise RuntimeError("fresh governed producer failed: " + str(run_root))

    generated_path = scratch_root / "governed-artifact.json"
    if not generated_path.is_file():
        raise RuntimeError("fresh governed producer emitted no scientific artifact")
    sidecar_path = scratch_root / "execution-sidecar.json"
    validate_sidecar(sidecar_path, generated_path, scratch_root)
    controls = {
        "baseline_missing_after_generation": must_refuse(
            "baseline_missing_after_generation",
            lambda: require_comparison(
                generated_path, run_root / "deliberately-absent-baseline.json"
            ),
        )
    }
    require_comparison(generated_path, args.baseline)
    baseline = json.loads(args.baseline.read_text())
    controls.update(mutation_controls(baseline))
    evidence = {
        "baseline_sha256": digest(args.baseline),
        "controls": controls,
        "execution_sidecar": str(sidecar_path),
        "generated_sha256": digest(generated_path),
        "producer_raw_rc": result.returncode,
        "scratch_root": str(scratch_root),
    }
    (run_root / "regression-evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    )
    print("generated_sha256=" + digest(generated_path))
    print("baseline_sha256=" + digest(args.baseline))
    print("fresh_scratch_root=" + str(scratch_root))
    print("PASS fresh governed generation, exact comparison, and controls")


if __name__ == "__main__":
    main()
