#!/usr/bin/env python3
"""Fresh Phase-5C governed regression with one portable provenance field."""

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


COMPILER_PATH = ("provenance", "toolchain", "compiler")
CANDIDATE_BASELINE_ALLOWED = {
    ("candidate_only",),
    ("classification",),
    ("governed_baseline",),
    COMPILER_PATH,
}
GOVERNED_BASELINE_ALLOWED = {COMPILER_PATH}
EXPECTED_CLASSIFICATION = (
    "GOVERNED REGRESSION BASELINE; human-ratified generic/free-gas "
    "Phase-5C coefficient scope"
)


class ComparisonFailure(RuntimeError):
    pass


def digest(data):
    return hashlib.sha256(data).hexdigest()


def compiler(artifact):
    try:
        value = artifact["provenance"]["toolchain"]["compiler"]
    except (KeyError, TypeError) as error:
        raise ComparisonFailure("compiler provenance is missing") from error
    if not isinstance(value, str) or not value.strip():
        raise ComparisonFailure("compiler provenance is empty")
    return value


def difference_paths(left, right):
    differences = []

    def visit(a, b, path=()):
        if type(a) is not type(b):
            differences.append(path)
        elif isinstance(a, dict):
            for key in sorted(set(a) | set(b)):
                if key not in a or key not in b:
                    differences.append(path + (key,))
                else:
                    visit(a[key], b[key], path + (key,))
        elif isinstance(a, list):
            if len(a) != len(b):
                differences.append(path + ("<length>",))
            for index, (a_value, b_value) in enumerate(zip(a, b)):
                visit(a_value, b_value, path + (index,))
        elif a != b:
            differences.append(path)

    visit(left, right)
    return differences


def compare(left, right, allowed):
    left_compiler = compiler(left)
    right_compiler = compiler(right)
    differences = difference_paths(left, right)
    unexpected = [path for path in differences if path not in allowed]
    if unexpected:
        raise ComparisonFailure(f"unexpected differences: {unexpected}")
    return {
        "compiler_left": left_compiler,
        "compiler_right": right_compiler,
        "compiler_same": left_compiler == right_compiler,
        "difference_paths": [list(path) for path in differences],
    }


def require_governed(name, artifact):
    if artifact.get("candidate_only") is not False:
        raise ComparisonFailure(f"{name} is still candidate-only")
    if artifact.get("governed_baseline") is not True:
        raise ComparisonFailure(f"{name} lacks governed-baseline classification")
    if artifact.get("classification") != EXPECTED_CLASSIFICATION:
        raise ComparisonFailure(f"{name} classification is not canonical")
    if artifact.get("predeclaration_sha") != (
        "a87f0212c2bd7bfba92db91dfac82447a6561334"
    ):
        raise ComparisonFailure(f"{name} predeclaration identity differs")
    if not artifact.get("provenance", {}).get("source_sha256"):
        raise ComparisonFailure(f"{name} source provenance is incomplete")
    compiler(artifact)


def must_fail(label, operation):
    try:
        operation()
    except ComparisonFailure:
        print(f"negative_{label}=PASS")
        return "PASS"
    raise ComparisonFailure(f"negative {label} mutation was accepted")


def run_controls(baseline):
    controls = {}

    numeric = copy.deepcopy(baseline)
    numeric["values"]["W"][0] = math.nextafter(
        numeric["values"]["W"][0], math.inf
    )
    controls["altered_scientific_field"] = must_fail(
        "altered_scientific_field",
        lambda: compare(numeric, baseline, GOVERNED_BASELINE_ALLOWED),
    )

    provenance = copy.deepcopy(baseline)
    source_name = sorted(provenance["provenance"]["source_sha256"])[0]
    provenance["provenance"]["source_sha256"][source_name] = "0" * 64
    controls["altered_scientific_provenance"] = must_fail(
        "altered_scientific_provenance",
        lambda: compare(provenance, baseline, GOVERNED_BASELINE_ALLOWED),
    )

    compiler_only = copy.deepcopy(baseline)
    compiler_only["provenance"]["toolchain"]["compiler"] += " portable-control"
    compiler_result = compare(
        compiler_only, baseline, GOVERNED_BASELINE_ALLOWED
    )
    if compiler_result["compiler_same"]:
        raise ComparisonFailure("compiler-only control did not report difference")
    controls["compiler_only"] = "PASS"
    print("positive_compiler_only=PASS")

    missing = copy.deepcopy(baseline)
    del missing["provenance"]["toolchain"]["compiler"]
    controls["missing_compiler"] = must_fail(
        "missing_compiler",
        lambda: compare(missing, baseline, GOVERNED_BASELINE_ALLOWED),
    )

    empty = copy.deepcopy(baseline)
    empty["provenance"]["toolchain"]["compiler"] = ""
    controls["empty_compiler"] = must_fail(
        "empty_compiler",
        lambda: compare(empty, baseline, GOVERNED_BASELINE_ALLOWED),
    )

    other_toolchain = copy.deepcopy(baseline)
    other_toolchain["provenance"]["toolchain"]["GSL"] += "-mutation"
    controls["other_toolchain_field"] = must_fail(
        "other_toolchain_field",
        lambda: compare(other_toolchain, baseline, GOVERNED_BASELINE_ALLOWED),
    )
    return controls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", required=True, type=Path)
    parser.add_argument("--characterization-executable", required=True, type=Path)
    parser.add_argument("--production-executable", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--reviewed-candidate", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    run_root = Path(tempfile.mkdtemp(prefix="regression-", dir=output_root))
    generated_root = run_root / "generated"
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    # Production completes before either expected artifact is read.
    command = [
        sys.executable,
        str(args.producer.resolve()),
        "--characterization-executable",
        str(args.characterization_executable.resolve()),
        "--production-executable",
        str(args.production_executable.resolve()),
        "--source-root",
        str(args.source_root.resolve()),
        "--output-dir",
        str(generated_root),
    ]
    result = subprocess.run(command, capture_output=True, text=True, env=environment)
    (run_root / "producer.log").write_text(result.stdout + result.stderr)
    print(result.stdout, end="")
    if result.returncode:
        raise SystemExit(f"Phase-5C producer failed; see {run_root / 'producer.log'}")

    generated_path = (
        generated_root / "phase5c_chemical_coefficients.json"
    ).resolve()
    baseline_path = args.baseline.resolve()
    candidate_path = args.reviewed_candidate.resolve()
    if generated_path in (baseline_path, candidate_path):
        raise ComparisonFailure("Refusing generated/expected self-comparison")

    generated_bytes = generated_path.read_bytes()
    baseline_bytes = baseline_path.read_bytes()
    candidate_bytes = candidate_path.read_bytes()
    generated = json.loads(generated_bytes)
    baseline = json.loads(baseline_bytes)
    candidate = json.loads(candidate_bytes)
    require_governed("generated artifact", generated)
    require_governed("governed baseline", baseline)

    governed_result = compare(generated, baseline, GOVERNED_BASELINE_ALLOWED)
    candidate_result = compare(candidate, baseline, CANDIDATE_BASELINE_ALLOWED)
    actual_candidate_paths = {tuple(path) for path in candidate_result["difference_paths"]}
    if actual_candidate_paths != CANDIDATE_BASELINE_ALLOWED:
        raise ComparisonFailure(
            f"candidate/baseline differences are not exact: {actual_candidate_paths}"
        )

    controls = run_controls(baseline)
    evidence = {
        "baseline_sha256": digest(baseline_bytes),
        "candidate_sha256": digest(candidate_bytes),
        "candidate_to_baseline": candidate_result,
        "controls": controls,
        "generated_sha256": digest(generated_bytes),
        "governed_run_to_baseline": governed_result,
        "portable_path": list(COMPILER_PATH),
    }
    (run_root / "comparison-evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    )
    print(f"baseline_compiler={governed_result['compiler_right']}")
    print(f"current_compiler={governed_result['compiler_left']}")
    print(f"compiler_same={str(governed_result['compiler_same']).lower()}")
    print(f"generated_sha256={digest(generated_bytes)}")
    print(f"baseline_sha256={digest(baseline_bytes)}")
    print(f"isolated_evidence={run_root}")


if __name__ == "__main__":
    main()
