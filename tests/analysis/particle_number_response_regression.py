"""Regenerate the governed Phase-5B response with narrow compiler portability."""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


COMPILER_PATH = "provenance.build.compiler"
COMPILER_SENTINEL = "<compiler-excluded-by-ADR-0011-portability-rule>"


class ComparisonFailure(RuntimeError):
    pass


def compiler(artifact, name):
    try:
        value = artifact["provenance"]["build"]["compiler"]
    except (KeyError, TypeError) as error:
        raise ComparisonFailure(f"{name} compiler provenance is missing") from error
    if not isinstance(value, str) or not value.strip():
        raise ComparisonFailure(f"{name} compiler provenance is empty")
    return value


def parse_artifact(data, name):
    try:
        artifact = json.loads(data)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ComparisonFailure(f"{name} is not valid JSON") from error
    if not isinstance(artifact, dict):
        raise ComparisonFailure(f"{name} JSON root is not an object")
    compiler(artifact, name)
    return artifact


def exact_parsed_bytes(artifact):
    return json.dumps(
        artifact,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()


def compare_artifacts(generated_bytes, baseline_bytes):
    generated = parse_artifact(generated_bytes, "generated artifact")
    baseline = parse_artifact(baseline_bytes, "governed baseline")
    generated_compiler = compiler(generated, "generated artifact")
    baseline_compiler = compiler(baseline, "governed baseline")
    compiler_same = generated_compiler == baseline_compiler

    if compiler_same:
        if generated_bytes != baseline_bytes:
            raise ComparisonFailure(
                "Same-compiler Phase-5B artifact is not raw-byte identical"
            )
    else:
        generated_projection = copy.deepcopy(generated)
        baseline_projection = copy.deepcopy(baseline)
        generated_projection["provenance"]["build"]["compiler"] = COMPILER_SENTINEL
        baseline_projection["provenance"]["build"]["compiler"] = COMPILER_SENTINEL
        if exact_parsed_bytes(generated_projection) != exact_parsed_bytes(
            baseline_projection
        ):
            raise ComparisonFailure(
                "Cross-compiler Phase-5B artifact differs outside " + COMPILER_PATH
            )

    return {
        "baseline_compiler": baseline_compiler,
        "compiler_same": compiler_same,
        "generated_compiler": generated_compiler,
        "portable_path": COMPILER_PATH,
    }


def encode(artifact):
    return (json.dumps(artifact, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def must_fail(label, operation):
    try:
        operation()
    except ComparisonFailure:
        print(f"negative_{label}=PASS")
        return "PASS"
    raise ComparisonFailure(f"negative {label} control was accepted")


def run_controls(baseline_bytes):
    baseline = parse_artifact(baseline_bytes, "control baseline")
    controls = {}

    identical = compare_artifacts(baseline_bytes, baseline_bytes)
    if not identical["compiler_same"]:
        raise ComparisonFailure("identical-artifact control reported compiler difference")
    controls["identical_artifact_identical_compiler"] = "PASS"
    print("positive_identical_artifact_identical_compiler=PASS")

    compiler_only = copy.deepcopy(baseline)
    compiler_only["provenance"]["build"]["compiler"] += " portable-control"
    compiler_only_result = compare_artifacts(encode(compiler_only), baseline_bytes)
    if compiler_only_result["compiler_same"]:
        raise ComparisonFailure("compiler-only control did not report compiler difference")
    controls["compiler_only"] = "PASS"
    print("positive_compiler_only=PASS")

    missing_generated = copy.deepcopy(baseline)
    del missing_generated["provenance"]["build"]["compiler"]
    controls["missing_generated_compiler"] = must_fail(
        "missing_generated_compiler",
        lambda: compare_artifacts(encode(missing_generated), baseline_bytes),
    )

    missing_baseline = copy.deepcopy(baseline)
    del missing_baseline["provenance"]["build"]["compiler"]
    controls["missing_baseline_compiler"] = must_fail(
        "missing_baseline_compiler",
        lambda: compare_artifacts(baseline_bytes, encode(missing_baseline)),
    )

    empty = copy.deepcopy(baseline)
    empty["provenance"]["build"]["compiler"] = ""
    controls["empty_compiler"] = must_fail(
        "empty_compiler",
        lambda: compare_artifacts(encode(empty), baseline_bytes),
    )

    numeric = copy.deepcopy(baseline)
    numeric["coefficients"][0]["K"] = math.nextafter(
        numeric["coefficients"][0]["K"], math.inf
    )
    controls["numeric_coefficient"] = must_fail(
        "numeric_coefficient",
        lambda: compare_artifacts(encode(numeric), baseline_bytes),
    )

    numerical_error = copy.deepcopy(baseline)
    numerical_error["coefficients"][0]["K_error"] = math.nextafter(
        numerical_error["coefficients"][0]["K_error"], math.inf
    )
    controls["numerical_error"] = must_fail(
        "numerical_error",
        lambda: compare_artifacts(encode(numerical_error), baseline_bytes),
    )

    eos_sha = copy.deepcopy(baseline)
    eos_sha["provenance"]["EOS_table_sha256"] = "0" * 64
    controls["eos_table_sha"] = must_fail(
        "eos_table_sha",
        lambda: compare_artifacts(encode(eos_sha), baseline_bytes),
    )

    architecture = copy.deepcopy(baseline)
    architecture["provenance"]["build"]["architecture"] += "-mutation"
    controls["architecture"] = must_fail(
        "architecture",
        lambda: compare_artifacts(encode(architecture), baseline_bytes),
    )

    configuration = copy.deepcopy(baseline)
    configuration["provenance"]["build"]["configuration"] += "-mutation"
    controls["configuration"] = must_fail(
        "configuration",
        lambda: compare_artifacts(encode(configuration), baseline_bytes),
    )

    other_provenance = copy.deepcopy(baseline)
    other_provenance["provenance"]["EOS_revision"] += "-mutation"
    controls["other_provenance"] = must_fail(
        "other_provenance",
        lambda: compare_artifacts(encode(other_provenance), baseline_bytes),
    )

    compiler_and_science = copy.deepcopy(baseline)
    compiler_and_science["provenance"]["build"]["compiler"] += " portable-control"
    compiler_and_science["coefficients"][0]["I_phys"] = math.nextafter(
        compiler_and_science["coefficients"][0]["I_phys"], math.inf
    )
    controls["compiler_plus_science"] = must_fail(
        "compiler_plus_science",
        lambda: compare_artifacts(encode(compiler_and_science), baseline_bytes),
    )
    return controls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    run_root = Path(tempfile.mkdtemp(prefix="regression-", dir=output_root))
    generated_root = run_root / "generated"
    command = [
        sys.executable,
        str(args.producer.resolve()),
        "--executable",
        str(args.executable.resolve()),
        "--source-root",
        str(args.source_root.resolve()),
        "--output-dir",
        str(generated_root),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    (run_root / "producer.log").write_text(result.stdout + result.stderr)
    print(result.stdout, end="")
    if result.returncode:
        raise SystemExit(f"Phase-5B producer failed; see {run_root / 'producer.log'}")

    generated = (generated_root / "structural-response.json").resolve()
    baseline = args.baseline.resolve()
    if generated == baseline:
        raise SystemExit(
            "Refusing a self-comparison: generated artifact is the governed baseline"
        )

    generated_bytes = generated.read_bytes()
    baseline_bytes = baseline.read_bytes()
    generated_sha = hashlib.sha256(generated_bytes).hexdigest()
    baseline_sha = hashlib.sha256(baseline_bytes).hexdigest()
    comparison = compare_artifacts(generated_bytes, baseline_bytes)
    controls = run_controls(baseline_bytes)
    evidence = {
        "baseline_sha256": baseline_sha,
        "comparison": comparison,
        "controls": controls,
        "generated_sha256": generated_sha,
    }
    (run_root / "comparison-evidence.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    )
    print(f"baseline_compiler={comparison['baseline_compiler']}")
    print(f"generated_compiler={comparison['generated_compiler']}")
    print(f"compiler_same={str(comparison['compiler_same']).lower()}")
    print(f"generated_sha256={generated_sha}")
    print(f"baseline_sha256={baseline_sha}")
    print(f"isolated_evidence={run_root}")


if __name__ == "__main__":
    main()
