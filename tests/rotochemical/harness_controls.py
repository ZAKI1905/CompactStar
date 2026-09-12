#!/usr/bin/env python3
"""Fast negative controls for the fresh-context orchestration boundary."""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
from compare_artifacts import validate_promotion
from fresh_context import run_command
from qualified_suite import prepare, verify_hashes


def must_refuse(label, action):
    try:
        action()
    except (RuntimeError, subprocess.CalledProcessError):
        print(label, "REFUSED")
    else:
        raise RuntimeError("negative control escaped: " + label)


def main():
    source = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="phase5d-harness-controls-") as name:
        scratch = Path(name)
        forged = scratch / "phase5b-entry.rc"
        forged.write_text("0\n")
        must_refuse(
            "forged_saved_rc",
            lambda: run_command(
                [sys.executable, "-c", "raise SystemExit(9)"],
                scratch / "actual-failure.log", source,
            ),
        )
        if (scratch / "actual-failure.rc").read_text().strip() != "9":
            raise RuntimeError("actual process return code was not recorded")

        empty_qualification = scratch / "missing-qualification"
        empty_qualification.mkdir()
        output = scratch / "qualification-consumer"
        output.mkdir()
        must_refuse(
            "missing_qualification", lambda: prepare(empty_qualification, output)
        )

        thermal = scratch / "eos.thermo"
        thermal.write_text("fresh bytes\n")
        import hashlib
        expected = {thermal: hashlib.sha256(thermal.read_bytes()).hexdigest()}
        thermal.write_text("altered bytes\n")
        must_refuse("thermal_source_mutation", lambda: verify_hashes(expected))

        # Unknown schema members fail closed. This minimal object is rejected
        # before any equality comparison or arbitrary ignore list can act.
        must_refuse("unknown_schema_field", lambda: validate_promotion({"extra": 1}))

        producer_sources = "\n".join(
            (source / path).read_text()
            for path in [
                "tests/rotochemical/fresh_context.py",
                "tests/rotochemical/produce_qualification.py",
                "tests/rotochemical/artifact_schema.py",
            ]
        )
        if "phase5d1_controlled_evolution_candidate.json" in producer_sources:
            raise RuntimeError("producer names the historical candidate")
        if "tests/baselines/phase5d" in producer_sources:
            raise RuntimeError("producer names a future Phase-5D baseline")
        print("historical_candidate_missing GENERATION_INDEPENDENT")
        print("promotion_candidate_missing GENERATION_INDEPENDENT")

    print("PASS live rc, missing qualification, thermal mutation, closed schema, input independence")


if __name__ == "__main__":
    main()
