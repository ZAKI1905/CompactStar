#!/usr/bin/env python3
"""Create the non-governed Phase-6A-1 candidate from passed evidence only."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

CLASSIFICATION = (
    "CONTROLLED MATHEMATICAL / ARCHITECTURE BNV CANDIDATE; "
    "NOT PHYSICAL BNV MODEL; NOT GOVERNED BASELINE; "
    "NOT OWNER-RATIFIED NUMERICAL RESULT"
)
CARDS = (
    ("RF-P0-TRANSIENT-v1", "P0", -1e-13, 1e6, 1025, True),
    ("CPL-P0-TRANSIENT-v1", "P0", -1e-13, 1e5, 2049, False),
    ("CPL-P1-TRANSIENT-v1", "P1", -1e-13, 1e5, 2049, False),
    ("CPL-P2-LINEAR-QSS-v1", "P2", -1e-12, 5e5, 8193, False),
)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        source = list(csv.DictReader(stream, delimiter="\t"))
    if not source:
        raise RuntimeError(f"empty candidate trajectory {path}")
    result = []
    for row in source:
        item = {}
        for name, value in row.items():
            if value is None:
                raise RuntimeError(f"malformed candidate row {path}")
            try:
                number = float(value)
                if not math.isfinite(number):
                    raise RuntimeError(f"nonfinite candidate value {path}:{name}")
                item[name] = number
            except ValueError:
                item[name] = value
        result.append(item)
    return result


def git(root: Path, *arguments: str) -> str:
    return subprocess.run(("git", *arguments), cwd=root, check=True,
                          text=True, capture_output=True).stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--entry-hashes", type=Path, required=True)
    parser.add_argument("--frozen-certificate", type=Path, required=True)
    parser.add_argument("--pretrajectory", type=Path, required=True)
    parser.add_argument("--predeclaration-sha", required=True)
    parser.add_argument("--implementation-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    if git(root, "rev-parse", "HEAD") != args.implementation_sha:
        raise RuntimeError("implementation SHA is not current HEAD")
    if git(root, "status", "--porcelain"):
        raise RuntimeError("candidate producer refuses a dirty implementation tree")
    validation = json.loads(args.validation.read_text())
    if any(validation.get(gate) != "PASS" for gate in ("BA11", "BA12", "BA13", "BA14", "BA16", "BA17")):
        raise RuntimeError("candidate producer requires BA11-BA17 validation PASS")
    pretrajectory = args.pretrajectory.read_text()
    if "PRETRAJECTORY PASS" not in pretrajectory:
        raise RuntimeError("candidate producer requires durable pretrajectory PASS")
    raw_hashes = {path.name: sha(path) for path in sorted(args.raw.glob("*")) if path.is_file()}
    trajectories = {}
    for identity, partition, drive, duration, checkpoints, reaction_free in CARDS:
        baseline_path = args.raw / f"{identity}.baseline.tsv"
        refined_path = args.raw / f"{identity}.refined.tsv"
        control_path = args.raw / f"{identity}.control.baseline.tsv"
        if any(not path.is_file() for path in (baseline_path, refined_path, control_path)):
            raise RuntimeError(f"incomplete candidate evidence for {identity}")
        baseline_rows = rows(baseline_path)
        if len(baseline_rows) != checkpoints:
            raise RuntimeError(f"candidate checkpoint count changed for {identity}")
        trajectories[identity] = {
            "run_card": {"partition": partition, "fractional_drive_per_year": drive,
                         "duration_year": duration, "checkpoints": checkpoints,
                         "reaction_free": reaction_free, "physical_rate": False,
                         "physical_model": False},
            "baseline_rows": baseline_rows,
            "refined_sha256": sha(refined_path),
            "matched_control_baseline_sha256": sha(control_path),
            "validation": validation["cards"][identity],
        }
    tracked = (
        list((root / "CompactStar/Physics/BNV").rglob("*"))
        + [root / "CompactStar/Analysis/EquilibriumBaryonTangent.hpp",
           root / "CompactStar/Analysis/src/EquilibriumBaryonTangent.cpp"]
        + list((root / "tests/bnv").glob("*"))
    )
    source_hashes = {str(path.relative_to(root)): sha(path) for path in sorted(tracked)
                     if path.is_file() and "__pycache__" not in str(path)}
    payload = {
        "schema_id": "compactstar.phase6a1.controlled-bnv-candidate.v1",
        "classification": CLASSIFICATION,
        "candidate_only": True, "governed_baseline": False,
        "owner_ratified_numerical_result": False,
        "physical_BNV_rate_selected": False, "physical_BNV_model_selected": False,
        "canonical_entry_sha": "961dfa0de6f76df71df4cb98edc8e1b35a5c21b1",
        "predeclaration_sha": args.predeclaration_sha,
        "implementation_sha": args.implementation_sha,
        "branch": "physics/phase6a1-controlled-bnv-implementation",
        "fixture": {"radial_resolution": 80000, "EOS_resolution": 8192,
                    "rho_c_g_cm3": 1.10e15, "initial_Tinf_K": 1e8,
                    "initial_eta_MeV": [0.0, 0.0], "spin": "OFF",
                    "run_purpose": "AnalyticControl"},
        "solver": {"name": "GSL RKF45", "baseline_rtol": 1e-7,
                   "baseline_atol": [1e-12, 1e-18, 1e-18],
                   "refined_rtol": 1e-9, "refined_atol": [1e-14, 1e-20, 1e-20]},
        "units": {"time": "s", "temperature": "K", "chemical_energy": "MeV",
                  "power": "erg s^-1", "particle_number": "count",
                  "particle_rate": "count s^-1"},
        "entry_hashes": json.loads(args.entry_hashes.read_text()),
        "frozen_certificate_sha256": sha(args.frozen_certificate),
        "pretrajectory_record_sha256": sha(args.pretrajectory),
        "source_hashes": source_hashes, "raw_evidence_hashes": raw_hashes,
        "validation": validation, "trajectories": trajectories,
        "forbidden_scope": {"A18": False, "superfluidity": False, "Regime_II": False,
                            "MixedStar_thermal_evolution": False, "sliding_background": False,
                            "variable_Z": False, "governed_baseline_created": False},
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    payload["candidate_payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"CANDIDATE_WRITTEN {args.output} payload_sha256 {payload['candidate_payload_sha256']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"STOP {error}", file=sys.stderr)
        raise SystemExit(1)
