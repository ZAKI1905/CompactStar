#!/usr/bin/env python3
"""Run-local Phase-5D protected-source manifest production and validation."""

import hashlib
import json
import re
import subprocess
from pathlib import Path

ENTRY = "f7116c1408c06f976527f86d4397ad6d4540dedf"
SPECIAL = {
    "tests/baselines/phase5b_structural_response.json":
        "7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa",
    "tests/baselines/phase5c_chemical_coefficients.json":
        "7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7",
    "docs/validation/phase5c_chemical_coefficients_candidate.json":
        "a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b",
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git(source_root, *args, binary=False):
    return subprocess.check_output(
        ["git", *args], cwd=source_root, text=not binary
    )


def protected_authority(source_root):
    source_root = Path(source_root).resolve()
    plan = source_root / "docs/validation/PHASE5D1_FROZEN_CONTEXT_PROVENANCE_PLAN.md"
    values = {
        path: value
        for value, path in re.findall(
            r"^\| `([0-9a-f]{64})` \| `([^`]+)` \|$", plan.read_text(), re.M
        )
    }
    if len(values) != 33:
        raise RuntimeError("frozen-context plan does not define exactly 33 protected paths")
    return values


def generate_manifest(source_root, output):
    """Hash the current tree; no prior manifest or candidate is an input."""
    source_root = Path(source_root).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise RuntimeError("entry manifest destination must be fresh")
    authority = protected_authority(source_root)
    protected = []
    for path, expected in sorted(authority.items()):
        actual = digest(source_root / path)
        if actual != expected:
            raise RuntimeError("protected hash mismatch: " + path)
        protected.append({"path": path, "expected": expected, "actual": actual})
    paths = _git(source_root, "ls-tree", "-r", "--name-only", ENTRY,
                 "tests/baselines").splitlines()
    if len(paths) != 10:
        raise RuntimeError("entry does not contain exactly ten governed baselines")
    baselines = []
    for path in sorted(paths):
        expected = hashlib.sha256(
            _git(source_root, "show", ENTRY + ":" + path, binary=True)
        ).hexdigest()
        actual = digest(source_root / path)
        if actual != expected:
            raise RuntimeError("governed baseline changed: " + path)
        baselines.append({"path": path, "sha256": expected})
    special = []
    for path, expected in sorted(SPECIAL.items()):
        actual = digest(source_root / path)
        if actual != expected:
            raise RuntimeError("protected artifact changed: " + path)
        special.append({"path": path, "expected": expected, "actual": actual})
    result = {
        "entry_sha": ENTRY,
        "protected": protected,
        "baselines": baselines,
        "special": special,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    check_manifest(source_root, output)
    return result


def check_manifest(source_root, entry):
    """Fail closed on topology and compare every byte to current source."""
    source_root = Path(source_root).resolve()
    entry = Path(entry).resolve()
    data = json.loads(entry.read_text())
    if set(data) != {"entry_sha", "protected", "baselines", "special"}:
        raise RuntimeError("entry manifest schema changed")
    if data["entry_sha"] != ENTRY:
        raise RuntimeError("entry SHA mismatch")
    authority = protected_authority(source_root)
    observed = {x["path"]: x["actual"] for x in data["protected"]}
    if len(data["protected"]) != 33 or observed != authority:
        raise RuntimeError("manifest differs from exact 33-path authority")
    if len(observed) != len(data["protected"]):
        raise RuntimeError("duplicate protected manifest path")
    for path, expected in authority.items():
        if digest(source_root / path) != expected:
            raise RuntimeError("protected hash mismatch: " + path)
    baseline_paths = _git(
        source_root, "ls-tree", "-r", "--name-only", ENTRY, "tests/baselines"
    ).splitlines()
    baseline_map = {x["path"]: x["sha256"] for x in data["baselines"]}
    if (len(data["baselines"]) != 10 or len(baseline_map) != 10
            or set(baseline_map) != set(baseline_paths)):
        raise RuntimeError("baseline manifest incomplete")
    for path in baseline_paths:
        expected = hashlib.sha256(
            _git(source_root, "show", ENTRY + ":" + path, binary=True)
        ).hexdigest()
        if baseline_map[path] != expected or digest(source_root / path) != expected:
            raise RuntimeError("baseline changed: " + path)
    special_map = {x["path"]: x["expected"] for x in data["special"]}
    if (len(data["special"]) != 3 or len(special_map) != 3
            or special_map != SPECIAL):
        raise RuntimeError("special artifact manifest mismatch")
    for path, expected in SPECIAL.items():
        if digest(source_root / path) != expected:
            raise RuntimeError("protected artifact changed: " + path)
    return data
