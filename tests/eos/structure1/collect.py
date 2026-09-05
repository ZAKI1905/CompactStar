#!/usr/bin/env python3
"""Collect Structure-1 generated diagnostics; never create a scientific baseline.

Run the two C++ executables serially with fresh output directories, then pass
those directories here. The compiler/test executables own all numerical outputs;
this standard-library-only tool packages values, provenance and SHA-256.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure", type=Path, required=True)
    ap.add_argument("--local", type=Path, required=True)
    ap.add_argument("--regression", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    with (args.structure / "stars.tsv").open() as f:
        stars = list(csv.DictReader(f, delimiter="\t"))
    numeric = [k for k in stars[0] if k not in ("group", "status")]
    for row in stars:
        for k in numeric:
            row[k] = float(row[k])
        assert row["status"] == "SURFACE_REACHED"
    midpoint = next(q for q in stars if q["group"] == "radial" and q["grid"] == 8192 and q["partition"] == 40000)
    finest = [q for q in stars if q["group"] == "radial" and q["grid"] == 8192]
    sequence = [q for q in stars if q["group"] == "rho-bin" and q["grid"] == 65]
    floors = [q for q in stars if q["group"] == "floor"]
    branch = [q for q in stars if q["group"] == "high-branch"]
    local = [line.split("\t") for line in (args.local / "local.tsv").read_text().splitlines()]
    sources = ["CompactStar/EOS/LocalThermodynamics.hpp", "CompactStar/EOS/src/LocalThermodynamics.cpp",
               "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp", "CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp",
               "tests/CMakeLists.txt", "tests/eos/rotochemical_trackr_freegas_barotrope.cpp",
               "tests/eos/rotochemical_trackr_freegas_structure.cpp"]
    sources += [str(p.relative_to(root)) for p in sorted((root / "tests/eos/structure1").glob("*")) if p.is_file()]
    artifacts = {}
    for name, folder in [("structure", args.structure), ("local", args.local), ("regression", args.regression)]:
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                artifacts[name + "/" + str(path.relative_to(folder))] = {"sha256": sha(path), "bytes": path.stat().st_size}
    regress = json.loads((args.regression / "runtimes.json").read_text())
    baseline = {str(p.relative_to(root)): sha(p) for p in sorted((root / "tests/baselines").glob("*.tsv"))}
    assert len(baseline) == 8
    for name, actual in baseline.items():
        expected = hashlib.sha256(subprocess.check_output(["git", "show", "f4ae22d971e25bdd74530aec184f3fe0c3440b95:" + name], cwd=root)).hexdigest()
        assert actual == expected
    protected = {}
    for name in ["CompactStar/Core/src/TOVSolver.cpp", "CompactStar/Core/TOVSolver.hpp", "CompactStar/Core/src/NStar.cpp", "CompactStar/Core/StarProfile.hpp"]:
        actual = (root / name).read_bytes()
        assert actual == subprocess.check_output(["git", "show", "f4ae22d971e25bdd74530aec184f3fe0c3440b95:" + name], cwd=root)
        protected[name] = hashlib.sha256(actual).hexdigest()
    tables = {}
    for path in sorted(args.structure.glob("*.tsv")):
        if not path.with_suffix(path.suffix + ".provenance.txt").exists():
            continue
        with path.open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        tables[path.name] = {"row_count": len(rows), "lower": rows[0], "upper": rows[-1],
                             "parameters": path.with_suffix(path.suffix + ".provenance.txt").read_text().splitlines()}
    result = {
        "role": "GENERATED_REPRODUCIBILITY_RECORD_NOT_SCIENTIFIC_BASELINE",
        "generator": "structure1-candidate-v1",
        "starting_commit": "f4ae22d971e25bdd74530aec184f3fe0c3440b95",
        "local_provider_ratification_commit": "933494d86daf2cf8965079ece49fabd66d9390e5",
        "collection_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "collection_worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=root)),
        "source_sha256": {name: sha(root / name) for name in sources},
        "artifact_directories": {"structure": str(args.structure.resolve()), "local": str(args.local.resolve()), "regression": str(args.regression.resolve())},
        "artifacts": artifacts,
        "midpoint": midpoint,
        "tables": tables,
        "rho_bin_image": {k: [min(q[k] for q in sequence), max(q[k] for q in sequence)] for k in ("M", "R0_est", "Rinf0_profile")},
        "midpoint_radial_spread": {k: max(q[k] for q in finest) - min(q[k] for q in finest) for k in ("M", "R_cut", "lapse_profile")},
        "finest_midpoint_oracle_disagreement": {"M": max(abs(q["M"] - q["oracle_M_cut"]) for q in finest), "R_cut_km": max(abs(q["R_cut"] - q["oracle_R_cut"]) for q in finest)},
        "final_central_density_max_error_g_cm3": max(abs(q["analytic_rho"] - q["requested_rho"]) for q in stars if q["group"] != "radial" or q["grid"] == 8192),
        "surface_sweep": floors,
        "high_branch": branch,
        "local_validation": local,
        "structure_summary": (args.structure / "summary.txt").read_text().splitlines(),
        "regression": regress,
        "baseline_sha256": baseline,
        "protected_sha256": protected,
        "new_baseline": False,
        "canonical_TOV_behavior_changed": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
