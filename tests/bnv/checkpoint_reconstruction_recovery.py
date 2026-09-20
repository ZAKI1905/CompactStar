#!/usr/bin/env python3
"""Recovery-only authentication, staging and accounting for the Phase-6 checkpoint
reconstruction candidate campaign.

The failed validation attempt (8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead) remains immutable.
This script never integrates an ODE. It authenticates the reused 478 rk8pd oracle solves from
stored bytes only, checks the corrected profile/EOS binding, verifies each replay stage, and
enforces the 956-solve accounting. Scientific adjudication is delegated unchanged to
checkpoint_reconstruction_verify.py `oracle` and `final`.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import struct
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

FAILED_VALIDATION_SHA = "8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead"
PREFLIGHT_SHA = "93e93c7f91a3cd8fced2f7a0961eda9c469c43fe"
CANONICAL_SHA = "bd697ffdc474863d7a39f42e17ad8e8dbf105e5d"
MATRIX_SHA = "32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616"
ORACLE_RESULT_SHA = "f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0"
FAILED_HARNESS_SHA = "82eea47fba62f1440a1b87abec95dc4ac6287756ad44b990564794b9cf703044"
EXPECTED_MAX_DO = 0.5794839113173227
EXPECTED_MAX_UO = 0.5793505315921852
EXPECTED_MAX_DIAG = 0.17336139714051704
PASSIVE_HASHES = {
    "schedule": "43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e",
    "accepted": "7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459",
    "brackets": "0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d",
    "internal": "fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8",
    "trajectory": "8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c",
    "steps": "912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb",
    "library": "b9b767dbc0114563e1d556e296b6d7fc9d680a9d90e8deae9b44357010dd6499",
}
# Production qualification literals (CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:89-92).
PROFILE_ROOT_FILES = {
    "profile.tsv": "e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890",
    "model.txt": "3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4",
    "freegas.tsv": "7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a",
}
PROFILE_ROOT_TREE_SHA = "233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72"
SCIENTIFIC_INPUT_HASHES = {
    "certificate": "7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453",
    "thermal": "1cdb98507958b50fbd1e1eff3eb3bfa4f48c790fef86a8740855884fc43b9d94",
    "frozen": "9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab",
    "coefficients": "3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6",
    "entry": "5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d",
    "pretrajectory": "1d4e77780fde474f24782879d6ac43fb5bd6d01f3919c359437c0c5f01d327a6",
}
CANDIDATE_STAGES = ("linear", "linear-repeat", "hermite", "hermite-repeat",
                    "replay1", "replay2", "replay1-repeat", "replay2-repeat")
INTEGRATION_STAGES = ("replay1", "replay2", "replay1-repeat", "replay2-repeat")
AUTHORIZED_NEW_INTEGRATIONS = 956
REUSED_ORACLE_INTEGRATIONS = 478
IDENTITY = (
    "run_card_identity", "partition_identity", "source_identity", "domain_identity",
    "revision_identity", "product_fate_identity", "actual_potential_provenance",
    "finite_T_weighting_class",
)
STATE = ("x_state", "eta_e_MeV", "eta_mu_MeV")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(path: Path) -> str:
    path = Path(path)
    if path.is_file():
        return sha256(path)
    digest = hashlib.sha256()
    for item in sorted(item for item in path.rglob("*") if item.is_file()):
        relative = item.relative_to(path).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(bytes.fromhex(sha256(item)))
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def bits(value: str | float) -> bytes:
    return struct.pack(">d", float(value))


def command(args: list[str], cwd: Path | None = None) -> str:
    return subprocess.run(args, cwd=cwd, check=True, text=True, capture_output=True).stdout.strip()


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    print(json.dumps(value, indent=2, sort_keys=True))


def read_matrix(path: Path) -> list[dict]:
    parsed = []
    for raw in rows(path):
        parsed.append({
            **raw,
            "observation": int(raw["observation_index"]),
            "left_index": int(raw["left_endpoint_index"]),
            "right_index": int(raw["right_endpoint_index"]),
            "category": raw["cstar_category"],
            "deep": raw["deep_interior"] == "1",
            "exact": raw["exact_endpoint"] == "1",
            "strict": raw["strict_interior"] == "1",
            "local_integrations": int(raw["local_integrations"]),
        })
    if len(parsed) != 240:
        raise RuntimeError("solve matrix row count changed")
    return parsed


def matrix_summary(matrix: list[dict]) -> dict:
    return {
        "positive_observations": len(matrix),
        "no_knot": sum(row["category"] == "A" for row in matrix),
        "one_knot": sum(row["category"] == "B" for row in matrix),
        "multiple_knots": sum(row["category"] == "C" for row in matrix),
        "exact_endpoint": sum(row["exact"] for row in matrix),
        "strict_interior": sum(row["strict"] for row in matrix),
        "deep_interior": sum(row["deep"] for row in matrix),
        "authorized_local_integrations": sum(row["local_integrations"] for row in matrix),
    }


EXPECTED_SUMMARY = {
    "positive_observations": 240, "no_knot": 237, "one_knot": 3, "multiple_knots": 0,
    "exact_endpoint": 1, "strict_interior": 239, "deep_interior": 81,
    "authorized_local_integrations": 1434,
}


def expected_platform() -> dict:
    compiler = command(["/usr/bin/clang++", "--version"]).splitlines()
    return {
        "macos_product_version": command(["sw_vers", "-productVersion"]),
        "macos_build_version": command(["sw_vers", "-buildVersion"]),
        "darwin_release": platform.release(),
        "architecture": platform.machine(),
        "compiler": compiler[0],
        "compiler_target": next(line for line in compiler if line.startswith("Target:")),
        "gsl": command(["/opt/local/bin/gsl-config", "--version"]),
    }


PLATFORM_EXPECTED = {
    "macos_product_version": "26.6.2", "macos_build_version": "25G83",
    "darwin_release": "25.6.0", "architecture": "arm64",
    "compiler": "Apple clang version 21.0.0 (clang-2100.3.34.2)",
    "compiler_target": "Target: arm64-apple-darwin25.6.0", "gsl": "2.7.1",
}


def solve_id(method: str, observation: int) -> str:
    return f"{method}-obs-{observation:03d}"


def authorized_solve_ids(matrix: list[dict]) -> list[str]:
    return [solve_id(method, row["observation"]) for method in INTEGRATION_STAGES
            for row in matrix if row["strict"]]


def method_batch_audit(root: Path, method: str, matrix: list[dict], integrate: bool) -> dict:
    """Authenticate one complete method batch against the solve matrix and its own ledger."""
    method_root = root / method
    problems = []
    if not (method_root / "COMPLETE").exists():
        problems.append("COMPLETE marker missing")
    if (root / f"{method}.error").exists():
        problems.append("method error file present")
    names = sorted(item.name for item in method_root.iterdir()) if method_root.exists() else []
    expected_names = sorted([f"obs-{i:03d}.tsv" for i in range(241)] + [f"obs-{i:03d}.meta.tsv" for i in range(241)]
                            + ["COMPLETE"] + (["solve_ledger.tsv"] if integrate else []))
    if names != expected_names:
        problems.append("method output manifest differs from the expected 241 result + 241 meta files")
    ledger_rows: list[dict[str, str]] = []
    ledger_hash = None
    if integrate and (method_root / "solve_ledger.tsv").exists():
        ledger_rows = rows(method_root / "solve_ledger.tsv")
        ledger_hash = sha256(method_root / "solve_ledger.tsv")
        ids = [row["solve_id"] for row in ledger_rows]
        if len(ids) != 239 or len(set(ids)) != 239:
            problems.append(f"ledger has {len(ids)} rows / {len(set(ids))} unique ids, expected 239")
        expected_ids = [solve_id(method, row["observation"]) for row in matrix if row["strict"]]
        if ids != expected_ids:
            problems.append("ledger solve ids differ from the strict-interior matrix rows")
        by_observation = {int(row["observation_index"]): row for row in ledger_rows}
        for row in matrix:
            if not row["strict"]:
                continue
            entry = by_observation.get(row["observation"])
            if entry is None:
                problems.append(f"observation {row['observation']} missing from ledger")
                continue
            if entry["method"] != method or int(entry["left_endpoint_index"]) != row["left_index"]:
                problems.append(f"observation {row['observation']} ledger method/left endpoint differ")
            if bits(entry["t_left_s"]) != bits(row["t_left_s"]) or bits(entry["target_s"]) != bits(row["t_obs_s"]):
                problems.append(f"observation {row['observation']} ledger left/target time differ from matrix")
            if entry["exit_status"] != "0":
                problems.append(f"observation {row['observation']} nonzero exit status")
            if sha256(method_root / f"obs-{row['observation']:03d}.tsv") != entry["result_sha256"]:
                problems.append(f"observation {row['observation']} result hash differs from ledger")
            if sha256(method_root / f"obs-{row['observation']:03d}.meta.tsv") != entry["meta_sha256"]:
                problems.append(f"observation {row['observation']} meta hash differs from ledger")
    identities = None
    for index in range(241):
        result_file = method_root / f"obs-{index:03d}.tsv"
        meta_file = method_root / f"obs-{index:03d}.meta.tsv"
        if not result_file.exists() or not meta_file.exists():
            continue
        data = rows(result_file)
        meta = rows(meta_file)
        if len(data) != 1 or len(meta) != 1:
            problems.append(f"observation {index} row count changed")
            continue
        row, m = data[0], meta[0]
        identity = tuple(row[name] for name in IDENTITY)
        if identities is None:
            identities = identity
        elif identity != identities:
            problems.append(f"observation {index} identity fields differ within the batch")
        if row["valid_through_sample"] != "1":
            problems.append(f"observation {index} not valid_through_sample")
        for name in STATE:
            if not math.isfinite(float(row[name])):
                problems.append(f"observation {index} nonfinite state {name}")
        if m["method"] != method or int(m["observation_index"]) != index:
            problems.append(f"observation {index} meta method/index differ")
        if index:
            matrix_row = matrix[index - 1]
            if bits(m["t_obs_s"]) != bits(matrix_row["t_obs_s"]) or bits(m["t_left_s"]) != bits(matrix_row["t_left_s"]) \
                    or bits(m["t_right_s"]) != bits(matrix_row["t_right_s"]):
                problems.append(f"observation {index} meta times differ from matrix")
            if m["category"] != matrix_row["category"] or (m["strict"] == "1") != matrix_row["strict"] \
                    or (m["exact"] == "1") != matrix_row["exact"]:
                problems.append(f"observation {index} meta category/flags differ from matrix")
            if integrate and matrix_row["strict"] and int(m["accepted"]) <= 0:
                problems.append(f"observation {index} integration recorded no accepted step")
            if integrate and not matrix_row["strict"] and int(m["accepted"]) != 0:
                problems.append(f"observation {index} exact endpoint consumed an integration")
        elif int(m["accepted"]) != 0:
            problems.append("initial observation consumed an integration")
    return {
        "method": method, "pass": not problems, "problems": problems,
        "ledger_rows": len(ledger_rows), "ledger_sha256": ledger_hash,
        "identity": dict(zip(IDENTITY, identities)) if identities else None,
        "tree_sha256": tree_sha256(method_root) if method_root.exists() else None,
    }


def scientific_inputs(args) -> dict:
    profile_root = Path(args.profile_root)
    inputs = {
        "profile_root": {
            "path": str(profile_root), "is_directory": profile_root.is_dir(),
            "files": {name: (sha256(profile_root / name) if (profile_root / name).is_file() else None)
                      for name in PROFILE_ROOT_FILES},
            "tree_sha256": tree_sha256(profile_root) if profile_root.is_dir() else None,
        },
    }
    for name in ("certificate", "thermal", "frozen", "coefficients", "entry", "pretrajectory"):
        inputs[name] = {"path": str(getattr(args, name)), "sha256": tree_sha256(Path(getattr(args, name)))}
    inputs["profile_root"]["pass"] = (inputs["profile_root"]["is_directory"]
                                      and inputs["profile_root"]["files"] == PROFILE_ROOT_FILES
                                      and inputs["profile_root"]["tree_sha256"] == PROFILE_ROOT_TREE_SHA)
    inputs["pass"] = inputs["profile_root"]["pass"] and all(
        inputs[name]["sha256"] == SCIENTIFIC_INPUT_HASHES[name] for name in SCIENTIFIC_INPUT_HASHES)
    return inputs


def passive_hashes(args) -> dict:
    return {name: sha256(Path(getattr(args, name))) for name in PASSIVE_HASHES}


def authenticate(args) -> None:
    repo = Path(args.repo)
    matrix = read_matrix(args.matrix)
    summary = matrix_summary(matrix)
    hashes = passive_hashes(args)
    inputs = scientific_inputs(args)
    forbidden = command([
        "git", "diff", "--name-only", PREFLIGHT_SHA, "--", "CompactStar", "tests/baselines",
        "CMakeLists.txt", "EOS", "data", "literature", "dependencies", "tests/rotochemical", "tests/analysis", "tests/eos",
    ], repo).splitlines()
    oracle_root = Path(args.oracle_root)
    oracle = {name: method_batch_audit(oracle_root, name, matrix, True) for name in ("oracle1", "oracle2")}
    oracle_ids = []
    for name in ("oracle1", "oracle2"):
        if (oracle_root / name / "solve_ledger.tsv").exists():
            oracle_ids.extend(row["solve_id"] for row in rows(oracle_root / name / "solve_ledger.tsv"))
    oracle_result = Path(args.oracle_result)
    platform_now = expected_platform()
    head = command(["git", "rev-parse", "HEAD"], repo)
    checks = {
        "matrix_hash": sha256(args.matrix) == MATRIX_SHA,
        "matrix_counts": summary == EXPECTED_SUMMARY,
        "passive_hashes": all(hashes[name] == PASSIVE_HASHES[name] for name in PASSIVE_HASHES),
        "scientific_inputs": inputs["pass"],
        "forbidden_paths_unchanged": forbidden == [],
        "platform_match": platform_now == PLATFORM_EXPECTED,
        "canonical_ref": command(["git", "rev-parse", "master"], repo) == CANONICAL_SHA
                         and command(["git", "rev-parse", "origin/master"], repo) == CANONICAL_SHA,
        "failed_validation_is_ancestor": subprocess.run(
            ["git", "merge-base", "--is-ancestor", FAILED_VALIDATION_SHA, "HEAD"], cwd=repo).returncode == 0,
        "oracle_batches": all(value["pass"] for value in oracle.values()),
        "oracle_ledger_count": len(oracle_ids) == REUSED_ORACLE_INTEGRATIONS and len(set(oracle_ids)) == REUSED_ORACLE_INTEGRATIONS,
        "oracle_phase_complete": (oracle_root / "PHASE_COMPLETE").exists(),
        "oracle_result_hash": oracle_result.exists() and sha256(oracle_result) == ORACLE_RESULT_SHA,
        "failed_harness_hash": Path(args.failed_harness).exists() and sha256(Path(args.failed_harness)) == FAILED_HARNESS_SHA,
    }
    result = {
        "classification": "CHECKPOINT-RECONSTRUCTION RECOVERY AUTHENTICATION; NO INTEGRATION",
        "head": head, "failed_validation_sha": FAILED_VALIDATION_SHA, "preflight_sha": PREFLIGHT_SHA,
        "matrix_sha256": sha256(args.matrix), "matrix_summary": summary,
        "passive_hashes": hashes, "scientific_inputs": inputs,
        "forbidden_path_diff": forbidden, "platform": platform_now,
        "oracle": oracle, "oracle_root_tree_sha256": tree_sha256(oracle_root),
        "oracle_result_sha256": sha256(oracle_result) if oracle_result.exists() else None,
        "reused_oracle_integrations": len(oracle_ids),
        "authorized_new_integrations": AUTHORIZED_NEW_INTEGRATIONS,
        "authorized_solve_ids_sha256": hashlib.sha256("\n".join(authorized_solve_ids(matrix)).encode() + b"\n").hexdigest(),
        "checks": checks, "pass": all(checks.values()), "main_trajectory_integrations": 0,
        "new_integrations_executed": 0,
    }
    if args.authorized_output:
        Path(args.authorized_output).write_text("\n".join(authorized_solve_ids(matrix)) + "\n")
    write_json(args.output, result)
    if not result["pass"]:
        raise SystemExit(2)


def requalify(args) -> None:
    """Recompute the oracle self-qualification from stored bytes through the unchanged verifier."""
    verifier = Path(args.verifier)
    subprocess.run([sys.executable, str(verifier), "oracle", "--matrix", str(args.matrix), "--root", str(args.oracle_root),
                    "--output", str(args.output), "--flag", str(args.flag)], check=True, text=True, capture_output=True)
    data = json.loads(Path(args.output).read_text())
    checks = {
        "result_pass": data["pass"] is True,
        "executed_local_integrations": data["executed_local_integrations"] == REUSED_ORACLE_INTEGRATIONS,
        "result_sha256_reproduced": sha256(Path(args.output)) == ORACLE_RESULT_SHA,
        "max_d_O_over_D_O1": bits(data["maximum_d_O_over_D_O1"]["value"]) == bits(EXPECTED_MAX_DO)
                              and data["maximum_d_O_over_D_O1"]["observation"] == 117
                              and data["maximum_d_O_over_D_O1"]["component"] == "x_state"
                              and data["maximum_d_O_over_D_O1"]["category"] == "B"
                              and data["maximum_d_O_over_D_O1"]["deep"] is False,
        "max_U_O_over_0p20F": bits(data["maximum_U_O_over_0p20F"]["value"]) == bits(EXPECTED_MAX_UO)
                               and data["maximum_U_O_over_0p20F"]["observation"] == 117,
        "max_diagnostic": bits(data["maximum_diagnostic_U_O_over_0p20F_P"]["value"]) == bits(EXPECTED_MAX_DIAG)
                          and data["maximum_diagnostic_U_O_over_0p20F_P"]["observable"] == "Pnet_erg_s"
                          and data["maximum_diagnostic_U_O_over_0p20F_P"]["observation"] == 117,
        "flag_written": Path(args.flag).exists() and "ORACLE QUALIFIED" in Path(args.flag).read_text(),
    }
    summary = {
        "classification": "ORACLE SELF-QUALIFICATION RECOMPUTED FROM STORED BYTES; NO INTEGRATION",
        "recomputed_result_sha256": sha256(Path(args.output)), "expected_result_sha256": ORACLE_RESULT_SHA,
        "maxima": {key: data[key] for key in ("maximum_d_O_over_D_O1", "maximum_U_O_over_0p20F", "maximum_diagnostic_U_O_over_0p20F_P")},
        "checks": checks, "pass": all(checks.values()),
    }
    write_json(Path(args.summary), summary)
    if not summary["pass"]:
        raise SystemExit(2)


def stage_check(args) -> None:
    matrix = read_matrix(args.matrix)
    method = args.method
    integrate = method in INTEGRATION_STAGES
    audit = method_batch_audit(Path(args.candidate_root), method, matrix, integrate)
    oracle2 = method_batch_audit(Path(args.oracle_root), "oracle2", matrix, True)
    ledger_rows = rows(Path(args.execution_ledger)) if Path(args.execution_ledger).exists() else []
    ids = [row["solve_id"] for row in ledger_rows]
    stage_ids = [row["solve_id"] for row in ledger_rows if row["method"] == method]
    checks = {
        "batch": audit["pass"],
        "identity_matches_oracle": audit["identity"] == oracle2["identity"],
        "execution_ledger_unique": len(ids) == len(set(ids)),
        "execution_ledger_within_cap": len(ids) <= AUTHORIZED_NEW_INTEGRATIONS,
        "stage_ledger_count": len(stage_ids) == (239 if integrate else 0),
        "stage_ledger_matches_method_ledger": (not integrate) or (
            stage_ids == [row["solve_id"] for row in rows(Path(args.candidate_root) / method / "solve_ledger.tsv")]),
    }
    result = {"classification": "REPLAY STAGE CHECK", "stage": method, "audit": audit,
              "execution_ledger_rows": len(ids), "checks": checks, "pass": all(checks.values())}
    write_json(Path(args.output), result)
    if not result["pass"]:
        raise SystemExit(2)


def determinism(args) -> None:
    matrix = read_matrix(args.matrix)
    root = Path(args.candidate_root)
    report = {}
    overall = True
    for method, repeat in (("replay1", "replay1-repeat"), ("replay2", "replay2-repeat")):
        unequal = []
        max_difference = 0.0
        hashes_equal = 0
        for index in range(241):
            a, b = root / method / f"obs-{index:03d}.tsv", root / repeat / f"obs-{index:03d}.tsv"
            if a.read_bytes() == b.read_bytes():
                hashes_equal += 1
            else:
                unequal.append(index)
                ra, rb = rows(a)[0], rows(b)[0]
                for name in STATE:
                    max_difference = max(max_difference, abs(float(ra[name]) - float(rb[name])))
        ledger_a = [row["result_sha256"] for row in rows(root / method / "solve_ledger.tsv")]
        ledger_b = [row["result_sha256"] for row in rows(root / repeat / "solve_ledger.tsv")]
        steps_unequal = [index for index in range(241) if any(
            rows(root / method / f"obs-{index:03d}.meta.tsv")[0][k] != rows(root / repeat / f"obs-{index:03d}.meta.tsv")[0][k]
            for k in ("accepted", "rejected", "rhs"))]
        passed = not unequal and ledger_a == ledger_b and not steps_unequal
        overall &= passed
        report[method] = {"comparisons": 241, "strict_interior_comparisons": 239, "unequal_result_count": len(unequal),
                          "unequal_observations": unequal, "unequal_step_count_observations": steps_unequal,
                          "max_state_difference": max_difference, "byte_identical_count": hashes_equal,
                          "ledger_result_hashes_identical": ledger_a == ledger_b, "pass": passed}
    result = {"classification": "REPLAY DETERMINISM", "tiers": report, "pass": overall}
    write_json(Path(args.output), result)
    if not overall:
        raise SystemExit(2)


def category_maxima(final: dict, matrix: list[dict], method: str) -> dict:
    observations = final["methods"][method]["observations"]
    out = {}
    for label, predicate in (("no_knot", lambda r: r["category"] == "A"), ("one_knot", lambda r: r["category"] == "B"),
                             ("exact_endpoint", lambda r: r["exact"]), ("strict_interior", lambda r: r["strict"]),
                             ("deep_interior", lambda r: r["deep"]), ("all", lambda r: True)):
        selected = [observations[str(r["observation"])] for r in matrix if predicate(r)]
        state_values = [v for item in selected for v in item["state"].values()]
        ledger_values = [v for item in selected for v in item["ledger"].values()]
        out[label] = {"count": len(selected), "max_state_utilization": max(state_values) if state_values else None,
                      "max_ledger_utilization": max(ledger_values) if ledger_values else None,
                      "failures": sum(not item["pass"] for item in selected)}
    worst_ledger = max(((v, name, index) for index, item in observations.items() for name, v in item["ledger"].items()),
                       key=lambda x: x[0])
    out["worst_ledger"] = {"utilization": worst_ledger[0], "observable": worst_ledger[1], "observation": int(worst_ledger[2]),
                           "t_obs_s": matrix[int(worst_ledger[2]) - 1]["t_obs_s"]}
    return out


def summarize(args) -> None:
    matrix = read_matrix(args.matrix)
    final = json.loads(Path(args.final).read_text())
    ledger_rows = rows(Path(args.execution_ledger))
    ids = [row["solve_id"] for row in ledger_rows]
    authorized = set(authorized_solve_ids(matrix))
    per_stage = {method: sum(row["method"] == method for row in ledger_rows) for method in INTEGRATION_STAGES}
    accounting = {
        "reused_oracle_integrations": REUSED_ORACLE_INTEGRATIONS,
        "new_integrations": len(ids), "new_integrations_unique": len(set(ids)),
        "new_integrations_authorized": all(item in authorized for item in ids) and len(set(ids)) == AUTHORIZED_NEW_INTEGRATIONS,
        "per_stage": per_stage, "complete_campaign": REUSED_ORACLE_INTEGRATIONS + len(ids),
        "verifier_executed_local_integrations": final["executed_local_integrations"],
        "pass": len(ids) == AUTHORIZED_NEW_INTEGRATIONS and len(set(ids)) == AUTHORIZED_NEW_INTEGRATIONS
                and final["executed_local_integrations"] == 1434 and all(v == 239 for v in per_stage.values()),
    }
    methods = {name: category_maxima(final, matrix, name) for name in ("linear", "hermite", "replay2")}
    for name in ("linear", "hermite", "replay2"):
        methods[name]["R20"] = final["methods"][name]["R20"]["metrics"]
        methods[name]["R20_pass"] = final["methods"][name]["R20"]["pass"]
        methods[name]["pass"] = final["methods"][name]["pass"]
        methods[name]["worst_utilization"] = final["methods"][name]["worst_utilization"]
        methods[name]["failure_count"] = final["methods"][name]["failure_count"]
    perf = final["performance"]
    replay_wall = sum(perf[name]["wall_sum_s"] for name in INTEGRATION_STAGES)
    replay_cpu = sum(perf[name]["cpu_user_sum_s"] + perf[name]["cpu_sys_sum_s"] for name in INTEGRATION_STAGES)
    all_wall = [float(rows(Path(args.candidate_root) / name / f"obs-{i:03d}.meta.tsv")[0]["solve_wall_s"])
                for name in INTEGRATION_STAGES for i in range(1, 241) if matrix[i - 1]["strict"]]
    ordered = sorted(all_wall)
    result = {
        "classification": "PHASE-6 CHECKPOINT-RECONSTRUCTION RECOVERY VALIDATION EVIDENCE; NOT BNV CANDIDATE; NOT GOVERNED BASELINE; NOT PHYSICAL RESULT",
        "accounting": accounting, "methods": methods,
        "hybrid_eligible": final["hybrid_eligible"], "hybrid": final["hybrid"],
        "selected_method": final["selected_method"], "selection_reason": final["selection_reason"],
        "replay_witness_max_utilization": final["methods"]["replay2"]["replay_witness_maximum_utilization"],
        "replay_performance": {
            "wall_sum_s": replay_wall, "cpu_sum_s": replay_cpu, "solves": len(all_wall),
            "median_solve_wall_s": statistics.median(all_wall), "p95_solve_wall_s": ordered[max(0, math.ceil(0.95 * len(ordered)) - 1)],
            "max_solve_wall_s": max(all_wall), "throughput_solves_per_wall_s": len(all_wall) / replay_wall,
            "per_stage": {name: perf[name] for name in INTEGRATION_STAGES},
            "estimated_8191_strict_interior_solve_wall_s_per_tier": {
                name: perf[name]["wall_sum_s"] * 8191 / 239 for name in ("replay1", "replay2")},
        },
        "final_pass": final["pass"], "immutable_inputs_pass": final["immutable_inputs_pass"],
        "historical_BA12": "FAIL", "historical_BA12R": "FAIL", "main_trajectory_integrations": 0,
        "bnv_candidate_created": False, "production_adr_created": False,
        "pass": final["pass"] and accounting["pass"],
    }
    write_json(Path(args.output), result)
    if not result["pass"]:
        raise SystemExit(2)


parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest="command", required=True)

p = sub.add_parser("authenticate")
p.add_argument("--repo", type=Path, required=True)
p.add_argument("--matrix", type=Path, required=True)
for name in ("schedule", "accepted", "brackets", "internal", "trajectory", "steps", "library"):
    p.add_argument(f"--{name}", dest=name, type=Path, required=True)
p.add_argument("--profile-root", dest="profile_root", type=Path, required=True)
for name in ("certificate", "thermal", "frozen", "coefficients", "entry", "pretrajectory"):
    p.add_argument(f"--{name}", dest=name, type=Path, required=True)
p.add_argument("--oracle-root", dest="oracle_root", type=Path, required=True)
p.add_argument("--oracle-result", dest="oracle_result", type=Path, required=True)
p.add_argument("--failed-harness", dest="failed_harness", type=Path, required=True)
p.add_argument("--authorized-output", dest="authorized_output", type=Path)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=authenticate)

p = sub.add_parser("requalify")
p.add_argument("--verifier", type=Path, required=True)
p.add_argument("--matrix", type=Path, required=True)
p.add_argument("--oracle-root", dest="oracle_root", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.add_argument("--flag", type=Path, required=True)
p.add_argument("--summary", type=Path, required=True)
p.set_defaults(function=requalify)

p = sub.add_parser("stage-check")
p.add_argument("--matrix", type=Path, required=True)
p.add_argument("--method", required=True, choices=CANDIDATE_STAGES)
p.add_argument("--candidate-root", dest="candidate_root", type=Path, required=True)
p.add_argument("--oracle-root", dest="oracle_root", type=Path, required=True)
p.add_argument("--execution-ledger", dest="execution_ledger", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=stage_check)

p = sub.add_parser("determinism")
p.add_argument("--matrix", type=Path, required=True)
p.add_argument("--candidate-root", dest="candidate_root", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=determinism)

p = sub.add_parser("summarize")
p.add_argument("--matrix", type=Path, required=True)
p.add_argument("--final", type=Path, required=True)
p.add_argument("--candidate-root", dest="candidate_root", type=Path, required=True)
p.add_argument("--execution-ledger", dest="execution_ledger", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=summarize)

arguments = parser.parse_args()
arguments.function(arguments)
