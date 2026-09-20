#!/usr/bin/env python3
"""Authenticate and adjudicate the bounded Phase-6 checkpoint reconstruction experiment."""

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

PREFLIGHT_SHA = "93e93c7f91a3cd8fced2f7a0961eda9c469c43fe"
PREDECLARATION_SHA = "94788a2f0941eb9f9510ffb810c7f89b9e3ae26d"
CANONICAL_SHA = "bd697ffdc474863d7a39f42e17ad8e8dbf105e5d"
MATRIX_SHA = "32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616"
PASSIVE_HASHES = {
    "schedule": "43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e",
    "accepted": "7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459",
    "brackets": "0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d",
    "internal": "fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8",
    "trajectory": "8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c",
    "steps": "912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb",
    "library": "b9b767dbc0114563e1d556e296b6d7fc9d680a9d90e8deae9b44357010dd6499",
}
STATE = ("x_state", "eta_e_MeV", "eta_mu_MeV")
ATOL_ULTRA = (1e-16, 1e-22, 1e-22)
ATOL_ORACLE1 = (1e-17, 1e-23, 1e-23)
ATOL_ORACLE2 = (1e-18, 1e-24, 1e-24)
LEDGER = (
    "P_dir_eq_erg_s", "P_dir_actual_erg_s", "LH_erg_s", "DeltaLnu_erg_s",
    "DeltaPbeta_erg_s", "Lnu_eq_erg_s", "Lnu_full_erg_s",
    "L_out_fluid_inf_erg_s", "Lgamma_erg_s", "Lother_erg_s", "Pnet_erg_s",
)
GROSS = (
    "P_dir_actual_erg_s", "LH_erg_s", "DeltaLnu_erg_s", "Lnu_eq_erg_s",
    "Lgamma_erg_s", "Lother_erg_s", "L_out_fluid_inf_erg_s",
)
IDENTITY = (
    "run_card_identity", "partition_identity", "source_identity", "domain_identity",
    "revision_identity", "product_fate_identity", "actual_potential_provenance",
    "finite_T_weighting_class",
)
MEV_TO_ERG = 1.602176634e-6
EPS = sys.float_info.epsilon


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def finite(value: str | float, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"nonfinite {label}")
    return result


def bits(value: str | float) -> bytes:
    return struct.pack(">d", float(value))


def exact(a: str | float, b: str | float) -> bool:
    return bits(a) == bits(b)


def command(args: list[str], cwd: Path | None = None) -> str:
    return subprocess.run(args, cwd=cwd, check=True, text=True, capture_output=True).stdout.strip()


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    print(json.dumps(value, indent=2, sort_keys=True))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    return ordered[max(0, math.ceil(fraction * len(ordered)) - 1)]


def read_matrix(path: Path) -> list[dict]:
    parsed = []
    for raw in rows(path):
        parsed.append({
            **raw,
            "observation": int(raw["observation_index"]),
            "left_index": int(raw["left_endpoint_index"]),
            "right_index": int(raw["right_endpoint_index"]),
            "t_obs": finite(raw["t_obs_s"], "matrix time"),
            "t_left": finite(raw["t_left_s"], "matrix left time"),
            "t_right": finite(raw["t_right_s"], "matrix right time"),
            "left": tuple(finite(raw[name], name) for name in ("x_left", "eta_e_left_MeV", "eta_mu_left_MeV")),
            "right": tuple(finite(raw[name], name) for name in ("x_right", "eta_e_right_MeV", "eta_mu_right_MeV")),
            "category": raw["cstar_category"],
            "deep": raw["deep_interior"] == "1",
            "exact": raw["exact_endpoint"] == "1",
            "strict": raw["strict_interior"] == "1",
            "local_integrations": int(raw["local_integrations"]),
        })
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


def preflight(args) -> None:
    matrix = read_matrix(args.matrix)
    summary = matrix_summary(matrix)
    hashes = {
        "matrix": sha256(args.matrix), "schedule": sha256(args.schedule),
        "accepted": sha256(args.accepted), "brackets": sha256(args.brackets),
        "internal": sha256(args.internal), "trajectory": sha256(args.trajectory),
        "steps": sha256(args.steps), "library": sha256(args.production_library),
    }
    expected_summary = {
        "positive_observations": 240, "no_knot": 237, "one_knot": 3,
        "multiple_knots": 0, "exact_endpoint": 1, "strict_interior": 239,
        "deep_interior": 81, "authorized_local_integrations": 1434,
    }
    repo = args.repo
    predeclared_unchanged = command([
        "git", "diff", "--name-only", PREDECLARATION_SHA, "--",
        "docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_VALIDATION.md",
        "docs/validation/phase6a1_checkpoint_reconstruction_solve_matrix.tsv",
    ], repo) == ""
    forbidden = command([
        "git", "diff", "--name-only", PREFLIGHT_SHA, "--", "CompactStar", "tests/baselines",
        "CMakeLists.txt", "EOS", "data", "literature",
    ], repo).splitlines()
    platform_now = expected_platform()
    platform_expected = {
        "macos_product_version": "26.6.2", "macos_build_version": "25G83",
        "darwin_release": "25.6.0", "architecture": "arm64",
        "compiler": "Apple clang version 21.0.0 (clang-2100.3.34.2)",
        "compiler_target": "Target: arm64-apple-darwin25.6.0", "gsl": "2.7.1",
    }
    checks = {
        "matrix_hash": hashes["matrix"] == MATRIX_SHA,
        "passive_hashes": all(hashes[name] == PASSIVE_HASHES[name] for name in PASSIVE_HASHES),
        "matrix_counts": summary == expected_summary,
        "predeclaration_commit": command(["git", "rev-parse", PREDECLARATION_SHA], repo) == PREDECLARATION_SHA,
        "predeclaration_parent": command(["git", "rev-parse", PREDECLARATION_SHA + "^"], repo) == PREFLIGHT_SHA,
        "predeclaration_unchanged": predeclared_unchanged,
        "forbidden_paths_unchanged": forbidden == [],
        "platform_match": platform_now == platform_expected,
        "canonical_ref": command(["git", "rev-parse", "master"], repo) == CANONICAL_SHA and
                         command(["git", "rev-parse", "origin/master"], repo) == CANONICAL_SHA,
    }
    result = {
        "classification": "PRE-RUN CHECKPOINT-RECONSTRUCTION AUTHENTICATION",
        "hashes": hashes, "matrix_summary": summary, "planned_decomposition": {
            "oracle1": 239, "oracle2": 239, "replay1": 239, "replay2": 239,
            "replay1_repeat": 239, "replay2_repeat": 239, "total": 1434,
        },
        "platform": platform_now, "expected_platform": platform_expected,
        "forbidden_path_diff": forbidden, "checks": checks, "pass": all(checks.values()),
        "main_trajectory_integrations": 0,
    }
    write_json(args.output, result)
    if not result["pass"]:
        raise SystemExit(2)


def result_path(root: Path, method: str, observation: int) -> Path:
    return root / method / f"obs-{observation:03d}.tsv"


def meta_path(root: Path, method: str, observation: int) -> Path:
    return root / method / f"obs-{observation:03d}.meta.tsv"


def read_result(root: Path, method: str, observation: int) -> dict[str, str]:
    data = rows(result_path(root, method, observation))
    if len(data) != 1:
        raise RuntimeError(f"{method} observation {observation}: diagnostic row count changed")
    return data[0]


def read_meta(root: Path, method: str, observation: int) -> dict[str, str]:
    data = rows(meta_path(root, method, observation))
    if len(data) != 1:
        raise RuntimeError(f"{method} observation {observation}: meta row count changed")
    return data[0]


def state(row: dict[str, str]) -> tuple[float, float, float]:
    return tuple(finite(row[name], name) for name in STATE)


def state_scales(matrix_row: dict, oracle1: tuple[float, ...], oracle2: tuple[float, ...], component: int) -> tuple[float, float, float, float]:
    magnitude = max(abs(matrix_row["left"][component]), abs(matrix_row["right"][component]),
                    abs(oracle1[component]), abs(oracle2[component]))
    d_ultra = ATOL_ULTRA[component] + 1e-11 * magnitude
    floor = max(d_ultra, 64 * math.ulp(magnitude))
    return magnitude, d_ultra, floor, 0.25 * floor


def gross_power(row: dict[str, str]) -> float:
    return max(1.0, sum(abs(finite(row[name], name)) for name in GROSS))


def oracle_diagnostic_uncertainty(o1: dict[str, str], o2: dict[str, str], name: str) -> tuple[float, float, float]:
    gp = gross_power(o2)
    a, b = finite(o1[name], name), finite(o2[name], name)
    magnitude = max(abs(a), abs(b))
    oracle_floor = max(1e-13 * gp, 64 * math.ulp(magnitude))
    uncertainty = 2 * max(abs(a - b), oracle_floor)
    common_floor = max(1e-11 * gp, 64 * math.ulp(magnitude))
    return uncertainty, common_floor, 0.25 * common_floor


def validate_identity(row: dict[str, str], reference: dict[str, str]) -> bool:
    return (row["valid_through_sample"] == "1" and
            all(row[name] == reference[name] for name in IDENTITY) and
            all(math.isfinite(finite(row[name], name)) for name in LEDGER + STATE))


def ledger_rows(path: Path) -> list[dict[str, str]]:
    data = rows(path)
    if len(data) != 239:
        raise RuntimeError(f"{path}: expected 239 integration ledger rows, got {len(data)}")
    return data


def oracle(args) -> None:
    matrix = read_matrix(args.matrix)
    if not (args.root / "oracle1" / "COMPLETE").exists() or not (args.root / "oracle2" / "COMPLETE").exists():
        raise RuntimeError("oracle method batch incomplete")
    ledgers = ledger_rows(args.root / "oracle1" / "solve_ledger.tsv") + ledger_rows(args.root / "oracle2" / "solve_ledger.tsv")
    state_pass = True
    diagnostic_pass = True
    identity_pass = True
    max_do = (-1.0, None)
    max_uo = (-1.0, None)
    max_diag = (-1.0, None)
    per_observation = []
    for matrix_row in matrix:
        index = matrix_row["observation"]
        r1, r2 = read_result(args.root, "oracle1", index), read_result(args.root, "oracle2", index)
        y1, y2 = state(r1), state(r2)
        detail = {"observation": index, "category": matrix_row["category"], "deep": matrix_row["deep"], "components": {}}
        for component, name in enumerate(STATE):
            _, _, floor, _ = state_scales(matrix_row, y1, y2, component)
            d_o = abs(y2[component] - y1[component])
            d_o1 = ATOL_ORACLE1[component] + 1e-12 * max(abs(y1[component]), abs(y2[component]))
            d_o2 = ATOL_ORACLE2[component] + 1e-13 * max(abs(y1[component]), abs(y2[component]))
            f_o = max(d_o2, 64 * math.ulp(max(abs(y1[component]), abs(y2[component]))))
            u_o = 0.0 if matrix_row["exact"] else 2 * max(d_o, f_o)
            ratio_do = d_o / d_o1
            ratio_uo = u_o / (0.20 * floor)
            passed = d_o <= d_o1 and u_o <= 0.20 * floor
            state_pass &= passed
            item = {"d_O": d_o, "D_O1": d_o1, "F_i": floor, "U_O": u_o,
                    "d_O_over_D_O1": ratio_do, "U_O_over_0p20F": ratio_uo, "pass": passed}
            detail["components"][name] = item
            candidate = {"observation": index, "component": name, "category": matrix_row["category"], "deep": matrix_row["deep"]}
            max_do = max(max_do, (ratio_do, candidate), key=lambda value: value[0])
            max_uo = max(max_uo, (ratio_uo, candidate), key=lambda value: value[0])
        for name in LEDGER:
            uncertainty, floor, _ = oracle_diagnostic_uncertainty(r1, r2, name)
            ratio = uncertainty / (0.20 * floor)
            passed = uncertainty <= 0.20 * floor
            diagnostic_pass &= passed
            max_diag = max(max_diag, (ratio, {"observation": index, "observable": name,
                                              "category": matrix_row["category"], "deep": matrix_row["deep"]}), key=lambda value: value[0])
        identity_pass &= validate_identity(r1, r2) and validate_identity(r2, r2)
        per_observation.append(detail)
    passed = state_pass and diagnostic_pass and identity_pass and len(ledgers) == 478
    result = {
        "classification": "TWO-LEVEL LOCAL RK8PD ORACLE QUALIFICATION",
        "oracle1": {"method": "rk8pd", "rtol": 1e-12, "atol": ATOL_ORACLE1},
        "oracle2": {"method": "rk8pd", "rtol": 1e-13, "atol": ATOL_ORACLE2},
        "executed_local_integrations": len(ledgers),
        "state_self_qualification_pass": state_pass,
        "diagnostic_self_qualification_pass": diagnostic_pass,
        "identity_currentness_pass": identity_pass,
        "maximum_d_O_over_D_O1": {"value": max_do[0], **max_do[1]},
        "maximum_U_O_over_0p20F": {"value": max_uo[0], **max_uo[1]},
        "maximum_diagnostic_U_O_over_0p20F_P": {"value": max_diag[0], **max_diag[1]},
        "observations": per_observation,
        "pass": passed,
    }
    write_json(args.output, result)
    if passed:
        payload = f"ORACLE QUALIFIED\npredeclaration={PREDECLARATION_SHA}\nresult_sha256={sha256(args.output)}\n"
        temporary = args.flag.with_name(args.flag.name + f".tmp.{os.getpid()}")
        temporary.write_text(payload)
        temporary.replace(args.flag)
    else:
        raise SystemExit(2)


def trapezoid(sequence: list[dict[str, str]], value) -> float:
    total = 0.0
    for left, right in zip(sequence, sequence[1:]):
        dt = finite(right["t_s"], "time") - finite(left["t_s"], "time")
        if not dt > 0:
            raise RuntimeError("nonincreasing reconstructed time grid")
        total += 0.5 * dt * (value(left) + value(right))
    return total


def r20(sequence: list[dict[str, str]]) -> dict[str, float]:
    first, last = sequence[0], sequence[-1]
    delta_eeq = MEV_TO_ERG * finite(first["mu_B_inf_MeV"], "mu_B") * (
        finite(last["B_count"], "B final") - finite(first["B_count"], "B initial"))
    delta_echem = MEV_TO_ERG * (finite(last["Echem_MeV"], "Echem final") - finite(first["Echem_MeV"], "Echem initial"))
    delta_uth = sum(0.5 * (finite(a["Cstar_erg_K"], "Cstar") + finite(b["Cstar_erg_K"], "Cstar")) *
                    (finite(b["Tinf_K"], "Tinf") - finite(a["Tinf_K"], "Tinf"))
                    for a, b in zip(sequence, sequence[1:]))
    luminosity = lambda row: sum(finite(row[name], name) for name in
                                 ("L_out_fluid_inf_erg_s", "Lnu_full_erg_s", "Lgamma_erg_s", "Lother_erg_s"))
    outgoing = trapezoid(sequence, luminosity)
    residual = delta_eeq + delta_echem + delta_uth + outgoing
    normalizer = max(1.0, abs(delta_uth), abs(delta_echem), trapezoid(sequence, lambda row: sum(
        abs(finite(row[name], name)) for name in ("Lnu_full_erg_s", "Lgamma_erg_s", "Lother_erg_s"))))
    return {"DeltaEeq_erg": delta_eeq, "DeltaEchem_erg": delta_echem, "DeltaUth_erg": delta_uth,
            "outgoing_erg": outgoing, "endpoint_sum_erg": delta_eeq + delta_echem + delta_uth,
            "R20_erg": residual, "N_R20_erg": normalizer, "normalized": abs(residual) / normalizer}


def method_sequence(root: Path, method: str) -> list[dict[str, str]]:
    return [read_result(root, method, index) for index in range(241)]


def deterministic(candidate_root: Path, method: str, repeat: str) -> dict:
    unequal_results = []
    unequal_steps = []
    endpoint_rhs_failures = []
    endpoint_rhs_values: dict[tuple[int, str], tuple[bytes, bytes, bytes]] = {}
    for index in range(241):
        if result_path(candidate_root, method, index).read_bytes() != result_path(candidate_root, repeat, index).read_bytes():
            unequal_results.append(index)
        a, b = read_meta(candidate_root, method, index), read_meta(candidate_root, repeat, index)
        if any(a[name] != b[name] for name in ("accepted", "rejected", "rhs")):
            unequal_steps.append(index)
        if method == "hermite" and index not in (0, 240):
            for side, endpoint_field in (("L", "left_endpoint_index"), ("R", "right_endpoint_index")):
                names = tuple(f"f{side}_{component}" for component in ("x", "eta_e", "eta_mu"))
                values_a = tuple(bits(a[name]) for name in names)
                values_b = tuple(bits(b[name]) for name in names)
                if values_a != values_b:
                    endpoint_rhs_failures.append({"observation": index, "side": side})
                key = (int(a[endpoint_field]), side)
                if key in endpoint_rhs_values and endpoint_rhs_values[key] != values_a:
                    endpoint_rhs_failures.append({"observation": index, "side": side, "reason": "same endpoint changed"})
                endpoint_rhs_values[key] = values_a
    return {"pass": not unequal_results and not unequal_steps and not endpoint_rhs_failures,
            "unequal_result_observations": unequal_results, "unequal_step_observations": unequal_steps,
            "endpoint_rhs_failures": endpoint_rhs_failures}


def candidate_observation(matrix_row: dict, candidate: dict[str, str], oracle1: dict[str, str], oracle2: dict[str, str]) -> dict:
    yc, y1, y2 = state(candidate), state(oracle1), state(oracle2)
    state_utilizations = {}
    passed = True
    for component, name in enumerate(STATE):
        _, _, floor, budget = state_scales(matrix_row, y1, y2, component)
        d_o = abs(y2[component] - y1[component])
        d_o2 = ATOL_ORACLE2[component] + 1e-13 * max(abs(y1[component]), abs(y2[component]))
        f_o = max(d_o2, 64 * math.ulp(max(abs(y1[component]), abs(y2[component]))))
        u_o = 0.0 if matrix_row["exact"] else 2 * max(d_o, f_o)
        utilization = (abs(yc[component] - y2[component]) + u_o) / budget
        state_utilizations[name] = utilization
        passed &= utilization <= 1.0
    ledger_utilizations = {}
    for name in LEDGER:
        uncertainty, _, budget = oracle_diagnostic_uncertainty(oracle1, oracle2, name)
        utilization = (abs(finite(candidate[name], name) - finite(oracle2[name], name)) + uncertainty) / budget
        ledger_utilizations[name] = utilization
        passed &= utilization <= 1.0
    identity = validate_identity(candidate, oracle2)
    passed &= identity
    pscale = max(1.0, abs(finite(candidate["P_dir_actual_erg_s"], "Pdir")), abs(finite(candidate["P_dir_eq_erg_s"], "Pdir eq")))
    residual_utilizations = {name: abs(finite(candidate[name], name)) / (64 * EPS * pscale) for name in
                             ("R18_residual_erg_s", "Ra_Rb_residual_erg_s", "Rb_Rc_residual_erg_s")}
    residual_pass = all(value <= 1.0 for value in residual_utilizations.values())
    passed &= residual_pass
    return {"pass": passed, "state": state_utilizations, "ledger": ledger_utilizations,
            "identity": identity, "residuals": residual_utilizations, "residual_pass": residual_pass}


def evaluate_method(matrix: list[dict], oracle_root: Path, candidate_root: Path, method: str,
                    deterministic_result: dict) -> dict:
    observation_results = {}
    failures = []
    state_values, ledger_values = [], []
    for matrix_row in matrix:
        index = matrix_row["observation"]
        result = candidate_observation(matrix_row, read_result(candidate_root, method, index),
                                       read_result(oracle_root, "oracle1", index), read_result(oracle_root, "oracle2", index))
        observation_results[index] = result
        state_values.extend(result["state"].values())
        ledger_values.extend(result["ledger"].values())
        if not result["pass"]:
            failures.append(index)
    candidate_sequence = method_sequence(candidate_root, method)
    oracle_sequence = method_sequence(oracle_root, "oracle2")
    candidate_r20, oracle_r20 = r20(candidate_sequence), r20(oracle_sequence)
    scale = oracle_r20["N_R20_erg"]
    r20_metrics = {
        "luminosity_difference_utilization": abs(candidate_r20["outgoing_erg"] - oracle_r20["outgoing_erg"]) / (5e-6 * scale),
        "thermal_difference_utilization": abs(candidate_r20["DeltaUth_erg"] - oracle_r20["DeltaUth_erg"]) / (5e-6 * scale),
        "endpoint_propagation_utilization": abs(candidate_r20["endpoint_sum_erg"] - oracle_r20["endpoint_sum_erg"]) / (5e-6 * scale),
        "R20_normalized_utilization": candidate_r20["normalized"] / 2e-4,
    }
    r20_pass = all(value <= 1.0 for value in r20_metrics.values())
    witness_pass = True
    witness_max = 0.0
    if method == "replay2":
        for matrix_row in matrix:
            index = matrix_row["observation"]
            y1 = state(read_result(candidate_root, "replay1", index))
            y2 = state(read_result(candidate_root, "replay2", index))
            o1 = state(read_result(oracle_root, "oracle1", index))
            o2 = state(read_result(oracle_root, "oracle2", index))
            for component in range(3):
                _, _, floor, _ = state_scales(matrix_row, o1, o2, component)
                utilization = abs(y2[component] - y1[component]) / floor
                witness_max = max(witness_max, utilization)
                witness_pass &= utilization <= 1.0
    worst = max(state_values + ledger_values + list(r20_metrics.values()) + ([witness_max] if method == "replay2" else [0.0]))
    passed = not failures and r20_pass and deterministic_result["pass"] and witness_pass
    return {
        "pass": passed, "failure_count": len(failures), "failure_observations": failures,
        "non_knot_failures": [index for index in failures if matrix[index - 1]["category"] == "A"],
        "knot_failures": [index for index in failures if matrix[index - 1]["category"] in ("B", "C")],
        "maximum_state_utilization": max(state_values), "median_state_utilization": statistics.median(state_values),
        "p90_state_utilization": percentile(state_values, 0.90), "p95_state_utilization": percentile(state_values, 0.95),
        "p99_state_utilization": percentile(state_values, 0.99), "maximum_ledger_utilization": max(ledger_values),
        "R20": {"candidate": candidate_r20, "oracle2": oracle_r20, "metrics": r20_metrics, "pass": r20_pass},
        "determinism": deterministic_result, "replay_witness_pass": witness_pass,
        "replay_witness_maximum_utilization": witness_max, "worst_utilization": worst,
        "observations": observation_results,
    }


def performance(root: Path, method: str) -> dict:
    data = [read_meta(root, method, index) for index in range(1, 241)]
    integrated = [row for row in data if int(row["accepted"]) > 0]
    wall = [finite(row["solve_wall_s"], "wall") for row in integrated]
    user = [finite(row["solve_cpu_user_s"], "user") for row in integrated]
    system = [finite(row["solve_cpu_sys_s"], "sys") for row in integrated]
    return {"local_solves": len(integrated), "wall_sum_s": sum(wall), "cpu_user_sum_s": sum(user),
            "cpu_sys_sum_s": sum(system), "median_wall_s": statistics.median(wall),
            "p95_wall_s": percentile(wall, 0.95), "maximum_wall_s": max(wall),
            "throughput_solves_per_wall_sum_s": len(wall) / sum(wall),
            "estimated_8192_wall_sum_s": sum(wall) * (8192 / 240)}


def immutable_hashes(args) -> dict:
    return {"matrix": sha256(args.matrix), "schedule": sha256(args.schedule), "accepted": sha256(args.accepted),
            "brackets": sha256(args.brackets), "internal": sha256(args.internal), "trajectory": sha256(args.trajectory),
            "steps": sha256(args.steps), "library": sha256(args.production_library)}


def final(args) -> None:
    matrix = read_matrix(args.matrix)
    methods = {"linear": "linear-repeat", "hermite": "hermite-repeat", "replay2": "replay2-repeat"}
    deterministic_results = {name: deterministic(args.candidate_root, name, repeat) for name, repeat in methods.items()}
    replay1_determinism = deterministic(args.candidate_root, "replay1", "replay1-repeat")
    method_results = {name: evaluate_method(matrix, args.oracle_root, args.candidate_root, name, deterministic_results[name])
                      for name in methods}
    method_results["replay2"]["replay1_determinism"] = replay1_determinism
    method_results["replay2"]["pass"] &= replay1_determinism["pass"]
    hermite_failures = method_results["hermite"]["failure_observations"]
    hybrid_eligible = (bool(hermite_failures) and not method_results["hermite"]["non_knot_failures"] and
                       all(method_results["replay2"]["observations"][index]["pass"] for index in hermite_failures))
    hybrid_result = None
    if hybrid_eligible:
        hybrid_rows = []
        for index in range(241):
            if index == 0:
                selected = "hermite"
            else:
                selected = "replay2" if matrix[index - 1]["category"] in ("B", "C") else "hermite"
            hybrid_rows.append(read_result(args.candidate_root, selected, index))
        oracle_rows = method_sequence(args.oracle_root, "oracle2")
        h_r20, o_r20 = r20(hybrid_rows), r20(oracle_rows)
        scale = o_r20["N_R20_erg"]
        metrics = {
            "luminosity_difference_utilization": abs(h_r20["outgoing_erg"] - o_r20["outgoing_erg"]) / (5e-6 * scale),
            "thermal_difference_utilization": abs(h_r20["DeltaUth_erg"] - o_r20["DeltaUth_erg"]) / (5e-6 * scale),
            "endpoint_propagation_utilization": abs(h_r20["endpoint_sum_erg"] - o_r20["endpoint_sum_erg"]) / (5e-6 * scale),
            "R20_normalized_utilization": h_r20["normalized"] / 2e-4,
        }
        observation_worst = max(
            max((method_results["replay2"] if row["category"] in ("B", "C") else method_results["hermite"])["observations"][row["observation"]][kind].values())
            for row in matrix for kind in ("state", "ledger"))
        hybrid_result = {"pass": all(value <= 1 for value in metrics.values()), "R20": metrics,
                         "worst_utilization": max(observation_worst, *metrics.values())}
        hybrid_eligible &= hybrid_result["pass"]
    passing_singles = [name for name, value in method_results.items() if value["pass"]]
    comfortable = [name for name in ("linear", "hermite", "replay2") if name in passing_singles and method_results[name]["worst_utilization"] <= 0.25]
    if comfortable:
        selected = comfortable[0]
        reason = "factor-four comfortably passing single; fixed cost order linear, Hermite, Replay-2"
    elif passing_singles:
        selected = min(passing_singles, key=lambda name: (method_results[name]["worst_utilization"], ("linear", "hermite", "replay2").index(name)))
        reason = "passing single with smallest worst budget utilization"
    elif hybrid_eligible:
        selected = "hybrid"
        reason = "no single passed; predeclared knot-only Hermite/replay hybrid eligible"
    else:
        selected = None
        reason = "no method satisfied the complete predeclared budget"
    ledgers = []
    for root, names in ((args.oracle_root, ("oracle1", "oracle2")),
                        (args.candidate_root, ("replay1", "replay1-repeat", "replay2", "replay2-repeat"))):
        for name in names:
            ledgers.extend(ledger_rows(root / name / "solve_ledger.tsv"))
    solve_ids = [row["solve_id"] for row in ledgers]
    hashes = immutable_hashes(args)
    immutable_pass = hashes["matrix"] == MATRIX_SHA and all(hashes[name] == PASSIVE_HASHES[name] for name in PASSIVE_HASHES)
    result = {
        "classification": "PHASE-6 CHECKPOINT-RECONSTRUCTION NUMERICAL VALIDATION EVIDENCE",
        "authorized_local_integrations": 1434, "executed_local_integrations": len(ledgers),
        "unique_solve_ids": len(set(solve_ids)), "solve_count_pass": len(ledgers) == len(set(solve_ids)) == 1434,
        "main_trajectory_integrations": 0, "methods": method_results,
        "hybrid_eligible": hybrid_eligible, "hybrid": hybrid_result,
        "selected_method": selected, "selection_reason": reason,
        "performance": {name: performance(root, name) for root, names in
                        ((args.oracle_root, ("oracle1", "oracle2")),
                         (args.candidate_root, ("replay1", "replay1-repeat", "replay2", "replay2-repeat"))) for name in names},
        "process_concurrency": 2, "immutable_hashes": hashes, "immutable_inputs_pass": immutable_pass,
        "historical_BA12": "FAIL", "historical_BA12R": "FAIL",
        "production_candidate_created": False, "production_ADR_created": False,
    }
    result["pass"] = bool(selected) and result["solve_count_pass"] and immutable_pass
    write_json(args.output, result)
    if not result["pass"]:
        raise SystemExit(2)


parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest="command", required=True)

p = sub.add_parser("preflight")
p.add_argument("--repo", type=Path, required=True)
p.add_argument("--matrix", type=Path, required=True)
for name in ("schedule", "accepted", "brackets", "internal", "trajectory", "steps"):
    p.add_argument(f"--{name.replace('_', '-')}", dest=name, type=Path, required=True)
p.add_argument("--production-library", type=Path, required=True)
p.add_argument("--predeclaration", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=preflight)

p = sub.add_parser("oracle")
p.add_argument("--matrix", type=Path, required=True)
p.add_argument("--root", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.add_argument("--flag", type=Path, required=True)
p.set_defaults(function=oracle)

p = sub.add_parser("final")
p.add_argument("--matrix", type=Path, required=True)
p.add_argument("--oracle-root", type=Path, required=True)
p.add_argument("--candidate-root", type=Path, required=True)
for name in ("schedule", "accepted", "brackets", "internal", "trajectory", "steps"):
    p.add_argument(f"--{name.replace('_', '-')}", dest=name, type=Path, required=True)
p.add_argument("--production-library", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=final)

arguments = parser.parse_args()
arguments.function(arguments)
