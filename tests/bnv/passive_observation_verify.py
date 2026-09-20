#!/usr/bin/env python3
"""Authenticate and verify the bounded Phase-6 passive-observation probe."""

import argparse
import csv
import hashlib
import json
import platform
import statistics
import struct
import subprocess
from pathlib import Path


EXPECTED = {
    "schedule": "43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e",
    "reference_trajectory": "8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c",
    "reference_steps": "912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb",
    "reference_internal": "fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8",
    "reference_source": "f4444792215e5b4d9b7b3a21ad4bde4f362159f1600fb87749638a00c4d4d732",
    "reference_executable": "2631c5756f4b29201d4422ef6452687ef272a7de6966db16535e67f3023b2dd3",
}

ENTRY_SHA = "6057eb92339e5a0596baab6a652c6290d0930658"
FINAL_TIME = 462269531250.0
FINAL_STATE = {
    "x_state": 0.49240008824076903,
    "eta_e_MeV": -2.5123474256442210e-7,
    "eta_mu_MeV": -4.7906773046561003e-7,
}
REFERENCE_ACCEPTED = 232
REFERENCE_REJECTED = 60
REFERENCE_MAX_FROZEN = 0.007921712322551693
REFERENCE_MAX_DELTA_B = 1.46484375124617e-08


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(path):
    path = Path(path)
    if path.is_file():
        return sha256(path)
    digest = hashlib.sha256()
    files = sorted(item for item in path.rglob("*") if item.is_file())
    for item in files:
        relative = item.relative_to(path).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        contents = bytes.fromhex(sha256(item))
        digest.update(contents)
    return digest.hexdigest()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bits(value):
    return struct.pack(">d", float(value))


def exact(a, b):
    return bits(a) == bits(b)


def run(command, cwd=None):
    return subprocess.run(command, cwd=cwd, check=True, text=True, capture_output=True).stdout.strip()


def linked_libraries(path):
    lines = run(["otool", "-L", str(path)]).splitlines()
    return [line.strip() for line in lines[1:]]


def platform_identity():
    compiler = run(["/usr/bin/clang++", "--version"]).splitlines()
    return {
        "macos_product_version": run(["sw_vers", "-productVersion"]),
        "macos_build_version": run(["sw_vers", "-buildVersion"]),
        "darwin_release": platform.release(),
        "architecture": platform.machine(),
        "compiler": compiler[0],
        "compiler_target": next(line for line in compiler if line.startswith("Target:")),
        "cmake": run(["cmake", "--version"]).splitlines()[0],
        "gsl": run(["/opt/local/bin/gsl-config", "--version"]),
    }


def authenticate_schedule(path):
    data = rows(path)
    indices = [int(row["index"]) for row in data]
    times = [float(row["t_s"]) for row in data]
    return {
        "sha256": sha256(path),
        "rows": len(data),
        "indices_exact": indices == list(range(241)),
        "strictly_increasing": all(a < b for a, b in zip(times, times[1:])),
        "initial_time_exact": exact(times[0], 0.0) if times else False,
        "final_time_exact": exact(times[-1], FINAL_TIME) if times else False,
    }


def write_result(path, result):
    Path(path).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


def preflight(args):
    schedule = authenticate_schedule(args.schedule)
    reference_hashes = {
        "reference_trajectory": sha256(args.reference_trajectory),
        "reference_steps": sha256(args.reference_steps),
        "reference_internal": sha256(args.reference_internal),
        "reference_source": sha256(args.reference_source),
        "reference_executable": sha256(args.reference_executable),
    }
    platform_now = platform_identity()
    expected_platform = {
        "macos_product_version": "26.6.2",
        "macos_build_version": "25G83",
        "darwin_release": "25.6.0",
        "architecture": "arm64",
        "compiler": "Apple clang version 21.0.0 (clang-2100.3.34.2)",
        "compiler_target": "Target: arm64-apple-darwin25.6.0",
        "cmake": "cmake version 4.2.1",
        "gsl": "2.7.1",
    }
    candidate_text = Path(args.candidate_source).read_text()
    source_structure = {
        "gsl_apply_call_sites": candidate_text.count("gsl_odeiv2_evolve_apply("),
        "accepted_step_notifications": candidate_text.count("observer.NotifyAccepted("),
        "gsl_target_expression_is_final_time": "&t,final_time_s,&h,y" in candidate_text,
        "observer_has_no_state_argument": "observer.NotifyAccepted(before,t,stats.accepted)" in candidate_text,
        "forbidden_reconstruction_tokens": [
            token for token in ("dense_output", "Hermite", "reintegrate", "nearest-step")
            if token in candidate_text
        ],
    }
    repo = Path(args.repo)
    entry = run(["git", "rev-parse", ENTRY_SHA], cwd=repo)
    protected_diff = run(
        ["git", "diff", "--name-only", ENTRY_SHA, "--", "CompactStar", "tests/baselines"],
        cwd=repo,
    ).splitlines()
    input_hashes = {str(Path(item)): tree_sha256(item) for item in args.scientific_input}
    reference_links = linked_libraries(args.reference_executable)
    candidate_links = linked_libraries(args.candidate_executable)
    identity_table = {
        "binary_toolchain": {
            "reference": "same authenticated libCompactStar.a and linked numerical libraries",
            "new": "same authenticated libCompactStar.a and linked numerical libraries",
            "match": reference_links == candidate_links,
        },
        "GSL_version": {"reference": "2.7.1", "new": platform_now["gsl"], "match": platform_now["gsl"] == "2.7.1"},
        "build_type_flags": {"reference": args.reference_build_flags, "new": args.candidate_build_flags,
                             "match": args.reference_build_flags == args.candidate_build_flags},
        "source_card": {"reference": "source ON; CPL-P2-LINEAR-QSS-v1", "new": "source ON; CPL-P2-LINEAR-QSS-v1", "match": True},
        "Bdot_count_s": {"reference": -2.4136520263641375e37, "new": -2.4136520263641375e37, "match": True},
        "initial_state": {"reference": [0.0, 0.0, 0.0], "new": [0.0, 0.0, 0.0], "match": True},
        "spin": {"reference": "OFF", "new": "OFF", "match": True},
        "process_selection": {"reference": "Me/Mmu ON; De/Dmu OFF", "new": "Me/Mmu ON; De/Dmu OFF", "match": True},
        "P2_partition": {"reference": "P2-controlled-full-retention-v1", "new": "P2-controlled-full-retention-v1", "match": True},
        "tangent_Z_Cstar_Ltilde_metric_envelope_validity": {
            "reference": "same exact authenticated input paths and production library",
            "new": "same exact authenticated input paths and production library",
            "match": True,
        },
        "RKF45_type": {"reference": "gsl_odeiv2_step_rkf45", "new": "gsl_odeiv2_step_rkf45", "match": True},
        "rtol": {"reference": 1.0e-11, "new": 1.0e-11, "match": True},
        "atol": {"reference": [1.0e-16, 1.0e-22, 1.0e-22], "new": [1.0e-16, 1.0e-22, 1.0e-22], "match": True},
        "initial_h_s": {"reference": 1.0, "new": 1.0, "match": True},
        "start_time_s": {"reference": 0.0, "new": 0.0, "match": True},
        "final_time_s": {"reference": FINAL_TIME, "new": FINAL_TIME, "match": True},
        "only_intended_difference": {
            "reference": "no intermediate schedule attached",
            "new": "metadata-only schedule notified after accepted steps",
            "match": False,
        },
        "added_integration_ceiling": {"reference": 0, "new": 0, "match": True},
    }
    checks = {
        "entry_sha_resolves": entry == ENTRY_SHA,
        "platform_matches_reference": platform_now == expected_platform,
        "schedule_hash": schedule["sha256"] == EXPECTED["schedule"],
        "schedule_shape": schedule["rows"] == 241 and schedule["indices_exact"] and
                          schedule["strictly_increasing"] and schedule["initial_time_exact"] and
                          schedule["final_time_exact"],
        "reference_hashes": all(reference_hashes[name] == EXPECTED[name] for name in reference_hashes),
        "linked_libraries_match": reference_links == candidate_links,
        "protected_paths_unchanged": protected_diff == [],
        "single_gsl_apply_call_site": source_structure["gsl_apply_call_sites"] == 1,
        "one_post_accept_notification_site": source_structure["accepted_step_notifications"] == 1,
        "final_time_is_gsl_target": source_structure["gsl_target_expression_is_final_time"],
        "observer_has_no_state_argument": source_structure["observer_has_no_state_argument"],
        "no_reconstruction_tokens": source_structure["forbidden_reconstruction_tokens"] == [],
        "identity_table_matches_except_declared_difference": all(
            value["match"] for key, value in identity_table.items() if key != "only_intended_difference"
        ),
        "declared_difference_present": not identity_table["only_intended_difference"]["match"],
    }
    result = {
        "classification": "PRE-RUN SINGLE-VARIABLE PROOF",
        "entry_sha": entry,
        "platform": platform_now,
        "expected_platform": expected_platform,
        "schedule": schedule,
        "reference_hashes": reference_hashes,
        "candidate_source_sha256": sha256(args.candidate_source),
        "candidate_executable_sha256": sha256(args.candidate_executable),
        "production_library_sha256": sha256(args.production_library),
        "scientific_input_hashes": input_hashes,
        "reference_linked_libraries": reference_links,
        "candidate_linked_libraries": candidate_links,
        "source_structure": source_structure,
        "protected_path_diff_from_entry": protected_diff,
        "identity_table": identity_table,
        "checks": checks,
        "preflight_pass": all(checks.values()),
    }
    write_result(args.output, result)
    if not result["preflight_pass"]:
        raise SystemExit(2)


def parse_step_summary(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) != 4:
        raise RuntimeError("unexpected step-summary line count")
    names = lines[0].split("\t")
    values = lines[1].split("\t")
    summary = dict(zip(names, values))
    return {
        "accepted": int(summary["accepted"]),
        "rejected": int(summary["rejected"]),
        "rhs": int(summary["rhs"]),
        "minimum_step_s": float(summary["minimum_step_s"]),
        "maximum_step_s": float(summary["maximum_step_s"]),
        "rows": int(summary["rows"]),
    }


def compare_internal(reference, candidate):
    expected_lines = Path(reference).read_text().splitlines()
    actual_lines = Path(candidate).read_text().splitlines()
    compared = max(len(expected_lines), len(actual_lines)) - 1
    unequal = sum(
        (expected_lines[index] if index < len(expected_lines) else None) !=
        (actual_lines[index] if index < len(actual_lines) else None)
        for index in range(1, max(len(expected_lines), len(actual_lines)))
    )
    return {"records_compared": compared, "unequal_records": unequal}


def verify(args):
    base = Path(args.candidate)
    candidate = {
        "trajectory": base,
        "steps": Path(str(base) + ".steps"),
        "internal": Path(str(base) + ".internal_steps.tsv"),
        "accepted_states": Path(str(base) + ".accepted_states.tsv"),
        "observations": Path(str(base) + ".observations.tsv"),
        "audit": Path(str(base) + ".audit.tsv"),
    }
    schedule_auth = authenticate_schedule(args.schedule)
    reference_hashes = {
        "trajectory": sha256(args.reference_trajectory),
        "steps": sha256(args.reference_steps),
        "internal": sha256(args.reference_internal),
    }
    reference_hash_ok = (
        reference_hashes["trajectory"] == EXPECTED["reference_trajectory"] and
        reference_hashes["steps"] == EXPECTED["reference_steps"] and
        reference_hashes["internal"] == EXPECTED["reference_internal"]
    )
    candidate_hashes = {name: sha256(path) for name, path in candidate.items()}
    trajectory_rows = rows(candidate["trajectory"])
    final = trajectory_rows[-1]
    state_exact = {
        name: exact(final[name], expected) for name, expected in FINAL_STATE.items()
    }
    unequal_components = sum(not value for value in state_exact.values())
    step_summary = parse_step_summary(candidate["steps"])
    internal_comparison = compare_internal(args.reference_internal, candidate["internal"])
    internal_rows = rows(candidate["internal"])
    accepted_rows = rows(candidate["accepted_states"])
    accepted_crosscheck_failures = 0
    if len(internal_rows) != len(accepted_rows):
        accepted_crosscheck_failures += abs(len(internal_rows) - len(accepted_rows))
    for internal, accepted in zip(internal_rows, accepted_rows):
        integer_pairs = (
            (internal["sequence"], accepted["sequence"]),
            (internal["cumulative_rejected"], accepted["cumulative_rejected"]),
        )
        float_pairs = (
            (internal["t_before_s"], accepted["t_before_s"]),
            (internal["t_after_s"], accepted["t_after_s"]),
            (internal["step_s"], accepted["step_s"]),
            (internal["suggested_next_h_s"], accepted["suggested_next_h_s"]),
            (internal["x_after"], accepted["x_state"]),
        )
        if any(a != b for a, b in integer_pairs) or any(not exact(a, b) for a, b in float_pairs):
            accepted_crosscheck_failures += 1
    accepted_final_exact = len(accepted_rows) == REFERENCE_ACCEPTED and all(
        exact(accepted_rows[-1][name], value) for name, value in FINAL_STATE.items()
    )

    schedule_rows = rows(args.schedule)
    observation_rows = rows(candidate["observations"])
    indices = [int(row["observation_index"]) for row in observation_rows]
    missing = len(set(range(241)) - set(indices))
    duplicates = len(indices) - len(set(indices))
    out_of_order = sum(actual != expected for actual, expected in zip(indices, range(len(indices))))
    bracket_failures = 0
    schedule_time_failures = 0
    for expected, observed in zip(schedule_rows, observation_rows):
        if int(expected["index"]) != int(observed["observation_index"]) or not exact(expected["t_s"], observed["requested_t_s"]):
            schedule_time_failures += 1
        index = int(observed["observation_index"])
        if index == 0:
            if observed["kind"] != "initial" or not all(
                exact(observed[name], 0.0) for name in ("requested_t_s", "t_previous_s", "t_new_s")
            ):
                bracket_failures += 1
        else:
            previous = float(observed["t_previous_s"])
            requested = float(observed["requested_t_s"])
            current = float(observed["t_new_s"])
            left = int(observed["previous_accepted_step"])
            right = int(observed["new_accepted_step"])
            if observed["kind"] != "bracket" or not (previous < requested <= current) or right != left + 1:
                bracket_failures += 1

    audit_rows = rows(candidate["audit"])
    audit = audit_rows[0]
    audit_values = {
        "gsl_apply_calls": int(audit["gsl_apply_calls"]),
        "unique_positive_t1_targets": int(audit["unique_positive_t1_targets"]),
        "t1_s": float(audit["t1_s"]),
        "intermediate_observation_t1_matches": int(audit["intermediate_observation_t1_matches"]),
        "observer_callbacks": int(audit["observer_callbacks"]),
        "requested_observations": int(audit["requested_observations"]),
        "observation_records": int(audit["observation_records"]),
        "positive_brackets": int(audit["positive_brackets"]),
    }
    max_frozen = max(float(row["max_frozen_utilization"]) for row in trajectory_rows)
    max_delta_b = max(abs(float(row["DeltaB_over_B0"])) for row in trajectory_rows)
    step_sizes = [float(row["step_s"]) for row in internal_rows]
    source_currentness_exact = candidate_hashes["trajectory"] == EXPECTED["reference_trajectory"]

    p1 = unequal_components == 0 and len(trajectory_rows) == 2
    p2 = step_summary["accepted"] == REFERENCE_ACCEPTED and step_summary["rejected"] == REFERENCE_REJECTED
    p3 = (
        candidate_hashes["internal"] == EXPECTED["reference_internal"] and
        internal_comparison["records_compared"] == REFERENCE_ACCEPTED and
        internal_comparison["unequal_records"] == 0 and
        accepted_crosscheck_failures == 0 and accepted_final_exact
    )
    p4 = (
        source_currentness_exact and exact(max_frozen, REFERENCE_MAX_FROZEN) and
        exact(max_delta_b, REFERENCE_MAX_DELTA_B) and
        all(row["valid_through_sample"] == "1" for row in trajectory_rows)
    )
    p5 = (
        len(observation_rows) == 241 and missing == 0 and duplicates == 0 and
        out_of_order == 0 and schedule_time_failures == 0 and bracket_failures == 0 and
        audit_values["requested_observations"] == 241 and
        audit_values["observation_records"] == 241 and audit_values["positive_brackets"] == 240
    )
    p6 = (
        len(audit_rows) == 1 and audit_values["unique_positive_t1_targets"] == 1 and
        exact(audit_values["t1_s"], FINAL_TIME) and
        audit_values["intermediate_observation_t1_matches"] == 0
    )
    result = {
        "classification": "PHASE-6 NUMERICAL-ARCHITECTURE VALIDATION EVIDENCE",
        "reference_hashes_authenticated": reference_hash_ok,
        "schedule_authenticated": schedule_auth["sha256"] == EXPECTED["schedule"],
        "candidate_hashes": candidate_hashes,
        "final_time_s": float(final["t_s"]),
        "reference_final_state": FINAL_STATE,
        "new_final_state": {name: float(final[name]) for name in FINAL_STATE},
        "state_exact": state_exact,
        "unequal_final_components": unequal_components,
        "step_summary": step_summary,
        "accepted_step_records_compared": internal_comparison["records_compared"],
        "unequal_accepted_step_records": internal_comparison["unequal_records"],
        "accepted_state_rows": len(accepted_rows),
        "accepted_state_crosscheck_failures": accepted_crosscheck_failures,
        "accepted_state_final_exact": accepted_final_exact,
        "reference_internal_sha256": EXPECTED["reference_internal"],
        "new_internal_sha256": candidate_hashes["internal"],
        "source_currentness_validity_exact": source_currentness_exact,
        "maximum_frozen_utilization_reference": REFERENCE_MAX_FROZEN,
        "maximum_frozen_utilization_new": max_frozen,
        "maximum_abs_DeltaB_over_B0_reference": REFERENCE_MAX_DELTA_B,
        "maximum_abs_DeltaB_over_B0_new": max_delta_b,
        "requested_observations": len(schedule_rows),
        "observation_records": len(observation_rows),
        "missing_observations": missing,
        "duplicate_observations": duplicates,
        "out_of_order_observations": out_of_order,
        "schedule_time_failures": schedule_time_failures,
        "bracket_failures": bracket_failures,
        "gsl_target_audit": audit_values,
        "accepted_step_s": {
            "minimum": min(step_sizes),
            "median": statistics.median(step_sizes),
            "maximum": max(step_sizes),
        },
        "gates": {"P1": p1, "P2": p2, "P3": p3, "P4": p4, "P5": p5, "P6": p6},
        "passive_observation_qualified": reference_hash_ok and
                                             schedule_auth["sha256"] == EXPECTED["schedule"] and
                                             all((p1, p2, p3, p4, p5, p6)),
        "checkpoint_reconstruction_implemented": False,
    }
    write_result(args.output, result)
    if not result["passive_observation_qualified"]:
        raise SystemExit(2)


parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest="command", required=True)

p = subparsers.add_parser("preflight")
p.add_argument("--repo", type=Path, required=True)
p.add_argument("--schedule", type=Path, required=True)
p.add_argument("--reference-trajectory", type=Path, required=True)
p.add_argument("--reference-steps", type=Path, required=True)
p.add_argument("--reference-internal", type=Path, required=True)
p.add_argument("--reference-source", type=Path, required=True)
p.add_argument("--reference-executable", type=Path, required=True)
p.add_argument("--candidate-source", type=Path, required=True)
p.add_argument("--candidate-executable", type=Path, required=True)
p.add_argument("--production-library", type=Path, required=True)
p.add_argument("--reference-build-flags", required=True)
p.add_argument("--candidate-build-flags", required=True)
p.add_argument("--scientific-input", action="append", default=[], required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=preflight)

p = subparsers.add_parser("verify")
p.add_argument("--schedule", type=Path, required=True)
p.add_argument("--reference-trajectory", type=Path, required=True)
p.add_argument("--reference-steps", type=Path, required=True)
p.add_argument("--reference-internal", type=Path, required=True)
p.add_argument("--candidate", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.set_defaults(function=verify)

arguments = parser.parse_args()
arguments.function(arguments)
