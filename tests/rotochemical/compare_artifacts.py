#!/usr/bin/env python3
"""Fail-closed comparator for Phase-5D historical, promotion, and baseline states."""

import argparse
import copy
import json
from pathlib import Path

from artifact_schema import ENVELOPE, SCHEMA_NAME, SCHEMA_VERSION

TOP_LEVEL = {
    "schema", "classification", "scientific_contract", "fixture",
    "source_provenance", "coefficient_authority", "thermal_authority",
    "spin_driver", "solver", "oracles", "trajectory_checkpoints",
    "trajectory_summary", "convergence", "quasi_steady",
    "thermal_sign_crossings", "numerical_diagnostics", "protected_upstream",
    "suite_results", "producer_provenance",
}
BASELINE_TRANSITIONS = {
    "classification": ("promotion_candidate", "governed"),
    "candidate_only": (True, False),
    "governed_baseline": (False, True),
}
SECTION_KEYS = {
    "schema": {"name", "version", "producer_contract"},
    "classification": {"classification", "candidate_only", "governed_baseline",
                       "benchmark_scope", "physical_spin_interpretation", "super_kepler"},
    "scientific_contract": {"normalizations", "enabled_processes", "disabled_processes",
                            "frozen_background", "dot_Z", "thermal_ledger",
                            "unchanged_caveat_ids"},
    "fixture": {"radial_resolution", "EOS_resolution", "rho_c_g_cm3", "initial",
                "state_layout", "mass_geometric_km", "radius_km"},
    "source_provenance": {"canonical_sha", "frozen_context_plan_sha",
                          "artifact_preparation_plan_sha", "science_implementation_sha",
                          "scientific_production_source_hashes", "qualification_goal_source"},
    "coefficient_authority": {"Ltilde", "W_I", "Z", "support_km",
                              "qualification_values"},
    "thermal_authority": {"source_hashes", "envelope_provenance", "formula_changed"},
    "spin_driver": {"B_G", "P0_s", "PPdot", "run_interval_yr",
                    "physical_spin_interpretation"},
    "solver": {"configuration", "step_statistics"},
    "oracles": {"spin_only_max_relative_error", "reaction_only_max_relative_error",
                "independent_mutation_families", "named_mutation_count", "lyapunov_pass",
                "same_ltilde_RE9_pass"},
    "trajectory_summary": {"final", "ranges", "ledger_residuals",
                           "power_convergence_absolute"},
    "convergence": {"ODE", "initial_conditions"},
    "thermal_sign_crossings": {"incremental_roots", "full_roots",
                               "aggregate_incremental", "aggregate_full", "method"},
    "numerical_diagnostics": {"radial_Ltilde", "jacobians", "physical_timescales",
                              "net_power_cancellation"},
    "protected_upstream": {"entry_sha", "protected", "baselines", "special",
                           "all_equal_at_release"},
    "producer_provenance": {"producer", "version",
                            "historical_candidate_read_during_generation",
                            "promotion_candidate_read_during_generation",
                            "phase5d_baseline_read_during_generation",
                            "logical_qualification_hashes"},
}


def strip_keys(value, keys):
    if isinstance(value, dict):
        return {
            key: strip_keys(item, keys)
            for key, item in value.items() if key not in keys
        }
    if isinstance(value, list):
        return [strip_keys(item, keys) for item in value]
    return value


def validate_artifact(artifact, expected_state):
    if set(artifact) != TOP_LEVEL:
        raise RuntimeError("promotion top-level schema differs")
    if artifact["schema"] != {
        "name": SCHEMA_NAME, "version": SCHEMA_VERSION,
        "producer_contract": "fresh-context-v1",
    }:
        raise RuntimeError("promotion schema identity differs")
    for section, keys in SECTION_KEYS.items():
        if not isinstance(artifact[section], dict) or set(artifact[section]) != keys:
            raise RuntimeError("closed schema differs in section: " + section)
    expected_suites = {
        "response", "component_tolerances", "phase5b_regression",
        "phase5c_regression", "coupled_oracles", "trajectory", "validator",
    }
    if set(artifact["suite_results"]) != expected_suites:
        raise RuntimeError("suite result inventory differs")
    if expected_state not in {"promotion_candidate", "governed"}:
        raise RuntimeError("unknown expected artifact state: " + expected_state)
    expected_classification = {
        "classification": expected_state,
        "candidate_only": expected_state == "promotion_candidate",
        "governed_baseline": expected_state == "governed",
        "benchmark_scope": "CONTROLLED_MATHEMATICAL_ARCHITECTURE",
        "physical_spin_interpretation": False, "super_kepler": True,
    }
    if artifact["classification"] != expected_classification:
        raise RuntimeError(expected_state + " classification differs")
    if artifact["thermal_authority"]["envelope_provenance"] != ENVELOPE:
        raise RuntimeError("future envelope provenance not corrected")
    if artifact["thermal_authority"]["formula_changed"] is not False:
        raise RuntimeError("envelope correction claims a formula change")
    flags = artifact["producer_provenance"]
    for name in [
        "historical_candidate_read_during_generation",
        "promotion_candidate_read_during_generation",
        "phase5d_baseline_read_during_generation",
    ]:
        if flags[name] is not False:
            raise RuntimeError("producer input boundary violated: " + name)
    for name, receipt in artifact["suite_results"].items():
        if receipt.get("raw_rc") != 0:
            raise RuntimeError("nonzero suite result: " + name)
        if receipt.get("failures", 0) != 0 or receipt.get("unexplained_skips", 0) != 0:
            raise RuntimeError("suite failures/skips: " + name)


def validate_promotion(artifact):
    validate_artifact(artifact, "promotion_candidate")


def validate_governed(artifact):
    validate_artifact(artifact, "governed")


def difference_inventory(left, right):
    differences = []

    def visit(a, b, path=()):
        if type(a) is not type(b):
            differences.append({"path": list(path), "left": a, "right": b})
        elif isinstance(a, dict):
            for key in sorted(set(a) | set(b)):
                if key not in a:
                    differences.append({
                        "path": list(path + (key,)), "left_missing": True,
                        "right": b[key],
                    })
                elif key not in b:
                    differences.append({
                        "path": list(path + (key,)), "left": a[key],
                        "right_missing": True,
                    })
                else:
                    visit(a[key], b[key], path + (key,))
        elif isinstance(a, list):
            if len(a) != len(b):
                differences.append({
                    "path": list(path + ("<length>",)),
                    "left": len(a), "right": len(b),
                })
            for index, (a_value, b_value) in enumerate(zip(a, b)):
                visit(a_value, b_value, path + (index,))
        elif a != b:
            differences.append({"path": list(path), "left": a, "right": b})

    visit(left, right)
    return differences


def _historical_payload(historical):
    return {
        "fixture": {
            "fixture": historical["fixture"],
            "normalizations": historical["normalizations"],
            "initial": historical["initial"],
            "spin": historical["spin"],
            "state_layout": historical["state_layout"],
        },
        "source": {
            "canonical": historical["canonical"],
            "frozen_context_plan_sha": historical["frozen_context_plan_sha"],
            "science_implementation_sha": historical["implementation_sha"],
            "production_hashes": {
                key: value for key, value in historical["source_hashes"].items()
                if key.startswith(("CompactStar/Physics/", "CompactStar/Analysis/"))
            },
        },
        "coefficients": {
            "Ltilde": historical["Ltilde"],
            "W_I": historical["semantic_W_I"],
            "Z": historical["semantic_Z"],
            "support_km": historical["Ltilde_support_km"],
        },
        "thermal_hashes": sorted(historical["thermal_source_hashes"].values()),
        "solver": historical["solver"],
        "step_statistics": historical["stiffness"]["step_statistics"],
        "oracles": {
            key: historical["oracle_results"][key]
            for key in ["spin_only_max_relative_error",
                        "reaction_only_max_relative_error",
                        "independent_mutation_families"]
        },
        "checkpoints": historical["checkpoints"],
        "trajectory_summary": {
            "final": historical["final"], "ranges": historical["ranges"],
            "ledger_residuals": historical["ledger_residuals"],
            "power_convergence_absolute": historical["power_convergence_absolute"],
        },
        "convergence": {
            "ODE": historical["ODE_convergence"],
            "initial_conditions": historical["initial_condition_convergence"],
        },
        "quasi_steady": strip_keys(historical["quasi_steady"], {"authority"}),
        "crossings": {
            "incremental_roots": strip_keys(
                historical["incremental_root_crossings"], {"meaning"}
            ),
            "full_roots": strip_keys(
                historical["full_root_crossings"], {"meaning"}
            ),
            "aggregate_incremental": strip_keys(
                historical["aggregate_incremental_power_crossings"], {"meaning"}
            ),
            "aggregate_full": strip_keys(
                historical["aggregate_full_power_crossings"], {"meaning"}
            ),
        },
        "diagnostics": {
            "radial_Ltilde": historical["physical_Ltilde_radial_comparison"],
            "jacobians": historical["stiffness"]["jacobians"],
            "physical_timescales": strip_keys(
                historical["stiffness"]["physical_timescales"], {"definition"}
            ),
            "net_power_cancellation": strip_keys(
                historical["net_power_cancellation_diagnostic"], {"interpretation"}
            ),
        },
        "protected": {
            "entry_sha": historical["protected_entry"]["entry_sha"],
            "protected": sorted([
                {key: item[key] for key in ["path", "expected", "actual"]}
                for item in historical["protected_entry"]["protected"]
            ], key=lambda item: item["path"]),
            "baselines": sorted([
                {key: item[key] for key in ["path", "sha256"]}
                for item in historical["protected_entry"]["baselines"]
            ], key=lambda item: item["path"]),
            "special": sorted([
                {key: item[key] for key in ["path", "expected", "actual"]}
                for item in historical["protected_entry"]["special"]
            ], key=lambda item: item["path"]),
        },
        "candidate_state": {
            "candidate_only": historical["candidate_only"],
            "governed_baseline": historical["governed_baseline"],
            "protected_baseline_installed": historical["protected_baseline_installed"],
        },
        "suite_science": {
            "raw_rc": historical["complete_validation"]["raw_rc"],
            "failures": historical["complete_validation"]["failures"],
            "skips": historical["complete_validation"]["skips"],
            "trajectory_rc": historical["all_six_trajectories_raw_rc"],
            "validator_rc": historical["trajectory_validator_raw_rc"],
        },
    }


def _promotion_payload(promotion):
    return {
        "fixture": {
            "fixture": {key: promotion["fixture"][key] for key in
                        ["radial_resolution", "EOS_resolution", "rho_c_g_cm3"]},
            "normalizations": promotion["scientific_contract"]["normalizations"],
            "initial": promotion["fixture"]["initial"],
            "spin": {key: promotion["spin_driver"][key]
                     for key in ["B_G", "P0_s", "PPdot"]},
            "state_layout": promotion["fixture"]["state_layout"],
        },
        "source": {
            "canonical": promotion["source_provenance"]["canonical_sha"],
            "frozen_context_plan_sha": promotion["source_provenance"]["frozen_context_plan_sha"],
            "science_implementation_sha": promotion["source_provenance"]["science_implementation_sha"],
            "production_hashes": promotion["source_provenance"]["scientific_production_source_hashes"],
        },
        "coefficients": {
            key: promotion["coefficient_authority"][key]
            for key in ["Ltilde", "W_I", "Z", "support_km"]
        },
        "thermal_hashes": sorted(promotion["thermal_authority"]["source_hashes"].values()),
        "solver": promotion["solver"]["configuration"],
        "step_statistics": promotion["solver"]["step_statistics"],
        "oracles": {
            key: promotion["oracles"][key]
            for key in ["spin_only_max_relative_error",
                        "reaction_only_max_relative_error",
                        "independent_mutation_families"]
        },
        "checkpoints": promotion["trajectory_checkpoints"],
        "trajectory_summary": promotion["trajectory_summary"],
        "convergence": promotion["convergence"],
        "quasi_steady": promotion["quasi_steady"],
        "crossings": {
            key: strip_keys(promotion["thermal_sign_crossings"][key], {"meaning"})
            for key in ["incremental_roots", "full_roots",
                        "aggregate_incremental", "aggregate_full"]
        },
        "diagnostics": {
            "radial_Ltilde": promotion["numerical_diagnostics"]["radial_Ltilde"],
            "jacobians": promotion["numerical_diagnostics"]["jacobians"],
            "physical_timescales": strip_keys(
                promotion["numerical_diagnostics"]["physical_timescales"], {"definition"}
            ),
            "net_power_cancellation": strip_keys(
                promotion["numerical_diagnostics"]["net_power_cancellation"], {"interpretation"}
            ),
        },
        "protected": {
            "entry_sha": promotion["protected_upstream"]["entry_sha"],
            "protected": sorted(promotion["protected_upstream"]["protected"],
                                key=lambda item: item["path"]),
            "baselines": sorted(promotion["protected_upstream"]["baselines"],
                                key=lambda item: item["path"]),
            "special": sorted(promotion["protected_upstream"]["special"],
                              key=lambda item: item["path"]),
        },
        "candidate_state": {
            "candidate_only": promotion["classification"]["candidate_only"],
            "governed_baseline": promotion["classification"]["governed_baseline"],
            "protected_baseline_installed": False,
        },
        "suite_science": {
            "raw_rc": 0, "failures": 0, "skips": 0,
            "trajectory_rc": promotion["suite_results"]["trajectory"]["raw_rc"],
            "validator_rc": promotion["suite_results"]["validator"]["raw_rc"],
        },
    }


def compare_historical(historical, promotion):
    validate_promotion(promotion)
    if "iron Potekhin1997 envelope" not in historical["thermal_authorities"]["photon"]:
        raise RuntimeError("historical envelope label no longer matches ratified evidence")
    if promotion["thermal_authority"]["envelope_provenance"] != ENVELOPE:
        raise RuntimeError("sole metadata correction differs")
    old = _historical_payload(historical)
    new = _promotion_payload(promotion)
    if old != new:
        for key in old:
            if old[key] != new.get(key):
                raise RuntimeError("historical scientific payload differs at class: " + key)
        raise RuntimeError("historical scientific payload differs")


def compare_baseline(candidate, baseline):
    validate_promotion(candidate)
    validate_governed(baseline)
    candidate_copy = copy.deepcopy(candidate)
    baseline_copy = copy.deepcopy(baseline)
    for key, transition in BASELINE_TRANSITIONS.items():
        if (candidate_copy["classification"].pop(key),
                baseline_copy["classification"].pop(key)) != transition:
            raise RuntimeError("unauthorized baseline classification transition: " + key)
    if candidate_copy != baseline_copy:
        raise RuntimeError("candidate-to-baseline changed nonclassification content")
    return difference_inventory(candidate, baseline)


def compare_governed(left, right):
    validate_governed(left)
    validate_governed(right)
    differences = difference_inventory(left, right)
    if differences:
        raise RuntimeError(
            "governed artifact differs at: "
            + ".".join(str(value) for value in differences[0]["path"])
        )
    return differences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["historical-promotion", "determinism",
                                         "candidate-baseline", "governed-baseline"])
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--difference-output", type=Path)
    args = parser.parse_args()
    inventory = []
    if args.mode == "determinism":
        if args.left.read_bytes() != args.right.read_bytes():
            raise RuntimeError("producer-authoritative artifact bytes differ")
    else:
        left = json.loads(args.left.read_text())
        right = json.loads(args.right.read_text())
        if args.mode == "historical-promotion":
            compare_historical(left, right)
        elif args.mode == "governed-baseline":
            inventory = compare_governed(left, right)
        else:
            inventory = compare_baseline(left, right)
    if args.difference_output:
        output = args.difference_output.resolve()
        if output.exists():
            raise RuntimeError("difference inventory output must be fresh")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps({
            "mode": args.mode,
            "differences": inventory,
            "unexpected_differences": [],
        }, indent=2, sort_keys=True) + "\n")
    print("PASS", args.mode)


if __name__ == "__main__":
    main()
