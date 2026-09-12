#!/usr/bin/env python3
"""Closed producer-authoritative Phase-5D promotion-candidate schema."""

import hashlib
import json
import re
from pathlib import Path

import numpy as np

SCHEMA_NAME = "compactstar.phase5d.controlled-evolution-promotion"
SCHEMA_VERSION = 1
SCIENCE_IMPLEMENTATION_SHA = "d3670f6d4e021def0483909b6d2fdeed1c6973a4"
CANONICAL_SHA = "d019ae390be4f5e3daba05039903485cb497e397"
FROZEN_CONTEXT_SHA = "f7116c1408c06f976527f86d4397ad6d4540dedf"
PLAN_SHA = "b53f26afeaa29e2e43744e27edda2cb5a94fa7c6"
ENVELOPE = "FR2005 Eq. (49) / PCY97 fully accreted"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _supports(console):
    mapping = {"2": "Me", "3": "Mmu"}
    result = {}
    for line in console.splitlines():
        if not line.startswith("FROZEN_LTILDE "):
            continue
        words = line.split()
        process = mapping[words[1]]
        positions = [i for i, value in enumerate(words) if value == "support_km"]
        result[process] = [
            [float(words[i + 1]), float(words[i + 2])] for i in positions
        ]
    if set(result) != {"Me", "Mmu"}:
        raise RuntimeError("missing Ltilde support metadata")
    return result


def _oracle_metrics(console):
    def value(label):
        match = re.search(label + r" PASS max_relative_error (\S+)", console)
        if not match:
            raise RuntimeError("missing oracle metric: " + label)
        return float(match.group(1))
    mutations = sorted(set(re.findall(r"^MUTATION (\S+) DETECTED", console, re.M)))
    if len(mutations) != 22:
        raise RuntimeError("expected 22 named current coupled mutation controls")
    return {
        "spin_only_max_relative_error": value("SPIN_ONLY"),
        "reaction_only_max_relative_error": value("REACTION_ONLY"),
        "independent_mutation_families": 16,
        "named_mutation_count": len(mutations),
        "lyapunov_pass": "ACTIVE_LYAPUNOV PASS" in console,
        "same_ltilde_RE9_pass": "COUPLED_RE9 PASS exact same-Ltilde equilibrium" in console,
    }


def _without(value, excluded):
    if isinstance(value, dict):
        return {key: _without(item, excluded) for key, item in value.items()
                if key not in excluded}
    if isinstance(value, list):
        return [_without(item, excluded) for item in value]
    return value


def build_artifact(raw, trajectory_root, qualification_root, entry_manifest,
                   oracle_root, suite_results):
    trajectory_root = Path(trajectory_root)
    qualification_root = Path(qualification_root)
    entry_manifest = Path(entry_manifest)
    oracle_root = Path(oracle_root)
    manifest = json.loads(entry_manifest.read_text())
    qualification = json.loads(
        (qualification_root / "qualification-evidence.json").read_text()
    )
    console = (trajectory_root / "console.log").read_text()
    oracle_console = (oracle_root / "console.log").read_text()
    profile_paths = list(
        qualification_root.glob("chemical-characterization/run-*/t8192-r80000")
    )
    if len(profile_paths) != 1:
        raise RuntimeError("qualification profile topology changed")
    profile = profile_paths[0]
    profile_data = np.loadtxt(profile / "profile.tsv", skiprows=1)
    baseline = np.genfromtxt(trajectory_root / "trajectory.tsv", names=True)
    full_power = baseline["LH"] - baseline["Lnu_full"]
    # A structured-array copy cannot acquire a field; construct the crossing
    # input explicitly while retaining the exact checkpoint coordinates.
    aggregate_full = []
    year = 365.25 * 86400
    for index in np.flatnonzero(full_power[:-1] * full_power[1:] < 0):
        fraction = float(-full_power[index] / (full_power[index + 1] - full_power[index]))
        aggregate_full.append({
            "time_bracket_yr": [float(baseline["t_s"][index] / year),
                                float(baseline["t_s"][index + 1] / year)],
            "power_bracket_erg_s": [float(full_power[index]),
                                    float(full_power[index + 1])],
            "estimated_time_yr_log_interpolation": float(
                np.exp((1 - fraction) * np.log(baseline["t_s"][index] / year)
                       + fraction * np.log(baseline["t_s"][index + 1] / year))
            ),
        })
    thermal_hashes = {
        Path(path).name: value
        for path, value in raw["run_provenance"]["inputs"].items()
        if "/thermal/eos." in path
    }
    if set(thermal_hashes) != {"eos.t", "eos.nb", "eos.yq", "eos.thermo"}:
        raise RuntimeError("thermal source inventory incomplete")
    production_hashes = {
        path: value for path, value in raw["source_hashes"].items()
        if path.startswith(("CompactStar/Physics/", "CompactStar/Analysis/"))
    }
    if not production_hashes:
        raise RuntimeError("scientific production source inventory missing")
    artifact = {
        "schema": {
            "name": SCHEMA_NAME,
            "version": SCHEMA_VERSION,
            "producer_contract": "fresh-context-v1",
        },
        "classification": {
            "classification": "promotion_candidate",
            "candidate_only": True,
            "governed_baseline": False,
            "benchmark_scope": "CONTROLLED_MATHEMATICAL_ARCHITECTURE",
            "physical_spin_interpretation": False,
            "super_kepler": True,
        },
        "scientific_contract": {
            "normalizations": raw["normalizations"],
            "enabled_processes": ["Me", "Mmu"],
            "disabled_processes": ["De", "Dmu"],
            "frozen_background": True,
            "dot_Z": False,
            "thermal_ledger": "Pnet=LH-DeltaLnu-Lnu_eq-Lgamma-Lother_neutrino",
            "unchanged_caveat_ids": [
                "LOTHER_PBF_ONLY_INTERFACE", "RKF45_STABILITY_BOUND",
                "HEAT_CAPACITY_CACHE_KINKS", "SOMMERFELD_FREE_GAS",
                "NO_DOT_Z", "NO_DU", "NO_SUPERFLUIDITY",
            ],
        },
        "fixture": {
            **raw["fixture"],
            "initial": raw["initial"],
            "state_layout": {
                "y0": "ln(Tinf/1e8 K)", "y1": "eta_npe_inf [MeV]",
                "y2": "eta_npmu_inf [MeV]",
                "typed_channels": ["Npe", "NpMu"],
                "public_deltaN_state": False,
            },
            "mass_geometric_km": float(profile_data[-1, 1]),
            "radius_km": float(profile_data[-1, 0]),
        },
        "source_provenance": {
            "canonical_sha": CANONICAL_SHA,
            "frozen_context_plan_sha": FROZEN_CONTEXT_SHA,
            "artifact_preparation_plan_sha": PLAN_SHA,
            "science_implementation_sha": SCIENCE_IMPLEMENTATION_SHA,
            "scientific_production_source_hashes": production_hashes,
            "qualification_goal_source": qualification["goals_source"],
        },
        "coefficient_authority": {
            "Ltilde": raw["Ltilde"],
            "W_I": raw["semantic_W_I"],
            "Z": raw["semantic_Z"],
            "support_km": _supports(console),
            "qualification_values": qualification["outputs"]["t8192-r80000"]["values"],
        },
        "thermal_authority": {
            "source_hashes": thermal_hashes,
            "envelope_provenance": ENVELOPE,
            "formula_changed": False,
        },
        "spin_driver": {
            **raw["spin"],
            "P0_s": 0.001,
            "run_interval_yr": [0.0, 1.0e10],
            "physical_spin_interpretation": False,
        },
        "solver": {
            "configuration": raw["solver"],
            "step_statistics": raw["stiffness"]["step_statistics"],
        },
        "oracles": _oracle_metrics(oracle_console),
        "trajectory_checkpoints": raw["checkpoints"],
        "trajectory_summary": {
            "final": raw["final"], "ranges": raw["ranges"],
            "ledger_residuals": raw["ledger_residuals"],
            "power_convergence_absolute": raw["power_convergence_absolute"],
        },
        "convergence": {
            "ODE": raw["ODE_convergence"],
            "initial_conditions": raw["initial_condition_convergence"],
        },
        "quasi_steady": _without(raw["quasi_steady"], {"authority"}),
        "thermal_sign_crossings": {
            "incremental_roots": raw["incremental_root_crossings"],
            "full_roots": raw["full_root_crossings"],
            "aggregate_incremental": raw["aggregate_incremental_power_crossings"],
            "aggregate_full": aggregate_full,
            "method": "checkpoint_bracket_log_interpolation",
        },
        "numerical_diagnostics": {
            "radial_Ltilde": raw["physical_Ltilde_radial_comparison"],
            "jacobians": raw["stiffness"]["jacobians"],
            "physical_timescales": _without(
                raw["stiffness"]["physical_timescales"], {"definition"}
            ),
            "net_power_cancellation": raw["net_power_cancellation_diagnostic"],
        },
        "protected_upstream": {
            "entry_sha": manifest["entry_sha"],
            "protected": manifest["protected"],
            "baselines": manifest["baselines"],
            "special": manifest["special"],
            "all_equal_at_release": raw["protected_gate"]["all_equal_entry"],
        },
        "suite_results": suite_results,
        "producer_provenance": {
            "producer": "tests/rotochemical/fresh_context.py",
            "version": 1,
            "historical_candidate_read_during_generation": False,
            "promotion_candidate_read_during_generation": False,
            "phase5d_baseline_read_during_generation": False,
            "logical_qualification_hashes": qualification["logical_hashes"],
        },
    }
    return artifact


def write_artifact(artifact, path):
    Path(path).write_text(
        json.dumps(artifact, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
