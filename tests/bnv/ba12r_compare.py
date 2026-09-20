#!/usr/bin/env python3
"""Predeclared three-level BA12R comparison for the P2 source/control pair only."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

sys.dont_write_bytecode = True

EPS = sys.float_info.epsilon
MEV_TO_ERG = 1.602176634e-6
EXPECTED_ROWS = 8193
CARD = "CPL-P2-LINEAR-QSS-v1"
HISTORICAL_AGGREGATE = "81314ddceba77eaaf29f493928560e93f3a6b45d59f302c11fbec0e799da0945"
HISTORICAL_HASHES = {
    f"{CARD}.baseline.tsv": "f331b5bbfd6bdfb0f8e93095052cf1edc712219180370d23383bed1c18b31093",
    f"{CARD}.refined.tsv": "1889336755a87377336e1b688e34032dab428304861bd2a63e142e27555b7a3a",
    f"{CARD}.control.baseline.tsv": "be807f432c81b7f62594005b593bdd5bcfd301ddfd9ed2837178831fc2dbb9b8",
    f"{CARD}.control.refined.tsv": "be807f432c81b7f62594005b593bdd5bcfd301ddfd9ed2837178831fc2dbb9b8",
}
STATE = ("x_state", "eta_e_MeV", "eta_mu_MeV")
ATOL_REFINED = (1e-14, 1e-20, 1e-20)
ATOL_ULTRA = (1e-16, 1e-22, 1e-22)
LEDGER = (
    "P_dir_eq_erg_s", "P_dir_actual_erg_s", "LH_erg_s", "DeltaLnu_erg_s",
    "DeltaPbeta_erg_s", "Lnu_eq_erg_s", "Lnu_full_erg_s",
    "L_out_fluid_inf_erg_s", "Lgamma_erg_s", "Lother_erg_s", "Pnet_erg_s",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(value: str | float, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"nonfinite {label}")
    return result


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        result = list(csv.DictReader(stream, delimiter="\t"))
    if len(result) != EXPECTED_ROWS:
        raise RuntimeError(f"{path}: {len(result)} rows, expected {EXPECTED_ROWS}")
    return result


def ordered_manifest_digest(root: Path) -> str:
    payload = b""
    historical_label = "build/phase6a1-controlled-bnv-debug/phase6a1-candidate-raw-authoritative2"
    for path in sorted(item for item in root.iterdir() if item.is_file()):
        # This exactly reproduces: find "$root" -type f -print0 | sort -z |
        # xargs -0 shasum -a 256 | shasum -a 256
        payload += f"{digest(path)}  {historical_label}/{path.name}\n".encode()
    return hashlib.sha256(payload).hexdigest()


def gross_power(row: dict[str, str]) -> float:
    return max(1.0, sum(abs(finite(row[name], name)) for name in (
        "P_dir_actual_erg_s", "LH_erg_s", "DeltaLnu_erg_s", "Lnu_eq_erg_s",
        "Lgamma_erg_s", "Lother_erg_s", "L_out_fluid_inf_erg_s",
    )))


def trapezoid(sequence: list[dict[str, str]], value) -> float:
    total = 0.0
    for left, right in zip(sequence, sequence[1:]):
        dt = finite(right["t_s"], "time") - finite(left["t_s"], "time")
        if not dt > 0:
            raise RuntimeError("nonincreasing time grid")
        total += 0.5 * dt * (value(left) + value(right))
    return total


def cumulative_trapezoid(sequence: list[dict[str, str]], value) -> list[float]:
    result = [0.0]
    for left, right in zip(sequence, sequence[1:]):
        dt = finite(right["t_s"], "time") - finite(left["t_s"], "time")
        result.append(result[-1] + 0.5 * dt * (value(left) + value(right)))
    return result


def thermal_energy(sequence: list[dict[str, str]]) -> float:
    total = 0.0
    for left, right in zip(sequence, sequence[1:]):
        total += 0.5 * (finite(left["Cstar_erg_K"], "Cstar") +
                        finite(right["Cstar_erg_K"], "Cstar")) * (
            finite(right["Tinf_K"], "Tinf") - finite(left["Tinf_K"], "Tinf"))
    return total


def cumulative_thermal_energy(sequence: list[dict[str, str]]) -> list[float]:
    result = [0.0]
    for left, right in zip(sequence, sequence[1:]):
        result.append(result[-1] + 0.5 * (
            finite(left["Cstar_erg_K"], "Cstar") + finite(right["Cstar_erg_K"], "Cstar")) * (
            finite(right["Tinf_K"], "Tinf") - finite(left["Tinf_K"], "Tinf")))
    return result


def r20(sequence: list[dict[str, str]]) -> dict[str, float]:
    first, last = sequence[0], sequence[-1]
    delta_eeq = MEV_TO_ERG * finite(first["mu_B_inf_MeV"], "mu_B") * (
        finite(last["B_count"], "B final") - finite(first["B_count"], "B initial"))
    delta_echem = MEV_TO_ERG * (
        finite(last["Echem_MeV"], "Echem final") - finite(first["Echem_MeV"], "Echem initial"))
    delta_uth = thermal_energy(sequence)
    luminosity = lambda row: sum(finite(row[name], name) for name in (
        "L_out_fluid_inf_erg_s", "Lnu_full_erg_s", "Lgamma_erg_s", "Lother_erg_s"))
    outgoing = trapezoid(sequence, luminosity)
    residual = delta_eeq + delta_echem + delta_uth + outgoing
    normalizer_integral = trapezoid(sequence, lambda row: sum(abs(finite(row[name], name)) for name in (
        "Lnu_full_erg_s", "Lgamma_erg_s", "Lother_erg_s")))
    normalizer = max(1.0, abs(delta_uth), abs(delta_echem), normalizer_integral)
    even = sequence[::2]
    if even[-1] is not sequence[-1]:
        raise RuntimeError("even R20 quadrature grid lost endpoint")
    quadrature = abs(outgoing - trapezoid(even, luminosity))
    thermal_quadrature = abs(delta_uth - thermal_energy(even))
    return {
        "DeltaEeq_erg": delta_eeq,
        "DeltaEchem_erg": delta_echem,
        "DeltaUth_erg": delta_uth,
        "R20_erg": residual,
        "N_R20_erg": normalizer,
        "R20_over_N": abs(residual) / normalizer,
        "luminosity_quadrature_over_N": quadrature / normalizer,
        "thermal_quadrature_over_N": thermal_quadrature / normalizer,
    }


def percentile95(values: list[int]) -> int:
    ordered = sorted(values)
    return ordered[max(0, math.ceil(0.95 * len(ordered)) - 1)]


def step_trend(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if len(lines) != EXPECTED_ROWS + 2:
        raise RuntimeError(f"{path}: step evidence row count differs")
    summary = dict(zip(lines[0].split("\t"), (finite(x, path.name) for x in lines[1].split("\t"))))
    names = lines[2].split("\t")
    if names != ["t_s", "cumulative_accepted", "cumulative_rejected", "last_step_s",
                 "minimum_step_since_previous_output_s"]:
        raise RuntimeError(f"{path}: BA12R step schema differs")
    output = [dict(zip(names, (finite(x, path.name) for x in line.split("\t")))) for line in lines[3:]]
    previous_a = previous_r = 0
    intervals = []
    for item in output:
        accepted = int(item["cumulative_accepted"])
        rejected = int(item["cumulative_rejected"])
        intervals.append({**item, "accepted_increment": accepted - previous_a,
                          "rejected_increment": rejected - previous_r})
        if intervals[-1]["accepted_increment"] <= 0 or intervals[-1]["rejected_increment"] < 0:
            raise RuntimeError(f"{path}: invalid step-count increment")
        previous_a, previous_r = accepted, rejected
    final_time = intervals[-1]["t_s"]
    decade = [item for item in intervals if item["t_s"] >= final_time / 10.0]
    midpoint = 0.5 * (final_time / 10.0 + final_time)
    first = [item["accepted_increment"] for item in decade if item["t_s"] <= midpoint]
    last = [item["accepted_increment"] for item in decade if item["t_s"] > midpoint]
    accepted = sum(item["accepted_increment"] for item in decade)
    rejected = sum(item["rejected_increment"] for item in decade)
    rejection_fraction = rejected / (accepted + rejected)
    minimum = min(item["minimum_step_since_previous_output_s"] for item in decade)
    spacing = final_time / (EXPECTED_ROWS - 1)
    first_median, last_median = statistics.median(first), statistics.median(last)
    first_p95, last_p95 = percentile95(first), percentile95(last)
    passed = (rejection_fraction <= 0.10 and last_median <= 2 * first_median and
              last_p95 <= 4 * first_p95 and minimum >= 1e-6 * spacing and
              summary["accepted"] < 100000 * (EXPECTED_ROWS - 1))
    return {
        "accepted_steps": int(summary["accepted"]),
        "rejected_steps": int(summary["rejected"]),
        "rhs_evaluations": int(summary["rhs"]),
        "overall_minimum_step_s": summary["minimum_step_s"],
        "final_decade_rejection_fraction": rejection_fraction,
        "final_decade_minimum_accepted_step_s": minimum,
        "fixed_output_spacing_s": spacing,
        "first_half_median_accepted_steps_per_interval": first_median,
        "last_half_median_accepted_steps_per_interval": last_median,
        "first_half_p95_accepted_steps_per_interval": first_p95,
        "last_half_p95_accepted_steps_per_interval": last_p95,
        "pass": passed,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-root", type=Path, required=True)
    parser.add_argument("--ultra-source", type=Path, required=True)
    parser.add_argument("--ultra-control", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    historical_root = args.historical_root.resolve()
    if ordered_manifest_digest(historical_root) != HISTORICAL_AGGREGATE:
        raise RuntimeError("historical 32-file aggregate authentication failed")
    for name, expected in HISTORICAL_HASHES.items():
        if digest(historical_root / name) != expected:
            raise RuntimeError(f"historical hash mismatch: {name}")

    ultra_paths = {"source": args.ultra_source.resolve(), "control": args.ultra_control.resolve()}
    evidence: dict[str, object] = {
        "schema_id": "compactstar.phase6a1.ba12r-validation-evidence.v1",
        "classification": "BA12R VALIDATION EVIDENCE; NOT BNV CANDIDATE; NOT GOVERNED BASELINE",
        "historical_aggregate_sha256": HISTORICAL_AGGREGATE,
        "historical_hashes": HISTORICAL_HASHES,
        "ultra_hashes": {},
        "state": {}, "ledger": {}, "endpoint": {}, "identity_budgets": {},
        "solver_trend": {}, "matched_source_minus_control": {},
    }
    levels: dict[str, dict[str, list[dict[str, str]]]] = {"source": {}, "control": {}}
    for mode in ("source", "control"):
        stem = CARD + (".control" if mode == "control" else "")
        levels[mode]["baseline"] = rows(historical_root / f"{stem}.baseline.tsv")
        levels[mode]["refined"] = rows(historical_root / f"{stem}.refined.tsv")
        levels[mode]["ultra"] = rows(ultra_paths[mode])
        evidence["ultra_hashes"][mode] = {
            "trajectory_sha256": digest(ultra_paths[mode]),
            "steps_sha256": digest(Path(str(ultra_paths[mode]) + ".steps")),
        }
        grids = [[row["t_s"] for row in levels[mode][level]] for level in ("baseline", "refined", "ultra")]
        if not grids[0] == grids[1] == grids[2]:
            raise RuntimeError(f"{mode}: output grid mismatch")
        for level, sequence in levels[mode].items():
            expected_identity = CARD + ("-MATCHED-CONTROL" if mode == "control" else "")
            for index, row in enumerate(sequence):
                if row["run_card_identity"] != expected_identity or row["valid_through_sample"] != "1":
                    raise RuntimeError(f"{mode} {level} row {index}: run-card/frozen mismatch")
                for key, value in row.items():
                    if key not in {"frozen_limiting_quantity", "run_card_identity", "partition_identity",
                                   "source_identity", "domain_identity", "revision_identity",
                                   "product_fate_identity", "actual_potential_provenance",
                                   "finite_T_weighting_class", "regime_e", "regime_mu",
                                   "sign_observable", "sign_classification", "terminal_fate_branch_ids",
                                   "terminal_fate_channel_ids", "terminal_fate_branch_weights"}:
                        finite(value, f"{mode} {level} row {index} {key}")
        final_depletion = finite(levels[mode]["ultra"][-1]["DeltaB_over_B0"], "final depletion")
        if mode == "source" and abs(final_depletion + 5e-7) > 64 * EPS:
            raise RuntimeError("source final depletion changed")
        if mode == "control" and final_depletion != 0.0:
            raise RuntimeError("control depletion is nonzero")

    state_max_stability = (-1.0, None)
    state_max_contraction = (-1.0, None)
    state_max_floor = (-1.0, None)
    resolvable = floor_limited = 0
    state_pass = True
    for mode in ("source", "control"):
        for index, triplet in enumerate(zip(levels[mode]["baseline"], levels[mode]["refined"], levels[mode]["ultra"])):
            baseline, refined, ultra = triplet
            for component, name in enumerate(STATE):
                values = tuple(finite(row[name], name) for row in triplet)
                d_br, d_ru = abs(values[0] - values[1]), abs(values[1] - values[2])
                magnitude = max(map(abs, values))
                d_refined = ATOL_REFINED[component] + 1e-9 * magnitude
                d_ultra = ATOL_ULTRA[component] + 1e-11 * magnitude
                # All three tables use max_digits10 decimal serialization and therefore
                # round-trip to the exact original binary64 values: declared Q = 0.
                floor = max(d_ultra, 64 * math.ulp(magnitude), 0.0)
                stability = d_ru / d_refined
                detail = {"mode": mode, "component": name, "index": index,
                          "time_s": finite(refined["t_s"], "time"), "d_BR": d_br,
                          "d_RU": d_ru, "D_R": d_refined, "D_U": d_ultra, "Q": 0.0, "F": floor}
                state_max_stability = max(state_max_stability, (stability, detail), key=lambda item: item[0])
                state_pass &= stability <= 1.0
                if d_br > 10 * floor:
                    resolvable += 1
                    ratio = d_ru / d_br
                    state_max_contraction = max(state_max_contraction, (ratio, detail), key=lambda item: item[0])
                    state_pass &= ratio <= 0.10
                else:
                    floor_limited += 1
                    utilization = d_ru / (10 * floor) if floor else (0.0 if d_ru == 0 else math.inf)
                    state_max_floor = max(state_max_floor, (utilization, detail), key=lambda item: item[0])
                    state_pass &= d_ru <= 10 * floor
    evidence["state"] = {
        "Q_definition": "0: max_digits10 serialization round-trips each binary64 value exactly",
        "maximum_d_RU_over_D_R": {"value": state_max_stability[0], **state_max_stability[1]},
        "maximum_resolvable_d_RU_over_d_BR": ({"value": state_max_contraction[0], **state_max_contraction[1]}
                                                   if state_max_contraction[1] else None),
        "maximum_floor_limited_utilization": ({"value": state_max_floor[0], **state_max_floor[1]}
                                                 if state_max_floor[1] else None),
        "resolvable_comparisons": resolvable,
        "floor_limited_comparisons": floor_limited,
        "pass": state_pass,
    }

    ledger_max_stability = (-1.0, None)
    ledger_max_contraction = (-1.0, None)
    ledger_max_floor = (-1.0, None)
    ledger_pass = True
    for mode in ("source", "control"):
        for index, triplet in enumerate(zip(levels[mode]["baseline"], levels[mode]["refined"], levels[mode]["ultra"])):
            scales = [gross_power(row) for row in triplet]
            scale = max(scales)
            for name in LEDGER:
                values = tuple(finite(row[name], name) for row in triplet)
                d_br, d_ru = abs(values[0] - values[1]), abs(values[1] - values[2])
                stability = d_ru / scale
                magnitude = max(map(abs, values))
                floor = max(1e-11 * scale, 64 * math.ulp(magnitude))
                detail = {"mode": mode, "observable": name, "index": index,
                          "time_s": finite(triplet[1]["t_s"], "time"), "d_BR": d_br,
                          "d_RU": d_ru, "G_P": scale, "floor": floor}
                ledger_max_stability = max(ledger_max_stability, (stability, detail), key=lambda item: item[0])
                ledger_pass &= stability <= 1e-9
                if d_br > 10 * floor:
                    ratio = d_ru / d_br
                    ledger_max_contraction = max(ledger_max_contraction, (ratio, detail), key=lambda item: item[0])
                    ledger_pass &= ratio <= 0.10
                else:
                    utilization = d_ru / floor if floor else (0.0 if d_ru == 0 else math.inf)
                    ledger_max_floor = max(ledger_max_floor, (utilization, detail), key=lambda item: item[0])
                    ledger_pass &= d_ru <= floor
    evidence["ledger"] = {
        "maximum_normalized_refined_ultra_difference": {"value": ledger_max_stability[0], **ledger_max_stability[1]},
        "maximum_resolvable_contraction": ({"value": ledger_max_contraction[0], **ledger_max_contraction[1]}
                                              if ledger_max_contraction[1] else None),
        "maximum_floor_utilization": ({"value": ledger_max_floor[0], **ledger_max_floor[1]}
                                         if ledger_max_floor[1] else None),
        "pass": ledger_pass,
    }

    endpoint_pass = True
    for mode in ("source", "control"):
        closures = {level: r20(levels[mode][level]) for level in ("refined", "ultra")}
        scale = max(closures["refined"]["N_R20_erg"], closures["ultra"]["N_R20_erg"])
        differences = {}
        for name in ("DeltaEeq_erg", "DeltaEchem_erg", "DeltaUth_erg"):
            value = abs(closures["refined"][name] - closures["ultra"][name]) / scale
            differences[name + "_difference_over_max_N"] = value
            endpoint_pass &= value <= 1e-9
        propagation = abs(sum(closures["refined"][name] for name in
                              ("DeltaEeq_erg", "DeltaEchem_erg", "DeltaUth_erg")) -
                          sum(closures["ultra"][name] for name in
                              ("DeltaEeq_erg", "DeltaEchem_erg", "DeltaUth_erg"))) / scale
        differences["endpoint_state_propagation_over_max_N"] = propagation
        endpoint_pass &= propagation <= 5e-5
        for closure in closures.values():
            endpoint_pass &= (closure["R20_over_N"] <= 2e-4 and
                              closure["luminosity_quadrature_over_N"] <= 5e-5 and
                              closure["thermal_quadrature_over_N"] <= 5e-5)
        evidence["endpoint"][mode] = {"tiers": closures, **differences}
    evidence["endpoint"]["pass"] = endpoint_pass

    identity_pass = True
    for mode in ("source", "control"):
        mode_result = {}
        for level in ("refined", "ultra"):
            maxima = {name: {"absolute": -1.0, "utilization": -1.0, "index": None, "time_s": None}
                      for name in ("R18_residual_erg_s", "Ra_Rb_residual_erg_s", "Rb_Rc_residual_erg_s")}
            for index, row in enumerate(levels[mode][level]):
                scale = max(1.0, abs(finite(row["P_dir_actual_erg_s"], "Pdir actual")),
                            abs(finite(row["P_dir_eq_erg_s"], "Pdir eq")))
                budget = 64 * EPS * scale
                for name in maxima:
                    value = abs(finite(row[name], name))
                    utilization = value / budget
                    if utilization > maxima[name]["utilization"]:
                        maxima[name] = {"absolute": value, "budget": budget,
                                        "utilization": utilization, "index": index,
                                        "time_s": finite(row["t_s"], "time")}
                    identity_pass &= value <= budget
            mode_result[level] = maxima
        evidence["identity_budgets"][mode] = mode_result
    evidence["identity_budgets"]["pass"] = identity_pass

    frozen_max = -1.0
    max_depletion = 0.0
    for mode in ("source", "control"):
        for level in ("refined", "ultra"):
            for row in levels[mode][level]:
                frozen_max = max(frozen_max, finite(row["max_frozen_utilization"], "frozen utilization"))
                max_depletion = max(max_depletion, abs(finite(row["DeltaB_over_B0"], "depletion")))
    frozen_pass = frozen_max <= 1.0 and max_depletion <= 1e-6
    evidence["frozen"] = {"maximum_utilization": frozen_max,
                          "maximum_abs_DeltaB_over_B0": max_depletion, "pass": frozen_pass}

    solver_pass = True
    for mode in ("source", "control"):
        trend = step_trend(Path(str(ultra_paths[mode]) + ".steps"))
        evidence["solver_trend"][mode] = trend
        solver_pass &= trend["pass"]
    evidence["solver_trend"]["pass"] = solver_pass

    matched_pass = True
    matched_max_stability = (-1.0, None)
    matched_max_contraction = (-1.0, None)
    matched_max_floor = (-1.0, None)
    source_minus_control: dict[str, dict[str, list[float]]] = {}
    r20_scales = {}
    for level in ("baseline", "refined", "ultra"):
        source, control = levels["source"][level], levels["control"][level]
        source_u, control_u = cumulative_thermal_energy(source), cumulative_thermal_energy(control)
        source_power = cumulative_trapezoid(source, lambda row: finite(row["Pnet_erg_s"], "Pnet"))
        control_power = cumulative_trapezoid(control, lambda row: finite(row["Pnet_erg_s"], "Pnet"))
        source_minus_control[level] = {
            "Tinf_K": [finite(a["Tinf_K"], "Tinf") - finite(b["Tinf_K"], "Tinf") for a, b in zip(source, control)],
            "Tsurface_inf_K": [finite(a["Tsurface_inf_K"], "Tsurface") - finite(b["Tsurface_inf_K"], "Tsurface") for a, b in zip(source, control)],
            "Lgamma_erg_s": [finite(a["Lgamma_erg_s"], "Lgamma") - finite(b["Lgamma_erg_s"], "Lgamma") for a, b in zip(source, control)],
            "Uth_erg": [a - b for a, b in zip(source_u, control_u)],
            "integrated_total_power_erg": [a - b for a, b in zip(source_power, control_power)],
        }
        r20_scales[level] = max(r20(source)["N_R20_erg"], r20(control)["N_R20_erg"])
    for name in ("Tinf_K", "Tsurface_inf_K", "Lgamma_erg_s", "Uth_erg", "integrated_total_power_erg"):
        for index in range(EXPECTED_ROWS):
            values = tuple(source_minus_control[level][name][index] for level in ("baseline", "refined", "ultra"))
            d_br, d_ru = abs(values[0] - values[1]), abs(values[1] - values[2])
            magnitude = max(map(abs, values))
            if name in ("Lgamma_erg_s",):
                scale = max(gross_power(levels["source"][level][index]) for level in ("baseline", "refined", "ultra"))
                stability_budget, floor = 1e-9 * scale, max(1e-11 * scale, 64 * math.ulp(magnitude))
            elif name in ("Uth_erg", "integrated_total_power_erg"):
                scale = max(r20_scales.values())
                stability_budget, floor = 1e-9 * scale, max(1e-11 * scale, 64 * math.ulp(magnitude))
            else:
                parents = [finite(levels[mode][level][index][name], name)
                           for mode in ("source", "control") for level in ("baseline", "refined", "ultra")]
                cancellation_floor = 64 * max(math.ulp(abs(value)) for value in parents)
                stability_budget = max(1e-9 * magnitude, cancellation_floor)
                floor = max(1e-11 * magnitude, cancellation_floor)
            utilization = d_ru / stability_budget if stability_budget else (0.0 if d_ru == 0 else math.inf)
            detail = {"observable": name, "index": index,
                      "time_s": finite(levels["source"]["refined"][index]["t_s"], "time"),
                      "d_BR": d_br, "d_RU": d_ru, "stability_budget": stability_budget, "floor": floor}
            matched_max_stability = max(matched_max_stability, (utilization, detail), key=lambda item: item[0])
            matched_pass &= d_ru <= stability_budget
            if d_br > 10 * floor:
                ratio = d_ru / d_br
                matched_max_contraction = max(matched_max_contraction, (ratio, detail), key=lambda item: item[0])
                matched_pass &= ratio <= 0.10
            else:
                floor_utilization = d_ru / floor if floor else (0.0 if d_ru == 0 else math.inf)
                matched_max_floor = max(matched_max_floor, (floor_utilization, detail), key=lambda item: item[0])
                matched_pass &= d_ru <= floor
    evidence["matched_source_minus_control"] = {
        "maximum_stability_utilization": {"value": matched_max_stability[0], **matched_max_stability[1]},
        "maximum_resolvable_contraction": ({"value": matched_max_contraction[0], **matched_max_contraction[1]}
                                              if matched_max_contraction[1] else None),
        "maximum_floor_utilization": ({"value": matched_max_floor[0], **matched_max_floor[1]}
                                         if matched_max_floor[1] else None),
        "heating_cooling_classification_made": False,
        "pass": matched_pass,
    }

    overall = state_pass and ledger_pass and endpoint_pass and identity_pass and frozen_pass and solver_pass and matched_pass
    evidence["BA12R"] = "PASS" if overall else "FAIL"
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"BA12R {evidence['BA12R']}")
    print(f"state_stability_max={state_max_stability[0]:.17g} state_contraction_max={state_max_contraction[0]:.17g}")
    print(f"ledger_normalized_max={ledger_max_stability[0]:.17g} endpoint_pass={endpoint_pass}")
    print(f"solver_pass={solver_pass} matched_pass={matched_pass} frozen_max={frozen_max:.17g}")
    return 0 if overall else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"STOP {error}", file=sys.stderr)
        raise SystemExit(1)
