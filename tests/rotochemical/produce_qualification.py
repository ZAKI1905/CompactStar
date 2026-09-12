#!/usr/bin/env python3
"""Produce the radial-80000/EOS-8192 qualification entirely from fresh outputs."""

import argparse
import copy
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_logged(command, log_path, cwd):
    with log_path.open("w") as log:
        result = subprocess.run(
            [str(x) for x in command], cwd=cwd, stdout=log,
            stderr=subprocess.STDOUT
        )
    (log_path.with_suffix(".rc")).write_text(str(result.returncode) + "\n")
    if result.returncode:
        raise RuntimeError(f"fresh command failed rc={result.returncode}: {command[0]}")
    return result.returncode


def read_rows(path):
    with Path(path).open() as stream:
        return [
            {
                key: int(value) if key in {"resolution", "species"} else float(value)
                for key, value in row.items()
            }
            for row in csv.DictReader(stream, delimiter="\t")
        ]


def require_pass(path):
    text = Path(path).read_text()
    if "STOP " in text or "PASS" not in text:
        raise RuntimeError("fresh Phase-5B evidence did not pass: " + str(path))
    return text


def parse_production(path):
    matrix_names = {
        "G", "E_G", "E_quadrature", "Q", "E_Q", "E_Schur_arithmetic",
        "Z", "E_Z", "E_global_solve", "E_Z_arithmetic",
    }
    values = {}
    provenance = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("PROVENANCE "):
            _, key, value = line.split(" ", 2)
            provenance[key] = value
        elif line.startswith("RESULT "):
            _, key, count, *raw = line.split()
            count = int(count)
            data = list(map(float, raw))
            if key in matrix_names or key.startswith((
                "G_ladder_", "E_local_", "E_center", "E_tail",
                "E_background", "E_refusal",
            )):
                values[key] = np.array(data).reshape(count, count).tolist()
            else:
                values[key] = data
    return values, provenance


def reconstruct_structural(output):
    rows = read_rows(output / "structural-m1m2/radial.tsv")
    by_key = {(row["resolution"], row["species"]): row for row in rows}
    finest = [by_key[(80000, species)] for species in range(4)]
    middle = [by_key[(40000, species)] for species in range(4)]
    for resolution in (20000, 40000, 80000):
        if any((resolution, species) not in by_key for species in range(4)):
            raise RuntimeError("incomplete structural resolution ladder")

    pb6_text = require_pass(output / "pb6.log")
    knot_relative = np.zeros(4)
    for species, value in re.findall(r"PB6 knot species=(\d+) relative=(\S+)", pb6_text):
        knot_relative[int(species)] = float(value)
    if np.any(knot_relative <= 0):
        raise RuntimeError("incomplete fresh PB6 knot evidence")
    a_error = np.array([abs(row["A"]) for row in finest]) * knot_relative
    a_b_error = a_error[0] + a_error[1]
    b_values = np.array([row["B"] for row in finest])
    b_b = finest[0]["B_B"]
    pb6 = a_error + np.abs(b_values) * a_b_error / abs(b_b)

    pb7_text = require_pass(output / "pb7.log")
    b_oracle = np.zeros(4)
    for species, value in re.findall(
        r"PB7 species=(\d+) sequence_B=\S+ homogeneous_B=(\S+)", pb7_text
    ):
        b_oracle[int(species)] = float(value)
    if np.any(b_oracle == 0):
        raise RuntimeError("incomplete fresh PB7 evidence")
    difference = np.abs(b_values - b_oracle)
    difference_b = difference[0] + difference[1]
    denominator = abs(b_b) - difference_b
    pb7 = abs(finest[0]["A_B"]) * (
        difference / denominator
        + np.abs(b_values) * difference_b / (abs(b_b) * denominator)
    )

    pb11_text = require_pass(output / "pb9-11.log")
    quotients = {}
    current_q = None
    delta_b_at_fine = None
    fine_q, coarse_q = 1.25e-7, 2.5e-7
    for line in pb11_text.splitlines():
        match = re.search(r"PB11 q=(\S+).*Delta_N_B=(\S+)", line)
        if match:
            current_q = float(match.group(1))
            if np.isclose(current_q, fine_q, rtol=0, atol=1e-22):
                delta_b_at_fine = float(match.group(2))
            continue
        match = re.search(r"PB11 species=(\d+) quotient=(\S+)", line)
        if match and current_q is not None:
            quotients[(current_q, int(match.group(1)))] = float(match.group(2))
    if delta_b_at_fine is None:
        raise RuntimeError("fresh PB11 fine-q baryon residual missing")
    k_values = np.array([row["K"] for row in finest])
    pb11_direct = np.array([
        abs(2 * quotients[(fine_q, species)]
            - quotients[(coarse_q, species)] - k_values[species])
        for species in range(4)
    ])
    residual_b = abs(delta_b_at_fine) / fine_q
    pb11_constraint = np.abs(b_values / b_b) * residual_b

    pb12_text = require_pass(output / "pb12.log")
    pb12_values = {
        (int(resolution), int(species)): float(value)
        for resolution, species, value in re.findall(
            r"PB12 EOS_resolution=(\d+) species=(\d+) K=(\S+)", pb12_text
        )
    }
    pb12 = np.array([
        abs(pb12_values[(8192, species)] - pb12_values[(4096, species)])
        for species in range(4)
    ])
    require_pass(output / "pb13.log")
    pb10_rows = sorted(
        read_rows(output / "structural-m1m2/pb10.tsv"),
        key=lambda row: row["species"],
    )
    pb10 = np.array([abs(row["Delta_K"]) for row in pb10_rows])
    direct_radial = np.abs(k_values - np.array([row["K"] for row in middle]))
    representation = np.maximum(direct_radial, pb6)
    v_k = representation + pb7 + pb11_direct + pb11_constraint + pb12 + pb10
    conversion = np.abs(np.array([row["I_phys"] for row in finest]) / k_values)
    v_i_all = v_k * conversion
    reconstructed_e_i = np.array(
        [row["K_numerical_error"] for row in finest]
    ) * conversion
    emitted_e_i = np.array([row["E_I_numerical"] for row in finest])
    if not np.allclose(
        reconstructed_e_i, emitted_e_i,
        rtol=4 * np.finfo(float).eps, atol=0,
    ):
        raise RuntimeError("fresh K-error to I-error conversion mismatch")
    return {
        "rows": rows,
        "pb10_rows": pb10_rows,
        "V_K_validation": v_k.tolist(),
        "V_I_validation_all_species": v_i_all.tolist(),
        "V_I_validation_consumed_e_mu": v_i_all[[2, 3]].tolist(),
        "E_I_numerical_all_species": emitted_e_i.tolist(),
        "E_I_numerical_consumed_e_mu": emitted_e_i[[2, 3]].tolist(),
        "conversion_I_per_K": conversion.tolist(),
        "ingredients": {
            "direct_radial": direct_radial.tolist(),
            "PB6": pb6.tolist(),
            "representation_max": representation.tolist(),
            "PB7": pb7.tolist(),
            "PB11_direct": pb11_direct.tolist(),
            "PB11_constraint": pb11_constraint.tolist(),
            "PB12": pb12.tolist(),
            "PB10": pb10.tolist(),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--trackr-executable", type=Path, required=True)
    parser.add_argument("--structural-executable", type=Path, required=True)
    parser.add_argument("--phase5b-executable", type=Path, required=True)
    parser.add_argument("--production-executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError("qualification output must be fresh")
    output.mkdir(parents=True)
    analysis = source_root / "tests/analysis"
    sys.path.insert(0, str(analysis))
    from chemical_production_evidence import (  # noqa: E402
        Model, certificates, integrate, matrix, tail, transport,
    )

    run_logged(
        [sys.executable, analysis / "chemical_trackr_budget.py",
         args.trackr_executable.resolve(), output / "chemical-characterization"],
        output / "chemical-characterization.log", source_root,
    )
    run_logged(
        [args.structural_executable.resolve(), output / "structural-m1m2"],
        output / "structural-m1m2.log", source_root,
    )
    for stage, label in [
        ("PB6", "pb6"), ("PB7", "pb7"), ("PB9-11", "pb9-11"),
        ("PB12", "pb12"), ("PB13", "pb13"),
    ]:
        run_logged(
            [args.phase5b_executable.resolve(), stage, output / label],
            output / (label + ".log"), source_root,
        )

    structural = reconstruct_structural(output)
    goals_path = source_root / "docs/validation/phase5c2_preproduction_evidence.json"
    goals = json.loads(goals_path.read_text())["goals"]
    fixed = {
        "goals": goals,
        "V_I_validation": structural["V_I_validation_consumed_e_mu"],
    }
    result_paths = list(
        (output / "chemical-characterization").glob("run-*/result.json")
    )
    if len(result_paths) != 1:
        raise RuntimeError("expected exactly one fresh chemical characterization")
    characterized_path = result_paths[0]
    characterized = json.loads(characterized_path.read_text())
    characterization_root = characterized_path.parent
    outputs = {}
    for fixture, resolution, mode in [
        ("t4096-r40000", 40000, "global"),
        ("t8192-r40000", 40000, "global"),
        ("t8192-r80000", 80000, "full"),
    ]:
        directory = characterization_root / fixture
        run_characterization = copy.deepcopy(characterized)
        if fixture != "t8192-r80000":
            data = np.loadtxt(directory / "profile.tsv", skiprows=1)
            model = Model(directory)
            run_characterization["tail"] = tail(
                data, model, matrix(integrate(data, model, 16)[0])
            )[0]
        certificate, _ = certificates(directory, run_characterization, goals)
        certificate_text = output / f"certificate-{fixture}.txt"
        certificate_json = output / f"certificate-{fixture}.json"
        transport(certificate, certificate_text, fixed)
        certificate_json.write_text(
            json.dumps(certificate, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        production_dir = output / f"production-{fixture}"
        production_log = output / f"production-{fixture}.log"
        run_logged(
            [args.production_executable.resolve(), directory, certificate_text,
             production_dir, resolution, mode],
            production_log, source_root,
        )
        values, provenance = parse_production(production_log)
        outputs[fixture] = {
            "mode": mode,
            "radial_resolution": resolution,
            "values": values,
            "provenance": provenance,
            "certificate_sha256": digest(certificate_json),
            "certificate_transport_sha256": digest(certificate_text),
            "log_sha256": digest(production_log),
        }

    lower = np.asarray(outputs["t4096-r40000"]["values"]["G"])
    middle = np.asarray(outputs["t8192-r40000"]["values"]["G"])
    fine = np.asarray(outputs["t8192-r80000"]["values"]["G"])
    crosscheck = 2 * (np.abs(middle - lower) + np.abs(fine - middle))
    component = np.asarray(outputs["t8192-r80000"]["values"]["E_background"])
    if not np.all(crosscheck <= component):
        raise RuntimeError("fresh background ladder exceeds transported error")
    fine_values = outputs["t8192-r80000"]["values"]
    if not np.array_equal(
        np.asarray(fine_values["E_I_numerical"]),
        np.asarray(structural["E_I_numerical_consumed_e_mu"]),
    ):
        raise RuntimeError("production E_I differs from fresh reconstruction")
    z = np.asarray(fine_values["Z"])
    ez = np.asarray(fine_values["E_Z"])
    physical_i = np.asarray(fine_values["I"])
    numerical_i = np.asarray(fine_values["E_I_numerical"])
    w_error = (
        np.abs(z) @ numerical_i + ez @ np.abs(physical_i)
        + ez @ numerical_i + np.asarray(fine_values["E_W_arithmetic"])
    )
    if not np.allclose(
        w_error, np.asarray(fine_values["E_W_numerical"]),
        rtol=4 * np.finfo(float).eps, atol=0,
    ):
        raise RuntimeError("W error reconstruction mismatch")
    evidence = {
        "schema": "phase5d-fresh-qualification-v1",
        "fixture": {"EOS_resolution": 8192, "radial_resolution": 80000,
                    "rho_c_g_cm3": 1.10e15},
        "freshness": {
            "baseline_or_candidate_read": False,
            "stored_V_I_read": False,
            "prior_generated_state_read": False,
        },
        "goals_source": {
            "logical_path": "docs/validation/phase5c2_preproduction_evidence.json",
            "sha256": digest(goals_path), "accessed_field": "goals",
        },
        "structural": structural,
        "outputs": outputs,
        "background_crosscheck": crosscheck.tolist(),
        "W_error_reconstructed": w_error.tolist(),
        "logical_hashes": {
            "profile": digest(characterization_root / "t8192-r80000/profile.tsv"),
            "freegas": digest(characterization_root / "t8192-r80000/freegas.tsv"),
            "model": digest(characterization_root / "t8192-r80000/model.txt"),
            "certificate": digest(output / "certificate-t8192-r80000.txt"),
        },
    }
    target = output / "qualification-evidence.json"
    target.write_text(
        json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("FRESH QUALIFICATION COMPLETE", digest(target), flush=True)


if __name__ == "__main__":
    main()
