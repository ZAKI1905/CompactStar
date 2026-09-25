#!/usr/bin/env python3
"""Authenticate and verify the bounded ADR-0017 production qualification."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import re
import struct
import subprocess
from pathlib import Path


EXPECTED = {
    "matrix": "32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616",
    "oracle_root": "21f1ff9eaf23b078f86aa2b10ec45be55f84a98b9aedc4b94f4a7c4bf999f866",
    "oracle1_tree": "89bc3cb0bde7e58f60e1e44fdee9332e920ae163870213fe79d7bb81ba65ecf9",
    "oracle2_tree": "f0bf0fb4e1e0acef58da1d770c0a71d87c4021264017f75cbfa8aa9da1a0af2d",
    "oracle1_ledger": "0fec6717e153ce19ecafda0fe0ed5c01d1be3934342ac9776f05ad2da769fbcc",
    "oracle2_ledger": "b1bc3feb6bce96b7e4f1feecf2660b1f88781941d66d25ae832e4103a3248f74",
    "oracle_result": "f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0",
    "profile_tree": "233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72",
    "profile.tsv": "e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890",
    "model.txt": "3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4",
    "freegas.tsv": "7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a",
    "certificate": "7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453",
    "thermal": "1cdb98507958b50fbd1e1eff3eb3bfa4f48c790fef86a8740855884fc43b9d94",
    "frozen": "9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab",
    "coefficients": "3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6",
    "entry": "5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d",
    "schedule": "43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e",
    "accepted": "7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459",
    "brackets": "0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d",
    "internal": "fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8",
    "trajectory": "8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c",
    "steps": "912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(path: Path) -> str:
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


def command(*args: str) -> str:
    return subprocess.check_output(args, text=True).strip()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bits(value: str | float) -> bytes:
    return struct.pack(">d", float(value))


def exact(a: str | float, b: str | float, message: str) -> None:
    require(bits(a) == bits(b), message + f": {a} != {b}")


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def authenticate(args: argparse.Namespace) -> None:
    oracle_root, matrix = Path(args.oracle_root), Path(args.matrix)
    profile, output = Path(args.profile), Path(args.output)
    require(not output.exists(), "fresh authentication output root required")
    output.mkdir(parents=True)
    observed = {
        "matrix": sha256(matrix),
        "oracle_root": tree_sha256(oracle_root),
        "oracle1_tree": tree_sha256(oracle_root / "oracle1"),
        "oracle2_tree": tree_sha256(oracle_root / "oracle2"),
        "oracle1_ledger": sha256(oracle_root / "oracle1" / "solve_ledger.tsv"),
        "oracle2_ledger": sha256(oracle_root / "oracle2" / "solve_ledger.tsv"),
        "oracle_result": sha256(Path(args.oracle_result)),
        "profile_tree": tree_sha256(profile),
        "profile.tsv": sha256(profile / "profile.tsv"),
        "model.txt": sha256(profile / "model.txt"),
        "freegas.tsv": sha256(profile / "freegas.tsv"),
        "certificate": sha256(Path(args.certificate)),
        "thermal": tree_sha256(Path(args.thermal)),
        "frozen": sha256(Path(args.frozen)),
        "coefficients": sha256(Path(args.coefficients)),
        "entry": sha256(Path(args.entry)),
    }
    require(observed == {key: EXPECTED[key] for key in observed}, "immutable/oracle hash mismatch")
    require(len(rows(matrix)) == 240, "solve matrix row count changed")
    for tier in ("oracle1", "oracle2"):
        batch = oracle_root / tier
        result_files = [p for p in batch.glob("obs-*.tsv") if re.fullmatch(r"obs-\d{3}\.tsv", p.name)]
        meta_files = list(batch.glob("obs-*.meta.tsv"))
        require(len(result_files) == 241 and len(meta_files) == 241, f"{tier} observation file count changed")
        ledger = rows(batch / "solve_ledger.tsv")
        require(len(ledger) == 239 and len({r["solve_id"] for r in ledger}) == 239, f"{tier} ledger count changed")
        result_paths = [batch / f"obs-{int(r['observation_index']):03d}.tsv" for r in ledger]
        meta_paths = [batch / f"obs-{int(r['observation_index']):03d}.meta.tsv" for r in ledger]
        require(all(path.is_file() for path in result_paths + meta_paths), f"{tier} result missing")
        require(all(sha256(path) == row["result_sha256"] for path, row in zip(result_paths, ledger)),
                f"{tier} result hash mismatch")
        require(all(sha256(path) == row["meta_sha256"] for path, row in zip(meta_paths, ledger)),
                f"{tier} meta hash mismatch")
    platform_evidence = {
        "macos": command("sw_vers", "-productVersion"),
        "build": command("sw_vers", "-buildVersion"),
        "darwin": command("uname", "-r"),
        "architecture": platform.machine(),
        "compiler": command("clang++", "--version").splitlines()[0],
        "gsl": command("gsl-config", "--version"),
        "cmake": command("cmake", "--version").splitlines()[0].split()[-1],
    }
    require(platform_evidence == {
        "macos": "26.6.2", "build": "25G83", "darwin": "25.6.0", "architecture": "arm64",
        "compiler": "Apple clang version 21.0.0 (clang-2100.3.34.2)", "gsl": "2.7.1", "cmake": "4.2.1",
    }, "qualification platform mismatch")
    evidence = {"pass": True, "hashes": observed, "platform": platform_evidence,
                "oracle_integrations_reused": 478, "new_integrations": 0}
    atomic_json(output / "authentication.json", evidence)
    (output / "oracle-qualified.flag").write_text("ADR0017_ORACLE_AUTHENTICATED\n")
    print("ADR0017 ORACLE AUTHENTICATION PASS")


FLOAT_MAP = {
    "P_dir_eq": "P_dir_eq_erg_s", "P_dir_actual": "P_dir_actual_erg_s", "LH": "LH_erg_s",
    "DeltaLnu": "DeltaLnu_erg_s", "DeltaPbeta": "DeltaPbeta_erg_s", "Lnu_eq": "Lnu_eq_erg_s",
    "Lnu_full": "Lnu_full_erg_s", "Lgamma": "Lgamma_erg_s", "Lother": "Lother_erg_s",
    "Pnet": "Pnet_erg_s", "mu_B": "mu_B_inf_MeV", "mu_n_actual": "mu_n_actual_inf_MeV",
    "sigma_e": "sigma_e_count_s", "sigma_mu": "sigma_mu_count_s", "Echem": "Echem_MeV",
    "Cstar": "Cstar_erg_K", "Tinf": "Tinf_K", "B": "B_count",
}


def oracle_row(root: Path, tier: str, observation: int) -> dict[str, str]:
    path = root / tier / f"obs-{observation:03d}.tsv"
    data = rows(path)
    require(len(data) == 1, f"{tier} observation {observation} row count")
    return data[0]


def compute_oracle_r20(oracle2: list[dict[str, str]]) -> dict[str, float]:
    mev_to_erg = 1.602176634e-6
    first, last = oracle2[0], oracle2[-1]
    delta_eq = mev_to_erg * float(first["mu_B_inf_MeV"]) * (float(last["B_count"]) - float(first["B_count"]))
    delta_chem = mev_to_erg * (float(last["Echem_MeV"]) - float(first["Echem_MeV"]))
    delta_uth = outgoing = outgoing_abs = 0.0
    for a, b in zip(oracle2, oracle2[1:]):
        dt = float(b["t_s"]) - float(a["t_s"])
        delta_uth += 0.5 * (float(a["Cstar_erg_K"]) + float(b["Cstar_erg_K"])) * (float(b["Tinf_K"]) - float(a["Tinf_K"]))
        pa = float(a["L_out_fluid_inf_erg_s"]) + float(a["Lnu_full_erg_s"])
        pa += float(a["Lgamma_erg_s"])
        pa += float(a["Lother_erg_s"])
        pb = float(b["L_out_fluid_inf_erg_s"]) + float(b["Lnu_full_erg_s"])
        pb += float(b["Lgamma_erg_s"])
        pb += float(b["Lother_erg_s"])
        outgoing += 0.5 * dt * (pa + pb)
        aa = abs(float(a["Lnu_full_erg_s"])) + abs(float(a["Lgamma_erg_s"]))
        aa += abs(float(a["Lother_erg_s"]))
        ab = abs(float(b["Lnu_full_erg_s"])) + abs(float(b["Lgamma_erg_s"]))
        ab += abs(float(b["Lother_erg_s"]))
        outgoing_abs += 0.5 * dt * (aa + ab)
    residual = delta_eq + delta_chem + delta_uth + outgoing
    normalizer = max(1.0, abs(delta_uth), abs(delta_chem), outgoing_abs)
    return {"R20_residual_erg": residual, "N_R20_erg": normalizer, "R20_normalized": residual / normalizer}


def verify(args: argparse.Namespace) -> None:
    output, oracle_root = Path(args.output), Path(args.oracle_root)
    production = rows(output / "checkpoints.tsv")
    require(len(production) == 241, "production checkpoint row count changed")
    main_hashes = {
        "schedule": sha256(output / "schedule.tsv"),
        "accepted": sha256(output / "main.accepted_states.tsv"),
        "brackets": sha256(output / "main.observations.tsv"),
        "internal": sha256(output / "main.internal_steps.tsv"),
        "trajectory": sha256(output / "main.tsv"),
        "steps": sha256(output / "main.steps"),
    }
    require(main_hashes == {key: EXPECTED[key] for key in main_hashes}, "production main Arm-E bytes differ")
    audit = rows(output / "main.audit.tsv")
    require(len(audit) == 1 and audit[0]["unique_positive_t1_targets"] == "1"
            and audit[0]["t1_s"] == "462269531250"
            and audit[0]["intermediate_observation_t1_matches"] == "0", "main target audit failed")

    oracle2_rows: list[dict[str, str]] = []
    exact_endpoint = strict = failures = knots = 0
    max_d = (-1.0, None, None)
    max_u = (-1.0, None, None)
    components = (("x", "x_state"), ("eta_e", "eta_e_MeV"), ("eta_mu", "eta_mu_MeV"))
    identity_names = ("run_card_identity", "source_identity", "domain_identity", "revision_identity",
                      "partition_identity", "product_fate_identity")
    for index, row in enumerate(production):
        require(int(row["observation_index"]) == index, "production observation ordering changed")
        o1, o2 = oracle_row(oracle_root, "oracle1", index), oracle_row(oracle_root, "oracle2", index)
        oracle2_rows.append(o2)
        exact(row["t_obs_s"], o2["t_s"], f"observation {index} target")
        if index > 0 and row["source"] == "MAIN_ENDPOINT":
            exact_endpoint += 1
            require(row["rk8pd_invocations"] == "0", "exact endpoint invoked rk8pd")
        if row["source"] == "RK8PD_RECONSTRUCTED":
            strict += 1
            require(row["rk8pd_invocations"] == "2", "strict interior did not invoke exactly O1/O2")
            for short, oracle_name in components:
                exact(row[f"{short}_O1"], o1[oracle_name], f"O1 state {index} {short}")
                exact(row[f"{short}_O2"], o2[oracle_name], f"O2 state {index} {short}")
                d_ratio = float(row[f"d_{short}"]) / float(row[f"D_O1_{short}"])
                u_ratio = float(row[f"U_{short}"]) / (0.20 * float(row[f"F_{short}"]))
                if d_ratio > max_d[0]: max_d = (d_ratio, index, short)
                if u_ratio > max_u[0]: max_u = (u_ratio, index, short)
        require(row["status"] == "QUALIFIED" and row["self_qualified"] == "1"
                and row["diagnostic_self_qualified"] == "1", f"checkpoint {index} unresolved")
        failures += row["status"] != "QUALIFIED"
        knots += index > 0 and row["category"] == "B"
        for short, oracle_name in FLOAT_MAP.items():
            exact(row[f"{short}_O1"], o1[oracle_name], f"O1 diagnostic {index} {short}")
            exact(row[f"{short}_O2"], o2[oracle_name], f"O2 diagnostic {index} {short}")
        for identity in identity_names:
            require(row[identity.removesuffix("_identity") if identity == "run_card_identity" else identity]
                    == o2[identity], f"identity mismatch {index} {identity}")

    require(exact_endpoint == 1 and strict == 239 and failures == 0 and knots == 3,
            "checkpoint classification counts changed")
    historical = json.loads(Path(args.oracle_result).read_text())
    exact(max_d[0], historical["maximum_d_O_over_D_O1"]["value"], "maximum d/D_O1")
    exact(max_u[0], historical["maximum_U_O_over_0p20F"]["value"], "maximum U/(0.20F)")
    require(max_d[1:] == (117, "x") and max_u[1:] == (117, "x"), "self-qualification maximum identity changed")

    performance = {row["key"]: row["value"] for row in rows(output / "performance.tsv")}
    oracle_r20 = compute_oracle_r20(oracle2_rows)
    for name, value in oracle_r20.items():
        exact(performance[name], value, name)
    uncertainty = float(performance["R20_reconstruction_uncertainty_erg"])
    require(math.isfinite(uncertainty) and uncertainty >= 0
            and uncertainty <= 5.0e-6 * float(performance["N_R20_erg"]), "R20 uncertainty budget failed")
    require(performance["main_run_count"] == "1" and performance["process_concurrency"] == "1",
            "run/concurrency accounting changed")
    result = {
        "pass": True,
        "classification": "ADR-0017 PRODUCTION IMPLEMENTATION BOUNDED QUALIFICATION",
        "main_hashes": main_hashes,
        "main_final_state_exact": True,
        "main_accepted": 232,
        "main_rejected": 60,
        "positive_gsl_t1_targets": 1,
        "exact_endpoint_observations": exact_endpoint,
        "strict_interior_qualified": strict,
        "strict_interior_failures": failures,
        "one_knot_qualified": knots,
        "max_d_O_over_D_O1": {"value": max_d[0], "observation": max_d[1], "component": max_d[2]},
        "max_U_O_over_0p20F": {"value": max_u[0], "observation": max_u[1], "component": max_u[2]},
        "production_oracle_state_bit_identity": True,
        "production_oracle_diagnostic_bit_identity": True,
        "R20_oracle_bit_identity": True,
        "R20": oracle_r20,
        "R20_reconstruction_uncertainty_erg": uncertainty,
        "performance": {key: float(value) for key, value in performance.items()},
        "new_main_integrations": 1,
        "new_rk8pd_integrations": 478,
        "BA12_rerun": False,
        "BA12R_rerun": False,
        "future_six_run_campaign": False,
    }
    atomic_json(Path(args.result), result)
    print("ADR0017 PRODUCTION VERIFICATION PASS")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    commands = result.add_subparsers(dest="command", required=True)
    auth = commands.add_parser("authenticate")
    for name in ("oracle_root", "oracle_result", "matrix", "profile", "certificate", "thermal", "frozen", "coefficients", "entry", "output"):
        auth.add_argument(f"--{name.replace('_', '-')}", required=True)
    auth.set_defaults(function=authenticate)
    final = commands.add_parser("verify")
    final.add_argument("--output", required=True)
    final.add_argument("--oracle-root", required=True)
    final.add_argument("--oracle-result", required=True)
    final.add_argument("--result", required=True)
    final.set_defaults(function=verify)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    try:
        arguments.function(arguments)
    except Exception as error:
        raise SystemExit(f"STOP {error}")
