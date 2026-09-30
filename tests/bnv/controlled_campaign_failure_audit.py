#!/usr/bin/env python3
"""Read-only audit of a stopped clean campaign; never invokes an ODE solver.

Accepts raw output or the deterministic gzip evidence archive. Recomputes the
ADR-0017 component budgets from the recorded O1/O2 and main bracket endpoints.
This verifies evidence arithmetic, not the numerical validity of a failed state.
"""
import argparse
import csv
import ctypes
import gzip
import io
import json
import math
from pathlib import Path
import sys

# The qualified AppleClang build contracts a + r*M in this translation unit.
# Reproduce that operation explicitly for exact binary64 budget comparison;
# no tolerance is added to either acceptance inequality.
fma = ctypes.CDLL(None).fma
fma.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
fma.restype = ctypes.c_double


def read(root, name):
    path = root / name
    return (path.read_text() if path.exists() else
            gzip.decompress(Path(str(path) + ".gz").read_bytes()).decode())


def rows(root, name):
    return list(csv.DictReader(io.StringIO(read(root, name)), delimiter="\t"))


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def audit(root):
    journal = json.loads(read(root, "execution.json"))
    need([r["name"] for r in journal["runs"]] == ["source-baseline", "control-baseline"],
         "unexpected executed run inventory")
    need(journal["runs"][0]["exit_code"] == 1 and
         journal["runs"][1]["cancelled_after_peer_failure"], "missing fail-closed stop")
    source = root / "source-baseline"
    main = rows(source, "main.accepted_states.tsv")
    observations = rows(source, "main.observations.tsv")
    checkpoints = rows(source, "checkpoints.tsv")
    failures = []
    for index, q in enumerate(checkpoints):
        need(int(q["observation_index"]) == index, "checkpoint order mismatch")
        if q["source"] == "MAIN_ENDPOINT":
            need(index == 0 and q["status"] == "QUALIFIED", "unexpected exact endpoint")
            continue
        bracket = observations[index]
        left = int(bracket["previous_accepted_step"])
        right = int(bracket["new_accepted_step"])
        need(left + 1 == right and left > 0, "invalid interior bracket")
        t = float(q["t_obs_s"])
        need(float(bracket["t_previous_s"]) < t < float(bracket["t_new_s"]),
             "checkpoint is not strictly interior")
        need(t == float(bracket["requested_t_s"]), "observation time mismatch")
        detail = {}
        for i, (component, state_name) in enumerate(zip(
                ("x", "eta_e", "eta_mu"), ("x_state", "eta_e_MeV", "eta_mu_MeV"))):
            o1, o2 = (float(q[component + suffix]) for suffix in ("_O1", "_O2"))
            mo = max(abs(o1), abs(o2))
            mi = max(mo, abs(float(main[left-1][state_name])),
                     abs(float(main[right-1][state_name])))
            d = abs(o2-o1)
            d1 = fma(1e-12, mo, (1e-17, 1e-23, 1e-23)[i])
            fo = max(fma(1e-13, mo, (1e-18, 1e-24, 1e-24)[i]), 64*math.ulp(mo))
            fi = max(fma(1e-11, mi, (1e-16, 1e-22, 1e-22)[i]), 64*math.ulp(mi))
            u = 2*max(d, fo)
            for prefix, value in (("d_", d), ("D_O1_", d1), ("F_", fi), ("U_", u)):
                need(float(q[prefix+component]) == value, "recorded budget differs: " + prefix+component)
            detail[component] = {"d_O": d, "D_O1": d1, "U_O": u, "F_i": fi,
                                 "d_over_D_O1": d/d1, "U_over_point2F": u/(.20*fi),
                                 "pass": d <= d1 and u <= .20*fi}
        passed = all(v["pass"] for v in detail.values())
        need(passed == (q["self_qualified"] == "1"), "qualification status differs")
        if not passed:
            need(q["status"] == "NUMERICALLY_UNRESOLVED" and
                 q["diagnostic_self_qualified"] == "0", "failed state published as qualified")
            failures.append({"index": index, "time_s": t, "time_year": t/(365.25*86400),
                             "category": q["category"], "knots": int(q["knots"]),
                             "left_accepted_ordinal": left, "right_accepted_ordinal": right,
                             "left_time_s": float(bracket["t_previous_s"]),
                             "right_time_s": float(bracket["t_new_s"]), "components": detail})
        else:
            need(q["status"] == "QUALIFIED" and q["diagnostic_self_qualified"] == "1",
                 "missing prior diagnostic qualification")
    need(len(failures) == 1 and failures[0]["index"] == len(checkpoints)-1,
         "execution continued after failure")
    partial = {}
    for mode in ("source", "control"):
        directory = root/(mode+"-baseline")
        table = rows(directory, "trajectory.tsv")
        complete = [r for r in table if all(v is not None for v in r.values())]
        need(table[:len(complete)] == complete, "incomplete nonterminal row")
        need(mode == "control" or len(complete) == len(checkpoints)-1,
             "unresolved source state serialized as diagnostic output")
        identity = {}
        for name in ("R18_residual_erg_s", "Ra_Rb_residual_erg_s", "Rb_Rc_residual_erg_s"):
            ratios = [abs(float(r[name]))/(64*sys.float_info.epsilon*max(
                1., abs(float(r["P_dir_actual_erg_s"])), abs(float(r["P_dir_eq_erg_s"]))))
                for r in complete]
            identity[name] = {"maximum_utilization": max(ratios), "partial_rows_pass": max(ratios) <= 1}
        main_audit = rows(directory, "main.audit.tsv")[0]
        need(main_audit["unique_positive_t1_targets"] == "1" and
             main_audit["intermediate_observation_t1_matches"] == "0", "nonpassive main")
        b0 = float(complete[0]["B_count"])
        bdot = float(complete[0]["Bdot_count_s"])
        depletion = abs((fma(bdot, float(main_audit["t1_s"]), b0)-b0)/b0)
        partial[mode] = {
            "complete_retained_diagnostic_rows": len(complete),
            "truncated_terminal_rows": len(table)-len(complete),
            "last_complete_time_s": float(complete[-1]["t_s"]),
            "partial_identity_checks": identity,
            "partial_max_abs_DeltaB_over_B0": max(abs(float(r["DeltaB_over_B0"])) for r in complete),
            "prescribed_full_horizon_max_abs_DeltaB_over_B0": depletion,
            "partial_max_frozen_utilization": max(float(r["max_frozen_utilization"]) for r in complete),
            "partial_all_valid": all(r["valid_through_sample"] == "1" for r in complete),
            "main_one_terminal_ceiling": True,
            "full_campaign_qualified": False,
        }
    return {"classification": "FAILED NUMERICAL EVIDENCE; NOT CANDIDATE",
            "budget_arithmetic": "explicit binary64 fma matching qualified AppleClang a+r*M",
            "read_only_arithmetic_audit_pass": True, "new_ODE_runs": 0,
            "source_checkpoint_rows": len(checkpoints), "failure": failures[0],
            "partial_diagnostics": partial, "final_state_published": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.campaign_root)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print("PASS evidence arithmetic audit; campaign remains NUMERICALLY_UNRESOLVED")
