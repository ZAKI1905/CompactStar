#!/usr/bin/env python3
"""Collect job 19171's startup refusal. Reads evidence; never launches a run."""
import argparse
import hashlib
import json
import pathlib
import re
import shutil


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=pathlib.Path)
    parser.add_argument("evidence", type=pathlib.Path)
    args = parser.parse_args()
    args.evidence.mkdir(parents=True, exist_ok=True)
    log = (args.run / "stage.log").read_text()
    refused = "/usr/lib64/ld-linux-x86-64.so.2"
    assert "ValueError: missing/mismatched authenticated dependency: " + refused in log
    assert not (args.run / "run").exists(), "Unexpected numerical output; inspect manually"
    usage = json.loads(log.split("RESOURCE_USAGE_JSON ", 1)[1])
    assert usage["exit_code"] == 1
    audits = {}
    for name in ("login", "compute"):
        audit = json.loads((args.run / (name + "-runtime-audit.json")).read_text())
        assert audit["checked_files"] == 508
        assert all(c["returncode"] == 0 for c in audit["commands"].values())
        audit.pop("files")
        # Address-space randomization is irrelevant to library path identity.
        audit["resolved_interpreter_paths"] = re.findall(
            r"(/[^\s]+)", audit["commands"]["interpreter_list"]["stdout"])
        audits[name] = audit
    assert audits["login"]["mismatches"] == []
    assert [r["path"] for r in audits["compute"]["mismatches"]] == [refused]
    assert audits["login"]["snapshot_sha256"] == audits["compute"]["snapshot_sha256"]
    assert audits["login"]["resolved_interpreter_paths"] == audits["compute"]["resolved_interpreter_paths"]
    assert refused not in audits["compute"]["resolved_interpreter_paths"]
    assessment = dict(
        status="BLOCKED_BEFORE_EVOLUTION", stage_A="STARTUP_AUTHENTICATION_REFUSED",
        stage_B="NOT_SUBMITTED_STAGE_A_DID_NOT_PASS", thermal_job=19171,
        read_only_runtime_audit_job=19172, numerical_runner_invoked=False,
        completed_trajectories=0, trajectory_samples=0, figures=[],
        scientific_gates={gate: "NOT_EVALUATED_NO_TRAJECTORY" for gate in (
            "state_accuracy", "reference_convergence", "accumulated_roundoff", "R18", "R20",
            "restart_agreement", "checkpoint_reconstruction", "actual_breakpoint_coverage",
            "source_control_consistency", "source_minus_control_signal")},
        resource_usage=usage, runtime_audits=audits,
        interpretation="The manifest additionally pins a login-host loader observed through ldd; "
                       "the ELF interpreter and its resolved shared libraries are authenticated shared paths. "
                       "No manifest change, bypass, or thermal retry was performed.")
    (args.evidence / "assessment.json").write_text(json.dumps(assessment, indent=2) + "\n")
    for name in ("stage.log", "allocation.txt", "slurm-19171.log",
                 "runtime-audit-allocation.txt", "runtime-audit-19172.log"):
        shutil.copyfile(args.run / name, args.evidence / name)
    files = sorted(p for p in args.run.rglob("*") if p.is_file())
    (args.evidence / "external-SHA256SUMS").write_text(
        "".join(digest(p) + "  " + str(p) + "\n" for p in files))
    print("Confirmed startup refusal: 0 trajectories; compute node matches 507/508 dependencies.")


if __name__ == "__main__":
    main()
