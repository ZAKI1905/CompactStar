#!/usr/bin/env python3
"""Run the predeclared clean P2 hierarchy after authenticated prerequisites.

Independent source/control processes run together at each tier. A failed process
stops its peer and prevents all later tiers. No retry or candidate production.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from adr0017_production_verify import EXPECTED, sha256, tree_sha256


def run(args):
    evidence = args.prerequisites.resolve()
    authentication = json.loads((evidence / "input-authentication.json").read_text())
    inputs = {k: Path(v) for k, v in authentication["paths"].items()}
    if not authentication["pass"]:
        raise RuntimeError("input authentication failed")
    for key, path in inputs.items():
        if tree_sha256(path) != EXPECTED[key]:
            raise RuntimeError("input bytes changed: " + key)
    prerequisites = json.loads((evidence / "prerequisites.json").read_text())
    required = {"source_projection", "matched_control", "frozen_monitor",
                "direct_energy", "thermal_ledger", "passive_contract"}
    if set(prerequisites) != required or any(v["exit_code"] for v in prerequisites.values()):
        raise RuntimeError("prerequisite execution incomplete or failed")
    regression = json.loads(args.regression.read_text())
    baseline = EXPECTED_BASELINE = "2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd"
    if (regression["producer_raw_rc"] != 0 or regression["generated_sha256"] != baseline
            or regression["baseline_sha256"] != baseline
            or any(v != "PASS" for v in regression["controls"].values())):
        raise RuntimeError("fresh governed regression did not pass")
    generated = Path(regression["scratch_root"]) / "governed-artifact.json"
    if sha256(generated) != EXPECTED_BASELINE:
        raise RuntimeError("fresh governed artifact changed")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    flag = root / "prerequisites.flag"
    flag.write_text("BNV_RESUME_PREREQUISITES_PASS\n")
    binary = args.executable.resolve()
    journal = {"executable": str(binary), "executable_sha256": sha256(binary),
               "input_hashes": authentication["hashes"],
               "prerequisites_sha256": sha256(evidence / "prerequisites.json"),
               "regression_sha256": sha256(args.regression), "runs": []}

    def save():
        (root / "execution.json").write_text(json.dumps(journal, indent=2) + "\n")

    save()
    environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    for tier in ("baseline", "refined", "ultra"):
        active = []
        for mode in ("source", "control"):
            name = mode + "-" + tier
            cmd = [binary, inputs["profile_tree"], inputs["certificate"], inputs["thermal"],
                   root / (name + "-context"), inputs["entry"], inputs["frozen"],
                   inputs["coefficients"], root / name, mode, tier, flag]
            record = {"name": name, "command": list(map(str, cmd)), "started_unix_s": time.time()}
            journal["runs"].append(record)
            log = (root / (name + ".log")).open("w")
            process = subprocess.Popen(cmd, stdout=log, stderr=log, env=environment)
            active.append((process, record, log, time.monotonic()))
        failed = False
        try:
            while active:
                for item in active[:]:
                    process, record, log, start = item
                    rc = process.poll()
                    if rc is None:
                        continue
                    record.update(exit_code=rc, wall_s=time.monotonic()-start)
                    log.close()
                    active.remove(item)
                    failed |= rc != 0
                    print(record["name"], "exit", rc, "wall_s", record["wall_s"], flush=True)
                    save()
                if failed:
                    for process, record, log, start in active:
                        process.terminate()
                        record.update(exit_code=process.wait(), wall_s=time.monotonic()-start,
                                      cancelled_after_peer_failure=True)
                        log.close()
                    active.clear()
                    save()
                    raise RuntimeError("campaign gate failed; later tiers not executed")
                if active:
                    time.sleep(.25)
        finally:
            for process, record, log, start in active:
                process.terminate()
                record.update(exit_code=process.wait(), wall_s=time.monotonic()-start,
                              interrupted=True)
                log.close()
            save()
    result = subprocess.run([sys.executable, "-B", str(Path(__file__).with_name("controlled_campaign_verify.py")),
                             "--campaign-root", str(root), "--output", str(root / "comparison.json")])
    journal["comparison_exit_code"] = result.returncode
    save()
    return result.returncode


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--prerequisites", type=Path, required=True)
    parser.add_argument("--regression", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    try:
        raise SystemExit(run(parser.parse_args()))
    except Exception as error:
        print("STOP", error, file=sys.stderr)
        raise SystemExit(1)
