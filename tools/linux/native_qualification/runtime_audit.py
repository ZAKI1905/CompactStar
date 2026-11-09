#!/usr/bin/env python3
"""Read-only inventory after a native snapshot refusal; never runs the solver."""
import argparse
import hashlib
import json
import pathlib
import platform
import re
import subprocess


def digest(path):
    value = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshots", type=pathlib.Path)
    args = parser.parse_args()
    snapshots = [json.loads((args.snapshots / (role + ".json")).read_text())
                 for role in ("source", "control")]
    assert snapshots[0]["files"] == snapshots[1]["files"]
    rows = []
    for item in snapshots[0]["files"]:
        row = dict(item)
        try:
            row["actual_sha256"] = digest(item["path"])
            row["match"] = row["actual_sha256"] == item["sha256"]
        except OSError as error:
            row.update(match=False, error=str(error))
        rows.append(row)
    binary = snapshots[0]["stage_binary"]
    commands = {}
    for command in (["readelf", "-l", binary], ["ldd", binary]):
        result = subprocess.run(command, text=True, capture_output=True)
        commands[command[0]] = dict(command=command, returncode=result.returncode,
                                   stdout=result.stdout, stderr=result.stderr)
    match = re.search(r"Requesting program interpreter: (.*?)\]", commands["readelf"]["stdout"])
    if match:
        command = [match.group(1), "--list", binary]
        result = subprocess.run(command, text=True, capture_output=True)
        commands["interpreter_list"] = dict(command=command, returncode=result.returncode,
                                            stdout=result.stdout, stderr=result.stderr)
    print(json.dumps(dict(host=platform.node(), os_release=pathlib.Path("/etc/os-release").read_text(),
                          snapshot_sha256={role: digest(args.snapshots / (role + ".json"))
                                           for role in ("source", "control")},
                          checked_files=len(rows), mismatches=[r for r in rows if not r["match"]],
                          commands=commands, files=rows), indent=2))


if __name__ == "__main__":
    main()
