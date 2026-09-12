#!/usr/bin/env python3
"""Execute one Phase-5D oracle or trajectory command from explicit fresh state."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np

from manifest import check_manifest, digest

QUALIFICATION_HASHES = {
    "freegas.tsv": "7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a",
    "profile.tsv": "e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890",
    "model.txt": "3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4",
    "certificate": "7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453",
}


def verify_hashes(expected):
    for path, value in expected.items():
        if not Path(path).is_file() or digest(path) != value:
            raise RuntimeError("generated source changed after construction: " + str(path))


def prepare(qualification, output):
    profiles = list(qualification.glob("chemical-characterization/run-*/t8192-r80000"))
    if len(profiles) != 1:
        raise RuntimeError("requires exactly one freshly generated radial80000 qualification")
    profile = profiles[0]
    certificate = qualification / "certificate-t8192-r80000.txt"
    for path, expected in [
        (profile / "freegas.tsv", QUALIFICATION_HASHES["freegas.tsv"]),
        (profile / "profile.tsv", QUALIFICATION_HASHES["profile.tsv"]),
        (profile / "model.txt", QUALIFICATION_HASHES["model.txt"]),
        (certificate, QUALIFICATION_HASHES["certificate"]),
    ]:
        if not path.is_file() or digest(path) != expected:
            raise RuntimeError("qualification input changed: " + str(path))
    thermal = output / "thermal"
    thermal.mkdir()
    data = np.loadtxt(profile / "profile.tsv", skiprows=1)
    meta = (profile / "model.txt").read_text().splitlines()
    masses = np.array([float(x.split()[0]) for x in meta[:4]])
    hc = float(meta[0].split()[1])
    data = data[data[:, 3] > 0]
    data = data[np.argsort(data[:, 3])]
    nb, index = np.unique(data[:, 3], return_index=True)
    density = data[index, 6:10]
    pf = hc * np.cbrt(3 * np.pi**2 * density)
    mu = np.hypot(pf, masses)
    slope = np.sum(pf * mu, axis=1) / (3 * hc**3 * nb)
    if not np.all(np.isfinite(slope) & (slope > 0)):
        raise RuntimeError("invalid entropy adapter")
    temperatures = [0, 1, 2, 4]
    charges = [0, 1]
    for name, values in [("t", temperatures), ("nb", nb), ("yq", charges)]:
        (thermal / ("eos." + name)).write_text(
            "1 " + str(len(values)) + "\n"
            + "\n".join(format(x, ".17g") for x in values) + "\n"
        )
    with (thermal / "eos.thermo").open("w") as stream:
        stream.write(f"{masses[0]:.17g} {masses[1]:.17g} 0\n")
        for it, temperature in enumerate(temperatures):
            for ib, value in enumerate(slope):
                for iy in range(len(charges)):
                    stream.write(
                        f"{it + 1} {ib + 1} {iy + 1} 0 "
                        f"{value * temperature:.17g} 0 0 0 0 0 0\n"
                    )
    return profile, certificate, thermal


def run_regressions(source_root, build_dir, output, suffix):
    codes = {}
    for name in [
        "phase5b_structural_response_regression",
        "phase5c_chemical_coefficient_regression",
    ]:
        log_path = output / f"{name}-{suffix}.log"
        command = [
            "ctest", "--test-dir", str(build_dir), "-R", "^" + name + "$",
            "--output-on-failure", "--no-tests=error",
        ]
        with log_path.open("w") as log:
            result = subprocess.run(
                command, cwd=source_root, stdout=log, stderr=subprocess.STDOUT
            )
        codes[name] = result.returncode
        (output / f"{name}-{suffix}.rc").write_text(str(result.returncode) + "\n")
        if result.returncode:
            raise RuntimeError(f"STOP BEFORE TRAJECTORY: {name} rc={result.returncode}")
    return codes


def pretrajectory_gate(source_root, build_dir, output, entry):
    check_manifest(source_root, entry)
    codes = run_regressions(source_root, build_dir, output, "pretrajectory")
    check_manifest(source_root, entry)
    result = {
        "protected_count": 33,
        "all_equal_entry": True,
        "regression_rc": codes,
        "entry_manifest_sha256": digest(entry),
    }
    (output / "pretrajectory-gate.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("mode", choices=["oracles", "trajectory"])
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--entry-manifest", type=Path, required=True)
    parser.add_argument("--qualification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--oracle-executable-sha256")
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    build_dir = args.build_dir.resolve()
    entry = args.entry_manifest.resolve()
    executable = args.executable.resolve()
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError("qualified-suite output must be fresh")
    output.mkdir(parents=True)
    check_manifest(source_root, entry)
    if args.mode == "trajectory":
        if not args.oracle_executable_sha256:
            raise RuntimeError("same-invocation oracle executable identity required")
        if args.oracle_executable_sha256 != digest(executable):
            raise RuntimeError("oracle executable differs from trajectory executable")
    elif args.oracle_executable_sha256:
        raise RuntimeError("oracle executable identity is a trajectory-only input")

    entry_codes = {}
    profile, certificate, thermal = prepare(args.qualification.resolve(), output)
    executed = output / "phase5d_evolution.executed"
    shutil.copy2(executable, executed)
    command = [
        str(executed), str(profile), str(certificate), str(thermal), str(output),
        args.mode, str(entry),
    ]
    inputs = [
        profile / "freegas.tsv", profile / "profile.tsv", profile / "model.txt",
        certificate, *sorted(thermal.glob("eos.*")),
    ]
    provenance = {
        "command": command,
        "source_root": str(source_root),
        "build_dir": str(build_dir),
        "output_root": str(output),
        "executable_sha256": digest(executable),
        "executed_image": str(executed),
        "executed_sha256": digest(executed),
        "inputs": {str(path): digest(path) for path in inputs},
        "cmake_cache_sha256": digest(build_dir / "CMakeCache.txt"),
        "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=source_root, text=True
        ).strip(),
        "entry_regression_rc": entry_codes,
    }
    (output / "run-provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    verify_hashes({path: digest(path) for path in inputs})
    check_manifest(source_root, entry)
    print("QUALIFIED COMMAND", json.dumps(command), flush=True)
    gate_error = None
    with (output / "console.log").open("w") as log:
        process = subprocess.Popen(
            command, cwd=source_root, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        for line in process.stdout:
            log.write(line)
            log.flush()
            if line.strip() == "PRE_TRAJECTORY_READY":
                try:
                    pretrajectory_gate(source_root, build_dir, output, entry)
                    process.stdin.write("PROTECTED_GOVERNED_GATES_PASS\n")
                    process.stdin.flush()
                except Exception as error:
                    gate_error = str(error)
                    process.stdin.close()
        returncode = process.wait()
    if gate_error:
        (output / "gate-error.txt").write_text(gate_error + "\n")
        returncode = returncode or 1
    (output / "raw-rc.txt").write_text(str(returncode) + "\n")
    print((output / "console.log").read_text(), flush=True)
    print("RAW EXECUTABLE RC", returncode, "EVIDENCE", output, flush=True)
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
