#!/usr/bin/env python3
"""Build current-profile Phase-5C input certificates for the fixed BA13 grid.

This is a validation-tool adapter only.  It reuses the accepted Phase-5C
certificate construction and immutable Phase-5C goals; it does not fit or
enlarge a frozen-background threshold.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True


def thermal_source(profile: Path):
    """Reproduce the governed Phase-5D free-gas entropy adapter on this star."""
    output = profile / "thermal"
    output.mkdir()
    data = np.loadtxt(profile / "profile.tsv", skiprows=1)
    meta = (profile / "model.txt").read_text().splitlines()
    masses = np.array([float(line.split()[0]) for line in meta[:4]])
    hc = float(meta[0].split()[1])
    data = data[data[:, 3] > 0]
    data = data[np.argsort(data[:, 3])]
    nb, index = np.unique(data[:, 3], return_index=True)
    density = data[index, 6:10]
    pf = hc * np.cbrt(3 * np.pi**2 * density)
    mu = np.hypot(pf, masses)
    slope = np.sum(pf * mu, axis=1) / (3 * hc**3 * nb)
    if not np.all(np.isfinite(slope) & (slope > 0)):
        raise RuntimeError("invalid current-star entropy adapter")
    temperatures = [0, 1, 2, 4]
    charges = [0, 1]
    for name, values in (("t", temperatures), ("nb", nb), ("yq", charges)):
        (output / ("eos." + name)).write_text(
            "1 " + str(len(values)) + "\n"
            + "\n".join(format(value, ".17g") for value in values) + "\n"
        )
    with (output / "eos.thermo").open("w") as stream:
        stream.write(f"{masses[0]:.17g} {masses[1]:.17g} 0\n")
        for it, temperature in enumerate(temperatures):
            for ib, value in enumerate(slope):
                for iy in range(len(charges)):
                    stream.write(
                        f"{it + 1} {ib + 1} {iy + 1} 0 "
                        f"{value * temperature:.17g} 0 0 0 0 0 0\n"
                    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--phase5-characterization", type=Path, required=True)
    args = parser.parse_args()
    source = args.source_root.resolve()
    sys.path.insert(0, str(source / "tests/analysis"))
    from chemical_production_evidence import certificates, transport

    fixed = json.loads(
        (source / "docs/validation/phase5c2_preproduction_evidence.json").read_text()
    )
    characterized = json.loads(args.phase5_characterization.read_text())
    profiles = args.profiles.resolve()
    directories = sorted(
        profiles.glob("target-*"), key=lambda path: int(path.name.split("-")[-1])
    )
    if len(directories) != 21:
        raise RuntimeError("frozen certificate requires exactly 21 current profiles")
    for index, directory in enumerate(directories):
        if directory.name != f"target-{index}":
            raise RuntimeError("noncanonical frozen-profile order")
        certificate, _ = certificates(directory, characterized, fixed["goals"])
        transport(certificate, directory / "certificate.txt", fixed)
        (directory / "certificate.json").write_text(
            json.dumps(certificate, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        thermal_source(directory)
        print(f"FROZEN_CHEMICAL_CERTIFICATE {index} PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
