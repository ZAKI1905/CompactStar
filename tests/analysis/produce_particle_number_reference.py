"""Canonical candidate-only Phase-5B structural reference producer.

Run only after PB1-PB14 and the mutation campaign pass. Every run solves all
stars in fresh output. No historical scratch coefficient or baseline is read.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--executable", required=True, type=Path)
parser.add_argument("--source-root", required=True, type=Path)
parser.add_argument("--output-dir", required=True, type=Path)
args = parser.parse_args()
root = args.source_root.resolve()
output = args.output_dir.resolve()
output.mkdir(parents=True, exist_ok=False)
data = output / "solve"
result = subprocess.run([str(args.executable.resolve()), "emit", str(data)],
                        capture_output=True, text=True)
(output / "producer.log").write_text(result.stdout + result.stderr)
if result.returncode:
    raise SystemExit(f"Producer refused; see {output / 'producer.log'}")

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rows(name, text_columns):
    with (data / name).open() as stream:
        return [{key: value if key in text_columns else float(value)
                 for key, value in row.items()}
                for row in csv.DictReader(stream, delimiter="\t")]

source_files = [
    "CompactStar/Analysis/ParticleNumberResponse.hpp",
    "CompactStar/Analysis/src/ParticleNumberResponse.cpp",
    "CompactStar/Core/NStar.hpp", "CompactStar/Core/src/NStar.cpp",
    "CompactStar/Core/StarProfile.hpp", "CompactStar/Core/RotationSolver.hpp",
    "CompactStar/Core/src/RotationSolver.cpp", "CompactStar/Core/src/TOVSolver.cpp",
    "CompactStar/RelativityUnits.hpp", "CompactStar/Geometry.hpp",
    "CompactStar/AngularVelocity.hpp", "CompactStar/Units.hpp",
    "CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp",
    "tests/eos/structure1/table.hpp", "tests/eos/structure1/local_oracle.hpp",
    "tests/eos/structure1/tov_oracle.hpp",
    "tests/analysis/phase5b_freegas_validation.cpp",
    "tests/analysis/particle_number_homogeneous.hpp",
    "tests/analysis/produce_particle_number_reference.py",
]
summary = dict(line.split("\t") for line in (data / "summary.tsv").read_text().splitlines())
artifact = {
    "schema": "ADR-0011-structural-candidate-v1",
    "status": "CANDIDATE ONLY; independent review and human ratification pending; not a governed baseline",
    "starting_canonical_sha": "a43d02227bf53c3242d3212f81dd71963804f3aa",
    "fixture": {"identity": "human-ratified Structure-1 midpoint; FR2005 common static state",
                "rho_c_g_cm3": 1.10e15, "EOS_intervals": 8192, "radial_resolution": 80000,
                "source_qualified_Mmax": False, "published_I_Omega_target": False},
    "domain": "WholeStar P=0, canonical finite cut with bounded remainder; no core cutoff",
    "surface": {"policy": "positive-source monotone pe comparison inequalities",
                "revision": "PB13-comparison-v1", "correction": "zero estimate plus explicit validated bound"},
    "units": {"N": "count", "A": "count km^2", "B": "count km^2, derivative with respect to epsilon_c in km^-2",
              "K": "count km^2", "I_geom": "K, conditional whole-star mapping", "I_phys": "count s^2",
              "q": "Omega_geom^2 in km^-2", "Omega_geom": "Omega_phys / c",
              "c_km_s": 299792.458, "equilibrium_rate": "2 I_phys Omega_phys Omega_dot_phys [count/s]"},
    "coefficients": rows("coefficients.tsv", {"species", "label"}),
    "summary": {key: float(value) for key, value in summary.items()},
    "numerical_errors": "reported quadrature/roundoff, central-step and explicit tail budgets; radial/EOS and independent-oracle convergence are separately recorded in the validation evidence",
    "sequence": {"coordinate": "ln(rho_c/rho_ref)", "steps": [.001, .0005, .00025],
                 "method": "five-point at achieved abscissae, compared with three-point",
                 "branch_rho_g_cm3": [1.095e15, 1.105e15], "radial_resolution": 80000,
                 "stencils": rows("stencil.tsv", {"species"})},
    "provenance": {
        "profile_identity": "runtime owning profile address plus ProfileVersion; addresses are not persistent or serialized",
        "roles": "central and sequence ordinals identify this execution's contributors only; they are not persistent object IDs",
        "contributors": rows("sources.tsv", {"role"}),
        "Hartle_inputs": "each contributing star's own first-order and monopole response, matching profile identity/version",
        "EOS_identity": "FR2005 equilibrium free gas, vacuum/pe/npe/npe-mu below Sigma ceiling",
        "EOS_revision": "canonical-A1-structure1-table-v1",
        "EOS_table_sha256": digest(data / "freegas.tsv"),
        "EOS_table_authority": (data / "freegas.tsv.provenance.txt").read_text(),
        "source_sha256": {name: digest(root / name) for name in source_files},
        "build": {"configuration": "Debug", "architecture": platform.machine(),
                  "compiler": subprocess.check_output(["/usr/bin/clang++", "--version"], text=True).splitlines()[0]},
    },
    "INV09": "INTENDED BUT UNVERIFIED", "INV11": "UNRESOLVED",
    "non_goals": ["Btilde", "chemical Z/W", "eta evolution", "weak-rate evolution", "BNV"],
}
target = output / "structural-response.json"
target.write_text(json.dumps(artifact, indent=2, sort_keys=True, allow_nan=False) + "\n")
print(digest(target), target)
