# Reproduce the D1 resolution evidence

The authorized calculation has completed once in job 19185. These commands are
reproducibility instructions, not authorization for another retained-data run.
All files in the original `../d1/` directory and original verifier remain intact.
See [PREDECLARATION.md](PREDECLARATION.md#L1) for method, criteria and assumptions.

Authenticate the manifests from the repository root:

```sh
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/d1/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/d1-resolution/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/evidence/external-SHA256SUMS
```

The reviewed source is ad4e91a; the preserved original D1 receipt is 871dbd7.
The predeclared resolution implementation is commit 7842938. The full execution
commit and source hashes are recorded in `evidence/job-19185/diagnostic/results.json`.
Authentication discovers the study from the original committed inventory,
checks all 354 external files afresh, checks the physical tables through their
committed authentication record, and verifies the original D1 bytes against
871dbd7. New scripts and predeclaration must match committed HEAD before a run.

The actual launch used the following environment and wrapper (the output root
already exists and must never be reused/overwritten):

```sh
export D1_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export D1_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/bnv-thermo-01-d1-resolution-20261109-r1
export D1_PYTHON=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python
sbatch --export=ALL --output="$D1_OUTPUT/slurm-%j.out" tools/linux/native_qualification/chemical_reference_resolution.sbatch
```

The wrapper runs all 17 analytic tests before the matrix. The interpreter has
no third-party numerical dependency. It requires Slurm and exactly one logical
CPU affinity and creates a new diagnostic directory exclusively. Resource
limits are one task/CPU, 4 GiB and 30 minutes. Time is local seconds from exact
binary64 t_left; eta and Echem use MeV; x is ln(Tinf/1e8 K).

K=8,16,32,64 is the number of steps per Hermite segment. JSON uses nominal
N=8K keys for the shared unchanged convergence gate, and explicitly records
actual total steps. The entire path and its tiny pieces remain unchanged.

Reporting can be reproduced without integrating anything:

```sh
python3 -B tools/linux/native_qualification/summarize_chemical_reference_resolution.py \
  docs/validation/results-first-bnv-pilot/echem-refinement/d1-resolution/evidence/job-19185/diagnostic/results.json \
  /tmp/compactstar-d1-resolution-report
```

Compare the generated `summary.json`, `TABLES.md`, `chemical-states.tsv`,
`chemical-energy.tsv`, and `original-differences.tsv` with the committed copies.
`original-failure-analysis.json` is produced by algebra on retained old outputs
and RHS evaluations, not by rerunning the original global-grid IVPs.

`evidence/capture.json` maps every copied raw file to its cluster location and
SHA-256. No raw result is hand-edited. `failure-details.json` is empty for this
PASS; the original D1 failure file remains in the original evidence directory.
No production execution is invoked by any of these scripts.
