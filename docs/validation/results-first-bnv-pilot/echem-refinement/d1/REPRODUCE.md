# D1 evidence reproduction

The retained result is an authentication failure, not a numerical reference
result. The original review must be supplied and authenticated before any
retained-interval calculations. Do not invent a replacement or hash a newly
written reconstruction of the review and treat it as independent provenance.

From the repository root, verify the package and external inputs without a
production solver:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/d1/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/evidence/external-SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/d1/evidence/external-SHA256SUMS
```

The complete stopped result is `evidence/job-19183/results.json`; individual
`evidence/point-*.json` records explicitly identify calculations not run. The
first launch failure and its exact executed scripts remain under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/bnv-thermo-01-d1-20261109-r1`.
The authenticated-input check is in the corresponding `-r2` directory. These
directories are preserved and must not be reused as new output roots.

The maintained wrapper is `tools/linux/native_qualification/chemical_reference.sbatch`.
Its fixed Slurm limits are one logical execution CPU, 4 GiB and 30 minutes. It
runs the analytic tests before invoking the diagnostic. Inputs are configured
through the following environment variables:

| Variable | Meaning |
| --- | --- |
| `D1_SOURCE` | Actual repository worktree descended from the reviewed commit |
| `D1_OUTPUT` | Fresh cluster receipt directory under `/mnt/sdd`; its `diagnostic` subdirectory must not exist |
| `D1_PYTHON` | `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python` |
| `D1_ADJUDICATION` | Actual supplied adjudication file, not a paraphrase |
| `D1_ADJUDICATION_SHA256` | Independently supplied/authenticated review SHA-256 |

Jobs 19182 and 19183 did not set the last two variables because those inputs
were unavailable. Authentication therefore could not pass. The wrapper's PATH
now includes the retained toolchain Git launcher, which repairs job 19182's
environment defect. Numerical definitions, tolerance, point set and RK4 matrix
were identical across the attempts.

The diagnostic has no option to skip authentication, relax the fidelity limit,
change the thermal band, choose another point list or alter the RK4 matrix.
When a proper review artifact becomes available, inspect Section D against the
implementation before executing; this package does not claim that comparison
has happened. Preserve the supplied review and its provenance in Git before
any further numerical work.
