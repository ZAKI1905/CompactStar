# D1 evidence reproduction

The final result is **C1 FAIL at points 1, 64 and 1440**, from the sole complete
retained-interval run, job 19184. No modified-parameter retry is authorized.
The earlier missing-review stop is historical and remains preserved.

From the repository root, verify the committed package and all original and
new external artifacts without executing a production solver:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/d1/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/evidence/external-SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/d1/evidence/external-SHA256SUMS
```

The complete result is `evidence/job-19184/results.json`. It includes the
110-digit references, all 156 RK4 endpoints, uncertainty components, comparisons
and input authentication chain. `evidence/summary.json` contains every failed
C1 comparison; `evidence/summary.tsv` provides the point readout; the current
`evidence/point-*.json` files are exact copies of job 19184's per-point outputs.
Job 19183's full receipt retains all 13 earlier NOT_RUN dispositions.

The original owner-supplied review is preserved byte-for-byte in
`review/BNV-THERMO-01_FABLE_ADJUDICATION.md`, with attachment provenance in
`review/provenance.json`. Its SHA-256 is
`46201727575270055d237178b2ad0433ddc22837c0cc1fa646a8d34eddcb0066`.
It was committed at `06066be` before numerical execution. No separately signed
reviewer digest is claimed.

The maintained wrapper is `tools/linux/native_qualification/chemical_reference.sbatch`.
Its fixed limits are one logical execution CPU, 4 GiB and 30 minutes. It runs the
analytic tests before invoking the diagnostic. Job 19184 used these values:

| Variable | Executed value |
| --- | --- |
| `D1_SOURCE` | `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot` |
| `D1_OUTPUT` | `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/bnv-thermo-01-d1-20261109-r3` |
| `D1_PYTHON` | `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python` |
| `D1_ADJUDICATION` | `D1_SOURCE` plus `/docs/validation/results-first-bnv-pilot/echem-refinement/d1/review/BNV-THERMO-01_FABLE_ADJUDICATION.md` |
| `D1_ADJUDICATION_SHA256` | `46201727575270055d237178b2ad0433ddc22837c0cc1fa646a8d34eddcb0066` |

The wrapper was submitted with `sbatch --parsable`, with its stdout path set to
`D1_OUTPUT/slurm-%j.log`. It invokes the maintained Python script with
`--output D1_OUTPUT/diagnostic --adjudication D1_ADJUDICATION
--adjudication-sha256 D1_ADJUDICATION_SHA256`. Retained runs must not be
reused as output roots. The same script rejects an existing output directory.

The `-r1` directory preserves the first launch failure and exact executed
scripts; `-r2` preserves the missing-review authentication stop; `-r3` preserves
the final calculation and exact executed scripts. All are under `/mnt/sdd` and
included in the new external manifest. The original 354-file study inventory
remains separate and unchanged. The script has no option to skip authentication,
relax fidelity, change the thermal band, choose another point list or alter the
RK4 matrix. No production solver is invoked.
