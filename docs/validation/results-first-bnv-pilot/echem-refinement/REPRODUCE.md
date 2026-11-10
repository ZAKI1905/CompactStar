# Reading and reproducing the handoff

Read `REPORT.md`, `PREDECLARATION.md`, `evidence/summary.json` and the full
`evidence/assessment.json` first. All statuses remain UNQUALIFIED. The report's
containing Git commit is the handoff revision; the executed binary's source
revision is recorded independently in `evidence/build-manifest.json`.

From the repository root, verify the committed bundle:

```sh
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/SHA256SUMS
```

On EKU, verify retained external evidence without executing any physics:

```sh
sha256sum -c docs/validation/results-first-bnv-pilot/echem-refinement/evidence/external-SHA256SUMS
```

The following optional commands recompute **offline assessments of existing
records**, not trajectories. They write fresh comparison files in /tmp; they do
not overwrite evidence. They are supplied for Fable, not executed again during
handoff. Python standard library only.

```sh
study=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/echem-refinement-20261109-r1
inputs=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery/run1/raw
gsl_source=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/source/gsl-2.7.1/ode-initval2
python3 -B tools/linux/native_qualification/assess_refinement.py \
  "$study/study/output" "$inputs/coefficients.tsv" /tmp/echem-assessment-review.json
cmp /tmp/echem-assessment-review.json \
  docs/validation/results-first-bnv-pilot/echem-refinement/evidence/assessment.json
python3 -B tools/linux/native_qualification/audit_chemical_roundoff.py \
  "$study/study/output" "$inputs/coefficients.tsv" "$inputs/tangent.tsv" \
  "$gsl_source" /tmp/echem-conditional-roundoff-review.json
cmp /tmp/echem-conditional-roundoff-review.json \
  docs/validation/results-first-bnv-pilot/echem-refinement/evidence/conditional-roundoff.json
```

`assess_refinement.py` computes configured gate results, observed local-witness
envelopes and 100-digit Echem algebra. Its `complete` flag means the declared
matrix is present, not scientifically qualified. Its Stage A gate is explicitly
false because complete accumulated-roundoff evidence is missing.

`audit_chemical_roundoff.py` is exploratory, unverified analysis conditional on
recorded thermal stages. Review its tableau extraction, context classification,
step reconstruction and Jacobian enclosure. Its JSON's excluded terms and false
Stage A gate are essential. Do not use its conditional numbers as a certified
coupled or original-t0 error bound.

No sbatch command is provided: handoff does not initiate another experiment.
The committed wrappers and build/invocation manifests preserve how the completed
study was run. Production settings remain unchanged.
