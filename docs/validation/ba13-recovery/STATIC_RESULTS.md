# Native BA13 static producer — candidate, not promoted

Status: **STATIC PRODUCER PASS; SCIENTIFIC/PLATFORM QUALIFICATION UNQUALIFIED**.
The original final fitter and raw samples are recovered. No replacement fitter was
selected or executed. The remaining blocker is numerical-error qualification,
not missing historical evidence. See proposed ADR-0021.

## Authentication and original execution

The transferred ZIP is 116506018 bytes, SHA-256
`990d3708a24c71e518d4bb59e6ad3c23ae47f337ba924007b083ca598e6b8e2b`.
A fresh cluster extraction passed ZIP CRC, safe/unique paths, all 351 checksums
and complete manifest coverage. All 21 original target/tangent records, 63 final
coefficient/thermal/profile shards, assembled tables, 21 profile/certificate sets,
and final JSON/TSV certificates are present. Original evidence was not modified.

The original stages are recovered C++ `stars`, recovered Python certificates,
recovered C++ `coefficient(s)`, then original Python `compare_candidate.py
merge-frozen` and `frozen`. Exact source identities are in
`recovered-identities.json`; successful historical commands are in
`archive-recovery/historical_commands.txt`. C++ retains the accepted ADR-0016
ownership relocation and otherwise compiles unchanged against current Linux owners.

The original archived final JSON/TSV hashes are respectively
`6643054d8cba36035af15a15615baf9920d9a56adf3f5211f13ea9f2745f4409` and
`9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab`.
The owner verified byte-identical fitter reproduction on the Mac. On this Linux
Python 3.11.9 runtime, replaying the unchanged fitter on those same Mac raw tables
passes but is not byte-identical: JSON `5ea4f4ed2a012ddb8f79bcaf70a27f7f088e2628e34c332efe372c25a7ba30bf`,
TSV `639cc4ce3bb169be6dcbc257470d3d78aca64e3ecc0aba83dbae1e1af6b7e313`.
Small arithmetic differences occur in fit coefficients/covariances and residuals;
a unique Python/libm/platform cause has not been isolated. No historical output
identity is claimed for the Linux replay. The independent refit validates the
archived Mac certificate within its separately stated arithmetic allowance.

## Native results

Two complete Linux Release runs regenerated 21 achieved target-B stars,
21 Phase-5C certificates/Phase-5D adapters, 21 coefficient rows, 273 thermal rows
and 43029 profile-metric rows. All five raw tables, all 231 profile/adapter files,
and both final certificate files are byte-identical between native runs. Radial
80000, EOS parameter 8192, the 21 target fractions, 13 thermal knots, active species,
source support and fixed fit indices 10–20 remain unchanged. Abscissae use achieved B.

| Check | Actual result |
| --- | --- |
| Maximum target residual / B0 | 5.3752496879551676e-15, PASS |
| Maximum bracket width / B0 | 4.88350583511898e-11, PASS; limit 5e-11 |
| All samples/quantities | 21 stars; 22 groups/star; 462 TSV rows |
| Historical residual/monotone-envelope criteria | PASS |
| Maximum native utilization | 0.54075347985134414; index 20, N_mu |
| Independent refit | PASS: all 10303 components, separate 80/112-digit normal equations |
| Historical Mac independent refit | PASS: all 10303 components |
| Mac/native support grid and active dimensions | Exact at all 43029 rows |
| Mac/native thermal knots | Exact at all 273 rows |
| Native depletion gate | Inside -9.9999997155751619e-7 accepted; outside -1.0000000283979224e-6 refused |
| M20 over-budget mutation | Refused before RHS |
| Native owning profile and B0 | Exact binding PASS |
| One-ULP profile mutation | Refused by existing owning-star guard |
| Changed source bytes / stale frozen token | Refused |
| Independent uncertainty/platform qualification | UNQUALIFIED; proposed ADR-0021 |

Native certificate JSON SHA-256:
`4212cad05bb4244ab75e80b8dfea023e6d0fcee69086c9993bd48d26a40e41fb`.
Native certificate TSV SHA-256:
`38bc85c62621b390ff2d7c8affee1f5ba4e9f0f44233047bb9959f2d685fd5ba`.
Both are under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery/run1/fit`.
The input hash binding is retained with the candidate. The binding test consumes
this native certificate with its native canonical star only; it invokes no trajectory.

## Uncertainty limitation and scientific disposition

The historical fitter assigns zero errors to radius, envelope and neutron/metric
profiles; thermal errors are roundoff multiples; coefficient order differences
can be zero; target chemical adapters inherit B0 characterization. The independent
refit reproduces the assigned error/covariance treatment, not its physical validity.
The historical envelope does not add pointwise numerical uncertainty explicitly;
its interpretation against preflight §11 is recorded for owner review, unchanged.

Cross-platform component differences are measured at each platform's achieved B,
with abscissa differences and linear-alignment diagnostics retained separately.
Several differences exceed the sum of declared errors: Cstar 214/273, N_mu 15/21,
surface gravity 21/21, for example. Small drift is not independent error validation.
The full report retains every component comparison; compact summaries are in Git.
No blanket tolerance, new floor, grid change or fit-window adjustment was made.

The original Python adapter's inherited-tail refusal is preserved. The separate
native adapter recomputes each actual-cut witness using unchanged Phase-5C tail
mathematics. Inherited B0 center/background estimates are explicitly identified,
not presented as per-target refinement evidence. No production owner changed.

## Cost, artifacts and next step

The two complete native static runs consumed 6711.735467 CPU seconds (1.864371 h),
including the refused original adapter. Earlier profile/currentness controls used
215.876689 CPU seconds. Post-archive fits, two independent refits, full comparison,
monitor and native frozen binding used 170.716227 CPU seconds. These measured
scientific/validation processes total about 1.972 CPU-hours; build, packaging and
previous reused Phase-5C characterization costs are separate. Peak producer child
RSS was 285260 KiB; earlier binding controls peaked at 312900 KiB. Summed concurrent
process wall time is not elapsed project time. Detailed receipts are retained.

The new candidate archive contains raw samples, final fit/covariance JSON, TSV,
bindings, independent fit records, comparisons and provenance. Its inventory/hash
is under `evidence/native-candidate-package.json`; both full native runs and the
original Mac ZIP remain outside Git. `RECIPE.md` is the complete static recipe.

The next bounded step is owner review of proposed ADR-0021, then independently
validated per-quantity/per-target error propagation and resolution evidence with
fixed nominal producer grids. This must precede further certificate-dependent
pilots. No certificate promotion, full BNV campaign, or nonlinear stability-ceiling
qualification occurred. Master, dependency masters, Mac baselines and historical
failure records remain unchanged.

Changes are dependency/build, validation-tool numerical-method, engineering,
documentation and candidate generated evidence. Authority: GOVERNANCE §§2–5,
preflight §11, INV-14/15, ADR-0015/0016/0017 and proposed ADR-0020. The original
scientific formulas are reproduced; qualification is not ratification.
