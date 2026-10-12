# Native static sample recovery — checkpoint, not certificate qualification

Two complete Linux Release runs regenerated all 21 achieved target-B stars,
21 Phase-5C certificates/Phase-5D adapters, coefficients, 273 thermal rows and
43029 profile-metric rows. All five raw tables and all 231 profile/adapter files
are byte-identical between runs. The target spacing, radial 80000 and EOS 8192
parameters, active species, support and 13 thermal knots were unchanged.

| Check | Result |
| --- | --- |
| Maximum achieved target residual / B0 | 5.3752496879551676e-15 |
| Maximum bracket width / B0 | 4.88350583511898e-11 |
| Required limit for each | 5e-11 |
| Static sample integrity and sampled support | PASS |
| Native exact-profile binding | PASS |
| One-ULP profile mutation | Refused by existing owning-star guard |
| Changed source-certificate bytes | Refused: StaleChemicalSourceBytes |
| Final 22-quantity fit, covariance, envelopes and thresholds | Pending historical archive audit |
| Independent refit and uncertainty-qualified Mac agreement | Pending final fit |
| Native depletion boundary and M20 over-budget refusal | Not run; require final native certificate |

The original Python adapter refused an inherited Debug B0 tail witness on the
new Release profile. That failure is preserved. The separate native adapter uses
the unchanged Phase-5C tail mathematics at each actual profile cut; the original
producer files retain their authenticated historical bytes. Inherited B0
center/background estimates are explicitly identified, not represented as new
per-target refinement evidence. See `tail-witness-mismatch.json` and each profile's
`tail-witness.json` for this limitation. Final uncertainty qualification must audit it.

New-work CPU cost for the two static runs is 6711.735467 seconds (1.864371 hours),
including the original refused adapter. The separate profile/currentness controls
used 215.876689 CPU seconds. Maximum child RSS was 285260 KiB for the producer
stages and 312900 KiB for the controls. Concurrent process wall times are not
elapsed campaign time. Prior authenticated Phase-5C characterization was reused
and its cost is excluded. Detailed command/Slurm/node/resource receipts are retained.

The completed sample archive is outside Git at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery/native-static-inputs.tar.gz`
(approximately 106 MiB), SHA-256
`d2c80ce474c7276d09a619ae9d5274a64a74cb6ac1cbf6d36b2281fefa01ad30`.
Its compressed stream was read back and its checksum verified. It contains the
first run's complete raw inputs, repeat hash manifest, integrity/resource receipts,
input-authentication receipts and build environment. Both full uncompressed runs
remain on the cluster. Git contains compact manifests and validation evidence.

The owner located completed historical BA13 evidence on Scriptorium and is
packaging it with hashes. Its cluster path/hash has not yet been supplied. The
next bounded step is to authenticate that package, establish the original fitting
semantics from its implementation/raw tables/final JSON, assemble the native
candidate with achieved B abscissae and fixed indices 10–20, independently refit,
and complete the native certificate boundary/refusal controls. No final native
frozen certificate is published, attached to a star, or promoted at this checkpoint.
No trajectory or nonlinear stability-ceiling qualification was run.

Scientific impact: production owners and protected Mac baselines are unchanged;
new candidate numerical inputs are generated on Linux, with the documented
actual-cut tail-witness adaptation. Authority and change classes are listed in
`RECOVERY_AUDIT.md`. This checkpoint is not scientific ratification.
