# Original BA13 evidence recovered 2026-10-12

This package preserves the completed original 21-star BA13 frozen-sensitivity evidence. Original files are byte copies; no simulation, source modification, Git commit, merge, or change to archived evidence was performed. New files in the package are recovery inventories, documentation, and verification records.

Validity: CONTROLLED MATHEMATICAL STATIC CERTIFICATE; NO BNV TRAJECTORY. Recovery authenticates historical bytes and completeness. It does not promote this certificate to a physical model, calibration, new numerical qualification, or owner-ratified scientific result. See historical-source/docs/validation/PHASE6A1_CONTROLLED_BNV_PRETRAJECTORY.md:40-79.

Archive root:
/Volumes/Scriptorium/Research/CompactStar/space_recovery/20261002T001142Z/builds/worktrees/CompactStar-phase6a1-controlled-bnv-implementation/build/phase6a1-controlled-bnv-debug

Original historical build path (now a symlink to that archive):
/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-controlled-bnv-implementation/build/phase6a1-controlled-bnv-debug

Original clean source checkout: /Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-controlled-bnv-implementation
Branch: physics/phase6a1-controlled-bnv-implementation
HEAD: e56e6e50040dbcd9dcbee1acecf58843f3dddf1c

## Authentication

Final JSON SHA-256: 6643054d8cba36035af15a15615baf9920d9a56adf3f5211f13ea9f2745f4409
Final TSV SHA-256: 9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab
Both match the user-specified hashes and the committed pretrajectory evidence exactly. All 347 frozen-run files (including stopped attempts) also match both original October 2 offload manifests. Additional selected original build metadata and characterization inputs match those manifests too. The full original 002-source.json and 002-destination.json are included under historical/ so inventories can be checked against their original records.

## Completeness

- 21 achieved target-B stars, indices 0..20, nominal grid zero through -1e-6.
- Every target residual and final bracket width <= tau_B_target = 3.8084532593594452e+46 counts.
- Maximum absolute residual = 5.4009617277691353e+42; maximum bracket width = 3.7197817215720943e+46.
- target_stars.tsv and tangent.tsv: 21 rows each.
- 21 coefficient shards and assembled coefficients.tsv: 21 rows.
- 21 thermal-surface shards and assembled thermal_surface.tsv: 273 rows, 13 exact knots per star.
- 21 profile-metric shards and assembled profile_metrics.tsv: 43,029 rows, 2,049 points per star.
- All three assembled tables are byte-identical to ordered concatenation of their 21 original shards with one header.
- 21 chemical certificate JSON files and 21 TXT certificates; all 21 radial profiles, model files, free-gas tables, window tables, and per-star thermal sources included (10 files per star).
- Final certificate.json: 21 samples and 22 quantity groups; certificate.tsv: 462 rows, 22 quantities per star. Outer half remains indices 10..20.
- Maximum recorded utilization: 0.54075347894390757, index 20, quantity N_mu.
- Required final missing files: NONE. See MISSING_FILES.txt and star_completeness.tsv.

## Original assembly and fitting recovered

The exact successful commands, working directory, UTC timestamps, exit codes, and outputs are preserved in historical_commands.txt and original session excerpts. Assembly completed at 2026-09-19T03:59:42.114Z; successful final fitting completed at 2026-09-19T04:00:27.370Z, both exit 0.

The recovered historical fitter is historical-source/tests/bnv/compare_candidate.py from commit 460ca713a6bd2633b4afbf3926e959b751cea0c6, Git blob 6790b8983d63fb908af43bf2c8e30e3a5fed7e94, SHA-256 9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370. It contains both merge-frozen and frozen modes. A verification-only run over the copied serialized original tables produced JSON and TSV with both exact historical hashes; outputs were written to a separate temporary directory. No stellar calculation was rerun.

The original earlier fit at 2026-09-19T03:59:50.098Z exited 1 with STOP rank-deficient frozen sensitivity fit. The preserved patch at 04:00:15.865Z changes the regression to centered sums before the successful 04:00:27.370Z fit. The final historical fitter contains that centered regression. Earlier failure is not presented as successful production.

## Generation stages and provenance distinctions

See generation_stages.tsv. The earlier target-tangent target table has all common values identical to the final target table; the final version additionally carries u_N_n, u_N_e, u_N_mu. Both original tables are preserved with their own hashes. Earlier and final tangent.tsv are exact byte matches.

The separate failed-attempt1, interrupted-sequential-attempt2 (exit 130), and domain-refusal-attempt3 directories are explicitly nonfinal. Their original logs are included; all intermediate files are hashed and inventoried, with excluded TSVs identified. The successful final shards live only under original/phase6a1-frozen-certificate-final. The committed pretrajectory record states stopped attempts contributed no row to the final certificate; the final shard-to-assembly byte checks establish which data were actually assembled.

## Package layout and verification

original/: original final evidence, earlier target tables, original run logs, selected original build metadata and inputs.
historical-source/: exact Git blobs of the historical fitter, adapters, validation source and governing records; these are recovered historical sources, not copied binaries.
historical/: original offload manifests and selected original command records.
source_inventory.tsv: full authenticated frozen-generation file inventory, exact archived and historical paths, SHA-256, sizes, metadata and inclusion/omission decisions.
directory_inventory.tsv: original frozen-generation directory inventory.
provenance.json, completeness.json, star_completeness.tsv, generation_stages.tsv: generated recovery records.
SHA256SUMS: SHA-256 for every packaged regular file except SHA256SUMS itself (standard self-exclusion).

After extraction, change into this package directory and run:
    shasum -a 256 -c SHA256SUMS

The delivered ZIP is independently extracted into a fresh temporary directory and every SHA256SUMS entry, file-set coverage, ZIP CRC and packaged historical inventory hash is checked. The resulting verification report and ZIP SHA-256 sidecar are delivered beside the ZIP because an archive cannot contain a stable hash of itself. See ZIP_VERIFICATION.json outside the ZIP.

No compiled binaries or unrelated simulation outputs are included. All original final radial profiles are included for independent verification. Metadata preservation means original content and recorded metadata, not a claim that ZIP reproduces original filesystem extended attributes or ownership.
