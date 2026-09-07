# Phase-5B-I — governed structural-response installation and INV-09 closeout

**Date:** 2026-09-06
**Disposition:** **PHASE-5B STRUCTURAL RESPONSE GOVERNED ARTIFACT INSTALLED AND VALIDATED.**

This record covers only the governed installation, regression enforcement, post-installation
validation, and canonical-integration preparation for the already implemented, independently
reviewed, and human-ratified ADR-0011 structural particle-number response. It changes no
production scientific source, accepted formula, domain, numerical tolerance, or EOS policy.

## 1. Authenticated authority and history

| Item | Authenticated value |
|---|---|
| Starting canonical `master` | `a43d02227bf53c3242d3212f81dd71963804f3aa` |
| Phase-5B implementation | `fe08c94ed5ae525a9cb78331c5dd69d9c617d591` |
| Human-ratification commit | `d67ac9be27677aa038a9ac7cdf5f262a32a008ef` |
| Integration branch | `integration/phase5b-global-particle-response` |
| Governed artifact | `tests/baselines/phase5b_structural_response.json` |
| Governed artifact SHA-256 | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` |
| Baseline/test commit | `ab8882255c6f939483068ec40b851b24d2f4575e` |
| Independent review | **PASS WITH EXPLICIT CLAIM-NARROWING CAVEATS** |
| Human authority | **RATIFIED for governed integration** |

The branch history was authenticated as the exact linear chain starting canonical master,
implementation, and ratification before work began. No earlier governed Phase-5B structural
artifact existed. The established governed-artifact hierarchy is `tests/baselines/`; the new
artifact is the ninth governed baseline and the first for the ADR-0011 response.

## 2. Fresh canonical production and installation

The sole canonical producer is
`tests/analysis/produce_particle_number_reference.py`, invoking the compiled
`phase5b_freegas_validation emit` path. The producer solves every contributing star in fresh
output and never reads a historical scratch coefficient or governed baseline.

| Production | Isolated output | SHA-256 | Exact disposition |
|---|---|---|---|
| Fresh generation 1 | `/private/tmp/compactstar-phase5bi.bWbmJ7/generation-1/structural-response.json` | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` | byte-identical to ratified candidate |
| Fresh generation 2 | `/private/tmp/compactstar-phase5bi.bWbmJ7/generation-2/structural-response.json` | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` | byte-identical to generation 1 and ratified candidate |
| Installed governed artifact | `tests/baselines/phase5b_structural_response.json` | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` | installed from fresh generation 1; byte-identical to generation 2 |
| Post-suite regeneration | `/private/tmp/compactstar-phase5bi.bWbmJ7/generation-3/structural-response.json` | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` | byte-identical to governed artifact |

All three independent producer invocations returned zero. Their equality authenticates
deterministic reproduction; it is not a substitute for the ratified independent scientific
review.

## 3. Regression enforcement

`tests/analysis/particle_number_response_regression.py` invokes the canonical producer through
the configured Python interpreter, creates a new isolated evidence directory, and compares the
fresh `structural-response.json` to the governed artifact byte for byte. It reports both hashes,
propagates a producer failure, and explicitly refuses a generated-path/baseline-path
self-comparison. `tests/CMakeLists.txt` registers this as
`phase5b_structural_response_regression`, serial, self-contained, scientific, and Phase-5B
labelled. The baseline is output only; it is never producer input.

The baseline/test commit contains exactly the governed artifact, this comparator, and its CTest
registration. Its tree `f6d8c61152b3a97637d5d5523afaed5b5a3f35a1` is exactly the tree tested before commit;
Git metadata caused no test-relevant change.

## 4. Post-installation validation

Configuration was Debug with AppleClang 17.0.0, Python 3.12.10, and authenticated external data
root `/Users/keeper/Documents/CompactStar/data/compose`. Every suite ran serially with
`--output-on-failure`; no failure was skipped, reclassified, or reinterpreted.

| Gate | Authenticated inventory | Result | Raw CTest rc |
|---|---:|---|---:|
| Focused `particle_number_analytic` plus `phase5b_*` | 9 | **9/9 PASS** | 0 |
| Complete data-free (`-LE external-data`) | 39 | **39/39 PASS** | 0 |
| Complete authenticated external-data suite | 62 | **62/62 PASS** | 0 |

The post-suite producer regeneration also returned the exact governed SHA and exact governed
bytes. `git diff --check` passed. The only pre-documentation changes were the three authorized
baseline/test paths; no unexpected generated or tracked path appeared.

## 5. Historical governed-baseline immutability

The eight pre-existing governed artifacts remained byte-identical before installation, after all
suites, and before commit:

| ID | Artifact | SHA-256 |
|---|---|---|
| C1 | `tests/baselines/passive_cooling_cmf_1p6_debug.tsv` | `8fef2314673fceb939f859612f4befe94117115d6d6b3ad0dcc59d1faa68c9f9` |
| C2 | `tests/baselines/grid_convergence_cmf_1p6_debug.tsv` | `b48519c3e948e9979a385d19facee2777d15955eeb8711b4bdd46b81fef74741` |
| C3 | `tests/baselines/grid_convergence_cmf_1p6_trajectory.tsv` | `d5b753932c0523e67a7f25b460c7494bec1a006a8d01c9e43124cb2e78f0720f` |
| C4 | `tests/baselines/hartle_I_dscmf1_debug.tsv` | `034ecddbd9bd847650429d7dc87d0331ec9e87aca3862ff87594e4bff5b707dd` |
| C5 | `tests/baselines/baryon_number_dscmf1_reference.tsv` | `90d607519cbdf3c4a0bf6ef50cc8fd22a8526b5db0354dc319e96854da29041d` |
| C6 | `tests/baselines/hartle_monopole_dscmf1_debug.tsv` | `caaa0ac0d3219cda0a9fb518b27688afc23c6cdad1ec76a2bcd7359614a8d4e8` |
| H7 | `tests/baselines/tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| H8 | `tests/baselines/tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

Historical-baseline diff: **NONE**. Production-science diff after the ratified implementation:
**NONE**. PB tolerances changed: **NO**.

## 6. Ratified qualifications retained

INV-09 closure does not widen the reviewed evidence. The controlling qualifications remain:

1. PB9's baryon identity is algebraically tautological and is not independent validation.
2. The mutation count combines production-discriminating mutations with local algebraic/null checks.
3. PB9/PB10 propagated budgets are conservative; the strongest evidence is PB10's independent
   `K` reconstruction, PB11 nonlinear closure, and the achieved charge residual.
4. `B_B conditioning = 1` is a cancellation ratio, not a general condition number.
5. PB8 is floor-limited for some species.
6. Generic declared physical jumps require authenticated adapter metadata; ordinary Track-R
   onsets are continuous.
7. The PB13 response bound is a conservative enclosure, not a precision estimate.
8. EOS and tail authority are adapter-owned.
9. PB6 is partition/refinement validation, not an independent physics oracle.

Further scope limits also remain: there is no authenticated free-gas core numerical
`I_Omega` benchmark and no source-qualified `M_max`; the Track-R fixture is the ratified
whole-star free-gas validation fixture, not either source claim. The explicit fixed-isobar
mapping machinery is governed only within its authenticated domain. Genuine EOS jumps remain
conditional on authenticated adapter metadata.

## 7. Closeout boundary

Upon fast-forward integration of this exact branch, ADR-0011's structural contract is
implemented, validated, independently reviewed, human-ratified, regression-protected, and
governed for the ordinary-`NStar` slow-rotation structural particle-number response:
domain-qualified PN1-PN8, whole-star fixed-baryon validation, the Track-R free-gas fixture, and
explicit fixed-isobar mapping machinery.

**INV-09 is VERIFIED / RESOLVED only for that structural layer. INV-11 remains UNRESOLVED.**
Corrected R2006 `Btilde`, reduced chemical `Z/W`, evolved eta, weak-rate evolution,
rotochemical evolution, realistic-EOS off-equilibrium extension, and BNV were not begun. The
next physics layer is a separately governed corrected R2006 chemical-susceptibility preflight.
