# ADR-0012 A1 governed baseline supersession and integration record

**UNIT-1I CLOSEOUT — 2026-09-06**

**C1-C6 GOVERNED SUPERSESSION VALIDATED.**

**H7/H8 RETAINED — NO SUPERSESSION.**

## Authority and authenticated history

- Canonical pre-integration master: `1d22dd1f5a0d1afa18c4cedebb36b28fdae49df4`.
- Ratified implementation (`UNIT1_SHA`): `b3ce4f1303dbab68b68a82614c944c269cefebdc`.
- Human-ratification commit (`RATIFICATION_SHA`):
  `432e64c4ed06108915184022584bfaf9dfec7219`.
- Baseline-only commit (`BASELINE_SHA`): `cca6570b742a2b21258bd3653710caa6da48e57f`.
- Supersession branch: `fix/relativistic-unit-boundary-a1-baselines`.
- Fresh worktree:
  `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-relativity-units-a1-baselines`.
- Human authority: `docs/validation/RELATIVISTIC_UNIT_BOUNDARY_RATIFICATION.md`.
- Independent-review authority: the independently reproduced review and its **PASS WITH
  NONBLOCKING FINDINGS** disposition, recorded in the ratification record, including fresh
  53-test and 30-test inventories, unit/detector checks, and repeated C5/C6 regeneration.

The authenticated ancestry at entry was exactly canonical master -> `UNIT1_SHA` ->
`RATIFICATION_SHA`. Canonical local `master`, `origin/master`, and live `refs/heads/master`
all named the canonical pre-integration master. Local, upstream, and live
`fix/relativistic-unit-boundary-a1` all named `RATIFICATION_SHA`.

## Producer environment and input authority

All six artifacts were freshly regenerated from the ratified code at `RATIFICATION_SHA` in
fresh scratch directories outside `tests/baselines/`. The authenticated environment was macOS
arm64, AppleClang 17.0.0, Debug configuration, GSL 2.7.1, Python 3.12.4, and NumPy 2.4.0. The
external-data authority supplied to every producer was
`/Users/keeper/Documents/CompactStar/data/compose`. Producer identity, argument order, governed
path, and expected bytes were taken from
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_CANDIDATES.json`; no cached candidate artifact was
used as a regeneration substitute.

## Fresh regeneration and exact-hash gate

| ID | Governed artifact | Authenticated canonical producer and command form | Historical SHA-256 | Fresh / ratified SHA-256 | Exact? |
|---|---|---|---|---|---|
| C1 | `tests/baselines/passive_cooling_cmf_1p6_debug.tsv` | `tests/passive_cooling_regression DATA_ROOT --emit OUTPUT` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` | `8fef2314673fceb939f859612f4befe94117115d6d6b3ad0dcc59d1faa68c9f9` | YES |
| C2 | `tests/baselines/grid_convergence_cmf_1p6_debug.tsv` | `tests/grid_convergence_cmf DATA_ROOT --emit-dir OUTPUT_DIR` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` | `b48519c3e948e9979a385d19facee2777d15955eeb8711b4bdd46b81fef74741` | YES |
| C3 | `tests/baselines/grid_convergence_cmf_1p6_trajectory.tsv` | `tests/grid_convergence_cmf DATA_ROOT --emit-dir OUTPUT_DIR` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` | `d5b753932c0523e67a7f25b460c7494bec1a006a8d01c9e43124cb2e78f0720f` | YES |
| C4 | `tests/baselines/hartle_I_dscmf1_debug.tsv` | `tests/hartle_moment_inertia_cmf DATA_ROOT --emit OUTPUT` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` | `034ecddbd9bd847650429d7dc87d0331ec9e87aca3862ff87594e4bff5b707dd` | YES |
| C5 | `tests/baselines/baryon_number_dscmf1_reference.tsv` | `tests/baryon_number_cmf DATA_ROOT --emit OUTPUT` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` | `90d607519cbdf3c4a0bf6ef50cc8fd22a8526b5db0354dc319e96854da29041d` | YES |
| C6 | `tests/baselines/hartle_monopole_dscmf1_debug.tsv` | `tests/hartle_monopole_regression DATA_ROOT --emit OUTPUT` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` | `caaa0ac0d3219cda0a9fb518b27688afc23c6cdad1ec76a2bcd7359614a8d4e8` | YES |

Only after all six fresh hashes matched were those emitted bytes installed. The baseline-only
diff contains exactly the six paths in the table: 100 insertions and 100 deletions across six
files. There was no production diff and no test-source diff relative to `RATIFICATION_SHA`.
No tolerance, test semantic, or production change was made.

## Retained governed artifacts

| ID | Governed artifact | Historical and retained SHA-256 | Result |
|---|---|---|---|
| H7 | `tests/baselines/tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` | byte-identical; no supersession |
| H8 | `tests/baselines/tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` | byte-identical; no supersession |

All eight historical hashes were authenticated before regeneration. After installation and
again after the complete suite, the governed set was exactly C1-C6, H7, H8. Each C1-C6 governed
file also compared byte-identically with its fresh producer output after the suite.

## Post-supersession validation

| Gate | Inventory | Raw result | Failures | Return code |
|---|---:|---:|---:|---:|
| Data-free / self-contained suite | 30 | 30 / 30 PASS | 0 | 0 |
| Complete external-data suite, serial | 53 | 53 / 53 PASS | 0 | 0 |

The complete suite was executed with one test job. The three historical Category-B comparator
failures disappeared under the governed successors; no A/B/C reinterpretation was used. Test
processes altered no tracked file. `git diff --check` passed and the worktree was clean at
`BASELINE_SHA` before this documentation-only closeout.

## Scientific closeout boundary

ADR-0012 A1 is **ACCEPTED, IMPLEMENTED, INDEPENDENTLY REVIEWED, HUMAN-RATIFIED, GOVERNED
BASELINES SUPERSEDED AND POST-SUPERSESSION VALIDATED** for ordinary-`NStar` A1 scope.

Upon the authorized fast-forward of this package to canonical master, ADR-0012 A1 is canonically
integrated and the ADR-0011 Q4 **RELATIVISTIC-UNIT-BOUNDARY PREREQUISITE** is closed. This closes
only that unit-boundary prerequisite of Q4; it does not complete ADR-0011 as a whole.

- PB7: **NOT COMPLETE** as a Phase-5B production-validation claim.
- INV-09: **INTENDED BUT UNVERIFIED / UNRESOLVED**.
- INV-11: **UNRESOLVED**.
- Source-qualified `M_max` selection semantics: **UNRESOLVED**.
- `MixedStar` unit reconciliation: **NONBLOCKING FUTURE DEBT** requiring separate authority.
- Phase-5B production implementation: **NOT YET BEGUN**.
- Btilde, paper Z/W, rotochemical evolution, and BNV: **NOT BEGUN**.

The independent-review nonblocking findings remain durable: `MixedStar` is outside A1 and
retains its historical conversion behavior; historical status prose is preserved as historical
context where explicitly labelled; and future major validation campaigns should retain durable
focused-selection provenance. None broadens this closeout.
