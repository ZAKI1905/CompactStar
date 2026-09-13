# Phase-5D-1 controlled rotochemical evolution — human-owner ratification

**Date:** 2026-09-12

> **PHASE-5D CONTROLLED NON-SUPERFLUID ROTOCHEMICAL EVOLUTION —
> IMPLEMENTED / CANDIDATE-VALIDATED / INDEPENDENTLY REVIEWED /
> HUMAN-RATIFIED FOR CONTROLLED FROZEN-V1 SCOPE —
> NOT YET GOVERNED-BASELINE PROMOTED OR CANONICALLY INTEGRATED.**

**Final disposition: A — PHASE-5D CONTROLLED EVOLUTION HUMAN-RATIFIED FOR
CONTROLLED FROZEN-V1 SCOPE — READY FOR GOVERNED REPRODUCIBILITY-HARNESS
REPAIR.**

This is a documentation/governance ratification. It changes no production
source, test, candidate artifact, baseline, EOS/data, literature, or trajectory
evidence. It does not install a governed Phase-5D baseline or integrate the
candidate into canonical `master`.

## 1. Authenticated authority

| Item | Authenticated value |
|---|---|
| Canonical `master` | `d019ae390be4f5e3daba05039903485cb497e397` |
| Candidate branch | `physics/phase5d-controlled-rotochemical-evolution` |
| Candidate entry / `PHASE5D1_CANDIDATE_SHA` | `3486b972f71f57e8351fa8320c1ffb250fcd5c42` |
| Coupled implementation / `PHASE5D1_EVOLUTION_SHA` | `d3670f6d4e021def0483909b6d2fdeed1c6973a4` |
| Candidate artifact | `docs/validation/phase5d1_controlled_evolution_candidate.json` |
| Candidate artifact SHA-256 | `47da15fa9e32be095a78d1b78ed17c3ffa14e5ef9be5528db3780a76afd0c079` |

The worktree, local branch, upstream tracking ref, and live branch ref were
authenticated at the candidate entry before this documentation change. The
evolution commit is the immediate parent of the candidate evidence commit. The
candidate artifact remained byte-identical throughout ratification. Canonical
local, `origin`, and live `master` remained at the SHA above. The implementation
and evidence lineage is recorded in
`PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_IMPLEMENTATION.md`; the accepted
contract is ADR-0014.

## 2. Independent-review evidence and owner disposition

The independent scientific review was the Opus review, not the implementation
team's earlier read-only specialist work. Its four specialist/review roles were:

- A — source physics;
- B — numerics;
- C — architecture/provenance; and
- D — lead red-team;

with lead regeneration/orchestration in addition. All scientific specialists
returned **PASS**. The lead disposition was:

> **B — PHASE-5D CONTROLLED EVOLUTION INDEPENDENT REVIEW PASS WITH
> NONBLOCKING FINDINGS — CANDIDATE READY FOR HUMAN-OWNER RATIFICATION WITH
> EXPLICIT CAVEATS.**

The totals are **0 BLOCKING, 0 MATERIAL, 7 NONBLOCKING, and 12 NOTES**. Fable
was **not needed**.

Independent regeneration authenticated the candidate artifact and matched all
scientific/equality-bearing regenerated fields. It reproduced all six
trajectories, step statistics, oracle errors, `Ltilde`, and thermal-source
hashes; confirmed all 33 protected paths unchanged; passed the Phase-5B and
Phase-5C regressions; and obtained 50/50 data-free and 73/73 authenticated
tests, with no failures or skips. These results are accepted as completed
independent review; this ratification does not independently reconstruct the
science.

The human owner **RATIFIES** the candidate for its explicitly declared
**CONTROLLED / MATHEMATICAL / ARCHITECTURE frozen-v1 scope**, retaining every
independent-review caveat below.

## 3. Exact claim ratified and realistic claim rejected

The owner ratifies exactly this controlled-model statement:

> **We have implemented a complete controlled, non-superfluid,
> modified-Urca, two-channel npe-mu rotochemical evolution on a fixed free-gas
> neutron-star background.**

For this statement, "complete" means that the candidate evolves
`T_infinity`, `eta_npe^infinity`, and `eta_npmu^infinity` with spin-down
compression drive, two-channel chemical relaxation, electron/muon cross
coupling through `Z`, nonequilibrium modified-Urca neutrino emission, chemical
heating, photon cooling, heat capacity, GR-weighted global coefficients, a
prescribed spin history, frozen `Z/W`, and a fixed stellar background.

The owner does **not** ratify this statement:

> We reproduced a realistic FR2005 neutron-star thermal history.

That statement remains **false / unsupported**. The ratified trajectory is not
a realistic Fernández–Reisenegger thermal reproduction and does not establish
the physical validity of the prescribed 1-ms spin history for this free-gas
fixture.

## 4. Ratified state, chemical dynamics, and frozen coefficients

The chemical state ordering is `(Npe, NpMu)`, stored in MeV, with
`eta^infinity = exp(nu) eta_local` and
`xi = eta^infinity/(k_B T_infinity)`. The convention `R_l > 0` denotes the net
reaction that destroys positive `eta_l`, and dissipation requires
`eta_l R_l >= 0`.

The ratified frozen-v1 chemical equation is

```text
dot eta = -Z R + 2 W Omega OmegaDot.
```

Both `Z` cross terms are retained. `dot Z = 0`; no time-dependent-`Z`
correction is included or ratified in this scope.

## 5. Ratified Urca and thermal ledger

For each controlled modified-Urca channel, the ratified ledger is

```text
Lnu_eq   = Ltilde T^8
Lnu_full = Ltilde F_M(xi) T^8
DeltaLnu = Ltilde [F_M(xi) - 1] T^8
R        = (Ltilde/k_B^erg) T^7 H_M(xi)
LH       = sum_l eta_l R_l
Pnet     = LH - DeltaLnu - Lnu_eq_controlled - Lgamma - Lnu_other.
```

Chemical heating crosses the MeV-to-erg boundary exactly once. The same
`Ltilde` authority feeds equilibrium cooling, full/incremental nonequilibrium
neutrino emission, the reaction rate, and chemical heating. Historical
placeholder modified-Urca cooling is not simultaneously active for the
controlled channels.

`F_D`, `H_D`, `F_M`, and `H_M` are ratified as implemented under ADR-0014. The
last modified-Urca `H_M` term uses `pi^8`. Its source status remains
**CONFIRMED PRINTED TYPO / INTERNAL SOURCE INCONSISTENCY; NO PUBLISHED ERRATUM
LOCATED**. This ratification does not upgrade that finding to a published
erratum.

## 6. Controlled trajectory ratified

The first trajectory covers `0` through `1e10 yr`. Its final state is
candidate-validated and independently reviewed:

| Quantity | Final value |
|---|---:|
| `Tinf` | approximately `1.0304414596e6 K` |
| `Tsurface_inf` | approximately `7.50828055e4 K` |
| `eta_npe` | approximately `0.0215104459 MeV` |
| `eta_npmu` | approximately `0.0242991323 MeV` |
| `xi_e` | approximately `242.2440954` |
| `xi_mu` | approximately `273.6494325` |

The ratified qualitative sequence is:

1. initial cooling;
2. spin-driven chemical-imbalance growth;
3. incremental beta cooling at small `xi`;
4. summed incremental beta power crossing zero near `3.36e5 yr`;
5. continued cooling to a temperature minimum near `1.585e6 yr`;
6. rotochemical reheating;
7. a temperature maximum of approximately `6.525e6 K` near `5.96e6 yr`;
8. slow late-time decline;
9. approach to/tracking of quasi-steady spin/reaction balance by the chemical
   imbalance; and
10. late thermal balance `LH - DeltaLnu approximately equals Lgamma`.

These are controlled benchmark results only. The summed beta-power crossings
must not be confused with single-channel source-function roots.

## 7. Analytic, asymptotic, and quasi-steady validation

The independently confirmed modified-Urca roots are:

- incremental `LH = DeltaLnu` root:
  `xi approximately 4.90971002892413`; and
- full heating-minus-full-neutrino root:
  `xi approximately 5.63371746764834`.

Their definitions are part of the ratified result. The large-`xi`
modified-Urca relation is
`DeltaLnu/LH -> 3/8` and net beta heating over `LH -> 5/8`. The late candidate
value `DeltaLnu/LH approximately 0.37571485` is consistent with finite-`xi`
corrections.

The unchanged predeclared `1e9..1e10 yr` window contains 35 eligible
checkpoints. The fitted slopes are `0.1413043725` for the electron channel and
`0.1416691560` for the muon channel, compared with
`1/7 = 0.1428571429`. This supports the controlled modified-Urca quasi-steady
solution.

## 8. Numerical-solver scope and cancellation-sensitive residual

Scaled RKF45 is ratified as adequate **for this controlled benchmark only**.
The baseline settings are `rtol=1e-7` and
`atol=(1e-12,1e-18,1e-18)`; the refined settings are `rtol=1e-9` and
`atol=(1e-14,1e-20,1e-20)`. The independently reviewed baseline/refined
comparisons satisfy the frozen temperature and eta acceptance criteria.

| Common-checkpoint maximum relative difference | Reviewed value |
|---|---:|
| `Tinf` | `1.7408584683e-5` |
| `eta_npe` | `1.5866651419e-7` |
| `eta_npmu` | `1.9479546822e-7` |

The baseline/refined runs used 8,828/10,415 accepted steps,
1,618/2,183 rejected steps, and 71,505/86,004 RHS evaluations, respectively.
This is numerical convergence and execution evidence for the controlled
fixture, not a general solver certification.

No general claim that RKF45 is adequate for future realistic rotochemical
evolution is ratified. The late evolution is stability-bound, with
`h|lambda|` approximately `3.5-3.7`, close to the relevant Fehlberg real-axis
stability boundary. Solver strategy must therefore be revisited before direct-
Urca or superfluid extensions if needed.

The approximately `1.02%` maximum relative baseline/refined difference in
`Pnet` occurs near the reheating temperature maximum, where `Pnet` crosses
zero. The absolute discrepancy is tiny relative to the gross energy ledger.
The independent review classifies this as **NONBLOCKING**. It is not evidence
of achieved 1% accuracy for the gross thermal terms.

## 9. Mandatory physical and provenance caveats

### 9.1 Super-Kepler prescribed spin history

The controlled spin history prescribes `P0=1 ms` for a free-gas fixture with
approximately `M=0.624 Msun` and `R=12.77 km`. Independent review found this
super-Kepler / physically inadmissible under empirical mass-shedding estimates
for roughly the first `~2 Gyr`.

The prescribed spin law is therefore accepted **only as a mathematical driver
for the frozen-`W` architecture benchmark**. No physical pulsar interpretation
of this trajectory is permitted. A future physically interpreted benchmark
requires, at minimum, a physically admissible initial spin such as `P0` of
order `>=2 ms` for this fixture, or a more realistic/compact stellar model with
an authenticated mass-shedding limit. The original `P0=1 ms` trajectory is not
rerun or changed; it remains ratified only for its declared mathematical
purpose.

### 9.2 Envelope provenance label

The implementation/metadata label `iron Potekhin1997` is incorrect. The actual
implemented envelope is the **FR2005 eq. (49) / PCY97 fully accreted-envelope
fit**. This is a label/provenance-description error only; no equation or
trajectory value changes in this ratification. The candidate must not be
described as an iron-envelope result. Before any governed Phase-5D baseline is
installed, the future governed metadata label must be corrected without
altering the ratified trajectory.

### 9.3 `Lother` interface limitation

`Lother=0` in controlled v1. The present implementation reads only the PBF
component rather than a fully configurable total "other neutrino"
contribution. This has zero numerical impact in v1 because those sources are
hard-disabled. The interface must be corrected/reviewed before configurable
other-neutrino physics is enabled; it is not fixed here.

### 9.4 Fixed-background thermal/EOS limitation

The candidate uses a free-gas fixed-background thermal adapter, Sommerfeld
entropy treatment, `m*=mu`, and a 160-node log-temperature heat-capacity
interpolation. This is accepted for the controlled mathematical benchmark. It
is not finite-temperature EOS authority and not realistic FR2005
thermodynamics. The documented small interpolation sawtooth and caching caveats
remain.

### 9.5 Mutation-evidence accounting

All required mutation families have detectors, so the mutation coverage is
adequate. The 22-case transformed-evaluation loop is not 22 independent
production mutant builds. Independent review identified approximately eight
detector-distinct independent groups and two near-aliases. No unqualified
claim of "16 independent mutation families" is retained. Mutation-count
accounting requires conservative wording.

### 9.6 Reproducibility-harness precondition

Before any governed Phase-5D baseline installation, the validation harness
must be repaired and demonstrated from a genuinely fresh repository/build
context. The candidate harness depends on gitignored/external build-state
details, including `ROOT/'build'` assumptions, retained entry-hash state,
trusted `-entry.rc` text, and a "self-contained" label that is not literally
self-contained. The successful independent regeneration means this does not
invalidate the candidate; it does block baseline promotion until repaired.

At minimum, the repair must remove dependence on gitignored build state and
hard-coded `ROOT/build` assumptions, regenerate rather than trust saved
`-entry.rc` provenance text, and correct the inaccurate self-contained label.

### 9.7 Candidate-JSON authority and crossing diagnostics

The candidate JSON contains author-added narrative/postprocessed keys that are
not generated directly by the in-repository validator. Independent review
checked their values. Before baseline promotion, governance must define which
fields are producer-authoritative and which are postprocessed narrative
metadata.

All crossing times are fixed-checkpoint brackets with log interpolation. They
are estimates, not event-localized ODE roots.

## 10. Frozen-v1 scope boundary

The ratified scope is limited to non-superfluid, modified-Urca-only electron
and muon channels on a fixed free-gas background, with frozen `Z/W`, a
prescribed external spin history, and the declared mathematical
normalizations.

Not included are `dot(Z)`, direct-Urca evolution, superfluidity, crust
chemistry, realistic A18, realistic Yakovlev/FR2005 rate normalization, BNV,
state-coupled spin, physically admissible spin histories, or chemical
dependence of the crust/envelope. No realistic extension is inferred from the
controlled result.

## 11. INV-11 owner disposition

After this owner ratification:

| Subpart | Owner disposition |
|---|---|
| INV-11b | **RESOLVED FOR CONTROLLED FROZEN-V1 SCOPE / NOT YET CANONICALLY INTEGRATED** |
| INV-11c | **RESOLVED FOR CONTROLLED FROZEN-V1 SCOPE / NOT YET CANONICALLY INTEGRATED** |
| INV-11d | **RESOLVED FOR CONTROLLED FROZEN-V1 SCOPE / NOT YET CANONICALLY INTEGRATED** |
| INV-11e | **RESOLVED FOR FROZEN-V1 COEFFICIENT SEMANTICS / NOT YET CANONICALLY INTEGRATED** |
| INV-11f | **RESOLVED FOR CONTROLLED FROZEN-V1 ODE/SOURCE COUPLING / NOT YET CANONICALLY INTEGRATED** |

Global INV-11 remains **UNRESOLVED** until canonical integration and because
the broader realistic extensions remain open: realistic A18/FR2005
normalization, `dot(Z)`, direct Urca, superfluidity, crust/envelope chemical
dependence, state-coupled spin, physically admissible spin histories, and BNV.

## 12. Promotion and integration boundary

No governed Phase-5D baseline is installed. No canonical integration or merge
is performed. Canonical `master` does not contain the Phase-5D implementation.
Realistic A18 work has not begun. BNV has not begun.

The next authorized work is not part of this ratification. Before any governed
baseline promotion, it must also decide how to correct the accreted-envelope
label in future governed metadata without changing the ratified scientific
trajectory.

## 13. Exactly one recommended next action

Perform a bounded Phase-5D reproducibility-harness and governed-artifact
preparation task before any baseline installation. Remove dependence on
gitignored build-state and hard-coded ROOT/build assumptions, make provenance
checks regenerate rather than trust saved -entry.rc text, define the exact
producer-authoritative candidate-artifact schema versus postprocessed
narrative metadata, and correct the future governed envelope provenance label
to FR2005 eq. (49) / PCY97 fully accreted without changing the ratified
scientific trajectory. Re-run the independent scientific equality-bearing
checks and Phase-5B/5C regressions after those non-scientific repairs. Only
then consider Phase-5D baseline promotion and canonical integration. Do not
change the ratified P0=1 ms benchmark trajectory; retain it explicitly as a
mathematical/super-Kepler control. Do not begin realistic A18 closure or BNV
during that repair task.

That next action is not begun by this record.

---

## 14. Governed integration follow-up — 2026-09-12

The separately authorized promotion installed the unchanged controlled
trajectory as a fresh-producer governed baseline. After a launch-level Python
CTest isolation repair outside all protected source identity, complete
validation passed 51/51 data-free and 75/75 authenticated tests, followed by a
no-cleanup governed Phase-5D regeneration. The branch is ready for canonical
fast-forward. This follow-up does not revise the ratification or its scientific
caveats. See
`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md`.
