# Phase-5C-2RAT corrected chemical coefficient ratification

**PHASE-5C CORRECTED CHEMICAL COEFFICIENTS —
IMPLEMENTED / CANDIDATE-VALIDATED / INDEPENDENTLY REVIEWED /
HUMAN-RATIFIED FOR GOVERNED GENERIC/FREE-GAS SCOPE.**

**CANONICAL INTEGRATION IS STILL REQUIRED BEFORE PHASE-5C IS CLOSED.**

## 1. Authority, identity, and disposition

| Authority | SHA / value |
|---|---|
| Canonical master | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Accepted Phase-5C preflight | `54ec7abac38fa0a32c5fb3a82e424b496361966a` |
| ADR-0013 ratification | `4780121f21010374da2eb50898e90a067795b6e5` |
| Reviewed numerical plan (`PHASE5C1P_SHA`) | `09d1b3c935919ec85f3c629607797ae87c42d8e5` |
| Ratified UQ semantics (`PHASE5C1R_SHA`) | `d3b102d2225c44bea581832b27cdc33eb5874e83` |
| Immutable production predeclaration (`PHASE5C2_PREDECLARATION_SHA`) | `a87f0212c2bd7bfba92db91dfac82447a6561334` |
| Production implementation candidate (`PHASE5C2_SHA`) | `4d78bf4000848ddecc2127daa2f2840872f266f5` |
| Candidate artifact SHA-256 | `a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b` |

The production candidate is the exact head of
`physics/phase5c-corrected-chemical-coefficients` in
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-chemical-coefficients`.
Its direct parent is the immutable predeclaration, whose direct parent is the ratified UQ
semantics. The implementation record preserves this sequencing and the fact that the first
production coefficient operation followed the predeclaration
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:9-24`).

The complete independent Opus review disposition supplied to the owner was:

> PHASE-5C CORRECTED CHEMICAL COEFFICIENTS
> INDEPENDENT REVIEW PASS WITH NONBLOCKING FINDINGS —
> CANDIDATE READY FOR HUMAN RATIFICATION WITH EXPLICIT CAVEATS.

The independent review found **0 BLOCKING**, **0 MATERIAL**, **4 NONBLOCKING**, and
**4 NOTES**. Fable was not needed.

The human-owner disposition is:

> **RATIFIED WITH EXPLICIT INDEPENDENT-REVIEW CAVEATS.**

The owner accepts the production central values, the ADR-0013 implementation,
`numerical_error` accounting, separate `validation_envelope` semantics, GC1-GC12 candidate
validation, GC14 provenance/staleness validation, the deterministic candidate artifact, and
the generic/free-gas scope. The owner does not ratify realistic A18 closure, GC13 as passed,
INV-11, eta evolution, weak rates, neutrino/heating evolution, or BNV.

## 2. Human-ratified central physics

The canonical source axes of `G_y` are `(N_n, N_e, N_mu)` and its units are `count / MeV`.
The human-ratified candidate is approximately

```text
G_y = [[2.904621518e55,  0,                 0],
       [0,                 2.201380364e53, -1.035411553e51],
       [0,                -1.035411553e51,  9.746440586e51]] count / MeV.
```

Global baryon reduction occurs **after** integration. The derived diagnostic is

```text
Q = [[ 2.184981542e53, -1.100609594e51],
     [-1.100609594e51,  9.743848458e51]] count / MeV.
```

`Q` remains derived and is not independent authority. The canonical chemical matrix is

```text
Z = [[4.579303181e-54, 5.172519910e-55],
     [5.172519910e-55, 1.026872798e-52]] MeV / count.
```

The human-ratified structural vector and spin drive are

```text
I_phys = [-1.1637998545112904e47,
          -1.484481985023382e46] count s^2,

W      = [-5.406177501704724e-7,
          -1.584571948010365e-6] MeV s^2,

W = Z I_phys.
```

These values reproduce the full-precision candidate entries
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:230-345`).
For supplied `Omega > 0` and `Omega_dot < 0`, `2 W Omega Omega_dot` has the accepted positive
chemical-imbalance-driving sign. No spin-down law is owned or implemented by this layer; its
only dynamic convenience consumes supplied scalar values
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:81-86`).

## 3. Independent-review confirmations

The owner records the independent review as confirming all of the following:

- production GC9 uses the actual injectable production integration kernel;
- an independent curved-GR oracle agrees at approximately `1e-17` relative scale;
- all six wrong lapse/volume routes remain separated by enormous margins;
- double-onset mapping is correct;
- refusal intervals are excluded and bounded;
- the neutron ULP/source enclosure contains the provider refusal;
- the chemical-tail analytic bound encloses both independent continuations;
- `G_y` was independently reconstructed from free-gas beta equilibrium;
- `Q` and `Z` were independently reconstructed;
- `E_Z` passed exhaustive 64-corner perturbation verification;
- `rho << 1` under the accepted inverse-perturbation criterion;
- `W` was independently reconstructed;
- the two-track `numerical_error` / `validation_envelope` semantics are preserved;
- `V_K` / `V_I` arithmetic independently reproduces;
- the old F2005 versus corrected R2006 gate is genuinely independent and strongly
  correction-sensitive;
- lifetime ownership is safe;
- the GC14 stale-input campaign passes;
- the candidate artifact regenerates byte-identically;
- focused 5/5, data-free 44/44, and complete 67/67 pass with raw rc=0 and no skips;
- all nine historical baselines, EOS/data, and literature are unchanged.

The production record independently preserves the kernel, error, lifetime, test, artifact, and
protected-inventory evidence
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:53-90`,
`:129-170`, `:188-228`).

## 4. Nonblocking findings and owner decisions

### NB-1 — paper-Z axis typing

`PaperZ` currently uses one `ImbalanceChannel` enum for both row and column. ADR-0013 Q3 asked
for semantically distinct output beta-channel and input-lepton orientation. The owner accepts
the current representation **for v1 with explicit limitation** because the current matrix is
symmetric, the supported mappings are one-to-one (`Npe <-> Electron`, `NpMu <-> Muon`), and
independent review found no numerical ambiguity or incorrect result and reproduced the matrix
and channel mapping.

This is a narrow implementation concession, not a precedent that one enum is always sufficient.
Before any extension in which output and input spaces differ, the matrix can become
non-symmetric, species/leptons/channels are added, or channel ordering is no longer one-to-one,
the public API **MUST** introduce distinct semantic types for `BetaChannel` output and
`LeptonInput` input, or an equivalently strong typed representation. No production API change is
required now, and current v1 numerical results are not weakened.

### NB-2 — local H error semantics

The current Track-R production path passes zero provider `h_error` into the local
susceptibility solve. For this candidate scope,
`ChargeNeutralNumberSusceptibility::NumericalError()` represents only local congruence,
factorization, and solve arithmetic uncertainty for the supplied local thermodynamic response.
It is **not** the total uncertainty of the EOS/provider model or background. Provider/background
characterization is accounted for separately at global level, especially through
`E_background`, which dominates the realized `G` uncertainty. The production record already
separates local binary-H arithmetic from independent provider, anchor, and background accounting
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:59-66`, `:88-90`).

Future prose and API documentation must not call the local `NumericalError` a total local
physical-model error. A future realistic EOS provider must supply nonzero provider uncertainty,
or govern it separately, where available. No production code change is required for v1.

### NB-3 — Python bytecode hygiene

Running `chemical_production_validation.py` can create
`tests/analysis/__pycache__/` because bytecode suppression/ignore coverage is incomplete. This is
non-scientific housekeeping only and affects no scientific output, candidate artifact, test
result, provenance, or baseline content. Neither tests nor `.gitignore` is changed by this
ratification. A subsequent integration task must begin from a clean integration worktree and may
address this only as a separately identified non-scientific housekeeping change if governance
permits, without altering the reviewed scientific candidate. This does not block ratification.

## 5. Review notes retained as caveats

### N-a — validation envelope dominated by the PB11 constraint

Approximately 92-95% of `V_K` / `V_I` is contributed by the finite-q PB11 baryon-constraint
term. It is an `O(q^2)` nonlinear finite-spin residual used as empirical validation sensitivity.
It is not a measured linear-coefficient error, certified truncation error, confidence interval,
or achieved-accuracy estimate. `V_I_validation` and `V_W_validation` remain described only as a
`validation_envelope`; `V_W_validation` must not later be quoted as achieved numerical accuracy.
The production record identifies the PB11 term as empirical finite-spin sensitivity and as the
dominant conservative contribution
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:42-51`).

### N-b — PB7 transfer micro-difference

The independent reviewer obtained a PB7-transfer value approximately `9e-7` relatively larger,
with approximately `5e-7` relative impact on `V_I`. This is immaterial and nonblocking because
PB7 is a minor envelope contributor and PB11 dominates `V_K`. Both W validation gates retain
large margins. The frozen envelope is not refitted or altered after production.

### N-c — GC9 M20 mutant construction

M20 uses an equivalent transformed injected susceptibility that reproduces the
inverted-proper-volume wrong integrand; it does not literally change the production weight
expression. Independent review confirmed mathematical equivalence for the tested route and the
intended wrong integrand. Together with M8 and the independent curved-GR oracle, it strongly pins
the production proper-volume exponent. It remains a valid mutation detector, but must not be
overstated as a literal mutation of the production weight expression.

### N-d — thin but valid numerical margins

Current refusal-geometry containment checking evaluates containing-cell edges. This is sufficient
for the current injected Track-R background because that representation is cellwise linear.
Future nonlinear background adapters must provide a stronger extrema contract or explicit
interior bound. The production certificate authenticates containing cells and whole-cell metric
extrema (`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:101-104`).

The tail old-M mutant's mass-upper radius construction is correctly separated from the old-M
route, but its absolute radius separation is small. Independent high-precision verification
confirmed the required inequality and positive margin. Future numerical changes must retain
outward-rounding and bound discipline. Neither point blocks v1.

## 6. Correction-sensitive validation

The independently reviewed absolute old-F2005 versus corrected-R2006 differences are

| Entry | Absolute difference (`MeV / count`) |
|---|---:|
| `Z_npe` | `6.9884708244e-56` |
| `Z_np` | `2.3428105032e-55` |
| `Z_npmu` | `3.7102536060e-54` |
| matrix Frobenius | `3.7256732162e-54` |

The minimum separation divided by combined numerical uncertainty is approximately `4.22e4`.
The production ledger preserves the full-precision separations, uncertainties, and margins
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:401-407`). The owner
accepts that the free-gas candidate proves that corrected R2006 electrostatic machinery is active
and numerically consequential. Free gas does not substitute for A18; GC13 remains
**SOURCE-LIMITED / BLOCKED**.

## 7. Two-track uncertainty ratification

`numerical_error` and `validation_envelope` remain distinct. The candidate values are

```text
E_W_numerical      ~= [1.591009671e-13, 1.484320882e-12] MeV s^2,
V_W_validation     ~= [2.263435255e-11, 5.882694394e-11] MeV s^2.
```

Both immutable predeclared componentwise gates passed. `V_W_validation` is not an error bar,
confidence interval, certified bound, formal truncation error, or achieved accuracy. The
controlling propagation and separation are recorded in ADR-0013
(`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:443-478`, `:507-540`) and the
candidate values are recorded in the production implementation
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:338-345`).

## 8. GC, test, artifact, and protected-inventory disposition

| Gate | Human-ratified candidate status |
|---|---|
| GC1 | PASS |
| GC2 | PASS |
| GC3 | PASS |
| GC4 | PASS |
| GC5 | PASS |
| GC6 | PASS |
| GC7 | PASS |
| GC8 | PASS |
| GC9 | PASS |
| GC9b | PASS |
| GC10 | PASS |
| GC11 | PASS — SOURCE TRACEABILITY ONLY |
| GC12 | PASS |
| GC13 | SOURCE-LIMITED / BLOCKED |
| GC14 | PASS |

The machine-readable disposition is retained at
`docs/validation/phase5c_chemical_coefficient_validation.json:4-21`. The independent review's
fresh prior inventories were focused **5/5 PASS**, data-free **44/44 PASS**, and complete
**67/67 PASS**, each raw rc=0 with no skips. This documentation-only ratification does not rerun
them; the production record preserves those inventories at
`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:188-205`.

The candidate artifact regenerated byte-identically in fresh producer directories at SHA-256
`a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b`. It remains a candidate,
not an installed governed baseline
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:195-203`, `:226-228`).
All nine governed historical baselines, fourteen authenticated external-data files, and twenty-eight
literature files retain their entry hashes
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:207-209`). This
ratification changes none of them and changes no EOS/data, production code, or test.

## 9. Ratified scope and unresolved boundary

The human-ratified candidate scope comprises the corrected R2006 local-neutral susceptibility
adapter, global `G_y` integration, global baryon reduction, chemical `Z`, structural spin drive
`W`, the generic/free-gas Track-R mathematical fixture, candidate numerical-error accounting,
validation-envelope methodology, lifetime/provenance semantics, GC1-GC12, and GC14.

The following are not ratified and are not implemented: GC13 realistic A18 closure; eta
storage/evolution; weak reaction rates; neutrino enhancement; rotochemical heating; thermal
coupling; time-dependent coefficients; superfluidity; and BNV. No downstream evolution is
authorized by this ratification.

INV-09 remains **VERIFIED / RESOLVED**, unchanged. INV-11 remains **UNRESOLVED**. Phase-5C
coefficient-level redshift and sign semantics are accepted, but that does not resolve evolved
chemical-state ownership, storage, coefficient evolution, or reaction/evolution coupling.

## 10. Ratification result

ADR-0013 remains **ACCEPTED**. Phase-5C corrected coefficients are **IMPLEMENTED /
CANDIDATE-VALIDATED / INDEPENDENTLY REVIEWED / HUMAN-RATIFIED** for the governed generic/free-gas
scope and are **NOT YET CANONICALLY INTEGRATED**. GC13 remains **SOURCE-LIMITED / BLOCKED** and
INV-11 remains **UNRESOLVED**. This record installs no baseline, performs no canonical merge, and
authorizes no eta evolution, weak rates, heating/cooling, realistic A18 closure, or BNV.

## 11. Subsequent canonical-integration status — 2026-09-10

Section 10 remains the historical disposition of this ratification gate. Subsequent governed
integration preserved this reviewed candidate byte-for-byte, installed a separately classified
baseline, and added a fresh-producer regression. The generic/free-gas coefficient scope is now
**CANONICALLY INTEGRATED / GOVERNED / CLOSED**, while every ratified caveat remains controlling.
GC13 remains **SOURCE-LIMITED / BLOCKED**, global INV-11 remains **UNRESOLVED**, and no eta
evolution, weak rates, heating/cooling, realistic A18 closure, Phase-5D, or BNV was implemented.
See `docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_INTEGRATION.md`.
