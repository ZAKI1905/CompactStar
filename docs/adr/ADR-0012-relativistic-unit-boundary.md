# ADR-0012: relativistic unit boundary for canonical ordinary-star backgrounds

**Status: PROPOSED — 2026-09-05**

**Owner decisions pending. No production correction or baseline supersession is authorized
by this PROPOSED record. ADR-0011 remains ACCEPTED and its Q4 implementation block remains.**

## Context and authority

The canonical TOV integrator uses GSL G,c and publishes a literal mass ratio to GSL's solar
mass. NStar multiplies that ratio by Zaki's nominal solar mass length, converts pressure and
energy with Zaki's modern-G factors, and carries the original GSL nuprime. The resulting
background is internally inconsistent: the full-star mass identity has a systematic
`2.56486e-4` relative residual and the nuprime identity up to `7.33756e-5`. Coherent conversion
of the same solved state closes those algebraic identities to roundoff and reduces the leading
homogeneous/complete-star particle derivative discrepancy from `9.1742e-4` to `1.3684e-6`.
This is scratch evidence, not a production PB7 pass
(`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`).

Governing authority: `GOVERNANCE.md:15`, `:39`, `:65`, `:88`; accepted canonical numerical
primitive (`docs/adr/ADR-0005-canonical-tov-numerical-primitive.md:208`); geometry
(`docs/adr/ADR-0004-proper-volume-and-metric-measure-ownership.md:29`); angular and monopole
contracts (`docs/adr/ADR-0006-hartle-first-order-physical-normalization.md:243`,
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:418`); measure/surface authority
(`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:61`,
`docs/adr/ADR-0009-tov-surface-event-and-termination.md:47`); unit prerequisite
(`docs/adr/ADR-0011-particle-number-structural-response.md:92`).

The eventual change is **scientific-semantic**, with a bounded architecture change and
independent numerical validation. A successful round trip or old baseline is not sufficient
scientific authority for a conversion. This document does not adjudicate a chemical convention.

## Five owner decisions

### Q1 — actual solve convention

**Proposed answer: retain the current canonical physical solve, GSL 2.7.1's
G=6.673e-8 cm^3 g^-1 s^-2 and c=2.99792458e10 cm/s.** This selects the numerical convention of
this calculation, not the ideal modern estimate of G for all future work. A later change to
modern G requires a separately authorized physical re-solve/revalidation.

Alternative: migrate canonical TOV to Zaki/CODATA G=6.67430e-8 CGS now, with a separately
chosen public mass convention. The audit's independent enthalpy re-solve quantifies that
alternative; it preserves the qualified printed Structure-1 bins but changes the physical
star. It is not required to reconcile the existing solved spacetime
(`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`).

### Q2 — authoritative geometric representation

**Proposed answer: all ordinary StarProfile m,eps,p and nuprime must represent the same
Q1 TOV solution, using one explicit conversion owner.** For CGS input and geometric km:

```
r_km = r_cm/1e5
m_km = G_TOV*m_grams/c_TOV^2/1e5
eps_km_minus2 = G_TOV*rho_cgs/c_TOV^2*1e10
p_km_minus2 = G_TOV*p_cgs/c_TOV^4*1e10
nuprime_km_inverse = nuprime_cm_inverse*1e5
```

The existing nu reconstruction matches to the corrected geometric surface mass. Geometry's
formulas, RotationSolver equations, fixed-central-energy normalization, energy measure,
finite-pressure surface and AngularVelocity c conversion are unchanged. All five TOV
identities, not just lambda or pprime versus carried nuprime, must be validated.
The derivation, current source paths and radial residuals are in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`.

### Q3 — public solar mass and geometric mass transport

**Proposed answer: keep current public `tp.m`, `SeqPoint::m` and NStar mass semantics as the
literal ratio `m_grams/(1.98892e33 g)`, explicitly named/documented as the historical GSL
normalization. Preserve it independently of the geometric mass in existing metadata.**
Public central e,p remain physical CGS. Copy the physical public inputs where they already
exist; otherwise use the matching inverse conversion. Do not reconstruct literal solar mass
by dividing geometric mass by Zaki's nominal solar length.

Zaki `SUN_M_KM` denotes `(GM_sun)^N/c^2` from IAU 2015 B3; a nominal gravitational-parameter
ratio is a different public convention. It may be useful as a future separately named output,
but no silent relabelling or denominator substitution is accepted here. The existing
`SeqPoint` and `StarProfile` fields suffice; no new TOVPoint schema is required for A1
(`CompactStar/Core/TOVSolver.hpp:316`; `CompactStar/Core/StarProfile.hpp:353`, `:381`;
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`).

### Q4 — minimum migration

**Proposed answer: Option A1, including the existing dual-metadata handling of Q3.**
Introduce the small proposed top-level `CompactStar/RelativityUnits.hpp` boundary owner;
update both NStar construction paths, FinalizeSurface/sequence and Mass accessor handling,
plus StarContext's physical-density inverse and TbDefinition's physical-density threshold.
Update only the necessary contract comments and independently built fixture adapters.

No canonical TOV numerical behavior, G/c, EOS, event, grid, search or derivative-policy change.
No RotationSolver or Geometry equation changes. No MixedStar, BNV or chemical/evolution
migration. Exact file-level MUST/MAY/MUST-NOT scope and U1–U14 are in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`.

A2 (additional authoritative geometric TOVPoint fields) is scientifically sound but larger;
choose it only if multi-convention transport is required now. D (explicit convention tags)
requires complete tuple consistency and fail-closed mixed-input handling; tags cannot make
today's mixed background valid. A global C migration without additional scientific benefit
is outside this minimum. These alternatives are not fallback behavior silently added to A1.

### Q5 — scientific correction, revalidation and reference supersession

**Proposed answer: upon acceptance, explicitly adjudicate the current mixed tuple as internally
inconsistent and authorize only Q1–Q4's correction under the requirements of GOVERNANCE §3.1,
with U1–U14 and scoped independent review/owner ratification before reference supersession.**
This is a correction with existing historical reference artifacts, not an assertion that no
prior baseline exists. The exception substitutes independent evidence for agreement with the
rejected behavior; it never waives validation or allows a baseline to precede validation.

Required disposition after independent validation:

- Supersede `passive_cooling_cmf_1p6_debug.tsv`, both `grid_convergence_cmf_1p6` artifacts,
  `hartle_I_dscmf1_debug.tsv`, `baryon_number_dscmf1_reference.tsv`, and
  `hartle_monopole_dscmf1_debug.tsv` through their canonical producers, with repeatable bytes
  and explicit historical/new hashes. No hand editing or wider regression tolerance.
- Revalidate `tov_dscmf1_reference.tsv` and `tov_path_equivalence_dscmf1.tsv`; preserve bytes
  if, as expected under A1, their physical/relational expectations remain unchanged.
- Preserve the old Phase-4/Structure-1/preflight records as historical evidence. Recheck the
  qualified common-state Structure-1 claim without promoting its unresolved Mmax semantics.
- Revalidate static identities first, then independent Hartle/measure/count/thermal results,
  immediately establish their validated superseding baselines, and complete the coherent PB7
  method comparison before declaring the separate ADR-0011 Q4 prerequisite complete.

The audit's baseline matrix, exact chronology, independent oracle definitions, pre-existing
PB7 2e-4 target and tolerance-basis requirements are in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`. Successful Q4 prerequisite work does
not itself close INV-09; accepted ADR-0011 PB1–PB14, independent review and human ratification
remain necessary. INV-11 is not addressed.

## Proposed GOVERNANCE §3.1 record — inactive until acceptance

1. **Named invalid behavior:** GSL-solved TOV mass/pressure/energy/nup represented with
   unmatched nominal-mass and modern-G profile factors.
2. **Why old output cannot govern the correction:** preserving it as physical truth enshrines
   failure of exact Einstein/TOV identities and the independent sequence derivative.
3. **Minimum correction:** Q1–Q4 A1 only, including the required inverse/public boundaries.
4. **Independent evidence:** high-precision unit algebra; exact analytic TOV; physical-RHS
   identities; independent enthalpy TOV/phase-space EOS; conservative first order;
   `(mhat,hhat)`/Stieltjes and published/analytic controls; complete-star versus homogeneous
   particle derivative; independently justified convergence and negative controls (U1–U13).
5. **Narrow scope:** ordinary NStar unit boundary only; no adjacent physics/method cleanup.
6. **Historical output disposition:** the six affected baseline artifacts remain in history
   and are identified as superseded only after validated correction. The two unaffected
   artifacts retain their role if revalidated. Full inventory is in the audit.
7. **Immediate subsequent baseline:** U14 canonical-producer supersession follows independent
   validation of the corrected quantities, with old/new hashes and exact repeatability.

This PROPOSED record does not invoke the exception. No production modification, baseline
replacement, Phase-5B implementation, Btilde, paper Z/W, evolution or BNV is implemented here.

## Acceptance and completion boundaries

All five decisions above are pending owner ratification. No question is implicitly accepted
by committing or pushing this document. No extension of historical validation scope is implied.
The audit recommends one minimum correction; subsequent implementation and revalidation remain
separate authorized work. ADR-0011 remains ACCEPTED; INV-09 remains INTENDED BUT UNVERIFIED;
INV-11 remains UNRESOLVED.

**Exactly one recommended next action: owner review and ratification of this proposed A1
repair/revalidation contract.** Do not begin implementation automatically.
