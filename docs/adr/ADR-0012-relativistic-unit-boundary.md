# ADR-0012: relativistic unit boundary for canonical ordinary-star backgrounds

**Status: ACCEPTED — 2026-09-05**

> **ADR-0012 ACCEPTED —**
> **A1 RELATIVISTIC UNIT-BOUNDARY CORRECTION AUTHORIZED;**
> **PRODUCTION CORRECTION AND REVALIDATION NOT YET PERFORMED.**

ADR-0011 remains **ACCEPTED**. Its Q4 prerequisite remains blocked until the A1 correction and
revalidation are complete. INV-09 remains **INTENDED BUT UNVERIFIED**; INV-11 remains
**UNRESOLVED**.

## Context and authority

The canonical TOV integrator uses GSL G,c and publishes a literal mass ratio to GSL's solar
mass. NStar multiplies that ratio by Zaki's nominal solar mass length, converts pressure and
energy with Zaki's modern-G factors, and carries the original GSL nuprime. The resulting
background is internally inconsistent: the full-star mass identity has a systematic
`2.56486e-4` relative residual and the nuprime identity up to `7.33756e-5`. Coherent conversion
of the same solved state closes those algebraic identities to roundoff and reduces the leading
homogeneous/complete-star particle derivative discrepancy from approximately `9.17e-4` to
approximately `1.37e-6`. This is pre-implementation scratch evidence, not a production PB7 pass
(`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`).

Governing authority: `GOVERNANCE.md:15`, `:39`, `:65`, `:88`; accepted canonical numerical
primitive (`docs/adr/ADR-0005-canonical-tov-numerical-primitive.md:208`); geometry
(`docs/adr/ADR-0004-proper-volume-and-metric-measure-ownership.md:29`); angular and monopole
contracts (`docs/adr/ADR-0006-hartle-first-order-physical-normalization.md:243`,
`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:418`); measure/surface authority
(`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:61`,
`docs/adr/ADR-0009-tov-surface-event-and-termination.md:47`); unit prerequisite
(`docs/adr/ADR-0011-particle-number-structural-response.md:92`).

The authorized future change is **scientific-semantic**, with a bounded architecture change and
independent numerical validation. A successful round trip or old baseline is not sufficient
scientific authority for a conversion. This decision does not adjudicate a chemical convention.

## Owner decision

The human project owner accepts Q1–Q5 below as normative.

### Q1 — actual solve convention

Retain the current canonical ordinary-star TOV numerical solve convention:

```
G_TOV = GSL 2.7.1 G = 6.673e-8 cm^3 g^-1 s^-2
c_TOV = 2.99792458e10 cm/s
```

This is the numerical convention defining the existing canonical TOV calculation. It is not a
declaration that this G value is the preferred modern physical constant for all future
calculations. Changing canonical TOV to modern/CODATA G would constitute a separate physical
re-solve and requires separate authorization.

### Q2 — authoritative geometric star representation

Every ordinary `StarProfile` geometric field derived from a canonical TOV star must represent
that **same solved spacetime** using the same `G_TOV` and `c_TOV`. For CGS input and geometric km:

```
r_km = r_cm / 1e5
m_km = G_TOV * m_grams / c_TOV^2 / 1e5
epsilon_km^-2 = G_TOV * rho_cgs / c_TOV^2 * 1e10
p_km^-2 = G_TOV * p_cgs / c_TOV^4 * 1e10
nuprime_km^-1 = nuprime_cm^-1 * 1e5
```

The corrected profile must satisfy one coherent geometric TOV system:

```
exp(-2 lambda) = 1 - 2m/r
dm/dr = 4 pi r^2 epsilon
dp/dr = -(epsilon+p)(m+4 pi r^3 p)/[r(r-2m)]
nu' = (m+4 pi r^3 p)/[r(r-2m)]
dp/dr = -(epsilon+p) nu'
```

No ordinary-star `StarProfile` tuple may mix constants from distinct relativistic conventions.
The existing nu reconstruction matches to the corrected geometric surface mass. TOVSolver
equations, Geometry equations, RotationSolver equations, fixed-central-energy normalization,
the energy measure, the finite-pressure surface, and AngularVelocity c conversion do not require
modification under A1.

### Q3 — public solar-mass semantics

Preserve the existing public ordinary-TOV mass-number semantics as the literal historical ratio:

```
M_public = m_grams / GSL_CONST_CGSM_SOLAR_MASS
GSL_CONST_CGSM_SOLAR_MASS = 1.98892e33 g
```

This public mass ratio is distinct from geometric mass `m_km`. Do not reconstruct `M_public`
from `m_km / Zaki::Physics::SUN_M_KM`. Zaki `SUN_M_KM` is based on the IAU nominal solar
gravitational parameter and belongs to a different convention.

The implementation must keep distinct (1) the physical/public literal mass ratio and (2) the
geometric mass length belonging to the solved TOV spacetime. No silent conversion between them
through a mismatched solar-length constant is permitted. Preserve public central e,p in physical
CGS and copy the physical public inputs where they already exist; otherwise use the matching
inverse conversion. The existing `SeqPoint` and `StarProfile` fields suffice. No new `TOVPoint`
schema is required for the minimum A1 repair.

### Q4 — minimum migration

Accept **Option A1**. Introduce one narrowly scoped relativistic conversion owner, conceptually
`CompactStar/RelativityUnits.hpp`, for the ordinary canonical TOV-to-geometric-profile boundary.
The repair updates only the places required to make forward, inverse, and public metadata
conversions coherent:

- both ordinary NStar construction paths;
- sequence/finalization handling as required;
- NStar mass/public access semantics where necessary;
- StarContext physical-density inverse conversion;
- TbDefinition physical-density threshold conversion; and
- narrowly required unit-contract comments.

No canonical TOV numerical behavior changes: no TOV G/c, EOS, event/surface, grid, search, or
derivative-policy change. No RotationSolver equation, Geometry equation, AngularVelocity
normalization, MixedStar, BNV, or chemistry/evolution migration is authorized. Option A2 is not
selected. Option B is not selected. Option C is not selected. Option D is not selected.

### Q5 — scientific correction, revalidation, and baseline supersession

The current mixed-convention ordinary-`StarProfile` tuple is adjudicated as scientifically and
internally inconsistent. The bounded A1 repair is authorized under `GOVERNANCE.md` §3.1. The
existing affected numerical baselines are historical evidence, not authority for preserving the
known defective representation.

Required chronology:

1. accepted ADR-0012;
2. bounded A1 production correction;
3. independent unit/static identity validation;
4. qualified Structure-1 recheck;
5. independent first-order Hartle revalidation;
6. independent monopole/measure revalidation;
7. baryon-number and relevant thermal validation;
8. coherent homogeneous-vs-complete-star PB7 prerequisite check;
9. independent review;
10. human ratification; and
11. only then supersede affected governed artifacts through canonical producers.

The six artifacts expected to require validated supersession are:

- `passive_cooling_cmf_1p6_debug.tsv`;
- `grid_convergence_cmf_1p6_debug.tsv`;
- `grid_convergence_cmf_1p6_trajectory.tsv`;
- `hartle_I_dscmf1_debug.tsv`;
- `baryon_number_dscmf1_reference.tsv`; and
- `hartle_monopole_dscmf1_debug.tsv`.

The following must be revalidated and retain their bytes if the A1 expectation is correct:

- `tov_dscmf1_reference.tsv`; and
- `tov_path_equivalence_dscmf1.tsv`.

No baseline may be hand edited, and no tolerance may be widened merely to preserve old output.
Old validation records and old artifact hashes remain preserved historically. Successful
unit-boundary repair does not itself resolve INV-09.

## Unit-authority boundary

This decision chooses consistency with the numerical TOV solve for the representation of **that
solved star**. It does not claim that GSL's historical G is more physically accurate than modern
CODATA G, and it does not forbid a future intentional migration to modern constants. Such a
migration would be a different scientific change requiring a new physical re-solve and
validation.

IAU nominal `GM_sun` remains valid as a reporting standard where explicitly named, but it must
not be silently used to reconstruct the geometric mass of a star solved with another
G/`M_sun` convention.

## Structure-1 boundary

Scratch coherent-A1 results preserve the human-ratified qualified Structure-1 claim. At the
printed `rho_c = 1.10e15 g cm^-3`, the coherent representation remains approximately:

```
M_public = 0.6236355691
R_0 = 12.7681549010 km
R_infinity = 13.80244 km
```

These values remain inside the published FR2005 `0.62 / 12.77 / 13.80` rounding bins. This is
**pre-implementation scratch evidence**. Structure-1 is not already revalidated under production
A1, and source-qualified `M_max` semantics remain unresolved.

## Phase-4 and PB7 boundary

The audit found a leading-species PB7 method discrepancy of approximately `9.17e-4` under the
current mixed convention and approximately `1.37e-6` under the coherent A1 scratch convention.
This demonstrates that A1 is sufficient in principle to remove the identified unit-boundary
blocker. It is not a production PB7 pass.

First-order and monopole numerical outputs are expected to change. Historical Phase-4
mathematical validation remains evidence for the equations, but the affected numerical baselines
must be independently revalidated and superseded after the production correction. No
RotationSolver equation change is authorized or required by this decision.

## Accepted revalidation plan: U1–U14

U1–U14 are the accepted plan for the eventual correction; they are **not fourteen completed
tests**:

- U1 constant/version authentication;
- U2 cgs-to-geometric algebra and round-trip fixtures;
- U3 geometric TOV identities;
- U4 TOVPoint-to-NStar path consistency;
- U5 analytic TOV oracle;
- U6 qualified Structure-1 recheck;
- U7 TOV/path/grid regressions;
- U8 first-order Hartle independent validation;
- U9 monopole independent validation;
- U10 baryon-number revalidation;
- U11 thermal impact/conversion validation;
- U12 coherent PB7 comparison;
- U13 negative mixed-convention refusal tests; and
- U14 validated baseline supersession/provenance.

The independent oracle definitions, error-budget requirements, negative controls, file-level
scope, and baseline matrix are normative as recorded in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:1`. They are not executed by ratification.

## GOVERNANCE §3.1 authorization record

1. **Named invalid behavior:** GSL-solved TOV mass/pressure/energy/nuprime represented with
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
   validation, review, and human ratification of the corrected quantities, with old/new hashes
   and exact repeatability.

Acceptance activates this bounded authorization; it does not perform the correction or permit a
baseline to precede validation. No production modification, baseline replacement, Phase-5B
implementation, Btilde, paper Z/W, evolution, or BNV is implemented here.

## Consequences and completion boundary

The unit convention and minimum repair contract are settled. Production correction,
independent revalidation, review, owner ratification, and governed artifact supersession remain
future work. ADR-0011 remains ACCEPTED, but its Q4 implementation prerequisite remains blocked
until the A1 correction and revalidation are complete. INV-09 remains **INTENDED BUT
UNVERIFIED**; INV-11 remains **UNRESOLVED**.

**Exactly one recommended next action:** open a fresh A1 implementation/revalidation branch from
canonical master, implement the bounded RelativityUnits conversion repair, execute U1–U14,
independently revalidate affected static/rotation/count/thermal quantities, and prepare but do
not automatically ratify or supersede affected baselines. Do not begin that action automatically.

---

## Ratification / status addendum — 2026-09-06

**UNIT-1 A1 IMPLEMENTATION AND INDEPENDENT REVIEW HUMAN-RATIFIED — 2026-09-06**

This dated addendum records later execution and review without rewriting the original accepted
Q1-Q5 history above. Full owner decision and caveats:
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_RATIFICATION.md:1`.

| Item | Status |
|---|---|
| Candidate | `b3ce4f1303dbab68b68a82614c944c269cefebdc` |
| Independent review disposition | **PASS WITH NONBLOCKING FINDINGS** |
| Human owner disposition | **RATIFIED** |
| Production correction | **RATIFIED** |
| U1-U14 candidate validation | **RATIFIED** |
| Six successor candidates | **RATIFIED FOR GOVERNED SUPERSESSION** |
| Governed baseline supersession | **NOT YET PERFORMED** |
| Canonical integration | **NOT YET PERFORMED** |
| ADR-0011 Q4 | **NUMERIC UNIT PREREQUISITE RATIFIED AS SATISFIED; CANONICAL CLOSURE PENDING SUPERSESSION/INTEGRATION** |
| INV-09 | **UNRESOLVED — INTENDED BUT UNVERIFIED** |
| INV-11 | **UNRESOLVED** |

The B3c-prime and surface-coordinate test-semantics adjudications are human-ratified as part of
Unit-1. C1-C6 are authorized successor contents for the next governed supersession task; H7/H8
remain byte-identical and require no supersession. No baseline is replaced and master is not
modified by this ratification.
