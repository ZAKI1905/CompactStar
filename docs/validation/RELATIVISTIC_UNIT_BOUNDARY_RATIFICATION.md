# ADR-0012 A1 relativistic-unit boundary — human ratification

**HUMAN-RATIFIED — 2026-09-06**

**THIS RATIFICATION DOES NOT ITSELF MODIFY GOVERNED BASELINES OR MASTER.**

The human project owner ratifies the ordinary-star ADR-0012 A1 correction at
`b3ce4f1303dbab68b68a82614c944c269cefebdc` (`fix: reconcile ordinary-star
relativistic units`). The independent review disposition is:

> **UNIT-1 INDEPENDENT REVIEW PASS WITH NONBLOCKING FINDINGS —
> READY FOR HUMAN RATIFICATION WITH EXPLICIT CAVEATS**

Owner disposition: **RATIFIED**. The review reported no BLOCKING findings, no MATERIAL
findings, three NONBLOCKING findings, and one NOTE. This decision accepts A1 and its
validated successor package. It does not install those successors, merge the candidate,
or authorize any downstream physics. ADR-0012's accepted decision and A1 scope are recorded
in `docs/adr/ADR-0012-relativistic-unit-boundary.md:38`; the completed candidate evidence is
recorded in `docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:68`.

## 1. Authenticated subject and review

- Canonical parent: `1d22dd1f5a0d1afa18c4cedebb36b28fdae49df4`.
- Candidate branch: `fix/relativistic-unit-boundary-a1`.
- Candidate worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-relativity-units-a1`.
- Ratified Unit-1 commit (`UNIT1_SHA`): `b3ce4f1303dbab68b68a82614c944c269cefebdc`.
- Independent review worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-unit1-a1-review`.
- Independent review disposition: **PASS WITH NONBLOCKING FINDINGS**.
- Human owner disposition: **RATIFIED**.

At ratification entry, the canonical checkout was clean; local, origin-tracking, and live
remote `master` all named the canonical parent. Local, upstream, and live remote candidate
all named `UNIT1_SHA`; its parent was exactly canonical master and it was exactly one commit
ahead and zero behind. All eight governed baselines retained their historical hashes and the
candidate baseline diff was empty. Candidate preservation and baseline evidence are recorded
in `docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:275` and
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:424`.

The reviewer independently reproduced the complete 53-test suite and the 30-test data-free
inventory, independently checked the unit constants and wrong-convention detector margins,
and independently regenerated C5 and C6 twice. The two C5 regenerations both hashed to
`90d607519cbdf3c4a0bf6ef50cc8fd22a8526b5db0354dc319e96854da29041d`; the two C6
regenerations both hashed to
`caaa0ac0d3219cda0a9fb518b27688afc23c6cdad1ec76a2bcd7359614a8d4e8`.

## 2. Ratified scientific claim

The canonical ordinary-star TOV numerical solution itself is unchanged. Before A1, the
`NStar`/`StarProfile` publication boundary mixed GSL-solved quantities with Zaki/IAU-derived
geometric conversion factors and therefore did not represent one self-consistent geometric
spacetime. The original audit demonstrates that mismatch and derives the coherent mapping in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:290` and
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:254`.

ADR-0012 A1 corrects only that representation boundary:

- every ordinary-`StarProfile` geometric quantity uses the same `G_TOV` and `c_TOV` as the
  canonical TOV solve;
- public mass remains the literal GSL solar-mass ratio;
- geometric mass is a separate kilometre-valued representation;
- forward and inverse physical/geometric boundaries are coherent;
- the TOVSolver, Hartle, Geometry, and EOS equations are unchanged.

Independent review confirmed the load-bearing mass-accumulation, hydrostatic, `nuprime`, and
pressure/`nuprime` identities close near machine precision after A1; the rejected mixed
convention fails at approximately `1e-4` scale. The candidate identity results and negative
controls are recorded in `docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:91` and
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:263`.

## 3. Ratified U1-U14 status

| Gate | Ratified disposition |
|---|---|
| U1 | PASS |
| U2 | PASS |
| U3 | PASS |
| U4 | PASS |
| U5 | PASS |
| U6 | PASS |
| U7 | PASS as candidate/revalidation |
| U8 | PASS |
| U9 | PASS |
| U10 | PASS |
| U11 | PASS |
| U12 | Numeric unit prerequisite satisfied in candidate/review |
| U13 | PASS |
| U14 | Validated successor package complete |

The detailed evidence and scope of every gate remain those in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:68`; ratification does not
broaden them. In particular, U12 does **not** mean PB7 is complete. ADR-0011 Q4's numeric-unit
prerequisite is human-ratified as satisfied by A1, but Q4 remains open until the subsequent
governed baseline-supersession and canonical-integration step establishes the accepted
canonical state. INV-09 remains **INTENDED BUT UNVERIFIED / unresolved** and INV-11 remains
**UNRESOLVED**. ADR-0011's Q4/PB7 boundary is retained in
`docs/adr/ADR-0011-particle-number-structural-response.md:92`.

## 4. Qualified Structure-1 result

The qualified FR2005 Structure-1 common-state result survives A1. At
`rho_c = 1.10e15 g cm^-3`, the ratified values remain approximately:

| Quantity | Value |
|---|---:|
| `M_public` | `0.6236355691` |
| `R_0` | `12.7681549010 km` |
| `R_infinity` | `13.8024398765 km` |

These map to the source's printed `0.62 / 12.77 / 13.80` bins. This is the already-qualified
common-state claim, not a source-qualified maximum-mass claim. `M_max` remains unresolved.
The A1 Structure-1 comparison is recorded in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:491`; the binding qualification remains
the Structure-1 ratification referenced by
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:138`.

## 5. Ratified test-semantics adjudications

### B3c

The historical `worst/best < 3` residual-ratio criterion characterized the old systematic
offset and was superseded by ADR-0012. The accepted B3c-prime criterion is retained:

```text
res_k <= 0.5 s_k
```

This is not a tolerance loosening: it rejects the old defective data and passes the corrected
data. The derivation, pre-A1 falsification, and A1 results are in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_B3C_ADJUDICATION.md:15` and
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_B3C_ADJUDICATION.md:45`.

### Surface coordinate

The original combined assertion was:

```cpp
Require(r == p.back().r && m == p.back().m * relativity_fixture::solar_km,
        "surface profile coordinate");
```

Exact surface-radius equality remains valid because the radius is copied directly. Exact
floating-point mass bit equality is not a physical or ADR invariant across independently
ordered arithmetic graphs. The ratified replacement retains two independently localizing
checks:

```cpp
Require(r == p.back().r,
        "surface profile radius");

const double m_expect =
    p.back().m * relativity_fixture::solar_km;

Require(std::fabs(m / m_expect - 1.0) <=
            16.0 * std::numeric_limits<double>::epsilon(),
        "surface profile geometric mass");
```

The accepted relative gate is exactly `16 * std::numeric_limits<double>::epsilon()`. Two
independently ordered expressions each contain at most about five rounded IEEE-754 operations,
giving a first-order separation estimate of approximately `10u`, or about `5 epsilon`; the
next conservative power-of-two bound is `16 epsilon`. The observed worst four-star difference
was one ULP and was not used to select the bound. Wrong mass conventions and unit powers remain
many orders of magnitude outside the gate. The high-precision oracle, four-star result,
operation-count derivation, and mutation table are in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_SURFACE_ASSERTION_ADJUDICATION.md:84`,
`:101`, `:120`, and `:143`.

## 6. Ratified successor package

The following exact contents are authorized for the **next** governed baseline-supersession
task. They are not installed by this ratification.

| ID | Artifact | Authorized successor SHA-256 |
|---|---|---|
| C1 | passive cooling | `8fef2314673fceb939f859612f4befe94117115d6d6b3ad0dcc59d1faa68c9f9` |
| C2 | grid debug | `b48519c3e948e9979a385d19facee2777d15955eeb8711b4bdd46b81fef74741` |
| C3 | grid trajectory | `d5b753932c0523e67a7f25b460c7494bec1a006a8d01c9e43124cb2e78f0720f` |
| C4 | Hartle I | `034ecddbd9bd847650429d7dc87d0331ec9e87aca3862ff87594e4bff5b707dd` |
| C5 | baryon number | `90d607519cbdf3c4a0bf6ef50cc8fd22a8526b5db0354dc319e96854da29041d` |
| C6 | Hartle monopole | `caaa0ac0d3219cda0a9fb518b27688afc23c6cdad1ec76a2bcd7359614a8d4e8` |

The package's producer, repeatability, validation, and numerical-delta evidence is recorded in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_CANDIDATES.json:1` and summarized in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:275`.

H7 and H8 require no supersession:

| ID | Artifact | Retained SHA-256 | Disposition |
|---|---|---|---|
| H7 | TOV reference | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` | Revalidated, byte-identical; retain |
| H8 | TOV path equivalence | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` | Revalidated, byte-identical; retain |

## 7. Full-suite review result

The independent review reproduced the complete suite:

- 53 total;
- 50 raw PASS;
- three historical-baseline mismatches;
- Category A = 50;
- Category B = 3;
- Category C = 0.

The three Category-B results are solely the unchanged historical comparisons for passive
cooling, Hartle monopole, and baryon number; each emitted output exactly matches its ratified
successor. The candidate's complete inventory and byte-authentication table are in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:297`.

## 8. Durable nonblocking findings

### NONBLOCKING-1 — MixedStar divergence

ADR-0012 explicitly excluded `MixedStar`. Ordinary `NStar` is coherent under the accepted GSL
relativistic convention, while `MixedStar` retains historical conversion behavior. This does
not block A1 ratification. Future roadmap debt: **audit and, if scientifically authorized,
reconcile the MixedStar relativistic unit boundary in a separate task.** Nothing in this
ratification implements or authorizes that task. The original exclusion is recorded in
`docs/adr/ADR-0012-relativistic-unit-boundary.md:121`.

### NONBLOCKING-2 — ADR-0012 historical status text

The original ADR header correctly remains historical and says the production correction and
revalidation had not yet been performed. It is not rewritten. A dated ratification/status
addendum now records the completed implementation, independent review, and human decision.

### NONBLOCKING-3 — focused-selection provenance

The historical Unit-1 focused 42-test selection was not fully reconstructible from the
committed tree because its selection list lived in build evidence. This does not undermine
ratification: the independent reviewer reproduced the complete 53-test suite, reproduced the
30-test data-free inventory, matched all 53 reported classifications, and performed separate
numerical checks. This remains validation-provenance debt. Future major validation campaigns
must make named focused selections durable in tracked test metadata or documentation. No test
infrastructure is added or redesigned here.

## 9. Note — tautological metric identity

U3's `exp(-2 lambda) = 1 - 2m/r` is structurally tautological because `lambda` is constructed
from `m` and `r`. It remains useful as an internal representation check but is not an
independent falsifier. The load-bearing independent identities are mass accumulation, the
hydrostatic equation, the `nuprime` equation, and the pressure/`nuprime` relation. They close
at approximately `1e-15` after A1 and fail strongly under the rejected mixed convention;
the identity roles were distinguished before implementation in
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:273` and
`docs/validation/RELATIVISTIC_UNIT_BOUNDARY_AUDIT.md:283`.

## 10. Exact remaining boundaries

- Governed baseline supersession: **NOT YET PERFORMED**.
- Canonical integration: **NOT YET PERFORMED**.
- ADR-0011 Q4: numeric unit prerequisite **RATIFIED AS SATISFIED**; canonical closure pending
  governed supersession/integration.
- PB7: **INCOMPLETE / NOT CLOSED**.
- INV-09: **INTENDED BUT UNVERIFIED / unresolved**.
- INV-11: **UNRESOLVED**.
- FR2005 source-qualified `M_max`: **UNRESOLVED**.
- Phase-5B, Btilde, paper Z/W, evolution, and BNV: **NOT BEGUN**.
- `MixedStar` reconciliation: separate future audit; **NOT BEGUN**.

**THIS RATIFICATION DOES NOT ITSELF MODIFY GOVERNED BASELINES OR MASTER.**

## 11. Exactly one recommended next action

Perform a separate governed baseline-supersession and canonical-integration task: regenerate
and install only C1-C6 through their canonical producers, preserve H7/H8, rerun the complete
53-test suite expecting 53/53 PASS after supersession, independently authenticate all new
governed hashes, then fast-forward master through `UNIT1_SHA` and `RATIFICATION_SHA` only if
every post-supersession check passes. Do not begin that action automatically.
