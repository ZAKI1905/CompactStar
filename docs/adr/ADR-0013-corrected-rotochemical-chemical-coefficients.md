# ADR-0013 — corrected rotochemical chemical coefficients

**Status:** PROPOSED
**Decision:** PENDING OWNER ADJUDICATION — no question is accepted by this draft.
**Date:** 2026-09-06
**Starting canonical SHA:** `49ab2b8c2881b6ef7b9309307d18cea51d557f72`
**Change class:** scientific-semantic and architecture proposal; documentation only.
**Evidence companion:** `docs/validation/PHASE5C0_CORRECTED_CHEMICAL_COEFFICIENT_PREFLIGHT.md`.

## Context and controlling authority

ADR-0010 is ACCEPTED and governs cold local neutral H_x, active-species domains and the
corrected R2006 interpretation. ADR-0011 is ACCEPTED and its structural particle-number
response is canonically integrated: Phase-5B COMPLETE / GOVERNED and INV-09 VERIFIED /
RESOLVED for ordinary NStar. INV-11 remains UNRESOLVED. This proposal concerns only the
chemical-coefficient layer between those existing capabilities and future secular evolution
(`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:357`;
`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:100`;
`docs/SCIENTIFIC_INVARIANTS.md:946`).

The governing primary source is
`/Users/keeper/Documents/CompactStar/literature/rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`,
SHA-256 `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8`,
journal569–571 / PDF2–4, especially (10)–(19). F2005 supplies the retained non-superfluid
framework and closing benchmark where not superseded; source hashes/pages, supersession and
scope are in the companion record, section2. All literature remains read-only.

This ADR does not amend the local EOS contract, structural PN1–PN8, or evolved state.
It proposes the missing ownership and scientific semantics, not production implementation.
The governing report/proposal boundary is `GOVERNANCE.md:64`; acceptance and validation are
separate (`GOVERNANCE.md:159`).

## Source mathematics proposed for adoption

In species order `(n,p,e,mu)`, q=(0,1,-1,-1). R2006 has
`mu^infinity=(mu+q psi)e^Phi`. In the chemical Cowling approximation, with intrinsic
`chi=partial n/partial mu` and `u=e^-Phi delta mu^infinity`, neutrality implies

```text
delta psi=(q^T chi u)/(q^T chi q),
C_projected=chi-chi q(q^T chi q)^-1 q^T chi.
```

The printed (11) appears to omit e^-Phi when its numerator uses delta-mu-infinity.
Classification is **INFERRED PRINTED/SOURCE OMISSION**; no authenticated erratum is asserted.
The inference follows by substitution into (10) and agrees with (13). Beta-reaction charge
cancellation is algebraic; it never licenses omitting electrostatics from the unconstrained
density response. `C_projected q=0`, and the full integrated corrected matrix is singular.
Its inverse is forbidden (companion sections3–5).

On the accepted neutral chart,

```text
x=(n_B,n_e,n_mu), g_x=(mu_n,-eta_npe,-eta_npmu),
y=(n_n,n_e,n_mu), g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu),
T=[[1,-1,-1],[0,1,0],[0,0,1]], y=T x,
H_x=T^T H_y T,
C_y=T H_x^-1 T^T=H_y^-1.
```

For a valid intrinsic extension, constrained quadratic minimization proves
`C_projected=S_x H_x^-1 S_x^T=S_y C_y S_y^T`, with the explicit species lifts in companion
section4. The neutral provider already incorporates the correction. No public delta-psi,
individual charged-potential reconstruction, or second projection is proposed.

With `nu=Phi`, r,m in km, and local C in fm^-3 MeV^-1, define

```text
G_y = 10^54 integral_D 4 pi r^2 e^-nu C_y /sqrt(1-2m/r) dr,
delta N_y = G_y delta g_y^infinity.
```

G_y is the reduced source Btilde in count/MeV; proper volume and one inverse lapse occur
exactly once. D is a declared connected diffusive reservoir, not an implicit reaction mask.
Use the canonical Geometry measure owner; chemistry owns its inverse lapse and count
conversion (companion sections6,9; `CompactStar/Geometry.hpp:1`).

For a globally baryon-conserving reservoir,

```text
b=(1,1,1), b^T delta N_y=0,
L=[[-1,-1],[1,0],[0,1]], delta N_y=L delta N_l,
eta^infinity=-L^T delta g_y^infinity,
Z=L^T G_y^-1 L,
eta^infinity=-Z delta N_l.
```

In channel order `(npe,np-mu)`, the canonical matrix is
`Z=[[Z_npe,Z_np],[Z_np,Z_npmu]]` in MeV/count. It is symmetric positive definite where G_y
has supported stable modes. Equivalently, after integration transform to G_x and take its
baryon Schur complement Q, then Z=Q^-1. **No local fixed-n_B reduction** is allowed; it
changes the global physics. Global channel reduction must follow actual mode support, not
fabricated absent-species rank (companion sections7,10).

For the whole-star structural mapping,

```text
W=Z (I_phys,e,I_phys,mu)^T                 [MeV s^2],
dot N_l^eq=+2 I_phys,l Omega Omega_dot,
delta N_l=N_l-N_l^eq,
dot eta^infinity=-Z R_l+2 W Omega Omega_dot  (frozen coefficients).
```

The last relation records the source sign, not an implemented evolution equation. R_l is
positive for net beta decay producing lepton l. I_phys is consumed through the existing
governed structural object; the chemical consumer never repeats the c^-2 conversion. The
whole-star versus core distinction and boundary exchange must be preserved (companion
sections8–9; `CompactStar/Analysis/src/ParticleNumberResponse.cpp:511`).

## Pending owner questions

Each recommendation is pending. Alternatives are genuine ownership/domain choices where
compatible with the source; source-inconsistent full inverses or missing corrections are
not presented as acceptable owner alternatives.

| ID | Owner decision | Recommendation | Alternative / consequence |
|---|---|---|---|
| Q1 | Canonical reduced global basis and object | Source y basis, `Analysis::GlobalChemicalNumberResponse`, with named axes and count/MeV G_y. | Canonical G_x with an exact y view is equivalent but adds source-mapping burden. A bare paper B class or full corrected inverse is not valid. |
| Q2 | Baryon-reduction owner and supported-space policy | Stellar chemical layer reduces **after** integration; retains G_y and a separate Z result; explicit smaller supported channel set if required. | A single combined immutable G/Z result can own both stages, provided the global order and support proof remain explicit. Pointwise baryon closure is invalid. |
| Q3 | Public chemical Z representation | One symmetric named matrix in `Analysis::ChemicalImbalanceResponse`; scalar paper accessors are derived views of that owner. | Three canonical named scalars with a derived matrix view are equivalent if there is still one authority and named orientation. Recommendation favors matrix action/provenance. |
| Q4 | W ownership and structural coupling | Separate immutable `Analysis::RotochemicalSpinDrive`, dependent on Z and complete governed `FixedBaryonNumberResponse`; whole-star path calls `WholeStarIPhysical()`. | Combined coefficient bundle is possible but forces structural invalidation on purely chemical consumers. Baseline JSON runtime coupling is forbidden. |
| Q5 | Chemical domains, active embeddings, thresholds and interfaces | First bounded whole-star Track-R coefficient implementation, with explicit 1D/2D/3D embeddings and certified finite-cut/onset handling; refuse ungoverned gaps/interfaces/core mapping. | Begin with an authenticated core model only after its reservoir, cutoff, boundary-flux and phase authority exist. No invented free-gas core cutoff or source-core claim. |
| Q6 | Correction-sensitive Track-R closure gate | Preserve A18+delta-v+UIX* closing EOS, require corrected Figure1/author coefficient comparison with authenticated source product and extraction uncertainty. | Authenticated correction-sensitive transient may supplement or replace that comparison if full evolution and source configurations exist. Quasi-steady temperature alone cannot qualify. |
| Q7 | Redshift facts versus evolved-state decisions | Fix coefficient eta-infinity semantics and e^-nu measure now; defer actual state layout, ordering/storage integration and INV-11 ratification. | Defer all chemical implementation until a broader evolution ADR, at the cost of coupling two otherwise separable contracts. No convenience-driven ChemState assignment. |
| Q8 | Provenance and numerical refusal | Immutable, dependency-complete results with RequireCurrent-style refusal and condition/error-enclosure gates; numerical cutoffs only after predeclared convergence/threshold evidence. | Explicit recomputation on every access avoids cached stale data but still needs full dependency/domain authentication and fails on scientifically unresolved inputs. |

## Recommended architecture and provenance details

Local provider ownership remains unchanged. Its proposed local adapter is
`CompactStar::ChargeNeutralNumberSusceptibility`, which performs only qualified active H
solves and basis embedding. Stellar G integration, global baryon reduction, and Z are in the
Analysis coefficient layer. Phase-5B owns structural response. Future evolution consumes
results and never recomputes local thermodynamics or structure implicitly (companion section16).

Recommended channel names are `BetaChannel::Npe` and `BetaChannel::NpMu`; G axes are
`NeutralNumberCoordinate::{Neutron,Electron,Muon}`. Proposed accessors are
`ResponseCountPerMeV(row,column)`, `ResponseMeVPerCount(output_channel,input_lepton)`,
`DriveMeVSecondsSquared(channel)`, and paper views `PaperZnpe()`, `PaperZnp()`,
`PaperZnpMu()`. They are proposed semantic names only, not current callable symbols.
No unqualified B/Z/W/MatrixB object or duplicate scalar/matrix authority is recommended.

Provenance includes profile identity/version, metric normalization, provider identity/revision
and bytes/constants, equilibrium match, source/basis/units, active branch support, reservoir,
partitions, thresholds/interfaces, surface/tail policy, numerical conditioning and uncertainties.
W adds every contributing structural source and its domain/spin normalization. RequireCurrent
must reject stale sequence profiles and EOS bytes as well as the main star. Lifetimes must
be owned or token-validated before dereference. A copied I array is not adequate runtime
provenance (`CompactStar/Analysis/src/ParticleNumberResponse.cpp:247`; companion section16).

Active embeddings are full T, npe `[[1,-1],[0,1],[0,0]]`, pe `[[0],[1],[0]]`; vacuum and
threshold objects never provide a fabricated H. Continuous onset limits are checked on C,
not divergent H. Positive-width provider refusal windows require certified numerical handling
or refusal. True phase/interface response requires separate authenticated metadata and a
chemical-interface law; structural jump terms do not automatically supply it (companion section10).

## Proposed validation and ratification gate

The companion section17 predeclares **GC1–GC14**: units/signs, neutral reconstruction, H
stability, analytic inverse, basis equivalence, independent intrinsic projection, charge-null
and no-second-projection control, active embeddings, independent GR integration, global baryon
reduction, source Z, structural W, correction-sensitive published comparison, and provenance/
domain/mutation coverage. Each specifies oracle, independence boundary, defect, metric,
tolerance source and negative controls. Section18 lists **M1–M18** with explicit detectors.

The exact coupled toy agrees by rational arithmetic and has the expected singular full
matrix. Whole-star free-gas diagnostic matrices and independent local/GR checks are recorded
in companion sections11–13. **These are not Z/W validation, production GC passes, or golden
targets.** Production numerical tolerances remain unratified until onset-aware convergence,
refusal-gap, provider/profile matching, conditioning and source-error evidence exists.

Two important detector qualifications are binding parts of the proposed validation claims:
repeating the identical charge projector is idempotent, and transposing symmetric Z changes
nothing. Numerical comparisons cannot detect those no-ops. Typed orientation/call-path audits
must cover the forbidden operation/semantic mixup; do not inflate numerical mutation counts
with equivalent mutants. Structural PB9/PB caveats are retained, not converted into new
independent chemistry evidence (companion sections5,17–18).

A18 arbitrary-composition, crust/core/phase authority and quantitative corrected benchmark
remain missing. R06 Figure2's high-mass panel says2.14 M_sun while its caption says2.13;
that benchmark choice needs source adjudication. Figure1 permits a coefficient benchmark
without that mass ambiguity, once its EOS and extraction are authenticated (companion
sections14–15). Missing realistic closure inputs do not authorize another EOS.

## Consequences, non-goals, and exact boundary

Accepting this ADR later would govern coefficient semantics and ownership; it would not
itself certify numbers or implement the layer. Independent scientific review and owner
ratification precede any separately authorized production implementation. An implementation
must predeclare its numerical/error policy within the accepted contract and demonstrate the
relevant GC gates; the realistic closing benchmark remains a separate required Track-R gate.

No production Btilde, Z, W, eta state/evolution, Urca rates, heating/cooling, superfluidity,
BNV, EOS product, regression baseline, or source-data modification is authorized by this draft.
There is no production or test code in this change and no scientific runtime output changes.

INV-11 remains UNRESOLVED: evolved channel/storage representation, redshifted state units,
reaction stoichiometry/net-rate sign, thermal/neutrino bookkeeping, coefficient variation and
solver coupling require later governance. Phase6 must permit generic externally supplied
particle-number sources separately from spin driving, retaining the unreduced neutral response
and explicit baryon assumptions; baryon-changing sources cannot be forced into a fixed-baryon
2-channel lift. No BNV equations or thermal assumptions are introduced (companion section19).

**ADR-0013 remains PROPOSED. No owner answer, merge, ratification, or production implementation
has occurred.**
