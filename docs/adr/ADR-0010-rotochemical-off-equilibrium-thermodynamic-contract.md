# ADR-0010 — Rotochemical off-equilibrium thermodynamic contract

| Field | Value |
|---|---|
| **Status** | **ACCEPTED** |
| **Date** | 2026-09-04 |
| **Proposal starting SHA** | `77a328676d83f515fe603cb62341d2efcd70ed78` |
| **Acceptance starting SHA** | `e4d1b8f4b30d10760bc0c894cee959040db6f9a3` |
| **Change class** | scientific-semantic and structural/architecture decision; documentation-only ratification in Phase 5A-1B |
| **Governing authority** | `GOVERNANCE.md`; ADR-0001; ADR-0002; ADR-0007; ADR-0008; ADR-0009; INV-01, INV-09, INV-11; `docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md` |
| **Evidence companions** | `docs/validation/PHASE5A1_THERMODYNAMIC_CONTRACT_DERIVATION.md`; `docs/validation/PHASE5A1A_ADR0010_INDEPENDENT_ADJUDICATION.md` |
| **Implementation state** | no thermodynamic EOS, derivative, interpolation, matrix, coefficient, rate, or evolution implementation is authorized or added |

## Context

Phase 4 is complete and its existing `HartleFirstOrderResponse` and
`HartleMonopoleResponse` objects are the ratified structural inputs for Phase 5
(`docs/validation/PHASE4_CLOSEOUT.md:64-176`). Phase 5A-0 established that the current
canonical cold CompactStar EOS path is a one-dimensional beta-equilibrium trajectory and cannot
provide the independent-composition response required by the corrected rotochemical formalism
(`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:10-31`,
`:293-324`, `:511-528`).

The architectural choice is therefore not whether to expose one more column from the existing
table. It is what local, off-equilibrium thermodynamic contract can support both a published
literature reproduction and future CompactStar research EOSs without making the EOS own stellar
integration or rotochemical evolution.

Phase 5A-1 drafted this ADR as **PROPOSED**. The Phase 5A-1A independent Astra adjudication then
checked the source mapping and the neutral-Hessian construction and required revisions R1-R7
(`docs/validation/PHASE5A1A_ADR0010_INDEPENDENT_ADJUDICATION.md:7-19`, `:468-614`). In Phase
5A-1B the project owner ratified the six decisions recorded below with those revisions. This
acceptance governs the contract; it does not certify an implementation or reproduce any
rotochemical result.

## Source authority

The scientific source root is `/Users/keeper/Documents/CompactStar/literature`. The exact sources
used are:

1. `/Users/keeper/Documents/CompactStar/literature/rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`
   — primary authority where electrostatic potential, local charge neutrality, the singular full
   response, or the reduced inverse is at issue: eqs. (2), (4), and (9)-(19), sections 2-5,
   PDF pp. 2-4 / journal pp. 569-571.
2. `/Users/keeper/Documents/CompactStar/literature/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`
   — primary non-superfluid reproduction benchmark: eqs. (7)-(8), PDF p. 3, and eqs. (9)-(13),
   PDF p. 4; eq. (24), PDF p. 6;
   eq. (30), PDF p. 7; section 3.1, PDF pp. 8-9; eq. (50), PDF p. 12; and eqs. (51)-(58),
   PDF p. 13. Its eqs. (9), (12), (13), and (54)-(56) are superseded as specified by 2006
   section 5.
3. `/Users/keeper/Documents/CompactStar/literature/rotochemical/1995-Reisenegger-Rotochemical-Heating.pdf`
   — supporting origin and coordinate-count example: section 2, paper pp. 4-5 / PDF pp. 14-15,
   especially eq. (3).
4. `/Users/keeper/Documents/CompactStar/literature/rotochemical/2020-Yanagi-NS-Therm-Thesis.pdf`
   — supporting review: section 4.1.1, thesis pp. 58-60 / PDF pp. 60-62, especially eq. (4.7).

Where the 2005 and 2006 treatments conflict, the 2006 paper governs. The gravitochemical paper is
not authority for this decision. The detailed equation ledger and the derivation from the source
variables to the candidate software bases are in the evidence companion §§2 and 5-6.

## Problem

For cold npe-mu matter, the local species densities obey

```text
n_B = n_n + n_p,
n_p = n_e + n_mu.
```

The first equation defines baryon density; the second is the physical local charge constraint of
2006 eq. (4). Thus four species densities have three independent local degrees of freedom. Global
baryon conservation is applied only after stellar integration (2006 eqs. (14)-(15) and footnote
4); it must not be misapplied as a pointwise condition that removes the local `n_B` direction.

Beta equilibrium then selects a one-dimensional curve within the three-dimensional local state
space:

```text
eta_npe  = mu_n - mu_p - mu_e  = 0,
eta_npmu = mu_n - mu_p - mu_mu = 0
```

where the relevant species are present. The current equilibrium curve contains only one
directional combination of a three-dimensional response and cannot identify the two transverse
composition directions (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:293-324`).

The electrostatic potential is indispensable in the rigorous derivation. Total redshifted
chemical potentials contain `q_i psi` (2006 eq. (2)); perturbations contain `q_i delta psi`
(eq. (10)); and local neutrality determines `delta psi` (eq. (11)). It cancels from the beta
imbalances only because each beta channel is electrically neutral (eq. (9)). Dropping it before
the projection reproduces the 2005 error.

## Scientific scope

The owner-specified Phase-5 order is locked:

1. reproduce standard non-superfluid rotochemical heating;
2. reproduce or check the corrected 2006 electrostatic formalism;
3. only then consider superfluidity or other refinements; and
4. keep BNV in Phase 6, blocked until the standard rotochemical benchmark is reproduced.

This accepted decision governs only the local thermodynamic contract and its interfaces to global stellar
susceptibility integration and the verified Phase-4 structural response. It does not implement
any of those layers.

## Track R vs Track P

| Track | Scientific purpose | Accepted entry and gate |
|---|---|---|
| **R — literature reproduction** | Reproduce the standard 2005 non-superfluid calculation with the 2006 correction where it governs. | Use analytic lepton/toy checks, the 2005 whole-star free-gas case, optional authenticated BPAL, and the A18 + delta-v + UIX* closing benchmark. Authenticate the exact source definition, phase/crust prescription, and comparison data before closure. Apply the 2006 correction from the outset and include a correction-sensitive coefficient or transient comparison; quasi-steady temperature alone cannot validate it. A missing closing source blocks closure; another EOS does not silently substitute. |
| **P — CompactStar production** | Apply already-validated machinery to DS(CMF)-1 and later research EOS families. | Canonical DS(CMF) identity requires the same authenticated model revision, parameter set, species, phases, mass/lepton conventions, and matched equilibrium background. Supply the neutral composition response through a matched arbitrary-composition product/provider; a two-dimensional hadronic surface plus validated source-compatible leptons is allowed. The current equilibrium file supplies the background only. Restrict v1 applicability to a verified npe-mu domain; additional responding species require separate governance. |

The accepted architecture is one generic consumer-facing contract with model-specific providers.
It does not require Track R and Track P to share an EOS model or underlying algorithm. It requires
them to share state, unit, derivative, constraint, provenance, and failure semantics.

## Required local state

### Physical state count

The charge-neutral density manifold has dimension three. A provider must expose one explicit
three-coordinate chart, its domain, and exact reconstruction of all four species densities. A
muon threshold or phase boundary is a domain boundary to report, not a reason to silently change
the meaning or dimension of returned derivatives.

At a smooth allowed state, the intrinsic zero-temperature energy differential is

```text
d epsilon = mu_n d n_n + mu_p d n_p + mu_e d n_e + mu_mu d n_mu.
```

On the candidate `(n_B,n_e,n_mu)` chart this becomes

```text
d epsilon = mu_n d n_B - eta_npe d n_e - eta_npmu d n_mu.
```

This follows by substituting the two physical relations above; it is not a new dynamical
assumption. It identifies the thermodynamically conjugate vector as
`(mu_n,-eta_npe,-eta_npmu)`.

## Candidate coordinate bases

| Candidate | Completeness and neutrality | Conditioning and units | Relation to corrected 2006 system | Provider compatibility | Assessment |
|---|---|---|---|---|---|
| **A. `(n_B,n_e,n_mu)`** | Complete three-density chart; `n_p=n_e+n_mu`, `n_n=n_B-n_e-n_mu` enforce neutrality exactly. | All inputs are `fm^-3`; no division by `n_B`; positivity domain is explicit. Species scales still differ and must be reported by numerical implementations. | Linearly equivalent to the source-reduced basis. Its conjugates are `(mu_n,-eta_npe,-eta_npmu)`, so observable imbalances are direct. | Natural `n_B` axis for analytic, tabulated, and RMF providers; requires two independent lepton-density directions. | **ACCEPTED canonical basis.** |
| **B. `(n_n,n_e,n_mu)`** | Complete; `n_p=n_e+n_mu`, `n_B=n_n+n_e+n_mu`. | All inputs are `fm^-3`; linear constraints. | Directly conjugate to the 2006 variables `(mu_n,mu_p+mu_e,mu_p+mu_mu)` used before eq. (14). | Natural for the paper and energy Hessians; less directly aligned with standard `n_B`-axis tables. | Equally physical; mandatory exact transform target for validation. |
| **C. `(n_B,Y_e,Y_mu)`** | Complete for `n_B>0`; `n_i=n_BY_i`. | Mixed units and a nonlinear transform; degenerates as `n_B` approaches zero and scales composition derivatives by `n_B`. | Equivalent away from `n_B=0`, but not the source's linear density chart. | Convenient table/query representation. | Allow as an adapter or data-storage chart, not recommended as the canonical derivative chart. |
| **D. full `(n_n,n_p,n_e,n_mu)`** | Contains an unphysical charged direction unless accompanied by a constraint/projection. | Uniform density units but redundant on the physical manifold. | Needed only as an intermediate reconstruction of 2006 eqs. (11)-(13). | Some microscopic EOSs may naturally provide it. | Not recommended as the public Phase-5 primitive; must never imply that the corrected full matrix is invertible. |

For A and B, define

```text
(n_n,n_e,n_mu)^T = T (n_B,n_e,n_mu)^T,
T = [[1,-1,-1],[0,1,0],[0,0,1]],

(n_n,n_p,n_e,n_mu)^T = S_x (n_B,n_e,n_mu)^T,
S_x = [[1,-1,-1],[0,1,1],[0,1,0],[0,0,1]].
```

The determinant is one, and the inverse is exact. Any conforming provider or adapter must produce
equivalent physical responses under this transformation. The owner accepts A as the canonical
software chart and requires exact response equivalence to B on the declared smooth npe-mu domain.
Thresholds or missing species must not silently change the coordinate meaning.

## Candidate derivative authorities

Let `K_ij = partial mu_i / partial n_j` and
`chi_ij = partial n_i / partial mu_j`, with the independent variables and smooth phase domain
made explicit.

| Candidate authority | Relation to sources | Singularity and symmetry | Numerical / architectural consequence | Assessment |
|---|---|---|---|---|
| **Full intrinsic `K` over four species** | Inverse orientation of 2005 eq. (10), before the 2006 charge projection. | A thermodynamic Hessian is symmetric where smooth; it is not itself the corrected paper `B`. | Requires an unphysical charged direction, a local solve/inverse, the 2006 projection, and later volume integration. | Broader than Phase 5 requires; easy to misuse. Not recommended as the minimum contract. |
| **Full intrinsic `chi` over four species** | Direct local object in 2005 eq. (10). | Symmetric by the identity cited after 2006 eq. (13). After applying eq. (11), the charge-projected response has a null charged/gauge mode. | Most literal source input, but requires non-neutral response data and invites forbidden inversion of the corrected full global matrix. | Permissible provider capability, not the accepted public primitive. |
| **Reduced constrained chemical Hessian** | Restriction of intrinsic thermodynamics to the local charge-neutral manifold; exactly transformable to the 2006 reduced variables. | Symmetric in conjugate coordinates where the neutral energy is `C^2`. Invertibility must be established on the declared active-species domain; strict stability is sufficient. | The EOS returns the most local physical derivative. A later local adapter obtains the reduced number susceptibility; the stellar layer integrates it. | **ACCEPTED derivative authority.** |
| **Reduced constrained number susceptibility** | Direct local integrand orientation for the reduced 2006 `Btilde`. | Symmetric inverse response only where the corresponding reduced Hessian is nonsingular; no smooth full-rank continuation through absent-species or phase boundaries is implied. | Simplifies the stellar integrand but makes every provider own the inverse-orientation result. | Derived by the adapter from the accepted Hessian; not the EOS-facing authority. |

The accepted semantic primitive is
`ChargeNeutralChemicalHessian`:

```text
H_ab = partial g_a / partial x_b,
x = (n_B, n_e, n_mu),
g = (mu_n, -eta_npe, -eta_npmu).
```

Its units are `MeV fm^3`. The inverse response, if later formed, has units
`fm^-3 MeV^-1`. This decision specifies neither how a provider obtains `H` nor how a later layer
solves for its inverse response.

The source-basis equivalent uses
`y=(n_n,n_e,n_mu)` and
`g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu)`. If `H_y` is its reduced Hessian, then
`H_x=T^T H_y T`. This makes the accepted basis exactly equivalent to the reduced system described
between 2006 eqs. (13) and (14).

This sufficiency statement concerns smooth cold bulk npe-mu response in the 2006
Cowling/local-neutrality approximation. The internal equilibration prescription must be
consistent across energy, conjugates, and Hessian; weak equilibrium is not imposed during
independent composition differentiation. Symmetry requires a twice differentiable energy
potential, and invertibility requires a nonsingular physical branch; strict stability is
sufficient, smoothness alone is not. At absent-species boundaries, critical points, and phase
transitions, report the domain restriction rather than fabricate a full-rank Hessian. A regular
Hessian does not supply moving-interface response. Source-specific phase and species treatment
must be governed before a global calculation crossing those boundaries.

**A full paper `B` inverse is forbidden as an API.** The corrected four-species global matrix of
2006 eq. (13) is singular. Only a qualified reduced physical response may be solved or inverted,
and that operation belongs to later governed work.

## Charge-neutrality/electrostatic ownership

| Responsibility | Accepted owner | Boundary |
|---|---|---|
| Validate or construct an allowed local charge-neutral state; return neutral conjugates and the selected local derivative primitive, with optional model-defined intrinsic species potentials | **Off-equilibrium EOS provider** | No GR volume integration, no redshift profile, no paper `B/Z/W`, and no electrostatic gauge variable. |
| Derive the source-basis reduced susceptibility from the constrained EOS Hessian and verify its equivalence to the 2006 charge projection | **Local rotochemical thermodynamics adapter** | `C_y=T H_x^{-1} T^T` already incorporates electrostatic elimination. Do not project it again. A full intrinsic toy/provider extension may be used privately for equivalence checks; neither that extension nor `delta psi` is required from a neutral provider. |
| Integrate the local reduced number response with `e^{-Phi} dV`, construct the corrected reduced global response, apply global baryon conservation, and later obtain paper `Z` coefficients | **Stellar rotochemical coefficient layer (Layer B)** | Must consume the star's governed metric/volume and provider response; must never ask the EOS to integrate a star. |
| Supply the verified fixed-`epsilon_c` rotation response and later baryon-conserving sequence particle response | **Phase-4 structural layer plus governed Phase-5 sequence reduction (Layer C)** | Existing Phase-4 objects remain unchanged; INV-09 work is separate. |
| Evolve redshifted imbalances, temperature, reactions, and heating | **Secular evolution layer (Layer D)** | Not part of this ADR's implementation and still blocked by INV-11 and later rate decisions. |

The provider's mandatory conjugates are local intrinsic neutral combinations. It does not add
electrostatic or redshift terms. If intrinsic species potentials are supplied, they must have a
declared model convention and reproduce those combinations. The neutral Hessian determines the
corrected bulk density response, not the individual charged-potential split or the electrostatic
potential profile. The electrostatic term cancels after enforcing neutrality, not by omitting it
from the unconstrained derivation.

In using 2006 eq. (11), retain the `e^{-Phi}` factor required by substitution into eq. (10) and
present in eq. (13); the printed eq. (11) appears to omit it. This is recorded as an **inferred
printed/source omission**, not as an authenticated published erratum. The source quotation is not
silently rewritten (`docs/validation/PHASE5A1A_ADR0010_INDEPENDENT_ADJUDICATION.md:194-219`).

## Required API/data semantics

This section is accepted semantics, not a C++ design.

| Capability | Phase-5 classification | Required semantics |
|---|---|---|
| Accepted local state coordinates and domain | **REQUIRED NOW** | Three independent charge-neutral coordinates; explicit species set, phase/threshold domain, and `fm^-3` units. |
| `n_B,n_n,n_p,n_e,n_mu` | **REQUIRED NOW** | Inputs plus exact dependent reconstruction; no stored fraction may masquerade as a number density (ADR-0001). |
| Neutral conjugates `g=(mu_n,-eta_npe,-eta_npmu)` | **REQUIRED NOW** | Local, in `MeV`, with one energy/rest-mass convention; consistent with the validation energy and `H`; no electrostatic or GR term. |
| Intrinsic `mu_n,mu_p,mu_e,mu_mu` | **OPTIONAL MODEL CAPABILITY** | Local, in `MeV`, with declared microscopic/component convention; if supplied, require `S_x^T mu=g`. The charged split cannot be inferred from neutral energy/H alone. |
| `eta_npe`, `eta_npmu` | **REQUIRED NOW, DERIVED** | `eta_npe=-g_e` and `eta_npmu=-g_mu`; if species potentials are supplied, their beta combinations must agree. No independently drifting authority. |
| Equilibrium state at a supplied `n_B` | **REQUIRED NOW** | Returns the model's equilibrium composition and active-species/threshold status; recovers the applicable zero-imbalance conditions. |
| `H_ab=partial g_a/partial x_b` | **REQUIRED NOW** | Zero temperature, with the other independent coordinates fixed; coordinate names, units, active smooth domain, phase validity, and failure status are part of every result. |
| Off-equilibrium energy density `epsilon(x)` | **REQUIRED FOR VALIDATION** | Prefer `MeV fm^-3`; supplies a thermodynamic potential for finite-perturbation and integrability checks. It is not required by the Phase-5A global coefficient runtime if the chemical response is supplied directly. |
| Off-equilibrium pressure `P(x)` | **FUTURE-GENERIC** | Desirable for a broader EOS interface and consistency checks, but not used to perturb the stellar structure in the 2005/2006 Cowling treatment. Equilibrium `P,epsilon` continue to come from the governed stellar-background EOS. |
| Temperature coordinate, entropy, second derivatives beyond the selected first response | **FUTURE-GENERIC** | No implicit finite-temperature generalization in a cold v1 contract. |
| Effective masses | **LATER REPRODUCTION / MICROPHYSICS** | Required for the 2005 heat-capacity/rate benchmark (2005 section 3.1 and eq. (50)), not for the local corrected susceptibility contract itself. |
| Provenance and capability metadata | **REQUIRED NOW** | EOS/model identifier, parameter/data revision, coordinate chart, particle/rest-mass convention, units, valid range, phase/threshold flags, derivative orientation, and lepton ownership. |

The consumer-facing provider returns one coherent npe-mu state. It may compose a hadronic model
with independently validated analytic free-Fermi electrons and muons, provided provenance names
every component, rest-mass conventions agree, lepton contributions are neither omitted nor
double counted, and all returned thermodynamics and derivatives are mutually consistent.

No algorithm is selected. An implementation may use analytic differentiation, automatic
differentiation, differentiated interpolation, finite perturbations, or a validated data product
only after a later governed task decides and validates that method.

## Validation requirements

The future implementation must follow this ladder. Numerical tolerances are intentionally absent
until source precision, conditioning, and convergence evidence exist.

| ID | Required check | Defect it falsifies |
|---|---|---|
| **V1** | Units, rest-mass convention, index order, and imbalance signs | Hidden `hbar c`, cgs/geometric, rest-mass, index, or sign error |
| **V2** | Exact charge-neutral reconstruction for every accepted state | Wrong dependent-density map or drift off the physical subspace |
| **V3** | Recover beta equilibrium and species-threshold conditions from `EquilibriumStateAt(n_B)` | Wrong equilibrium anchor, chemical-potential convention, or active-species rule |
| **V4** | Where analytic free-Fermi lepton components or individual lepton capabilities are supplied, check their chemical potentials and derivatives against analytic cold free-Fermi results | Lepton mass, Fermi momentum, unit, or derivative error |
| **V5** | Analytic toy npe and npe-mu energy functionals with known reduced Hessian and susceptibility | Wrong held-fixed variables, orientation, indices, or inversion target |
| **V6** | Maxwell/Hessian symmetry where the source is smooth | Non-integrable or internally inconsistent thermodynamics |
| **V7** | Shrinking finite perturbations versus the returned linear response within one smooth phase | Incorrect derivative values or signs; lack of convergence |
| **V8** | Exact A/B response equivalence; for valid fraction charts include the nonlinear chain-rule term in energy-Hessian transformations | Basis dependence, density/fraction mixing, or omitted second-coordinate derivatives |
| **V9** | In an independent full-intrinsic toy fixture, compare the 2006 projected response with `S_x H_x^{-1} S_x^T`; then verify full global charge null and proton-row identity | Missing or double electrostatic projection; a null identity alone is insufficient to validate response magnitudes |
| **V10** | Establish reduced symmetry, rank, stability, and conditioning on the declared active domain; embed absent-species responses explicitly and verify global mode support before inversion | Hidden gauge mode, fabricated absent-species rank, instability, or invalid threshold continuation |
| **V11** | Reproduce published corrected 2006 `Z` information or its stated changes once authenticated numerical reference data exist | Wrong GR integration, source-basis transform, or global baryon reduction |
| **V12** | After Layers B-D exist, reproduce a published 2005/2006 non-superfluid thermal benchmark | End-to-end coefficient, drive, reaction, redshift, or evolution error |

V9-V12 belong to later layers; listing them here prevents a local API from being accepted without
the downstream falsifiers it must support.

## Alternatives considered

### Expose only the current beta-equilibrium trajectory

Rejected as a contract candidate. Derivatives along that curve provide one directional
combination and cannot recover transverse composition response
(`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:293-324`).

### Make the EOS return the already integrated paper `B_ij`

Rejected. It would make a local material model own the star's GR volume integral and profile
provenance, violate the established layer boundary, and encourage a full inverse that 2006 proves
singular.

### Make a full unconstrained four-species `K` or `chi` mandatory

Not recommended for the minimum interface. Such a provider can be useful, but it requires an
unphysical charged direction that neither a charge-neutral table nor the Phase-5 consumer needs.
The public contract is the physical constrained response; a microscopic provider may retain the
full object internally and demonstrate equivalence to the 2006 projected response.

### Use fractions as the canonical derivative basis

Not recommended. Fractions remain useful data coordinates and ADR-0001 profile semantics, but
their nonlinear, mixed-unit transform obscures derivative units and becomes singular at
`n_B=0`. V8 must nevertheless prove equivalence where the chart is valid.

### Use separate model-specific Phase-5 interfaces

Rejected as the target architecture. Track R and Track P may use different sources and algorithms,
but their consumer-facing semantics should be common so a reproduced literature EOS and a CMF
provider are tested by the same physical invariants.

### Require a temperature-aware interface immediately

Deferred. The standard 2005 susceptibility coefficients are built from a cold equilibrium EOS;
the current task does not establish that a temperature coordinate improves the first contract
enough to justify its added state dimension.

## Non-goals

- No EOS derivative, interpolation, table, model, or matrix algorithm is implemented.
- No corrected paper `B`, reduced global matrix, paper `Z/W`, structural `A_i/B_i/Z_i`, or
  baryon-conserving sequence response is constructed.
- No evolved chemical-state ordering/redshift convention is accepted; INV-11 remains unresolved
  beyond the local neutral-conjugate semantics governed here.
- No rotochemical ODE, reaction rate, neutrino modification, heating, superfluidity, BNV, or
  baseline is added.
- No current EOS, `StarProfile`, Phase-4 response, test, or scientific data product is changed.
- No APR, BPAL, CMF, or CompOSE data are downloaded, fabricated, or selected as accepted
  authority.

## Migration implications

This accepted ADR requires a new local thermodynamic provider boundary rather than modification
of the current `TOVSolver` background API in place. The current
`P(n_B),epsilon(n_B),Y_i(n_B)` path remains authoritative for nonrotating structure; a matching
off-equilibrium provider supplies the transverse thermodynamics and proves equilibrium-point
identity.

Track R providers follow the staged ladder and close on the authenticated A18 + delta-v + UIX*
benchmark. Canonical Track P requires the source-matched composition provider defined above. A
distinct analytic/RMF model is a separate owner-governed research EOS and does not retain
canonical DS(CMF) identity by assertion. Current
`SigmaOmegaRho`, `Fermi_Gas_Many`, and dormant `SigmaOmegaRho_npemu` code is evidence only and is
not activated (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:161-188`, `:461-479`).

The global stellar layer would consume local reduced response values at the nonrotating profile
state, integrate them with the governed geometry, and preserve source/profile provenance. The
Phase-4 objects remain the only structural rotation inputs. Particle-number measures and the
baryon-conserving sequence reduction remain governed separately by INV-09 and ADR-0008.

## Owner-ratified decisions

The project owner accepts all six decisions below in Phase 5A-1B.

| ID | Ratified answer | Consequence |
|---|---|---|
| **Q1 — Track R acceptance benchmark** | Track R closes on quantitative reproduction of A18 + delta-v + UIX*. Analytic/toy checks and the 2005 whole-star free-gas case precede it; an authenticated BPAL intermediate is optional. The 2006 correction governs its supersession scope. Closure requires a correction-sensitive coefficient or transient comparison; quasi-steady temperature alone is insufficient. CMF cannot substitute. | Missing authenticated source definition, phase/crust prescription, or comparison data blocks Track R closure. |
| **Q2 — canonical local basis** | Use `x=(n_B,n_e,n_mu)`, with `n_p=n_e+n_mu` and `n_n=n_B-n_e-n_mu`. Require exact equivalence to `y=(n_n,n_e,n_mu)` on the declared smooth npe-mu domain. | Thresholds and absent species are explicit domain boundaries and never silently change coordinate meaning. |
| **Q3 — derivative authority** | The EOS-facing primitive is `H_ab=partial g_a/partial x_b`, where `g=(mu_n,-eta_npe,-eta_npmu)`, at zero temperature with the other independent coordinates fixed. The constrained Hessian already includes electrostatic elimination. The adapter may derive the reduced susceptibility but must not project it again. No full paper-`B` inverse is exposed. Individual charged intrinsic potentials and `delta psi` are optional model-specific capabilities only. | Applies only on a smooth, nonsingular active-species domain. Phase transitions, critical points, absent-species boundaries, and moving interfaces require later explicit treatment. |
| **Q4 — lepton boundary** | Expose one coherent consumer-facing npe-mu provider. It may internally compose a hadronic EOS with independently validated analytic free-Fermi electrons and muons when provenance, rest masses, interactions, and cross derivatives are consistent and lepton terms are neither omitted nor double counted. | A single concrete monolithic model is not required, but one mutually consistent outward thermodynamic state is. |
| **Q5 — Track P data authority** | Canonical DS(CMF) identity requires an authenticated source-matched arbitrary-composition provider/product: model revision, parameter set, particle content, phase treatment, mass and lepton conventions, and equilibrium-background agreement. A matched hadronic composition surface plus validated leptons is allowed. The current beta-equilibrium table, distinct finite-T `Y_q` table, dormant RMF code, or a family name alone are not equivalent. Thermodynamic-response v1 is restricted to a verified npe-mu domain; additional responding species require separate governance. | If no matched source exists, canonical Track P stops for owner adjudication rather than substituting another EOS. |
| **Q6 — temperature scope** | Thermodynamic response v1 is explicitly cold-only. This does not imply zero-temperature thermal evolution: finite-temperature weak rates, heat capacities, and neutrino processes remain later inputs. | Finite-temperature thermodynamic response requires a separate contract revision naming the thermodynamic potential and held-fixed variables. |

## Decision

**ACCEPTED — 2026-09-04**

The project owner ratifies Q1-Q6 exactly as recorded in the preceding table, including the
Phase 5A-1A revisions R1-R7. The accepted contract governs the cold charge-neutral local
thermodynamic provider boundary and its validation obligations. It does not select a derivative,
interpolation, or inversion algorithm; authorize a realistic EOS substitution; implement global
paper `B/Z/W`; resolve INV-09 or the remaining INV-11 evolution-state convention; or authorize
rotochemical ODEs, reaction rates, superfluidity, or BNV.

## Provenance

The historical sequence is preserved:

1. **Phase 5A-0 — capability audit:** established that the current equilibrium EOS path lacks the
   required independent-composition response
   (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:10-31`, `:511-528`).
2. **Phase 5A-1 — proposed contract:** drafted ADR-0010 as PROPOSED from
   `77a328676d83f515fe603cb62341d2efcd70ed78`, with Decision `PENDING OWNER ADJUDICATION`.
3. **Phase 5A-1A — independent Astra adjudication:** independently checked the local 2005/2006 papers
   and proved equivalence of the neutral Hessian to the corrected reduced bulk response. Its
   proposed revisions and exact Q1-Q6 answers are recorded in
   `docs/validation/PHASE5A1A_ADR0010_INDEPENDENT_ADJUDICATION.md:7-19`, `:468-614`. That review
   remained historical evidence and did not exercise owner authority or authorize implementation.
4. **Phase 5A-1B — owner acceptance:** the project owner accepted Q1-Q6 and R1-R7 on 2026-09-04
   from the authenticated starting commit `e4d1b8f4b30d10760bc0c894cee959040db6f9a3`.

The four PDFs named under Source authority were read directly from the shared read-only library,
with the 2006 paper governing its stated correction scope. Current-code and data-product facts
were checked against the Phase 5A-0 audit and current source. No gravitochemical paper or untracked
scientific source was used as authority.

Phase 5A-1B changes governance documentation only. It does not modify production source, tests,
baselines, EOS data, or the shared literature library; implement the accepted provider; or begin
stellar coefficient integration, evolution, or BNV.
