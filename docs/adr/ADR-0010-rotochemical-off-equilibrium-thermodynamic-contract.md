# ADR-0010 — Rotochemical off-equilibrium thermodynamic contract

| Field | Value |
|---|---|
| **Status** | **PROPOSED** |
| **Date** | 2026-09-04 |
| **Starting SHA** | `77a328676d83f515fe603cb62341d2efcd70ed78` |
| **Change class** | scientific-semantic and structural/architecture proposal; documentation only in Phase 5A-1 |
| **Governing authority** | `GOVERNANCE.md`; ADR-0001; ADR-0002; ADR-0007; ADR-0008; ADR-0009; INV-01, INV-09, INV-11; `docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md` |
| **Evidence companion** | `docs/validation/PHASE5A1_THERMODYNAMIC_CONTRACT_DERIVATION.md` |
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

This ADR is deliberately **PROPOSED**. The cited papers constrain the physical subspace and the
observable response, but do not force one software coordinate basis, one derivative orientation,
one lepton-composition boundary, one data product, or one temperature-generalization policy.

## Source authority

The scientific source root is `/Users/keeper/Documents/CompactStar/literature`. The exact sources
used are:

1. `/Users/keeper/Documents/CompactStar/literature/rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`
   — primary authority where electrostatic potential, local charge neutrality, the singular full
   response, or the reduced inverse is at issue: eqs. (2), (4), and (9)-(19), sections 2-5,
   PDF pp. 2-4 / journal pp. 569-571.
2. `/Users/keeper/Documents/CompactStar/literature/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`
   — primary non-superfluid reproduction benchmark: eqs. (7)-(13), PDF p. 4; eq. (24), PDF p. 6;
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

This proposal governs only the local thermodynamic contract and its interfaces to global stellar
susceptibility integration and the verified Phase-4 structural response. It does not implement
any of those layers.

## Track R vs Track P

| Track | Scientific purpose | Proposed entry and acceptance logic |
|---|---|---|
| **R — literature reproduction** | Reproduce the standard 2005 non-superfluid calculation with the 2006 correction where it governs. | Use a staged validation ladder: analytic free-gas/toy cases, then a published realistic EOS, ending with the paper's reference A18 + delta-v + UIX* benchmark if the required primary EOS definition and numerical data can be authenticated. A current CMF curve is not a substitute. |
| **P — CompactStar production** | Apply already-validated machinery to DS(CMF)-1 and later research EOS families. | Supply the same local contract through a source-matched arbitrary-composition EOS or data product. The current equilibrium DS(CMF)-1 file remains a stellar-background source only. |

The proposed architecture is one generic consumer-facing contract with model-specific providers.
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
| **A. `(n_B,n_e,n_mu)`** | Complete three-density chart; `n_p=n_e+n_mu`, `n_n=n_B-n_e-n_mu` enforce neutrality exactly. | All inputs are `fm^-3`; no division by `n_B`; positivity domain is explicit. Species scales still differ and must be reported by numerical implementations. | Linearly equivalent to the source-reduced basis. Its conjugates are `(mu_n,-eta_npe,-eta_npmu)`, so observable imbalances are direct. | Natural `n_B` axis for analytic, tabulated, and RMF providers; requires two independent lepton-density directions. | **Recommended, pending Q2.** |
| **B. `(n_n,n_e,n_mu)`** | Complete; `n_p=n_e+n_mu`, `n_B=n_n+n_e+n_mu`. | All inputs are `fm^-3`; linear constraints. | Directly conjugate to the 2006 variables `(mu_n,mu_p+mu_e,mu_p+mu_mu)` used before eq. (14). | Natural for the paper and energy Hessians; less directly aligned with standard `n_B`-axis tables. | Equally physical; mandatory exact transform target for validation. |
| **C. `(n_B,Y_e,Y_mu)`** | Complete for `n_B>0`; `n_i=n_BY_i`. | Mixed units and a nonlinear transform; degenerates as `n_B` approaches zero and scales composition derivatives by `n_B`. | Equivalent away from `n_B=0`, but not the source's linear density chart. | Convenient table/query representation. | Allow as an adapter or data-storage chart, not recommended as the canonical derivative chart. |
| **D. full `(n_n,n_p,n_e,n_mu)`** | Contains an unphysical charged direction unless accompanied by a constraint/projection. | Uniform density units but redundant on the physical manifold. | Needed only as an intermediate reconstruction of 2006 eqs. (11)-(13). | Some microscopic EOSs may naturally provide it. | Not recommended as the public Phase-5 primitive; must never imply that the corrected full matrix is invertible. |

For A and B, define

```text
(n_n,n_e,n_mu)^T = T (n_B,n_e,n_mu)^T,
T = [[1,-1,-1],[0,1,0],[0,0,1]].
```

The determinant is one, and the inverse is exact. Any conforming provider or adapter must produce
equivalent physical responses under this transformation. The canonical choice remains an owner
decision because the papers establish equivalence, not a software preference.

## Candidate derivative authorities

Let `K_ij = partial mu_i / partial n_j` and
`chi_ij = partial n_i / partial mu_j`, with the independent variables and smooth phase domain
made explicit.

| Candidate authority | Relation to sources | Singularity and symmetry | Numerical / architectural consequence | Assessment |
|---|---|---|---|---|
| **Full intrinsic `K` over four species** | Inverse orientation of 2005 eq. (10), before the 2006 charge projection. | A thermodynamic Hessian is symmetric where smooth; it is not itself the corrected paper `B`. | Requires an unphysical charged direction, a local solve/inverse, the 2006 projection, and later volume integration. | Broader than Phase 5 requires; easy to misuse. Not recommended as the minimum contract. |
| **Full intrinsic `chi` over four species** | Direct local object in 2005 eq. (10). | Symmetric by the identity cited after 2006 eq. (13). After applying eq. (11), the charge-projected response has a null charged/gauge mode. | Most literal source input, but requires non-neutral response data and invites forbidden inversion of the corrected full global matrix. | Permissible provider capability, not the proposed public primitive. |
| **Reduced constrained chemical Hessian** | Restriction of intrinsic thermodynamics to the local charge-neutral manifold; exactly transformable to the 2006 reduced variables. | Symmetric in conjugate coordinates and nonsingular/physically stable inside a smooth stable phase; threshold and phase-boundary failures must be explicit. | The EOS returns the most local physical derivative. A later local adapter obtains the reduced number susceptibility; the stellar layer integrates it. | **Recommended, pending Q3.** |
| **Reduced constrained number susceptibility** | Direct local integrand orientation for the reduced 2006 `Btilde`. | Symmetric and nonsingular on the same physical domain; it is the inverse response of the reduced chemical Hessian, not the inverse of the singular full paper `B`. | Simplifies the stellar integrand but makes every provider own the inverse-orientation result. | Source-aligned alternative; owner must choose whether it or the Hessian is authoritative. |

The recommended semantic primitive is provisionally named
`ChargeNeutralChemicalHessian`:

```text
H_ab = partial g_a / partial x_b,
x = (n_B, n_e, n_mu),
g = (mu_n, -eta_npe, -eta_npmu).
```

Its units are `MeV fm^3`. The inverse response, if later formed, has units
`fm^-3 MeV^-1`. This proposal specifies neither how a provider obtains `H` nor how a later layer
solves for its inverse response.

The source-basis equivalent uses
`y=(n_n,n_e,n_mu)` and
`g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu)`. If `H_y` is its reduced Hessian, then
`H_x=T^T H_y T`. This makes the proposed basis exactly equivalent to the reduced system described
between 2006 eqs. (13) and (14).

**A full paper `B` inverse is forbidden as an API.** The corrected four-species global matrix of
2006 eq. (13) is singular. Only a qualified reduced physical response may be solved or inverted,
and that operation belongs to later governed work.

## Charge-neutrality/electrostatic ownership

| Responsibility | Proposed owner | Boundary |
|---|---|---|
| Validate or construct an allowed local charge-neutral state; return intrinsic species chemical potentials and the selected local derivative primitive | **Off-equilibrium EOS provider** | No GR volume integration, no redshift profile, no paper `B/Z/W`, and no electrostatic gauge variable. |
| Apply the source-defined charge-neutral basis/projection, demonstrate equivalence to 2006 eqs. (11)-(13), and transform between the canonical and source bases | **Local rotochemical thermodynamics adapter** | Owns the electrostatic **correction/projection**, not a new EOS and not a stellar integral. `delta psi` may be represented for validation but is not an observable state output. |
| Integrate the local reduced number response with `e^{-Phi} dV`, construct the corrected reduced global response, apply global baryon conservation, and later obtain paper `Z` coefficients | **Stellar rotochemical coefficient layer (Layer B)** | Must consume the star's governed metric/volume and provider response; must never ask the EOS to integrate a star. |
| Supply the verified fixed-`epsilon_c` rotation response and later baryon-conserving sequence particle response | **Phase-4 structural layer plus governed Phase-5 sequence reduction (Layer C)** | Existing Phase-4 objects remain unchanged; INV-09 work is separate. |
| Evolve redshifted imbalances, temperature, reactions, and heating | **Secular evolution layer (Layer D)** | Not part of this ADR's implementation and still blocked by INV-11 and later rate decisions. |

The EOS returns **intrinsic local** chemical potentials. It must not add `q_i psi`, choose an
electrostatic zero, or silently redshift them. The local adapter carries the rigorous 2006
projection. The observable imbalances are gauge invariant because their charge sums vanish, not
because the electrostatic term may be omitted from the derivation.

## Required API/data semantics

This section is proposed semantics, not a C++ design.

| Capability | Phase-5 classification | Required semantics |
|---|---|---|
| Accepted local state coordinates and domain | **REQUIRED NOW** | Three independent charge-neutral coordinates; explicit species set, phase/threshold domain, and `fm^-3` units. |
| `n_B,n_n,n_p,n_e,n_mu` | **REQUIRED NOW** | Inputs plus exact dependent reconstruction; no stored fraction may masquerade as a number density (ADR-0001). |
| Intrinsic `mu_n,mu_p,mu_e,mu_mu` | **REQUIRED NOW** | Local, in `MeV`; consistent rest-mass convention declared; no electrostatic or GR redshift term. |
| `eta_npe`, `eta_npmu` | **REQUIRED NOW, DERIVED** | Derived from returned intrinsic chemical potentials with the stated signs; not an independently drifting stored authority. |
| Equilibrium state at a supplied `n_B` | **REQUIRED NOW** | Returns the model's equilibrium composition and active-species/threshold status; recovers the applicable zero-imbalance conditions. |
| Selected reduced local derivative authority | **REQUIRED NOW** | Coordinate names, held-fixed semantics, orientation, units, phase validity, and failure status are part of every result. |
| Off-equilibrium energy density `epsilon(x)` | **REQUIRED FOR VALIDATION** | Prefer `MeV fm^-3`; supplies a thermodynamic potential for finite-perturbation and integrability checks. It is not required by the Phase-5A global coefficient runtime if the chemical response is supplied directly. |
| Off-equilibrium pressure `P(x)` | **FUTURE-GENERIC** | Desirable for a broader EOS interface and consistency checks, but not used to perturb the stellar structure in the 2005/2006 Cowling treatment. Equilibrium `P,epsilon` continue to come from the governed stellar-background EOS. |
| Temperature coordinate, entropy, second derivatives beyond the selected first response | **FUTURE-GENERIC** | No implicit finite-temperature generalization in a cold v1 contract. |
| Effective masses | **LATER REPRODUCTION / MICROPHYSICS** | Required for the 2005 heat-capacity/rate benchmark (2005 section 3.1 and eq. (50)), not for the local corrected susceptibility contract itself. |
| Provenance and capability metadata | **REQUIRED NOW** | EOS/model identifier, parameter/data revision, coordinate chart, particle/rest-mass convention, units, valid range, phase/threshold flags, derivative orientation, and lepton ownership. |

The consumer-facing provider should return one coherent npe-mu state even if its implementation
composes a hadronic model with analytic free-Fermi leptons. Whether that composite boundary is
mandatory is owner question Q4.

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
| **V4** | Electron and muon chemical potentials and derivatives against analytic cold free-Fermi results | Lepton mass, Fermi momentum, unit, or derivative error |
| **V5** | Analytic toy npe and npe-mu energy functionals with known reduced Hessian and susceptibility | Wrong held-fixed variables, orientation, indices, or inversion target |
| **V6** | Maxwell/Hessian symmetry where the source is smooth | Non-integrable or internally inconsistent thermodynamics |
| **V7** | Shrinking finite perturbations versus the returned linear response within one smooth phase | Incorrect derivative values or signs; lack of convergence |
| **V8** | Exact response equivalence under A/B and, for `n_B>0`, C coordinate transforms | Hidden basis dependence or fraction/density mixing |
| **V9** | Test-side construction of corrected full paper `B` has the 2006 charge/null mode and row identity | Missing electrostatic projection or accidental use of the 2005 matrix |
| **V10** | The reduced `3 x 3` physical system has the expected rank, symmetry/stability, and documented conditioning across its valid domain | Wrong reduction, threshold handling, or unstable physical response |
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
The public contract should be the physical constrained response; a microscopic provider may retain
the full object internally and demonstrate the required projection.

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
- No chemical-state ordering/redshift convention is accepted; INV-11 remains unresolved.
- No rotochemical ODE, reaction rate, neutrino modification, heating, superfluidity, BNV, or
  baseline is added.
- No current EOS, `StarProfile`, Phase-4 response, test, or scientific data product is changed.
- No APR, BPAL, CMF, or CompOSE data are downloaded, fabricated, or selected as accepted
  authority.

## Migration implications

If accepted in a later owner action, this ADR would require a new local thermodynamic provider
boundary rather than modification of the current `TOVSolver` background API in place. The current
`P(n_B),epsilon(n_B),Y_i(n_B)` path remains authoritative for nonrotating structure; a matching
off-equilibrium provider supplies the transverse thermodynamics and proves equilibrium-point
identity.

Track R providers would begin with independent analytic fixtures and later an authenticated
published realistic EOS. Track P would require an explicitly source-matched multidimensional
DS(CMF)-family product or a separately governed analytic/model provider. Current
`SigmaOmegaRho`, `Fermi_Gas_Many`, and dormant `SigmaOmegaRho_npemu` code is evidence only and is
not activated (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:161-188`, `:461-479`).

The global stellar layer would consume local reduced response values at the nonrotating profile
state, integrate them with the governed geometry, and preserve source/profile provenance. The
Phase-4 objects remain the only structural rotation inputs. Particle-number measures and the
baryon-conserving sequence reduction remain governed separately by INV-09 and ADR-0008.

## Open owner questions

These six questions are the substantive choices not forced by the source.

| ID | Owner question | Recommended answer | Consequence |
|---|---|---|---|
| **Q1 — Track R acceptance benchmark** | Must Track R close on the paper's reference A18 + delta-v + UIX* EOS, or may another published realistic EOS close it? | Require the staged ladder and use A18 + delta-v + UIX* as the closing quantitative benchmark; BPAL may be an intermediate, not a silent substitute. | Requires acquisition/authentication of the primary APR EOS definition, exact fit/table and benchmark data before Track R can close; CMF work cannot claim literature reproduction. |
| **Q2 — canonical local basis** | Should the generic contract use A `(n_B,n_e,n_mu)`, B `(n_n,n_e,n_mu)`, or another source-equivalent chart? | A, with exact B transformation as a mandatory validation. | Aligns the API with standard EOS `n_B` axes and makes `-eta` thermodynamic conjugates; adapters must preserve source-basis equivalence. |
| **Q3 — derivative authority** | Should the EOS own the reduced chemical Hessian, reduced number susceptibility, or a full unconstrained object? | Own the reduced charge-neutral chemical Hessian; derive the reduced susceptibility in the local rotochemical adapter. Never expose a full paper-`B` inverse. | Gives one local, gauge-free, symmetry-testable primitive; a later governed task must choose and validate the solve/differentiation algorithm. |
| **Q4 — lepton boundary** | Must leptons be supplied by the same concrete model as baryons, or may a unified provider compose a hadronic EOS with analytic free-Fermi leptons? | Require one unified consumer-facing npe-mu provider but permit internally composable analytic lepton components. | Consumers see one consistent state; Track R gains independent analytic lepton falsifiers; provenance must identify both components. |
| **Q5 — Track P data authority** | Is Track P to obtain a DS(CMF)-matched multidimensional data product, or to adopt and validate an analytic/RMF model as a distinct research EOS? | Prefer a source-matched DS(CMF)-family arbitrary-composition product; if it is unavailable, stop for explicit owner adjudication rather than treating the current finite-T `Y_q` table or dormant RMF code as equivalent. | Determines data acquisition, equilibrium-background matching, and whether Track P results retain the canonical EOS identity. |
| **Q6 — temperature scope** | Is the first contract explicitly cold-only, or temperature-aware from the start? | Cold-only v1 with an explicit future extension point. | Matches the standard benchmark and limits validation dimensionality; finite-temperature response requires a later contract revision rather than an implicit extra held-fixed variable. |

## Decision

**PENDING OWNER ADJUDICATION**

The recommendations above are not accepted authority. No implementation may treat Q1-Q6 as
settled until the project owner records the decisions and changes this ADR's status in a separate
action.

## Provenance

Drafted by one Codex agent in Phase 5A-1 from the exact starting commit
`77a328676d83f515fe603cb62341d2efcd70ed78`. The four PDFs named under Source authority were read
directly from the shared read-only library, with the 2006 paper governing its stated correction
scope. Current-code and data-product facts were rechecked against the Phase 5A-0 audit and current
source. No memory, gravitochemical paper, web substitute, or untracked scientific source was used
as authority.

This proposal changes documentation only. It does not modify production source, tests, baselines,
EOS data, or the shared literature library, and it does not begin BNV.
