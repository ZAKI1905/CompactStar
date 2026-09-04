# Phase 5A-1 — thermodynamic contract derivation and evidence

**Date:** 2026-09-04

**Starting SHA:** `77a328676d83f515fe603cb62341d2efcd70ed78`

**Branch:** `analysis/phase5a-eos-thermodynamic-audit`

**Worktree:** `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-eos-audit`

**Change class:** source-backed scientific-semantic and structural/architecture governance;
documentation only.

**Result:** `docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md` is a
**PROPOSED** contract with **PENDING OWNER ADJUDICATION**. Nothing in this record is an accepted
scientific decision or an implementation authorization.

## 1. Entry authentication and inherited gates

At entry the worktree was clean at exactly
`77a328676d83f515fe603cb62341d2efcd70ed78`; local, cached-upstream, and live-remote branch heads
were equal. Local `master`, cached `origin/master`, and live remote `master` were all exactly
`59526c56122dbbb8c0a8ef808bdf627453d99c3a`. No old Phase-4 worktree was reused.

The inherited gates were read before drafting:

- Phase 4 is complete and its structural interface is ratified
  (`docs/validation/PHASE4_CLOSEOUT.md:1-9`, `:64-176`).
- Phase 5A-0 is complete as a documentation audit and found a minimum gap of class D: a new
  governed off-equilibrium thermodynamic extension or data product is required
  (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:10-31`, `:511-528`).
- No Phase-5 physics was implemented by Phase 5A-0
  (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:559-565`).
- INV-09 remains unresolved for particle-number response and fixed-baryon sequence reduction;
  INV-11 remains unresolved for the evolved chemical-imbalance convention
  (`docs/SCIENTIFIC_INVARIANTS.md:721-805`, `:816-836`).

This task does not change any of those statuses.

## 2. Source ledger

The library root `/Users/keeper/Documents/CompactStar/literature`, its `README.md`, and its
`catalog.tsv` were read first. The four PDFs below were then read directly. The 2006 paper governs
over 2005 wherever electrostatic potential, local charge neutrality, or inversion of the full
thermodynamic response conflicts.

| Requirement or statement used here | Exact PDF filename | Equation / page / section | Role |
|---|---|---|---|
| The npe and np-mu imbalance definitions | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (7)-(8), PDF p. 4, section 2.1 | benchmark |
| 2005 treated intrinsic redshifted chemical-potential deviations as uniform without the electrostatic completion | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (9), PDF p. 4, section 2.1 | superseded where corrected |
| The local linear response is `delta n_i=sum_j (partial n_i/partial mu_j) delta mu_j`, evaluated at beta equilibrium | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (10), PDF p. 4, section 2.1 | benchmark requirement |
| Integrated number response, original paper `B_ij`, and naive full inverse | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (11)-(13), PDF p. 4, section 2.1 | eqs. (12)-(13) corrected by 2006 |
| Global reaction and changing-equilibrium particle numbers | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (14)-(17), PDF pp. 4-5, section 2.1 | benchmark |
| Rotating scalar particle number contains composition, metric/volume, frame-dragging, and boundary pieces | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (24), PDF p. 6, section 2.2 | Layers B/C |
| Equilibrium particle-number spin drive `I_{Omega,i}` | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (30), PDF p. 7, section 2.2 | Layers B/C |
| Quantitative work needs equilibrium structure/composition, the partials in eq. (10), and effective masses; energy versus baryon density and proton fraction can supply the partials | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | section 3.1, PDF pp. 8-9 | benchmark EOS requirements |
| The three EOS sets are two APR models, five BPAL models, and a noninteracting Fermi gas; A18 + delta-v + UIX* is the reference EOS | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | section 3.1, PDF pp. 8-9 | reproduction authority |
| Effective masses enter the non-superfluid heat capacity | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (50), PDF p. 12, section 3.4 | later reproduction input |
| Non-superfluid evolution and original `Z/W` definitions | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (51)-(58), PDF p. 13, section 3.6 | benchmark; eqs. (54)-(56) corrected |
| Total redshifted chemical potential includes intrinsic `mu_i` and `q_i psi` | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (2), PDF p. 2 / journal p. 569, section 2 | corrected authority |
| Local charge neutrality is imposed on densities | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (4), PDF p. 2 / journal p. 569, section 2 | corrected authority |
| Electrostatic contributions cancel in observable beta-channel imbalances | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (9), PDF p. 3 / journal p. 570, section 3 | corrected authority |
| Perturbing the total redshifted potential retains the electrostatic term; the metric perturbation is separately neglected in Cowling | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (10) and following text, PDF p. 3 / journal p. 570, section 4 | corrected authority |
| Neutrality determines `delta psi` | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (11), PDF p. 3 / journal p. 570, section 4 | corrected authority |
| The corrected integrated response contains the charge-projection term | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eqs. (12)-(13), PDF p. 3 / journal p. 570, section 4 | corrected authority |
| The corrected full four-species paper `B` is singular; `B_pj=B_ej+B_muj`; a charge/gauge perturbation changes no distribution | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | text after eq. (13), PDF p. 3 / journal p. 570, section 4 | corrected authority |
| Define `tilde mu_n=delta mu_n` and `tilde mu_l=delta mu_p+delta mu_l`, remove the proton row/column, and invert the symmetric `3 x 3` reduced matrix | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | text between eqs. (13)-(14), PDF p. 3 / journal p. 570, section 4 | corrected authority |
| Global baryon conservation reduces the integrated departures further | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eqs. (14)-(15) and footnote 4, PDF p. 3 / journal p. 570, section 4 | corrected authority |
| Corrected `Z_np`, `Z_npe`, `Z_npmu`, and `W` | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eqs. (16)-(19), PDF p. 3 / journal p. 570, section 4 | corrected authority |
| Only 2005 eqs. (9), (12), (13), and (54)-(56) need formal replacement; the quasi-steady state is unchanged and transient changes are small | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | sections 5-6, PDF p. 4 / journal p. 571 | corrected scope |
| Charge-neutral npe matter has one composition parameter `x=n_p/n` and equilibrium minimizes energy with respect to it | `1995-Reisenegger-Rotochemical-Heating.pdf` | section 2, paper pp. 4-5 / PDF pp. 14-15, eq. (3) | supporting coordinate example |
| Later review states that the corrected charge identity makes full `B_ij` noninvertible and requires a submatrix | `2020-Yanagi-NS-Therm-Thesis.pdf` | section 4.1.1, thesis pp. 58-60 / PDF pp. 60-62, eq. (4.7) and surrounding text | supporting review |

No source equation was reconstructed from memory. No gravitochemical paper was used as primary
authority.

## 3. Locked Phase-5 scientific order

The owner intent supplied for Phase 5A-1 fixes the order:

1. reproduce **standard non-superfluid rotochemical heating**;
2. reproduce/check the corrected 2006 electrostatic treatment;
3. only afterward extend to superfluidity or other refinements; and
4. leave BNV in Phase 6, blocked until the standard rotochemical benchmark is reproduced.

This prevents a DS(CMF)-1 result from being relabelled as a reproduction of the 2005 calculation
and prevents a later BNV generalization from being built on an unvalidated standard mechanism.

## 4. Current capability baseline

Phase 5A-0's code-and-data inventory remains the factual baseline:

- `TOVSolver` reads three fixed EOS columns plus generic extras, then constructs Steffen splines
  of `epsilon(P)`, `n_B(P)`, and each extra
  (`CompactStar/Core/src/TOVSolver.cpp:597-715`).
- Its relevant public values are `GetEDens(P)`, `GetRho(P)`, generic/misnamed `GetRho_i(P)`, and
  the barotropic `GetEDensDeriv(P)`; no chemical-potential or composition-Jacobian API exists
  (`CompactStar/Core/src/TOVSolver.cpp:1094-1104`, `:1551-1583`).
- ADR-0001 makes the carried extras dimensionless `Y_i=n_i/n_B`, not species densities
  (`docs/SCIENTIFIC_INVARIANTS.md:64-111`).
- The raw cold DS(CMF)-1 product has one `T` point, one placeholder `Y_q` point, and 1,191
  `n_B` points; its collective `Q3-Q5` potentials lie only on that equilibrium path and are not
  exposed by the production reader
  (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:201-225`).
- The separate DNS CMF table has `(T,n_B,Y_q)` and raw `Q1-Q7`, but current
  `CompOSE_Thermo` stores/exposes only `Q2`; it has no independent electron/muon split and is a
  different product from the canonical cold structure EOS
  (`CompactStar/EOS/src/CompOSE_Thermo.cpp:278-356`;
  `docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:227-242`).
- The compiled `Model`, RMF, particle, free-gas, and CompOSE classes do not provide a tested,
  generic, arbitrary-composition npe-mu state and derivative contract
  (`CompactStar/EOS/CMakeLists.txt:1-39`;
  `docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:147-188`).
- Dormant `SigmaOmegaRho_npemu` and the old rotochemical code are historical candidates, not
  authority or a current capability
  (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:461-479`).

Therefore neither equilibrium-path chemical potentials nor `d mu_i/d n_B` along the current curve
can identify the required transverse response. If `x_eq(n_B)` is the equilibrium curve, then an
along-curve derivative has the form

```text
d mu_i/d n_B = sum_j (partial mu_i/partial x_j) (d x_j/d n_B),
```

which is one contraction of a three-dimensional Jacobian. Infinitely many transverse Jacobians
have the same contraction. This is a rank/identifiability gap, not an interpolation gap.

## 5. Local coordinate count and exact basis equivalence

### 5.1 Count the physical degrees of freedom

Use species order `(n,p,e,mu)`. There are four local densities. Local neutrality, 2006 eq. (4),
provides one independent relation:

```text
n_p - n_e - n_mu = 0.
```

`n_B=n_n+n_p` defines a useful coordinate but does not fix its value. Therefore the local
charge-neutral manifold has dimension `4-1=3`.

The two beta-equilibrium conditions `eta_npe=0` and `eta_npmu=0` select a one-dimensional
equilibrium trajectory inside that manifold where both lepton species are active. Small
rotochemical departures need the `n_B` direction plus two independent composition directions.

Global baryon conservation is different. The 2006 paper explicitly says it is only globally
conserved (footnote 4) and applies it to integrated `delta N_i` in eq. (14). Imposing
`delta n_B=0` locally would erase the compression direction and would not be source equivalent.

### 5.2 Source-reduced basis

Choose

```text
y = (n_n,n_e,n_mu)^T
```

and map it into the four species densities with

```text
n = S_y y,

S_y = [[1,0,0],
       [0,1,1],
       [0,1,0],
       [0,0,1]].
```

Thus `n_p=n_e+n_mu`. For intrinsic chemical-potential vector
`mu=(mu_n,mu_p,mu_e,mu_mu)^T`, the energy differential restricted to the neutral manifold is

```text
d epsilon = mu^T S_y d y = g_y^T d y,
g_y = S_y^T mu = (mu_n,mu_p+mu_e,mu_p+mu_mu)^T.
```

The perturbations of `g_y` are exactly the intrinsic, local versions of the three combinations
defined between 2006 eqs. (13) and (14). Because `S_y^T q=0` for
`q=(0,+1,-1,-1)^T`, the electrostatic contribution cancels only after projection onto this
physical basis.

### 5.3 Baryon-density basis

Choose instead

```text
x = (n_B,n_e,n_mu)^T,
n_n = n_B-n_e-n_mu,
n_p = n_e+n_mu.
```

The two reduced density bases obey

```text
y = T x,
T = [[1,-1,-1],
     [0, 1, 0],
     [0, 0, 1]],
det(T)=1.
```

Their conjugates transform as

```text
g_x = T^T g_y
    = (mu_n,mu_p+mu_e-mu_n,mu_p+mu_mu-mu_n)^T
    = (mu_n,-eta_npe,-eta_npmu)^T.
```

Hence

```text
d epsilon = mu_n d n_B - eta_npe d n_e - eta_npmu d n_mu.
```

This makes candidate A physically transparent while retaining exact source equivalence. The
equivalence is mathematical; choosing A as the canonical software chart remains Q2 for the owner.

### 5.4 Fraction basis

For `z=(n_B,Y_e,Y_mu)` with `n_e=n_BY_e` and `n_mu=n_BY_mu`, the conjugates become

```text
(mu_n-Y_e eta_npe-Y_mu eta_npmu,
 -n_B eta_npe,
 -n_B eta_npmu).
```

This chart is equivalent only for `n_B>0`; it is nonlinear, has mixed derivative units, and
degenerates at zero density. It is useful for storage and adapters but is not recommended as the
canonical derivative chart.

## 6. Derivative objects, electrostatic projection, and the null mode

### 6.1 Full intrinsic response

In a smooth phase, define the intrinsic four-species objects

```text
K = partial mu / partial n,
chi = partial n / partial mu = K^{-1}.
```

`chi` is the local object appearing in 2005 eq. (10). Thermodynamic integrability gives the
symmetry identity cited after 2006 eq. (13). The full intrinsic `K` or `chi` need not itself have
the 2006 gauge null; the null appears after the electrostatic/neutrality projection. Conflating
these objects with the corrected global paper `B` is a notation and physics error.

### 6.2 Corrected neutral response

Equation (10) of 2006 gives the intrinsic perturbation in terms of the total redshifted
perturbation and `delta psi`. Equation (11) selects `delta psi` so that `q^T delta n=0`. In matrix
form the resulting local response is

```text
chi_CN = chi - (chi q) (q^T chi q)^{-1} (q^T chi).
```

This is the matrix content of the bracket in 2006 eq. (13). It obeys

```text
chi_CN q = 0,
q^T chi_CN = 0,
rank(chi_CN) <= 3.
```

After integration with `e^{-Phi}dV`, the full four-species paper `B` is therefore singular and
has the row identity stated after eq. (13). Its null vector represents an electrostatic-potential
zero-point change, not a physical composition response.

### 6.3 Reduced physical Hessian and susceptibility

Restrict the intrinsic energy Hessian to the source basis:

```text
H_y = S_y^T K S_y = partial g_y/partial y.
```

For a stable smooth phase, the corresponding local reduced number response is

```text
chi_y = H_y^{-1} = partial y/partial g_y.
```

It represents the same physical response as `chi_CN`:

```text
chi_CN = S_y chi_y S_y^T.
```

This identity is also a future validation target; no production matrix is formed here. In the
recommended baryon-density basis,

```text
H_x = T^T H_y T,
chi_x = T^{-1} chi_y T^{-T}.
```

The corrected reduced global response is then semantically

```text
Btilde = integral_core dV e^{-Phi} chi_y,
```

with the source's `(n_n,n_e,n_mu)` particle-number ordering. It is Layer B, not an EOS output.
Its inverse can enter corrected `Z` only after global integration. This is why neither
`integral(H)`, nor `inverse(integral(full B))`, nor an equilibrium-curve derivative is a valid
substitute.

### 6.4 Candidate derivative-authority assessment

| Authority | Source fit | Full charged null risk | Symmetry/stability | Later solve needed | Proposed disposition |
|---|---|---|---|---|---|
| Full four-species `K` | inverse orientation of 2005 eq. (10) | no null until projected, but exposes an unphysical charged direction | direct Hessian tests | local inverse plus projection | optional provider capability; not minimum API |
| Full four-species `chi` | exact 2005 eq. (10) input | corrected projection is rank three; full integrated corrected `B` singular | direct susceptibility symmetry | projection, reduction, global inverse | optional provider capability; must not expose a full paper-`B` inverse |
| Reduced `H_y` or `H_x` | exactly equivalent to 2006 reduced physical variables | charged direction absent | direct Hessian symmetry and stability | local reduced solve to obtain susceptibility | **recommended authority, pending Q3** |
| Reduced `chi_y` or `chi_x` | direct reduced global-integrand orientation | charged direction absent | inverse-response symmetry and stability | no local orientation conversion; later global inverse | viable source-aligned alternative, pending Q3 |

The recommendation is a semantic reduced Hessian, not a numerical differentiation or inversion
choice.

## 7. Proposed minimum generic contract

### 7.1 Required now

For any accepted cold local state:

| Item | Semantics and units |
|---|---|
| Independent coordinates | provisionally `(n_B,n_e,n_mu)` in `fm^-3`; explicit chart and held-fixed meaning |
| All species densities | `n_B,n_n,n_p,n_e,n_mu` in `fm^-3`, with exact dependent reconstruction |
| Intrinsic chemical potentials | `mu_n,mu_p,mu_e,mu_mu` in `MeV`; local, no `q_i psi`, no redshift, rest-mass convention declared |
| Imbalances | derived `eta_npe=mu_n-mu_p-mu_e`, `eta_npmu=mu_n-mu_p-mu_mu` in `MeV` |
| Equilibrium anchor | equilibrium state at `n_B`, including active-species/threshold status and applicable beta-equilibrium conditions |
| Local derivative | selected reduced object, orientation, basis, held-fixed variables, `MeV fm^3` for `H` or `fm^-3 MeV^-1` for `chi` |
| Domain/failure | valid density/composition range, phase identity, differentiability, threshold status; fail closed rather than silently cross a boundary |
| Provenance | EOS/model and revision, parameter/data source, component provenance, species and rest-mass convention, units, coordinate basis, lepton ownership |

At a muon threshold the provider must not pretend that a smooth four-active-species derivative
exists through the boundary. Below threshold, `n_mu=0` and the appropriate equilibrium
complementarity/appearance condition replaces a false requirement that both channel imbalances
vanish for an absent species. The exact threshold contract remains a later source-specific item.

### 7.2 Required for validation, not Phase-5A coefficient runtime

An off-equilibrium energy density `epsilon(x)` in `MeV fm^-3` is required to validate that the
chemical potentials and reduced Hessian are derivatives of one thermodynamic potential. The
2005/2006 global coefficients can consume a directly supplied local response without evaluating
off-equilibrium structure, so `epsilon(x)` need not be on their runtime hot path.

The validation relation in the recommended basis is

```text
d epsilon = g_x^T d x,
g_x=(mu_n,-eta_npe,-eta_npmu).
```

### 7.3 Future-generic, not minimum Phase 5A

Off-equilibrium `P(x)`, an explicit temperature coordinate, entropy, higher derivatives, and a
general extra-species state are valuable extensions but are not required to construct the cold
linear response used by the standard benchmark. The 2005 paper explicitly notes in section 2.2,
PDF p. 5, that away from beta equilibrium the pressure-density relation also depends on local
abundances; the 2006 calculation separately adopts the Cowling approximation after eq. (10),
PDF p. 3 / journal p. 570, for the small composition perturbation. Thus off-equilibrium pressure
is real thermodynamic information, but is not needed to re-solve the background structure when
forming the paper's linear susceptibility coefficients. At zero temperature, a provider that
supplies both `epsilon` and `P` can be checked against the corresponding thermodynamic identity;
a separate `P(x)` return is not made mandatory by this proposal.

Effective masses are outside this local susceptibility contract. They are nevertheless required
later for quantitative reproduction of the 2005 heat-capacity and reaction inputs (2005 section
3.1 and eq. (50)).

### 7.4 Unit chain

| Object | Proposed local/source unit | Later conversion boundary |
|---|---|---|
| `n_B,n_i` | `fm^-3` | stellar count integrals convert proper volume from `km^3` to `fm^3`; no such integral is implemented here |
| `Y_i` | dimensionless | `n_i=Y_i n_B` only on a state where both quantities share provenance (ADR-0001) |
| intrinsic `mu_i`, local `eta` | `MeV` | redshifted `eta^infinity=e^Phi eta` belongs to the later chemical/evolution contract |
| off-equilibrium `epsilon`, optional `P` | `MeV fm^-3` | the existing stellar background stores `P,epsilon` in geometric `km^-2`; conversion remains at the EOS/structure boundary |
| reduced chemical Hessian `H` | `MeV fm^3` | none before the local reduced response solve |
| reduced number susceptibility `chi` | `fm^-3 MeV^-1` | multiply by proper volume in `fm^3` and dimensionless `e^-Phi` |
| corrected reduced global `Btilde` | `MeV^-1` | its reduced inverse and paper `Z` coefficients have `MeV` units |

Natural-unit quantities in the current particle/RMF classes (`fm^-1`, derivative labels in
`fm^2`) are not interchangeable with this proposed unit contract without an explicit `hbar c`
conversion (`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:492-509`).

## 8. Ownership boundary

| Layer | Owns | Must not own |
|---|---|---|
| **A1 — off-equilibrium EOS provider** | charge-neutral local state, intrinsic `mu_i`, equilibrium anchor, selected local derivative, units/domain/provenance | electrostatic gauge, GR redshift, star integration, paper `B/Z/W`, ODE state |
| **A2 — local rotochemical thermodynamics adapter** | exact basis transforms; the 2006 electrostatic/charge-neutral projection; conversion from selected Hessian orientation to reduced susceptibility; source notation mapping | EOS model evaluation or global volume integration |
| **B — stellar rotochemical coefficients** | `e^{-Phi}dV` integration, corrected reduced global response, global baryon constraint, later paper `Z`; profile/EOS provenance | inventing local derivatives or asking the EOS to integrate the star |
| **C — rotational structure and fixed-baryon reduction** | existing Phase-4 response; later measure-complete particle response and sequence reduction under INV-09 | local EOS susceptibility or secular evolution |
| **D — secular evolution** | redshifted `eta` state, reaction restoration, heating/cooling, ODEs | local derivative reconstruction or structural re-integration |

This separation preserves the electrostatic physics without making `psi` a false observable EOS
field. A validation harness may explicitly reconstruct `delta psi` using 2006 eq. (11); the
production consumer needs the gauge-invariant reduced response.

## 9. EOS and data-product option comparison

| Option | Equilibrium only? | Independent composition? | Chemical potentials? | Derivative route possible? | Source-validated / compiled now? | Track R | Track P |
|---|---|---|---|---|---|---|---|
| Current cold DS(CMF)-1 CompOSE product | Yes: one cold beta-equilibrium trajectory | No; one placeholder `Y_q` point | Raw collective `Q3-Q5` only along the path; production API discards them | No transverse derivative is identifiable | Background TOV/rotation path validated and compiled; not an off-equilibrium source | No | Background only |
| Current finite-T DNS CMF product `(T,n_B,Y_q)` | No in strong charge fraction | One `Y_q` direction, but no independent `e/mu` split; different product | Raw `Q3-Q5`; current API stores only `Q2` | A future differentiated table could be possible, but no method or matching contract is validated | Entropy/heat-capacity path compiled/tested; off-equilibrium EOS response not validated | No | Candidate component only, insufficient as-is |
| A new CompOSE multidimensional product | Depends on requested table | Standard grid can supply `T,n_B,Y_q`; a composite lepton treatment is still needed for two lepton directions | Potentially collective chemical potentials and energy/pressure | Table, analytic, or hybrid route possible; none selected | Not present or authenticated for DS(CMF)-1 | Only if it reproduces a literature EOS | Plausible if source-matched and validated |
| Current `Model` / `SigmaOmegaRho` classes | Existing public flow solves equilibrium | Solver variables exist, but not a public arbitrary-composition thermodynamic state | Equilibrium rows include `mu_i`; derivative-labelled fields have no governed fixed-variable meaning | Underlying RMF model could in principle support analytic/AD/numerical derivatives | Compiled; not production-used or Phase-5 tested (`CompactStar/EOS/src/SigmaOmegaRho.cpp:82-253`, `:917-960`) | Not APR/BPAL without a new source implementation | Possible distinct research provider after governance/validation |
| `Fermi_Gas_Many` | Solves a beta-equilibrium free-gas composition | No generic off-equilibrium API; current pool is developmental | Low-level mutable particle potentials | Analytic free-gas response is possible, but current class is not validated authority | Compiled, dormant/developmental (`CompactStar/EOS/src/Fermi_Gas_Many.cpp:35-175`, `:300-370`) | Historical starting point only; a clean independent fixture is preferred | No realistic production EOS |
| Dormant `SigmaOmegaRho_npemu` | Equilibrium solver | No governed arbitrary-composition contract | Candidate internal values | Candidate Jacobian is an equilibrium-solver Jacobian, not the required response | Not in EOS/Physics CMake lists; stale include; untested | No | Historical evidence only |
| Noninteracting Fermi gas from the literature benchmark | Paper uses it as one EOS set | Analytic physical state can be defined | Analytic | Analytic derivatives available in principle | Published benchmark source, but no authenticated CompactStar implementation | Strong first falsifier / simple published star | Not realistic production |
| BPAL family | Published analytic nuclear EOS family | Energy as density/composition function is intended by the 2005 method | Derivable from the authenticated model | Analytic derivatives possible in principle | Primary source/parameters not in the shared library; no current implementation | Realistic intermediate candidate | Different from canonical CMF |
| APR A18 + delta-v + UIX* | 2005 reference EOS; equilibrium background plus off-path energy fit required | Energy versus baryon density/proton fraction is the stated route | Derivable from the authenticated fit/model | Analytic derivatives possible in principle | Primary source/fit/data not in the shared library; no current implementation | **Proposed closing benchmark** | Different from canonical CMF |

“Possible in principle” is not “available.” No option is promoted to authority by this table.

## 10. Literature reproduction proposal

### 10.1 Options assessed

| Option | Fidelity | Difficulty / independence | Use as falsifier | Published-output reach |
|---|---|---|---|---|
| **R1 free Fermi gas first** | Low realism, but one of the 2005 EOS sets | Lowest; analytic and independent | Excellent for units, signs, constraints, derivatives | Can reproduce the paper's free-gas structural case, not the reference realistic curves alone |
| **R2 one BPAL EOS first** | Realistic analytic intermediate | Requires exact BPAL variant and primary parameters; more independent than CMF | Strong test of interacting composition response | Can target one published 2005 family member |
| **R3 A18 + delta-v + UIX* first** | Highest direct fidelity to the paper's reference EOS | Highest source/data and phase-transition burden | Strong final benchmark but poor first debugger | Required for the reference 2005/2006 coefficient and thermal comparisons |
| **R4 staged ladder** | Progresses from exact local checks to the reference calculation | More work overall, but localizes failures | Strongest | Supports both simple and realistic published benchmarks |

### 10.2 Recommendation

Recommend **R4**, pending owner adjudication:

1. analytic relativistic free-electron and free-muon components;
2. an analytic charge-neutral toy npe and npe-mu energy functional with exact reduced matrices;
3. the 2005 noninteracting Fermi-gas stellar EOS case;
4. one authenticated BPAL variant if it supplies a useful intermediate before APR;
5. the reference A18 + delta-v + UIX* EOS with the 2006 corrected response; and only then
6. Track P on a source-matched CompactStar CMF provider.

The toy cases are validation fixtures, not substitutes for a published realistic result. The
current DS(CMF)-1 trajectory cannot close Track R.

### 10.3 Missing primary sources and numerical data

The shared library does not currently contain the following material needed for a defensible
quantitative reproduction. Nothing was downloaded or fabricated in Phase 5A-1.

1. **APR/A18 + delta-v + UIX***: Akmal, Pandharipande & Ravenhall, *Phys. Rev. C* **58**, 1804
   (1998), including the exact table/analytic-fit definition used by Fernández & Reisenegger,
   its valid density range, the normal-to-pion-condensed phase handling, and enough data to
   reproduce equilibrium composition and the transverse energy derivatives.
2. **Crust completion for the 2005 star models**: Pethick, Ravenhall & Lorenz, *Nucl. Phys. A*
   **584**, 675 (1995), and Haensel & Pichon, *A&A* **283**, 313 (1994), plus the exact joins and
   core-integration boundary used by the benchmark.
3. **BPAL**: Prakash, Ainsworth & Lattimer, *Phys. Rev. Lett.* **61**, 2518 (1988), and the
   Prakash et al. *Physics Reports* **280**, 1 (1997) labelling/parameter definitions for the
   exact five variants used in 2005.
4. **Effective-mass prescription**: the exact APR/PAL effective-mass source used for 2005
   eq. (50), cited there via Page et al. (2004), plus any reaction-rate inputs needed to reproduce
   the later thermal curves.
5. **Machine-readable benchmark outputs or author code** for the 2005 `I_Omega`, `Z`, `W`, and
   temperature/imbalance curves, and the 2006 corrected `Z` curves. The papers provide figures
   and qualitative/fractional statements but not a complete numerical fixture. If such data do
   not exist, an owner-approved digitization and uncertainty protocol will be required before
   setting acceptance tolerances.
6. **Free-gas reference details** if exact replication of the paper's complete free-gas star is
   required, including the cited Shapiro & Teukolsky formulation and its composition/maximum-density
   convention.

The exact implementation provenance used by the 2005 authors is a scientific input. Recreating a
plausible APR or BPAL formula without authenticating that lineage would not be a reproduction.

## 11. Validation ladder and falsified defects

No tolerance is assigned here. Each later implementation must derive tolerances from analytic
precision, source tabulation/figure uncertainty, conditioning, and measured convergence before
the check is run.

| ID | Future experiment | Pass meaning | Defect falsified |
|---|---|---|---|
| **V1 — units/signs** | Check declared `fm^-3`, `MeV`, `MeV fm^-3`, `MeV fm^3`, `fm^-3 MeV^-1`; rest-mass inclusion; species order; `eta_npe`/`eta_npmu` signs; local versus redshifted status. | Every field and transform is dimensionally and semantically explicit. | Hidden natural-unit, cgs/geometric, redshift, rest-mass, species-index, or sign mismatch. |
| **V2 — exact neutrality** | Generate allowed states in each supported chart and reconstruct all species. Check `n_p=n_e+n_mu` by construction and positivity/domain failure. | No accepted state leaves the physical charge-neutral manifold. | Constraint drift, wrong dependent species, or invalid accepted state. |
| **V3 — beta equilibrium** | Evaluate `EquilibriumStateAt(n_B)` across the domain; check active-channel imbalance equalities and threshold/absent-species conditions; match the provider's governed equilibrium background. | The off-equilibrium provider and stellar-background EOS meet at the same physical curve. | Wrong equilibrium anchor, mass convention, model mismatch, or threshold rule. |
| **V4 — free leptons** | Compare electron/muon `mu_l(n_l)` and analytic first derivatives to the exact zero-temperature relativistic Fermi-gas expressions, including threshold behavior. | Lepton component is independently right. | Missing `hbar c`, wrong rest mass/Fermi momentum, wrong units, or derivative orientation. |
| **V5 — analytic toy EOS** | Use smooth npe and npe-mu energy functionals with closed-form `mu`, `H`, `chi`, equilibrium, and basis transforms. | The complete local contract reproduces an independent exact oracle. | Held-fixed ambiguity, wrong indices, incorrect constraint, or wrong solve target. |
| **V6 — Maxwell/symmetry** | Check `H=H^T` (or corresponding `chi` symmetry) in every smooth phase and compare mixed derivatives from independent routes. | Returned response derives from one thermodynamic potential where theory requires. | Non-integrable interpolation/model response or mismatched differentiation paths. |
| **V7 — finite perturbation** | At fixed declared coordinates and within one phase, compare finite `delta g` with the linear response under shrinking perturbations and establish convergence. | The returned Jacobian is the local derivative it claims to be. | Wrong numerical derivative, fixed-variable convention, sign, or unreported nonsmoothness. |
| **V8 — basis equivalence** | Transform the same state/perturbation among `(n_B,n_e,n_mu)`, `(n_n,n_e,n_mu)`, and valid fraction coordinates; compare `delta n`, `delta epsilon`, imbalances, `H`, and `chi`. | Software coordinates do not change the physics. | Nonlinear-transform omission, density/fraction unit mix, or source-basis mismatch. |
| **V9 — corrected full null mode** | In a test-only toy star/provider, construct the 2006 eq. (13) charge-projected full paper `B`; check symmetry, the charge null vector, and `B_pj=B_ej+B_muj`. | The electrostatic correction and local neutrality are present. | Silent reversion to 2005 eq. (12), incorrect charge vector, or illegal full inverse. |
| **V10 — reduced system** | Remove the source proton row/column or use the exact reduced basis; establish the expected `3 x 3` rank, symmetry/stability, conditioning record, and equivalence to the full projected action over the valid domain. | Only physical modes remain and the later solve is defined. | Wrong reduction, hidden gauge mode, threshold crossing, or unstable data. |
| **V11 — corrected published coefficients** | With an authenticated published EOS and numerical reference, reproduce corrected 2006 `Z` values or the explicitly stated changes: `Z_np` nearly unchanged and `Z_npe/Z_npmu` changes reaching roughly 10-20% at the largest masses (2006 section 5/Figure 1). | Local response, GR integration, reduced inverse, and baryon relation agree with the corrected benchmark. | Wrong redshift/volume, matrix projection/reduction, source mapping, or particle-number ordering. |
| **V12 — non-superfluid evolution benchmark** | Only after Layers B-D exist, reproduce selected 2005/2006 temperature and imbalance curves for the authenticated realistic EOS and initial/spin parameters. | The standard mechanism works end to end before superfluidity or BNV. | Coefficient, spin drive, reaction, heating/neutrino, redshift, integration, or ODE defect. |

V1-V8 primarily validate Layer A. V9-V11 validate the Layer A-to-B interface. V12 is deliberately
deferred to Layer D.

## 12. Scientific invariants independent of algorithm

Any later implementation method must satisfy:

- exact/tolerance-accounted local charge neutrality on every accepted state;
- explicit local intrinsic versus total/redshifted chemical-potential semantics;
- recovery of the matching beta-equilibrium curve and active-species thresholds;
- thermodynamic integrability and symmetry where the EOS is smooth;
- physical stability/rank properties on the declared phase domain;
- derivative convergence under independent finite perturbations;
- exact coordinate-basis transformation equivalence;
- the expected charged null mode only in the corrected full projected response;
- a nonsingular reduced physical system in its declared domain;
- explicit failure or phase metadata at thresholds and non-differentiable transitions; and
- data/model, units, component, and stellar-profile provenance sufficient to prevent Track R and
  Track P from being silently mixed.

These requirements do not select finite differences, automatic differentiation, analytic
differentiation, spline differentiation, a matrix library, a solver, or an interpolation scheme.

## 13. Notation hygiene

The repository's historical structural `A_i`, `B_i`, and `Z_i` are not the papers' thermodynamic
`B_ij`, `Z_np`, `Z_npe`, `Z_npmu`, `W_npe`, or `W_npmu`
(`docs/validation/PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:132-144`, `:461-479`).

Proposed descriptive names, subject to later API design, are:

| Role | Descriptive name |
|---|---|
| Allowed local density state | `ChargeNeutralCompositionState` |
| Intrinsic local chemical potentials | `IntrinsicChemicalPotentials` |
| Recommended local derivative | `ChargeNeutralChemicalHessian` |
| Its inverse-orientation response | `ReducedNumberSusceptibility` |
| Corrected, volume-integrated reduced response | `IntegratedThermodynamicResponse` |
| Paper `Z` coefficients | `RotochemicalImbalanceCoupling` |
| Equilibrium particle-number spin drive | `SpinCompositionDrive` |
| Fixed-central-density / sequence particle response | `FixedCentralDensityParticleResponse` / `SequenceParticleResponse` |

New code must not introduce an unqualified bare `B`, `B_i`, `B_ij`, or `Z_i` for these new
objects. Documentation may retain paper notation only when prefixed by “paper” or otherwise
qualified by its source role.

## 14. Owner questions and recommendations

| ID | Unresolved decision | Recommendation | Architectural/scientific consequence |
|---|---|---|---|
| **Q1** | Which realistic published EOS closes Track R? | Make A18 + delta-v + UIX* the closing benchmark; use free gas/toy first and BPAL only as an authenticated intermediate. | Requires primary APR/crust/effective-mass and numerical benchmark acquisition; prevents CMF substitution. |
| **Q2** | Which canonical local coordinates? | `(n_B,n_e,n_mu)`, with exact transformation to the 2006 `(n_n,n_e,n_mu)` basis as V8. | Natural standard EOS axis and direct `-eta` conjugates; source equivalence remains testable. |
| **Q3** | Which local derivative orientation does the EOS own? | Reduced charge-neutral chemical Hessian, with the local adapter deriving the reduced number susceptibility. | Gauge-free and symmetry-testable; later algorithm governance is required for differentiation/solve. |
| **Q4** | How are leptons owned? | One unified consumer-facing npe-mu provider, allowed to compose an independently validated analytic free-Fermi lepton component with a hadronic provider. | One coherent state for consumers; component provenance and cross-derivative semantics must be explicit. |
| **Q5** | What supplies Track P transverse thermodynamics? | Prefer a source-matched DS(CMF)-family multidimensional product; if unavailable, stop before treating a distinct RMF model as the canonical EOS. | Determines whether Track P preserves DS(CMF) identity or becomes a separately governed EOS study. |
| **Q6** | Cold-only or temperature-aware first contract? | Cold-only v1 with an explicit future extension point. | Matches the standard benchmark and avoids an ungoverned held-fixed temperature/entropy variable. |

These are recommendations only. ADR-0010 remains PROPOSED and its Decision is
`PENDING OWNER ADJUDICATION`.

## 15. Remaining unresolved scientific/data questions

The six owner questions above are the high-leverage architecture choices. The following details
remain later source-specific requirements rather than extra owner ballots:

- exact rest-mass and collective-potential conventions for mapping any CompOSE `Q3-Q5` data to
  individual intrinsic `mu_i`;
- the active-species and density-domain rule for applying an npe-mu contract to a DS(CMF)-1
  background that contains hyperons, deltas, or quarks at some densities;
- the muon-appearance and phase-boundary differentiability/complementarity contract;
- representation of first-order phase transitions in local response and in later
  measure-complete global particle/susceptibility integrals;
- the exact core boundary over which each paper coefficient is integrated;
- acquisition and authentication of the missing APR/BPAL/crust/effective-mass sources and
  machine-readable coefficient/evolution benchmarks; and
- INV-11's eventual redshifted evolution-state units, ordering, and owner. Local `eta` semantics
  proposed here do not accept that separate evolution contract.

No ambiguity above is resolved by guessing. A source-specific implementation task must stop if
its chosen provider cannot state these semantics.

## 16. Migration and next governed boundary

If the owner accepts ADR-0010 later, the next task should design and validate one provider
implementation plan against the accepted Q1-Q6. It must separately select the data/model source
and numerical derivative algorithm; this record does not do so.

The current cold EOS, `StarProfile`, Phase-4 rotation responses, production files, tests,
baselines, and external data/literature remain unchanged in Phase 5A-1. The proposed contract does
not authorize Layer B coefficient construction, Layer C particle-number reduction, or Layer D
evolution.

## 17. Scope closure

**NO THERMODYNAMIC EOS IMPLEMENTATION IN PHASE 5A-1**

No EOS derivative, interpolation, matrix inversion, rotochemical coefficient, ODE, reaction rate,
heating term, superfluid refinement, BNV path, test, scientific baseline, or data product was
implemented. ADR-0010 is **PROPOSED**, not ACCEPTED.
