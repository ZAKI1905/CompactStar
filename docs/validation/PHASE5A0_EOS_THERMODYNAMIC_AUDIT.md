# Phase 5A-0 EOS / thermodynamic capability audit

**Date:** 2026-09-04

**Starting SHA:** `59526c56122dbbb8c0a8ef808bdf627453d99c3a`

**Branch:** `analysis/phase5a-eos-thermodynamic-audit`

**Worktree:** `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-eos-audit`

**Scope:** source-backed inventory only; no EOS, rotochemical, evolution, test, baseline, or
literature implementation change.

## 1. Result first

CompactStar's production stellar path knows a **one-dimensional, cold, beta-equilibrium
background trajectory**. For the canonical DS(CMF)-1 stars it can evaluate and carry

- `epsilon(P)`, `n_B(P)`, and the dimensionless equilibrium fractions `Y_i(P)`;
- the barotropic derivative `d epsilon / dP`; and
- the TOV metric and the Phase-4 fixed-central-density rotational response.

It cannot evaluate the independent-composition derivatives required by the corrected
rotochemical formalism. In particular, no current production API supplies either
`partial mu_i / partial n_j` at fixed independent remaining state variables or
`partial n_i / partial mu_j`. The cold table's raw CompOSE files contain collective chemical
potentials only **along the beta-equilibrium path**; the production `.eos` reader discards those
chemical potentials. A separate finite-temperature CMF data product has an independent `Y_q`
axis, but the current API retains only entropy `Q2`, and the available files do not establish the
independent electron/muon composition coordinates needed for the corrected npe-mu system.

**Minimum-gap classification: D — a new, governed thermodynamic EOS extension/data product is
scientifically required.** This classification does not choose the EOS model, derivative method,
interpolation, reduced-coordinate basis, or inversion algorithm. Those are owner/source decisions
for a later task.

## 2. Start and scope authentication

The canonical checkout was authenticated at
`/Users/keeper/Documents/CompactStar/repo/CompactStar` on `master`. Local `HEAD`, cached
`origin/master`, and live `origin/master` were all
`59526c56122dbbb8c0a8ef808bdf627453d99c3a`; the tree was clean. The dedicated branch and worktree
above were created from that exact commit and authenticated clean before investigation.

The governing Phase-4 record states that the Phase-4 prerequisite is discharged but particle-number
response, baryon-conserving reduction, and chemical evolution remain Phase-5 work
(`docs/validation/PHASE4_CLOSEOUT.md:201-217`). The exact starting state was:

> PHASE 4 ROTATION CORRECTNESS COMPLETE — PHASE-5 STRUCTURAL INTERFACE RATIFIED
>
> PHASE 5 NOT YET BEGUN

This audit enters only Phase 5A-0 documentation/reconnaissance. It does not change the Phase-4
status and does not activate Phase-5 physics.

## 3. Source authority and source ledger

The shared library root is `/Users/keeper/Documents/CompactStar/literature`. Its `README.md` and
`catalog.tsv` were read before the papers. All catalog paths are relative to that root and resolve;
all 22 entries in `SHA256SUMS.txt` verified. The library is external to Git and was not modified.

Authority order used here:

1. `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`
   — **primary authority** for electrostatic perturbations, local charge neutrality, the singular
   full thermodynamic matrix, and the reduced inverse.
2. `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` — **primary non-superfluid
   reproduction benchmark**, except where corrected by the 2006 paper.
3. `1995-Reisenegger-Rotochemical-Heating.pdf` and
   `2020-Yanagi-NS-Therm-Thesis.pdf` — supporting material only.

The gravitochemical paper was not used as authority for the correction. The following ledger ties
every literature-derived requirement used by this audit to an exact local PDF, location, and role.

| Requirement or claim | Exact PDF filename | Equation / location | Source role |
|---|---|---|---|
| Chemical imbalances are `eta_npe = delta mu_n - delta mu_p - delta mu_e` and the analogous muon channel | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (7)-(8), PDF p. 4, section 2.1 | 2005 benchmark |
| Redshifted intrinsic chemical-potential perturbations are spatially uniform in the 2005 treatment | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (9), PDF p. 4, section 2.1 | 2005 benchmark; electrostatic completion superseded |
| Linear local response uses `delta n_i = sum_j (partial n_i / partial mu_j) delta mu_j`, evaluated at beta equilibrium | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (10), PDF p. 4, section 2.1 | 2005 benchmark |
| Integrated particle-number response and the original thermodynamic matrix | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (11)-(12), PDF p. 4, section 2.1 | 2005 benchmark; matrix definition superseded |
| The 2005 calculation naively inverts the full matrix | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (13), PDF p. 4, section 2.1 | superseded by 2006 |
| Reaction and changing-equilibrium particle numbers drive `delta N_i` | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (14)-(17), PDF pp. 4-5, section 2.1 | 2005 benchmark |
| Rotating-star scalar particle count and its proper-volume/metric/composition terms | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (24), PDF p. 6, section 2.2 | 2005 benchmark; structural interface |
| Equilibrium particle-number driving coefficient `I_{Omega,i}` | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (30), PDF p. 7, section 2.2 | 2005 benchmark; Layers B/C |
| Equilibrium EOS, composition, adiabatic index, and off-path energy versus baryon density and proton fraction were separate numerical inputs | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | section 3.1, PDF pp. 8-9 | 2005 benchmark |
| Effective masses enter the specific-heat/reaction benchmark, not the corrected susceptibility identity itself | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eq. (50), PDF p. 12, section 3.4, and section 3.1 | 2005 benchmark |
| Standard non-superfluid evolution equations and original `Z`, `W` definitions | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | eqs. (51)-(58), PDF pp. 13-14, sections 3.5-3.6 | ODE benchmark; eqs. (54)-(56) superseded |
| Diffusive equilibrium is governed by total redshifted chemical potential including `q_i psi` | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (2), PDF p. 2 / journal p. 568, section 2 | corrected 2006 authority |
| Local charge neutrality is imposed on the densities | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (4), PDF p. 2 / journal p. 568, section 2 | corrected 2006 authority |
| Observable beta-channel imbalances do not contain the electrostatic potential because charges cancel | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (9), PDF p. 4 / journal p. 570, section 3 | corrected 2006 authority |
| Perturbation of total redshifted chemical potential includes intrinsic, electrostatic, and metric terms | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (10), PDF p. 4 / journal p. 570, section 4 | corrected 2006 authority |
| Charge neutrality determines the electrostatic-potential perturbation | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eq. (11), PDF p. 4 / journal p. 570, section 4 | corrected 2006 authority |
| Corrected integrated response uses the charge-projected `B_ij` | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eqs. (12)-(13), PDF p. 4 / journal p. 570, section 4 | corrected 2006 authority |
| The corrected full four-species matrix is singular, with a charge/gauge null mode | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | text immediately after eq. (13), PDF pp. 4-5 / journal pp. 570-571, section 4 | corrected 2006 authority |
| Proton row/column is eliminated and a symmetric invertible `3 x 3` reduced matrix over `n,e,mu` is inverted | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | text between eqs. (13)-(14), PDF p. 5 / journal p. 571, section 4 | corrected 2006 authority |
| Global baryon conservation further reduces particle-number perturbations to observable imbalances | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eqs. (14)-(15), PDF p. 5 / journal p. 571, section 4 | corrected 2006 authority |
| Corrected `Z_np`, `Z_npe`, `Z_npmu` and `W` expressions | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | eqs. (16)-(19), PDF p. 5 / journal p. 571, section 4 | corrected 2006 authority |
| 2006 replaces 2005 eqs. (9), (12), (13), (54)-(56); quasi-steady results remain unchanged and transient changes are small | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | sections 5-6, PDF pp. 5-6 / journal pp. 571-572 | corrected 2006 authority |
| Original rotochemical mechanism uses an independent composition coordinate under charge neutrality | `1995-Reisenegger-Rotochemical-Heating.pdf` | section 2, PDF p. 4: `x=n_p/n`, with `n_p=n_e` | supporting |
| Later review confirms the charge-neutral identity makes the full `B_ij` noninvertible and requires a submatrix | `2020-Yanagi-NS-Therm-Thesis.pdf` | eq. (4.7) and surrounding text, PDF pp. 58-60, section 4.1.1 | supporting |

The CompOSE variable interpretation was checked independently against the official
[CompOSE Reference Manual v3.00](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf):
section 3.4, eqs. (3.15), (3.26), and section 4.2, eqs. (4.5)-(4.9), (4.18), manual pp. 14, 20,
33, and 36. It is a data-format authority, not a rotochemical authority.

## 4. Corrected thermodynamic requirements before any ODE

To avoid the repository's historical naming collision, this report uses:

- **paper `B_ij`** for the global, integrated thermodynamic susceptibility of 2005/2006;
- **paper `Btilde`** for the corrected reduced `3 x 3` matrix of 2006; and
- **legacy structural `B_i`** only for the repository candidate's equilibrium-sequence derivative
  `(partial N_i / partial epsilon_c)_Omega`.

They are not interchangeable.

| Needed object | Role | Primary source | Layer dependency |
|---|---|---|---|
| Equilibrium local `n_i(r)` | Background species populations and scalar counts | 2005 eqs. (10), (20)-(24) | A, then B |
| Intrinsic local `mu_i` and total redshifted `mu_i^infinity` | Defines diffusive equilibrium and observable imbalances | 2006 eqs. (2), (9)-(10) | A |
| Local susceptibility `chi_ij = partial n_i / partial mu_j` at equilibrium, with the correct fixed independent variables | Maps small chemical-potential departures to density departures | 2005 eq. (10); corrected projection in 2006 eqs. (11)-(13) | A |
| Equivalently, a governed inverse response `K_ij = partial mu_i / partial n_j` on an independent coordinate basis | May encode the same local response if coordinates and constraints are explicit | equivalent form of 2005 eq. (10), subject to 2006 charge reduction | A |
| Local charge constraint and charges `q_i` | Determines `delta psi` and removes the unphysical charged mode | 2006 eqs. (4), (11) | A |
| Corrected paper `B_ij` | Volume-integrated, redshifted, charge-projected response | 2006 eqs. (12)-(13) | A + B |
| Null/singular charged mode | Forbids inversion of the full four-species paper `B_ij` | 2006 text after eq. (13) | A + B |
| Reduced paper `Btilde` and `Btilde^{-1}` | Physical nonsingular response in variables `tilde mu_n=delta mu_n` and `tilde mu_l=delta mu_p+delta mu_l`, after eliminating the proton row/column | 2006 text between eqs. (13)-(14) | A + B |
| Baryon-number constraint | Removes the remaining globally conserved baryon-number perturbation in the imbalance relation | 2006 eqs. (14)-(15) and footnote 4 | B |
| Paper `Z_np`, `Z_npe`, `Z_npmu` | Converts particle-number departures into observable imbalances | 2006 eqs. (15)-(18) | A + B |
| `I_{Omega,i}` / equilibrium particle-number rotational response | Describes changing equilibrium populations under spin down | 2005 eq. (30) | B + C + fixed-baryon reduction |
| Paper `W_npe`, `W_npmu` | Couples the equilibrium structural response to `2 Omega dotOmega` | 2006 eq. (19), retaining the role of 2005 eqs. (57)-(58) | A + B + C + fixed-baryon reduction |

Only after these quantities and their conventions exist can Layer D form the `T_infinity`,
`eta_npe`, and `eta_npmu` ODEs. This audit does not form them.

## 5. Current EOS architecture map

### 5.1 Production path and interpolation

`CompactStar::Model` is the abstract analytic/table-generator base. Its public contract is
`EDens(rho)`, `Press(rho)`, `EOSRow(rho)`, `EOSHeader()`, `FindEOS()`, and `ExportEOS()`
(`CompactStar/EOS/Model.hpp:54-85`). The default standardized output is energy-density-like mass
density in `g cm^-3`, pressure in `dyne cm^-2`, and baryon density in `fm^-3`
(`CompactStar/EOS/src/Model.cpp:83-96`). All concrete source files listed below are compiled into
the library (`CompactStar/EOS/CMakeLists.txt:1-37`), but compilation is not evidence of production
use.

The actual canonical star path is `TOVSolver`, not a polymorphic `Model` call. `ImportEOS()` reads a
standardized tab-separated file: the first three columns become `epsilon`, `P`, `n_B`, and all
remaining columns are copied verbatim as generic extras
(`CompactStar/Core/src/TOVSolver.cpp:597-640`). It constructs monotone GSL Steffen splines for
`epsilon(P)`, `n_B(P)`, and every extra (`:667-715`). The public queries are:

| Class / API | Input | Output and units | Composition domain | Production / tests / status |
|---|---|---|---|---|
| `TOVSolver::GetEDens(P)` | `P`, `dyne cm^-2` | table `epsilon`, `g cm^-3` | equilibrium path | Production; extensively used and tested (`CompactStar/Core/src/TOVSolver.cpp:1054-1059`) |
| `TOVSolver::GetRho(P)` | `P`, `dyne cm^-2` | `n_B`, `fm^-3` | equilibrium path | Production and tested (`:1551-1556`) |
| `TOVSolver::GetRho_i(P)` | `P`, `dyne cm^-2` | generic extra values; for canonical DS file these are dimensionless `Y_i`, despite the name | equilibrium path | Production and tested through stars/thermal code (`:1565-1583`) |
| `TOVSolver::GetEDensDeriv(P)` | `P`, `dyne cm^-2` | geometric `d epsilon / dP`, dimensionless | barotropic equilibrium path | Production, EOS-owned and tested (`:1069-1104`; `tests/eos/eos_derivative_cmf.cpp:100-170`) |

The TOV solve samples exactly those queries into each `TOVPoint`
(`CompactStar/Core/src/TOVSolver.cpp:2580-2592`). There is no `mu_i`, Fermi momentum, effective
mass, arbitrary-composition state, composition Jacobian, or susceptibility in `TOVSolver`.

### 5.2 Concrete EOS and thermodynamic classes

| Class | Relevant APIs; inputs -> outputs | Thermodynamic meaning / units | Equilibrium vs arbitrary composition | Use / test / disposition |
|---|---|---|---|---|
| `Fermi_Gas`, `Polytrope`, `CoulombLattice`, `SigmaOmega` | `EDens(n_B)`, `Press(n_B)` | analytic barotropes; standardized export converts to the base header's cgs/fm units | one density coordinate; no composition API | Compiled; no canonical Phase-5 star/test caller found; generic but dormant for this purpose |
| `Fermi_Gas_Many` | equilibrium solver followed by `EDens(n_B)`, `Press(n_B)`, `EOSRow(n_B)` | free-gas energy/pressure and equilibrium fractions; internal momenta/chemical potentials use natural `fm^-1` conventions | explicitly solves beta equilibrium and charge/baryon relations, not arbitrary composition (`CompactStar/EOS/src/Fermi_Gas_Many.cpp:51-156,203-279`) | Compiled; no caller/test found; developmental/dormant |
| `SigmaOmegaRho` | `SetRho(n_B)`, `SolveEq[Jacob]()`, `EDens`, `Press`, `EOSRow` | RMF equilibrium; row exports cgs `epsilon/P`, `n_B`, fractions, `mu_i` in MeV, and derivative-labeled columns in natural `fm^2` | `SetXRhos` reconstructs dependent baryon/lepton densities from conservation and `Eq` imposes equilibrium (`CompactStar/EOS/src/SigmaOmegaRho.cpp:82-145,148-253`) | Compiled, but no canonical test/production caller; dormant generator |
| `SigmaOmegaRho::Jacob` | solver vector of scalar field plus selected fractions -> residual Jacobian | derivatives of the **nonlinear equilibrium residuals**, not the paper susceptibility at fixed other physical densities (`:365-486`) | equilibrium solve only | Untested for Phase 5; not a rotochemical derivative API |
| `SigmaOmegaRho::EOSRow` | `n_B`, `fm^-3` -> equilibrium row | includes `Y_i`, `mu_p,mu_n,mu_e,mu_mu` in MeV and columns labelled `d_Mu_* (fm^2)` (`:917-960`) | values/derivatives generated on the solved equilibrium state; no documented fixed-variable contract | Compiled but unused/untested; cannot be promoted to corrected susceptibility authority |
| `SigmaOmegaRho_nstar` | same shape as `SigmaOmegaRho`, with a project-specific additional neutron-star species | same natural/cgs export conventions | equilibrium solver only | Compiled; project-specific, no caller/test found; dormant |
| `Particle`, `Baryon`, `Lepton` | `SetRho`, `GetRho`, `k`, `Set_Mu`, `Mu`, `Set_Rho_From_Mu`, `DMu_DRho`, plus baryon `Mu(gsig)`, `Dmu(...)` | low-level mutable natural-unit primitives; density `fm^-3`, momentum/chemical potential internally `fm^-1`, free-particle `dmu/dn` in `fm^2` | individual values can be set, but no coherent public arbitrary-composition EOS state or constrained derivative basis exists (`CompactStar/EOS/Particle.hpp:53-166`; `src/Particle.cpp:63-153`; `src/Baryon.cpp:179-192,256-290`) | Compiled as generator internals; no direct Phase-5 production/test contract |
| `CompOSE_EOS` | raw `ImportGrid/ImportThermo/ImportCompo/ImportMicro`; standardized `ImportEOS`; data getters; `GetFermiE` | can build standardized `epsilon/P/n_B/Y_i` and, if separate auxiliary files exist, effective masses/potentials and derived Fermi energies (`CompactStar/EOS/CompOSE_EOS.hpp:74-229`) | follows supplied grids; no public arbitrary-composition chemical-potential or susceptibility evaluator | Compiled; used by some BNV/LightDM code, not by canonical star loading; no Phase-5 thermodynamic tests |
| `CompOSE_Thermo` | `Q2(T,n_B,Y_q)`, `dQ2dT`, heat-capacity helpers | `T` MeV, `n_B` fm^-3, `Y_q` dimensionless; returns entropy/heat-capacity quantities (`CompactStar/EOS/CompOSE_Thermo.hpp:7-28,93-182`) | arbitrary supplied `Y_q` for **entropy only**; does not impose beta equilibrium | Production thermal support and tested, but not an EOS pressure/chemical-potential surface |

`CompOSE_EOS::ImportThermo` consumes pressure `Q1` and energy `Q7` but not `Q3-Q5` chemical
potentials (`CompactStar/EOS/src/CompOSE_EOS.cpp:270-317`). `CompOSE_Thermo` parses all `Q1-Q7`
fields but stores only `Q2` (`CompactStar/EOS/src/CompOSE_Thermo.cpp:278-355`). Thus chemical
potentials physically present in raw `eos.thermo` are an **API exposure gap**, while missing
independent composition directions are a **physics/data gap**.

`CompOSE_EOS::ImportEOS` additionally expects `_m_eff.micro`, `_V.micro`, and `_U.micro` files
(`CompactStar/EOS/src/CompOSE_EOS.cpp:1217-1229`). Its electron and muon Fermi energies are
analytic free-Fermi expressions from `Y_i n_B`; baryon expressions additionally require the
effective-mass/potential tables (`:1242-1252,1293-1308,1321-1343`). Those auxiliary files do not
exist in either current CMF directory, and this route is not the canonical star loader.

## 6. EOS files actually used and their columns

Canonical full tests receive external tables through `COMPACTSTAR_EOS_DATA_ROOT`
(`tests/CMakeLists.txt:179-218`). TOV, baryon-number, and rotation tests use
`DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos`; thermal/cooling tests pair that cold structure with
`DNS-CMF-Hadronic-with-electrons/eos.{t,nb,yq,thermo}` (for example
`tests/thermal/heat_capacity_real_star.cpp:21-22,90-112`).

### 6.1 Cold DS(CMF)-1 with crust

External directory:
`/Users/keeper/Documents/CompactStar/data/compose/DS-CMF-1-with-crust`

The standardized production file has 1,191 data rows plus a header:

`e(g/cm^3), p(dyne/cm^2), rho(1/fm^3), 500, 501, 502, 0, 1, 11, 10, 100, 112, 111, 110, 121, 120, 23, 22, 21, 20`.

The numeric labels map to electron `0`, muon `1`, neutron `10`, proton `11`, hyperons/deltas, and
quarks (`CompactStar/EOS/src/CompOSE_EOS.cpp:23-44`). Inspection of representative rows and the
nucleonic-core residuals confirms the ADR-0001 interpretation that these extras are fractions:
where the row is nucleonic,
`Y_n+Y_p=1` and `Y_p=Y_e+Y_mu` to file precision. At `n_B=0.16 fm^-3`, for example,
`Y_e=0.046448922`, `Y_mu=0.0045602062`, `Y_p=0.051009128`, and `Y_n=0.94899087`.

The raw axes have exactly one `T=0` point, one placeholder `Y_q=0` point, and 1,191 `n_B` points
from `1e-7` to `3.03 fm^-3`; therefore the raw cold product is also one-dimensional. Its
`eos.thermo` physically contains `Q1-Q7`. Under CompOSE v3.00 eq. (4.18), `Q3-Q5` are the
collective baryon, charge, and lepton chemical potentials normalized by the neutron mass. These
permit reconstruction of **equilibrium-path** species chemical potentials using each particle's
conserved charges (CompOSE eq. (3.26)). At `n_B=0.16 fm^-3`, the raw row has
`Q3=0.046910179`, `Q4=-0.12688433`, `Q5=0`, and header neutron mass `938.9187125 MeV`.

It does **not** provide an independent composition coordinate, an `eos.compo`, an `eos.micro`, or
the `_m_eff.micro`, `_V.micro`, `_U.micro` products expected by `CompOSE_EOS::ImportEOS`. The
standardized production `.eos` contains no `mu_i`, Fermi momentum, or effective-mass column.

### 6.2 Finite-temperature DNS CMF with electrons

External directory:
`/Users/keeper/Documents/CompactStar/data/compose/DNS-CMF-Hadronic-with-electrons`

The axes contain 81 temperatures (`0` to `160 MeV`), 301 baryon densities (`0.01` to
`3.01 fm^-3`), and 54 charge fractions (`0` to `0.53`), yielding 1,316,574 thermodynamic rows.
The raw `eos.thermo` has `Q1-Q7`, including collective `Q3-Q5`, on this `(T,n_B,Y_q)` grid.
This is genuine independent-`Y_q` table information, but:

- current `CompOSE_Thermo` exposes only `Q2`/heat capacity;
- the directory has no `eos.compo`, `eos.micro`, or auxiliary standardized micro files;
- it is a distinct CMF data product from the cold DS(CMF)-1 structure table; and
- no independent electron-versus-muon composition coordinate is present.

Accordingly, it is **not authenticated as an arbitrary-composition npe-mu susceptibility EOS**.
Whether a future governed extension may use it is an owner/source decision.

### 6.3 Five-way data classification

| Quantity | A. physically present | B. current API exposes | C. StarProfile carries | D. exactly reconstructible | E. unavailable |
|---|---|---|---|---|---|
| `n_B` | cold `.eos` and raw axis | yes, `GetRho(P)` | yes | yes | no |
| equilibrium `Y_i` | cold `.eos` | yes, generic `GetRho_i(P)` | yes | yes | no |
| equilibrium `n_i` | not stored directly | no typed API | no | yes, exactly `Y_i n_B` on the carried background | no |
| equilibrium-path `mu_i` | reconstructible from cold raw `Q3-Q5` | no | no | yes from raw file plus charge conventions; not from production `.eos`/profile alone | no |
| arbitrary-composition `P`, `epsilon` | finite-T raw `Q1,Q7` on `(T,n_B,Y_q)`, different product | no | no | not for the canonical cold model or full npe-mu coordinate set | ambiguous/insufficient |
| arbitrary-composition individual `n_i`, `mu_i` | not established in current files | no | no | no | yes for required full state |
| effective masses | absent in current CMF directories | no | no | no | yes |
| `partial mu_i/partial n_j` | no | no | no | no | yes |
| `partial n_i/partial mu_j` | no | no | no | no | yes |

## 7. Built-star / StarProfile exposure

`StarProfile` defines core columns `r`, `m`, `nu'`, `P`, `epsilon`, `n_B`, `nu`, and `lambda`
(`CompactStar/Core/StarProfile.hpp:223-249`). `NStar::BuildFromTOV` gives them these exact current
units and labels:

| Field | Current stored unit / semantics | Evidence |
|---|---|---|
| radius `r` | km | `CompactStar/Core/src/NStar.cpp:101-104,178-180` |
| enclosed mass `m` | geometric km | `:105-106,182-184` |
| metric derivative `nu'` | km^-1 | `:108-109,186-187` |
| pressure `P` | geometric km^-2 | `:111-112,189-192` |
| energy density `epsilon` | geometric km^-2 | `:114-115,194-197` |
| baryon density `n_B` | fm^-3 | `:117-118,199-200` |
| `nu`, `lambda` | dimensionless metric exponents | `:120-124,202-214` |
| `d epsilon/dP` | dimensionless, separate EOS-owned column | `CompactStar/Core/StarProfile.hpp:313-331,893-913` |

The species registry and accessors are generic label-to-column mappings
(`CompactStar/Core/StarProfile.hpp:710-864`). They carry only the raw `Y_i` values copied from the
TOV extras; they do not distinguish number density, chemical potential, Fermi momentum, or
effective mass at the type level. No `mu_i`, `p_Fi`, or `m_i*` column exists. The comments/export
header still misleadingly say “species densities” (`CompactStar/Core/StarProfile.hpp:291-310`;
`CompactStar/Core/src/StarProfile.cpp:66-75`), but ADR-0001 is binding:

> `Y_i = n_i/n_B` is the stored primitive; consumers form `n_i=Y_i n_B`
> (`docs/adr/ADR-0001-species-profile-semantics.md:87-103,125-149`).

The production direct-Urca consumer follows that rule explicitly, multiplying the neutron, proton,
and electron fractions by `n_B` before computing Fermi momenta
(`CompactStar/Physics/Evolution/src/StarContext.cpp:493-498,550-597`). The charge-fraction cache
also correctly sums strong-sector `q_i Y_i` (`:697-750`). The only identified wrong consumer is the
uncompiled `RotochemicalCache`, discussed in section 12; this task does not repair it.

## 8. Equilibrium versus off-equilibrium determination

For canonical ordinary stars, the state actually available to production code is

`P(n_B), epsilon(n_B), Y_i(n_B)`

on a one-dimensional beta-equilibrium path, represented internally as splines against `P`. The raw
cold table adds reconstructible `mu_i(n_B)` along that same path, but current production code does
not expose it.

No compiled, tested, generic implementation was found that evaluates

`epsilon(n_B, Y_p, Y_e, Y_mu, ...), P(n_B, composition), mu_i(n_B, composition)`

with the required independent variables. The apparent alternatives do not change that result:

- `SigmaOmegaRho`, `_nstar`, and `Fermi_Gas_Many` solve equilibrium constraints internally; their
  solver coordinates/Jacobians are not a documented off-equilibrium thermodynamic API.
- `Particle` setters are mutable low-level primitives, not a consistent interacting EOS surface.
- `CompOSE_Thermo` is arbitrary in `Y_q` only for entropy/heat capacity.
- The finite-T DNS raw table has a `Y_q` direction, but not an authenticated independent
  electron/muon composition direction or current chemical-potential API, and it is not the
  canonical cold model.
- `SigmaOmegaRho_npemu` is historical candidate code outside the EOS CMake lists, contains a stale
  include, and still solves equilibrium; it is not a capability.

Therefore, derivatives **along** the beta-equilibrium curve cannot substitute for the corrected
susceptibility. Along the curve,

`d mu_i/d n_B = sum_j (partial mu_i/partial n_j) (d n_j/d n_B)`

is one constrained directional combination. It contains no information about independent
composition directions transverse to the curve and cannot recover the full local Jacobian or its
charge-reduced physical subspace. This is an identification/rank problem, not an interpolation
detail.

## 9. The 2005 treatment and the 2006 correction

### 9.1 What 2005 did

The 2005 benchmark linearized intrinsic chemical potentials and densities (eq. (10)), integrated
the response into a global paper `B_ij` (eqs. (11)-(12)), and inverted the full matrix (eq. (13)).
It then combined the inverse with reaction particle-number changes and rotationally changing
equilibrium populations to obtain the `Z` and `W` coefficients used in the non-superfluid evolution
equations (eqs. (51)-(58)). Its eq. (9) redshifted only the intrinsic chemical-potential
perturbation.

### 9.2 Why that full inverse is not rigorous

The 2006 authority includes the electrostatic potential in total redshifted chemical potentials
(eqs. (2), (10)) and imposes local charge neutrality (eq. (4)). Solving the neutrality constraint
for the electrostatic perturbation (eq. (11)) adds the charge-projection term missing from 2005
eq. (12), producing corrected paper `B_ij` in eq. (13).

That projection builds charge conservation into the integrated response. For the four-species
`n,p,e,mu` system the rows are linearly dependent (`B_pj=B_ej+B_muj` in the paper's unit-charge
convention). A perturbation along the charge vector merely changes the electrostatic zero point and
does not change particle distributions. The full corrected matrix is therefore singular, so
inverting it would attempt to invert a gauge/null direction.

### 9.3 What 2006 changes

The corrected treatment defines three physical combinations,
`tilde mu_n=delta mu_n` and `tilde mu_l=delta mu_p+delta mu_l` for `l=e,mu`, eliminates the proton
row and column, and inverts the resulting symmetric `3 x 3` paper `Btilde`. Proton density and
particle number follow from local charge neutrality. Global baryon conservation then eliminates
the conserved baryon-number perturbation, yielding the observable imbalance relations in eqs.
(14)-(15) and corrected `Z` coefficients in eqs. (16)-(18). Equation (19) writes the corresponding
`W` coefficients.

The physical imbalances remain invariant because the proton and lepton charges cancel the
electrostatic term (2006 eq. (9)). The 2005 reaction-rate, thermal-balance, structural-compression,
and ODE framework remains a useful reproduction benchmark. The 2006 paper states narrowly that
2005 eqs. (9), (12), (13), and (54)-(56) must be replaced; the coefficient values can change, while
the old-star quasi-steady solutions are unchanged and the reported transient differences are
small. This audit does not independently rederive or extend that conclusion.

## 10. Explicit susceptibility answers

**Q1. Can current CompactStar evaluate `partial mu_i / partial n_j` at fixed independent values of
the other appropriate variables?**

No. Low-level `Dmu` methods and the `SigmaOmegaRho::Jacob` do not document or implement that fixed-
physical-variable contract. Their equilibrium-solver coordinates are not the paper's independent
density coordinates.

**Q2. Can it evaluate `partial n_i / partial mu_j` required by the paper?**

No. There is no such production API, table column, tested implementation, or candidate with a
corrected charge-projected contract.

**Q3. Does an existing arbitrary-composition EOS contain enough information to calculate the
derivatives correctly?**

Not on authenticated current evidence. The finite-T DNS table has `T,n_B,Y_q` and collective
chemical potentials, which is potentially useful, but it is a different model/data product and has
no independent electron-versus-muon coordinate or composition/micro table. The current code
exposes only `Q2`. It cannot be certified as sufficient for the corrected npe-mu susceptibility
without owner/source adjudication and a governed data contract.

**Q4. Are beta-equilibrium curve derivatives sufficient?**

No. They provide only one directional derivative after beta equilibrium and charge constraints
have already tied all fractions to `n_B`. The corrected local response requires enough independent
directions to build the physical charge-reduced susceptibility; transverse information cannot be
reconstructed from a one-dimensional curve.

**Q5. How are electron and muon chemical potentials handled?**

For the canonical star path, neither is exposed. The raw cold CompOSE `Q3-Q5` allow their
equilibrium-path values to be reconstructed using the CompOSE charge convention, but the
standardized `.eos`, `TOVSolver`, and `StarProfile` omit them. Dormant RMF/free-particle code treats
leptons analytically as free Fermi gases (`Particle.cpp:123-153` and
`CompOSE_EOS.cpp:1242-1245,1293-1299,1326-1332`); that is not wired into canonical stars and is not
a governed Phase-5 choice.

**Q6. Is local charge neutrality imposed, and where?**

The cold EOS data already encode it: `Y_p=Y_e+Y_mu` on nucleonic rows. The TOV loader does not
impose or validate neutrality; it copies extras. `StarContext` reconstructs the strong charge
fraction for thermal use but does not create the equilibrium. Dormant RMF solvers impose dependent
charged densities inside their equilibrium solve. Perturbative neutrality and `delta psi` are not
implemented anywhere.

**Q7. Is there an independent charged-composition degree of freedom?**

Not in the canonical cold table: the single-point `Y_q` axis and beta-equilibrium fractions show
that composition has already been reduced to one trajectory. The separate finite-T DNS table has
an independent strong `Y_q`, but no independent electron/muon split and no current API beyond
entropy.

**Q8. Does code already compute paper `B_ij`, `Z`, `W`, susceptibility, or imbalances?**

No corrected paper `B_ij`, reduced `Btilde`, susceptibility, or governed chemical imbalance exists.
The uncompiled candidate uses names `B_i` and `Z_i` for structural particle-number derivatives,
not the paper thermodynamic matrix/coefficients, and its channel mapping is unfinished. No `W_npe`
or `W_npmu` implementation was found.

## 11. Capability matrix

`YES*` means usable only for the equilibrium background role stated, not for the independent-
composition susceptibility.

| REQUIREMENT | SOURCE AUTHORITY | CURRENT EOS FILE HAS IT? | CURRENT EOS API EXPOSES IT? | STARPROFILE CARRIES IT? | EXACTLY RECONSTRUCTIBLE? | EQUILIBRIUM-ONLY? | ARBITRARY-COMPOSITION? | DERIVATIVE AVAILABLE? | TESTED? | STATUS |
|---|---|---|---|---|---|---|---|---|---|---|
| `n_B` background | 2005 §2.2 | yes | yes | yes | yes | yes | no | `dn_B/dP` spline derivative not exposed; values tested | yes | READY |
| equilibrium `Y_i` | 2005 §3.1 | yes | yes, generically/misnamed | yes | yes | yes | no | along-path spline exists internally, no typed derivative | production consumers tested | EQUILIBRIUM-PATH-ONLY |
| equilibrium `n_i=Y_i n_B` | 2005 eqs. (10), (24); ADR-0001 | indirect | no typed API | indirect | yes, exactly on background | yes | no | no thermodynamic partial | DU reconstruction tested | READY |
| equilibrium-path `mu_i` | 2005 eqs. (7)-(10); 2006 eq. (2) | raw `Q3-Q5`, not standardized `.eos` | no | no | yes from raw collective potentials and charges | yes | no | no | no | AVAILABLE-BUT-NOT-EXPOSED |
| arbitrary-composition `P,epsilon` matching canonical cold EOS | 2005 §2.2/§3.1 | no; separate DNS product is not matched authority | no | no | no | n/a | no certified capability | no | no | MISSING |
| separate DNS `(T,n_B,Y_q)` thermodynamic surface | CompOSE format; no rotochemical authority decision | raw `Q1-Q7` | only `Q2` | no | not as the required npe-mu state | no | `Y_q` only | only `dQ2/dT` | entropy path tested | AMBIGUOUS |
| independent individual `n_i,mu_i` state | 2005 eq. (10) | no | no | no | no | n/a | no | no | no | MISSING |
| local `partial mu_i/partial n_j` | equivalent susceptibility form subject to 2006 constraint | no | no | no | no | n/a | no | no | no | MISSING |
| local `partial n_i/partial mu_j` | 2005 eq. (10) | no | no | no | no | n/a | no | no | no | MISSING |
| electrostatic/charge-projected local response | 2006 eqs. (10)-(13) | no | no | no | no | evaluated at equilibrium but needs independent response | no | no | no | MISSING |
| local neutrality of equilibrium background | 2006 eq. (4) | yes, encoded | not enforced; values exposed | fractions carry evidence | yes/checkable | yes | DNS `Y_q` only | n/a | indirectly in composition consumers | READY |
| perturbative neutrality and `delta psi` | 2006 eqs. (10)-(11) | no | no | no | no | linearized at equilibrium | no | no | no | MISSING |
| corrected global paper `B_ij` | 2006 eqs. (12)-(13) | no | no | no | no | based on equilibrium susceptibility | no | prerequisite absent | no | MISSING |
| singular/null-mode identification | 2006 text after eq. (13) | literature fact | no representation | no | conceptually known, not numerically constructible | n/a | n/a | n/a | no | DERIVABLE-WITH-GOVERNED-EXTENSION |
| reduced paper `Btilde` / inverse | 2006 between eqs. (13)-(14) | no | no | no | no | based on equilibrium susceptibility | no | prerequisite absent | no | MISSING |
| paper `Z_np`, `Z_npe`, `Z_npmu` | 2006 eqs. (15)-(18) | no | no | no | no | global background coefficient | no | prerequisite absent | no | MISSING |
| `I_{Omega,i}` / fixed-baryon particle-number drive | 2005 eq. (30) | equilibrium composition exists | Phase-4 structural fields exist; particle response API absent | structural fields yes, result no | not yet: measure-complete species response/reduction absent | background yes | no | structural derivative interface only | Phase-4 structure tested, coefficient not | DERIVABLE-WITH-GOVERNED-EXTENSION |
| paper `W_npe`, `W_npmu` | 2006 eq. (19) | no | no | no | no | global coefficient | no | prerequisites absent | no | MISSING |
| effective masses for benchmark rates/heat | 2005 §3.1, eq. (50) | no current auxiliary files | dormant API expects them | no | no | n/a | n/a | n/a | no | MISSING |
| `d epsilon/dP` for Phase-4 structure | Phase-4 authority, not paper susceptibility | yes via barotrope | yes | yes | yes | yes | no | yes | yes | READY |

No row is marked `READY` merely because a similarly named value exists. In particular, the ready
background `n_i` does not imply a ready susceptibility.

## 12. Existing rotochemical candidate code and history

Searches covered `rotochemical`, chemical imbalance variants, `eta_npe`, `eta_npmu`, beta
equilibrium, chemical potential, spin-down heating, `Z_npe`, `Z_npmu`, `W_npe`, `W_npmu`, `B_ij`,
susceptibility, and naming variants in the current tree and Git history.

| Candidate | Lineage / status | Assumptions and defects | Disposition |
|---|---|---|---|
| `Physics/Evolution/RotochemicalCache.{hpp,cpp}` | introduced by `675b4a9` (“Add Hartle O(Omega^2) second-order solver and rotochemical heating framework”); source absent from `Physics/Evolution/CMakeLists.txt:22-32`; `Build()` has no caller/test | computes enclosed counts and structural `A_i/B_i/Z_i`; passes stored `Y_i` as though it were `n_i` (`src/RotochemicalCache.cpp:25-44,47-104,147-171`); predates ADR-0001 and Phase-4 correction; omits required measure/normalization terms recorded by INV-09 | historical candidate evidence only; not authority or capability |
| `Physics/Driver/Chem/Rotochemical.{hpp,cpp}` | same `675b4a9` lineage; header installed but Chem source list is empty (`Physics/Driver/Chem/CMakeLists.txt:1-11`) | comments cite 2005; writes per-species structural `Z_i` while acknowledging unresolved reaction-channel mapping (`src/Rotochemical.cpp:98-123`); no corrected paper `Btilde`/`Z` | historical scaffold; uncompiled/untested |
| `Physics/State/ChemState.hpp` | compiled generic state container | says units/meaning belong elsewhere; no governed redshift frame or channel ordering; INV-11 records no convention (`docs/SCIENTIFIC_INVARIANTS.md:820-836`) | generic storage only, not scientific authority |
| `Microphysics/Rates/Urca.hpp` | header-only abstract rate interface under an old `ChemicalHeating` namespace; no implementation/CMake source found | accepts a vector named `eta`, but defines no npe-mu convention, susceptibility, or coefficient | dormant interface evidence only |
| `Physics/SigmaOmegaRho_npemu.{hpp,cpp}` | earlier history includes `d07ac2a` and repository reorganization `dfb4443`; current source not in an EOS/Physics CMake source list | source labels itself incomplete/stale, contains an obsolete include, and solves beta equilibrium; no arbitrary-composition contract | historical/project-specific candidate |

The repository's structural `B_i`/`Z_i` names in these candidates are **not** the 2006 paper's
thermodynamic paper `B_ij`/corrected `Z` objects. Existing documentation already classifies the
candidate as uncompiled and nonconformant (`docs/architecture/CURRENT_ARCHITECTURE.md:192,334-344`;
`docs/SCIENTIFIC_INVARIANTS.md:712-799`). Nothing was revived or modified.

## 13. Four physics layers and current interfaces

| Layer | Scope | Current state | Boundary exposed to next layer |
|---|---|---|---|
| **A — local EOS thermodynamics** | `n_i`, `mu_i`, susceptibilities/Jacobians, electrostatic/charge constraint | equilibrium `n_B,Y_i` only; raw equilibrium chemical potentials not exposed; required off-path response missing | Not ready for corrected paper `B_ij` construction |
| **B — nonrotating stellar integration** | redshifted/global paper `B_ij`, reduced `Btilde`, paper `Z`; scalar particle-number integrals | metric/proper-volume infrastructure exists; no susceptibility or corrected integrals | Blocked by Layer A and unresolved measure-complete particle response |
| **C — rotational structural driving** | verified Phase-4 Hartle response and fixed-baryon sequence reduction | fixed-`epsilon_c` Phase-4 response verified/ratified; particle-number and baryon-conserving reduction remain INV-09 work (`docs/validation/PHASE4_CLOSEOUT.md:201-215`) | Structural fields ready; `I_Omega/W` not ready |
| **D — secular evolution** | `T_infinity(t)`, `eta_npe(t)`, `eta_npmu(t)`, reactions/heating | passive thermal path exists; chemical conventions/rates/ODEs absent | Out of scope; not implemented here |

This audit is primarily Layer A, with factual interface checks into B and C. It makes no Layer-D
scientific decision.

## 14. Unit and variable audit

| Quantity | Current/source unit | Conversion / warning |
|---|---|---|
| `n_B` | `fm^-3` in CompOSE, standardized `.eos`, TOV, and StarProfile | a physical count integral converts `fm^-3` to `km^-3` by `10^54`; this factor is not centralized in `Units.hpp` (`CompactStar/Units.hpp:23-34`) |
| `Y_i=n_i/n_B` | dimensionless | never treat as density; `n_i=Y_i n_B` in `fm^-3` |
| `n_i` | required locally in `fm^-3` | derived exactly only for the equilibrium background; no independent off-path `n_i` API |
| intrinsic `mu_i` | papers and CompOSE reconstructed values: MeV; legacy particle internals: `fm^-1` | legacy conversion uses the external `MEV_2_INV_FM`; no typed unit contract |
| `P` | standardized table/TOV input: `dyne cm^-2`; StarProfile: geometric `km^-2`; raw CompOSE `Q1=P/n_B` in MeV | conversion occurs in `NStar.cpp:189-192`; do not combine raw and standardized quantities silently |
| `epsilon` | standardized table/TOV input: mass-density-style `g cm^-3`; StarProfile: geometric `km^-2`; raw CompOSE `Q7=epsilon/(n_B m_n)-1` | conversion occurs in `NStar.cpp:194-197`; `MeV fm^-3 -> erg cm^-3` is `1.602176634e33` (`CompactStar/Units.hpp:58-60`) |
| `d epsilon/dP` | dimensionless after the production cgs-to-geometric factor | this barotropic derivative is not a composition susceptibility (`TOVSolver.hpp:800-826`) |
| `partial mu_i/partial n_j` | `MeV fm^3` in the requested local convention | legacy RMF headers label natural-unit values `fm^2`; multiplication by `hbar c` is required to express `MeV fm^3` |
| `partial n_i/partial mu_j` | `fm^-3 MeV^-1` | reciprocal-matrix units only after an independent, nonsingular coordinate choice |
| corrected global paper `B_ij` | `MeV^-1` when `dV` converts to `fm^3` and `mu` is MeV | redshift/charge factors are dimensionless; full four-species matrix is singular |
| reduced inverse / paper `Z` | MeV in the same convention | must come from `Btilde^{-1}`, not legacy structural `Z_i` counts |

`Units.hpp` intentionally is not a general units library and explicitly leaves geometric EOS and
`fm^-3 <-> km^-3` conversions elsewhere (`CompactStar/Units.hpp:17-39`). Any later extension must
make these boundaries explicit before numerical work.

## 15. Minimum-gap classification

**D. A new thermodynamic EOS extension/data product is scientifically required.**

Why the other choices are rejected:

- **A** is false: the corrected susceptibility and reduced system are absent.
- **B** is only partly true for equilibrium-path chemical potentials. Exposing raw `Q3-Q5` would
  not create missing transverse composition information.
- **C** is not established: the dormant RMF solvers are equilibrium generators, while the finite-T
  DNS table is a distinct, incomplete-for-npe-mu candidate data product rather than a certified
  existing off-equilibrium EOS implementation.
- **E** is not required to answer the capability question: current capability is decisively
  insufficient. Owner/source adjudication is still required to choose the future scientific
  extension.

This is a gap classification, not a design. No finite differencing, automatic differentiation,
analytic-derivative choice, matrix inversion method, independent-coordinate basis, interpolation
scheme, or EOS replacement is selected here.

## 16. Explicit unresolved questions

1. Which arbitrary-composition EOS model/data product is authoritative for the same ordinary-star
   backgrounds used by Phase 5: a DS(CMF)-1-compatible extension, a governed RMF implementation,
   another published table, or something else?
2. What independent thermodynamic coordinates and held-fixed convention define the local
   susceptibility before the 2006 charge projection?
3. How are electron and muon lepton-number degrees of freedom represented independently, and are
   they free-Fermi analytic contributions or EOS-table-owned quantities?
4. Which species set is in scope when hyperons/deltas/quarks appear in the DS(CMF)-1 composition,
   versus the four-species npe-mu formalism of the benchmark?
5. What are the authoritative rest-mass and collective-chemical-potential conventions for mapping
   CompOSE `Q3-Q5` to individual `mu_i`?
6. How will phase transitions/composition jumps be represented in local susceptibility and in the
   measure-complete global integrals?
7. What exact core boundary/species-domain rule governs each global paper `B_ij` integral?
8. How will the corrected thermodynamic paper `Z/W` names coexist with the repository's historical
   structural `B_i/Z_i` names without semantic collision?
9. INV-11 remains unresolved: chemical-imbalance ordering, redshift frame, and units require an ADR
   before any evolution implementation.

## 17. Recommended next governed task

Create a source/owner-adjudicated **Phase 5A-1 thermodynamic EOS requirements and data-contract
decision record**. It should select the authoritative arbitrary-composition EOS/data source and
independent variables, specify the 2006 charge-reduced susceptibility contract and units, and define
validation fixtures—without yet implementing derivatives, matrix inversion, or evolution.

## 18. Scope closure

**NO ROTOCHEMICAL PHYSICS IMPLEMENTED IN PHASE 5A-0**

No production source, test, baseline, EOS data, or literature file was changed. Phase-4 status is
unchanged. BNV and Layer-D secular evolution were not begun.
