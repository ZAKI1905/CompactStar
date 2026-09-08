# Phase-5D-0 — Secular rotochemical evolution scientific preflight

**Class:** documentation-only scientific preflight. No production source, test, baseline, EOS/data,
CMake, literature byte, or accepted-ADR semantic change is made by this task.
**Companion decision:** `docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md` (PROPOSED).
**Disposition:** post-Phase-5D-0RR material-closure revision; see §33.4–§35.

> **PHASE-5D MATERIAL CLOSURE COMPLETE — R1/R2/R3 CLOSED —
> READY FOR FINAL BOUNDED INDEPENDENT RE-REVIEW.**

---

## 1. Authenticated entry

| Item | Value | Check |
|---|---|---|
| Canonical master (local) | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` | matches required |
| Canonical master (`origin/master`) | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` | equal to local |
| Phase-5C production implementation | `4d78bf4000848ddecc2127daa2f2840872f266f5` | ancestor of ratification tip |
| Phase-5C immutable predeclaration | `a87f0212c2bd7bfba92db91dfac82447a6561334` | ancestor of ratification tip |
| Phase-5C ratification tip (local) | `27727016856a6a25a46e447c70e380722ea8ddbf` | matches required |
| Phase-5C ratification tip (`origin/physics/phase5c-corrected-chemical-coefficients`) | `27727016856a6a25a46e447c70e380722ea8ddbf` | equal to local |
| `master` is ancestor of ratification tip | yes | `git merge-base --is-ancestor` |

Recorded ancestry of the branch point (newest first):

```text
2772701 docs: ratify corrected chemical coefficients          <- branch point / PHASE5C2_RATIFICATION_SHA
4d78bf4 feat: implement corrected chemical coefficients        <- PHASE5C2_SHA
a87f021 docs: predeclare chemical coefficient acceptance       <- immutable predeclaration
d3b102d docs: ratify structural uncertainty semantics
09d1b3c test: plan corrected chemical coefficients
4780121 docs: ratify corrected rotochemical coefficients
54ec7ab docs: preflight corrected rotochemical coefficients
49ab2b8 docs: close particle-number structural response        <- canonical master
```

Branch `analysis/phase5d-secular-rotochemical-evolution-preflight` did not exist locally or on
`origin` before this task; the worktree path
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-rotochemical-preflight` did not
exist. Both were created fresh from `27727016856a6a25a46e447c70e380722ea8ddbf`; the worktree was
clean at creation. Phase-5C is **not** canonically integrated, and this preflight does not
integrate it.

---

## 2. Literature authority

Library root `/Users/keeper/Documents/CompactStar/literature` was read-only throughout. No byte was
added, removed, or modified. Hashes recomputed in this task:

| Short name | File under `literature/rotochemical/` | SHA-256 | Declared in task? | Role used here |
|---|---|---|---|---|
| **FR2005** | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88` | yes, matches | **Primary** non-superfluid formalism; eqs. (1)–(87), App. A |
| **R2006** | `2006-Reisenegger-...-Electrostatic-Potential-Perturbations.pdf` | `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8` | yes, matches | **Primary** corrected electrostatic/chemical coefficients; eqs. (10)–(19) |
| **R1995** | `1995-Reisenegger-Rotochemical-Heating.pdf` | `9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff` | yes, matches | **Primary** origin of the imbalance functions; Appendix eqs. (22)–(34), pp. 14–15 |
| **Y2020** | `2020-Yanagi-NS-Therm-Thesis.pdf` | `69590539c275fa679a5521a9c5abedd9fdc58718b1554827786c8f707bc618cc` | yes, matches | **Supporting** cross-check only; §2.4.2, §4.1, eqs. (2.75)–(2.76), (4.3)–(4.20) |
| **R1997** | `1997-Reisenegger-Constraining-Dense-Matter-Superfluidity-...pdf` | `19a10133511aefc05ece33d7454cba60def03de85d9980c55fd7229058c2d08b` | not declared | Supporting physical context; `x>1/9` direct-Urca statement |
| **JRF2006** | `2006-Jofre-Reisenegger-Fernandez-Gravitochemical-Heating-Model.pdf` | `2dd5444d19cebae12509fe4ecb7dac31957d332e2131813894a10cceac403109` | not declared | **Supporting** generic-external-source precedent; eqs. (4)–(8) |

All four hashes declared in the task statement match exactly. The two undeclared companions are
catalogued in `literature/catalog.tsv` and were consulted only in their catalogued supporting role.

**Supersession discipline applied.** Where R2006 and FR2005 disagree, R2006 governs (coefficient
definitions `Z_np`, `Z_npe`, `Z_npmu`, `W_npl`, and the `B_ij` integrand). Everywhere else FR2005
governs. Y2020 was used **only** to cross-check quantities the primaries state independently, and
in exactly one place (§9.3, the FR2005 eq. (37) printed typo/internal inconsistency) as a third witness to a
conclusion already established by independent derivation. Y2020 is never allowed to originate a
convention. Its §2.4.2 branch normalization is flagged in §10.4 as a *supporting-source dependency*,
not adopted.

**Sign-convention hazard documented in the sources themselves.** Y2020 footnote 3 (p. 61) records
that its `I^N_{M,Gamma}` has *the opposite sign* to the corresponding phase-space integrals in two of
its own references. R1995 footnote 3 (p. 14) records that its `eta` convention is opposite to that of
Sawyer (1989) and Haensel (1992). The literature therefore genuinely contains both signs. §7 fixes
one convention and authenticates it against the source *source-equations*, not against isolated
formulas.

---

## 3. Source equation crosswalk

Page numbers are FR2005 preprint page numbers as printed (`– n –`).

| Quantity | FR2005 | R2006 | R1995 | Y2020 | CompactStar name |
|---|---|---|---|---|---|
| Redshifted temperature `T_inf = T(r) e^{Phi}` | eq. (1), p. 3 | §3, p. 569 | — | §4.1 | `ThermalState::Tinf()` (INV-10) |
| Thermal balance | eq. (2), p. 3 | eq. (8), p. 569 | — | eq. (4.3) | ADR-0002 governing invariant |
| Heat capacity `C = int dV c_V` | eq. (3), p. 3 | — | — | — | `StarContext::HeatCapacityStar_Tinf` |
| `L_H^inf = int dV Q_H e^{2Phi}` | eq. (4), p. 3 | — | — | eq. (4.5) | *not implemented* |
| `L_nu^inf = int dV Q_nu e^{2Phi}` | eq. (5), p. 3 | — | — | — | `NeutrinoCooling` (structure only) |
| `L_gamma^inf = 4 pi sigma R_inf^2 (T_s^inf)^4` | eq. (6), p. 3 | — | — | — | `PhotonCooling` |
| `eta_npe`, `eta_npmu` definitions | eqs. (7)–(8), p. 3 | §3, p. 569 | — | §4.1 fn. 1 | ADR-0010 local `g_x` |
| `delta mu_i^inf = delta mu_i e^{Phi}` | eq. (9), p. 4 | eq. (2)+(10), p. 569–570 | — | — | ADR-0013 `eta^inf = e^nu eta_local` |
| `delta N_i = sum_j B_ij delta mu_j^inf` | eq. (11), p. 4 | eq. (12), p. 570 | — | eq. (4.6) | `GlobalChemicalNumberResponse` (`G_y`) |
| `B_ij = int dV (dn_i/dmu_j) e^{-Phi}` | eq. (12), p. 4 | **eq. (13)**, p. 570 (corrected) | — | — | `G_y` integrand (ADR-0013 §3.2) |
| `Ndot_i = int dV e^{Phi} sum_a DeltaGamma_ia` | eq. (15), p. 4 | eq. (6), p. 569 | — | eq. (4.8) | *not implemented* |
| `Ndot_i^eq = 2 Omega Omegadot I_{Omega,i}` | eq. (30), p. 7 | eq. (5), p. 569 | — | eq. (4.9) | `RotochemicalSpinDrive::IPhysical()` |
| `Q_a = Q_a^eq F_*(eta/kT)` | eq. (32), p. 10 | — | eqs. (22),(29),(31) | eq. (4.15),(4.19) | *not implemented* |
| `DeltaGamma_a = (1/kT) Q_a^eq H_*(eta/kT)` | eq. (33), p. 10 | — | eqs. (23),(30),(32) | eq. (4.15),(4.20) | *not implemented* |
| `F_D, H_D, F_M, H_M` | eqs. (34)–(37), p. 10 | — | eqs. (29)–(32), p. 15 | eqs. (4.19)–(4.20) | *not implemented* |
| `Q_H = sum_a DeltaGamma_a eta_a` | eq. (38), p. 10 | — | — | eq. (4.4),(4.5) | *not implemented* |
| Sign convention `DeltaGamma = Gamma_{A->B} - Gamma_{B->A}`, `eta = delta mu(A) - delta mu(B)` | eqs. (39)–(40), p. 10 | — | — | eq. (4.8) | §7 below |
| `Q_a^eq = S_a(n) T^q` | eq. (41), p. 10 | — | eqs. (33)–(34), p. 15 | eqs. (2.75)–(2.76) | `NeutrinoCoolingCachePayload` (placeholder `S`) |
| `xi = eta/(kT) = eta^inf/(k T_inf)` | eq. (42), p. 11 | — | — | eq. (4.18) | *not implemented* |
| `Ltilde_a = integral_{D_a} 4 pi r^2 e^Lambda S_a e^{(2-q)Phi} dr` | eq. (43), p. 11 | — | — | — | `K_DU_erg_s_K6`, `K_MU_erg_s_K8` (same structure; canonical unit `erg s^-1 K^-q`) |
| `L_a^inf = Ltilde_a F_*(xi) T_inf^q` | eq. (44), p. 11 | — | — | — | *not implemented* |
| `integral_{D_a} DeltaGamma_a e^Phi dV = [Ltilde_a/k_B^(erg)] H_*(xi) T_inf^{q-1}` | eq. (45), p. 11 | — | — | — | *not implemented* |
| `L_{H,a}^inf = Ltilde_a xi H_*(xi) T_inf^q` | eq. (46), p. 11 | — | — | — | *not implemented* |
| `M_*(xi) = xi H_*(xi) - F_*(xi)` | eq. (47), p. 11 | — | p. 8 cooling eqs. | — | *not implemented* |
| `L_H - L_nu = Ltilde M_*(xi) T_inf^q` | eq. (48), p. 11 | — | — | — | *not implemented* |
| Envelope `T_s^4(T_inf)` | eq. (49), p. 11 | — | — | — | `EnvelopePotekhin1997/2003` |
| `Ctilde = C/T_inf` | eq. (50), p. 12 | — | — | — | ADR-0002 (equivalent) |
| Thermal ODE | eq. (51), p. 12 | eq. (8) | — | eq. (4.3) | §16 |
| `etadot_npe`, `etadot_npmu` | eqs. (52)–(53), p. 13 | §5 (unchanged form) | — | eqs. (4.10)–(4.11) | §6 |
| `Z_np, Z_npe, Z_npmu` | eqs. (54)–(56), p. 13 (**superseded**) | **eqs. (16)–(18)**, p. 570 | — | eq. (4.7) | `ChemicalImbalanceResponse` |
| `W_npe, W_npmu` | eqs. (57)–(58), p. 13 | **eq. (19)**, p. 570 | — | eq. (4.12) | `RotochemicalSpinDrive` |
| `eta^inf = -Z delta N_l` | (implied by (52)–(58)) | **eq. (15)**, p. 570 | — | **eq. (4.7)** | ADR-0013 §3.3 |
| Large-`xi` asymptotics `C_H, C_M` | eqs. (59)–(60), p. 14 | — | — | — | §9.5 |
| Quasi-steady `eta_qs`, `L_gamma,qs` | eqs. (61)–(66), p. 14–15 | §3, p. 570 | §5 | §4.1 | §22 |
| `tau_eq`, `A` | eqs. (79)–(83), p. 17 | — | — | — | §23, §24 |
| Generic external source (`Gdot` analogue) | — | — | — | — | JRF2006 eqs. (5)–(8) → §28 |

The source crosswalk fixes the physical integrand, but source traceability is not an executable
detector for an implementation lapse-power error. Future gate RE10b therefore exercises
`GlobalUrcaChannelCoefficient` directly against an independent analytic/high-precision quadrature
of the displayed `Ltilde_a` integral; it does not reuse a precomputed or injected `Ltilde`. The
global-integrand mutations M10/M11, M12, M30, and M33–M37 point to this executable future oracle.

---

## 4. Current-code architecture audit

Audited at `27727016856a6a25a46e447c70e380722ea8ddbf`.

### 4.1 What is live and reusable unchanged

| Component | File | Finding |
|---|---|---|
| Evolved thermal DOF | `CompactStar/Physics/State/ThermalState.hpp:150-200` | `values_[0] = ln(T_inf / T_ref)`, `T_ref = 1e8 K`; `Tinf()`, `SetTinf()` accessors. Exactly the FR2005 `T_inf` in Kelvin. **Reusable unchanged.** |
| Heat-capacity ownership | ADR-0002 §Decision; `StarContext::HeatCapacityStar_Tinf` (`StarContext.hpp:152`) | `C_*(T_inf) = int c_V(T_inf e^{-nu}) 4 pi r^2 e^Lambda dr`. Identical to FR2005 eqs. (3)+(50). ADR-0002 *already* writes its governing invariant as `C_* dT_inf/dt = -L_nu - L_gamma + L_H + ...`. **No ADR-0002 change is needed to add rotochemical terms.** |
| Photon cooling | `PhotonCooling.hpp:287-297` | `DependsOn = {Thermal}`, `Updates = {Thermal}`. No chemical coupling. **Unchanged** (§17). |
| Envelope models | `Driver/Thermal/Boundary/EnvelopePotekhin1997.hpp`, `EnvelopePotekhin2003.hpp`, `IEnvelope.hpp`, `SurfaceGravity.hpp`, `TbDefinition.hpp` | Potekhin 1997 is FR2005's own envelope family (FR2005 eq. 49). **Reusable unchanged.** |
| Evolution core | `EvolutionSystem.hpp/.cpp`, `RHSAccumulator.hpp`, `StateLayout.hpp`, `StatePacking.hpp`, `StateVector.hpp` | Tag-blocked flat state vector; `StateLayout::Configure(state, {tags})` selects active blocks; drivers accumulate via `RHSAccumulator::AddTo(tag, i, value)`. **The engine already supports a clean vector state.** |
| Chemical state block | `CompactStar/Physics/State/ChemState.hpp` | Exists, `Resize(N)`, `Eta(i)`, `PackTo/UnpackFrom`. **Storage is reusable; its semantics are undefined (INV-11).** |
| Spin state / driver | `SpinState.hpp:175` (`Omega()` rad/s), `Driver/Spin/src/MagneticDipole.cpp` | `dOmega/dt = -K |Omega|^n sign(Omega)`. Spin is already an independent driver. **Reusable unchanged.** |
| Boltzmann constant | `NeutrinoCooling_Details.cpp:213-220` | `MEV_PER_K = Zaki::Physics::K_BOLTZ_EV * 1e-6` = `8.61733326214518e-11 MeV/K`, single authority (`docs/validation/PHASE3C_BOLTZMANN_AUTHORITY.md`). **No new constant is needed for `xi`.** |
| Chemical coefficients | `CompactStar/Analysis/ChemicalResponse.hpp` | `GlobalChemicalNumberResponse` (`G_y`), `ChemicalImbalanceResponse` (`Z`), `RotochemicalSpinDrive` (`W`, `IPhysical()`, `Evaluate(Omega, Omegadot)` returning `2 W Omega Omegadot`). **Reusable unchanged; `Evaluate` is exactly the ratified spin drive.** |

### 4.2 What must be extended

| Component | Why |
|---|---|
| `NeutrinoCoolingCachePayload` | Holds two **lumped** coefficients `K_DU_erg_s_K6`, `K_MU_erg_s_K8`. The chemical ODE needs the **electron/muon-resolved** set `{Ltilde_De, Ltilde_Dmu, Ltilde_Me, Ltilde_Mmu}` (FR2005 eqs. 52–53). Must gain channel resolution; the equilibrium driver then consumes the sum (§13). |
| `NeutrinoCooling_Details::ComputeDerived` | Currently `L_nu = K_DU T^6 + K_MU T^8`. The controlled benchmark must construct equilibrium cooling and its nonequilibrium extension from the same declared channel-resolved `Ltilde`; exact reduction at `xi=0` is a same-coefficient identity, not a promise to match the historical placeholder coefficients (§13). |
| `StarContext` direct-Urca support | `BuildDirectUrcaMaskCache_` (`StarContext.cpp:475-600`) tests only `kF_n <= kF_p + kF_e`. A muon criterion does not exist. Its last-index representation can also sweep a closed inner region if a future allowed outer shell appears; support/applicability must become explicit (§15). |
| `EvolutionConfig::n_eta` | `= 0` in `RunBuilder.cpp:39` and every `main/Test` program. Must become 2 with a named channel order. |
| Driver ordering | `EvolutionSystem::operator()` (`EvolutionSystem.cpp:119-128`) runs drivers in **registration order** with no dependency sort, although `IDriver::DependsOn()` exists. A rotochemical driver must not depend on that order to obtain `Omegadot` (§18). |

### 4.3 Dead scaffolding — present in the tree, **not** in the build

None of the following is compiled. All predate the governance regime and **must not be treated as
prior art**:

| Path | Build status | Defect |
|---|---|---|
| `CompactStar/Microphysics/Rates/Urca.hpp` | `Microphysics/CMakeLists.txt:5` — `# add_subdirectory(Rates)` | File is named `Urca.hpp` but its content is a `RateChannels.hpp` interface in a nonexistent `CompactStar::ChemicalHeating` namespace. No implementation. Not included anywhere. |
| `CompactStar/Physics/Driver/Chem/src/Rotochemical.cpp` | `Driver/Chem/CMakeLists.txt` — sources list is **empty** | Implements `d(eta_i)/dt += Z_i * 2 Omega Omegadot` calling `Z_i = dN_i/dOmega^2` — this conflates `W` with `Z`, uses per-species rather than per-channel indices, and has no relaxation term. **Scientifically wrong under the ratified convention.** Recorded as such in INV-11. |
| `CompactStar/Physics/Driver/Chem/WeakRestoration.hpp` | header list only | **0 bytes.** |
| `CompactStar/Physics/Driver/Thermal/HeatingFromChem.hpp` | `Driver/Thermal/CMakeLists.txt` — `# HeatingFromChem.hpp` commented out | No `.cpp`. |
| `CompactStar/Physics/Evolution/RotochemicalCache.hpp` | header only | Recorded as INV-01 nonconformant in `docs/SCIENTIFIC_INVARIANTS.md`. |

### 4.4 Where should `eta` live?

**Recommendation: a rotochemical module supplying additional state semantics and RHS terms; not the
generic evolution core.** Evidence: the core (`StateLayout`, `RHSAccumulator`, `StatePacking`) is
already fully generic over tagged blocks and needs no rotochemical knowledge; `ChemState` already
provides untyped storage; and every rotochemical quantity (`Z`, `W`, `Ltilde`, `F_*`, `H_*`,
channel order, `xi`) is physics that must carry Phase-5C provenance. Putting channel semantics into
the core would force the core to depend on `CompactStar/Analysis`. The channel *names and order*
belong to the rotochemical module and are simply honoured by `ChemState` indices.

---

## 5. Canonical evolved chemical state

### 5.1 Decision

Store the **redshifted chemical imbalances** directly:

```text
y_chem = ( eta_npe^inf , eta_npmu^inf )        [MeV], channel order (Npe, NpMu)
```

**Not** `delta N_e`, `delta N_mu`.

### 5.2 Authentication

FR2005 §2.1, p. 4: *"Being uniform in the star, these redshifted chemical imbalances are the ideal
variables to quantify the departure from chemical equilibrium and follow its time evolution"*
(R2006 §3, p. 570, states this verbatim of `eta_npl^inf`). FR2005 integrates eqs. (52)–(53) in
`eta^inf`; R2006 §5, p. 571, states that the corrected formalism changes *only* the coefficient
values *"[i]n their dynamical system for the time evolution of the variables `T^inf` and
`eta_npl^inf`"*; Y2020 eqs. (4.10)–(4.11) evolve `eta_e^inf`, `eta_mu^inf`. All three primaries and
the supporting source evolve `eta^inf`.

Decisive additional reasons against `delta N`:

1. `xi_a = eta_a^inf/(k_B^(MeV) T_inf)` is needed **every RHS evaluation**. From `delta N` it costs a
   matrix multiply by `Z` each time, and `Z` is exactly the ill-conditioned object
   (`G_condition = 2995.7` for the governed free-gas fixture).
2. `eta^inf` is **spatially uniform** (diffusive equilibrium) and hence a genuine scalar state;
   `delta N` is a derived global integral.
3. `eta^inf` is the quantity with a source benchmark (FR2005 eqs. 64–65; R2006 Fig. 2).
4. The `delta N` route only becomes preferable when baryon number is *not* conserved — which is
   the future-BNV case and is handled by the unreduced `G_y` seam (§28), not by changing v1.

### 5.3 Full specification

| Property | Value | Authority |
|---|---|---|
| Symbol | `eta_npe^inf`, `eta_npmu^inf` | FR2005 eqs. (7)–(9),(16)–(17) |
| Units | **MeV** | ADR-0013 §3.3 (`Z` in `MeV/count`) |
| Local definitions | `eta_npe = mu_n - mu_p - mu_e`, `eta_npmu = mu_n - mu_p - mu_mu`, each the *deviation* `delta mu_n - delta mu_p - delta mu_l` | FR2005 eqs. (7)–(8); ADR-0010 |
| Redshift relation | `eta^inf = e^nu eta_local`, `nu = Phi` | FR2005 eq. (9); ADR-0013 §3.4 (already governed) |
| Sign convention | `eta^inf = -Z delta N_l`, `Z` symmetric positive definite in `MeV/count` | **R2006 eq. (15)**; Y2020 eq. (4.7); ADR-0013 §3.3 |
| Channel order | index 0 = `Npe`, index 1 = `NpMu` | ADR-0013 §3.3 (`BetaChannel::Npe`, `BetaChannel::NpMu`) |
| Initial condition | `eta_npe^inf(0) = eta_npmu^inf(0) = 0` (caller-supplied, default) | FR2005 Fig. 4 caption; R2006 Fig. 2 caption |
| Local reconstruction | `eta_local(r) = eta^inf e^{-nu(r)}` | FR2005 eq. (9) |
| Validity domain | `\|eta\| << mu_i^eq` (linear-response `B_ij`); connected diffusive chemical domain declared by `G_y` | FR2005 §2.1, p. 4; ADR-0013 §3.2 |
| Ownership | a rotochemical evolution module; `ChemState` provides storage only | §4.4 |

The ratified `eta^inf = -Z delta N_lepton` convention is **confirmed, not merely assumed** — see §6.2
for the independent re-derivation and §6.3 for the three-witness check.

---

## 6. Full chemical ODE derivation

### 6.1 Statement

```text
d/dt ( eta_npe^inf  )  = - [ Z_npe  Z_np   ] ( R_e  )  +  2 ( W_npe  ) Omega Omegadot
       ( eta_npmu^inf )      [ Z_np   Z_npmu ] ( R_mu )         ( W_npmu )
```

with the **global net lepton-production rates**

```text
R_l  ==  integral_D dV e^{Phi} sum_{a in l-channels} DeltaGamma_a        [count / s]
      =  sum_{a in l} [Ltilde_a/k_B^(erg)] H_a(xi_l) T_inf^(q_a-1)
```

Written out, exactly as asked in §6 of the task statement:

```text
etadot_npe^inf   = - Z_npe  * DeltaGamma_npe  - Z_np    * DeltaGamma_npmu + 2 W_npe  Omega Omegadot
etadot_npmu^inf  = - Z_np   * DeltaGamma_npe  - Z_npmu  * DeltaGamma_npmu + 2 W_npmu Omega Omegadot
```

where `DeltaGamma_npe == R_e` and `DeltaGamma_npmu == R_mu` are the **global, `e^Phi`-weighted,
channel-summed** net reaction rates in `count/s` — *not* local rate densities. The schematic
`etadot^inf = -Z R + 2 W Omega Omegadot` of ADR-0013 §3.4 is thereby confirmed with `R = (R_e, R_mu)^T`.

Here `D` is the declared chemical/reaction domain consistent with the associated `G_y`, and each
channel integral is restricted to its declared support/applicability subset `D_a subseteq D`; no
implicit core boundary is permitted. Units check:

```text
[Ltilde_a/k_B^(erg)] T_inf^(q_a-1)
  = (erg s^-1 K^-q)/(erg K^-1) K^(q-1) = s^-1 == count/s
[Z][R] = (MeV/count)(count/s) = MeV/s
[W][Omega][Omegadot] = (MeV s^2)(s^-1)(s^-2) = MeV/s.
```

### 6.2 Derivation (not assumed)

Start from the two independently authenticated pieces.

**(a) Kinematics.** R2006 eq. (15) / Y2020 eq. (4.7):
`eta_npe^inf = -Z_npe delta N_e - Z_np delta N_mu`, `eta_npmu^inf = -Z_np delta N_e - Z_npmu delta N_mu`,
i.e. `eta^inf = -Z delta N_l` with `Z` symmetric. `Z` is time-independent for frozen coefficients
(§19), so

```text
etadot^inf = - Z  d(delta N_l)/dt .
```

**(b) Number balance.** FR2005 eqs. (14)–(15) and (30); R2006 eq. (7); Y2020 eqs. (4.8)–(4.9):

```text
delta N_i  =  N_i - N_i^eq
d(delta N_l)/dt = Ndot_l - Ndot_l^eq
                = int dV e^Phi sum_a DeltaGamma_{l,a}  -  2 Omega Omegadot I_{Omega,l}
                = R_l - 2 Omega Omegadot I_{Omega,l} .
```

**(c) Combine.**

```text
etadot^inf = - Z R + 2 Omega Omegadot ( Z I_Omega ) .
```

**(d) Identify `W`.** R2006 eq. (19): `W_npe = Z_npe I_{Omega,e} + Z_np I_{Omega,mu}`,
`W_npmu = Z_np I_{Omega,e} + Z_npmu I_{Omega,mu}`, i.e. exactly `W = Z I_Omega`. Substituting gives
the boxed system. **This closes the derivation without assuming FR2005 eqs. (52)–(53).**

**Consistency with the superseded FR2005 form.** FR2005 eqs. (57)–(58) write
`W_npe = (Z_npe - Z_np) I_{Omega,e} + Z_np I_{Omega,p}`. With charge neutrality
`I_{Omega,p} = I_{Omega,e} + I_{Omega,mu}` (FR2005 p. 14) this is algebraically identical to
`W = Z I_Omega`. R2006 §5, p. 570, states this explicitly: the `W_npl` *"also remain correctly
expressed by equations (57) and (58), although they are more easily written as"* eq. (19).
**No discrepancy.**

**Numerical confirmation against the governed Phase-5C artifact.** From
`docs/validation/phase5c_chemical_coefficients_candidate.json` (free-gas Structure-1, whole star):

```text
Z      = [[4.5793031807026964e-54, 5.1725199102788050e-55],
          [5.1725199102788054e-55, 1.0268727975139168e-52]]   MeV/count
I      = [-1.1637998545112904e+47, -1.4844819850233820e+46]   count s^2
W      = [-5.4061775017047240e-07, -1.5845719480103649e-06]   MeV s^2
```

Recomputing `Z I` gives `[-5.40617750e-07, -1.58457195e-06]`, and solving `Z x = W` recovers
`x = [-1.16379985e+47, -1.48448199e+46]`. `W = Z I_Omega` holds in the shipped artifact to full
working precision. Scratch: `scratchpad/calc/evolve.py`.

### 6.3 Three-witness check on the `eta = -Z delta N` sign

| Witness | Statement |
|---|---|
| Independent derivation (this preflight) | Substituting `eta^inf = -Z delta N` and `W = Z I_Omega` into FR2005 eqs. (52)–(53) makes the `2 W Omega Omegadot` terms cancel identically against `-Z * (-2 Omega Omegadot I_Omega)`, leaving `etadot^inf = -Z d(delta N)/dt`. Any other sign fails this cancellation. |
| R2006 eq. (15), p. 570 | Printed relation (its glyphs are degraded in text extraction; the *structure* `eta_npl = ∓Z_.. delta N_e ∓ Z_.. delta N_mu` is unambiguous and its sign is fixed by the cancellation above). |
| **Y2020 eq. (4.7), p. 59** | Clean text: `eta_e^inf = -Z_npe delta N_e - Z_np delta N_mu`, `eta_mu^inf = -Z_np delta N_e - Z_npmu delta N_mu`. |

The already-ratified ADR-0013 §3.3 convention is therefore **source-confirmed**, not merely
internally consistent.

### 6.4 `Z` cross-coupling map

```text
Z = [[Z_npe , Z_np  ],        row  = OUTPUT beta channel  (Npe, NpMu)
     [Z_np  , Z_npmu]]        col  = INPUT  lepton number (e,   mu  )
```

`Z_np` is the **shared** off-diagonal: electron-channel reactions feed the muon-channel imbalance and
vice versa, through the shared neutron/proton reservoir. Dropping it is mutation M7 (§32).
`Z_npe = Btilde^-1_nn - 2 Btilde^-1_ne + Btilde^-1_ee`,
`Z_npmu = Btilde^-1_nn - 2 Btilde^-1_nmu + Btilde^-1_mumu`,
`Z_np = Btilde^-1_nn - Btilde^-1_ne - Btilde^-1_nmu + Btilde^-1_emu`
(R2006 eqs. (16)–(18)); in CompactStar these are read-only views onto the one stored matrix
(ADR-0013 §3.3), never independently stored.

### 6.5 `W` sign and the spin-down direction

`I_{Omega,l} = (dN_l^eq/dOmega^2)_A < 0` (FR2005 p. 14: *"since `Omega Omegadot` and `I_{Omega,i}`
are negative"*); `Z` is positive definite; hence `W = Z I_Omega < 0` componentwise. The governed
free-gas `W = (-5.406e-7, -1.585e-6) MeV s^2` confirms both components negative.

For a spinning-down star, `Omega > 0` and `Omegadot < 0`, so `Omega Omegadot < 0` and

```text
2 W Omega Omegadot  >  0    (both channels)
```

i.e. **spin-down drives both redshifted imbalances positive**, exactly as FR2005 states on p. 13
(*"equations (52) and (53) have a positive term ... which makes the chemical imbalances grow"*) and
p. 14 (*"Both chemical imbalances at quasi-equilibrium are positive"*).

---

## 7. Reaction-rate sign convention

### 7.1 Canonical definition

For a beta channel `a` with lepton `l`, take the forward direction to be **neutron decay**:

```text
A -> B :    n  ( + N_1 )  ->  p  ( + N_2 )  +  l  +  nubar_l
B -> A :    p  ( + N_2 )  +  l  ->  n  ( + N_1 )  +  nu_l          (lepton capture)

DeltaGamma_a  ==  Gamma_{A->B} - Gamma_{B->A}         [count cm^-3 s^-1]
eta_a         ==  delta mu(A) - delta mu(B) = delta mu_n - delta mu_p - delta mu_l
```

This is FR2005 eqs. (39)–(40) applied with `A = {n, N_1}`, `B = {p, N_2, l}`, which reproduces
FR2005 eqs. (7)–(8) for `eta` — so the pairing is fixed by the sources, not chosen.

### 7.2 What `DeltaGamma_a > 0` means

`DeltaGamma_a > 0` means **net neutron decay**: neutrons are destroyed, protons and leptons `l` are
created. Authenticated by Y2020 eq. (4.8), which is explicit about every sign:

```text
Ndot_n = - sum_l sum_N int dV DeltaGamma_{M,Nl} e^Phi          (neutrons destroyed)
Ndot_p = + sum_l sum_N int dV DeltaGamma_{M,Nl} e^Phi          (protons created)
Ndot_l = +       sum_N int dV DeltaGamma_{M,Nl} e^Phi          (leptons l created)
```

Consistency with FR2005 eq. (15) (`Ndot_i = int dV e^Phi sum_a DeltaGamma_{i,a}`, with
`DeltaGamma_{i,a}` the net creation rate of species `i`) is exact: `DeltaGamma_{l,a} = +DeltaGamma_a`,
`DeltaGamma_{p,a} = +DeltaGamma_a`, `DeltaGamma_{n,a} = -DeltaGamma_a`. Baryon number and charge
are conserved term by term, as FR2005 p. 4 asserts.

### 7.3 Positive `DeltaGamma` destroys positive `eta` — proof

Yes. Chain:

1. `DeltaGamma_a = [1/(k_B^(erg)T)] Q_a^eq H_*(xi_a)` with `H_*` **odd** and `x H_*(x) >= 0` (§9.4), so
   `sign(DeltaGamma_a) = sign(eta_a)`.
2. `R_l > 0` for `eta_l > 0`, so `d(delta N_l)/dt` receives a positive reaction contribution.
3. `eta^inf = -Z delta N_l` with `Z` positive definite, so the reaction part of the RHS is
   `etadot^inf|_react = -Z R`.
4. Consider the positive-definite quadratic form `V = eta^T Z^-1 eta`. Then

```text
dV/dt |_react  =  2 eta^T Z^-1 ( -Z R )  =  -2 eta . R  =  -2 sum_l eta_l R_l  <=  0 ,
```

The derivative is nonpositive. Strict decrease requires that every nonzero imbalance direction
being considered couple to at least one physically active channel with positive normalization and
nonzero dissipative response. If all applicable channel normalizations for an uncoupled direction
vanish (`sum_{a in l} Ltilde_a = 0`), that **uncoupled dead imbalance** may freeze: `V` is then only
nonincreasing and the reaction operator is negative-semidefinite. With nonzero `Z` cross-coupling,
a dead individual reaction channel does not imply that its `eta` component remains at its initial
value. Tested at both signs of `eta` in §21 (limit F).

### 7.4 Thermodynamic dissipation

`eta_a DeltaGamma_a >= 0` pointwise. In the chosen mixed-unit contract,
```text
C_(MeV->erg) eta_a[MeV] DeltaGamma_a[count cm^-3 s^-1]
  = Q_a^eq[erg cm^-3 s^-1] xi_a H_*(xi_a) >= 0.
```
The positive conversion does not change the sign.
Numerically verified at `x = 0.3, 1, 5, 50` for both `H_D` and `H_M` (all strictly positive;
`x H(x) = 0` only at `x = 0`). In the source's single-energy-unit notation,
`Q_H = sum_a DeltaGamma_a eta_a >= 0` (FR2005 eq. 38);
the chosen implementation contract applies `C_(MeV->erg)` once at the global thermal boundary.
Chemical heating can never be negative. Y2020 §4.1 gives the thermodynamic reason: in the isolated limit
`T^inf dS = L_H^inf dt`, and the second law forces `L_H^inf >= 0`.

### 7.5 No hidden sign flip

The three places a flip could hide are pinned as follows.

| Interface | Convention | Anchor |
|---|---|---|
| paper `DeltaGamma` | `Gamma_{n->p+l+nubar} - Gamma_{p+l->n+nu}`, positive for net neutron decay | FR2005 eqs. (39)–(40); Y2020 eq. (4.8) |
| `H_*` | odd, `H_*(x) > 0` for `x > 0`, `H_*(0) = 0` | §9; FR2005 p. 10 (*"reaction rates are enhanced in the direction which restores equilibrium, since the functions `H_*` are odd"*) |
| chemical ODE | `etadot = -Z R + 2 W Omega Omegadot` with `R_l` the **global lepton-creation** rate | §6 |

**Warning recorded for implementers.** Y2020 footnote 3, p. 61, states its `I^N_{M,Gamma}` has the
*opposite* sign to two of its own cited references, and R1995 footnote 3, p. 14, states its `eta` is
opposite to Sawyer (1989) and Haensel (1992). Any future rate source must be re-anchored to the
table above before use; the sign in a formula copied out of context is not trustworthy.

---

## 8. Redshifted imbalance / temperature ratio

### 8.1 Derivation

From FR2005 eq. (1), `T(r) e^{nu(r)} = T_inf`; from FR2005 eq. (9), `eta_a(r) e^{nu(r)} = eta_a^inf`.
Therefore

```text
xi_a(r)  ==  eta_a(r) / (k_B^(MeV) T(r))
         =  [eta_a^inf e^{-nu(r)}] / (k_B^(MeV) T_inf e^{-nu(r)})
         =  eta_a^inf / (k_B^(MeV) T_inf)      —  independent of r.
```

**Authenticated.** FR2005 eq. (42), p. 11, states exactly `xi_a == eta_a/(kT) = eta_a^inf/(k T_inf)`;
Y2020 eq. (4.18) defines `xi_l == eta_l/T` in units with `k_B = 1`. The cancellation is *exact* and
requires **both** redshift conventions to use the same `e^{nu}`. This spatial constancy is what
allows `F_*(xi)` and `H_*(xi)` to be pulled out of the stellar integrals (FR2005 eqs. 44–45) —
it is not a convenience, it is the structural reason the whole `Ltilde` factorization works.

### 8.2 Units and constant ownership

```text
xi_a  =  eta_a^inf [MeV] / ( k_B^(MeV) [MeV/K] * T_inf [K] )  (dimensionless)
k_B^(MeV) = Zaki::Physics::K_BOLTZ_EV * 1e-6                    MeV/K
C_(MeV->erg) = Units::MEV_FM3_TO_ERG_CM3 / 10^39               erg/MeV
k_B^(erg) = C_(MeV->erg) k_B^(MeV)                              erg/K
```

Audited: this is already the single Boltzmann authority
(`NeutrinoCooling_Details.cpp:213-220`, `docs/validation/PHASE3C_BOLTZMANN_AUTHORITY.md`), which
replaced two divergent local literals in Phase 3C. `C_(MeV->erg)` is derived from the governed
energy-density conversion and the exact `1 fm^-3 = 10^39 cm^-3`; it is not an independent energy
literal. **No new Boltzmann literal may be introduced.** The MeV/K view owns `xi`; the derived
erg/K view owns the reaction-rate formula against the erg-normalized `Ltilde`.

---

## 9. Exact non-superfluid Urca functions

### 9.1 Independent derivation from the phase-space integrals

The functions were **derived**, not copied. R1995 Appendix, p. 14–15, defines

```text
F(x) = int_0^inf dy  y^3 P(y-x) / (1 + exp(y-x))
G(x) = int_0^inf dy  y^2 P(y-x) / (1 + exp(y-x))
F_+(x) = F(x) + F(-x)              (even by construction)
G_-(x) = G(x) - G(-x)              (odd  by construction)
P_D(y) = pi^2 + y^2                        (direct Urca)
P_M(y) = 9 pi^4 + 10 pi^2 y^2 + y^4        (modified Urca)
```

and (FR2005 eqs. 32–33 restated) `F_*(x) = F_+(x)/F_+(0)`, `H_*(x) = G_-(x)/F_+(0)`.

Evaluating with the Bernoulli identity R1995 eq. (28),
`int_0^inf dy y^{2n-1}/(1+e^y) = (1 - 2^{1-2n}) (2 pi)^{2n} |B_{2n}| / (4n)`:

```text
int y  /(1+e^y) = pi^2/12       int y^3/(1+e^y) = 7 pi^4/120
int y^5/(1+e^y) = 31 pi^6/252   int y^7/(1+e^y) = 127 pi^8/240

F_+(0)_D = 2 * ( pi^2 * 7pi^4/120 + 31pi^6/252 )                  =  457 pi^6 / 1260
F_+(0)_M = 2 * ( 9pi^4 * 7pi^4/120 + 10pi^2 * 31pi^6/252
                 + 127 pi^8/240 )                                  = 11513 pi^8 / 2520
```

**The paper constants 457 and 11513 are reproduced exactly by this derivation.** Carrying the same
algebra through the finite integrals `int_0^x (y-x)^3 P(y) dy` and `int_0^x (y-x)^2 P(y) dy`
yields the four polynomials below with every coefficient exact.

**Independent numerical confirmation.** `scratchpad/calc/urca_derive.py` evaluates `F_+(x)/F_+(0)`
and `G_-(x)/F_+(0)` by 40-digit quadrature at `x = 0.5, 1, 2, 4, 8` and compares to the polynomials:
**ratio = 1.0 to 18 significant digits at all five points for all four functions.** Five points
over-determine a degree-6 and a degree-8 polynomial, so the identification is exact. The wrong
modified-Urca weight `P_M = (9/4)pi^4 + 10 pi^2 y^2 + y^4` was also tested and *fails* (ratios
0.66–0.99), confirming `9 pi^4`.

### 9.2 The four functions

```text
F_D(xi) = 1 + 1071 xi^2/(457 pi^2) + 315 xi^4/(457 pi^4) + 21 xi^6/(457 pi^6)

H_D(xi) =      714 xi  /(457 pi^2) + 420 xi^3/(457 pi^4) + 42 xi^5/(457 pi^6)

F_M(xi) = 1 + 22020 xi^2/(11513 pi^2) + 5670 xi^4/(11513 pi^4)
            +   420 xi^6/(11513 pi^6) +    9 xi^8/(11513 pi^8)

H_M(xi) =     14680 xi  /(11513 pi^2) + 7560 xi^3/(11513 pi^4)
            +   840 xi^5/(11513 pi^6) +   24 xi^7/(11513 pi^8)
```

`F_*` multiplies the **neutrino emissivity** (FR2005 eq. 32); `H_*` multiplies the **net reaction
rate** (FR2005 eq. 33). They are not interchangeable, and the `1/kT` prefactor belongs to
`DeltaGamma` only.

### 9.3 FR2005 eq. (37) — confirmed printed typo / internal source inconsistency

**CONFIRMED PRINTED TYPO / INTERNAL SOURCE INCONSISTENCY; NO PUBLISHED ERRATUM LOCATED.** As printed,
FR2005 eq. (37) gives the last `H_M` term as `24 xi^7 / (11513 pi^6)`. The proposed normative
implementation formula uses `11513 pi^8`. Independent evidence:

1. **Literal source:** FR2005 prints `pi^6` in eq. (37).
2. **Exact independent phase-space/Fermi-convolution derivation** (§9.1): the `xi^7` coefficient is `(1/105)(2520/11513)/pi^8`, i.e.
   `24/(11513 pi^8)`.
3. **Numerical quadrature** (§9.1): the `pi^8` form matches `G_-(x)/F_+(0)` to 18 digits; the `pi^6`
   form does not.
4. **FR2005's later coefficient structure**, eq. (60), p. 14: `M_M(x) ~= 15 x^8/(11513 pi^8)`. Since `M_M = xi H_M - F_M` and
   `F_M`'s leading term is `9 x^8/(11513 pi^8)`, the leading `xi H_M` term must be
   `24 x^8/(11513 pi^8)` — which requires `pi^8` in eq. (37). The `pi^6` reading gives
   `24 pi^2 - 9 != 15`.
5. **R1995 source-variable form**, eq. (32), gives the consistent final coefficient after
   `u=xi/pi`.
6. **Y2020 eq. (4.20)**, p. 62, prints `24 xi^7/(11513 pi^8)`.

FR2005's own downstream results (eqs. 59–66, the `5/8` heating fraction, Fig. 2) all use the correct
`pi^8`. A bounded search found no published erratum. The correction remains proposed and normative;
it is not classified as an authenticated erratum.

### 9.4 Verified properties

| Property | Result |
|---|---|
| `F_D(0) = F_M(0) = 1` | exact (constant term) |
| `H_D(0) = H_M(0) = 0` | exact (no constant term) |
| `F_*` even | exact — only even powers; numerically `\|F(3) - F(-3)\| = 0` |
| `H_*` odd | exact — only odd powers; numerically `\|H(3) + H(-3)\| = 0` |
| `xi H_*(xi) >= 0` | all coefficients of `H_*` positive ⇒ `xi H_*(xi) = sum c_k xi^{2k}` with `c_k > 0`; verified at `xi = 0.3, 1, 5, 50` |
| `F_* > 0` | all coefficients positive |

### 9.5 Small- and large-`xi` behaviour

**Small `xi`** (linear response):

```text
H_D(xi) -> [714/(457 pi^2)] xi   = 0.158300492605 xi
H_M(xi) -> [14680/(11513 pi^2)] xi = 0.129192649689 xi
F_D(xi) -> 1 + 0.237450738908 xi^2
F_M(xi) -> 1 + 0.193788974534 xi^2
```

So the reaction term is **linear** in `eta` at small imbalance, and the neutrino enhancement is
**quadratic** — the chemical relaxation is a linear (matrix-exponential) problem near equilibrium
(§21 limit D).

**Large `\|xi\|`:**

```text
F_D -> 21 xi^6/(457 pi^6)     H_D -> 42 xi^5/(457 pi^6)      M_D -> 21 xi^6/(457 pi^6)
F_M ->  9 xi^8/(11513 pi^8)   H_M -> 24 xi^7/(11513 pi^8)    M_M -> 15 xi^8/(11513 pi^8)

heating fraction  C_M/C_H  =  21/42 = 1/2   (direct)      15/24 = 5/8   (modified)
neutrino fraction 1 - C_M/C_H =  21/42 = 1/2 (direct)     9/24 = 3/8 (modified)
```

matching FR2005 p. 11 verbatim (*"a fixed fraction of the energy released is emitted as neutrinos,
3/8 and 1/2 for modified and direct Urca"*) and FR2005 eqs. (59)–(60)
(`C_H = 24/(11513 pi^8)`, `C_M = 15/(11513 pi^8)`).

### 9.6 Two-source coefficient check

| Function | FR2005 eqs. (34)–(37) | R1995 eqs. (29)–(32), `u = xi/pi` | Y2020 eqs. (4.19)–(4.20) | Agree? |
|---|---|---|---|---|
| `F_D` | `1071, 315, 21 / 457` | `1071, 315, 21 / 457` | — | yes |
| `H_D` | `714, 420, 42 / 457` | `714, 420, 42 / 457` | — | yes |
| `F_M` | `22020, 5670, 420, 9 / 11513` | `22020, 5670, 420, 9 / 11513` | `22020, 5670, 420, 9 / 11513` | yes |
| `H_M` | `14680, 7560, 840, 24 / 11513` (last `pi` exponent mis-set) | `14680, 7560, 840, 24 / 11513` | `14680, 7560, 840, 24 / 11513`, `pi^8` | yes |

R1995 uses `u = eta/(pi k T)`; FR2005 and Y2020 use `xi = eta/(kT) = pi u`. The presentations are
identical under that substitution. **No primary-authority disagreement.**

Cross-check on the combined functions: R1995's own cooling equations (p. 8) read
`Tdot_8 ∝ T_8^5 (21u^6 + 105u^4 - 357u^2 - 457)/457` (direct) and
`Tdot_8 ∝ T_8^7 (15u^8 + 420u^6 + 1890u^4 - 7340u^2 - 11513)/11513` (modified). Forming
`M_*(xi) = xi H_*(xi) - F_*(xi)` from §9.2 with `u = xi/pi` gives **exactly** those two polynomials.
This is an independent primary-source confirmation of both families simultaneously.

---

## 10. Rate normalization

### 10.1 The single normalization identity

Because `xi` is spatially constant (§8) and `Q_a^eq = S_a(n) T^q` factorizes (FR2005 eq. 41), **one**
coefficient governs four different global quantities:

```text
Ltilde_a  ==  integral_{D_a} 4 pi r^2 e^Lambda S_a(n) e^{(2-q) Phi} dr     [erg s^-1 K^-q]

L_{nu,a}^inf (T_inf, xi_a)      =  Ltilde_a F_*(xi_a)      T_inf^q             (FR2005 eq. 44)
R_a = integral_{D_a} DeltaGamma_a e^Phi dV
                                  = [Ltilde_a/k_B^(erg)] H_*(xi_a) T_inf^{q-1} [count/s]
L_{H,a}^inf                     =  Ltilde_a xi_a H_*(xi_a) T_inf^q             (FR2005 eq. 46)
L_{nu,a,eq}^inf                 =  Ltilde_a                T_inf^q             (xi_a = 0)
```

Derivation of eq. (45), reproduced because it is the load-bearing one:

```text
integral_{D_a} DeltaGamma_a e^Phi dV
  = integral_{D_a} dV e^Phi (1/(k_B^(erg) T_local)) S_a T_local^q H_*(xi_a)
  = [1/k_B^(erg)] H_*(xi_a) T_inf^{q-1} integral_{D_a} dV S_a e^{(2-q)Phi}
  = [Ltilde_a/k_B^(erg)] H_*(xi_a) T_inf^{q-1}
```

`H_*(xi_a)` leaves the integral **only** because `xi_a` is spatially constant. `q = 6` (direct),
`q = 8` (modified). Dimensionally,

```text
(erg s^-1 K^-q)/(erg K^-1) K^(q-1) = s^-1 == count/s.
```

The same channel normalization feeds equilibrium `L_nu`, nonequilibrium `L_nu`/`DeltaL_nu`, the
net reaction rate, and chemical-heating bookkeeping. The chemical module's canonical power is
`P_H,chem^inf = sum_l eta_l^inf[MeV] R_l[count/s]` in `MeV/s`; it is converted exactly once by
`C_(MeV->erg)` when contributed to the thermal RHS in `erg/s`. Nonequilibrium neutrino luminosity
remains `erg/s`. Algebraically,
`C_(MeV->erg) eta R = Ltilde xi H T_inf^q` because
`k_B^(erg)=C_(MeV->erg)k_B^(MeV)` and `eta=xi k_B^(MeV)T_inf`.

**This identity is the whole no-double-counting architecture** (§13): equilibrium cooling and every
rotochemical correction are the *same* `Ltilde_a` evaluated at different `xi`.

### 10.2 Audit of the current CompactStar equilibrium normalization

`CompactStar/Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp:60-207` builds

```text
K_MU_erg_s_K8 = Q0_MU * (1e-9)^8 * KM3_TO_CM3 * int_0^{R}   rho15 e^{-8 nu} (4 pi r^2 e^Lambda) dr
K_DU_erg_s_K6 = Q0_DU * (1e-9)^6 * KM3_TO_CM3 * int_0^{r_DU} rho15 e^{-6 nu} (4 pi r^2 e^Lambda) dr
L_nu^inf      = K_DU T_inf^6 + K_MU T_inf^8
```

with (`:101-103`)

```c
// Placeholder normalizations (must match your emissivity model).
constexpr double Q0_DU = 1.0e27; // erg cm^-3 s^-1 * rho15 * T9^6
constexpr double Q0_MU = 1.0e21; // erg cm^-3 s^-1 * rho15 * T9^8
```

| Aspect | Finding |
|---|---|
| **Structure** | `K_n = int dV A(r) e^{(2-n)nu}` with `L = K_n T_inf^n` — **structurally identical to FR2005 eq. (43)**. The geometric/redshift half is already right, including `e^{(2-q)Phi}` and the proper-volume factor `4 pi r^2 e^Lambda`. |
| **Normalization** | `Q = Q0 * (rho/1e15 g cm^-3) * T9^n`, with `Q0` explicitly labelled *"Placeholder"*. `rho` is the **mass-energy density** (`StarContext::MassDensity_gcm3`), not the source-prescribed `(x_eq n / n_0)^{1/3}` of R1995 eqs. (33)–(34). **Not source-authoritative.** |
| **Electron/muon split** | **Absent.** One lumped DU and one lumped MU coefficient. |
| **Neutron/proton branch split** | **Absent.** |
| **Composition dependence** | **Absent** (only mass density). |
| **Redshift convention** | Correct: exactly `e^{(2-n)nu}`, one proper-volume factor. |
| **Support domain** | DU restricted to `[0, durca_last]`, MU over the whole profile. |
| **Governance status** | `docs/SCIENTIFIC_INVARIANTS.md` explicitly records: *"Placeholder emissivities and arbitrary normalizations are recorded, not blessed."* |

### 10.3 Can the existing coefficient be reused?

**Structurally: yes, and it must be.** The cached-`K` mechanism, its geometry, its redshift powers,
its profile-versioned cache, and its `L = K T^n` evaluation are exactly the `Ltilde` mechanism and
should be extended, not duplicated.

**Numerically: no.** The stored `Q0_DU`, `Q0_MU` are declared placeholders with the wrong density
dependence and no lepton resolution. They cannot be reused as a source-faithful `S_a(n)`.

**Conversion, if a source-faithful `S_a` is supplied.** Nothing else is needed:

```text
S_a(n)  [erg cm^-3 s^-1 K^-q]  ==  Q_a^eq(n, T) / T^q
Ltilde_a = integral_{D_a} 4 pi r^2 e^Lambda S_a(n) e^{(2-q)nu} dr
             (declared reaction/chemical domain; km^3 -> cm^3 once)
```
and then §10.1 supplies `DeltaGamma`, `L_nu(eta)`, and `L_H` with **no further microphysical input**.
The only missing ingredient is `S_a(n)` itself, per channel.

### 10.4 What `S_a(n)` authority exists

| Source | Provides | Status |
|---|---|---|
| **R1995 eqs. (33)–(34)**, p. 15 (primary; citing Haensel 1992) | `eps_d(T,0) ~= 4.3e21 (x_eq n/n_0)^{1/3} T_8^6`, `eps_m(T,0) ~= 3.5e13 (x_eq n/n_0)^{1/3} T_8^8 erg cm^-3 s^-1`, `n_0 = 0.16 fm^-3` | **Historical/supporting normalization only.** It is npe-lumped and the source calls it *"somewhat uncertain"*. It is not sufficient as the definitive FR2005 realistic channel normalization. |
| **FR2005 §3.2**, p. 11 | Declares `alpha = De, Dmu, Me, Mmu`, each *"adding the contributions of the neutron and proton branches (Yakovlev et al. 2001)"* | **Delegates** the branch- and lepton-resolved normalization to Yakovlev et al. (2001). |
| **Yakovlev et al. 2001 (YKGH2001)**, Phys. Rep. 354, 1, arXiv:astro-ph/0012122 | Modified-Urca neutron-branch normalization, proton-branch factor, electron/muon substitutions, threshold structure, and effective-mass dependence | **Publicly available and sufficient in form** to define the channel-normalization architecture. Not installed/authenticated here. The exact `alpha_n` choice matching the intended FR2005 reproduction remains unresolved: an FM79-style constant and later density-dependent/OPE prescriptions are plausible alternatives. No convention is selected here. |
| **FR2005 §3.4**, p. 12 | Effective masses *"can be obtained analytically for the APR and PAL EOSs (see, e.g., Page et al. 2004)"*; **for the noninteracting Fermi gas `m*_i = mu_i/c^2`** | Free-gas effective masses are **fully specified by a primary source**. Realistic ones are not. |
| **Y2020 eqs. (2.75)–(2.76)**, p. 29 (supporting) | Branch- and lepton-resolved modified Urca: `Q^(0)_{M,nl} = 8.05e21 v_{F,l} (m*_n/m_n)^3 (m*_p/m_p) (p_{F,p}/k_0) T_9^8 alpha beta`, `Q^(0)_{M,pl} = Q^(0)_{M,nl} (m*_p/m*_n)^2 (p_{F,l}+3p_{F,p}-p_{F,n})^2/(8 p_{F,l} p_{F,p}) Theta(...)`, `alpha = 1.76 - 0.63 (n_0/n_n)^{2/3}`, `beta = 0.68` | **SUPPORTING SOURCE, NOT ADOPTED.** This reproduces the Yakovlev-family formulas FR2005 delegates to. It may be used only as a *declared benchmark model input* with that classification recorded, never as FR2005 reproduction authority. |

**Conclusion.** R1995 is useful historical support but is not the definitive FR2005 realistic
normalization. YKGH2001 is sufficient in form for the architecture, while authentication/install
and the exact `alpha_n` prescription remain realistic-source blockers (§27). A controlled free-gas
benchmark may instead use declared positive coefficients with that non-realistic classification.

---

## 11. Global GR integrals

Metric `ds^2 = -e^{2 Phi} dt^2 + e^{2 Lambda} dr^2 + r^2 dOmega^2`, `nu == Phi`,
`e^Lambda = (1 - 2m/r)^{-1/2}`, proper volume `dV = 4 pi r^2 e^Lambda dr`.

Every normative reaction/luminosity coefficient integral is `integral_D` or
`integral_{D_a}`. `D` is the declared chemical/reaction domain consistent with the associated
`G_y`; `D_a subseteq D` is each channel's support/applicability subset. For the free-gas fixture,
`D` is exactly the connected whole-star, source-valid Phase-5C domain from the centre through the
authenticated `npemu`/`npe`/`pe` branch partition to the governed physical neutron-onset/vacuum
boundary, including the accepted refusal-window and tail treatment. No arbitrary core, crust, or
saturation-density cutoff is introduced.

Every factor below is derived from proper time, proper volume and energy redshift; none is copied.

- **Proper time.** A local rate per unit proper time `tau` contributes, per unit coordinate time `t`,
  a factor `dtau/dt = e^{Phi}`. Local rates therefore acquire **one `e^{+Phi}`**.
- **Energy redshift.** A quantum of local energy `E` arrives at infinity with `E_inf = e^{Phi} E`.
  Energy fluxes therefore acquire **a second `e^{+Phi}`**.
- **Redshifted potential.** A local intensive potential reconstructed from a redshifted one carries
  `X_local = X^inf e^{-Phi}`, giving **`e^{-Phi}`**.
- **Local temperature.** `T_local = T_inf e^{-Phi}`, so an integrand `∝ T_local^q` contributes
  `e^{-q Phi}`.

| Quantity | Local integrand | Proper-volume factor | Time-redshift | Energy-redshift | Other | **Net factor** | Source |
|---|---|---|---|---|---|---|---|
| Heat capacity `C_*(T_inf)` | `c_V(T_local)` | `4 pi r^2 e^Lambda` | `e^{+Phi}` (energy at infinity) | — | `c_V ∝ T_local = T_inf e^{-Phi}` | **`e^{0}`** (net 1), with `T_local` inside | FR2005 eqs. (3),(50); ADR-0002 |
| Global net reaction rate `Ndot_i` | `DeltaGamma_{i,a}` | `4 pi r^2 e^Lambda` | `e^{+Phi}` | — | — | **`e^{+Phi}`** | FR2005 eq. (15); R2006 eq. (6); Y2020 eq. (4.8) |
| Neutrino luminosity `L_nu^inf` | `Q_nu` | `4 pi r^2 e^Lambda` | `e^{+Phi}` | `e^{+Phi}` | — | **`e^{+2Phi}`** | FR2005 eq. (5) |
| Chemical heating `L_H^inf` | `Q_H = sum_a DeltaGamma_a eta_a^local` | `4 pi r^2 e^Lambda` | `e^{+Phi}` | `e^{+Phi}` | `eta_local = eta^inf e^{-Phi}` | **`e^{+2Phi}` on `Q_H`, equivalently `e^{+Phi}` on `DeltaGamma` with `eta^inf` outside** | FR2005 eqs. (4),(38); Y2020 eq. (4.5) |
| Chemical response `G_y` / `B_ij` | `dn_i/dmu_j` | `4 pi r^2 e^Lambda` | — | — | `delta mu_local = delta mu^inf e^{-Phi}` | **`e^{-Phi}`** | FR2005 eq. (12); R2006 eq. (13); ADR-0013 §3.2 |
| Luminosity coefficient `Ltilde_a` | `S_a(n)` | `4 pi r^2 e^Lambda` | `e^{+Phi}` | `e^{+Phi}` | `T_local^q = T_inf^q e^{-q Phi}` | **`e^{(2-q)Phi}`** | FR2005 eq. (43) |
| Photon luminosity `L_gamma^inf` | surface | — | — | — | `4 pi sigma R_inf^2 (T_s^inf)^4`, `R_inf = R e^{-Phi_s}`, `T_s^inf = T_s e^{Phi_s}` | `= 4 pi sigma R^2 T_s^4 e^{2 Phi_s}` | FR2005 eq. (6) |

**The heating identity.** Because `eta^inf` is spatially uniform, the `e^{2Phi}` on `Q_H` and the
`e^{-Phi}` inside `eta_local` collapse:

```text
P_H,chem^inf [MeV/s]
        = sum_a eta_a^inf[MeV] integral_{D_a} dV e^{Phi} DeltaGamma_a[local count/(volume time)]
        = sum_a eta_a^inf R_a .
L_H^inf [erg/s] = C_(MeV->erg) P_H,chem^inf .
```

The canonical chemical module owns `P_H,chem^inf` in MeV/s and converts it exactly once at the
thermal-luminosity boundary. The governed km^3-to-cm^3 conversion is already part of the coefficient
integration and is not repeated here. This reproduces FR2005 eq. (46) in erg/s on substituting eq. (45).

**Warning: the `e^{-Phi}` of `G_y` and the `e^{+Phi}` of the reaction integral are different
operations on different integrands.** ADR-0013 §3.2 already says so (*"The positive lapse used for
reaction-rate conversion is a different future operation and does not alter the coefficient
integral"*). Confusing them is mutation M10/M11 (§32).

---

## 12. Thermal ledger

### 12.1 Two different questions

The sources contain **two** differences that must never be confused:

```text
FULL         :  L_H^inf - L_nu^inf(T, eta)                 = Ltilde_a M_*(xi_a) T_inf^q ,  M_* = xi H_* - F_*
INCREMENTAL  :  L_H^inf - [ L_nu^inf(T,eta) - L_nu^inf(T,0) ]
                                                            = Ltilde_a [ M_*(xi_a) + 1 ] T_inf^q
```

`M_*` (FR2005 eq. 47) has the equilibrium cooling **already subsumed in it** — the constant term of
`F_*` is `1`, and FR2005 p. 11 says so directly: *"With `xi = 0`, the constant term in `F_*` gives
the conventional cooling case."* So FR2005 eq. (48) answers *"heating minus **all** neutrino
emission"*, whereas the incremental question is *"what does non-zero `eta` add, relative to the same
star cooling passively at the same `T_inf`?"* They differ by exactly `Ltilde_a T_inf^q`, i.e. by the
whole equilibrium cooling luminosity. **Both are correct; they are not interchangeable.**

### 12.2 Complete ledger

For each channel `a` with lepton `l` and exponent `q`:

```text
L_{nu,a}^inf(T_inf, xi)   =  Ltilde_a F_*(xi_a) T_inf^q                        (>= Ltilde_a T_inf^q)
L_{nu,a,eq}^inf(T_inf)    =  Ltilde_a T_inf^q                                  (xi = 0)
DeltaL_{nu,a}^inf         =  Ltilde_a [F_*(xi_a) - 1] T_inf^q                  (>= 0, even in xi)
L_{H,a}^inf               =  Ltilde_a xi_a H_*(xi_a) T_inf^q                   (>= 0, even in xi)
DeltaP_{beta,a}           =  L_{H,a}^inf - DeltaL_{nu,a}^inf
                          =  Ltilde_a [ xi_a H_*(xi_a) - F_*(xi_a) + 1 ] T_inf^q
                          =  Ltilde_a [ M_*(xi_a) + 1 ] T_inf^q
```

Total: `L_H^inf = sum_a L_{H,a}^inf`, `DeltaL_nu^inf = sum_a DeltaL_{nu,a}^inf`.

Note `L_H` and `DeltaL_nu` are both **even** in `xi` (since `H_*` is odd and `F_*` even), so the
ledger is insensitive to the sign of the imbalance, as FR2005 p. 11 states of `M_*`.

### 12.3 Chosen thermal-RHS architecture

**Incremental.** The thermal RHS gains, per rotochemical channel,

```text
+ L_{H,a}^inf   and   - DeltaL_{nu,a}^inf ,
```

leaving the equilibrium `NeutrinoCooling` path as the sole owner of the equilibrium term for **all**
channels, rotochemical or not. In the controlled benchmark that path is instantiated from the same
declared benchmark `Ltilde_a` used by the rotochemical extension, not from the historical placeholder
normalization. The default production normalization and its historical baselines remain unchanged;
source-authoritative replacement is a later realistic-physics task. Justification in §13.

---

## 13. No-double-counting architecture

### 13.1 The hazard

CompactStar already subtracts an equilibrium Urca luminosity. Adding a rotochemical
`L_nu(T, eta)` naively would count the equilibrium part twice.

### 13.2 The decision

**Incremental, with a shared coefficient.** Three requirements, all mandatory:

**(R1) Incremental form.** For each rotochemical channel add
`DeltaL_{nu,a} = L_{nu,a}(T,eta) - L_{nu,a}(T,0)` and `L_{H,a}`, never the full `L_nu(T,eta)`. At
`eta = 0`, `F_*(0) = 1` and `H_*(0) = 0` give `DeltaL_nu = 0` and `L_H = 0` **identically**, so the
extended channel reduces exactly to its own declared equilibrium luminosity. This is a same-
coefficient structural identity, not a requirement to match today's historical placeholder cooling
bit-for-bit.

**(R2) Shared coefficient.** The equilibrium term and the correction must be built from **the same**
`Ltilde_a`, structurally:

```text
NeutrinoCooling equilibrium   :  L_nu,eq  = sum_a Ltilde_a T_inf^{q_a}
Rotochemical correction       :  DeltaL_nu = sum_a Ltilde_a [F_*(xi_a) - 1] T_inf^{q_a}
Rotochemical heating          :  L_H       = sum_a Ltilde_a xi_a H_*(xi_a) T_inf^{q_a}
Chemical ODE rate             :  R_l       = sum_{a in l} [Ltilde_a/k_B^(erg)] H_*(xi_a) T_inf^{q_a - 1}
```

If a *second*, independently normalized coefficient set were introduced for the correction, the
enhancement `F_*` — which is physically a **ratio** to the star's own equilibrium emissivity — would
be applied to a different emissivity than the one being subtracted. `eta = 0` would still give
`DeltaL_nu = 0`, so the double-counting test would *pass* while the physics was wrong. **This is why
(R2) is mandatory and why a parallel duplicate module is forbidden.**

**(R3) Single conversion boundary.** The chemical module forms
`P_H,chem^inf=sum_l eta_l^inf R_l` in MeV/s and converts once to `L_H^inf` in erg/s at the thermal
RHS/luminosity boundary. Nonequilibrium neutrino luminosity remains erg/s. For the controlled
free-gas benchmark, equilibrium cooling must be constructed from the same declared benchmark
`Ltilde_a`; the historical `Q0`/`K` placeholders are neither source authority nor a validation
oracle.

### 13.3 Consequence for the existing code

`NeutrinoCoolingCachePayload` must be extended from two lumped coefficients to the channel-resolved
set `{Ltilde_De, Ltilde_Dmu, Ltilde_Me, Ltilde_Mmu}`, with the existing driver consuming the sums
`K_DU = Ltilde_De + Ltilde_Dmu` and `K_MU = Ltilde_Me + Ltilde_Mmu`. That change:

- is **required** by the chemical ODE independently (FR2005 eqs. 52–53 need separate `H_M(xi_npe)`
  and `H_M(xi_npmu)`);
- requires the controlled benchmark to instantiate its equilibrium Urca term from the **same
  declared benchmark `Ltilde_a`** used for `F_*`, `H_*`, rates, and heating. The existing historical
  placeholder `NeutrinoCooling` normalization is not the controlled benchmark's equilibrium
  coefficient;
- **changes the numerical behaviour of default production equilibrium cooling** only if the
  placeholder `Q0` values are replaced. That is a separate governed realistic-physics change with
  its own baseline consequences and is **not** authorized by ADR-0014. This clarification changes no
  historical baseline;
- Channels with no rotochemical treatment (PBF and any future non-beta channel) stay entirely inside
  the equilibrium driver and receive no correction.

---

## 14. Sign-crossing roots

Computed independently from the §9.2 polynomials (`scratchpad/calc/roots.py`, mpmath, 30 digits).

### 14.1 Root definitions

```text
INCREMENTAL root  :  smallest xi > 0 with  M_*(xi) + 1 = 0 ,  i.e.
                     L_H^inf = L_nu^inf(T,eta) - L_nu^inf(T,0)
                     — equilibrium neutrino emission IS subtracted.

FULL root         :  smallest xi > 0 with  M_*(xi) = 0 ,  i.e.
                     L_H^inf = L_nu^inf(T,eta)
                     — equilibrium neutrino emission is NOT subtracted.
```

### 14.2 Results

| Quantity | Direct Urca | Modified Urca |
|---|---|---|
| **Incremental root `xi`** | **`4.7870134733369`** | **`4.90971002892413`** |
| Incremental root, exact form | `xi = pi sqrt((sqrt(93) - 5)/2)` | `xi = pi sqrt(v)`, `v` the positive root of `3 v^3 + 84 v^2 + 378 v - 1468 = 0`, `v = 2.44237272219924` |
| Defining polynomial (`u = xi/pi`) | `21 u^4 + 105 u^2 - 357 = 0`, i.e. `u^4 + 5u^2 - 17 = 0` | `15 u^6 + 420 u^4 + 1890 u^2 - 7340 = 0` |
| **Full root `xi`** | **`5.4585315948676`** | **`5.63371746764834`** |
| Defining polynomial (`u = xi/pi`) | `21 u^6 + 105 u^4 - 357 u^2 - 457 = 0` | `15 u^8 + 420 u^6 + 1890 u^4 - 7340 u^2 - 11513 = 0` |
| Maximum incremental cooling at | `xi = 3.49729392477`, `M_D = -1.5277746` | `xi = 3.61254068813`, `M_M = -1.467659` |

Residuals `M_*(root)` and `M_*(root)+1` evaluate to `<= 1.6e-30` and `0.0` respectively at 30 digits.

### 14.3 Adjudication of the previous project estimate

The earlier rough project figures **4.79 (direct)** and **4.91 (modified)** are hereby identified as
the **INCREMENTAL** roots (`4.78701`, `4.90971`) — i.e. equilibrium neutrino subtraction *is*
included. They were **not** the FR2005 eq. (48) roots. The full (`M_* = 0`) roots are `5.459` and
`5.634`, consistent with FR2005 p. 11: *"completely balancing the cooling for `|xi| ~ 5.5`"*.
FR2005's *"reaching a maximum cooling at `|xi| ~ 3.5`"* likewise matches the computed `3.497`/`3.613`.
**Both project numbers were right; their meaning was unlabelled.** ADR-0014 makes the label
mandatory.

### 14.4 Test-oracle status

All four roots are **exact closed-form or exact-polynomial** quantities depending on nothing but the
imbalance functions. They are durable oracles (RE6). The `4.787` root is fully closed-form.

---

## 15. Direct / modified Urca channel support

### 15.1 Current direct-Urca threshold ownership

`StarContext::BuildDirectUrcaMaskCache_` (`CompactStar/Physics/Evolution/src/StarContext.cpp:475-600`):

| Aspect | Finding |
|---|---|
| Criterion | `kF_n <= kF_p + kF_e` with `kF = (3 pi^2 n)^{1/3}`, `n_i = Y_i n_B` |
| Composition source | profile species columns `"10"` (n), `"11"` (p), `"0"` (e) — **composition-source authoritative** (INV-01), not a hardcoded fraction |
| Lepton coverage | **electron only.** No `kF_n <= kF_p + kF_mu` test exists anywhere in the tree. |
| Domain guards | `nB_min = 1e-6 fm^-3`, `n_min = 1e-12 fm^-3`, documented in-source as *"numerical/semantic guards, not DU thresholds"*; they exist to stop the degenerate `0 <= 0 + kF_e` false positive in an electron-only crust |
| Boundary semantics | A last-index representation is consumed as `[0,last]`. The builder appears to assume the first allowed region begins at the centre; a future physical outer allowed shell could therefore sweep a closed inner region unless support is represented explicitly. |
| Consumers | `NeutrinoCooling_Details.cpp:83,188-197` — the DU coefficient integral runs over `[0, durca_last]` only |
| Cache | invalidated with the profile version (ADR-0003) |

**Assessment.** The triangle criterion is source-consistent (R1997 p. 4: *"`x > 1/9` for direct
Urca reactions"* is the npe form), but triangle support and applicability under the degenerate-Urca
model are distinct. `nB_min` is a numerical/semantic guard, **not** the DU threshold. A future
adapter must represent support explicitly, add the muon criterion with `kF_mu`, and assert the
profile is ordered innermost first before using any ordering-dependent integration.

### 15.2 Free-gas fixture: is direct Urca open?

Computed for the governed FR2005 free-gas EOS (`scratchpad/calc/freegas_du.py`), solving npe-mu beta
equilibrium with `m_n, m_p, m_e, m_mu, hbar c` taken from the Phase-5C provenance record:

| `n_B` [fm^-3] | `rho` [g/cm^3] | `x_p` | `kF_n` | `kF_p + kF_e` | DU(e)? | `kF_p + kF_mu` | DU(mu)? | muons? |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 1.72e14 | 0.00327 | 283.0 | 84.1 | no | 42.1 | no | no |
| 0.30 | 5.30e14 | 0.00816 | 407.6 | 164.6 | no | 82.3 | no | no |
| **0.6049** | **1.100e15** | **0.01566** | **513.6** | **252.8** | **no** | **193.4** | **no** | **yes** |
| 1.00 | 1.87e15 | 0.02663 | 605.0 | 344.0 | no | 304.7 | no | yes |
| 3.00 | 6.16e15 | 0.05680 | 863.5 | 617.1 | no | 596.2 | no | yes |

Muon DU triangle support does not open where muons are present. Electron DU triangle support is
closed through the dense stellar interior by a wide margin at the centre: `x_p <= 0.0157` against
`1/9 = 0.111`. A very-low-density electron-DU kinematic sliver nevertheless exists near neutron
onset, approximately `n_B ~ 7.36e-9` to `6.67e-8 fm^-3`, with neutron Fermi temperature near the
upper edge `T_F,n ~ 3.5e7 K`. The sliver is outside the intended source-validity regime used to
motivate the controlled benchmark at early/high-temperature epochs, but its degeneracy status
changes as the benchmark cools and is **not** a static exclusion rule. **Muons are present** in the
interior (`mu_e = 123.6 MeV` at the centre against `m_mu = 105.66 MeV`; onset near
`n_B ~ 0.45 fm^-3`).

Central state `n_B = 0.604925 fm^-3` at `rho_c = 1.10e15 g/cm^3`. A scoping TOV solve on this EOS
gives `M = 0.62362 M_sun`, agreeing with the governed Structure-1 `M_cut ~= 0.623635569 M_sun`
(`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1.md:371`) to five digits, and with FR2005
Table 1's free-gas `M_max = 0.62 M_sun`.

**Controlled-benchmark consequence.** The free-gas Phase-5D architecture benchmark is statically
configured **MODIFIED-URCA-ONLY**: its enabled-process set is `{Me, Mmu}` and excludes `{De, Dmu}`.
Therefore `D_De = empty` and `D_Dmu = empty` for this benchmark by the declared benchmark
process-configuration contract, not because the triangle condition never fires, not because of a
temperature-dependent degeneracy predicate, and not because of `nB_min`. Triangle support, declared
benchmark process support, and general future physical DU applicability are distinct. The
low-density triangle-open sliver remains a negative control: a future implementation must prove
that it does not activate DU when that process is disabled or outside the declared channel-support
contract, and that an outer allowed shell cannot make `[0,last]` integrate a closed inner region.

Future realistic or source-authoritative evolution must obtain DU support from a governed
applicability provider capable of representing kinematics, composition, degeneracy/model-validity
requirements, and possibly disconnected support intervals. The controlled benchmark's disabled-DU
configuration is not a universal physical rule.

### 15.3 Moving thresholds in secular evolution

For **frozen** v1 coefficients (§19) the background model and the benchmark's declared process set
are fixed. No temperature-dependent support predicate is introduced and there is no moving-threshold
term. Two future-only caveats, recorded so they cannot be introduced silently:

1. If the background is ever allowed to evolve (`Z(t)`), the DU boundary radius becomes `r_DU(t)`
   and `d Ltilde_D/dt` acquires a boundary term. Out of scope.
2. `eta != 0` shifts the *local* chemical potentials but **not** the DU kinematic threshold, which is
   a Fermi-momentum condition on the equilibrium background. Future physical applicability still
   requires the governed provider specified in §15.2; the current benchmark configuration is not a
   universal DU rule.

---

## 16. Modified-Urca branch accounting

### 16.1 What FR2005 requires

FR2005 §3.2, p. 11: *"We denote modified Urca reactions with electrons and muons by `alpha = Me` and
`alpha = Mmu` respectively, in each case adding the contributions of the neutron and proton branches
(Yakovlev et al. 2001). Direct Urca with electrons or muons is denoted by `alpha = De` and `Dmu`."*

So FR2005's channel set is **four** channels (`De, Dmu, Me, Mmu`), each already **summed over the
neutron and proton branches**. The branches are an *internal decomposition of `Ltilde_a`*, not
separate evolution channels: `Ltilde_Me = Ltilde_Me^{(n branch)} + Ltilde_Me^{(p branch)}`.

### 16.2 Minimum honest v1 scope

| Requirement | Verdict |
|---|---|
| Electron/muon channels separate | **Mandatory.** FR2005 eqs. (52)–(53) evaluate `H_M` at `xi_npe` and at `xi_npmu` separately. Collapsing them is mutation M5/M6 and destroys the two-channel physics. |
| Neutron/proton branches separate as *evolution* channels | **Not required.** Both branches share the same `eta_l` and the same `xi_l`; they differ only in `S_a(n)` and in the proton-branch kinematic threshold `p_{F,l} + 3 p_{F,p} > p_{F,n}`. Summing them inside `Ltilde_{M,l}` is exactly what FR2005 does. |
| Direct Urca channels | **Present in the contract, dormant in the free-gas v1** because the benchmark's static enabled-process set excludes `De` and `Dmu`, so their declared benchmark domains are empty (§15.2). The slots and their `q = 6` exponent must exist, but the free-gas star exercises neither. |

**Minimum v1: four declared channels `{De, Dmu, Me, Mmu}`; two live (`Me`, `Mmu`) for the free-gas
fixture; branch summation internal to each `Ltilde`.**

### 16.3 Current code

The current equilibrium cooling distinguishes **neither** the lepton nor the branch. It must gain the
lepton split (mandatory) and should gain the branch split inside `Ltilde_{M,l}` when a
branch-resolved `S_a` is supplied.

---

## 17. Heat-capacity coupling and photon cooling

### 17.1 Heat capacity — unchanged

ADR-0002 fixes exactly one denominator, `C_*(T_inf)`, and already writes its governing invariant as

```text
C_*(T_inf) dT_inf/dt  =  - L_nu,inf  -  L_gamma,inf  +  L_H,inf  +  ...
```

with the explicit note *"Terms may be added to the right-hand side — chemical heating `L_H,∞`, BNV
heating, and others — without altering this ADR."*
**Phase-5D requires no change to ADR-0002 and none is proposed.** The FR2005 definition
`C = int dV c_V(T_local)` (eqs. 3, 50) is identical to `C_*(T_inf)` including the `T_local = T_inf e^{-nu}`
substitution and the proper-volume measure.

### 17.2 Thermal ODE in the variable the solver actually evolves

The evolved variable is `x = ln(T_inf / T_ref)`, `T_ref = 1e8 K` (`ThermalState.hpp:150-200`).
Therefore the rotochemical extension contributes

```text
dx/dt  =  [ - L_nu,eq^inf(T_inf)  -  L_gamma^inf(T_inf)
            + sum_a Ltilde_a ( xi_a H_*(xi_a) - F_*(xi_a) + 1 ) T_inf^{q_a} ]
          / ( T_inf * C_*(T_inf) )
```

where the first two terms are **exactly** what the existing `NeutrinoCooling` and `PhotonCooling`
drivers already contribute, and the bracketed sum is the new rotochemical contribution
`sum_a DeltaP_{beta,a}` of §12.2. Equivalently, isolating the new driver's own contribution:

```text
(dx/dt)_rotochem  =  [ L_H^inf - DeltaL_nu^inf ] / ( T_inf C_*(T_inf) )
                  =  sum_a Ltilde_a [ M_*(xi_a) + 1 ] T_inf^{q_a} / ( T_inf C_*(T_inf) ) .
```

Consistency with FR2005 eq. (51): substituting `L_nu,eq = sum_a Ltilde_a T_inf^{q_a}` and
`C_* = Ctilde T_inf` recovers
`Tdot_inf = (1/Ctilde)[ sum_a M_*(xi_a) Ltilde_a T_inf^{q_a - 1} - L_gamma^inf / T_inf ]`, which is
FR2005 eq. (51) term for term (their `T_inf^5` and `T_inf^7` are `T_inf^{q-1}` for `q = 6, 8`).
**The incremental and full formulations give the identical total RHS**; they differ only in which
driver owns which term. Unit check: `[erg/s]/([K][erg/K]) = 1/s`. Correct for `dx/dt`.

### 17.3 Photon cooling — unchanged (v1 contract)

`PhotonCooling` declares `DependsOn = {Thermal}` only (`PhotonCooling.hpp:287-292`) and takes
`T_inf -> T_s` through the envelope. **No chemical-imbalance correction to the envelope is
introduced.** Authority: FR2005 §3.3 applies the Potekhin et al. (1997) envelope as a pure function
of `T_inf` and surface gravity, with no `eta` dependence anywhere; chemical/reaction quantities are
integrated over the declared reaction/chemical domain `D` and channel subsets `D_a`, while the
surface envelope is not part of those integrals. FR2005 §3.5's realistic model uses its own governed
core domain, but that is not an implicit domain rule for the controlled free-gas fixture. Recorded
as a v1 assumption: *photon cooling is affected by rotochemical heating only through `T_inf`.* Any
future envelope `eta` dependence requires its own source authority and its own ADR.

---

## 18. Spin ownership

### 18.1 Requirement

The rotochemical module must not hard-wire a spin-down law. It needs `Omega(t)` and `Omegadot(t)`,
nothing more.

### 18.2 Audit of the three available routes

| Route | Mechanism | Verdict |
|---|---|---|
| (a) Recompute the torque inside the rotochemical driver | read `SpinState::Omega()` and apply a braking law | **Rejected.** Duplicates `MagneticDipole`, silently diverges from it, and hard-wires one law. |
| (b) Read the spin RHS out of the accumulator | `RHSAccumulator::Peek(StateTag::Spin, 0)` after the spin driver has run | **Rejected as the contract.** `EvolutionSystem::operator()` (`EvolutionSystem.cpp:119-128`) iterates drivers in **registration order** with no topological sort, even though `IDriver::DependsOn()` exists. Correctness would silently depend on registration order — a defect the accumulator's `Clear()`-then-accumulate design does not protect against. |
| (c) **An explicit spin-history interface supplied through `DriverContext`** | `ISpinHistory { double Omega(t, Y) const; double OmegaDot(t, Y) const; }` | **Recommended.** |

### 18.3 Recommendation

A small `ISpinHistory` (or equivalently named) read-only interface, obtained from the
`DriverContext`, with two implementations in the first cut:

- a **prescribed** history (analytic `Omega(t)`, `Omegadot(t)`, or tabulated observed pulsar timing) —
  this is what the deterministic v1 benchmark uses, and it makes the benchmark reproducible without
  any spin DOF at all;
- a **state-coupled** adapter that evaluates the *same* torque object the spin driver uses, so both
  cannot diverge.

`RotochemicalSpinDrive::Evaluate(omega, omega_dot)` already exists
(`ChemicalResponse.cpp:824-835`) and returns `2 W Omega Omegadot` per channel — the rotochemical
module simply feeds it whatever the spin history reports. This keeps the machinery reusable for
observed pulsar timing, magnetic-dipole braking, gravitational-wave spin-down, and future exotic
torques, and it keeps `MagneticDipole` an *example driver*, not a physics dependency.

**If option (b) is ever chosen instead, a driver dependency sort becomes mandatory, not optional.**

---

## 19. Frozen vs time-dependent coefficients

### 19.1 What the sources assume

FR2005 §2.1, p. 4, is explicit: *"Since the `B_ij` do not depend on time, we invert and take the time
derivative of equation (11)"*. The `B_ij` (and hence `Z`) are evaluated on the **nonrotating
background star** (FR2005 eq. 12 integrates over the unperturbed core; R2006 eq. (6) says *"the
volume of the unperturbed, nonrotating stellar core"*). `I_{Omega,i}` is likewise *"a constant,
depending on the structure of the unperturbed, nonrotating star"* (R2006 §3, p. 569), and FR2005
eq. (19) treats `(dN/dOmega^2)_{P,rho_c}` as independent of `Omega^2`.

**Frozen `Z` and `W` evaluated on the nonrotating background is what FR2005 and R2006 actually do.**
It is the correct v1 contract, not an approximation invented here.

### 19.2 The extra term if `Z` were ever allowed to vary

From `eta^inf = -Z delta N_l`, without the frozen assumption:

```text
etadot^inf = - Zdot delta N_l  -  Z d(delta N_l)/dt
           = + Zdot Z^-1 eta^inf  -  Z R  +  2 W Omega Omegadot
```

i.e. an **extra term `+ Zdot Z^-1 eta^inf`**, and `W = Z I_Omega` would additionally acquire the
time dependence of both factors. Note the sign: it is `+Zdot Z^-1 eta`, not `-`, because
`delta N_l = -Z^-1 eta^inf`.

**Not implemented, not authorized.** Recorded here specifically so that a future agent cannot enable
time-dependent coefficients without noticing that a term is missing. Silently making `Z` time
dependent while keeping the v1 RHS is mutation M18 (§32).

### 19.3 Provenance consequence

Because `Z`, `W`, `I_Omega` and `Ltilde` are all frozen, they are **inputs** to an evolution run and
must be captured in the run provenance with their Phase-5C identity, exactly as ADR-0013 §6 requires.
A global reaction coefficient must additionally retain its declared chemical/reaction domain `D`
and channel support/applicability domain `D_a`. A run whose coefficients cannot be traced to a
specific `ChemicalImbalanceResponse` / `RotochemicalSpinDrive` result must refuse.

---

## 20. Initial conditions

### 20.1 What FR2005 actually does

FR2005 Fig. 4 caption, p. 29: *"with initial conditions `T_inf = 10^8 K` and null chemical imbalances
at `t = 0`"*. R2006 Fig. 2 caption, p. 571: *"with initial conditions of vanishing chemical
imbalances and `T^inf = 10^8 K` at `t = 0`"*, `B = 10^8 G`, `P_0 = 1 ms`. FR2005 Fig. 5 explores
`eta^inf(0) = 0, kT_inf, 10 kT_inf, 20 kT_inf` and `T_inf(0)` over `10^4`–`10^9 K`, and concludes
(p. 15) *"the value of the temperature at the quasi-equilibrium state and the time required to arrive
do not depend on initial conditions"*.

### 20.2 Contract

| Item | v1 value | Status |
|---|---|---|
| `eta_npe^inf(0)`, `eta_npmu^inf(0)` | `0` | **Conventional, not required.** FR2005 Fig. 5 (right) demonstrates the quasi-steady state is reached from non-zero initial imbalances too. Physically motivated: a young, hot star equilibrates fast (Y2020 §4.1: *"For a very young NS, the Urca reaction is fast enough so that the chemical equilibrium is maintained"*). |
| `T_inf(0)` | caller input; `1e8 K` reproduces FR2005/R2006 | **Caller input, never a physics constant.** |
| `Omega(0)` / `P_0` | caller input; `P_0 = 1 ms` reproduces FR2005 Fig. 4 / R2006 Fig. 2 | **Caller input.** |
| `B` (if a dipole driver is used) | caller input; `1e8 G` reproduces FR2005 Fig. 4 | **Caller input, and owned by the spin driver, not the rotochemical module.** |

Initial chemical equilibrium is therefore **conventional and default, not mandatory**; the contract
must accept arbitrary `eta(0)` so that FR2005 Fig. 5 (right) can be reproduced as a test.

---

## 21. Analytic limits

All nine limits from the task statement, with their status. `Ltilde_a > 0`, `Z` symmetric positive
definite, `W = Z I_Omega`.

| # | Limit | Statement | Status |
|---|---|---|---|
| **A** | `eta = 0` | `xi = 0`, `H_*(0) = 0` ⇒ `R_l = 0` exactly ⇒ `DeltaGamma = 0`. | Exact; structural zero |
| **B** | no spin-down (`Omegadot = 0`) | `etadot = -Z R`. With frozen symmetric SPD `Z`, `dV/dt = -2 sum_l eta_l R_l <= 0`. Strict decrease requires an active positive-normalization channel with nonzero dissipative response for every nonzero direction considered; a direction with all coupled normalizations zero may freeze (§7.3). | **Proved nonincrease**; strict only under all-active qualification |
| **C** | no reactions (`Ltilde_a = 0`) | `etadot^inf = 2 W Omega Omegadot` exactly; integrating with `W` frozen, `eta_l^inf(t) = W_l [Omega^2(t) - Omega_0^2]`, i.e. `|W_l|(Omega_0^2 - Omega^2)` for spin-down — **FR2005 eq. (77)**. | Exact; closed form |
| **D** | small `xi` | `R_l -> [1/k_B^(erg)] [Ltilde_{D,l} 0.158300492605 + Ltilde_{M,l} 0.129192649689 T_inf^2] xi_l T_inf^5`, linear in `eta`; substituting `xi=eta/[k_B^(MeV)T]` gives the governed mixed-unit linear matrix rate. | Exact |
| **E** | large `\|xi\|` | `R_l -> [1/k_B^(erg)][Ltilde_{D,l} 42 xi^5/(457 pi^6) T^5 + Ltilde_{M,l} 24 xi^7/(11513 pi^8) T^7]`; `DeltaP_beta -> Ltilde_a M_*(xi) T^q`. | Exact power laws |
| **F** | heating positivity | `C_(MeV->erg) eta_a[MeV] DeltaGamma_a = Q_a^eq[erg...] xi_a H_*(xi_a) >= 0` pointwise, both signs of `eta`. ⇒ `L_H^inf >= 0`. | **Proved**; verified at `xi = ±0.3, ±1, ±5, ±50` |
| **G** | equilibrium `DeltaL_nu` | `F_*(0) = 1` ⇒ `DeltaL_nu = Ltilde_a [F_*(0) - 1] T^q = 0` exactly. | Exact; structural zero |
| **H** | zero imbalance, same coefficient | For a channel whose equilibrium and rotochemical coefficient are the same declared `Ltilde`, `H(0)=0`, `F(0)=1`, `R=DeltaGamma=0`, `DeltaL_nu=0`, and `L_H=0`; it reduces exactly to its own `L_nu,eq=Ltilde T^q`. | Exact same-coefficient identity; historical placeholder is not the oracle |
| **I** | zero `W` | `2 W Omega Omegadot = 0` ⇒ no spin-generated disequilibrium; combined with (B), `eta(0) = 0` implies `eta(t) = 0` for all `t`. | Exact |

---

## 22. Quasi-steady state

### 22.1 Derivation

Set `etadot^inf = 0`, `Tdot_inf = 0`. In the pure-modified-Urca, `xi >> 1` regime (FR2005 eqs.
59–60, error `< 1%` at `xi ~ 200`):

```text
H_M(xi) ~= C_H xi^7 ,  C_H = 24/(11513 pi^8)
M_M(xi) ~= C_M xi^8 ,  C_M = 15/(11513 pi^8)
```

`etadot = 0` gives (FR2005 eqs. 62–63)

```text
Z_npe Ltilde_Me (eta_npe^inf)^7 + Z_np   Ltilde_Mmu (eta_npmu^inf)^7
  = [2 k_B^(erg) (k_B^(MeV))^7/C_H] W_npe Omega Omegadot
Z_np Ltilde_Me (eta_npe^inf)^7 + Z_npmu Ltilde_Mmu (eta_npmu^inf)^7
  = [2 k_B^(erg) (k_B^(MeV))^7/C_H] W_npmu Omega Omegadot
```

whose solution — using `W = Z I_Omega`, so that `Z` cancels exactly — is FR2005 eqs. (64)–(65):

```text
eta_npe^inf  = k_B^(MeV) [ 2 k_B^(erg) I_{Omega,e}  / (C_H Ltilde_Me ) ]^{1/7} (Omega Omegadot)^{1/7}
eta_npmu^inf = k_B^(MeV) [ 2 k_B^(erg) I_{Omega,mu} / (C_H Ltilde_Mmu) ]^{1/7} (Omega Omegadot)^{1/7}
```

Both are positive because `I_{Omega,l} < 0` and `Omega Omegadot < 0`. **The cancellation of `Z` is
exactly the `W = Z I_Omega` identity of §6.2** — the quasi-steady imbalances depend on `I_Omega`, not
on `Z`. This is why R2006 §3, p. 570, can state that *"the quasi–steady state reached by NSs with
slow spin-down evolution ... is strictly unaffected by these corrections"*: the electrostatic
correction changes `Z`, and `Z` cancels.

Then `Tdot = 0` with FR2005 eq. (61) gives eq. (66):

```text
L_gamma,eq^inf = C_M (2 k_B^(erg) / C_H)^{8/7} [ (I_{Omega,e}^8 /Ltilde_Me )^{1/7}
                                       + (I_{Omega,mu}^8/Ltilde_Mmu)^{1/7} ] |Omega Omegadot|^{8/7}
```

**independent of envelope model**, because at `xi >> 1` reactions are imbalance-dominated and the
`T_inf` dependence drops out (FR2005 p. 15).

### 22.2 Source scalings, in our convention

| Quantity | FR2005 | Our variables |
|---|---|---|
| `eta_qs` | eq. (64)/(65) | `∝ (I_{Omega,l} / Ltilde_{M,l})^{1/7} \|Omega Omegadot\|^{1/7}` — depends on `I_Omega`, **not** on `Z` |
| `L_gamma,qs` | eq. (66) ⇒ (67) | `~= 10^{30-31} (Pdot_{-20}/P_ms^3)^{8/7} erg/s` |
| `T_s,qs^inf` | eq. (68) | `~= (2-3) x 10^5 (Pdot_{-20}/P_ms^3)^{2/7} K` |
| `L_gamma,qs / Edot` | eq. (69) | `~ (0.3-3) x 10^{-5} (Pdot_{-20}/P_ms^3)^{1/7}` |
| `xi_qs` | eq. (70) | `∝ \|Omega Omegadot\|^{(alpha-8)/(7 alpha)}`, `alpha = 2.42` (Potekhin envelope exponent) |
| `tau_eq` | eq. (81) | `~ 1.6 x 10^7 (P_ms^3/Pdot_{-20})^{6/7} yr` |
| `A = tau_eq/tau_sd` | eqs. (79),(83) | Expressed with `k_B^(erg)(k_B^(MeV))^7` under the seventh root, consistent with the mixed-unit rate convention; no bare `k^8` implementation formula is permitted. |

### 22.3 Numerical validation of the quasi-steady map (scoping)

`scratchpad/calc/evolve.py` integrates the full system (§25) with the **governed** free-gas `Z`, `W`
and a scoping `Ltilde`, then compares the endpoint to FR2005 eqs. (64)–(65):

```text
t = 1e10 yr, B = 1e8 G, P_0 = 1 ms
integrated : eta_npe = 2.2467e-02 MeV   eta_npmu = 2.5750e-02 MeV   T_inf = 1.5573e6 K
FR2005 (64)/(65) : eta_npe = 2.2500e-02 MeV   eta_npmu = 2.5782e-02 MeV
agreement  : 0.15 %  and  0.12 %
xi at quasi-steady : xi_npe = 167.4 , xi_npmu = 191.9      (FR2005 Fig. 4 quotes xi ~ 200)
```

**This is a high-value durable validation target (RE14): the asymptote is a closed-form function of
`I_Omega` and `Ltilde` alone.**

---

## 23. ODE stiffness and solver strategy

### 23.1 Measured

Full three-variable system (`ln T_inf`, `eta_npe`, `eta_npmu`), free-gas star,
`t = 0` to `1e10 yr`, `B = 1e8 G`, `P_0 = 1 ms` (`scratchpad/calc/stiff.py`):

| Method | rtol | atol | steps | RHS evals | Jacobians | wall (s) | endpoint `T_inf` | endpoint `eta_npe` |
|---|---|---|---|---|---|---|---|---|
| RK45 (explicit) | 1e-6 | 1e-10 | 4829 | 33260 | 0 | 0.71 | 1.5573e6 | 2.2467e-2 |
| RK45 (explicit) | 1e-9 | 1e-12/1e-24 | 6420 | 38522 | 0 | 0.82 | 1.5573e6 | 2.2467e-2 |
| DOP853 (explicit) | 1e-9 | 1e-12/1e-24 | 2778 | 48221 | 0 | 1.00 | 1.5573e6 | 2.2467e-2 |
| **LSODA (auto-switching)** | 1e-9 | 1e-12/1e-24 | **561** | **1105** | 33 | 0.02 | 1.5573e6 | 2.2467e-2 |
| BDF | 1e-9 | 1e-12/1e-24 | 686 | 2029 | 43 | 0.09 | 1.5573e6 | 2.2467e-2 |
| Radau | 1e-9 | 1e-12/1e-24 | 644 | 5325 | 104 | 0.16 | 1.5573e6 | 2.2467e-2 |

All methods agree on the endpoint to five digits. LSODA's non-zero Jacobian count shows it switched
to its stiff branch.

### 23.2 Finding

**The system is moderately stiff but comfortably tractable with the existing explicit RKF45.** An
explicit method costs roughly **20–30× more RHS evaluations** than a stiff method, but succeeds. The
repository default `max_internal_steps = 10000` per output interval is above the ~5–6 × 10^3 total
steps required.

Two caveats, recorded:

1. This measurement is for **modified Urca only** — the free-gas case. With direct Urca open
   (`q = 6`, much larger `Ltilde`, `xi^5` rather than `xi^7` in the rate) the relaxation timescale
   drops by orders of magnitude and stiffness will be far worse. FR2005 Fig. 6 shows exactly that
   regime (a "metastable" quasi-equilibrium for the electron channel while the muon channel is still
   growing) — a genuinely two-timescale problem.
2. `EvolutionConfig` (`EvolutionConfig.hpp:170-190`) records that `MSBDF` is *rejected* because GSL
   requires a Jacobian CompactStar does not supply, and its own comment already anticipates this:
   *"A future rotochemical/chemical system may need a stiff method and a real Jacobian; that is NOT
   decided here."*

### 23.3 Recommendations (no solver change in this task)

| Item | Recommendation |
|---|---|
| State variables | Evolve `eta^inf` in **MeV**, not `xi`. `xi = eta/(k_B^(MeV) T)` couples the chemical and thermal errors, has no fixed scale, and becomes singular as `T -> 0`; `eta` is bounded, monotone in the spin-driven phase, and is the quantity with a source benchmark. |
| Absolute tolerances | The current scalar `atol = 1e-10` is applied to **every** component of a heterogeneous vector (`ln T` dimensionless, `Omega ~ 10^3 rad/s`, `eta ~ 10^-3..10^-2 MeV`). A **per-component `atol`** is recommended: `~1e-12` for `ln T`, `~1e-18 MeV` for `eta` (twelve decades below the quasi-steady scale). Not required for v1 correctness — the free-gas run converges with the scalar default — but it removes an obvious accuracy trap. |
| Relative tolerance | `rtol = 1e-6` reproduces the endpoint to five digits; `1e-9` for validation runs. |
| Stepper | Keep `RKF45` for v1. Revisit **before** enabling direct Urca or superfluidity. |
| Event handling | The two interesting events (`T_inf` minimum; `DeltaP_beta` sign change at `xi = 4.787`/`4.910`) are **diagnostics**, not integration events. No root-finding in the integrator is required. |
| Positivity | `eta` legitimately changes sign (a spun-**up** star), so no positivity clamp may be applied to `eta`. `ln T_inf` already guarantees `T_inf > 0`. |

---

## 24. Energy conservation and thermodynamic sanity

### 24.1 The ledger

Per unit time, measured at infinity, for a star with frozen background and fixed baryon number:

```text
  spin-down power             Edot_spin = -I_star Omega Omegadot          (> 0, supplied by the torque)
      |
      +--> stored as chemical free energy at rate     Pdot_chem, driven by 2 W Omega Omegadot
      |
  chemical free energy  U_chem
      |
      +--> released by beta reactions as              P_H,chem^inf = sum_a eta_a^inf R_a [MeV/s] >= 0
           and converted once at thermal boundary     L_H^inf = C_(MeV->erg) P_H,chem^inf [erg/s]
                |
                +--> escapes as extra neutrinos       DeltaL_nu^inf = sum_a Ltilde_a [F_*(xi_a)-1] T^q  >= 0
                |
                +--> deposited as heat                DeltaP_beta = L_H^inf - DeltaL_nu^inf   (either sign)
                          |
                          +--> raises T_inf, then escapes as photons (L_gamma^inf) and
                               as equilibrium neutrinos (L_nu,eq^inf)
```

### 24.2 Signs of every term

| Term | Sign | Reason |
|---|---|---|
| `L_H^inf` | **always `>= 0`** | `C_(MeV->erg) eta_a[MeV] R_a = Ltilde_a xi_a H_*(xi_a)T^q >= 0` (§7.4). The single positive unit conversion preserves the dissipative sign. |
| `DeltaL_nu^inf` | **always `>= 0`** | `F_*(xi) >= F_*(0) = 1`; all `F_*` coefficients positive. FR2005 p. 10: *"A finite `eta_a` of either sign enhances neutrino emission due to the even nature of the functions `F_*`."* |
| `L_nu,eq^inf` | `> 0` | ordinary cooling |
| `L_gamma^inf` | `> 0` | ordinary cooling |
| **`DeltaP_beta`** | **either sign** | `= Ltilde_a [M_*(xi_a) + 1] T^q`; negative for `0 < xi < 4.787` (D) / `4.910` (M), positive beyond |
| `2 W Omega Omegadot` | `> 0` for spin-down | §6.5 |

### 24.3 Why "chemical heating is always dissipative" and "the incremental beta effect can be negative" are both true

They answer different questions and there is no contradiction:

- **`L_H` is the energy actually released from the chemical reservoir into the matter.** It is a
  dissipation rate, `sum_a eta_a DeltaGamma_a`, a product of a generalized force with its conjugate
  flux, and the second law forces it non-negative. It cannot be negative at any `xi`.
- **`DeltaP_beta` is the *net incremental* effect on the thermal budget** of running the same star
  out of equilibrium rather than in it. Non-zero `eta` does two things at once: it releases `L_H`
  *and* it opens extra neutrino phase space (`DeltaL_nu`). At modest imbalance the second dominates —
  the star loses more energy to enhanced neutrino emission than it gains from chemical heating — so
  the *net* effect on `T_inf` is **extra cooling**. Only past `xi ~ 4.8`–`4.9` does the heating win.

Quantitatively at the maximum: `M_* + 1` reaches `-0.528` (direct, at `xi = 3.497`) and `-0.468`
(modified, at `xi = 3.613`), i.e. the incremental effect is at worst about half an extra equilibrium
cooling luminosity. In the large-`xi` limit a fixed fraction of the released energy escapes as
neutrinos — `1/2` (direct), `3/8` (modified) — and the rest stays behind, which is FR2005 p. 11.

**This ledger is the control case for a future BNV analysis.** Any future claim of BNV heating must
be stated against the same four-way split (chemical release / neutrino escape / photon escape /
thermal deposition) and must say explicitly which of the two differences it means.

---

## 25. Free-gas end-to-end benchmark design

### 25.1 Purpose

Validate **architecture, signs, energy bookkeeping and ODE correctness** end to end. **This is not a
source reproduction**, and no claim about realistic nuclear matter follows from it. Because the
benchmark injects declared `Ltilde` values, it validates state evolution, reaction/heating algebra,
thermal coupling, and spin coupling, but does **not** by itself validate the stellar construction of
`Ltilde`; future gate RE10b validates that separate `GlobalUrcaChannelCoefficient` layer.

### 25.2 Fixture

| Item | Value | Provenance |
|---|---|---|
| EOS | FR2005 noninteracting `npe-mu` Fermi gas, whole star | `track-r-fernandez-reisenegger-2005-free-gas-local` provider |
| Structure | Track-R Structure-1 midpoint, `rho_c = 1.10 x 10^15 g/cm^3` | `docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1.md` |
| `M`, `R` | `0.6236 M_sun`, `~12.77 km` | Structure-1; matches FR2005 Table 1 free-gas row |
| `Z` | `[[4.5793031807026964e-54, 5.1725199102788050e-55], [5.1725199102788054e-55, 1.0268727975139168e-52]]` MeV/count | `phase5c_chemical_coefficients_candidate.json` (**candidate, not yet canonically integrated**) |
| `W` | `[-5.4061775017047240e-07, -1.5845719480103649e-06]` MeV s^2 | same |
| `I_Omega` | `[-1.1637998545112904e+47, -1.4844819850233820e+46]` count s^2 | same |
| Active Urca channels | Static enabled-process set `{Me, Mmu}`; `{De, Dmu}` disabled, so `D_De = D_Dmu = empty` for this controlled benchmark | §15.2 |
| Muon support | present for `n_B >~ 0.45 fm^-3` | §15.2 |
| Heat capacity | `StarContext::HeatCapacityStar_Tinf` (ADR-0002) | existing |
| Envelope | existing Potekhin model | existing |

### 25.3 Run specification

| Item | Value |
|---|---|
| `T_inf(0)` | `1e8 K` (FR2005 Fig. 4 / R2006 Fig. 2 convention) |
| `eta_npe^inf(0) = eta_npmu^inf(0)` | `0` |
| Spin history | prescribed magnetic-dipole braking, `B = 1e8 G`, `P_0 = 1 ms`, `P(t) = sqrt(P_0^2 + 2 (P Pdot) t)`, `P Pdot = (B/3.2e19)^2` — supplied through `ISpinHistory`, **not** owned by the rotochemical module |
| Duration | `0` to `1e10 yr` (quasi-steady is reached and tracked well before the end) |
| Reaction normalization | Positive **declared mathematical/architecture benchmark coefficients**, provenance-recorded and explicitly not realistic source normalizations; the same values feed equilibrium cooling, `F`, `H`, and chemical heating (§10.4, §27) |
| Process configuration | Static MODIFIED-URCA-ONLY; no temperature-dependent support predicate is introduced |

### 25.4 Expected qualitative behaviour (from the scoping integration)

| `t` [yr] | `T_inf` [K] | `eta_npe` [MeV] | `eta_npmu` [MeV] | `xi_npe` | `xi_npmu` |
|---|---|---|---|---|---|
| 1e2 | 9.997e7 | 1.31e-6 | 3.85e-6 | 0.00 | 0.00 |
| 1e4 | 9.716e7 | 1.20e-4 | 3.47e-4 | 0.01 | 0.04 |
| 1e6 | 1.719e7 | 1.21e-2 | 3.34e-2 | 8.18 | 22.54 |
| 1e7 | 9.739e6 | 3.85e-2 | 4.44e-2 | 45.9 | 52.9 |
| 1e8 | 9.272e6 | 3.80e-2 | 4.38e-2 | 47.6 | 54.8 |
| 1e10 | 1.557e6 | 2.25e-2 | 2.58e-2 | 167.4 | 191.9 |

Phases: (i) `eta` grows from zero driven purely by `2 W Omega Omegadot` while `T_inf` cools
conventionally; (ii) `xi` crosses the incremental sign-crossing roots — measured at
`xi_npe = 4.906` (`t ~ 8.0e5 yr`) and `xi_npmu = 4.86` (`t ~ 4.8e5 yr`), against the predicted
`4.90971` (grid-resolution limited); (iii) the star settles into quasi-steady, thereafter tracking
`|Omega Omegadot|`.

### 25.5 Test observables

| ID | Observable | Oracle |
|---|---|---|
| B1 | `eta(t)` for `Ltilde = 0` | `eta_l(t) = W_l [Omega^2(t) - Omega_0^2]`, exact (§21 C) |
| B2 | `eta(t)` for `Omegadot = 0`, `eta(0) != 0` | `eta^T Z^-1 eta` nonincreasing; strictly decreasing only for all-active dissipative directions; an uncoupled dead imbalance may freeze, while `Z` cross-coupling can still move an `eta` component whose own reaction channel is dead (§21 B) |
| B3 | thermal RHS at `eta = 0` | Exact same-coefficient channel identity: the extension reduces to its own `L_nu,eq=Ltilde T^q`; the historical placeholder baseline is not an oracle (§21 H) |
| B4 | `L_H`, `DeltaL_nu` at every step | both `>= 0` (§24.2) |
| B5 | `DeltaP_beta` sign change | at `xi = 4.90971` for the modified-Urca channel |
| B6 | quasi-steady endpoint | FR2005 eqs. (64)–(65); scoping agreement 0.15% |
| B7 | `L_gamma,qs` scaling | `∝ \|Omega Omegadot\|^{8/7}` (FR2005 eq. 66) |
| B8 | initial-condition independence | endpoint invariant under `eta(0) in {0, kT, 10kT, 20kT}` and `T_inf(0) in {1e7..1e9 K}` — FR2005 Fig. 5, **and this claim is EOS-independent, so the free gas can test it honestly** |
| B9 | `xi` spatial invariance | `eta^inf/(k_B^(MeV) T_inf)` equals `eta_local(r)/(k_B^(MeV) T_local(r))` at every sampled radius |
| B10 | `W = Z I_Omega` | recomputation matches the shipped artifact (already verified, §6.2) |

These evolution observables consume the declared coefficients and therefore cannot substitute for
RE10b's independent validation of the stellar `Ltilde` construction.

---

## 26. A18 feasibility audit (bounded; no implementation)

### 26.1 Scope and honesty statement

**CONFIRMED FROM BOUNDED PRIMARY-SOURCE PREPRINT INSPECTION:** APR 1998 provides enough form to
reconstruct the core arbitrary-composition A18+delta-v+UIX* functional `E(n_B,x_p)`. The inspected
public preprint remains discovery/feasibility evidence only: no source byte was installed, no
journal-version authority was authenticated, and no A18 implementation began.

### 26.2 What the realistic reproduction actually requires

FR2005 §3.1, p. 9, states the requirement precisely:

> *"In addition to the usual variables given in most published EOSs (pressure, energy density,
> chemical composition, and adiabatic index, all evaluated in chemical equilibrium), we need to know
> the partial derivatives of equation (10) (which are easy to compute if we know the energy density
> or energy per baryon of interacting particles as a function of baryon number density and proton
> fraction) and the baryon effective masses for computing reaction rates and heat capacities."*

So the requirement is exactly:

1. `E(n_B, x_p)` — energy per baryon at **arbitrary composition**, to give `dn_i/dmu_j` off the
   beta-equilibrium line;
2. baryon **effective masses**;
3. the beta-equilibrium EOS itself, for structure.

### 26.3 Public beta-equilibrium tables — sufficient?

**No.** A beta-equilibrium table gives `E`, `P`, `Y_i` **along a one-parameter curve** in
`(n_B, x_p)` space. The off-equilibrium Hessian `dn_i/dmu_j` requires the second derivatives of `E`
**transverse to** that curve. No amount of interpolation along the curve recovers them. This is the
same distinction ADR-0010/ADR-0013 already draw between the equilibrium anchor and the
charge-neutral susceptibility. A CompOSE-style `APRP_1998` table is therefore **necessary but not
sufficient**.

### 26.4 What FR2005 tells us is reconstructable

| Piece | Evidence | Status |
|---|---|---|
| `E(n_B, x_p)` for A18+δv+UIX* | Bounded APR primary-preprint inspection plus FR2005 §3.1 | **Sufficient in form:** effective Hamiltonian/effective masses, proton-fraction interpolation, Appendix-A fit functions, and Table-XII `p1...p21` for A18+delta-v+UIX* LDP and HDP are present. |
| Baryon effective masses | FR2005 §3.4, p. 12: *"The latter can be obtained analytically for the APR and PAL EOSs (see, e.g., Page et al. 2004)"* | Reconstructable, but the cited route is **Page et al. (2004)**, which is **not in the library**. |
| Pion-condensed phase | FR2005 §3.1, p. 9: A18+δv+UIX* has *"a phase transition associated with the appearance of a neutral pion condensate at a density `~4 x 10^14 g cm^-3`"* | The two-phase (normal / pion-condensed) structure is part of the APR fit. |
| Phase-transition construction | FR2005 §3.1, p. 9: *"we assumed a **Maxwell transition** ... resulting in an energy-density jump of **6.6%**"* | The convention and jump are specified **in form** by FR2005. Exact phase matching/reconstruction remains a realistic blocker. |
| Causality | FR2005 §3.1: *"becomes non-causal at densities greater than `2 x 10^15 g cm^-3`, which is the central density of a star of `2.14 M_sun`"*; `M_max = 2.19 M_sun` | **Fully specified by FR2005.** Bounds the usable mass range. |
| Crust | FR2005 §3.1: Pethick, Ravenhall & Lorentz (1995) inner crust, Haensel & Pichon (1994) outer crust | **Two further sources, neither in the library.** |
| Core/crust domain | FR2005 §3.5: `B_ij` integrated **over the core only**; crust processes neglected | Specified by FR2005. |
| Weak-rate normalization | FR2005 §3.2: branches summed *"(Yakovlev et al. 2001)"* | **Not in the library** (§10.4). |
| Envelope | FR2005 §3.3: Potekhin et al. (1997) accreted envelope, eq. (49) given explicitly | Formula given in FR2005; a Potekhin implementation already exists in-tree. |

**DISCOVERY ONLY (recorded, not installed or adopted as repository authority).** The bounded search confirmed the primary's
bibliographic identity — Akmal, Pandharipande & Ravenhall, *Equation of state of nucleon matter and
neutron star structure*, Phys. Rev. C **58**, 1804–1828 (1998), preprint `arXiv:nucl-th/9804027` —
and the inspected preprint supplies the form described above. A CompOSE `APRP_1998` table remains a
beta-equilibrium table and is insufficient per §26.3. No discovered byte enters the library.

### 26.5 Verdict

| Question | Answer |
|---|---|
| Is the arbitrary-composition `E(n_B, x_p)` reconstructable from the primary source? | **Confirmed in form by bounded primary-source preprint inspection**, subject to authenticated journal authority and exact FR2005 construction. |
| Reconstructable pieces in form | Effective Hamiltonian/effective masses, proton-fraction interpolation, Appendix-A fit functions, Table-XII `p1...p21` for A18+delta-v+UIX* LDP/HDP; plus FR2005's Maxwell-transition convention. |
| Still missing for realistic closure | Authenticated APR journal authority; exact FR2005 Maxwell/phase construction; crust joins; authenticated YKGH2001 rate authority; effective-mass interpretation; `alpha_n` choice; direct/MU support authority; Page/crust/envelope authority as required; benchmark arrays/digitization. |
| Is a public beta-equilibrium table sufficient? | **No** (§26.3) |
| Digitization required? | **Not for the EOS** if the analytic fit is confirmed. **Yes, probably, for the R2006 Fig. 1 coefficient benchmark and the FR2005 Fig. 4/6 transient benchmarks**, unless author arrays are obtained — as ADR-0013 Q6 already anticipated. |
| Action in this task | None. No source byte installation and no A18 implementation are authorized. |

---

## 27. Realistic-source blocker ledger

Two **disjoint** blocker sets. They must not be conflated.

### 27.1 (A) First free-gas evolution — blockers

| # | Blocker | Severity | Resolution |
|---|---|---|---|
| A1 | Per-channel equilibrium Urca normalization with an electron/muon split | **Not a blocker to the mathematical benchmark** | Use positive declared benchmark coefficients with recorded provenance, classified as not realistic source normalizations. The same coefficients must feed equilibrium cooling, `F`, `H`, and chemical heating. R1995 is historical/supporting only; YKGH2001 is sufficient in form but its `alpha_n` choice remains unresolved for realistic reproduction. |
| A2 | Channel-resolved `Ltilde` does not exist in `NeutrinoCoolingCachePayload` | Implementation, not source | §13.3 |
| A3 | Explicit DU applicability/support representation and muon criterion absent | Implementation, not source | §15 — not exercised because the benchmark's static process set disables DU; future tests must cover the triangle-open sliver as a disabled/out-of-contract negative control and the disconnected outer-shell sweep hazard |
| A4 | Phase-5C `Z`/`W` are a ratified **candidate**, not canonically integrated | Process | Canonical integration of Phase-5C is a prerequisite for a *production* benchmark, not for the contract |
| A5 | `n_eta = 0` everywhere; no spin-history interface | Implementation | §4.2, §18 |

**Nothing in (A) is a scientific stop condition for the controlled benchmark.** It can validate ODE
wiring, signs, energy bookkeeping, quasi-steady scaling, spin coupling, and thermal coupling. It
cannot validate the FR2005 absolute temperature/history.

### 27.2 (B) Realistic FR2005 reproduction — blockers

| # | Blocker | Status |
|---|---|---|
| B1 | Authenticated APR journal authority and arbitrary-composition reconstruction | Preprint sufficient in form only; exact FR2005 Maxwell/phase construction and crust joins remain |
| B2 | Authenticated YKGH2001 branch/lepton rate normalization and exact `alpha_n` choice | Public and sufficient in form; not installed/authenticated; `alpha_n` deliberately undecided |
| B3 | Effective-mass interpretation and required Page authority | Unresolved |
| B4 | Direct/MU support authority; PRL inner crust; Haensel-Pichon outer crust; Page/crust/envelope authority as required | Unresolved |
| B5 | R2006 Fig. 1 / FR2005 Fig. 3, 4, 6 numerical arrays | Author arrays preferred; governed digitization otherwise (ADR-0013 Q6) |
| B6 | Superfluid gaps | Explicitly out of scope for the non-superfluid programme |

**(B) is unchanged from ADR-0013's finding** that realistic Track-R closure is blocked on
authenticated A18 authority. Phase-5D adds B2 and B3 as *new* rate-specific dependencies that the
coefficient layer did not have.

---

## 28. Future BNV-compatible seam

### 28.1 The structural requirement

BNV changes baryon number. The v1 two-channel space is built by the **fixed-baryon lift**
(ADR-0013 §3.3)

```text
L = [[-1,-1],[1,0],[0,1]] ,   delta N_y = L delta N_l ,   Z = L^T G_y^-1 L
```

which *encodes* `delta N_n = -(delta N_e + delta N_mu)`. A baryon-number-violating source does not
satisfy that constraint and **cannot be represented in the reduced space at all**. Forcing it
through `L` would silently project away the baryon-violating component.

### 28.2 The seam

ADR-0013 §3.2 already retains `G_y` in the **unreduced** source basis `y = (N_n, N_e, N_mu)` as
*"the canonical unreduced authority ... retained as the source-aligned unreduced physical authority
for later non-fixed-baryon sources."* That is the seam, and Phase-5D must not remove it. Concretely:

```text
generic source :   d(delta N_y)/dt |_ext  =  Sigma_y(t, state)          [count/s], 3-vector
reduced path   :   Sigma_y = L Sigma_l  =>  etadot^inf |_ext = -Z Sigma_l     (fixed baryon number)
general path   :   delta g_y^inf = -G_y^-1 delta N_y  ,  eta^inf = -L^T delta g_y^inf
                   — requires the unreduced G_y, not Z
```

Spin-down is one adapter into `Sigma_l`; `RotochemicalSpinDrive` supplies
`Sigma_l = -2 Omega Omegadot I_Omega`. ADR-0013 §3.4 already says `W` *"is a spin-drive adapter, not
a universal chemical source"* and that future baryon-changing sources *"must not be forced into the
fixed-baryon lift `L`."* **Phase-5D inherits that constraint verbatim.**

### 28.3 Source precedent for a generic external source

**JRF2006** (supporting) is the existence proof that the FR2005 formalism accommodates a *different*
external parameter: it replaces `2 Omega Omegadot I_{Omega,i}` with `Gdot I_{G,i}`, where
`I_{G,i} = (dN_i^eq/dG)_A`, and obtains the same ODE structure (its eqs. 5–8, with
`C_npl = (Z_npl - Z_np) I_{G,l} + Z_np I_{G,p}` — structurally identical to FR2005 eqs. 57–58).

**Critical limitation, recorded:** JRF2006's source is still **at constant total baryon number `A`**
(its `I_{G,i}` is defined at fixed `A`, exactly as `I_{Omega,i}` is). So JRF2006 authorizes a generic
**fixed-baryon** external source seam — covering spin-down, `Gdot`, and any other background-parameter
drift — but it does **not** authorize the baryon-violating case. That still needs the unreduced
`G_y` path and its own source derivation.

**No BNV rate is derived here. No BNV heating is claimed. No BNV code is written.**

---

## 29. Existing neutrino-module double-counting audit

Every current Urca contribution in the tree:

| Contribution | File:line | Equilibrium only? | Local/global | e/mu split | D/M split | Normalization | Redshift | Support domain |
|---|---|---|---|---|---|---|---|---|
| `L_nu^DU = K_DU T_inf^6` | `NeutrinoCooling_Details.cpp:953-954`, built at `:188-197` | **Yes** (`xi = 0` implicitly) | Global, cached per profile version | **No** | Yes | `Q0_DU = 1e27 * rho15 * T9^6`, **declared placeholder** | `e^{(2-6)nu} = e^{-4nu}` inside `K`; `4 pi r^2 e^Lambda` proper volume — **correct** | `[0, durca_last]` from `StarContext::DirectUrcaLastAllowedIndex()` |
| `L_nu^MU = K_MU T_inf^8` | `NeutrinoCooling_Details.cpp:957-958`, built at `:184-185` | **Yes** | Global, cached | **No** | Yes | `Q0_MU = 1e21 * rho15 * T9^8`, **declared placeholder** | `e^{-6nu}` inside `K`; correct | whole profile `[0, N-1]` |
| `L_nu^PBF` | `NeutrinoCooling_Details.cpp:960-964` | Hook, `K_PBF = 0` | — | — | — | none | — | none |
| `Microphysics/Rates/Urca.hpp` | not built | — | — | — | — | — | — | **dead** |

**There is exactly one equilibrium Urca path, it is global, it is `xi`-independent, its redshift
convention is correct, and its normalization is an acknowledged placeholder.**

### 29.1 Implementation-ready no-double-counting recipe

1. Extend `NeutrinoCoolingCachePayload` to hold `Ltilde_De, Ltilde_Dmu, Ltilde_Me, Ltilde_Mmu`
   (`erg s^-1 K^-6` and `erg s^-1 K^-8`), each built with the existing integrand structure
   `integral_{D_a} 4 pi r^2 e^Lambda S_a(n) e^{(2-q)nu} dr`, over the channel's declared
   support/applicability subset `D_a subseteq D`.
2. `NeutrinoCooling` computes `L_nu,eq = (Ltilde_De + Ltilde_Dmu) T_inf^6 + (Ltilde_Me + Ltilde_Mmu) T_inf^8`
   from the same declared channel coefficients used by the extension. `NeutrinoCooling` never sees
   `eta`. The existing historical placeholder `NeutrinoCooling` normalization is not the controlled
   benchmark's equilibrium coefficient or a validation oracle. Replacing the default production
   normalization with source-authoritative values remains a later realistic-physics task and this
   clarification changes no historical baseline.
3. A new rotochemical thermal contribution reads the **same** payload plus `ChemState` and adds
   `sum_a Ltilde_a [xi_a H_*(xi_a) - F_*(xi_a) + 1] T_inf^{q_a} / (T_inf C_*)` to `dx/dt`.
4. A new rotochemical chemical contribution reads the same payload and adds
   `-Z R + 2 W Omega Omegadot` to the `Chem` block.
5. **Invariant, testable:** with `eta = 0`, `H(0)=0` and `F(0)=1` make `R`, `DeltaGamma`,
   `DeltaL_nu`, and `L_H` identically zero, so each extended channel reduces exactly to its own
   declared `L_nu,eq=Ltilde T^q`. This is the RE9 same-coefficient identity, not a historical-
   baseline comparison.
6. **Invariant, testable:** `Ltilde` appears in exactly one place. A grep-level test that no second
   Urca normalization constant exists is a legitimate contract test (M-series, §32).

---

## 30. INV-11 subparts

| Subpart | Scope | Status before | Would ADR-0014 resolve it? |
|---|---|---|---|
| **INV-11a** — redshift / coefficient semantics | ADR-0013 coefficient-object scope: `eta^inf = e^nu eta_local`; one `e^{-nu}` in `G_y`; `Z` acts on redshifted imbalance; `W` units/sign. Secular extension: evolved-`eta` ownership, `xi`, global reaction-rate and neutrino-luminosity integrals, and global heating coupling. | **PARTIALLY RESOLVED UPSTREAM** for static `G_y`/`Z`/`W` coefficient-object semantics only | **PROPOSED COMPLETION / EXTENSION HERE** for the secular-evolution layer; not globally resolved |
| **INV-11b** — evolved `eta` state ownership | variable, units, ordering, storage, initial condition, domain | **UNRESOLVED** | **Proposed resolution** (§5) |
| **INV-11c** — reaction sign/index convention | `DeltaGamma` direction, `H_*` parity, channel indices, `Z` orientation | **UNRESOLVED** | **Proposed resolution** (§7, §6.4) |
| **INV-11d** — thermal energy ledger/no double counting | `L_H`, `L_nu(eta)`, `DeltaL_nu`, conversion boundary | **UNRESOLVED** | **Proposed resolution** (§12, §13, §24) |
| **INV-11e** — frozen coefficient lifetime/update policy | frozen `Z`, `W`, `Ltilde`; `Zdot` term if relaxed | **UNRESOLVED** | **Proposed resolution for frozen v1**; time-dependent case recorded and forbidden |
| **INV-11f** — ODE/source coupling | state vector, spin ownership, solver, tolerances, ordering | **UNRESOLVED** | **Partially proposed**; numerical solver implementation and validation remain future (§23) |

**INV-11 must NOT be marked globally resolved.** ADR-0014 is PROPOSED, not accepted; nothing is
implemented; INV-11e's time-dependent branch and INV-11f's stiff-solver branch remain open by
design. ADR-0013's accepted scope does not itself govern evolved-state ownership, `xi` use in
evolution, global reaction-rate integrals, global neutrino-luminosity integrals, or global heating
coupling.

---

## 31. RE validation ladder

Classification: **IND** = independent analytic/numerical oracle; **STR** = source traceability;
**CON** = contract; **CNV** = convergence; **SL** = source-limited.

| ID | Gate | Class | Oracle |
|---|---|---|---|
| **RE1** | Units and channel ordering | CON | `eta` in MeV; index 0 = `Npe`, 1 = `NpMu`; `Z` in MeV/count; `W` in MeV s^2; `Ltilde` in erg s^-1 K^-q; `k_B^(MeV)` and `k_B^(erg)` are derived views of one repository Boltzmann authority through the governed energy conversion |
| **RE2** | `eta` redshift and `xi` spatial invariance | IND | `eta_local(r) = eta^inf e^{-nu}`; `eta^inf/(k_B^(MeV) T_inf) = eta_local(r)/(k_B^(MeV) T_local(r))` at every sampled radius |
| **RE3** | Source polynomial coefficients | STR | The four polynomials of §9.2, coefficient by coefficient; §9.3 classifies the `H_M` exponent as a confirmed printed typo/internal inconsistency, with no published erratum claimed |
| **RE4** | Parity and normalization | IND | `F_*(0)=1`, `H_*(0)=0`, `F_*` even, `H_*` odd, `xi H_*(xi) >= 0`, `F_* >= 1` |
| **RE5** | Small-`xi` limit | IND | `H_D -> 0.158300492605 xi`, `H_M -> 0.129192649689 xi`, `F -> 1 + O(xi^2)` with the §9.5 coefficients |
| **RE6** | Large-`xi` asymptotics and sign-crossing roots | IND | Leading powers; `C_H = 24/(11513 pi^8)`, `C_M = 15/(11513 pi^8)`; heating fractions `1/2`, `5/8`; roots `4.7870134733369`, `4.90971002892413` (incremental) and `5.4585315948676`, `5.63371746764834` (full) |
| **RE7** | Reaction sign, `eta·DeltaGamma >= 0`, and Lyapunov qualification | IND | Both signs of `eta`; for frozen SPD `Z`, `Vdot=-2 sum eta_l R_l<=0`; strict decrease requires an active positive-normalization dissipative channel in every nonzero direction; an uncoupled dead imbalance may freeze, while cross-coupling can move an `eta` component whose own reaction channel is dead |
| **RE8** | Neutrino equilibrium limit | CON | `eta = 0 ⇒ DeltaL_nu = 0` and `L_H = 0` **identically** (structural zero, not tolerance-based) |
| **RE9** | Same-coefficient equilibrium/no-double-counting identity | CON | For each channel using the same declared `Ltilde` for equilibrium and extension: `eta=0 => H=0, F=1, R=DeltaGamma=DeltaL_nu=L_H=0`, and the extended channel equals its own `L_nu,eq=Ltilde T^q`; historical placeholders are not the oracle; exactly one coefficient authority |
| **RE10** | Spin-only analytic `eta` source | IND | `Ltilde = 0 ⇒ eta_l(t) = W_l [Omega^2(t) - Omega_0^2]` (FR2005 eq. 77) |
| **RE10b** | **GLOBAL URCA COEFFICIENT GR-INTEGRAND ORACLE** | **IND** | Exercise `GlobalUrcaChannelCoefficient` itself on an explicit synthetic radial background with declared nonconstant `nu(r)`, `lambda(r)`, local `S_a(r)`, support domain `D_a`, and fixed `q_a`. The independent expected value is `Ltilde_a = integral_{D_a} 4 pi r^2 e^{lambda(r)} S_a(r) e^{(2-q_a)nu(r)} dr`, with `4 pi` and every ordinary unit prefactor declared, including the exact length-volume conversion when the fixture coordinate is not in cm. Compute it outside the production integration kernel by a closed-form fixture where practical and an independent high-precision quadrature. The gate must reject `e^{(1-q)nu}`, `e^{(3-q)nu}`, `e^{-(2-q)nu}`, omitted or inverted `e^lambda`, direct `q=8` / modified `q=6`, and a wrong `D_a`. A precomputed or injected `Ltilde` may not be the expected value. |
| **RE11** | Reaction-only relaxation | IND | `Omegadot = 0 => V` nonincreasing; strict decay to zero only under the RE7 all-active qualification; an uncoupled dead imbalance may freeze; small-`xi` active-channel rate matches RE5 |
| **RE12** | Coupled toy analytic solution | IND | Linearized two-channel system at two distinct constant `T_inf` values, using `k_B^(erg)` in the rate normalization and `k_B^(MeV)` in `xi`, including `Z_np` cross-coupling; the two-temperature comparison independently detects `T_inf^q` substituted for `T_inf^(q-1)` in `R_l` |
| **RE13** | Free-gas end-to-end numerical evolution | CNV | §25.5 observables B1–B10; tolerance-convergent under `rtol`, `atol`, and profile refinement. Its declared-`Ltilde` input validates state evolution, reaction/heating algebra, thermal coupling, and spin coupling, but not stellar `Ltilde` construction; RE10b owns that layer. |
| **RE14** | Quasi-steady asymptote | IND | FR2005 eqs. (64)–(65) with `Z` cancelling; scoping agreement 0.15% |
| **RE15** | FR2005 source benchmark | **SL** | R2006 Fig. 1 coefficients; FR2005 Fig. 4/6 transients; FR2005 eqs. (67)–(69) brackets. **BLOCKED** on §27.2 |
| **RE16** | Provenance and staleness | CON | Every evolution result carries `Z`/`W`/`Ltilde`/profile/spin-history identity plus the chemical/reaction domain `D` and channel support/applicability domain `D_a`; a changed dependency refuses before scientific access (ADR-0013 §6) |
| **RE17** | Direct/modified support domains | CON | Triangle support, declared benchmark process support, and future physical applicability are distinct. The low-density electron-DU triangle sliver must not activate the kernel when DU is disabled/outside the declared channel-support contract; an outer allowed shell must not sweep a closed inner region; any ordering-dependent integration asserts the profile is innermost first; muon support is separately represented (§15.2). No temperature-dependent support predicate is introduced into frozen v1. |
| **RE18** | Frozen-coefficient contract | CON | `Z`, `W`, `Ltilde` byte-identical at every RHS evaluation of a run; time-dependent `Z` refuses |

RE1–RE14 (including RE10b, which uses its own synthetic coefficient fixture) and RE16–RE18 are
achievable with controlled fixtures. **RE15 alone is source-limited.**

---

## 32. Mutation plan

Mutations a future test suite must kill. Distinct mutations are separated from algebraic aliases.

### 32.1 Genuinely distinct mutations

M11 is retained in §32.2 as the historical label for the same integrated lapse-power mutation as
M10; it is deliberately not counted as an additional independent mutation.

| ID | Mutation | Killed by |
|---|---|---|
| M1 | Flip the sign of `eta` (store `+Z delta N`) | RE7, RE10 (`eta` would grow negative under spin-down), RE14 |
| M2 | Flip the sign of `DeltaGamma` (use `-H_*`) | RE7 (`eta·DeltaGamma < 0`), RE11 (`eta` diverges instead of relaxing) |
| M3 | Flip the sign of `Omegadot` | RE10 (`eta` decreases under spin-down), RE14 |
| M4 | Flip the sign of `W` | RE10, RE14 |
| M5 | Swap the `npe`/`npmu` **state** channels | RE13/B6 (wrong `eta_e`/`eta_mu` ratio), RE14 |
| M6 | Transpose the semantic channel map (`Z` row/column roles swapped) | RE12 (cross-coupling eigenvectors), RE14; **not** killed by symmetry since `Z` is symmetric — needs the `Ltilde` asymmetry |
| M7 | Drop the `Z_np` cross-coupling (diagonal `Z`) | RE12, RE13/B6 |
| M8 | Use `eta_local` with `T_inf` in `xi` | RE2 (`xi` no longer spatially constant) |
| M9 | Double or omit the redshift in `xi` (`e^{2nu}` or `e^{0}`) | RE2 |
| M10 | In the shared integrated-`Ltilde` construction, omit one of the two `e^{+nu}` factors, yielding `e^{(1-q)nu}` instead of `e^{(2-q)nu}` | RE10b direct `GlobalUrcaChannelCoefficient` oracle. M11 is the algebraically equivalent historical layer label and is not counted separately (§32.2). |
| M12 | Substitute the `G_y` inverse-lapse route into the shared coefficient construction, yielding `e^{-q nu}` instead of `e^{(2-q)nu}` | RE10b; explicitly distinguished from `G_y` by ADR-0013 §3.2 |
| M13 | Double-count equilibrium Urca (add full `L_nu(T,eta)` **and** keep the same-coefficient equilibrium channel) | RE8/RE9 (the declared equilibrium term appears twice) |
| M14 | Subtract the **full** `L_nu(T,eta)` instead of `DeltaL_nu` in the incremental channel | RE9; also detected by the sign-crossing root moving from `4.910` to `5.634` (RE6) |
| M15 | Omit chemical heating `L_H` entirely | RE6 (no sign crossing at all), RE14 (no quasi-steady) |
| M16 | Use `eta_local` in the global heating integral without the `e^{-Phi}` correction | RE13, RE14 |
| M17 | Direct-Urca support mismatch (rotochemical applicable DU support != cooling applicable DU support) | RE17 |
| M18 | Make `Z` time-dependent without adding the `+Zdot Z^-1 eta` term | RE18 (coefficient byte-identity), plus energy-ledger drift |
| M19 | Freeze `eta` (never integrate the `Chem` block) | RE10, RE13 |
| M20 | Muon/electron **normalization** swap (`Ltilde_Me <-> Ltilde_Mmu`) | RE13/B6, RE14 — distinct from M5 because it swaps the *rate* coefficient, not the state slot |
| M21 | Introduce a second, independently normalized `Ltilde` for the correction | RE9's "exactly one authority" check. **Not killed by RE8**, since `eta = 0` still gives zero — this is precisely why RE9 must check authority uniqueness structurally, not only the `eta = 0` limit |
| M22 | Use `H_*` where `F_*` belongs (or vice versa) | RE4 (parity), RE8 (`F_*(0)=1` vs `H_*(0)=0`) |
| M23 | Use `pi^6` in the last `H_M` term (the FR2005 printed typo) | RE3, RE6 (large-`xi` heating fraction becomes `(24 pi^2 - 9)/(24 pi^2) != 5/8`) |
| M24 | Apply a positivity clamp to `eta` | RE7 with `eta(0) < 0`; physically legitimate for a spun-up star |
| M25 | Obtain `Omegadot` by re-deriving the torque inside the rotochemical driver | RE13 with a prescribed spin history that differs from the dipole law |
| M26 | Omit `k_B^(erg)` from the rate normalization | RE1 dimensional check; RE12 |
| M27 | Double `k_B^(erg)` in the rate normalization | RE12; RE14 |
| M28 | Use the wrong Boltzmann energy unit: MeV/K against erg-normalized `Ltilde`, or erg/K against a MeV-normalized coefficient | RE1 dimensional/authority check; RE12 |
| M29 | Use `T_inf^q` rather than `T_inf^(q-1)` in `R_l` | RE12 independent two-temperature analytic solution |
| M30 | Integrate reaction `Ltilde` over a domain different from the associated `G_y`/declared chemical domain | RE10b explicit `D_a` restriction and wrong-domain mutant |
| M31 | Activate the electron-DU kernel in the low-density triangle-open sliver even though DU is disabled/outside the declared channel-support contract | RE17 applicability/process-configuration negative control |
| M32 | Select an outer DU-allowed shell by last index and incorrectly integrate the closed inner region | RE17 explicit-support and innermost-first negative control |
| M33 | Add one extra lapse factor in `GlobalUrcaChannelCoefficient`, using `e^{(3-q)nu}` | RE10b |
| M34 | Use the wrong sign in the global coefficient lapse exponent, `e^{-(2-q)nu}` | RE10b |
| M35 | Omit the proper-volume curvature factor `e^lambda` | RE10b |
| M36 | Invert the proper-volume curvature factor, using `e^{-lambda}` | RE10b |
| M37 | Use the wrong channel exponent in coefficient construction: direct with `q=8` or modified with `q=6` | RE10b |

### 32.2 Algebraic aliases — NOT independent mutations

Recorded so a future suite does not over-count its own coverage:

- **M1 ≡ M2 ≡ M4** in the *reaction-only* limit: all three flip the sign of the same product. They
  are distinguished **only** when both the reaction and spin terms are active; a reaction-only test
  cannot separate them.
- **M3 ≡ M4** in the spin-only limit (both flip `2 W Omega Omegadot`). Separated only by a run in
  which `Omegadot` changes sign, or by an independent check of `W = Z I_Omega` (RE13/B10).
- **M5 ≡ M20** whenever `Ltilde_Me = Ltilde_Mmu`. On the free-gas fixture they are well separated:
  the scoping calculation gives `Ltilde_Mmu/Ltilde_Me ~ 0.05`.
- **M6** is invisible in any test that uses a symmetric `Z` with equal channels; it requires the
  asymmetric `Ltilde` and unequal `I_Omega` of the real fixture.
- **M9 (`e^{0}`) ≡ M8** when the star is evaluated at a single radius. A multi-radius `xi` check is
  required.
- **M10 ≡ M11** in the chosen shared, integrated single-`Ltilde` representation. Calling the missing
  `e^{+nu}` factor “time dilation” (M10) or “the second energy-redshift factor” (M11) produces the
  same mutated `GlobalUrcaChannelCoefficient` integrand `e^{(1-q)nu}`. They are one independent
  mutation detection, killed by RE10b; M11 is retained as a historical label and does not receive
  separate coverage credit.
- **M26, M27, and M28 are distinct unit/normalization faults**, not algebraic duplicates: omitted,
  factor-two, and energy-unit mismatch failures require separate future controls. M28 groups its two
  reciprocal cross-unit directions to avoid inflating the mutation count with algebraic duplicates.

---

## 33. Exact implementation recommendation

### 33.1 Smallest scientifically complete Phase-5D-1

| # | Object | Responsibility |
|---|---|---|
| 1 | `ChemicalImbalanceState` (immutable value) | Two-channel `eta^inf` in MeV, named order `(Npe, NpMu)`, with `Xi(k_B^(MeV), T_inf)` and `Local(nu)` accessors. Backed by `ChemState` storage in the ODE vector. |
| 2 | `UrcaImbalanceFunctions` (pure, stateless) | `F_D, H_D, F_M, H_M, M_D, M_M` exactly as §9.2, plus the four root constants of §14.2 as named compile-time oracles. No stellar or state dependence. |
| 3 | `UrcaChannelMicrophysics` / `UrcaChannelNormalization` | Process, lepton, nucleon branch, effective masses, matrix-element and alpha/beta factors, local equilibrium `S_a(r) [erg cm^-3 s^-1 K^-q_a]`, support/applicability metadata, and source provenance. No stellar integration. |
| 4 | `GlobalUrcaChannelCoefficient` | Stellar GR integration of provider (3) over declared `D_a subseteq D` only, yielding channel-resolved `Ltilde_a [erg s^-1 K^-q_a]` with `D`/`D_a` and profile-versioned provenance. No F/H or eta. Its lapse, proper-volume, exponent, and domain semantics are independently validated by future gate RE10b. |
| 5 | `RotochemicalReactionResponse` | Combines (2)+(4)+`T_inf`+`eta^inf` to provide `R_l`, `DeltaL_nu^inf`, and the MeV/s chemical-heating contribution with its one thermal-boundary conversion. |
| 6 | `SecularEvolutionRHS` | Combines channels and adds `-Z R + 2 W Omega Omegadot` to `Chem` and the thermal contribution to `Thermal`; consumes `Z`, `W`, and `ISpinHistory`. |
| 7 | `ISpinHistory` | `Omega(t, Y)`, `OmegaDot(t, Y)`. Prescribed and state-coupled implementations. |
| 8 | Equilibrium cooling adapter | Consumes the same (4) channel coefficients for `L_nu,eq`; historical placeholders are not an oracle. |
| 9 | Free-gas benchmark | §25. |

**Normative separation:** local/source microphysics normalization, global stellar GR integration,
pure dimensionless F/H mathematics, the reaction response, and the secular RHS are five separate
layers. No object may own more than one. In particular, the normalization provider cannot integrate
the star; the global coefficient cannot choose microphysics or evaluate F/H; and
`UrcaImbalanceFunctions` knows neither stars nor state.

### 33.2 Explicitly reused unchanged

`ThermalState`; `C_*(T_inf)` and ADR-0002; `PhotonCooling` and the envelope models; `SpinState` and
`MagneticDipole`; `EvolutionSystem`, `RHSAccumulator`, `StateLayout`, `StatePacking`, `GSLIntegrator`
with `RKF45`; `ChemState` storage; the `k_B` authority; `ChemicalImbalanceResponse`; and
`RotochemicalSpinDrive` including `Evaluate(Omega, Omegadot)`. Only the current fixture's static
MODIFIED-URCA-ONLY process configuration is reused. Any enabled, nonempty, or outer-shell DU adapter
first requires explicit support representation plus the §15/RE17 ordering and applicability
controls.

### 33.3 Explicitly excluded

Superfluidity (**not included**); A18 implementation (**not begun**); BNV (**not begun**);
time-dependent `Z`/`W`; envelope `eta` dependence; crustal processes; hyperons; solver replacement;
source-authoritative replacement of the default production placeholder equilibrium normalization
(a separate governed realistic-physics change);
any change to ADR-0002, ADR-0010, ADR-0011, ADR-0012 or ADR-0013 semantics.

### 33.4 Disposition

> **PHASE-5D MATERIAL CLOSURE COMPLETE — R1/R2/R3 CLOSED —
> READY FOR FINAL BOUNDED INDEPENDENT RE-REVIEW.**

Disposition **A** for this bounded author revision. R1 now rests on the benchmark's static
MODIFIED-URCA-ONLY process set; R2 now has a dedicated future `GlobalUrcaChannelCoefficient` oracle;
and R3 restores INV-11a to partial upstream resolution with a proposed secular extension here.
ADR-0014 remains PROPOSED, and this preflight is **REVIEWED / REVISED / AWAITING FINAL BOUNDED
RE-REVIEW**. The separately classified APR/YKGH discoveries do not close the realistic blockers.

No new unresolved scientific ambiguity was introduced by this material-closure revision. The
FR2005 eq. (37) issue remains classified as a confirmed printed typo/internal source inconsistency,
with no published erratum located (§9.3).

### 33.5 Scratch evidence (not committed)

`scratchpad/calc/urca_derive.py` (phase-space derivation of `F_*`, `H_*`);
`roots.py` (sign-crossing roots, parity, limits);
`freegas_du.py` (free-gas composition and DU threshold);
`freegas_star.py`, `freegas_evol.py` (scoping TOV and `Ltilde`/`Ctilde`);
`evolve.py`, `stiff.py` (scoping evolution, quasi-steady check, stiffness).
None is production code; none is committed; none is a validated numerical result.

---

## 34. Post-independent-review revision — 2026-09-07

The owner-supplied Phase-5D-0R Opus review returned prior disposition C,
`PHASE-5D CONTRACT REVIEW BLOCKED — ADR-0014 REVISION REQUIRED BEFORE RATIFICATION`, with zero
blocking scientific derivation failures and eight material documentation/contract findings. This
revision does not reopen any independently confirmed load-bearing physics. The table below preserves
the revision record as it stood at `PHASE5D0_REVISION_SHA`; Phase-5D-0RR later found the recorded
M-2 closure partial, as corrected in §35.

| Finding | Original defect | Revision | Remaining blocker | Status |
|---|---|---|---|---|
| M-1 | Reaction normalization mixed an MeV/K Boltzmann view with an unspecified luminosity energy unit | Canonical `Ltilde` is `erg s^-1 K^-q`; rates use derived `k_B^(erg)`; `eta R` is canonical MeV/s and converts once at the thermal boundary; dimensional proof and mutations added | None for text closure; implementation and validation remain future | **CLOSED BY TEXT REVISION** |
| M-2 | “DU closed everywhere” conflated triangle support with model applicability | Low-density electron triangle sliver, keV-scale non-degeneracy, empty applicable support, `nB_min` semantics, explicit-support/outer-shell hazard, and negative controls recorded | None for controlled fixture; future adapter/test work remains | **CLOSED BY TEXT REVISION** |
| M-3 | Unqualified global strict Lyapunov decay | `Vdot=-2 sum eta_l R_l<=0`; strictness now requires an active positive-normalization dissipative channel in every nonzero direction; an uncoupled dead imbalance may freeze | None for text closure | **CLOSED BY TEXT REVISION** |
| M-4 | RE9 promised broad bit-for-bit agreement with historical passive cooling | RE9 is the exact same-coefficient identity; benchmark equilibrium cooling uses the same declared `Ltilde`; historical placeholders are not the oracle | None for text closure | **CLOSED BY TEXT REVISION** |
| M-5 | Normative global integrals used implicit core/domain notation | All reaction/coefficient integrals use declared `D`/`D_a`, matched to `G_y` and support; free-gas whole-star domain is stated without arbitrary cutoff | Realistic core/support domain remains governed-source work | **CLOSED BY TEXT REVISION** |
| M-6 | Printed `H_M` issue was called an authenticated erratum | Classified as confirmed printed typo/internal source inconsistency with no published erratum located; independent evidence listed; `pi^8` remains proposed normative | None for text closure | **CLOSED BY TEXT REVISION** |
| M-7 | One proposed object owned microphysical normalization and stellar integration | Five normative layers now separate local/source microphysics, global GR integration, pure F/H functions, reaction response, and secular RHS | Implementation remains future | **CLOSED BY TEXT REVISION** |
| M-8 | R1995 was overstated and YKGH2001/APR discovery understated | R1995 qualified; YKGH2001 public and sufficient in form with `alpha_n` unresolved; APR feasibility upgraded from bounded primary-preprint inspection; realistic blockers retained | Authenticated realistic authorities/construction and benchmark data remain blocked | **CLOSED BY TEXT REVISION** |

The first controlled free-gas secular evolution remains possible as a mathematical/architecture
benchmark using declared positive electron and muon coefficients. The same coefficients must feed
equilibrium cooling, the `F` correction, the `H` rate, and chemical heating. DU is dormant because
the benchmark's static enabled-process set excludes `De` and `Dmu`, so `D_De = D_Dmu = empty`. It
can validate ODE wiring, signs, energy bookkeeping, quasi-steady scaling, spin coupling, and thermal
coupling, but not the stellar construction of `Ltilde` or the FR2005 absolute temperature/history.

Realistic FR2005 reproduction remains separately blocked on authenticated APR journal authority and
reconstruction, exact phase/crust construction, authenticated YKGH2001 normalization, the `alpha_n`
decision, effective masses, direct/MU support authority, required Page/crust/envelope authority, and
published benchmark arrays or governed digitization. Global INV-11 remains unresolved: INV-11a is
only **PARTIALLY RESOLVED UPSTREAM** for coefficient-object semantics and ADR-0014 proposes its
secular completion/extension; INV-11b–e have proposed resolutions and INV-11f is partially
proposed. Nothing becomes accepted until human-owner ratification.

---

## 35. Post-Phase-5D-0RR material-closure revision — 2026-09-08

The independent Phase-5D-0RR re-review returned zero blocking findings and three material findings.
This bounded author revision changes only the contract and validation-plan text needed to close
those findings, plus the directly adjacent enumerated cleanups. It does not implement the future
oracle, alter production or tests, acquire a source, select `alpha_n` or realistic `Ltilde` values,
begin A18 or BNV, or reopen independently confirmed physics.

| Item | Review defect | Corrected contract | Remaining future implementation obligation | Status |
|---|---|---|---|---|
| R1 — DU-support applicability | The low-density electron-DU triangle sliver was called non-degenerate unconditionally, making a temperature-dependent observation the purported static exclusion rule. | The controlled benchmark is statically MODIFIED-URCA-ONLY: `{Me, Mmu}` enabled and `{De, Dmu}` disabled, so `D_De = D_Dmu = empty`. Triangle support, benchmark process support, and future physical applicability are distinct. `nB_min` is only a numerical/semantic guard. The sliver remains a temperature-aware factual observation and negative control. | Before realistic DU evolution, implement a governed applicability provider covering kinematics, composition, degeneracy/model validity, and possibly disconnected support intervals; exercise the sliver and outer-shell hazards. No temperature-dependent support predicate is introduced into frozen v1. | **CLOSED BY TEXT REVISION** |
| R2 — global-`Ltilde` GR-factor detector | M11 pointed to gates that could not independently detect a wrong lapse power, while the declared-`Ltilde` benchmark did not exercise stellar coefficient construction. | RE10b directly exercises `GlobalUrcaChannelCoefficient` against an independent closed-form where practical plus high-precision quadrature value for `integral_{D_a} 4 pi r^2 e^lambda S_a(r) e^{(2-q_a)nu(r)} dr`, with all ordinary unit prefactors declared. It covers missing/extra/wrong-sign lapse powers, omitted/inverted proper volume, wrong `q`, and wrong domain. M10 and M11 are recorded as one algebraic mutation in the shared integrated-`Ltilde` representation. | Implement RE10b independently of the production integration kernel. The present revision declares but does not implement the oracle. | **CLOSED BY TEXT REVISION** |
| R3 — INV-11a governance status | INV-11a was incorrectly marked already resolved upstream. | INV-11a is **PARTIALLY RESOLVED UPSTREAM; PROPOSED COMPLETION / EXTENSION HERE**. ADR-0013 covers static `G_y`/`Z`/`W` coefficient-object semantics only; ADR-0014 proposes evolved-`eta`, `xi`, global-rate, global-luminosity, and heating-coupling semantics. Global INV-11 remains UNRESOLVED. | Final bounded independent re-review and human-owner ratification remain required; no accepted ADR, invariant closure, or implementation is implied. | **CLOSED BY TEXT REVISION** |

All three Phase-5D-0RR material findings are **CLOSED BY TEXT REVISION**. ADR-0014 remains
**PROPOSED**, and this preflight remains **REVIEWED / REVISED / AWAITING FINAL BOUNDED RE-REVIEW**.
