# ADR-0014 — secular rotochemical evolution contract

**Status:** PROPOSED
**Decision:** NOT ACCEPTED — awaiting independent scientific review and owner ratification
**Date:** 2026-09-07
**Starting canonical SHA:** `49ab2b8c2881b6ef7b9309307d18cea51d557f72`
**Branched from (Phase-5C human-ratified candidate):** `27727016856a6a25a46e447c70e380722ea8ddbf`
**Change class:** scientific-semantic and architectural contract; documentation-only proposal.
**Evidence companion:** `docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md`
**Implementation state:** no production evolution object, weak rate, imbalance function, chemical or
thermal RHS term, superfluid treatment, A18 model, BNV source, test, baseline, EOS/data, or
literature byte is created, modified, or authorized by this proposal.

> **PHASE-5D SECULAR ROTOCHEMICAL EVOLUTION CONTRACT — PROPOSED.
> NOT RATIFIED. NOT IMPLEMENTED.**

---

## 1. Context and authority

ADR-0002 (ACCEPTED) fixes one thermal denominator `C_*(T_inf)` and already admits chemical heating on
its right-hand side. ADR-0003 (ACCEPTED) governs profile-derived cache provenance. ADR-0004
(ACCEPTED) governs the proper-volume measure. ADR-0010 (ACCEPTED) governs the local cold
charge-neutral conjugates. ADR-0011 (ACCEPTED) governs the structural particle-number response;
INV-09 is VERIFIED/RESOLVED in that scope. ADR-0012 (ACCEPTED) governs the relativistic unit
boundary. ADR-0013 (ACCEPTED) governs `G_y`, `Z`, `W`, coefficient redshift semantics, and the
active-branch contract; its Phase-5C implementation is human-ratified but **not yet canonically
integrated**. INV-11 remains UNRESOLVED.

This ADR governs **only** the layer between the ADR-0013 coefficients and a first controlled secular
evolution: the evolved chemical state, the reaction-rate convention, the global GR integrals, the
thermal ledger, the frozen-coefficient contract, spin ownership, and the first free-gas benchmark.

**Primary sources** (SHA-256 verified in this task; library read-only throughout):

| Short name | SHA-256 | Role |
|---|---|---|
| FR2005 | `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88` | Primary non-superfluid formalism |
| R2006 | `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8` | Primary corrected coefficients; **supersedes FR2005 eqs. (9), (12), (13), (54)–(56)** |
| R1995 | `9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff` | Primary origin of the imbalance functions |
| Y2020 | `69590539c275fa679a5521a9c5abedd9fdc58718b1554827786c8f707bc618cc` | **Supporting cross-check only** — may never originate a convention |
| JRF2006 | `2dd5444d19cebae12509fe4ecb7dac31957d332e2131813894a10cceac403109` | Supporting precedent for a generic fixed-baryon external source |
| R1997 | `19a10133511aefc05ece33d7454cba60def03de85d9980c55fd7229058c2d08b` | Supporting physical context |

---

## 2. Scope

**Covered.** Evolved chemical state variables and ordering; `eta` sign, units and redshift; the
dimensionless imbalance `xi`; the two-channel chemical ODE; the reaction-rate sign convention; global
GR redshift factors; the non-superfluid Urca imbalance functions and their exact coefficients;
neutrino/heating semantics and the no-double-counting architecture; heat-capacity coupling; photon
cooling invariance; frozen-coefficient lifetime; spin-driver ownership; dependency direction;
provenance; the first free-gas implementation scope; exclusions; the validation ladder; and the
unresolved realistic-source blockers.

**Not covered.** Production implementation; realistic A18 closure; replacement of the existing
placeholder equilibrium Urca normalization; time-dependent coefficients or evolving backgrounds;
superfluidity; BNV; crustal processes; hyperons; envelope `eta` dependence; solver replacement.

---

## 3. Ratified-on-acceptance contract

### 3.1 Evolved state

```text
y_chem = ( eta_npe^inf , eta_npmu^inf )        [MeV]
index 0 = BetaChannel::Npe    index 1 = BetaChannel::NpMu
```

The redshifted imbalances are the evolved variables; `delta N_e`, `delta N_mu` are **not** evolved.
Authority: FR2005 §2.1 and eqs. (52)–(53); R2006 §3, p. 570 (*"the ideal variables to quantify the
departure from chemical equilibrium and follow its time evolution"*) and §5, p. 571; Y2020
eqs. (4.10)–(4.11).

```text
eta_npe   = delta mu_n - delta mu_p - delta mu_e
eta_npmu  = delta mu_n - delta mu_p - delta mu_mu
eta^inf   = e^nu eta_local ,     nu = Phi          (ADR-0013 §3.4, unchanged)
eta^inf   = - Z delta N_l ,      Z symmetric positive definite, [MeV/count]
Z = [[Z_npe, Z_np], [Z_np, Z_npmu]]      row = output channel, column = input lepton
```

Authority for the sign: **R2006 eq. (15)**; Y2020 eq. (4.7); and the independent derivation in
preflight §6.2, where any other sign breaks the exact cancellation of the spin term against
`W = Z I_Omega`. Default initial condition `eta(0) = 0` is **conventional, not required**: the
contract must accept arbitrary `eta(0)` so FR2005 Fig. 5 (right) is reproducible.

Validity domain: `|eta| << mu_i^eq`, on the connected diffusive chemical domain declared by `G_y`.
`ChemState` provides storage only; channel semantics belong to the rotochemical module.

### 3.2 Dimensionless imbalance

```text
xi_a = eta_a^inf / (k_B T_inf) = eta_a(r) / (k_B T(r))        — spatially constant
k_B  = Zaki::Physics::K_BOLTZ_EV * 1e-6 = 8.61733326214518e-11 MeV/K   (single existing authority)
```

Authority: FR2005 eq. (42); Y2020 eq. (4.18). Spatial constancy follows exactly from
`T_inf = T e^nu` **and** `eta^inf = eta e^nu` using the same `nu`, and is the structural reason the
`Ltilde` factorization works. **No new MeV/K constant may be introduced.**

### 3.3 Reaction-rate convention

Forward direction is **neutron decay**:

```text
A -> B :  n (+N_1) -> p (+N_2) + l + nubar_l
DeltaGamma_a = Gamma_{A->B} - Gamma_{B->A}          [count cm^-3 s^-1]
eta_a        = delta mu(A) - delta mu(B)
```

`DeltaGamma_a > 0` means **net neutron decay**: neutrons destroyed, protons and leptons `l` created.
Authority: FR2005 eqs. (39)–(40); Y2020 eq. (4.8) (which states every species sign explicitly).

Required consequences, all provable and testable:

```text
DeltaGamma_a = (1/(k_B T)) Q_a^eq H_*(xi_a)      H_* odd,  H_*(0) = 0
eta_a DeltaGamma_a >= 0   pointwise, both signs of eta       (thermodynamic dissipation)
V = eta^T Z^-1 eta  is a strict Lyapunov function of the reaction sub-system
```

**Warning, normative:** the literature contains both signs. Y2020 fn. 3, p. 61, and R1995 fn. 3,
p. 14, each record a sign opposite to cited references. Any future rate source must be re-anchored to
this table before use.

### 3.4 Chemical ODE

```text
etadot_npe^inf  = - Z_npe R_e - Z_np    R_mu + 2 W_npe  Omega Omegadot
etadot_npmu^inf = - Z_np  R_e - Z_npmu  R_mu + 2 W_npmu Omega Omegadot

R_l = int_core dV e^{Phi} sum_{a in l} DeltaGamma_a
    = (1/k_B) [ Ltilde_{D,l} H_D(xi_l) T_inf^5 + Ltilde_{M,l} H_M(xi_l) T_inf^7 ]     [count/s]

W = Z I_Omega ,     I_{Omega,l} = (dN_l^eq/dOmega^2)_A
```

`R_l` is the **global, `e^{+Phi}`-weighted, channel-summed** net lepton-creation rate — not a local
density. The ADR-0013 §3.4 schematic `etadot^inf = -Z R + 2 W Omega Omegadot` is thereby given its
exact content. Authority: FR2005 eqs. (14)–(15), (30), (45), (52)–(53); R2006 eqs. (5)–(7), (15),
(19); Y2020 eqs. (4.8)–(4.12).

For `Omega > 0`, `Omegadot < 0`: `W < 0` componentwise, so `2 W Omega Omegadot > 0` and **spin-down
drives both imbalances positive** (FR2005 pp. 13–14).

`Z_np` is a **shared** off-diagonal coupling the two channels through the neutron/proton reservoir.
Dropping it is forbidden.

### 3.5 Global GR redshift factors

| Quantity | Net factor over `dV = 4 pi r^2 e^Lambda dr` | Authority |
|---|---|---|
| Heat capacity `C_*(T_inf)` | `e^0`, with `T_local = T_inf e^{-nu}` inside `c_V` | FR2005 eqs. (3),(50); ADR-0002 |
| Global net reaction rate | **`e^{+nu}`** (proper-time dilation) | FR2005 eq. (15); R2006 eq. (6); Y2020 eq. (4.8) |
| Neutrino luminosity at infinity | **`e^{+2nu}`** (time dilation × energy redshift) | FR2005 eq. (5) |
| Chemical heating at infinity | **`e^{+2nu}`** on `Q_H`; equivalently `e^{+nu}` on `DeltaGamma` with `eta^inf` outside | FR2005 eqs. (4),(38); Y2020 eq. (4.5) |
| Chemical response `G_y` | **`e^{-nu}`** (redshifted potential) | FR2005 eq. (12); R2006 eq. (13); ADR-0013 §3.2 |
| Luminosity coefficient `Ltilde_a` | **`e^{(2-q)nu}`** | FR2005 eq. (43) |

Because `eta^inf` is spatially uniform:

```text
L_H^inf = sum_a eta_a^inf R_a          (exactly FR2005 eq. 46)
```

**The `e^{-nu}` of `G_y` and the `e^{+nu}` of the reaction integral are different operations on
different integrands and may never be interchanged** (already stated in ADR-0013 §3.2).

### 3.6 Non-superfluid Urca imbalance functions

```text
F_D(xi) = 1 + 1071 xi^2/(457 pi^2) + 315 xi^4/(457 pi^4) + 21 xi^6/(457 pi^6)
H_D(xi) =      714 xi  /(457 pi^2) + 420 xi^3/(457 pi^4) + 42 xi^5/(457 pi^6)
F_M(xi) = 1 + 22020 xi^2/(11513 pi^2) + 5670 xi^4/(11513 pi^4)
            +   420 xi^6/(11513 pi^6) +    9 xi^8/(11513 pi^8)
H_M(xi) =     14680 xi  /(11513 pi^2) + 7560 xi^3/(11513 pi^4)
            +   840 xi^5/(11513 pi^6) +   24 xi^7/(11513 pi^8)
M_*(xi) = xi H_*(xi) - F_*(xi)
```

`F_*` multiplies the **neutrino emissivity**; `H_*` multiplies the **net reaction rate**. Authority:
FR2005 eqs. (34)–(37); R1995 eqs. (29)–(32) (with `u = xi/pi`); Y2020 eqs. (4.19)–(4.20).

**Erratum, normative.** FR2005 eq. (37) prints the last `H_M` denominator as `11513 pi^6`. **The
correct denominator is `11513 pi^8`.** Established by independent derivation from R1995's phase-space
integrals, by 40-digit quadrature, by FR2005's own eq. (60), and by Y2020 eq. (4.20)
(preflight §9.3). This is a typesetting error, not a source disagreement. Implementing `pi^6` is
mutation M23.

Required properties: `F_*(0)=1`; `H_*(0)=0`; `F_*` even; `H_*` odd; `xi H_*(xi) >= 0`; `F_* >= 1`.

Large-`xi`: `C_H = 24/(11513 pi^8)`, `C_M = 15/(11513 pi^8)`; heating fraction `1/2` (direct) and
`5/8` (modified), neutrino fraction `1/2` and `3/8`.

**Sign-crossing roots** (exact; durable test oracles), with the **root definition mandatory** at every
citation:

| Definition | Direct Urca | Modified Urca |
|---|---|---|
| **Incremental**, `M_*(xi) + 1 = 0` — equilibrium neutrino emission **is** subtracted, `L_H = DeltaL_nu` | `4.7870134733369` (`= pi sqrt((sqrt(93)-5)/2)`) | `4.90971002892413` |
| **Full**, `M_*(xi) = 0` — equilibrium neutrino emission **not** subtracted, `L_H = L_nu(T,eta)` | `5.4585315948676` | `5.63371746764834` |
| Maximum incremental cooling | `xi = 3.49729392477` | `xi = 3.61254068813` |

A bare number without its definition is non-compliant.

### 3.7 Thermal ledger and no-double-counting

```text
L_{nu,a}^inf     = Ltilde_a F_*(xi_a) T_inf^{q_a}
DeltaL_{nu,a}^inf = Ltilde_a [F_*(xi_a) - 1] T_inf^{q_a}        >= 0
L_{H,a}^inf      = Ltilde_a xi_a H_*(xi_a) T_inf^{q_a}          >= 0
DeltaP_{beta,a}  = L_{H,a}^inf - DeltaL_{nu,a}^inf
                 = Ltilde_a [M_*(xi_a) + 1] T_inf^{q_a}          (either sign)
```

**Architecture: INCREMENTAL, with a shared coefficient.** Two mandatory requirements:

- **(R1)** For each rotochemical channel the thermal RHS adds `+L_{H,a}` and `-DeltaL_{nu,a}`, never
  the full `L_nu(T,eta)`. The existing equilibrium `NeutrinoCooling` term is retained intact for all
  channels. At `eta = 0` both new terms are **identically zero** (structural, not tolerance-based),
  so the thermal RHS reduces bit-for-bit to today's passive cooling.
- **(R2)** The equilibrium term and every correction are built from **the same** `Ltilde_a`.
  `NeutrinoCoolingCachePayload` gains the channel-resolved set
  `{Ltilde_De, Ltilde_Dmu, Ltilde_Me, Ltilde_Mmu}`; the equilibrium driver consumes the sums.
  **A second, independently normalized coefficient set is forbidden** — it would still pass the
  `eta = 0` test while applying the enhancement ratio `F_*` to the wrong emissivity.

`Ltilde_a` is the **single normalization authority** for equilibrium cooling, `DeltaL_nu`, `L_H`, and
`R_l`, because `xi` is spatially constant (preflight §10.1).

Whether and when to replace the existing placeholder normalization constants
(`Q0_DU = 1e27`, `Q0_MU = 1e21`, `NeutrinoCooling_Details.cpp:101-103`) with a source-authoritative
`S_a(n)` is a **separate governed change with baseline consequences** and is **not authorized here**.

### 3.8 Thermal ODE in the evolved variable

The evolved thermal DOF is `x = ln(T_inf/T_ref)`, `T_ref = 1e8 K`. The rotochemical driver adds

```text
(dx/dt)_rotochem = sum_a Ltilde_a [ M_*(xi_a) + 1 ] T_inf^{q_a} / ( T_inf * C_*(T_inf) )
```

**ADR-0002 is unchanged and no change to it is proposed.** Its governing invariant already reads
`C_* dT_inf/dt = -L_nu - L_gamma + L_H + ...`. Substituting `L_nu,eq = sum_a Ltilde_a T_inf^{q_a}`
and `C_* = Ctilde T_inf` recovers FR2005 eq. (51) term for term.

### 3.9 Photon cooling — unchanged

Photon cooling is affected by rotochemical heating **only through `T_inf`**. No chemical-imbalance
correction to the envelope is introduced. Authority: FR2005 §3.3 applies the Potekhin et al. (1997)
envelope as a function of `T_inf` and surface gravity alone, and the chemical imbalance is a core
quantity while `B_ij` is integrated over the core only (FR2005 §3.5). Recorded as a **v1 contract
assumption**; any future envelope `eta` dependence requires its own source authority and its own ADR.

### 3.10 Direct- and modified-Urca channel support

Four declared channels `{De, Dmu, Me, Mmu}`, `q = 6` for direct and `q = 8` for modified. Neutron and
proton branches are summed **inside** `Ltilde_{M,l}`; they are not separate evolution channels
(FR2005 §3.2). Electron and muon channels **must** stay separate: the chemical ODE evaluates `H_*` at
`xi_npe` and `xi_npmu` independently.

Direct-Urca support is owned by `StarContext::DirectUrcaLastAllowedIndex()`
(`kF_n <= kF_p + kF_e`, composition-source authoritative, contiguous core region). **Phase-5D reuses
it unchanged and introduces no competing threshold policy.** The absent muon support domain must be
added as a *second instance of the same criterion* with `kF_mu` substituted. For frozen v1
coefficients the support domains are static; there is no moving-threshold term, and `eta != 0` does
not shift the kinematic threshold.

### 3.11 Frozen coefficients (v1)

`Z`, `W`, `I_Omega` and `Ltilde` are **frozen**, evaluated on the nonrotating background star, and
byte-identical at every RHS evaluation of a run. Authority: FR2005 §2.1 (*"the `B_ij` do not depend
on time"*) and eq. (19); R2006 §3 (*"a constant, depending on the structure of the unperturbed,
nonrotating star"*).

If `Z` were ever allowed to vary, the ODE acquires an **extra term**:

```text
etadot^inf = + Zdot Z^-1 eta^inf - Z R + 2 W Omega Omegadot
```

(note the `+` sign, from `delta N_l = -Z^-1 eta^inf`). **Not implemented and not authorized.**
Enabling time-dependent coefficients without this term is mutation M18 and must refuse.

### 3.12 Spin ownership

The rotochemical module consumes `Omega(t)` and `Omegadot(t)` through an explicit read-only
spin-history interface supplied by the `DriverContext`. It **must not** re-derive a torque and
**must not** depend on driver registration order to read the spin RHS: `EvolutionSystem::operator()`
executes drivers in registration order with no dependency sort. A prescribed-history implementation
(analytic or observed pulsar timing) and a state-coupled adapter are both required. Magnetic-dipole
braking is an **example driver**, never a rotochemical dependency. `RotochemicalSpinDrive::Evaluate`
already returns `2 W Omega Omegadot` per channel and is consumed unchanged.

If the accumulator-reading route is ever adopted instead, a driver dependency sort becomes
**mandatory, not optional**.

### 3.13 Dependency direction

```text
UrcaImbalanceFunctions        (pure mathematics; knows nothing about stars or eta)
        |
UrcaChannelLuminosityCoefficients   (microphysical S_a + stellar integration + support domains)
        |
RotochemicalReactionResponse / RotochemicalThermalPower   (global rates and powers)
        |
RotochemicalEvolutionRHS      (secular coupling; consumes Z, W, ISpinHistory)
```

No object may own more than one of these layers. The rotochemical module depends on
`CompactStar/Analysis` (for `Z`, `W`); the generic evolution core must **not** acquire a dependency
on rotochemical physics.

### 3.14 Provenance and lifetime

Every evolution result retains: the `ChemicalImbalanceResponse` and `RotochemicalSpinDrive` identity
and revision; the `Ltilde` coefficient identity and its normalization classification (source-authority
vs declared benchmark input); the `StarProfile` and geometry identity; the spin-history identity; the
channel order; the solver identity and tolerances; and initial conditions. A changed dependency
**refuses before scientific access**. This inherits ADR-0013 §6 verbatim and adds nothing weaker.

### 3.15 First implementation scope

Objects (1)–(8) of preflight §33.1. Reused unchanged: `ThermalState`, `C_*(T_inf)`, `PhotonCooling`
and envelopes, `SpinState`/`MagneticDipole`, the evolution core and `GSLIntegrator` with `RKF45`,
`ChemState` storage, the `k_B` authority, the direct-Urca support logic, and the Phase-5C `Z`/`W`
objects.

**First benchmark:** the governed Track-R free-gas Structure-1 star
(`rho_c = 1.10e15 g/cm^3`, `M = 0.6236 M_sun`), for which **direct Urca is closed everywhere**
(`x_p <= 0.0157` against the `1/9` threshold) and **muons are present** in the core — i.e. a pure
modified-Urca, two-lepton-channel run, exactly FR2005's analytically tractable case. Run
`T_inf(0) = 1e8 K`, `eta(0) = 0`, prescribed dipole spin history `B = 1e8 G`, `P_0 = 1 ms`, to
`1e10 yr`. Observables B1–B10 of preflight §25.5.

**The free-gas benchmark validates architecture, signs, energy bookkeeping and ODE correctness. It is
NOT an FR2005 source reproduction and no claim about realistic nuclear matter follows from it.**

### 3.16 Numerical contract

Evolve `eta^inf` in MeV, **not** `xi`. Keep `RKF45` for v1: the free-gas system is moderately stiff
but tractable explicitly (measured 4829–6420 steps and ~33–38k RHS evaluations to `1e10 yr`, against
561 steps for an auto-switching stiff method — a 20–30× cost, not a failure). Per-component absolute
tolerances are recommended (`~1e-12` for `ln T_inf`, `~1e-18 MeV` for `eta`) in place of the current
scalar `atol = 1e-10` applied across a heterogeneous state vector. **No positivity clamp may be
applied to `eta`** — it legitimately changes sign for a spun-up star. **The solver decision must be
revisited before direct Urca or superfluidity is enabled**, where two-timescale behaviour is expected
(FR2005 Fig. 6).

### 3.17 Future BNV seam

`G_y` **must remain available in the unreduced source basis** `(N_n, N_e, N_mu)`. A
baryon-number-violating source cannot be represented in the fixed-baryon two-channel space, because
the lift `L = [[-1,-1],[1,0],[0,1]]` encodes `delta N_n = -(delta N_e + delta N_mu)`. The seam is a
generic external source `Sigma_y(t, state)` on `delta N_y`; the reduced path `Sigma_y = L Sigma_l`
applies only at fixed baryon number. JRF2006 is the source precedent that the formalism accepts a
different external parameter (`Gdot` in place of `2 Omega Omegadot`) — **but its source is also at
fixed `A`, so it authorizes only the fixed-baryon seam.**

**No BNV rate is derived, no BNV heating is claimed, and no BNV code is authorized.**

---

## 4. Exclusions

Superfluidity; A18 implementation; BNV implementation; time-dependent `Z`/`W`/background; replacement
of the placeholder equilibrium Urca normalization; envelope `eta` dependence; crustal processes;
hyperons; solver replacement; and any semantic change to ADR-0002, ADR-0004, ADR-0010, ADR-0011,
ADR-0012 or ADR-0013.

---

## 5. Validation ladder

`RE1`–`RE18` as specified in `docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md`
§31, classified INDEPENDENT / SOURCE TRACEABILITY / CONTRACT / CONVERGENCE / SOURCE-LIMITED.
**RE1–RE14 and RE16–RE18 are achievable with the free-gas fixture; RE15 (FR2005 source benchmark) is
SOURCE-LIMITED and blocked.** The mutation inventory `M1`–`M25`, with its explicit list of algebraic
aliases that must not be counted as independent coverage, is preflight §32.

---

## 6. Unresolved realistic-source blockers

**(A) First free-gas evolution.** One genuine source gap: the **muon-channel equilibrium Urca
normalization**. R1995 eq. (34) supplies primary authority for the electron modified-Urca channel;
FR2005 §3.2 delegates the branch- and lepton-resolved normalization to Yakovlev et al. (2001), which
is not in the shared library. **Resolution:** the v1 free-gas normalization is a **declared benchmark
input** with recorded provenance, explicitly classified as *not an FR2005 reproduction*. Every
benchmark oracle except the absolute quasi-steady value holds for any positive `Ltilde`. The
remaining items (channel-resolved `Ltilde`, muon DU support domain, `n_eta = 0`, spin-history
interface, Phase-5C canonical integration) are implementation or process, not source.

**(B) Realistic FR2005 reproduction.** Blocked on: APR (1998) arbitrary-composition `E(n_B, x_p)`
for A18+δv+UIX* (**not in the library and not inspected in this task**); Yakovlev et al. (2001) rate
normalization; Page et al. (2004) effective masses; PRL (1995) inner crust; Haensel & Pichon (1994)
outer crust; and R2006 Fig. 1 / FR2005 Fig. 3–4–6 numerical arrays or governed digitization. A public
beta-equilibrium table is **necessary but not sufficient** — it constrains only a one-parameter curve
and cannot supply the transverse second derivatives `dn_i/dmu_j`.

**(A) and (B) are disjoint and must never be conflated.**

---

## 7. INV-11 disposition

| Subpart | Would this ADR resolve it? |
|---|---|
| INV-11a coefficient redshift semantics | **Completed** (was partial under ADR-0013 Q7) |
| INV-11b evolved `eta` state ownership | **Yes** |
| INV-11c reaction-rate sign / index convention | **Yes** |
| INV-11d thermal energy ledger | **Yes** |
| INV-11e coefficient lifetime / update policy | **Frozen v1 only**; the time-dependent branch is recorded and forbidden, not resolved |
| INV-11f ODE / source coupling | **Partially** — contract fixed; the stiff-solver decision is deferred |

**INV-11 is NOT marked resolved.** This ADR is PROPOSED, nothing is implemented, and INV-11e's
time-dependent branch and INV-11f's solver branch remain open by design. INV-11 stays UNRESOLVED and
fail-closed until this ADR is accepted **and** its implementation is validated.

---

## 8. Consequences

**If accepted.** A future implementer has no remaining freedom on any sign, redshift factor,
channel ordering, unit, ownership boundary, or double-counting question. The reusable production
objects, their dependency direction, the validation ladder and the mutation inventory are fixed. The
free-gas benchmark is executable without A18.

**If rejected or amended.** No code changes, because none exists. Only this document and the
preflight record change.

**Costs.** `NeutrinoCoolingCachePayload` must gain channel resolution, touching a currently-passing
driver. `EvolutionConfig::n_eta` must become non-zero in the benchmark path. A spin-history interface
must be added. None of these is authorized by this proposal.

---

## 9. Acceptance requirements

This ADR may be accepted only after an independent scientific review verifies, at minimum:

1. the reaction-rate sign convention and the `eta DeltaGamma >= 0` dissipation claim (§3.3);
2. every global GR redshift factor in §3.5, derived rather than checked against this document;
3. the four imbalance polynomials and the `pi^8` erratum (§3.6);
4. the four sign-crossing roots **and their definitions** (§3.6);
5. the rate-normalization identity and the single-`Ltilde` requirement (§3.7);
6. that the incremental architecture cannot double-count the existing equilibrium cooling, including
   the R2-without-R1 and R1-without-R2 failure modes;
7. the frozen-coefficient semantics and the `+Zdot Z^-1 eta` term (§3.11);
8. the free-gas benchmark contract, including the claim that direct Urca is closed and muons are
   present (§3.15).

**Do not begin BNV before the standard rotochemical evolution machinery is validated.**

---

## 10. Status

**PROPOSED. NOT ACCEPTED. NOT IMPLEMENTED.** Acceptance is not implementation validation.
