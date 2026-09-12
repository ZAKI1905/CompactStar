# ADR-0014 — secular rotochemical evolution contract

**Status:** ACCEPTED / HUMAN-RATIFIED / CANONICALLY INTEGRATED
**Decision:** ACCEPTED WITH RETAINED FINAL INDEPENDENT-REVIEW CAVEATS
**Date:** 2026-09-07
**Starting canonical SHA:** `49ab2b8c2881b6ef7b9309307d18cea51d557f72`
**Branched from (Phase-5C human-ratified candidate):** `27727016856a6a25a46e447c70e380722ea8ddbf`
**Change class:** scientific-semantic and architectural contract; documentation-only ratification.
**Evidence companion:** `docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md`
**Implementation state:** no production evolution object, weak rate, imbalance function, chemical or
thermal RHS term, superfluid treatment, A18 model, BNV source, test, baseline, EOS/data, or
literature byte is created or modified by this ratification. Phase-5C canonical integration is now
complete; controlled non-superfluid v1 implementation is authorized only in a separate governed
production task.

> **PHASE-5D SECULAR ROTOCHEMICAL EVOLUTION CONTRACT — SCIENTIFIC PREFLIGHT COMPLETE /
> INDEPENDENTLY REVIEWED / HUMAN-RATIFIED FOR CONTROLLED NON-SUPERFLUID V1 SCOPE —
> PRODUCTION IMPLEMENTATION NOT YET BEGUN.**

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

## 3. Accepted contract

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
xi_a = eta_a^inf / (k_B^(MeV) T_inf) = eta_a(r) / (k_B^(MeV) T(r))  — spatially constant
k_B^(MeV) = Zaki::Physics::K_BOLTZ_EV * 1e-6                         MeV/K
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
DeltaGamma_a = (1/(k_B^(erg) T)) Q_a^eq H_*(xi_a)  H_* odd, H_*(0) = 0
eta_a DeltaGamma_a >= 0   pointwise, both signs of eta       (thermodynamic dissipation)
V = eta^T Z^-1 eta  is nonincreasing for the reaction sub-system
```

For `Omegadot = 0` with frozen symmetric positive-definite `Z`,

```text
dV/dt = -2 sum_l eta_l R_l <= 0.
```

Strict decrease is claimed only for a nonzero imbalance direction for which at least one
physically active reaction channel has positive normalization and nonzero dissipative response.
If every applicable channel coupled to an imbalance direction has zero normalization
(`sum_{a in l} Ltilde_a = 0`), an **uncoupled dead imbalance** may freeze; with nonzero `Z`
cross-coupling, a dead individual reaction channel does not imply that its `eta` component remains
at its initial value. The Lyapunov result is negative-semidefinite/nonincreasing, not globally
strict.

**Warning, normative:** the literature contains both signs. Y2020 fn. 3, p. 61, and R1995 fn. 3,
p. 14, each record a sign opposite to cited references. Any future rate source must be re-anchored to
this table before use.

### 3.4 Chemical ODE

```text
etadot_npe^inf  = - Z_npe R_e - Z_np    R_mu + 2 W_npe  Omega Omegadot
etadot_npmu^inf = - Z_np  R_e - Z_npmu  R_mu + 2 W_npmu Omega Omegadot

R_l = integral_D dV e^{Phi} sum_{a in l} DeltaGamma_a
    = sum_{a in l} [Ltilde_a / k_B^(erg)] H_a(xi_l) T_inf^(q_a-1)     [count/s]

W = Z I_Omega ,     I_{Omega,l} = (dN_l^eq/dOmega^2)_A
```

`R_l` is the **global, `e^{+Phi}`-weighted, channel-summed** net lepton-creation rate — not a local
density. `D` is the declared chemical/reaction domain defined in §3.10 and must be consistent with
the associated `G_y` and channel-support contract; no implicit core boundary is permitted. The
ADR-0013 §3.4 schematic `etadot^inf = -Z R + 2 W Omega Omegadot` is thereby given its exact content.
Authority: FR2005 eqs. (14)–(15), (30), (45), (52)–(53); R2006 eqs. (5)–(7), (15), (19); Y2020
eqs. (4.8)–(4.12).

The canonical normalization and unit boundary are:

```text
Ltilde_a                         [erg s^-1 K^-q_a]
k_B^(erg)                       [erg K^-1]
[Ltilde_a/k_B^(erg)] T_inf^(q_a-1)
  = (erg s^-1 K^-q)/(erg K^-1) K^(q-1) = s^-1 == count/s
Z R                             = (MeV/count)(count/s) = MeV/s.
```

`k_B^(erg)` must be derived from the repository's single Boltzmann authority
`Zaki::Physics::K_BOLTZ_EV` and the governed energy-unit conversion:
`C_(MeV->erg) = Units::MEV_FM3_TO_ERG_CM3 / 10^39` and
`k_B^(erg) = C_(MeV->erg) [Zaki::Physics::K_BOLTZ_EV * 10^-6]`. It must not be a new literal.
The MeV/K view remains the one used in `xi`; an MeV/K value may not divide an erg-normalized
`Ltilde`.

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

Because `eta^inf` is spatially uniform, define the canonical chemical-module power and its thermal
boundary conversion separately:

```text
P_H,chem^inf [MeV/s] = sum_a eta_a^inf[MeV] R_a[count/s]
L_H^inf [erg/s]      = C_(MeV->erg) P_H,chem^inf
                      = Ltilde_a xi_a H_*(xi_a) T_inf^q summed over a.
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

**CONFIRMED PRINTED TYPO / INTERNAL SOURCE INCONSISTENCY; NO PUBLISHED ERRATUM LOCATED.** FR2005
eq. (37) literally prints the last `H_M` denominator as `11513 pi^6`. **The proposed normative
implementation formula uses `11513 pi^8`.** The evidence is the literal FR2005 printing; the later
FR2005 coefficient structure in eq. (60), which requires `pi^8`; the R1995 source-variable form;
the independent phase-space/Fermi-convolution derivation and quadrature; and Y2020 eq. (4.20), which
prints `pi^8` (preflight §9.3). A bounded search located no published erratum. Implementing `pi^6`
is mutation M23.

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
P_{H,chem,a}^inf = eta_a^inf R_a                               [MeV/s]
L_{H,a}^inf      = C_(MeV->erg) P_{H,chem,a}^inf
                 = Ltilde_a xi_a H_*(xi_a) T_inf^{q_a}          [erg/s] >= 0
DeltaP_{beta,a}  = L_{H,a}^inf - DeltaL_{nu,a}^inf
                 = Ltilde_a [M_*(xi_a) + 1] T_inf^{q_a}          (either sign)
```

**Architecture: INCREMENTAL, with a shared coefficient.** `Ltilde_a` is stored/exposed in
`erg s^-1 K^-q_a`, matching the existing thermal luminosity interface. Three mandatory
requirements:

- **(R1)** For each rotochemical channel the thermal RHS adds `+L_{H,a}` and `-DeltaL_{nu,a}`, never
  the full `L_nu(T,eta)`. At `eta = 0`, `H(0)=0` and `F(0)=1` make `R`, `DeltaL_nu`, and `L_H`
  identically zero, and the extended channel reduces exactly to **its own declared** equilibrium
  luminosity `L_nu,eq = Ltilde_a T_inf^q`. This is a same-coefficient identity, not a requirement to
  reproduce the historical placeholder cooling baseline bit-for-bit.
- **(R2)** The equilibrium term and every correction are built from **the same** `Ltilde_a`.
  `NeutrinoCoolingCachePayload` gains the channel-resolved set
  `{Ltilde_De, Ltilde_Dmu, Ltilde_Me, Ltilde_Mmu}`; the equilibrium driver consumes the sums.
  **A second, independently normalized coefficient set is forbidden** — it would still pass the
  `eta = 0` test while applying the enhancement ratio `F_*` to the wrong emissivity.
- **(R3)** The canonical chemical power is formed internally as
  `P_H,chem^inf = sum_l eta_l^inf[MeV] R_l[count/s]` in `MeV/s` (with count dimensionless). It is
  converted exactly once through the
  governed energy-unit authority when contributed at the thermal-RHS luminosity boundary in
  `erg/s`. Nonequilibrium neutrino luminosity remains `erg/s`. An erg-native rearrangement is not a
  second public ownership route.

`Ltilde_a` is the **single normalization authority** for equilibrium cooling, nonequilibrium
`L_nu`/`DeltaL_nu`, `R_l`, and chemical-heating bookkeeping, because `xi` is spatially constant
(preflight §10.1). The dimensional equivalence of the MeV-native power and the displayed erg-native
formula follows from `k_B^(erg) = C_(MeV->erg) k_B^(MeV)` and
`eta^inf = xi k_B^(MeV) T_inf`: converting `eta R` once gives
`C_(MeV->erg) eta R = Ltilde xi H T_inf^q`.

Whether and when to replace the existing placeholder normalization constants
(`Q0_DU = 1e27`, `Q0_MU = 1e21`, `NeutrinoCooling_Details.cpp:101-103`) with a source-authoritative
`S_a(n)` is a **separate governed change with baseline consequences** and is **not authorized here**.
The existing historical placeholder `NeutrinoCooling` normalization is not the controlled
benchmark's equilibrium coefficient. The controlled benchmark instantiates its equilibrium Urca
contribution from the **same declared benchmark `Ltilde_a`** used for `F_*`, `H_*`, reaction rates,
and heating. Source-authoritative replacement of the **default production normalization** is a
later realistic-physics task; this documentation contract changes no historical baseline.

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
envelope as a function of `T_inf` and surface gravity alone. Chemical/reaction integration uses the
explicitly declared domain `D`; no envelope contribution is implied. Recorded as a **v1 contract
assumption**; any future envelope `eta` dependence requires its own source authority and its own ADR.

### 3.10 Direct- and modified-Urca channel support

Four declared channels `{De, Dmu, Me, Mmu}`, `q = 6` for direct and `q = 8` for modified. Neutron and
proton branches are summed **inside** `Ltilde_{M,l}`; they are not separate evolution channels
(FR2005 §3.2). Electron and muon channels **must** stay separate: the chemical ODE evaluates `H_*` at
`xi_npe` and `xi_npmu` independently.

`D` is the chemical/reaction domain used by the associated `G_y`. For the controlled free-gas
fixture it is exactly the connected whole-star, source-valid Phase-5C chemical domain: centre,
authenticated `npemu`/`npe`/`pe` branch partition, and the governed physical neutron-onset/vacuum
boundary with its accepted refusal-window and tail treatment. No saturation-density, crust, or
other arbitrary cutoff is introduced. Each channel has an explicit support/applicability subset
`D_a subseteq D`.

The controlled free-gas Phase-5D architecture benchmark is intentionally configured
**MODIFIED-URCA-ONLY**. Its static declared enabled-process set is `{Me, Mmu}` and excludes `{De,
Dmu}`. Therefore, for this benchmark, `D_De = empty` and `D_Dmu = empty` by the benchmark
process-configuration contract, not by inference from a triangle condition, a degeneracy threshold,
or `nB_min`.

Triangle support, declared benchmark process support, and general future physical DU applicability
are three distinct concepts. In the authenticated fixture, the muon-DU triangle does not open where
muons are present. The electron-DU triangle does open in a very-low-density outer sliver near neutron
onset, approximately `n_B = 7.36e-9 ... 6.67e-8 fm^-3`, where `E_F,n` is only of order keV. That
sliver is outside the intended source-validity regime used to motivate the controlled benchmark at
early/high-temperature epochs, but its degeneracy status changes with temperature and is **not** the
static exclusion rule for the complete evolution. The existing `nB_min` remains only a
numerical/semantic guard: it is not a physical DU threshold, a degeneracy threshold, or the reason
`D_De` is empty.

Future realistic or source-authoritative evolution must obtain DU support from a governed
applicability provider capable of representing kinematics, composition, degeneracy/model-validity
requirements, and possibly disconnected support intervals. The controlled benchmark's disabled-DU
configuration is not a universal physical rule and does not claim that DU is physically impossible
everywhere.

Future support representation must be explicit rather than only a last index. The low-density
triangle-open electron sliver is a required **negative control**: a future implementation must show
that triangle opening does not activate a DU channel when the process is disabled or outside its
declared channel-support contract. The current
`BuildDirectUrcaMaskCache_` scans for an allowed region while the luminosity integral consumes
`[0,last]`; a future outer allowed shell could therefore sweep a closed inner region. This is
nonblocking for the controlled fixture but remains a required disconnected-support negative
control, together with an explicit innermost-first profile-order assertion. The absent muon
criterion remains a second instance with `kF_mu` substituted. Frozen v1 benchmark process support
is static; no temperature-dependent support predicate is introduced, and `eta != 0` does not move
the equilibrium-background kinematic threshold.

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
UrcaChannelMicrophysics / UrcaChannelNormalization
        |  local S_a(r), process/lepton/branch, effective masses, matrix elements,
        |  alpha/beta factors, applicability metadata, source provenance
        |  S_a(r) [erg cm^-3 s^-1 K^-q_a]
        v
GlobalUrcaChannelCoefficient
        |  stellar GR integration over declared D_a only -> Ltilde_a
        +-------------------------------+
                                        v
UrcaImbalanceFunctions ----------------> RotochemicalReactionResponse
  (pure dimensionless F/H only)           (Ltilde + F/H + T_inf + eta_inf
                                           -> R, DeltaL_nu, L_H)
                                        |
                                        v
Secular RHS                    (combines channels with thermal and spin sources)
```

The separation between local/source microphysics normalization and stellar GR integration is
normative; one object may not own both. No object may own more than one of these layers. The
rotochemical module depends on
`CompactStar/Analysis` (for `Z`, `W`); the generic evolution core must **not** acquire a dependency
on rotochemical physics.

### 3.14 Provenance and lifetime

Every evolution result retains: the `ChemicalImbalanceResponse` and `RotochemicalSpinDrive` identity
and revision; the `Ltilde` coefficient identity and its normalization classification (source-authority
vs declared benchmark input); the chemical/reaction domain `D` and channel support/applicability
domain `D_a`; the `StarProfile` and geometry identity; the spin-history identity; the channel order;
the solver identity and tolerances; and initial conditions. A changed dependency
**refuses before scientific access**. This inherits ADR-0013 §6 verbatim and adds nothing weaker.

### 3.15 First implementation scope

Objects (1)–(9) of preflight §33.1. Reused unchanged: `ThermalState`, `C_*(T_inf)`, `PhotonCooling`
and envelopes, `SpinState`/`MagneticDipole`, the evolution core and `GSLIntegrator` with `RKF45`,
`ChemState` storage, the `k_B` authority, and the Phase-5C `Z`/`W` objects. For the current fixture
only, the static MODIFIED-URCA-ONLY benchmark configuration is reused. Before any enabled, nonempty,
or outer-shell DU adapter, the last-index representation must be replaced by explicit support and
the ordering/applicability negative controls of §3.10 must pass.

**First benchmark:** the governed Track-R free-gas Structure-1 star
(`rho_c = 1.10e15 g/cm^3`, `M = 0.6236 M_sun`), whose static enabled-process set is `{Me, Mmu}` and
therefore has `D_De = D_Dmu = empty` as specified in §3.10; **muons are present** in the interior.
This is a modified-Urca, two-lepton-channel run. Run
`T_inf(0) = 1e8 K`, `eta(0) = 0`, prescribed dipole spin history `B = 1e8 G`, `P_0 = 1 ms`, to
`1e10 yr`. Observables B1–B10 of preflight §25.5.

Its positive electron and muon normalizations are declared benchmark coefficients, not realistic
source normalizations, and the same `Ltilde` values must feed equilibrium cooling, the `F`
correction, the `H` reaction rate, and chemical heating. Electron and muon channels remain distinct;
no realistic astrophysical normalization claim is made. **This mathematical/architecture benchmark
validates ODE wiring, signs, energy bookkeeping, quasi-steady scaling, spin coupling, and thermal
coupling. Because it consumes declared `Ltilde` values, it does not by itself validate the stellar
construction of `Ltilde` or the `GlobalUrcaChannelCoefficient` GR integrand. That layer is validated
separately by RE10b. It cannot validate the FR2005 absolute temperature/history.**

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

`RE1`–`RE18` plus `RE10b` as specified in
`docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md` §31, classified INDEPENDENT /
SOURCE TRACEABILITY / CONTRACT / CONVERGENCE / SOURCE-LIMITED. RE10b is the dedicated **INDEPENDENT
ANALYTIC / NUMERICAL ORACLE** for the `GlobalUrcaChannelCoefficient` GR integrand and does not use an
injected/precomputed `Ltilde` as its expected value. **RE1–RE14 (including RE10b) and RE16–RE18 are
achievable with controlled fixtures; RE15 (FR2005 source benchmark) is SOURCE-LIMITED and blocked.**
The mutation inventory `M1`–`M37`, with its explicit list of algebraic
aliases that must not be counted as independent coverage, is preflight §32.

---

## 6. Unresolved realistic-source blockers

**(A) First free-gas evolution.** R1995 eq. (34) is historical/supporting normalization: it is
npe-lumped, the source itself calls it *"somewhat uncertain"*, and it is not sufficient as the
definitive FR2005 realistic normalization. The v1 free-gas run may instead use positive declared
benchmark coefficients with recorded provenance, explicitly classified as *not realistic source
normalizations* and used consistently in all four channel roles. Yakovlev et al. (2001), Phys. Rep.
354, 1, arXiv:astro-ph/0012122, is publicly available and sufficient **in form** to define the
channel-normalization architecture: modified-Urca neutron normalization, proton-branch factor,
electron/muon substitutions, threshold structure, and effective-mass dependence. It is not installed
or authenticated in the library here. The exact `alpha_n` prescription matching the intended
FR2005 reproduction remains unresolved (for example, an FM79-style constant versus later density-
dependent/OPE choices); this ADR chooses none.

**(B) Realistic FR2005 reproduction.** Bounded primary-source preprint inspection confirms that the
APR 1998 publication provides enough **form** to reconstruct the A18+δv+UIX* arbitrary-composition
core functional `E(n_B,x_p)`: effective Hamiltonian/effective masses, proton-fraction interpolation,
Appendix-A fit functions, and Table-XII `p1...p21` for LDP and HDP. This is discovery/feasibility
evidence only; no source byte is installed and no A18 implementation begins. Realistic closure
remains blocked on authenticated journal-version authority; exact FR2005 Maxwell/phase construction;
crust joins; authenticated YKGH2001 normalization; effective-mass interpretation; `alpha_n` choice;
direct/MU support authority; Page/crust/envelope authority as required; and published benchmark
arrays or governed digitization. A beta-equilibrium table remains insufficient for transverse
second derivatives.

**(A) and (B) are disjoint and must never be conflated.**

---

## 7. INV-11 disposition

| Subpart | Ratified disposition |
|---|---|
| INV-11a redshift / coefficient semantics | **PARTIALLY RESOLVED UPSTREAM; ACCEPTED EXTENSION / CLARIFICATION HERE.** ADR-0013 resolves coefficient-object semantics for `G_y`, `Z`, and `W`; this ADR accepts their secular extension to evolved `eta`, `xi`, global reaction-rate and neutrino-luminosity integrals, and global heating coupling. |
| INV-11b evolved `eta` state ownership | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| INV-11c reaction-rate sign / index convention | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| INV-11d thermal energy ledger / no double counting | **CONTRACT RESOLVED / IMPLEMENTATION PENDING** |
| INV-11e frozen coefficient lifetime / update policy | **CONTRACT RESOLVED / IMPLEMENTATION PENDING**; time-dependent branch recorded and forbidden in v1 |
| INV-11f ODE / source coupling | **UNRESOLVED / IMPLEMENTATION + VALIDATION PENDING** |

**INV-11 is NOT marked resolved.** The accepted contract resolves INV-11b–e at contract level only;
nothing is implemented, and INV-11e's time-dependent branch and INV-11f's implementation/validation
branch remain open by design. Global INV-11 stays UNRESOLVED.

---

## 8. Consequences

**Acceptance effect.** A future implementer has no remaining freedom on any sign, redshift factor,
channel ordering, unit, ownership boundary, or double-counting question. The reusable production
objects, their dependency direction, the validation ladder and the mutation inventory are fixed. The
free-gas benchmark is executable without A18.

**If rejected or amended.** No code changes, because none exists. Only this document and the
preflight record change.

**Costs.** `NeutrinoCoolingCachePayload` must gain channel resolution, touching a currently-passing
driver. `EvolutionConfig::n_eta` must become non-zero in the benchmark path. A spin-history interface
must be added. None of these is implemented by this ratification; implementation remains a
separately governed task after the accepted Phase-5C coefficient history is canonically integrated.

---

## 9. Satisfied acceptance requirements

The final bounded independent scientific review and human-owner disposition satisfy the requirement
that review verify, at minimum:

1. the reaction-rate sign convention, unit boundary, and `eta DeltaGamma >= 0` dissipation claim
   with the qualified nonincrease/strictness statement (§3.3–§3.4);
2. every global GR redshift factor in §3.5, derived rather than checked against this document;
   RE10b must independently exercise the `GlobalUrcaChannelCoefficient` integrand rather than an
   injected `Ltilde`;
3. the four imbalance polynomials and the proposed `pi^8` correction, without claiming a published
   erratum (§3.6);
4. the four sign-crossing roots **and their definitions** (§3.6);
5. the rate-normalization identity and the single-`Ltilde` requirement (§3.7);
6. the same-coefficient equilibrium/no-double-counting identity, without using the historical
   placeholder baseline as an oracle;
7. the frozen-coefficient semantics and the `+Zdot Z^-1 eta` term (§3.11);
8. the free-gas benchmark contract, including its static MODIFIED-URCA-ONLY enabled-process set,
   the temperature-dependent degeneracy status of the low-density electron-DU triangle sliver,
   empty benchmark DU domains, explicit support representation hazards, and muon presence (§3.10,
   §3.15).

**Do not begin BNV before the standard rotochemical evolution machinery is validated.**

---

## 10. Status

**ACCEPTED / HUMAN-RATIFIED. NOT IMPLEMENTED.** Acceptance is not implementation validation.

---

## 11. Post-independent-review revision — 2026-09-07

The independent Phase-5D-0R Opus adjudication returned disposition C with zero blocking scientific
derivation failures and material documentation/contract findings M-1 through M-8. This revision
changes only the affected contract text. It does not reopen or alter the independently confirmed
signs, redshifts, imbalance functions, roots, frozen-coefficient equation, spin ownership, thermal
variable, or unreduced `G_y` seam.

| Original proposal wording | Review finding | Corrected proposal wording |
|---|---|---|
| Rate formula mixed MeV/K `k_B` with an unspecified luminosity energy unit | M-1 | `Ltilde` is erg-based; `R` uses derived `k_B^(erg)`; chemical power is MeV/s and converts once at the thermal boundary |
| DU described as closed everywhere | M-2 | Kinematic low-density electron sliver is recorded; applicable degenerate-DU support is empty under the eligibility contract |
| Global strict Lyapunov claim | M-3 | `V` is nonincreasing; strictness requires an active dissipative channel in each nonzero direction |
| RE9 implied equality to the historical passive-cooling baseline | M-4 | RE9 is the exact same-coefficient channel identity only |
| Reaction integral used implicit `core` notation | M-5 | All normative global coefficient/rate integrals use declared `D`/`D_a` consistent with `G_y` and support |
| Printed `H_M` issue called an erratum | M-6 | Confirmed printed typo/internal inconsistency; no published erratum located; `pi^8` remains proposed normative formula |
| One proposed object owned microphysics and stellar integration | M-7 | Local normalization and global GR integration are separate normative layers |
| Source/blocker ledger understated YKGH2001/APR discovery and overstated R1995 | M-8 | R1995 qualified; YKGH2001 sufficient in form with `alpha_n` unresolved; APR reconstruction feasibility confirmed in form, realistic closure still blocked |

At `PHASE5D0_REVISION_SHA`, this table recorded all eight findings as **CLOSED BY TEXT REVISION**.
Phase-5D-0RR subsequently found M-2 only partially closed; that residual and the two additional
material findings are closed by the bounded revision below. Nothing becomes accepted until final
bounded independent re-review and human-owner ratification.

---

## 12. Post-Phase-5D-0RR material-closure revision — 2026-09-08

The independent Phase-5D-0RR re-review returned zero blocking findings and three material findings.
This bounded revision changes only their contract text and the directly adjacent clarifications
identified by that review. It does not implement the future validation oracle and does not reopen
the independently confirmed Phase-5D physics.

| Item | Review defect | Corrected contract | Remaining future implementation obligation | Status |
|---|---|---|---|---|
| R1 — DU-support applicability | The low-density electron-DU triangle sliver was called non-degenerate unconditionally, making a temperature-dependent observation the purported static exclusion rule. | The controlled benchmark is statically MODIFIED-URCA-ONLY: `{Me, Mmu}` enabled, `{De, Dmu}` disabled, so `D_De = D_Dmu = empty`. Triangle support, benchmark process support, and future physical applicability are distinct; the sliver remains a temperature-aware factual observation and negative control, not the exclusion rule. | Implement a governed future applicability provider with kinematics, composition, degeneracy/model-validity, disconnected-domain support, and the sliver/outer-shell negative controls before any realistic DU evolution. | **CLOSED BY TEXT REVISION** |
| R2 — global-`Ltilde` GR-factor detector | The declared-`Ltilde` architecture benchmark and its cited gates could not independently detect a wrong lapse power in stellar coefficient construction. | RE10b directly exercises `GlobalUrcaChannelCoefficient` against an independent analytic/high-precision quadrature value for `integral_{D_a} 4 pi r^2 e^lambda S_a(r) e^{(2-q_a)nu(r)} dr`, including declared unit prefactors and lapse, proper-volume, exponent, and domain mutants. | Implement RE10b independently of the production integration kernel; the present revision specifies but does not implement it. | **CLOSED BY TEXT REVISION** |
| R3 — INV-11a governance status | INV-11a was incorrectly stated as already resolved upstream. | INV-11a is **PARTIALLY RESOLVED UPSTREAM; PROPOSED COMPLETION / EXTENSION HERE**. ADR-0013 governs only static `G_y`/`Z`/`W` coefficient-object semantics; ADR-0014 proposes the secular-evolution extension. Global INV-11 remains UNRESOLVED. | Final bounded review and human-owner ratification remain required; no global invariant closure or implementation is implied. | **CLOSED BY TEXT REVISION** |

ADR-0014 remains **PROPOSED. NOT ACCEPTED. NOT IMPLEMENTED.**

---

## 13. Human-owner ratification — 2026-09-09

The bounded final independent re-review returned disposition B:

> **PHASE-5D FINAL BOUNDED RE-REVIEW PASS WITH NONBLOCKING FINDINGS — R1/R2/R3 CLOSED —
> ADR-0014 READY FOR HUMAN-OWNER RATIFICATION WITH EXPLICIT CAVEATS.**

Review totals were **0 BLOCKING, 0 MATERIAL, 12 NONBLOCKING, and 7 NOTES**; Fable was not needed.
The owner ratifies this ADR with the retained review caveats. The full authority, accepted contract,
caveats, source blockers, and INV-11 subpart disposition are recorded in
`docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_RATIFICATION.md`.

This ratification is based on the durable closure at
`5f04b5ef7cefc7ceb0d73fb0b3927bfbb28508be` (`PHASE5D0_CLOSURE_SHA`) plus the owner-supplied final
RR2 disposition. It does not reconstruct an unavailable external report. The historical proposal,
review, and revision language in sections 11–12 remains as the ledger of those earlier states; it
does not override the current accepted status in the document header and section 10.

ADR-0014 is **ACCEPTED / HUMAN-RATIFIED** for the controlled non-superfluid v1 scope. Production
implementation has not begun. Phase-5C canonical integration remains the next repository dependency;
realistic FR2005/A18 closure remains source-limited and blocked; global INV-11 remains unresolved;
and BNV has not begun.

---

## 14. Canonical integration addendum — 2026-09-10

The already-reviewed Phase-5D history was reconciled with canonical Phase-5C integration by the
governed merge recorded in
`docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md`. The historical dependency
statement in section 13 records the state at ratification; that dependency is now satisfied.

ADR-0014 is **ACCEPTED / HUMAN-RATIFIED / CANONICALLY INTEGRATED**. Its scientific contract is
unchanged. Phase-5C remains implemented, governed-regression protected, and closed for the
generic/free-gas coefficient scope. Phase-5B and Phase-5C compiler-portability semantics remain
unchanged. Phase-5D production implementation has **NOT BEGUN**. Global INV-11 remains
**UNRESOLVED**; realistic FR2005/A18 remains **SOURCE-LIMITED / BLOCKED**; BNV has **NOT BEGUN**.

---

## 15. Phase-5D-1 implementation-ratification addendum — 2026-09-12

The accepted equations and ownership contract above are unchanged. The
controlled non-superfluid implementation now exists on
`physics/phase5d-controlled-rotochemical-evolution`, with implementation commit
`d3670f6d4e021def0483909b6d2fdeed1c6973a4` and candidate evidence commit
`3486b972f71f57e8351fa8320c1ffb250fcd5c42`. The Opus independent review returned
**B — PASS WITH NONBLOCKING FINDINGS**, with 0 blocking, 0 material, 7
nonblocking, and 12 notes; Fable was not needed. The human owner ratifies the
implementation only for the controlled mathematical/architecture frozen-v1
scope documented in
`docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_RATIFICATION.md`.

The prescribed `P0=1 ms` spin history is super-Kepler / physically inadmissible
for the approximately `0.624 Msun`, `12.77 km` free-gas fixture for roughly the
first `~2 Gyr` under empirical mass-shedding estimates. It is accepted only as
a mathematical frozen-`W` driver; no physical pulsar interpretation is
permitted. The implemented envelope formula is FR2005 eq. (49) / PCY97's fully
accreted-envelope fit even though historical implementation metadata labels it
`iron Potekhin1997`; that label is not physically correct and must be repaired
before governed baseline promotion without changing the ratified trajectory.

Scaled RKF45 is accepted as adequate only for this controlled benchmark. The
late evolution is stability-bound with `h|lambda|` approximately `3.5-3.7`, so
solver strategy must be reconsidered before direct-Urca or superfluid
extensions if needed. No general realistic-evolution adequacy claim follows.

INV-11b/c/d are owner-resolved for controlled frozen-v1, INV-11e for frozen-v1
coefficient semantics, and INV-11f for controlled frozen-v1 ODE/source
coupling, all **not yet canonically integrated**. Global INV-11 remains
**UNRESOLVED**. Realistic Fernández–Reisenegger/A18 reproduction and
normalization remain **SOURCE-LIMITED / BLOCKED**; no realistic A18 work or BNV
has begun. No governed Phase-5D baseline is installed, and canonical `master`
does not yet contain this implementation.
