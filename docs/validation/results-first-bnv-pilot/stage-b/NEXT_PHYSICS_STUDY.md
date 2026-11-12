# Proposed study: competition between BNV chemical driving and rotochemical heating

**PROPOSED — owner consideration. Revised 2026-11-12 from the draft at
`0e88ec7e5ee981fcc2d93c85cf8a17d53d6d15fd`. No implementation, trajectory or
run is authorized by this document, and none was performed for it.**

The earlier draft made the retained-energy fraction the central variable. This
revision makes the **competition between BNV-induced chemical disequilibrium and
spin-down (rotochemical) disequilibrium** the central objective. The retained-
energy partition becomes a secondary axis (§6, §7.5). The alternatives the owner
must choose between are in
[ADR-0025](../../../adr/ADR-0025-bnv-retained-energy-study.md).

The Stage A/B P2 results are fixed-input mathematical benchmarks. Their
equations, numerical evidence and conditional qualification are not altered here:
[CONDITIONAL_ACCEPTANCE.md](CONDITIONAL_ACCEPTANCE.md), [CLOSEOUT.md](CLOSEOUT.md)
and `review-evidence/` are unchanged, and CLOSEOUT.md's one-line description of
this file refers to the earlier draft. Stage B ratification and approval of any
new physics remain separate owner decisions.

## 1. Summary

1. **The combined equation is already accepted.** ADR-0014 §3.4 and ADR-0015 A3
   give one chemical ODE with both drivers. Only its combined use, its diagnostics
   and its energy ledger with imposed spin are new (§2).
2. **The two mechanisms are not additive heaters.** They add as *count-rate
   drives* into the same two imbalances, and heating is a nonlinear function of
   the sum (§3).
3. **A neutron sink opposes spin-down.** Spin-down raises the equilibrium lepton
   content (η > 0, neutron decays). Neutron loss lowers it (η < 0, captures). On
   the governed Structure-1 coefficients the two drive vectors are antiparallel to
   within 16%, so they can nearly cancel in both channels at once (§3, §4).
4. **One number controls the competition:**
   `beta_l = D t_l / (2 |I_Omega,l| Omega |Omegadot|)`. On Structure-1,
   `beta_e = 1` at `Gamma ≈ 1.2e-10 yr^-1 × (Pdot_-20 / P_ms^3)`, where
   `Gamma = D/B`. For a 5 ms, 1e-20 pulsar that is ≈ 1e-12 yr^-1, the Stage B
   rate (§4). This is a fixture ratio, not an astrophysical prediction.
5. **Negligible direct heating is not generic.** For an event whose products all
   leave the fluid, the deposit is the Fermi-hole energy `mu_n,actual − E_n`.
   Competition requires that deposit to be below about `|t_l·eta*|/beta`, which
   is ≈ 0.5–1 keV per destroyed baryon on Structure-1. Sea-uniform decays deposit
   tens of MeV (§6).
6. **Recommended first calculation:** a static quasi-steady-state map of
   `(eta_e, eta_mu, T)` over the two drives at zero direct deposition, followed
   by four short matched coupled runs that verify it (§7, §9, §10).

## 2. Governing coupled equations

### 2.1 Conventions and authority

Infinity-coordinate time; redshifted `T`, `eta` and potentials; energies per
baryon in MeV; global powers in erg/s; `c = MeVToErg`. State basis
`y = (N_n, N_e, N_mu)`, `N_p = N_e + N_mu`, `b = (1,1,1)`, `P` selects e/μ,
`eta = −Z ell` with `ell = N_l − N_l^eq`, `R` positive for neutron decay.
These are ADR-0015 §3 and ADR-0014 §§3.1–3.4, unchanged.

| Object | Definition | Authority |
|---|---|---|
| `Z` | symmetric positive 2×2, MeV/count | ADR-0013; ADR-0014 §3.1 |
| `I_Omega` | `(dN_l^eq/dOmega^2)_B`, count s² | ADR-0014 §3.4 |
| `W` | `Z I_Omega` | ADR-0014 §3.4 |
| `t` | `(dN_y^eq/dB)_Omega`, `b·t = 1` | ADR-0015 §§6–7 |
| `R_l` | `(Ltilde_l/k_B) T^7 H_M(xi_l)`, `xi_l = eta_l/(k_B T)` | ADR-0014 §§3.4, 3.6 |

Fernández & Reisenegger 2005 (FR05) is the published formalism. ADR-0014 is
the project's governed reading of it: R2006 supersedes FR05 eqs. (9), (12),
(13), (54)–(56), the `H_M` denominator is `pi^8`, and the literature sign
warning of ADR-0014 §3.3 applies. Nothing below re-derives those.

### 2.2 Moving equilibrium reference with two drivers

The equilibrium composition is a function of two slowly changing quantities,
`N_y^eq(B, Omega^2)`. Its total rate of change is

```text
dN_l^eq/dt = t_l Bdot + I_Omega,l d(Omega^2)/dt
           = −D t_l + 2 I_Omega,l Omega Omegadot,        D = −Bdot > 0.      (1)
```

The BNV source also changes the actual composition directly. For a charge-
consistent raw source `S_y` with `b·S_y = Bdot`, the moving-reference source is
`sigma = P(S_y − t Bdot) = S_l + D t_l`, where `S_l = (S_e, S_mu)` is the raw
lepton source (ADR-0015 §7). Spin-down changes no particle number, so its whole
effect is the reference motion.

### 2.3 Chemical evolution and the net drive

```text
ell_dot = R + s,                 eta_dot = −Z (R + s)                        (2)
s = sigma + s_Omega              [count/s]
sigma   = S_l + D t_l            (BNV; = D t_l for a pure neutron sink)
s_Omega = −2 I_Omega Omega Omegadot   (spin; componentwise < 0 for spin-down,
                                       because I_Omega < 0 and Omega Omegadot < 0)
```

Equation (2) is ADR-0015 A3 with `Zdot = 0`, and reduces to ADR-0014 §3.4 at
`sigma = 0` and to the Stage B equation at `Omegadot = 0`. `−Z s_Omega` equals
`2 W Omega Omegadot`, the governed spin term. **Both mechanisms act only through
the single vector `s`.** Reactions see `eta`, not its origin.

Signs for spin-down and a neutron sink:

| Driver | Effect on `N_l^eq` | Sign of `s_l` | Sign of driven `eta_l` | Net reaction |
|---|---|---|---|---|
| Spin-down | rises (compression) | `s_Omega,l < 0` | `eta_l > 0` | neutron decay, `R > 0` |
| Neutron sink | falls (`t_l > 0`) | `sigma_l > 0` | `eta_l < 0` | capture, `R < 0` |

The first row is ADR-0014 §3.4; the second is ADR-0015 §7. They are opposite.

Stoichiometry decides the sign. Per destroyed baryon, with `nu = S_y/D`:

| Process class | `nu` | `sigma/D = (e, mu)` | Relation to spin-down |
|---|---|---|---|
| Neutron sink (one- or two-neutron) | `(−1,0,0)` | `(t_e, t_mu)` | opposes, both channels |
| Proton + electron sink | `(0,−1,0)` | `(t_e − 1, t_mu)` | reinforces e strongly; opposes μ weakly |
| Proton + muon sink | `(0,0,−1)` | `(t_e, t_mu − 1)` | reinforces μ strongly; opposes e weakly |

On Structure-1, `t_l = (0.0309, 0.00338)`, so a lepton-removing process drives
the electron channel about 31 times harder per baryon than a neutron sink, with
the opposite sign. The neutron sink is the weakest chemical driver. The first
experiment keeps it because it is the only source already implemented.

### 2.4 Reaction model and relaxation

Unchanged from ADR-0014 §§3.6–3.7: `L_H = c eta·R ≥ 0`,
`DeltaLnu_l = Ltilde_l T^8 [F_M(xi_l) − 1] ≥ 0`, one `Ltilde` for every role.

- Linear regime (`|xi| << 1`): `R ≈ Lambda(T) eta`, `Lambda` diagonal,
  `∝ T^6`. Then `eta_dot = −Z Lambda eta − Z s`; relaxation rates are the
  eigenvalues of `Z Lambda`.
- Nonlinear regime (`|xi| >> 1`): `R_l ∝ eta_l^7`, independent of `T`.
  Local relaxation rate near a driven state is `7 Z_ll |s_l| / |eta_l|`.
- Quasi-steady state (QSS): `R_l(T, eta_l*) = −s_l` in each channel. It is
  independent of `Z` and requires a positive `Ltilde_l` (ADR-0015 §10).

### 2.5 Energy ledger

```text
Echem = (1/2) eta^T Z^-1 eta                                   [MeV]
A_B     = −c eta·sigma                      BNV transfer into the reservoir
A_Omega = +c 2 Omega Omegadot eta·I_Omega = −c eta·s_Omega      spin work into it
c Echem_dot = A_B + A_Omega − L_H                                            (3)

mu_n,actual = mu_B + t_l·eta
P_actual = c (mu_n,actual − E_esc,fluid) D          (neutron sink, no inflow)
P_0      = P_actual + A_B = c mu_B D − L_out,fluid
Cstar Tdot = P_actual + L_H − DeltaLnu − Lnu,eq − Lgamma                     (4)

d/dt [ U + c Echem ] = P_actual + A_B + A_Omega − DeltaLnu − Lnu,eq − Lgamma (5)
```

Equation (3) is ADR-0015 A5 with frozen coefficients. Equation (4) is ADR-0015
A6 unchanged: **spin adds no term to the thermal equation.** Equation (5) is
their sum and is the conservation law this study must satisfy.

`A_Omega > 0` when spin-down acts on `eta > 0`: compression work is stored.
`A_B > 0` when a neutron sink acts on `eta < 0`. When the two drivers oppose,
each drains what the other stores: `A_B < 0` on a spin-built reservoir, and
`A_Omega < 0` on a BNV-built one. At QSS, `L_H = A_B + A_Omega`.

Where the drained energy goes depends on the partition. With all energy escaping
at the actual potential (`E_esc = mu_n,actual`, `P_actual = 0`), each event
carries `mu_B + t_l·eta` out of the star. For `eta > 0` that exceeds `mu_B`: the
sink exports spin-supplied chemical energy that Urca reactions would otherwise
have turned into heat and neutrinos.

**Double-counting refusals.**

1. The thermal right-hand side contains `P_actual + L_H − DeltaLnu`. It never
   contains `A_B`, `A_Omega`, `−Echem_dot` or `Omegadot J` (ADR-0015 §§9, 19).
2. `A_Omega` enters the chemical reservoir only. Rotational energy loss is not
   heat.
3. Fermi-hole energy is already inside `mu_n,actual − E_esc,fluid`.
4. Escape energy is counted once. `E_esc = mu_n,actual` already includes the
   chemical excess `t_l·eta`; do not also subtract `A_B` from the thermal
   equation.
5. A prescribed observed `Omegadot` already contains any BNV-induced spin
   change. Do not add a second BNV torque term to the chemical drive (§8.4).
6. The interaction contrast of §7.4 is a diagnostic, not a third heat source.

### 2.6 What is accepted and what is derived here

| Statement | Status |
|---|---|
| Eq. (2), combined chemical ODE | Accepted contract: ADR-0015 A3, ADR-0014 §3.4 |
| Eq. (3), reservoir derivative with spin term | Accepted: ADR-0015 A5 |
| Eq. (4), thermal equation | Accepted: ADR-0015 A6 |
| Eq. (5), thermal + chemical subsystem ledger with spin | **Derived here; unratified candidate** |
| Global closure including `E_eq(B, Omega)` | **Not claimed** (ADR-0015 §9) |

The global rotating-star closure cannot be tested at chemical precision in the
frozen-v1 model. `mu_B` and `t` are evaluated on the nonrotating background, so
`(dE_eq/dB)_Omega` carries an unowned `O(Omega^2/Omega_K^2)` relative error.
That is tens of MeV per baryon, against chemical energies below 1 keV. The
thermal and chemical evolution do not depend on it: with `P_actual` defined as a
residual, only the subsystem ledger (5) enters. A global statement needs an
owner for `Omegadot J` work and rotating first-law conjugates.

Two sign-conditional statements in ADR-0015 assume a neutron-sink-built state
and **do not hold when spin-down makes `eta > 0`**: the bound on `P_dir` that
requires `delta mu_n ≤ 0` (§5), and "nonpositive while storage fills, exactly
`−DeltaLnu` in reached QSS" (§10). With spin, the correction relative to `P_0`
is `−c Echem_dot + A_Omega − DeltaLnu`, which is `A_Omega − DeltaLnu` at QSS
and can have either sign.

## 3. Physical interpretation

Spin-down compresses the star, and the equilibrium lepton number rises. The
actual lepton number lags, so neutrons are over-abundant and decay. Neutron loss
decompresses and dilutes the star, and the equilibrium lepton number falls. The
actual lepton number lags the other way, so leptons are over-abundant and are
captured. A star subject to both is driven by the difference.

The coupling is nonlinear in three distinct ways.

1. **Energies are quadratic even when kinetics are linear.** At fixed `T` and
   small `xi`, `eta = eta_RC + eta_BNV` exactly. But
   `L_H = c eta^T Lambda eta` and `Echem` contain the cross terms
   `2c eta_RC^T Lambda eta_BNV` and `eta_RC^T Z^-1 eta_BNV`, which are negative
   for opposing drives.
2. **Rates are nonlinear.** At large `xi`, `R ∝ eta^7`, so `L_H ∝ |s|^(8/7)`.
3. **Temperature feeds back.** `Lambda ∝ T^6`, and the net weak response changes
   sign at `|xi| ≈ 4.9097` (ADR-0014 §3.6). A partial cancellation that lowers
   `|xi|` can turn net weak heating into net weak cooling during a transient.

Cancellation is per channel and never complete. The electron and muon nulls sit
at different `beta` because `t_l` and `I_Omega,l` are not exactly parallel (§4).
Between them the two imbalances have opposite signs.

## 4. Dimensionless measures

### 4.1 BNV versus spin-down driving

```text
beta_l = −sigma_l / s_Omega,l = D t_l / (2 |I_Omega,l| Omega |Omegadot|)      (6)
s_l    = s_Omega,l (1 − beta_l)
```

for a neutron sink. `beta_l < 0` denotes reinforcement (lepton-removing
stoichiometry). `beta` depends only on stellar structure (`t`, `I_Omega`) and the
two drives. It does not depend on reaction rates, superfluidity or temperature.

With `Gamma = D/B`, characteristic age `tau_c = P/(2 Pdot)` and
`Omega Omegadot = −4 pi^2 Pdot/P^3`:

```text
beta_l = Gamma / Gamma_crit,l,
Gamma_crit,l = 2 |I_Omega,l| Omega |Omegadot| / (t_l B)
             = 8 pi^2 |I_Omega,l| Pdot / (t_l B P^3)                          (7)
```

`Gamma_crit` depends on spin only through `Pdot/P^3`, the same combination that
sets the FR05 stationary luminosity (FR05 eqs. 66–68).

### 4.2 Direct versus chemical thermal power

Let `q_dep = P_actual/(cD)` be the direct deposit per destroyed baryon and
`a = |t_l·eta|` the chemical work per destroyed baryon.

```text
delta     = P_actual / L_H                      direct vs gross chemical heating
delta_net = P_actual / (L_H − DeltaLnu)         direct vs net weak response
delta_RC  = P_actual / L_H^RC = beta q_dep / a_RC                             (8)
```

`L_H^RC` is the spin-only QSS heating and `a_RC = sum_l t_l |eta_l^RC|`. The
last equality uses `s_Omega,l = −D t_l / beta_l` and one common `beta`. Ratios
with a vanishing denominator are undefined; report raw powers there.

BNV modifies rotochemical heating materially, without being swamped by its own
direct heating, when

```text
beta ≳ 0.1     and     delta_RC ≲ 1     <=>     q_dep ≲ a_RC / beta.          (9)
```

| Regime | Condition | Thermal character |
|---|---|---|
| I | `|beta| << 1`, `delta_RC << 1` | Conventional rotochemical heating |
| II | `|beta| << 1`, `delta_RC ≳ 1` | Direct BNV heating over a rotochemical background; chemistry irrelevant |
| III | `|beta| ≳ 0.1`, `delta_RC ≲ 1` | **BNV-modified rotochemical heating (target)** |
| IV | `|beta| ≳ 1`, `delta_RC ≳ 1` | BNV dominates both; chemistry is an `O(a/q_dep)` correction |

### 4.3 Reachability and relevance

- `rho_build = tau_build / tau_evol`, with `tau_build,l = |eta_l*| / |(Z s)_l|`
  and `tau_evol` the shortest of `Omega/|Omegadot|`, `D/|Ddot|`, `T/|Tdot|` and
  the age. QSS claims need `rho_build << 1` in the applicable regime.
- Depletion needed to build a BNV-driven state from zero:
  `Gamma tau_build ≈ |eta_l*| / ((Z t)_l B)`. It depends on `Gamma` only through
  `eta* ∝ Gamma^(1/7)`.
- Thermal relevance: `(L_H − DeltaLnu) / (Lnu,eq + Lgamma)`.
- Spin-feedback size: `2 kappa |I_Omega,l| Omega^2 / (t_l B)`, with
  `kappa = (Omegadot/Omega)_BNV / Gamma` (§8.4).

### 4.4 Values on the governed Structure-1 fixture

Inputs, all existing governed or retained values:

| Quantity | Value | Source |
|---|---|---|
| `t_e, t_mu` | 0.030852171, 0.0033777438 | `review-evidence/outputs/longterm.json:1` |
| `mu_B^inf` | 868.4354349 MeV | same |
| `I_Omega,e`, `I_Omega,mu` | −1.16380e47, −1.48448e46 count s² | `tests/baselines/phase5d1_controlled_evolution.json:1` |
| `Z` | [[4.5793e-54, 5.1725e-55],[5.1725e-55, 1.02687e-52]] MeV/count | same |
| `Ltilde_Me`, `Ltilde_Mmu` | 1.81291e-32, 9.85050e-34 erg s⁻¹ K⁻⁸ | same |
| `N_n, N_p, N_e, N_mu` | 7.5683e56, 4.8595e54, 4.8119e54, 4.7577e52 | `tests/baselines/phase5b_structural_response.json:1` |
| `B` | 7.6169e56 | `N_n + N_p`, as in `tests/native/pair.hpp:20` |

Derived ratios (arithmetic on the table above; no solve, no integration).
Constants used: `k_B = 8.617333262e-11 MeV/K`, `1 MeV = 1.602176634e-6 erg`,
Julian year `3.15576e7 s`, `Pdot_-20 = Pdot/1e-20`, `P_ms = P/1 ms`:

| Measure | Electron | Muon |
|---|---|---|
| `t_l / |I_Omega,l|` | 2.651e-49 | 2.275e-49 (ratio 1.165) |
| `beta_l / [Gamma tau_c (P/1 ms)^2]` | 5.115 | 4.390 |
| `Gamma_crit,l / (Pdot_-20/P_ms^3)` | 1.23e-10 yr⁻¹ | 1.44e-10 yr⁻¹ |
| `d ln N_l^eq / d ln B` | 4.88 | 54.1 |
| `|I_Omega,l| Omega^2 / N_l` at P = 1, 5, 10 ms | 0.95, 0.038, 0.0095 | 12.3, 0.49, 0.12 |
| Spin feedback `2|I|Omega^2/(tB)` at 1 ms | 0.39 | 0.46 |

Consequences on this fixture:

- `beta_e = 1` at `Gamma ≈ 1e-12 yr^-1` for `P = 5 ms`, `Pdot = 1e-20`. At
  that point depletion over one characteristic age is
  `Gamma tau_c ≈ 0.2 (1 ms/P)^2 ≈ 0.8%`.
- The mixed-sign window is `1 < beta_e < 1.165`.
- The muon channel is near threshold. Its equilibrium population responds 54
  times faster than `B` and by order unity to rotation faster than about 10 ms.
  Muon results on this star are not representative.
- An empirical mass-shedding estimate gives `P_K ≈ 1.5 ms` for this star,
  consistent with ADR-0014 §15's super-Kepler classification of the governed
  `P0 = 1 ms` history.

Large-`xi` QSS illustrations, from `|s_l| = (Ltilde_l/k_B) [24/(11513 pi^8)]
(eta_l/k_B)^7` and `L_H = c sum_l |eta_l s_l|`. These are closed-form estimates
to fix scales, not the proposed map:

| Case (P, Pdot, Gamma) | `eta_e, eta_mu` (keV) | `L_H` (erg/s) | `N_LH` (§7.4) |
|---|---|---|---|
| 5 ms, 1e-20, spin only | +19.0, +21.5 | 2.6e28 | — |
| 5 ms, 1e-20, BNV only at 1e-12 yr⁻¹ | −19.1, −21.1 | 2.5e28 | — |
| 5 ms, 1e-20, both at 1e-12 yr⁻¹ | −10.2, +16.1 | 4.7e26 | −0.99 |
| 5 ms, 1e-20, both at 3e-13 yr⁻¹ | +18.1, +20.6 | 1.7e28 | −0.47 |
| 2.5 ms, 1e-19, both at 1e-11 yr⁻¹ | +34.9, +39.5 | 3.3e30 | −0.21 |

Further scales from the same arithmetic at 5 ms, 1e-20, 1e-12 yr⁻¹:

- Chemical work per destroyed baryon `a ≈ 0.66 keV`, so `a/mu_B ≈ 7.6e-7`.
- Build time from zero: 1.7e8 yr (e), 7.6e7 yr (μ) for either driver alone.
  Build depletion: 1.7e-4 (e), 7.6e-5 (μ). One e-fold near QSS costs one seventh
  of that.
- Stored chemical energy at the spin-only state ≈ 6.6e43 erg, roughly a hundred
  times the thermal content at `T ~ 1e6 K` (using `Cstar/T ≈ 1.3e30 erg/K²`
  from `evidence/final-summary.json:1`).
- Interior QSS temperature: a crude power-law extrapolation of the Stage B
  endpoint `Lgamma(T)` puts `L ~ 1e28–1e29 erg/s` at `T ~ 3e5–1e6 K`, at or
  below the 1e6 K certificate floor. The envelope owner must replace this
  estimate.

## 5. Quasi-steady-state scalings and distinguishing signatures

In the modified-Urca large-`xi` regime, `eta_l* ∝ |s_l|^(1/7)` and
`L_H ∝ |s_l|^(8/7)` per channel, with net retained fraction 5/8 (ADR-0014
§3.6). With photon-dominated cooling, `T_s ∝ L^(1/4)`. FR05 eqs. (64)–(68) are
the `sigma = 0` case. Replacing `s_Omega` by `s` gives, for a neutron sink:

```text
L_gamma,QSS ∝ sum_l (|I_Omega,l|^8 / Ltilde_l)^(1/7) |1 − beta_l|^(8/7)
              |Omega Omegadot|^(8/7)
T_s(BNV + spin) / T_s(spin) ≈ |1 − beta|^(2/7)      (one-channel form)      (10)
tau_eq ∝ |s|^(−6/7)                                  (FR05 eq. 81 with s)
```

These are scaling predictions to be tested by the map of §7. Use the full
polynomials near a null, where `|xi|` need not be large.

Qualitative signatures that distinguish BNV-modified from conventional
rotochemical heating:

1. **Luminosity is no longer a function of spin-down alone.** FR05's stationary
   photon luminosity depends only on `Omega Omegadot`. With BNV it depends on
   `Pdot/P^3` and `Gamma` through `1 − beta_l`.
2. **A neutron sink changes the temperature by `|1 − beta|^(2/7)`**: colder
   than the spin-only star for `0 < beta < 2`, coldest at `beta = 1`, and
   hotter for `beta > 2`. Because `beta ∝ Gamma P^3/Pdot`, at fixed `Gamma`
   the fastest spin-down stars are unchanged, intermediate ones are colder, and
   the slowest are hotter than the spin-only prediction.
3. **A cold gap at `beta ≈ 1`**, between a spin-dominated branch
   `T_s ∝ (Pdot/P^3)^(2/7)` and a BNV-dominated branch `T_s ∝ Gamma^(2/7)`
   that is independent of spin-down. A direct-heating floor is also spin-
   independent but scales as `Gamma^(1/4)` and carries no imbalance.
4. **Sign reversal of the imbalance.** `eta < 0` (net capture) for `beta > 1`,
   against `eta > 0` for conventional rotochemical heating.
5. **A mixed-sign window** between the electron and muon nulls.
6. **Asymmetric non-additivity.** For opposing drives the interaction contrast
   reaches −1. For reinforcing drives it is at most about +0.10 in the large-`xi`
   regime (`2^(8/7) = 2.208` against 2) and +1 in the linear regime.
7. **Stoichiometric discrimination.** A lepton-removing process heats the star
   above the FR05 value, with about 31 times the leverage per baryon.
8. **Slower equilibration near cancellation.** `tau_eq ∝ |s|^(−6/7)` grows, so
   stars near `beta = 1` are more likely to be out of QSS and history-dependent.
9. **Energy export.** Under cancellation, spin-supplied chemical work leaves
   with the BNV products instead of appearing as heat.

None of these is an observational claim. FR05 reports a predicted temperature
for PSR J0437−4715 about 20% below the inferred one. A neutron sink would move
the prediction further below only for `0 < beta < 2`, and above it for
`beta > 2`; no value of `beta` is inferred for that pulsar here.

## 6. When can direct deposition be negligible?

### 6.1 The ledger constraint

For an event that removes a neutron of energy `E_n` and whose products deposit
nothing further, `E_esc,fluid = E_n` and the deposit is the hole energy:

```text
q = mu_n,actual − E_n ≥ 0        (cold, spontaneous removal; ADR-0015 §5)    (11)
```

This holds whether the products are visible, invisible, escaping or retained in
a hidden sector. Product fate changes the global ledger and the X reservoir, not
(11). Zero direct deposition therefore means **events drawn from the Fermi
surface**, not merely invisible products. ADR-0015 §5 already states that
invisible is not generally zero heat.

Condition (9) sets the scale: `q_dep ≲ a_RC/beta`, about 0.5–1 keV per
destroyed baryon on Structure-1 at `beta ~ 1`, a fraction ≲ 1e-6 of `mu_B`.

### 6.2 In-medium energy and momentum conservation

- **Neutron spectrum.** In relativistic mean-field matter
  `E(k) = sqrt(m*^2 + |k − Sigma|^2) + Sigma^0` (Berryman, Gardner & Zakeri
  2023, eq. 4.8, PDF p. 25). `m*` falls and `Sigma^0` rises with density, so
  `E(0) > m_B` at high density and decays to final states heavier than the
  vacuum baryon become possible in dense cores (text after their eq. 4.8, and
  their §5 discussion of Fig. 9).
- **Which neutrons decay.** For `n → chi + light`, the invariant
  `E(k)^2 − |k|^2` must exceed the squared threshold mass. With a repulsive
  vector self-energy it increases with `k`, so a threshold can admit only the
  top of the Fermi sea. Away from threshold the whole Fermi sphere contributes,
  weighted by the momentum-dependent width (their eq. 5.10, PDF p. 33).
- **Escape.** A massive product escapes only if `e^Phi E_local ≥ m c^2`
  (ADR-0015 §§5, 13). The Killing energy available per removed neutron is at
  most `mu_n,actual^inf = mu_B^inf + t_l·eta`, which is 868.4 MeV on
  Structure-1 and lower for more compact stars. The model window quoted by
  Berryman et al. (their eq. 2.10, PDF p. 12) is
  `937.993 MeV < m_chi < 1077.84 MeV`. A dark baryon in that window cannot
  escape promptly from this star; it is retained unless dark-sector decay or
  annihilation to light states removes it. That removal is the working
  assumption of their analysis (their §3, PDF p. 15).

### 6.3 Classification

| Class | Mechanism | Mean deposit per event | Negligible? |
|---|---|---|---|
| M0 | Mathematical limit `E_esc = mu_n,actual` (existing P0 fixture) | 0 by construction | Benchmark only; no process claimed |
| M1 | Sea-uniform disappearance: decay or collisional conversion from anywhere in the Fermi sea, all products leave the fluid | Tens of MeV. ADR-0015 §5: `2 E_F,kin/5` for a uniform nonrelativistic sea. McKeen, Pospelov & Raj 2021: volume-averaged `<mu_n> = 68 MeV` (PDF p. 3); 2020: `Delta E ≈ 100 MeV` adopted (PDF p. 6) | **No.** `delta_RC ~ 1e4–1e5 × beta` |
| M2 | Threshold-gated in-medium decay: only a shell of width `Delta(r)` below the Fermi surface is allowed | Up to `Delta`, which varies by ~100 MeV across the star | Only if the threshold is barely open at the centre. Fine-tuned; `D` then depends on compression |
| M3 | Fermi-surface-selective process: thermally activated, or otherwise restricted by Pauli blocking to within `k_B T` of the Fermi surface. No published example is identified here | `O(k_B T)`, either sign | **Candidate, unsourced.** Marginal: `k_B T/a ≈ 0.05–0.4` at `T = 3e5–3e6 K`; needs the finite-`T` terms and a rate derivation |
| M4 | Retained product (cannot escape; §6.2) | Same hole energy as M1 for the ordinary fluid | No for sea-uniform events. Adds an X reservoir (ADR-0015 Regime II) unless removed |
| M5 | Light invisible final states allowed in vacuum (neutrinos, light dark states) | Hole energy as M1 | No. Also `tau(n → inv) > 9.0e29 yr` (SNO+ 2022) makes `beta` negligible by many orders |

Collisional conversion deposits more than the hole: McKeen, Pospelov & Raj 2021
count the converted-neutron hole, the spectator hole and the recoil (PDF p. 2).

The resonant baryon–dark-baryon conversion of Berryman et al. 2023 is **not** a
Fermi-surface mechanism. Their resonance condition (eq. 5.13, PDF p. 38) holds
for at most one momentum magnitude at a given density, in a shell of width
`~epsilon_Bchi` inside the Fermi sphere, and they conclude that it is not a
phenomenologically relevant production mechanism (PDF p. 40). Its deposit would
be small only where the resonant momentum happens to equal the Fermi momentum.
It is not a candidate for class M3.

### 6.4 Conclusions for the study

- The zero-deposition four-case experiment measures the chemical coupling
  cleanly and bounds the purely chemical effect. It is a diagnostic limit.
- A physical claim requires a process in class M2 or M3 with a derived rate,
  stoichiometry, spectrum and fate, and `q_dep ≲ a_RC/beta` demonstrated.
- For class M1, regime III is empty: where `beta` matters, direct heating is
  `1e4–1e5` times the rotochemical power; where direct heating is small, `beta`
  is negligible. On this reading, thermal limits on such processes are far
  tighter than binary-timing limits, in line with McKeen, Pospelov & Raj.
  Berryman et al. 2023 leave the net thermal effect open (§11); whether Urca
  cooling changes this conclusion is one of the questions the study answers.
- An arbitrary retention fraction is not a physical decay. Report `q_dep` in keV
  per event beside `a`, and name the class.
- Superfluid gaps would raise `eta*` and therefore `a`, relaxing (9) in
  proportion. That is an expectation, not a result.

Finite-temperature status. The smooth-weight floor
`(pi^2/6)(k_B T)^2/E_F,kin` is about 1e-7 keV at `T ~ 1e6 K`, far below `a`.
The ADR-0015 §15 statement that the P0 tier is not a resolved thermal-sign
experiment was made at 1e8 K and the old rates, where `a` was 0.065 eV. In the
cold QSS regime the smooth-weight class clears that floor. A Fermi-surface-
concentrated source does not: its `O(k_B T)` residual is a tenth or more of `a`
and must be owned before any sign claim.

## 7. Matched four-case numerical experiment

### 7.1 Cases

| Case | Spin drive | BNV drive | Direct deposition | Existing analogue |
|---|---|---|---|---|
| 1 `(0,0)` | off | off | — | Stage B exact-zero control |
| 2 `(1,0)` | on | off | — | Phase-5D governed rotochemical run |
| 3 `(0,1)` | off | on | zero (M0) | P0 cold-cancellation fixture |
| 4 `(1,1)` | on | on | zero (M0) | none; combination currently refused |

### 7.2 Matching rules

- Identical star, `Z`, `W`, `I_Omega`, `t`, `Ltilde`, `Cstar`, envelope, solver
  and tolerances; identical `T(0)` and `eta(0)`.
- Cases 2 and 4 use the **same prescribed `Omega(t), Omegadot(t)`**. This is
  observational matching: two stars with the same timing solution. Matching the
  external torque instead changes `Omegadot` by `kappa Gamma Omega`; its
  chemical effect relative to `sigma` is `0.39 kappa (1 ms/P)^2` on Structure-1,
  the same order as terms frozen-v1 already drops (§8.4).
- Cases 1 and 3 use the existing static zero-spin history. In frozen-v1 a
  rotating but non-decelerating star is identical.
- Cases 3 and 4 use the same `D(t)` and the same partition.
- Choose `(Pdot/P^3, Gamma)` to place the pair at declared `beta_e`, e.g.
  0.3, 1 and 3, and keep the spin sub-Kepler with `|I_Omega,l|Omega^2/N_l`
  reported.

### 7.3 Observables

Per case and per time or map point: `T_inf`, `T_s^inf`, `Lgamma`;
`eta_e`, `eta_mu`, `xi_e`, `xi_mu`; `R_l`; gross weak heating `L_H`; excess
neutrino emission `DeltaLnu`; net weak response `L_H − DeltaLnu`; net chemical
thermal response relative to `P_0`, `−c Echem_dot + A_Omega − DeltaLnu`;
`A_B`, `A_Omega`; `Echem` and `Echem_dot`; `P_actual`, `P_0`;
`Lnu,eq`; the residual of ledger (5); `beta_l`, `delta`, `delta_RC`,
`rho_build`; depletion and frozen-validity utilization; `|I_Omega,l|Omega^2/N_l`.

### 7.4 Quantifying the nonlinear interaction

For any observable `X` with case values `X_00, X_10, X_01, X_11`:

```text
I_X = X_11 − X_10 − X_01 + X_00                    interaction contrast
N_X = I_X / ( |X_10 − X_00| + |X_01 − X_00| )      normalized              (12)
```

- `N_X = 0` means additive; `N_X = −1` means complete cancellation;
  `N_X > 0` means reinforcement beyond the sum.
- Report the four marginal effects together with `I_X`: spin added first
  (`X_10 − X_00`) and last (`X_11 − X_01`), BNV added first and last. Their
  disagreement is `I_X`. No single number assigns heating to one mechanism.
- If one symmetric attribution is wanted, use the order average
  `phi_RC = [(X_10 − X_00) + (X_11 − X_01)]/2` and likewise `phi_BNV`; they sum
  to `X_11 − X_00`. Label it as a convention.
- Apply (12) first to quantities that are additive when the physics is:
  `eta_l` at fixed `T` in the linear regime (expected `N = 0`), then powers
  `L_H`, `L_H − DeltaLnu`, `Lgamma` and integrated energies. Temperature is a
  nonlinear function of power, so `N_T` is secondary and must not be read as an
  interaction of mechanisms.
- Reference curves for one channel: linear regime
  `N_LH = −2 beta/(1 + beta^2)`; large-`xi` regime
  `N_LH = (|1 − beta|^(8/7) − 1 − |beta|^(8/7))/(1 + |beta|^(8/7))`.
- BNV *enhances*, *suppresses* or *dominates* the rotochemical response
  according to `X_11 − X_10 > 0`, `< 0`, or `|X_11 − X_10| > |X_10 − X_00|`, each
  stated with its observable and time window (ADR-0015 §14).

### 7.5 Sensitivity axes

1. `Gamma` at fixed spin history: scan `beta_e` from 0.03 to 30.
2. Spin-down rate at fixed `Gamma`: scan `Pdot/P^3`.
3. Deposit per event at the `beta` values of interest:
   `q_dep/a ∈ {0, 0.1, 1, 10}` and one M1-like value, through the actual-
   potential partition of ADR-0025 option P-A.
4. Long-term QSS: follow the QSS locus as `Omega Omegadot` declines along a
   spin history, with `rho_build` reported at every point.

## 8. Stellar regimes, input adequacy and limitations

### 8.1 Target regime

Old neutron stars and millisecond pulsars: photon-dominated cooling,
`T_inf ~ 1e5–1e7 K`, `|xi| >> 1`, ages comparable to or exceeding the build
times of §4.4. This is where rotochemical heating matters (FR05) and where a
cold QSS makes the comparison insensitive to initial conditions.

### 8.2 Adequacy of existing inputs

| Input | Adequate for a first controlled comparison? | Limit |
|---|---|---|
| Structure-1 free gas, 0.6236 M☉, 12.77 km | Yes, for signs, `beta`, exponents and `N_X` | Not a millisecond-pulsar mass or EOS; absolute temperatures meaningless |
| `Z`, `W`, `I_Omega`, `t` (governed) | Yes | Nonrotating background; muon channel near threshold |
| Modified-Urca `F_M`, `H_M` | Yes | Non-superfluid; no direct Urca |
| `Ltilde` benchmark values | Yes for ratios. FR05 eq. 66 gives `L_gamma ∝ Ltilde^(−1/7)`, so `T_s ∝ Ltilde^(−1/28)` | Declared benchmark, not a source normalization (ADR-0014 §3.7) |
| Envelope and `Cstar` | Yes within their certified range | BNV certificate spans 1e6–1e9 K only |
| Prescribed dipole history `P0 = 1 ms` | Only as a mathematical driver | Super-Kepler for ~2 Gyr on this star (ADR-0014 §15) |
| Uniform neutron sink, P0/P2 partitions | Yes | Test fixtures; abstract source |

### 8.3 Limitations by topic

- **Superfluidity.** Gaps suppress Urca rates until the imbalance is comparable
  to the gap, so `eta*`, the stored energy and the QSS luminosity all change.
  FR05 notes superfluidity should raise its predicted temperature. `beta` is
  unaffected; `a` and condition (9) are.
- **Direct Urca.** Large-`xi` exponents become `L_H ∝ |s|^(6/5)` and
  `T_s ∝ |s|^(3/10)`, with smaller `eta*` and a two-timescale system. The solver
  decision must be revisited (ADR-0014 §3.16).
- **Slow rotation.** Coefficients are `O(Omega^0)` and the drive `O(Omega^2)`.
  Order-unity `|I_Omega,l|Omega^2/N_l` marks breakdown; the muon channel reaches
  it near 5 ms on this star.
- **Frozen structure.** `Z`, `t`, `I_Omega`, `Ltilde` drift with depletion and
  with spin. The cross-derivative `d^2 N^eq/dB dOmega^2` is omitted. Evolving
  coefficients need the `+Zdot Z^-1 eta` term (ADR-0014 §3.11).
- **Finite temperature.** §6.4.
- **In-medium kinematics.** §6.2; free-gas Structure-1 has no `m*` or `Sigma^0`.

### 8.4 BNV-related angular momentum

Berryman, Gardner & Zakeri 2022 take BNV as spherically symmetric, so it does
not change the angular momentum but does change the moment of inertia:
`Idot/I = Pdot_s/P_s ≈ ±Bdot/B` (their eq. 29), sign depending on mass. So
`kappa = O(1)`. Products that carry angular momentum add to `kappa`. Under
observational matching this is already inside the measured `Omegadot`. It
matters only for the torque-matched counterfactual and for inferring a magnetic
field.

### 8.5 Existing gates that the target regime touches

- **Depletion gate.** `CompactStar/Physics/BNV/src/FrozenBnvValidity.cpp:66` refuses depletion beyond 1e-6.
  Building a BNV-driven QSS from zero needs about 1e-4; one e-fold near QSS
  needs about 1e-5. This is a limit of the Stage A/B certificate, not a physical
  breakdown: the estimated coefficient drift at 2e-4 depletion is about 1%
  (from `d ln N_mu/d ln B ≈ 54`). That estimate is not a certificate.
- **Temperature gate.** `CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:108` requires
  1e6–1e9 K. Realistic millisecond-pulsar drives may put the Structure-1 QSS
  below it, and case 1 cools through it.
- **Spin with BNV.** Refused at construction and at every evaluation (§9).

### 8.6 First experiment versus publishable prediction

| Requirement | First controlled experiment | Publishable astrophysical prediction |
|---|---|---|
| Star | Structure-1 | ≥ 1.4 M☉ with an off-equilibrium EOS giving `Z`, `t`, `I_Omega` |
| Reactions | Non-superfluid modified Urca, benchmark `Ltilde` | Source normalization; gaps; direct Urca where open |
| BNV source | Abstract uniform neutron sink | Process-specific rate, stoichiometry, radial profile |
| Deposit | M0, plus a `q_dep` sweep | Derived spectrum, hole energy, finite-`T` terms, fate |
| Spin | Prescribed, sub-Kepler, observationally matched | Torque-matched with `kappa`; second-order rotation |
| Background | Frozen; QSS snapshots | Sliding background with evolving coefficients |
| Ledger | Subsystem ledger (5) | Rotating global first law with `Omegadot J` owner |
| Result | Ratios, exponents, `N_X`, regime map in `(beta, delta_RC)` | Temperatures and limits for named pulsars |

Realistic A18 closure remains source-limited (ADR-0014 §6). This study does not
unblock it.

## 9. Prioritized investigation plan

Each stage needs its own owner approval. None requires a new broad numerical-
certification campaign: validation is by regression to existing governed corners,
ledger closure, and step- or tolerance-halving.

**S1 — Static QSS map (recommended first; no production-source change).**
For each `(Pdot/P^3, Gamma)`:
solve `R_l(T, eta_l) = −s_l` for `l = e, mu` and the thermal balance
`P_actual + L_H − DeltaLnu − Lnu,eq − Lgamma(T) = 0` at `P_actual = 0`.
The four cases are the corners `s_Omega = 0` and `D = 0` of the same map.
A standalone tool reads the governed coefficient artifacts; it is three
monotone scalar equations per point. Outputs: every observable of §7.3 that
is defined at QSS, `N_X`, the regime map, and the gate status of every point.
Checks: FR05-form scaling at `Gamma = 0`; exponent 8/7 and 2/7 fits; null
positions against (6); Jacobian stability; the 1e6 K floor.
`Lgamma(T)` comes from the governed envelope relation (ADR-0014 §3.9 and §16:
FR2005 eq. 49 / PCY97 fully accreted). Its validity below 1e6 K has not been
established here. Until it is, restrict reported points to drives whose QSS lies
at or above 1e6 K, or flag colder points as outside the certified range.
**S1 is the stage that yields temperature comparisons.**

**S2 — Minimal combined-drive enablement and four short coupled runs.**
The right-hand side already sums both drives: the ordinary evaluator carries
`2 W Omega Omegadot` (`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:141`) and the BNV context
subtracts `Z sigma` (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:139`). After approval,
the changes are:

1. Admit spin with a nonzero BNV source. It is refused in the constructor
   (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:35`, `:63`) and at each evaluation
   (`:114`).
2. Add a parameterized sub-Kepler prescribed spin history.
   `PrescribedDipoleHistory` is hard-coded (`CompactStar/Physics/Rotochemical/PrescribedSpinHistory.hpp:26`) and
   is required by type (`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:76`, `:105`). Native
   Structure-1 inputs are restricted to analytic controls (`:90`).
3. Add `A_Omega` to the reservoir diagnostic and the residual of ledger (5).
   `Echem_dot_total` omits the spin term
   (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:164`).
4. Promote the P0 partition from a test fixture
   (`tests/bnv/energy_partition_fixtures.hpp:28`) to a run-card partition.

Runs: the four cases from one hot initial state, inside the existing 1e-6
depletion and 1e6–1e9 K gates. They test the drive slopes
`eta_dot = −Z s` and their cancellation, the linear-regime superposition, the
storage cross term and ledger (5). Regressions: case 2 against the Phase-5D
baseline, case 1 against the Stage B control, case 3 against P0 fixtures.
**S2 does not yield a temperature comparison.** Inside the 1e-6 depletion gate
the chemical temperature contrasts are at Stage B's microkelvin scale, below what
tolerance halving in double precision resolves. S2 verifies imbalance slopes,
powers and the ledger only.

**S3 — Reaching and tracking QSS in time (owner choice of route).**
Near-QSS initial data from S1, then either a narrowly scoped frozen-certificate
extension to depletion of a few 1e-4 with the existing sensitivity producer
(`tests/bnv/produce_frozen_certificates.py:1`), or short windows only. Add the
`q_dep` sweep (ADR-0025 option P-A) and the `Gamma` and spin-rate scans. Time-
dependent temperature comparisons come from this stage, after the gate decision.

**S4 — Toward a prediction.** One item at a time from §8.6, starting with
whichever S1–S3 show to be controlling.

Resources. S1 is a parameter grid and is trivially parallel. S2–S3 are
independent trajectories; run them concurrently on EKU according to measured
need and scheduler availability. No CPU or walltime ceiling is proposed.

## 10. Recommendation

**The simplest physically meaningful calculation is the static quasi-steady-
state map of S1 at zero direct deposition on the existing Structure-1
coefficients**: two decoupled channel balances `R_l(T, eta_l) = −s_l` and one
thermal balance, over the plane of spin drive `Pdot/P^3` and baryon-loss rate
`Gamma`. The four requested cases are its corners. It needs no new physics in
the production code and no time integration.

It establishes when BNV-induced disequilibrium materially changes a
rotochemically heated star, in three statements to be confirmed or refuted:

1. The change is controlled by `beta_l = Gamma/Gamma_crit,l`, with
   `Gamma_crit,l = 8 pi^2 |I_Omega,l| Pdot/(t_l B P^3)`. It is material for
   `|beta| ≳ 0.1` and reverses the imbalance at `beta > 1`.
2. For a neutron sink the photon luminosity follows `|1 − beta_l|^(8/7)` per
   channel, with the interaction contrast `N_X` near −1 at `beta ≈ 1`.
3. The chemical effect survives the process's own direct heating only if
   `q_dep ≲ a_RC/beta`, about a keV per destroyed baryon on this star.

The four coupled runs of S2 then verify that the time-dependent system has the
signs, slopes and ledger the map assumes. Statement 3 is why negligible direct
deposition must be treated as a hypothesis about the microscopic process and not
as a default.

## 11. Sources and their status

The source PDFs are tracked in [docs/references](../../../references/README.md),
supplied by the owner on 2026-11-12, with hashes in
`docs/references/SHA256SUMS`. An earlier state of this section rested on web
renderings; every equation, table value and quoted statement cited in this
document has since been checked against the PDF text at the pages below. The
check covered the cited passages, not a full reading of each paper.

Paths are under `docs/references/`. Page numbers are PDF pages.

| Short name | File | Checked items |
|---|---|---|
| FR2005 | `rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | Abstract (p. 1): quasi-equilibrium luminosity depends only on spin-down power; J0437−4715 prediction 20% low; superfluidity should raise it. Eqs. (64)–(68), p. 15, including the 5/8 and 3/8 split. Eq. (81), p. 17 |
| BGZ2022 | `bnv/zakeri/2022-Berryman-Gardner-Zakeri-Neutron-Stars-with-Baryon-Number-Violation-Probing-Dark-Sectors.pdf` | Eq. (29) and its paragraph, p. 18. Table 2 last row, p. 20 |
| BGZ2023 | `bnv/zakeri/2023-Berryman-Gardner-Zakeri-Macroscopic-Limits-on-Neutron-Star-Baryon-Loss.pdf` | Eq. (2.10), p. 12. §3 passages, p. 15. Table I, p. 22. Eq. (4.8), p. 25. Eq. (5.10), p. 33. Eqs. (5.13), (5.17) and the thin-shell statement, pp. 38–40. §6 passage, p. 49. §8 passage, p. 56 |
| MPR2021 | `bnv/2021-McKeen-Pospelov-Raj-Smoke-and-Mirrors-Neutron-Star-Internal-Heating.pdf` | Three heat contributions and eq. (9), p. 2. Eq. (11) and `<mu_n> = 68 MeV`, p. 3 |
| MPR2020 | `bnv/2020-McKeen-Pospelov-Raj-Cosmological-and-Astrophysical-Probes-of-Dark-Baryons.pdf` | `Delta E ≈ 100 MeV` per decay and applicability to other exotic modes, p. 6. Eq. (20), p. 8 |
| SNOplus2022 | `bnv/2022-SNOplus-Improved-Search-for-Invisible-Modes-of-Nucleon-Decay.pdf` | `tau(n → inv) > 9.0e29 yr`, abstract p. 2 and Table V p. 6 |
| RFJ2006 | `rotochemical/2006-Reisenegger-Fernandez-Jofre-Internal-Heating-and-Thermal-Emission-from-Old-Neutron-Stars.pdf` | Not cited for any statement here. Skimmed only |

FR2005 is byte-identical to the macOS-catalogue file: its SHA-256 equals the
value in ADR-0014 §1.

**Still missing.** `R2006`, the electrostatic-potential formalism paper that
supersedes FR05 eqs. (9), (12), (13), (54)–(56), was not among the supplied
files. The supplied 2006 paper (RFJ2006, astro-ph/0610955) is a different one.
R2006 is used here only through ADR-0014, so no statement in this document
depends on reading it; it is still needed in the tracked set. The macOS
`catalog.tsv` and literature README remain unavailable on EKU.

**Corrections made after reading the PDFs.**

1. The baryon–dark-baryon resonance of BGZ2023 was wrongly listed as a Fermi-
   surface mechanism. It selects one momentum inside the Fermi sphere and the
   authors find it irrelevant. §6.3 is corrected.
2. The β-disequilibrium and Urca passage is in BGZ2023 §6 (p. 49), not §8.
3. Section references for the in-medium statements are corrected in §6.2.

Limits on `Gamma`: 3.7e-13, 7e-12 and 1.1e-12 yr⁻¹ (BGZ2022 Table 2:
J0737−3039A/B, B1913+16, J1713+0747) and 1.8e-10, 2.0e-11, 4.0e-13 yr⁻¹
(BGZ2023 Table I: J0348+0432, J1614−2230, J0737−3039A/B), all 2σ. They show
that `Gamma ~ 1e-13–1e-12 yr^-1` is not excluded by timing for those systems.
They are not limits on the pulsars or the fixture star discussed here.

**What BGZ2023 says about thermal effects, and how this study treats it.**
The earlier draft of this section reported an apparent tension. Reading the PDF
and the owner's clarification remove it. The three passages are:

- p. 15: the authors "expect the deviations from a degenerate state at
  β–equilibrium due to BNV should be negligibly small" and "leave a more
  detailed study of possible thermal effects on neutron stars from BNV to future
  work". The context is the justification for a cold degenerate equation of
  state.
- p. 49: "BNV decays lead to β-disequilibrium, which leads to neutrino cooling
  via (direct and modified) Urca processes, which impact how the energy released
  in the decays is deposited back into the SM fluid. We cannot say at this time
  whether a constraint from thermal heating is more stringent or not."
- p. 56: the `O(E_F)` energy deposited on neutron removal "does not characterize
  a net gain of thermal energy by the system. For that we need to include the
  cooling from Urca processes", because the thermal evolution "is determined by
  the total rate of heating or cooling".

The owner, a co-author, stated on 2026-11-12 that these passages record an open
question: the paper did not calculate whether the Urca cooling induced by BNV
can compete with the direct heating, and it makes no claim about the size or
sign of the net effect. **They are not to be relied on as results.** Questions
about that paper's intent go to the owner.

This study is consistent with them and does not use them as evidence:

- The imbalances found here are tens of keV against a baryon chemical potential
  of 868 MeV, a fraction of about 2e-5. That is negligible for the equation of
  state, as p. 15 expects, and it is still the quantity that sets rotochemical
  heating in old pulsars (FR2005).
- The net thermal effect of hole heating together with the induced Urca response
  is exactly the quantity in ledger (4)–(5). §6 class M1 gives this study's own
  estimate, on the fixture, that direct hole heating exceeds the chemical
  response by `1e4–1e5`. That estimate is to be tested by the map of §9, not
  assumed.

No trajectory, root solve, predictor or qualification campaign was run for this
proposal. The numbers in §4.4 are closed-form arithmetic on existing governed
values and are not retained evidence.
