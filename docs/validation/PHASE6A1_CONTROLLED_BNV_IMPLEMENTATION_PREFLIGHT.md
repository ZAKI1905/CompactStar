# Phase-6A-1 controlled abstract-BNV implementation preflight

**Status:** CANDIDATE IMPLEMENTATION CONTRACT / DOCUMENTATION ONLY / NO BNV
IMPLEMENTATION / NO BNV TRAJECTORY / NO BNV BASELINE.

**Review correction status:** **PHASE-6A-1 REVIEW CORRECTIONS APPLIED —
ACTUAL-POTENTIAL LEDGER AND ZERO-SPIN QUALIFICATION CLOSED — CANDIDATE READY
FOR BOUNDED INDEPENDENT RE-REVIEW.** The independent review of
`f9b9029600bab47c61a527672da45516bbe218a6` found **0 BLOCKING / 2 MATERIAL /
5 NONBLOCKING / 7 NOTE**. This document remains a candidate implementation
preflight: it is not owner-accepted, canonically integrated, or implementation
authority.

This plan governs the first abstract ordinary-neutron-disappearance experiment on
the already-governed Structure-1 / Phase-5D free-gas machinery. It does not select
a physical BNV rate, model `n -> chi gamma`, begin A18, add superfluidity, add
Regime-II or MixedStar thermal physics, or authorize a trajectory. Every numerical
trajectory remains blocked by the pretrajectory gates in section 15.

## 1. Authenticated entry and change boundary

`PHASE6A1_ENTRY_SHA = 15a224c804605877103cf8464d48a89781097f74`.

Before this file was written, the canonical checkout was clean and all four
identities were equal:

```text
canonical checkout HEAD =
local master            =
origin/master           =
live refs/heads/master  =
15a224c804605877103cf8464d48a89781097f74
```

The fresh planning branch and worktree are:

```text
analysis/phase6a1-controlled-bnv-implementation-preflight
/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-bnv-implementation-preflight
```

The review-correction entry is
`PHASE6A1_REVIEWED_SHA = f9b9029600bab47c61a527672da45516bbe218a6`.
The permanent diff for this correction is restricted to this document. No production
source, test, baseline, EOS/data file, literature file, or Phase-5 artifact may
change. No ADR-0016 is proposed: ADR-0015 already decides the scientific and
architectural seam; the remaining choices are implementation planning within that
accepted seam.

## 2. Authority read and authenticated upstream artifacts

The following were read completely in authority order: `AGENTS.md`,
`GOVERNANCE.md`, ADR-0015, the Phase-6A-0 ratification, preflight and Cowling
diagnostic, ADR-0011/0013/0014, the Phase-5B integration and ratification records,
the Phase-5C integration and ratification records, the Phase-5D implementation,
provenance-preserving final report, independent-review/ratification and integration
records, `SCIENTIFIC_INVARIANTS.md`, `MODERNIZATION_ROADMAP.md`, and
`CURRENT_ARCHITECTURE.md`.

ADR-0015 is **ACCEPTED / HUMAN-RATIFIED / CANONICALLY INTEGRATED**. Existing BNV
code is architectural evidence only. It cannot override that contract. The current
architecture records that the legacy channels are compiled but unexercised and do
not implement ADR-0015
(`docs/architecture/CURRENT_ARCHITECTURE.md:1`-`11`).

The eleven governed files under `tests/baselines/` were authenticated at entry:

| Governed artifact | SHA-256 |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `90d607519cbdf3c4a0bf6ef50cc8fd22a8526b5db0354dc319e96854da29041d` |
| `grid_convergence_cmf_1p6_debug.tsv` | `b48519c3e948e9979a385d19facee2777d15955eeb8711b4bdd46b81fef74741` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `d5b753932c0523e67a7f25b460c7494bec1a006a8d01c9e43124cb2e78f0720f` |
| `hartle_I_dscmf1_debug.tsv` | `034ecddbd9bd847650429d7dc87d0331ec9e87aca3862ff87594e4bff5b707dd` |
| `hartle_monopole_dscmf1_debug.tsv` | `caaa0ac0d3219cda0a9fb518b27688afc23c6cdad1ec76a2bcd7359614a8d4e8` |
| `passive_cooling_cmf_1p6_debug.tsv` | `8fef2314673fceb939f859612f4befe94117115d6d6b3ad0dcc59d1faa68c9f9` |
| `phase5b_structural_response.json` | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` |
| `phase5c_chemical_coefficients.json` | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| `phase5d1_controlled_evolution.json` | `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

The Phase-5D manifest defines exactly 33 protected upstream paths and checks the ten
pre-Phase-5D baselines plus the Phase-5B/5C special hashes
(`tests/rotochemical/manifest.py:31`-`:82`). The eleventh, Phase-5D, hash above is
the canonically governed controlled benchmark. All 33 protected paths and all
eleven baseline bytes are immutable inputs to the future implementation.

The independent-review findings are reconciled as follows:

| Finding | Correction location |
|---|---|
| Material C-1 actual-potential direct ledger | sections 5, 8, 12, 14-17; R18/R20 and M21 are nonzero-`eta` requirements |
| Material C-2 spin-off qualification | sections 5, 10, 14-15; identity-bearing zero spin, `AnalyticControl` transfer, BA10a/BA10b |
| N-1 error-aware depletion | section 11 and BA13 |
| N-2 one consistent relativistic P1 convention | section 8 and BA7 |
| N-3 governed `Z` uncertainty | sections 7 and 11; `u_Z=E_Z`, no invented `Z` validation envelope |
| N-4 drive feasibility | section 13; infeasible high-`xi` root target removed |
| N-5 reached linear QSS | O16 and BA16 |
| Remaining notes | O17/BA17 B1 bound, frozen-`Z` deferral, distinct conversion sites, trial-state validity enforcement, legacy R-b wording, uniform-sink normalization, `E_F,kin` notation and architecture citations |

## 3. Governing mathematical contract

The ordinary charge-neutral source basis is

```text
y = (N_n, N_e, N_mu),       b = (1,1,1)^T,       B = b^T N_y,

    [-1 -1]
L = [ 1  0],                b^T L = 0,
    [ 0  1]

t = (partial N_y^eq / partial B)_Omega,           b^T t = 1,
Bdot = b^T S_y,
Sigma_y = S_y - t Bdot,
sigma = P Sigma_y = (Sigma_e,Sigma_mu)^T,
Sigma_y = L sigma.
```

There is exactly one physical chemical seam: `{Bdot, sigma}`. The chemical RHS is

```text
eta_dot = -Z(R + sigma) + 2 W Omega OmegaDot
```

for frozen coefficients, with `+ Zdot Z^-1 eta` only in a future consistently
sliding-coefficient implementation. `G_y` remains the upstream authority used to
construct `Z` and may construct diagnostic `k`; neither raw `G_y S_y` nor
`k = G_y b/(b^T G_y b)` may drive physical BNV chemistry.

For the first neutron-sink fixture,

```text
S_y = (Bdot,0,0)^T,    Bdot < 0,
sigma_e  = -t_e Bdot = t_e |Bdot| > 0,
sigma_mu = -t_mu Bdot = t_mu |Bdot| > 0.
```

With spin off and `eta=0`, both channels of `eta_dot=-Z sigma` must have the
ratified negative sign. Structure-1 long-digit values are arithmetic oracles only:

```text
t = (0.9657700849496014,
     0.030852171225661786,
     0.0033777438247248118)

absolute propagated t budgets =
    (1.4428870e-7, 4.3530477e-9, 2.1785265e-9)

(eta_dot_e,eta_dot_mu)/|Bdot| =
    (-1.4302859054377823e-55,
     -3.628096720520591e-55) MeV/count.
```

They support only the rounded fixture claims `(-1.4303e-55,-3.6281e-55)
MeV/count`; they are not physical precision. Provenance and the ratio-error formula
are recorded in the accepted preflight
(`docs/validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md:403`-`:435`).

## 4. Current-code audit and disposition

The audit covered every tracked BNV-named file, `StateTag::BNV`, the Phase-5D
pipeline, legacy thermal/heating code, BNV sequence/rate classes, and MixedStar BNV
analysis. CMake compiles all legacy Microphysics BNV sources
(`CompactStar/Microphysics/BNV/CMakeLists.txt:1`-`:10`) but repository search finds
no `main/` or `tests/` caller of a legacy BNV channel. The Chem BNV driver is empty
and has no compiled source (`CompactStar/Physics/Driver/Chem/CMakeLists.txt:1`-`:11`).

| Existing object/path | Classification | Finding and Phase-6A-1 rule |
|---|---|---|
| `Physics/State/BNVState.*`, `StateTag::BNV`, `StateVector::GetBNV` | **DO_NOT_USE_FOR_PHASE6** | The compiled state mixes a model-specific `eta_I` with a cached spin-down limit and describes both as ODE-compatible (`BNVState.hpp:13`-`:30`, `:84`-`:139`). It has no governed `B`, source, `t`, `sigma`, energy or currentness semantics. Leave it untouched; do not activate its tag. |
| `Physics/BNV.hpp` | **RETIRE_LATER** | Declares spin-down/heating-bound utilities but has no tracked definition or caller. It is not the ADR-0015 namespace contract (`Physics/BNV.hpp:60`-`:92`). The new layer may occupy `Physics/BNV/`; this header must not be silently treated as implemented authority. |
| `Driver/Chem/BNVSource.hpp`, `Driver/Spin/BNVSpinTorque.hpp` | **RETIRE_LATER** | Empty scaffolding. Do not fill these ambiguous files; use the explicit paths in section 14. |
| `Microphysics/BNV/Analysis/BNV_Analysis.*` | **DO_NOT_USE_FOR_PHASE6** | Integrates species-dependent static rate factors, then `Evolve` hard-codes `Gamma_BNV=1e-10/yr` and channel rescalings (`BNV_Analysis.cpp:49`-`:96`, `:217`-`:256`). It has no moving-reference source, matched control or thermal ledger. |
| `Microphysics/BNV/Analysis/BNV_Sequence.*` | **DO_NOT_USE_FOR_PHASE6** | Owns a physical-looking default rate `1e-10/yr`, exponential `B(t)`, mass cutoffs and a separate spin ODE (`BNV_Sequence.hpp:128`-`:149`, `:210`-`:217`, `:340`-`:362`; source `:580`-`:594`). It would mix changing-background/sequence physics into the first oracle campaign. |
| `Microphysics/BNV/Analysis/Decay_Analysis.*` | **RETIRE_LATER / DO_NOT_USE_FOR_PHASE6** | The header labels it “old and obsolete” and embeds particle masses/mixing (`Decay_Analysis.hpp:27`-`:37`, `:80`-`:120`). |
| `Microphysics/BNV/Internal/BNV_Chi.*` | **DO_NOT_USE_FOR_PHASE6** | Abstract base is already specialized to chi production, physical masses, mixing, EOS and pulsar rate limits (`BNV_Chi.hpp:27`-`:44`, `:161`-`:208`, `:246`-`:367`). It cannot own the abstract source contract. |
| `Channels/BNV_B_Chi_Transition.*`, `BNV_B_Chi_Combo.*`, `BNV_B_Psi_Pion.*` | **DO_NOT_USE_FOR_PHASE6** | Physical channel/mass/matrix-element implementations. They may be historical references in a later physical specialization only. |
| `Channels/BNV_B_Chi_Photon.*` | **DO_NOT_USE_FOR_PHASE6** | It is explicitly `B -> chi + photon`; it calculates a Fermi-hole term (`BNV_B_Chi_Photon.cpp:522`-`:617`) and combines it with a photon term (`:1722`-`:1751`, `:1881`-`:1898`). Hole plus independently deposited photon energy is a permitted R-b representation when one event, energy zero and terminal-fate ledger make it internally consistent; that decomposition is not intrinsically wrong. This class still cannot be reused because it is model-specific, owns incompatible event semantics and local literal conversions (`:646`-`:658`), has no governed moving-reference chemistry or Phase-6 provenance/currentness, and supplies no generic direct-energy interface. |
| `Extensions/MixedStar/DarkCore_Analysis.*`, `Core/MixedStar.*` | **DO_NOT_USE_FOR_PHASE6** | Separate two-fluid/dark-core sequence analysis with a fixed chi mass and a hard-coded rate in its limit analysis (`DarkCore_Analysis.cpp:16`-`:35`, `:82`-`:138`, `:266`-`:302`). MixedStar thermal physics is out of scope. |
| `Physics/SigmaOmegaRho_npemu.*` | **RETIRE_LATER** | Its only BNV claim is an “incomplete” file comment; it is an EOS model, not an ADR-0015 source/thermal implementation (`SigmaOmegaRho_npemu.hpp:27`-`:37`, `:54`-`:129`). |
| `Analysis::EquilibriumSequenceNumberDerivative` | **REUSE** | This is the exact Phase-5B sequence-derivative owner. It solves independent canonical neighbor stars, retains metadata/currentness and produces `dN_i/d epsilon_c` with numerical errors (`ParticleNumberResponse.cpp:349`-`:432`). A typed derived view constructs `t`; no recomputation through `G_y`. |
| `FrozenRotochemicalRunContext` | **WRAP; LEAVE SOURCE UNTOUCHED** | It owns current semantic `Z/W/Ltilde`, thermal bytes, spin owner, star/geometry and currentness (`FrozenRotochemicalRunContext.hpp:61`-`:129`). The BNV context composes it. The spin-off campaign deliberately uses the `AnalyticControl` base and re-imposes every relevant omitted production qualification gate through the wrapper, as specified in section 10. |
| `SecularEvolutionDriver` | **REUSE AS O1 ORACLE; LEAVE UNTOUCHED** | Its exact Phase-5D RHS updates one thermal and two chemical slots (`SecularEvolutionDriver.hpp:6`-`:29`). A new sibling driver owns the BNV-augmented RHS; the two drivers are never registered together. |
| `ChemicalImbalanceState`, `RotochemicalReactionResponse`, `RotochemicalThermalPower` | **REUSE** | They preserve typed channel order and the beta ledger. The existing thermal boundary converts MeV/s to erg/s once (`ChemicalImbalanceState.hpp:12`-`:35`; `RotochemicalReactionResponse.hpp:48`-`:54`). |
| `ScaledRKF45` | **REUSE; LEAVE UNTOUCHED** | It already requires the exact three-slot thermal/chemical layout and component-scaled tolerances (`ScaledRKF45.hpp:12`-`:18`, `:29`-`:58`). The wrapped Phase-5D context remains its currentness input; the BNV driver checks the BNV wrapper before and after every evaluation. |

No existing BNV production object satisfies ADR-0015. No audited legacy BNV code
currently routes raw `S_y` through `G_y`, because it has no governed chemical
response path at all; the risk is killed prospectively by BA5/M1. No existing BNV
object distinguishes `E_esc,fluid`, `E_esc,star`, `E_X` and terminal fate, or owns a
matched no-BNV control. The legacy photon class is the concrete double-count hazard.

## 5. Ownership answers A-L

| Question | Sole proposed owner |
|---|---|
| **A. `B(t)` / `Bdot`** | One immutable `OrdinaryMatterBnvHistory` interface returns an atomic sample containing `B`, `Bdot` and `S_y`. The prescribed history—not the ODE state, driver or observer—owns their mutual consistency. |
| **B. Actual ordinary-matter `S_y`** | The same history delegates to one typed `IOrdinaryMatterBnvSource`; its atomic sample is the only source consumed by projection and energy ledgers. The neutron-sink fixture is a test/validation specialization. |
| **C. `t`** | New `Analysis::EquilibriumBaryonTangent`, constructed only from a current `EquilibriumSequenceNumberDerivative` on the identical whole-star domain and sequence state. |
| **D. `sigma`** | New `Physics::BNV::MovingReferenceSource`; exactly one constructor evaluates `P(S_y-t Bdot)` and all identities/refusals. No caller may supply bare `sigma`. |
| **E. Chemical consumer** | New `ControlledBnvSecularDriver`, a sibling of the untouched Phase-5D driver, adds `-Z sigma` to the Phase-5D `eta_dot`. It consumes only typed `sigma`, never raw `S_y`. |
| **F. Direct partition / `P_dir`** | New `BnvDirectEnergyLedger`, using a separate `IDirectBnvEnergyPartition` and the same source event measure. It consumes the current typed chemical state, `mu_B^infinity(B0)`, governed `t`, source stoichiometry and product fate to reconstruct actual potentials; it owns GR integration and the exactly-once direct-power MeV-to-erg boundary. |
| **G. Product fate** | Separate immutable `ProductFateLedger`, keyed by source event/channel ID, with terminal weighted branches. Direct energy may read it but cannot invent or alter it. |
| **H. Matched control** | The validation run factory owns two distinct controls: BA10a is the governed spin-on `ControlledTrajectory` regression; BA10b uses the same identity-bearing `StaticZeroSpinHistory`, `AnalyticControl` base context and re-imposed BNV-wrapper qualification gates as the spin-off target, replacing only the BNV history/direct bundle by exact zero. |
| **I. Frozen validity/depletion** | `FrozenBnvValidityMonitor`, owned by the wrapped run context and parameterized by a pretrajectory `FrozenSensitivityCertificate`. It checks before/inside every RHS evaluation, including trial states, and at every emitted checkpoint; it fails closed before any out-of-domain point can be reported. |
| **J. Diagnostics/schema** | `BnvDiagnosticSnapshot` is assembled by the wrapped context from semantic owners; the candidate producer serializes only that schema. Neither observer nor serializer recomputes physics. |
| **K. Phase-5D object treatment** | Wrap `FrozenRotochemicalRunContext`; leave it, `SecularEvolutionDriver`, `ChemicalImbalanceState`, `RotochemicalReactionResponse`, `FrozenThermalSource`, prescribed-spin abstractions and `ScaledRKF45` unchanged. The spin-off wrapper deliberately constructs the base with `RunPurpose::AnalyticControl`, then re-imposes every relevant qualification gate omitted by that purpose. Add a sibling BNV driver; do not stack it with `SecularEvolutionDriver`. |
| **L. Generic vs fixture-specific** | History/source, domain/provenance, tangent, projection, product fate, energy ledger, validity monitor, driver and diagnostics are model-independent. Neutron-only `S_y`, P0/P1/P2, spin-off history, free-gas event weighting and Structure-1 oracle values are fixture-specific and live under `tests/bnv/`. |

## 6. `t` production contract

`EquilibriumBaryonTangent` consumes a shared current
`EquilibriumSequenceNumberDerivative`; it neither solves a new sequence nor reads
`G_y`. At `Omega=0`, it selects the exact named species axes and computes

```text
B_B,canonical = sum_baryons b_i B_i = B_n + B_p,
t_n = B_n/B_B,canonical,
t_e = B_e/B_B,canonical,
t_mu = B_mu/B_B,canonical.
```

`B_B,canonical` is the ADR-0011 / Phase-5B whole-star baryon derivative; the
reduced-basis expression `B_n+B_e+B_mu` is retained only as a charge-closure
cross-check. The two denominators must agree within the governed combined baryon
and charge numerical budget. Production `t` always uses the canonical baryon
denominator (`docs/adr/ADR-0011-particle-number-structural-response.md:46`-`47`,
`:156`-`:167`; `CompactStar/Analysis/src/ParticleNumberResponse.cpp:439`-`:467`).

The species map is explicit, never positional. The object retains:

- the complete Phase-5B `NumberMetadata`, contributing-star ownership, EOS bytes,
  central state, whole-star domain/surface policy, branch, stencil/step ladder and
  currentness dependencies;
- `Omega=0`, `B0`, sequence-state identity and a serialized domain identity;
- dimensionless `count/count` units;
- raw ratios, propagated absolute errors
  `(error_Bi+|t_i|error_BB)/(|B_B|-error_BB)`, denominator conditioning, and the
  raw closure residual;
- a closed representation storing `t_e,t_mu` and deriving
  `t_n=1-t_e-t_mu`. The closure adjustment must lie within the raw propagated
  budget; the accessor then satisfies `b^Tt=1` to at most two floating-point ulps.

Construction fails unless the domain is the same whole-star ordinary domain as the
source; the species and charge identities are exact; `B_B` passes the Phase-5B
conditioning rule; every dependency is current; and

```text
|b^T t_raw - 1| <= tau_t,
tau_t = sum_i error(t_i) + 32 epsilon_machine sum_i |t_i|.
```

For the Structure-1 oracle, `tau_t` is approximately `1.51e-7`. Every accessor
calls `RequireCurrent()`. A profile version, provider/revision, EOS bytes,
sequence-star identity, domain, central state or source-policy change makes the
object stale and causes refusal. A frozen object is explicitly `t(B0)`; a sliding
object would require a new sequence-state sample at each `B` and is not implemented
under this plan.

## 7. Generic source and the sole source-to-`sigma` transform

An `OrdinaryMatterBnvHistorySample` contains:

```text
epoch [s]
B [count]
Bdot_declared [count/s]
S_y = (S_n,S_e,S_mu) [count/s]
source/channel identity and revision
whole-star domain identity and sequence-state identity
event-measure identity
product-fate identity
direct-energy-partition identity
currentness token
```

Each event-measure branch also carries a dimensionless ordinary stoichiometric
vector `Delta y_event` and a nonnegative coordinate-time event measure. The atomic
sample must independently reproduce

```text
S_y = sum_branches Delta y_event dR_event^infinity
```

within `tau_S`. Neutron disappearance has `Delta y_event=(-1,0,0)` and therefore
positive event rate `R_event=-Bdot`. This keeps event energy positive while retaining
the signed ordinary source.

The generic API is not neutron-specific. If a microscopic adapter provides a
proton source, it must also satisfy `S_p=S_e+S_mu`; otherwise the charge-neutral
ordinary representation refuses it. The atomic sample requires

```text
tau_S = 32 epsilon_machine max(1 count/s,
                               |Bdot|, sum_i |S_i|),
|Bdot_declared - b^T S_y| <= tau_S.
```

For the first fixture the source is exactly `(Bdot,0,0)`, so the mismatch is zero.
The magnitude remains symbolic in this preflight. If the prescribed global
neutron sink is realized as a uniform proper sink proportional to neutron density,
its purely mathematical normalization is

```text
gamma = |Bdot| / integral_D e^Phi n_n dV.
```

This follows the ADR-0015 infinity-coordinate count-rate convention. `gamma` is a
fixture normalization only; it is not a physical decay rate.

`MovingReferenceSource` is the only class allowed to construct `sigma`:

```text
Sigma_y = S_y - t Bdot;
sigma = (Sigma_e,Sigma_mu);
lift = L sigma.
```

Its independent acceptance budgets are

```text
tau_Sigma,i = 64 epsilon_machine max(1 count/s,
                   |S_i|, |t_i Bdot|, |Sigma_i|)
              + |Bdot| error(t_i),

tau_b = tau_S + |Bdot| tau_t
        + 64 epsilon_machine max(1 count/s,
                                 sum_i |S_i|, |Bdot|),

|b^T Sigma_y| <= tau_b,
|Sigma_y - L sigma|_infinity <= max_i tau_Sigma,i.
```

It fails closed before returning a value if source/tangent domains, central-star
identity, revision, `B0`, or sequence state differ; either object is stale; source
charge closure fails; `Bdot` closure fails; tangent closure exceeds `tau_t`; or
either lift identity exceeds its budget. It does not accept precomputed `sigma`,
`G_y`, `k`, a bare coefficient array or an unqualified source callback.

Fixture claims that combine `Z` and `t` use, for output row `a`, the conservative
bound

```text
u_drive,a = sum_j [ |Z_aj| u_t,j
                   + u_Z,aj (|t_j|+u_t,j) ]
            + 64 epsilon_machine sum_j |Z_aj t_j|,
```

where `u_Z,aj = E_Z,aj` is the actual governed Phase-5C entrywise
`numerical_error.Z` for that entry. Phase-5C does not govern a structural
validation envelope for `Z`; none is invented or added here. Any future
Cowling/physical-model uncertainty remains a separately labelled uncertainty, not
`numerical_error` (`docs/validation/PHASE5C1R_STRUCTURAL_UNCERTAINTY_SEMANTICS_RATIFICATION.md:226`-`:246`;
`tests/baselines/phase5c_chemical_coefficients.json:121952`-`:121995`). No printed
oracle digit is used as a tighter bound.

For the governed Structure-1 artifact, the entrywise `E_Z` matrix is

```text
[[4.499198710915918e-61, 4.285424700647402e-61],
 [4.285424702045283e-61, 2.908175030741565e-59]] MeV/count.
```

These are numerical-error inputs to propagation, not physical-model uncertainty or
a validation envelope.

## 8. Direct energy and product fate

The direct-energy interface is a separate semantic owner from the source-to-`sigma`
transform, but it is **not independent of the chemical state**. It consumes the
source's immutable event measure, current typed `eta`, frozen/reference
`mu_B^infinity(B0)`, governed `t`, event stoichiometry, and a partition/fate object
with the same star, domain, revision, metric and energy-zero identities. It
constructs the actual ordinary-matter conjugate vector from ADR-0015:

```text
g_actual^infinity
  = mu_B^infinity b - (I - b t^T) P^T eta,

mu_n,actual^infinity
  = mu_B^infinity + t_e eta_npe + t_mu eta_npmu.
```

Thus `mu_n,actual` is state dependent even though the structural background and
`mu_B^infinity` are frozen. This correction is never reconstructed through a
changing-`B` `G_y` response; `G_y` remains upstream authority for `Z` only
(`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:212`-`:249`).

The production seam is generic. For an ordinary-matter event with signed count
change `Delta y_event`, the rest-mass-inclusive infinity-frame event ledger is

```text
Q_dir,event
  = -Delta y_event^T g_actual^infinity
    - E_esc,fluid^infinity
    + E_in,explicit^infinity.
```

`E_in,explicit` exists only when an independently owned external-inflow object is
present. For neutron disappearance, `Delta y_event=(-1,0,0)` and, without inflow,

```text
Q_dir,event = mu_n,actual^infinity - E_esc,fluid^infinity,
E_esc,fluid = E_esc,star + E_X.
```

All terms include rest mass and use one common energy zero. The product-fate
ledger uses terminal branches `PROMPT_ESCAPE`, `SM_THERMALIZATION`, `BOUND_INERT`
and `BOUND_INTERACTING`; branch weights are finite, nonnegative, sum to one within
`32 epsilon_machine`, and an event/channel ID may occur exactly once. Mixed fates
are separate weighted branches. Intermediate particles cannot be booked as a
second terminal fate.

If `dR_event^infinity=e^Phi Gamma_event dV` is the coordinate-time event measure,
the global power at infinity is

```text
P_dir,infinity [MeV/s]
  = integral e^Phi Q_dir,event,local dR_event^infinity
  = integral e^(2Phi) Gamma_event Q_dir,event,local dV.
```

The event measure independently integrates to the ordinary `S_y` sample. There
are two semantically distinct conversion sites in the future architecture:

1. beta chemical power is converted by existing
   `RotochemicalThermalPower::From`; and
2. direct BNV event power is converted by `BnvDirectEnergyLedger` after global
   MeV/s integration.

Each individual quantity crosses MeV-to-erg exactly once, using
`Rotochemical::MeVToErg` from `ChemicalImbalanceState.hpp:14`. What is forbidden
is omitted or double conversion of the same quantity. Partition implementations
return MeV/event only; sources return count/s only; no lower direct-energy layer
may return erg/s. The direct output is typed `DirectPowerErgPerSecond`
(`CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp:42`-`:47`).

The ledger exposes both equilibrium and actual direct power and must satisfy at
nonzero `eta`

```text
P_dir(actual) = P_dir(eq) + C_(MeV->erg) eta^T sigma.       (R18)
```

`P_dir(eq)` uses the equilibrium-sequence potential; `P_dir(actual)` uses the
current disequilibrated potential. BA7 exercises P2 at synthetic nonzero `eta`
and a second generic non-neutron stoichiometric fixture. A zero-`eta` comparison
cannot discharge this oracle (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:233`-`:250`).

Three representations of one event must agree: **R-a** is the open-system source
form, **R-b** is hole plus independently deposited event energy, and **R-c** is the
inclusive `-Delta y_event^T g_actual-E_esc,fluid` residual. BA8 requires the same
energy zero, event, event rate and terminal product bookkeeping for all three.
This is the explicit guard against extra hole heat, product-energy duplication and
mixed rest-mass conventions.

There is no independent Fermi-hole luminosity, generic PdV/gravity heat,
`-Echem_dot` heat, or unspecified positive efficiency. Mutation gates reject each.

### P0, P1 and P2 mathematical fixtures

| Fixture | Exact definition | Allowed claim |
|---|---|---|
| **P0** | At every cold event, `E_esc,fluid=mu_n,actual` in the same local energy zero. | `Q_dir,event,cold=0` and `P_dir,cold=0` exactly. It is a moving-reference/chemical transient oracle, not a physical channel and not a resolved finite-temperature sign result. |
| **P1** | At each radius, occupied neutron momentum is drawn from the normalized uniform occupied sea, `w(p)=3p^2/p_F^3`, `0<=p<=p_F`; `E_esc,fluid=E_n(p)` uses the same rest-mass-inclusive relativistic single-particle convention as the governed Track-R provider. Radial events use the source event measure. | The exact relativistic R10 average below is the production oracle. `(2/5)E_F,kin` is only its nonrelativistic-limit secondary check. The weighting is mathematical and is never generalized to a physical matrix element. |
| **P2** | `E_esc,fluid=E_esc,star=E_X=0` for every event; terminal fate is full retention in the mathematical ordinary-sector ledger. | `Q_dir,event=mu_n,actual`. It is the maximum-retention oracle only under the declared nonnegative escaping-energy convention and current cold fixture; it is not a universal physical bound after adding incoming energy or another energy zero. |

The first implementation deliberately omits finite-temperature direct-energy
corrections. Each run records `finite_T_direct_terms_included=false`, its weighting
class and an omitted-power bound. For smooth toy weights:

```text
P0 event floor = (pi^2/6) (k_B T)^2/E_F,kin,
P1 event floor = (pi^2/3) (k_B T)^2/E_F,kin.
```

A source concentrated within `O(k_B T)` of the Fermi surface instead receives an
`O(k_B T)` bound. “eV/event” is forbidden as a universal floor. A thermal sign may
be classified only when the named observable's complete uncertainty band excludes
zero after this omitted term is propagated through the same event measure.

For P1 the one permitted production convention is the relativistic uniform-sea
average

```text
E_F = sqrt(p_F^2 + m_n^2),
<E> = 3/[8 p_F^3]
      {p_F E_F (2 p_F^2 + m_n^2) - m_n^4 asinh(p_F/m_n)},
<mu-E> = E_F - <E>.                                      (R10)
```

Independent quadrature must match this expression within its predeclared
quadrature/roundoff budget. The nonrelativistic `2E_F,kin/5` limit is checked only
on a separately chosen `p_F/m_n << 1` fixture; it is not compared to the
relativistic production value at `1e-12`
(`docs/validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md:352`-`:365`).

## 9. Background and state-vector decision

Both structural policies were evaluated:

- **Sliding background:** scientifically possible later, but would require a
  current sequence-state owner at every `B`, recomputed `t,Z,C_star,Ltilde`, metric,
  photon/envelope data, `mu_n` event weights, a nonzero `Zdot Z^-1 eta` term and new
  ODE/refinement/energy-storage validation. It would mix changing-background
  physics into the first ADR-0015 projection test.
- **Frozen background:** tests the moving-reference and direct-energy contracts
  with the fewest new owners. Drift is bounded by an explicit certificate and a
  hard depletion stop.

**Decision: frozen background.** Operationally, `B(t)` is the cumulative ordinary
baryon count declared by the prescribed history,

```text
B(t)=B0+integral_0^t Bdot(t') dt',
```

and is used only for source consistency, depletion/currentness diagnostics and
the stop condition. It does not select a new star or alter frozen `Z,t,C_star,
Ltilde`, metric/structure, `mu_B^infinity(B0)`, equilibrium
`mu_n^infinity(B0)`, support or surface data. Those are frozen
structural/reference quantities. The current disequilibrated potential is not
frozen away:

```text
mu_n,actual^infinity
  = mu_B^infinity(B0) + t_e eta_npe + t_mu eta_npmu
```

within the frozen quadratic contract. The depletion certificate measures drift
of the structural/reference quantities; every direct-energy evaluation still
applies the current `eta` correction.

`B` is **not** a new ODE degree of freedom in the first implementation. The existing
Phase-5D state remains exactly

```text
(ln(Tinf/1e8 K), eta_npe,infinity [MeV], eta_npmu,infinity [MeV]).
```

The prescribed history supplies analytic `B(t),Bdot,S_y` atomically; BA1 checks
`B(t)-B0-integral Bdot dt=0` for constant and independently quadrature-checked
nonconstant test histories. This avoids activating the ambiguous `BNVState` while
preserving exact provenance and source consistency. Adding a B state is reconsidered
only when the source depends on accumulated B in a way that cannot be represented
by a current prescribed history, or when the background slides.

## 10. Matched control, chemical energy and sign observables

The first decisive BNV experiment is spin off. Current Phase-5D qualification does
not permit that through `RunPurpose::ControlledTrajectory`: production construction
requires `PrescribedDipoleHistory`. Existing zero-spin analytic tests instead use
`RunPurpose::AnalyticControl`, which omits several literal production-qualification
checks (`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:45`-`:103`;
`tests/rotochemical/coupled_oracles.hpp:53`-`:67`, `:93`-`:94`).

The future implementation therefore adds an identity-bearing
`StaticZeroSpinHistory` under the existing `PrescribedSpinHistory` abstraction. It
returns identically

```text
Omega = 0,
OmegaDot = 0,
identity = "static zero spin Omega=0 rad s^-1 OmegaDot=0 rad s^-2; Phase-6A-1 controlled abstract BNV".
```

Its provenance/currentness token is retained by the wrapper. Anonymous lambdas and
untracked fixture state are forbidden. `FrozenControlledBnvRunContext` constructs
the underlying Phase-5D context with `RunPurpose::AnalyticControl`, then fails
closed unless it re-imposes these omitted qualifications at its own boundary:

1. the exact `StaticZeroSpinHistory` identity and current provenance;
2. enabled process selection exactly `{Me,Mmu}`, with `De,Dmu` disabled;
3. metric identity exactly
   `"qualified Structure-1 radial80000 canonical nu/lambda"`;
4. normalization identities exactly
   `"predeclared mathematical benchmark SMe=1e-51 erg cm^-3 s^-1 K^-8"` and
   `"predeclared mathematical benchmark SMmu=2e-51 erg cm^-3 s^-1 K^-8"`.

Those literals are the strings production compares today
(`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:84`-`:103`;
`tests/rotochemical/coupled.hpp:11`-`:17`). The underlying context is intentionally
entered through `AnalyticControl` to permit zero spin. Every relevant production
qualification omitted by that purpose is transferred and re-imposed explicitly by
the BNV wrapper, so the net qualified experiment does not silently weaken the
governed Phase-5D fixture contract.

Two controls have different jobs:

- **BA10a — governed spin-on zero-BNV regression** uses
  `RunPurpose::ControlledTrajectory` with the existing
  `PrescribedDipoleHistory`, zero BNV source, and the governed Phase-5D controlled
  trajectory. It checks RHS identity, state order, `Z/W/Ltilde`, thermal ledger
  and the predeclared scaled-state criterion.
- **BA10b — spin-off zero-BNV matched control** uses
  `RunPurpose::AnalyticControl`, the same `StaticZeroSpinHistory` and the same
  re-imposed wrapper gates as the BNV target. It replaces only the BNV source and
  direct-energy bundle by exact zero. The untouched `SecularEvolutionDriver` under
  that same zero-spin history is the ordinary chemical/thermal reference path.

Thus every target still has a **PASSIVE SAME-INITIAL-B0 NO-BNV STAR**, while no
spin-off run is incorrectly compared to the spin-on governed trajectory as though
their physical histories were identical.

The diagnostic energy remains

```text
Echem = 1/2 eta^T Z^-1 eta,
Echem_dot = -eta^T(R+sigma)        (frozen Z, spin off),
BNV filling = -eta^T sigma,
beta release = eta^T R.
```

The implementation reports these components and checks the identity, but never
adds `-Echem_dot` to the thermal RHS. Beta thermal power remains
`DeltaP_beta=eta^T R-DeltaLnu`.

Diagnostics remain separate:

```text
LH, DeltaLnu, DeltaPbeta, Lnu_eq, Lnu_full,
P_dir, Lgamma, Pnet.
```

For the governed non-superfluid functions, the modified-Urca incremental sign root
is `|xi|=4.9097100289`; direct Urca is `|xi|=4.7870134733`. DU is disabled in the
first fixture unless a separately governed configuration enables it; this plan does
not. No run is called quasi-steady unless a measured `tau_relax(T,xi)` is much less
than the relevant source/thermal/depletion evolution time. Reports distinguish
linear QSS, freeze-out, drive-dominated transient and large-`|xi|` QSS.

No “heating”, “cooling” or “near zero” label may come from one RHS term. It must name
one of `DeltaTinf(t)`, `DeltaTsurface_inf(t)`, `DeltaLgamma(t)`, `DeltaU_th(t)` or
`integral_[t0,t1] DeltaP dt`, the exact time/window, the matched control and the
combined numerical/finite-T/partition/frozen-background error band. If the band
contains zero, the label is **SIGN UNRESOLVED**.

## 11. Frozen-background and depletion certificate

No trajectory may cross the first failed condition below. The absolute hard ceiling
for the first controlled campaign is

```text
|DeltaB|/B0 <= 1.0e-6.
```

This is a validity ceiling, not a selected drive or claim that all frozen quantities
are accurate at that depletion. Before any trajectory, a test-only sensitivity tool
must solve 21 independent equilibrium stars on the same authenticated stable branch,
uniformly spaced in `DeltaB/B0` over `[-1.0e-6,0]`, including `B0`. Target `B` is
obtained by a bracketed canonical sequence solve, not by altering the TOV equations.
For every target, both the final bracket width and target residual must satisfy

```text
|B_solved-B_target| <= tau_B,target,
DeltaB_bracket <= tau_B,target,
tau_B,target = max(u_B,solved,
                   64 epsilon_machine max(|B_target|,|B_solved|)).
```

`u_B,solved` is the propagated Phase-5B whole-baryon numerical error for that
solved star. This rule is fixed before the certificate is produced and cannot be
tuned to obtain a desired sensitivity result.
At every point it recomputes each quantity through its governed owner. The resulting
`FrozenSensitivityCertificate` retains all source bytes, star/profile revisions,
domain identities, numerical errors and the complete sample table.

For a scalar nonzero quantity, drift is

```text
D_X(B) = |X(B)-X(B0)| / max(|X(B0)|, absolute_scale_X).
```

Matrices use a row-scaled infinity norm; profiles use the maximum scaled nodewise
norm on the source support; exact-zero/near-zero quantities use an absolute scale.
The certificate forms a conservative monotone envelope by taking the maximum of all
sampled drifts up to each depletion point. Its smoothness test is error aware. For
each scalar/component or declared matrix/profile norm, fit
`X(x)=a+s x`, `x=DeltaB/B0`, by uncertainty-weighted least squares on the outer
certified half `x in [-1e-6,-5e-7]`. Propagate each point's governed numerical
error and the fit covariance to a residual uncertainty `u_res,i`. Let `T_X` be the
predeclared frozen threshold in the table below, expressed in the same norm. Require

```text
|X_i-(a+s x_i)| <= 3 u_res,i + 0.10 T_X
```

at every fitted point, and include the largest accepted residual as a conservative
nonlinearity term in the monotone envelope. The `0.10 T_X` term is the fixed allowed
physical linearization error; the full drift plus numerical and nonlinearity terms
must still remain below `T_X`. If the fit is rank deficient, any uncertainty is
nonfinite, or any criterion fails, the certificate refuses the full window. Raw
one-sided-slope percentages are forbidden when the drift is comparable to numerical
error. No implementation agent may enlarge, refit selectively or repair the window
after observing a trajectory.

The following limits are predeclared. `u_X` is the governed numerical error
propagated into the same norm. A condition passes only if its measured envelope is
at most the listed threshold:

| Frozen quantity | Sensitivity/recomputation owner | Exact Phase-6A-1 threshold |
|---|---|---|
| `t_n,t_e,t_mu` | `EquilibriumSequenceNumberDerivative` -> `EquilibriumBaryonTangent` | componentwise `abs(delta t_i) <= max(5 u_ti, 1e-4 abs(t_i))` and both lepton lower bounds remain positive |
| all four `Z` entries | current `GlobalChemicalNumberResponse` / `ChemicalImbalanceResponse` | row-scaled infinity drift `<= max(5 E_Z,row,1e-4)` using only governed `numerical_error.Z`; channel order unchanged; no nonexistent `Z` validation envelope |
| `C_star(T)` | `StarContext::HeatCapacityStar_Tinf` with current thermal/geometry owners | maximum relative drift over every sampled trajectory-temperature knot `<=1e-4` plus `5u_C` |
| enabled `Ltilde` coefficients | `GlobalUrcaChannelCoefficient` | each nonzero coefficient relative drift `<=1e-4` plus `5u_L`; enabled/disabled process mask identical |
| equilibrium/reference `mu_B^infinity(B)` and, where required, equilibrium `mu_n^infinity(B,r)` plus P0/P1/P2 structural event averages | Track-R local thermodynamic provider on each equilibrium star plus direct-energy integrator | source-support profile/event-average drift `<=1e-4` plus `5u_mu`; current `eta`-dependent actual-potential correction is evaluated separately and never omitted |
| `N_n,N_e,N_mu` | Phase-5B `ParticleNumbers` | each relative drift from `B0` `<=1e-4`; no small-muon exception |
| species support/domain | Track-R active-chart/EOS adapter | exact same active species and source support; no onset, threshold or domain crossing |
| metric/structure (`r,m,nu,lambda`) | canonical TOV/profile/geometry owners | source-support row-scaled profile drift `<=1e-4` plus `5u_structure` |
| radius, surface gravity, envelope quantities and `Tsurface_inf(Tinf)` | canonical surface and existing envelope owners | each scalar/table-knot drift `<=1e-4` plus `5u_surface` |

The Structure-1 estimate `d ln N_mu/d ln B approximately 54` means the hard
`1e-6` ceiling predicts an already material `approximately 5.4e-5` relative muon
drift. Therefore `N_mu`, `t_mu`, the muon row/column of `Z`, muon support and Mmu
`Ltilde` are measured separately and receive no shared proxy or relaxed threshold.

At runtime the monitor reports `DeltaB/B0`, each direct history ratio
`DeltaN_i/N_i` when available, and every certificate-envelope utilization
`U_X=drift_bound/threshold`. The frozen experiment is valid only while

```text
|DeltaB|/B0 <= 1e-6,
max_X U_X <= 1,
all source/domain/currentness identities remain exact.
```

Enforcement occurs before/inside every RHS evaluation, including integrator trial
states, and again at every emitted checkpoint/sample. A trial-state violation fails
closed immediately; no reported point may exist beyond the first failed check.
The untouched `ScaledRKF45` has internal accepted-state validation but no public
BNV accepted-step callback (`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp:38`-`:57`),
so this plan does not claim accepted-step-only BNV monitoring. Such semantics would
require a separate future integrator seam. A “valid” trajectory cannot be truncated
after the fact to hide a missed failure; M20 proves the RHS/checkpoint stop is active.

## 12. Analytic and bookkeeping oracles

The future implementation must expose these named tests before any BNV trajectory:

For the closed spin-off frozen experiment, BA11 uses the explicit ADR-0015 R20
finite-interval residual

```text
R20_[ti,tf] =
    [E_eq + C_(MeV->erg) E_chem + U_th]_f
  - [E_eq + C_(MeV->erg) E_chem + U_th]_i
  + integral_ti^tf dt [
        L_out,fluid^infinity
      + L_nu,full^infinity
      + L_gamma^infinity
      + L_other^infinity ],

Eeq_dot = C_(MeV->erg) mu_B^infinity Bdot,
L_out,fluid^infinity = L_esc,star^infinity + J_X^infinity.
```

An explicitly owned external input, if ever present, is subtracted as
`- integral P_in^infinity dt`; it may not be hidden in direct heat. For this
controlled closed fixture it is absent. The acceptance is fixed before any
trajectory:

```text
|R20_[ti,tf]| /
max(1 erg,
    |Delta E_eq|,
    |C_(MeV->erg) Delta E_chem|,
    |Delta U_th|,
    integral dt sum |terminal luminosity components|)
<= 2e-4,
```

with the independently predeclared quadrature error and omitted finite-temperature
floor also reported; neither may be tuned after the run. `Echem_dot` remains a
reservoir diagnostic and is never inserted as a thermal source
(`docs/validation/PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md:554`-`:603`).

| Oracle | Input | Required result and tolerance |
|---|---|---|
| **O1 ZERO SOURCE** | `S_y=0`, zero direct partition, otherwise exact Phase-5D context/state | New driver RHS equals the untouched `SecularEvolutionDriver` RHS bit for bit; all BNV diagnostics are exact zero. |
| **O2 PHYSICAL SLIDING NULL** | Synthetic `S_y=t Bdot` constructed from the same typed tangent; spin off, `eta(0)=0`, frozen coefficients | `sigma=0` within `tau_Sigma`; `eta_dot=0` within propagated `Z*tau_Sigma`; reaction-free integrated `eta(t)=0` within the component ODE tolerance. A future sliding background must reproduce this with current `t(B)`. |
| **O3 NEUTRON SINK SIGN** | `(Bdot,0,0)`, `Bdot<0`, Structure-1 tangent/Z, spin off, `eta=0` | `sigma_e,sigma_mu>0` after subtracting their propagated lower errors; both `eta_dot` upper error bounds are `<0`; rounded ratios match `(-1.4303e-55,-3.6281e-55) MeV/count` within combined Phase-5B/5C errors. |
| **O4 LINEAR TRANSIENT** | weak reactions disabled, frozen `Z`, arbitrary piecewise source, spin off | Independent quadrature gives `eta(t)=eta(0)-Z integral sigma dt`; for constant sigma, direct matrix multiplication gives `eta(0)-Z sigma t`. RHS comparison uses 64 ulps plus propagated coefficient/source error; integrated comparison has scaled ODE norm `<=1`. |
| **O5 BARYON NEUTRALITY** | generic finite synthetic source/tangent vectors | `abs(b^T(S_y-tBdot))<=tau_b`. |
| **O6 EXACT LIFT** | same as O5 | Componentwise `abs(L sigma-(S_y-tBdot))<=tau_Sigma`. |
| **O7 P0 COLD DIRECT NULL** | P0 at `T=0` convention | Event residual and integrated direct power are exact floating-point zero. Finite-T omission metadata/floor remain nonzero diagnostics. |
| **O8 P2 DIRECT POWER** | P2 event samples | Each residual equals `mu_n,actual` exactly; independent GR quadrature agrees within `max(1e-12 relative,10 times reported quadrature error)`. |
| **O9 NO DOUBLE COUNT** | fault adds a separate hole term to the complete ledger | BA9 comparison differs beyond its budget and the ledger/event-ID ownership check refuses it. |
| **O10 NO PdV HEAT** | fault registers synthetic PdV/gravity power | No authorized thermal-ledger slot exists; construction refuses before RHS accumulation. |
| **O11 NO ECHEM-AS-HEAT** | fault adds `-Echem_dot` to `Pnet` | BA8 independent ledger and BA9 mutation both fail. |
| **O12 UNIT BOUNDARY** | known MeV/s fixture and omitted/doubled conversion faults | Nominal erg/s divided by MeV/s equals `MeVToErg` within 8 ulps; both faults are detected. |
| **O13 R18 ACTUAL/EQUILIBRIUM** | P2 at synthetic nonzero two-channel `eta`, plus a second non-neutron stoichiometric event; common source/fate | `P_dir(actual)-P_dir(eq)-C_(MeV->erg)eta^T sigma=0` within the sum of direct quadrature, `t/Z/source` and 64-ulp arithmetic budgets. A zero-`eta` fixture does not count. |
| **O14 R-a/R-b/R-c EQUIVALENCE** | one neutron event with a nonzero hole and independently deposited component, plus product-fate branch | all three representations agree within `max(64 ulps,10 times independent quadrature error)` under the same event/rate/energy zero/fate; adding R-b pieces to R-c is detected. |
| **O15 R20 FINITE-INTERVAL CLOSURE** | independently reconstructed BA11 primitive outputs | the normalized residual defined above is `<=2e-4`, and remains inside the declared quadrature plus finite-T omission budget. |
| **O16 REACHED LINEAR QSS** | spin off, frozen coefficients, owner-accepted mathematical drive within the depletion ceiling | componentwise `abs(R_l+sigma_l)/max(abs(sigma_l),R_resolution_l)<=0.05` over the named terminal interval and measured `tau_relax<=0.10` times its elapsed forcing duration. Otherwise classify as transient/freeze-out, not QSS. |
| **O17 B1 NON-SUPERFLUID BOUND** | each enabled modified-Urca process at fixed `T`, both eta signs and independent polynomial extrema | `-DeltaP_beta,l <= 0.467659 Lnu_eq,M,l` with violation allowance no larger than `1e-10 max(Lnu_eq,M,l,1 erg/s)`. This is not a superfluid claim; separately qualified DU would use its governed bound. |

The chemical-energy check independently verifies

```text
Echem_dot = (-eta^T R) + (-eta^T sigma)
```

against a centered finite difference of `Echem` in a reaction/source-only RHS test.
The finite-difference discrepancy must be at most
`max(1e-10 relative,10 times the reported roundoff/step estimate)`. It is a
bookkeeping oracle, never a luminosity.

## 13. Mathematical drive selection policy

No physical BNV rate and no mathematical drive value is selected or ratified here.
The implementation task is divided into pretrajectory infrastructure/gates and a
later candidate-production step. After BA1-BA10b and BA13's sensitivity certificate
pass, but before BA11 or any trajectory, a committed/reviewable run card must state
the exact `Bdot/B0` values, durations, initial state, event support and P0/P1/P2
choices.

Values are chosen only from dimensionless mathematical targets:

- a reaction-free linear-transient case with maximum predicted `|eta|` at least
  100 times the analytic/numerical resolution and with `|DeltaB|/B0<=1e-7`;
- a coupled small-imbalance case whose predicted `max |xi|<=0.5`;
- a reachable linear-QSS case satisfying O16 while remaining within the certified
  frozen-depletion domain; and
- an optional drive-dominated transient that remains below the frozen ceiling,
  without a sign-root or large-`|xi|` requirement.

The reaction-free feasibility estimate is

```text
eta approximately -(Z sigma) Delta t_run,
|eta_l| approximately |(Z t_lepton)_l| |DeltaB|.
```

On the governed fixture, the reviewed arithmetic gives

```text
(Z t_lepton) B0 approximately (108.9,276.3) MeV.
```

Therefore `|DeltaB|/B0<=1e-6` caps the electron imbalance near
`1.1e-4 MeV` before reactions, and at `T=1e8 K` gives approximately
`|xi_e|<=0.0126`. The modified-Urca incremental root near `4.9097` is
unreachable at that temperature inside this campaign's depletion ceiling. It is
removed as a Phase-6A-1 target; root crossing belongs to a separately qualified
lower-temperature, wider-depletion or sliding-background experiment.

For each candidate strength, the run card computes the predicted depletion,
`|Z sigma|` build time, weak-reaction time and thermal resolution before acceptance.
It refuses strengths that violate the frozen envelope, fail to separate the target
signal from combined errors by a factor of ten, or map to a particle mass, matrix
element or physical lifetime. If no value meets all conditions, BA11 is blocked;
the agent does not relax the conditions. Exact values require separate owner
acceptance before trajectory generation.

A small-`xi` reachable-QSS run may use the reviewed order-of-magnitude guide of no
more than roughly `1e-11 yr^-1` fractional loss at `T approximately 1e8 K`, but that
number is not installed here as a drive and has no physical-rate meaning. The exact
mathematical value remains a future owner-accepted run-card field.

The first Phase-6A-1 implementation is **frozen `Z`**. ADR-0015's variable-`Z`
formula remains governed, but executable `Zdot Z^-1 eta` chain-rule validation is
explicitly deferred to the future sliding-background phase. This is a scope
deferral, not an unresolved ownership ambiguity.

## 14. Exact implementation map

All paths below are proposed future changes. None exists or changes in this
preflight.

| Class | Proposed path | New/modified | Role / semantic owner | Dependencies and state | RHS / units | Currentness, tests and controls |
|---|---|---|---|---|---|---|
| **PRODUCTION** | `CompactStar/Analysis/EquilibriumBaryonTangent.hpp/.cpp` | new | sole typed `t(B0)` derived view | current Phase-5B sequence derivative; no ODE state | dimensionless | BA2; stale/domain/axis/closure refusals; M2/M17/M18 |
| **PRODUCTION** | `CompactStar/Physics/BNV/OrdinaryMatterSource.hpp` | new | abstract atomic history/source/event-measure types | typed domain/revision/currentness tokens | `B` count, `Bdot,S_i` count/s | BA1; charge/Bdot/history consistency; M16/M19 |
| **PRODUCTION** | `CompactStar/Physics/BNV/MovingReferenceSource.hpp/.cpp` | new | only `sigma=P(S-tBdot)` owner | atomic source sample + typed `t` | `sigma` count/s; no RHS write | BA3-BA5; M1-M5/M17-M19 |
| **PRODUCTION** | `CompactStar/Physics/BNV/ProductFate.hpp` | new | terminal fate branches and once-only event identity | source event/channel identity | weights, energies remain in partition | BA1/BA7/BA9; M15/M16 |
| **PRODUCTION** | `CompactStar/Physics/BNV/DirectEnergyLedger.hpp/.cpp` | new | generic actual-potential reconstruction, direct partition, GR integration and direct-power unit boundary | source event measure/stoichiometry, product fate, metric, `mu_B^infinity(B0)`, current typed `eta`, governed `t` | returns typed `P_dir(eq/actual)` erg/s; no chemical RHS | BA7-BA9/O7-O15; M8-M10/M13-M16/M21 |
| **PRODUCTION / CONTROLLED-EXPERIMENT** | `CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp/.cpp` | new | identity-bearing spin-off history under `PrescribedSpinHistory` | owned currentness token and exact identity string | `Omega=OmegaDot=0` | BA10b construction/currentness and identity mutations |
| **PRODUCTION** | `CompactStar/Physics/BNV/FrozenBnvValidity.hpp/.cpp` | new | sensitivity certificate and runtime stop | frozen authorities and prescribed `B` | diagnostics/refusal only | BA13; M17/M18/M20 |
| **PRODUCTION** | `CompactStar/Physics/BNV/FrozenControlledBnvRunContext.hpp/.cpp` | new | wrapper/composition root; transfers the omitted `AnalyticControl` qualifications and assembles semantic snapshot | shared const Phase-5D context plus exact zero-spin owner, process/metric/normalization literals, `t`, history, projector, fate, partition, validity | evaluates full diagnostics; same 3-state vector | every dependency retained/current; BA1-BA14/BA16-BA17; qualification mutations fail closed |
| **PRODUCTION** | `CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp/.cpp` | new | sole augmented chemical/thermal RHS owner | wrapper; Thermal(1)+Chem(2); never co-registered with Phase-5D driver | adds `-Z sigma` to eta and `+P_dir(actual)/(Cstar*Tinf)` through full ledger | BA4-BA12; all RHS mutations |
| **PRODUCTION** | `CompactStar/Physics/BNV/BnvDiagnostics.hpp` | new | immutable diagnostic schema assembled from owners | wrapper evaluation only | fields/units in section 17 | BA14; serializer cannot recompute |
| **BUILD** | `CompactStar/Physics/BNV/CMakeLists.txt`, `CompactStar/Physics/CMakeLists.txt` | new/modified | compile/install only the new layer | no legacy BNV dependency | none | build and dependency audit |
| **TEST-ONLY** | `tests/bnv/controlled_neutron_sink_fixture.hpp` | new | neutron-only history and spin-off source fixture | Structure-1, symbolic/predeclared drive | exact `(Bdot,0,0)` | O1-O6; M1-M7/M17-M20 |
| **TEST-ONLY** | `tests/bnv/energy_partition_fixtures.hpp` | new | P0/P1/P2, R18/R-a-b-c and independent event quadrature | direct ledger inputs at zero and nonzero `eta` | MeV/event until direct boundary | O7-O15; M8-M16/M21 |
| **TEST-ONLY** | `tests/bnv/source_projection.cpp`, `tests/bnv/direct_energy.cpp`, `tests/bnv/thermal_ledger.cpp` | new | BA1-BA9 executable oracles | production public APIs | no artifact | analytic/refusal/mutation coverage |
| **VALIDATION TOOL** | `tests/bnv/frozen_sensitivity.cpp` | new | produces the 21-star sensitivity certificate before trajectories | governed sequence/EOS/chemical/thermal owners | static calculations only | BA13; no BNV evolution |
| **TEST-ONLY** | `tests/bnv/matched_control.cpp`, `tests/bnv/coupled_trajectory.cpp`, `tests/bnv/ode_refinement.cpp` | new | BA10a spin-on governed regression, BA10b spin-off matched control, BA11-BA12 target/refinement | purpose-qualified run factories, exact spin owners, same relevant state/output grid | trajectory only after gates/run card | M11/M12/M20/M21 and coupled checks |
| **TEST-ONLY** | `tests/bnv/qss_beta_bounds.cpp` | new | BA16 reached linear-QSS and BA17 non-superfluid B1 diagnostic | accepted mathematical run card and independent R1995/FR2005 polynomials | diagnostics only | O16/O17; no physical-rate claim |
| **VALIDATION TOOL** | `tests/bnv/produce_candidate.py`, `tests/bnv/compare_candidate.py` | new | deterministic candidate production/comparison | completed applicable BA1-BA17 evidence | JSON only | no baseline promotion |
| **BUILD** | `tests/CMakeLists.txt` | modified | registers named BA tests and qualified producer | existing authenticated-data guard | none | absence of data fails/excludes exactly as predeclared |
| **DOCUMENTATION** | `docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION.md` | new in implementation task | evidence/status report | exact commands, return codes, hashes | none | cannot claim ratification/integration |
| **CANDIDATE ARTIFACT** | `docs/validation/phase6a1_controlled_bnv_candidate.json` | new only after gates | non-governed candidate | schema/provenance below | diagnostic/trajectory values | produced only after applicable BA1-BA17 preconditions |

The legacy `Physics/BNV.hpp`, `BNVState`, empty drivers, Microphysics BNV sources,
MixedStar and every Phase-5 source remain untouched by the initial implementation
map. Any future need to alter a Phase-5 protected file is a plan mismatch and stops
the task for review.

## 15. Validation ladder BA1-BA17

The tolerances below are part of this contract. `epsilon_machine` is IEEE-754 binary64
epsilon. “Independent” means the oracle does not call the production function under
test. Dependency order is `BA1-BA10a/BA10b -> BA13 static certificate -> accepted
run card -> BA11-BA12/BA16 -> BA13 runtime-stop controls -> BA14-BA15/BA17`.
A failed prerequisite
blocks every dependent gate and no candidate value is retained; the numbers name
the governed subjects rather than permitting BA11 to precede the static BA13 gate.

| Gate | Input | Expected result / numerical tolerance | Independent oracle / independence level | Failure disposition |
|---|---|---|---|---|
| **BA1 — source semantics** | exact zero, neutron-only, charge-balanced mixed, nonconstant analytic histories; mismatched `Bdot`, charge, event IDs and stale revisions | atomic sample; `abs(Bdot-b^TS)<=tau_S`; analytic `B-B0-integral Bdot=0` for constant history and `<=max(64 ulps,10 quadrature errors)` otherwise; all invalid cases refuse before publication | direct scalar sums and analytic/independent Gauss-Legendre time integral | stop; source API or fixture is invalid; M16/M19 must fire |
| **BA2 — `t` provenance** | current Phase-5B Structure-1 derivative and mutations of every retained dependency | oracle values lie within their propagated absolute budgets; raw closure `<=tau_t`; closed accessor within two ulps; exact domain/state/currentness retained; every mutation refuses | direct extraction from governed Phase-5B JSON and independent ratio/error calculation | stop; no `t`, projection or trajectory; M2/M17/M18 must fire |
| **BA3 — moving-reference projection** | at least 100 deterministic synthetic finite source/tangent cases spanning signs/scales plus domain/currentness faults | O5/O6 within `tau_b/tau_Sigma`; nominal source unchanged; every mismatch refuses | standalone matrix/vector implementation of `S-tBdot`, `P` and `L` | stop; M3-M5/M18/M19 must fire |
| **BA4 — neutron-sink signs** | Structure-1 `t,Z`, spin off, `eta=0`, symbolic nonzero negative `Bdot` scaled out | positive `sigma_e,mu` with positive lower bounds; negative `eta_dot_e,mu` with negative upper bounds; rounded drive coefficients agree within combined propagated Phase-5B/5C errors | direct high-precision matrix product from authenticated JSON values | stop; neutron-sink specialization/channel order invalid; M3/M5-M7 must fire |
| **BA5 — sliding null / raw-G refusal** | `S=tBdot`, both signs and at least three scales; fault routes raw source via `G_y` or substitutes diagnostic `k` | O2 zero within projection/ODE budgets; raw-G and k faults reproduce a nonzero result and are detected; no production BNV dependency accepts `G_y` | analytic zero and accepted R14 negative oracle from ADR-0015 preflight | stop; physical seam is violated; M1-M4 must fire |
| **BA6 — reaction-free transient** | weak reactions disabled; constant and piecewise sigma; nonzero two-channel initial eta; spin off | O4 RHS to 64 ulps plus input errors; integrated scaled state error `<=1` using `D_i=atol_i+rtol abs(y_i)` | independent matrix multiplication and quadrature, never the driver | stop; M5-M7 must fire |
| **BA7 — P0/P1/P2 plus R18 direct energy** | same synthetic event measure/metric and all three partitions; P2 and a generic event at nonzero `eta`; invalid energy zeros/fates | P0/P2 exact identities; relativistic P1 R10 matches within `max(1e-12 relative,10qerr)`; its separately nonrelativistic fixture approaches `2E_F,kin/5`; O13 R18 closes within propagated direct/source/tangent plus 64-ulp budget; fate/energy closure and finite-T floors present | analytic R10 expression, standalone quadrature, direct matrix evaluation of actual `g`, and GR sum share no production-ledger code | stop; M8/M13-M16/M21 must fire |
| **BA8 — thermal/open-system ledger plus R-a/R-b/R-c** | arbitrary signed eta, rates, sigma, direct power, photons/neutrinos; one event represented three ways; beta sign-root and chemical finite-difference cases | exact named decomposition; `DeltaPbeta=LH-DeltaLnu`; `Pnet=Pdir(actual)+DeltaPbeta-Lnu_eq-Lgamma-Lother`; O14 representations agree within `max(64 ulps,10qerr)`; Echem identity within section-12 bound; roots agree to `5e-10` absolute in `xi` | independently coded scalar/event ledger and accepted imbalance polynomials; no production partition routine reused | stop; M8-M16/M21 must fire |
| **BA9 — double-count/unit/product-fate mutations** | every forbidden extra term/energy/fate/unit/actual-potential mutation | nominal output passes; every M8-M16/M21 changes or refuses beyond BA7/BA8 budget; zero aliases are tested at nonzero discriminating fixtures | fault-injection adapter outside production | stop; no coupled run |
| **BA10a — governed spin-on zero-BNV regression** | `ControlledTrajectory`, governed `PrescribedDipoleHistory`, zero BNV source, same governed Phase-5D inputs | new-wrapper/driver RHS identity where expected; exact state order and same `Z/W/Ltilde` and thermal ledger; repeat bytes identical and scaled state difference from governed Phase-5D trajectory `<=1` | untouched `SecularEvolutionDriver` and governed Phase-5D baseline producer on the same spin history | stop; Phase-6 wrapper perturbs ordinary Phase-5D behavior |
| **BA10b — spin-off zero-BNV matched control** | `AnalyticControl`, same `StaticZeroSpinHistory`, exact transferred process/metric/normalization gates; source/direct bundle OFF versus ON | wrapper construction refuses every relaxed-gate mutation; zero-BNV RHS equals untouched `SecularEvolutionDriver` under the same zero-spin owner; target/control differ only in BNV bundle | untouched Phase-5D driver on the same analytic context; identity/provenance mutations independent of BNV source | stop; zero-spin qualification or matched-control isolation failed; M11/M12 also exercised |
| **BA11 — coupled spin-off BNV+beta and R20 closure** | owner-accepted predeclared mathematical run card, spin-off target plus BA10b control, primitive ledger outputs | completes without refusal; all prerequisite identities hold at every output/RHS check; no frozen utilization exceeds one; sign labels include observable/window/error; explicit O15 R20 normalized residual `<=2e-4` with all terms and omitted floors reported | independently reconstructed `[E_eq+C Echem+U_th]` endpoints, terminal-luminosity quadrature and matched-control differences from primitive columns | stop; discard candidate; do not tune drive or residual after seeing result; M9-M12/M21 must fire |
| **BA12 — ODE refinement** | identical BA11 run with baseline Phase-5D tolerances `(rtol=1e-7,atol=(1e-12,1e-18,1e-18))` and refined `(1e-9,(1e-14,1e-20,1e-20))` | at every common output and for each state, `abs(y_base-y_refined)/(atol_base+rtol_base max(abs(y_base),abs(y_refined)))<=1`; ledger residual still `<=2e-4`; no step-collapse trend in final decade | refined run and exact shared checkpoints | stop; no candidate artifact |
| **BA13 — depletion/frozen budget** | 21-star certificate with target-`B` bracket evidence, then forced histories ending just below and above each stop | target brackets meet `tau_B,target`; uncertainty-weighted outer-half fits meet `3u_res+0.10T_X`; full monotone envelopes remain within thresholds; first trial/checkpoint over-limit refuses before serialization; runtime max utilization `<=1` | independently parsed sample/error/covariance table, independent refit and forced-threshold test | stop; M20 must fire; do not reduce or retune evidence window after a trajectory |
| **BA14 — output/diagnostics** | one RHS snapshot, matched pair and eligible trajectory | all required fields, units, identities, error/utilization fields and sign metadata present; serializer values equal owner snapshot bit for bit; independent recomputation passes BA8/BA11 budgets | schema validator and direct recomputation from primitive columns | stop; artifact incomplete or unauditable |
| **BA15 — governed regression protection** | clean source plus complete future diff | 33 protected paths exact, all eleven baseline hashes exact, Phase-5B/5C/5D focused regressions green, complete data-free and authenticated suites return zero, `git diff --check` clean; only predeclared implementation/test/docs/candidate paths changed | fresh-producer regressions and exact SHA-256 manifest, not saved return-code files | stop; restore/replan; no candidate acceptance |
| **BA16 — reached linear QSS** | spin off, frozen `Z`, owner-accepted purely mathematical small-`xi` drive inside certificate | O16 component balance `<=0.05` and measured `tau_relax/elapsed<=0.10` on named interval; otherwise no QSS label | rates and sigma recomputed independently from primitive columns and accepted polynomials | stop QSS claim; trajectory may remain a transient candidate if all other gates pass |
| **BA17 — non-superfluid instantaneous B1 bound** | each enabled Me/Mmu process at fixed T, both signs and extrema; separately qualified DU only if ever enabled | O17 MU inequality with `1e-10 max(Lnu_eq,1 erg/s)` allowance; process attribution retained | independent polynomial evaluation/maximization, not production thermal ledger | stop beta-bound diagnostic; no superfluid or generic claim |

The only values deliberately deferred are the mathematical drive magnitudes. Their
selection procedure and acceptance criteria are fixed in section 13; the exact
values must be committed and owner-accepted before BA11/BA16 starts. No tolerance
in BA1-BA17 may be chosen after observing a trajectory.

## 16. Mutation/falsifier matrix

Each mutation is implemented outside production against a nonzero discriminating
fixture. A mutation “fires” only when the named test fails or explicitly refuses for
the intended reason.

| Mutation | Required detector(s) |
|---|---|
| **M1 raw `G_y S_y` route** | `BA5_sliding_null_raw_G_refusal`: physical sliding source becomes nonzero under the mutant. |
| **M2 `k` substituted for `t`** | `BA2_t_oracle` and `BA5_sliding_null_k_refusal`; accepted fixture `t_e/k_e` and `t_mu/k_mu` distinctions are load-bearing. |
| **M3 wrong sign on `t Bdot`** | `BA3_projection` and `BA4_neutron_sign`. |
| **M4 omit `t Bdot`** | `BA3_projection` and exact BA5 sliding null. |
| **M5 omit one sigma channel** | BA3 exact lift, BA4 two-channel sign and BA6 transient. |
| **M6 transpose/cross-`Z` error** | BA4 high-precision product and BA6 independent full matrix transient; both cross entries nonzero. |
| **M7 wrong eta channel ordering** | BA4 named ratios and BA6 asymmetric eta/sigma fixture. |
| **M8 Fermi-hole added twice** | BA7 P0/P1 identities, BA8 R-a/R-b/R-c equivalence and `BA9_double_hole`; event/partition ID duplicate refuses. |
| **M9 `Echem_dot` added as heat** | BA8 independent `Pnet` and `BA9_echem_heat`. |
| **M10 PdV/gravity added as heat** | construction refusal plus `BA9_pdv_heat`. |
| **M11 `DeltaLnu` omitted** | BA8 root/ledger, BA10b nonzero-eta matched control and BA11 R20 closure. |
| **M12 equilibrium neutrinos double counted** | BA8 full ledger, O1/BA10a/BA10b zero-source reductions and BA11 R20 closure. |
| **M13 MeV-to-erg omitted** | O12 known-unit fixture and BA7 global power. |
| **M14 MeV-to-erg doubled** | O12 and BA8 `Pnet`; mutant ratio differs by `MeVToErg`. |
| **M15 `E_esc,fluid` confused with `E_esc,star`** | BA7 branch with nonzero retained `E_X`; closure `Efluid=Estar+EX` fails. |
| **M16 product fate double-booked** | BA1 duplicate event ID and BA7 fate-weight/terminal-energy closure. |
| **M17 stale `t`** | BA2 mutation of profile/revision/source bytes and BA13 runtime revocation. |
| **M18 source/`t` domain mismatch** | BA3 constructor refusal before any sigma publication. |
| **M19 `Bdot != b^T S_y`** | BA1 atomic sample refusal and BA3 no-call assertion. |
| **M20 ignored frozen-budget violation** | BA13 forced just-over-threshold history; absence of pre-serialization refusal fails the test. |
| **M21 equilibrium `mu_n` used at nonzero `eta`** | BA7/O13 nonzero-eta R18 identity and BA11/O15 finite-interval R20 closure. A zero-eta fixture is forbidden as its sole detector. |

In addition, a source dependency scan must show that the production BNV subtree has
no include or member dependency on `GlobalChemicalNumberResponse::Values()` or any
diagnostic-`k` constructor. `Z` is consumed only through the wrapped Phase-5C/5D
semantic owner.

## 17. Output schema

The future candidate uses schema ID
`compactstar.phase6a1.controlled-bnv-candidate.v1`. Every trajectory row contains
the following typed fields; none is inferred by column position:

| Group | Fields and units |
|---|---|
| epoch/background | `t_s`, `B_count`, `Bdot_count_s`, `DeltaB_over_B0` |
| raw source | `S_n_count_s`, `S_e_count_s`, `S_mu_count_s`, source/event/domain/revision IDs |
| projection | `sigma_e_count_s`, `sigma_mu_count_s`, `bSigma_residual_count_s`, `lift_residual_count_s`, `t_n/e/mu`, `t_error_n/e/mu` |
| chemical state | `eta_e_MeV`, `eta_mu_MeV`, `xi_e`, `xi_mu`, `R_e_count_s`, `R_mu_count_s`, `eta_dot_from_sigma_e/mu_MeV_s`, `eta_dot_from_beta_e/mu_MeV_s` |
| chemical energy | `Echem_MeV`, `Echem_dot_reaction_MeV_s=-eta^T R`, `Echem_dot_source_MeV_s=-eta^T sigma`, `Echem_dot_total_MeV_s` |
| actual potential | `mu_B_inf_MeV`, `mu_n_actual_inf_MeV`, actual-potential provenance and `g_actual` reconstruction residual |
| direct/product | `partition_id`, terminal fate branch IDs/weights, `Eesc_fluid/Star/X_MeV` summaries, `P_dir_eq_erg_s`, `P_dir_actual_erg_s`, R18 residual, R-a/R-b/R-c residuals, `L_out_fluid_inf_erg_s`, `L_esc_star_inf_erg_s`, `J_X_inf_erg_s`, finite-T weighting class and omitted-floor `P_erg_s` |
| beta/thermal | `LH_erg_s`, `DeltaLnu_erg_s`, `DeltaPbeta_erg_s`, `Lnu_eq_erg_s`, `Lnu_full_erg_s`, `Lgamma_erg_s`, `Lother_erg_s`, `Pnet_erg_s` |
| temperature | `Tinf_K`, `Tsurface_inf_K`, `DeltaTinf_K`, `DeltaTsurface_inf_K`, `DeltaLgamma_erg_s`, `DeltaU_th_erg` |
| regime/QSS | per-channel `tau_relax_s`, relevant evolution time, ratio, and enum `LINEAR_QSS/FREEZE_OUT/DRIVE_TRANSIENT/LARGE_XI_QSS/NOT_CLASSIFIED` |
| frozen validity | `DeltaN_n/e/mu_over_N_i`, every named drift bound/threshold/utilization, `max_frozen_utilization`, `valid_through_sample` |
| sign statement | exact observable enum, point or `[t0,t1]`, central value, lower/upper error, classification `HEATING/COOLING/NEAR_ZERO/SIGN_UNRESOLVED` |

Run-level metadata records the entry/implementation commit, dirty-state refusal,
compiler/build identity, all eleven baseline hashes, 33-path manifest hash,
Phase-5 authority identities, source/tangent/partition/fate/sensitivity/run-card
hashes, `RunPurpose`, exact spin-owner identity, transferred qualification literals,
state layout, solver tolerances, initial state, output grid, physical-rate
flag `false`, physical-model flag `false`, candidate/governed status and every gate's
raw return code. The primitive endpoint energies and terminal luminosity columns
must independently reconstruct R20. `Lnu_full` must equal `Lnu_eq+DeltaLnu` within 16 ulps and remains a
diagnostic; `Pnet` subtracts the two components once, not `Lnu_full` in addition.

## 18. Candidate artifact and comparison policy

No governed BNV baseline exists and none is created here. After all pretrajectory
gates and an accepted exact run card, the implementation task may create only:

```text
docs/validation/phase6a1_controlled_bnv_candidate.json
```

Its classification is `CANDIDATE / NOT GOVERNED / NOT PHYSICAL BNV MODEL`. JSON is
canonicalized with sorted keys and binary64 round-trip `max_digits10` values. The
artifact contains the schema ID, complete metadata and named row objects. It embeds
SHA-256 values for source tree/commit, run card, inputs, sensitivity certificate,
all upstream governed artifacts and its producer version; it cannot embed or compare
against itself.

Same-context repeat production must be byte-identical. Scientific comparisons are
field-aware: exact strings/enums/IDs; exact zero/identity fields where specified;
BA1-BA10b analytic budgets; BA11 R20 residual `2e-4`; BA12 component-scaled ODE
norm one; BA13 utilization one; and BA16/BA17 diagnostic criteria. Printed digits are serialization precision, not
physical precision. Promotion would require a separate independent review, owner
ratification and governed integration task with a fresh non-self-comparison producer.

## 19. Answers to the twenty review questions

1. **Can the experiment avoid changing the Phase-5D state vector?** Yes. It reuses
   Thermal(1)+Chem(2); `B` is a prescribed/history diagnostic.
2. **Prescribed or evolved `B`?** Prescribed atomically with `Bdot,S_y`; not an ODE
   state for the first frozen campaign.
3. **Exact `t` owner?** `Analysis::EquilibriumBaryonTangent`, derived only from the
   current Phase-5B `EquilibriumSequenceNumberDerivative` on the same whole-star
   `Omega=0` sequence state.
4. **Exact `sigma` owner?** `Physics::BNV::MovingReferenceSource`, once only.
5. **Exact direct-energy owner?** `BnvDirectEnergyLedger`, separate from chemical
   projection but consuming current typed `eta`, frozen `mu_B^infinity`, governed
   `t`, event stoichiometry and fate to construct `g_actual`; it owns the generic
   event ledger and direct-power GR/unit boundary.
6. **Exact product-fate owner?** Immutable `ProductFateLedger`, terminal branches
   keyed once to source events.
7. **Does existing BNV production satisfy ADR-0015?** No.
8. **Which BNV code must not be reused?** `BNVState`, old BNV header/empty drivers,
   all legacy Microphysics BNV analysis/sequence/channel classes, and MixedStar
   BNV analysis; especially `BNV_B_Chi_Photon` heating.
9. **Where is MeV-to-erg performed?** At two semantically distinct boundaries:
   beta chemical power remains in `RotochemicalThermalPower::From`, while direct
   BNV event power crosses once in `BnvDirectEnergyLedger` after global MeV/s
   integration. The same quantity may never cross twice.
10. **How is P0 represented at finite T?** Exact cold null plus explicit omitted
    finite-T flag, weighting class and propagated floor; never “physical zero heat.”
11. **How is P1 weighted?** Normalized uniform occupied relativistic Track-R
    free-gas sea `3p^2/p_F^3`, with rest-mass-inclusive `E_n(p)` and the exact R10
    average. `2E_F,kin/5` is a separate nonrelativistic-limit check only.
12. **How is P2 bounded?** It is full retention, `Eesc=0`, and a maximum only under
    the declared nonnegative escaping-energy/cold fixture convention.
13. **How is the control constructed?** BA10a separately protects the governed
    spin-on `ControlledTrajectory`. BA10b constructs target and no-BNV control with
    the same identity-bearing `StaticZeroSpinHistory`, `AnalyticControl` base and
    explicitly transferred process/metric/normalization qualification gates,
    replacing only the BNV history/direct bundle by exact zero.
14. **What does sliding null prove?** A physical pure movement along the equilibrium
    sequence has `S=tBdot`, hence no baryon-neutral departure and no chemical drive.
    Raw-G/k routes fail this exact physical negative control.
15. **How is raw-G/k killed?** No production dependency/API accepts them; BA5 injects
    both mutations and requires nonzero wrong results to be detected.
16. **How is double counting killed?** One event ID, one terminal fate, one complete
    direct ledger, one thermal RHS owner, separate Echem diagnostics, and M8-M16.
17. **What exact depletion condition invalidates the run?** The first of
    `|DeltaB|/B0>1e-6`, any certified utilization `>1`, any identity/currentness
    mismatch, or support/domain change; refusal occurs before the offending sample.
18. **Which diagnostics determine thermal sign?** Matched-control
    `DeltaTinf`, `DeltaTsurface_inf`, `DeltaLgamma`, `DeltaU_th` or integrated
    `DeltaP`, with an exact time/window and complete error band—not a single RHS term.
19. **What must hold before any trajectory?** BA1-BA10b and static BA13 certificate
    pass; all protected hashes match; exact mathematical drives/run card are
    committed and owner-accepted; no physical-rate/model flag; fresh output root;
    control/target provenance identical except BNV bundle.
20. **What must hold before a candidate is acceptable?** BA1-BA17 all applicable gates pass, every
    mutation fires, target/control/refinement complete within budgets, output schema
    and hashes are complete, frozen ceiling never fails, repeat bytes match, all
    Phase-5 regressions remain exact, and status remains candidate/not governed.

## 20. Unresolved questions, blockers and exact next action

The two material review findings are closed in this candidate: direct energy now
uses actual state-dependent potentials and the spin-off campaign explicitly
transfers the Phase-5D qualifications omitted by `AnalyticControl`. No new ADR is
required. No ownership ambiguity remains for bounded re-review. The following are
deliberately unresolved but nonblocking:

- exact mathematical drive values and durations, which must be selected by section
  13 and owner-accepted before trajectory generation;
- the measured 21-star error-aware sensitivity certificate and exact target-B
  brackets; failure blocks trajectories rather than authorizing a looser budget;
- physical event weighting, products, rates and finite-temperature direct terms;
- sliding-background ownership and executable `Zdot Z^-1 eta` validation, explicitly
  deferred beyond frozen Phase-6A-1;
- realistic A18, superfluidity, Regime-II/MixedStar thermal evolution and any
  physical BNV specialization.

There is no present blocker to bounded re-review. ADR-0015 is sufficient;
ADR-0016 is not needed. This document does not authorize implementation until the
bounded re-review passes and the owner explicitly accepts the plan.

**Exact recommended next action:** Run a BOUNDED fresh-context Claude Opus 5 XHIGH
re-review of ONLY the delta from
`f9b9029600bab47c61a527672da45516bbe218a6` to
`PHASE6A1_CORRECTION_SHA`, plus the immediately surrounding paragraphs governing
actual versus equilibrium chemical potentials; R18 direct-power identity;
R-a/R-b/R-c; R20 finite-interval closure; zero-spin `RunPurpose` and transferred
qualification gates; BA10a/BA10b; depletion certificate; P1 convention; and the
QSS/B1 additions. The reviewer should not redo the entire ADR-0015 derivation.
Required pass condition: **0 BLOCKING / 0 MATERIAL** before return to the owner for
explicit acceptance of the Phase-6A-1 implementation plan. Do not begin that
re-review automatically.
