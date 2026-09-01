# Passive-cooling regression baseline — Phase 2B-1

> **STATUS: PASSIVE COOLING BASELINE BLOCKED.**
> The authenticated live thermal configuration **cannot execute**. `Config::stepper` defaults to
> `MSBDF`, an implicit BDF method, while `GSLIntegrator::Integrate` hard-codes
> `sys.jacobian = nullptr`. GSL calls a null function pointer and the process dies with `SIGSEGV`.
> No golden values were frozen, and **`GOVERNANCE.md` §3.1 condition 7 is NOT yet satisfied.**

| Field | Value |
|---|---|
| **Purpose** | Freeze the coherent thermal equation `C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞` as a regression reference |
| **Governance role** | `GOVERNANCE.md` §3.1 **condition 7** — the baseline owed immediately after the Phase-2A-3 correction |
| **Source commit** | `56aee7c5a132ffb4d922cac160fda917363ef8e7` |
| **Build configuration** | `Debug` |
| **Toolchain** | AppleClang 17.0.0.17000604, CMake 4.2.1, macOS 15.6.1, arm64 |
| **Harness** | `tests/thermal/passive_cooling_regression.cpp` — built, **not registered as a CTest** |

## 1. §3.1 condition-7 chain

1. **Authorizing ADR** — ADR-0002.
2. **Correction made (Phase 2A-3)** — `−L_γ/1e40` replaced by `−L_γ/C_⋆(T∞)`.
3. **Independent evidence already completed** — ADR-0002 §V1, **V1 VERIFIED**
   (`docs/validation/HEAT_CAPACITY_V1.md`).
4. **Condition being discharged here** — the regression baseline.
   **NOT DISCHARGED.** It remains outstanding, and no unrelated modernization item may begin
   before it is closed.

## 2. The blocker

### Symptom

Running the harness in the authenticated live configuration terminates with `SIGSEGV` (exit 139)
during the first integration segment. Under `lldb`:

```
stop reason = EXC_BAD_ACCESS (code=1, address=0x0)
frame #0: 0x0000000000000000
```

Execution jumps to address `0x0` — a call through a **null function pointer**.

### Cause

| Fact | Evidence |
|---|---|
| The default stepper is an implicit BDF method | `EvolutionConfig.hpp:162` — `StepperType stepper = StepperType::MSBDF;` |
| It maps to GSL's `msbdf` | `GSLIntegrator.cpp:59` — `return gsl_odeiv2_step_msbdf;` |
| No Jacobian is ever supplied | `GSLIntegrator.cpp:334` — `sys.jacobian = nullptr; // analytic Jacobian not implemented` |
| The live program never overrides the stepper | `main/Test/spin_therm_evol_2_main.cpp` contains no `cfg.stepper` assignment |

`gsl_odeiv2_step_msbdf` is implicit and **requires** a Jacobian; GSL invokes `sys.jacobian`
unconditionally. With it null, the call dispatches to `0x0`.

### Scope — this is not confined to the test harness

**Any run using the default stepper crashes**, including the governed end-to-end program
`spin_therm_evol_2_main.cpp`, which §4 of the baseline brief names as the provenance source.

Confirmed by a controlled probe: setting `cfg.stepper = RKF45` (explicit, needs no Jacobian) in
the harness lets the full `100 yr → 1 Myr` integration complete through all nine checkpoints. The
probe was removed; the committed harness uses the authenticated default.

### Pre-existing, not introduced by Phase 2A-3

`git log -S` shows the `MSBDF` default and the `sys.jacobian = nullptr` line have coexisted since
the integrator was added (`ec136fd`, `8586741`, `d07ac2a`). Phase 2A-3 touched only
`PhotonCooling`, `NeutrinoCooling_Details` and documentation — no integrator or config file.

### Not repaired here

The baseline brief freezes production (§18, §20, §25.4). Choosing between *(a)* changing the
default stepper, *(b)* supplying a Jacobian, and *(c)* using an explicit stepper for this workload
is a **numerical-method** decision under `GOVERNANCE.md` §2 — it changes results and needs
convergence evidence and its own classification. It is not an engineering tidy-up, and it must not
be smuggled into a baseline task.

## 3. Second defect found — thermal-only harness impossible

Brief §5 asks for a thermal-only harness (Configuration B) if spin can be shown not to feed the
thermal RHS. **Configuration B cannot run.**

`EvolutionSystem::operator()` (`EvolutionSystem.cpp:103-112`) unconditionally calls
`m_state.GetSpin()` inside a block whose only consumer is commented-out logging:

```cpp
{
    const auto &th = m_state.GetThermal();
    const auto &sp = m_state.GetSpin();
    // Z_LOG_INFO("DEBUG t=" ... );
}
```

A `StateVector` without a registered `Spin` block therefore throws
`StateVector::Get: requested tag 'Spin' is not registered.` The evolution system **cannot run a
purely thermal state at all**, purely as a vestige of removed debug code.

§5 forbids modifying production to force equivalence, so this was **not** repaired. The harness
uses **Configuration A** — the authenticated live spin + thermal wiring.

### Decoupling is instead provable by the energy identity

A stronger check than trajectory comparison is built into the harness: at every checkpoint it
verifies

```
d ln T∞/dt  ==  −(L_ν,∞ + L_γ,∞) / (C_⋆(T∞) · T∞)
```

If the spin block contributed anything to the thermal RHS, this identity would fail. The harness
also asserts that `NeutrinoCooling_Details::C_eff_erg_K` equals `PhotonCooling`'s
`C_star_erg_K` to `1e-12`, so an ADR-0002 Pattern-A violation is caught directly. **Neither check
has been exercised end-to-end**, because the integration cannot run.

## 4. What was completed

- **EOS provenance, both representations** (brief §24) — recorded below.
- **Configuration authenticated from source** at `56aee7c`, not from documentation.
- **Harness written** against production APIs only: `NStar::SolveTOV_Profile`, `StarProfile`,
  `StarContext`, `GeometryCache`, `CompOSE_Thermo`, `PhotonCooling`, `NeutrinoCooling`,
  `EvolutionSystem`, `GSLIntegrator`. No mocked star, no mocked heat capacity, no simplified
  cooling law, no historical result file as input.
- **TOV path verified working** — it reaches `1.6 M☉` before the integrator fails.

### Authenticated configuration (frozen in the harness)

| Parameter | Value | Source |
|---|---|---|
| Structural EoS | `DS(CMF)-1_with_crust` (processed `.eos`) | `spin_therm_evol_2_main.cpp:116` |
| Finite-`T` thermo | `CMF-1_general` = official CompOSE *CMF hadronic (with electrons)* | `:150` |
| Target mass | `1.6 M☉` | `:123` |
| Initial `T∞` | `1.0e9 K` | `:199` |
| Interval | `100 yr → 1e6 yr` | `:324-325` |
| Stepper | **`MSBDF` (default, never overridden)** | `EvolutionConfig.hpp:162` |
| `rtol` / `atol` | `1e-6` / `1e-10` | `:173-174` |
| Cadence | `LogTime`, 150 samples/decade | `:179`, `:182` |
| Envelope | `EnvelopeTbTs`, Potekhin-1997 Iron, `ξ=0`, `ρ_b=1e10` | `:188`, `:243`; `PhotonCooling.hpp:247` |
| Photon | `radiating_fraction=1`, `global_scale=1`, `C_⋆` from `StarContext` | `:244-249` |
| Neutrino | DU on, MU on, PBF off, `global_scale=1` | `:253-256` |
| Spin | braking index 3, `K=1e-15`, `use_moment_of_inertia=false` | `:231-233` |

### EOS provenance — dual representation (brief §24)

CompactStar consumes two forms, and both are recorded:

**Processed** (consumed by TOV via `NStar::SolveTOV_Profile`):

| SHA256 | File |
|---|---|
| `5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd` | `DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos` |

**Raw CompOSE** (the cold table the processed form derives from):

| SHA256 | File |
|---|---|
| `416444999ccac569e2c9b34808888949c36d759f30cce25dab0d42c13e900ce3` | `DS-CMF-1-with-crust/eos.thermo` |
| `d9c8e78c2fcf37fe770fecfc2d3a211d840a28299821a56c77e66f9ff74edef8` | `DS-CMF-1-with-crust/eos.nb` |
| `1a37b9563c40962b203e7bca1aa3b41e8c8b1427953df68095a51dd2cc17ff96` | `DS-CMF-1-with-crust/eos.t` |
| `1a37b9563c40962b203e7bca1aa3b41e8c8b1427953df68095a51dd2cc17ff96` | `DS-CMF-1-with-crust/eos.yq` |
| `4e69b9193e0f075584239d818e1b459791da4d12427531914c86cdd209c898a8` | `DS-CMF-1-with-crust/README` |
| `8b4472405295655cf530572af7edd7448efd0393d0bf9ad86ead3e4c87228c90` | `DS-CMF-1-with-crust/eos.init` |

**Raw CompOSE finite-`T`** (consumed directly by `CompOSE_Thermo`):

| SHA256 | File |
|---|---|
| `a456fb8595208ddf3119350a856fbf2b906c0a0e19bb7c716571748d0aa0724b` | `DNS-CMF-Hadronic-with-electrons/eos.thermo` |
| `3f79dbcc6f8b519696377f89ebc86464bc55cd61d9e2459f6e21e2d9e00f380d` | `DNS-CMF-Hadronic-with-electrons/eos.nb` |
| `2e4c6ec1feb85b16d0ee7036dce183782a9f681577e79c72315171069aa8513d` | `DNS-CMF-Hadronic-with-electrons/eos.t` |
| `d98fcd2f7752039c552c2ef2d04ab485b75db47a61f8ae1740875b54bf9824fd` | `DNS-CMF-Hadronic-with-electrons/eos.yq` |
| `48af68b1b4f6727252ae0051fc35c6445240241d654566973a068d20e1f35222` | `DNS-CMF-Hadronic-with-electrons/README` |
| `412f6739c769df650b3238a6e4b6f0d0f2d7d4a5df43e7d16c80e913bcaddbfb` | `DNS-CMF-Hadronic-with-electrons/eos.init` |

Supplied via `COMPACTSTAR_EOS_DATA_ROOT`. **No EOS data is committed** and no absolute path appears
in committed source. The EOS converter itself was not validated or redesigned (brief §24).

## 5. What could NOT be completed

Every remaining brief item depends on a runnable integration:

| Item | Status |
|---|---|
| Golden checkpoint values | ⛔ not produced |
| Energy-equation identity at each checkpoint | ⛔ not exercised |
| Neutrino/photon dominance ratios | ⛔ not measured |
| Repeat-run determinism | ⛔ not measured |
| Tighter-integrator comparison | ⛔ not measured |
| Cadence-sensitivity check | ⛔ not measured |
| Final tolerance policy | ⛔ not fixed — provisional values in source, unaccepted |
| Regression-detection perturbation proof | ⛔ not demonstrated |
| Runtime measurement | ⛔ n/a |
| CTest registration | ⛔ deliberately withheld |

**No tolerance was accepted and no golden file was written.** The `kTolState` / `kTolLumin`
constants in the harness are placeholders that have not been justified by any measurement, and
must not be treated as an approved policy.

## 6. Standing statements

> **This is a regression baseline, not a validation of the placeholder neutrino emissivity
> normalizations.** `Q0_DU = 1.0e27` and `Q0_MU = 1.0e21` are self-labelled placeholders in
> `NeutrinoCooling_Details.cpp:100-102`. A scientifically justified change to them is expected to
> move the baseline and must be separately reviewed and re-baselined.

> **All passive-cooling outputs predating Phase 2A-3 remain superseded as validation references.**
> They are retained as historical artifacts, not deleted, and must not be used as expected values.

## 7. Required next step

Repair the integrator configuration as its own **numerical-method** change under
`GOVERNANCE.md` §2, then freeze the baseline in the same increment. The decision to make is
between supplying a Jacobian for `MSBDF`, changing the default stepper, and selecting an explicit
stepper for this workload — with convergence evidence for whichever is chosen. `GOVERNANCE.md`
§3.1 condition 7 stays open until that lands.
