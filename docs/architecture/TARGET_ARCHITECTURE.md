# CompactStar — Target Architecture

> **STATUS: NON-NORMATIVE** (`GOVERNANCE.md` authority rank 8). This is intended design.
> It describes where the project is going, **not what it does today**.
>
> For present behavior see [`CURRENT_ARCHITECTURE.md`](CURRENT_ARCHITECTURE.md), which wins on
> any disagreement about the present.
>
> **Provenance.** Distilled from the untracked root `ARCHITECTURE.md` (written against `d91c31b`)
> and `notes.txt` (the owner's original design sketch). The useful intended design is preserved
> here; every claim is re-labeled against authenticated status.

## Labels

**IMPLEMENTED** · **PARTIAL** · **SCAFFOLDED** (exists, not reachable or not compiled) ·
**ABSENT** (planned, no code) · **PROPOSED** (design idea only)

---

## 1. Layered design

| Layer | Intent | Status |
|---|---|---|
| Data primitives (Zaki) | DataColumn / DataSet / Constants | **IMPLEMENTED** (vendored binary; see roadmap Phase 1) |
| EOS | Model hierarchy, CompOSE ingestion, RMF composition solvers | **IMPLEMENTED** |
| Core structure | TOV, profiles, rotation, sequences | **PARTIAL** — see §2 |
| Physics state | Spin / Thermal / Chem / BNV as ODE states | **PARTIAL** — Chem and BNV unexercised |
| Evolution | Layout, packing, caching, RHS assembly, integration | **IMPLEMENTED** |
| Drivers | Pure functions: read state + context, accumulate RHS | **PARTIAL** |
| Microphysics | Reusable rates: Urca, PBF, BNV channels | **ABSENT** for weak rates |
| Observers / diagnostics | Time series, diagnostics, checkpoints | **PARTIAL** |
| Extensions | LightDM, dark-core analysis | **IMPLEMENTED** |

---

## 2. Core structure

| Element | Status |
|---|---|
| `StarProfile` as the single canonical star representation | **IMPLEMENTED** for `NStar` |
| Version-gated invalidation via `m_version` + `EditScope` | **IMPLEMENTED** |
| `MixedStar` migrated to `StarProfile` | **ABSENT** — still `ds_vis`/`ds_dar` with private int indices |
| MixedStar master-grid totals (r, M, p, ε on one grid) | **IMPLEMENTED** by `3639d71` — requires re-audit |
| Single canonical TOV integration path | **ABSENT** — two live paths |
| Hartle O(Ω): ω̄, J, I | **PARTIAL** — live, untested, normalization unresolved (INV-07) |
| Hartle O(Ω²): m₀, p₀, ξ₀, δM | **SCAFFOLDED · CANDIDATE** — unreachable; dropped j² factor; incomplete δM |
| `HartleResult` reachable from `NStar` | **ABSENT** — `rot_solver` is private |
| Profile-backed interpolation at the true RHS radius | **IMPLEMENTED** by `3639d71` — requires re-audit |

---

## 3. Evolution framework

| Element | Status |
|---|---|
| Tag-keyed `StateVector`, flat `y[]` via `StateLayout`/`StatePacking` | **IMPLEMENTED** |
| Additive `RHSAccumulator` — drivers never talk to each other | **IMPLEMENTED** |
| `DriverContext` as the only channel for static context | **IMPLEMENTED** |
| Single uniform cache-invalidation rule | **ABSENT** — five rules coexist (INV-12) |
| `GeometryCache` version-gated | **ABSENT** — no gate at all |
| Driver dependency DAG resolving `DependsOn()`/`Updates()` | **ABSENT** — planned as `Evolution/Graph` in `notes.txt`; `Coupling.hpp` is 0 bytes |
| Diagnostics share driver physics via `*_Details` | **IMPLEMENTED** — genuinely prevents drift |
| Populated `UnitContract` per producer | **SCAFFOLDED** — always empty |
| `CheckpointObserver` | **ABSENT** — 0-byte file |

---

## 4. Rotochemical heating — target pipeline

The scientific objective: reproduce standard non-superfluid rotochemical heating
(Fernández & Reisenegger 2005) before extending to BNV.

```
Ω ──► spin-down driver ──► Ω̇                                   IMPLEMENTED
      ↓
   structural coefficients  A_i, B_i → Z_i = A_i − B_i(A_B/B_B)  SCAFFOLDED
      ↓
   η_npe , η_npμ                                                 ABSENT
      ↓
   out-of-equilibrium weak rates  ΔΓ(η,T),  F/H(ξ = η/k_BT)      ABSENT
      ├─► chemical restoration                                   ABSENT
      ├─► modified neutrino luminosity                           ABSENT
      └─► heating                                                ABSENT
      ↓
   T∞                                                            IMPLEMENTED (passive only)
```

| Component | Status |
|---|---|
| Enclosed species numbers `N_i` | **SCAFFOLDED** — not compiled; must construct `n_i = Y_i n_B` per ADR-0001 (currently does not) |
| `A_i = (∂N_i/∂Ω²)|_{ε_c}` | **PARTIAL** — never divided by Ω² |
| `B_i = (∂N_i/∂ε_c)|_Ω` | **PARTIAL** — perturbed grids × unperturbed volume weight |
| `Z_i` reduction | **IMPLEMENTED in form** — correct, but not compiled |
| Spin-down forcing `2ΩΩ̇` | **SCAFFOLDED** |
| η definitions, ordering, redshift frame | **ABSENT** (INV-11) |
| Out-of-equilibrium weak rates | **ABSENT** — no function anywhere takes η |
| `WeakRestoration` | **ABSENT** — 0-byte file |
| `HeatingFromChem` | **SCAFFOLDED** — header only |
| Single-source Γ preventing double counting | **PROPOSED** |
| F&R 2005 regression | **ABSENT** |

**The out-of-equilibrium weak rates are the terminal blocker.** Everything downstream of η
depends on them, and none of it exists.

---

## 5. Design principles (intended, and worth keeping)

1. **Generic machinery over bespoke conditions.** No reaction, campaign, or species name is
   hardcoded into generic numerical infrastructure.
2. **One authoritative owner per quantity.** A profile belongs to the structural representation;
   a solver returns a solver result; a cache owns derived coefficients; a driver contributes an
   RHS; diagnostics never independently re-derive physics.
3. **Pure-function drivers.** `AccumulateRHS` is const; contributions are additive.
4. **Explicit context, no hidden global state.**
5. **Fail closed on ambiguity** (`GOVERNANCE.md` §3).

Principle 2 is currently violated by proper volume (INV-04) and the duplicated TOV paths, and —
**in the source only** — by heat capacity. Heat-capacity *ownership* is no longer ambiguous:
**ADR-0002 (ACCEPTED 2026-08-31)** names one physical owner, `C_⋆(T∞)`, for the thermal degree of
freedom (INV-15). What remains is a live **source nonconformance** — `PhotonCooling` still divides
by a driver-local constant — scheduled as roadmap Phase 2A. A governed convention that the code
does not yet meet is a different thing from an ungoverned ambiguity, and only the second is a
fail-closed condition.

ADR-0002 deliberately leaves open *where* the division by `C_⋆(T∞)` happens: each driver dividing
by a shared `C_⋆` (preserving principle 3 above) or a centralized thermal-balance owner consuming
power contributions. That choice is a target-architecture question and is recorded as an
anticipated ADR, not decided here.

---

## 6. Planned but never built

From `notes.txt`, still absent: `Pulsar` as a first-class component, `Evolution/Graph` (the
dependency DAG), `EOSHooks`, `Rates/PairBreaking`, `Rates/Emissivities`, `LogObserver`,
`Utilities/`.

These are recorded so the gap between intent and reality stays visible. Their absence is not a
defect — only their absence *while being described as present* would be.
