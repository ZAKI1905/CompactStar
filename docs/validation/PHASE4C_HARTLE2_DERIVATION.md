# Phase 4C-G — Hartle O(Ω²) monopole (l = 0) system: primary-source derivation and governance record

> **FORMAL STATUS: `PHASE-4C O(OMEGA^2) MONOPOLE ADR PROPOSED — OWNER ADJUDICATION REQUIRED`**
>
> This is a **primary-source derivation / governance** record (increment 4C-G). **No production
> source, test, baseline or CMake file was changed.** The O(Ω²) candidate (`ODE_Hartle2_N_Fast`,
> `SolveHartle2_N`, the two MixedStar stubs, `HartleResult`) is byte-identical to `bb073c8` and to
> `df859b5`. Nothing here validates any production O(Ω²) number; no candidate output was
> executed, baselined or used as an oracle. `ADR-0007` is drafted **PROPOSED** — it authorizes
> nothing until the owner adjudicates it.

| Field | Value |
|---|---|
| **Starting HEAD** | `bb073c8ed0c7ce15f3a8c960e9f76173bde51a39` (`test: validate physical Hartle first-order response`) — `origin/physics/rotation-correctness` equal, tree clean, **4 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Branch / worktree** | `physics/rotation-correctness` at `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotation` |
| **Change class** | documentation / governance (`GOVERNANCE.md` §2); every claim about the repository carries a `path:line` citation at `bb073c8`; every claim about the physics carries a page/equation citation to the primary source |
| **Authority read, in order** | `GOVERNANCE.md` (§1–§5, §3.1 verbatim); ADR-0001…ADR-0006 (all ACCEPTED; ADR-0006 Decision and binding clarification); `docs/SCIENTIFIC_INVARIANTS.md` INV-05, INV-06, INV-07, INV-08, INV-09, INV-14; `docs/MODERNIZATION_ROADMAP.md` Phase 4 (incl. the ratified exit criterion); `docs/architecture/CURRENT_ARCHITECTURE.md`; `docs/validation/PHASE4_ROTATION_ENTRY.md`; `PHASE4A_FIRST_ORDER_NORMALIZATION.md`; `PHASE4B_FIRST_ORDER_PHYSICS.md`; `HARTLE_MOMENT_INERTIA.md`; `TOV_REFERENCE.md` §5; `AGENTS.md` |
| **Primary literature** | **J. B. Hartle, "Slowly Rotating Relativistic Stars. I. Equations of Structure", ApJ 150, 1005–1029 (1967)** — read in full from the journal scan (§2) |
| **Secondary literature** | J. B. Hartle & K. S. Thorne, "Slowly Rotating Relativistic Stars. II. Models for Neutron Stars and Supermassive Stars", ApJ 153, 807–834 (1968) — §II read from the scan (conventions, `r` vs `R`, boundary conditions, `δM`, binding energy, the printed HW / Vγ tables); O. Benhar, V. Ferrari, L. Gualtieri, S. Marassi (2005), arXiv:gr-qc/0504068 (a modern transcription in different variable names, used only to confirm two misprints) |
| **Inherited, binding** | first order is physically normalized (ADR-0006 ACCEPTED, Phase 4A) and independently verified (Phase 4B): the canonical first-order response is `I`, `s(r) = ω̄/Ω`, `s'(r) = ω̄'/Ω`, seed-free. **No arbitrary first-order seed enters anything proposed here.** |

**Where the task's twenty required items live in this record**

| Required item | Section |
|---|---|
| 1 literature authorities | §2 |
| 2 exact metric / matter conventions | §3 |
| 3 variable definitions | §4 |
| 4 unit table | §5 |
| 5 l = 0 equations term by term | §6 |
| 6 normalization by Ω² | §7 |
| 7 centre series | §9 |
| 8 fixed-ε_c condition | §8 |
| 9 EOS derivative authority | §10 |
| 10 surface-convention audit | §11 |
| 11 exterior matching | §12 |
| 12 δM | §12 |
| 13 ξ₀ | §13 |
| 14 l = 0 / l = 2 scope | §14 |
| 15 Phase-5 structural response interface | §15 |
| 16 candidate discrepancy table | §16 (+ §17 hypothesis classification) |
| 17 §3.1 applicability | §18 |
| 18 implementation dependency map | §20 |
| 19 validation plan | §21 (+ §22 detectors) |
| 20 owner questions | §23 |

---

## 1. Authentication — state and inherited baseline

```
worktree  /Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotation
branch    physics/rotation-correctness      HEAD = bb073c8ed0c7ce15f3a8c960e9f76173bde51a39
upstream  origin/physics/rotation-correctness = bb073c8   (0 ahead / 0 behind upstream)
master    df859b5a73c4cac0c115f240744d89ce9f830b8d       (4 ahead / 0 behind)
status    clean
```

Pre-task suite, run on this exact tree before any file was written (`Debug`, AppleClang 17,
CMake 4.2, GSL 2.7.1, macOS arm64 — `docs/build/MACOS_BUILD.md`):

| Configuration | Registered | Result | Wall time |
|---|---|---|---|
| Full — `-DCOMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose` | 23 | **23/23 PASS** | 211.40 s |
| Self-contained — no data root (`build-selfcontained/`) | 12 | **12/12 PASS** | 13.78 s |

Seven durable artifacts, `shasum -a 256`, identical to `PHASE3_CLOSEOUT.md` §1 and to the 4A-0,
4A and 4B records (no regeneration):

| Artifact | SHA-256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |
| `tov_path_equivalence_dscmf1.tsv` | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |

Next ADR number authenticated: `docs/adr/` holds `ADR-0001` … `ADR-0006` and the template; the
index (`docs/adr/README.md`) lists ADR-0006 as the last allocated number. **ADR-0007** is next.

## 2. Literature authorities — primary access authenticated

**How the primary source was read.** The NASA ADS journal scan of ApJ 150, 1005–1029 (1967) and
of ApJ 153, 807–834 (1968) were obtained as PDFs (the scan's OCR text layer is unusable — every
text-extraction attempt returned garbage), rendered page by page, and read **visually**. Every
equation cited below was read from the page image; the four equations that matter most, and
the three places where a sign or a factor was ambiguous at page resolution, were re-rendered at
220 dpi and re-read. Nothing below is reconstructed from memory or from the candidate's
comments. The 4A-0 transcription (`PHASE4_ROTATION_ENTRY.md` §10.1) is treated as a hypothesis
and classified in §17.

**2.1 Primary — Hartle (1967), ApJ 150, 1005 ("H67"). Citations used in this record.**

| Page | Eq. | Content |
|---|---|---|
| 1007 | (1)–(2) | slow-rotation condition `Ω² ≪ (c/R)²(GM/Rc²)`, i.e. `RΩ ≪ c` |
| 1008–1009 | (6)–(9) | the constant-density label coordinate `R`: `Θ = θ`, `ρ[r(R,Θ),Θ] = ρ(R) = ρ⁽⁰⁾(R)`; `r = R + ξ(R,Θ) + O(Ω⁴)`; **"if the rotating configuration is chosen to have the same central density as the non-rotating configuration … ξ vanishes at R = 0"** (p. 1009, below eq. 8) |
| 1009 | (11)–(13) | `ξ = Σ ξ_l P_l`; only `l = 0, 2` survive |
| 1010 | (14)–(19) | Newtonian `l = 0` system; `δM = 4π∫(−ξ₀ dρ/dR)R²dR = M⁽²⁾(a)` (18); centre `p* → (1/3)Ω²R²` (19) |
| 1011 | (20)–(21) | Newtonian `ξ₀(a) = (a²/GM)p*(a)` |
| 1012 | (25)–(29) | background metric `ds² = −e^{ν}dt² + e^{λ}dr² + r²dΩ²`, `e^{−λ} = 1 − 2M(r)/r`; TOV (28); `dM/dr = 4πr²E`, `dν/dr = −(2/(E+P))dP/dr` (29) |
| 1012–1014 | (30)–(41) | stationary axisymmetric metric (30); `L = ω + O(Ω³)` (31); field equation `R_φ^t = 8πT_φ^t` (32) with the identity (33)–(35); `ω̄ = Ω − ω` (38); `j(r) = exp[−(ν+λ)/2]` (40) |
| 1015 | (46)–(49) | `(1/r⁴)(r⁴jω̄')' + (4/r)j'ω̄ = 0` (46); exterior `ω̄ = Ω − 2J/r³` (47); `J = IΩ` (48); `GI/c² = −(2/3)∫₀^a r³ j' (ω̄/Ω) dr` (49) |
| 1016 | (61)–(62) | `j > 0`; `−dj/dr = 4πr e^{−ν}(E+P)/j > 0` |
| 1017 | (66)–(70) | second-order metric: `ds² = −e^{ν}(1+2h)dt² + e^{λ}[1 + 2m/(r−2M)]dr² + r²(1+2k)[dθ² + sin²θ(dφ − ωdt)²] + O(Ω³)`; `h = h₀ + h₂P₂`, `m = m₀ + m₂P₂`, `k = k₀ + k₂P₂`; gauge `k₀ = 0` |
| 1018 | (71)–(81) | `ℰ[r(R,Θ),Θ] = E(R)`, `r = R + ξ + O(Ω⁴)`; `ΔG = δG + ξ ∂_R G⁽⁰⁾` (73)–(78); **(2d-order term in ℰ) = −ξ dE/dR, (2d-order term in 𝒫) = −ξ dP/dR** (79)–(80); `ξ = ξ₀ + ξ₂P₂` (81) |
| 1019 | (83)–(91) | integral of hydrostatic equilibrium `μ_c = ((ℰ+𝒫)/u^t)exp(−∫dℰ/(ℰ+𝒫))`; `μ_c = μ[1 + γ + O(Ω⁴)]`; `γ = p* + h − (1/2)e^{−ν}ω̄²R²sin²Θ` (86); **`p*(R,Θ) = (1/2)ν_R ξ = −ξ (1/(E+P)) dP/dR`** (87); `p₀* = −ξ₀(1/(E+P))dP/dR` (88); `l = 0`: `γ = p₀* + h₀ − (1/3)e^{−ν}R²ω̄²` (90); `l = 2`: `0 = p₂* + h₂ + (1/3)e^{−ν}R²ω̄²` (91) |
| 1020 | (92)–(98) | `g_tt → −(1−2M/r)` defines `M`; `(ΔG_t^t)_{l=0}`, `(ΔG_R^R)_{l=0}` (93)–(94); `(ΔT_t^t)_{l=0} = (2/3)(j²)_R Ωω̄R` (95); `(ΔT_R^R)_{l=0} = 0` (96); **`dm₀/dR = 4πR²(dE/dP)(E+P)p₀* + (1/12)j²R⁴(ω̄_R)² − (1/3)R³(j²)_R ω̄²`** (97); **`dh₀/dR − (m₀R²/(R−2M)²)(8πP + 1/R²) = 4π(E+P)R²p₀*/(R−2M) − (1/12)R⁴j²(ω̄_R)²/(R−2M)`** (98) |
| 1021 | (99)–(107) | `p₀* = −ξ₀(1/(E+P))dP/dR` (99); **`−dp₀*/dR + (1/12)(R⁴/(R−2M))j²(ω̄_R)² + (1/3)(d/dR)[R³j²ω̄²/(R−2M)] = 4π((E+P)R²/(R−2M))p₀* + (m₀R²/(R−2M)²)(8πP + 1/R²)`** (100); "the linearity of these equations … means that from the solution for one value of Ω the solution for other values may be obtained by scaling" (p. 1021); Newtonian limit (101)–(104); exterior **`m₀ = δM − J²/r³`** (105), `h₀ = −δM/(r−2M) + J²/(r³(r−2M))` (106), **`δM = m₀(a) + J²/a³`** (107) |
| 1022 | (108)–(113) | **centre: `p₀* → (1/3)(j_c ω̄_c)²R²`, `m₀ → (4π/15)(E_c+P_c)[(dE/dP)_c + 2](j_c ω̄_c)²R⁵`**, "these boundary conditions guarantee that the central density of the rotating and non-rotating configurations are the same" (108); baryon number `A = ∫d³x (−g)^{1/2} u^t 𝔑` (109); `d𝔑/dℰ = 𝔑/(ℰ+𝒫)` (110); `−E_B = M − μA` (111); `−δE_B = δM − μδA` (112); `ε = ℰ − μ𝔑` (113); "**Only l = 0 perturbations will survive in this expression**" (p. 1022, on δA) |
| 1023 | (114)–(117) | `δE_B = −J²/a³ + ∫4πR²B(R)dR` (114); `B(R)` (115); surface `r = a + ξ₀(a) + ξ₂(a)P₂` (116); `ξ₀(a) = p₀*(a)a(a − M)/M` as printed (117) — see 2.3 |
| 1023–1027 | (118)–(154) | the `l = 2` sector (Clairaut generalization); `Q = J²/M + 16AM³/5` (138, corrected in HT68 fn 5) |
| 1028–1029 | (A.1)–(A.15) | second-order Einstein-tensor components; `(2d order G_t^t)_{l=0}` (A.9), `(2d order G_r^r)_{l=0}` (A.11) |

**2.2 Secondary — Hartle & Thorne (1968), ApJ 153, 807 ("HT68"), §II.**

| Page | Eq. | Content |
|---|---|---|
| 808 | (1)–(3d) | one-parameter EOS `E = E(P)`, `N = N(P)`; TOV; `dν/dr = −2(E+P)^{−1}dP/dr` |
| 809 | (3e) | non-rotating baryon number `A = ∫₀^R N(r)[1 − 2M(r)/r]^{−1/2}4πr²dr` — **identical in form to INV-14** |
| 809 | (4)–(8) | second-order metric (same as H67 66); `u^t = e^{−ν/2}[1 + (1/2)r²sin²θ(Ω−ω)²e^{−ν} − h₀ − h₂P₂]` (5); **`P + (E+P)(p₀* + p₂*P₂) ≡ P + ΔP`, `E + (E+P)(dE/dP)(p₀* + p₂*P₂) ≡ E + ΔE`, `N + (E+P)(dN/dP)(p₀* + p₂*P₂) ≡ N + ΔN`** "in a reference frame momentarily moving with the fluid" at `(r,θ)` (7a–c); "`p₀*` and `p₂*` are dimensionless functions of `r`, proportional to `Ω²`" |
| 810 | (9)–(14) | first order; `j = e^{−ν/2}[1 − 2M/r]^{1/2}` (10); `J = (1/6)R⁴(dω̄/dr)_R`, `Ω = ω̄(R) + 2J/R³` (13); **`ω̄_new = ω̄_old(Ω_new/Ω_old)`** (14) — the INV-07 / ADR-0006 rescaling, in the primary author's own words |
| 810 | (15a)–(15b) | the `l = 0` system written in `r` — **byte-for-byte the same terms as H67 (97) and (100)**; footnote 4: "in the limit of slow rotation the equations are the same when written in terms of `r` as when written in terms of `R`"; "**integrated outward, with the boundary conditions that both `m₀` and `p₀*` vanish at the origin. With these boundary conditions, the rotating star will have the same central density as the non-rotating one**" |
| 811 | (15c)–(17b) | `m₀ = δM − J²/r³` outside (15c); `M + δM = M + m₀(R) + J²/R³` (16); `h₀` outside (17a) and **inside `h₀ = −p₀* + (1/3)r²e^{−ν}ω̄² + h_{0c}`**, `h_{0c}` fixed by continuity at the surface (17b) |
| 811 | (18)–(20) | `E_B = μA − M`; `ε = E − μN`; `δE_B = −J²/R³ + ∫4πr²B dr` with `B(r)` (20) — **the `dε/dP` term carries a minus sign** |
| 812 | (24b), (25b) | **`ξ₀ ≡ δr = −p₀*(E+P)/(dP/dr)`**; mean radius `r̄* = r + ξ₀(r)`, `R̄ = R + δR` |
| 813 | (27) | `dN/N = dE/(E+P)` (adiabatic compression, isentropic) |
| 814–819 | Tables 1–6 | printed HW and Vγ EOS tables; non-rotating and slowly rotating models (`δM/M`, `δR/R`, `ω_c/Ω`, `ω_s/Ω`, `R_g/R`, `e_s`, `Q/MR²` at `Ω² = M/R³`, "accurate to about 1 per cent") |

**2.3 Two misprints in the primary source, resolved.** Neither affects the governed system.

| Place | As printed | Correct form | Evidence |
|---|---|---|---|
| H67 (117), p. 1023 | `ξ₀(a) = p₀*(a) a(a − M)/M` | `ξ₀(a) = p₀*(a) a(a − 2M)/M` | follows from H67's own (99) with the exterior `ν_R = 2M/[r(r−2M)]` from (27); HT68 (24b); Benhar et al. (2005) eq. (51) |
| H67 (115), p. 1023 | `+ (dε/dP)(1 − 2M/R)^{−1/2}` inside the `p₀*` bracket | `− (dε/dP)(1 − 2M/R)^{−1/2}` | HT68 (20); the independent derivation of `δA` in §14–§15 reproduces HT68 (20) term by term, sign included |

**2.4 Naming hazard in the secondary literature.** Benhar et al. write the `l = 0` system in a
variable they call `δp₀` which multiplies `(ε+P)` in their `m₀` equation and enters their
`ξ₀(R) = δp₀(R)R(R−2M)/M` — i.e. **their `δp₀` is Hartle's dimensionless `p₀*`, not the Eulerian
pressure perturbation.** HT68 (7a) likewise labels the *Eulerian* change at fixed `(r,θ)` as
"`ΔP`". Both are consistent with H67 once the definitions of §4 are fixed; the repository must
not import either name without the mapping.

## 3. Conventions — metric, matter, coordinates

| Object | Hartle (H67) | CompactStar (INV-03, INV-02) | Mapping |
|---|---|---|---|
| background metric | `ds² = −e^{ν_H}dt² + e^{λ_H}dr² + r²dΩ²` (25) | `g_tt = −e^{2ν}`, `g_rr = e^{2λ}` (`StarProfile.hpp:246-247`) | `ν_H = 2ν`, `λ_H = 2λ` |
| `g_rr` | `e^{−λ_H} = 1 − 2M(r)/r` (26) | `e^{2λ} = 1/(1 − 2m/r)`; `m` in km (`Mass`, `NStar.cpp:105`) | `M_H ≡ m` |
| TOV | `−dP/dr = (E+P)(M + 4πr³P)/(r(r−2M))` (28) | `dp/dr = −(ε+p)ν'`, `ν' = (m + 4πr³p)/(r(r−2m))` (`MetricNuPrime`, km⁻¹, `NStar.cpp:187`) | `dν_H/dr = −2(E+P)^{−1}dP/dr` (29) ⇔ `ν' = −(ε+p)^{−1}dp/dr` |
| `j` | `j = exp[−(ν_H+λ_H)/2]` (40) `= e^{−ν_H/2}(1 − 2M/r)^{1/2}` (HT68 10) | `j = e^{−(ν+λ)} = e^{−ν}(1 − 2m/r)^{1/2}` | identical function |
| `j²`, `dj²/dr` | `−dj/dr = 4πr e^{−ν_H}(E+P)/j` (62) | `j² = e^{−2ν}(1 − 2m/r)`; **`dj²/dr = −8πr(ε+p)e^{−2ν} = −8πr²(ε+p)j²/(r−2m)`**; `j²/(r−2m) = e^{−2ν}/r` | identities used throughout §6–§7 |
| first order | `ω̄ = Ω − ω` (38); eq. (46) | `RotationSolver.hpp:72`; `ODE_N_Fast`; `EQUATION MATCH` (`HARTLE_MOMENT_INERTIA.md` §4) | validated (2B-4B, 4B) |
| second-order `l = 0` metric | `−e^{ν_H}(1 + 2h₀)dt²`, `e^{λ_H}[1 + 2m₀/(r−2M)]dr²`, `k₀ = 0` (66)–(70) | `m₀` [km] is the O(Ω²) perturbation of the enclosed-mass function; `h₀` dimensionless; gauge `k₀ = 0` adopted | — |
| matter labelling | `R`: `ℰ[r(R,Θ),Θ] = E(R)`, `r = R + ξ(R,Θ)` (71)–(72) | integration on the profile grid `r`; coefficient functions evaluated at `r` (HT68 fn 4: same to O(Ω²)) | `R`-vs-`r` difference is O(Ω⁴) in every coefficient |
| units | `G = c = 1` | geometric km (INV-02); `Ω_geom = Ω_phys/c` (ADR-0006, sole owner `AngularVelocity::GeomKmInverse`) | — |

**Regime.** The expansion is in `Ω²` with the requirement H67 (1)–(2); HT68 tabulate at
`Ω² = M/R³` and call it "approximately the critical angular velocity at which rotational
shedding will occur … an upper bound". ADR-0006 P10's diagnostic `Ω/Ω_K`, `Ω_K = (M/R³)^{1/2}`,
is exactly HT68's scale. Phase 4B's disclaimer carries forward unchanged: correctness of the
O(Ω²) coefficients is not accuracy of the slow-rotation truncation — at `Ω/Ω_K ≈ 0.4–0.6`
(716 Hz on the CMF sequence) the O(Ω²) terms are `10–20 %` effects and O(Ω⁴) is not negligible.

## 4. Variable definitions — what `p₀*` is

**4.1 Derivation from the primary text.** H67 (83)–(86): the first integral of hydrostatic
equilibrium for uniform rotation and a one-parameter EOS is `μ_c = ((ℰ+𝒫)/u^t)exp(−∫dℰ/(ℰ+𝒫))`
= const. Expanding to O(Ω²) with `u^t` from (36) gives `μ_c = μ[1 + γ]` with

```
γ = p*(R,Θ) + h(R,Θ) − (1/2) e^{−ν_H(R)} ω̄(R)² R² sin²Θ                        H67 (86)
p*(R,Θ) ≡ (1/2) ν_{H,R} ξ(R,Θ) = −ξ (1/(E+P)) dP/dR                                  H67 (87)
p₀*(R)  = −ξ₀(R) (1/(E+P)) dP/dR ,     p₂*(R) = −ξ₂(R) (1/(E+P)) dP/dR              H67 (88)–(89), (99)
```

Hartle names `p*` the "relativistic pressure perturbation factor". It is **dimensionless**
(`ξ₀ [km] × ν' [km⁻¹]`), and by (87) it is, equivalently, `p₀* = ν' ξ₀` in CompactStar's `ν`
(`(1/2)ν_{H,R} = ν'`).

**4.2 The Eulerian perturbations.** H67 (79)–(80): the second-order term in `𝒫` at fixed
coordinate point `(r,θ)` is `−ξ dP/dR`; with (87) that is `(E+P)p*`. HT68 (7a)–(7c) write the
same thing for pressure, energy density and baryon density. Hence, for the monopole:

```
δp₀ ≡ p(r) − P(r)  = (ε+p) p₀*                 [km⁻²]      Eulerian pressure perturbation
δε₀ ≡ ε(r) − E(r)  = (ε+p) (dε/dp) p₀*         [km⁻²]      Eulerian energy-density perturbation
δn₀ ≡ n(r) − N(r)  = (ε+p) (dn/dp) p₀*         [fm⁻³]      Eulerian number-density perturbation (any barotropic n(p))
```

**4.3 The Lagrangian perturbation.** By construction of `R` (H67 (6), (9), (71)) the pressure,
energy density and number density at the *displaced* point `r = R + ξ₀` equal their
non-rotating values at `R`: the Lagrangian monopole changes vanish identically,

```
Δp₀ ≡ δp₀ + ξ₀ dP/dr = (ε+p)p₀* − p₀*(ε+p) = 0 ,     Δε₀ = 0 ,     Δn₀ = 0 .
```

So `ξ₀` is *defined* by `Δp₀ = 0`; it is not an additional condition that can be imposed at a
surface.

**4.4 The displacement.** `ξ₀(r) = −δp₀/(dp/dr) = p₀*/ν'` [km] (H67 (99), HT68 (24b)): the
outward displacement of the isobar (equivalently the constant-density surface) that sits at `r`
in the non-rotating star. Sign: `dp/dr < 0` and `p₀* > 0` ⇒ `ξ₀ > 0` (the rotating star is
puffed up at fixed `ε_c`).

**4.5 Physical reading.** For a barotrope `dh/h = dp/(ε+p)` with `h = (ε+p)/n` the enthalpy per
baryon, so `p₀* = δp₀/(ε+p) = δh/h|_{Euler}`: **`p₀*` is the Eulerian fractional enthalpy
perturbation**, and (90) `γ = p₀* + h₀ − (1/3)e^{−2ν}r²ω̄² = const` is the statement that the
redshifted injection energy stays uniform through O(Ω²). That constancy is an internal
conservation identity the implementation can be checked against (§21, line G′).

**4.6 The four quantities side by side.**

| Quantity | Definition | Dimension (geometric) | `r → 0` | at the surface |
|---|---|---|---|---|
| `p₀*` | `−ξ₀ (dP/dr)/(ε+p) = ν'ξ₀ = δp₀/(ε+p)` | 1 | `(1/3)j_c²ω̄_c²r²` (§9) | finite, generally `≠ 0` |
| `δp₀` (Eulerian) | `(ε+p)p₀*` | km⁻² | `∝ r²` | `→ (ε_s + p_s)p₀*(R_*)`; `→ 0` only if `ε + p → 0` there |
| `Δp₀` (Lagrangian) | `δp₀ + ξ₀ dp/dr` | km⁻² | `0` | `0` (identically) |
| `ξ₀` | `p₀*/ν'` | km | `∝ r` (§13) | finite, generally `≠ 0` (§13) |

## 5. Unit table

Geometric units, `G = c = 1`, lengths in km (INV-02); `Ω_geom = Ω_phys/c` with
`c = LIGHT_C_KM_S` owned by `AngularVelocity` (ADR-0006 P2). "Normalized" means divided by
`Ω_geom²` (§7).

| Symbol | Meaning | Geometric unit | Normalized (`÷ Ω_geom²`) |
|---|---|---|---|
| `r`, `R_*`, `m`, `M` | radius, surface radius (production last node), enclosed mass, total mass | km | — |
| `ε`, `p` | energy density, pressure (profile `eps(km^-2)`, `p(km^-2)`) | km⁻² | — |
| `n_B`, `n_i = Y_i n_B` | baryon / species number density (ADR-0001) | fm⁻³ (×`1e54` → km⁻³ count density, INV-14) | — |
| `ν`, `λ`, `h₀`, `k` | metric functions | 1 | `ĥ₀`: km² |
| `ν'`, `λ'` | radial derivatives | km⁻¹ | — |
| `j`, `j²` | `e^{−ν−λ}` | 1 | — |
| `dε/dp` | EOS derivative `= 1/c_s²` | 1 | — |
| `ω̄`, `ω`, `Ω` | frame-dragging variables | km⁻¹ | `s = ω̄/Ω`: 1 |
| `ω̄'` | | km⁻² | `s' = ω̄'/Ω`: km⁻¹ |
| `J` | angular momentum | km² | `J/Ω = I`: km³ |
| `I` | moment of inertia | km³ | — |
| `m₀` | O(Ω²) enclosed-mass perturbation | km | **`m̂₀`: km³** |
| `p₀*` | pressure perturbation factor | 1 | **`p̂₀*`: km²** |
| `δp₀` | Eulerian pressure perturbation | km⁻² | **`δp̂₀`: 1** |
| `ξ₀` | isobar displacement | km | **`ξ̂₀`: km³** |
| `δM` | total-mass change at fixed `ε_c` | km | **`δM̂`: km³** |
| `δA`, `δN_i` | baryon / species count change | 1 | **`δÂ`, `δN̂_i`: km²** |
| `Ω_phys` | physical angular velocity | s⁻¹ (`= c·Ω_geom`) | — |

Physical values follow by one multiplication: `Q_phys = Q̂ · Ω_geom²`. Conversion of `δM` to
`M☉` uses `SUN_M_KM` (the same unadjudicated solar-mass authority ADR-0006 defers; it does not
enter any geometric result). Per-`Ω_phys²` coefficients are `Q̂/c²`.

## 6. The `l = 0` equations, term by term

The two Einstein equations Hartle solves are `G_t^t = 8πT_t^t` and `G_R^R = 8πT_R^R` for
`l = 0` (H67 p. 1020), with (A.9), (A.11) for the left-hand sides, (74)–(75) for the `ξ`
shifts and (95)–(96) for the sources. Their combination with the differentiated integral (90)
gives (97) and (100), which HT68 (15a)–(15b) restate in `r`. CompactStar form, with every metric
factor explicit (`j² = e^{−2ν}(1 − 2m/r)`; `j²/(r−2m) = e^{−2ν}/r`; `dj²/dr = −8πr(ε+p)e^{−2ν}`):

```
dm₀/dr  = 4π r² (ε+p) (dε/dp) p₀*  +  (1/12) j² r⁴ (ω̄')²  −  (1/3) r³ (dj²/dr) ω̄²                            (H1) ≡ H67 (97)

dp₀*/dr = − m₀ (1 + 8π r² p)/(r − 2m)²  −  4π (ε+p) r² p₀*/(r − 2m)
          + (1/12) r⁴ j² (ω̄')²/(r − 2m)  +  (1/3) (d/dr)[ r³ j² ω̄²/(r − 2m) ]                                  (H2) ≡ H67 (100)

dh₀/dr  = + m₀ (1 + 8π r² p)/(r − 2m)²  +  4π (ε+p) r² p₀*/(r − 2m)  −  (1/12) r⁴ j² (ω̄')²/(r − 2m)            (H3) ≡ H67 (98)

p₀* + h₀ − (1/3) e^{−2ν} r² ω̄² = γ  (constant)                                                                  (H4) ≡ H67 (90)
```

(H2) is (H3) with `dh₀/dr` eliminated through (H4); (H3) is therefore redundant for the
solution but available as an internal check (§21). `m₀(1+8πr²p)/(r−2m)² ≡ m₀r²(8πp + 1/r²)/(r−2m)²`
as printed.

**6.1 Equation (H1) — `dm₀/dr`, dimension of every term: 1 (dimensionless).**

| Literature term (H67 97) | Symbol definition | CompactStar background needed | Dimension | Candidate counterpart (`RotationSolver.cpp`) | Candidate verdict |
|---|---|---|---|---|---|
| `4πR²(dE/dP)(E+P)p₀*` | EOS derivative × Eulerian `δp₀` | `r`, `ε`, `p` (profile); `dε/dp` (**EOS**, §10) | km²·km⁻²·1·1 = 1 ✓ | `4π r² dEdP p0` (`:1170`) with `p0 ≡ δp` and profile-FD `dEdP` (`:1254-1277`) | form matches for `p0 ≡ δp`; derivative source wrong authority |
| `+(1/12) j² R⁴ (ω̄_R)²` | rotational-energy source | `e^{−2ν}` (`MetricNu`), `1 − 2m/r`, `ω̄'` | km⁴·km⁻⁴ = 1 ✓ | `(1/12) r⁴ (ω̄')²/(1 − 2m/r)` (`:1204`) | **wrong**: reciprocal of the `(1 − 2m/r)` factor, `e^{−2ν}` absent |
| `−(1/3) R³ (j²)_R ω̄²` `= +(8π/3) r⁴ (ε+p) e^{−2ν} ω̄²` | frame-dragging source | `ε`, `p`, `e^{−2ν}`, `ω̄` | km³·km⁻¹·km⁻² = 1 ✓ | `−(4π/3) r³ (ε+p) ω̄²/(1 − 2m/r)` (`:1204`) | **wrong**: sign, factor 2, one power of `r` (dimension km⁻¹), metric factor |

**6.2 Equation (H2) — `dp₀*/dr`, dimension of every term: km⁻¹.**

| Literature term (H67 100) | Symbol definition | CompactStar background needed | Dimension | Candidate counterpart | Candidate verdict |
|---|---|---|---|---|---|
| `−(m₀R²/(R−2M)²)(8πP + 1/R²)` | gravitational pull of the mass perturbation | `m₀`, `r`, `m`, `p` | km·1/km² = km⁻¹ ✓ | `−m₀(ε+p)/(r²(r−2m))` (`:1184`) | **wrong**: dimension km⁻⁵·km = km⁻⁴ ≠ km⁻¹ (an extra `1/r`), spurious `(ε+p)`, GR factor `(1+8πr²p)/(r−2m)` absent |
| `−4π(E+P)R²p₀*/(R−2M)` | self-coupling | `ε`, `p`, `r`, `m` | km⁻²·km²·1/km = km⁻¹ ✓ | `−p0[4π(ε+p)r + m/r²]/(r−2m)` (`:1184`) | **wrong**: dimension km⁻² (extra `1/r`); the `m/r²` piece belongs to the `δp` formulation's `(1+dε/dp)ν'` term and is incomplete there too |
| `+(1/12)(R⁴/(R−2M))j²(ω̄_R)²` `= (1/12) r³ e^{−2ν} (ω̄')²` | rotational-energy source | `e^{−2ν}`, `ω̄'` | km³·km⁻⁴ = km⁻¹ ✓ | `(1/12) r (1 − 2m/r)(ω̄')²` (`:1216`) | **wrong** structure (`r`, `(1−2m/r)` instead of `r³e^{−2ν}`); dimension happens to be km⁻³·km = km⁻² … not km⁻¹ |
| `+(1/3)(d/dR)[R³j²ω̄²/(R−2M)]` `= (1/3)(d/dr)[r² e^{−2ν} ω̄²]` | centrifugal / frame-dragging source | `e^{−2ν}`, `ν'`, `ω̄`, `ω̄'` | `d/dr`[km²·km⁻² = 1] = km⁻¹ ✓ | `(1/3) r ν' ω̄²` (`:1216`) | **wrong**: not the derivative term; dimension km·km⁻¹·km⁻² = km⁻² |

**6.3 Dimensional audit — result.** Every term of (H1) is dimensionless and every term of (H2)
is km⁻¹ **in the governed system**; the audit passes. In the candidate, three of the seven
terms are dimensionally inconsistent with their own equation; under `GOVERNANCE.md` §3 that
alone fails closed, independent of any physics argument.

**6.4 (H3) and (H4).** `dh₀/dr` is km⁻¹ term by term (same three terms as (H2) with the signs
of the first two reversed and the derivative term absent). (H4) is a dimensionless identity.
`h₀` is fixed only up to the constant `γ` from the interior; HT68 (17a)–(17b) fix it by
continuity with the exterior solution — this is the one place in the `l = 0` sector where a
*surface* matching determines a *constant*, and it concerns `h₀`, never `p₀*` (§8).

## 7. Normalization by `Ω²`

(H1)–(H2) are linear in `(m₀, p₀*)` with sources quadratic in `ω̄`. Writing the verified
first-order response `ω̄ = Ω s(r)`, `ω̄' = Ω s'(r)` (`HartleFirstOrderResponse`, seed-free:
`RotationSolver.hpp:159`, `RotationSolver.cpp:786-802`) and defining

```
m̂₀ ≡ m₀/Ω²  [km³] ,     p̂₀* ≡ p₀*/Ω²  [km²] ,     Ω ≡ Ω_geom
```

the system becomes, with the derivative term expanded analytically
(`r³j²/(r−2m) = r²e^{−2ν}` ⇒ `(1/3)(d/dr)[r²e^{−2ν}s²] = (2/3) r e^{−2ν} s (s + r s' − r ν' s)`):

```
dm̂₀/dr  = 4π r² (ε+p) (dε/dp) p̂₀*  +  (1/12) r⁴ e^{−2ν} (1 − 2m/r) s'²  +  (8π/3) r⁴ (ε+p) e^{−2ν} s²          (N1)

dp̂₀*/dr = − m̂₀ (1 + 8π r² p)/(r − 2m)²  −  4π (ε+p) r² p̂₀*/(r − 2m)
           + (1/12) r³ e^{−2ν} s'²  +  (2/3) r e^{−2ν} s ( s + r s' − r ν' s )                                  (N2)
```

Dimensions: (N1) km² term by term (`dm̂₀/dr`), (N2) km term by term (`dp̂₀*/dr`) ✓. **No
seed appears:** `s`, `s'` are the published seed-free arrays, `e^{−2ν}` is the profile's
`MetricNu`, `ν'` its `MetricNuPrime`, `1 − 2m/r` its `Mass`/`Radius`, `dε/dp` the EOS's. Because
the fixed-`ε_c` solution (§8) is the unique regular solution of a linear system whose only
inhomogeneity scales as `Ω²`, `(m̂₀, p̂₀*)` are **independent of `Ω`** — exactly as `I` is at first
order — and every physical O(Ω²) quantity is `Q = Q̂ Ω_geom²`. Derived normalized fields:

```
δp̂₀ = (ε+p) p̂₀*                      [1]
ξ̂₀  = p̂₀*/ν'                         [km³]
δM̂  = m̂₀(R_*⁻) + 4π R_*² ε(R_*) ξ̂₀(R_*) + I²/R_*³      [km³]        (§12)
ĥ₀  = γ̂ − p̂₀* + (1/3) r² e^{−2ν} s²    [km²], γ̂ from the exterior matching (§12; optional field)
```

ADR-0006 P9 offered two routes to a seed-free second order — (a) from `ω̄_phys` at a requested
`Ω`, or (b) coefficients per `Ω_geom²` from the normalized response. Route (b) is the one that
makes the O(Ω²) product a property of the star, computed once, materialized at any `Ω` by a
multiplication (§19); route (a) would re-solve the ODE per `Ω` for a result that is provably a
multiple of the same coefficient.

## 8. The fixed-central-density family and the centre condition

**8.1 Solution space.** The homogeneous part of (N1)–(N2) (sources off) is a linear
second-order system with two independent solutions near `r = 0`:

- one with `m̂₀(0) ≠ 0`: then `dp̂₀*/dr ≈ −m̂₀/r²` diverges — **irregular**, excluded;
- one with `p̂₀*(0) = c₀`, `m̂₀ ≈ (4π/3)(ε_c+p_c)(dε/dp)_c c₀ r³` — **regular**. Since
  `δε_c = (ε_c+p_c)(dε/dp)_c c₀` (§4.2), this solution is the derivative of the non-rotating
  TOV solution with respect to its central density: **the regular homogeneous solution is a
  shift along the non-rotating sequence.**

Regularity therefore leaves exactly one free constant, `c₀ = p̂₀*(0)`, and it is fixed by
*which* non-rotating star the rotating one is compared with.

**8.2 Fixed `ε_c`.** H67 requires `ξ(0) = 0` "if the rotating configuration is chosen to have
the same central density as the non-rotating configuration" (p. 1009), states the centre
conditions (108) with "these boundary conditions guarantee that the central density of the
rotating and non-rotating configurations are the same" (p. 1022), and HT68 §II f repeat it:
"both `m₀` and `p₀*` vanish at the origin. With these boundary conditions, the rotating star
will have the same central density as the non-rotating one". In the chosen variable:

```
fixed ε_c  ⇔  p̂₀*(0) = 0   and   m̂₀(0) = 0        (no homogeneous admixture)
```

Together with regularity this determines the solution **uniquely; no surface condition exists,
and none may be imposed.** The candidate's `p0(R) = 0` shooting (§16) selects `c₀ ≠ 0` — a
different member of the family (the central density is shifted by
`δε_c = (ε_c+p_c)(dε/dp)_c c₀`, measured at `−0.47 p_c`… in the 4A-0 record) — and is thereby
not the fixed-`ε_c` response INV-09 needs. Where a surface matching *does* determine a constant
in the primary formalism it is `h_{0c}` in the `(m₀, h₀)` formulation (HT68 17a–b; §6.4) — a
metric normalization, not a matter condition. The conflation of the two is the likeliest origin
of the candidate's condition.

**8.3 A (fixed `ε_c`) versus B (fixed baryon number) — kept separate.** Phase 5 compares the
star with itself as it spins down at constant total baryon number `A`. That member of the
family is `(particular) + α·(homogeneous)` with `α = −δA_part/δA_hom` — a *linear combination
formed after both solutions exist*, i.e. INV-09's `Z_i = A_i − B_i(A_B/B_B)` in the language of
sequence derivatives. Phase 4 delivers the particular (fixed-`ε_c`) solution; the homogeneous
solution is available from the same solver at zero extra derivation and is offered as an
optional field (§19, owner question Q8) because it *is* the sequence derivative `∂/∂ε_c`
INV-09 calls `B_i`. **Constant `A` is never imposed inside the Hartle solve.**

## 9. Central regularity series and the numerical start

**9.1 First-order input near the centre.** From the production first-order equation
(`ω̄'' = −4ω̄'/r + [16π(ε+p)/(1−2m/r)]ω̄ + …`, `HARTLE_MOMENT_INERTIA.md` §3): with `ω̄ = ω̄_c(1 + a r²)`,
`2a + 8a = 16π(ε_c+p_c)` ⇒

```
s  = s_c [ 1 + (8π/5)(ε_c+p_c) r² + O(r⁴) ] ,      s' = (16π/5)(ε_c+p_c) s_c r + O(r³) ,      j² = j_c² [1 + O(r²)] .
```

**9.2 Leading powers (derived from (N1)–(N2), then compared with H67 (108)).** Try
`p̂₀* = a₂r² + …`, `m̂₀ = b₅r⁵ + …`, with `m ≈ (4π/3)ε_c r³`, `ν' ≈ (4π/3)(ε_c+3p_c) r`:

- (N2): `2a₂r = (2/3) r j_c² s_c² + O(r³)` (the `m̂₀`, `p̂₀*` and `s'²` terms are O(r³) or O(r⁵))
  ⇒ **`a₂ = (1/3) j_c² s_c²`**.
- (N1): `5b₅r⁴ = 4πr²(ε_c+p_c)(dε/dp)_c a₂ r² + (8π/3) r⁴ (ε_c+p_c) j_c² s_c² + O(r⁶)`
  ⇒ **`b₅ = (4π/15)(ε_c+p_c)[(dε/dp)_c + 2] j_c² s_c²`**.

Both reproduce H67 (108) exactly (with `ω̄_c → Ω s_c`):

```
p̂₀* = (1/3) j_c² s_c² r² + O(r⁴)             m̂₀ = (4π/15)(ε_c+p_c)[(dε/dp)_c + 2] j_c² s_c² r⁵ + O(r⁷)
δp̂₀ = (ε_c+p_c)(1/3) j_c² s_c² r² + O(r⁴)     ξ̂₀ = p̂₀*/ν' = j_c² s_c² r/[4π(ε_c+3p_c)] + O(r³)   (→ 0 at the centre, as H67 (8) requires)
```

**9.3 Numerical start at `r₀ = 1e-5 km`** (INV-05: production profiles begin at `r₀`, never at 0;
the first-order solve starts there with `ω̄'(r₀) = 0`).

| Option | Start values at `r₀` | Error committed | Assessment |
|---|---|---|---|
| A — exact limiting values | `m̂₀ = 0`, `p̂₀* = 0` | the homogeneous solution with `p̂₀*(r₀) = −(1/3)j_c²s_c²r₀²` is excited: relative error on `p̂₀*(R_*)` of order `(r₀/R)² × O(1) ≈ 6e-13`; the central density is perturbed by `δε_c/ε_c ≈ (dε/dp)_c (1+p_c/ε_c) × (1/3)j_c²s_c²r₀²/1 ≈ 1e-12` | acceptable, must be recorded; this is what the candidate does (`:1287`) — an *approximation*, not a defect |
| **B — regular series at `r₀`** | `p̂₀*(r₀) = (1/3)j(r₀)²s(r₀)²r₀²`, `m̂₀(r₀) = (4π/15)(ε+p)(r₀)[(dε/dp)(p_c) + 2] j(r₀)² s(r₀)² r₀⁵` | O(r₀⁴/R⁴) ≈ 4e-25 | **recommended (Q4)**: two one-line expressions from quantities the solver already holds; makes "fixed `ε_c`" exact to rounding; removes a documented approximation instead of inheriting it |

Neither option needs the `r < 1e-10` guard the candidate carries (`:1158`), which is never
reached on a production grid.

## 10. EOS derivative authority

**10.1 Which derivative.** H67 (97) and HT68 (7b), (15a) use `dE/dP`: the derivative of the
one-parameter equation of state `E = E(P)` (HT68 (1)) — the **barotropic** `dε/dp = 1/c_s²`
(`c_s²` the adiabatic sound speed of the cold EOS; adiabatic index `Γ = (ε+p)/p · dp/dε`). Not
`dρ/dp` (rest-mass density; HT68 (27) relates them through `dN/N = dE/(E+P)` for isentropic
matter) and not a profile derivative. Phase 5 additionally needs `dn_i/dp = d(Y_i n_B)/dp` on the
same footing (§15). Incompressible matter has `dε/dp = 0` (§17: the 4A-0 correction of the old
"→ ∞" wording stands); the causal limit is `dε/dp = 1`.

**10.2 Audit of the EOS layer at `bb073c8`** (every path checked; a repository-wide search for
`gsl_spline_eval_deriv` / `gsl_interp_eval_deriv` finds nothing).

| Fact | Citation |
|---|---|
| The EOS is tabulated as `EOSTable{eps, pre, rho, rho_i[], extra_labels}` | `TOVSolver.hpp:153-166`, instances `:470-471` |
| **One** interpolation type for every EOS spline: `TOV_gsl_interp_type = gsl_interp_steffen` (monotone, C¹) | `TOVSolver.hpp:488`; splines built at `TOVSolver.cpp:595,596,622` (`eps(p)`, `rho(p)`, `rho_i(p)`) |
| All evaluations go through `SafeSplineEval` (clamps to the table range, warns) → `gsl_spline_eval_e` | `TOVSolver.cpp:40-92`; `GetEDens` `hpp:739`/`cpp:969`; `GetRho` `hpp:754`/`cpp:1437`; `GetRho_i` `hpp:770`/`cpp:1451` |
| **No derivative facility exists**: no `dε/dp`, `dp/dε`, `c_s²`, `Γ`, and no `dn_i/dp` | whole of `CompactStar/EOS/**`, `Core/TOVSolver.*`, `Core/NStar.*`, `Core/StarProfile.*` |
| The only `dε/dp` in the tree is the candidate's centred finite difference on the **profile** with a `1.0` fallback | `RotationSolver.cpp:1254-1277`; consumed at `:1147`, `:1170` |
| The only `dn_i/dp` in the tree is a profile finite difference in an **uncompiled** file | `Physics/Evolution/src/RotochemicalCache.cpp:72-91` (not in any CMake list; INV-09) |
| A generic spline-derivative utility exists and is used for **sequence** derivatives (`dB/dε_c`, `dI/dε_c`) | `Zaki::Vector::DataSet::Derivative`; `Core/src/StarBuilder.cpp:256,259` |
| Stale comments: the spline block still describes "cubic spline with natural boundary conditions"; the cutoff comment says "`1e-5`" where the code uses `1e-15` | `TOVSolver.hpp:494-498`; `TOVSolver.hpp:549-554` vs `TOVSolver.cpp:1206` |
| Species columns are fractions `Y_i`; `n_i = Y_i n_B` | ADR-0001; `StarProfile.hpp:45` comment still says "density" |

**10.3 Options.**

| | Source of `dε/dp` | Assessment |
|---|---|---|
| **(i)** | **derivative of the same Steffen interpolant `ε(p)` that built the star** (`gsl_spline_eval_deriv` on `visi_eps_p_spline`), evaluated by the EOS/TOV layer and carried onto the profile as a column at build time | consistent with the profile by construction (INV-13 spirit: the derivative of the function actually integrated); continuous (Steffen is C¹), piecewise-smooth, O(h²) accurate between knots; no new interpolant, no new authority. **Recommended (Q5)** |
| (ii) | a tabulated `c_s²` column from the EOS source (CompOSE optionally provides one) | not imported today (`Hidden_ImportEOS_Vis` reads `eps, pre, rho, species` only: `TOVSolver.cpp:419-528`); would be a second authority that can disagree with (i) at the table's resolution. Use as a **cross-check** in 4D (line J), not as the source |
| (iii) | profile finite differences (the candidate) | equal to (i) by the chain rule on a barotrope only where `dp/dr ≠ 0`; amplifies grid noise (range `3.8–4.9e3` on the CMF star, 4A-0 §11C); `1.0` fallback silently substitutes the causal limit. **Rejected as authority** |
| (iv) | re-interpolating `ε(p)` with a smoother scheme (cspline/akima) for the derivative only | makes `dε/dp` the derivative of a *different* function from the one the star was built with; needs INV-13 re-governed. Rejected |

**10.4 Prerequisite surfaced.** There is no authoritative EOS-derivative API today. Implementing
option (i) requires one — an EOS/TOV-layer method with explicit units (inside `TOVSolver` the
table is cgs, INV-02: `dε/dp` in `(g cm⁻³)/(dyn cm⁻²) = s² cm⁻²`, made dimensionless by `×c²`
before it leaves the TOV boundary) and a profile column so that `RotationSolver` never touches
the EOS — plus a way for point-constructed stars (`NStar(std::vector<TOVPoint>)`, used by every
analytic test) to supply the column or fail closed. **This is an ADR-0007 implementation
prerequisite (Q5), not a Phase-4C-G change.** Profile finite differencing must not become the
scientific authority merely because the candidate does it.

## 11. Surface-convention audit (mandatory)

**11.1 Facts.** Production terminates the TOV integration when
`p ≤ PressureCutoff() = max(1e-15·p_c, eos_tab.pre[0])` (`TOVSolver.cpp:1204-1208`; loop breaks
at `:2512`, RHS guard `:1372/1384`). For DS(CMF)-1 the table floor dominates:
`p_floor = 3.351885e25 dyn cm⁻²`, `ρ_floor = 1.658808e8 g cm⁻³`, `n_B = 1e-7 fm⁻³`, `p_surf/p_c = 3.2e-10`
(`TOV_REFERENCE.md:249-251`, `:319`). The omitted layer below the floor is `0.094 / 0.064 / 0.054 km`
at `1.0 / 1.4 / 1.6 M☉` and explains the `0.20–0.35 %` radius residual against CompOSE
(`TOV_REFERENCE.md:264-278`). INV-06 defines the surface as that last node. Write **`R_*`** for
the production surface radius and `ε_*`, `p_*` for the floor values there; the physical `p = 0`
surface is `R_* + ΔR` with `ΔR/R_* ≈ 4–7 × 10⁻³`. Scale of the floor matter: `ε_*/⟨ε⟩ ≈ 6e-7`
(`⟨ε⟩ ≈ 2.7e14 g cm⁻³` for the 1.4 M☉ star), `(ε_*+p_*) R_*² ≈ 2e-8` in geometric units.

**11.2 Which `l = 0` relations require `p = 0`, and which tolerate the floor node.**

| Relation | Needs | At the floor node `R_*` | Effect |
|---|---|---|---|
| exterior `ω̄ = Ω − 2J/r³`, `J = r⁴ω̄'/6`, `I` | vacuum (`ε = p = 0`, `j = 1`) | `j' ∝ (ε+p)`: matching error `O((ε_*+p_*)R_*²) ≈ 2e-8` | already validated (2B-4B, 4B); **unaffected** |
| exterior `m₀ = δM − J²/r³` (H67 105) | vacuum | in vacuum (N1) reduces to `dm₀/dr = (1/12)r⁴(6J/r⁴)² = 3J²/r⁴`, so `m₀ + J²/r³` is **constant at every radius where the source is negligible** — including `R_*` — to `O(ε_*/⟨ε⟩)` plus the omitted layer's own contribution `≲ (ε_layer/⟨ε⟩)(dε/dp)_layer(ΔR/R) ≈ 1e-7 × 1e3 × 5e-3 ≈ 5e-7` | **`δM` adequate as-is** (§12) |
| surface mass shell `4πR_*²ε_*ξ₀(R_*)` (§12.2) | the density at the last node | `≈ 3(ε_*/⟨ε⟩)(δR/R)/(δM/M) ≲ 1e-5` of `δM` | negligible, **computed explicitly, never assumed zero** |
| `p₀*(r)`, `m₀(r)`, `ξ₀(r)`, `δp₀(r)` in the interior | nothing at the surface (initial-value problem from the centre) | exact | unaffected |
| `ξ₀(R_*) = p₀*(R_*)/ν'(R_*)` | an isobar | **it is the displacement of the `p = p_floor` isobar**; `p₀*` and `ν'` both vary on the scale `R` through the omitted layer, so the displacement of the true `p = 0` surface differs by `O(ΔR/R) ≈ 4–7e-3` **relative** | bounded systematic on `δR`, comparable to the existing radius residual; recorded, not corrected |
| `h_{0c}` (HT68 17b) | exterior `h₀` at the matching radius | same `O(ε_*)` argument as `δM` | adequate (optional field) |
| `δA`, `δN_i` boundary term `4πR_*²e^{λ}n(R_*)ξ₀(R_*)` (§15) | the number density at the last node | `n_B = 1e-7 fm⁻³` ⇒ `≲ 1e-6` of the interior integral | negligible, **computed explicitly** |
| omitted-layer contribution to `δA`, `δN_i` | — | `≲ (n_*/⟨n⟩)(dε/dp)_layer(ΔR/R) ≈ 1e-8` | negligible |

**11.3 Classification: `SURFACE ADEQUATE AS-IS`**, with a **semantics contract** rather than a
correction: every O(Ω²) surface quantity is defined at `R_*`, the production termination node;
the exterior matching identities are exact in vacuum and hold at `R_*` to better than `1e-6`
relative; `ξ₀(R_*)` is the displacement of the floor isobar and carries an `O(ΔR/R)` relative
systematic with respect to the `p = 0` surface; the mass-shell and number-boundary terms are
evaluated from the actual last-node `ε_*`, `n_*`, so that a future change of the surface
convention changes numbers, never formulas. No extrapolation or surface reconstruction is
required; INV-06 does not need further governance before 4C-I; the dependency order (§20) is
unchanged. **The EOS-floor node is not identified with the physical surface anywhere in the
proposed contract — it is named, and its consequences are bounded.**

## 12. Exterior matching and `δM`

**12.1 Matching.** Outside the matter (`ε = p = 0`, `j = 1`, `ω̄ = Ω − 2J/r³`) (N1) integrates to
`m₀ = δM − J²/r³` (H67 105) and (N2)/(H3) to `h₀ = −δM/(r−2M) + J²/[r³(r−2M)]` (H67 106),
where `δM` is by definition the change of the Schwarzschild mass read from `g_tt → −(1 − 2(M+δM)/r)`
(H67 92): **`δM = m₀(a) + J²/a³`** (H67 107; HT68 16) is the quantity CompactStar must call
`delta_M`. The candidate's `delta_M = m0[-1]` (`:1384`) omits `J²/R³`.

**12.2 Surface mass shell.** (N1)'s first term is `4πr²(ε+p)(dε/dp)p₀* = −4πr²ξ₀ dε/dr`
(by (99) and the chain rule). Where the density is discontinuous at the surface —
`ε(R_*⁻) = ε_* → 0` — this term carries a delta function and `m₀` jumps by `4πR_*²ε_*ξ₀(R_*)`:
the mass of the displaced surface shell (the same term is the *entire* Newtonian `δM` of a
homogeneous star, H67 (18)). On production stars `ε_*` is the floor density and the shell is
`≲ 1e-5` of `δM` (§11.2); on the constant-density analytic star it is the dominant term. The
governed definition therefore keeps it explicit:

```
δM̂ = m̂₀(R_*⁻) + 4π R_*² ε_* ξ̂₀(R_*) + I²/R_*³           [km³] ,      δM = δM̂ Ω_geom²   [km]
```

using `J/Ω = I` (`J²/R³ = I²Ω²/R³`) — no raw `J`, no seed. Physical `δM/M`, `δM` in `M☉`
(via `SUN_M_KM`) follow by multiplication. Sign and size: HT68 Tables 5–6 give `δM/M` positive,
`0.05–0.5` at `Ω² = M/R³`; the candidate's `δM/M = −1.6` at a raw seed equivalent to
`Ω/Ω_K ≈ 0.4` (4A-0 §12) is therefore unphysical in sign and magnitude — consistent with the
hypothesis classification in §17.

## 13. `ξ₀` and the surface displacement

- **Definition.** `ξ₀ = −δp₀/(dp/dr) = p₀*/ν'` [km] (H67 99; HT68 24b) — the outward
  displacement of the isobar at `r`. **Eulerian ↔ Lagrangian:** `δp₀` is the Eulerian change,
  `ξ₀` is the displacement that makes the Lagrangian change vanish (§4.3); the two are not
  independent, and "`δp₀(R) = 0`" is a statement about a *particular* member of the family, not
  a boundary condition of the formalism.
- **Centre.** `ξ̂₀ = j_c²s_c² r/[4π(ε_c+3p_c)] + O(r³)` → 0 (§9.2), as H67 (8) requires for the
  fixed-`ε_c` family.
- **Surface.** `ξ₀(R_*) = p₀*(R_*)/ν'(R_*)` with `ν'(R_*) = (m + 4πr³p)/(r(r−2m))|_{R_*}` from the
  profile (`→ M/[R_*(R_*−2M)]` as `p_* → 0`, the corrected H67 (117)); generally **nonzero and
  positive** (`δR/R = 0.15–0.33` at `Ω² = M/R³` in HT68 Tables 5–6). Newtonian homogeneous-star
  limit (§21 line F): `ξ̂₀(a) → a⁴/(3M)`.
- **Adjudication of the candidate's condition.** Forcing `δp₀(R) = 0` forces `ξ₀(R) = 0`
  (measured `−5e-10 km`, `1e-16 km`) — no monopole surface displacement — and moves the central
  density (`p0_c = −0.47 p_c, −0.52 p_c`). Against the primary text (§8.2) this is
  **REFUTED as Hartle's condition** and **CONFIRMED as a different, unphysical-for-Phase-5
  member of the one-parameter family**. It is not a matter of convention.
- **Semantics at `R_*`.** Per §11.3, `ξ₀(R_*)` is the floor-isobar displacement, with an
  `O(ΔR/R)` relative systematic against the `p = 0` surface; `δR/R ≡ ξ₀(R_*)/R_*` is reported
  with that caveat.

## 14. `l = 0` versus `l = 2` — decided by the angular integration

**14.1 The baryon integral to O(Ω²).** H67 (109): `A = ∫d³x (−g)^{1/2} u^t 𝔑`. With the metric
(66)–(70) (`k₀ = 0`):

```
−g            = e^{ν_H+λ_H}(1+2h)[1 + 2m/(r−2M)](1+2k)² r⁴ sin²θ          (the ω² cross terms cancel exactly in g_tt g_φφ − g_tφ²)
(−g)^{1/2}    = e^{(ν_H+λ_H)/2}(1+h)(1 + m/(r−2M))(1+2k) r² sinθ + O(Ω⁴)
u^t           = e^{−ν_H/2}[1 − h + (1/2)e^{−ν_H} r² sin²θ ω̄²]              HT68 (5)
(−g)^{1/2}u^t = e^{λ_H/2} r² sinθ [1 + m/(r−2M) + 2k + (1/2) e^{−ν_H} r² sin²θ ω̄²]     — h cancels identically
```

In the label coordinates `r = R + ξ(R,Θ)`, `θ = Θ`, the number density is `𝔑 = N(R)` **exactly**
(a barotrope's `n` is a function of `ε`, and `R` labels constant-`ε` surfaces), the Jacobian is
`(R+ξ)²(1 + ∂_Rξ)`, and `e^{λ_H(r)/2} = e^{λ(R)}(1 + λ'ξ)` (CompactStar `λ`). Hence

```
A = ∫dΩ ∫dR N e^{λ} R² [ 1 + 2ξ/R + ∂_Rξ + λ'ξ + m/(R−2M) + 2k + (1/2)e^{−2ν}R² sin²Θ ω̄² ] + O(Ω⁴) .
```

**Angular integration.** Every `l = 2` piece — `ξ₂P₂`, `m₂P₂`, `k₂P₂` — integrates to zero
(`∫P₂ dΩ = 0`); `k₀ = 0`; `∫sin²Θ dΩ = 8π/3`. **Only `l = 0` survives**, exactly as H67 states
on p. 1022. The `ξ₀` terms combine into a total derivative plus `−N'ξ₀`:

```
N e^{λ} R² [2ξ₀/R + ξ₀' + λ'ξ₀] = (d/dR)[N e^{λ} R² ξ₀] − (dN/dR) e^{λ} R² ξ₀ ,     −(dN/dR)ξ₀ = (dN/dp)(ε+p) p₀* = (dN/dp) δp₀ .
```

Therefore, per unit `Ω²` and on the production star:

```
δÂ = 4π ∫₀^{R_*} dr r² e^{λ} { (dN/dp) δp̂₀ + N [ m̂₀/(r − 2m) + (1/3) r² e^{−2ν} s² ] }  +  4π R_*² e^{λ(R_*)} N(R_*) ξ̂₀(R_*)        (A2)
```

and the unperturbed `A = 4π∫r²e^{λ}N dr` is INV-14 / HT68 (3e). **Cross-check against the
primary source:** forming `μδA − δM` with `N = (E−ε)/μ` and (N1) reproduces HT68 (20) term by
term — `(E+P)p₀*{(dE/dP)[(1−2M/r)^{−1/2} − 1] − (dε/dP)(1−2M/r)^{−1/2}}`,
`(E−ε)(1−2M/r)^{−3/2}[m₀/r + (1/3)j²r²ω̄²]`, `−(1/4πr²)[…]`, `−J²/R³` — with the minus sign on
`dε/dP` (H67 (115) as printed has `+`; §2.3). The derivation is thereby anchored to the primary
formalism, not to the secondary paper.

**14.2 Classification.** For **every scalar count** (`A`, `N_i`, integrated composition) at
O(Ω²): **A — `l = 0` is sufficient**; the `l = 2` sector contributes neither through the volume
integrals nor through the boundary deformation (`ξ₂(R_*)P₂` integrates to zero on the sphere).
The same holds for `M` (H67 §VII) and `I`. The `l = 2` sector is required for the **shape**
(`ξ₂`, ellipticity `ε(R) = −3ξ₂/2R`), the quadrupole moment `Q` (H67 138 as corrected by HT68
fn 5: `Q = J²/M + 8KM³/5`), and `h₂, v₂, m₂` — none of which Phase 5 consumes.

**14.3 The ratified exit criterion.** `MODERNIZATION_ROADMAP.md` Phase 4: "O(Ω) and O(Ω²)
validated and reachable. The candidate status of `675b4a9` is resolved — ratified or replaced."
The candidate is **`l = 0` only** (`ODE_Hartle2_N_Fast` integrates `(m0, p0)`, nothing else;
`RotationSolver.hpp:203-217`), so *replacing* it with the governed `l = 0` system discharges
"the candidate status … resolved". Whether "O(Ω²) validated" is read as the monopole sector
(what the candidate ever claimed) or as the full `l = 0 + l = 2` Hartle–Thorne structure
(shape, `Q`) is **an owner decision** (ADR-0007 Q11, §23 Q6): the recommendation is to read it as
the monopole sector for Phase 4C/4D — sufficient for Phase 5 and for the candidate's
replacement — and to schedule `l = 2` as a separate, later ADR. Nothing here narrows the roadmap
silently: the question is surfaced.

## 15. Phase-5 structural particle-number response interface

**15.1 Every O(Ω²) contribution to a scalar count** (from (A2), per species `i` with
`n_i = Y_i n_B` per ADR-0001 and `N → n_i`):

| # | Contribution | Term in (A2) | Origin | Species-dependent? |
|---|---|---|---|---|
| 1 | composition / density perturbation | `(dn_i/dp) δp̂₀` | Eulerian `δn_i = (ε+p)(dn_i/dp)p₀*` (HT68 7c) | yes — needs `dn_i/dp` from the EOS |
| 2 | proper-volume (radial metric) perturbation | `n_i m̂₀/(r − 2m)` | `g_rr → e^{2λ}[1 + 2m₀/(r−2m)]`, `(−g)^{1/2}` | no — multiplies `n_i` |
| 3 | time-dilation of the rotating fluid | `n_i (1/3) r² e^{−2ν} s²` | `u^t`'s `(1/2)e^{−ν_H}r²sin²θω̄²` after `∫sin²Θ` | no |
| 4 | lapse perturbation `h₀` | — | cancels exactly between `(−g)^{1/2}` and `u^t` | none |
| 5 | moving isobar / boundary | `4πR_*²e^{λ}n_i(R_*)ξ̂₀(R_*)` | integration by parts of the `ξ₀` Jacobian terms; the bulk part became #1 | yes (through `n_i(R_*)`) — `≲ 1e-6` on production stars, non-zero in principle |
| — | `l = 2` (`ξ₂, m₂, k₂, h₂, p₂*`) | 0 | `∫P₂ dΩ = 0` | — |

So `∫(dn_i/dp)δp₀ dV` alone is **not** complete: terms 2, 3 and 5 are of the same order and
carry no EOS derivative. (INV-09's defect list — `A_i` computed from `∫(dn_i/dp)p₀ w_V dr` only,
and never divided by `Ω²` — is thereby confirmed to be two defects, not one.)

**15.2 What Phase 4 must expose** (per node on the profile grid, seed-free, per `Ω_geom²`, fixed
`ε_c`): `m̂₀(r)`, `p̂₀*(r)` (or `δp̂₀ = (ε+p)p̂₀*`), `ξ̂₀(r)`, and — from first order — `s(r)`; plus the
surface scalars `ξ̂₀(R_*)`, `δM̂`, and `R_*`, `ε_*`, `n_B(R_*)`, `e^{λ(R_*)}` for the boundary
terms. Phase 5 supplies `n_i(r) = Y_i n_B`, `dn_i/dp` (EOS-owned, §10) and the proper-volume
weight `w_V = 4πr²e^{λ}` (`GeometryCache::WV`, ADR-0004) and evaluates

```
A_i ≡ (∂N_i/∂Ω_geom²)|_{ε_c} = δN̂_i = ∫₀^{R_*} w_V { (dn_i/dp) δp̂₀ + n_i [ m̂₀/(r−2m) + (1/3) r² e^{−2ν} s² ] } dr + w_V(R_*) n_i(R_*) ξ̂₀(R_*)       [km²]
```

(`× 1e54` for a dimensionless count with `n_i` in fm⁻³, INV-14; `÷ c²` for a per-`Ω_phys²`
coefficient). `B_i = (∂N_i/∂ε_c)|_Ω` and `Z_i = A_i − B_i(A_B/B_B)` remain Phase 5 (INV-09); the
optional homogeneous field of §8.3/§19 gives `B_i` from the same solver, and the identity
"homogeneous `δM̂_hom` ∝ `dM/dε_c` from the TOV sequence" is a 4D consistency line (§21, K).
The ADR establishes the **structural response fields Phase 5 may trust**; it decides nothing
about weak-reaction rates or the chemical-imbalance convention (INV-11).

## 16. Candidate archaeology — discrepancy table (written only after §3–§15)

The candidate (`ODE_Hartle2_N_Fast` `RotationSolver.cpp:1140-1226`; `SolveHartle2_N` `:1242-1386`;
`HartleResult` `RotationSolver.hpp:203-217`) is compared with the governed system. **This table
confers no authority on the candidate; a matching term is archaeology, not ratification**
(`GOVERNANCE.md` §5).

| Feature | Governed derivation (§) | Candidate (`RotationSolver.cpp` at `bb073c8`) | Verdict |
|---|---|---|---|
| cited authority | H67 (97), (100), (107), (108); HT68 (15a–c), (16) | comment "Hartle (1967), Eqs. (30)–(33)" (`:1129`), "Eq. (33)" (`:1203`), "Eq. (30)" (`:1211`); "`[FIX: confirm exact from textbook]`" (`:1135`); "`???`" (`:1207`) | **wrong citations**: H67 (30)–(33) are the general stationary metric and the first-order `R_φ^t` identity, not the `l = 0` system |
| integrated variable | Hartle's dimensionless `p₀*` (§4) | `p0` annotated `delta_p(r) [km^-2]` (`hpp:207`), i.e. Eulerian `δp₀` by intent | different variable; its ODE is not the correct equation for **either** (§6) |
| dimensions | (H1) dimensionless, (H2) km⁻¹ term by term (§6.3) | three of seven terms dimensionally inconsistent (`:1184`, `:1204` 2nd term, `:1216` both terms) | **fails closed** (`GOVERNANCE.md` §3) |
| `dm₀/dr` term 1 | `4πr²(ε+p)(dε/dp)p₀*` | `4πr²(dEdP)p0` (`:1170`) | form-equivalent for `p0 ≡ δp₀`; derivative source rejected (§10) |
| `dm₀/dr` term 2 | `(1/12)e^{−2ν}(1−2m/r)r⁴ω̄'²` | `(1/12)r⁴ω̄'²/(1−2m/r)` (`:1204`) | **wrong** (`j²` inverted, `e^{−2ν}` absent) |
| `dm₀/dr` term 3 | `+(8π/3)r⁴(ε+p)e^{−2ν}ω̄²` | `−(4π/3)r³(ε+p)ω̄²/(1−2m/r)` (`:1204`) | **wrong** (sign, factor 2, power of `r`, metric factor) |
| `dp₀*/dr` homogeneous terms | `−m₀(1+8πr²p)/(r−2m)²`, `−4π(ε+p)r²p₀*/(r−2m)` | `−p0[4π(ε+p)r + m/r²]/(r−2m) − m0(ε+p)/(r²(r−2m))` (`:1184`) | **wrong** (both dimensionally; GR factors absent) |
| `dp₀*/dr` source terms | `(1/12)r³e^{−2ν}ω̄'² + (1/3)(d/dr)[r²e^{−2ν}ω̄²]` | `(1/12)r(1−2m/r)ω̄'² + (1/3)rν'ω̄²` (`:1216`) | **wrong** (structure; second term dimensionally) |
| `j²` | `e^{−2ν}(1−2m/r)`, `e^{−2ν}` from `MetricNu` | comment "we don't have nu directly in the fast cache" (`:1194`); `1/(1−2m/r)` used instead | **false premise** (profile carries `MetricNu`, `StarProfile.hpp:246`) |
| `dε/dp` | EOS derivative, EOS-owned (§10) | profile centred differences, `1.0` fallback (`:1254-1277`) | **rejected as authority**; fallback is the causal limit, not incompressible |
| centre condition | `m̂₀(0) = p̂₀*(0) = 0`; series start at `r₀` recommended (§9) | `{0,0}` particular (`:1287`) — acceptable as an approximation — **plus** `{0,1}` homogeneous (`:1321`) admixed by shooting | admixture **wrong** for fixed `ε_c` (§8) |
| surface condition | none exists for `p₀*` (§8.2) | `p0(R) = 0` shooting (`:1349-1355`), `p0_c = −p0_part(R)/p0_hom(R)` | **REFUTED** (§13); forces `ξ₀(R) = 0`, shifts `ε_c` |
| `ξ₀` | `p₀*/ν' = −δp₀/(dp/dr)` | `−p0/(dp/dr)`, `dp/dr = −(ε+p)ν'` (`:1368-1380`) | form correct for `p0 ≡ δp₀`; value wrong because `p0` is |
| `δM` | `m̂₀(R_*⁻) + 4πR_*²ε_*ξ̂₀(R_*) + I²/R_*³` (§12) | `m0[-1]` (`:1384`) | **incomplete** (`J²/R³`, shell absent); sign/magnitude unphysical (4A-0 §12) |
| normalization | per `Ω_geom²` from `s`, `s'` (§7) | reads `stored_omega_bar_`, `stored_domega_bar_` (`:1300-1301`) written at the raw seed (`:753-754`); every output `∝ seed²` (4A-0 §12) | **violates ADR-0006 P9** |
| superposition mechanism | the homogeneous solution *is* the sequence derivative and may be exposed (§8.3) | used only as a shooting device | mechanism sound, purpose wrong |
| reachability | new seed-free product through `NStar` (§19) | `SolveHartle2_N`, `GetHartleResult` public, zero callers, executable from user code (4A-0 §13) | disposition §19 |
| MixedStar second order | out of scope (ADR-0004 §0-Q2 track) | empty stubs (`:1230-1236`, `:1390-1393`) | remove with the candidate (§19) |

## 17. The 4A-0 findings, classified against the primary source

| 4A-0 hypothesis (`PHASE4_ROTATION_ENTRY.md` §10–§11) | Classification | Primary evidence |
|---|---|---|
| the transcribed system (H1)–(H2) with `dj²/dr = −8πr²(ε+p)j²/(r−2m)` | **CONFIRMED** | H67 (97), (100), (62); HT68 (15a–b) |
| `p₀*` is dimensionless; `δp₀ = (ε+p)p₀*` | **CONFIRMED** | H67 (87)–(88), (99); HT68 (7a) |
| `ξ₀ = −δp₀/(dp/dr) = p₀*/ν'` | **CONFIRMED** | H67 (99) (+ (117) corrected, §2.3); HT68 (24b) |
| `δM = m₀(R) + J²/R³` | **CONFIRMED** (and completed by the explicit shell term) | H67 (105)–(107); HT68 (15c)–(16); §12.2 |
| centre series `p₀* ≈ (1/3)j_c²ω̄_c²r²`, `m₀ ≈ (4π/15)(ε_c+p_c)[(dε/dp)_c+2]j_c²ω̄_c²r⁵` | **CONFIRMED** (re-derived §9.2) | H67 (108) |
| fixed-`ε_c` family: `m₀(0) = p₀*(0) = 0`, no surface shooting | **CONFIRMED** | H67 p. 1009, p. 1022; HT68 §II f |
| `j²` factor defect in `S_m` (reciprocal, `e^{−2ν}` absent) | **CONFIRMED** | H67 (97) |
| `S_m` second term: sign / factor 2 / power of `r` / metric factor | **CONFIRMED** | H67 (97) with (62) |
| homogeneous `p0`, `m0` coefficients off by `1/r`; GR factors absent | **CONFIRMED** | H67 (100) |
| `S_p` terms structurally wrong | **CONFIRMED** | H67 (100) |
| `δp(R) = 0` shooting is not Hartle's condition; forces `ξ₀(R) = 0`; shifts `ε_c` | **CONFIRMED** | §8.2, §13 |
| candidate's `p0` is Eulerian `δp` by intent | **CONFIRMED** (code evidence; convention-independent once §4 is fixed) | `hpp:207`, `:1170`, `:1377` |
| `1.0` fallback is the causal limit; incompressible ⇒ `dε/dp → 0` | **CONFIRMED** (physics; HT68 (27) consistent) | — |
| `dε/dp` should be the EOS derivative | **CONFIRMED** | H67 (97) "`dE/dP`"; HT68 (1), (7b) |
| (T2), the `δp₀` form of (H2) | **CONFIRMED** (algebraic consequence of (H2) and TOV) | — |
| literal-zero start: relative truncation `≈ 6e-13`, acceptable | **PARTLY CORRECT** — the bound is confirmed, but the series start is now recommended over the approximation (§9.3) | H67 (108) |
| `δM/M = −1.6` unphysical; physical `δM > 0`, small | **CONFIRMED** | HT68 Tables 5–6 |
| `l = 0` suffices for `M` and for `N_i` | **CONFIRMED** (derived §14.1) | H67 §VII, p. 1022 |
| equation-number citations in the candidate ("(30)–(33)") | **REFUTED** | H67 pp. 1012–1013 vs 1020–1021 |
| the `R`-versus-`r` coordinate question (not raised in 4A-0) | **CONVENTION-DEPENDENT**, resolved: same to O(Ω²) | HT68 fn 4 |
| "no available `e^{−2ν}`" premise in the candidate | **REFUTED** | `StarProfile.hpp:246` |

## 18. `GOVERNANCE.md` §3.1 — applicability, condition by condition

Replacing the candidate is a scientific-semantic change with **no** prior baseline (no
`SolveHartle2_N` output has ever been baselined; none exists under `tests/baselines/`). The
normal rule "creating the baseline is the task" would freeze the candidate's numbers — the
exact situation §3.1 exists for.

| §3.1 condition | Status after ADR-0007 acceptance | How ADR-0007 records it |
|---|---|---|
| 1 an ACCEPTED ADR identifies a specific current behaviour as invalid / inconsistent | satisfiable — ADR-0007 names `ODE_Hartle2_N_Fast` / `SolveHartle2_N` and the defects of §16 (dimensional inconsistency, wrong sources, non-Hartle condition, incomplete `δM`, seed² normalization) | §1 "invalid / superseded behaviour" |
| 2 capturing it as golden would enshrine the rejected behaviour | yes — `δM/M = −1.6`, `ξ₀(R) = 0`, `ε_c` shifted, seed-quadratic | §2 "why capturing is forbidden" |
| 3 the ADR explicitly identifies the minimum correction | the governed `l = 0` system (N1)–(N2), fixed-`ε_c` centre condition with series start, EOS-owned `dε/dp`, `δM` with shell and `I²/R³`, exposure per `Ω_geom²`; **no `l = 2`, no rotochemical, no MixedStar** | §3 "minimum correction" |
| 4 independent verification that does not depend on the superseded output | §21 lines A–K: analytic limits, series, an independent solver in different variables, conservation identity (H4), published HT68 tables, convergence, EOS-derivative sensitivity | §4 |
| 5 narrowly scoped to the governed defect | the replacement touches only the O(Ω²) path and its result type; first order, `I`, the seven goldens stay bitwise | §5 |
| 6 records which historical outputs are no longer reference results | every value ever produced by `SolveHartle2_N` / `GetHartleResult` (including the diagnostic numbers in `PHASE4_ROTATION_ENTRY.md` §12) is **not** a reference | §6 |
| 7 a baseline is created immediately after correction + independent validation | 4D ends by baselining the governed monopole response on DS(CMF)-1 (proposed `hartle_monopole_dscmf1_debug.tsv`) | §7 |

**Conclusion: §3.1 is applicable and every condition can be discharged; it is invoked by
ADR-0007, not by this record.** It is not invoked casually: without it, condition 6 of §3 would
force a golden of numbers already adjudicated invalid.

## 19. Candidate API disposition, result type, materialization

**19.1 Disposition of the public candidate API** (`SolveHartle2_N`, `GetHartleResult`,
`HartleResult`, `ODE_Hartle2_N_Fast`, `ODE_Hartle2_Mixed_Fast`, `SolveHartle2_Mixed`,
`include_m0p0_source_`, `fast_dEdP*`, `fast_nu`; zero callers).

| | Option | Assessment |
|---|---|---|
| **A** | **atomic replacement** in 4C-I: the governed implementation and its new result type land in the same commit that deletes the candidate functions and `HartleResult` | zero callers ⇒ zero churn; the name `HartleResult` (which has meant three different things since `675b4a9`) retires, so no reader can mistake old semantics for new; the §3.1 record and the deletion are one change. **Recommended (Q13)** |
| B | fail closed now (make `SolveHartle2_N` return `valid = false` / throw) until replacement | requires a production change **in this task**, which is forbidden; remains available as a one-line 4C-I-0 step if implementation is delayed after acceptance |
| C | retain the API name, replace internals | invites continuity of a struct whose fields (`p0 [km⁻²]`, `p0_c`, `delta_M` without `J²/R³`) do not exist in the governed contract; rejected |

Until 4C-I lands, the candidate stays exactly as it is: publicly callable, zero callers, marked
`UNVERIFIED SCIENTIFIC CANDIDATE` in source (`RotationSolver.hpp:187-202`) and in every document.

**19.2 Result type (concept; names indicative).** Seed-free coefficients per `Ω_geom²` at fixed
`ε_c`, one canonical geometric representation (ADR-0006 Q3 precedent), no physical/geometric
duplication:

```
struct HartleMonopoleResponse {                    // fixed eps_c, per Omega_geom^2, seed-free
  DataColumn m0_over_Omega2;        // [km^3]   m̂₀(r)
  DataColumn p0star_over_Omega2;    // [km^2]   p̂₀*(r)            (canonical variable, Q1)
  DataColumn delta_p0_over_Omega2;  // [1]      (eps+p) p̂₀*         (derived; Phase-5 convenience)
  DataColumn xi0_over_Omega2;       // [km^3]   p̂₀*/nu'
  double deltaM_over_Omega2;        // [km^3]   m̂₀(R*) + 4π R*² ε* ξ̂₀(R*) + I²/R*³
  double shell_over_Omega2;         // [km^3]   the explicit surface-shell term (expected ≲ 1e-5 δM̂ on EOS-floor stars)
  double xi0_R_over_Omega2;         // [km^3]   floor-isobar displacement (§11.3 semantics)
  const DataColumn *r_grid;  bool valid;
  PhysicalMonopoleRotation At(AngularVelocity) const;   // Q_phys = Q̂ · Ω_geom²; no new solve
};
```

Optional (owner question Q8): a parallel `HartleMonopoleHomogeneous` (per unit `p₀*(0)`) — the
TOV sequence derivative from the same solver, for Phase 5's `B_i` and for the 4D identity K.
`h₀` (needs the exterior constant) is derivable algebraically and is **not** required for `A_i`;
it is left out unless Phase 5 shows a need. Nothing from the candidate sits behind these
accessors.

**19.3 Physical materialization.** Option **B** — coefficients plus `At(AngularVelocity)` —
parallel to `HartleFirstOrderResponse::At` (`RotationSolver.cpp:1083-1122`): `m₀ = m̂₀Ω²` [km],
`p₀* = p̂₀*Ω²` [1], `δp₀ = δp̂₀Ω²` [km⁻²], `ξ₀ = ξ̂₀Ω²` [km], `δM = δM̂Ω²` [km], `δR/R`, `Ω/Ω_K`,
from **one** coefficient source, by multiplication, with the ADR-0006 binding clarification
extended verbatim: no `NStar` acquires an O(Ω²) physical perturbation without an explicit
`AngularVelocity`. Zero spin materializes zeros.

## 20. Implementation dependency map (proposed; ratified item list unchanged)

```
4C-G  this record + ADR-0007 PROPOSED                         (documentation)                 ← done
  │
  ▼  owner adjudicates ADR-0007 (Q1–Q13)                                                     ← next
  │
4C-I-0  EOS-derivative authority: EOS/TOV-layer dε/dp (and dn_i/dp) from the Steffen
        interpolant; profile column; fail-closed path for point-built stars; unit contract     (engineering + numerical rationale; own tests)
  │
4C-I    governed l = 0 solver consuming HartleFirstOrderResponse (s, s') + profile columns
        + EOS derivative; HartleMonopoleResponse; NStar accessor; atomic candidate removal    (scientific-semantic under §3.1; I and goldens bitwise)
  │
4D      independent validation A–K (§21), detectors (§22), then the first monopole baseline   (validation; §3.1 condition 7)
  │
4E      Phase-5 structural fields: δN̂_i machinery on n_i = Y_i n_B with EOS dn_i/dp           (structural)
  │
Phase 5 (INV-09 Z_i, INV-11 …)
```

Surface convention: **no prerequisite** (§11.3). `l = 2`: separate later ADR if the owner reads
the exit criterion that way (Q11). MixedStar: separate track (ADR-0004 §0-Q2).

## 21. Validation plan for 4D (plan only; nothing implemented, no candidate oracle)

| Line | Check | Analytic / external | Predeclared bound (basis) |
|---|---|---|---|
| A | constant-density Schwarzschild interior (the 4B fixture): `dε/dp = 0` inside ⇒ `m̂₀` from the rotational sources only; the whole matter `δM` sits in the explicit shell term; `p̂₀*`, `ξ̂₀` node-by-node against line C | analytic background, numerical reference | to be predeclared at 4D entry from the reference solver's measured floor (4B precedent: `1e-7` analytic) |
| B | regular-centre series (§9.2): `p̂₀*/r² → (1/3)j_c²s_c²`, `m̂₀/r⁵ → b₅`, `ξ̂₀/r → …` as `r → 0` | analytic | relative deviation `∝ r²` with the fitted coefficient of the correct sign; agreement at the first ten nodes `≤ 1e-8` |
| C | independent solver in **different variables**: the Eulerian pair `(m₀, δp₀)` via (T1)–(T2) (`PHASE4_ROTATION_ENTRY.md` §10.3) and/or the `(m₀, h₀)` pair via (H1)+(H3) with the exterior constant, built on `hartle_reference.hpp`'s background, never on production helpers | analytic + CMF | analytic `1e-7`; CMF `1e-4` (4B precedent), to be re-predeclared at 4D entry |
| D | `Ω²`-coefficient invariance under the internal first-order seed over `[1e-3, 1e3]` (via `RotationSolverTestSeam`) | analytic | `≤ 1e-10` (ADR-0006 §7 item 1; integrator abstol `1e-10`) |
| E | materialized quantities quadratic in the requested `Ω`: `Q(2Ω) = 4Q(Ω)` node by node; zero spin ⇒ zeros | analytic | `≤ 1e-14` (arithmetic) |
| F | Newtonian limit on the homogeneous star as `M/R → 0`: `δM̂ → R³`, `3Mξ̂₀(R)/R⁴ → 1`, `p̂₀* → r²/3` interior, monotone in `M/R` (4B Exp E method) | analytic | deviation `∝ M/R`, sequence over four compactnesses, sign and monotonicity asserted; bound on the extrapolated intercept `≤ 5e-3` |
| G | exterior identity: `m̂₀ + I²/r³` constant across the last nodes where `(ε+p)r² < 1e-6`; `δM̂` independent of the matching node | analytic + CMF | `≤ 1e-6` relative (§11.2 bound) |
| G′ | conservation identity (H4): integrate (H3) alongside and assert `p̂₀* + ĥ₀ − (1/3)r²e^{−2ν}s²` constant | analytic + CMF | `≤ 1e-9` relative to `p̂₀*(R_*)` |
| H | published results: HT68 Table 5 (`δM/M`, `δR/R`, `ω_c/Ω`, `ω_s/Ω`, `R_g/R` at `Ω² = M/R³`, "accurate to about 1 per cent") on the HW EOS printed as HT68 Table 1 — requires transcribing the table into test data; Chandrasekhar & Miller (1974, MNRAS 167, 63) for the homogeneous star — **not authenticated in this task**, to be authenticated before use | external (small, transcribed) | `≤ 2e-2` (HT68's stated accuracy plus interpolation-scheme differences: HT68 used logarithmic interpolation, CompactStar Steffen) |
| I | radial-resolution convergence of `δM̂`, `ξ̂₀(R_*)`, `p̂₀*(R_*)` on `N ∈ {1001, 2001, 4001, 8001}` (analytic) and the CMF resolution ladder of `GRID_CONVERGENCE.md` | analytic + CMF | order consistent with the profile's interpolation order (INV-13); Richardson residual reported |
| J | EOS-derivative sensitivity: Steffen `dε/dp` vs a tabulated `c_s²` column (if the source provides one) vs the profile finite difference — `δM̂` spread reported; the governed source must sit inside the spread and the FD noise must be visible | CMF | diagnostic; predeclared expectation: spread `≤ 1e-3` on `δM̂` |
| K | sequence-derivative identity: the homogeneous solution's `δM̂_hom/p₀*(0)` versus `dM/dε_c` from a TOV central-density sweep (`StarBuilder`'s existing sequence-derivative machinery, `StarBuilder.cpp:256-259`) | CMF | `≤ 1e-3` (finite-difference sweep accuracy) |

Then, and only then, the first monopole baseline (§18 condition 7).

## 22. Detector plan (later; no production mutation in this task)

| Detector | Mutation | Must be caught by |
|---|---|---|
| M1 | drop `e^{−2ν}` from one source term (`j² → 1−2m/r`) | A, C (interior profile), F (wrong Newtonian intercept survives? no — F is insensitive; A/C carry it) |
| M2 | flip the sign of the `(8π/3)r⁴(ε+p)e^{−2ν}s²` term | A, C, B (centre coefficient `b₅` changes sign structure) |
| M3 | omit `I²/R_*³` in `δM̂` | G (constancy of `m̂₀ + I²/r³` fails), H |
| M4 | omit the surface-shell term | A (homogeneous star: `δM̂` loses its dominant term), F |
| M5 | substitute profile finite differences for the EOS derivative | J (spread), I (convergence order degrades), B (centre coefficient noise) |
| M6 | impose `δp₀(R_*) = 0` (the candidate's condition) | B (`p̂₀*(0) ≠ 0`), K, H (`δR/R` vanishes) |
| M7 | start with `{0,0}` at `r₀` after the ADR mandates the series start | B at the predeclared `1e-8` level only if `r₀` is enlarged in the detector (`6e-13` is below the bound at `1e-5 km`) — record as a **weak** detector |
| M8 | let the raw first-order seed leak (use `ω̄_raw` instead of `Ω s`) | D, E |
| M9 | drop the `(1 + 8πr²p)` GR factor | C, H (compact stars), F insensitive |

Each detector maps to a distinct line of §21; a detector that no line catches is a coverage gap
to be closed before 4D is declared complete.

## 23. Owner questions (genuine decisions only)

| | Question | Recommendation, from the primary equations and downstream needs |
|---|---|---|
| **Q1** | Canonical integrated pressure variable: Hartle's `p₀*`, the Eulerian `δp₀`, or the `(m₀, h₀)` pair? | **`p₀*`** — the primary equations are most direct in it (no `(1+dε/dp)ν'` coupling in the homogeneous part), it is dimensionless with a literal-zero centre condition and a clean `r²` series, `δp₀ = (ε+p)p₀*` is a well-conditioned *product* everywhere (the inverse division is ill-conditioned at the surface), `ξ₀ = p₀*/ν'` is a ratio of profile quantities, and (H4) gives a conservation check. `δp₀` stays a derived field; `(m₀, h₀)` is the independent 4D formulation |
| **Q2** | Fixed `ε_c` as the governed family? | **Yes** — `m̂₀(0) = p̂₀*(0) = 0`, no surface shooting (§8); the homogeneous (sequence-derivative) solution optionally exposed for Phase 5 (Q8) |
| **Q3** | Store coefficients per `Ω_geom²` (seed-free), materialize by multiplication? | **Yes** (§7, §19.3) — ADR-0006 P9 route (b) |
| **Q4** | EOS owns `dε/dp` (and `dn_i/dp`)? | **Yes**: derivative of the same Steffen interpolant, evaluated in the EOS/TOV layer, carried on the profile; profile finite differences rejected; API is a 4C-I-0 prerequisite (§10) |
| **Q5** | Does the EOS-floor surface require a Phase-4 correction? | **No** — `SURFACE ADEQUATE AS-IS` with the `R_*` semantics contract, explicit shell and boundary terms, and the recorded `O(ΔR/R)` systematic on `ξ₀(R_*)` (§11) |
| **Q6** | Is `l = 0` sufficient for the governed Phase-5 deliverable, with `l = 2` separately deferred — and does the ratified Phase-4 exit criterion "O(Ω²) validated" mean the monopole sector? | **`l = 0` is sufficient for every scalar count and for `M`, `I`** (§14, derived). Recommendation: read the exit criterion as the monopole sector (the candidate never contained more), schedule `l = 2` (shape, `Q`) as a separate later ADR. **Owner must confirm the reading** |
| **Q7** | Disposition of the public candidate API | **A — atomic replacement** in 4C-I (§19.1) |
| **Q8** | Expose the homogeneous (sequence-derivative) solution alongside the fixed-`ε_c` response? | **Yes, as an optional second struct** — it costs one linear solve, gives Phase 5's `B_i` from the same machinery, and provides identity K |

## 24. Documents synchronized in this task (documentation class)

| File | Change |
|---|---|
| `docs/validation/PHASE4C_HARTLE2_DERIVATION.md` | **new** — this record |
| `docs/adr/ADR-0007-hartle-second-order-monopole-response.md` | **new — PROPOSED** (Q1–Q13; §3.1 record; not accepted) |
| `docs/SCIENTIFIC_INVARIANTS.md` | INV-08 headline and status → `UNVERIFIED SCIENTIFIC CANDIDATE — GOVERNED REPLACEMENT ADR PROPOSED (ADR-0007)`; Phase-4C-G note (primary source authenticated, hypotheses classified, surface audit, derived `δN_i` formula); INV-09 note; INV-06 note (surface semantics for O(Ω²); stale `1e-5` comment); INV-14 cross-reference; summary table |
| `docs/architecture/CURRENT_ARCHITECTURE.md` | evidence-scope row; O(Ω²) row; §6 non-claims |
| `docs/MODERNIZATION_ROADMAP.md` | Phase-4 status and 4C row (4C-G done, adjudication next); ratified bullet annotated; exit-criterion reading surfaced; dependency summary |
| `docs/adr/README.md` | ADR-0007 PROPOSED index row; anticipated row resolved |

No accepted ADR decision text was altered. Nothing under `CompactStar/**`, `tests/**`,
`tests/baselines/**` or any `CMakeLists.txt` changed.

## 25. Exact non-scope

No production change; no repair of the candidate; no candidate execution or baseline; no
Phase-4D test; no rotochemical work; no `l = 2` derivation beyond the scope determination of
§14; no MixedStar work; no EOS-derivative API (surfaced as a prerequisite only); no change to
INV-06's convention; no adjudication of the solar-mass / `G` authority; no ADR acceptance.

## 26. Status

> **`PHASE-4C O(OMEGA^2) MONOPOLE ADR PROPOSED — OWNER ADJUDICATION REQUIRED`**

Primary-source access was authenticated (§2). The governed `l = 0` system, its variables,
normalization, centre condition, EOS ownership, surface semantics, exterior matching and
Phase-5 interface are derived and dimensionally audited (§3–§15). The surface convention needs
no prior governance (§11.3). The candidate is byte-identical and remains an unverified
candidate (§16). ADR-0007 is drafted PROPOSED.

**Exact next action.** The owner reviews and adjudicates ADR-0007 (Q1–Q13) **before any
O(Ω²) production change.**
