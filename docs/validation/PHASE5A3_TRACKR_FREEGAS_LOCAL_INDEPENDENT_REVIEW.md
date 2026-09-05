# Phase 5A-3 — Independent review of the Track-R free-gas local thermodynamics

> **Status (2026-09-04): PHASE 5A-3 ACCEPTABLE WITH NONBLOCKING FOLLOW-UP —
> READY FOR CANONICAL INTEGRATION.**
>
> No blocking production defect was found. The local cold ideal `npe-mu`
> thermodynamics of `6cc1bbc` is source-faithful, analytically correct, and
> independently reproduced to double round-off by calculations that share no
> code with the repository. Ten nonblocking findings are ranked in section 16;
> **R1 (no muon-free `npe` branch) must be resolved before whole-star Track-R
> reproduction**, and the phrase "WHOLE-STAR FREE-GAS REPRODUCTION READY" in
> the implementation record overstates readiness because of it.
>
> This review changed no production source, no test, and no baseline.

---

## 1. Reviewed SHA and authentication

| Item | Authenticated value |
|---|---|
| Reviewed HEAD | `6cc1bbc3f134f9c322899967e864c78fe48166fb` |
| Branch | `physics/phase5a-trackr-freegas-local` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-trackr-freegas` |
| Canonical master | `8658fee45a77c483c2fbab6b250b9be6d5b4acdf` |
| Relationship | **1 ahead / 0 behind** `master` (`git rev-list --left-right --count master...HEAD` = `0 1`) |
| Working tree | clean, including untracked files (`git status --porcelain` empty) |
| Local / upstream | `HEAD` = `origin/physics/phase5a-trackr-freegas-local` = `6cc1bbc` |
| Live remote | `git ls-remote origin`: branch `6cc1bbc`, `master` `8658fee` — both equal to local after `git fetch --all --prune` |
| Governing authority | `GOVERNANCE.md` (rank 1); accepted **ADR-0010** (rank 2); `docs/SCIENTIFIC_INVARIANTS.md` (rank 3); `docs/architecture/CURRENT_ARCHITECTURE.md` (rank 6) |
| Change class of this review | **documentation** — one new validation record; no production, test, baseline, or literature change |

`AGENTS.md` §1 preflight was performed before any reading conclusion was drawn.
`GOVERNANCE.md` §3.1 is **not** invoked by the reviewed commit and is not
invoked by this review.

## 2. Source ledger

All PDFs were read directly from the shared read-only library
`/Users/keeper/Documents/CompactStar/literature/rotochemical/`. Page numbers are
PDF pages.

| Source | Location read | What it establishes |
|---|---|---|
| `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | p. 2, §1 | "we have chosen the simplest possible core composition, namely neutrons, protons, electrons, and muons (npeµ matter)". Superfluidity ignored. |
| Same | pp. 8–9, §3.1 | Third EOS set is "a non-interacting Fermi gas (see, e.g., Shapiro & Teukolsky 1983) **for the whole star**" — i.e. no separate crust EOS for this case. "The values for the Fermi gas EOS listed in Table 1 correspond to a central density slightly below the appearance of `Σ−` hyperons." |
| Same | p. 12, §3.4, eq. (50) | "for the noninteracting Fermi gas we use `m*_i = µ_i/c²` … The sum is performed over all particle species and the integral is done over the region where each particle species exists." |
| Same | p. 16, §4.3, eqs. (71)–(74) | "The first particles to appear after muons are probably the `Σ−` and `Λ0` hyperons"; `η_nnΣp = 2µ_n − µ_Σ − µ_p` (71); `η_nΛ = µ_n − µ_Λ` (72); reactions `n+n ⇌ Σ−+p` (73), `n+Λ0 ⇌ Σ−+p` (74). Hyperons were **not** included in the calculations. |
| Same | p. 16, §4.3, eq. (75) | `Σ− → n+e+ν̄_e`, with "associated chemical imbalance `η_Σne = η_npe − η_nnΣp`" — i.e. `η_Σne = µ_Σ − µ_n − µ_e`. |
| Same | p. 26, Fig. 1 (left) | Whole-star spin-down compression profile; a "Fermi gas" curve exists at fixed central pressure `4.5e34 dyn cm^-2`, normalized by `Ω_K`. No local EOS convention. |
| Same | p. 30, Table 1 + footnotes | Free-gas row: `M_max=0.62 M_sun`, `ρ_c=1.10e15 g cm^-3`, `R=12.77 km`, `R∞=13.80 km`, `P_K=0.98 ms`. Footnote c: mass is "before appearance of `Σ−` hyperons"; footnote d: the free-gas `P_K` was **adopted** from Haensel, Salgado & Bonazzola (1995) for a **pure neutron gas**, not computed for this EOS. |
| Same | p. 33, Fig. 8 | Whole-star quasi-equilibrium `T_s∞` versus `R∞`, a "Fermi gas" curve among the EOS set. Belongs to Layers B–D. |
| `2006-Reisenegger-…-Electrostatic-Potential-Perturbations.pdf` | p. 2 / journal 569, eqs. (2), (4) | Total redshifted potential `µ_i^∞=[µ_i(r)+q_iψ(r)]e^{Φ}`; local charge neutrality `Σ n_i q_i = 0` replaces hydrostatic charge balance. |
| Same | p. 3 / journal 570, eqs. (9)–(19) | `η_npl^∞ = µ_n^∞−µ_p^∞−µ_l^∞`, with the electrostatic potential cancelling "due to the equal and opposite charges of protons and leptons" — so `η_npl` is built from **intrinsic** potentials. Eq. (11) fixes `δψ` from local neutrality; eq. (13) is singular (`B_pj=B_ej+B_µj`); the reduction eliminates the proton row/column to a symmetric invertible 3×3 `B̃`. |
| `1995-Reisenegger-Rotochemical-Heating.pdf` | pp. 14–15 / preprint §2, eq. (3) | At fixed baryon density the equilibrium composition minimizes the total energy per baryon **including the electron contribution**: `∂E/∂x = µ_p + µ_e − µ_n = 0`. |

Two consequences that this review relies on, both established from the sources
and not from the implementation record:

1. **FR2005 supplies no local numerical free-gas table.** Every free-gas number
   the paper prints (Table 1, Figs. 1 and 8) is a whole-star quantity.
2. **The Σ− threshold used in the implementation is the source's own threshold.**
   FR2005 eq. (71) with `µ_Σ(n_Σ=0)=M_Σ` and FR2005 eq. (75)'s channel condition
   `µ_Σ = µ_n + µ_e` are algebraically the same condition once `η_npe=0`, because
   `η_Σne = η_npe − η_nnΣp`. Section 11 verifies they give the identical density
   to 59 significant digits.

## 3. Exact scientific claim established by `6cc1bbc`

Against the four-way classification:

| Class | Verdict |
|---|---|
| A — correct generic ideal-Fermi formulas only | **Understates.** The commit does more than supply generic formulas. |
| **B — source-faithful implementation of the LOCAL FR2005 free-gas model** | **ESTABLISHED, with one scope qualifier** (below). |
| C — reproduction of published FR2005 local numerical data | **NOT established, and not claimed.** No such data exist in the authority set. |
| D — reproduction of the published FR2005 whole-star free-gas benchmark | **NOT established, and explicitly disclaimed.** |

**The qualifier.** What is implemented is the **active `npe-mu` branch** of the
local FR2005 free-gas model — the interval `0.45698… < n_B < 0.61736…` fm^-3.
The muon-free `npe` branch of the *same* local model, which FR2005 uses over
most of the free-gas star, is not reachable through any public entry point:
`EquilibriumStateAt` throws below muon onset, and `Evaluate` rejects `n_mu = 0`
outright (verified directly — section 16, R1). The honest statement of the
claim is therefore:

> **The local cold ideal `npe-mu` FR2005 free-gas model is faithfully and
> correctly implemented on its declared active branch.**

**Overclaim check.** The code and the metadata do not overclaim. The
implementation record's §9 correctly states "NO DIRECT LOCAL NUMERICAL SOURCE
TABLE EXISTS IN THE CURRENT AUTHORITY SET" and the record repeats
"WHOLE-STAR FREE-GAS REPRODUCTION NOT YET PERFORMED". The one overclaim is the
record's own status line and roadmap wording
(`docs/validation/PHASE5A3_TRACKR_FREEGAS_LOCAL_PROVIDER.md:3-4`;
`docs/MODERNIZATION_ROADMAP.md:583-587`): **"WHOLE-STAR FREE-GAS REPRODUCTION
READY"** / "Phase 5A-3 COMPLETE" is not supported while the provider cannot
return the `npe` equilibrium composition that most of the free-gas star
occupies. See R1.

## 4. Particle-model source fidelity

Established from FR2005 directly, then compared with the implementation:

| Source assumption | FR2005 location | Implementation | Verdict |
|---|---|---|---|
| Species: n, p, e, µ | p. 2, §1 | `neutron_`, `proton_`, `electron_`, `muon_` (`TrackRFreeGasThermodynamics.hpp:71-74`) | **FAITHFUL** |
| All species noninteracting cold Fermi gases | p. 9 §3.1 ("non-interacting Fermi gas for the whole star") reinforced by p. 12 §3.4 (`m*_i=µ_i/c²` for **every** species) | all four are `ColdRelativisticIdealFermion` | **FAITHFUL** |
| Fully relativistic, no non-relativistic nucleon approximation | implied by `m*_i=µ_i/c²` and by S&T's relativistic ideal gas | `µ=sqrt(M²+p_F²)`, exact relativistic energy integral, for all four species | **FAITHFUL** |
| Beta equilibrium `µ_n−µ_p=µ_e` (and `=µ_µ` where muons exist) | R1995 eq. (3); R2006 eq. (9) | `EquilibriumResidualAtCommonLeptonMu` roots `µ_n−µ_p−λ` with `λ=µ_e=µ_µ` | **FAITHFUL** |
| Local charge neutrality `n_p=n_e+n_µ` | R2006 eq. (4) | imposed by construction in `MakeChargeNeutralCompositionState`; never clipped | **FAITHFUL** |
| Rest masses retained in the energy | implicit in a relativistic ideal gas; Table 1 `ρ_c` is total mass-energy | `epsilon_i` includes `M_i n_i`; metadata declares "total species energy including … rest masses" | **FAITHFUL** |
| Muon appearance criterion | not stated explicitly; standard threshold `µ_e = m_µ` at `n_µ=0` | `ComputeMuonOnsetBaryonDensityFm3` | **FAITHFUL** (section 10) |
| High-density termination at `Σ−` appearance | p. 9 §3.1 and Table 1 footnote c | `ComputeSigmaMinusOnsetBaryonDensityFm3` via eq. (71) | **FAITHFUL** (section 11) |
| Role of `Σ−` | p. 16 §4.3: first particle after muons; eq. (71); hyperons **not** included | used only to locate a domain endpoint; no hyperon energy, density, or susceptibility exists in the code | **FAITHFUL** |
| Other species ignored below the endpoint | p. 16 §4.3 names `Σ−` **and** `Λ0` | only `Σ−` is checked | **FAITHFUL — and independently confirmed sufficient**: see below |

**`Λ0` is not missed.** FR2005 eq. (72) gives `η_nΛ = µ_n − µ_Λ`, so a free `Λ0`
first appears when `µ_n = M_Λ = 1115.683 MeV`. Solving the same free-gas
equilibrium independently gives that at `n_B = 0.98300978247 fm^-3`, far above
the `Σ−` endpoint `0.61735520797 fm^-3`, where `µ_n = 1072.418 MeV < M_Λ`. The
`Σ−` endpoint is therefore genuinely the first hyperon threshold of this model,
and the implementation's single-threshold treatment is complete for its domain.
`Σ0` (`µ_Σ0=µ_n`, `M=1192.6`) and `Σ+` (`µ_Σ+=µ_n−2µ_e≈822 MeV`, `M=1189.4`) are
also far above threshold.

## 5. Mass-convention assessment

FR2005 states no particle masses; it cites Shapiro & Teukolsky (1983) for the
noninteracting Fermi gas. The repository therefore supplies the convention:

```
M_n = 939.56542052 MeV   M_p = 938.27208816 MeV
M_e = 0.51099895 MeV     M_mu = 105.6583755 MeV
M_Sigma- = 1197.449 MeV  hbar c = 197.32698045930249 MeV fm
```

(read from the compiled `Zaki::Physics` authority at full precision;
`hbar c = 1/MEV_2_INV_FM`). These are modern CODATA/PDG values; FR2005-era
values differ at the `1e-5`–`1e-4` relative level.

**Quantified impact on an eventual Table-1 reproduction** (independent
recomputation of both onsets under alternative conventions):

| Convention change | `n_B,µ-onset` shift | `n_B,Σ` shift |
|---|---:|---:|
| `M_Σ` 1197.449 → 1197.3 (1983-era value) | 0 | **−9.95e-04 rel** |
| `M_Σ` within its PDG uncertainty (±0.030 MeV) | 0 | ±2.00e-04 rel |
| `M_n,M_p` → 1983-era (939.5731, 938.2796) | +8.4e-06 rel | −4.1e-05 rel |
| `M_µ` → 105.6594 | +1.65e-05 rel | −4.6e-07 rel |
| `M_n = M_p` (degenerate-nucleon toy) | +1.83e-02 rel | +6.63e-03 rel |

Sensitivity: `dn_B,Σ/dM_Σ = 4.12e-3 fm^-3 MeV^-1`.

**Assessment.** Modern constants are **source-faithful enough** for the
quantitative Table-1 reproduction. The largest realistic convention ambiguity is
the `Σ−` mass, worth `1e-3` in the endpoint density. Table 1 quotes
`M_max` to two significant figures (`0.62`) and `ρ_c` to three (`1.10`), and
footnote c makes the endpoint itself approximate ("slightly below"). A `1e-3`
endpoint shift cannot move `M_max` at that precision. The one convention that
*would* matter — collapsing `M_n` to `M_p` — is not used. **No mass-convention
blocker.** The convention is declared in metadata
(`rest_mass_convention`), satisfying ADR-0010's provenance row.

**Independent external corroboration of the whole mass/threshold convention.**
Evaluating the free-gas beta-equilibrium energy density at the derived `Σ−`
endpoint gives

```
n_B = 0.617355207967 fm^-3  ->  eps = 630.3746 MeV fm^-3  ->  rho = 1.1237e15 g cm^-3
```

FR2005 Table 1 lists the free-gas central density as `1.10e15 g cm^-3`, described
in §3.1 as "slightly below the appearance of `Σ−` hyperons" — **2.1 % below** the
independently computed threshold density, in the correct direction. This is a
source-side check that the implementation did not use and did not assume (the
record explicitly declined to reinterpret Table 1's `g cm^-3` as `n_B`). It
independently corroborates the species set, the rest-mass convention, the
threshold condition, and the `Σ−` mass simultaneously. Repeating it for an
`npe`-only free gas gives `1.1160e15 g cm^-3`, so the check confirms the model
and the endpoint but does not by itself discriminate muon inclusion.

## 6. Independent species formulas

Derived from first principles (spin degeneracy `g=2`, `T=0`), then compared
with `CompactStar/EOS/src/LocalThermodynamics.cpp:104-186`:

```
n      = (2/(2 pi hbar)^3) * (4 pi/3) p_F^3 = p_F^3 / (3 pi^2 hbar^3)
  =>  p_F = hbar c (3 pi^2 n)^(1/3)                                  [MeV]
mu     = sqrt(M^2 + p_F^2)                                            [MeV]
eps    = (1/(pi^2 (hbar c)^3)) Int_0^{p_F} k^2 sqrt(M^2+k^2) dk
       = M^4 [ x sqrt(1+x^2)(1+2x^2) - asinh(x) ] / (8 pi^2 (hbar c)^3),  x=p_F/M
dmu/dn = (dmu/dp_F)(dp_F/dn) = (p_F/mu)(p_F/(3n)) = p_F^2/(3 mu n)   [MeV fm^3]
n(mu)  = [ sqrt((mu-M)(mu+M)) / (hbar c) ]^3 / (3 pi^2)
```

Every one of these is exactly what the production code computes, including the
closed form of the energy bracket (the `1/8` and the `asinh`), the `hbar c`
placement, and the `hypot`/`sqrt((mu-M)(mu+M))` cancellation-safe forms.

| Property | Production behavior | Verdict |
|---|---|---|
| Spin degeneracy | `2` (through `p_F^3/(3 pi^2 hbar^3)`) | correct |
| `hbar c` units | `1/MEV_2_INV_FM` = 197.32698045930249 MeV fm; `p_F` and `mu` in MeV, `eps` in MeV fm^-3, `dmu/dn` in MeV fm^3 | correct and dimensionally consistent |
| Rest-energy convention | included in `mu` and in `eps` (`eps -> M n` as `n -> 0`) | correct, and declared |
| Zero density | `p_F=0`, `mu=M`, `eps=0`, `dmu/dn = std::nullopt`; scalar accessor throws | correct: `dmu/dn ~ n^{-1/3}` genuinely diverges |
| Small-`x` branch (`x<1e-2`) | series `x^3(8/3 + x^2(4/5 + x^2(-1/7 + x^2(1/18 - 5x^2/176))))` | **verified exactly equal to the Taylor expansion through `O(x^11)`.** Independent 60-digit evaluation at the switch point `x=1e-2` gives relative deviation `6.31e-23` from the closed form, scaling as `x^13` — pure truncation, no coefficient error |
| Relativistic limits | `mu -> p_F`, `dmu/dn -> p_F/(3n)` as `x -> inf`; `mu -> M + p_F^2/2M` as `x -> 0` | correct |
| `NumberDensityForChemicalPotentialFm3` | throws for `mu < M`; returns exactly `0.0` at `mu = M`; exact inverse otherwise | correct; fails closed, no clamp |

**Generalization `ColdRelativisticFreeLepton -> ColdRelativisticIdealFermion`.**
The diff (`git diff 8658fee 6cc1bbc -- CompactStar/EOS/src/LocalThermodynamics.cpp`)
shows the numeric body of `Evaluate` is **byte-for-byte unchanged**; only the
enclosing class/enum/struct names and the error strings changed, and
`ColdRelativisticFreeLepton::Evaluate` now delegates to
`ColdRelativisticIdealFermion::Electron()/Muon()` and copies the six fields.
The two lepton factories bind identical constants before and after. Electron and
muon behavior is therefore preserved **by construction**, and RFG1 confirms
bit-identity at `n = {0, 1e-18, 1e-6, 0.01} fm^-3` for both leptons — including
the zero-density `std::nullopt`.

Two caveats on that evidence, neither a defect (R9): the bit-identity assertion
is tautological because the legacy path now *is* the new path, and it does not
compare `number_density_fm3` or `rest_mass_energy_MeV`. The genuinely
load-bearing lepton evidence is RFG1's comparison against the independent
32,768-panel Simpson oracle, plus this review's independent 60-digit
reproduction.

**Platform note (R8).** `sizeof(long double) == 8` and `__LDBL_DIG__ == 15` on
this arm64 macOS toolchain (measured). The `long double` intermediates carry **no
extra precision on the validation platform**; the small-`x` series is therefore
genuinely load-bearing rather than belt-and-braces, and all reported `~1e-15`
agreements are pure double round-off.

## 7. Independent derivation of the neutral conjugates `g`

Starting from the charge-neutral potential with `n_p = n_e+n_mu` and
`n_n = n_B-n_e-n_mu` substituted **before** differentiation:

```
eps_CN(n_B,n_e,n_mu) = eps_n(n_B-n_e-n_mu) + eps_p(n_e+n_mu)
                     + eps_e(n_e) + eps_mu(n_mu)
```

Using `d eps_i/d n_i = mu_i`:

```
d eps_CN/d n_B  = mu_n * (+1)                              =  mu_n
d eps_CN/d n_e  = mu_n*(-1) + mu_p*(+1) + mu_e*(+1)        = -(mu_n-mu_p-mu_e) = -eta_npe
d eps_CN/d n_mu = mu_n*(-1) + mu_p*(+1) + mu_mu*(+1)       = -(mu_n-mu_p-mu_mu)= -eta_npmu
```

so `d eps_CN = mu_n d n_B - eta_npe d n_e - eta_npmu d n_mu`, exactly as ADR-0010
requires, and `g = (mu_n, -eta_npe, -eta_npmu)`.

**Production output matches exactly.** `TrackRFreeGasThermodynamics.cpp:249-254`
forms `eta_npe = mu_n-mu_p-mu_e`, `eta_npmu = mu_n-mu_p-mu_mu`, and returns
`{mu_n, -eta_npe, -eta_npmu}`; `NeutralConjugates::EtaNpeMeV()` returns
`-value_MeV[1]`, restoring `eta_npe` with the accepted sign
(`LocalThermodynamics.hpp:85-92`). Signs, units (`MeV`), and index order all
agree with ADR-0010 and with R2006 eq. (9) (intrinsic potentials only, no
electrostatic or redshift term).

**Independent numerical confirmation** at the off-equilibrium state
`x = (0.52, 0.030, 0.011) fm^-3`, computed in 50-digit decimal arithmetic sharing
no code with the repository:

| Quantity | Provider | Independent | Rel. dev. |
|---|---|---|---:|
| `eps` | `528.42025446273863` | `528.42025446273874` | `2.9e-16` |
| `g_0 = mu_n` | `1054.0056932138762` | `1054.0056932138762` | `6.4e-18` |
| `g_1 = -eta_npe` | `97.275574846119724` | `97.27557484611971` | `1.6e-16` |
| `g_2 = -eta_npmu` | `79.624664008286175` | `79.62466400828616` | `1.8e-16` |

(`eta_npe = -97.2756 MeV`, `eta_npmu = -79.6247 MeV` — the fixture is genuinely
far off equilibrium, so `Evaluate` demonstrably does not impose beta
equilibrium.)

Additionally, the optional `IntrinsicChemicalPotentialsMeV` was checked against
ADR-0010's `S_x^T mu = g` requirement at four states spanning
`n_B = 0.30 … 0.615 fm^-3`: **worst relative deviation `0.000e+00`** (bit-exact).
This API has no committed test (R2).

## 8. Independent derivation of the Hessian `H`

Differentiating `g` again, with `D_i = dmu_i/dn_i`:

```
d g_0/d n_B = D_n            d g_0/d n_e  = -D_n           d g_0/d n_mu = -D_n
d g_1/d n_B = -D_n           d g_1/d n_e  = D_n+D_p+D_e    d g_1/d n_mu = D_n+D_p
d g_2/d n_B = -D_n           d g_2/d n_e  = D_n+D_p        d g_2/d n_mu = D_n+D_p+D_mu
```

which is exactly

```
H = [ D_n,   -D_n,           -D_n          ]
    [ -D_n,  D_n+D_p+D_e,    D_n+D_p       ]
    [ -D_n,  D_n+D_p,        D_n+D_p+D_mu  ]     [MeV fm^3]
```

matching `TrackRFreeGasThermodynamics.cpp:255-259` entry for entry, sign for
sign. This derivation was performed from the potential, **not** read off the
checked-in matrix.

| Property | Independent verdict |
|---|---|
| Every sign | correct (three `-D_n` in row 0, three in column 0, `+` elsewhere) |
| Units | `MeV fm^3` — `d(MeV)/d(fm^-3)`; consistent with ADR-0010 |
| Symmetry | exact: `H_01=H_10=H_02=H_20=-D_n`, `H_12=H_21=D_n+D_p`. Measured asymmetry `0` on 36 states (RFG7) |
| `H = S_x^T K S_x` | verified term by term with `S_x=[[1,-1,-1],[0,1,1],[0,1,0],[0,0,1]]`, `K=diag(D_n,D_p,D_e,D_mu)`: `(S_x^T K S_x)_{ab} = sum_i S_{ia} D_i S_{ib}` reproduces all nine entries |
| Positive definiteness | `v^T H v = sum_i D_i (S_x v)_i^2`. `S_x` has **full column rank 3** (`a c_0 + b c_1 + c c_2 = 0` forces `b=c=0` from the e- and µ-rows, then `a=0`). With `D_i>0` for `n_i>0`, `H` is PD **everywhere on the declared `Evaluate` domain**, not merely on the tested grid |

**Independent numerical confirmation** at `x = (0.52, 0.030, 0.011)`:
`D_n,D_p,D_e,D_mu = 150.629797191, 374.638252893, 2107.60782932, 3246.76519876
MeV fm^3`; the nine provider entries agree with the independently assembled
matrix to a **maximum relative deviation of `2.6e-16`**.

## 9. Equilibrium solver assessment

The beta-equilibrium problem at fixed `n_B` is: find `(n_e,n_mu)` with
`mu_n-mu_p=mu_e` and `mu_n-mu_p=mu_mu`, subject to `n_p=n_e+n_mu`,
`n_n=n_B-n_p`. Because the two conditions force `mu_e = mu_mu ≡ λ`, and because
each free lepton density is a strictly monotone function of its own chemical
potential, the two-dimensional problem collapses **exactly** to one scalar
unknown `λ`. **The chosen root variable is sufficient**; nothing is lost.

```
n_e(λ)  = [sqrt(λ²-M_e²)/hbar c]³/(3 pi²)
n_mu(λ) = [sqrt(λ²-M_mu²)/hbar c]³/(3 pi²)   (0 for λ <= M_mu)
F(λ)    = mu_n(n_B-n_e-n_mu) - mu_p(n_e+n_mu) - λ
```

| Question | Independent finding |
|---|---|
| Is the bracket valid? | **Yes.** `λ_low = M_mu`: then `n_mu=0` and `n_p=n_e(M_mu)`, and `F(M_mu)` is exactly the muon-onset residual — zero at onset and **strictly positive above it**, which is precisely the interval `SolveActiveEquilibriumUnchecked` requires. `λ_high = mu_n(n_B)-M_p`: since `n_n<n_B` and `n_p>0`, `F(λ_high) < mu_n(n_B)-M_p-λ_high = 0` strictly. Both are established analytically, not assumed. |
| Is the residual monotone? | **Yes, strictly decreasing.** `dn_e/dλ>0` and `dn_mu/dλ>=0`, so `n_p` increases, `n_n` decreases, `mu_n` decreases, `mu_p` increases, and `-λ` decreases: every term of `dF/dλ` is `<0`. |
| Is the root unique? | **Yes**, by strict monotonicity on the bracket. |
| Does bisection fail closed if unbracketed? | **Yes.** `BisectDecreasingRoot` throws unless `F(lo)>=0` and `F(hi)<=0` and both are finite; `EquilibriumResidualAtCommonLeptonMu` throws if `n_B <= n_p`; `ChargeDensityAtCommonLeptonMuFm3` throws below `M_e`; non-finite midpoints throw; 256 iterations without convergence throws. |
| Stopping criteria reasonable? | **Yes.** Residual `<= 2e-12 MeV` on a `~125 MeV` scale (`1.6e-14` relative), or a bracket width of `64 eps * max(1,|lo|,|hi|) ~ 1.8e-12 MeV`, reached in ~44 bisections — far inside the 256 cap. Bisection on a strictly monotone residual cannot converge falsely. |
| Any arbitrary density floor or unbracketed solve? | **None found.** No clamp, no epsilon floor, no fallback root, no secant/Newton escape hatch. |
| Does `Evaluate` ever impose equilibrium? | **No.** `Evaluate` never calls the solver; verified numerically (`eta_npe = -97.28 MeV` at the probe state). |

**Independent verification of the solved compositions.** Six equilibrium states
were solved with a completely separate 50-digit implementation and compared with
the provider:

| `n_B` | max rel. deviation over `(n_n, n_e, n_mu)` |
|---:|---:|
| 0.46 | `3.7e-12` (in `n_mu`) |
| 0.48 | `7.8e-14` |
| 0.50 | `3.9e-13` |
| 0.55 | `3.1e-14` |
| 0.60 | `7.0e-14` |
| 0.617 | `9.3e-14` |

`n_n` agrees to `<= 6.3e-16` at every point. The larger `n_mu` deviations are
the expected threshold amplification: `n_mu ∝ (λ-M_mu)^{3/2}`, so a `2e-12 MeV`
residual tolerance yields `~1.5 δλ/(λ-M_mu)` relative error in `n_mu`, which is
`~7e-12` at `n_B=0.46` (0.4 MeV above threshold) and shrinks rapidly. This is
solver tolerance, not a formula error.

**Assessment: the equilibrium solver is VALID.**

## 10. Independently calculated muon onset

Derivation. At the threshold the muon Fermi sea is empty, so `mu_mu = M_mu`;
beta equilibrium forces `mu_e = mu_n - mu_p = M_mu`; and `n_mu=0` gives
`n_p = n_e`. Hence the closed construction

```
n_e     = n(mu = M_mu; M_e)                    = 5.18461054657707e-03 fm^-3
mu_p    = mu(n_e; M_p)                         = 944.202278441171 MeV
mu_n    = mu_p + M_mu                          = 1049.86065394117 MeV
n_n     = n(mu = mu_n; M_n)                    = 0.451800194865842 fm^-3
n_B,mu  = n_n + n_e
```

**Independent value (60-digit arithmetic):**

```
n_B,mu-onset = 0.4569848054124190326  fm^-3
```

**Repository value:** `0.45698480541241987 fm^-3`
**Agreement:** absolute `8.4e-16`, **relative `1.8e-15`** — pure double round-off
(`~15 ulp` accumulated through `cbrt`, `sqrt`, and the reconstruction, with no
`long double` headroom on arm64). **AGREES.**

**Behavioral correctness at the threshold** (probed directly against the built
library):

| Behavior | Observed | Verdict |
|---|---|---|
| `npe` branch below onset (`n_mu=0`, `mu_e<M_mu`) | RFG10 independently solves it and measures `mu_e = 105.0358 MeV < M_mu`. It exists physically. | correct physics |
| Threshold at `n_mu=0` | `EquilibriumStateAt(onset)` throws | correct: `D_mu` is undefined there |
| Active branch above onset | `EquilibriumStateAt(onset*(1+1e-6))` returns `n_mu>0` | correct |
| Singular `dmu_mu/dn` at zero | `Evaluate({0.9*onset, 0.04, 0.0})` throws; `ColdRelativisticIdealFermion::Muon().Evaluate(0)` returns `nullopt` | correct — `D_mu ~ n^{-1/3} -> inf` |
| No smoothing | `Evaluate({0.9*onset,0.04,1e-30})` succeeds and `H(2,2)-H(2,1)` equals the exact `D_mu(1e-30)` to `2e-14` | **no floor, no regularization** |
| Off-equilibrium `Evaluate` below onset | `Evaluate({0.30, 0.010, 0.004})` succeeds (verified directly) | correct: `Evaluate` has its own, wider domain |

**Is failing closed at/below onset consistent with the accepted contract?**
For `EquilibriumStateAt` **as a paired input to a fixed-dimension 3×3 Hessian**,
yes: ADR-0010 requires that at "absent-species boundaries … report the domain
restriction rather than fabricate a full-rank Hessian", and throwing is a
report, not a silent change of coordinate meaning. **But it is stricter than the
ADR requires.** `ChargeNeutralCompositionState` can represent `n_mu = 0`
exactly, `MakeChargeNeutralCompositionState` accepts it, and the Hessian is a
*separate* call — so returning the `npe` equilibrium composition would violate
nothing. ADR-0010's REQUIRED-NOW row asks the provider to return "the model's
equilibrium composition **and active-species/threshold status**"; the current
`ChargeNeutralCompositionState` return type carries no status channel at all, so
the only available report is an exception. This is exactly the "phase/threshold/
failure status" item the Phase 5A-2 review classified as **FUTURE HARDENING now;
MUST FIX BEFORE TRACK-R crosses such domains**
(`docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_INDEPENDENT_REVIEW.md:404`).
Track-R now *does* need to cross that domain — see R1. This is **not** a blocker
for the local increment as scoped, and this review does not redesign it.

## 11. Independently assessed `Σ−` endpoint

**What FR2005 actually says.** §3.1 (p. 9) and Table 1 footnote c (p. 30) say the
free-gas entry is taken at a central density *slightly below* `Σ−` appearance.
§4.3 (p. 16) supplies the imbalance `η_nnΣp = 2µ_n − µ_Σ − µ_p` (eq. 71) for the
weak reaction `n+n ⇌ Σ−+p` (eq. 73), and the `Σ−` direct-Urca channel
`Σ− → n+e+ν̄_e` (eq. 75) with `η_Σne = η_npe − η_nnΣp`. Hyperons are not in the
calculations; the paper gives no `Σ−` mass and no threshold density.

**The physical threshold for a free gas.** First appearance means `n_Σ = 0`, so
for an ideal `Σ−` gas `µ_Σ(0) = M_Σ`. Setting `η_nnΣp = 0` gives
`2µ_n − µ_p = M_Σ`. Along the already beta-equilibrated `npe-mu` branch
(`η_npe = 0`) this is *identical* to the channel condition
`µ_Σ = µ_n + µ_e`, because `η_Σne = η_npe − η_nnΣp`.

**Independent verification of that equivalence.** Solving the free-gas
equilibrium and rooting the two residuals separately in 60-digit arithmetic:

```
eq. (71)  2 mu_n - mu_p = M_Sigma   ->  n_B = 0.61735520796652986888639657543168059625650749407468821577880
eq. (75)  mu_n + mu_e   = M_Sigma   ->  n_B = 0.61735520796652986888639657543168059625650749407468821577880
difference                              =  0E-59
```

They are the same density to every digit carried. The implementation's use of
eq. (71) is therefore **the source's own condition**, not a substituted
convention.

**Independent endpoint value:** `n_B,Σ = 0.6173552079665298689 fm^-3`
**Repository value:** `0.61735520796653498 fm^-3`
**Agreement:** absolute `5.1e-15`, **relative `8.3e-15`** — well inside the
`2e-11 MeV` residual tolerance of the outer bisection (which corresponds to
`~1.8e-13 fm^-3`). The committed test's own second solve reports
`0.61735520796652943`, `4.4e-16` from the independent value. **AGREES.**

**Itemized answers to the task's checklist:**

| Item | Finding |
|---|---|
| `Σ−` rest mass used | `Zaki::Physics::SIGMA_MINUS_M_MEV = 1197.449 MeV` (PDG). Not supplied by FR2005 — a repository convention, declared in the source comment but **not** in the provider metadata (R10-adjacent) |
| Chemical-equilibrium condition used | FR2005 eq. (71) with `µ_Σ(n_Σ=0)=M_Σ`, evaluated on the beta-equilibrated `npe-mu` branch |
| Is `p_F,Σ = 0` at onset? | **Yes** — that is exactly the content of `µ_Σ = M_Σ`. The code introduces no `Σ−` density, energy, or derivative |
| Equivalent to `µ_Σ = µ_n + µ_e`? | **Yes, exactly**, on the equilibrium branch — verified to 59 digits above |
| Does the paper determine the endpoint uniquely? | **Not by itself.** FR2005 fixes the *condition* (eq. 71 / eq. 75) and the *species* (`Σ−`), but supplies neither `M_Σ` nor a threshold density. The endpoint is uniquely determined only once the mass convention is added |
| Convention introduced beyond the source? | **One: the numerical `M_Σ`.** It is a standard PDG value, it is the same authority used for every other mass, and its `±0.030 MeV` uncertainty moves the endpoint by `2e-4` relative. **Flagged, not a defect.** |
| Effect on future Table-1 reproduction | Bounded and small. Using the 1983-era `M_Σ=1197.3` shifts the endpoint by `-1.0e-3` relative; Table 1 reports `M_max` to two significant figures. The independent `ρ(n_B,Σ) = 1.1237e15 g cm^-3` versus Table 1's `1.10e15` (2.1 % below, correct direction) is strong evidence the convention reproduces the paper's own endpoint |

**No hyperon physics is implemented**, as required.

## 12. RFG1–RFG11 coverage assessment

Every line was re-read and re-run. The focused executable reproduces every
number in the implementation record exactly (section 17).

| ID | Intended falsifier | Oracle | Independent? | Could a plausible wrong implementation pass? |
|---|---|---|---|---|
| **RFG1** | wrong mass/`hbar c`/`p_F`/`mu`/`eps`/`D` per species; refactor regression | 32,768-panel long-double Simpson of `∫t²√(1+t²)dt`, plus separate `p_F`, `mu`, reciprocal-`D` formulas, at `x={1e-3,0.1,1,10}` for **all four** species | **Partly.** The energy quadrature is a genuinely different algorithm from the closed form/series. But `p_F(n)`, `mu(p_F)` and `D` are formula-identical to production | **Yes, for one class:** a wrong spin degeneracy in `n ↔ p_F` would be invisible (R3). Everything else is covered. Legacy-lepton bit-identity is tautological by delegation and omits two fields (R9) |
| **RFG2** | wrong neutral energy assembly | independent four-species Simpson sum at three off-equilibrium states; max rel. err `3.86e-16` | **Yes for the sum and the species reconstruction** | A *common* prefactor error in `eps` would cancel here — but **not** in RFG4, which pins `d eps/d n = mu` against an independently computed `mu`. Net: **no** |
| **RFG3** | wrong `eta` signs, wrong index order, accessor sign flip | independently computed `mu_n,mu_p,mu_e,mu_mu`, then `eta = mu_n-mu_p-mu_l`; **public accessors called directly**; fixture asserted `|eta|>1 MeV` | **Yes** | **No.** Both nonzero `eta` signs are independently fixed |
| **RFG4** | wrong `g`, wrong held-fixed variables, wrong energy normalization | centered differences of the **independent** energy quadrature along all three canonical axes; four step halvings, monotone convergence required; final errors `3.1e-7 / 4.7e-7 / 1.8e-7 MeV` | **Yes** | **No.** All three canonical energy derivatives are independently pinned. This is the load-bearing thermodynamic-consistency check |
| **RFG5** | wrong `H` entries | matrix assembled from independently computed `D_n,D_p,D_e,D_mu` — **but with the 3×3 structure hard-coded identically to production**; max abs err `0` | **Partly** — the `D_i` are independent, the *structure* is not | **Yes in principle** for a structural sign error — but RFG11 and RFG6 close that gap. Documentation nit only (R4) |
| **RFG6** | wrong `H` column, wrong held-fixed coordinate | centered differences of the provider's `g` along each column; four halvings; per-column final errors `8.1e-7 / 2.5e-5 / 2.8e-5 MeV fm^3` | Structurally independent of the hard-coded matrix | **No.** Every column is perturbatively checked, and chained with RFG4 (`g` independently pinned) this validates `H` end to end |
| **RFG7** | asymmetry, indefiniteness, ill conditioning | Sylvester minors + `κ_∞` on a 36-state grid; asymmetry `0`; min minors `140.14 / 237340.6 / 5.553e8`; `κ_∞ <= 51.44` | Self-computed from provider output | **Coverage gap (R5):** the grid spans only `n_B ∈ [onset+0.08Δ, onset+0.88Δ]` and `n_mu/n_B ∈ [0.008,0.04]` — it does **not** cover the sub-onset part of the declared `Evaluate` domain, nor `n_mu → 0⁺` where `κ` diverges. PD itself is analytically guaranteed everywhere (section 8), so this is a documentation/coverage issue, not a soundness one |
| **RFG8** | wrong equilibrium root; only one beta channel satisfied; equilibrium not an energy minimum | exact `n_p == n_e+n_mu`; **both** `|eta_npe|` and `|eta_npmu|` `<= 1.15e-12 MeV` at four densities; symmetric energy rise in three composition directions, min rise `1.39e-7 MeV fm^-3` | Chains through RFG3-validated `eta` | **No.** Both beta channels **are** tested independently, and the minimum property adds a second, distinct falsifier |
| **RFG9** | equilibrium-solver test not load-bearing | **none at runtime** — as committed it is an unconditional `std::cout` line. The mutation (`mu_n-mu_p-1.01λ`) and its SHA-256-attested revert are recorded in the implementation record | Not a runtime check | The *stated* claim is honest ("exercised separately and reverted"). Independent runtime falsification of the solver is nevertheless supplied by RFG8 and by RFG10's fully independent Σ-onset solve (R6) |
| **RFG10** | smoothing, muon floor, wrong onset, wrong domain, missing provenance | independent closed-form muon onset (`2e-14`); independent `npe` sub-threshold solve (`mu_e=105.0358<M_mu`); independent Σ-onset solve using a **test-local** equilibrium solver (`2e-11`); six domain probes; `n_mu=1e-30` unfloored via `H(2,2)-H(2,1)`; eight metadata assertions | **Yes** — the strongest single test in the file | **No.** Genuinely proves no smoothing, no floor, and correct fail-closed domain behavior at both ends |
| **RFG11** | fabricated Hessian structure; hidden projector | explicit `S_x` from ADR-0010 contracted with independently computed diagonal `K`, at three **physical density-dependent** states; max err `0` | **Yes** — structurally independent of the hard-coded matrix | **No.** This is the check that makes RFG5's non-independence harmless |

**Material missing falsifiers** (both nonblocking, both verified correct by this
review): **R2** — `IntrinsicChemicalPotentialsMeV` is a public production API with
**zero** callers in the whole repository and no test, yet ADR-0010 imposes
`S_x^T mu = g` on it. **R3** — the `n ↔ p_F` spin-degeneracy normalization is not
independently falsified anywhere in RFG1–RFG11.

## 13. Available whole-star FR2005 source targets for the next stage

Authentication of the boundary claim: **NO DIRECT LOCAL NUMERICAL SOURCE TABLE
EXISTS IN THE CURRENT AUTHORITY SET — CONFIRMED.** Pages 2, 8–9, 12, 16, 26, 30
and 33 of FR2005 were read in full; every free-gas number is a whole-star
quantity. Planning inventory only — **nothing below was computed in this review
beyond the two consistency checks explicitly marked as such.**

| # | Quantity | Source location | Published value | Type | Convention / ambiguity |
|---|---|---|---|---|---|
| T1 | `M_max`, free gas | Table 1, p. 30 | `0.62 M_sun` | exact tabular, 2 s.f. | Footnote c: **maximum mass before `Σ−` appearance** — a truncated sequence, not the true OV maximum of the EOS |
| T2 | `ρ_c`, free gas | Table 1, p. 30 | `1.10e15 g cm^-3` | exact tabular, 3 s.f. | Total mass-energy density; §3.1 "slightly below the appearance of `Σ−`". **Independent consistency check performed in this review:** the free-gas `Σ−` threshold sits at `1.1237e15 g cm^-3`, so `1.10e15` is 2.1 % below it — consistent, correct direction |
| T3 | `R`, free gas | Table 1, p. 30 | `12.77 km` | exact tabular, 4 s.f. | Schwarzschild coordinate radius of the same configuration |
| T4 | `R∞`, free gas | Table 1, p. 30 | `13.80 km` | exact tabular, 4 s.f. | `R∞ = R/sqrt(1-2GM/Rc²)`; consistency with T1/T3 is itself a check |
| T5 | `P_K`, free gas | Table 1, p. 30 | `0.98 ms` | exact tabular, but **adopted** | Footnote d + §3.1: taken from Haensel, Salgado & Bonazzola (1995) **for a pure neutron gas**, *not* computed for this `npe-mu` EOS and *not* from the Lasota formula. **Do not treat as an output of the reproduction.** |
| T6 | Spin-down compression profile, free gas | Fig. 1 left, p. 26 | curve `~2.0–2.35` over `r/r_core ∈ [0,1]` | **figure-derived only** | Fixed central pressure `P_0 = 4.5e34 dyn cm^-2` (not fixed mass); ordinate normalized by `Ω_K`; digitization uncertainty |
| T7 | Quasi-equilibrium `T_s∞` for PSR J0437-4715, free gas | Fig. 8, p. 33 | `log T_s∞ ≈ 4.67–4.70` near `R∞ ≈ 13.2–13.8 km` | **figure-derived only** | Requires Layers B–D (rates, heat capacity, envelope, evolution) and the `M_PSR = 1.58 ± 0.18 M_sun` bold-line band; needs the 2006 correction where it governs |
| T8 | `m*_i = µ_i/c²` for the free gas | §3.4, eq. (50), p. 12 | — | prescription, not a number | Required input for the free-gas heat capacity `C̃`; the species sum runs only over regions where each species exists — i.e. it **needs the muon-free region too** |
| T9 | Crust prescription for the free gas | §3.1, p. 9 | — | prescription | The Fermi gas is used **for the whole star**; the PRL95/HP94 crusts apply to the APR and BPAL sets only. A crustless free-gas star is the correct reproduction target |
| T10 | Global scalings `L∞_γ,eq`, `T_s,eq` | eqs. (67)–(69), p. 15 | `L ≈ 10^{30-31}(...)^{8/7}`, `T_s,eq ≈ (2–3)e5 K` | ranges over **all** EOSs | Not free-gas-specific; usable only as an order-of-magnitude sanity band |

Ordering note for Track-R: ADR-0010 Q1 makes the free-gas case a **precursor**,
not the closing benchmark; closure is on A18 + δυ + UIX* with a
correction-sensitive comparison.

## 14. Architecture / API assessment

| Check | Finding | Class |
|---|---|---|
| Conforms to ADR-0010 | **Yes.** Canonical chart `x=(n_B,n_e,n_mu)`; exact `n_p=n_e+n_mu`, `n_n=n_B-n_p`; `g=(mu_n,-eta_npe,-eta_npmu)` in MeV; `H_ab=∂g_a/∂x_b` in `MeV fm^3` at fixed remaining coordinates; energy density supplied for validation; `EquilibriumStateAt` present; provenance metadata present | **NO ISSUE** |
| `Evaluate()` genuinely off equilibrium | **Yes.** It never calls the solver; measured `eta_npe = -97.28 MeV`, `eta_npmu = -79.62 MeV` at the probe state | **NO ISSUE** |
| Equilibrium recovery separate | **Yes** — distinct method, distinct (narrower) domain | **NO ISSUE** |
| Metadata honesty | Model id, revision, particle content, chart, `T=0` scope, rest-mass convention, lepton ownership ("exactly once"), and a **numeric** smooth-domain string carrying both onsets to 17 digits; plus typed `MuonOnsetBaryonDensityFm3()` / `SigmaMinusOnsetBaryonDensityFm3()`. This closes the Phase 5A-2 "provenance for realistic Track R" item. **Not** carried: derivative orientation and explicit unit fields (in the type/comments only), and the `Σ−` mass convention | **FUTURE HARDENING** (generic-metadata item, consistent with the 5A-2 deferral) |
| No stellar integration in the local model | **Confirmed.** The only includes are `LocalThermodynamics.hpp`, `<array>`, five std headers, and `Zaki/Physics/Constants.hpp` | **NO ISSUE** |
| No electrostatic projection | **Confirmed** by source inspection; the strings `psi`, `electrostatic`, `project` appear only inside negative-statement comments | **NO ISSUE** |
| No paper-`B`/`Z`/`W` object, no inverse | **Confirmed.** No matrix inverse exists in production; `Inverse3` lives only in the test | **NO ISSUE** |
| No coupling to TOV/rotation/evolution | **Confirmed.** `TrackRFreeGasThermodynamicProvider` has zero callers outside its own translation unit; registered only as an EOS library source | **NO ISSUE** |
| `ColdRelativisticFreeLepton -> ColdRelativisticIdealFermion` generalization | Numeric body unchanged; leptons delegate. **Latent hazard:** the legacy wrapper now ignores its own `rest_mass_energy_MeV_`/`hbar_c_MeV_fm_` members and rebuilds from `kind_`. Unreachable today (private ctor, two factories with identical constants), but a future third factory would silently diverge from its own accessors | **FUTURE HARDENING** (R9) |
| Muon-free `npe` branch unreachable | Neither `EquilibriumStateAt` nor `Evaluate` will produce or accept `n_mu = 0` | **MUST FIX BEFORE WHOLE-STAR REPRODUCTION** (R1) |
| `IntrinsicChemicalPotentialsMeV` untested | Public API, zero callers, ADR-0010 imposes `S_x^T mu = g` on it | **FUTURE HARDENING** (R2) |
| `Σ−` onset bracket robustness | The outer bracket starts at `nextafter(muon_onset, +inf)`, where the **inner** solve's bracket is valid by `~1 ulp`; `BisectDecreasingRoot` throws on `lower_value < 0`. A 1-ulp rounding change in `(n_n0+n_e0)-n_e0` could make the **constructor** throw and the provider unconstructible. Also, bracketing evaluates `SolveActiveEquilibriumUnchecked` up to `2 × muon_onset ≈ 0.914 fm^-3`, well past the source ceiling (internal only) | **FUTURE HARDENING** (R7) |
| `CURRENT_ARCHITECTURE.md` component table | A header banner was added, but the "EOS — Phase 5A local thermodynamics" table (`:178-181`) gained **no row** for `ColdRelativisticIdealFermion` or `TrackRFreeGasThermodynamicProvider`, and no `GOVERNANCE.md` §5 candidate marking. The implementation record does carry the candidate statement | **FUTURE HARDENING** (R10) |

**No BLOCKING LOCAL-PROVIDER DEFECT was identified.**

## 15. Scope-integrity result

`git diff --name-only 8658fee 6cc1bbc` touches exactly eleven files, all under
`CompactStar/EOS/`, `docs/`, or `tests/`. Nothing outside those trees changed.

| Prohibited action | Result |
|---|---|
| Alters `TOVSolver` behavior | **NO** |
| Alters `StarProfile` | **NO** |
| Alters `RotationSolver` | **NO** |
| Activates legacy `Rotochemical` / `RotochemicalCache` / `ChemState` code | **NO** — `Physics/Evolution/CMakeLists.txt` and `Driver/Chem/CMakeLists.txt` are untouched |
| Implements whole-star free gas | **NO** |
| Implements APR / BPAL | **NO** |
| Implements DS(CMF) off-equilibrium physics | **NO** |
| Constructs stellar `Btilde` | **NO** |
| Constructs `Z` / `W` | **NO** |
| Implements reaction rates or evolution | **NO** |
| Begins BNV | **NO** |
| Integrates stellar susceptibilities | **NO** |

**Eight governed baseline artifacts:** `git diff 8658fee 6cc1bbc -- tests/baselines/`
is **empty**, and the live `shasum -a 256` of all eight files reproduces the
implementation record's table byte for byte:

```
7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b  baryon_number_dscmf1_reference.tsv
2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb  grid_convergence_cmf_1p6_debug.tsv
7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f  grid_convergence_cmf_1p6_trajectory.tsv
a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9  hartle_I_dscmf1_debug.tsv
bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52  hartle_monopole_dscmf1_debug.tsv
afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de  passive_cooling_cmf_1p6_debug.tsv
3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a  tov_dscmf1_reference.tsv
5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678  tov_path_equivalence_dscmf1.tsv
```

**SCOPE INTEGRITY: CLEAN.**

## 16. Ranked findings

| # | Finding | Class | Smallest required response |
|---|---|---|---|
| **R1** | **The muon-free `npe` equilibrium branch is unreachable.** `EquilibriumStateAt` throws for `n_B <= 0.45698 fm^-3` and `Evaluate` rejects `n_mu = 0`, so the provider can supply neither composition nor thermodynamics for muon-free matter. For the FR2005 free-gas star this is most of the star: muon onset is `n_B = 0.457 fm^-3` (`ρ = 8.20e14 g cm^-3`) while Table 1's free-gas central density is `ρ_c = 1.10e15 g cm^-3` — muons occupy only the inner core. FR2005 §3.4 also requires the species sum "over the region where each particle species exists", i.e. the muon-free region explicitly. Consequently the record's **"WHOLE-STAR FREE-GAS REPRODUCTION READY"** and the roadmap's "Phase 5A-3 COMPLETE" overstate readiness | **MUST FIX BEFORE WHOLE-STAR REPRODUCTION** (not a local-provider defect; the declared domain is honest and ADR-0010 permits declaring one) | Governed follow-up that (a) exposes the `npe` equilibrium branch and its thermodynamics, and (b) supplies the per-result active-species/threshold status the Phase 5A-2 review already deferred (`PHASE5A2_…_INDEPENDENT_REVIEW.md:404`), **without** fabricating a full-rank 3×3 Hessian at `n_mu=0`. Qualify the status wording in the same change |
| **R2** | `IntrinsicChemicalPotentialsMeV` is a public production API with **zero callers anywhere in the repository** and no RFG coverage, yet ADR-0010 requires `S_x^T mu = g` of any supplied intrinsic potentials. *This review verified it holds bit-exactly (worst deviation `0.000e+00`) at four states spanning `n_B = 0.30–0.615 fm^-3`* | **FUTURE HARDENING** — material missing falsifier | Add one `S_x^T mu == g` assertion (plus its domain-failure probes) to the RFG ladder |
| **R3** | The `n ↔ p_F` spin-degeneracy normalization is not independently falsified inside RFG1–RFG11: the test oracle reuses `p_F = hbar c (3 pi² n)^{1/3}` and the same energy prefactor. RFG4 pins the energy normalization relative to `mu` but not the degeneracy. *This review pins it externally, from first principles and via the Table-1 `ρ_c` consistency check* | **FUTURE HARDENING** | One oracle that derives `n` from the phase-space integral `2/(2 pi hbar)^3 ∫ d³p` rather than the closed form |
| **R4** | RFG5's "independent" Hessian oracle hard-codes the same 3×3 structure as production; only the `D_i` are independent. Structural independence is in fact supplied by RFG11 (`S_x^T K S_x`) and RFG6 (column-wise finite differences) | **NO ISSUE** for coverage; documentation nit | Note in the record that RFG11, not RFG5, is the structural falsifier |
| **R5** | RFG7's positive-definiteness/conditioning grid covers `n_B ∈ [onset+0.08Δ, onset+0.88Δ]` and `n_mu/n_B ∈ [0.008,0.04]` only — not the sub-onset part of the **declared `Evaluate` domain**, and not `n_mu → 0⁺` where `κ_∞` diverges. PD itself is analytically guaranteed on the whole domain (`H=S_x^T K S_x`, `S_x` full column rank, `D_i>0`) | **FUTURE HARDENING** | Either extend the grid to the declared domain or state explicitly that `κ_∞ <= 51.4` is a grid statement and PD is analytic |
| **R6** | RFG9 as committed is an unconditional print, not a runtime falsifier; the mutation evidence lives only in the record. Runtime falsification of the solver **is** supplied by RFG8 and RFG10 | **NO ISSUE** (the printed message is honest) | Optional: rename the line so it cannot be read as a runtime check |
| **R7** | `Σ−`-onset bracket robustness: `lower = nextafter(muon_onset,+inf)` leaves the **inner** equilibrium bracket valid by `~1 ulp`, and `BisectDecreasingRoot` throws if `F(lo) < 0` — so a 1-ulp rounding change could make the **constructor** throw. Bracketing also evaluates the model up to `~0.914 fm^-3`, past the source ceiling (internal only) | **FUTURE HARDENING** (robustness) | Start the outer bracket a small relative margin above muon onset rather than one ulp |
| **R8** | `sizeof(long double)==8`, `__LDBL_DIG__==15` on arm64 macOS (measured): the `long double` intermediates give **no** extra precision on the validation platform. The small-`x` series is therefore load-bearing — and is verified exactly correct through `O(x^11)` (`8/3, 4/5, -1/7, 1/18, -5/176`), with `x^13` truncation `6.3e-23` relative at the `x=1e-2` switch | **NO ISSUE** — portability/documentation note | Record the platform fact so the series is never treated as redundant |
| **R9** | `ColdRelativisticFreeLepton::Evaluate` now ignores its own `rest_mass_energy_MeV_`/`hbar_c_MeV_fm_` and rebuilds from `kind_` — unreachable today, latent divergence hazard. RFG1's bit-identity check is tautological by delegation and omits `number_density_fm3` / `rest_mass_energy_MeV` | **FUTURE HARDENING** | Delegate through the stored members, or retire the wrapper once no caller needs it |
| **R10** | `docs/architecture/CURRENT_ARCHITECTURE.md` §2 "EOS — Phase 5A local thermodynamics" table gained no row for `ColdRelativisticIdealFermion` or `TrackRFreeGasThermodynamicProvider`, and no `GOVERNANCE.md` §5 candidate marking (the implementation record does carry the candidate statement). `GOVERNANCE.md` §2 asks for `CURRENT_ARCHITECTURE.md` to be updated in the same change for structural work | **FUTURE HARDENING** — documentation | Add the two rows with status and candidate marking |

**No finding is a BLOCKING LOCAL-PROVIDER DEFECT.** No production formula, sign,
unit, domain rule, or convention was found to be wrong.

## 17. Validation actually run in this review

Per the review validation policy, the ~11-minute full suite was **not** rerun.

| Action | Result |
|---|---|
| Configure + build `build-selfcontained` at `6cc1bbc` | success (exit 0) |
| Focused `rotochemical_trackr_freegas_local` | **RFG1–RFG11 PASS**, 0.25 s |
| Focused `rotochemical_local_thermodynamics` | **V1–V10 PASS**, 0.04 s |
| Self-contained registered-test inventory (`ctest -N`) | **24** — matches the record's 24/24 |
| Full registered-test inventory (`ctest -N`, `COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`) | **47** — matches the record's 47/47 |
| Complete self-contained suite | **not rerun** (policy); the record's 24/24 in 92.23 s stands |
| Complete full suite | **not rerun** (policy); the record's 47/47 in 684.10 s stands |
| Eight governed baselines | re-hashed live; all eight match |

Every measured number in the implementation record's RFG table was reproduced
**exactly** by the focused run, including
`muon onset = 0.45698480541241987`, `Sigma-minus onset = 0.61735520796653498`,
`RFG1 max energy rel. err = 3.0231150144907791e-15`,
`RFG2 3.8594513925747193e-16`,
`RFG3 eta_npe = -92.634662115730947 / eta_npmu = -72.631533467752774`,
`RFG4 4.7153037030511769e-07`, `RFG5 0`,
`RFG6 2.7530119496077532e-05`,
`RFG7 asymmetry 0`, minors `140.14373125628921 / 237340.63311412412 / 555320814.43648303`, `κ_∞ = 51.4353890897873`,
`RFG8 |eta| = 1.1510792319313623e-12`, min rise `1.3918548802394071e-07`,
`RFG10 mu_e = 105.03583860356072`, independent Σ onset `0.61735520796652943`,
`RFG11 0`.

Independent scratch calculations (60-digit decimal, no shared code) are recorded
in sections 5–11. Scratch programs and scripts were kept in the session
scratchpad and are **not** in the tree.

## 18. Final disposition

**PHASE 5A-3 ACCEPTABLE WITH NONBLOCKING FOLLOW-UP —
READY FOR CANONICAL INTEGRATION**

No production source or test was modified by this review. The only permanent
change is this file.

Per `GOVERNANCE.md` §5, the reviewed code remains an **unverified scientific
candidate** until a human with domain authority ratifies it; this independent
review is evidence, not ratification.

**Exactly one recommended next action:** fast-forward the complete Phase 5A-3
Track-R free-gas local provider and this independent review to canonical
`master` before beginning whole-star reproduction. Do not begin the whole-star
work automatically; open **R1** as the first governed task of that stage.
