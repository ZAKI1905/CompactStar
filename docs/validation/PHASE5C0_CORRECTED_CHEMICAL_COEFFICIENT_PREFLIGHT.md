# Phase-5C-0 — corrected chemical-coefficient scientific preflight

**Status:** PHASE-5C CORRECTED CHEMICAL-COEFFICIENT PREFLIGHT COMPLETE — PROPOSED CONTRACT REQUIRES INDEPENDENT REVIEW AND OWNER RATIFICATION

**Date:** 2026-09-06 (local task date; work continued after midnight).
**Change class:** documentation; scientific-semantic and architecture **proposal only**.
**Starting SHA:** `49ab2b8c2881b6ef7b9309307d18cea51d557f72`.
**Branch:** `analysis/phase5c-corrected-chemical-coefficients-preflight`.
**Worktree:** `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-chemical-preflight`.
**Companion:** `docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md`, **PROPOSED**.

This derives a coherent proposed cold bulk coefficient contract. It does not implement or
validate production chemical coefficients. All new numbers below are scratch diagnostics,
not reference targets, baselines, source-core reproductions, or realistic-EOS claims.
Independent review and owner ratification have not happened.

## 1. Authentication, authority, and exact scope

Canonical `/Users/keeper/Documents/CompactStar/repo/CompactStar` was clean, including
untracked status, on master. HEAD, local master, cached origin/master, and live remote master
were the starting SHA above. The proposed local branch, remote branch, and worktree path were
absent before creation. A fresh branch/worktree was created from that exact master. Other
registered worktrees were preserved. Before documentation writes, `git log --all -- <both new
paths>` was empty and neither path existed in any registered worktree; there was no competing
file history to adjudicate. Evidence: `/private/tmp/compactstar-phase5c0-preflight/pre-document-authentication.txt:1`.

The complete AGENTS, GOVERNANCE, ADR-0010/0011, Phase-5A audit/derivation/adjudication and
all Phase-5A implementation/review/ratification records, Phase-5B integration and
ratification, invariant register, roadmap, and current architecture were read. Source-level
APIs were authenticated afresh. Authority order and documentation requirements are
`GOVERNANCE.md:17`, `GOVERNANCE.md:45`; proposed decisions are permitted as the report of an
unresolved contract (`AGENTS.md:81`, `GOVERNANCE.md:64`). No accepted contract was amended.
No actual architecture ownership moved, so CURRENT_ARCHITECTURE remains byte-unchanged.

Current controlling status is Phase-5B **COMPLETE / GOVERNED**, INV-09 **VERIFIED / RESOLVED**
for ordinary-NStar structural response, ADR-0010 **ACCEPTED**, INV-11 **UNRESOLVED**.
Older experiment-time status strings in documents and the immutable structural JSON remain
historical; they do not override the integration record (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:100`,
`docs/SCIENTIFIC_INVARIANTS.md:739`, `docs/SCIENTIFIC_INVARIANTS.md:946`). The source-qualified
free-gas M_max question remains unresolved. MixedStar is separate.

Permanent allowlist consists of this record and proposed ADR-0013 only. No production C++,
test, baseline, EOS data, literature, build configuration, or evolution code is changed.
No existing science suite was rerun or baseline regenerated for this documentation task.
The fresh scratch build and standalone probes below are explicitly outside Git.

## 2. Primary-source authentication and supersession

Root `/Users/keeper/Documents/CompactStar/literature` was read-only. All 22 entries in its
SHA256SUMS manifest verified at closure. Filenames below are relative to `rotochemical/`.

| ID | Exact PDF | SHA-256 | Authority and read locations |
|---|---|---|---|
| R06 | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8` | Primary corrected authority. PDF p.2 / journal 569: sections 2–3, (2),(4)–(8); PDF p.3 / journal 570: (9)–(19), footnote 4; PDF p.4 / journal 571: section 5, Figures 1–2, conclusions. All three pages rendered and visually inspected. |
| F05 | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88` | Non-superfluid benchmark. PDF pp.3–5: (7)–(17); pp.6–8: (24),(30),(31); pp.8–13: all section 3 including (50)–(58); pp.21–23: appendices A–B. |
| R95 | `1995-Reisenegger-Rotochemical-Heating.pdf` | `9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff` | Supporting neutral composition and diffusion discussion, text pp.4–6 / PDF pp.14–16. Its delta-mu sign is opposite to the eta used here. |
| Y20 | `2020-Yanagi-NS-Therm-Thesis.pdf` | `69590539c275fa679a5521a9c5abedd9fdc58718b1554827786c8f707bc618cc` | Supporting section 4.1.1, thesis pp.58–60 / PDF pp.60–62, especially (4.6)–(4.12). Its summary does not replace R06 electrostatic authority. |

These locators are based on the actual PDFs. Some older Phase-5A-0 page numbers differ and
must not be copied; no historical file was rewritten. No gravitochemical source is used.

| F05 item | R06 disposition |
|---|---|
| (9), uniform redshifted **intrinsic charged** potential deviations | Replace by total electrochemical potential, R06 (2),(10); neutral combinations remain gauge-independent. |
| (12), unprojected integrated intrinsic susceptibility | Replace by charge-projected R06 (13). |
| (13), full four-species inverse | Forbidden; invert physical reduced response only, R06 prose after (13). |
| (54)–(56), old scalar Z | Replace by R06 (16)–(18). |
| (57),(58), W relations | Algebraic form survives using corrected Z and I_p=I_e+I_mu; R06 (19) gives the direct lepton form. |
| (14)–(17), (24),(30),(31), reaction sign and structural conventions | Retain with corrected potential semantics, explicit domains, and governed Phase-5B mapping. |
| Final quasi-steady temperature agreement | Survives correction; cannot discriminate old and corrected coefficient mathematics (R06 section 5). |

A bounded web erratum search and the primary [arXiv record](https://arxiv.org/abs/astro-ph/0606322)
(v1, journal DOI 10.1086/506601) yielded no authenticated erratum used here. This is not a
claim that no erratum exists. The journal PDF remains the authenticated numerical/source
object. The eq.(11) issue below is **INFERRED PRINTED/SOURCE OMISSION**.

## 3. Independent electrostatic derivation and physical null mode

Species order is `(n,p,e,mu)`, charges in units of proton charge are `q=(0,1,-1,-1)^T`.
Intrinsic chemical potentials include rest energy and arise from the local energy. R06 (2)
is `mu_i^infinity=(mu_i+q_i psi)e^Phi`, constant in a connected diffusive reservoir.
Linearizing (10) gives

```text
delta mu^infinity = e^Phi [delta mu + q delta psi + (mu+q psi) delta Phi].
```

R06 separately neglects delta-Phi for the small chemical perturbation (Cowling). The background
Phi is retained. Set `u=e^-Phi delta mu^infinity` and intrinsic `chi=partial n/partial mu`.
Then `delta n=chi(u-q delta psi)`. Enforcing `q^T delta n=0` yields

```text
delta psi = (q^T chi u)/(q^T chi q)
          = e^-Phi (q^T chi delta mu^infinity)/(q^T chi q),
C_projected = chi - chi q (q^T chi q)^-1 q^T chi,
delta n = e^-Phi C_projected delta mu^infinity.
```

The denominator is positive for an SPD intrinsic toy. This full-intrinsic representation is
an independent derivation/validation capability, not a requirement on a neutral provider.
The journal's printed (11) lacks the e^-Phi present in the result obtained by substituting
(10); its (13) agrees with this derivation. Verdict: **INFERRED PRINTED/SOURCE OMISSION**,
not a silently corrected quotation or a published erratum.

The beta stoichiometric vectors `(1,-1,-1,0)` and `(1,-1,0,-1)` annihilate q.
Consequently `eta_l^infinity=e^Phi(mu_n-mu_p-mu_l)` has no psi. A precise qualification to
“cancels only after neutrality”: cancellation in the **neutral reaction combination** is an
algebraic charge identity even before solving the density constraint. The density response
must nevertheless enforce neutrality before eliminating psi. Cancellation from eta does not
justify setting delta-psi=0 in the unconstrained density equations. Doing that produces the
F05 defect `delta n=chi u`, generally with `q^T delta n != 0`.

Symmetry of chi gives `C_projected q=q^T C_projected=0`, and therefore
`row_p=row_e+row_mu` and the matching column identity. With SPD intrinsic chi, rank is exactly
3 and the only null direction is q. After any positive scalar GR weighting and integration,
the full paper response retains this charge/gauge null mode; absent physical support can
lower rank further. **A full four-species corrected inverse is mathematically forbidden.**
A pseudoinverse is not the source's reduced physical inverse and is not silently permitted.

## 4. Neutral coordinates, Hessian, and reduced equivalence

ADR-0010's local x chart and conjugates follow directly from the energy differential:

```text
x=(n_B,n_e,n_mu)^T,  y=(n_n,n_e,n_mu)^T,
n_p=n_e+n_mu,  n_n=n_B-n_e-n_mu,
S_y=[[1,0,0],[0,1,1],[0,1,0],[0,0,1]],
T=[[1,-1,-1],[0,1,0],[0,0,1]],
S_x=S_y T=[[1,-1,-1],[0,1,1],[0,1,0],[0,0,1]],
y=T x, U=T^-1=[[1,1,1],[0,1,0],[0,0,1]].

d epsilon = mu^T S_y dy = g_y^T dy = (T^T g_y)^T dx,
g_y=(mu_n,mu_p+mu_e,mu_p+mu_mu)^T,
g_x=T^T g_y=(mu_n,-eta_npe,-eta_npmu)^T.
H_x=partial g_x/partial x=T^T H_y T,
C_x=H_x^-1,
C_y=T C_x T^T=H_y^-1.
```

The chart is linear; no nonlinear-coordinate Hessian terms are missing. H is taken at cold
beta equilibrium but with the **other independent composition coordinates held fixed**,
not differentiated along the equilibrium curve. Rest-mass linear terms vanish from second
derivatives but fix equilibrium/onset/conjugate conventions. A barotropic derivative cannot
supply H (`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:157`;
`CompactStar/EOS/LocalThermodynamics.hpp:98`, `:300`).

For an independent intrinsic extension `K=chi^-1`, restrict the energy: `H_x=S_x^T K S_x`.
For any forcing u, minimize `1/2 delta n^T K delta n-u^T delta n` subject to
`q^T delta n=0`. The Lagrange-multiplier solution is `C_projected u`; parameterizing the same
unique minimizer as `delta n=S_x delta x` gives `S_x H_x^-1 S_x^T u`. Hence

```text
C_projected=S_x H_x^-1 S_x^T=S_y C_y S_y^T.
```

This proves equivalence for every forcing, not just the charge-null identity. Individual
charged potentials, K, chi, and psi cannot be reconstructed uniquely from a neutral H.
The public neutral path is **already corrected**, and never projects C_y again.

## 5. Exact coupled intrinsic toy and wrong routes

Scratch energy is `epsilon(n)=1/2 n^T K n` about an arbitrary positive neutral reference
state, with arbitrary consistent unit scales. Adding linear rest energies does not alter its
response. The explicit matrix has charged cross-couplings and is strictly diagonally dominant
with positive diagonal, hence SPD:

```text
K = [[4,1,1,0],[1,5,0,1],[1,0,6,1],[0,1,1,7]],
H_x = [[4,-2,-3],[-2,11,8],[-3,8,16]],
C_projected = (1/381) [[105,-18,-21,3],
                      [-18,43,29,14],
                      [-21,29,55,-26],
                      [3,14,-26,40]].
```

Python Fraction elimination independently computes both routes. Equality residual **exactly
0**, symmetry residual **0**, charge-null residual **0**, proton identity residual **0**,
rank **3**; full inverse fails at an exact zero pivot. This is exact rational evidence, not
a floating-point near-singularity test (`/private/tmp/compactstar-phase5c0-preflight/toy.txt:1`).

| Wrong route | Scratch result / required detector |
|---|---|
| Omit electrostatics | `chi q=(-18,170,-98,-113)^T/719 !=0`; exact gauge forcing produces nonzero charged density. |
| Apply the same projector again | With `P=I-chi q(q^T chi q)^-1 q^T`, **P C_projected=C_projected exactly**. A numerical-output detector cannot detect this idempotent no-op. Reapplying the susceptibility subtraction instead encounters `q^T C_projected q=0` and must refuse. C_y has no four-species charge coordinate. Future typed API/call-count or dependency audit detects the forbidden extra operation; do not count the null mutation as independent validation. |
| Invert full corrected matrix | Exact rational zero pivot; full inverse refused. |
| Impose baryon conservation locally | Independent two-zone example below changes the exact result. |
| Flip eta sign | For positive lepton excess `(1,2)`, the two-zone result requires `eta=-(21,29)/53`; a plus sign produces its negative and fails the exact signed-vector oracle. This is an algebraic detector, not executed production mutation. |

For the local/global baryon counterexample choose two SPD local C_x matrices with unit
positive volume weights:

```text
C1=[[3,1,0],[1,2,0],[0,0,1]], C2=[[2,0,1],[0,1,0],[1,0,3]].
Schur(C1+C2) = [[14,-1],[-1,19]]/5,
Schur(C1)+Schur(C2) = diag(8/3,7/2),
Z_correct=[[19,1],[1,14]]/53,
Z_wrong=diag(3/8,2/7).
```

Thus integration and the baryon Schur reduction do not commute. A one-zone fixture could
miss this defect. Identities built from the same H/projection are **CONTRACT / CONSISTENCY**;
this coupled intrinsic energy and two-zone constraint example provide independent analytic
oracles for the future consumer. They are not realistic nuclear-matter validation.

## 6. GR integral, dimensions, and redshift

Use the nonrotating reference star of R06. CompactStar's `nu` is Phi itself, not 2Phi:
`g_tt=-e^(2nu)`, `dV=4 pi r^2 (1-2m/r)^-1/2 dr` in km^3. The canonical mathematical measure
owner is `CompactStar::Geometry`; consumers own chemical factors and unit conversion
(`docs/SCIENTIFIC_INVARIANTS.md:186`, `:231`; `CompactStar/Geometry.hpp:1`).

For a declared connected diffusive chemical domain D, `delta g_y(r)=e^-nu(r) delta g_y^infinity`.
Integrating the number excess at fixed Cowling background therefore gives

```text
delta N_y = G_y delta g_y^infinity,
G_y = 10^54 integral_D [4 pi r^2 e^-nu / sqrt(1-2m/r)] C_y(r) dr.
```

Rows label `(N_n,N_e,N_mu)`; columns label `(g_n^infinity,g_e^infinity,g_mu^infinity)`.
G_y corresponds to R06's reduced Btilde (delete proton row/column from the projected full
response). The full matrix is `S_y G_y S_y^T`. There is exactly **one inverse lapse**, due to
local energy versus energy at infinity, and **no extra time lapse in a particle count**.
Rotation/velocity number perturbations belong to the separately governed Phase-5B driver.
They are not added a second time to this static chemical susceptibility integral.

Particle count is dimensionless physically; the ledger retains “count” as a bookkeeping label.
Local published units suppress this label in density and its derivative. All rows use the
same particle-number normalization; “count” does not introduce a new conversion constant.

| Quantity | Unit with count bookkeeping | Factors / semantic boundary |
|---|---|---|
| x,y,n_i | count fm^-3 | ADR-0001 profiles carry Y_i; recover n_i=Y_i n_B. |
| g_x,g_y,mu_i, local eta | MeV | Intrinsic neutral combinations, rest-energy convention declared. |
| H_x,H_y | MeV fm^3 / count (normally MeV fm^3) | partial g / partial density; 1D/2D/3D only on active chart. |
| C_x,C_y | count fm^-3 MeV^-1 | inverse active H; no lapse locally. |
| dV | km^3 | Proper volume; multiply by (10^18 fm/km)^3=10^54 once. |
| G_y; G_x; baryon-constrained number response Q | count MeV^-1 | G integral has e^-nu; Q is a global Schur complement. |
| G_y^-1; chemical Z | MeV / count | No additional redshift or speed-of-light conversion. |
| delta N_y, delta N_l | count | actual minus corresponding equilibrium reservoir counts. |
| eta^infinity | MeV | e^nu eta_local; uniform under the stated diffusion assumptions. |
| structural K, I_geom | count km^2 | Whole-star conditional mapping I_geom=K; core mapping differs. |
| I_phys | count s^2 | I_geom/c^2 once, c=299792.458 km/s, existing AngularVelocity owner. |
| Omega, Omega_dot | s^-1, s^-2 | Physical spin; radians dimensionless. |
| chemical W=Z I_phys,l | MeV s^2 | 2 W Omega Omega_dot has MeV/s. |
| net local reaction density rate DeltaGamma | count fm^-3 s^-1 (or declared cm^-3 s^-1) | Local proper time. No rate implemented. |
| reaction number rate R_l^infinity | count/s | integral e^nu DeltaGamma dV; 10^54 for fm or 10^15 for cm. This **positive** lapse belongs to rate time conversion, not G. |
| paper Figure-1 Z | erg/count, plotted in 10^-60 erg | 1 MeV=1.602176634e-6 erg; preserve source plotting scale. |

Later energy power has an additional energy-redshift factor; neither its implementation nor
its partition into heat/neutrinos is part of this contract. No coefficient unit remains ambiguous.

## 7. Global baryon reduction and corrected chemical Z

Local neutrality does not fix local n_B. Globally the closed modeled reservoir obeys
`delta N_B=delta N_n+delta N_e+delta N_mu=b^T delta N_y=0`, with `b=(1,1,1)^T`.
The common baryon chemical-potential direction is delta-g_y proportional to b. Let

```text
L=[[-1,-1],[1,0],[0,1]],
delta N_l=(delta N_e,delta N_mu)^T,
delta N_y=L delta N_l,
eta^infinity=-L^T delta g_y^infinity,
Z=L^T G_y^-1 L,
eta^infinity=-Z delta N_l.
```

This is R06 (14)–(18). In order `(npe,np-mu)`, Z is the named symmetric 2x2 matrix
`[[Z_npe,Z_np],[Z_np,Z_npmu]]`. With `M=G_y^-1`, its entries are

```text
Z_npe = M_nn - 2 M_ne + M_ee,
Z_npmu = M_nn - 2 M_nmu + M_mumu,
Z_np = M_nn - M_ne - M_nmu + M_emu.
```

Strict positivity of G_y on its supported 3D space makes Z SPD. The cross entry need not be
positive for every EOS. Symmetry relates row/output imbalance and column/input lepton excess;
it is not permission to disregard named channel orientation.

Independent reduction route: `G_x=U G_y U^T=[[a,h^T],[h,D]]`. With conjugates
`delta g_x^infinity=(alpha,-eta^infinity)`, global delta-N_B=0 gives
`alpha=h^T eta^infinity/a`. Therefore

```text
Q = D - h h^T/a,
delta N_l = -Q eta^infinity,
Z=Q^-1.
```

Q has dimension 2x2 and count/MeV units, while G_y has dimension 3x3. Both are regular if the
physical modes have sufficient support and strict stability. This route provides a useful
conditioning diagnostic, not a license to apply the Schur complement at each radius. R06
footnote 4 explicitly distinguishes global baryon conservation from local neutrality.

## 8. W, signs, and the governed structural boundary

The current semantic consumer is `CompactStar::Analysis::FixedBaryonNumberResponse`, retaining
its complete metadata, live sources, contributing stars, and validity. Whole-star W consumes
its `WholeStarIPhysical()`, after `RequireCurrent()`, and maps explicit species labels to e/mu.
For checking an actual supplied spin, `WholeStarEquilibriumNumberRate(omega,omega_dot)` is the
higher-level number-rate action. No K/c conversion is duplicated by the chemical consumer.
The API is authenticated at `CompactStar/Analysis/ParticleNumberResponse.hpp:135` and
`CompactStar/Analysis/src/ParticleNumberResponse.cpp:511`. Production must **not** parse
`tests/baselines/phase5b_structural_response.json`.

The baseline SHA is `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`.
It is validation authority; its preserved candidate-era status strings are superseded by the
integration record, not edited here. This scratch W uses its authenticated named numerical
values, without claiming to have constructed a runtime provenance-compatible chemical object.

| Sign definition | Consequence |
|---|---|
| delta-N_l=N_l-N_l^eq | Excess means actual minus equilibrium, R06 (7). |
| eta_l=mu_n-mu_p-mu_l | eta^infinity=-Z delta-N_l. |
| R_l>0 denotes net neutron beta decay producing lepton l | delta-dot-N_l=R_l-dot-N_l^eq. |
| Phase-5B dot-N_l^eq=+2 I_phys,l Omega Omega_dot | Explicit source sign, F05 (30), R06 (5). |
| Fixed background coefficients | dot-eta^infinity=-Z R+2(Z I_phys,l) Omega Omega_dot. |
| W=Z (I_phys,e,I_phys,mu)^T | Matches R06 (19); both source signs reconciled, no analogy-only inference. |

The differential equation is a source-level sign derivation, **not an ODE implementation**.
Frozen coefficients are the approximation of the source; changing coefficients would add a
term from differentiating Z and requires later evolution governance. For this scratch free gas
I_e and I_mu are negative and both W entries are negative, so physical spin-down with
Omega>0, Omega_dot<0 gives positive driving in both eta channels. This entrywise sign is not
postulated for all EOSs. The 1995 delta-mu convention must be explicitly negated before use.

`FixedIsobarIGeometric(outer_Y,outer_Q_B,inner_Y,inner_Q_B)` carries the distinct PN8 boundary
mapping (`CompactStar/Analysis/src/ParticleNumberResponse.cpp:517`). There is no existing
high-level physical core driver with all chemical-reservoir assumptions. Initial whole-star
chemical work must refuse a core/shell W request until that separate adapter is governed.
Equal domain names alone do not establish equal reservoir or boundary-flux semantics.

## 9. Chemical, structural, reaction, and support domains

| Object | Source/model domain | Required qualification |
|---|---|---|
| Background TOV and total baryon constraint | Whole star | Physical P=0 or canonical finite cut with explicitly bounded tail. |
| F05 (12), corrected R06 counterpart | Core in F05 section 3.5; R06 (12),(13) do not print a new domain | R06 retains the framework; its prose “star” supplies no new numerical core boundary. Preserve core scope for source reproduction. |
| F05 (30) I_Omega | Core composition drive | Governed Phase-5B distinguishes whole-star K/I and fixed-isobar core mapping; no invented cutoff. |
| Net reaction count | Core, and for each process its actually active subregion | R06 (6),(7), F05 (15),(43); direct-Urca domain may be smaller than susceptibility domain. |
| Specific heat, F05 (50) | Each free species' support | Not identical to core chemistry or a direct-Urca region. |
| Track-R free-gas diagnostic in this task | All physical pe/npe/npemu support, finite cut plus assessed remainder | Whole-star mathematical validation fixture, **not** a source-core coefficient. |

No source-authenticated free-gas core cutoff exists in the current authority. None is invented.
Whole-star free-gas G/Z/W can test the mathematics and numerical coupling. It cannot establish
a published core I_Omega, a corrected A18 coefficient, or source-qualified M_max. A core
reservoir's own baryon closure/exchange prescription and structurally compatible drive must
be explicit before a source-core calculation; global whole-star conservation cannot be
silently substituted for core conservation. The later thermal benchmark requires authenticated
core/crust joins, free-particle and reaction support, phase treatment, and transport assumptions.

## 10. Active-species embeddings and interfaces

Use the actual `ActiveLocalThermodynamicEvaluation` variant, not a padded H
(`CompactStar/EOS/LocalThermodynamics.hpp:277`, `:319`). For each active chart z,
`C_y=E H_z^-1 E^T` with

| Branch | z, active H dimension | E into `(n_n,n_e,n_mu)` |
|---|---|---|
| npemu | `(n_B,n_e,n_mu)`, 3 | T from section 4 |
| npe | `(n_B,n_e)`, 2 | `[[1,-1],[0,1],[0,0]]` |
| pe | `(n_B)`, 1 | `[[0],[1],[0]]` |
| vacuum | value-only boundary, no H | Zero susceptibility as the one-sided pe limit; not an inverse of a 0D/padded matrix. |

Thus no local fictitious absent-muon or absent-neutron susceptibility is inverted. In pe the
active conjugate is mu_p+mu_e, not mu_n. Beta channels absent locally are not independently
constrained there. Their global response can still be supported elsewhere in the connected
reservoir. Global rank is determined by the union of local ranges: every nonzero global
forcing must have nonzero response on a set of positive measure. If muons are absent everywhere,
reduce to the supported `(N_n,N_e)` space and one beta channel explicitly; all-pe matter has
only one count mode and no active beta channel. Do not retain an invertible 2-channel Z by fiat.

At continuous cold free-Fermi onset, newly appearing density is proportional to the positive
chemical-potential gap to the 3/2 power and its susceptibility vanishes as the square root.
The embedded susceptibility, not the divergent H, has the proper lower-branch limit. The
moving continuous onset supplies no first-order density-jump atom. It still produces a
nonsmooth integrand and demands onset-aware partition/convergence control.

Exact provider threshold objects are value-only. Quadrature must avoid endpoint H requests or
use a separately certified one-sided C limit. The positive-width near-neutron-onset response
refusal is different: nodal sampling that misses it does not close the gap. Production requires
a validated bounding/limit adapter with its own error budget or must refuse. No density floor,
state substitution, extrapolated full H, or silent omitted cell is accepted. In source,
`n_n >= 2^30 [n_B-nextafter(n_B,0)]` and separate bracket guards remain unchanged
(`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:460`).

A genuine first-order phase transition requires authenticated phase/interface metadata and a
chemical perturbation/motion law. The bulk C integral alone need not include its chemical
interface response. F05 (31) and Appendix B govern the structural equilibrium drive's jump,
not an automatically reusable chemical interface term. No jump is inferred from steepness.
Crossing an ungoverned interface refuses; no APR Maxwell interface implementation is supplied.

## 11. Fresh Track-R scratch method and independence boundaries

All work is under `/private/tmp/compactstar-phase5c0-preflight/`. Fresh Debug CMake build of
current CompactStar used AppleClang 17, Python 3.12, GSL 2.7.1 and installed Zaki archive.
`probe.cpp` builds the canonical Structure-1 table with the existing test fixture generator,
uses `SingleStarSolveToTOVPoints(1.10e15,...)`, requires SURFACE_REACHED, constructs NStar,
and calls the **actual** provider `EquilibriumAt(n_B)` at every profile node. H entries,
active dimensions and states are dumped without changing production. H inversions and GR
integration are separate NumPy/SciPy scratch code (`analyze.py`). Typed absence is retained;
zero slots in the dump are serialization padding and are never included in a local inverse.

The 8192-interval table SHA is exactly
`7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a`, matching the governed
Phase-5B fixture. Its upper density is .616 fm^-3 and floor 1e-14 fm^-3; the actual pressure
cut sets n_cut=6.281380783102176e-12 fm^-3. Center n_B=.60492524140311366 fm^-3,
R_cut=12.766174760730493 km, M_geom=.92093273643866647 km. This is the ratified midpoint
fixture, not a newly selected extremum. Table construction is `tests/eos/structure1/table.hpp:49`;
source and scratch hashes are `/private/tmp/compactstar-phase5c0-preflight/source-hashes.json:1`.

Primary scratch quadrature is trapezoidal on the complete returned radial partition, plus the
regular-center `4 pi r0^3 e^-nu_c C_c/3` contribution (relative norm 1.29e-18). Radial targets
20000/40000/80000 were varied at fixed table 8192, and table intervals 4096/8192/16384 at
fixed radial 80000. Simpson comparison is a numerical consistency check on the same sampled
integrand, not independent physics validation.

`oracle.py` instead solves free-gas equilibrium in the common lepton potential, obtains
intrinsic susceptibilities directly from phase space
`chi_i=mu_i p_Fi/(pi^2 (hbar c)^3)`, projects four-species charge and selects the n/e/mu
submatrix. It never calls the production H, inverse, equilibrium solver, or geometry helper.
It integrates with 8/16/32-point Gauss rules, splits at the physical onsets in an explicitly
linear n_B(r) background, and evaluates the analytic local response at quadrature points.
It shares the current star profile and authenticated particle constants, so it is independent
of local implementation and quadrature, **not** an independent TOV/EOS-background validation.
The analytic oracle can traverse a production response-refusal window only as **scratch
source evidence**; it is not a newly authorized runtime fallback.

A separate manufactured metric oracle uses constant C, m=0, and
`nu=-ln(1+a r/R)`, a=.4. Its exact integral is `4 pi R^3(1/3+a/4) C`, and Simpson yields
relative error 2.22e-16. This is a test of integration/redshift, not an asserted TOV solution.
An independent curved-metric exact-star oracle remains a future GC9 requirement.

## 12. Numerical diagnostics, conditioning, and convergence

At table8192/radial80000 there are 20259 nodes, zero sampled provider refusals, and:

| Branch | Sampled radial support (km) | Nodes | Eigenvalue range H (MeV fm^3) | kappa_2(H) range |
|---|---|---:|---|---|
| npemu | 1e-5 to 3.27866515662 | 9416 | 133.0813 to 735974.95 | 78.3802 to 4873.2302 |
| npe | 3.27954015650 to 12.64291381887 | 10702 | 151.0583 to 4.9313012e7 | 51.0469 to 9463.5087 |
| pe | 12.64378881875 to 12.76617476073 | 141 | 4.9738772e7 to 1.2873618e9 | 1 to 1 |

The small intervals between sampled branch endpoints contain the continuous onsets; they
are not omitted physical shells. No exact boundary H was queried. Separate near-onset probes
at relative offsets 10^-2 through 10^-12 reproduced npe refusals on the positive side at
10^-7,10^-8,10^-10,10^-12 and reached kappa_2(H)=4.1789e7 above muon onset. Thus the nodal
maximum 9463.51 is **not** a universal condition bound. Physical H conditioning diverges at
the appropriate limiting boundary. Evidence: `numeric.json:1`, `oracle.json:1`, and
`local-probe.txt:1` under the scratch root; full C_x/C_y node dumps are in each run's matrices.npz.

Rounded main diagnostic matrices (exact free-gas neutron/charged cross entries in G_y are
zero; floating residuals smaller than 4.7e36 count/MeV are displayed as zero):

```text
G_y [count/MeV] =
[[2.904619168725205e55, 0,                     0],
 [0,                    2.201380346798831e53, -1.035413409841426e51],
 [0,                   -1.035413409841426e51,  9.746459033521813e51]]

Q [count/MeV] =
[[ 2.184981511700372e53, -1.100611626042510e51],
 [-1.100611626042510e51,  9.743866893671463e51]]

Z [MeV/count] =
[[4.579303248926041e-54, 5.172519750054969e-55],
 [5.172519750054969e-55, 1.026870855745314e-52]]

I_phys,(e,mu) [count s^2] = (-1.1637998545112904e47, -1.4844819850233820e46)
W [MeV s^2] = (-5.406177578724549e-7, -1.584569063625169e-6).
```

kappa_2(G_y)=2981.7377879; kappa_2(Q)=kappa_2(Z)=22.4381301. G symmetry relative Frobenius
residual is 2.224e-19; x/y inverse-transform residual is at most 3.929e-16; Schur and lifted
inverse Z agree to 5.64e-18 relative Frobenius norm. These are **CONTRACT / CONSISTENCY**.
They do not make the integral independently correct. The source y basis has small roundoff
cross entries because neutral H was inverted before transforming; no ad hoc symmetrization
or projection was applied to erase them.

| Change from table8192/radial80000 | Relative Frobenius change G | Relative Frobenius change Z |
|---|---:|---:|
| radial20000 | 1.0863e-6 | 2.0694e-5 |
| radial40000 | 1.7508e-7 | 1.6151e-7 |
| table4096 | 2.2684e-9 | 4.2212e-8 |
| table16384 | 2.5963e-10 | 5.8581e-9 |

Radial differences are nonmonotone; no Richardson order is claimed. Simpson versus trapezoid
Z differs by 1.8709e-6 on the finest radial grid. Independent equilibrium/projection gives
local C agreement at most 2.199e-11; its independent integrated Z is
`[[4.579303180700017e-54,5.172519910355882e-55],[5.172519910355882e-55,1.026872797719940e-52]]`,
1.8892e-6 from the provider/trapezoid diagnostic. Independent Gauss G differences 8→16 and
16→32 are 5.526e-10 and 9.452e-11 in norm, respectively. These latter numbers only measure
Gauss convergence on the stipulated interpolated background; they do not bound star error.

The threshold-containing neutron and muon cells have G fractions 9.465e-7 and 3.295e-5;
omitting the whole cells changes Z by 8.843e-9 and 1.107e-5. These are **omission sensitivity
diagnostics**, not rigorous bounds on the narrower response-refusal window. The main star's
local-baryon wrong route changes Z by 6.998e-3. Omitting the inverse lapse changes Z by
.22065, adding another inverse lapse by .18078, and omitting proper volume by .016951.
These scratch mutants are readily distinguishable; no production mutation campaign ran.

### 12.1 Surface and remainder

The finite cut is not P=0. For cold monotone pe matter let
`h_cut=ln[(mu_p+mu_e)/(m_p+m_e)]=1.307364207e-5`. Positive density and pressure, m>=M_cut,
and the TOV enthalpy equation imply the comparison bounds

```text
R_upper=2 M_cut/[1-(1-2 M_cut/R_cut) exp(2 h_cut)],
M_upper=M_cut+(4 pi/3) epsilon_cut (R_upper^3-R_cut^3),
0 <= Delta G_ee <= 10^54 (4 pi/3)(R_upper^3-R_cut^3)
                    e^-nu_cut C_pe,cut /sqrt(1-2 M_upper/R_cut).
```

The last inequality uses monotonically increasing C_pe with density and the minimum metric
factor over the tail; it is a proposed analytic **static chemical-tail** bound, separate from
Phase-5B's structural tail bound. epsilon_cut=7.80499e-15 km^-2; R_upper=12.7681549034 km,
M_upper=.9209327364386981 km, Delta-G_ee upper=3.68181e45 count/MeV. Only the source e/e
entry has support in this tail. Adding that endpoint rank-one bound changes Z by at most
7.49e-10 relative norm for this fixed background. Matching the continued pe sphere's lapse
must also be accounted for: the possible shell-mass correction is bounded at order
`(M_upper-M_cut)/(R_cut(1-2M_upper/R_cut)) < 3e-15`; constant enclosed-mass enthalpy evolution
otherwise preserves the interior lapse normalization. Center and tail estimates do not bound
background interpolation or the near-threshold response gap. This proposal is not independent
ratification of a production P=0 coefficient. The scalar pressure cutoff and TOV solve remain
unchanged.

## 13. Conditioning and tolerance proposal

Use symmetric positive solves/factorizations on active spaces, with the mathematical inverse
notation above. No solver/library algorithm is selected or implemented by this preflight.

| Solve | Size/rank and stability | Refusal and evidence |
|---|---|---|
| Intrinsic toy chi=K^-1; charge scalar | 4 full rank SPD; scalar q^T chi q>0 | Private validation only. Invalid intrinsic model or zero denominator refuses. |
| Active local H | 1,2,3, full rank on declared smooth chart | Nonsymmetric beyond independently established numerical error, nonfinite, nonpositive eigenvalue, rank loss, absent dimension or unknown phase refuses. |
| Full projected response | 4, rank at most3, charge null | **Never solve/invert**; no pseudoinverse escape. |
| Integrated G_y | 3 if all physical modes supported | Require support proof, SPD and resolvable smallest eigenvalue; no fabricated muon support. |
| Global Schur Q / Z | 2 supported beta channels | Must remain SPD and sufficiently conditioned after propagated integration error. Lower rank requires a separately named reduced channel set. |

Proposed acceptance method, not a chosen production numeric cutoff: for each scaled matrix
use kappa_2 from its symmetric eigenspectrum; record the basis and scaling. With dimension n,
IEEE epsilon u, use `gamma_n=n*u/(1-n*u)` to organize operation-count roundoff estimates.
Backward residual `||H C-I||` and forward error amplified by kappa must both be reported.
For matrix uncertainty E, require `rho=||G^-1 E||<1`; inverse perturbation is bounded by
`||G^-1|| rho/(1-rho)`. Propagate this through L and the Schur operation, using absolute
component budgets for exact-zero or cancellation-sensitive entries. Do not use relative
errors against an analytically zero cross term. Positive-eigenvalue tests must exceed their
error enclosures; a bare positive floating eigenvalue is insufficient.

Exact rational toy identities have zero tolerance. The manufactured cubic integrand's
Simpson truncation is exactly zero; only a separately counted arithmetic budget applies.
No literal universal tolerance is inferred from the displayed scratch agreement. On the
sampled smooth star, kappa*u is roughly 2.1e-12 locally at the largest sampled H condition and
6.6e-13 for G; these are first-order sensitivity scales, not certified error bounds. The
near-muon sample already raises the local scale to about 9.3e-9. Scalar pe kappa=1 says
nothing about absolute magnitude, parameter/equilibrium error, or boundary smoothness.

**Still unresolved before a production tolerance:** onset-aware radial partitions over
multiple node phases; a proved bound for the finite response-refusal gap; local equilibrium
anchor error and table/profile consistency; independent curved-GR oracle; EOS/radial and
surface error propagation into the small global eigenmode; conditioning over a predeclared
mass/density/domain set; componentwise W uncertainty including structural K_errors/c^2;
source curve extraction uncertainty and realistic phase response. Predeclare budgets from
these inputs before testing a future implementation. The nonmonotone current grid differences
cannot select a universal order or a convenient acceptance number.

## 14. R2006 correction-sensitive benchmark inventory

| Candidate | Classification | Available information / limitation |
|---|---|---|
| R06 (13)–(19) | FORMULA ONLY | Exact algebraic corrected coefficients; source formula validation is available now. |
| R06 Figure 1 | FIGURE-DERIVED | A18+delta-v+UIX*, mass x axis 1–2 M_sun linear; Z y axis 0–2.5 in 10^-60 erg linear. Three labelled coefficients, old dashed versus corrected solid. Source pixels are available; no extraction was performed. |
| R06 section 5 change statement | EXACT TEXTUAL (approximate numerical precision) | Z_np hardly changes; Z_npe and Z_npmu fractional changes reach approximately 10–20% at the largest masses. Exact text does not make 10–20% a precision target at a specified mass. |
| R06 Figure 2 | FIGURE-DERIVED | log10 t in years versus log10 T_infinity in K and log10(eta_infinity/k) in K; corrected solid/old dashed. Left mass1.47 M_sun. Right **panel says2.14; caption says2.13**. Do not choose one silently. |
| Figure 2 initial/spin conditions | EXACT TEXTUAL | Caption: zero initial imbalances, T_infinity=1e8 K at t=0, magnetic field1e8 G, P0=1 ms, dipole braking. Left modified Urca only; high-mass electron direct Urca, muon modified Urca (section5). |
| Corrected coefficient tables | EXACT TABULAR: none located | No exact tabular Z fixture in R06 or authenticated author numerical data in the current library. |
| Author arrays/code and exact high-mass configuration | UNAVAILABLE in authenticated materials | Do not reverse-infer them from plotted line centers or substitute F05 parameters. |
| Quasi-steady temperature alone | Available but correction-insensitive | Cannot close ADR-0010's correction-sensitive gate. Transient differences and Figure1 coefficient changes can. |

The presence of reproducible source pixels means a correction-sensitive benchmark is
**potentially reproducible**, not wholly unavailable. It is not yet an executable quantitative
benchmark: no governed extraction or authenticated matching EOS exists. A Figure1 comparison
can avoid the Figure2 high-mass ambiguity; it does not avoid missing EOS authority.

Required extraction protocol: authenticate R06 SHA above, render journal571/PDF4 at at least
two resolutions, record crop/panel/axes/ticks and pixel mapping, identify each coefficient and
old/new curve, and preselect masses away from unresolved overlaps. Save both line edges,
line-width and pixel quantization uncertainties, axis-fit residuals, and curve ambiguity.
Use linear transforms for Figure1 and log10 transforms for Figure2, converting eta/k to MeV
only with the declared k. Require an **independent second extraction**, with no access to the
first picked points, and compare the traces before combining uncertainties. Record
interpolation uncertainty separately, including any interpolation between masses/times; do
not treat interpolated pixels as independent samples. Select comparison metrics and bounds
from those envelopes before running a coefficient implementation. No digitized value may be
labelled exact. No figure-derived targets were created in this task. Figure2 high-mass
comparison additionally requires source/owner adjudication of 2.13 versus 2.14.

## 15. A18 + delta-v + UIX* closing authority

Inventory was checked against the current literature catalog and filenames under the shared
literature/data roots; no APR/Akmal/A18/UIX*/BPAL source product was located. This is a
bounded local inventory, not a claim of global nonexistence. F05 section3.1 identifies a
many-body analytic fit, Maxwell joining to the neutral-pion phase near4e14 g/cm^3 with6.6%
energy-density jump, and crust sources; it does not supply a complete authenticated runtime
product. A18+delta-v+UIX* is the accepted closing EOS, not CMF or BPAL
(`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:357`).

| Input | Status | Exact missing authority / next-stage role |
|---|---|---|
| Model identity and published methodology | AVAILABLE NOW | F05 section3.1 and R06 Figure1 specify A18+delta-v+UIX*, source references and comparison scope. |
| Equilibrium EOS/composition used by authors | MISSING authenticated product | Exact APR revision/fit coefficients/tables, units, valid domains, equilibrium reconstruction and match to source. A generic downloaded APR barotrope would not suffice. |
| Arbitrary-composition nuclear response | MISSING | Source-authenticated energy versus density/proton fraction, its derivatives and phase-specific off-equilibrium conventions, plus consistent independently validated leptons. This is the direct realistic chemical-layer blocker. |
| Potential online/equilibrium APR substitutes | AVAILABLE BUT UNAUTHENTICATED as benchmark substitutes | Existence of an APR-labelled product is not lineage equivalence; none was acquired/promoted here. |
| Crust and joins/core cutoff | MISSING | Exact Pethick–Ravenhall–Lorenz1995 inner crust, Haensel–Pichon1994 outer crust, joins, free-particle support and authors' core boundary. |
| Phase treatment | AVAILABLE NOW qualitatively; MISSING quantitative off-equilibrium law | F05 Maxwell choice and6.6% jump are known; phase-specific response/interface displacement and composition jump data are not authenticated. |
| Effective masses | NOT YET NEEDED for cold bulk coefficient identity; MISSING for later thermal closure | Source-specific APR/Page et al. prescription, not bare-mass substitution. |
| Corrected Z comparison | Source pixels/formula AVAILABLE NOW; exact arrays UNAVAILABLE | Governed Figure1 extraction or authenticated author data, matching EOS/mass. |
| F05/R06 thermal reproduction | NOT EXECUTABLE NOW | Above inputs plus exact weak-rate normalization/control functions, effective masses, envelope, initial/spin setup, phase assumptions, curve reference; future evolution governance and implementation. |
| R06 high-mass transient configuration | Conflicting source labels | Caption2.13 versus panel2.14 M_sun; unresolved benchmark-only source finding. |

Free gas validates coordinate orientation, charge correction, GR factors, constraints, units,
signs, supported-mode conditioning and structural coupling mechanics. It cannot validate
realistic nuclear interactions, published A18 Z, final thermal evolution, or a missing source
core I_Omega. Track R closure remains blocked; preflight mathematics is still draftable.

## 16. Proposed owners, names, immutable results and provenance

Names follow existing descriptive value/result classes and explicit method semantics
(`CompactStar/EOS/LocalThermodynamics.hpp:98`, `CompactStar/Analysis/ParticleNumberResponse.hpp:85`).
They are **recommendations**, not declarations or promises of an implemented API.

| Layer | Proposed owner / exact descriptive names | Boundary |
|---|---|---|
| A | Existing `ILocalThermodynamicProvider`, `ChargeNeutralChemicalHessian` and lower charts | Own local cold energy/conjugates/H, identity, equilibrium and branch; no stellar integration. |
| A adapter | `CompactStar::ChargeNeutralNumberSusceptibility` | Derive active C_y and named embedding, retain H/domain provenance; no public psi/projector. |
| B | `CompactStar::Analysis::GlobalChemicalNumberResponse` | Own G_y, source basis, GR integration and error/domain provenance; use Geometry measure. |
| B2 | `CompactStar::Analysis::ChemicalImbalanceResponse` | Own canonical symmetric Z and global baryon reduction, retain originating G result. |
| C | Existing `FixedBaryonNumberResponse` | Structural K and domain-qualified I/rates stay Phase-5B owned. |
| B+C | `CompactStar::Analysis::RotochemicalSpinDrive` | Separate immutable W, referencing chemical Z and complete governed structural response. |
| D | Future secular evolution owner | Consume coefficients and named sources, choose states/rates/thermal accounting later; never rebuild EOS or structure implicitly. |

Recommend named `BetaChannel::Npe`, `BetaChannel::NpMu`; source axes
`NeutralNumberCoordinate::{Neutron,Electron,Muon}`. Recommended accessors:
`ResponseCountPerMeV(row,column)` on G, `ResponseMeVPerCount(output_channel,input_lepton)`
on Z, and `DriveMeVSecondsSquared(channel)` on W. Source-facing `PaperZnpe()`, `PaperZnp()`,
`PaperZnpMu()` are views of the **one** canonical Z matrix, never duplicate scalar storage.
Both named matrix and scalar views are useful; ownership is not split between them.
No bare B/MatrixB/Z/W class is proposed. Historical structural Z_i and current structural
A_i/B_i/K_i remain distinct from chemical Z.

Keep G, Z, and W as separate immutable dependency-linked results: thermodynamic Z is usable
without a spin response; W adds structural dependence and invalidation. Every public scientific
read/action validates currency. A result must carry star/profile identity and version, metric
identity and normalization, schema/units, provider identity/revision/data bytes and component
constants, equilibrium-anchor match, active branch map, basis/orientation, domain definition,
reservoir/boundary policy, integration partitions and onset locations, interface/tail policies,
source convention (R06 correction/Cowling), numerical method version, eigenvalues/conditioning,
and error budgets. Provider identity strings alone are insufficient mutable-data provenance.

W additionally retains the **entire** structural provenance dependency set, including the
central and sequence source profiles, EOS table bytes, Hartle snapshots, surface/domain and
physical/geometric spin convention. A copied I vector is insufficient. Existing structural
`RequireCurrent` checks all sources and table contents
(`CompactStar/Analysis/src/ParticleNumberResponse.cpp:247`). Own lifetimes or an explicit
validated lifetime token are mandatory; a dangling raw pointer cannot safely be validated.
Any changed dependency, invalidated profile, mismatched provider/background, stale sequence
star, changed bytes with unchanged label, or core/whole-star mismatch refuses before returning
numbers. Rebuilding is an explicit new result, never lazy stale science or an undocumented
fall-through to historical RotochemicalCache.

## 17. Predeclared GC validation ladder — fourteen required gates

These are **proposed tests**, not fourteen achieved passes. Each row states a production target,
independent oracle, defect, independence limit, metric/tolerance source, and negative controls.
“C/C” means CONTRACT / CONSISTENCY, not INDEPENDENT PHYSICS VALIDATION. Scratch coverage in
sections5 and11–13 informs the plan but does not run or certify future production tests.

| ID / production target | Independent oracle and independence boundary | Defect; metric and tolerance derivation | Negative controls |
|---|---|---|---|
| GC1 units/index/sign/rest mass | Source energy differential, SI/fm/km and physical c literals, signed named forcing vectors; independent of consumer conversions | Wrong units, lapse, rest-energy/eta convention; exact symbolic dimensions and arithmetic error budget; no output-fitted bound | 10^54 missing/double, erg/MeV mix, sign/c^2 mutations |
| GC2 neutral reconstruction | Independent species lift from positive free coordinates and energy differential; reconstruction identity alone C/C | Wrong dependent species/order; exact rational fixtures, input-domain invalidity, float roundoff derived from operations | Wrong proton map, negative composition, density/fraction substitution |
| GC3 local H stability | Independently differentiated coupled energy and shrinking finite perturbations on declared branches, not provider H against itself | Wrong held-fixed variables/cross terms, instability; perturbation convergence plus truncation/roundoff enclosure and positive eigenvalue margin | Cross-term/sign corruption, nonsymmetry, indefinite fixture |
| GC4 local inverse action | Exact rational solution and high-precision independently minimized quadratic energy; H*C identity only C/C | Orientation/solve error; exact rational answers, kappa/backward-error bounds in floating cases | Invert H elementwise, swap a nonconjugate axis, near-singular fixture |
| GC5 x/y equivalence | Independent energy-work response to nontrivial named perturbations in both charts; algebraic T checks C/C | T versus inverse/transpose misuse; exact rational equality and operation-count float bounds | Wrong T, treat fractions as densities, reverse index map |
| GC6 corrected response | Full intrinsic coupled energy + constrained minimization/projection versus neutral consumer, no consumer projector reused | F05 uncorrected response; exact toy equality and independent arbitrary-composition perturbation limit | Omit charge subtraction, wrong charge, omit charged cross-coupling |
| GC7 null/rank/no extra projection | Exact four-species rank/null proof plus typed path/call audit; null identity alone C/C | Illegal full inverse/extra operation; exact singular refusal and zero public projector calls on neutral path | Full inverse/pseudoinverse, repeated-projector call (idempotent output explicitly not a numeric detector) |
| GC8 active embeddings | Independent species phase-space limits, one-sided supported response, no production threshold evaluator reused | Padded H, false mu support, boundary continuation; analytic onset exponent/enclosure and refined threshold partitions | Add absent row, floor density, skip finite refusal window, undeclared phase jump |
| GC9 GR integral | Manufactured variable-lapse exact integral and independent curved-metric exact-star quadrature; separate full background validation | Wrong proper volume/lapse/count conversion; analytic remainder plus independent radial/EOS/center/tail enclosures | Omit e^-nu, extra lapse, omit metric factor, change nu to2nu |
| GC10 global baryon solve | Two-zone exact constrained energy minimization and independent global Schur route; b^T L=0 alone C/C | Local versus global constraint, unsupported mode, poor conditioning; exact toy Z, eigenvalue uncertainty and inverse perturbation bounds | Local Schur-before-integral, delete n_B locally, almost unsupported muon region |
| GC11 source Z | R06 (14)–(18) evaluated from independently established G and signed excess vectors; formula substitution alone C/C | Incorrect coefficients/cross term/sign; exact formula fixture plus independent source-model GC13 | Old F05 scalar formula, wrong minus sign, missing offdiagonal |
| GC12 W coupling | Governed independently validated Phase-5B semantic result + independent physical spin/rate and signed number-excess perturbation; W=ZI alone C/C | Wrong channels/sign/units; propagated Z and I uncertainties plus finite signed action; retains ratified PB caveats | Wrong Omega_dot, K in place of I, second c^-2, omit/swap lepton |
| GC13 corrected published benchmark | Authenticated A18 product and independent Figure1 two-extractor envelope or author arrays; neither fitted to output | Missing correction masked by steady-state agreement; curve/relative-change metric bounded by source extraction, model match and numerical errors | Old F05 formulation must be separable from corrected envelope; steady temperature alone fails gate |
| GC14 provenance/domain/mutations | Deliberate mutations of **every** dependency, including same-version different profile and sequence source; independent expected-refusal table | Stale numbers, domain mixing, hidden state reuse; exact refusal/no-number-on-failure, all declared cases counted by detector type | Stale bytes/revision, foreign metric, core with whole-star I, dangling lifetime token, unsupported phase, mismatched channel labels |

Do not credit PB9's baryon identity, G symmetry by symmetric summation, a Z transpose no-op,
or W=ZI checked against the same multiplication as independent physics evidence. Phase-5B
integration retains all nine ratified qualifications (PB9 tautology, mixed mutation types,
conservative PB9/PB10 budgets, B_B ratio not condition number, PB8 floor, metadata-only true
jumps, conservative PB13 enclosure, adapter EOS/tail authority, PB6 refinement-only evidence)
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:100`). Its strongest independent
structural evidence remains PB10 reconstruction and PB11 nonlinear closure; this task does not
re-credit or rerun that evidence as chemical validation.

## 18. Explicit mutation plan — eighteen cases

Future execution must report each attempted mutation, actual failure, and independence class.
No aggregate “all mutations fired” claim is made here.

| ID | Mutation | Specific detector / current scope |
|---|---|---|
| M1 | Omit electrostatic correction | Exact toy chi*q nonzero and projected-action mismatch; scratch demonstrated. |
| M2 | Second projection | Neutral-path call audit rejects any such call; rederived subtraction zero denominator refuses. Identical P twice is numerically invisible, scratch exact no-op demonstrated. |
| M3 | Full singular inverse or hidden pseudoinverse | Rank/null precondition and unsupported API/type; exact zero-pivot scratch refusal. |
| M4 | Local baryon constraint | Exact two-zone Z, plus free-gas .6998% norm discrepancy. |
| M5 | Flip beta imbalance sign | Named excess vector requires eta=-(21,29)/53 in the two-zone fixture. |
| M6 | Omit inverse lapse | Manufactured exact integral and independent GR oracle; scratch Z changes22.1%. |
| M7 | Extra inverse lapse / use2nu | Same exact oracle and source factor count; extra inverse lapse scratch changes18.1%. |
| M8 | Omit proper volume | Curved-metric oracle and independent formula; scratch changes1.695%. |
| M9 | Flip Omega_dot | Supplied signed physical spin-down rate must yield positive drive for this fixture; compare source excess-number equation, not abs(rate). |
| M10 | Use K as I without c^-2 | Independent SI spin perturbation and unit ledger catches factor c^2. |
| M11 | Apply c^-2 twice | Same physical number-rate action catches extra1/c^2. |
| M12 | Omit electron channel | Separate named electron-only and muon-only forcing with nonzero cross entries; full W comparison. |
| M13 | Swap e/mu | Unequal diagonal entries and unequal named I inputs; signed channel-specific forcing, not norm alone. |
| M14 | Use old F05 Z | Intrinsic toy wrong-route comparison and eventual corrected A18 Figure1 gate; freegas alone not published validation. |
| M15 | Transpose chemical matrix | **Z=Z^T, so pure numerical transpose is an equivalent mutation, not detectable.** Audit row/output and column/input label schema and nonsymmetric rectangular lift/source-map action. Reject semantically reversed types; test nonsymmetric mapping harness separately without calling that harness a physical Z. |
| M16 | Core chemical response with whole-star structural driver | Exact domain/reservoir/boundary-policy mismatch refuses before multiplication. |
| M17 | Stale structural/provider provenance | Change each central/sequence profile, EOS bytes, provider revision, geometry normalization individually; RequireCurrent must refuse. |
| M18 | Invent absent-species support / cross ungoverned interface | Active variant/rank/support metadata and one-sided analytic oracle refuse; no density floor or phase-by-steepness. |

The five deliberately wrong routes in the request are thus either numerically falsified or
explicitly assigned structural refusal/call-path detectors where output mathematics cannot
detect an equivalent no-op. This qualification prevents the same logical overclaim seen in
older validation counts.

## 19. INV-11 and Phase-6 boundary

Coefficient-level facts can be fixed without choosing the evolved state: named beta
combinations and sign, eta^infinity=e^nu eta_local in MeV, spatial constancy under diffusion,
one inverse lapse in G, global baryon reduction, Z action, physical W units and source sign.
These source-level statements are proposed contract semantics, not INV-11 closure.

INV-11 still needs owner-ratified evolved-channel ordering and named/indexed representation,
redshifted storage and units/conversion boundary, reaction stoichiometric mapping and net-rate
sign, thermal-energy accounting, neutrino versus retained-heat partition, dependence on
changing coefficients/background, and solver/driver coupling. No ChemState slot or final
storage layout is chosen (`docs/SCIENTIFIC_INVARIANTS.md:946`). Weak rates, heating/cooling,
chemical/thermal evolution and BNV remain outside this task.

Phase6 architecture note only: future externally supplied particle-number sources must be
separable from structural/equilibrium spin driving and may feed the same chemical-response
machinery. Do not hard-code “chemical source = spin-down only.” Retain the unreduced neutral
G response and explicit conservation assumptions: a source changing total baryon number
cannot be forced into the two-column fixed-baryon lift L. Its later mapping requires its own
conservation/source contract. This preserves the ability to investigate net heating **or**
net cooling without assuming either. No BNV formula, rate, energy partition or implementation
is added.

## 20. Owner decisions and implementation gates

ADR-0013 presents **eight** owner decisions, all pending:

1. Q1: canonical global reduced source basis and G_y object, with full corrected inverse forbidden.
2. Q2: global baryon reduction after integration, supported mode/rank policy.
3. Q3: one named symmetric chemical Z matrix with scalar paper views, owned by the coefficient layer.
4. Q4: separate immutable W using the governed structural semantic object and its full provenance.
5. Q5: explicit chemical reservoirs/domains, active embeddings, tail/onset/interface refusal policy.
6. Q6: A18 closing benchmark and correction-sensitive source/extraction gate; no EOS substitution.
7. Q7: source redshift semantics now while evolved storage/INV-11 remain later decisions.
8. Q8: dependency-complete provenance and error-informed conditioning/refusal contract.

Remaining implementation gates: independent scientific review, owner ratification and a
separate bounded implementation authorization; predeclared numerical algorithms/budgets with
threshold/refusal-gap and profile/provider matching evidence. Realistic Track-R closure further
requires authenticated APR/off-equilibrium, crust/phase/core and benchmark data. Figure2 high
mass requires its own source adjudication. These do not prevent the coherent proposed bulk
contract. No convenient threshold, core cutoff, realistic EOS replacement or state-storage
choice was used to conceal a gap.

## 21. Scratch reproduction, audit trail, and closure

Read-only source extraction used pdftotext; pdftoppm rendered R06 PDF2–4 at scale1800.
Numerical commands, with `worktree` set to the absolute branch path and `scratch` set to the
absolute scratch root above:

```sh
cmake -S "$worktree" -B "$scratch/build" -DPython3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3 -DCMAKE_BUILD_TYPE=Debug
cmake --build "$scratch/build" --target CompactStar -j1
/usr/bin/clang++ -std=c++17 -g -I"$scratch/build/generated/include" -I"$worktree" -isystem /opt/local/include -isystem "$worktree/dependencies/include" "$scratch/probe.cpp" "$scratch/build/libCompactStar.a" "$worktree/dependencies/lib/Zaki/Darwin/arm64/libZaki.a" "$worktree/dependencies/lib/Confind/Darwin/arm64/libConfind.a" /opt/local/lib/libgsl.dylib /opt/local/lib/libgslcblas.dylib /Users/keeper/miniforge3/lib/libpython3.12.dylib /opt/local/lib/libomp.dylib -lz -Wl,-rpath,/Users/keeper/miniforge3/lib -o "$scratch/probe"
"$scratch/probe" "$scratch/run-8192-80000" 8192 80000
"$scratch/probe" "$scratch/run-8192-40000" 8192 40000
"$scratch/probe" "$scratch/run-8192-20000" 8192 20000
"$scratch/probe" "$scratch/run-4096-80000" 4096 80000
"$scratch/probe" "$scratch/run-16384-80000" 16384 80000
/usr/bin/clang++ -std=c++17 -I"$scratch/build/generated/include" -I"$worktree" -isystem /opt/local/include -isystem "$worktree/dependencies/include" "$scratch/local-probe.cpp" "$scratch/build/libCompactStar.a" "$worktree/dependencies/lib/Zaki/Darwin/arm64/libZaki.a" -o "$scratch/local-probe"
"$scratch/local-probe" > "$scratch/local-probe.txt"
/Users/keeper/miniforge3/bin/python3 "$scratch/analyze.py" > "$scratch/analysis.log"
/Users/keeper/miniforge3/bin/python3 "$scratch/oracle.py" > "$scratch/oracle.log"
```

Probes require fresh output directories. To rerun, copy the authenticated scratch scripts to a
fresh root and change their root constant explicitly; do not overwrite evidence. NumPy2.3.1,
SciPy1.16.0 and Python Fraction were used. No SymPy was installed. Initial manual linking
failed on compress, then the executable failed on the Python dynamic-library rpath; adding
`-lz` linkage and the explicit rpath resolved both. These were scratch build failures, not
scientific passes. A manifest check initially used the library directory instead of its
parent and failed to open its prefixed paths; rerunning from `/Users/keeper/Documents/CompactStar`
verified22/22. Initial live Git authentication hit sandbox DNS restrictions and succeeded with
a reviewed network escalation. No literature or repo source was altered to resolve these.

| Scratch artifact | SHA-256 |
|---|---|
| probe.cpp | `e3e04f02e8977f7d1aac506e3be0f511f6f62b3635482165b174f22673400d6a` |
| local-probe.cpp | `5d4a26b146cab69d2bc527f7dbdffaf53641c35f12964d339a3a4c96cd9a93ce` |
| analyze.py | `de132df305475f966144f876d680b890beb11252beef3e41b992f4ae9e89553a` |
| oracle.py | `f78a070cbdb761a63252b19f536037a65b49b07ea31b649be5c20ab9351eaca2` |
| numeric.json | `e4bb42934193c379d1a3bd9006d73c49b5ec22ccf4268406850711c6dd6f07e5` |
| oracle.json | `0c7c3ad8447adadf702e813d6b5e6410d1ebe76bc78eeb7ea374926b9aaaca18` |
| toy.txt | `003cd55cad85de643acd533cfc3de3a50e1315b171c473cb31a6853a0589f7f4` |

No production result changes. G/Z/W diagnostics are not golden targets. No production
Btilde, Z, W, eta evolution, weak rates, heating/cooling, or BNV was implemented. No baseline,
data, literature or test implementation changed. ADR-0013 is PROPOSED; INV-09 stays
VERIFIED / RESOLVED in its governed scope; INV-11 stays UNRESOLVED; BNV is NOT BEGUN.
Commit/push verification is reported after the commit; this document does not invent its own
commit SHA. No merge or owner ratification is performed.
