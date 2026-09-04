# Phase 5A-1A — Independent adjudication of ADR-0010


**Date:** 2026-09-04

**Status:** INDEPENDENT REVIEW / PROPOSED OWNER DECISION; not ratification

**Reviewed HEAD:** `6c8b3686cf54b81b19c220bd5d0299f6cd2bb7d0`

**Branch:** `analysis/phase5a-eos-thermodynamic-audit`

**Worktree:** `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-eos-audit`

**Change class:** documentation; scientific and architecture recommendations only.

**Scientific impact:** no change to any number produced by the code.

## 1. Decision proposed to the owner

**B. ADR-0010 is ready after the specific revisions in §7.**

The central construction is correct: the charge-neutral Hessian in
`x=(n_B,n_e,n_mu)` is a complete **local linear-response primitive for the reduced
2006 bulk npe-mu system**, on its smooth, nonsingular domain. No charged-direction
susceptibility or electrostatic state variable is missing from that response.
It is not a complete specification of individual charged intrinsic chemical potentials,
the electrostatic potential profile, or phase-interface response. The ADR needs to make
those distinctions explicit rather than place a second electrostatic projection after an
already constrained Hessian. The derivation and limits are independently established in §4.

| Question | Disposition | Reason and evidence |
|---|---|---|
| Q1 — Track R closing benchmark | **ACCEPT RECOMMENDATION** | A18 + delta-v + UIX* is the 2005 reference and the EOS used in the 2006 corrected coefficient/evolution figures. It is the appropriate closing target for the requested reproduction. BPAL can reproduce a labelled family member but cannot silently replace that target; CMF is a separate application. [S05 §3.1; S06 §5, Figs. 1–2; §5 below.] |
| Q2 — canonical local basis | **ACCEPT RECOMMENDATION** | `(n_B,n_e,n_mu)` is complete and has direct imbalance conjugates. Its source-equivalent `(n_n,n_e,n_mu)` transformation is exact. Choosing the former is an architecture judgment supported by the existing density axis, not a uniquely required law of physics. [S06 §4; §4.1–4.2; `CompactStar/Core/src/TOVSolver.cpp:628-639`.] |
| Q3 — derivative authority | **REVISE RECOMMENDATION** | Retain precisely the proposed Hessian. Revise its sufficiency claim, intrinsic-potential requirements, domain qualifications, and adapter wording: inversion of the constrained Hessian already includes the 2006 correction. No second projection or full charged extension is required. [S06 eqs. (10)–(18); §4; revisions R2–R5.] |
| Q4 — lepton boundary | **ACCEPT RECOMMENDATION** | A unified thermodynamic provider may compose source-compatible free leptons and a hadronic EOS. One total potential, component provenance, and consistent cross derivatives matter; a single concrete model class does not. Free leptons must not be added twice or replace model interactions silently. [S05 §3.1; §4.7.] |
| Q5 — Track P data authority | **REVISE RECOMMENDATION** | Make source matching mandatory for a claim of canonical DS(CMF) identity, not merely preferable. Family naming alone is insufficient. State the npe-mu species-domain gate and allow a matched hadronic composition surface plus validated leptons, rather than require a monolithic three-density table. [Phase-5A0 audit `:199-243`, `:534-547`; §4.7; revision R6.] |
| Q6 — temperature scope | **ACCEPT RECOMMENDATION** | Cold response v1 matches the benchmark's degenerate thermodynamics. Later thermal evolution still has nonzero temperature in heat capacities and rates. Finite-T response needs a new potential and held-fixed convention. [S05 §3.1, §3.4, eq. (50); §4.8.] |

No seventh owner preference is required to settle these recommendations for the stated goal.
The papers do not uniquely choose Q1's project acceptance target, Q2's chart, Q3's derivative
orientation, Q4's software boundary, or Q6's scope. Those are reasoned recommendations under the
owner's stated reproduction and architecture goals; they remain the owner's decisions. A choice
to change that goal or replace the canonical EOS would require a separate owner decision.

**Owner acceptance is scientifically defensible after R1–R7 are incorporated.** Accepting the
current text unchanged is not recommended. Missing APR inputs block quantitative reproduction,
not acceptance of this bounded local contract. Phase-interface treatment, expanded species,
INV-09, and INV-11 remain gates on subsequent work, not facts silently settled here.

## 2. Authentication, authority, and review method

Entry path, branch, HEAD, and `git status --porcelain=v1 -uall` matched the request exactly;
the worktree was clean. Seven worktrees existed. Local and live remote heads of this branch
were both `6c8b3686cf54b81b19c220bd5d0299f6cd2bb7d0`; local and remote `master` were
`59526c56122dbbb8c0a8ef808bdf627453d99c3a`. All local and cached remote refs were checked for
commits touching the ADR, companion derivation, and new deliverable: none had a file-touching
commit outside reviewed HEAD's ancestry. Live remote branch tips matched the cached tips.
No other worktree was changed. These are entry observations, not assertions about future refs.

Authority consulted:

- `AGENTS.md:11-27`, `:135-148`: checkout authentication, governing hierarchy, local library.
- `GOVERNANCE.md:14-39`, `:43-55`, `:61-79`: authority, documentation evidence, fail-closed scope.
- ADR-0001's accepted species convention (`docs/adr/ADR-0001-species-profile-semantics.md:3-13`);
  ADR-0008 Q11's measure-complete particle response
  (`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:78`).
- Phase 4's ratified inputs and remaining limits
  (`docs/validation/PHASE4_CLOSEOUT.md:64-88`, `:158-181`, `:201-217`).
- INV-09 and INV-11 remain unresolved
  (`docs/SCIENTIFIC_INVARIANTS.md:712-725`, `:820-836`);
  candidate status remains descriptive, not authority
  (`docs/architecture/CURRENT_ARCHITECTURE.md:192-202`, `:334-344`).
- The complete proposed ADR and its derivation were reviewed; the Phase-5A0 inventory was
  consulted, with focused source checks of the current EOS readers and build list.
  Numerical/validation policy slots in governance are explicitly not binding unwritten
  documents (`GOVERNANCE.md:35-39`).

Repository citations in this review, including replacement locations in §7, refer to the
**reviewed HEAD**, so small companion corrections do not silently move their evidence anchors.
Paper references use the exact local PDF, equation/section, and PDF page in §3. Algebra marked
as independent below is derived here from those equations, not claimed as a quotation.
One agent performed this review. No prior-run memory supplied scientific evidence. PDF text was
extracted outside the worktree; the key 2006 equation page was also rendered and visually read.
No production calculation, matrix inversion routine, or differentiation algorithm was implemented.

## 3. Scientific source ledger

All local paper paths below are relative to the read-only root
`/Users/keeper/Documents/CompactStar/literature/rotochemical/`.
Library roles were read in `../catalog.tsv:2-4,7` and `../README.md:13-16` before derivation.

| ID | Exact filename | Directly checked locations and role |
|---|---|---|
| S05 | `2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | Primary non-superfluid benchmark. Eqs. (7)–(8) on PDF p. 3, continued explanation and eqs. (9)–(13) on p. 4; eqs. (24), (30) on pp. 6–7; §3.1 pp. 8–9; §3.2 p. 10; §3.4/eq. (50) p. 12; eqs. (51)–(58) p. 13; Table 1 p. 30. |
| S06 | `2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | Correcting authority. Eqs. (2), (4), PDF p. 2 / journal p. 569; eqs. (9)–(19), surrounding reduction, and footnote 4, PDF p. 3 / journal p. 570; §5 and Figs. 1–2, PDF p. 4 / journal p. 571. Appendix p. 5 illustrates the separate Cowling approximation. |
| S95 | `1995-Reisenegger-Rotochemical-Heating.pdf` | Supporting conceptual source. §2, paper pp. 4–5 / PDF pp. 14–15; npe coordinate count and energy-minimization eq. (3). Its imbalance sign is opposite to the `eta` used here. |
| S20 | `2020-Yanagi-NS-Therm-Thesis.pdf` | Supporting review only. §4.1.1, thesis pp. 58–60 / PDF pp. 60–62, especially the charge-row identity and singularity discussion before eq. (4.7), PDF p. 61. Its earlier shorthand for redshifted intrinsic potentials does not override S06. |

SHA256 of the local PDFs read:

```text
S05 f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88
S06 a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8
S95 9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff
S20 69590539c275fa679a5521a9c5abedd9fdc58718b1554827786c8f707bc618cc
```

A narrow external bibliographic check located the primary
[APR article, Akmal, Pandharipande & Ravenhall, Phys. Rev. C 58, 1804 (1998)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.58.1804)
and the [S06 author/preprint record](https://arxiv.org/abs/astro-ph/0606322).
These are acquisition/provenance pointers; the local S05/S06 PDFs settle the adjudication.
No external EOS dataset was acquired or authenticated as a reproduction input.

## 4. Independent Q3 check

### 4.1 State count and conjugates

Assume cold, neutrino-transparent bulk matter with only `n,p,e,mu` as independent particle
species. Work within one smooth phase, with the same internal equilibration prescription for
all returned thermodynamic quantities. Mean fields or other fast internal variables must be
relaxed consistently; weak beta equilibrium must **not** be reimposed during a composition
perturbation. S06 eq. (4) imposes local neutrality to the accuracy of its bulk approximation;
it does not claim the exact Maxwell solution has identically zero charge density everywhere.

In species order `(n,p,e,mu)`, set `q=(0,1,-1,-1)^T`, measuring charge in units of the
proton charge. Four species densities minus `q^T n=0` leave **three** coordinates.
`n_B=n_n+n_p` defines a coordinate, not a second imposed constraint. At fixed `n_B` there
are two composition directions; across the local provider's state space there are three.
Global baryon conservation removes one integrated degree of freedom later, as S06 footnote 4
explicitly distinguishes. Active beta equilibrium selects a one-dimensional curve.

Define the constant tangent maps

```text
y = (n_n,n_e,n_mu)^T = T x,       x = (n_B,n_e,n_mu)^T,

T = [[1,-1,-1],                 U = T^{-1} = [[1,1,1],
     [0, 1, 0],                              [0,1,0],
     [0, 0, 1]],                             [0,0,1]],

n = S y,                       S_x = S T,
S = [[1,0,0],                  S_x = [[1,-1,-1],
     [0,1,1],                         [0, 1, 1],
     [0,1,0],                         [0, 1, 0],
     [0,0,1]].                        [0, 0, 1]].
```

Both maps have rank three; `S^T q=S_x^T q=0`, and `det(T)=1`.
Pull back `d epsilon=mu^T dn` to the neutral manifold:

```text
g_y = S^T mu = (mu_n, mu_p+mu_e, mu_p+mu_mu)^T,
g_x = T^T g_y = (mu_n, -eta_npe, -eta_npmu)^T,
d epsilon_CN = mu_n dn_B - eta_npe dn_e - eta_npmu dn_mu.
```

These are exactly the local versions of S06's three reduced conjugates between eqs. (13)
and (14). The minus signs follow from neutron destruction in the two beta channels.
The first-law argument does not require an off-neutral EOS to be publicly evaluable.

### 4.2 Symmetry and coordinate transformation

On a `C^2` neutral energy surface, at `T=0` with the other two independent densities held
fixed, define `H_x=partial g_x/partial x`. Equality of mixed derivatives requires

```text
H_Be = partial mu_n/partial n_e = -partial eta_npe/partial n_B = H_eB,
H_Bmu = partial mu_n/partial n_mu = -partial eta_npmu/partial n_B = H_muB,
H_emu = -partial eta_npe/partial n_mu = -partial eta_npmu/partial n_e = H_mue.
```

Symmetry belongs to these **conjugate** coordinates. An arbitrary array of species-potential
partials in some constrained chart is not automatically symmetric. Smoothness alone does not
prove positive definiteness or invertibility; those require a nonsingular physical branch,
with strict stability a sufficient condition. Critical points, coexistence, and active-species
boundaries must be reported. At `n_mu=0`, for example, the cold free-muon chemical derivative
is not a finite smooth continuation of the interior Hessian.

Because `T` is constant, the exact transformations are

```text
H_x = T^T H_y T,                 H_y = U^T H_x U,
C_x = H_x^{-1},                 C_y = T C_x T^T = H_y^{-1}.
```

`H` has units `MeV fm^3`; `C` has units `fm^-3 MeV^-1`.
For nonlinear fraction charts `x=x(z)`, the energy Hessian instead transforms as
`H_z=J^T H_x J + sum_a g_x,a partial^2 x_a/partial z partial z`.
A congruence-only transformation off equilibrium is generally wrong; V8 must account for
this term if fraction-chart energy Hessians are compared. This is an exact chain-rule
identity, not a selected derivative algorithm.

### 4.3 Electrostatic elimination, including a printed redshift omission

S06 eq. (10), with its separately justified Cowling assumption `delta Phi=0`, gives

```text
u(r) = e^{-Phi(r)} delta mu^infinity,
delta mu = u - q delta psi,
delta n = chi (u - q delta psi),        chi = partial n/partial mu.
```

Here `q delta psi` denotes electrostatic energy in the chosen charge units. Where the
unconstrained map exists and `q^T chi q` is nonzero, local neutrality fixes its multiplier:

```text
delta psi = (q^T chi u)/(q^T chi q),
delta n = C_CN u,
C_CN = chi - chi q (q^T chi q)^{-1} q^T chi.                 (A)
```

**Source note:** the journal PDF visibly prints S06 eq. (11) with
`delta mu_m^infinity` in its numerator and no `e^{-Phi}`. Substitution into the line
immediately preceding that equation, or into eq. (10), requires that factor. Equivalently,
the printed right-hand side would give `e^Phi delta psi`. This review identifies an
**inferred printed omission**, not an authenticated erratum. Equation (13) includes the
correct overall `e^{-Phi}` and agrees with (A). No change to the projected bracket is needed.
The companion now records this caveat explicitly.

The field did not vanish physically: neutrality solved for its effect on species forces.
Multiplying by `S^T` then cancels it from neutral conjugates:
`delta g_y=S^T u`. Dropping `delta psi` from the full four-species response before this
elimination would give `chi u`, generally violating local neutrality. Working directly
with the correctly constrained potential performs the elimination by construction.

### 4.4 Proof of exact sufficiency, without assuming the answer

For comparison to (A), suppose a differentiable full intrinsic extension exists with
`K=chi^{-1}`. The projected equations imply

```text
delta n = S delta y,
K S delta y = u - q delta psi,
S^T K S delta y = S^T u,
H_y = S^T K S.
```

If `H_y` is nonsingular, the last equation has a unique physical solution,
`delta y=H_y^{-1} S^T u`. Therefore

```text
C_CN = S H_y^{-1} S^T = S_x H_x^{-1} S_x^T.                (B)
```

This proves equality to the corrected S06 response, not merely matching dimension or a
plausible analogy. Conversely, solving the constrained equation leaves a residual
`u-K S delta y` orthogonal to every neutral tangent; that residual is proportional to `q`
and is precisely the electrostatic multiplier. Thus no additional bulk density response
can be hidden in that direction. A provider can define `epsilon_CN`, `g_x`, and `H_x`
directly without ever constructing `K` or `chi`; the full extension is used only for this
source-equivalence proof.

Deleting proton row/column in (B) selects the identity rows of `S`, so the resulting
principal submatrix is exactly `H_y^{-1}`. It is **not** the neutron/electron/muon
principal submatrix of the unprojected `chi`. This distinction implements the paper's
reduction correctly.

### 4.5 Global singularity and legitimate reduced inversion

For the fixed nonrotating core and frame used by S06 eqs. (12)–(13),

```text
B_full = integral_core e^{-Phi} C_CN dV,
Btilde = integral_core e^{-Phi} C_y dV,
B_full = S Btilde S^T,
B_full q = 0,                    B_pj = B_ej + B_muj.
```

The same constant charge vector is null at every radius, so integration cannot restore the
missing rank. Its conjugate displacement changes only the electrostatic zero/external
potential and no density distribution (S06 text after eq. (13)). A numerical nonzero fourth
eigenvalue would not establish a physical fourth mode. The full corrected matrix has no
ordinary inverse.

On a regular strictly stable three-species-coordinate domain, `C_y` is positive definite.
With positive `e^{-Phi}dV`, for every nonzero constant vector `v`,
`v^T Btilde v = integral e^{-Phi} v^T C_y v dV > 0`; hence the integrated reduced response
is invertible. If species are absent in part of the core, their domain-specific contribution
must be embedded with explicit active coordinates; global rank requires supported response
for every retained mode. If no muons exist anywhere, a three-mode global inversion is not
legitimate. This cannot be repaired by an invented tiny muon density.

Only **after** that integral is inverted does the global condition
`delta N_n=-delta N_e-delta N_mu` give S06 eqs. (14)–(18). In particular,
`eta_npl^infinity=(e_n-e_l)^T Btilde^{-1} delta N_y` gives the paper's signs and
`Z_npe=Btilde^{-1}_nn-2 Btilde^{-1}_ne+Btilde^{-1}_ee`, with the corresponding muon
and cross combinations. Neither imposing `delta n_B=0` pointwise nor integrating `H`
is equivalent. These are identities establishing contract sufficiency, not coefficients
computed or implemented in this task.

### 4.6 What is, and is not, determined

The derivative primitive is sufficient for **all local bulk density response needed by the
corrected reduced coefficient construction**. The background metric/volume, equilibrium
profile, phase-domain information, and later structural/rate inputs are separate required
inputs, already assigned to other layers. `H` alone is not the energy potential or its
first derivative: the contract still needs `g_x`, an equilibrium anchor, and an energy
potential available for validation.

There is one important information limit. On the neutral manifold, `mu` and `mu+q f(x)`
have the same `g_x` for any scalar function `f`. To see why neutral data cannot resolve
this, two full energy extensions differing by `(q^T n) f(n)` have identical restricted
energy, tangent gradient, and tangent Hessian, but different individual charged intrinsic
potentials on that manifold. This is an ambiguity of the off-neutral extension; it does
**not** license changing the intrinsic chemical potentials of a fixed microscopic model
by an arbitrary electromagnetic gauge transformation.

Consequently, neutral `epsilon_CN/g_x/H_x` cannot by themselves recover the individual
`mu_p,mu_e,mu_mu` split, or the source's model-dependent `delta psi(r)`. Those are unnecessary
for (B), `Btilde`, the beta imbalances, and corrected `Z/W` construction. The electrostatic
field would be a separate physical prediction; it must not be advertised as recovered from
`H`. If explicitly required later, one model-defined charged intrinsic potential function
on the neutral chart (including its tangent derivative for linear response), together with
an electrostatic reference/boundary condition, is sufficient extra information to reconstruct
the split and multiplier. For example,
`q_p delta psi=e^{-Phi} delta mu_p^infinity-delta mu_p(x)`.
A public four-dimensional charged EOS is still not necessary for that limited extension.

The minimal response contract is therefore the neutral state/domain, `g_x`, `H_x`, equilibrium
anchor, energy potential for validation, units and provenance. Individual intrinsic species
potentials can be model-supplied capabilities, checked against `S_x^T mu=g_x`; they must
not be claimed as deductions from a neutral energy surface alone. The existing proposal's
mandatory four-potential output is a stronger provider requirement, not a consequence of Q3.

### 4.7 Q4/Q5: composable leptons and the actual data dimension

A source that assumes independent cold leptons can supply

```text
epsilon_CN(x) = epsilon_h(n_n=n_B-n_e-n_mu, n_p=n_e+n_mu)
                + epsilon_e(n_e) + epsilon_mu(n_mu).
```

For a hadronic Hessian in `(n_n,n_p)`, the baryonic tangent map is
`L=[[1,-1,-1],[0,1,1]]`. The total Hessian is
`H_x=L^T H_h L + diag(0, dmu_e/dn_e, dmu_mu/dn_mu)`.
Thus a **two-dimensional hadronic** composition product plus independently variable analytic
leptons can span the full three-dimensional neutral state. Lack of an explicit lepton-split
axis in a hadronic table is not, by itself, a no-go. S05 §3.1 expressly identifies energy
versus baryon density and proton fraction as adequate interacting-particle input.

For such a separable source, `dmu_l/dn_l=(hbar c)^2(3 pi^2 n_l)^(2/3)/(3 mu_l n_l)`
with `mu_l=sqrt((m_l c^2)^2+(hbar c)^2(3 pi^2 n_l)^(2/3))`. This analytic check is valid
at positive density with the declared rest-mass convention. It is not a prescription to
replace a model's interacting leptons. Hadronic terms generate nonzero constrained cross
entries even if the free-lepton sectors have no direct interaction.

A table already containing electrons can only be decomposed/recomposed if its electron
contribution and other equilibrated species are authenticated; subtracting a guessed free
component would change the EOS. The current distinct DNS finite-T product is **not thereby
certified** as a DS(CMF)-1 replacement (`PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:227-243`).
Its reader currently stores `Q2` only (`CompactStar/EOS/src/CompOSE_Thermo.cpp:323-355`).
The current barotropic star reader likewise cannot provide the missing surface
(`CompactStar/Core/src/TOVSolver.cpp:628-639`, `:1551-1583`).

For canonical Track P, source matching means the same model/parameter revision, particle
content, mass convention, phase treatment, and equilibrium background within subsequently
validated source precision. A matching family name, or even agreement on one equilibrium
curve, does not uniquely determine transverse derivatives. The cold DS(CMF) file contains
extra-species columns (`PHASE5A0_EOS_THERMODYNAMIC_AUDIT.md:204-213`): an npe-mu-only contract
cannot cover regions with additional responding hyperon/delta/quark species by omitting them.
Restrict applicability explicitly to a verified npe-mu domain; any extension or elimination
of extra degrees of freedom requires a separately governed equilibration/response contract.

### 4.8 Q6 and nonsmooth-domain limits

Cold v1 fixes thermodynamic response at zero temperature, as appropriate to S05's degenerate
EOS treatment; it does not assert a zero-temperature thermal evolution. S05 eq. (50) and
§3.2 separately use finite-temperature heat capacities and weak rates. At finite temperature
an isothermal Hessian derives from Helmholtz free energy, whereas other constraints require
other specified thermodynamic derivatives. A temperature argument silently appended to the
cold energy Hessian would not establish either contract.

A regular `H` is not a phase-transition or moving-interface response. S05 §3.1 chooses a
Maxwell join for A18 + delta-v + UIX*, with a 6.6% energy-density jump; §2.2 and its footnote
4 also account for composition discontinuities in the spin drive. Smooth-phase local
sufficiency does not justify integrating across such a jump as though every response were
an ordinary function. Source-matched phase information, domain treatment, and any applicable
interface contributions must be established before an APR global benchmark. This is a
later gate, consistent with ADR-0008 Q11, not an invitation to select an interface or
numerical method here.

### 4.9 Independent exact arithmetic check

A scratch-only rational calculation used the analytic fixture

```text
K = diag(2,3,5,7),
H_y = [[2,0,0],[0,8,3],[0,3,10]],
H_x = [[2,-2,-2],[-2,10,5],[-2,5,12]],
C_y = [[1/2,0,0],[0,10/71,-3/71],[0,-3/71,8/71]].
```

`K` values are arbitrary positive test units, not EOS data. Seven exact checks passed:
neutral tangent, both basis maps, `H_y C_y=I`, `H_x(U C_y U^T)=I`, equality of (A) and
`S C_y S^T`, charge null, and proton-row identity. The calculation used rational addition
and multiplication of these analytically supplied matrices and the analytic diagonal `chi`.
It implemented no differentiation or matrix inversion routine. The general proof in §4.4,
not this single fixture, establishes sufficiency.

## 5. Track R ladder and honest reproduction claims

Accept the proposed order, with the following interpretation:

1. **Analytic leptons and toy neutral EOSs:** exact local checks of units, conjugates,
   signs, mixed derivatives, charge projection, and basis transforms. These are validation
   fixtures, not literature reproduction.
2. **S05 free-gas source case:** tests the connection from local response to stellar
   structure. S05 §3.1 uses free gas for the **whole star**, not a free core with the APR
   crust. Table 1's free-gas endpoint is limited by Sigma-minus appearance; it is not an
   unconstrained extension of an npe-mu sequence. Its listed Kepler period is adopted
   from a pure-neutron-gas calculation, as §3.1 explains.
3. **Optional authenticated BPAL variant:** useful interacting intermediate if the exact
   parameters are available at modest cost. BPAL11, 21, 31, 32, and 33 have distinct source
   labels (S05 §3.1/Table 1); reproducing one does not require acquiring or implementing all
   five. Skipping BPAL is acceptable; making it mandatory would add no scientific necessity.
4. **A18 + delta-v + UIX* closing benchmark:** reproduce selected S05 reference outputs
   with S06 controlling the corrected response. Build the correction into physical response
   checks from the start. An uncorrected S05 branch, if later used for comparison, must be
   labelled historical, not temporarily accepted physical authority.
5. **DS(CMF) production after Track R closure:** use the validated generic contract with
   an authenticated matched provider. Successful CMF curves do not retroactively establish
   the literature benchmark.

A quasi-steady temperature alone is a weak discriminator: S06 §3 shows it is unchanged by
the correction. Closure must include a correction-sensitive coefficient or transient
comparison as well as the final thermal benchmark. S06 §5/Fig. 1 reports nearly unchanged
`Z_np` and changes reaching roughly 10–20% in `Z_npe/Z_npmu` at high mass; those qualitative
statements are useful checks but not precise numerical fixtures. S06 Fig. 2 supplies explicit
1.47 and 2.13 solar-mass cases, initial `T^infinity=10^8 K`, zero imbalances, `B=10^8 G`, and
`P_0=1 ms`. Those are candidate source targets, not newly chosen numerical acceptance tests.

A faithful reproduction also labels historical domain limitations. S05 §3.1 states the
reference becomes acausal above approximately `2 x 10^15 g cm^-3`, corresponding to a
2.14 solar-mass model, below its formal maximum mass. Reproducing a historical acausal
point would not validate it for production. Do not silently causalize or smooth APR and
still call it the same benchmark.

## 6. Missing acquisition and authentication gates

The local catalog has S05/S06 and a broad Prakash review, but no dedicated APR, exact BPAL,
crust, or effective-mass primary package listed below. Neither a paper title nor the current
CMF products supply those missing inputs. This review found bibliographic pointers, not an
authenticated APR implementation lineage.

| Needed before the relevant claim | Required material and what must be established | Source basis |
|---|---|---|
| Quantitative APR thermodynamics | Primary APR 1998 definition and the exact A18 + delta-v + UIX* analytic fit/parameters or equivalent composition-resolved product used by S05/S06; energy normalization, neutron/proton masses, composition interpolation, density validity, branches, and equilibrium/transverse response. A beta-equilibrium `P(epsilon)` table or an unspecified modern “APR” fit is insufficient. | S05 §3.1; [APR primary record](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.58.1804). |
| Quantitative APR star/drive | Pethick, Ravenhall & Lorenz, Nucl. Phys. A 584, 675 (1995), and Haensel & Pichon, A&A 283, 313 (1994); actual crust data/definitions, joins, core boundary, and S05 Maxwell phase prescription. Authenticate numerical realization rather than infer it from a label. | S05 §3.1; §2.2 eq. (30) and footnote 4. S05 spells the former author's surname “Lorentz” in its citation. |
| Optional BPAL comparison | Prakash, Ainsworth & Lattimer, Phys. Rev. Lett. 61, 2518 (1988), and the Prakash et al. Physics Reports 280, 1 (1997) labels/parameters for the **chosen** variant. Optional if BPAL is skipped. | S05 §3.1 and references. |
| Quantitative heat-capacity/rate/evolution reproduction | Source-consistent effective masses (S05 points to Page et al. 2004), Urca rates and reduction functions, source envelope prescription, spin law, initial conditions, and all normalizations used in the selected figure. These cannot be inferred from `H`. | S05 §§3.2–3.6, especially eq. (50); S06 Fig. 2. |
| Honest numerical comparison | Existing tables plus sufficient independent reference values for the selected composition/structure, drive, corrected coefficients, and evolution outputs. Prefer author tables/code; if unavailable, declared figure digitization with quantified uncertainty and approved acceptance tolerances is a legitimate alternative. Table 1 supplies some published structure numbers, so “no numerical references exist” would be false. | S05 Table 1 and figures; S06 Figs. 1–2. |
| Exact free-gas source comparison | Source conventions and the Shapiro & Teukolsky formulation cited in S05, composition/endpoint rule, and whole-star setup. An analytic gas fixture may be validated before claiming exact source replication. | S05 §3.1 and Table 1. |

Author code is useful, not logically mandatory if the published definition and reference
uncertainty suffice. Conversely, acquiring only the APR paper does not prove that an
implementation uses the precise fit, branch joining, and conventions of S05/S06. The eventual
reproduction record must authenticate that chain. No missing APR/BPAL data were invented,
downloaded into the library, or selected as authority here.

## 7. Exact proposed final answers and ADR edits

The following text is proposed for an owner-reviewed ADR revision. **It is not applied to
ADR-0010 by this task.** Its status and Decision section remain PROPOSED / PENDING OWNER
ADJUDICATION. The owner alone can accept the revised ADR in a separate action.

### 7.1 Exact proposed answers to Q1–Q6

**Q1:** “Track R closes on quantitative reproduction of the A18 + delta-v + UIX* reference
calculation, using the 2006 correction wherever it supersedes 2005. Analytic lepton/toy checks
and the source free-gas case precede it; an authenticated BPAL variant is optional. BPAL and
CMF cannot silently replace the closing target. Missing authenticated inputs block Track R
closure. Include a correction-sensitive comparison, not quasi-steady temperature alone.”

**Q2:** “Use `(n_B,n_e,n_mu)` as the canonical local density chart, with exact reconstruction
of dependent species and mandatory response equivalence to `(n_n,n_e,n_mu)` under the stated
linear transformation. The chart applies only on its declared npe-mu phase/species domain;
thresholds do not silently change coordinate meaning.”

**Q3:** “The EOS-facing derivative authority is
`H_ab=partial g_a/partial x_b`, for `x=(n_B,n_e,n_mu)` and
`g=(mu_n,-eta_npe,-eta_npmu)`, at zero temperature with the other independent coordinates
held fixed. Together with the neutral state, conjugates, equilibrium anchor, validation energy,
and domain/provenance, it completely specifies the smooth local reduced response used by
2006. The adapter derives `C_y=T H_x^{-1} T^T`; neutrality and the electrostatic elimination
are already included, with no second projection. Only qualified nonsingular reduced responses
may be solved. Individual charged intrinsic potentials are optional model-supplied capabilities,
not deductions from the neutral energy/Hessian. No full `K`, `chi`, paper-`B` inverse, or
electrostatic-potential output is required.”

**Q4:** “Require one coherent consumer-facing npe-mu thermodynamic provider; permit internal
composition of a hadronic EOS with independently validated free-Fermi leptons when consistent
with the source model. Component ownership, rest masses, interactions, and cross derivatives
must be declared, with no omitted or double-counted lepton contribution.”

**Q5:** “Preserving canonical DS(CMF) identity requires an authenticated source-matched
arbitrary-composition provider/product, including model revision, species, phase, mass and
lepton conventions, and equilibrium-background matching. A matched hadronic composition
surface plus validated leptons is allowed. The current beta-equilibrium table, distinct finite-T
Y_q table, and dormant RMF code are not automatically equivalent. Restrict v1 use to a verified
npe-mu domain; additional responding species require a separately governed extension. If the
matched source is unavailable, stop canonical Track P work for owner adjudication rather than
substitute a different EOS.”

**Q6:** “V1 is explicitly cold-only for the thermodynamic response. Nonzero-temperature rates
and heat capacities remain separate later inputs. Finite-temperature response is deferred to
an explicit contract revision specifying the potential, thermal constraint, and validation.”

### 7.2 Exact companion edits required in ADR-0010

Locations refer to ADR-0010 at the reviewed HEAD. These are the complete substantive edits
recommended by this adjudication, in addition to replacing Q1–Q6's Recommended answer cells
with §7.1. Preserve historical provenance and owner-pending status.

**R1 — source precision and Track R gate.** In Source authority item 2 (`:42`), replace
“eqs. (7)-(13), PDF p. 4” with “eqs. (7)-(8), PDF p. 3, and eqs. (9)-(13), PDF p. 4”.
Replace Track R's Proposed entry and acceptance logic cell (`:105`) with:

> Use analytic lepton/toy checks, the 2005 whole-star free-gas case, optional authenticated
> BPAL, and the A18 + delta-v + UIX* closing benchmark. Authenticate the exact source
> definition, phase/crust prescription, and comparison data before closure. Apply the 2006
> correction from the outset and include a correction-sensitive coefficient or transient
> comparison; quasi-steady temperature alone cannot validate it. A missing closing source
> blocks closure; another EOS does not silently substitute.

**R2 — specify the sufficiency domain.** Append after the Hessian transformation (`:187`):

> This sufficiency statement concerns smooth cold bulk npe-mu response in the 2006
> Cowling/local-neutrality approximation. The internal equilibration prescription must be
> consistent across energy, conjugates, and Hessian; weak equilibrium is not imposed during
> independent composition differentiation. Symmetry requires a twice differentiable energy
> potential, and invertibility requires a nonsingular physical branch; strict stability is
> sufficient, smoothness alone is not. At absent-species boundaries, critical points, and
> phase transitions, report the domain restriction rather than fabricate a full-rank Hessian.
> A regular Hessian does not supply moving-interface response. Source-specific phase and
> species treatment must be governed before a global calculation crossing those boundaries.

In the derivative-authority table (`:167-168`), replace the two Singularity and symmetry
cells respectively with:

> Symmetric in conjugate coordinates where the neutral energy is C^2. Invertibility must be
> established on the declared active-species domain; strict stability is sufficient.

> Symmetric inverse response only where the corresponding reduced Hessian is nonsingular;
> no smooth full-rank continuation through absent-species or phase boundaries is implied.

**R3 — remove ambiguous double-projection ownership.** Replace the Local rotochemical
thermodynamics adapter row (`:198`) with:

> Derive the source-basis reduced susceptibility from the constrained EOS Hessian and verify
> its equivalence to the 2006 charge projection | Local rotochemical thermodynamics adapter |
> `C_y=T H_x^{-1} T^T` already incorporates electrostatic elimination. Do not project it
> again. A full intrinsic toy/provider extension may be used privately for equivalence checks;
> neither that extension nor `delta psi` is required from a neutral provider.

Replace the paragraph at `:203-206` with:

> The provider's mandatory conjugates are local intrinsic neutral combinations. It does not
> add electrostatic or redshift terms. If intrinsic species potentials are supplied, they
> must have a declared model convention and reproduce those combinations. The neutral
> Hessian determines the corrected bulk density response, not the individual charged-potential
> split or the electrostatic potential profile. The electrostatic term cancels after
> enforcing neutrality, not by omitting it from the unconstrained derivation. In using 2006
> eq. (11), retain the `e^{-Phi}` factor required by eq. (10) and present in eq. (13); the
> printed eq. (11) appears to omit it. This is an inferred source omission, not a cited erratum.

In the EOS provider row (`:197`), replace “return intrinsic species chemical potentials and
the selected local derivative primitive” with “return neutral conjugates and the selected
local derivative primitive, with optional model-defined intrinsic species potentials”.

**R4 — make the minimal potential authority self-consistent.** Replace the Intrinsic
chemical potentials and Imbalances rows at `:216-217` with these three rows:

| Capability | Phase-5 classification | Required semantics |
|---|---|---|
| Neutral conjugates `g=(mu_n,-eta_npe,-eta_npmu)` | REQUIRED NOW | Local MeV, one energy/rest-mass convention, consistent with the validation energy and `H`; no electrostatic or GR term. |
| Intrinsic `mu_n,mu_p,mu_e,mu_mu` | OPTIONAL MODEL CAPABILITY | Local MeV with declared microscopic/component convention; if supplied, require `S_x^T mu=g`. The charged split cannot be inferred from neutral energy/H alone. |
| `eta_npe`, `eta_npmu` | REQUIRED NOW, DERIVED | `eta_npe=-g_e`, `eta_npmu=-g_mu`; if species potentials are supplied, their beta combinations must agree. No independently drifting authority. |

**R5 — sharpen falsifiers.** Replace V8, V9, and V10 (`:248-250`) with:

| ID | Required check | Defect it falsifies |
|---|---|---|
| V8 | Exact A/B response equivalence; for valid fraction charts include the nonlinear chain-rule term in energy-Hessian transformations. | Basis dependence, density/fraction mixing, omitted second-coordinate derivatives. |
| V9 | In an independent full-intrinsic toy fixture, compare the 2006 projected response with `S_x H_x^{-1} S_x^T`; then verify full global charge null and proton-row identity. | Missing or double electrostatic projection. A null identity alone is insufficient to validate response magnitudes. |
| V10 | Establish reduced symmetry, rank, stability and conditioning on the declared active domain; embed absent-species responses explicitly and verify global mode support before inversion. | Hidden gauge mode, fabricated absent-species rank, instability, or invalid threshold continuation. |

**R6 — tighten production identity without overconstraining storage.** Replace Track P's
Proposed entry and acceptance logic cell (`:106`) with:

> Canonical DS(CMF) identity requires the same authenticated model/parameter revision,
> species, phases, mass/lepton conventions, and matched equilibrium background. Supply the
> neutral composition response through a matched arbitrary-composition product/provider;
> a two-dimensional hadronic surface plus validated source-compatible leptons is allowed.
> The current equilibrium file supplies the background only. Restrict v1 applicability to
> a verified npe-mu domain; additional responding species require separate governance.

Replace the first two sentences beginning “Track R providers” in Migration implications
(`:316-318`) with:

> Track R providers follow the staged ladder and close on the authenticated A18 + delta-v +
> UIX* benchmark. Canonical Track P requires the source-matched composition provider defined
> above. A distinct analytic/RMF model is a separate owner-governed research EOS and does
> not retain canonical DS(CMF) identity by assertion.

**R7 — record the evidence without exercising owner authority.** Add this record to the
Evidence companion field (`:10`), and append under Provenance:

> Phase 5A-1A independently checked the local 2005/2006 papers and proved equivalence of
> the neutral Hessian to the corrected reduced bulk response. Its proposed revisions and
> exact Q1–Q6 answers are recorded in
> `docs/validation/PHASE5A1A_ADR0010_INDEPENDENT_ADJUDICATION.md`. This review does not
> constitute owner acceptance or implementation authorization.

## 8. Scope and validation closure

The only edits in this task are this adjudication record and narrow factual corrections in
`PHASE5A1_THERMODYNAMIC_CONTRACT_DERIVATION.md`: the PDF page for S05 eqs. (7)–(8), the
inferred S06 eq. (11) redshift omission, and the invertibility/strict-domain qualifications.
ADR-0010 itself is byte-unchanged and remains PROPOSED. Its suggested changes are review text
in §7, not exercised owner authority.

Validation performed: checkout/ancestry/live-tip authentication; direct local paper reading
and visual verification of S06 p. 3; the seven exact rational identities in §4.9; source/link
and document-structure checks; `git diff --check`; scope checks against reviewed HEAD; and
before/after SHA256 comparison of all regular files in the shared literature root. The initial staged whitespace
check reported `trailing whitespace` on the six Markdown hard-break lines in the metadata
block (then lines 3–8); those spaces were removed. The final checks passed.
No production/test/baseline/build file changed. No literature bytes changed. No scientific
runtime tests were run because no runtime behavior was modified, and these algebraic checks
do not claim numerical EOS or end-to-end reproduction validation.

No APR/BPAL data were fabricated. No derivative algorithm, matrix inversion implementation,
rotochemical coefficient, evolution equation, superfluid extension, or BNV work was added.
No next phase begins under this record. Commit/push state is reported separately after the
record is committed, so the document does not claim a self-referential commit hash.
