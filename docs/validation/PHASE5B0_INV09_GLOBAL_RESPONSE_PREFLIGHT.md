# Phase 5B-0: global structural particle-number response preflight

**PHASE 5B-0 PREFLIGHT COMPLETE. ADR-0011 ACCEPTED.**

> **INV-09 IMPLEMENTATION CONTRACT DEFINED —**
> **PHASE-5B PRODUCTION IMPLEMENTATION BLOCKED ON UNIT-BOUNDARY RECONCILIATION.**

This is a source-backed documentation audit with disposable numerical prototypes. The human
owner accepts the structural formula, fixture, domain and sequence-ownership contract in
`docs/adr/ADR-0011-particle-number-structural-response.md:1`. Acceptance does not authorize
production implementation or accept a numerical coefficient. INV-09 remains **INTENDED BUT
UNVERIFIED**; INV-11 remains **UNRESOLVED**. The demonstrated GSL/Zaki unit-convention mismatch
must be reconciled in a separate governed task before Phase-5B production implementation.

NO CORRECTED PAPER BTILDE, PAPER Z/W, EVOLVED CHEMICAL STATE,
WEAK-RATE MODIFICATION, ROTOCHEMICAL EVOLUTION, SUPERFLUIDITY,
APR/BPAL, DS(CMF) OFF-EQUILIBRIUM PHYSICS, OR BNV IS IMPLEMENTED
BY PHASE 5B-0.

## 0. Owner adjudication — accepted contract, implementation blocked

The owner decisions are authoritative and are recorded in full in ADR-0011:

- **Q1:** the primary implementation/validation fixture is the human-ratified Structure-1
  midpoint, `rho_c=1.10e15 g cm^-3`. FR2005 Figure 1 at
  `P_c=4.5e34 dyn cm^-2` is secondary source-comparison material only. The midpoint is not a
  published free-gas `I_Omega` benchmark.
- **Q2:** the fixed-total-baryon constraint and INV-09 free-gas validation use whole-star
  counts, with `N_B=sum_baryons b_i N_i` and `N_B=N_n+N_p` for npe-mu. Generic response objects
  are domain-qualified. A core/fixed-isobar driver requires an explicit boundary and the PN7/PN8
  boundary-flux mapping; no unpublished cutoff is inferred.
- **Q3:** production `B_i` is owned by complete canonical-star differentiation, preferably in
  `x=ln(rho_c/rho_ref)`, using multi-scale symmetric steps and a higher-order estimate. Every
  neighboring star owns its grid, metric, composition and surface. The regular homogeneous
  solution is the mandatory PB7 independent oracle only; no public homogeneous Phase-4 API is
  created.
- **Q4:** a separate governed GSL/Zaki unit-boundary reconciliation is a prerequisite to any
  Phase-5B production implementation. The convention-matched scratch result is diagnostic only;
  PB7 is blocked until one convention is accepted, repaired as necessary and independently
  revalidated.

PB1–PB14 are the accepted validation plan, not achieved validation. Existing numerical budgets
remain proposed implementation acceptance targets. INV-09 cannot close without all PB1–PB14,
independent review and human ratification.

## 1. Starting authority, checkout and change class

Starting SHA: `c1be10a374ef893290583855e8fc4df3a815e66e`.
Canonical checkout: `/Users/keeper/Documents/CompactStar/repo/CompactStar`, branch `master`.
Canonical status was clean; HEAD, local master, origin/master and live
`git ls-remote origin refs/heads/master` were equal to that SHA. The initial restricted-network
query failed DNS; the authorized network query succeeded. No remote reference was changed.
All other local branch tips are contained in master; the all-ref log excluding master had no
unmerged changes to the requested documentation paths.

Fresh branch: `analysis/phase5b-inv09-global-response-preflight`.
Fresh worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5b-inv09-preflight`.
It was created at the exact starting SHA and authenticated clean. One agent; no delegation.

Class: **documentation**, including source derivations and scientific/numerical/architecture
recommendations, not implementation. Authority: `GOVERNANCE.md:39`, `GOVERNANCE.md:65`,
ADR-0001 species semantics (`docs/adr/ADR-0001-species-profile-semantics.md:80`), ADR-0003
provenance (`docs/adr/ADR-0003-profile-cache-provenance-and-invalidation.md:242`), ADR-0005
canonical solve (`docs/adr/ADR-0005-canonical-tov-numerical-primitive.md:208`), ADR-0007
fixed-center response and deferred sequence owner
(`docs/adr/ADR-0007-hartle-second-order-monopole-response.md:418`), ADR-0008 Q11 measure
completeness (`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md:78`), ADR-0009
finite-cutoff surface (`docs/adr/ADR-0009-tov-surface-event-and-termination.md:47`), and
ADR-0010 layer boundaries (`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:273`).
Accepted ADR decisions are untouched. Current response authority is the later ratification,
not stale candidate warnings in some headers: `docs/validation/PHASE4_CLOSEOUT.md:64`.

## 2. Direct source ledger

Library root: `/Users/keeper/Documents/CompactStar/literature`. The six ledger sources were consulted directly;
text was extracted outside the tree and equation/figure pages were visually inspected.
The complete `SHA256SUMS.txt` check passed from the CompactStar parent directory, including
catalog and README. Originals were not changed. Catalog roles give R2006 precedence for
its correction scope. Online bibliographic cross-checks agree with the local papers:
[FR2005](https://arxiv.org/abs/astro-ph/0502116),
[R2006](https://arxiv.org/abs/astro-ph/0606322).

| Key and exact relative PDF | SHA256 | Sections/equations used |
|---|---|---|
| F: `rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88` | pp. 3–8, §§2.1–2.3, (9)–(31); §3.1 pp. 8–9; §3.5 p. 12; §3.6 (52)–(58); Fig. 1 p. 26, Fig. 3 p. 28, Table 1 p. 30; Appendix B interface continuity |
| R: `rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf` | `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8` | §§3–5, pp. 569–571; (6)–(7), (10)–(19), Figs. 1–2 |
| H: `rotation/1967-Hartle-I.pdf` | `ed263946e9bc13842399b5c9e9c2eae31823e7323bc81b456fb5174697cefc35` | pp. 1019–1023; (87)–(100), (105)–(114), especially current integral (109), center (108) and angular reduction preceding (114) |
| HT: `rotation/1968-Hartle-II.pdf` | `4584cb2b299b84ed101fb8543abe9352be915ba9318d3d790837af54cbc7583e` | §II, (7), (15)–(16), supporting the accepted Phase-4 variables; not a new numerical target |
| R95: `rotochemical/1995-Reisenegger-Rotochemical-Heating.pdf` | `9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff` | §2 (4)–(5), §3.1: Lagrangian compression and diffusion versus reaction-volume distinction |
| Y: `rotochemical/2020-Yanagi-NS-Therm-Thesis.pdf` | `69590539c275fa679a5521a9c5abedd9fdc58718b1554827786c8f707bc618cc` | §4.1.1, pp. 58–61, (4.6)–(4.10): corrected submatrix treatment and structural/chemical distinction; supporting only |

The remaining rotation-library PDFs were checksum-authenticated; no new result is attributed
to an unread section of them. Accepted Phase-4 derivations/revalidation/baseline/closeout were
consulted with their historical-versus-current status preserved. In particular the number
formula is `docs/validation/PHASE4C_HARTLE2_DERIVATION.md:610`; the measure principle is
`docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md:271`; the accepted convergence
clarification is `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md:586`; the
baseline is regression evidence only (`docs/validation/PHASE4D_MONOPOLE_BASELINE.md:42`).

**Source typography:** F(22) actually prints `(-g)^(-1/2)` on the left, whereas F(20),
the right side of (22), H(109), and F(24) require `sqrt(-g)`. This is an inferred printed
exponent error, not an asserted published erratum. It does not change the independently
derived current measure. R(11)'s inferred missing redshift factor is already governed by
ADR-0010; nothing here uses that printed equation to invert a chemical response.

## 3. Source notation map

Counts are dimensionless; the word count is retained to make dimensions intelligible.
Here `q=Omega_geom^2`, `e_c=epsilon_c,geom`, `x=ln(rho_c/rho_ref)` and
`rho_c=epsilon_c,physical/c^2`. Superscript `D` always specifies a domain.

| Source symbol | Equation | Meaning / held fixed | Units | Local/global | CompactStar analogue / status |
|---|---|---|---|---|---|
| `N_i`, `N_i^eq`, `delta N_i` | F(11),(14),(30); R(6),(12) | Species counts, equilibrium counts, excess counts in the declared diffusion/core domain | count | global over D | New number functional needed; legacy cache invalid |
| `N` | F(20),(21) | Enclosed **baryon** count at a pressure surface, not species count | count | cumulative global | Generic cumulative baryon integral proposed |
| `A` (our `N_B`) | F §2.2, (18) | Whole-star total baryon count, held fixed in spin-down; surface P=0 | count | whole star | Background nB exists; complete response absent |
| `(partial N/partial Omega^2)_(P,rho_c)` | F(19),(24) | Fixed-isobar enclosed baryon rotation derivative at fixed center | geometric count km²; physical count s² | cumulative | `A_B(P)` from governed fields, scratch only |
| `(partial N/partial rho_c)_(P,Omega^2)` | F(26),(27) | Equilibrium-sequence derivative at fixed isobar | count cm³/g in mass-equivalent CGS | cumulative | `B_B(P)`; no public production owner |
| `(partial P/partial Omega^2)_(N,A)` | F(18) | Lagrangian compression at fixed enclosed and total baryon counts | geometric dimensionless; physical pressure s² | local shell, globally constrained | Cumulative reduced response, absent |
| `I_{Omega,i}` | F(30), R(7) | Integral of equilibrium composition change driven by compression, before susceptibility | count s² for physical Omega | global over core | Domain-qualified structural driver proposed; no production result |
| corrected `B_ij`, reduced `Btilde` | R(12),(13), text before (14) | Global chemical susceptibility at background/Cowling structure; charge projection included | count/energy | global | ADR-0010 later Layer B; forbidden here |
| `Z_np`, `Z_npe`, `Z_npmu` | R(16)–(18), superseding F(54)–(56) | Chemical imbalance per excess-particle response after constraints | energy/count | global | Not the historical structural `Z_i`; forbidden here |
| `W_npe`, `W_npmu` | F(57),(58); R(19) | Chemical susceptibility combined with structural driver | energy s² | global | Later consumer, not computed |

Use production names such as `ParticleNumbers`, `FixedCentralEnergyNumberResponse`,
`EquilibriumSequenceNumberDerivative`, `FixedBaryonNumberResponse` and
`CoreCompositionSpinDownDriver`. Fields should say `dN_dOmegaGeom2`, `dN_dCentralEnergyDensity`,
`dCentralEnergyDensity_dOmegaGeom2` and carry the domain. Avoid unqualified `A`, `B`, `Z` in
public APIs. Mathematical `A_i`, `B_i`, `K_i` below are explicitly structural; `K_i` is never
paper Z. Historical FR `A` also collides with structural `A_i` and with later heating constants.

## 4. Benchmark-star determination

F Table 1 is a static maximum-mass-labelled inventory, not an `I_Omega` table. F §3.1 says
the free-gas row lies slightly below Sigma-minus appearance. The qualified common-state
reproduction is ratified; the source's exact retained model/maximum-selection rule is not:
`docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_RATIFICATION.md:34`, `:80`.

F Fig. 1 **does** publish a free-gas compression curve at the separate input
`P_c=4.5e34 dyn/cm²`, with ordinate `-(Omega_K²/P)(partial P/partial Omega²)_(N,A)` and
abscissa `r/r_core`. It is not the Table-1 midpoint. Independent direct free-gas inversion
here gives approximately `rho_c=1.0377003368e15 g/cm³`, `n_B,c=0.5722204221 fm^-3`.
The paper does not give a numerical free-gas `r_core`/core-cut prescription sufficient to
silently define the domain of that plotted curve. The free-gas `P_K=0.98 ms` is adopted
from a pure-neutron calculation, not a new mass-shedding calculation of this npe-mu model
(F §3.1 and Table-1 footnotes). Record that convention if the figure is digitized later.

F Fig. 3's `I_{Omega,e}` and `I_{Omega,mu}` curves are for **A18 + delta-v + UIX*** as a
function of mass, not free gas and not Table 1. R Figs. 1–2 provide chemical/evolution
comparisons for that realistic EOS, not an additional free-gas structural coefficient.
There is no uniquely tabulated free-gas species-driver benchmark in these sources.

**Accepted primary implementation/validation fixture:** the already ratified Structure-1 midpoint,
requested and achieved
`rho_c=1.10e15 g/cm³`, `n_B,c≈0.6049252413 fm^-3`; finest imported `P_c=4.91565287895e34
 dyn/cm²`. The owner selected this source-compatible static fixture without calling it a
published `I_Omega` benchmark. The distinct Figure-1 pressure configuration remains secondary
source-comparison material until its domain/core-radius and `Omega_K` conventions are authenticated.
Neither source mass fitting nor a maximum-mass search was performed.

## 5. Current production and historical candidate audit

Current path: `SingleStarSolveToTOVPoints` -> `NStar`/`StarProfile` -> automatic normalized
first-order response -> explicit `ComputeHartleMonopoleResponse`. The source of n_i is
`Y_i*n_B`. The source of the energy measure is the profile. The Track-R value-only
`BarotropeAt` path is distinct from the local Hessian path and crosses continuous onsets
without weakening N-1. Evidence: `CompactStar/Core/src/TOVSolver.cpp:2511`,
`CompactStar/Core/src/NStar.cpp:176`, `CompactStar/Core/StarProfile.hpp:636`,
`CompactStar/Core/NStar.hpp:411`, `CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:326`,
`CompactStar/EOS/src/LocalThermodynamics.cpp:153`.

| Historical routine/surface | Classification | Authenticated defects / disposition |
|---|---|---|
| `RotochemicalCache::ComputeEnclosedNumber` | FORMULA CORRECT BUT IMPLEMENTATION INVALID in its call path | Trapezoid form and 1e54 factor are conceptually correct; caller passes Y, not n; no array-size/domain/provenance checks; omits unresolved center/surface policy (`CompactStar/Physics/Evolution/src/RotochemicalCache.cpp:25`, `:147`) |
| `ComputeStructuralDerivative` | DANGEROUS / MUST NOT ACTIVATE | Only nodal dn/dp times old p0; drops metric, velocity and surface; no Omega² normalization; absolute dp guard returns fabricated zeros; N<3 returns zero. Uses deleted `HartleResult`; misleading assertion that other terms are subdominant; citation to F Eq. A.14 is not the number formula (`:47`) |
| `Build` sequence finite difference | REUSABLE CONCEPT ONLY | Complete nearby stars are sound inputs; implementation uses central star geometry on other grids, possible out-of-bounds access, raw fractions, no central-step validation or achieved-state/provenance checks (`:107`–`:165`) |
| `Build` reduction | FORMULA CORRECT BUT IMPLEMENTATION INVALID | Algebra correct only for compatible domain/units; near-zero B_B falls back to A instead of refusing; valid can remain stale on early return; missing species skipped/zeroed; no conditioning or conservation checks (`:167`–`:179`) |
| `Rotochemical::AccumulateRHS` | DANGEROUS / MUST NOT ACTIVATE | Number driver written into energy-valued eta slots without susceptibility; species/channel mismatch; physical/geometric conversion absent; depends on driver order and earlier accumulated spin RHS (`CompactStar/Physics/Driver/Chem/src/Rotochemical.cpp:48`, `:102`) |
| old `HartleResult`, shooting/raw-seed path | OBSOLETE | Deleted by Phase 4; dangling candidate references cannot serve as a compatibility contract (`docs/validation/PHASE4_CLOSEOUT.md:33`) |
| `ChemState` pack/unpack/storage | REUSABLE CONCEPT ONLY for eventual state plumbing | Compiled storage, no accepted evolved chemical semantics or source-to-slot mapping; not a structural oracle (`CompactStar/Physics/State/ChemState.hpp:45`, `CompactStar/Physics/State/src/ChemState.cpp:9`) |
| `BNVSource`, `WeakRestoration`, other empty Chem placeholders | OBSOLETE as evidence of operational physics | Header/install scaffolding; no implementation authority (`CompactStar/Physics/Driver/Chem/CMakeLists.txt:1`) |

No legacy routine qualifies as a numerical oracle. Full tracked-code/caller search found
`Build`, `ComputeEnclosedNumber` and `ComputeStructuralDerivative` referenced only in the
candidate itself; the driver only has its own declaration/definition and cache accessors.
`Physics/Evolution/CMakeLists.txt:23` excludes the cache source; `Driver/Chem/CMakeLists.txt:9`
has an empty source list. ChemState is compiled through `Physics/State/CMakeLists.txt:14` and
used by StateVector, but that does not activate chemistry. Current-architecture evidence:
`docs/architecture/CURRENT_ARCHITECTURE.md:218`, `:371`. Nothing was added to CMake.

## 6. Nonrotating number functional and integration domain

For a comoving species current `j_i^mu=n_i u^mu`, its instantaneous count on a t=constant
slice is the oriented hypersurface integral, equivalently
`N_i=integral_D n_i u^t sqrt(-g) d^3x` (H109; F20 for baryons). Species currents can have
reaction sources; this does not change the instantaneous count measure. For the background,
`u^t=e^-nu`, `sqrt(-g)=e^(nu+lambda) r² sin(theta)`, hence

```
C = 10^54 (fm^3/km^3 conversion in the numerical density-times-volume product)
w(r) = 4 pi r^2 / sqrt(1-2m/r)
N_i^D = C integral_D w(r) n_i(r) dr,        n_i = Y_i n_B [fm^-3].                 (PN1)
```

There is no lapse/redshift factor in a count. r is the areal radial coordinate of the
nonrotating background. P=0 is the whole-star physical boundary in F §2.2. Current profiles
stop at `p_cut=max(1e-15 p_c,p_floor)`, not P=0; center starts at positive r and needs the
regular missing-center integral (ADR-0009; `TOVSolver.cpp:2537`, `:2550`).

**Domain is load-bearing.** F(12),(15),(30) say **core**; §2.1 explains free-particle diffusion
regions, and §3.5 explicitly neglects crust processes. The driver is not restricted to the
subset where a particular reaction is fastest (R95 §3.1). F's **total baryon constraint is
whole-star**, including matter outside that core. A core species count with a moving boundary
cannot be substituted for a whole-star species count without the boundary-flux term in §13.
The whole-star free-gas EOS has no nuclear crust splice (F §3.1), but that statement alone does
not supply an unpublished numerical core cutoff for F Fig. 1. This audit therefore labels its
scratch counts **whole free-gas star**, not authenticated published core counts.

For neutral npe-mu matter in the same domain, `N_B=N_n+N_p` and `N_p=N_e+N_mu`. For a future
more general species set the whole-star count is `N_B=sum b_i N_i`, with declared baryon
charges b_i; do not silently exclude baryons outside a selected reaction/diffusion domain.

## 7. Fixed-central-energy rotational derivative: independent derivation

Define `A_i^D=(partial N_i^D/partial q)_(e_c)` at q=0 with the boundary prescription fixed.
Use Eulerian Hartle variables `mhat=m0/q`, `phat=p0*/q`, `xihat=xi0/q=phat/nu'`, and
`s=omega_bar/Omega_geom`. In the particle-current product the lapse factors `(1+h0)` and
`(1-h0)` cancel. Radial metric contributes `mhat/(r-2m)`; angular averaging of
`(1/2)r² sin²(theta)e^-2nu s²` gives `r²e^-2nu s²/3`. Displacement of equilibrium isobars
gives `delta n_i/q=-xihat dn_i/dr`. A moving outer isobar adds `w n_i xihat` at its boundary.
The l=2 fields and boundary displacement multiply P2 whose spherical average is zero;
products of two l=2 perturbations enter only at q². These steps independently recover
H(109)–(114), F(20)–(24), and the accepted Phase-4 formula cited in §2.

For a single outer isobar R, the complete answer is

```
A_i = C { - integral_[0,R) w xihat dn_i
          + integral_0^R w n_i [ mhat/(r-2m) + r^2 exp(-2nu) s^2/3 ] dr
          + w(R) n_i(R-) xihat(R) }.                                          (PN2)
```

For a shell domain subtract the analogous inner-boundary term. If an imposed boundary has
displacement different from xihat, replace that boundary term with its actual displacement.
On smooth backgrounds the first term equals `integral w (dn_i/dp)(e+p) phat dr`.
This verifies, rather than changes, `PHASE4C_HARTLE2_DERIVATION.md:610` and ADR-0007 P13,
with ADR-0008 Q11's mandatory measure representation. No formula discrepancy was found.

## 8. Canonical dn_i measure, thresholds and independent oracle

The measure is the outward signed Lebesgue–Stieltjes measure of the **species number density**,
not dY alone: `-w xihat dn_i`. At a declared interface with inside/outside values n-/n+,
its contribution is `w_t xihat_t (n- - n+)`. At a finite-density terminal boundary the
n- -> 0 atom gives the last term of PN2 exactly once. Do not both append a vacuum density
row and add that shell again. No steepness detector is allowed (ADR-0008 Q3–Q7,Q11).

Recommended mandatory representation: explicitly reconstruct n_i at each profile node from
Y_i n_B, then use those n_i values and a declared radial partition as the source of truth.
On each segment integrate `-w xihat Delta n_i/Delta r`. Require every profile node to be an
integration boundary; optional immutable EOS-composition knots can refine the representation.
This extends the accepted energy-measure pattern; it is a **proposed particle representation**,
not an assertion that Q11 already selected its complete C++ implementation. Interpolating Y
and nB separately produces a different sub-cell polynomial from interpolating their nodal
product; the representation choice and its convergence error must be explicit.

Track-R onsets, from direct free-Fermi phase space and the ratified active-set ladder:

| Onset | Density behavior | Measure consequence |
|---|---|---|
| vacuum -> pe | n_p=n_e -> 0 continuously, nonrelativistic n proportional to enthalpy^(3/2) near zero | No physical vacuum atom; finite p_cut artificial terminal atom is bookkeeping only |
| pe -> npe | n_n proportional to (mu_n-m_n)^(3/2) on the allowed side; existing charged densities continuous | Continuous absolutely-continuous measure with nonanalytic higher derivatives; no neutron jump |
| npe -> npe-mu | n_mu proportional to (mu_e-m_mu)^(3/2); neutrality changes proton/electron response continuously | Continuous measure, no muon atom |

The equilibrium energy slope has a square-root cusp near neutron appearance; no full inactive
chemical Hessian is needed to integrate n_i values. The local N-1 refusal remains binding
for its own response API (`TrackRFreeGasThermodynamics.cpp:136`), not for value-only structure.

Independent measure oracle: pull the current back by r_phys=r+q xihat, or integrate by parts
`-integral w xihat dn_i + [w xihat n_i] = integral n_i d(w xihat)`. This differentiates the
geometry/displacement instead of the species density. Scratch 5-point segment quadrature
agreed with this oracle to `6.9e-15` relative at worst. A separate manufactured flat-current
fixture and a shrinking continuous transition versus a true atom are in §16. Future tests
must include atoms declared by fixture metadata; the current importer cannot represent
constant-pressure jumps (`TOVSolver.cpp:665`). Do not invent transition metadata now.

## 9. Surface and zero-pressure tail

F §2.2 uses P=0 for total A. PN1/PN2 for a whole-star implementation therefore require either
a P=0 estimator or a demonstrated error bound on a finite-cutoff result. ADR-0009 is unchanged.
At the scratch midpoint the cutoff is `4.91565287895e19 dyn/cm²` and the tail is wholly pe.
An independently coded direct free-gas enthalpy continuation gives the following diagnostic
values (counts and coefficients below include the factor 1e54):

| Quantity | n | p | e | mu |
|---|---:|---:|---:|---:|
| Omitted count, cutoff -> vacuum | 0 | 1.09619e43 | 1.09619e43 | 0 |
| Conservative count bound | 0 | 2.76e43 | 2.76e43 | 0 |
| Omitted **interior** A tail, count km² | 0 | 6.69685e49 | 6.69685e49 | 0 |
| Already included cutoff atom, count km² | 0 | 6.69287e49 | 6.69287e49 | 0 |
| Net tail replacement estimate, count km² | 0 | 3.98237e46 | 3.98237e46 | 0 |
| Absolute bound on tail minus cutoff atom, count km² | 0 | 1.35e50 | 1.35e50 | 0 |
| d(tail count)/dx from symmetric steps | 0 | 4.17984e42 | 4.17984e42 | 0 |

The number bound is `<5.7e-12` of the charged counts. Even the **triangle** bound on tail
replacement is `<7.4e-8` of A_e, below the **predeclared 1e-6** tail budget. It changes the
baryon-reduction central shift fraction by at most about `4e-10`; multiplying the converged
B-tail derivative by |dx/dq| adds only `1.22e46 count km²`. The conservative reduced charged
coefficient bound remains below `1.4e-8` relative. No zero-pressure correction is numerically
required at that proposed tolerance for this star, although reporting/bounding it is required.

Bound construction (not just radius-smallness): with H_s=1.30746556e-5 and monotone pe matter,
use `R_u=2m_s/[1-(1-2m_s/R_s)exp(2H_s)]`,
`m_u=m_s+4pi epsilon_s(R_u³-R_s³)/3`, and the maximum proper-volume factor W_u.
Then `N_tail <= C W_u n_s (R_u-R_s)`. Positive-source comparison of the tail Hartle equations
bounds s<=1, s', phat, xihat=phat/nu' and mhat: the computed enclosing values are
`s'_u=0.00710210 km^-1`, `phat_u=31.80252 km²`, `xihat_u=4817.989 km³`,
`mhat_u=473.1278 km³`. Thus density-measure tail <= C W_u xihat_u n_s, plus the number
bound times the bounded metric/velocity factor. Add the absolute cutoff atom for the
conservative replacement bound. These comparison inequalities and numerical envelopes are
implemented in scratch `tail_envelope.py:1`; the smooth central-step tail comparison is
`tail_bounds.py:1`. Species absent throughout this tail have exactly zero contributions.

The continuation from the stored-profile geometry has Delta R≈1.980439 m. The independent
all-GSL direct P=0 solve gives R0=12.7681549010 km; continuation from the stored-profile
mass gives 12.7681551998 km. This tiny difference is **not** silently called a new canonical
surface: it reflects the convention mismatch in §15. Rounded conservative bounds above
cover its much smaller effect. Future generic tails, self-bound surfaces and core boundaries
need their own EOS/domain treatment; this free-gas estimate is not a generic extrapolator.

## 10. Equilibrium-sequence derivative

`B_i=(partial N_i/partial e_c)_(q=0)` means the derivative of complete independently solved
nonrotating stars of the **same EOS and domain prescription**, including metric, threshold
and boundary motion. It is not the R2006 susceptibility. F(26),(27) explicitly provide the
analogous cumulative enclosed-baryon derivative at fixed P. For smooth backgrounds,

```
B_i = C { integral w [ (partial n_i/partial e_c)_r
                       + n_i (partial m/partial e_c)_r/(r-2m) ] dr
          + w n_i (partial R/partial e_c) } .                                  (PN3)
```

Interface terms are included in the derivative/measure; F(27) gives R_e=-P_e/P_r at a
fixed isobar. For the moving numerical p_cut(e_c), R_e=(p_cut,e-P_e)/P_r instead.
A P=0 count functional or explicitly differentiated tail eliminates that numerical-boundary
ambiguity. Scratch differentiates x=ln(rho_c/rho_ref): `B_x=e_c B_i`, exactly for constant
unit conversion. Both e_c and x must never be stored under an ambiguous `epsilon` name.

## 11. Production-method alternatives and numerical comparison

| Method | Exactness and moving features | Convergence/cost/provenance | Recommendation |
|---|---|---|---|
| A: regular homogeneous linearized TOV/Hartle | Continuum-equivalent to sequence derivative only when background and equations use consistent conventions; requires energy/species measures and interface/surface motion | One linear solve, no subtractive finite differences; center slope needed; same profile provenance; realistic sharp EOS requires partition/jump contract | Independent oracle; public ownership explicitly deferred by ADR-0007 |
| B: complete-star symmetric differentiation/local fit | Direct derivative of the canonical complete number functional; each star supplies its own metric/grid; threshold/surface motion follows the solved family | At least two step scales and a higher-order estimate; count errors amplified by 1/step; multiple solves and composite provenance; applicable to tabulated EOS with branch checks | **Accepted production owner** under ADR-0011, after the Q4 prerequisite |
| C: automatic differentiation/sensitivity of the actual canonical solver | Can differentiate the numerical primitive and surface event, but differentiability of branching/root/interpolation must be governed | Larger invasive ownership change; would modify protected canonical surfaces | Defer; not the smallest increment |

A vs B was actually exercised. On an analytic incompressible GR star, epsilon_c is constant
and cannot parameterize the sequence, so use central enthalpy H_c (and say so). For epsilon
=0.0004 km^-2, n=0.2 fm^-3, H_c=0.1, analytic `N/1e54=761.97584877` and
`(dN/dH_c)/1e54=8975.90134924`. The homogeneous calculation agrees to `1.34e-13` relative;
symmetric steps 0.004/0.002/0.001/0.0005 give
8974.965844 / 8975.667569 / 8975.842910 / 8975.886740, with quadratic step convergence.
It would be wrong to fabricate a B_epsilon for this incompressible fixture (`toy.py:1`).

On the Track-R star the **uncalibrated** homogeneous calculation fails the predeclared
`2e-4` particle-sequence target for n, p and e. At radial_res 10000/20000/40000 its B_x,n/1e54
is 132.635515 / 132.635713 / 132.635473, versus complete-star 132.757267.
This is not cured by radial refinement. The source/code diagnosis and calibrated control
are in §15. Choosing method B avoids importing this homogeneous discrepancy into B, but
does not retrospectively certify all current structural fields to a stricter physical
unit-consistency standard. ADR-0011 therefore blocks production implementation until a separate
unit-boundary reconciliation is accepted and completed.

## 12. Fixed-total-baryon reduction

At q=0 expand `N_i(e_c+delta e_c,q)=N_i0+B_i delta e_c+A_i q+O(q²)`.
Set `delta N_B=0`, with **whole-star** `A_B=sum b_i A_i`, `B_B=sum b_i B_i`:

```
(de_c/dq)_(N_B) = -A_B/B_B,             B_B must be nonzero and well-conditioned. (PN4)
K_i = (dN_i/dq)_(N_B) = A_i - B_i A_B/B_B.                                    (PN5)
sum b_i K_i = 0;            sum charge_i K_i = 0.                             (PN6)
```

For npe-mu, `K_n+K_p=0` and `K_p-K_e-K_mu=0`. No sign is chosen by convention: PN4 follows
from the expansion. F(18) is the corresponding cumulative chain rule, with denominator
`partial N/partial P` for compression. At an actual sequence turning point B_B can vanish;
F(28) discusses the divergence. An absolute `1e-30` fallback to A is forbidden. Require
uncertainty in B_B small compared with |B_B| and a branch-labelled, finite shift; do not
select an unqualified maximum-mass configuration.

## 13. Exact source I_Omega mapping and the core boundary term

Let `Q_B(P)=A_B(P)+B_B(P) de_c/dq` be the enclosed-baryon response on a fixed-P isobar.
F(18),(25) give

```
Pi(P) = (partial P/partial q)_(N,N_B) = -Q_B(P)/(partial N/partial P).
I_geom,i^core = integral_core dN (dY_i^eq/dP) Pi
             = -C integral_core Qbar_B(r) dY_i^eq(r),                         (PN7)
```

where Qbar_B=Q_B/C if the displayed outer C is used; never multiply a counted Q_B by C
again. This is the Stieltjes extension of F(30),(31), with signed interface increments.
It contains **only structural response**. R(7) preserves `dot N_i^eq=2 I_Omega,i Omega dot Omega`.
There is no additional e^nu multiplying the count driver; Omega and time are measured at
infinity. F §2.2 uses G=c=1 during the derivation, while the plotted I units are s².

For a whole star with P=0 and fixed total baryons, integration by parts gives
`I_geom,i^whole=K_i^whole`, because `Q_B(P=0)=0`. For a core ending on a fixed-pressure
boundary with fraction Y_i,c, the actual relation is

```
I_geom,i^core = K_i^(core, fixed-isobar boundary) - Y_i,c Q_B(P_core).           (PN8)
```

The second term subtracts particles carried by motion of that boundary. Thus merely reducing
core-integrated A_i and B_i does **not** guarantee the source driver, and its baryon sum need
not vanish until the boundary term is handled. The all-star baryon constraint is still used
in Q_B. A reaction-volume cutoff is a different object altogether.

For the ideal free-gas pe exterior, Y_p=Y_e=1 and Y_n=Y_mu=0, so dY_i=0 there. A domain
including all composition-changing layers gives the same PN7 as whole-star integration.
This establishes a useful conditional equivalence; it does not prove that FR2005 used an
unpublished core cutoff with that property. The production contract must either establish
that domain or expose cumulative Q_B and domain-qualified PN7/PN8. The scratch whole-star
K values in §18 are **not** asserted to reproduce the paper's free-gas core I_Omega.

## 14. Omega normalization and dimensions

C++ **proposed** number-response storage uses counted geometric coefficients. Existing Hartle
field units are authenticated at `CompactStar/Core/RotationSolver.hpp:283` and
`docs/validation/PHASE4_CLOSEOUT.md:85`. Conversion owner:
`CompactStar/AngularVelocity.hpp:150`, using `c=299792.458 km/s`.

| Quantity | Stored/proposed unit | Physical/source unit and conversion |
|---|---|---|
| n_i | fm^-3, reconstructed | multiply by 1e54 for km^-3 in km-volume integrals |
| N_i, N_B | count | no conversion/redshift |
| mhat, xihat | km³ per geometric Omega² normalization | m0, xi0 = hat times q, km |
| phat, delta phat | km², dimensionless respectively | p0*=phat q; delta p0=delta phat q in km^-2 |
| A_i, K_i, I_geom,i | count km² | I_phys=I_geom/c², count s² |
| B_i=dN/de_c | count km², e_c in km^-2 | B_rho=B_i*(de_c/drho_c), count cm³/g; actual conversion constant must be named |
| B_x=dN/dln rho_c | count | B_i=B_x/e_c |
| de_c/dq | dimensionless | d rho_c/dOmega_phys² = (de_c/dq)/[(de_c/drho_c)c²] |
| dx/dq | km² | dx/dOmega_phys²=(dx/dq)/c² |
| dot N_i^eq | count/s | `2*(I_geom/c²)*Omega_phys*dotOmega_phys`; equivalently `2 I_geom Omega_geom dotOmega_geom` |
| chemical B/Btilde | count/erg or count/MeV | later global thermodynamics, not sequence B_i |
| chemical Z; W | erg/count; erg s² | later energy-valued coefficients, not structural K |

The geometric and physical versions are not independent state fields. Use an explicitly
named conversion/materialization method through the existing AngularVelocity owner. A hidden
factor c² would be a factor 8.987551787e10 error. Radians are dimensionless. Every finite-spin
scratch q is explicit: 1e-6 km^-2 corresponds to 299.792458 rad/s, with subsequent q halvings
reducing Omega by sqrt(2).

## 15. dEdP adjudication and the separate convention mismatch

**Every current RotationSolver use was traced:**

| Path | Use of profile dEdP | Relevance of outer neutron-threshold error |
|---|---|---|
| First-order Hartle | No read | Exactly none |
| Monopole workspace | Requires complete column; gathers/interpolates it (`RotationSolver.cpp:1137`, `:1177`, `:1219`) | Presence/finite-value gate, not an interior source |
| Fixed-e_c center series | b0.dedp in r0^5 mass start (`:1481`) | Only if the center itself is near the feature; this midpoint center is far above it |
| Monopole interior | eps_slope measure, explicitly not dedp (`:1280`–`:1313`) | Exactly none at fixed profile values |
| Derived xi, pressure, surface | Uses governed fields, nu', values (`:1552`–`:1635`) | No additional derivative read |
| Homogeneous sequence candidate | Not a public current solver; scratch center normalization uses central dedp, interior energy measure | Outer error irrelevant; center accuracy and unit consistency matter |
| Particle A / complete-star B | Scratch measure / whole-star count differentiation | No derivative-column source; legacy nodal dn/dp is a different forbidden approximation |

Scratch perturbed **only TOVPoint.dedp** by +0.410% on all nodes with
`1e-10<n_B<1e-6 fm^-3`, covering neutron appearance, then reconstructed the star normally.
The derivative changed at 66 profile nodes. First-order and every exported monopole field
were bit-identical; A and derived K changed by **zero**.
A deliberately broader +0.410% all-node control changed the minute center mhat start but not
the final integrated A at stored precision. This does not license deleting the central
pointwise derivative or using a center at onset. Source tracing is the exact irrelevance
proof; the perturbation is a mutation-sensitive input-dependency check.

The **table** ladder also changes EOS values and the central pressure inversion, so its
changes must not be attributed to dEdP alone. Scratch maxima relative to the 8192 table:

| Table resolution | I km³ | deltaMhat km³ | max relative A_i | max relative B_x,i | max relative K_i |
|---|---:|---:|---:|---:|---:|
| 1024 | 31.4395102752 | 473.593234178 | 2.0263e-4 | 2.5443e-3 | 7.5313e-5 |
| 2048 | 31.4397154030 | 473.602565224 | 2.8878e-7 | 6.5203e-6 | 1.1304e-7 |
| 4096 | 31.4397160601 | 473.602589915 | 4.5585e-8 | 2.9228e-6 | 2.2386e-8 |
| 8192 | 31.4397161287 | 473.602592803 | reference for this comparison only | reference | reference |

At common fractional-radius sample points, 1024->8192 maximum norm-relative field differences
are s `2.55e-6`, mhat `1.99e-5`, phat `1.09e-5`, xihat `2.77e-5`; 4096->8192 are
`7.3e-10 / 6.9e-9 / 3.2e-9 / 8.1e-9`. Center dedp values are
13.2603500 / 13.2921041 / 13.2919769 / 13.2920153; center phat is about 2.91504584e-11 km².
The coarse sequence derivatives do not meet the independent-method target; this is retained,
not hidden by the stable reduced result. Nested threshold-distance grids are already in the
Structure-1 generator (`tests/eos/structure1/table.hpp:48`). Further localized derivative
correction is not required for the selected measure-plus-complete-star path.

**Classification of the carried 0.410% item: A, irrelevant as a derivative-column error to
this governed path at this star.** Table/value convergence is a separate mandatory check.
No analytic slope substitution, table-policy change, or production correction is made.

**Separate new finding: homogeneous/background convention mismatch.** The raw homogeneous
B_x errors versus complete-star differentiation are
`-9.1742e-4 / -2.1428e-4 / -2.3443e-4 / -3.0208e-5` for n/p/e/mu. GSL TOV uses
G=6.673e-8 CGS and M_sun=1.98892e33 g (`TOVSolver.cpp:1511`, `:2595`). NStar stores mass
through Zaki SUN_M_KM and pressure/energy through Zaki conversions (`NStar.cpp:180`–`:196`).
The loaded constants are:

```
Zaki G = 6.6743e-11 SI; SUN_M_KM = 1.4766250380501249 km
alpha = stored epsilon / epsilon_in_GSL_geometric_units = 1.0001948149258202
beta  = stored m / m_in_GSL_geometric_units             = 0.9999382793775072
```

Meanwhile nu' is carried from the GSL solve (`NStar.cpp:186`). Consequently the stored
background does not exactly satisfy a single conventional geometric TOV system. This
extends the previously recorded solar-mass/apparent-radius debt into a demonstrated
sequence-response error; it does not change the accepted static printed-bin result.

**Diagnostic control, not production repair:** convert only scratch homogeneous-equation
inputs to GSL geometric epsilon,p,m, retain the carried nu', then convert the resulting
mass derivative back to the stored-profile convention for the number metric. B_x becomes
132.75713829 / 4.705393366 / 4.241071799 / 0.4643215671 in units 1e54; disagreement falls to
`9.74e-7 / 6.36e-7 / 7.05e-7 / 5.82e-9`. This isolates the cause without changing either
canonical solver. Using the raw homogeneous B instead of complete-star B changes K by up
to `8.72e-4` relative. ADR-0011 requires a separate governed unit-boundary reconciliation
before Phase-5B production implementation. The scratch control is diagnostic only and does not
constitute the accepted PB7 oracle. Passing the looser proposed K budget does
not cure a failed independent B gate. No existing ADR or ratification is rewritten.

## 16. Scratch A results, partition and finite-spin checks

All prototypes/results are outside the tree at `/private/tmp/phase5b0-inv09` and are
**not baselines**. Fresh Debug production library built from the authenticated worktree;
CMake 4.2.1, AppleClang, C++17, existing GSL/Zaki dependencies. No source/test/build file in
the worktree was edited. The standalone link initially needed explicit zlib (`-lz`);
a scratch Python syntax error was corrected before results were generated. No failed
scratch output is represented as successful evidence.

Predeclared scratch targets were written in `protocol.txt` before number-result inspection:
1e-8 same-representation count quadrature; 1e-12 conservation residual scaled by the sum
of absolute summands;
2e-4 independent sequence agreement (existing ADR-0008 sequence standard, proposed extension);
1e-3 per nonzero species response stability with absolute floor 1e-8 sum|A|;
1e-6 relative tail contribution. These are proposed implementation budgets, not owner acceptance.
The raw B charge identity fails that original species-sum budget (§17). Explicit propagated
stencil error in the future PB10 contract is a proposed refinement, not a retroactive pass.

At table8192/radial_res40000 (10132 actual output nodes), all values in the following table
are divided by 1e54; A columns consequently have numeric units count km²/1e54:

| Species | N | density measure | radial metric | rotation time-dilation | terminal atom | complete A |
|---|---:|---:|---:|---:|---:|---:|
| n | 756.82753225 | 374027.869908 | 10201.681073 | 12204.832651 | 0 | 396434.383632 |
| p | 4.859479123 | 1755.804060 | 36.010307 | 48.456840 | 6.69287e-5 | 1840.271274 |
| e | 4.811901927 | 1744.641333 | 36.001257 | 48.400748 | 6.69287e-5 | 1829.043404 |
| mu | 0.04757719621 | 11.162726516 | 0.009050558 | 0.056092424 | 0 | 11.227869498 |
| B=n+p | 761.687011373 | 375783.673967 | 10237.691380 | 12253.289492 | 6.69287e-5 | 398274.654906 |

The independent direct-vacuum count oracle gives
756.83114853 / 4.859497616 / 4.811920397 / 0.04757721970, sharing neither production EOS
interpolant nor radial TOV ODE. The small difference is dominated by profile metric convention
and sampling, not the tail. The direct solver uses 96-point phase-space quadrature and
DOP853 in enthalpy; tightening rtol 2e-11->2e-12 changes neutron N by 2.6e-8 in these units.

A radial partition sweep 10000/20000/40000 changes the worst species A relative to the
finest by `5.14e-6 / 1.02e-6 / 0`. Both neutron and muon appearance are crossed; the table
ladder separately refines their distance grids. On the finest profile, even a corrected
**naive nodal dn/dp density term** differs by `9.83e-6 / 1.37e-6 / 1.38e-6 / 3.81e-7`
relative in complete A. Its small error on this smooth model does not override Q11 or
validate true jumps. The legacy expression additionally omits the large metric/velocity terms.

Independent nonlinear pullback materialization uses r_phys=r+q xihat, its Jacobian,
advected n_i, displaced background mass plus q mhat, and the angularly integrated fluid
Lorentz factor. It does **not** define N(q)=N(0)+qA. Relative errors of
`[N(q)-N(0)]/q` versus A_n at q=1e-6,5e-7,2.5e-7,1.25e-7,6.25e-8 are
`2.2436e-4 / 1.1217e-4 / 5.6085e-5 / 2.8042e-5 / 1.4021e-5`; the other species converge
with the same q order. This validates current assembly and Omega² normalization; it is
not an independent full rotating-star solve or evidence for physical O(Omega^4) fields.

Manufactured current fixture: flat geometry, n=1-r/R, xi=k r, R=2,k=0.3 gives
N=8.377580410, A=7.539822369 and N(q)=(1+kq)^3 N. Difference quotients at
q=.01/.005/.0025/.00125 are 7.562464455/7.551137757/7.545478649/7.542650155.
A unit density drop at r=.73 has exact measure contribution 1.466559539; continuous ramps
of width .1/.01/.001/.0001 give 1.473439627/1.466628340/1.466560227/1.466559546.
These independently check sign, Jacobian, discontinuity limit and quadratic spin dependence.
They are kinematic analytic fixtures, not a newly validated physical Hartle model.

## 17. Scratch B results

Complete canonical stars were solved at x=±0.004,±0.002,±0.001,±0.0005 for **each** table
resolution. Achieved central states were checked against requested states; no clamp or
inverse-mass fit occurred. Each star had its own grid, metric and n_i values. Count changes
from the P=0 tail were independently estimated (§9). For the finest table:

| Symmetric x step | B_x,n/1e54 | B_x,p/1e54 | B_x,e/1e54 | B_x,mu/1e54 |
|---|---:|---:|---:|---:|
| .004 | 132.757006706 | 4.705408908 | 4.241068672 | .464340236 |
| .002 | 132.757156808 | 4.705396837 | 4.241070784 | .464326054 |
| .001 | 132.757235557 | 4.705396074 | 4.241073269 | .464322806 |
| .0005 | 132.757259498 | 4.705396285 | 4.241074406 | .464321879 |
| 5-point, h=.0005 | 132.757267479 | 4.705396355 | 4.241074785 | .464321570 |

Five-point formula: `[8(N(h)-N(-h))-(N(2h)-N(-2h))]/(12h)`.
`B_x,B/1e54=137.462663834`; divide each B_x by
`e_c=8.16877629603e-4 km^-2` for the stored geometric B_i.
The direct free-gas vacuum sequence gives B_x values
132.760282809 / 4.705419114 / 4.241097355 / .464321759, within `2.28e-5` relative of the
canonical finite differences, including their differing metric convention. Its symmetric-step
ladder is retained in `direct.json`; no single step is called validation. Raw homogeneous
comparison **fails** as recorded in §11/§15; calibrated comparison succeeds as a diagnostic.

Raw floating-point charge residuals (in the displayed 1e54 scaling) are N `-1.07e-14`,
A `6.83e-13`, and B_x `1.48e-11`. B_x magnifies count-quadrature rounding by 1/h; it exceeds
an unamplified species-only 1e-12 relative identity budget. Production must use compensated
accumulation and explicit propagated stencil error, and/or reconstruct dependent counts
from governed baryon/charge identities while retaining independent raw checks. This audit
reports the residual rather than asserting an exact arithmetic zero.

## 18. Scratch fixed-baryon reduction and finite-spin closure

With the independent complete-star B estimate:

```
dx/dq = -2897.329673369 km^2
de_c/dq = -2.36676379576
K/1e54 = [11792.8132103, -11792.8132103, -10458.7484176, -1334.06479270] count km^2/1e54
I_phys (whole-star conditional mapping) =
[+1.31212742795e47, -1.31212742795e47, -1.16369270130e47, -1.48434726638e46] count s^2
```

The raw baryon residual K_n+K_p is `-3.64e-12` in 1e54 count km² units; raw charge residual
K_p-K_e-K_mu is `-4.28e-8`, explained by the B stencil's amplified count rounding (§17).
Both are negligible relative to the **assembly contribution** scale, but neither is hidden
by manually setting a reported result to zero. Exact dependent reconstruction is an additional
future production gate, not evidence that this audit closed INV-09.

Fresh complete stars at x(q)=(dx/dq)q were actually solved; their own governed Hartle responses
were materialized through the nonlinear current pullback. These are not central-star
coefficients algebraically rescaled to force a result:

| q [km^-2] | Omega_phys [rad/s] | fractional Delta N_B | error in Delta N_n/q versus K_n | error for mu |
|---|---:|---:|---:|---:|
| 1e-6 | 299.792458 | 8.99133e-7 | 5.61778e-2 | -8.82227e-3 |
| 5e-7 | 211.985280 | 2.24802e-7 | 2.80897e-2 | -4.41708e-3 |
| 2.5e-7 | 149.896229 | 5.61483e-8 | 1.40314e-2 | -2.20957e-3 |
| 1.25e-7 | 105.992640 | 1.41284e-8 | 7.06295e-3 | -1.10495e-3 |

Delta N_B contracts about fourfold per q halving: **O(q²)=O(Omega^4)**. Species quotient
errors contract about twofold. The neutron reduced coefficient is a cancellation of much
larger terms, making its finite-q relative error larger than the raw A error. This is the
expected asymptotic trend, not a claim that the smallest finite spin already meets the
1e-3 extrapolated coefficient tolerance. F1/F2 algebra and F3/F4 order are supported in
scratch; acceptance still requires the full ladder, independent review and human ratification.

## 19. Published numerical target inventory

| Category | Source item | What it can validate / limitation |
|---|---|---|
| EXACT TABULAR TARGET | F Table 1 common-state free-gas M,R,R_infinity at printed rho_c | Static input/structure, already qualified-ratified; no particle-response numbers |
| FIGURE-DERIVED TARGET | F Fig. 1 free-gas compression at P_c=4.5e34 | Structural compression profile, only after core-radius/domain and adopted Omega_K conventions are fixed |
| FIGURE-DERIVED TARGET | F Fig. 3 middle-right I_Omega,e and I_Omega,mu versus mass, ordinate 1e47 s² | Genuine pre-W structural target for APR, not an available free-gas target; realistic EOS remains outside this task |
| SOURCE FORMULA ONLY | F(18)–(31), H(109)–(114), R(7) | Count measure, fixed-baryon mapping, sign, normalization and interfaces |
| NO NUMERICAL TARGET | Free-gas tabulated species A/B/K or I_Omega | None found in F/R; scratch cannot become a published target by analogy |
| Outside structural scope | R Fig. 1 corrected Z; R Fig. 2 transients; F W curves | Later chemical-layer validation, forbidden to build here |

If later authorized, digitization must retain PDF hash/page, exact curve/EOS, axis transforms,
line-width/pixel uncertainty, independent second extraction, source configuration and domain,
and interpolation uncertainty distinct from solver tolerance. Do not infer intermediate
free-gas coefficients from chemical W or fit a star to a curve. No digitization was performed.

## 20. PB1–PB14 implementation validation contract

Let E_count and E_A denote independent quadrature/representation error estimates, and let
S_i=|A_i|+|B_i A_B/B_B|. PB1–PB14 are the accepted validation plan, not achieved validation.
Identity tolerances must include the propagated stencil error, not just machine epsilon.
The existing numerical budgets remain proposed implementation acceptance targets fixed before
the implementation's results; they are not permission to widen a failed bound. Near-zero K uses an absolute scale, and an
ill-conditioned B_B refuses. Convergence envelopes may be nonmonotone in node placement,
as already adjudicated for ADR-0008; arbitrary strict monotonicity is not substituted.

| ID | Independent oracle | Metric / predeclared tolerance logic | Exact defect falsified |
|---|---|---|---|
| PB1 | Analytic proper-volume counts; direct phase-space/enthalpy TOV | Same representation <=1e-8; cross-method <= bounded interpolation/convention budget; raw and reconstructed baryon/charge identities | Y used as n, wrong metric/count units/domain |
| PB2 | Current pullback and angular integration on analytic fixtures | Exact algebra; numerical IBP equivalence <=1e-10 | Missing lapse cancellation, 1/3 velocity factor, l=2 contamination |
| PB3 | Piecewise analytic densities with declared jumps plus shrinking continuous ramps | Exact signed atom, segment total measure; <=1e-10 same representation; convergent weight envelope | Nodal-only dn/dp, invented onset atom, steepness detector, duplicate terminal atom |
| PB4 | Analytic/toy physical Hartle model, independent current quadrature | <=1e-8 analytic regime plus demonstrated radial convergence | Missing metric/density/velocity/surface A term |
| PB5 | Explicit-spin nonlinear current materialization | Quotient error proportional to q before floor; exact physical/geometric conversion | Seed² contamination, missing Omega² division or c² |
| PB6 | Radial and independent EOS-knot/threshold partitions | <=1e-3 per nonzero species, absolute floor 1e-8 sum|A|; vanishing cell-weight error envelope | Lost sharp continuous composition, node-placement pathology |
| PB7 — **BLOCKED on Q4** | Complete-star sequence versus independently derived homogeneous/sensitivity solution | <=2e-4 after **one declared convention**; unresolved convention mismatch is a failure | Wrong held-fixed variables, solver/background inconsistency |
| PB8 | Symmetric step ladder and 5-point/local-polynomial derivative | Adjacent estimates agree within propagated count/root/step budget; <=2e-4 nonzero B target | Single-step cancellation, wrong achieved center, branch/clamp crossings |
| PB9 | Independent baryon sum and explicit expansion PN4 | Zero within compensated arithmetic + propagated ingredient uncertainty; conditioning gate | Wrong sign, wrong baryon species/domain, B_B fallback |
| PB10 | Independent raw charged-species integrals plus neutral reconstruction | Same propagated arithmetic/derivative budget, <=1e-12 of assembly contributions where achievable | Species ordering/charge map or inconsistent quadratures |
| PB11 | Newly solved shifted-center stars plus independent finite-current materialization | Delta N_B proportional to q²; Delta N_i/q -> K_i linearly in q before numerical floor | Fixed-e_c mislabeled fixed N_B, omitted shift/boundary or trivial linearized self-test |
| PB12 | dEdP-only perturbation, EOS-value ladder, center-location sweep | Exact zero where dependency absent; end-to-end K <=1e-3; PB7 still independent | Blaming column error for value error, threshold-center misuse, hidden pointwise consumer |
| PB13 | Direct-vacuum/tail comparison inequalities and central-step tail derivatives | <=1e-6 of count/response scales, or retain a bounded explicit correction; include shifted cutoff | R_cut silently called P=0, ignored species boundary motion |
| PB14 | F(18),(24),(26),(30),(31), R(7), conditional figure extraction | Exact symbol/sign/units/domain; published comparison tolerance set by source precision, not fitted output | Whole-star/core driver confusion, boundary flux omission, chemical Z/W substituted for I |

Additional required falsifiers under these rows: mutate n_i->Y_i, omit each A term, omit or
double the terminal atom, change c², use the central star's metric on a nearby grid, flip
PN4's sign, accept stale profile/response provenance, and introduce B_B near zero. All must
fail without being retained. The scratch campaign did not execute this entire production
mutation ladder or an independent physical rotating-star oracle; PB1–PB14 are the accepted
validation plan, **not fourteen claimed passes**. PB7 cannot execute as an accepted oracle until
the unit-boundary reconciliation is accepted and completed.

## 21. Accepted ownership/API contract

ADR-0011 accepts generic, domain-qualified `ParticleNumbers`,
`FixedCentralEnergyNumberResponse`, `EquilibriumSequenceNumberDerivative` and
`FixedBaryonNumberResponse` objects. Use one generic Core/Analysis particle-number functional and a generic structural
analysis helper consuming **the owning star's current** HartleFirstOrderResponse and
HartleMonopoleResponse. A sequence-analysis owner orchestrates complete canonical TOV stars
and differentiates that same number functional; it never supplies one star's geometry to
another. A thin later rotochemical structural adapter maps the domain-qualified cumulative
baryon response to F(30). Do not add these results to BNV or resurrect RotochemicalCache.

A result carries domain/surface semantics, species/baryon/charge metadata, e_c definition,
q normalization, units, central-step/partition error estimates, achieved central state,
validity/refusal reason and profile provenance. Runtime provenance uses profile identity +
Version and the existing lifetime rule (ADR-0003), including matching both Hartle inputs.
A sequence result also depends on every contributing star, EOS revision, stencil and domain
policy. Prefer eager value results with no second lazy cache; do not invent persistent IDs.
Optional immutable measure metadata must be profile-attached and version-changing. A generic
layer must not know free-gas masses, neutron/muon thresholds, N-1 workarounds or Track-R
BarotropeAt; model-specific validation/tails are separately named adapters/oracles.

ADR-0011 accepts complete canonical-star differentiation as the production owner, including
multi-scale symmetric steps, a higher-order estimate, achieved central states and an explicit
materialized derivative with respect to `epsilon_c`. The regular homogeneous method remains an
independent PB7 oracle only; no public homogeneous Phase-4 API is created. Exact source-file/API
placement remains an implementation design detail subordinate to this contract. No public
header/source is created by this audit, and implementation is blocked on Q4.

## 22. Exact INV-09 closure criteria and chemical boundary

INV-09 may be changed to VERIFIED/RESOLVED only after all of the following exist together:
correct PN1 species/domain semantics; implemented measure-complete PN2 including declared
interfaces and terminal policy; independently validated PN3 sequence derivative; conditioned
PN4/PN5 with whole-star baryon constraint; explicit Omega conversion; PB1–PB14 including
nontrivial finite-spin closure; raw and dependent baryon/charge checks; threshold, center,
surface and central-step convergence; dependency-complete provenance and stale-input
refusal; resolved source/core-domain and unit-convention questions; no legacy candidate
activation; independent scientific review; and separate human-owner ratification.
A later regression artifact follows validated implementation under applicable governance;
these scratch values must not be frozen as a reference in advance.

Phase 5B stops at domain-qualified I_Omega / particle-number drivers and their ingredients.
R §3 and §5 leave the compression/number driving quantities unchanged; R §4 corrects the
chemical-potential-to-number map, replacing F(9),(12),(13),(54)–(56). In intrinsic notation
its corrected susceptibility is `integral e^-Phi [chi-chi q(q^T chi q)^-1 q^T chi] dV`,
symmetric with charge null mode, and the full four-species matrix is singular (R13).
The later local neutral-Hessian adapter and reduced stellar coefficient layer are governed
by ADR-0010. They require this phase's **structural driver**, matching EOS/star/domain and
physical Omega units; they alone later combine it with corrected chemical susceptibilities.
No Btilde inverse, chemical Z, W, evolved eta ordering, rate, or evolution is formed here.
INV-11 remains unresolved regardless of any future INV-09 structural closure.

## 23. Owner decisions and remaining blocker

ADR-0011 records the owner's answers to Q1–Q4: the Structure-1 midpoint fixture; whole-star
fixed-baryon counts with domain-qualified source mapping; complete canonical-star differentiation
as the production `B_i` owner; and a mandatory separate GSL/Zaki unit-boundary reconciliation.
The first three define the implementation contract. The fourth blocks implementation until its
own governed audit, minimum correction and revalidation are accepted and completed.

No core boundary is supplied by this adjudication, so no numerical FR2005 free-gas core
`I_Omega` reproduction claim is available. INV-09 remains **INTENDED BUT UNVERIFIED** and
INV-11 remains **UNRESOLVED**.

## 24. Blocked implementation scope, evidence custody and verification

**Implementation scope after the separate Q4 prerequisite is completed:** generic nonrotating/cumulative
number functional; fixed-e_c number measure from existing governed responses; one owner for
complete-star sequence differentiation; conditioned whole-baryon reduction; domain-qualified
F(30) mapping; explicit units/provenance/errors; PB1–PB14 and an independent review. Production
TOV/NStar/RotationSolver changes needed for the unit-boundary decision belong only to the
separately authorized and validated prerequisite task. Do not smuggle them into particle-number
integration. **No Phase-5B production implementation is authorized now.**

Evidence directory `/private/tmp/phase5b0-inv09` contains `export.cpp` (canonical-star/response
export and derivative-column controls), `analyze.py` (measure, homogeneous, nonlinear current),
`direct.py` (independent phase-space EOS, enthalpy TOV/counts, tail), `calibrate.py`, `toy.py`,
`extra.py`, `tail_bounds.py`, `tail_envelope.py`, `protocol.txt`, tables/profiles, build/run logs,
`results.json`, `direct.json`, `closure.json`, `calibrated.json`, `toy.json`, and `SHA256SUMS.txt`.
These files are disposable audit evidence outside Git, never runtime inputs or baselines.
Reruns must use a fresh output directory because table generation refuses overwrites.

Reproduction: configure the exact starting checkout out-of-source in Debug with the existing
GSL/Zaki toolchain and `/Users/keeper/miniforge3/bin/python3`; build only CompactStar; compile
the standalone export with `compile.py` and explicit `-lz`; generate the 4x3 table/radial
ladder plus sequence stars; execute the listed Python scripts; solve the explicit shifted
x(q) stars before evaluating `finite` on them. The executable consumers use public NStar
responses, never baseline TSVs or raw private seeds. No production CTest suite was rerun:
this is documentation with scratch numerics, and no production/test behavior changed.

The pre-adjudication preflight verification passed: `git diff --check`; an exact four-document
tracked/untracked allowlist; unchanged bytes for every tracked accepted ADR and the complete
INV-11 section; zero production/test/baseline/data diff; and a repeated full literature checksum
check. Both new files were checked separately for trailing whitespace because ordinary
`git diff --check` does not inspect untracked files. Canonical and remote integration state is
authenticated separately when the accepted documentation commit is created and fast-forwarded.

The scratch checksum manifest verified all 21 listed source/result/metadata files when run
from its own directory. Its SHA256 is
`ba704893ae8167d2f56bd965a109902ab6b2a0e976ebb4a4787abb06f99f0454`.
The only permanent files are this record, roadmap/invariant pointers, and accepted ADR-0011.
The owner adjudication supplies a clear structural contract and authorizes documentation
integration. It does not authorize Phase-5B production implementation.

**Exactly one recommended next action:** open a fresh unit-boundary reconciliation audit from
canonical master to identify and adjudicate one consistent GSL/Zaki geometric convention across
TOVSolver -> NStar/StarProfile -> RotationSolver, quantify the minimal repair, and define the
required static/Phase-4/Structure-1 revalidation before any Phase-5B production implementation.
Do not begin it automatically.

## 25. Explicit non-goals and disposition

No production code, tests, baseline, EOS dataset, literature file, previously accepted ADR,
canonical TOV/NStar/RotationSolver behavior or evolved state was modified. ADR-0011 is accepted
as a structural contract, while implementation remains blocked. No Btilde, paper Z/W, weak-rate
change, evolution, superfluidity, APR/BPAL implementation, DS(CMF) off-equilibrium physics,
BNV or legacy activation occurred. The preflight establishes the equations and accepted
validation/ownership contract; it does not close INV-09 or INV-11 and does not waive the
unit-boundary prerequisite.
