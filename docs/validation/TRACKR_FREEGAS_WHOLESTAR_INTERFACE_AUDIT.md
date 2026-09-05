# Track-R Structure-0 — FR2005 free-gas static whole-star interface audit

> **TRACK-R FREE-GAS WHOLE-STAR STRUCTURE INTERFACE AUDITED — IMPLEMENTATION READY**
>
> This is an implementation/validation plan, not a reproduced star or a scientific
> baseline. The ratified local model remains unchanged. Readiness means that the
> next bounded task can begin with the value-only barotrope prerequisite below;
> it does not mean that today's public provider is already a total TOV EOS.

## 1. Authentication, authority, and scope

- Date: 2026-09-05; one agent, serial scientific diagnostics.
- Starting SHA: `a1a1ea5233c156d69410a0d3d61750f3da57e458`.
- Canonical checkout: `/Users/keeper/Documents/CompactStar/repo/CompactStar`.
  Clean `master`; HEAD = local master = origin/master = live
  `git ls-remote origin refs/heads/master` at the starting SHA.
- Fresh branch: `analysis/trackr-freegas-wholestar-interface`.
- Fresh worktree:
  `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-trackr-freegas-wholestar`.
  Created at precisely the starting SHA; existing thermodynamic worktrees were not reused.
  `git log --all --not <starting-SHA> --` the two authorized documentation paths
  returned no divergent commits. Other worktrees were inventoried before creation.
- Change class: **documentation** (`GOVERNANCE.md:43-57`). No number produced by
  production code changes. No pre-baseline correctness exception is invoked.

Governing ownership is ADR-0005, especially its canonical-primitive contract and
completed retirement of the duplicate loop
(`docs/adr/ADR-0005-canonical-tov-numerical-primitive.md:206-253`, `:566-614`),
ADR-0009's accepted surface decisions
(`docs/adr/ADR-0009-tov-surface-event-and-termination.md:48-68`), and ADR-0010's
cold local provider/stellar-background boundary
(`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:227-247`,
`:334-353`). INV-01/02/03/06/13 govern fractions, unit boundaries, metric, surface,
and interpolation (`docs/SCIENTIFIC_INVARIANTS.md:56-146`, `:148-174`,
`:276-327`, `:930-984`). Phase-4 objects are ratified only within their recorded
scope (`docs/validation/PHASE4_CLOSEOUT.md:64-176`); they do not supply global
particle response. Local ratification is current authority
(`docs/validation/PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:27-46`, `:94-135`).
The implementation record and independent review were read, especially the
N-1/P-7 distinction (`docs/validation/PHASE5A5_TRACKR_PE_BRANCH.md:193-221`;
`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:403-469`, `:682-726`).

Current and target architecture and roadmap were consulted. Historical blanket
claims in target architecture are not present-code authority. In particular,
the target document's two-live-TOV statement is superseded by accepted ADR-0005
and current architecture section 3
(`docs/architecture/TARGET_ARCHITECTURE.md:3-6`, `:36-49`;
`docs/architecture/CURRENT_ARCHITECTURE.md:253-297`). No broad documentation
cleanup is part of this audit.

## 2. Direct source ledger

All PDF page numbers below are **one-based file pages**. FR2005's printed
preprint page numbers agree with these; they are not journal pagination.
The shared library was read only; its README and catalog identify FR2005 as
the non-superfluid reproduction authority and R2006 as the corrected chemical
response authority (`/Users/keeper/Documents/CompactStar/literature/README.md:14-17`;
`/Users/keeper/Documents/CompactStar/literature/catalog.tsv:2-4`).

| Key | Exact source and locators read | Use in this audit |
|---|---|---|
| F05 | `/Users/keeper/Documents/CompactStar/literature/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`; section 3.1, PDF/printed pp. 8-9; Table 1, p. 30, caption, headings and **all** footnotes a-d | Authoritative EOS setup and numerical row. Rendered pp. 9 and 30 checked against text extraction. |
| F05-metric | Same PDF, section 2.1, p. 3, eq. (6) and following definitions; section 2.3, p. 6, eqs. (20)-(21) and metric definition | Coordinate radius, effective radius, local baryon density versus relativistic energy density; `G=c=1` in the structure discussion. |
| F05-context | Same PDF, section 2.3, p. 7, discussion of Fig. 1; Fig. 1 caption, p. 26; section 4.3, p. 16, eq. (71); references pp. 24-25 | Fig. 1 uses fixed central pressure `4.5e34 dyn cm^-2`; it does not identify Table 1's central input. Eq. (71) supports `2 mu_n = mu_Sigma + mu_p` at hyperon equilibrium. |
| R95 | `/Users/keeper/Documents/CompactStar/literature/rotochemical/1995-Reisenegger-Rotochemical-Heating.pdf`; section 2, PDF pp. 14-15 / printed pp. 4-5, eq. (3) | Electron-inclusive energy minimization and `n_p=n_e` on npe matter; **not** authority for FR2005's relativistic density convention. Eq. (3) visually checked despite Type-3-font rendering warnings. |
| R06 | `/Users/keeper/Documents/CompactStar/literature/rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`; introduction and section 2, PDF pp. 1-2 / journal pp. 568-569, eqs. (2), (4) | Supports intrinsic/local-neutral chemical semantics; its electrostatic-response correction does not redefine Table 1's static target. No B-matrix construction follows. |

Authenticated SHA-256:

```text
F05 f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88
R95 9af85e37c7a52fd5b704c0ba07cc0ad89741d23b049df31cb6867d501d91d0ff
R06 a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8
```

FR2005 directly settles the important cited-source qualifications: section 3.1
names Shapiro & Teukolsky (1983) for the ideal gas; Oppenheimer & Volkoff (1939)
for relativistic structure; Lasota et al. (1996) for the other rows' empirical
Kepler estimates; and Haensel, Salgado & Bonazzola (1995) for the adopted
pure-neutron period. These are citations **by FR2005**, not claims that their
complete contents were independently read here. An attempted ADS retrieval of
the last paper was unavailable. It is unnecessary to infer its provenance:
FR2005 states it explicitly. No numerical information from an inaccessible
source is used. The online [FR2005 record](https://arxiv.org/abs/astro-ph/0502116)
confirms bibliographic identity; the authenticated local PDF governs this ledger.

## 3. Exact Table-1 row and target classification

The p. 30 row label is `Fermi gas`. Exact heading/value transcription:

| Printed heading | Exact printed value | Unit | Significant figures | Meaning and acceptance role |
|---|---:|---|---:|---|
| `M_max` | `0.62` with superscript c | `M_sun` | 2 | Gravitational mass of the listed nonrotating free-gas configuration; **target**, subject to footnote c below. |
| `rho_c` | `1.10` | `10^15 g cm^-3` | 3 | Central energy density expressed as mass density; **central-state constraint/target**. It is not an independent predicted output if used to select the central input. |
| `R` | `12.77` | `km` | 4 | Coordinate/circumferential surface radius of that configuration; **target**. |
| `R_infinity` | `13.80` | `km` | 4 | Effective radius seen at infinity for the same configuration; **derived target**, correlated with M and R. |
| `P_K` with heading superscript d | `0.98` | `ms` | 2 | **Context only; excluded.** Adopted pure-neutron-gas rotation result, not a static npe-mu free-gas TOV prediction. |

All five numbers and their printed digits come from F05 Table 1, p. 30, not
earlier reviews. Its footnote c reads:

> Corresponds to maximum mass before appearance of Σ− hyperons

Footnotes a and b qualify the two APR rows (tabulation endpoint and noncausal
model); neither modifies the free-gas row. Footnote d identifies an empirical
calculation except for the last, adopted value. F05 section 3.1, p. 9, makes
the latter specifically a **“pure neutron gas”** calculation and explains why
the empirical formula was not used. Do not include `P_K` in S10.

The mass/radius values are presented as the paper's free-gas structural
calculation, with no attribution to another structure calculation. The central
density is the model parameter accompanying those outputs. The source does
not publish its full-precision input, constants inventory, scan grid, surface
algorithm, or numerical error budget. Printed precision is not proof of
computational accuracy.

## 4. Meaning of M_max and the central configuration

**Interpretation B is authenticated.** F05 section 3.1, p. 9, says the free-gas
central density is **“slightly below”** hyperon appearance; Table 1 footnote c
qualifies the apparent maximum. This is a source-domain-truncated mass label,
not evidence for a stationary point `dM/d rho_c=0` of an indefinitely extended
npe-mu sequence. There is no support for an unrestricted mathematical maximum,
or for replacing the model by a pure-neutron gas.

The paper does **not** publish exact central `n_B`, the last sequence step, or
an exact distance below the threshold. On an open increasing sequence there
is a supremum but no largest real-valued sub-threshold central density.
“Highest valid model” therefore requires a stated numerical sampling rule;
it cannot silently mean the authors' model or an arbitrary `nextafter` value.
Even monotonic mass growth over the relevant high-density branch is a
proposition for S8 to check, not a TOV result of this audit.

Using the current ratified masses and ideal-gas equations, independent local
calculations give the following (no stellar equation was integrated):

| Local state | n_B [fm^-3] | epsilon [MeV fm^-3] | epsilon/c^2 [g cm^-3] | P [MeV fm^-3] |
|---|---:|---:|---:|---:|
| Published central-density midpoint | 0.604925241328851 | 617.054746418490 | 1.10000000000000e15 | 30.6810920899594 |
| Sigma appearance, mathematical one-sided local limit | 0.617355207966530 | 630.374639464484 | 1.12374486613325e15 | 31.6881733964115 |

The provider's stored ceiling is `0.61735520796653498 fm^-3`; evaluating
strictly **below** it gives `rho=1.12374486613326e15`. The independent root differs
by about `5.2e-15 fm^-3`, compatible with the existing root precision, not a
new ceiling convention. Production rejects requests at/above its stored
ceiling (`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:189-218`, `:230-246`).

The Table-1 midpoint is about 2.01% below the ceiling in baryon density and
2.11% below it in energy density. Its inferred rounding interval maps to
`n_B in [0.602305396323278, 0.607544233003639) fm^-3`; the ceiling is outside
it. Thus `rho_c=1.10e15` and `rho_Sigma=1.123745e15` are conceptually
consistent with a selected sub-Sigma model, but are **not the same state
rounded differently**. The midpoint inversion is our source-constrained
choice, not recovery of an unpublished author input. Fig. 1's `P_c=4.5e34`
also is not that midpoint (`P_c=4.9156528852e34 dyn cm^-2`); no cross-identification
of the figure and table configurations is justified (F05-context).

## 5. Density, pressure, and constant conventions

F05 section 3.1, p. 9, explicitly calls `rho_c` **“central energy density”**.
With Table 1's `g cm^-3` unit it is total local energy density divided by
`c^2`, including rest energies. Section 2.3, p. 6, separately defines the
rest-frame **baryon number** density n and uses `G=c=1`. This rules out treating
`rho_c` as `m_n n_B`, `m_u n_B`, or a sum of baryon rest masses.

For a local energy density E in MeV fm^-3:

```text
E_cgs [erg cm^-3] = E * Units::MEV_FM3_TO_ERG_CM3
rho [g cm^-3]    = E_cgs / c_cgs^2
P_cgs [dyn cm^-2]= P [MeV fm^-3] * Units::MEV_FM3_TO_ERG_CM3

Units::MEV_FM3_TO_ERG_CM3 = 1.602176634e33
c_cgs = GSL_CONST_CGSM_SPEED_OF_LIGHT = 2.99792458e10 cm/s
rho / E = 1.7826619216278975e12
```

This is a derivation from **existing** conversion/constant owners, not a new
literal to add in production (`CompactStar/Units.hpp:56-62`;
`CompactStar/Core/src/TOVSolver.cpp:1505-1519`). The equivalent Zaki chain is
`E * MEV_2_INV_FM * INV_FM4_2_G_CM3`; its measured factor is
`1.7826619216278979e12`, a rounding difference of about `2.7e-16` relative.
For pressure the analogous chain uses `INV_FM4_2_Dyn_CM2`. Reuse the governed
owners and record conversion order in the later generator.

At the Sigma limit, `m_n n_B/c^2` in the same units would give about
`1.03402521388e15 g cm^-3`, and the baryon rest-mass sum about
`1.03400237807e15`; both differ materially from the total `1.12374486613e15`.
At the table midpoint those alternatives are `1.01320592096e15` and
`1.01318407792e15`, versus the required `1.10e15`. These independent consistency
checks reinforce the explicit source definition; they do not choose it by fitting.

There is existing solar-mass conversion debt, **not resolved here**:
the raw TOV output divides grams by GSL's solar mass, whereas `NStar` maps that
mass to km with `Zaki::Physics::SUN_M_KM`
(`CompactStar/Core/src/TOVSolver.cpp:2591-2592`;
`CompactStar/Core/src/NStar.cpp:182-197`;
`CompactStar/Units.hpp:30-35`). The measured km-per-solar-mass values are
`1.4767161818921164` (GSL G,M_sun,c) and `1.4766250380501249` (Zaki), differing
by `6.17e-5` relative. Retain today's conversion boundaries, label oracle
comparisons by convention, and do not “fix” this debt to improve a table match.

## 6. Equilibrium pressure from the active local model

The model is cold, ideal spin-1/2 n,p,e,mu matter with each rest energy included
once (`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:88-99`;
`CompactStar/EOS/src/LocalThermodynamics.cpp:65-145`). Use `k_i=p_F,i c` in MeV,
`b=hbar*c` in MeV fm, `m_i` as rest energy in MeV:

```text
n_i       = k_i^3 / (3 pi^2 b^3)
mu_i      = sqrt(m_i^2 + k_i^2)
epsilon_i = [1/(pi^2 b^3)] integral_0^k_i q^2 sqrt(m_i^2+q^2) dq
P_i       = [1/(3 pi^2 b^3)] integral_0^k_i q^4/sqrt(m_i^2+q^2) dq
P_i       = mu_i*n_i - epsilon_i                 (integration by parts)
P         = sum_i P_i = sum_i mu_i*n_i - epsilon
```

Absent species contribute exactly zero energy/pressure/density. No inactive
derivative is needed. With charge neutrality and baryon accounting,

```text
n_p = n_e+n_mu, n_B=n_n+n_p
sum mu_i*n_i = mu_n*n_B - eta_npe*n_e - eta_npmu*n_mu.
```

For **npe**, `n_mu=0` and only `eta_npe=0` is an active equilibrium equation.
For **npe-mu**, both imbalances vanish. Therefore on either equilibrium branch
`d epsilon_eq/dn_B = mu_n` and `P=n_B*mu_n-epsilon`.
This uses the intrinsic differential and exact constrained coordinates of
ADR-0010 (`docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:126-148`).
It assumes cold chemical equilibrium, neutrality, the stated baryon inventory,
and consistent rest-energy conventions. It is not a general off-equilibrium identity.

For **p-e**, there is no active neutron equilibrium equation:
`epsilon=epsilon_p(n_B)+epsilon_e(n_B)`,
`h_pe=mu_p+mu_e=d epsilon/dn_B`, and `P=n_B*h_pe-epsilon`.
Do not replace `h_pe` by the inactive neutron rest mass below onset.
At neutron onset `h_pe=m_n`; both formulas have the same limit. At muon onset
the added muon term vanishes. Vacuum has `epsilon=P=0`, with
`h_pe -> m_p+m_e` (`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:326-370`).

These are exact mathematical identities. **Numerical implementation should
evaluate stable species pressures**, not subtract two almost equal rest-energy
quantities near vacuum. For `k/m << 1`, the same integral gives
`P_i = n_i*k_i^2/(5m_i) * [1 - 5(k_i/m_i)^2/14 + ...]`.
Use a validated integral/series or cancellation-safe closed form, with truncation
and switching error checked. Today's energy accessor does not provide a pressure
primitive. Even it calls `Evaluate`, which computes a species derivative for
positive n (`CompactStar/EOS/src/LocalThermodynamics.cpp:106-177`): it avoids
the neutral Hessian, but is not strictly a derivative-free implementation.

## 7. P-7 and the smallest equilibrium interface

**Q1 — existing capability is partial, not a total API.**
`EquilibriumAt(n_B)` calls `SolveNpeEquilibrium` and then `EvaluateNpe`; the
N-1 `2^30`-ULP guard throws **before** assembling the energy, so the composite
result cannot supply epsilon or pressure in its refusal window
(`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:373-404`, `:431-451`).
`EquilibriumStateAt` omits that neutral-Hessian request (`:304-320`). A consumer
can use its composition plus the public ideal-fermion value primitives to obtain
epsilon and sum species pressures without H. This is executable today wherever
the **composition root** is certified, including the P-7 probes below.

Live scratch checks against unchanged source:

| Relative n_B above neutron onset | Composition + species values | Composite EquilibriumAt |
|---:|---|---|
| 1e-6 | available | available |
| 1e-7 | available; epsilon=6.91023228979755e-6, P=1.89648759299609e-9, MeV fm^-3 | N-1 response refusal |
| 1e-9 | available | N-1 response refusal |
| 1e-10 | available | N-1 response refusal |
| 1e-12 | composition bracket unresolved | composition bracket unresolved |

The older composition guards also affect immediate neutron/muon neighborhoods
(`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:139-174`, `:250-300`),
and the stored muon threshold has the carried finite-precision N-6 caveat
(`docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:710-715`). Thus
“just use composition” does not establish global totality. The all-species
`IntrinsicChemicalPotentialsMeV` accessor is not a shortcut: it requires every
species positive (`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:458-472`).

**Q2 — recommend a distinct value-only equilibrium capability**, on the
source-specific provider or a thin provider-owned companion, before TOV work.
It implements the existing equilibrium-background role, not a replacement for
ADR-0010's chemical-response authority. No new generic TOV interface is needed
for the recommended table route. Do not make all other local providers implement
a new pure virtual function merely for this source case.

**Q3 — proposed result semantics:**

```text
EquilibriumBarotropePoint
  n_B [fm^-3]
  epsilon, pressure [MeV fm^-3], total rest-energy convention
  n_n,n_p,n_e,n_mu [fm^-3] with baryon and charge accounting
  physical active-set identity and exact threshold/vacuum identity
  provider/model revision, rest-mass/lepton ownership, source-domain provenance
  value validity / numerical resolution status, separate from response availability
```

This is a proposed semantic shape, not an accepted C++ declaration. A vacuum
alternative must preserve the existing positive-density composition invariant;
fractions are defined only for `n_B>0`. No Hessian member or dummy response
dimension is necessary. A separate qualified **barotropic slope** capability
may provide `d epsilon/dP` or `c_s^2` at positive pressure; its vacuum limit
must be typed rather than stored as a fake finite derivative.

Use the ratified species values and branch physics. In threshold-sensitive
value calculations, prefer a stable chemical-potential/Fermi-momentum
parameterization with branch-aware inversion, as independently checked in
section 9, rather than cancellation in `n_B-n_e`. Demonstrate accuracy in the
P-7 window and the narrower composition-resolution bands. If a requested
composition still cannot be certified, return explicit value/composition
resolution status; **do not claim totality by catching an exception and
substituting a p-e state, a floor, or a neighbor**. A generated barotrope can
cover a numerically unresolvable query interval only using independently
certified bracketing values, exact threshold metadata, and an explicit bounded
interpolation error, not an unexplained hole or silent branch relabeling.
The next task must demonstrate this coverage before integrating a star.

N-1 remains unchanged and still rejects inaccurate **chemical responses**.
No caller may infer an available H from available values. This preserves the
ratification's mandatory P-7 planning boundary
(`docs/validation/PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:118-125`).

## 8. Actual ordinary-TOV interface and caller map

The consumed object is **`Core::EOSTable` owned by concrete `Core::TOVSolver`**,
populated from a tab-separated file. There is no abstract material EOS pointer
or virtual/callable epsilon(P) injection in this path. `CompactStar::Model`
is an abstract **table producer** (`EDens`, `Press`, `EOSRow`, `FindEOS`,
`ExportEOS`), not the object queried by the TOV RHS
(`CompactStar/EOS/Model.hpp:54-85`; `CompactStar/EOS/src/Model.cpp:24-61`,
`:84-97`). `CompOSE_EOS` is another producer, deriving from `Core::Prog`,
not an interchangeable TOV callback (`CompactStar/EOS/CompOSE_EOS.hpp:74`;
`CompactStar/EOS/src/CompOSE_EOS.cpp:429-430`). Do not activate the older
unratified free-gas/RMF producers to bypass the ratified local model.

Actual paths, traced through call bodies:

```text
main/Test/ns_build_main.cpp:81
main/Test/spin_therm_evol_2_main.cpp:124
  -> NStar::SolveTOV_Profile(file, target_M)
     -> new TOVSolver -> ImportEOS(file, true)
     -> SolveToProfile(target_M) -> SingleStarSolveToTOVPoints(ec)
     -> NStar::BuildFromTOV(points, labels)

main/Test/spin_therm_evol_main.cpp:57-67
main/Test/tov_debug_main.cpp:165
  -> ImportEOS -> Solve(Axis)
     -> SingleStarSolveToTOVPoints(Axis[idx])
     -> Append -> SurfaceIsReached -> FinalizeSurface
     -> Sequence::Add -> Reset -> ExportSequence

tests/core/tov_reference_cmf.cpp:205-214
  -> ImportEOS -> SingleStarSolveToTOVPoints(ec, points)
```

Implementations: `CompactStar/Core/src/NStar.cpp:1254-1318`;
`CompactStar/Core/src/TOVSolver.cpp:1807-1854`, `:2516-2800`.
The example `polytrope`, `coulatt`, `sig_omega`, and `Table_5-8_Glenn` paths also
import files before `Solve` (`main/Examples/polytrope.cpp:64-66`,
`main/Examples/coulatt.cpp:48-50`, `main/Examples/sig_omega.cpp:62-64`,
`main/Examples/Table_5-8_Glenn.cpp:48-63`). Names and commented-out calls were
not counted as alternate radial implementations.

The downstream `StarBuilder::BuildFromSequence` call in
`main/Test/spin_therm_evol_main.cpp:71` reads exported sequence/profile files,
selects a mass index and may blend neighboring profiles; it does not solve a
new TOV star (`CompactStar/Core/src/StarBuilder.cpp:93-212`). Do not use a
blended profile as the first Table-1 structure result.

| Requirement | Actual boundary / owner | Capability classification for Track R |
|---|---|---|
| Independent variable | Pressure for every visible EOS spline; TOV integrates in radius in cm | **NEEDS ADAPTER** from n_B/chemical-potential parameter to increasing P |
| epsilon(P) | File column 1 is epsilon/c^2 in g cm^-3; `GetEDens` reads its spline | Local epsilon **ALREADY PROVIDED DIRECTLY** outside response failures; continuous value path **NEEDS NEW PRODUCTION API**; TOV bridge **NEEDS TABULATION** |
| P | File column 2 in dyn cm^-2, not MeV fm^-3 or km^-2 | **DERIVABLE ANALYTICALLY**, stable pressure capability **NEEDS ADAPTER/API** |
| n_B(P) | File column 3, fm^-3; `GetRho` is baryon number density despite its name | Composition supplies n_B **ALREADY PROVIDED DIRECTLY**; **NEEDS TABULATION** for this consumer |
| Species columns | Optional extra columns; physical content is dimensionless Y_i, not n_i | **NOT NEEDED FOR STRUCTURE** RHS; proposed table should include them for coherent NStar composition/provenance. **NEEDS ADAPTER** dividing n_i by n_B at positive rows |
| d epsilon/dP | `GetEDensDeriv`: same energy interpolant differentiated, multiplied by c^2 | **NOT NEEDED FOR STRUCTURE** equations; nevertheless **required finite by current primitive publication**; analytic reference **DERIVABLE ANALYTICALLY**, current route supplies it from the table |
| Central state | Public primitive accepts epsilon_c/c^2 in g cm^-3; clamps then Brent-inverts the same epsilon(P) | **NEEDS ADAPTER** from requested n_B; verify achieved central epsilon/P |
| Off-equilibrium H | Never requested by the TOV path | **NOT NEEDED FOR STRUCTURE** |

Import stores columns without unit conversion or fraction normalization; header
strings do not enforce units. It constructs separate `gsl_interp_steffen` splines
for epsilon, n_B and every extra column, with a shared pressure accelerator
(`CompactStar/Core/TOVSolver.hpp:154-169`, `:508-545`;
`CompactStar/Core/src/TOVSolver.cpp:603-718`). Missing/malformed input handling
is not a complete scientific validator: validate finite aligned rows, units,
strictly increasing pressure, positive densities, labels and fractions **before**
import. Null species splines otherwise produce zeros (`:1567-1582`).

The RHS uses `GetEDens` only, with GSL G and c; it neither uses n_B, fractions
nor a thermodynamic derivative (`CompactStar/Core/src/TOVSolver.cpp:1490-1524`).
However, output-node assembly explicitly evaluates and checks **finite dedp**
as well as energy/n_B/species (`:2580-2592`). The older descriptive claim that
TOV does not require the derivative must therefore be narrowed to the
**equations**, not the current successful-return contract
(`docs/architecture/CURRENT_ARCHITECTURE.md:208`).

`TOVPoint` uses radius km, mass solar masses, pressure dyn cm^-2, energy g cm^-3,
nu-prime cm^-1, n_B fm^-3 and fractions; `NStar` converts to profile r,m in km,
epsilon,P in km^-2, nu-prime km^-1
(`CompactStar/Core/TOVSolver.hpp:317-342`;
`CompactStar/Core/src/NStar.cpp:178-208`, `:767-779`).
Both NStar builders also compute total baryon number and first-order I as
existing side effects (`CompactStar/Core/src/NStar.cpp:368-393`, `:660-706`).
They are not new scientific deliverables or validation targets. Prefer raw
primitive M/R for static sequence exploration; use ordinary NStar construction
only for the required path/interface comparison and existing metric carrier.

## 9. Equilibrium slope and threshold regularity

Define a **scalar ideal-gas compressibility**, not a chemical Hessian or stellar
susceptibility:

```text
s_i = dn_i/dmu_i = mu_i*k_i/(pi^2 b^3)       (positive-density species)
s_i -> 0                                   at first appearance
L = s_e + s_mu                              (s_mu=0 on npe)
C = dn_B/dmu_n = s_n + s_p*L/(s_p+L)
d epsilon/dP = mu_n*C/n_B
c_s^2/c^2 = dP/d epsilon = n_B/(mu_n*C)
```

Derivation: choose the common lepton potential t, set `n_p=n_e(t)+n_mu(t)`,
`mu_n=mu_p(n_p)+t`. Then `dn_p=L dt`, `dmu_p=L dt/s_p`,
`dn_n=s_n dmu_n`, yielding C above. `dP=n_B dmu_n` and
`d epsilon=mu_n dn_B` complete the result. On p-e use
`C_pe=1/(D_p+D_e)=s_p*s_e/(s_p+s_e)` and replace mu_n by `h_pe`.
At neutron onset `s_n=0` the formulas coincide; at muon onset `s_mu=0`
they reduce to npe. Zero s_i here is the true one-sided derivative of
**density with respect to potential**; it is not a fabricated finite `D_i`
or padded chemical response. No inversion of the ADR-0010 H is performed.

At a positive-density threshold the newly appearing species satisfies
`n_new proportional to (mu_new-m_new)^(3/2)`. Consequently:

- epsilon(n_B) is C2, generally not C3; its second derivative has a square-root
  cusp in its variation, with the next derivative singular.
- P(n_B), epsilon(P), n_B(P), and composition along the equilibrium barotrope
  are C1, generally not C2. Their first derivatives are continuous but have
  nonanalytic approach to the limiting value. There is no latent heat or
  pressure/energy jump and no density delta-function shell at these onsets.
- In the vacuum limit `P proportional to n_B^(5/3)`,
  `epsilon approximately (m_p+m_e)n_B`; hence `epsilon(P) proportional to P^(3/5)`.
  It is C0 at P=0 and its first derivative diverges there;
  `dP/d epsilon -> 0`. Do not supply finite dedp at vacuum.

These regularity conclusions are derived for the ratified ideal-gas branch
model, not inferred from the changing dimension of a software response type.
Independent 80-digit quadrature checks at the exact mathematical onsets gave:

| Threshold | epsilon [MeV fm^-3] | P [MeV fm^-3] | c_s^2/c^2 |
|---|---:|---:|---:|
| neutron | 6.91023159858475e-6 | 1.89648750263179e-9 | 3.87437706631595e-4 |
| muon | 460.159248800949 | 19.6111178505123 | 0.0657776313936194 |

For relative offsets `1e-8,1e-16,1e-24` in common lepton potential, epsilon and
P from both sides approach these values. At the smallest offset, maximum
relative epsilon/P discrepancies are respectively `3.56e-24/5.02e-24` at
neutron onset and `1.77e-24/2.74e-24` at muon onset. The slope converges much
more slowly at neutron onset: relative dedp differences `3.02`, `3.02e-4`,
`3.02e-8` across those offsets. This is physical stiffness, not a jump. Three
interior branch probes compared the analytic slope with independently
differentiated energy/pressure quadratures to about `1.1e-79` relative.
The integral pressure identities agree to better than `4e-77` at the recorded
threshold/central probes. Those are high-precision algebra checks, not promised
double-precision production accuracy.

Monotonicity follows from C>0 and positive h,n_B. Causality also follows:
the relaxed-composition stiffness is no larger than the fixed-fraction
stiffness, and `n_i D_i=k_i^2/(3mu_i)<=mu_i/3`. Thus
`0 < dP/d epsilon <= 1/3` on the positive-density equilibrium branches.
S4 must test the **implemented representation**, since interpolation can still
violate this property. No chemical Hessian availability enters this proof.

For the **table route**, the production derivative must remain the analytic
derivative of **the same interpolated epsilon(P)**, as it is today
(`CompactStar/Core/src/TOVSolver.cpp:1070-1103`). The analytic expression above
is the independent reference for its convergence, not a second derivative
silently attached to a different EOS. GSL's [Steffen specification](https://www.gnu.org/software/gsl/doc/html/interp.html)
is piecewise cubic and C1, with monotonicity between supplied points. It is not
a natural C2 spline despite stale nearby comments in `TOVSolver.hpp:534-538`.

## 10. Architecture alternatives and recommendation

| Strategy | Fidelity/error and thresholds | Reuse/invasiveness | Validation/provenance and verdict |
|---|---|---|---|
| A. Direct analytic pressure/density query | Eliminates table error; still needs stable equilibrium inversion and correct one-sided limits. Vacuum remains singular in epsilon(P). | No supported injection into today's TOV; requires a new EOS-consumer contract and coordinated derivative/domain/surface work. | Attractive later and as an independent oracle; excessive production change for the first reproduction. |
| B. Governed generated table from a value-only ratified-provider capability | Sampling/interpolation error measurable independently. Include threshold rows and refine both sides; do not insert crust or smooth physical onsets. | Reuses `ImportEOS`, Steffen, canonical primitive, derivative owner and ADR-0009 unchanged. | **Recommended first route.** Generator/config/constants/source hashes and accuracy record travel together. |
| C. Thin analytic adapter implementing an existing EOS abstraction | A `Model` subclass can generate rows, but cannot be queried polymorphically by TOV. | As a Model producer it reduces to B. As a supposed analytic TOV callback the necessary abstraction does not exist. | Do not present a producer base class as an available runtime interface. No need to subclass merely to reuse an unsuitable uniform log grid. |
| D. Provider-owned values plus a small dedicated exporter using existing file schema | Same physical equations as B, with explicit active-set, endpoint and unit checks rather than a new broad hierarchy. | Small additive local capability and a generator; no TOV change. | Preferred concrete realization of B; independent analytic EOS reserved for S12. |

Generate in a monotone physical parameter/n_B with logarithmic low-density
coverage and localized threshold refinement; export in strictly increasing P.
Predeclare accuracy checks at off-grid points. Include physical onset rows
once, with exact inactive zeros and branch metadata outside the numeric file.
Sample fractions accurately, including nonmonotone individual fractions; do
not assume sum/neutrality identities survive separate cubic interpolants exactly.
Check those identities off-grid too. Export enough digits for round-trip values
(at least `max_digits10` for double); do not inherit an uninspected table-writer
default. Store no regression TSV as a runtime EOS authority.

Finite-resolution interpolation approximates a C1 threshold; that is acceptable
only with explicit threshold nodes and convergence against the analytic
barotrope, especially its rapidly changing slope. Smoothing/reclassifying an
onset over an arbitrary width, interpolating across an unbounded unresolved
interval, or asserting cubic convergence without measurement is unsafe.
Do not insert a vacuum row in today's import merely to advertise zero pressure:
the primitive requires a positive floor and positive p_cut
(`CompactStar/Core/src/TOVSolver.cpp:2538-2551`).

No new architectural owner is selected: local matter values remain with the
provider, unit/schema export with its adapter, interpolation and structure
with TOV, and stellar postprocessing with NStar. ADR-0010 already separates
equilibrium backgrounds from off-equilibrium response. Accordingly this plan
requires **no new ADR or blocking owner question**. If future evidence instead
requires an analytic TOV injection or changing the accepted surface definition,
stop that implementation at a proposed owner question; this audit does not
authorize either alternative.

## 11. Surface: accepted p_cut versus the source's zero-pressure limit

F05 section 3.1 applies the free gas to the **“whole star”**, with no separate
crust supplement. The ratified ideal gas has P,epsilon -> 0. The accepted
CompactStar surface is instead exactly

```text
p_cut = max(1e-15 * P_c, eos_tab.pre.front())     [dyn cm^-2]
```

The number `1e-15` is fixed, not a configurable tolerance. There is no universal
absolute default: it depends on central pressure and the imported floor. The
implementation is `CompactStar/Core/src/TOVSolver.cpp:1329-1332`; the header's
`1e-5` comment at `CompactStar/Core/TOVSolver.hpp:589-594` is stale. The
DS(CMF) floor `3.351885e25 dyn cm^-2` is **not** a Track-R prescription
(`docs/SCIENTIFIC_INVARIANTS.md:289-297`). At the Table-1 central-density
midpoint, the smallest permitted cutoff, once the generated floor is lower,
is approximately `4.9156528852e19 dyn cm^-2`; its local p-e density is
`6.2813807792e-12 fm^-3`.

The event locator integrates the terminal segment in pressure coordinates,
`dr/dP=1/f_P`, `dm/dP=f_m/f_P`, using the same RHS; the final node is exactly
at p_cut and successful publication requires `SURFACE_REACHED`. Trial-only
epsilon is evaluated at max(P,P_floor), retaining trial P in pressure terms
(`CompactStar/Core/src/TOVSolver.cpp:1490-1524`, `:2479-2513`, `:2571-2596`).
This is not a true-P=0 treatment. At P=0 both the finite-density domain guard
and finite dedp contract fail, and dr/dP is singular. A separate enthalpy
coordinate reaches the analytic zero-enthalpy surface cleanly; implementing
it as a production surface would reopen ADR-0009 and is not recommended here.

**Estimate, not a computed star:** define dimensionless local enthalpy
`H(P)=integral_0^P dP/(epsilon+P)=ln[h(P)/(m_p+m_e)]` in the p-e tail.
At the minimal cutoff H=1.3073642065e-5. With the printed M,R as scale inputs,

```text
delta R approximately H_cut * R^2 * (1-2m/R) / m
                  = 0.0019948 km, about 2 metres.
```

This is 40% of the radius rounding half-width, so the default cannot honestly
be declared negligible without further evidence. Reducing the table floor
below `1e-15 P_c` will reach a **plateau**, not demonstrate convergence to
P=0. A 2-metre missing layer can affect whether a model rounds to 12.77.

Required S7 falsifier, at fixed achieved central state:

1. Generate nested, already EOS-converged tables with floors giving effective
   `p_cut/P_c` near `1e-11,1e-12,1e-13,1e-14,1e-15`, then lower the table floor
   further to demonstrate the fixed-ratio plateau. Record actual p_cut,
   endpoint P and completion status, M,R,R_infinity at every step. Resolve
   interpolation error separately; a moving last radial sample is not a surface.
2. Establish the remaining tail using independent enthalpy/pressure integration
   to P=0 **or controlled analytic bounds**, with error substantially below the
   publication interval. Keep canonical R(p_cut) and inferred R(0) distinctly
   named; do not modify the stored star or relabel it as a vacuum-surface star.
3. For an analytic bound in consistent geometric units, let r0,m0,epsilon0,P0
   be the cutoff values. Monotonicity gives the constant-mass upper radius
   `ru = 2m0/[1-(1-2m0/r0)*exp(2H0)]`, provided its denominator is positive.
   Then `dm_u = (4pi/3)*epsilon0*(ru^3-r0^3)` bounds the tail mass. With
   `mu=m0+dm_u`, use
   `g_u=(mu+4pi*ru^3*P0)/[r0*(r0-2mu)]` if r0>2mu, giving
   `r0+H0/g_u <= R(0) <= ru`. These inequalities follow directly from
   `-dH/dr=(m+4pi*r^3*P)/[r(r-2m)]`; they bound rather than assume the
   thin-tail approximation. Require the bracket width and corresponding mass/
   R_infinity uncertainty to meet the internal budget, or use the independent
   tail integrator. Propagate raw-GSL/profile conversion conventions separately.
4. Check that the zero-pressure estimate/bounds agree across the cutoff family.
   The asymptotic omitted radius scales as P_cut^(2/5), not linearly in P_cut.
   Any extrapolation must be checked against the enthalpy bound or independent
   oracle, not fitted solely to force 12.77. If those checks fail, stop; do not
   lower the production coefficient or change ADR-0009 automatically.

## 12. Central sequence and source-comparison rule

**Recommend D, constrained by C:** reproduce the simultaneous printed ranges
without claiming an exact unpublished central state. First invert
`rho_c=1.10e15` using the fixed local EOS; report that predeclared midpoint
star and all discrepancies. Independently sample/refine the full printed
central-density interval `[1.095,1.105)e15`, and show the image of that interval
in (M,R,R_infinity). If only a subset agrees, report the subset and the
midpoint outcome. One **same** central state must satisfy all output intervals;
separate matches at separate densities are not a reproduction. Using rho to
select the star is an input constraint, not a fourth independent prediction.
Never tune particle masses, degeneracies, interactions, threshold, surface
cutoff, or interpolation to fit these numbers.

Additionally scan the relevant high-density branch below Sigma to check mass
behavior and describe the domain truncation. Do not claim the unrestricted
sequence maximum, a dynamically established stability boundary, or that an
arbitrarily selected last grid point is the authors' last model. A finite
sample's largest mass is labeled as such. If the printed box has no common
intersection after convergence, report non-reproduction/source ambiguity and
seek scientific adjudication; do not widen tolerances automatically.

Use fixed-central-density calls to the canonical primitive, not
`SolveToProfile(0.62)`: the latter searches 25 coarse logarithmic samples,
uses the first rising mass bracket, and stops at absolute `1e-4 M_sun`
(`CompactStar/Core/src/TOVSolver.cpp:2654-2715`, `:2742-2800`). A target-mass
fit would hide a central-density discrepancy. Very low-density free-gas
configurations may exceed the default 70-km integration radius; failure there
is not evidence of an absent physical sequence. Bound the structural scan to
the high-density target branch, and classify any broader exploratory failures.

The current primitive clamps to `[10 epsilon_min, 0.999 epsilon_max]` by
default, then performs Brent pressure inversion with `p_of_e_prec=1e-4`
(`CompactStar/Core/TOVSolver.hpp:657-669`;
`CompactStar/Core/src/TOVSolver.cpp:1173-1246`, `:2538-2549`). Therefore:

- Validate requested and **achieved** central state; no silent clamp is acceptable
  as the requested source configuration. Record `GetInitPress()` and the first
  point's energy, and bound inversion error using the analytic EOS.
- Put the generated upper row below Sigma, but high enough that
  `0.999 epsilon_max > 1.105e15 g cm^-3`. This is possible with the existing
  `rho_Sigma approximately 1.123745e15`; no super-Sigma padding is needed.
- The exact Sigma supremum is inaccessible with a sub-Sigma table plus this
  clamp. Near-ceiling diagnostic stars must state the actual ceiling margin;
  reducing sequence spacing does not remove it. A supremum solve is not the
  Table-1 acceptance rule and does not justify changing the clamp.
- Pressure-inversion and ODE tolerances are not public runtime controls here.
  Check achieved-state residuals and use an independent high-accuracy oracle;
  if default inversion error prevents the internal target, stop for a separately
  authorized numerical-control change. Do not call `GenTestSequence` as a
  hidden tolerance setter: it mutates `p_of_e_prec` (`:3191`, `:3228`).

## 13. R_infinity and available metric information

F05 section 2.1, p. 3, defines `R_infinity = R exp(-Phi_s)` after eq. (6).
Schwarzschild exterior matching gives

```text
R_infinity = R / sqrt(1 - 2GM/(R c^2)).
```

With the printed M=0.62 and R=12.77, existing Zaki profile conventions give
`13.7974183735 km`; the raw GSL combination gives `13.7974896494 km`.
Both round to **13.80 km**, so the printed values are consistent. The roughly
0.071-metre difference between these checks is existing constant-convention
debt, not numerical integration error or evidence of a different FR2005 radius.

`NStar::BuildFromTOV` sets the surface lapse `z_surf=exp(nu(R))`, where
`nu(R)=0.5 ln(1-2m/R)`; `StarProfile::RadiusSurface()` and
`ExpNuSurface()` already supply `R/exp(nu_s)`
(`CompactStar/Core/src/NStar.cpp:385-393`, `:905-928`;
`CompactStar/Core/StarProfile.hpp:354-367`). Here `z_surf` is the lapse,
**not** the astronomical redshift z. No dedicated R_infinity accessor was
found on this audited path; no duplicate metric implementation is needed.
The `Append/FinalizeSurface` sequence path still leaves mirror surface scalars
zero (ADR-0005 M1), so use its actual final radius/metric columns or the
BuildFromTOV path, not division by the unset mirror lapse
(`docs/architecture/CURRENT_ARCHITECTURE.md:299-315`;
`CompactStar/Core/src/NStar.cpp:658-708`).

## 14. Publication precision versus numerical convergence

Assuming ordinary rounding to nearest at the printed decimal place, the
comparison bins are:

| Quantity | Printed | Inferred rounding interval | Role/ambiguity |
|---|---:|---|---|
| M/M_sun | 0.62 | [0.615, 0.625) | Truncated source-domain configuration; two significant figures. |
| rho_c [g cm^-3] | 1.10e15 | [1.095e15, 1.105e15) | Unknown precise author input; constrain central state. |
| R [km] | 12.77 | [12.765, 12.775) | Source zero-pressure limit versus computed finite-cutoff radius must be distinguished. |
| R_infinity [km] | 13.80 | [13.795, 13.805) | Correlated with M,R and metric convention. |

Endpoint inclusion depends on the unpublished rounding rule; do not declare
success on an exact tie. These are **precision-derived bins, not published
statistical confidence intervals or guaranteed author error bars**. If a
converged result misses, retain the residual and uncertainty; lack of a source
error budget does not license broadening a bin. Require a common central-state
solution whose full computed uncertainty fits inside all applicable bins.

Proposed internal error budget, to be established before S10: aggregate numerical
uncertainty no larger than **1% of each half-bin**, i.e. `5e-5 M_sun`,
`5e10 g cm^-3` achieved-central-density uncertainty, and `5e-5 km` (5 cm)
for each radius. This is an explicit engineering separation from publication
precision, **not a source-prescribed tolerance**. Divide/track this aggregate
budget among equilibrium root/pressure evaluation, inversion, EOS interpolation,
radial integration/output partition, central scan and surface-tail estimation;
passing each separately at the full budget does not establish the sum bound.
The finite-cutoff **bias** (about 2 m by the scale estimate) must be bounded or
corrected in the separate comparison estimate; it is not an internal error
to conceal inside this 5-cm budget. Likewise known constant-convention offsets
must be recorded separately, not absorbed into a convergence tolerance.

For radial output partitions retain ADR-0009's existing tighter checks
`relative delta M <=1e-9`, `relative delta R/lapse <=1e-8` and surface pressure
residual `<=1e-8`, using e.g. 2500,5000,10000,20000,40000 samples. The adaptive
RK8PD driver's current tolerances are `1e-10` absolute/relative, initial radial
step 0.1 cm, and center 1 cm; those settings are not error certificates
(`CompactStar/Core/src/TOVSolver.cpp:2550-2558`, `:2499-2509`;
`CompactStar/Core/TOVSolver.hpp:599-604`). Require EOS and sequence refinement
on nested levels with measured error, not presumed cubic/second-order scaling.
Compare same achieved central states; record numerical floor or nonmonotonic
convergence instead of picking the level closest to the paper.

## 15. Proposed independent validation ladder (not executed as structure tests)

| ID | Required evidence | Defect falsified |
|---|---|---|
| S1 | Independent phase-space species pressures versus sum(mu_i n_i)-epsilon and n_B*h-epsilon on each active branch; small-momentum series and roundoff probes; pressure agrees where H deliberately refuses | Missing rest masses/leptons, factor of 3 or degeneracy, wrong p-e conjugate, catastrophic cancellation, hidden Hessian dependency |
| S2 | Analytic equilibrium slope above versus independently differentiated quadratures; compare the imported interpolant's own derivative at off-grid points and under refinement | Chemical Hessian confused with barotropic slope, missing equilibrium chain rule or c^2 conversion, derivative from a different interpolant |
| S3 | One-sided epsilon/P/slope limits at both exact species onsets, vacuum limits, physical classification plus P-7 and composition-resolution probes; prove any interpolation bridge's error | Artificial EOS gap/jump, floors, smoothing, threshold branch substitution or nonexistent inactive derivative |
| S4 | Positive P and epsilon, strictly increasing P(n_B), positive finite slope on P>0, 0<dP/d epsilon<=1/3; baryon/charge/fraction identities, including off-grid | Reordered/duplicate table rows, negative pressure, interpolation overshoot, acausality or inconsistent composition |
| S5 | Same generated table and achieved epsilon_c through primitive plus NStar(points,labels) and supported sequence/finalization path; same raw M/R and metric with documented mirror-scalar distinction | Unit/schema/import mismatch, different radial owner, lost species/metric or stale sequence state |
| S6 | ADR-0009 radial partition matrix, complete event status, independent tighter integration check; center truncation bounded or independently reduced in oracle | Radius set by sampling, incomplete-star publication, unresolved adaptive/center error |
| S7 | Section 11 effective-cutoff sweep, fixed-ratio plateau, controlled zero-pressure tail bounds or independent enthalpy oracle, uncertainty propagation | Truncated radius falsely called P=0, missing no-crust envelope, convergence claimed from a fixed-cutoff plateau |
| S8 | Refine central density interval and nearby high-density mass scan; verify achieved state and clamp margin, invariant common matching subset | Rounded mass reverse-fit, last-grid-point “maximum”, hidden central clamp/inversion error, mesh-dependent acceptance |
| S9 | R/lapse identity at finite cutoff and separate zero-pressure estimate; distinguish raw GSL from profile Zaki convention | Redshift versus lapse confusion, sign/factor-of-two error, unit mismatch or duplicate metric authority |
| S10 | A single converged central state inside the rho bin has M,R(0),R_infinity(0) uncertainty intervals inside their printed bins; report midpoint and subset | Matching different observables at different states, overstated publication precision or EOS tuning |
| S11 | Reject n_B>=stored Sigma ceiling before generation/solve; no upper padding; audit requested, clamped and achieved central densities and all table rows | Out-of-source-domain star, exact-threshold substitution or extrapolated hyperon-free state |
| S12 | Small test-side independent static TOV oracle, preferably dimensionless enthalpy with direct phase-space EOS values, independently coded RHS and center handling, no production splines/ODE calls; compare raw cgs M/R at fixed central state, then known profile mapping | Common algorithm, units or surface defect that self-comparison/regression cannot reveal |

S5 is path-equivalence evidence, not independent physics validation. S12 must
not merely wrap the canonical ODE or interpolate the same generated table.
For each major falsifier, introduce a test-side wrong factor, sign, branch,
surface or unit fixture sufficient to demonstrate sensitivity without changing
ratified production sources. No test or baseline is added by this audit.

## 16. Exact next implementation scope and stop conditions

The next bounded **static structure** task should:

1. Add and independently validate the value-only equilibrium capability and
   stable pressure evaluation described in sections 6-7. Prove threshold and
   P-7 value coverage without weakening N-1 or ADR-0010. Existing local
   response tests remain unchanged in meaning.
2. Add a deterministic, provenance-bearing table generator/adapter for the
   existing TOV file schema, explicit positive surface floor and sub-Sigma
   upper domain. Validate values, fractions, slope convergence and imported
   units before running structure. Generated EOS artifacts are outputs with
   recorded provenance, not newly accepted scientific baselines.
3. Reuse the unmodified canonical TOV primitive for fixed-central-state
   static M/R checks and the source interval search. Add only the test-side
   independent structure/surface oracle and bounded comparison estimator
   needed for S1-S12. Keep finite-p_cut canonical profiles separate from
   zero-pressure comparison estimates and verify the existing NStar metric.
4. Report the common published-bin reproduction or explicit mismatch, all
   convergence evidence and remaining uncertainty. No baseline or ratification
   is conferred automatically by a match.

The next task must receive its own exact authorized file list and validation
requirements. It must stop before any needed TOV behavior/ownership change,
cutoff-coefficient change, above-Sigma extrapolation, unconstrained central
fit, silent response-policy weakening, or broadened source error allowance.
The missing precise author central state is a disclosed comparison limitation,
not a blocking owner decision under the interval-based rule. No unresolved
owner question blocks **beginning** this plan; unsuccessful later evidence
would require a new scientific decision rather than an automatic workaround.

Explicit non-targets: adopted `P_K=0.98 ms`; Fig. 1 rotational compression;
Fig. 8 thermal result; particle-number integrals/response; Btilde, paper Z/W;
Urca rates; chemical evolution; APR/BPAL; DS(CMF) off-equilibrium physics;
superfluidity and BNV. Existing NStar construction side effects are not
validation or acceptance of any of those quantities.

## 17. Audit evidence and delivery boundary

Fresh local diagnostics were compiled directly from the unchanged two local
thermodynamics translation units, using Apple clang 17.0.0, C++17, and the
vendored arm64 Zaki archive. Archive SHA-256:
`3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef`.
Independent mpmath 1.3.0 calculations used 80 decimal digits, the **exact
binary values of the existing production mass/hbar-c constants**, independently
integrated phase-space energy/pressure, and a chemical-potential parameterization.
Their Sigma evaluation is a local one-sided analytic limit, not a production
request outside the source domain. No TOV solver was run and no whole-star
oracle was implemented during this audit.

External scratch evidence (not repository artifacts):
`/Users/keeper/.codex/diagnostics/trackr-structure0-20260905/`:
`values.cpp`, `values.txt`, `independent.py`, `independent.txt`,
`continuity.py`, `continuity.txt`, and `SHA256SUMS.json`.
These provide fresh numerical evidence for sections 4-7, 9, 11 and 13.
The first compile attempt used a nonexistent Homebrew GSL include path; it
failed before execution. Retrying with `/opt/local/include`, authenticated
from canonical `build/CMakeCache.txt:341`, succeeded. No installed project
dependency was changed; mpmath was installed only in a scratch virtual environment.
The PDF rendering warnings did not prevent visual inspection of R95 eq. (3).

Repository validation for this documentation-only task is `git diff --check`,
source/line-reference checks, and exact changed-path/protected-tree checks.
Production, tests, baselines, data and shared literature remain unchanged.
All 22 entries in the shared literature checksum manifest were verified from
`/Users/keeper/Documents/CompactStar`; the three directly read PDF hashes also
match the entry authentication. `literature-integrity.json` records the check.
Existing CTest suites are **not rerun**; historical passes in local ratification
are corroborated records, not fresh results of this audit.

Only this audit and a concise roadmap status pointer are changed. Commit:
`docs: audit Track-R free-gas whole-star interface`; non-force push of the new
branch; no merge. The delivery commit and remote equality are reported after
commit so this document does not pretend to contain its own hash.

NO WHOLE-STAR FREE-GAS REPRODUCTION, STELLAR SUSCEPTIBILITY,
PAPER B/Z/W, ROTOCHEMICAL EVOLUTION, APR/BPAL, DS(CMF)
OFF-EQUILIBRIUM PHYSICS, SUPERFLUIDITY, OR BNV IS IMPLEMENTED
BY THIS AUDIT.

**TRACK-R FREE-GAS WHOLE-STAR STRUCTURE INTERFACE AUDITED — IMPLEMENTATION READY**

**Exactly one recommended next action:** open the bounded Track-R Structure-1
implementation/validation task specified in section 16. Do not begin it automatically.
