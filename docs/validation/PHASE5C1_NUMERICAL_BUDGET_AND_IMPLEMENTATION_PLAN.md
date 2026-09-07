# Phase-5C-1P — numerical budget and implementation plan

**Test/oracle/planning work only. Production G/Z/W is NOT IMPLEMENTED.**
This record proposes numerical implementation policy under accepted ADR-0013;
it does not ratify new physical semantics or claim that precursor tests validate
a production chemical path. Final execution and disposition are recorded below.

## 1. Authenticated entry and authority

| Item | Authenticated value |
|---|---|
| Canonical checkout | `/Users/keeper/Documents/CompactStar/repo/CompactStar` |
| Canonical local/origin/live master | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Accepted preflight | `54ec7abac38fa0a32c5fb3a82e424b496361966a` |
| Accepted ADR ratification / planning starting SHA | `4780121f21010374da2eb50898e90a067795b6e5` |
| Ratified branch, local/upstream/live | `analysis/phase5c-corrected-chemical-coefficients-preflight`, equal to ratification SHA |
| Planning branch | `analysis/phase5c-numerical-budget-planning` |
| Fresh worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-numerical-budget` |
| Parents | ratification → preflight → canonical master, exactly as above |
| Freshness | branch absent locally/remotely; worktree path absent before creation; fresh HEAD and clean tree authenticated |

Entry evidence included `git status --porcelain=v2`, `git worktree list --porcelain`,
`git branch -vv`, and `git log --graph --oneline --decorate -8`, plus live
`git ls-remote` comparisons. Canonical master remains intentionally unmerged.

Authority: `AGENTS.md`, `GOVERNANCE.md`, ADR-0010, accepted ADR-0013,
`PHASE5C0_CORRECTED_CHEMICAL_COEFFICIENT_PREFLIGHT.md`, its ratification record,
ADR-0011, `PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md`, scientific invariants,
roadmap, architecture, and the Phase-5A-2/3/4/5 implementation, correction,
independent-review and local-ratification records. Historical status statements
in those records do not override the accepted entry status.

The owner-supplied `PHASE5C0_INDEPENDENT_REVIEW.md` was recovered from the
attachment in the conversation **Rotochemical Heating Implementation** and read.
Its curved-star derivation and items 35, 58–63, 76, 87–98 informed the experiments.
Attachment SHA-256:
`98811606dd849ca34312dcfe172b853324de5d00097181c37c1e9864a5e3d9e4`.
The attachment is evidence, not a replacement authority. The shared literature
`README.md` and `catalog.tsv` authenticate R2006 as the corrected-charge authority.
No literature or source-data bytes were changed or acquired.

The accepted `G_y` basis, integration-before-baryon-reduction, redshift, Z/W,
active dimensions, full singular response, no pseudoinverse, no second projection,
and A18 benchmark identity remain fixed (ADR-0013 §§3–7). No scientific defect
in those contracts was found by this planning work.

## 2. Quadrature decision

| Candidate | Branch/onset handling | Reproducibility / cost | Error and authority assessment |
|---|---|---|---|
| A: raw profile trapezoid | Integrand kink lies inside a cell unless explicitly split; endpoint H may be unavailable | Fixed nodes, lowest cost | Useful diagnostic; cannot own response integration merely because it follows stored nodes |
| B: composite Simpson | Uniform-grid Simpson only on admissible smooth equal-spacing triples; nonuniform SciPy variant is separately labelled | Fixed nodes, low cost | Variable last step and threshold-straddling triples defeat an assumed fourth order; no automatic acceptance by name |
| C: fixed segment GL | Respects every profile interpolation knot, but physical onsets can remain inside cells | Fixed order and node policy, moderate cost | Smooth-cell high order is useful; an unsplit square-root onset invalidates its smooth remainder hypothesis |
| **D: onset-aware segment GL** | Authenticated branch events split cells; one-sided maps regularize continuous onsets; active embeddings are explicit | Fixed deterministic order/partition, reproducible accumulation; cost proportional to cells × order | **Recommended owner**: separates profile, local response and quadrature authorities, supports falsifiable refinement and missing-measure accounting |
| E: adaptive interpolated response | Requires branch-labelled interpolant with declared one-sided limits and all knots | Error-driven nodes, more bookkeeping and variable cost; deterministic only with fixed implementation/tie policy | Independent comparison; error estimator controls its declared interpolant, not interpolation bias or hidden onsets |

D is selected for its explicit domain and interpolation boundaries, not because
it wins a single numerical comparison. The experiment retains all five routes.
E deliberately integrates a *linear nodal response*; its difference from D is a
representation spread, not an estimate of D's quadrature truncation. B's
nonuniform experiment is not advertised as classical composite Simpson.

### Deterministic production policy proposed now

1. Freeze the exact radial profile knots, interpolation identities and domain;
   never replace the governing EOS/Geometry interpolation through this layer.
   The current test declares linear interpolation of stored `r,m,nu,nB` and
   evaluates the analytic model between nodes. Canonical TOV construction still
   uses its existing EOS interpolation. These are distinct authorities.
2. Obtain ordered event **brackets** from an authenticated branch-partition
   adapter. The generic contract is an adapter returning the complete event
   inventory, left/right status, active axes, certified bracket, and physical
   support/refusal metadata. A generic integrator must not discover onsets by
   response steepness. An adapter can implement authenticated branch-transition
   bracketing; lack of a complete inventory is a refusal.
3. Track-R's test adapter uses provider-declared neutron/muon onset densities,
   its six-way status semantics, independent phase-space checks, and root/bracket
   mapping into each declared background representation. Vacuum is the terminal
   boundary. Onset uncertainty is retained; a stored floating onset is not an
   exact real number.
4. Split at profile knots and each event bracket boundary. An isolated exact
   event has zero measure: request values/status only, **never H** on
   `MuonThresholdEvaluation`, `NeutronThresholdEvaluation`, or vacuum. Interior
   nodes must authenticate the expected branch. A finite unavailable interval
   is not an isolated zero-measure event.
5. Default segment rule: GL16. Compare GL8/16/32 on the identical sealed partition;
   add bisection of all segments and GL32/64 when needed. Sum in radial order
   with compensated or pairwise accumulation, explicitly recorded. On the
   segment adjacent to a continuous square-root onset, map
   `r=r_onset±length*t²`, `0<t<1`; include `2*length*t` Jacobian. This is a
   coordinate map, not smoothing or a density floor. General exponents or jumps
   require their adapter's source-backed map/one-sided policy.
6. Estimate truncation from independent order and partition refinements;
   require contraction or a demonstrated roundoff/provider-noise plateau.
   Do not infer Richardson order from nonmonotone differences. The test-side
   provisional enclosure uses twice the sum of adjacent differences plus the
   independent-method spread and accumulated roundoff. It is a numerical
   estimate, not a theorem. Persistent disagreement triggers refinement/refusal.
7. Acceptance uses requested **componentwise absolute accuracy goals**, declared
   before inspecting the final result, and the total propagated budget. The
   coefficient layer has no universal default relative tolerance. For exact
   structural zeros, only absolute budgets apply. A user-requested goal smaller
   than the available provider/background authority causes refusal.
8. Integrate the center with a documented regular-center expansion/enclosure;
   finite pressure cut is not vacuum. Require a separate positive tail bound,
   lapse-normalization correction, and background uncertainty. No skipped cell,
   arbitrary endpoint substitution, unreported threshold H, or density floor.

## 3. Durable GC9 curved-GR oracle

Expected mathematics is isolated in
`tests/analysis/chemical_curved_reference.py`. The numerical adapter and explicit
wrong routes are in `tests/analysis/chemical_curved_gc9.py`. The expected module
must never import a production integration helper or the adapter under test.
When production exists, redirect the numerical adapter to its public integration
path while preserving this pre-existing expected mathematics.

Classification: **INDEPENDENT ANALYTIC ORACLE**. Manufactured exact star:
`R=12 km`, `M=9/5 km`, `m(r)=M(r/R)^3`,
`exp(nu)=[3 sqrt(1-2M/R)-sqrt(1-2Mr²/R³)]/2`.
Both the lapse and proper-volume factor vary. Susceptibilities are scalar `1`
and `1+(r/R)²`; multiplying either by a declared constant symmetric matrix gives
the corresponding matrix reference. No production Geometry integral supplies
any expected value.

Let `k=2M/R³`, `a=3 sqrt(1-2M/R)`, `theta=asin(sqrt(k)R)`.
For constant C the independent primitive is

```text
8*pi/k^(3/2) * [a*theta + sin(theta)
 + (1-a²)*2/sqrt(a²-1)*atan(sqrt((a+1)/(a-1))*tan(theta/2))].
```

It follows by `r=sin(theta)/sqrt(k)` and polynomial division of
`sin²(theta)/(a-cos(theta))`. Independent 70- and 100-digit theta quadrature
checks the primitive and defines the variable-shape reference. Double radial
GL16/32/64 is compared against it. The predeclared roundoff envelope is
`gamma_4096`, with binary64 `u=eps/2`, plus the observed order difference;
it is specific to this smooth analytic fixture, not a production tolerance.

The controls separately evaluate M6 omit lapse; M7a extra inverse lapse;
M7b interpret supplied `nu` as `2Phi` and therefore use `exp(-nu/2)`; M8 omit
curvature; M19 use `exp(+nu)`; M20 invert curvature. The alternative convention
error `nu→2nu` is numerically identical to M7a and is not counted twice.
Every separation must exceed twice the independent oracle numerical envelope;
no observed percentage is used to set a threshold. Final measurements appear
in the execution section.

## 4. Refusal-window measure

`chemical_trackr_fixture.cpp` brackets live provider availability on the upper
neutron side and **both** muon sides, querying typed values at exact onsets.
`chemical_trackr_budget.py` independently derives model response bounds at
70-digit precision. The source guard is `n_n >= 2^30*downward_ulp(nB)`
(`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:465`); the accepted-root
residual and the 64-epsilon endpoint guards are distinct contributions.

The neutron upper enclosure starts from the source ULP guard, adds the root
density uncertainty `(5e-11 + subtraction margin)/(Dp+De)_min` and two ULPs,
then independently reconstructs the equilibrium charged density from the
neutron chemical potential. `(Dp+De)_min` is evaluated at an explicitly wider
density endpoint where monotonicity gives the lower derivative. The test verifies
the source-derived upper endpoint encloses the measured provider boundary.
Muon enclosures invert the signed endpoint residual around
`mu_n=mu_p(n_e(m_mu))+m_mu`, using a common conservative bound covering both
source guards. These are source/model bounds under the declared floating-point
error model, not merely a list of sampled refusals.

All equilibrium species densities and intrinsic `chi_i=mu_i*pFi/(pi² hc³)`
increase with the common potential. Thus each missing diagonal of neutral C is
bounded by its maximum intrinsic susceptibility; `|C_e,mu| <= min(chi_e,chi_mu)`.
This bounds the **entire unavailable response**. Free-gas `ne` and `nmu` entries
remain structural zeros. The nn-only effect is reported separately and must
not be promoted to a bound on the whole unavailable H bundle.

Two radial mappings are recorded: the sealed profile's linear nB map, and the
independent cold differential TOV relation

```text
|dr/dh| = r(r-2m)/(m+4*pi*r³P) <= r(r-2m_min)/m_min,
h = log(mu_B / reference_mass),
|Delta G_ij| <= 10^54 * max(weight) * max(|C_ij|) * Delta r.
```

Whole containing-cell extrema bound the geometric factors. The larger radial
enclosure is used. This prevents a narrow *interpolated* onset interval from
being mislabeled as a physical-model bound. Propagate the resulting matrix E
through global baryon reduction and Z, including nn-only and all-component
effects. The fixture acceptance allocation is the same predeclared numerical
separation budget as the global test, not the review's quoted final numbers.

**Generic future rule:** a finite refused interval intersecting physical support
is accepted only with a complete authenticated width bracket, a model-backed
susceptibility majorant, computed G/Q/Z error and proof it fits the already
declared global accuracy goal. Otherwise **REFUSE**. The test-side model adapter
supplies mathematics inside its declared free-gas domain; this does not authorize
a generic integrator to interpolate over missing provider output.

## 5. Finite-cut chemical tail and background separation

The fresh canonical fixture uses `structure1::generate`, central mass-energy
density `1.10e15 g/cm³`, and table/radial refinements `(4096,40000)`,
`(8192,40000)`, `(8192,80000)`. It must report `SURFACE_REACHED`; this is the
canonical positive-pressure cut, not P=0. The fine table's SHA-256 is
`7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a`.
Freshly generated profile/table hashes, model revision, masses and hc are output
with each test run. Review/cached G/Z values are never inputs.

In the pe-only outer support, `C_pe=chi_p*chi_e/(chi_p+chi_e)` is monotone
increasing in density: both chi functions increase and the harmonic combination
has positive partial derivatives. Only the y-e axis has tail support. With
`h_cut=log[(mu_p+mu_e)/(m_p+m_e)]`, positive pressure/density, outward decreasing
enthalpy, and `r>2m`, the comparison construction is

```text
R_u = 2 M_c / [1-(1-2M_c/R_c)*exp(2h_cut)]
V_u = 4*pi/3 * (R_u³-R_c³)
delta M_u = epsilon_cut_geom * V_u; M_u = M_c+delta M_u
E_tail,ee = 10^54 * V_u * exp(-nu_cut)*C_pe,cut/sqrt(1-2M_u/R_c).
```

These are analytic bounds conditional on their stated monotone cold-pe and
background endpoint assumptions. The small difference of cubes is evaluated
as `(Ru-Rc)(Ru²+Ru*Rc+Rc²)`. They do not bound interior TOV interpolation error.

The direct test integrates independent enthalpy TOV from the **same cut** to
exact h=0 with DOP853, evolving the tiny shell mass as a separate increment.
It uses phase-space momentum quadrature for epsilon/P and an analytic inversion
of `sqrt(mp²+p²)+sqrt(me²+p²)`. A second 60-digit exterior-shell integral uses
the independent constant-M radial solution and a separate susceptibility
expression; its neglected self-gravity has its own mass-based allowance.
The direct result, its refinement error and this independent check must lie
inside the analytic tail bound.

The lapse-normalization allowance is separately bounded by
`delta M_u/[R_c(1-2M_u/R_c)]`. The already normalized finite-cut lapse is not
silently treated as the exact-vacuum lapse. For full background characterization,
`chemical_enthalpy_reference.hpp` extends the pre-existing test-only
`structure1::EnthalpyOracle` with an independent integral, using its phase-space
EOS and RKF45, no production TOV/profile/Geometry integration. Two ODE/center
controls provide its numerical spread. The difference between this whole star
and the finite-cut profile includes background/interpolation differences and
cannot be assigned wholly to the tail. The shell-mass bound is **not** a bound
on the whole background.

## 6. Executable free-gas correction-separation gate

`chemical_trackr_budget.py` integrates on the same freshly authenticated fixture:

* **Old:** four unprojected intrinsic diagonal integrals `B_n,B_p,B_e,B_mu`;
  `Z_npe=1/Bn+1/Bp+1/Be`, `Z_np=1/Bn+1/Bp`,
  `Z_npmu=1/Bn+1/Bp+1/Bmu`. No projector is called on this route.
* **Corrected:** independently solve the equilibrium common-potential problem,
  build the analytic charge-neutral response, integrate in y, and solve
  `Z=L^T solve(G,L)`. Separately compare actual provider active-H solves.

The old error propagation follows its **four intrinsic axes**, including their
own tail/refusal bounds. It does not reuse a corrected response helper or simply
copy corrected error estimates. The gate requires each of the three named
absolute separations and a matrix-norm separation to exceed the sum of both
propagated uncertainty estimates. Relative percentages are output characterization
only. The Frobenius report explicitly states its denominator (old matrix).

This proves that the correction machinery is active and consequential in the
test mathematics. It does not prove production correctness or reproduce
published A18 coefficients. Under accepted ADR-0013 §7.1 it is an intermediate
GC13 rung; **GC13 realistic closure remains blocked**.

## 7. GC1–GC14 acceptance matrix

Statuses describe this planning increment. `EXACT` means exact arithmetic for
the declared fixture; `ANALYTIC-BOUNDED` means a mathematical bound with explicit
assumptions; `NUMERICALLY-BOUNDED` means a declared numerical estimate with
refinement/independent checks; `SOURCE-LIMITED` and `UNRESOLVED` are not passes.
All floating production gates remain to be exercised against their future API.

| Gate / target | Oracle / classification | Primary error | Secondary error | Conditioning | Source uncertainty | Acceptance / threshold status | Remaining blocker |
|---|---|---|---|---|---|---|---|
| GC1 units/sign/index | Exact units and signed eta fixture; CONTRACT plus signed falsifier | Arithmetic | Unit conversion | None in rational test | Rest-mass/convention identity | EXACT now; gamma/absolute unit closure predeclared | Future typed API wiring |
| GC2 neutral reconstruction | Rational species reconstruction; CONTRACT | Arithmetic | Axis mapping | None | Declared chart | EXACT equality; wrong reconstruction differs | Future API wiring |
| GC3 coupled-energy H | Independently substituted polynomial and finite differences; ANALYTIC ORACLE | Differentiation/step | Local evaluation roundoff | H active modes | Declared toy, no realism | EXACT rational central differences now; future numerical h ladder | Production local adapter comparison |
| GC4 inverse action | Exact rational inverse and residual; CONTRACT, not independent physics alone | Solve backward error | E_H | Declared scaled H | Inherited provider | EXACT fixture; future rho_H enclosure and forward goal | Actual factorization diagnostics |
| GC5 x/y equivalence | Independently differentiated y polynomial; ANALYTIC ORACLE plus identity | Coordinate map | Solve | Congruence scaling | Chart fixed | EXACT independent response; componentwise floating envelope later | Future API wiring |
| GC6 intrinsic vs neutral | Coupled K and independently supplied projected rational matrix; ANALYTIC ORACLE | Local solve | Projection arithmetic | Positive active K/H | Toy declared | EXACT rational equality; no shared projector between expected and neutral route | Production adapter comparison |
| GC7 singular/null/no extra projection | Rank/null plus test-only refusal/type harness; CONTRACT/HARNESS | Illegal inverse/projection | Numerical null residual | Full singularity is required | None | EXACT null/singular rejection; future compile-type and call-path tests | Production API types; idempotent repeat cannot be a numeric mutant |
| GC8 active/onset | Exact embedded npe/pe, phase-space exponents; ANALYTIC ORACLE; live status brackets | Onset partition | Refusal measure | H grows; C vanishes in appearing mode | Authenticated event/support map | ANALYTIC-BOUNDED assumptions; absolute embedding and missing-measure budgets | Future partition adapter and node-level checks |
| GC9a GR measure | Independent exact curved star and variable C; ANALYTIC ORACLE | Rule truncation | Arithmetic/reference | Smooth positive fixture | Manufactured exact | ANALYTIC/NUMERICALLY-BOUNDED; gamma plus independent refinement | Redirect actual production path |
| GC9b convergence | Five rules, fixed partition/order/background ladders, tail/refusal | Quadrature | Provider/profile/onset/tail | Propagate to G/Q/Z | Model-specific adapters | NUMERICALLY-BOUNDED method predeclared; requested absolute goals before run | Realized production error amounts must be measured |
| GC10 global baryon reduction | Exact two-zone independent closed form; ANALYTIC ORACLE | Schur arithmetic | E_G | a denominator, lambda_min(Q) | Toy declared | EXACT noncommutation now; future enclosed Q and global-only API | Production Schur path |
| GC11 named R2006 Z | Independent supplied symmetric M and source expansions; SOURCE TRACEABILITY / CONTRACT | Index/algebra | Global solve | Q/G enclosure | Formula authority only | EXACT fixture; cannot count as independent physics | Production named-axis/accessor checks |
| GC12 W/source action | Declared signed I/Z and c conversion; CONTRACT plus signed/unit falsifiers | Z/I uncertainty | Multiply/units | Z action, component cancellation | Phase-5B qualifications retained | EXACT signed action; future E_W vs requested goal | Real structural-owner connection and provenance |
| GC13 correction/source | Independent old/free-gas route; separate published A18 benchmark | Free-gas budget / extraction | Configuration match | G/Q/Z plus source envelope | **A18 unavailable** | Free gas NUMERICALLY-BOUNDED precursor; realistic SOURCE-LIMITED | A18 lineage and arbitrary-composition/phase data, author arrays or governed extraction |
| GC14 provenance/domain | Expected refusal inventory; PROVENANCE/MUTATION HARNESS | Stale/lifetime logic | Missing numerical metadata | Diagnostics must be retained | Source-byte identity | EXACT expected refusal semantics predeclared; runtime UNRESOLVED until API | Lifetime owner bundle and all-dependency mutations |

GC1–8,10–12 precursors are in `chemical_exact_oracles.py`; the exact two-zone
counterexample has `Q_global=[[14,-1],[-1,19]]/5`, whereas
`sum Q_local=diag(8/3,7/2)`. The independently supplied inverse is
`Z_global=[[19,1],[1,14]]/53`. Expected local/global sides do not share a
reduction helper. GC12 uses Z=`[[2,1],[1,3]]`, I=`(-3,-5) count s²`,
Omega=`2 s^-1`, Omega-dot=`-1 s^-2`; W=`(-11,-18) MeV s²` and the source
action is `(+44,+72) MeV/s`. Flipped spin derivative, omitted/doubled c^-2,
omitted e/mu, and swapped channels are explicit wrong routes.

### Preserved V1–V12 crosswalk

V1→GC1; V2→GC2; V3→inherited Phase-5A equilibrium plus GC8;
V4→inherited lepton validation, not new credit; V5→GC3–6;
V6→GC3; V7→GC3 and inherited local finite-perturbation validation;
V8→GC5; V9→GC6–7; V10→GC8 and GC10 conditioning;
V11→GC9–11 and GC13. **V12 is not discharged by chemical coefficients** and
remains the later end-to-end published thermal/evolution benchmark.
This preserves ADR-0013:296–313 rather than replacing ADR-0010's ladder.

## 8. Local and global numerical semantics

Local H must be symmetric, in the declared active chart, with dimension 1/2/3,
strict physical support and positive/stable active modes. Absent rows are never
padded into an H solve. Compute `H X=I` or requested response actions using a
stable symmetric factorization; explicit inverses are unnecessary. Library
selection remains replaceable. Report symmetry residual, minimum/maximum
eigenvalues, backward residual, scaling matrix, condition number in that scaling,
and a forward enclosure from both provider and solve errors.

Large kappa near continuous onset is expected. The appearing susceptibility
vanishes while an H mode diverges; a universal kappa cutoff would reject the
correct limit. Use actual absolute response error and support uncertainty.
If a certified response cannot be obtained, use the explicit interval-bound
contract or refuse. Never regularize the physics to make a matrix well conditioned.

G is a supported symmetric positive response accumulated in y. Require an
authenticated positive-measure support/rank proof from the union of local active
subspaces, symmetry residual within absolute arithmetic error, and
`lambda_min(G)>||E_G||_2` on that supported space. An absent global muon channel
is explicitly reduced; an all-pe star supports no baryon-conserving beta channel.
Do not infer rank by determinant or invent an absent response.

Q is formed **after** the complete global G integral and x transformation.
Require `a>E_a` and positive supported Q eigenvalues exceeding E_Q; support is
either the two named beta channels or an explicitly declared reduced channel set.
For Z prefer `L^T solve(G,L)` or a symmetric supported Q solve; neither permits
an inverse/pseudoinverse of the full corrected 4×4 response. Verify both routes
against independent oracles during implementation, not by shared expected code.

Require `rho=||G^-1 DeltaG||_2<1`, certified through an uncertainty majorant,
or an equivalently justified perturbation criterion. A small residual alone does
not establish forward accuracy. There is **NO universal kappa threshold**.

### Structural-zero register

Track-R free gas authenticates G_ne and G_nmu as identically zero: neutron
intrinsic decoupling and q_n=0 prevent charged projection from coupling them.
Each declaration records axes, model/revision, derivation/source, domain, and
absolute rounding/provider enclosure. Tests inspect raw computed entries against
absolute budgets. Do not zero entries or symmetrize a production result after
calculation. For a different interacting model, these declarations are absent
unless independently authenticated. The test's analytic reference has exact
zeros; provider-H numerical residuals are measured independently.

## 9. Provider/background matching

A valid match includes model/EOS identity and revision; exact source-data bytes;
particle masses, hc and unit conventions; lepton inclusion; composition and
equilibrium chemical conditions; total energy including rest energy; pressure
where the provider exposes it; active phase/domain and interface semantics;
and the requested/achieved central state. Equal names or a generic "free gas"
label are insufficient. The existing Track-R fixtures authenticate all constants,
the owning table builder and table hashes.

At each anchor compare model epsilon, pressure and composition to the actual
profile. Report energy/P relative error where nonzero, species absolute density
error and density-scaled error (never relative error against an absent species).
The test records a separate response effect from profile-composition versus
model-equilibrium C. It also compares provider active-H response to the independent
model. Finite interpolation mismatch is budgeted, not renamed exact agreement.
No universal acceptance tolerance is inferred from this one table. Future generic
matching requires source-backed interpolation/error certificates or independently
refined characterization sufficient for the requested coefficient goal. Identity,
phase, conventions, or unsupported-domain mismatch causes refusal; it cannot be
hidden by a larger numerical tolerance.

## 10. Exact proposed provenance fields — API planning only

The future immutable result metadata should use these explicit fields (names
are proposed, not implemented):

```text
authority.adr_id, authority.ratification_sha, authority.source_hashes
model.id, model.revision, model.constants, model.lepton_convention
background.profile_identity, background.version, background.content_sha256
background.eos_identity, background.eos_content_sha256, background.interpolation_id
background.central_requested, background.central_achieved, background.domain
axes.local_chart, axes.active_species, axes.global_y, axes.output_beta, axes.input_lepton
quadrature.rule_id, quadrature.implementation_revision, quadrature.orders_attempted
quadrature.accepted_order, quadrature.mapping_id, quadrature.partition_sha256
quadrature.partition_radii_km, quadrature.segment_branch_ids, quadrature.node_count
quadrature.refinement_history, quadrature.accumulation_rule, quadrature.comparison_route
onsets[].id, onsets[].provider_status_left, onsets[].provider_status_right
onsets[].density_bracket_fm3, onsets[].radial_bracket_km, onsets[].mapping_authority
onsets[].source_interval_error, onsets[].value_only, onsets[].query_count
refusals[].reason, refusals[].interval, refusals[].width_authority
refusals[].susceptibility_majorant, refusals[].E_G, refusals[].E_Q, refusals[].E_Z
tail.policy_id, tail.cut_pressure, tail.cut_radius_km, tail.h_cut
tail.R_upper_km, tail.M_upper_km, tail.support_axes, tail.susceptibility_bound
tail.lapse_normalization_bound, tail.bound_authority, tail.direct_validation_id
structural_zeros[].axes, structural_zeros[].derivation_id, structural_zeros[].absolute_budget
local_diagnostics.scaling, local_diagnostics.symmetry_max, local_diagnostics.eigenvalue_ranges
local_diagnostics.backward_residual_max, local_diagnostics.condition_range
local_diagnostics.forward_error_max, local_diagnostics.locations_of_extrema
global_diagnostics.G_support_proof, global_diagnostics.Q_support_proof
global_diagnostics.scaling, global_diagnostics.eigenvalues, global_diagnostics.symmetry
global_diagnostics.rho_G, global_diagnostics.rho_Q, global_diagnostics.solve_residuals
errors.norm_id, errors.componentwise_axes, errors.requested_absolute_goals
errors.components, errors.E_G, errors.E_Q, errors.E_Z, errors.E_W
errors.classifications, errors.assumptions, errors.unmeasured_components
lifetime.owner_bundle_id, lifetime.covered_dependency_ids, currency.validation_handle_id
structural.complete_provenance, structural.sequence_source_ids, structural.I_phys_units
```

The result stores realized details, not merely a policy name or input tolerances.
All attempted refinements, refusals and changed precision belong in provenance.

## 11. Lifetime-safety plan

Audit: `NumberInput` and `NumberProvenance` in
`CompactStar/Analysis/ParticleNumberResponse.hpp` retain raw star/profile and
Hartle response pointers while the EOS has shared ownership. `RequireCurrent()`
in `CompactStar/Analysis/src/ParticleNumberResponse.cpp` dereferences those live
objects; it cannot safely validate an already dangling pointer.
`EquilibriumSequenceNumberDerivative` retains contributing stars, but copying its
provenance into `FixedBaryonNumberResponse` does not automatically transfer that
sequence owner's lifetime. Retaining only a shared pointer to the final structural
response is therefore insufficient.

Minimum future policy: an explicit owner-held dependency bundle with shared
ownership of the central NStar, provider, EOS sources, complete structural result,
and **every contributing sequence star/owner**. Validate that each raw address
reachable from retained Phase-5B provenance is covered by the bundle **before
storing or dereferencing it**. Star ownership must cover its Hartle response
objects; replacement is detected by existing currency checks. Where external
ownership cannot be transferred, a validated safe lifetime handle is required
and failure to obtain it causes construction refusal. An identity hash alone
does not make raw-pointer validation safe.
The ownership bundle must be established while the structural sources are alive
and retained from their construction/capture; matching a newly reused address
cannot retroactively authenticate an expired dependency. The chemical factory
must refuse a legacy result whose original lifetime coverage cannot be proved.

Retain immutable identity snapshots plus live safe validation handles. Every
value/error accessor revalidates source currency. Shared ownership prevents
destruction, not mutation; version/content/response-identity checks remain
mandatory. W retains the complete Z and structural dependency closure; copied
I numbers alone are forbidden. GC14 must mutate central and sequence profiles,
same-version different profiles, EOS bytes, response identities, provider revision,
tail/onset policy, owner coverage, and stale/expired handles. Expected outcomes
are explicit refusals, including before any dangling dereference. No Phase-5B
redesign or raw-pointer edit is made here.

## 12. Proposed files and APIs — no implementation

Repository convention: public headers under `CompactStar/<subsystem>/`, sources
under its `src/`, namespace `CompactStar` for local EOS interfaces and
`CompactStar::Analysis` for number responses. Follow these existing conventions.
Names below are plans, not new production declarations.

| Future owner / placement | Inputs / compute contract | Immutable outputs / units / domain | Provenance, currency, error, refusal |
|---|---|---|---|
| `CompactStar/EOS/ChargeNeutralNumberSusceptibility.hpp` and `EOS/src/ChargeNeutralNumberSusceptibility.cpp` | Owned provider, authenticated active evaluation/chart, requested solve goal; `Compute(...)` returns result or typed refusal | Active response and explicitly embedded C_y; fm^-3/MeV; active dimension 1/2/3; named y axes; no threshold H | Provider identity/domain; H diagnostics and E_C; refuse invalid symmetry/support, unavailable finite interval without adapter, uncertain solve |
| `CompactStar/Analysis/GlobalChemicalNumberResponse.hpp` and `Analysis/src/GlobalChemicalNumberResponse.cpp` | Owned NStar/provider bundle, matched EOS identity, partition adapter, tail certificate, absolute error goal; `Compute(...)` | G_y in count/MeV, supported axes/rank and optional full singular view only; immutable metadata/E_G | Complete realized quadrature/onset/tail/background provenance; `RequireCurrent()` on values/errors; refuse unmatched source, missing measure, unsupported rank, unmet goal |
| `CompactStar/Analysis/ChemicalImbalanceResponse.hpp` and matching `Analysis/src/` | Owned current global G result, explicit supported beta channel set, solve goal; `Compute(...)` | One immutable Z matrix in MeV/count; output `BetaChannel::{Npe,NpMu}`, input distinct `Lepton::{Electron,Muon}`; eta=-Z deltaN | Global Schur/solve provenance and E_Q/E_Z; read-only `PaperZnpe()`, `PaperZnp()`, `PaperZnpMu()` views, never separate storage; refuse uncertain rank/inversion or unsupported requested channel |
| `CompactStar/Analysis/RotochemicalSpinDrive.hpp` and matching `Analysis/src/` | Owned current Z, complete `FixedBaryonNumberResponse` and lifetime bundle covering all sequence dependencies; request physical I via existing unit owner | W=Z I_phys, MeV s², named beta channels; no evolved eta or spin ODE | Seal both lineages and E_I/E_W; validate both on access; refuse incomplete owners, units, stale source, unsupported channels or unmet error goal |

Proposed typed refusal reasons: `UnmatchedModel`, `SourceBytesChanged`,
`UnsupportedDomain`, `IncompletePartitionAuthority`, `ValueOnlyThreshold`,
`UnboundedRefusalInterval`, `UnboundedTail`, `UncertainActiveSupport`,
`UnstableLocalMode`, `UnreliableLocalSolve`, `UnreliableGlobalSolve`,
`AccuracyGoalUnmet`, `LifetimeCoverageMissing`, `StaleDependency`.
No production structs or methods are added by this task.

Proposed callable surface, written here as API notation only:

```text
ChargeNeutralNumberSusceptibility::Compute(
    OwnedLocalProvider provider, ActiveLocalThermodynamicEvaluation evaluation,
    ProviderErrorCertificate provider_error, LocalAccuracyGoal absolute_goal)
  -> Result<ChargeNeutralNumberSusceptibility, ChemicalRefusal>

GlobalChemicalNumberResponse::Compute(
    ChemicalDependencyOwners owners, BackgroundProviderMatch match,
    ChemicalPartitionAuthority partition, ChemicalTailCertificate tail,
    GlobalAccuracyGoal absolute_goal)
  -> Result<GlobalChemicalNumberResponse, ChemicalRefusal>
  G(NeutralNumberAxis row, NeutralNumberAxis column) -> count/MeV
  ErrorG(NeutralNumberAxis row, NeutralNumberAxis column) -> count/MeV

ChemicalImbalanceResponse::Compute(
    shared_ptr<const GlobalChemicalNumberResponse> global,
    SupportedBetaChannels channels, ChemicalAccuracyGoal absolute_goal)
  -> Result<ChemicalImbalanceResponse, ChemicalRefusal>
  Response(BetaChannel output, Lepton input) -> MeV/count
  ErrorResponse(BetaChannel output, Lepton input) -> MeV/count
  PaperZnpe(), PaperZnp(), PaperZnpMu() -> read-only views

RotochemicalSpinDrive::Compute(
    shared_ptr<const ChemicalImbalanceResponse> chemical,
    shared_ptr<const FixedBaryonNumberResponse> structural,
    ChemicalDependencyOwners complete_owners, SpinDriveAccuracyGoal absolute_goal)
  -> Result<RotochemicalSpinDrive, ChemicalRefusal>
  Drive(BetaChannel channel), ErrorDrive(BetaChannel channel) -> MeV*s²

All results: Metadata() -> immutable provenance; RequireCurrent() -> success/refusal.
All value/error accessors call RequireCurrent(); construction is through Compute.
```

The input certificate names denote the exact responsibilities in §§9–11 and 14,
not pre-existing repository types. `Result` can use the repository's eventual
chosen error transport; typed refusal semantics are mandatory. Private immutable
storage and factory-only publication prevent partially computed valid objects.
The local result owns only active C plus an explicitly derived embedding; the
global result owns G, Z owns its one named matrix, and W owns its one named vector.
No duplicate independent Q/Btilde authority or evolved state is introduced.

## 13. Future algorithm sequence

| Step | Required semantics | Replaceable algorithm |
|---|---|---|
| A H→C | Authenticated active chart, symmetric stable modes, bounded response solve | Stable factorization/library and equilibration choice, recorded |
| B embed | Exact declared 1/2/3-dimensional embedding into y; no padded H | Static typed maps or equivalent verified algebra |
| C quadrature | Complete onset/refusal/partition authority and error control | Recommended mapped GL16 with independent refinement; documented equivalent rule only after validation |
| D accumulate G | e^-nu, proper volume, 10^54, center/tail accounting; global before reduction | Compensated/pairwise sums, deterministic implementation |
| E baryon Schur | Transform complete G to x, global a/h/D reduction on supported set | Explicit small symmetric Schur or constrained factorization |
| F Z | eta=-Z deltaN; reliable supported solve, named axes | Q solve or equivalent L^T G solve |
| G W | W=Z I_phys, one existing c^-2 unit owner, no channel omission | Stable matrix-vector multiplication |
| H seal | Immutable values/errors and full realized provenance/lifetimes | Value object/result factory details |
| I currency | Safe live validation of every dependency before access | Version/content handles compatible with existing owners |

## 14. Complete error ledger and propagation

Use spectral 2-norm for support/perturbation tests, componentwise nonnegative
absolute matrices for propagation, Frobenius norm only where explicitly reported
for correction separation. Spectral/Frobenius denominators and scaling must be
named; they are not interchangeable. Structural-zero and cancellation-sensitive
entries always retain componentwise absolute budgets. Correlations may only
reduce a bound when explicitly proved; otherwise add conservative components.

| Component | Meaning / propagation entry | Status |
|---|---|---|
| E_H_provider | Provider model/differentiation/interpolation uncertainty on active H | Method PREDECLARED NOW; actual certificate MUST BE MEASURED BEFORE IMPLEMENTATION ACCEPTANCE |
| E_local_solve | Backward residual, factorization/arithmetic, scaled forward response error | Method now; actual per-node values during implementation |
| E_equilibrium_anchor | Matched-source composition/epsilon/P interpolation effects | Test characterization now; production matching certificate/goal measured before acceptance |
| E_profile_background | TOV/profile/EOS representation plus independent background spread | Independent and refinement characterization now; production realized amount measured |
| E_quadrature | Order/partition/refinement and independent-method discrepancy | Rule/estimator predeclared now; realized amount measured |
| E_onset_partition | Bracket width × one-sided response majorants plus mapping uncertainty | Method now; each event certificate measured |
| E_refusal_window | Missing full susceptibility measure, source/model interval bound | Track-R method/test now; adapter-specific amount measured |
| E_surface_tail | Chemical tail, center remainder and separate lapse normalization | Track-R analytic bound/test now; each result endpoint certificate measured |
| E_roundoff | Unit conversion/evaluation/accumulation operation envelopes | gamma method now; actual operation counts/platform recorded |
| E_global_solve | Global transform/Schur/factorization backward error | Perturbation method now; actual residuals measured |
| E_structural_I | Existing physical I uncertainty plus all structural qualifications | Inherited authority; actual matching I error MUST be supplied for W |
| E_source_benchmark | Published arrays/extraction, EOS/configuration/phase uncertainty | SOURCE-LIMITED for A18; cannot assign invented numbers |

Local: with `H+DeltaH`, use a certified `rho_H<1` on the active space and
backward-solve error. A useful norm bound is
`||Delta C|| <= ||H^-1|| rho_H/(1-rho_H)`; embedding amplifies by at most
`||T_active||²`. Keep componentwise enclosures where a norm hides small modes.
If direct source/phase-space bounds are stronger near onset, report their
independent authority instead of forcing the normwise inverse estimate.

G: for integrand `w C`, bound

```text
E_G <= integral (|w| E_C + E_w |C| + E_w E_C)
     + E_quadrature + E_onset_partition + E_refusal_window
     + E_surface_tail + E_roundoff,
```

where E_C includes provider/local solve/anchor contributions and E_w includes
the matched background/lapse/geometry error. Do not double-count a component
without labeling conservative overlap. The current experiment deliberately
retains some overlap between full-background spread and explicit tail bounds;
it does not cancel them by subtracting uncertain numbers.

Q: propagate `E_x=|U| E_G |U|^T + E_transform` to blocks `(Ea,Eh,ED)`.
Provided `a>Ea`, an entrywise enclosure is

```text
E_Q <= ED + (|h| Eh^T + Eh |h|^T + Eh Eh^T)/(a-Ea)
          + |h h^T| Ea/[a(a-Ea)] + E_Schur_arithmetic.
```

Z: include global backward error as a matrix perturbation. With
`B=|G^-1|E_G`, require `||B||_2<1`; then

```text
E_inverse <= (I-B)^-1 B |G^-1|,
E_Z <= |L|^T E_inverse |L| + E_Z_arithmetic.
```

This nonnegative Neumann majorant preserves a tiny nn-only effect; a blanket
`kappa*relative_E_G` can obscure it. An independently validated Q-based
majorant is equivalent and must exceed the corresponding propagated uncertainty
before rank/precision claims. Negative entries in a numerically evaluated error
majorant indicate numerical failure, not negative uncertainty.

W: `E_W <= |Z| E_I + E_Z |I| + E_Z E_I + E_multiply` with I in count s².
For prescribed Omega and Omega-dot the signed action has error
`2|Omega Omega_dot| E_W`, plus spin-input errors if uncertain. No evolution is
implemented. Source benchmark error belongs to the comparison envelope; it is
not permission to enlarge numerical tolerances until a benchmark agrees.

Budget *logic* is predeclared now. Real production error magnitudes and requested
absolute accuracy goals must be fixed/measured before implementation acceptance.
Unmeasured input error is not zero; inability to certify a needed term is refusal.

## 15. GC13 source blocker and extraction plan

GC13 realistic status: **DESIGNED / BLOCKED ON AUTHORITY**. Required missing
authority remains authenticated A18+delta-v+UIX* equilibrium/composition lineage,
arbitrary-composition response, phase/crust/core/interface semantics, matching
configuration and author arrays or governed Figure-1 extraction. No A18 data was
acquired and no figure was digitized here. This blocks later realistic closure,
not free-gas numerical planning.

Carry forward: authenticate PDF/hash/version; close axis calibration against
labels/ticks; declare curve identity for all solid/dashed named pairs; exclude
unresolved crossings; demonstrate rasterization independence; use an independent
second extractor; retain a full uncertainty envelope; round-trip all identities
against source text; predeclare low-mass comparison points near 1.0–1.2 Msun
where curves are identifiable. Figure 2's 2.13/2.14 ambiguity remains an unresolved
source issue and is not silently adjudicated by this task.

## 16. Materialized versus future tests and execution

Three new registered CTests: `chemical_exact_oracles`, `chemical_curved_gc9`,
`chemical_trackr_budget`. They cover exact precursor mathematics, the independent
curved oracle and six wrong routes, and fresh Track-R quadrature/refusal/tail/
correction/background experiments. The C++ producer and independent enthalpy
reference are test-only. They do not instantiate future production classes.
Python dependencies are NumPy/SciPy/mpmath; absence fails visibly rather than
skipping these gates. Fresh results are JSON evidence outside the source tree.

Still future: production API redirection for all applicable oracles; actual
factorization/partition refinements and per-node provenance; full GC14 type,
currency and lifetime mutation suite; physical structural-owner W wiring;
authenticated A18 GC13; later V12 evolution benchmark. Exact identities are not
relabelled independent physics. No existing baseline is generated or superseded.

### Final numerical measurements

Final fresh Track-R evidence:
`/private/tmp/compactstar-phase5c1p.WwogFj/build/tests/phase5c-evidence/run-ds_kph5e/result.json`,
SHA-256 `722ac4162bad21799ed9fc813698665d3c9bd86cc26b75c4dc646a2c019b8e88`.
The fine exported profile (including provider H and actual profile anchor columns)
has SHA-256 `e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890`.
These are evidence hashes, not governed baselines. The test creates fresh outputs
and performs independent checks on every execution.

| Curved GC9 quantity | Constant C | Variable C=1+(r/R)² |
|---|---:|---:|
| Independent integral, before dimensional susceptibility/count multiplier | 9985.8765807897892591843 | 16030.732681068735338006 |
| Correct-route relative error | 1.822e-16 | 3.404e-16 |
| M6 omit lapse, signed relative separation | -0.1968245400 | -0.1933247435 |
| M7a extra inverse lapse | +0.2459600844 | +0.2404804151 |
| M7b nu as 2Phi | -0.1038801987 | -0.1019230964 |
| M8 omit volume | -0.0963708403 | -0.1033704333 |
| M19 sign-flipped lapse | -0.3544461369 | -0.3488482373 |
| M20 inverted volume | -0.1816021674 | -0.1943482469 |

Both M19 and M20 are durable, explicitly executed controls. The fixed arithmetic
envelope is `4.548e-13` relative, independently of the observed separations.

| Quadrature diagnostic | Measured value |
|---|---:|
| Mapped GL8→16 relative spectral G difference | 3.168e-12 |
| Mapped GL16→32 relative spectral G difference | 9.054e-14 |
| Bisected-partition GL16 vs original GL32, relative G | 3.355e-14 |
| Raw trapezoid vs mapped GL32, relative spectral Z | 1.891e-6 |
| Nonuniform Simpson vs mapped GL32, relative spectral Z | 3.764e-6 |
| Unsplit segment GL16 vs mapped GL32, relative spectral Z | 7.666e-8 |
| Mapped GL16 vs mapped GL32, relative spectral Z | 2.890e-15 |
| Adaptive linear-response representation vs mapped GL32, relative spectral Z | 1.959e-6 |
| Fine mapped GL32 realized node count | 648320 |
| Independent whole-star vs mapped profile Z | 1.770e-8 |

The last row is background/model representation characterization, not quadrature
error on the same interpolant. The separate table/radial G spreads are
`2.371e-9` and `1.242e-6`. No common convergence order is inferred from these
different refinements. In particular, accurate quadrature does not erase radial
profile interpolation error near an onset.

| Refusal result | Neutron upper side | Muon lower side | Muon upper side |
|---|---:|---:|---:|
| Source-enclosed density width fm^-3 | 1.26118e-15 | 4.94382e-13 | 4.92717e-13 |
| Differential-TOV radial-width bound km | 2.89734e-9 | 7.08600e-12 | 7.06213e-12 |
| Missing G_nn bound count/MeV | 4.95798e41 | 7.75641e42 | 7.73029e42 |
| Relative G_nn bound | 1.70693e-14 | 2.67037e-13 | 2.66137e-13 |
| nn-only relative spectral Z bound | 1.14453e-17 | 1.79054e-16 | 1.78451e-16 |
| Full missing-response relative spectral Z bound | 2.80578e-14 | 3.59248e-14 | 3.58042e-14 |

The neutron bracket measured live is last refused
`7.356730162240422e-9`, first available `7.356730162240423e-9 fm^-3`;
source-derived upper enclosure `7.356730164911233e-9`. Its downward-ULP neutron
guard is `8.881784197001252e-16 fm^-3`. The neutron missing charged diagonal
also has `E_Gee=1.37898e41 count/MeV`; omitting it would understate the Z budget.
The smaller linear-profile neutron radial width is `2.49454e-11 km`, which is
explicitly not used as the physical-model bound. The review's approximate nn-only
`9.8e-18` and this conservative spectral `1.145e-17` use different margins/norm
accounting; neither number is a tolerance. Full E_Q/E_Z matrices are emitted.

| Tail / background quantity | Measured value |
|---|---:|
| h_cut | 1.3073642069997383e-5 |
| Positive cut pressure km^-2 | 4.0608740658e-20 |
| R_upper km | 12.768154903424017 |
| M_upper-M_cut km | 3.1656856447e-14 |
| Direct same-cut tail G_ee count/MeV | 2.4255976742539e45 |
| Direct refinement allowance count/MeV | 1.0788043227e37 |
| Independent 60-digit exterior-shell integral count/MeV | 2.42559767427363e45 |
| Analytic tail bound count/MeV | 3.68181084450387e45 |
| Direct continuation enclosed | YES |
| Separate lapse-normalization relative bound | 2.8978359742e-15 |
| Tail-only relative Z bound | 7.48947e-10 |
| Independent full-background relative mass difference | approximately 1.5e-10 |

**Review characterization clarification:** the review's approximately `2.685e45`
is reproduced as `2.6843522e45` by subtracting the raw-profile trapezoidal G_ee
from the independent exact-vacuum star. It includes quadrature and background
differences. After changing the profile quadrature that subtraction is
`9.1579241e44`; the same-cut direct tail remains `2.4255976743e45`. These distinct
quantities must not be conflated. The accepted positive tail-bound construction
is confirmed; no accepted physical contract is changed by this clarification.

Maximum actual profile/provider anchor mismatches: energy `1.764e-10` relative,
pressure `2.428e-9` relative, largest density-scaled species mismatch `3.268e-10`.
Independent model/provider C spread on the available nodes is `4.830e-14`;
maximum sampled unscaled kappa(H) is `9463.51`, backward absolute residual
`5.190e-16`. These are sampled diagnostics, not universal provider-domain bounds.

| Correction gate | Relative change from old | Absolute separation MeV/count | Combined uncertainty MeV/count |
|---|---:|---:|---:|
| Z_npe | +1.5497499% | 6.9884708e-56 | 1.6974009e-60 |
| Z_np | +82.7933249% | 2.3428105e-55 | 1.3124176e-60 |
| Z_npmu | +3.7486008% | 3.7102536e-54 | 5.9449453e-59 |

Frobenius separation is `3.7602485%` of the old matrix, or `3.6244800%` of the
corrected matrix. Minimum named separation/uncertainty is `41171.6`; the
Frobenius separation/uncertainty is `62613.6`. All pass without fitting percentage
targets. The full entrywise G/Q/Z and old-intrinsic error matrices are retained
in JSON; the certified perturbation majorant is `rho=3.5262283e-6 < 1`.
Refusal/tail bounds are small within this predeclared *fixture separation* budget.
No A18 claim follows.

Final focused CTest inventory: **3/3 PASS**, 98.19 s, rc=0; Track-R test 97.86 s.
Existing affected data-free inventory: **39/39 PASS**, 374.46 s, rc=0. Runs
were sequential and serial. Together they cover the complete 42-test data-free
inventory with no failures or skipped tests. The existing 39 were selected with
`-LE external-data -E '^chemical_'` after the three focused new tests passed;
the new tests were not omitted from validation or redundantly rerun in that phase.
No unrelated external-data campaign was run. Test changes did not affect
external-data registration or implementations.

Exact existing inventory, all PASS:

```text
particle_number_analytic
phase5b_PB1
phase5b_PB6
phase5b_PB7
phase5b_PB9-11
phase5b_PB12
phase5b_PB13
phase5b_contracts
phase5b_structural_response_regression
compactstar_library_smoke
heat_capacity_v1
tov_reference_analytic
eos_derivative_contract
rotochemical_local_thermodynamics
rotochemical_trackr_freegas_local
rotochemical_trackr_npe
rotochemical_trackr_pe
hartle_monopole_contract
hartle_monopole_physics_analytic
hartle_monopole_measure_contract
hartle_monopole_published
cache_contract
cache_thermal_contract
hartle_moment_inertia_analytic
hartle_normalization_contract
hartle_first_order_physics_analytic
evolution_stepper_contract
photon_cooling_conformance
proper_volume_contract
geometry_cache_measure_contract
tov_surface_sweep_hw
tov_surface_audit_hw
tov_surface_derivatives_hw
tov_surface_consumers_hw
tov_surface_contracts_hw
rotochemical_trackr_freegas_barotrope
rotochemical_trackr_freegas_structure
relativistic_unit_boundary
relativistic_unit_background
```

Build/evidence root: `/private/tmp/compactstar-phase5c1p.WwogFj` (fresh Debug
configuration). Focused and existing logs are `focused-final.log` and
`existing-datafree-final.log`; `inventory-datafree.txt` records the 42-test
inventory. The existing governed Phase-5B fresh-producer regression passed.
`git diff --check` passed. The exact allowed permanent inventory contains eight
files: this record, `tests/CMakeLists.txt`, and the six `chemical_*` test files.
Production/baseline/EOS-data/literature diffs: **NONE**.

## 17. Scope and status boundary

ADR-0013 remains **ACCEPTED**. Production G, Z, W, eta evolution, weak rates,
heating/cooling changes and BNV: **NOT IMPLEMENTED / NOT BEGUN by this task**.
GC1–GC14: **PREDECLARED / TEST ORACLES PARTIALLY MATERIALIZED**.
GC13 realistic benchmark: **BLOCKED ON A18 SOURCE AUTHORITY**.
INV-09 remains **VERIFIED / RESOLVED only within its ratified structural scope**;
all nine qualifications in `PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:98`
remain in force, including the lack of an authenticated free-gas core I_Omega
benchmark and a source-qualified M_max claim. INV-11 remains **UNRESOLVED**.

Permanent scope is tests, test registration and this record only. Production,
baseline, EOS/data, literature, Geometry, RelativityUnits, TOV and Phase-5B
implementation changes must remain **NONE**. No merge or production work is
authorized by this record. A successful disposition A is still conditional on
independent review before the separately governed production task.

## 18. Final disposition and remaining obligations

Phase-5C numerical planning: **COMPLETE**. The recommended owner, durable curved
oracle and M19/M20 controls, source/model refusal bounds, same-cut tail
continuation, correction-separation gate, GC acceptance matrix, error ledger,
lifetime policy and exact API plan are all materialized. No unresolved numerical
semantic choice remains that prevents an independently reviewed production task.

**A. PHASE-5C NUMERICAL BUDGET / IMPLEMENTATION PLAN COMPLETE —
PRODUCTION IMPLEMENTATION MAY PROCEED AFTER INDEPENDENT REVIEW**

Open obligations are explicit rather than assigned zero uncertainty:

1. Independent review must assess this candidate numerical policy, including
   conditional floating-point/source bounds and the two characterization
   clarifications (same-cut tail versus cross-route difference; full refusal
   response versus nn-only effect).
2. Actual production node-level errors, requested absolute accuracy goals,
   support diagnostics, complete GC14 lifetime/currency mutations and oracle
   redirection must be measured/executed before production acceptance. Their
   methodology is predeclared here; those production results do not yet exist.
3. GC13/A18 source authority, later V12, INV-11 and the retained Phase-5B
   qualifications remain open in their existing scopes. No new realistic,
   evolutionary or BNV claim is made.

The commit uses exactly `test: plan corrected chemical coefficients`.
`PHASE5C1P_SHA` denotes the commit containing this final record, reported by Git
after commit rather than self-referentially embedded here. The branch is to be
pushed non-force with upstream; canonical master is not merged.

Exactly one recommended next action:

Run an independent Claude Opus 5 XHIGH review of PHASE5C1P_SHA, focusing on
the selected onset-aware quadrature policy, curved-GR GC9 oracle independence,
refusal-window and surface/tail bounds, free-gas old-vs-corrected
correction-separation gate, conditioning/error propagation, lifetime/provenance
plan, and whether the proposed GC1-GC14 budgets are sufficiently predeclared
to authorize production G/Z/W implementation.

That review is not started automatically.

## 19. Post-independent-adjudication addendum — two-track structural uncertainty, 2026-09-07

**Disposition:** the owner accepts the completed independent Phase-5B to Phase-5C uncertainty
adjudication:

> **END-TO-END K/I VALIDATION ENVELOPE IS SCIENTIFICALLY SUFFICIENT FOR PHASE-5C
> GENERIC/FREE-GAS CANDIDATE IMPLEMENTATION — REVISE N3 SEMANTICS BEFORE PRODUCTION OUTPUT.**

This addendum supersedes only the earlier interpretation that every inherited structural-I
contribution must be a complete deterministic certified bound. It changes no accepted Q1-Q8
mathematics, quadrature policy, `G_y/Q/Z/W` ownership, N1-N9 requirement, or Phase-5C-owned
numerical fail-closed rule. It also changes no Phase-5B central value, formula, implementation,
test, tolerance, or governed baseline. Phase-5B's implementation record already states that
`Errors()` does not universally bound EOS/profile interpolation and that independent/refinement
effects are recorded separately (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_IMPLEMENTATION.md:114`).

### 19.1 Terminology and replacement ledger

The terms below are binding and noninterchangeable:

| Term | Ratified meaning |
|---|---|
| `numerical_error` | Propagated uncertainty of the declared discrete representation/computation, including arithmetic, quadrature on that representation, finite-difference/stencil estimate, local solve residual, roundoff, or a proved remainder explicitly included by the computation. |
| `certified_bound` | A mathematically demonstrated enclosure under explicit hypotheses. Use only when such proof exists; do not infer it from an error field. |
| `validation_envelope` | A conservative predeclared empirical envelope from independent numerical/analytic discrepancies and controlled refinement/representation variants. It is not a probability distribution, confidence interval, formal truncation remainder, or mathematically certified continuum bound. |

The former single ledger row `E_structural_I` (`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:557`)
is replaced semantically by two separately stored and propagated quantities:

| Component | Definition / role | Acceptance status |
|---|---|---|
| `E_I_numerical,i` | `K_error_i/c^2`, for `i in {e,mu}`, using the existing governed unit owner exactly once | Inherited Phase-5B `NUMERICAL_ERROR`; required for numerical W propagation |
| `V_I_validation,i` | Componentwise end-to-end K/I `validation_envelope` assembled under sections 19.2-19.5 | Required and frozen before candidate acceptance; not interval certification |

`RotochemicalSpinDrive` directly consumes only `(I_phys,e,I_phys,mu)`, with
`I_phys,i=K_i/c^2` (`docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md:182`).
Neutron and proton remain necessary construction/validation evidence but are not direct entries in
`W=ZI`. A complete componentwise certified interval for `A_i` and `B_i` separately is not a
downstream prerequisite; it could destroy correlated cancellation information in K. This does
not reduce or discard A/B, neutron, or proton evidence.

The sentence “unmeasured input error is not zero; inability to certify a needed term is refusal”
(`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:611`) is scoped as follows:

- every required Phase-5C-owned `numerical_error` term for G, Q, or Z, and every explicitly
  required `certified_bound`, remains measured/bounded and fail-closed; inability to measure a
  required owned numerical term causes refusal;
- every required inherited structural `validation_envelope` ingredient must be measured before
  candidate acceptance, and a missing ingredient also causes refusal then, but no formal
  deterministic interval certification is required for that empirical envelope.

### 19.2 Required K/I evidence and classifications

For consumed species e and mu, the future implementation task records the following durable raw
evidence. Classifications are primary evidentiary roles; dependence and any corroborating role
must also be documented.

| Required evidence | Classification | Binding handling |
|---|---|---|
| Existing Phase-5B stored `K_error` | `NUMERICAL_ERROR` | Convert once to `E_I_numerical`; do not add again inside `V_I_validation` |
| PB11 finite-q, `q -> 0` direct K discrepancy | `VALIDATION_ENVELOPE_INGREDIENT` | Retain raw per-species Richardson/direct comparison and its method dependence |
| PB11 fixed-baryon central shift and `Delta N_B ~ q^2` | `FALSIFIER` | Keep distinct from the Richardson K comparison; retain both with dependence stated |
| PB7 independent-background homogeneous/sensitivity B discrepancy | `VALIDATION_ENVELOPE_INGREDIENT` | Transfer into consumed K with explicit sensitivity/amplification; document correlation with central-shift/B_B effects |
| PB12 K-level EOS/table-resolution variation | `VALIDATION_ENVELOPE_INGREDIENT` | Retain direct K variation at every declared setting |
| M1 direct K-level radial-resolution ladder | `VALIDATION_ENVELOPE_INGREDIENT` | New required measurement under the governed fixture/semantics before candidate acceptance |
| PB6 partition/knot variation | `VALIDATION_ENVELOPE_INGREDIENT` | Transfer A-level evidence to K with an explicit sensitivity factor; treat same-class radial/knot duplicates by a conservative max/envelope |
| M2 PB10 raw per-species direct K discrepancies | `VALIDATION_ENVELOPE_INGREDIENT` | New durable test-side reporting requirement; PASS/FAIL alone is insufficient |
| Separately reviewed independent K comparison | `FALSIFIER` | Corroboration; count as an envelope amount only if its numerical path is demonstrably distinct and its use is predeclared |
| `W=ZI` algebra/sign/unit fixture alone | `CONSISTENCY_ONLY` | Retain as contract checking; it is not candidate structural validation |

M1 records the adjudicated gap: PB6 primarily varies radial resolution at A level and PB12
varies EOS/table resolution at K level, but no direct K-level radial-resolution ladder presently
exists. The subsequent implementation task must measure K directly across a predeclared radial
ladder using the same governed fixture and physical semantics. This is implementation-time
validation, not new Phase-5B physics and not a prerequisite to begin writing candidate G/Z/W
code. It is mandatory before candidate acceptance and is not run in this documentation task.
The existing PB6 role is partition/refinement validation rather than an independent physics
oracle (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:225`).

M2 records that PB10 validates per-species K against the propagated budget but does not durably
record every raw per-species discrepancy needed for the envelope. The implementation task must
record each discrepancy, which may be done test-side without changing Phase-5B production source.
PB10 is not modified in this documentation task. Its existing strongest evidentiary role is the
independent reconstruction of all K values (`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_RATIFICATION.md:103`).

M3 classifies PB6 knot/profile shifts as raw representation sensitivity and
`VALIDATION_ENVELOPE_INGREDIENT` evidence, not violations of a certified bound. Stored `A_error`
does not claim to contain that reconstruction class. Therefore the large neutron/muon ratios to
`A_error` are retained but not interpreted as a broken bound. For p/e, the existing tail term may
dominate stored `A_error`, explaining a ratio below one. No variant is erased or hidden.

### 19.3 Ladder, floor, and refusal policy

1. For a contracting ladder with no demonstrated asymptotic order, use the measured finest-level
   discrepancy/envelope contribution; do not Richardson extrapolate.
2. For a nonmonotone or sign-alternating ladder consistent with a numerical floor, record that
   floor as a validation-envelope ingredient and do not invent an order.
3. Truly noncontracting or unbounded behavior with no stable envelope causes candidate refusal.
4. Preserve every raw variant measurement in durable evidence, including inconvenient variants.
5. No safety factor may be fitted after seeing candidate W.

These rules complement the pre-existing plan's prohibition on inferring Richardson order from
nonmonotone differences (`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:93`).

### 19.4 Double-counting policy

- Phase-5B K-error/tail components already included in `K_error` are not added again as
  independent validation-envelope terms.
- PB11's Richardson/direct K comparison and fixed-baryon central-shift closure constrain distinct
  aspects and may both be retained, with dependence stated.
- PB7 B discrepancy and central-shift/B_B effects may be correlated. Absent proved covariance or
  cancellation, combine conservatively and document possible overlap.
- PB6 radial and knot variants within the same representation class use a conservative maximum or
  envelope rule, not blind summation of duplicates.
- An independent-review discrepancy may falsify or corroborate; the same numerical pathway is not
  counted twice under different labels.

### 19.5 W propagation and dual predeclared acceptance

Retain the existing numerical propagation, with its inherited structural term renamed:

```text
E_W_numerical
  <= |Z| E_I_numerical
   + E_Z |I|
   + E_Z E_I_numerical
   + E_W_arithmetic.
```

`E_W_numerical` is `numerical_error` and contains no validation-envelope quantity. Separately,

```text
V_W_validation <= |Z| V_I_validation.
```

`V_W_validation` is the currently inherited structural `validation_envelope`. Any future
validation envelope on Z requires separate governance before addition.

For every required W component, candidate acceptance requires both

```text
E_W_numerical <= G_W_numerical
V_W_validation <= G_W_validation.
```

Both goals are componentwise, absolute, derived from pre-production authority, predeclared before
the first production G/Z/W output, and immutable afterward. Either failure causes refusal without
relaxation. Use `AccuracyGoalUnmet` for the numerical goal and a separate classification such as
`StructuralValidationEnvelopeUnmet` for the inherited structural gate. Neither numerical goal is
selected by this documentation task.

### 19.6 GC12 revision

GC12 remains a candidate validation gate and now validates both:

1. W numerical accuracy through `E_I_numerical`, `E_Z`, arithmetic, and the predeclared
   `G_W_numerical`; and
2. inherited structural validation stability through `V_I_validation`, including the new direct
   K-level radial ladder and raw PB10 per-species K discrepancies, propagated to
   `V_W_validation` and tested against `G_W_validation`.

The existing exact `W=ZI` algebra/sign/unit fixture remains `CONTRACT / CONSISTENCY_ONLY`, not
candidate validation (`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:291`).
GC12 cannot pass candidate acceptance until both tracks and both goals exist.

### 19.7 Baseline, realistic Track R, and status boundary

A later governed Phase-5C regression artifact may be installed under this two-track model after
candidate acceptance. It certifies reproducible regeneration of accepted deterministic bytes,
not a formal interval-certified continuum solution. Phase-5B's governed installation likewise
distinguishes byte reproducibility from independent scientific review
(`docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md:44`).

The separation method transfers to later realistic Track R; free-gas numerical envelope values
do not. A18 needs its own adapter authority, EOS/table resolution, phase/interface evidence,
radial resolution, correction-sensitive source benchmark, and validation envelope. GC13 remains
`SOURCE-LIMITED / BLOCKED` on A18 authority
(`docs/validation/PHASE5C1_NUMERICAL_BUDGET_AND_IMPLEMENTATION_PLAN.md:615`). Publication-level
claims may report a clearly labelled numerical-error ledger and validation envelope, but may not
call the envelope a certified bound, confidence interval, or formal truncation error without
separate proof.

At this ratification **NO PRODUCTION G_y/Z/W RESULT EXISTS**. The K-level radial ladder is
**PREDECLARED / NOT YET RUN**; PB10 raw per-species envelope values are
**PREDECLARED / NOT YET RECORDED**. INV-09 remains **VERIFIED / RESOLVED** with all nine
qualifications; INV-11 remains **UNRESOLVED**; GC13/A18 remains blocked on source authority; BNV
is not begun. No production code, test, baseline, EOS/data, literature, build file, or numerical
result changes in this addendum.
