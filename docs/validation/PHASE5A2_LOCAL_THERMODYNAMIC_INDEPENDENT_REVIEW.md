# Phase 5A-2R — Independent local thermodynamic review

## 1. Disposition and authentication

**C. PHASE 5A-2 REVIEW FOUND BLOCKING DEFECT —
CORRECTION REQUIRED BEFORE INTEGRATION**

The integration blocker is inadequate load-bearing validation, principally V8's
self-confirming coordinate assertions. No incorrect production lepton formula or
implemented neutral-state/conjugate convention was found. This is not a request
to redesign ADR-0010 or to change production physics. The smallest correction is
test/evidence work under a separate governed task; none is performed here.

| Item | Authenticated value |
|---|---|
| Reviewed SHA | `bdcaa4ba48ef1a8d20da940fde769785d0c71444` |
| Starting branch | `physics/phase5a-local-thermodynamics` |
| Worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotochemical-local-thermo` |
| Master and origin/master | `1e6dcd8c5cdb28c9d6d959e42e9e9745d54aeac8` |
| Initial branch/upstream/live remote | All equal to the reviewed SHA |
| Relationship to master | One ahead, zero behind; master is the parent |
| Initial tree | Clean, including untracked-file check |
| Authority | ADR-0010, ACCEPTED; no owner decisions reopened |
| Review execution | One agent; independent inspection, focused executable, scratch calculations |
| Repository change class | Documentation only; this review record only |

The original implementation record and the Phase 5A-1A adjudication remain
historical evidence, unmodified. This review does not claim owner scientific
ratification. Line references below refer to the reviewed Git object, not a
future corrected implementation.

## 2. Evidence ledger

Abbreviations identify exact repository paths for subsequent `path:line` evidence:

| ID | File and load-bearing locations |
|---|---|
| ADR | `docs/adr/ADR-0010-rotochemical-off-equilibrium-thermodynamic-contract.md:5` (accepted); `:178` (H); `:192` (basis); `:212` (ownership); `:238` (semantics); `:267` (validation) |
| ADJ | `docs/validation/PHASE5A1A_ADR0010_INDEPENDENT_ADJUDICATION.md:195` (electrostatic elimination); `:227` (sufficiency); `:384` (rational fixture) |
| IMPL | `docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_IMPLEMENTATION.md:43` (state); `:82` (toy); `:103` (basis claims); `:130` (coverage claims); `:156` (regression evidence) |
| API | `CompactStar/EOS/LocalThermodynamics.hpp:34` (coordinates); `:50` (state); `:89` (conjugates); `:99` (H); `:111` (evaluation); `:120` (metadata); `:140` (provider); `:165` (lepton result); `:179` (lepton API) |
| SRC | `CompactStar/EOS/src/LocalThermodynamics.cpp:15` (finite checks); `:21` (energy); `:37` (state); `:65` (constants); `:84` (lepton evaluation); `:134` (scalar wrappers) |
| TEST | `tests/eos/rotochemical_local_thermodynamics.cpp:34` (coefficients); `:142` (inverse); `:201` (matrix-form oracle); `:227` (toy); `:315` (lepton oracle); `:333` onward (V1–V10) |
| BUILD | `CompactStar/EOS/CMakeLists.txt:7`, `:29`; `tests/CMakeLists.txt:55` |

`AGENTS.md` and `GOVERNANCE.md` were read. The requested ADR, ADJ, IMPL,
production files, grouped test, and both CMake files were inspected directly.
Relevant invariant context remains `docs/SCIENTIFIC_INVARIANTS.md:712` (INV-09,
INTENDED BUT UNVERIFIED) and `:828` (INV-11, UNRESOLVED).

For the source check, the shared library README and catalog were consulted and
PDF page 3 (journal p. 570, section 4, equations (9)–(19), especially (10)–(14))
was read directly from:

`/Users/keeper/Documents/CompactStar/literature/rotochemical/2006-Reisenegger-Rotochemical-Heating-of-Neutron-Stars-Rigorous-Formalism-with-Electrostatic-Potential-Perturbations.pdf`.

That corrected source authenticates the elimination and reduced-system
interpretation; 2005 was not needed for a new source adjudication in this review.
The accepted benchmark authority remains
`2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`, subordinate to 2006 on
the specified correction. No literature file was changed. The printed 2006
equation (11) caveat remains an **inferred printed/source omission** of the
`e^{-Phi}` factor, not an authenticated published erratum; see ADR:229 and
ADJ:213. No source quotation is silently corrected here.

## 3. Exact public API and ADR conformance

| Public API | Reviewed semantics and disposition |
|---|---|
| `ChargeNeutralCoordinates::{n_B_fm3,n_e_fm3,n_mu_fm3}` | Exactly x=(n_B,n_e,n_mu), densities in fm^-3; not fractions. Mutable input values are validated at construction. |
| `MakeChargeNeutralCompositionState(coords)` | Rejects nonfinite inputs, n_B<=0, negative leptons, nonfinite proton sum, negative/nonfinite reconstructed neutron density. Computes n_p=n_e+n_mu and n_n=n_B-n_p; never clips. |
| `ChargeNeutralCompositionState::{Baryon,Neutron,Proton,Electron,Muon}DensityFm3()` and `Coordinates()` | Value-owned immutable validated densities; constructor private. General physical state may have an absent species; derivative providers must further restrict their active domain. |
| `NeutralConjugates::value_MeV`, `MuNMeV()`, `EtaNpeMeV()`, `EtaNpmuMeV()` | g=(mu_n,-eta_npe,-eta_npmu), in MeV; accessors return g_0,-g_1,-g_2. Correct actual source, but the public eta accessors lack a direct falsifier (F3). |
| `ChargeNeutralChemicalHessian::value_MeV_fm3`, `operator()(row,column)` | 3x3 H_ab=partial g_a/partial x_b, fixed other independent coordinates, MeV fm^3. Bounds-checked accessor. No inverse or full charged object. |
| `LocalThermodynamicEvaluation` | State, energy density MeV fm^-3, g, H. Energy includes the provider's declared rest-energy convention; ADR:247 requires it for validation, not a future coefficient-only runtime. Packaging it here satisfies this increment, not a demand for pressure. |
| `LocalThermodynamicProviderMetadata` | Eight strings: model_id, model_revision, particle_content, coordinate_chart, temperature_scope, rest_mass_convention, lepton_ownership, smooth_domain. No typed units/phase/capability schema or per-result status structure. See section 8. |
| `ILocalThermodynamicProvider::Metadata() const noexcept` | Borrowed const reference; concrete toy owns its metadata. Reference cannot outlive provider. |
| `ILocalThermodynamicProvider::Evaluate(coords) const` | Arbitrary allowed neutral composition, not forced beta equilibrium. Returns value-owned evaluation. Concrete implementation is responsible for domain and output consistency. |
| `ILocalThermodynamicProvider::EquilibriumStateAt(n_B) const` | Appropriate accepted obligation, **not an accidental extra requirement**: ADR:245 explicitly says REQUIRED NOW, and ADR:271 explicitly names this method. Restricting its supported smooth domain is legitimate. Threshold/phase handling remains an explicit future gate, not implicit total coverage. |
| Virtual provider destructor | Correct polymorphic destruction; no ownership factory or star lifetime coupling. |
| `ColdRelativisticFreeLepton::{Electron,Muon}()` | Reusable component factories using authoritative Zaki constants. |
| `Kind()`, `Name()`, `RestMassEnergyMeV()`, `HbarCMeVFm()` | Species/constant identity; not a generic provider's mandatory charged-potential split. |
| `Evaluate(n)`, `ChemicalPotentialMeV(n)`, `EnergyDensityMeVFm3(n)`, `ChemicalPotentialDerivativeMeVFm3(n)` | Cold spin-1/2 free-gas values, including rest mass; optional derivative at zero. Scalar wrappers delegate to complete evaluation. |

The generic provider requires neither individual charged intrinsic potentials
nor pressure. No temperature coordinate, hidden equilibrium restriction,
electrostatic gauge variable, projector, paper-B inverse, stellar integral, or
global response object is present (API:12, API:134; complete SRC inspection).
No production inverse is present at all. The toy's additional active-domain
checks are test-side (TEST:288), not changes to background EOS behavior.

## 4. Independent free-lepton derivation and numerical check

Let M=m_l c^2 in MeV, hbar-c in MeV fm, p=p_F c in MeV, and x=p/M.
Counting the two spin states in a filled momentum sphere independently gives

```text
n = p^3 / [3 pi^2 (hbar-c)^3],
epsilon = 1/[pi^2 (hbar-c)^3] integral_0^p k^2 sqrt(M^2+k^2) dk,
mu = d epsilon/dn = sqrt(M^2+p^2),
dn/dmu = mu p/[pi^2 (hbar-c)^3],
dmu/dn = pi^2 (hbar-c)^3/(mu p) = p^2/(3 mu n).
```

This route derives the reciprocal derivative from n(mu), rather than copying
SRC:113's expression. Dimensions are respectively fm^-3, MeV fm^-3, MeV,
fm^-3 MeV^-1 and MeV fm^3. No cgs or geometric-km conversion belongs here.
Integrating the momentum integral gives

```text
epsilon = [p mu (2p^2+M^2) - M^4 log((p+mu)/M)] /
          [8 pi^2 (hbar-c)^3].
```

The logarithm equals asinh(x); hence the production closed form is correct.
Independently expanding `8 integral_0^x t^2 sqrt(1+t^2) dt` gives

```text
8x^3/3 + 4x^5/5 - x^7/7 + x^9/18 - 5x^11/176 + 7x^13/416 + ... .
```

All five implemented coefficients and their signs are correct (SRC:21).
The first omitted term is relative O(x^10), below ordinary double roundoff at
x<0.01. This series is analytic evaluation, not threshold smoothing.

Limits provide independent normalization checks:

- Nonrelativistic: epsilon/(nM)=1+3x^2/10-3x^4/56+O(x^6),
  mu/M=1+x^2/2+O(x^4), dmu/dn proportional to n^(-1/3).
- Ultrarelativistic: mu approaches p, epsilon approaches 3np/4,
  dmu/dn approaches p/(3n).
- At zero density: epsilon=0, p=0, mu=M. The positive-density derivative
  diverges as n approaches zero; it must not become a finite zero-density H entry.

Constants are taken from `Zaki::Physics::{ELECTRON_M_MEV,MUON_M_MEV,MEV_2_INV_FM}`
(SRC:65; `dependencies/include/Zaki/Physics/Constants.hpp:119`, `:158`, `:163`).
The reciprocal MeV-to-inverse-fm convention is dimensionally correct. Existing
V4's literal checks also guard gross constant drift; this review did not replace
repository constants with a new authority.

### Scratch experiment and results

A separate C++ probe linked the existing `build/libCompactStar.a` and Zaki
archive. It called the actual production lepton evaluation; no production or
test source was copied and altered. A Python check used 80-digit Decimal
arithmetic for the momentum-space antiderivative and reciprocal `dn/dmu`, plus
SciPy quadrature of the independently normalized integral

`epsilon/(nM)=3 integral_0^1 u^2 sqrt(1+x^2 u^2) du`.

For **each** electron and muon, sampled target x values were
`1e-8,1e-5,1e-3,0.009999,0.01,0.010001,0.1,1,10,1e4,1e8`.
Actual rounded input n and the exposed constant values were used in the oracle,
not the rounded target x. True relative errors, not `max(1,epsilon)` scaling,
were measured.

| Quantity | Maximum sampled relative error |
|---|---:|
| p | 9.54e-17 |
| mu | 1.70e-16 |
| epsilon | 6.95e-13 |
| dmu/dn | 3.23e-16 |
| Independent quadrature versus Decimal energy | 4.45e-16 |

The energy maximum occurs just above the series switch, x=0.010001. On this
Apple arm64 build `numeric_limits<long double>::digits=53`; the direct branch
has small cancellation despite the `long double` spelling. Classification:
**BENIGN NUMERICAL**, not a wrong integral, series coefficient, or rest mass.
No production change is required by this measurement. It does deserve a
branch-aware regression bound, not a claim of machine-relative energy accuracy
everywhere.

Both species independently returned correct zero-density values with absent
derivatives and rejected negative, NaN, and infinite input. Source inspection
and the focused V10 run confirm that the zero-density scalar derivative throws.
Nonfinite computed results throw (SRC:120); there is no clipping. This does not
promise relative accuracy for subnormal densities or every representable double.

Scratch files were isolated under `/tmp/compactstar-5a2r.Smkf7D/` (`probe.cpp`,
`check.py`, probe executable), not repository deliverables. The equations,
sampling grid, precision and results above are the durable reproduction recipe.

## 5. Independent toy derivation

The stated test potential (IMPL:82; TEST:34, TEST:240) is

```text
d=(b,e,u)=x-(0.20,0.04,0.02),
epsilon = 200 + 950b + b^2 + 5e^2 + 6u^2 - 2be - 2bu + 5eu + (10/6)b^3.
```

Direct polynomial differentiation, without using either test helper, gives

```text
g = (950+2b-2e-2u+5b^2, -2b+10e+5u, -2b+5e+12u),
H = [[2+10b,-2,-2],[-2,10,5],[-2,5,12]].
```

These agree with the concrete toy and its matrix-form oracle. The coefficients
of H have units MeV fm^3; the cubic coefficient has units MeV fm^6. The
arbitrary linear energy term and offset are TEST physics, not neutron/proton
mass predictions. Nonzero cross derivatives test more than a diagonal model.

Solving the independent two-equation system g_e=g_mu=0 gives
`e=(14/95)b`, `u=(10/95)b`. This reproduces `EquilibriumStateAt` exactly
(TEST:263). Over 0.15<=n_B<=0.25 the leading Sylvester quantities are
`2+10b`, `16+100b`, `142+950b`, whose minima are 1.5, 11 and 94.5.
They are strictly positive on the whole declared interval, not merely the
three sampled points. Positive active densities are additionally checked.

The exact linearization residual is `(5 delta_b^2,0,0)`, independent of the
composition perturbations. Thus V7's direction with delta_b=h must give
`5h^2`, explaining its reported factor-four decrease. The observed minimum
order, 1.9999999944, is consistent with this polynomial and floating arithmetic.

At V1's x=(0.21,0.043,0.021), direct differentiation predicts
`g_e=0.015 MeV`, `g_mu=0.007 MeV`, hence
`eta_npe=-0.015 MeV`, `eta_npmu=-0.007 MeV`.
The actual public accessors have these signs. However, the existing test
constructs its optional charged split *from g* (TEST:279) and never calls those
eta accessors. That identity is not an independent sign oracle; see F3.

## 6. V1–V10 load-bearing assessment

PASS below means the focused executable's observed result. Coverage judgments
are separate; a printed PASS does not establish every claimed falsifier.

| Test / location | Intended defect and actual coverage | Oracle class / limitation |
|---|---|---|
| V1 / TEST:333 | PASS. Checks the optional split's algebra and a metadata units substring. Incorrect public eta accessor signs would not fail: neither accessor is called. A split derived from g cannot independently establish g's physical sign. | Contract/declaration check; materially incomplete sign falsifier (F3). |
| V2 / TEST:351 | PASS. Three states satisfy reconstructed charge/baryon identities exactly. Wrong dependent sums can fail. Preserving each input/getter and rejecting invalid inputs are separate concerns; V10 covers representative invalid inputs. | Contract checks, not exhaustive state-space proof; source construction is straightforward and correct. |
| V3 / TEST:369 | PASS. Three equilibrium anchors, largest composition-g residual 2.78e-17 MeV. A wrong nonzero composition residual fails. Does not cover absent-species equilibrium branches, intentionally outside this toy domain. | Analytic equilibrium/contract check; adequate scoped anchor, not threshold EOS validation. |
| V4 / TEST:385 | PASS. Independent expression form (`pow`, `sqrt`, `log` versus `cbrt`, `hypot`, `asinh`) but shared physical formula and exposed constants. Wrong leading factors at the sampled densities fail. No positive sample reaches the production small-x branch. | Analytic implementation comparison, partially independent, not a separate integral derivation. Missing production branch falsifier (F2). Review's integral/Decimal checks independently confirm current physics. |
| V5 / TEST:427 | PASS. Expanded polynomial provider versus matrix-form energy/g/H at one state; detects missing cross entries or unilateral formula changes. Does not numerically differentiate the energy. Correlated energy/gradient convention changes can escape an oracle sharing the same coefficients. | Analytic oracle with some shared assumptions, not a complete integrability test. Independent derivation in section 5 confirms actual formulas. |
| V6 / TEST:447 | PASS, zero asymmetry. Asymmetric mixed entries fail. A symmetric wrong Hessian or jointly reversed signs can pass symmetry. | Necessary identity only; not sufficient thermodynamic consistency evidence. |
| V7 / TEST:460 | PASS, O(h^2) remainder in one mixed direction. Wrong H acting on that direction generally gives O(h) and fails. Directions invisible to that vector are not tested. It does not hold two coordinates fixed in separate column checks. | Directional linear-response/convergence test, not a full derivative-orientation oracle; see section 6.2. |
| V8 / TEST:490 | PASS. Linear transform is a round trip from H_x through the same supplied U/T; inverse equality is another algebraic identity. Nonlinear term is added and immediately subtracted without an independent transformed-energy Hessian. | Equivalence/algebra test, materially self-confirming for coordinate correctness (F1). |
| V9 / TEST:573 | PASS, largest reported residual 1.39e-17. Diagonal intrinsic chi is supplied independently of provider evaluation and test-side cofactor inversion. Missing projection, wrong charge signs or wrong response magnitude generally fail this fixture. | Genuine independent local analytic equivalence oracle. Not a generic double-projection detector; see section 7. |
| V10 / TEST:609 | PASS, positive minors, finite kappa_inf, invalid and absent-muon cases rejected; an exactly singular fixture throws. No adversarial ill-conditioned nonsingular fixture or scale-invariance test of the inverse guard. | Domain/rank contract plus condition diagnostic. Adequate evidence for this stable toy, not a validated inverse policy; see section 6.3. |

### 6.1 F1: independent coordinate oracle is missing — MATERIAL, integration-blocking

At TEST:496, H_y is calculated as U^T H_x U; it is not obtained from the
energy expressed in y. With the **wrong** T=U=identity, both the round-trip H
assertion and inverse-response assertion still hold. The nonlinear test at
TEST:517 sets `full=naive`, adds g_e/g_mu, then asserts the differences equal
exactly those added terms. Even a wrong J would leave those assertions true.
These are concrete counterexamples to IMPL's claim that V8 falsifies a wrong
basis map. No repository mutation was needed to establish the algebraic holes.

An independent expected y Hessian, derived from the same potential expressed
in y, is

```text
H_y = [[2,0,0],[0,8,3],[0,3,10]] + 10b * ones(3,3).
```

At V8's point x=(0.218,0.049,0.026), it is
`[[2.18,.18,.18],[.18,8.18,3.18],[.18,3.18,10.18]]`.
The actual T gives zero sampled matrix discrepancy. Identity T would miss this
independent H_y by 2.18 MeV fm^3 while passing the current round-trip test.

For z=(B,a,c)=(n_B,Y_e,Y_mu), independently differentiating the expanded
epsilon(B,Ba,Bc) gives

```text
H_zz[0,0] = 2-4a-4c+10a^2+10ac+12c^2+10(B-.20),
H_zz[0,1] = 2B(-2+10a+5c)-.10,
H_zz[0,2] = 2B(-2+5a+12c)-.04,
H_zz[1,1] = 10B^2, H_zz[1,2] = 5B^2, H_zz[2,2] = 12B^2,
with symmetric partners.
```

The independent Hessian agrees with `J^T H_x J + sum_i g_i d2x_i/dz2`
to 6.67e-16 in scratch arithmetic. At the sampled point the required extra
mixed terms are 0.084 and 0.081; the naive transformation is genuinely wrong.
Those numbers come from differentiating the transformed potential, not
subtracting a term just added. Units in this fraction chart are mixed: the
B-B entry has MeV fm^3, B-fraction entries MeV, fraction-fraction entries
MeV fm^-3. This is an independent validation fixture, not a new provider API.

**Required correction:** compare both transformed Hessian and response
matrices against independently differentiated y/z fixtures, and demonstrate
that a wrong coordinate map or omitted/wrong chain-rule term fails. Keep the
accepted x basis and production interface unchanged.

### 6.2 Additional material missing falsifiers

**F2 — V4 misses the small-x production branch (MATERIAL).** The smallest
sampled muon x is about 0.124; electrons are more relativistic. Every positive
V4 sample therefore bypasses `x<0.01` (SRC:25). A wrong series coefficient
could leave V1–V10 green, including the zero-density case which returns before
the series. Required correction: exercise both species below and on both sides
of the branch using an independent integral/series oracle, meaningful relative
energy bounds and limiting/energy-gradient consistency checks. The verified
current series should not be changed merely to fill this testing gap.

**F3 — public imbalance accessor signs are untested (MATERIAL).** Reversing
`EtaNpeMeV()` or `EtaNpmuMeV()` in API:94–95 would not affect this executable:
it reads `value_MeV` instead. Required correction: directly check these public
accessors against fixed independently predicted nonzero eta values, and retain
the optional split consistency check as a separate contract identity.

**Additional V5–V7 limitation:** a symmetric Hessian error `alpha vv^T`, with
v=(0,2,3), has zero action on V7's direction (1,.3,-.2). Scratch evaluation
gives only roundoff residuals (~2.22e-16) for a nonzero composition-block error.
V5/V9 can catch such an error at their particular oracle points; this is not
a claim that it necessarily escapes the entire executable. Separate
held-fixed-coordinate perturbations and an energy-gradient check would remove
the directional/shared-formula weakness without broadening production scope.

### 6.3 Conditioning and fail-closed coverage

V10 actually measures `||H||_inf ||H^-1||_inf`, maximum 23.9259259259.
The three sampled Sylvester minima are confirmed over the full toy interval
by section 5's polynomial proof. Invalid input and the absent-muon derivative
boundary are genuinely exercised. Electron zero/invalid boundaries were also
checked independently in this review.

`Inverse3` (TEST:142) rejects a nonfinite determinant or `abs(det)<1e-14`.
That is not a scale-aware rank/conditioning criterion. For example,
diag(1e16,1,1) has determinant 1e16 and kappa_inf=1e16 but passes this guard.
Conversely scaling a well-conditioned identity can make its determinant cross
the fixed threshold. This helper is test-only and its present accepted inputs
are well-conditioned: **no production inversion defect is asserted**.
Classification: **FUTURE HARDENING**, and **MUST FIX BEFORE TRACK-R** if anyone
proposes to reuse it as a physical inverse/conditioning policy. The <100 fixture
bound is a diagnostic guard, not an authenticated general scientific tolerance.
No production inverse or algorithm should be introduced to address this review.

## 7. Independent V9 projection derivation

This is an independent check of the accepted local algebra, not a new physics
decision. At fixed background redshift let `u=e^{-Phi} delta mu_infinity`.
The 2006 source equation (10) gives

```text
delta n = chi (u-q delta psi),       q=(0,1,-1,-1)^T,
q^T delta n=0.
```

For positive intrinsic chi, `q^T chi q>0`. Solving the constraint yields
`delta psi=(q^T chi u)/(q^T chi q)` and hence

```text
C_CN = chi-chi q (q^T chi q)^(-1) q^T chi.
```

This is the local bracket in source equation (13), without a volume integral.
The potential's effect has been eliminated, not assumed absent. Let

```text
S_x = [[1,-1,-1],[0,1,1],[0,1,0],[0,0,1]],  q^T S_x=0.
```

For the independent full-intrinsic toy K=chi^-1, substitute delta n=S_x delta x
into `K delta n=u-q delta psi` and multiply by S_x^T. This gives
`H_x delta x=S_x^T u`, H_x=S_x^T K S_x. Uniqueness on the positive reduced
subspace proves `C_CN=S_x H_x^-1 S_x^T`, without defining a full paper-B inverse.
Consequently `C_CN q=0` and row(p)=row(e)+row(mu). The full rank is three;
the reduced rank is three. These are not statements about a constructed
stellar/global response in this increment.

For K=diag(2,3,5,7), the independent exact results are

```text
H_x^-1 = [[95/142,7/71,5/71],[7/71,10/71,-3/71],[5/71,-3/71,8/71]],
C_CN = [[1/2,0,0,0],
        [0,12/71,7/71,5/71],
        [0,7/71,10/71,-3/71],
        [0,5/71,-3/71,8/71]].
```

Scratch Python Fraction arithmetic confirmed H_x H_x^-1=I, the two C_CN
routes, the charge null, and the proton-row identity exactly. A separate
NumPy calculation agreed to 1.39e-17 on amplitudes. The checked-in V9's
independent diagonal chi and separate cofactor inverse are not circular:
matching coefficients define the same declared physical fixture, but the
oracle does not call the provider to obtain chi. The fixture matches the
provider at its reference point, where the cubic H correction vanishes;
it does not establish equivalence for a constant intrinsic K at other states.

**Detector-claim limit (DOCUMENTATION ONLY):** V9 cannot universally detect
"double projection." The linear projector
`P=I-chi q q^T/(q^T chi q)` is idempotent, and P C_CN=C_CN (scratch residual
2.78e-17). Repeating that mathematically identical projector leaves amplitudes
unchanged; reusing the subtraction formula with already projected C_CN instead
would have a zero denominator and is invalid. V9 validates one correct local
elimination and its response magnitude. Source inspection—not this single
test—establishes that production contains no second projection. The broad
claim in IMPL:130 and its V9 row should be narrowed in a future correction
record without rewriting this review's historical target.

## 8. API and future-provider hazards

| Topic | Classification | Assessment / bounded consequence |
|---|---|---|
| Equilibrium recovery on base provider | **NO ISSUE** | Explicitly required by accepted ADR:245 and :271. Does not force numerical root finding or require support outside the declared domain. No basis for removing it in this review. |
| Mandatory energy in this evaluation bundle | **NO ISSUE** for this increment | Required by the implementation task and ADR validation. It must be the same cold potential whose gradient/H are returned, including consistent linear rest-energy terms. A later coefficient-only consumer need not create a broader pressure interface. |
| Raw strings and public value aggregates | **FUTURE HARDENING** | Fixed type/comments define x, g and H; strings do not enforce units or consistency and result arrays can be modified by consumers. No wrong current value found. Avoid asserting compile-time dimensional safety. |
| Metadata lifetime and provider ownership | **NO ISSUE** | Virtual destructor; toy-owned stable metadata reference; value-owned returned state. Consumers must not retain the metadata reference after provider destruction. No hidden ownership framework is needed. |
| Provenance for realistic Track R | **MUST FIX BEFORE TRACK-R** provider acceptance | Populate/authenticate model revision, component/particle/rest-mass conventions and domain. Current toy strings are honest but not a realistic EOS provenance record. This can be addressed without inventing a new factory or changing accepted decisions. |
| Canonical DS(CMF) identity | **MUST FIX BEFORE TRACK-P** | Eight unrestricted strings alone do not authenticate model revision/parameter set/phase/mass/lepton conventions or background agreement. A real provider must supply the accepted source-matching evidence; not a claim that a typed metadata framework is required now. |
| Phase/threshold/failure status | **FUTURE HARDENING** now; **MUST FIX BEFORE TRACK-R** crosses such domains | Current toy fails closed with exceptions, and metadata states its positive active domain. Generic physical state can represent zero density, but evaluation has no explicit per-result phase/active-status tag. ADR:245–251 requires that information; document or extend its representation before mixed active domains, not by pretending H remains smooth. |
| Hadronic EOS plus analytic leptons | **NO ISSUE** in interface shape | A concrete coherent provider can own components, reconstruct hadronic densities, and sum epsilon/g/H with consistent conventions. No mandatory charged split prevents composition. Such a composed provider and its no-double-counting tests are not implemented or validated here. |

No **BLOCKING CONTRACT DEFECT** in the implemented public formulas/state
construction was identified. The integration blockers are the material
validation gaps, not hypothetical convenience features. This review makes no
choice about realistic EOS differentiation, interpolation, finite temperature,
or matrix-solve algorithms.

## 9. Scope integrity and regression provenance

`git diff --name-status master reviewed_SHA` shows exactly nine changed paths:
the two new LocalThermodynamics production files, EOS CMake registration,
the new grouped test and test registration, CURRENT_ARCHITECTURE, roadmap,
invariants, and IMPL. The complete production delta outside the two new files
is only their EOS source/header registration plus a final newline (BUILD).
The test CMake delta only adds the new grouped executable/labels. A reference
search finds no production consumer of the new header beyond its own source.

Therefore this commit does not change TOVSolver, NStar, StarProfile,
RotationSolver, existing EOS background behavior, cooling or evolution. It
does not activate RotochemicalCache/ChemState, implement APR or DS(CMF)
off-equilibrium physics, integrate a star, construct paper Z/W, implement
rotational particle-number or baryon-sequence response, or begin BNV.
Historical candidate code elsewhere in the repository is not new implementation.
INV-09/INV-11 remain unresolved in their governed senses; this review changes
neither status.

All eight working-tree baseline SHA256s were independently compared with
`git show 1e6dcd8c5cdb28c9d6d959e42e9e9745d54aeac8:<path>` byte hashes.
Every pair agrees, and the hashes also equal IMPL's eight recorded values:

| Artifact under `tests/baselines/` | SHA256 |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `7b036942f2ae599ace3b2fc8b9a7d91f6d11b5899ab7e8a88d2bb4ee6686493b` |
| `grid_convergence_cmf_1p6_debug.tsv` | `2c68e2f7e871192e00322bcc12b6d8b13eb7fdbda4a6d8d7050f26f0271ce5eb` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `7c84557742ec0ec747756d118b2837fb54095a94c10194d4ccc2d0a778f5b04f` |
| `hartle_I_dscmf1_debug.tsv` | `a21d4c3f6d89322cb4cca7a073e0e142f180e529332d704b7b16a145eae741c9` |
| `hartle_monopole_dscmf1_debug.tsv` | `bd49e5a091ebcc59f7c4899422200181d4e71ecf552284840454d01aac4b8d52` |
| `passive_cooling_cmf_1p6_debug.tsv` | `afcad5d078fe458bdb441278ce56caa2d025becc6b489d00e9baad91a101c0de` |
| `tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

The focused target was built with `cmake --build build --target
rotochemical_local_thermodynamics -j1`, then executed directly as
`./build/tests/rotochemical_local_thermodynamics`. It passes V1–V10, including
V7's measured order and V9's 1.39e-17 residuals. Direct invocation preserved
the existing complete-suite LastTest logs.

No broader suite was rerun. Read-only `ctest --test-dir build-selfcontained -N`
and `ctest --test-dir build -N` authenticate 23 and 46 current registrations.
The committed IMPL:164–165 evidence reports 23/23 in 91.14 s and 46/46 in
686.34 s. Existing logs corroborate 23 and 46 `Test Passed` records, zero
failures, and sequential Sep 04 windows 20:08–20:10 and 20:10–20:21 EDT.
Summed rounded individual test times are 91.15 and 686.34 s; the 0.01 s
difference from self-contained wall time is not a failure. The full cache
names `/Users/keeper/Documents/CompactStar/data/compose`; the self-contained
cache has no external EOS root. These logs corroborate the committed report;
they are local build artifacts, not independently signed run attestations.
This review does not present their old runs as newly executed evidence.

## 10. Ranked findings and correction boundary

| Finding | Severity / integration consequence | Exact smallest required response |
|---|---|---|
| F1: V8 uses self-confirming basis and nonlinear-Hessian assertions | **MATERIAL; blocks claimed V1–V10 verification and integration** | Independent transformed-potential Hessian/response fixtures with wrong-map/omitted-chain falsifiers; no new production transform API. |
| F2: production lepton small-x energy branch has no V4 coverage | **MATERIAL; include in correction before integration** | Independent branch/limit energy checks for both leptons with justified scaling and consistency checks. Current formulas independently verified; no physics repair inferred. |
| F3: mandatory public eta accessors not exercised | **MATERIAL; include in correction before integration** | Direct nonzero independent eta/accessor sign tests, separate from the g-derived optional-potential identity. |
| V5–V7 shared/directional limitations | **FUTURE HARDENING** | Column-wise fixed-coordinate and energy-gradient checks are a small compatible strengthening, not a new algorithm mandate. |
| V10 inverse helper is not a conditioning policy | **FUTURE HARDENING**, before any reuse for physical solves | Keep test-only; do not advertise generic ill-conditioned fail-closed behavior. No new production inverse required here. |
| Series-switch energy cancellation | **BENIGN NUMERICAL** | Record branch-aware accuracy; measured maximum ~7e-13 relative. No current physics defect. |
| Universal double-projection detector claim | **DOCUMENTATION ONLY** | Narrow the claim to what local V9 and separate source inspection actually establish. |

No production/test/data files were changed to address any finding. No original
validation record was rewritten. The review's sole repository addition is this
document; production and test diffs from the reviewed SHA are empty.

The findings do not authorize realistic free-gas stars, APR, DS(CMF), stellar
susceptibilities, Z/W, reaction/evolution physics, or BNV. They do not reopen
the accepted thermodynamic Hessian or equilibrium-state obligations.

**C. PHASE 5A-2 REVIEW FOUND BLOCKING DEFECT —
CORRECTION REQUIRED BEFORE INTEGRATION**

**Exactly one next action:** perform the smallest governed correction before integration.
