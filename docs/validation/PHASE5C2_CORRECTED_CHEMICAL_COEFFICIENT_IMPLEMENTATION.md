# Phase-5C-2 corrected chemical coefficient implementation

**PHASE-5C CORRECTED CHEMICAL COEFFICIENTS IMPLEMENTED AND
CANDIDATE-VALIDATED FOR GOVERNED GENERIC/FREE-GAS SCOPE —
INDEPENDENT REVIEW REQUIRED**

This AI-authored implementation remains an unratified scientific candidate.

## Authority and sequencing

Canonical master: `49ab2b8c2881b6ef7b9309307d18cea51d557f72`.
Reviewed plan: `09d1b3c935919ec85f3c629607797ae87c42d8e5`.
Starting implementation authority: `d3b102d2225c44bea581832b27cdc33eb5874e83`.
The existing clean production branch was authenticated at the reviewed plan and
fast-forwarded with `git merge --ff-only` to its sole direct child, the ratified
uncertainty semantics. No new branch/worktree, rebase, reset or merge commit.

`PHASE5C2_PREDECLARATION_SHA` is
`a87f0212c2bd7bfba92db91dfac82447a6561334`, with message
`docs: predeclare chemical coefficient acceptance`. Its structural measurements,
goals, policy, raw machine-readable ledger and measurement code preceded every
production chemical calculation, including synthetic calculations. The first
production coefficient operation was the synthetic contract executable after
that commit. The predeclaration has not been amended or rewritten.

The controlling semantics distinguish propagated `numerical_error`, an analytic
`certified_bound` under explicit hypotheses, and an empirical
`validation_envelope`. No confidence interval, probability interpretation, or
certified continuum/truncation claim is assigned to the structural envelope.

## Pre-production structural evidence

The immutable [measurement record](PHASE5C2_PREPRODUCTION_STRUCTURAL_ENVELOPE.md)
and [predeclaration](PHASE5C2_PRODUCTION_ACCEPTANCE_PREDECLARATION.md) are controlling.
The complete 20000/40000/80000 radial ladder recomputed A, B, fixed-baryon reduction
and K for all four species. The last gaps contract in magnitude with changing
signs; no radial Richardson order is asserted. The e/mu radial envelope terms
are `6.092770162828933e50` and `4.5385432256275946e50 count km^2`.
PB10 reporting preserves each original acceptance condition and reports all four
raw species discrepancies without changing Phase-5B formulas.

The conservative componentwise V_K values are `4.386661132887441e53` and
`4.927777951659812e52 count km^2`. V_I values are `4.880818755395441e42` and
`5.482892414134075e41 count s^2`. Inherited E_I values are
`2.082871228648889e40` and `9.660006381851622e39 count s^2`.
The ledger excludes K_error/PB13 from V, takes the maximum within overlapping
radial/representation evidence, records PB7/PB11 dependence while conservatively
adding the selected reduction terms, and assigns no second credit to prior-review
values. The finite-q PB11 central constraint is an empirical finite-spin
sensitivity ingredient under the explicit resume instruction, not a continuum
error estimate. This conservative term dominates V.

## Production objects and numerical method

The four objects live in `CompactStar/Analysis/ChemicalResponse.hpp` and
`src/ChemicalResponse.cpp`, compiled and installed by the Analysis target.
They add no downstream state or evolution interface.

* `ChargeNeutralNumberSusceptibility` consumes actual 3D npemu, 2D npe or 1D pe
  neutral H. Vacuum and both exact thresholds are value-only and refuse H use.
  It performs the integer congruence to the active y chart before scaled
  Cholesky, algebraically equivalent to T H_x^-1 T^T. This avoids cancellation
  of exactly decoupled free-gas blocks. Congruence arithmetic, solve residual,
  and input-H errors enter a componentwise positive inverse-perturbation
  majorant. No padded Hessian, four-species projector, pseudoinverse or second
  projection API exists. High condition alone is not a refusal.
* `GlobalChemicalNumberResponse` accepts injected metric and local-C evaluators.
  The actual production kernel integrates `1e54 4 pi r^2 exp(-nu) C_y /
  sqrt(1-2m/r)` with one inverse lapse and one proper-volume factor. Default
  GL16 is checked against GL8, GL32, and bisected GL32 on the sealed partition.
  Adjacent refinements must contract or lie inside the local arithmetic floor.
  The estimator retains twice all adjacent and bisection differences; it is
  numerical characterization, not a certified quadrature remainder. Ordered
  Neumaier accumulation, all ladder matrices and node counts are retained. Source-declared structural zeros are retained explicitly and checked at every local evaluation; a contradictory response refuses.
* `ChemicalImbalanceResponse` derives Q diagnostically by global Schur reduction
  only after G integration. Denominator uncertainty, product cross terms and
  Schur arithmetic are explicit. Z is obtained by supported SPD solves against
  the named beta directions. Q is not stored as authority; paper scalars are
  views into Z. The solve residual and E_G jointly enter the full Neumann
  denominator, including their cross terms. Supported uncertain modes refuse. Curvature-factor arithmetic includes its amplification by the distance from 1-2m/r=0 and refuses an unresolved geometric denominator.
* `RotochemicalSpinDrive` requires a complete, owned, current Phase-5B response.
  It calls `WholeStarIPhysical()` and maps inherited K errors through the
  governed AngularVelocity unit owner once. It retains W, E_W and V_W separately.
  E_W includes |Z|E_I + E_Z|I| + E_Z E_I + arithmetic. V_W is |Z|V_I, rounded
  outward. Both frozen absolute goals must pass. Its only dynamic convenience
  evaluates `2 W Omega Omega_dot` for supplied scalars.

The local binary H is the declared matrix for the local solve. Independent
provider, equilibrium-anchor and background characterization remains separately
charged to E_G; no provider accuracy is inferred from an exact binary solve.

## N1-N9 disposition

| Requirement | Implementation and test |
|---|---|
| N1 | Injected production kernel; GC9 numerical side calls that kernel, independent high-precision reference unchanged. |
| N2 | Both-ended continuous onset segments split at the midpoint before endpoint mapping. Analytic sqrt(x(1-x)) fixture checks pi/8, support, measure and partition. A genuine discontinuity refuses. |
| N3 | Ratified separate numerical_error and validation_envelope; immutable dual W goals. |
| N4 | Explicit global-solve, Schur, Z and W arithmetic; local congruence and solve arithmetic; complete input/solve perturbation denominator and cross terms. |
| N5 | Production verifies positive pe bootstrap, mass-upper, radius-upper, lapse and response inequalities. Known cut mass is used only in the boundary factor; total mass UPPER enters the radius numerator. Direct same-cut continuation and independent high-precision shell are enclosed. Incorrect mass/radius and cut mismatch tests refuse. |
| N6 | Certificate stores containing cells, interval, branch/source availability, whole-cell metric extrema and additive matrix. Production authenticates cell indices/edges and checks all contained metric nodes. Invalid containment, cell index, source authority and geometry refuse. The source-qualified linear-profile extrema theorem, not three sampled points, supplies interval authority. |
| N7 | Source-derived high-precision neutron endpoint uses a fixed 2^30-density-ULP source guard, source root-residual allowance, outward conversion, and fixed representational radius margin. Empirical provider brackets must be contained. No padding is fitted to production coefficients. |
| N8 | Highly conditioned resolved SPD passes; uncertain smallest mode, indefinite H, nonsymmetric H and resolved-G/unresolved-Q fixtures refuse explicitly. |
| N9 | Both finite refusal edges are in the actual partition; no excluded quadrature node is evaluated. Analytic missing-volume fixture and an evaluator that throws inside the interval check support and accounting. |

## Center, onset and tail certificates

The regular center uses cubic mass, constant leading density and nu inside the
first profile node, plus the separately characterized missing-term error.
Finite cuts remain nonvacuum: the omitted pe tail contributes only an error
matrix. The tail theorem uses positive pressure/energy, monotone pe energy,
authenticated source/cut and no-horizon hypotheses:

`R_0 <= R_cut / (1 - h_upper R_cut / M_cut) = R_boot`;

`M_total <= M_cut + (4 pi/3) epsilon_upper (R_boot_upper^3 - R_cut^3)`;

`R_0 <= 2 M_total_upper / [1 - (1 - 2 M_cut/R_cut) exp(2 h_upper)]`.

The first displayed radius expression defines an upper enclosure for the physical
vacuum radius; `R_boot_upper` is its outward-rounded evaluation. The final
inequality follows from monotonicity of
`nu - log(1-2m/r)/2`. The detailed artifact retains the numerical certificate,
source hypotheses and independent continuation results. Refusal intervals are
excluded and bounded using source-derived susceptibility, geometry extrema and
both representation/physical width enclosures. Their contribution is additive
and is never supplied as a quadrature value.

## Lifetime and scientific identity

Global results retain shared ownership of every central/sequence NStar and the
provider, plus actual source bytes and profile identity/version snapshots.
Ownership coverage and ordered central/sequence source identity are checked before dereferencing any Phase-5B raw pointer. Initial spin-drive construction also authenticates canonical species charges/baryon numbers, q convention, units, and common EOS authority.
Z owns its G dependency; W owns Z and the complete structural response. Length-prefixed structural snapshots cover material scalar and metadata fields without concatenation ambiguity. G, Z and
W cannot be reassigned after construction, so replacing a dependency through a
caller-retained mutable shared pointer is impossible. Scientific configuration
is retained by value and authenticated revision tokens; modifying a fresh input
request does not mutate an already-owned result.

Tests cover live validity, caller-handle release, owner removal, actual central
and sequence profile-version changes, a genuinely destroyed foreign source,
profile-owner substitution, swaps between distinct live owned sources, central and
sequence structural-source substitution, structural response fields, provider
revision/bytes, same-path different EOS/certificate/provider-copy bytes, domain, partition,
onset, tail, accuracy, metric, basis, constants and lifetime token. The actual
source files are restored after deliberate test mutations. Persistent provenance
contains hashes and semantic identifiers, never pointer addresses. The artifact also records the actual compiler string, linked GSL version, C++ language level, platform, profile versions and exported particle constants.

## Independent gate classifications

GC1 units/index/sign and GC2 neutral reconstruction remain contract checks.
GC3 finite perturbations and GC6 full-intrinsic corrected rational oracle retain
independent Fraction mathematics and check the actual local production response.
GC4 exercises production solves and refusal fixtures. GC5 checks independent x/y
polynomials against production. GC7 checks the unavailable inverse/projector API
at compile time and singularity in independent rational mathematics. GC8 checks
actual active branches, value-only boundaries and onset/refusal behavior.
GC9 uses the actual production integration kernel against independent 70/100-digit
curved references; all six metric mutants are separated using the unchanged
predeclared analytic criteria. GC9b requires actual GL/partition and table/profile
ladders, independent background, center/refusal/tail and arithmetic accounting.
GC10 compares production global reduction with the exact two-zone oracle; local
Schur reduction differs. GC11 checks R2006 source orientation from independently
supplied synthetic G/M and remains source traceability only. GC12 checks both W
goals and signed/unit/species/channel mutants. GC13 remains SOURCE-LIMITED /
BLOCKED; free gas is not A18. GC14 is executable ownership/staleness validation.

The old F2005 path remains independent intrinsic test mathematics. Named and
matrix separations are compared with combined old/corrected numerical uncertainty;
no historical percentage is a hard-coded acceptance threshold.

## Intermediate failures and their disposition

No raw failure is converted into a pass. An initial target-build invocation before
CMake reconfiguration had rc=2; reconfiguration and subsequent complete builds
returned rc=0. A later support-axis guard edit had a compile-only member-name error (rc=2); it was corrected before final testing. The first synthetic production GC9 call refused its fixed absolute
goal on a single broad integration interval; deterministic eight-cell geometric
refinement resolved the error estimator without changing the analytic goal.
The subsequent independent GC9 reference and mutant checks returned rc=0.

An exploratory reproducibility driver returned rc=1 because source files were
still being edited during that run. All four producer subprocesses returned
rc=0 and numerical values were byte-for-byte equivalent after JSON parsing, but
provenance hashes differed. That run is excluded from reproducibility acceptance.
Final validation uses a sealed source manifest and independently regenerated
artifacts. The predeclaration and all acceptance goals remained unchanged.

## Validation and candidate evidence

Final focused suite: **5/5 PASS, raw rc=0**. Complete data-free suite:
**44/44 PASS, raw rc=0**. Full authenticated repository suite with external
CompOSE data: **67/67 PASS, raw rc=0**. All ran serially on the same sealed
source tree, in that order. No failure in these final suites was absorbed.

Each stage generated a new characterized profile/EOS fixture and recomputed the
complete production response in fresh directories. All three full candidate JSON
files were byte-identical, SHA-256:

`a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b`.

The durable [validation record](phase5c_chemical_coefficient_validation.json)
contains all inventories, raw return codes/console hashes, independent artifact
hashes, actual production table/profile ladder, and protected-file manifests.
Every GC1-GC12 candidate gate and GC14 passed; GC9b passed; GC13 remains explicitly
SOURCE-LIMITED / BLOCKED. GC11 remains source traceability only. Both W goals pass.

All **9 governed baselines**, **14 authenticated external-data files**, and
**28 literature files** retain their entry SHA-256 values. The immutable
predeclaration files and measurement source remain byte-identical to their commit.
Among preexisting tracked files, changes are restricted to Analysis/test build
registration and the exact/curved test redirection. TOV, Geometry, RelativityUnits,
rotation equations and Phase-5B structural mathematics have no diff.

Builds used Apple LLVM 17.0.0 (clang-1700.6.4.2), C++17, Debug/assertions enabled,
GSL 2.7.1, Darwin arm64. The candidate records the toolchain reported by its actual
executable. Existing unrelated C++20-extension warnings were not rewritten.

Reproduction: configure Debug with BUILD_TESTING=ON and a Python environment
providing numpy/scipy/mpmath; run `ctest -L phase5c --output-on-failure -j1`.
For the full suite, configure COMPACTSTAR_EOS_DATA_ROOT to the authenticated
CompOSE estate and run unfiltered `ctest --output-on-failure -j1`. A separate
configuration with an empty COMPACTSTAR_EOS_DATA_ROOT supplies the 44-test
data-free inventory. The production CTest fixture regenerates its own planning
characterization and its own canonical structural/chemical calculations.

No candidate artifact was installed as a governed baseline. No open candidate
validation failure remains within the declared scope. Scientific ratification,
realistic A18 closure and downstream evolution remain outside this result.

## Governed candidate values

Full precision, the sealed partition, actual constants/toolchain and every numerical component are retained in [the candidate JSON](phase5c_chemical_coefficients_candidate.json). It is not a governed baseline.

G (count / MeV):

```json
[
  [
    2.904621518418346e+55,
    0.0,
    0.0
  ],
  [
    0.0,
    2.201380364481809e+53,
    -1.0354115533873709e+51
  ],
  [
    0.0,
    -1.0354115533873709e+51,
    9.746440586307364e+51
  ]
]
```

E_G (count / MeV):

```json
[
  [
    1.0219404103274857e+50,
    1.8349389199423292e+37,
    1.128861506525617e+36
  ],
  [
    1.3796979262051092e+37,
    1.5824359678370865e+46,
    2.9270347720524113e+44
  ],
  [
    5.4936893113781254e+35,
    2.927034771968926e+44,
    2.745206319837477e+45
  ]
]
```

Q (count / MeV; diagnostic only):

```json
[
  [
    2.184981541996094e+53,
    -1.1006095941901028e+51
  ],
  [
    -1.1006095941901028e+51,
    9.743848458410244e+51
  ]
]
```

E_Q (count / MeV):

```json
[
  [
    2.17914539749605e+46,
    5.478831452061045e+44
  ],
  [
    5.478831451947084e+44,
    2.7560651427082615e+45
  ]
]
```

Z (MeV / count):

```json
[
  [
    4.5793031807026964e-54,
    5.172519910278805e-55
  ],
  [
    5.1725199102788054e-55,
    1.0268727975139168e-52
  ]
]
```

E_Z (MeV / count):

```json
[
  [
    4.499198710915918e-61,
    4.285424700647402e-61
  ],
  [
    4.285424702045283e-61,
    2.908175030741565e-59
  ]
]
```

| Quantity | Npe / electron entry | NpMu / muon entry | Units |
|---|---:|---:|---|
| I | -1.1637998545112904e+47 | -1.484481985023382e+46 | count s^2 |
| E_I_numerical | 2.082871228648889e+40 | 9.6600063818516217e+39 | count s^2 |
| V_I_validation | 4.8808187553954412e+42 | 5.4828924141340749e+41 | count s^2 |
| W | -5.4061775017047245e-07 | -1.5845719480103649e-06 | MeV s^2 |
| E_W_numerical | 1.59100967128208e-13 | 1.4843208824397632e-12 | MeV s^2 |
| V_W_validation | 2.2634352552795975e-11 | 5.8826943936771424e-11 | MeV s^2 |

Support is Neutron/Electron/Muon, rank 3. Eigenvalues are `[9.741345082771019e+51, 2.2014313195171725e+53, 2.904621518418346e+55]` count/MeV. The infinity-norm condition estimate is `2995.7010369732106`; the inverse perturbation rho is `3.5183255579973513e-06`. Default GL16 uses `324176` nodes; the full order/partition characterization uses `2431320` nodes.

Frozen absolute goals:

```json
{
  "G": [
    [
      2e+50,
      1e+40,
      1e+40
    ],
    [
      1e+40,
      4e+46,
      7e+44
    ],
    [
      1e+40,
      7e+44,
      6e+45
    ]
  ],
  "Q": [
    [
      5e+46,
      2e+45
    ],
    [
      2e+45,
      8e+45
    ]
  ],
  "W_numerical": [
    1e-12,
    4e-12
  ],
  "W_validation": [
    1e-10,
    3e-10
  ],
  "Z": [
    [
      2e-60,
      2e-60
    ],
    [
      2e-60,
      9e-59
    ]
  ]
}
```

| Corrected vs old entry | Absolute separation (MeV/count) | Combined numerical uncertainty | Margin |
|---|---:|---:|---:|
| Z_npe | 6.9884708244366951e-56 | 1.6552041997265152e-60 | 42221.200414978288 |
| Z_np | 2.3428105032109998e-55 | 1.3019697243209341e-60 | 179943.54703086012 |
| Z_npmu | 3.7102536059678841e-54 | 5.8848648881359181e-59 | 63047.388113325724 |

Matrix Frobenius separation is `3.7256732162286055e-54` MeV/count, with margin `63253.453652792203`. Minimum named margin is `42221.200414978288`. Free gas does not substitute for A18.

The pe tail has R_cut `12.766174760730493` km, M_cut `0.92093273643866647` km, total-mass upper `0.92093273643870355` km, and R_upper `12.768154903424522` km. Its Gee enclosure is `3.681810845442163e+45` count/MeV; the independent same-cut and high-precision shell values are `2.4255976742539016e+45` and `2.4255976742736301e+45` count/MeV, both inside the enclosure.

## Limitations and scope

The whole-star Structure-1 midpoint at rho_c=1.10e15 g/cm^3 is the only primary
physical fixture. The result is generic/free-gas and unratified. INV-09 remains
at its previously ratified/integrated status. INV-11 is unresolved. GC13 realistic
closure is blocked. No ChemState, eta ODE, rates, Urca imbalance functions,
heating/cooling, thermal coupling, spin-down ownership, time-dependent Z, A18 or
BNV is added. Existing TOV, Geometry, RelativityUnits, rotation equations and
Phase-5B structural mathematics are unchanged.

## Phase-5C-2RAT status addendum — 2026-09-07

The historical implementation evidence above is unchanged. Independent Opus review is complete
with disposition B, **PASS WITH NONBLOCKING FINDINGS — candidate ready for human ratification
with explicit caveats**, with 0 blocking and 0 material findings. Human-owner ratification is
complete for the governed generic/free-gas candidate scope at implementation commit
`4d78bf4000848ddecc2127daa2f2840872f266f5`, subject to every explicit caveat in
`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_RATIFICATION.md`.

Canonical integration remains pending. The deterministic candidate artifact is not installed as
a governed baseline by this ratification; no production numerical evidence above is changed.
GC13 remains **SOURCE-LIMITED / BLOCKED**, INV-11 remains **UNRESOLVED**, and no eta evolution,
weak rates, heating/cooling, realistic A18 closure, or BNV is authorized or begun.
