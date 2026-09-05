# Track-R Structure-1 — FR2005 static free-gas structure

**Scientific candidate; independent static validation completed; human ratification is not implied.**
The predeclared `rho_c=1.10e15 g cm^-3` state satisfies the simultaneous FR2005
Table-1 mass, zero-pressure radius and apparent-radius rounding bins. This claim
uses the source-domain and central-state policy below. It does not identify the
authors' unpublished exact input or their retained sequence mesh.

## 1. Authentication, authority and classification

Starting canonical `master`, `origin/master`, and the live remote were all
`f4ae22d971e25bdd74530aec184f3fe0c3440b95`, with a clean canonical checkout.
The fresh worktree is
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-trackr-freegas-structure1`,
branch `physics/trackr-freegas-wholestar-structure1`, created at that exact SHA.
Existing branches touching the local provider contained the ratified ancestors;
no alternate implementation was merged. One agent performed this task.

Classes: EOS architecture/additive capability, numerical method (value inversion,
pressure evaluation and table sampling), test/build plumbing, documentation, and
generated reproducibility artifacts. Authority: GOVERNANCE, accepted ADR-0005
(canonical ordinary TOV), ADR-0009 (surface event and unchanged cutoff), ADR-0010
(local response versus equilibrium background), and
`PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:27-46`. The ratified provider implementation
is `933494d86daf2cf8965079ece49fabd66d9390e5`; its independent review is
`097546f441ec4497cc426a5bb7051d53c2d59da7`. This is an additive equilibrium
capability, not a refactor of unbaselined TOV physics, and invokes no pre-baseline
correctness exception. No scientific baseline is created.

## 2. Direct source authentication

The read-only primary PDF is
`/Users/keeper/Documents/CompactStar/literature/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`,
SHA-256 `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88`.
The PDF itself was read for §3.1, pp. 8–9; Table 1 and all footnotes, p. 30;
radius/metric definitions, pp. 3 and 6; and §4.3, p. 16. The shared catalog gives
it the primary non-superfluid reproduction role. Supporting direct reads were
R1995, PDF p. 15 (printed p. 5, electron-inclusive energy minimization), and
R2006, pp. 1–2 (intrinsic/electrostatic potentials and neutrality). R2006's
correction governs its response scope; it does not introduce a new static
free-gas EOS here.

| Printed free-gas row | Value | Disposition |
|---|---:|---|
| M_max | 0.62 M_sun | Static comparison; footnote c restricts the retained model before Sigma-minus appearance |
| rho_c | 1.10 × 10^15 g cm^-3 | Central mass-energy density, not baryon rest-mass density |
| R | 12.77 km | Coordinate radius at P=0 |
| R_infinity | 13.80 km | R exp(-Phi_surface), at P=0 |
| P_K | 0.98 ms | Excluded: adopted pure-neutron-gas rotation result |

Section 3.1 applies the noninteracting gas to the **whole star**, without a crust
splice. It describes this particular central state as slightly below Sigma
appearance. No rotation, compression figure, thermal figure, or coefficient
comparison enters S10.

## 3. Architecture and production changes

Only four production files receive additive declarations/definitions:

- `CompactStar/EOS/LocalThermodynamics.hpp` and `src/LocalThermodynamics.cpp`:
  `ColdIdealFermionValues` and `ColdRelativisticIdealFermion::Values`.
- `CompactStar/EOS/TrackRFreeGasThermodynamics.hpp` and
  `src/TrackRFreeGasThermodynamics.cpp`: `FreeGasBarotropePoint` and
  `TrackRFreeGasThermodynamicProvider::BarotropeAt`.

The smallest source-specific result fits the existing provider ownership; no
new pure virtual obligation is imposed on unrelated EOS providers. The generator
is a test/run adapter in `tests/eos/structure1/table.hpp:52`, not a new production
TOV implementation. Two self-contained executables own the local and stellar
falsifiers. `local_oracle.hpp` and `tov_oracle.hpp` are test-only independent
calculations. `collect.py` packages diagnostic outputs and source/artifact hashes.

## 4. Value-only API and source coverage

`FreeGasBarotropePoint` contains n_B, total epsilon, pressure, actual species
number densities in order (n,p,e,mu), physical domain identity, model id/revision,
and `values_resolved`. It contains no Hessian. Vacuum has zero densities, energy
and pressure; no vacuum fractions or finite vacuum slope are invented.
Positive-density fractions are formed only when exporting `Y_i=n_i/n_B`.
Numerically unrepresentable values fail explicitly; `values_resolved` never
promises chemical-response availability
(`CompactStar/EOS/TrackRFreeGasThermodynamics.hpp:48-64`, `:92-93`).

The value path uses stable independent-density parameterizations rather than
requesting an active thermodynamic evaluation:

- pe is charge-neutral n_p=n_e=n_B, with absent neutrons and muons.
- npe solves in delta=mu_n-m_n. The electron potential increment is evaluated
  algebraically relative to neutron appearance, using `log1p`/`expm1` for its
  density increment. Newborn n_n is evaluated directly from delta, rather than
  recovered from n_B-n_e.
- npe-mu inverts the monotone total density of a common lepton potential and
  evaluates intrinsic species values directly.
- Exact neutron/muon thresholds retain their value-only boundary identity.

The stored classification thresholds are unchanged. Even the immediate
`nextafter` neighbors return structure values in the tests. At ultra-close
thresholds, individual floating-point values may tie or differ by roundoff;
strict mathematical monotonicity is not a claim of distinct double values at
every adjacent representable density. Generated rows are explicitly required to
have strictly increasing P and epsilon. Neither species floors nor threshold
smoothing are used (`TrackRFreeGasThermodynamics.cpp:325-400`).

## 5. Pressure and stability

For each spin-1/2 ideal species, with k=p_F c, b=hbar c, and m the rest energy,

```
n = k^3/(3 pi^2 b^3)
epsilon = (pi^2 b^3)^-1 integral_0^k q^2 sqrt(m^2+q^2) dq
P = (3 pi^2 b^3)^-1 integral_0^k q^4/sqrt(m^2+q^2) dq
  = m^4/[24 pi^2 b^3] [x sqrt(1+x^2)(2x^2-3)+3 asinh(x)].
```

For x<0.1, production integrates the convergent binomial series for
`1/sqrt(1+x^2)` through order 12. Its leading term is `n k^2/(5m)`. Pressure
therefore avoids a rest-energy subtraction near vacuum. Energy uses the
existing stable energy-bracket owner. The old response-bearing `Evaluate`
implementation is unchanged (`LocalThermodynamics.cpp:153-180`).

S1 independently compares 96-node phase-space quadratures to production values;
checks `P=sum(mu_i n_i)-epsilon` through the equilibrium enthalpy identity where
subtraction is well conditioned; and checks each species' nonrelativistic limit.
At equilibrium, `h=mu_n` for npe/npe-mu and `h=mu_p+mu_e` for pe. Thus
`P=n_B h-epsilon=n_B d epsilon_eq/dn_B-epsilon` on each smooth branch.

## 6. Exact constants and units

No physical constant is fitted or introduced into production as a competing
owner. The species factories use Zaki; the table boundary uses Units.hpp and the
canonical GSL speed of light. Binary-double values used in this build are:

| Constant | Value |
|---|---:|
| m_n [MeV] | 939.56542051999998 |
| m_p [MeV] | 938.27208815999995 |
| m_e [MeV] | 0.51099894999999995 |
| m_mu [MeV] | 105.65837550000001 |
| m_Sigma-minus [MeV] | 1197.4490000000001 |
| MEV_2_INV_FM | 0.0050677307161563949 |
| hbar c [MeV fm] | evaluated as 1 / MEV_2_INV_FM |
| MEV_FM3_TO_ERG_CM3 | 1.602176634e33 |
| c_cgs [cm/s] | 29979245800 |
| G_cgs [cm^3 g^-1 s^-2] | 6.6730000000000003e-8 |
| raw GSL solar mass [g] | 1.9889199999999999e33 |
| Zaki SUN_M_KM [km] | 1.4766250380501249 |
| GSL G M_sun/c^2 [km] | 1.4767161818921164 |

The retained profile conversion owners additionally use `INV_FM4_2_INV_KM2=0.00026122803039748728`, `INV_FM4_2_Dyn_CM2=3.1615267734966909e+35`, `INV_FM4_2_G_CM3=351767294174610.88`. These are existing Zaki values, not new table constants (`CompactStar/Core/src/NStar.cpp:188-201`).

The imported energy column is **epsilon/c^2**, in g cm^-3:

```
rho = epsilon_MeV_fm3 * MEV_FM3_TO_ERG_CM3 / c_cgs^2
    = epsilon_MeV_fm3 * 1.7826619216278975e12
P_cgs = P_MeV_fm3 * MEV_FM3_TO_ERG_CM3.
```

The oracle independently uses the literal conversion defined by MeV and fm,
then `epsilon_geo=epsilon_cgs G/c^4 * 1e10` in km^-2. Its mass conversion uses
the raw GSL solar mass. The profile conversion follows existing Zaki ownership;
these conventions are kept separate (§19).

## 7. Table schema and provenance

The exact tab-separated header is

```
eps(g/cm^3)  pre(dyn/cm^2)  rho(1/fm^3)  10  11  0  1
```

The first three numerical columns are mass-energy density, pressure and n_B;
extra columns are Y_n,Y_p,Y_e,Y_mu (not n_i). Species ids match current consumers
(`CompactStar/Physics/Evolution/src/StarContext.cpp:495-497`). The importer
accepts a single header, so provenance is in a sidecar, not extra comment rows.
Rows use 17 significant digits and deterministic ordered sampling. Existing
files are never overwritten. No canonical EOS file is modified.

Sidecars record generator v1, the starting source commit (the sidecar's
`source_commit` field), ratified provider commit/model, thresholds, resolution,
positive floor, upper density and conversion. The accompanying generated JSON
also hashes the actual candidate source/test files: the starting commit alone
is **not** asserted to contain the new generator. It records collection HEAD
and dirty status honestly. The generator parameters plus exact source hashes
and artifact SHA-256 identify the uncommitted candidate used for validation.
No table is promoted above the analytic local provider as source authority.

## 8. Threshold-aware nested sampling

Default lower n_B is `1e-14 fm^-3`, giving
`P_floor=1.0850284294719746e15 dyn cm^-2`. Upper n_B is `.616 fm^-3`.
The base grid is logarithmic in n_B. Exact thresholds and logarithmic distances
from both sides are inserted. Base and distance increments halve with each
refinement; the neutron distance range also extends closer to appearance.
The threshold grids are nested. Muon distances extend to relative 1e-8;
neutron distances extend from 1e-8 to 1e-11 over the ladder.

| Base intervals | Total rows | max relative epsilon(P) error | max relative n_B(P) error | max absolute fraction error | max relative interpolant slope error |
|---:|---:|---:|---:|---:|---:|
| 1024 | 1255 | 4.1981e-5 | 4.1981e-5 | 6.5629e-6 | 0.12808 |
| 2048 | 2535 | 1.0256e-5 | 1.0256e-5 | 6.2762e-7 | 0.041902 |
| 4096 | 5127 | 2.5336e-6 | 2.5336e-6 | 2.3705e-7 | 0.013329 |
| 8192 | 10375 | 6.2958e-7 | 6.2958e-7 | 2.8784e-8 | 0.004096 |

These are off-grid direct-quadrature comparisons at every geometric cell
midpoint, not values at the spline's own knots. The largest slope discrepancy
is localized near neutron appearance; it is not a uniform high-precision
chemical-response validation. The stellar convergence and independent TOV
comparisons bound its effect on structure. Steffen remains the sole production
interpolation owner. An attempted over-dense threshold grid produced indistinct
floating-point P/epsilon rows and was rejected; the retained nested grid resolves
the threshold without that failure. No cutoff coefficient or interpolator was
changed (`tests/eos/structure1/table.hpp:52-104`).

## 9. Thresholds and P-7/N-1 separation

Stored neutron, muon and Sigma thresholds remain, respectively,
`7.3567289037328318e-9`, `.45698480541241987`, and
`.61735520796653498 fm^-3`. Values are continuous at each species appearance.
At relative neutron offsets 1e-7, 1e-9 and 1e-10, the new value-only path supplies
composition, epsilon and P while the old `EquilibriumAt` still throws the
N-1 response-resolution exception. The generator includes rows in this window.

The wider relative 1e-8 threshold bridges were independently sampled at 201
points each. Neutron-bridge maxima are approximately `3.0e-13` relative energy,
`5.17e-10` relative pressure and `3.66e-10` absolute fraction; muon errors are
at approximately 1e-14 or less. These bounds address interval coverage without
inventing absent-species Hessians. The original N-1 `2^30`-ULP guard and all
response methods remain unchanged (`TrackRFreeGasThermodynamics.cpp:451-482`).

At adjacent muon doubles a pressure difference of order `2e-14 MeV fm^-3` can
have a roundoff sign. This is explicitly below the independently checked bridge
error; these adjacent queries are not inserted as unordered spline knots. The
source threshold's carried N-6 rounding caveat remains. Available structure
values do not ratify a response in that interval.

## 10. Equilibrium slope, regularity and causality

Define `s_i=dn_i/dmu_i=mu_i k_i/(pi^2 b^3)`, with the genuine zero-density
limit s_i=0 at appearance. This is a scalar along the equilibrium sequence,
not a fabricated inverse inactive-species chemical Hessian. Let L=s_e+s_mu.
For npe/npe-mu,

```
C = dn_B/dmu_n = s_n + s_p L/(s_p+L)
d epsilon/dn_B = mu_n
dP/dn_B = n_B/C
d epsilon/dP = mu_n C/n_B.
```

For pe, replace C by `s_p s_e/(s_p+s_e)` and mu_n by h=mu_p+mu_e.
At neutron appearance s_n→0, so the two one-sided slopes agree. At muon
appearance s_mu→0 and the npe/npe-mu slopes agree. The numerical probes show
the approach; neutron convergence is unusually slow (square-root behavior),
so the algebraic limit must not be inferred from a single finite offset.
Vacuum has P proportional to n_B^(5/3), epsilon proportional to n_B, and a
divergent d epsilon/dP; it is not a finite-slope table knot.

The ideal-gas fixed-composition stiffness satisfies
`n_i dmu_i/dn_i=k_i^2/(3mu_i) <= mu_i/3`. Allowing equilibrium composition to
relax decreases the second variation at fixed n_B (the positive composition
Hessian's Schur complement). Consequently `0<dP/d epsilon<=1/3` for this
stable ideal mixture, also through the one-sided threshold limits. S4 checks
this over the off-grid source-domain sweep. The sparse diagnostic maximum is
0.0758988308951; the proof, not that sparse maximum, supplies the global bound.
Independent differentiated quadratures agree with the analytic slope to
`4.74e-6` over the finite-difference probes. This tolerance includes the
finite-step/roundoff limitation; it is not the stellar error budget.

## 11. Current canonical production path

The current code was traced directly. Equilibrium import remains
`TOVSolver::ImportEOS` → `Hidden_ImportEOS_Vis` (TSV columns, no unit inference)
→ separate Steffen epsilon(P), n_B(P), and fraction splines
(`CompactStar/Core/src/TOVSolver.cpp:506-728`). `GetEDens`, `GetRho`, and
`GetRho_i` consume those splines. `GetEDensDeriv` differentiates its own
energy spline and converts the result to dimensionless d epsilon/dP
(`:1056`, `:1096`, `:1553`, `:1567`). No analytic slope is injected into it.

Fixed-central-energy calls use `SingleStarSolveToTOVPoints` exclusively
(`:2516-2607`). It retains the central clamp, Brent pressure inversion,
1-cm start, canonical RK8PD controls, output partition and pressure-coordinate
surface locator. The driver records requested rho, requested analytic epsilon,
achieved first-profile rho, initial pressure, independent analytic rho inferred
from that pressure, cutoff and completion status. A request outside the adapter's
unclamped domain is rejected before solving. No target-mass solver is called.

`NStar::BuildFromTOV` publishes the canonical points; the supported `Solve(Axis)`
sequence path delegates to the same primitive and then uses Append/finalization
(`:1807-1854`). S5 compares every resulting radial/species/metric column
bit-for-bit at the same state. Existing baryon-number and first-order-I
postprocessing side effects are retained without making new scientific claims
about them. The known sequence mirror-scalar zeros are not used as observables.

## 12. Upper support and clamp authentication

| Quantity | Value |
|---|---:|
| Upper table n_B | 0.616 fm^-3 |
| Stored Sigma n_B minus upper row | 0.0013552079665349881 fm^-3 |
| Upper table rho | 1.1211542516959692e15 g cm^-3 |
| 0.999 upper rho | 1.1200330974442732e15 g cm^-3 |
| Clamp ceiling minus 1.105e15 | 1.503309744427325e13 g cm^-3 |

Thus the entire printed interval is reachable without super-Sigma padding.
No generated row or accepted central request is at/above Sigma, and no exact
Sigma-supremum star is claimed. The diagnostic scan stops at `.615 fm^-3`,
inside the unchanged clamp as well as the source domain.

## 13. Independent static TOV oracle

`tests/eos/structure1/tov_oracle.hpp` owns a test-only enthalpy formulation.
It queries direct analytic equilibrium composition and independent 96-node
phase-space quadratures. It calls neither the production spline nor canonical
TOV ODE, and does not wrap the production solve.

Use H=ln(h/(m_p+m_e)), x=H_c-H, u=r^2, v=m/r^3 in geometric km units:

```
A = (1-2uv)/(v+4 pi P)
du/dx = 2A
dv/dx = (4 pi epsilon-3v) A/u.
```

Independent center initialization at delta=x is
`u=3 delta/[2 pi (epsilon_c+3P_c)]`,
`v=4 pi epsilon_c/3 - (4 pi/5)(epsilon_c+P_c)(d epsilon/dP)_c delta`.
The regular-center expansion is refined with delta, separately from production's
1-cm initialization. A GSL RKF45 scheme (distinct from canonical RK8PD) integrates
to the matched finite-pressure H or vacuum H=0.

| Oracle tolerance | center delta/H_c | M [M_sun] | R_0 [km] |
|---:|---:|---:|---:|
| 1e-9 | 1e-7 | 0.623635588620 | 12.768154909480 |
| 1e-10 | 1e-8 | 0.623635569676 | 12.768154901783 |
| 1e-11 | 1e-9 | 0.623635569290 | 12.768154901038 |
| 1e-12 | 1e-10 | 0.623635569137 | 12.768154900998 |

For canonical comparisons, the oracle uses the **achieved initial pressure** to
recover the same analytic central state; it does not quietly compare different
central states. The finest midpoint finite-surface radius difference is below
`6e-10 km`, and its mass difference is below `3e-11 M_sun`. Every printed-bin
sample also passes the independent finite-surface comparison.

## 14. Radial convergence — S6

The governed ladder is 2500, 5000, 10000, 20000 and 40000 radial samples.
Each table resolution is tested at all five partitions at fixed requested rho.
All complete with exactly the canonical surface pressure. Pairwise partition
bounds remain the unwidened ADR-0009 values: relative mass <=1e-9,
relative radius/lapse <=1e-8 and surface-pressure residual <=1e-8.
The finest table's observed radius spread is below 1e-9 km, far inside that
contract. Sampling partitions do not redefine the star.

## 15. Cutoff sweep — S7

The effective ratios below are the actual cutoff/achieved-P_c ratios, rounded
for display. Every row has SURFACE_REACHED; raw diagnostics include floor,
cutoff, lapse, R_infinity,cut and achieved central state.

| Effective ratio | R_cut [km] | P=0 radius bounds [km] |
|---:|---:|---:|
| 1e-11 | 12.700695700 | [12.767768438, 12.768154902] |
| 1e-12 | 12.738870129 | [12.768082072, 12.768154901] |
| 1e-13 | 12.755989254 | [12.768142333, 12.768154902] |
| 1e-14 | 12.763219789 | [12.768152832, 12.768154900] |
| 1e-15 | 12.766174760 | [12.768154569, 12.768154903] |
| floor/P_c=1e-16; cutoff/P_c=1e-15 | 12.766174759 | [12.768154568, 12.768154902] |
| floor/P_c=1e-18; cutoff/P_c=1e-15 | 12.766174757 | [12.768154567, 12.768154901] |

M_cut stays near `.623635569 M_sun`; the missing tail mass is negligible at
the allocated mass precision. Lowering the floor below the canonical coefficient
produces the expected plateau. `R_cut` is never relabeled as the P=0 radius.

## 16. Independent P=0 tail bound

At cutoff (r,m), direct EOS enthalpy H_s bounds the exterior tail without a
production spline. In consistent geometric units,

```
r_upper = 2m / [1-(1-2m/r) exp(2H_s)]
m_upper = m + (4pi/3) epsilon_s (r_upper^3-r^3)
g_upper = [m_upper+4pi r_upper^3 P_s]/[r(r-2m_upper)]
r_lower = r + H_s/g_upper.
```

These follow from `-dH/dr=(m+4pi r^3P)/[r(r-2m)]`, positive mass growth and
monotone pressure/energy in the tail. Both denominators are checked positive.
For the finest midpoint, the bracket width is `3.330e-7 km` (0.333 mm),
with half-width 0.167 mm. The independent full enthalpy star lies inside the
bracket, allowing for the separately bounded canonical integration error.
The source-comparison estimate is the bracket midpoint. The residual tail is
about **1.980 m**, not an adjustable surface correction.

## 17. Predeclared midpoint

| Quantity | Result |
|---|---:|
| requested rho_c | 1.1000000000000000e15 g cm^-3 |
| n_B,c from value-only inversion | 0.6049252413288515 fm^-3 |
| requested epsilon_c | 617.05474641848934 MeV fm^-3 |
| achieved first-profile rho | 1.1000000000000000e15 g cm^-3 |
| independent analytic rho at achieved P_c | 1.099999999073754e15 g cm^-3 |
| achieved P_c | 4.9156528789505761e34 dyn cm^-2 |
| p_cut | 4.9156528789505761e19 dyn cm^-2 |
| M_cut | 0.623635569063 M_sun |
| R_cut | 12.7661747607 km |
| surface profile lapse | 0.925057781353 |
| R_infinity,cut (profile convention) | 13.8004079508 km |
| R_0 tail-bracket midpoint | 12.7681547369 km |
| R_infinity,0 (profile convention) | 13.8023679117 km |
| independent direct-vacuum R_0, requested midpoint | 12.7681549010 km |

A separate **90-digit mpmath** phase-space inversion, using the actual binary
constant values printed by the compiled provider, gives
`n_B,c=0.6049252413288510735586358360115556 fm^-3` and
`P_c=30.6810920899593572322004677385559 MeV fm^-3`.
The corresponding energy is `617.054746418489747041366820889524 MeV fm^-3`.
This independently recomputes the midpoint rather than importing the audit's
number. The supplemental source/check script, compiled-provider dump and JSON
are preserved under the final regression artifact directory's `precision/`.
At the three P-7 probes, fresh high-precision neutron-density relative errors
are `1.40e-9`, `1.67e-7`, and `1.78e-6`; energy and pressure agree to better than
`7e-16` there. This makes the newborn-species comparison more discriminating
than the fraction-normalized bridge check alone. The increased relative error
close to onset reflects the retained threshold's finite precision; it does not
license the chemical Hessian, whose N-1 guard remains intact.

The small difference between requested and achieved analytic rho is measured
rather than hidden by the imported interpolant's exact central value. It is
approximately `9.3e5 g cm^-3` here and well below the density budget.

## 18. Entire printed interval and nearby upper branch — S8

The independently recomputed input interval is approximately
`n_B in [0.602305396323278,0.607544233003639) fm^-3`.
Nested 17-, 33- and 65-state scans cover the complete half-open rho interval;
the upper diagnostic uses `nextafter(1.105e15,1.095e15)`.

| Image of rho interval | Minimum | Maximum |
|---|---:|---:|
| M [M_sun] | 0.623147391953 | 0.624120126175 |
| R_0 [km] | 12.7547527962 | 12.7816286996 |
| R_infinity,0 [km], profile convention | 13.7900057567 | 13.8147955958 |

M increases and R_0 decreases on every sample. Endpoint images are stable under
refinement. At inserted central midpoints, the radius error of interpolation
from the preceding mesh falls from `1.419e-7` to `3.559e-8 km`.
Fourteen final-mesh states between `1.0990625e15` and
`1.10109375e15 g cm^-3` pass all bins with the allocated numerical margins.
This is a **tested discrete common-state subset**, not a claim to have computed
the maximal continuous matching interval. The predeclared midpoint itself
already supplies a common state; no mass inversion was needed.

A separate 17-state scan over n_B=.59 through .615 fm^-3 shows increasing mass;
its finite-sample largest mass is `.625480949269 M_sun` at .615 fm^-3.
This lies above the printed mass bin, while still below Sigma. Therefore the
published row is not the largest member of **our** retained mesh. The authors'
exact mesh/last retained state is unpublished. The comparison establishes the
printed static row at a common source-compatible central state, not a continuous
Sigma-limit mass, an unrestricted mathematical maximum, or the authors' exact
retention algorithm.

## 19. Apparent radius and constant conventions — S9

`NStar::BuildFromTOV` sets the surface lapse by the existing profile metric;
S9 compares its R/lapse to the explicit Schwarzschild identity. At P=0 use
R_0 and the corresponding exterior metric, not the finite-cutoff lapse unchanged.
For the midpoint, the raw-GSL exterior formula gives
`R_infinity,0=13.8024397139 km`, while the profile convention gives
`13.8023679117 km`. The difference is `7.180e-5 km`, or 7.18 cm.
Both remain in the printed bin.

This is an existing **deterministic convention offset**, not numerical scatter
or an unknown error bar. The independent TOV falsifier uses raw GSL mass/metric
conventions consistently; profile publication and the principal R_infinity
comparison use existing Zaki profile ownership. No constant owner is changed,
and cross-convention agreement to 5 cm is explicitly not claimed.

## 20. S1–S12 validation matrix

| ID | Exact falsifier and independent reference | Tolerance/convergence basis | Observed result | Defect falsified |
|---|---|---|---|---|
| S1 | Species pressure versus independent phase-space quadrature; equilibrium enthalpy identity; four species' small-k limits | quadrature/value 1e-10 energy, 1e-8 pressure/identity; reliable-subtraction restriction explicit | PASS; energy 2.10e-11, pressure 2.90e-12, identity 3.92e-12 relative | degeneracy, rest-mass, omitted lepton, factor-three and cancellation defects |
| S2 | Analytic equilibrium compliance versus differentiated independent quadratures; imported Steffen derivative at every cell midpoint | finite-step slope 1e-5; strictly decreasing nested errors; final localized slope envelope .006 | PASS; differentiated slope 4.74e-6; table slope .1281 → .004096 | confusing H with equilibrium slope; untested interpolant derivative |
| S3 | Vacuum; exact/nextafter domains; 201-point threshold bridges; one-sided slope approach; P-7 table rows with H refusal | continuity errors bounded in §9; positive finite values; shrinking one-sided slope discrepancy | PASS; no artificial density interval gap or response-guard weakening | floors, branch substitution, missing onset, fake Hessian |
| S4 | Positivity/order in generated rows; off-grid quadrature causality and charge/baryon/fraction identities | analytic <=1/3 proof; finest identity error <1e-7 | PASS; observed fraction-identity error 1.15e-13 or less | acausal/negative EOS, n_i versus Y_i confusion |
| S5 | Canonical points → bulk NStar versus actual supported Solve(Axis)/finalization | every radial/species/metric column bit-identical; omit unsupported mirror accessors | PASS, both sequence captures identical | competing path physics, finalization/species drift |
| S6 | Five governed partitions at each of four tables | unchanged mass 1e-9, R/lapse 1e-8, pressure 1e-8 relative | PASS; finest R spread <1e-9 km | node-defined radius, incomplete surface, partition-dependent star |
| S7 | Seven floor levels; direct-EOS enthalpy bounds; independent vacuum oracle | effective cutoff plateau; final bracket width <1e-5 km | PASS; width 3.330e-7 km; 1.980-m tail | relabeling R_cut as R_0 or tuning surface coefficient |
| S8 | Nested complete rho interval and separate sub-Sigma branch | midpoint insertions <1e-5 km; endpoint image stable; no inverse mass | PASS; last interpolation error 3.56e-8 km | different central states for different observables, endpoint-as-exact-author-model |
| S9 | Actual profile lapse versus explicit identity; separately reported raw convention | identity 1e-14; fixed convention offset reported | PASS; 7.18-cm offset retained | cutoff/P=0 mixing, hidden constant inconsistency |
| S10 | One common predeclared midpoint with all three uncertainty intervals inside original bins | mass margin 1e-7 M_sun, radius margins 1e-5 km; no bin widening | PASS; midpoint and 14 final-mesh states | plausible-looking isolated observable matches |
| S11 | Invalid upper generator requests and value/central/sequence source checks; requested/achieved/clamp audit | strict stored Sigma ceiling and positive upper-clamp margin | PASS; upper margin .0013552 fm^-3; clamp margin 1.5033e13 g cm^-3 | super-Sigma padding or a silently clamped source state |
| S12 | Independently coded enthalpy RHS, center series, RKF45 and phase-space EOS | oracle tolerance/center ladder; finest finite-star mass 1e-7 M_sun and R 1e-6 km gates | PASS; midpoint R discrepancy <6e-10 km; full rho-bin comparisons pass | shared production spline/ODE error, self-confirming structure tests |

Executable assertions are in `tests/eos/rotochemical_trackr_freegas_barotrope.cpp:8-170`
and `tests/eos/rotochemical_trackr_freegas_structure.cpp:40-249`; the oracle source
and generated JSON provide the detailed independent evidence.

## 21. Aggregate numerical budget

This is an engineering bound conditional on the named constant convention and
source model, not a source error bar. The following conservative radius
allocations sum to `9e-6 km`, below the `1e-5 km` comparison margin and the
requested `5e-5 km` (5-cm) target. They are not twelve allocations each equal
to the entire target.

| Component | Radius allowance [km] | Evidence |
|---|---:|---|
| equilibrium values and pressure evaluation together | 1e-6 | independent quadratures and matched direct EOS |
| EOS sampling/interpolation | 2e-6 | nested stars and direct-EOS oracle |
| central pressure inversion/state recovery | 3e-6 | finest/sequence achieved-rho audit; high-branch maximum below 6.1e8 g cm^-3 |
| radial integration and oracle center/integrator | 1e-6 | partition and independent tolerance ladders |
| central-state sampling | 5e-7 | measured 3.56e-8-km interpolation residual |
| remaining surface tail | 5e-7 | half-width 1.665e-7 km |
| arithmetic and metric propagation within one convention | 1e-6 | direct identity and separately audited fixed conversion |

Mass is conservatively bounded by `1e-7 M_sun` and achieved central density by
`1e9 g cm^-3` on the retained fine tables/scans, below the requested `5e-5 M_sun`
and `5e10 g cm^-3`. Coarse-table diagnostics are convergence evidence, not the
accepted fine-grid uncertainty: the 1024 table's analytic central shift is
`2.565e10 g cm^-3`. At the midpoint dR_infinity/dR is about .99, and the mass
allowance contributes less than 2e-7 km, leaving the propagated R_infinity
allowance below the same 1e-5-km margin. The 7.18-cm convention difference is
reported separately; it is not charged to or hidden in this numerical allowance.
These are supported numerical/convergence envelopes, not interval arithmetic
certification of every floating-point operation.

## 22. Simultaneous publication-bin comparison — S10

The bins are ordinary-rounding precision bins, **not published error bars**:
M in [.615,.625), rho in [1.095e15,1.105e15), R_0 in [12.765,12.775), and
R_infinity,0 in [13.795,13.805), in the units above. The midpoint and its stated
numerical uncertainty intervals lie strictly inside every bin. Rho is a
state-selection constraint, not a prediction after selecting the state.
Both metric conventions pass the apparent-radius bin; none of the bins or
physical constants was widened/tuned. The accepted comparison is therefore:

**TRACK-R STRUCTURE-1 COMPLETE —
FR2005 FREE-GAS STATIC TABLE-1 STRUCTURE REPRODUCED**

This claim remains an AI-authored scientific candidate for human review. It
is narrower than reproducing an unpublished exact source star, the continuous
Sigma-limit mass, or any rotation/chemical-evolution benchmark.

## 23. Complete regression and protected surfaces

The two focused Structure-1 tests pass. A final packaging build/rerun also passes
after an explanatory header-comment correction; no compiled numerical logic
changed after the complete suites. The complete **data-free configuration**
contains 28 tests and passes 28/28 in 94.86 s CTest (94.866 s wall). A label-only
probe selected only six tests; its result is not used as the complete suite.
The initially started full run was interrupted to correct that ordering.
The accepted order is complete data-free suite, then full authenticated suite,
both `-j1`, never concurrent. The full configuration passes **51/51** in **670.28 s** CTest
(**670.286 s** wall), after the complete 28-test data-free run. Both final
inventories, complete stdout and elapsed times are included in the generated
record. All established local and Phase-4 tests pass.

Twelve unique external data files were authenticated against HEAT_CAPACITY_V1
and TOV_REFERENCE, including the official eos.mr. All eight governed baseline
files match their starting SHA-256 values; the generated record lists all eight.
TOVSolver header/source, NStar and StarProfile protected files are byte-identical
to the starting commit. No ninth baseline, EOS data mutation, or literature
mutation occurs. V1–V10, RFG1–RFG11, R1-V1–R1-V10, PE-V1–PE-V13 and all established
Phase-4 tests remain required; no tests are removed or reclassified.

## 24. Reproduction and artifact policy

Configure the normal Debug build with the existing GSL/Zaki toolchain and
Python/NumPy owner. A data-free build omits COMPACTSTAR_EOS_DATA_ROOT; the full
build uses the authenticated external root. Run all CTest tests serially in each
configuration. For a reviewable artifact set, run the two new executables with
fresh explicit output directories, then invoke:

```
python3 tests/eos/structure1/collect.py --structure <star-output-dir> \
  --local <local-output-dir> --regression <completed-regression-dir> \
  --output docs/validation/trackr_structure1_results.json
```

Tables and detailed diagnostics remain in generated build/run directories.
Only the compact generated values/checksum record is committed beside this
narrative. Paths, source hashes, generator parameters, table SHA-256 and exact
suite evidence are retained. Generated tables are reproducibility artifacts,
not ratified scientific baselines. Post-commit branch/remote authentication is
reported in the delivery record; no merge is part of this task.

## 25. Remaining findings and claim boundaries

- The derivative error near neutron appearance converges but is larger than
  smooth bulk interpolation errors. This task validates static structure, not
  a future response/integral consumer's derivative accuracy requirements.
- Adjacent-double threshold values can be indistinguishable or show last-bit
  roundoff. The stored local thresholds and response refusal policies remain
  authoritative; no response-resolution finding is silently closed.
- The raw GSL/profile Zaki convention offset is retained explicitly.
- The mass continues to rise on the sampled sub-Sigma upper branch; the printed
  row does not identify an unrestricted maximum or the authors' exact mesh.
- INV-09 and INV-11 remain unresolved. The analytic matter source is human
  ratified; this new structure implementation/results are not auto-ratified.

## 26. Explicit non-goals

NO STELLAR ROTOCHEMICAL SUSCEPTIBILITY, PAPER B/Z/W,
PARTICLE-NUMBER SPIN-DOWN RESPONSE, ROTOCHEMICAL EVOLUTION,
APR/BPAL, DS(CMF) OFF-EQUILIBRIUM PHYSICS, SUPERFLUIDITY,
OR BNV IS IMPLEMENTED IN TRACK-R STRUCTURE-1.

No canonical TOV equations, controls, clamp, surface ownership, or ADR-0009
coefficient change. No new production integrator, inverse-mass fit, realistic
EOS substitution, source-domain extension, new baseline, merge or human
ratification is included.
