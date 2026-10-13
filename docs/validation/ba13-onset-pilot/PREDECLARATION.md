# BA13 onset pilot — prospective predeclaration

UNVERIFIED SCIENTIFIC CANDIDATE. Validation-only numerical-method, documentation,
and generated-evidence changes. Entry 4f3cf97462ca030f715e26808938f781a2a40945;
branch candidate/linux-ba13-onset-pilot-20261013. GOVERNANCE §§2–5;
INV-02/03/14/15; accepted thermal owner ADR-0002; proposed ADR-0021 Option A.
Interpretation (2), frozen right-hand thresholds and joint deterministic input
errors/nonlinear remainder on the left, remains a proposal for owner adjudication.

## Authentication and blindness

Read-only review manifest verifies all 11 files. The entry-bound STUDY_COMPLETE
manifest verifies 1022 files and source hashes at dcd14d2. All other local branches
are ancestors; the source worktree is clean. See evidence/authentication.json.
Neither the review nor prior evidence may be written. Original 80000/8192
production settings and target B from ba13-recovery/run1/raw/target_stars.tsv
are immutable. Protected production paths, fitter and certificate are hashed.

Discovery: reviewer ν results at targets 0,10,20 on 20k,40k,80k,160k,320k;
all prior production/EOS/chemical evidence. Reanalysis of e2048/e4096 nominal
profiles is discovery, not prospective evidence. Targets 5/15 are secondary
transfer tests on authenticated existing profiles, not fully blind stars.
Actual holdouts: every newly solved star below. Their numerical outcomes must
not be accessed until this document AND numerical predictions.json are committed.
Freeze discovery calculations on Slurm first; commit tools before holdout jobs.

## Mathematical hypotheses and limitations

Cold charge-neutral equilibrium gives ε+P=nB μB and dP=nB dμB. TOV gives
P'=−(ε+P)ν', ν'=(m+4πr³P)/(r(r−2m)), m'=4πr²ε in geometric units.
Thus d(log μB+ν)/dr=0. In pe matter μB=μp+μe, and in npe(μ) μB=μn.
Neutron onset: μp+μe=mn at nn=0; muon onset: μe=mμ. These are constrained
energy-minimum conditions, not an extrapolation of a negative neutron density.
Single charged-species redshifted chemical potentials are NOT assumed constant;
the electrostatic contributions cancel in the neutral baryon combination.
References: FR2005, https://arxiv.org/pdf/astro-ph/0502116 eqs. 1,3;
R2006, https://arxiv.org/pdf/astro-ph/0606322 (neutral electrochemical combinations).
Local Mac literature catalog is absent on this cluster; consulted the primary
papers and authenticated provider/source instead. No realistic EOS claim follows.

Separate Steffen ε(P),n(P) interpolation is not thermodynamically exact. Constancy
must be checked, not imposed as proof of nodal accuracy. Same-profile quadrature
agreement does not independently validate the star, EOS, central B or metric.

## Full case list, execution, resource limits

Targets 0,10,20, identical original absolute B targets:
- A: r640000, EOS8192 (3 stars).
- B: r80000, EOS16384c/32768c/65536c, depth cap11 (9 stars);
  r80000 EOS32768c cap12 (3 stars).
- C: r80001…80005, EOS8192 (15 stars); r80000 reused authenticated nominal.
- D: 5 and 15, retained 20k/40k/80k/160k/320k (zero solves, secondary).
- E: r640000 EOS16384c cap11, target20 (1 star).
Total 31 new production-owner solves. No further resolutions or 21-star expansion.
Table-only checks at 2048/4096/8192 for identity, and all capped members including
32768 depth12; high-precision checks described below. One independent analytic-EOS
structure solver at each of three B targets, two prescribed tolerances 2e-10 and
2e-12, is a cross-check, not an additional production refinement level.

All numerical calculations use Slurm, single-core, at most 2 GiB. Build separate
and reported. Numerical CPU ceiling 3600 seconds including discovery, failed
attempts, high precision and scoring. Sequential controller caps total numerical
wall allowance after subtracting prior measured CPU; per subprocess <=600 seconds.
Reserve 100 seconds below CPU cap for reporting/scoring; limit BLAS/OMP to 1.
Estimated 900 seconds solves + 900 estimators + 600 independent checks =2400;
stop if remaining conservative allowance cannot cover the next task. Preserve
refusals, partials, receipts. No retry of a numerical failure with new parameters.
Fixed cases may complete despite a scored scientific failure to localize it;
no post-failure retuning or extra levels. This differs explicitly from review §9's
stop-on-first-estimator-failure recommendation; the fixed total ceiling still binds.

## EOS construction and admissibility (frozen before outcomes)

For resolution N use the original log-density base, same endpoints and exact
threshold nodes. Neutron distances d=10^(−j/(N/128)), j=N/128…D(N)N/128,
D(N)=min(8+floor(log2(N/1024)),11); muon D=8 unchanged. Cap12 experiment uses
min(original D,12). Nothing shifts or perturbs physical values. At N<=8192 this
is the exact historical construction, required bitwise table identity.
At N=32768 cap11 the last adjacent relative n separation is approximately
9e-14, versus approximately 9e-16 at historical depth13; it resolves distance
scales while avoiding representationally wasteful clustering. The minimum
unresolved density interval is 1e-11 relative, far narrower than radial and EOS
interpolation scales, and is explicitly probed by cap12. This is a physical scale
argument and prediction, not an ordering guarantee.
The set union of original base/onset grids is preserved (intentional shared nodes);
no P/ε duplicates may be removed. Every adjacent pressure and energy pair must be
strictly increasing both in provider units and serialized cgs units. Refuse if
not; retain the failed table audit. A representational failure means INCONCLUSIVE
for EOS estimator availability and falsifies the review's universal admissibility
prediction. It never authorizes weakening guards.
High precision: 70-digit free-gas momentum integrals/equilibrium at all historic
refused pairs and nearest pairs on either side of neutron onset in every new table;
verify strict order before and after binary64 rounding, inspect ULP gaps. Exact
binary input constants, not rounded decimal masses, define this oracle. Also
check analytic onset and ε(rho_b) inversion and ε+P=nμ, dP=n dμ independently.

## Estimators and frozen scoring rules

L: rν locates log(μ∞/mn) between metric nodes, μ∞ anchored at first pe node.
Lν=1e-36 ∫[0,rν] exp(−6ν)4πr²/sqrt(1−2m/r)dr. Compare production nB edge,
absolute L, paired L_i/L_0−1, and edge motion. Edge sensitivity
s=1e-36 exp(−6νe)4πre²/(sqrt(1−2me/re)Lν).

C: C/T∞=kB² conversion ∫ exp(−ν) Σ pF μ/(3 hc³) dV, split at neutron/muon
onsets with squared-distance quadrature. Compare nodal production, cache-removed
nodal trapezoid, nB/composition reconstruction and ν reconstruction. Cache factor
F(T)=((1−w)Ta+wTb)/T, w=log(T/Ta)/log(Tb/Ta), cache 160 log nodes 1e-5…1 MeV.
All 13 historical knots 10^(6+0.25k) K required. F is a common multiplicative
absolute bias; never count it as independent noise at each star.

Envelope: solve ε=1e10*c²/(MeV fm−3 to erg cm−3) by independent free gas;
μb defines νb=log(μ∞/μb). Tb=T∞exp(−νb); compare actual production energy-grid
index, log-density interpolation, and physical boundary. Same-profile identity
is algebraic evidence only. Propagate into original iron Ts∝g14^(1/4)Tb^0.605
and Ts∞=Ts exp(νR), for all 13 knots.

P1 absolute rν,Lν,Cν and Z entries: for retained 20,40,80,160,320k increments
require last three above 128eps*scale, same sign, last two ratios in (0,0.5).
q=max(last two ratios). Eligible prediction X640=X320+(X320−X160)/3;
sharp acceptance |X640−X320|<=2q|X320−X160|. Also report original loose band
2|d|/(1−q). Conditional tail estimate U=2q|d|/(1−q) at 320k, and propagated
q*U at 640k only IF future contraction <=q; no rigorous global enclosure.
P2 paired L,C: use same paired eligibility (arithmetic floor => INCONCLUSIVE for
order); |D640−D320|<=max(2q|D320−D160|,1e-7), AND conditional paired U<=1e-7.
If paired ladder is below floor, report stability band <=1e-7 separately, never
call its convergence validated. Physical thresholds stay 1e-4, not retuned.
P3 at640k and each jitter: |(Dprod−Dν)−s0*((re_prod,i−re_prod,0)−(reν,i−reν,0))|
<=0.02*|Dprod−Dν|. Tb paired snapping discrepancy must match linearized
Δν_b−Δν_node within 2%; exact exponential identity separately labeled tautological.
P4 secondary: D_i predicted=(x_i/x20)D20, x from frozen physical target Bs;
|observed−predicted|<=1e-9+1e-3|predicted| on each retained resolution, L and C.
P5 EOS: low-N exact identities, strict admissibility; freeze predictions on
2048/4096/8192 discovery with P1 eligibility and formula, recursively predict
16384c and32768c. 32768 sharp score uses fixed discovery q and measured 16384
increment, not refitted q. Report65536 consistency. Depth12 vs11 relative
change <=max(128eps,1e-12) for Lν,Cν,rν,R,M,Nμ. If eligibility fails, INCONCLUSIVE.
P6 cross additivity: |X640,c−X640,8192−X80,c+X80,8192| <=
0.1 max(|X640,8192−X80,8192|,|X80,c−X80,8192|), Lν,Cν,rν,R,M,Nμ.
P7 all six jitter levels: ν estimators L,C,edge relative range <=1e-9 absolute;
paired drift range <=1e-10. P3 applies at every phase. Tb snapped paired variation
<=1e-10 within groups sharing the same production index at BOTH targets.
Any numerical miss is FAIL; arithmetic floor, missing case, budget stop, or
unrepresentable EOS is INCONCLUSIVE. Report every score separately; overall PASS
requires all mandatory scores eligible and passing plus independent checks below.

## Required independent cross-checks and error accounting

- EOS high-precision order, onset threshold, thermodynamic identities; profile
  μ∞ constancy from independent species energies and (ε+P)/n; metric ODE identity.
- Independent analytic-EOS TOV using different coordinate/integration, at fixed B,
  two tolerances. Report solver/tolerance, support location, absolute and paired
  spread separately. A difference is not by itself a rigorous bound.
- Independent pressure-based support/base location (avoid nB-linear edge), and
  anchor-node sensitivity. Check physical edge lies within node bracket.
- Quadrature order24/48 on three r640000 profiles (same-background check only).
- Preserve inherited chemical background/tail diagnostics; compute ν G,Z at three
  r640000 profiles and production-style G,Z there only. No tangent-stencil solves;
  infer the three-target structure contribution from Nμ radial paired changes,
  explicitly not a full tangent bound.
- B residual/bracket and physical target bytes, common deterministic errors,
  EOS conditioning, source/binary/manifest identity, thermal-cache factor all logged.
- Independent oracle discrepancies must fit the sum of applicable conditional
  estimates and measured solver/EOS errors to support disposition A. If any
  missing common systematic cannot be bounded, disposition B. Material inconsistent
  equations/physical onset => D; production numerical correction need => separate C
  proposal, never an edit to production. Historical covariance stays historical.

No BNV trajectory, certificate promotion, owner/threshold/baseline/EOS-data edit,
master merge, automatic ADR ratification, or 21-star expansion. Deliver incremental
bundle based on 4f3cf97 for the office Mac unless established non-force push workflow
is available. No direct push credentials assumed.
