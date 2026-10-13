# BA13 onset-aware blind pilot — disposition B

**The proposed estimator failed prospective falsification. The dominant onset and
snapping mechanisms survived, but their proposed numerical bounds did not.**
No validated estimator suitable for immediate 21-star qualification results.
Further bounded evidence is required. This is an unratified validation-only
numerical-method investigation, not a production correction or certificate.

## Evidence and blindness

The review's 11 files and the prior study's 1022 retained artifacts authenticate.
The read-only review, historical fitter, frozen certificate, accepted owners,
EOS source/data, thresholds and baselines remain unchanged. Initial clean entry:
`4f3cf97462ca030f715e26808938f781a2a40945`. Descendant branch:
`candidate/linux-ba13-onset-pilot-20261013`.

Full predeclaration: `8bb72e2`; numeric radial/secondary predictions: `e657efb`;
EOS predictions and executed estimator tools: `650871e`, all before Slurm 14690.
The scorer was implemented while the fixed cases ran, from committed rules;
no bands were changed. Discovery and secondary data are explicitly distinguished:
reviewed 20k–320k and old EOS profiles are discovery; targets 5/15 are secondary
transfer tests on authenticated existing profiles, not fully blind stars.

All 31 new stars completed: three 640k, nine capped-EOS, three depth12, fifteen
jitter, one 640k×16384c target 20. Six analytic-EOS/common-B oracle solves also
completed at the two declared tolerances. All numerical work ran on Slurm with
one CPU. Total is about **795.5 CPU-seconds (0.221 CPU-hour)** including conservative
startup charge and post-outcome diagnostics; build about 9.44 CPU-seconds separate.
Exact accounting and peak RSS are in evidence/resource-summary.json and the
compressed resources ledger. The 3600-second ceiling was not approached.

## Frozen outcomes

Source: evidence/scores.json (193 individual score records; counts include controls,
discovery nominal pinning, and secondary tests, not 193 independent holdouts).

| Test | Outcome |
|---|---|
| P1 absolute radial edge,L,C and three Z entries | 18/18 PASS |
| P2 paired convergence/tail eligibility | 8 INCONCLUSIVE; all four separate stability screens pass |
| P3 edge-motion and Tb snapping mechanism | 28/28 PASS, including four nominal discovery controls |
| P4 targets 5/15, L/C, five levels | 20/20 PASS, secondary |
| P5 16384c predictions | 9/9 PASS |
| P5 withheld 32768c | 7 PASS, **2 FAIL**, both target 20 L/C |
| P5 65536c consistency | 6 PASS, **3 FAIL**, target 20 L/C/edge |
| P5 depth12 and table identity/order | 18/18 and7/7 PASS |
| P6 radial×EOS additivity | 6/6 PASS |
| P7 phase | **9 FAIL** (L/edge); C passes; both envelope group screens pass |
| High-precision pair ordering / independent oracle completion | 32/32 and6/6 PASS |

The two 32768c misses are small but remain failures: L increment 7.6568e-44
against 7.3837e-44; C/T increment 5.6388e18 against 4.9572e18. The subsequent
65536c increments also violate the fixed contraction forecasts. These increments
are about 4e-12 relative, where EOS changes compete with background/solver and
floating-point differences. That is a candidate explanation, not permission to
replace the failed limits by an arithmetic floor or retune q.

C radial increment ratios are 0.25005–0.25054; Z ratios 0.25015–0.25119.
The edge ratios are 0.27259–0.27602 and L ratios 0.28663–0.29202: compatible
with the declared sharp bands, but already depart from an exact quarter ratio.
The paired discovery ladders do not meet the frozen order-eligibility rules.
Their small observed increments establish stability on these cases, not paired
convergence order or an asymptotic physical bound.

## Physical checks and what failed

The independent derivation in EQUATIONS_AND_LIMITS.md establishes the neutral
chemical-potential identity from the first law plus TOV, including the pe branch.
The 70-digit independent EOS confirms ε+P=nμ and dP=n dμ, analytic neutron onset
n=7.3567289037328326e-9 fm^-3, and envelope-base μb=939.6306205322022 MeV.
The exported production onset differs by one density ULP. Physical branch and
support brackets are consistent; no inconsistent EOS/metric equation was found.
Nominal μ∞ varies about 4.01e-10 across a profile; this residual is measured rather
than erased by imposing the reconstruction. Provider enthalpy/species μ agree
near 1.3e-14. Metric secants are diagnostic approximations, not equation residual bounds.

The independent enthalpy solver uses direct momentum-integral EOS values, its own
common-B shooting, and no production table or radial profile. At 640k its tighter
run differs by about 1.8e-10 in L and1.2e-10 in C, relative. However changing its
prescribed tolerance moves L by 7.05e-10 and C by 2.59e-9. Thus fine-run agreement
is useful independent evidence, **not an independently certified reference at
1e-10**. Large nearly common oracle errors also show why paired agreement alone
cannot establish absolute accuracy. Original target-B residuals, brackets,
central-density changes and oracle alignment are retained in the JSON evidence.

The fixed 1e-9 absolute phase band fails for all three L and edge ladders:
L ranges 3.002–3.131e-9; edge ranges 1.242–1.295e-9. The paired L range reaches
2.925e-10 against 1e-10. C stays inside both phase bands.
The analytic inverse-interpolation error is

    r_linear − r_true ≈ −ν''/(2ν') θ(1−θ) h².

Using ν',ν'' from TOV at the physical onset explains the phase-dependent edge
error. Subtracting that leading term **post hoc, in a separate diagnostic only**
reduces L phase ranges to 3.73–4.35e-11 and edge ranges to 2.50–2.79e-11.
See evidence/geometry-diagnostic.json. P7 remains FAIL. This is discovery for a
future curvature-aware estimator, not a validated correction. Nearby-anchor
μ∞ changes are about 1.4–1.8e-14 in nominal jitter profiles, far smaller than
the edge phase variation; central/radius changes are also recorded separately.

## Frozen-relative versus absolute accuracy

At nominal settings the target 20 L drift is −1.22664e-6 by production edge
reconstruction, +1.80923e-6 by ν reconstruction, and +1.80921e-6 by the independent
analytic-EOS star. The corrected sign is supported independently. Production L
is about9.29e-5 high at target 0 relative to the same-profile ν estimate, while
the ν estimate itself is about 1.71e-8 high relative to the tighter oracle.
The production paired discrepancy is about 3.04e-6, roughly 3% of a 1e-4 physical
threshold. These are observed discrepancies, not a certificate enclosure.

All 13 thermal knots are evaluated. C cache bias reaches 6.55536e-4 absolute and
cancels in equal-knot ratios; the cache law reproduces to about8e-15. After
removing it, nominal nodal C differs from ν C by 1.82e-7 at target 0, with a
paired target 20 discrepancy1.10e-8. Separate 640k order 24/48, nodal trapezoid,
linear-species and nB-equilibrium reconstructions are retained. Agreement of
quadrature orders validates integration of that reconstructed background only.

The independently determined envelope base shows nominal Tb offsets 4.12–4.23e-6
and surface-temperature offsets 2.49–2.56e-6 relative. At target 20, paired Tb
error is 1.0713e-7 and paired local/redshifted Ts error about 6.4811e-8. At T∞=1e6 K,
target 0 Tb is 1081985.435 K at the selected node versus 1081980.977 K from the
physical-boundary reconstruction. Absolute and paired errors are not interchangeable.
The two literal repeated-index phase groups did provide comparisons and passed;
none was awarded a vacuous PASS.

The capped family preserves exact 2048/4096/8192 tables and strict pressure/energy
ordering in both provider and serialized cgs arithmetic. High precision verifies
all audited dangerous pairs before and after binary64 conversion. Historical
refused pairs have distinct correctly rounded outputs: loss occurs in numerical
intermediates, not unavoidable equality of the physical output values. The
historical 32768 refusal remains a refusal. Every-node extended-precision ULP
coverage and a proved interpolation/conditioning enclosure remain missing.

## Justified uncertainty and disposition

There is **no new unconditional global physical error bound**. The successful
radial forecasts justify only conditional reconstruction estimates on this pilot:
at 640k approximately 8.1e-10 relative for L and 4.23e-10 for C, assuming continued
frozen-q contraction and excluding common EOS/solver error. They cannot cover
the failed phase experiment, missing EOS floor, or absolute oracle uncertainty.
All measured finite-case spreads, errors and conditional estimates are retained;
none replaces historical errors, thresholds or covariance.

The three-target Z check supports the onset-coordinate explanation of its radial
profile dependence. G conditioning is about 2982; inherited background/tail
characterizations are preserved. Common G errors cannot simply be subtracted
before a nonlinear inverse; their changing Jacobian and remainder also matter.
Nmu paired radial changes supply a three-target secant structure diagnostic only,
not a completed tangent-stencil uncertainty. No extra tangent solves were made.

**Primary disposition B: further bounded evidence. A is not supported yet.**
No material inconsistency in the accepted physical equations was found, so this
pilot does not establish disposition D. The known production edge sign error and
absolute thermal/envelope biases **do warrant a separately governed numerical-method
amendment** (C as a production-change proposal). This pilot does not prove that a
production correction is logically necessary before frozen-relative qualification;
a validated error enclosure could instead suffice. No correction was implemented.

All non-exact BA13 classes remain unqualified for full physical uncertainty:
Ltilde, Cstar/thermal/surface/envelope/gravity, Z and chemical background/tail,
particle numbers and tangents, radius/structure, μ/P0/P1 and normalized profiles.
The inherited exact-P2 result remains separate; no new certificate field is issued.
Historical fit covariance remains historical. Inter-sample/inner-half remainder,
common deterministic systematics and all 21-star coverage are unresolved.
ADR-0021 Option A and interpretation (2) with frozen right-hand thresholds remain
proposals requiring owner adjudication, not ratified policy. **Expansion to 21 stars
is not justified on this estimator.**

## Costed continuation, not executed

Predeclare a new three-target experiment for the curvature term and EOS/solver
floor, retaining all present FAILs. Candidate ceiling 0.5 CPU-hour: eighteen new
jitter stars 80006…80011 at <=15s each (270s), three 1280000 stars at <=100s (300s),
reconstruction/independent solver checks <=330s, whole 32768c high-precision audit
and table-interpolation conditioning <=600s, scoring/overhead <=150s: total 1650s,
150s reserve. Benchmark the high-precision audit first and stop if its extrapolated
cost exceeds its allocation. These are estimates, not authorization or guaranteed
bounds. No 21-star expansion until fresh prospective validation and an independent
common-error enclosure succeed. Background/tangent/nonlinear remainder completion
needs its own subsequent scope and cost.

Changes are confined to the new validation tools/tests and this evidence directory.
Authority: GOVERNANCE §§2–5, INV-02/03/14/15, accepted ADR-0002/0012; proposed
ADR-0021 does not acquire normative status. No BNV trajectory, production change,
threshold/baseline edit, master merge, certificate promotion or direct remote push.
