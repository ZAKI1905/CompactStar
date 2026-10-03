# Phase-6 breakpoint candidate: successor validation predeclaration

Status: proposed successor criteria, written before new scoring, 2026-10-03.
Authority: owner's bounded implementation instruction; ADR-0017 Amendment 1 candidate.

1. Authenticate the imported independent review and immutable fixture inputs. Reuse
   valid historical reproduction/order/perturbation evidence; exclude invalid files
   according to evidence/obs32-independent-review/validity.json. Never repair run3.
2. Preserve exact old obs-32 and Arm-E reproduction separately. New main histories
   need not match the old step sequence. Use common physical times for new comparisons.
3. Exercise analytic coupled piecewise IVPs in both directions, knot starts, direction
   changes, stage excursions, independent event/endpoint convergence, missed crossings,
   wrong cells, stale rollback derivatives and work exhaustion.
4. Replay all 37 archived source/control crossing intervals at BASELINE, REFINED,
   ULTRA. Bound windows: archived source step20 to 1e12 s; archived control step49
   to its review window endpoint; final-decade source stability window of at most
   100 main steps; original 0..462269531250 s case. No full campaign integration.
5. Accuracy: per component absolute error plus qualified reference uncertainty <=
   atol_tier + rtol_tier max(abs(candidate),abs(reference)). Reference uses own events,
   tighter rk8pd with explicit step subdivision, and independent RKF45 agreement.
   Reference uncertainty is at least twice the refinement discrepancy and a binary64
   roundoff floor. Do not infer accuracy from mere O1/O2 identity.
6. Report original BA12R c1/c2/c3 and ledger gates unchanged. Classify boundary-limited
   segments from clipping/event counts and matched cell intervals before inspecting
   state differences. Successor acceptance requires accuracy in item5, state and
   ledger budgets, and explicit within-cell refinement/independent-method agreement.
7. Test schedule invariance bitwise; exact and nextafter-adjacent boundary outputs;
   smooth intervals from identical initial states; reconstructed diagnostics/R18/R20.
   Reconstruction allocation stays 0.20 F_i; R20 accuracy stays 2e-4 with the existing
   reconstruction allocation. No extra event points in the quadrature grid.
8. Tail: report accepted/rejected/RHS/event work, minimum/maximum step, elapsed time,
   finite progress and reference errors. Legacy rejection gate remains separately
   reported; no relaxed threshold. An unexplained stability/accuracy failure blocks.
9. Evidence: isolated task output, compact rows, checked writes, flush/close and
   read-back schema/row/time/completeness validation before publication. Inject a
   controlled write failure. Require available space > expected output plus 1 GiB;
   expected compact evidence <100 MiB and disposable build <2 GiB.

The empirically tested propagation envelope applies only to these bounded intervals,
fixture/domain and methods. It is not a universal event theorem or campaign approval.
Any failed budget must remain visible; ordinary code/build/test repairs are authorized.

Implementation detail predeclared before ledger scoring: the retained BA12R ledger
accuracy scale is 1e-9 G_P for REFINED/ULTRA, with original 1e-11 G_P floor and
0.10 contraction reported separately. Successor comparison to a qualified reference
uses the same 1e-9 G_P accuracy scale; state accuracy retains each tier's own scale.
The independent reference uncertainty must fit 0.05 D_ULTRA. U_event_i is the
full coupled endpoint envelope 2 max(|O2_i-H_i|, |H_i-I_i|)+64 ulp(|O2_i|), where H
is tighter half-step rk8pd and I tighter quarter-step RKF45, each with independent
event solves. Units are dimensionless for x and MeV for each eta. Since differences
are measured at the observation, the endpoint propagation factor is one; intermediate
coupling is integrated explicitly, not approximated by a scalar thermal decay.
This empirical envelope is tested against further step refinement and independent
endpoint/reference errors in the declared domain. Diagnostic propagation evaluates
all eight corners of the coupled state uncertainty box with a factor-two margin;
it must agree conservatively with independent diagnostics in the bounded tests.

Engineering runtime correction: profiling the retained full-currentness calls
identified character-by-character copying in ParticleNumberResponse.cpp::contents
and ChemicalResponse.cpp::bytes. Replace only these file-reading loops with checked
block reads. Every call still reads
and compares all bytes. Empty/binary files, exact/adjacent block lengths, content
mutation and missing-file refusal are tested. This is an explicit engineering
exception to source-byte non-change in shared Analysis code; all Phase-5 mathematics,
solver behavior and historical baselines remain unchanged. Historical Arm E and
obs-32 must reproduce exactly after this change. It does not address the deferred
Phase-5 knot-integration issue or weaken currentness to metadata checks.
