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
