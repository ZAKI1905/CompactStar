# RF-BNV-P2-20261108 — proposed bounded-pilot decision

**PROPOSED; NOT ACCEPTED. Owner, date, approval reference and reviewed delivery
commit: UNSET. No registry entry or operational certificate is authorized by this
document's creation.** Change class: documentation of a proposed scientific-scope
decision. Authority: GOVERNANCE §§1–5; ADR-0015 §§9,14–16; ADR-0017 Amendment 1;
accepted implementation preflight §11; ADR-0021 remains proposed.

The single next owner action is to accept, reject or amend the following package.
It is deliberately a decision about one P2 benchmark, not full BA13 or the full
campaign. Evidence source is `bb0f9c5af8ebaf86b5a59aece5551bb9dafbd264` and the
independent review copied unchanged as `evidence/INDEPENDENT_REVIEW.md`, SHA-256
`606c9cb0f4ea005b7439aeddf0a3fbeeacca2e7a972e34d74ddf32dccd34d412`.
The review's recommendation is not owner acceptance.

## Exact proposed acceptance text

“For RF-BNV-P2-NATIVE-BOUNDED-v1 only, I select ADR-0021 Option A and interpretation
(2), and accept A–E in `../ba13-final-decision/OWNER_DECISION.md:14` with the
following explicit amendments and limitations:

1. **A — interval:** retain the full depletion interval [-1e-6,0], all historical
   thresholds, grids, norms and fit windows. Close the attained endpoint sliver
   using C-CONT; do not reduce the certification interval to fit this pilot.
2. **B — evidence:** accept standard S and each of C-REF, C-TAN, C-CONT, C-INPUT,
   C-SUPPORT, C-CUT and C-MUN, with the quantitative statements in the authenticated
   review Appendix A (`evidence/INDEPENDENT_REVIEW.md:536`). The qualification
   wording is **qualified under named validated numerical estimates and itemized
   conditionals; not a rigorous enclosure**. Generic acceptance of B alone is
   insufficient. No unknown error is replaced by zero.
3. **C — object:** certify physical equilibrium response including production
   error. Define ln(P_s) linear in fractional depletion through the retained cuts
   (reviewed slope 8.3839; deviations <=2.1e-12). The tangent object is the
   derivative along that surface law. Keep the current production tangent and
   charge its fixed-cut/path discrepancy through C-TAN, without cancellation or
   silently modifying t. C-MUN applies only to the eta=0 equilibrium reference;
   the actual eta-dependent potential remains evaluated separately.
4. **D — restricted metric scope:** accept the P2 consumer inventory in
   `METRIC_CONSUMERS.md`. Replace only the independent nodewise metric claim by
   coverage through the listed functionals. Do not claim full BA13, absolute
   stellar accuracy, or nodewise accuracy. Retain all 22 historical runtime rows,
   including metric_structure, and their unchanged refusal behavior. A separate
   authenticated sidecar reports conditional utilization; the runtime continues
   to consume the unchanged historical utilization. This authorization does not
   turn its metric row into a nodewise PASS.
5. **E — reference:** accept quantity-specific achieved-B transfer with matched
   surface pressures, the existing prospective holdouts and the named conditions.
   Two holdouts are not a proof of a continuous supremum.
6. **Linux pilot scope:** permit the reviewed native candidate and authenticated
   existing Linux dependency closure for this bounded mathematical experiment
   only, with no whole-library/platform equivalence, ADR-0020 blanket ratification,
   historical Mac pin replacement, baseline promotion or master merge. Authorize
   an explicit native run-identity binding candidate after sidecar completion;
   it must preserve exact owning-profile, byte-currentness and refusal checks.
   The native artifact hashes in `evidence/preflight.json` identify the candidate,
   not self-issued scientific authority.
7. **Release sequence:** authorize the <=60 CPU-second mechanical static sidecar
   completion recommended by the review, then a short acceptance check. Record
   my selections in a separate registry commit before the assessor consumes them.
   Preserve all prior FAIL/INCONCLUSIVE/UNQUALIFIED records. Only after that check,
   exact native context binding and the affected prerequisite tests pass may the
   separately declared paired pilot in `EXPERIMENT.md` run. Failure of any gate
   stops it. This accepts no trajectory in advance and authorizes no 500 kyr run.”

## Quantitative condition index (all pending)

The exact statements and retained per-row values are in Appendix A of the copied
review; they are incorporated above, not newly inferred from available headroom.

| Selection | Quantitative content requiring explicit acceptance |
|---|---|
| S | >=3 levels, signed consistent increments, order within 0.5 of formal; finest tail or sensitivity with a convergence law; subordinate iteration/roundoff; independent corroboration where available |
| C-REF | Charged GL4/EOS/centre estimates; stage/arithmetic remainder <= finest increment; shared residual <=10 times the specified agreement; special Ltilde_Me X1–X2 value 9.1e-12 |
| C-TAN | Finest ratio tail + twice narrow-stencil tail + EOS/centre + moving-cut residual, all absolute; denominator relative error <=1e-14 |
| C-CONT | Per-group K from retained conditional intervals; positive B_H >=0.999 minimum sampled secant; fixed topology; endpoint derivative/sliver allowance |
| C-INPUT | Per-gap coordinate errors <=5.1e-12; charged B0 normalization; all-row fit factor 5.185 and declared second-order remainder; historical weights frozen |
| C-SUPPORT | No onset crossing under the accepted monotonicity condition |
| C-CUT | Declared surface law and <=2.1e-12 log-pressure residual propagated through surface sensitivity |
| C-MUN | Equilibrium redshifted mu_n = mu_B only at eta=0 |

The review's failure multipliers and scope limitations remain part of the record
(review §§3.9,5 and Appendix A). They are not new acceptance thresholds. In
particular, same-knot Cstar error cancellation does not establish absolute thermal
accuracy. The sidecar must propagate all these terms jointly and keep numerical
and scientific qualification distinct.

## Recording the answer

An explicit acceptance of **RF-BNV-P2-20261108 as proposed**, naming the reviewed
delivery commit, selects the entire text above. Alternatively list amended or
rejected selections. Record owner/date/approval reference verbatim; do not infer
them from this task's permission to investigate. The existing registry remains
empty until then. If declined, retain the historical reproduction claim and stop;
do not begin another numerical study by default.

There is also an identified implementation gate, not a missing scientific study:
the native thermal context is not yet wired to an accepted identity. All four
native input hashes differ from the fixed Mac pins and native Bdot differs from
the old literal. `METRIC_CONSUMERS.md` records the exact locations. Acceptance
above permits preparing that narrow candidate; it cannot make the existing
executable runnable by declaration.
