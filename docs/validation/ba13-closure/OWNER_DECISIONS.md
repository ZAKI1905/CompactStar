# Owner decisions for ADR-0021 — no acceptance implied

The proposed default is Option A, interpretation (2): retain the historical
fitter and fixed right-hand thresholds; propagate correlated numerical error,
fit and interval remainder on the left. The original covariance supplies no
physical uncertainty. No certificate is promoted. The present result is E:
insufficient complete evidence, even though measured terms are small.
See [REPORT.md](REPORT.md) and [dispositions.json](evidence/dispositions.json).

| Decision | Proposed default and justification | Effect on status | Alternatives and implications |
|---|---|---|---|
| 1. Explicit inter-sample remainder | Require a named full-window modulus. K*dx²/8 follows from linear-interpolation error if |Q_xx|<=K; current factor 10 times sampled curvature is conditional. | Interval rows remain UNQUALIFIED pending an accepted, validated construction. | A strictly sampled certificate is narrower than a continuous frozen window. An analytic enclosure is stronger but costs further mathematical work. |
| 2. Accepted evidence classes | Admit established bounds and validated numerical estimates under a cited prospective standard. A conditional estimate needs named explicit acceptance. Propagation algebra cannot upgrade an uncertain input. | M6 enforces this. Accepting the vocabulary alone does not fill missing stage, arithmetic, EOS/center or derivative remainders. | Requiring fully rigorous enclosures increases work. Accepting a named conditional sidecar changes the scope of the claim and must remain visible. |
| 3. Certified object | Exact physical equilibrium response plus production numerical uncertainty, consistent with preflight §11's full-drift requirement. Sampled production representation is an error term. | No fictitious production jump is charged to physical curvature. Positive B_H/domain continuation and physical curvature still need eligible evidence. | Certifying intermediate production output additionally requires adaptive-controller and moving-surface coverage. Equal sampled indices and large margins alone are insufficient. |
| 4. Record format | A validation sidecar bound to the immutable historical certificate/source hashes. It can report complete uncertainty without modifying historical fit semantics. | Present sidecar is an assessment with explicit refusals, not a runtime certificate. | A successor runtime record requires a distinct identity, monitor-binding validation and separate authorization. No successor fitter is presently needed. |
| 5. Transfer across targets | Require a quantity-specific sensitivity/enclosure or an explicitly accepted continuation hypothesis. Independent endpoint addition removes a signed-correlation assumption but does not prove the endpoint magnitudes. | M2 remains UNQUALIFIED. All interior nodes agree; the moving final node must be treated separately. The inherited chemical background magnitude also needs acceptance/validation at each target. | Direct 21-star solving would violate this task's scope and still would not prove continuous coverage. A bounded static sensitivity argument or limited prospective withheld-target check under an accepted standard is smaller. |

All five are scientific scope/evidence decisions. The implementation makes none
on the owner's behalf. Accepting these defaults would still require the
specific missing numerical inputs named in the report before any new PASS.
Changing a historical threshold, fit window, source model or certificate value
is not an alternative within this task.
