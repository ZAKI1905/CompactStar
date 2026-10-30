# ADR-0021 proposed owner ratification record

**Status: PROPOSED — NOT ACCEPTED.** This is the exact proposed record for owner
review, not a signature inferred from permission to investigate. Authority and
historical ambiguity: `docs/adr/ADR-0021-ba13-static-uncertainty-completion.md:1`;
original norms: `docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:691`.
Candidate source and evidence are bound by this directory's manifests and final
handoff commit. Owner name, approval date, approval evidence and reviewed commit
are all **UNSET**. The acceptance registry remains empty.

## Text proposed for explicit ratification

“I select ADR-0021 Option A and interpretation (2), subject to the following:

A. **Full interval.** BA13 evaluates the physical equilibrium response throughout
[-1e-6,0] in fractional baryon depletion, including nonlinear and numerical
remainders between all 21 stars and any attained-coordinate endpoint sliver.
The original right-hand thresholds, norms and positivity/support requirements
remain unchanged.

B. **Evidence standard.** I permit analytically established bounds and independently
validated numerical estimates with declared error models, scope, numerical
uncertainty and falsification tests. I do not require every estimate to be a
rigorous interval enclosure. A conditional estimate is eligible only when I
explicitly accept its named quantitative assumptions and residual limitations.
Acceptance of this standard alone accepts no conditional value or unknown bias.
Passing tests, observed differences and illustrative allowances are not by
themselves error bounds.

C. **Certified object.** The object is the physical equilibrium response on the
declared depletion interval, including the numerical error of the historical
production samples. Continuous physical response is distinguished from historical
adaptive mesh, support-edge pinning and envelope-node selection. I do not require
continuous reproduction of those discrete numerical rules. The equilibrium
mu_n_profile identity applies only to the equilibrium/reference redshifted field,
not the separate eta-dependent actual chemical potential.

D. **Record.** An authenticated validation sidecar may supplement the unchanged
historical certificate, with exact certificate/fitter identities and provenance
for every contribution. It cannot silently overwrite historical FAIL,
INCONCLUSIVE or UNQUALIFIED records. A scientific qualification or operational
certificate attachment requires separate review and explicit authorization.

E. **Target transfer.** All-target comparison uses independently checked X2
interpolation at achieved B with correctly matched surface pressures and
quantity-specific continuation evidence. Three reference stars or two passing
holdouts do not automatically establish a continuous bound. Tangents, chemical
rows, thermal knots and metric profiles must retain their actual scopes.”

## Clarification versus new requirements

A and C resolve the existing preflight §11 phrase “full drift plus numerical and
nonlinearity terms” and the object of the drift. They do not tighten its thresholds.
B is a substantive evidence-policy selection; generic acceptance must be separated
from itemized assumption acceptance. D authorizes a new record format, not altered
physics. E selects a validation method; independent error accounting follows from
A/B, while a particular extra holdout, factor-ten continuation envelope, category
budget, 10%-of-share refinement stop or observer order gate is an investigation
criterion, **not a new historical acceptance requirement**.

The 12 inherited category shares may guide work, but are not additional physical
limits from preflight §11. Their violation alone is not disposition C. Requiring
an independent solver for every tangent would also be an added requirement;
well-supported derivative remainders under B can suffice without one.

## Separate conditional selections (all pending)

- C-REF: quantitative joint bound/estimate for remaining X2 common error, including
  stage, EOS/inversion, centre/surface and arithmetic over the actual quantities.
  No value is proposed as established merely from the illustrative 1e-11 ledger.
- C-TAN: a quantitative enclosure of the remaining joint numerator/denominator
  derivative residual and stencil/moving-cut remainder in absolute tangent units.
- C-CONT: quantity-specific continuous curvature/positive-B/topology and cut-law
  hypotheses, their validation domain and endpoint completion. The inherited
  factor-ten estimates remain conditional.
- C-INPUT: bounded achieved-coordinate, normalizer, fit design/weight and surface
  pressure-input perturbations, including nonlinear response.

The final matrix gives required accuracy regions for these selections; those
regions are not empirical estimates. Full metric_structure additionally requires
nodewise data and independently validated interior interpolation. Acceptance of
A–E alone cannot fill that missing measurement. A narrower claim excluding this
or any other group requires a separate explicit D decision and cannot be called
full BA13 qualification.

## Recording an actual owner decision

Append the reviewed full commit and manifest hashes, owner/date, an exact quoted
approval reference, A–E accepted/rejected/amended, and each conditional identifier
with its accepted quantitative statement and limitations. Leave unselected items
pending. Do not populate authority fields or mark this record ACCEPTED until
that explicit authorization exists. No signature or approval is supplied here.
