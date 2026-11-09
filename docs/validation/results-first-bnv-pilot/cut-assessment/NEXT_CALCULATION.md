# Single smallest remaining calculation

**One bounded finite-cut surface-transfer witness, with its first variation and
fixed-B remapping.** This is a local static calculation, not another equilibrium
sequence, an interior solver comparison, or a trajectory. Its necessity is the
missing finite transfer identified in DERIVATION.md, not threshold headroom.
No new stellar calculation was launched in this assessment.

The reuseable inputs are the retained central/offset states at 0, 10, 20,
their existing numerical allowances, the original cut law, the analytic pe
terminal equations and the already accepted C-CONT slope/curvature constraints.
All are sufficient to formulate the following narrowly specified calculation;
they are not already a computed transfer witness.

1. Construct the terminal flow between two nearby pressure cuts at fixed Hc,
   transporting R, M, B, Ne, Cstar's raw integral, the ee susceptibility integral
   and K. The neutron/muon core integrals are unchanged during this local flow.
   Evaluate the exact nonlinear thermal relation and chemical inverse afterward.
   Include the state/parameter uncertainty of the retained starting states.
2. Differentiate this local map once. Combine it with the existing along-law
   family and solve the scalar fixed-B matching equation with a separated
   denominator. Use the on-law identity in DERIVATION.md to propagate C-CONT's
   count curvature into the path-tangent bound; do not require the rejected
   mixed-count Richardson estimate as an input.
3. Produce an explicit remainder/enclosure on the whole cut strip
   |delta ln Ps|<=2.1e-12 and the full x interval. Include baseline and final
   endpoint endcaps: fixed-B remapping can require an extremely small extension
   of the on-law family. Use the retained central-offset neighborhoods to
   establish that extension with a supported remainder; do not silently extend
   C-CONT's domain. If that extension cannot be established, report exactly that
   failure before considering any new interior calculation.

The prospective acceptance criteria are:

- A finite observable-wise U_Q, baseline U_Q(0), path-tangent U_T and endpoint
  contribution, with no unknown term replaced by zero. Each must follow from
  an analytic majorant/interval remainder or the accepted S convergence
  standard; sensitivity products and tiny matching residuals alone do not pass.
- Positive pressure, pe-only terminal topology, D>0, separated B_h/fixed-B
  mapping and chemical inverse denominators throughout the local tube. A
  domain check is required for both signs of the perturbation.
- Agreement with the independently checked first-variation identities in the
  zero-perturbation limit. For any numerical rather than analytic remainder,
  at least three levels, signed consistent increments, formal order +/-0.5,
  and subordinate arithmetic/iteration errors, as required by S. Precision
  checks must concern the new local map, not just decimal reformatting of old
  stellar outputs.
- Propagate the resulting bounds through the original paired difference,
  all-row fit operator, nonlinear baseline normalization, and existing ledger.
  Qualification requires every resulting budget below its unchanged historical
  threshold. That inequality is an acceptance check, not the source of U_Q.
- Preserve the 42 historical failures, 353 expectation violations, original
  matrix, accepted C-TAN/C-INPUT terms, and metric consumer limitation. No
  operational certificate or trajectory authorization follows automatically.

This witness can reuse all retained interiors. New global equilibrium solves
have **not** been shown indispensable. Enlarged arbitrary cut offsets and an
ordinary double-precision re-shoot at the literal 2.1e-12 residual are not the
minimal test: the induced Hc shift is approximately 1.2e-26, below the resolution
of the retained central coordinate. Do not start a 21-star refinement campaign.
