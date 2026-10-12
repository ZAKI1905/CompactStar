# Legacy T2 remains scientifically unresolved

The imported review §6.3 measures Linux/Mac neutron, Lambda and Sigma-minus
lifetime differences 2.2%, 2.5%, 4.3%; the within-Mac Debug/Release control is
0.59%, 0.60%, 0.66%. These are material; no equivalence pass is claimed.

Call-path separation, not a filename inference:
* `Core/src/TaskManager.cpp:581` converts each constant-B contour to Curve2D,
  intersects it with critical_curve, requires exactly one intersection, and
  bisects at that point. The exported B_conts then feed the legacy analysis
  (`Microphysics/BNV/Analysis/src/BNV_Analysis.cpp:217`, Evolve).
* Controlled Phase-6 constructs an `EquilibriumBaryonTangent` from
  `ParticleNumbers` and `EquilibriumSequenceNumberDerivative` at the Structure-1
  midpoint (`tests/bnv/campaign_fixture.hpp:21`). Its source amplitude is a
  declared fractional drive times B0/year, not a lifetime read from T2.
* `tests/bnv/breakpoint_fixture.hpp:31` creates
  `UniformProperNeutronSinkHistory`, product-fate partition and
  `FrozenControlledBnvRunContext`; its coefficient input is mu_B_inf from the
  native-bound coefficient table. `ControlledBnvSecularDriver` consumes this
  explicitly supplied history. The matched control has an exact-zero history.
* Source/include scans of `CompactStar/Physics/BNV` and `tests/bnv` find no
  TaskManager, MicroBNVAna, BNV_Analysis, B_conts or BNV_tau consumer. Symbol
  inspection of the linked readiness executable is recorded separately.

Thus this controlled fixed-background Phase-6 scope does not consume the legacy
T2 lifetime pipeline. Its unresolved conditioning is not a prerequisite for that
limited qualification. This statement does not extend to realistic lifetime-driven
campaigns, other CompactStar clients or the whole library.

Specific follow-up: independently refine the 19x19 sequence and critical/constant-B
curves, measure intersection crossing angle/Jacobian and root displacement under
both curve uncertainty envelopes, bracket near-tangent and multiple-root cases,
and propagate death-point uncertainty through each species lifetime integral.
Use held-out central-density ranges and the within-Mac control. Establish a
lifetime uncertainty budget before proposing any cross-platform acceptance.
Do not change radius deduplication or saddle handling in this correction.
