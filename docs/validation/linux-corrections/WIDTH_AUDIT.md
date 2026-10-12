# Per-quantity extended-precision audit (candidate)

Only the inverse-density method changes width. No wholesale replacement.

| Site | Quantity and conditioning | Linux validation / width consequence |
|---|---|---|
| EOS/src/LocalThermodynamics.cpp:21,121 | Energy bracket and energy, momentum, mu, dmu/dn. Small-x analytic series avoids subtracting its leading x^3 term; direct branch subtracts x*sqrt(1+x²)*(1+2x²)-asinh(x). | Retain extended intermediates. Local thermodynamic/PE independent quadrature tests plus Decimal oracle below exercise tiny and relativistic densities. x86 has 64 mantissa bits; arm64 typically 53. Report actual widths, not assumed equivalence. |
| EOS/src/LocalThermodynamics.cpp:208 | Inverse density p=sqrt((mu-m)(mu+m)), n=(p/hbarc)^3/(3*pi²). | Candidate binary64 inverse plus explicit fma onset; exact-input 80/120-digit tests exercise both production entrypoints. Muon threshold and generated metadata may move; fresh tables required. |
| Analysis/src/ChemicalResponse.cpp:455 | Pe tail geometry denominator, mass/radius/lapse enclosures and absolute response bound. Subtractions can amplify roundoff when denominators approach zero. | Retain long double and fail-closed positivity/enclosure checks. chemical_production_contract checks valid bound and radius/old-mass/cut mutants; chemical_production_validation checks physical tail certificate and refinement. A width change can alter a marginal inequality: actual Mac checks required, no blanket tolerance exemption. |
| Physics/Rotochemical/src/GlobalUrcaChannelCoefficient.cpp:44 | Sum of nonnegative quadrature terms formed in binary64. No subtractive cancellation in this accumulator. | Retain long double. For N terms, accumulation error bounded by gamma_(N-1) times sum(abs(term)) with u set by actual mantissa width, plus term formation/quadrature errors and final cast. phase5d response analytic oracles cover polynomial radial integrals and refusal controls. Width simulation is diagnostic only; actual Mac execution pending. |

Additional held-out local-EOS oracle rules (frozen before scoring): exact binary
n,m and rounded hbar-c member; 120-digit Decimal cbrt/sqrt/asinh with independent
phase-space closed form. Momentum/mu/derivative bound 32 epsilon times magnitude.
Energy direct branch uses 32 epsilon * prefactor * (|a|+|b|), a and b the two
subtracted analytic terms. This is an absolute conditioning-aware bound, not
relative error at cancellation. Small-x series uses 32 epsilon * |energy| plus
prefactor * 7*x^13/416, the first omitted alternating-series term (x<0.01).
The safety factor covers input/intermediate rounding and transcendental errors;
these tests measure binary64-level accuracy, not a correctly-rounded guarantee.
Test densities: 1e-18, 1e-12, 1e-6, 0.03125, 0.5 for electrons and neutrons.

Mac contraction remains its previously recorded implicit AppleClang default.
Linux -ffp-contract=off is unchanged. Explicit std::fma is deliberate at one
expression. Standalone default/off AppleClang probes and actual production tests
are separate handoff requirements. No simulation can discharge that requirement.

Known Mac record (not a fresh execution): cmake/dependency-lock.json and
phase6-campaign-readiness/environment.json:38 record Apple clang 21.0.0
(clang-2100.3.34.2), target arm64-apple-darwin25.6.0, CMake 4.2.1. The prior
CompactStar compile configuration omits an explicit contraction flag; dependency
HistoricalMath translation units separately disable contraction/builtins/LTO.
This correction changes none of those Mac settings.
