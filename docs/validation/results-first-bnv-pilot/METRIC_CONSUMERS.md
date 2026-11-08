# P2 frozen-metric consumer audit

Source: `bb0f9c5af8ebaf86b5a59aece5551bb9dafbd264`. Code reading only; no new
stellar evaluation. This conclusion is specific to spin-off, frozen Structure-1,
modified-Urca-only P2 with the current thermal/envelope owners. It supports a
proposed scope decision, not a certificate acceptance.

**Conclusion:** an independent nodewise BA13 metric reference is unnecessary for
this restricted P2 observable if the owner accepts qualification through the
following consumer functionals and their joint error propagation. The literal
preflight §11 still requires it for full BA13. Metric arrays are actually read
inside the thermal path; claiming they are never read would be wrong.

| Path and source location | Actual metric use | Quantity carrying it into this pilot |
|---|---|---|
| `tests/bnv/breakpoint_fixture.hpp:9`, `campaign_fixture.hpp:12` | Assemble common ordinary context, numbers and sequence tangent | N_n/N_e/N_mu, t_n/t_e/t_mu, both Z rows; W is multiplied by exactly zero spin |
| `tests/bnv/controlled_neutron_sink_fixture.hpp:49` | Integrate exp(nu) n_n dV using r,m,nu; store mean event metadata and gamma normalization | Total count loss is prescribed exactly. P2 does not use mean radius/density/lapse to compute deposited power. Gamma is a normalization diagnostic, not a microscopic rate or a newly certified observable |
| `tests/bnv/energy_partition_fixtures.hpp:39`, `coupled_trajectory.cpp:130` | P2 returns zero escaping/product energy; its local-energy check is disabled, as in the existing fixture | mu_B and t supply actual mu_n; P2_average is the exact escaping-energy null. The conservative finite-T floor uses the existing 1 MeV scale, not a new metric functional |
| `CompactStar/Physics/BNV/src/DirectEnergyLedger.cpp:77` | Optional local/infinity consistency check; P2 does not activate it | Existing R18/energy checks retained; no nodewise comparison |
| `CompactStar/Physics/Rotochemical/src/GlobalUrcaChannelCoefficient.cpp:55` | Global integral with exp(lambda+(2-q)nu) | Ltilde_Me, Ltilde_Mmu; full neutrinos, signed rates and reaction heating use these same coefficients |
| `CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:130` | Calls actual photon/thermal owners at current T | Cstar(T), local surface temperature, redshifted surface temperature, radius and surface gravity |
| `CompactStar/Physics/Driver/Thermal/src/PhotonCooling_Details.cpp:170` | Envelope-base redshift and surface g; surface R and exp(2nu); heat capacity integral/cache | envelope(T), Tsurface_inf(T), radius, surface_gravity, Cstar(T), all explicitly sampled by `tests/bnv/frozen_sensitivity.cpp:284` using this same owner |
| `CompactStar/Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp:856` | Disabled ordinary DU/MU/PBF options return zero before their integrated emission path | No additional neutrino metric functional; controlled beta owner supplies MU |
| `CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:45` | Metric identity string, owner/currentness checks, then validity before every RHS | Provenance and refusal, not a nodewise accuracy estimator |
| `CompactStar/Physics/BNV/src/FrozenBnvValidity.cpp:10` | Requires all 22 scalar group names, interpolates their historic utilizations | Keep metric_structure row and every stop; no nodewise array consumed by monitor |

Photon closure is explicit: with fixed emission options,
Lgamma_inf = 4*pi*sigma*R^2*Tsurface_local^2*Tsurface_inf^2, because
Tsurface_inf = exp(nu_surface)*Tsurface_local. Thus no extra surface-lapse
observable is needed. Bounds must propagate the radius and both temperature
errors jointly; each row's 1e-4 limit is not automatically a 1e-4 luminosity bound.
Likewise Cstar/Ltilde errors are model/frozen uncertainties, separate from ODE
accuracy. Thermal knots, between-knot interpolation and its absolute cache bias
retain their existing limitations (review §5; readiness predeclaration:1).

P1 has an additional integrated sea-energy owner at
`tests/bnv/energy_partition_fixtures.hpp:64`; it is not instantiated here. No claim
about a new spatial source, P1 physical decay, evolving metric, superfluidity or
additional radiation channel follows. A future consumer outside this inventory
must reopen the scope decision, not inherit it automatically.

## Native binding is a separate prerequisite

`FrozenRotochemicalRunContext.hpp:85` hardcodes four Mac hashes for profile,
model, chemical certificate and EOS, even for AnalyticControl. The fresh read-only
audit verifies all four differ from the native target-0 artifacts. The subsequent
check at line 91 refuses them. `tests/bnv/ba13_native_certificate_binding.cpp:1`
explicitly constructs no thermal context, RHS or trajectory; its retained PASS
cannot discharge this gate. `tests/rotochemical/fixture.hpp:60` additionally
requires entry-for-entry equality with the freshly constructed owning star.

`tests/bnv/breakpoint_fixture.hpp:28` requires the historical absolute Bdot
`-2.4136520263641375e37` count/s. Native B0 and the unchanged fractional drive
instead give `-2.413652026364133e37`. This is not permission to round or weaken
equality. The separately named native pilot must bind its own exact derived Bdot
and unchanged fractional drive; the original run card and literal remain intact.

Minimal subsequent implementation: an explicit, owner-authorized native identity
record for the new pilot, exact hashes and semantic ownership, distinct from the
historical Mac identity; source and control share it. Missing/foreign/mixed/stale
identities must refuse. No arbitrary caller-supplied hash may self-authorize a
fixture. This audit changes none of these guards. Evidence: `evidence/preflight.json`.
