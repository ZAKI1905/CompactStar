# Minimal physical follow-on (recommendation only)

Recommend **n -> chi + gamma**, restricted to a parameter region where chi
escapes and the photon thermalizes, after these assumptions are checked. It reuses
the neutron-source projection and energy ledger. Genuine total BNV depends on
chi's quantum numbers; ordinary neutron disappearance alone does not establish it.
The microscopic starting point is Berryman, Gardner & Zakeri, PRD109,023021,
[arXiv:2305.13377v2, §5 and Appendix G](https://arxiv.org/html/2305.13377v2#S5),
which treats in-medium radiative dark baryon decays. This is a recommendation,
not an accepted rate or model.

Required input is the differential proper in-medium rate over neutron momentum
and final-state phase space, with m_chi, coupling/mixing, dispersion relations,
occupations/blocking and photon medium effects consistent with the selected EOS.
No vacuum lifetime or uniform sink normalization substitutes for that input.

Under the stated fate assumptions, escaping energy is E_chi_inf; the direct
ordinary-fluid contribution per event is mu_n_actual_inf minus E_chi_inf. Retained
photon energy is already inside this ledger; do not add another hole/photon heat.
Check exp(Phi)*E_chi_local >= m_chi*c^2 and escape opacity, and quantify any bound
fraction. Accumulation invalidating Regime I requires additional physics, not an
unrecorded retention efficiency (ADR-0015 §§9,11–13).

Minimal implementation: a process-specific source/partition/fate adapter sharing
one event identity, differential-rate integration, and escape/thermalization
validation. Reuse the controlled thermal and chemical owners. The legacy
`CompactStar/Microphysics/BNV/Channels/BNV_B_Chi_Photon.hpp` is a candidate reference
to audit, not already-qualified thermal integration. First observational target:
redshifted surface-temperature/luminosity excess versus age and coupling; a late
temperature floor only if the evolution actually reaches it. No numeric prediction
or mass/coupling is chosen here.

The Mac source library catalog and PDFs are unavailable on this Linux host; the
linked primary paper's HTML was consulted. No claim is made to have read a local
PDF or to have validated its rate implementation.
