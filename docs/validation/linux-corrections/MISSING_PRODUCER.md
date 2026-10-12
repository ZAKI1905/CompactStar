# One specific missing-source/input request

Please supply the original `tests/bnv/produce_frozen_certificates.py`, recorded
SHA-256 `fb609dbd549763800cc47e962e4f12333835916ea299058388dcb6c1c592ad79`,
with its invoked stellar/coefficient sample producer source, exact build command,
and the retained 21 raw sample inputs/outputs or the complete regeneration recipe.
Include EOS/profile construction, achieved-B grid, derivative/coefficient authorities,
numerical settings, certificate tests and uncertainty calculations.

Evidence: PHASE6A1_CONTROLLED_BNV_RECOVERY_IMPLEMENTATION.md:109 records the
script hash; :152 records the retained 21-star refit and max utilization
0.54075347894390757 (N_mu, index 20). No tracked script of that name exists in
any supplied ref's history (`git log --all -- tests/bnv/produce_frozen_certificates.py`).
The CompactStar bundle contains the four already-present branch tips d8e7208,
9a5c1ec, 5433fc2 and 0c3f9269. Both incoming ZIP inventories lack the producer.
The six transferred Phase-6 input artifacts contain final frozen.tsv and
coefficients.tsv, not their complete producer/source/sample authority.
`mac-recovery.json` records a build-directory location only, not producer identity.

Native Phase-6 certificate construction is BLOCKED on this specific input;
no guessed producer, Mac sensitivity transfer, exact-profile guard bypass or
new 21-star baseline was attempted. Native profile diagnostics and static
coefficient refinement remain independently executable. Certificate-dependent
source/control/event pilots remain BLOCKED; their code and prior failure evidence
are preserved. Resume from these corrections when the source is supplied.
