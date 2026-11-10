# Accepted Stage A restart-selection amendment

The owner explicitly approved the narrowly scoped amendment proposed at
`3557da9861274033f642ca08f690b4fb59105797` in the current task:

> Replace the requirement to select the first accepted endpoint at or beyond
> half-time with the latest accepted noninitial endpoint at or before half-time,
> provided it leaves a valid, nonzero restart tail.

Authorization includes applying the exact patch, selector and binary64 round-trip
tests, and authentication/necessary republication of runtime snapshots. It does
not authorize another Stage A trajectory. Integration algorithms, controls,
tolerances, physical inputs and numerical acceptance criteria remain unchanged.
No Stage B, stellar solves, BA13 reopening or master merge is authorized.

Authority: the owner's explicit instruction; ADR-0017 §§4–6 and Amendment 1
(`docs/adr/ADR-0017-phase6-passive-checkpoint-reconstruction.md:108`;
`docs/validation/PHASE6A1_PIECEWISE_RECONSTRUCTION_PROPOSAL.md:16`). Class:
numerical-method/test-protocol, engineering, documentation and generated evidence
under GOVERNANCE §2. No physics/solver method or ownership change occurs.

The implementation is exactly the four-file patch retained at
`../restart-coverage/proposed-selector.patch`. That proposal's dated pending
status remains historical; this acceptance supersedes only its authorization
status. Its tests do not qualify trajectory accuracy. Previous D1 and Stage A
failure evidence, original snapshot identities and their raw outputs are retained.

Preparation checks are bounded, one logical CPU, 4 GiB, 30 minutes on EKU. They
build the harness, run contracts and standalone selector/serialization tests,
and evaluate the twelve existing static RHS states from the retained archive.
They execute no stellar or thermal integration. A future source-only retry must
be separately authorized against the new committed preparation receipt.
