# CompactStar External-Dependency Migration Acceptance

**Status:** Owner accepted

**Authority:** Qualified local-Mac dependency stack only

**Accepted CompactStar production candidate:** `c22e300ac802ab73602dba98ce8ff89c903a1843`

**Accepted dependency stack:**

- ZakiLib 2.0.1 production candidate: `e263a6e180c5c417198e7778bd21fc9c0a32dc33`
- CONFIND 2.0: `b0cbd510fd3fd0c772fa50499cd749287cb39e7b`

Canonical Zaki master now contains the accepted production candidate plus its acceptance documentation only. Canonical CONFIND remains unchanged.

The owner accepts ADR-0019 and the package-based dependency architecture. The historical vendored Zaki and CONFIND archives remain preserved scientific oracles and are inactive build inputs. The plotting, compiled-Python, and obsolete OpenMP dependency removal is accepted.

The accepted qualification established:

- TaskManager T1 exact in Debug and Release;
- all 38 deterministic TaskManager T2 outputs exact in Debug and Release;
- final neutron, Lambda, and Sigma- BNV lifetime outputs exact;
- Phase-5B, Phase-5C, Phase-5D, Phase-5D1, and Phase-5D coupled-oracle gates passed; and
- the ADR-0017 scientific result remained exact.

The combined independent review of ZakiLib 2.0.1 and the CompactStar migration reported:

- `BLOCKING = 0`
- `MATERIAL = 0`

Accepted nonblocking caveats remain deferred. This acceptance does not establish Linux, cluster, or multithreaded TaskManager authority.

This acceptance commit changes documentation only. CompactStar production, tests, CMake, baselines, and ADR-0019 remain byte-identical to the accepted production candidate above.
