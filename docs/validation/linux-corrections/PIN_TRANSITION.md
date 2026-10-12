# Proposed upstream frozen-context extension

Authority: frozen-context provenance plan §§5,8; ADR-0020 proposed amendment;
owner 2026-10-11 correction instruction. This is not a historical plan edit.
`protected-successor.json` enumerates all 33 protected paths, all ten historical
baselines and three special artifacts. `transitions` records old pins, exact new
bytes, origin, behavior analysis and acceptance limitation for each changed path.
The four already-stale paths and two newly authorized numerical-method paths and the generated-provenance path are
separately explained. No paths removed; no artifact or baseline hash refreshed.

The compiler-selection edit is retained: an absolute CXX records the actual Linux
compiler while keeping /usr/bin/clang++ on Darwin. Moving it outside the protected
producer would obscure which compiler it records without improving preservation.

ADR-0017 Amendment 1 accepted bounded numerical architecture at 0c3f9269;
PHASE6_BREAKPOINT_RATIFICATION.md explicitly declines blanket ancestry ratification.
Thus successful file contents and scientific formula preservation are supported,
but no claim of prior acceptance of these individual upstream source pins is made.
The current owner instruction authorizes the candidate correction, not ratification.

Run tests/linux/protected_successor.py for historical and proposed results,
including topology, forged-record and real-file mutation controls. Historical
phase5d_protected_manifest and coupled_oracles stay unchanged and must still fail.
The successor checks are provenance qualification, not a pass of coupled physics.

The final metadata audit also identified tests/eos/structure1/table.hpp as a
protected path containing hardcoded historical commits labeled as current source.
Its seventh transition explicitly renames them as historical origins and records
the proposed arithmetic contract. No table-grid or equation change is hidden in
that metadata fix. The pre-metadata six-mismatch receipt is retained alongside
the final seven-mismatch receipt.
