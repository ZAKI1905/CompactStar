# ADR-0024: Native stellar owner transport

Status: candidate implementation under explicit owner authorization to regenerate
missing inputs and complete native assembly; no historical scientific promotion.

ADR-0023 documented missing retained owner data. The subsequent owner request
explicitly authorizes one unchanged fifteen-star sequence. Transport the returned
Phase-5B objects losslessly, preserving full owner lifetime and source-currentness
semantics (ADR-0016), rather than replacing them with tangent summaries.

A narrow friend transport seam in NStar, StarProfile, RotationSolver and
NumberResult permits tests/native/stellar_archive.hpp to serialize the existing
results. No numerical producer implementation changes. The format stores radial
and EOS-derivative columns, species layout, scalar/version state, both published
Hartle responses, sequence results and full metadata. Solver scratch workspaces
are not runtime dependencies and are not restored. Load never solves or seals stale
provenance: it authenticates a required SHA256, validates shapes/versions, relocates
addresses to newly owned objects, checks original source bytes and invokes existing
RequireCurrent. Owners and EOS snapshots remain alive through shared ownership.
The historical builder remains available unchanged.

Alternative of storing only tangent ratios is rejected by ADR-0016: it loses
canonical derivatives, denominator and errors. Re-solving during load is excluded
by the restricted authorization. See native-regeneration/PREDECLARATION.md for
recipe, validation and stopping rules; tests/native/regenerate.cpp for the sole
producer call. Structural and engineering change; candidate until owner review.
