# ADR-0023: Retained native paired-context reconstruction

Status: **PROPOSED; no alternative adopted or guard changed.**

The owner requested native RF-BNV-P2 assembly from retained Linux artifacts,
without stellar solves, thermal trajectories, or weakened guards. The starting
revision is `afeb0418ab697829f3fec5722022f70db65f4b5b`. This record documents a
missing-input stop under GOVERNANCE §3 and AGENTS.md §8, not a new BA13 decision.

The native data are authenticated but do not contain a serialized Phase-5B
derivative owner. The target-0 sequence directories have no files. The native
preparation fixture likewise retained no central or neighboring-star files.
See `../validation/results-first-bnv-pilot/native-assembly/evidence/input-inventory.json:1`
and `../validation/results-first-bnv-pilot/native-assembly/REPORT.md:1`.

`tests/bnv/frozen_sensitivity.cpp:338` writes B0, closed t, propagated t errors,
raw closure and closure budget. It does not serialize the four canonical
derivatives, their errors, canonical baryon denominator, three/five-point
ladder, achieved central states, realized tails, or source owners.
`Analysis::EquilibriumSequenceNumberDerivative::Compute` retains those in memory
(`CompactStar/Analysis/src/ParticleNumberResponse.cpp:359`).

ADR-0016 §§1.1,3–5 requires the adapter to remain a typed Phase-5B consumer,
with its existing errors, denominator, provenance and currentness semantics.
The aggregate tangent does not determine the missing derivative denominator:
scaling all derivatives by the same positive factor leaves t unchanged.
It therefore cannot reconstruct the existing object or independently satisfy
its conditioning check. Inventing a `NumberResult` with matching ratios is not
deserialization (`CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp:41`).

Alternatives requiring an evidence-based decision:

1. Locate a complete retained native owner archive, authenticate its lineage,
   and implement a lossless loader preserving all current checks. This follows
   the requested scope if such an archive exists. No new solve is needed.
2. Propose an immutable artifact-backed authority with explicitly justified
   equivalence to every displaced live-owner guard. This changes the ownership
   boundary and needs a separate decision and evidence; acceptance of the BA13
   sidecar alone does not supply that evidence. No such adapter is implemented
   or approved here.

A new stellar solve, a foreign coefficient table, a fabricated denominator,
or promotion of the conditional sidecar to the historical runtime certificate
is outside this task. Nothing in this proposal authorizes those actions.

The historical builder and scientific/numerical sources remain byte-identical.
The next work remains paired assembly and the controlled heating calculation;
no new static qualification framework is proposed.
