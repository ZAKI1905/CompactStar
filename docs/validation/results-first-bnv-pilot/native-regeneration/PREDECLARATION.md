# Authorized native input regeneration (candidate)

Authority: owner request “Regenerate Missing Stellar Inputs and Complete Native P2
Assembly”, starting aa6239b963b7a48946f585af158482659d8fc59f. This supersedes the
previous prohibition on this single regeneration only. No thermal trajectories,
BA13 recalculation, certificate promotion, or model study is authorized.

Worktree: /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot;
branch candidate/results-first-bnv-pilot-20261108. Clean starting tree; all 15
local branch tips are ancestors of starting HEAD; nine retained worktrees.

The only producer invocation is EquilibriumSequenceNumberDerivative::Compute
(CompactStar/Analysis/src/ParticleNumberResponse.cpp:359). Its implementation,
TOV, Hartle, EOS/provider implementations, and historical builders remain unchanged.
Recipe: canonical Structure-1 freegas EOS8192, rho_ref=1.10e15 g/cm^3,
branch (1.095e15,1.105e15), radial80000, h=(.001,.0005,.00025),
k=(-2,-1,0,1,2), rho=rho_ref*exp(k*h). Fifteen contributors, eight distinct
noncentral configurations, three central contributors. PB13-comparison-v1 tails
are evaluated on each owning star by the existing callback; no borrowed tails.
Native EOS SHA256 f6b9c4068a430cc356b54e39c21a0ca17e17b8c2b35d17b836638298d1112d7f.
All retained input checksums are inherited verbatim from native-assembly/evidence/input-inventory.json.

The retained central export lacks full owner data. Reuse it as exact validation
and certificate authority, and use contributor ordinal2 from this sequence as the
complete runtime central owner. No additional central solve. Each returned
contributor carries its complete profile (including species fractions, nu-prime,
EOS d-epsilon/d-p), first-order and monopole Hartle dependencies. Persist all
four derivatives, errors, 3/5-point ladders, achieved energies, realized tails,
and source provenance. No aggregate-tangent substitution.

Allocation: phyall/phyn004, one task, four allocated CPUs for build (numerics
single-thread), 8 GiB, 15 minutes. Hard Slurm timeout, no automatic retries.
Persist capture immediately. A failed loader or downstream assembly must be
repaired against that capture; it does not authorize another producer call.

Before launch: SHA256-authenticate all native inputs, producer translation units,
provider source and dependency lock. Require committed clean source at launch.
Stop on any mismatch, producer refusal, nonfinite/incomplete owner, stale Hartle
provenance, central mismatch, derivative/tangent disagreement, missing certificate,
or coefficient/table mismatch. Do not change the model to obtain a pass.
Output checks: producer RequireCurrent; exactly15 owners and complete schemas;
lossless 17-digit binary64 text roundtrip and SHA256 identity; independent reload;
exact central retained-profile matching; existing chemical/fixed-B/tangent and
runtime guards; applicable static source/control and energy-ledger checks.

Change classes: engineering, structural/serialization, build, documentation,
generated evidence. Governing ADR-0003/0006/0007/0008/0011 ownership and Hartle,
ADR-0013/0014 chemistry, ADR-0015 open-system energy, ADR-0016 typed derivative
consumer, ADR-0017 checkpoint infrastructure. No numerical algorithm or physical
assumption changes. AI-authored transport remains a candidate, not ratification.
