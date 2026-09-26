# Phase-6 EKU cluster numerical-platform qualification preflight

<!-- Intended repository path: docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md -->

## 0. Status

**Status:** PROPOSED — OWNER ACCEPTANCE REQUIRED BEFORE CLUSTER MUTATION/JOB SUBMISSION.

**Change class:** documentation (planning). It cites `GOVERNANCE.md:36-57` for change classes
and `AGENTS.md` §8 for fail-closed conditions.

**Disposition proposed by this preflight:** see §24.

This preflight authorizes nothing by itself. It makes no change to production source, tests,
CMake, baselines, EOS/data, or literature. It runs no ODE, no BA12R, no six-trajectory campaign,
no Slurm job, no remote build, and no remote package install.

| Authority | Identity |
|---|---|
| `PHASE6_CLUSTER_PREFLIGHT_ENTRY_SHA` | `232565a32303a4953f3f516d1d5286b6663f8f99` |
| live `refs/heads/master` (read with `git ls-remote https://github.com/ZAKI1905/CompactStar.git`, 2026-09-25) | `232565a32303a4953f3f516d1d5286b6663f8f99` |
| ADR-0015 / ADR-0016 / ADR-0017 | ACCEPTED / ACCEPTED / ACCEPTED (human-ratified) |
| ADR-0017 production implementation | QUALIFIED / OWNER-ACCEPTED / CANONICALLY INTEGRATED (`docs/validation/PHASE6_ADR0017_PRODUCTION_ACCEPTANCE.md:5-7,132-148`) |
| historical BA12 / BA12R | FAIL / FAIL, permanently |
| clean six-run BA12R campaign | NOT AUTHORIZED |
| EKU cluster numerical platform | NOT QUALIFIED |

### 0.1 Execution-context note, recorded as a finding

This preflight was prepared by an agent session running **on the EKU login node `ekucn1.eku.edu`**
(user `zaki`, uid 1038), not on the Mac. The canonical Mac checkout
`/Users/keeper/Documents/CompactStar/repo/CompactStar` does not exist on this host. Cloning,
branching, or creating a worktree on the cluster is forbidden by the task. Therefore:

- canonical entry was authenticated **against live GitHub** (`git ls-remote` over HTTPS). The
  governing documents were read by streaming single files from
  `raw.githubusercontent.com/ZAKI1905/CompactStar/232565a…`. No clone was made, and no
  source tree was written to the cluster;
- the local branch `docs/phase6-cluster-qualification-preflight` and its worktree were **not
  created** from this session. This document was written to agent session scratch only;
- the Mac-side placement, `git diff --check`, commit, and push must be done by the owner, or by a
  Mac-hosted agent, exactly as in §23.

Before this preflight is committed, the canonical Mac checkout's clean state and its
HEAD/local/origin equality must be verified on the Mac (§23). Live `master` was equal to the
entry SHA when this document was prepared.

---

## 1. Authenticated read-only cluster discovery (2026-09-25)

Every fact below was obtained with inspection commands only. Compute-node facts that could
not be read live without `srun` (which this task forbids) are marked **SECONDARY**. Those come
from existing PWNSpall evidence files on `/mnt/sdd`, cited by path. They must be re-authenticated
live in CQ0.

### 1.1 Login node

| Item | Value | Source |
|---|---|---|
| hostname | `ekucn1.eku.edu` | `hostname -f` |
| OS / kernel | Rocky Linux 9.7 (Blue Onyx) / `5.14.0-570.42.2.el9_6.x86_64` | `/etc/os-release`, `uname -a` |
| CPU | AMD EPYC 7232P 8-Core (Zen 2, family 23 model 49), 1 socket, 8 cores, 16 threads, 1 NUMA node | `lscpu` |
| RAM / swap | 15 GiB (356 MiB free, 6.4 GiB available) / 7.7 GiB (6.2 GiB used) | `free -h` |
| root filesystem | 70 G xfs, **95 % used** | `df -hT` |
| glibc | 2.34; `libm.so.6` SHA-256 `d4d4aca358a3704f3e7e7bf3cc782dde7e629af2a7e4adc197d74782af9e0844` | `getconf`, `sha256sum` |
| role | Slurm controller (`SlurmctldHost=ekucn1`) **and** NFS server for `/home`, `/mnt/sdd`, `/mnt/hdd`, `/opt` | `scontrol show config`, `/etc/exports` |

**Consequence:** the login node is small, already swapping, and uses a *different CPU model* from
the compute nodes. No CompactStar build, test, or scientific executable may run on it. Only light
Git and inspection operations may.

### 1.2 Compute nodes and Slurm

| Item | Value | Source |
|---|---|---|
| Slurm | `slurm 22.05.9`; `ClusterName=warewulf` | `sinfo --version`, `scontrol show config` |
| nodes | `phyn001`–`phyn004` | `sinfo -N -l` |
| topology per node | 2 sockets × 24 cores × 2 threads = 96 logical CPUs (48 physical cores) | `scontrol show node`, `sinfo -N -l` (S:C:T = 2:24:2) |
| RealMemory per node | `125000` MB | `scontrol show node` |
| TmpDisk | `0` (Warewulf-provisioned, likely stateless; node-local `/tmp` capacity UNKNOWN) | `sinfo -N -l` |
| compute OS kernel | `5.14.0-427.20.1.el9_4.0.1.x86_64` (Rocky 9.4 image; differs from login 9.7) | `scontrol show node` |
| CPU model | AMD EPYC 7352 24-Core (Zen 2), recorded for phyn002/003/004; **SECONDARY** | `/mnt/sdd/zaki/PWNSpall/qualification/cluster-runtime-04/prelaunch/node-probe-prelaunch-51f9991aa678-phyn00{2,3,4}.json` |
| phyn001 CPU model | not recorded anywhere; same S:C:T/memory. **UNKNOWN — CQ0 must authenticate** | — |
| NUMA topology on compute | **UNKNOWN** (needs `srun lscpu`); expected 2+ NUMA nodes | — |
| `/mnt/sdd` on compute | `192.168.13.1:/mnt/sdd` mounted `nfs4` at `/mnt/sdd` on phyn001 (2026-08-21); **SECONDARY** | `/mnt/sdd/zaki/PWNSpall/qualification/cluster-scientific-pilot-01/mount-identities/phyn001.json` |
| `/opt/spack` on compute | exported `ro`. PWNSpall compute-node Slurm build jobs used `/opt/spack/.../gcc-14.1.0` and `cmake-3.29.4`; **SECONDARY** | `/etc/exports`; `/mnt/sdd/zaki/PWNSpall/qualification/cluster-runtime-01/slurm/root-build-2407.out` |

| Partition | Nodes | State (2026-09-25 09:45) | Default | MaxTime / DefaultTime | Memory policy | OverSubscribe | Preemption |
|---|---|---|---|---|---|---|---|
| `phy01` | phyn[001-002] | both idle | no | UNLIMITED / NONE | DefMemPerNode=MaxMemPerNode=UNLIMITED | NO | OFF |
| `phy02` | phyn[003-004] | phyn003 **drained** ("Kill task failed", root, 2026-09-18T15:38:05); phyn004 idle | no | UNLIMITED / NONE | same | NO | OFF |
| `phyall` | phyn[001-004] | 3 idle, 1 drained | **YES** | UNLIMITED / NONE | same | NO | OFF |

Scheduler policy (`scontrol show config`, `/etc/slurm/cgroup.conf`):

- `SelectType=select/cons_tres`, `SelectTypeParameters=CR_CORE`. Cores are the consumable
  resource and memory is **not**, so a `--mem` request is used only for node selection.
  Allocation granularity is one physical core, which is two logical CPUs.
- `TaskPlugin=task/affinity,task/cgroup`, `ProctrackType=proctrack/cgroup`. However,
  `cgroup.conf` has **`ConstrainCores=no` and `ConstrainRAMSpace=no`**. Slurm does not confine a
  job to its cores or its memory. Process isolation is therefore procedural: explicit
  `OMP_NUM_THREADS=1`, `--cpu-bind=cores`, and headroom-based packing.
- **Accounting is disabled** (`AccountingStorageType=accounting_storage/none`,
  `JobAcctGatherType=jobacct_gather/none`; `sacct` reports that accounting storage is disabled).
  `sacct` cannot provide provenance, so all resource usage must be captured inside each job (§17).
- `MinJobAge=300 s`: finished job records disappear from `scontrol` after 5 min.
- `MaxArraySize=1001`, `MaxJobCount=10000`, `ReturnToService=1`, `KillWait=30 s`,
  `PreemptType=preempt/none`, `sched/backfill`.

### 1.3 Filesystems

| Mount | Type / backing | Size / free | Exported to compute? | Suitability |
|---|---|---|---|---|
| `/mnt/sdd` | xfs on local LV `ssdvg-ssd_lv` of `ekucn1` | 894 G / 591 G free | yes, `rw,sync` to 192.168.13.0/24 | **persistent project root (selected)** |
| `/home` | xfs local LV | 874 G / 830 G free | yes, `rw,sync` | usable, but the owner's precedent places projects on `/mnt/sdd` |
| `/mnt/hdd` | xfs local LV | 9.1 T, **95 % used** | yes | rejected (nearly full) |
| `/mnt/synology` | NFS v3 from 192.168.13.56 | 104 T / 30 T free | **no** (the export line is commented out) | login-only; could serve as an owner-directed evidence archive, not for jobs |
| `/opt` | login root filesystem | — | yes, `ro` | Spack toolchains |
| node-local `/tmp` | UNKNOWN (TmpDisk=0) | UNKNOWN | — | CQ0 measures it; not relied on |

**Risk recorded:** `/mnt/sdd` is a single local volume served by the login node. A login-node
reboot or NFS outage makes it unavailable to every running job, and no backup policy was
found. Durable evidence must therefore also be copied back to the Mac and committed to Git (§16).

### 1.4 Software discovery

Environment Modules 5.3.0. The Spack module tree is `/opt/spack/share/spack/modules/linux-rocky9-zen2`
(282 modules). No module is loaded by default.

| Dependency | Mac reference | Cluster availability (exact) | Classification |
|---|---|---|---|
| C++ compiler | Apple clang 21.0.0 (`clang-2100.3.34.2`) | system `/usr/bin/g++` GCC 11.5.0 (login; compute image version UNKNOWN). Spack `gcc/14.1.0-gcc-11.4.1-yd5t5aj` (prefix `/opt/spack/linux-rocky9-zen2/gcc-11.4.1/gcc-14.1.0-yd5t5aj4bggpy5swlajrqvsiohzyyfbh`). Spack `llvm/12.0.1`, `llvm/17.0.6` | COMPATIBLE_DIFFERENT_VERSION |
| Clang/LLVM | AppleClang 21 | Spack LLVM 12.0.1, 17.0.6. System: only `llvm-libs-20.1.8` (no clang driver) | COMPATIBLE_DIFFERENT_VERSION (not selected) |
| CMake | 4.2.1 | Spack `cmake/3.29.4-gcc-11.4.1-gybrkgt` only. No system cmake | COMPATIBLE_DIFFERENT_VERSION (top level requires ≥3.10) |
| GSL | **2.7.1** | Spack `gsl/2.6-gcc-11.4.1-ikklbhc` only (`libgsl.so.25`, static `libgsl.a` present). No system GSL | **MISSING (2.7.1)**; 2.6 is COMPATIBLE_DIFFERENT_VERSION, but equivalence of its ODE sources to 2.7.1 is UNKNOWN (a live source diff attempt timed out) |
| Ninja | not recorded | Spack `ninja/1.10.2`, `ninja/1.12.0` | available |
| OpenMP runtime | AppleClang `-Xclang -fopenmp` (libomp) | GCC `libgomp` (system `libgomp-11.5.0`; GCC 14.1.0 bundled libgomp) | COMPATIBLE_DIFFERENT_VERSION |
| Python 3 interpreter + dev headers | project `.venv` Python (version not recorded in Phase-6 evidence) | system 3.9.25 (no `python3-devel`), Spack `python/3.8.7`, `python/3.11.9` (has `Python.h`) | COMPATIBLE_DIFFERENT_VERSION |
| **NumPy** (`find_package(Python3 … NumPy REQUIRED)`, `CMakeLists.txt`) | present on Mac | **absent** from system, Spack 3.11.9, and PWNSpall's `python-3.11.9` (PWNSpall's own ROOT build logged "Could NOT find Python3 (missing: … NumPy)") | **MISSING** |
| **`libZaki.a` / `libConfind.a`** (vendored, hard-set paths) | `dependencies/lib/{Zaki,Confind}/Darwin/arm64/` | Only `Darwin/arm64` and `Darwin/x86_64` exist in the tree. **No Linux build anywhere** on `/home/zaki` or `/mnt/sdd` | **MISSING — hard configure blocker** (§2) |
| Git | — | 2.47.3 | available |
| Slurm | — | 22.05.9 | available |
| `/usr/bin/time` (for RSS capture) | — | present | available |

---

## 2. Cluster blockers found by discovery (headline)

**B1 — The canonical SHA cannot configure on Linux.** `CMakeLists.txt` sets, unconditionally,

```cmake
set(ZAKI_LIB    ${DEP_LIB_DIR}/Zaki/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}/libZaki.a)
set(CONFIND_LIB ${DEP_LIB_DIR}/Confind/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}/libConfind.a)
```

and stops with `FATAL_ERROR "Missing Zaki static library"` when the file is absent. On the cluster
the path resolves to `Linux/x86_64`, which does not exist. Because these are plain `set()` calls,
not `CACHE` variables, `-DZAKI_LIB=…` cannot override them. `docs/build/MACOS_BUILD.md:10` states
that "No claim of Linux or Windows support is made or implied."

**B2 — Zaki source availability.** `ZAKI1905/Confind` is public (`master` =
`89c5d9b731534e4289d9f686549d9f0ac178e567`). `https://github.com/ZAKI1905/Zaki.git` returned
an authentication prompt, meaning it is private or has another name. **The owner must identify the
exact Zaki source repository and SHA** that produced the Darwin `libZaki.a`. Without it, a Linux
`libZaki.a` cannot be tied to the same source.

**B3 — NumPy is required at configure time and is absent** on every cluster Python. The Python
tests additionally need `mpmath` and `scipy` (§13.1). The C++ scientific paths load libpython
only as a link-time dependency; matplotlibcpp plotting is gated off (`gen_plots=false`) on every
registered path.

**B4 — GSL 2.7.1 is absent** (only Spack 2.6).

**B5 — The Phase-5 governed comparators and producers assume the same platform.** Every double
is compared bit-exactly. `architecture`, `platform`, `GSL` version, and hashes of regenerated
intermediates are equality-bearing. The 5B producer invokes the hard-coded `/usr/bin/clang++` to
record its compiler string. The 5D producer hard-codes SHA-256 values of regenerated
`freegas.tsv`/`profile.tsv`/`model.txt`/certificate (§12). Both portability ratifications cover
**compiler version on one Apple arm64 host only** and forbid generalizing to OS, architecture, or
GSL (`docs/validation/PHASE5C2_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md:76-78`,
`docs/validation/PHASE5B_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md:44-59`). The governed
comparators therefore **cannot and must not** be used as cross-platform comparators, and must
not be weakened.

**B6 — The ADR-0017 qualification verifier compares bit-exactly** against historical Mac oracle
bytes (`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:249-258`). It also hard-requires
macOS 26.6.2 / Darwin 25.6.0 / arm64 / Apple clang 21.0.0 / GSL 2.7.1 / CMake 4.2.1
(`tests/bnv/adr0017_production_verify.py:131-143`). A separate cluster-mode comparator
implementing §11 is needed, and the accepted verifier stays unchanged.

### 2.1 Governance consequence (PROPOSED-ADR alternatives; not decided here)

Building on Linux is a **dependency/build change that alters platform support**. Under
`GOVERNANCE.md:44-57` that requires recorded versions, reproducible build instructions, and **an
ADR**. Per `AGENTS.md` §8, the alternatives are recorded without selecting one on the owner's
behalf:

| Option | Mechanism | Consequence |
|---|---|---|
| L-A | Commit Linux `libZaki.a`/`libConfind.a` under `dependencies/lib/{Zaki,Confind}/Linux/x86_64/` | No CMake change. Adds opaque generated binaries tied to one cluster toolchain. They must be built from pinned Zaki/Confind SHAs with the cluster compiler. This is a generated-artifact class change needing its own ADR. The canonical SHA changes |
| L-B | CMake change: make the Zaki/Confind library paths `CACHE FILEPATH` overridable, and optionally make `NumPy`/`Development` conditional. Supply out-of-tree, versioned, hashed Linux builds | Keeps binaries out of Git. It is a small build-system change needing ADR plus regression that the Mac build is byte-identical in effective flags. The canonical SHA changes |
| L-C | Build Zaki/Confind from pinned source inside the CompactStar build (ExternalProject or vendored source) | Most reproducible, but has the largest build-system change |

The agent's **recommendation is L-B**, because it has the smallest tracked change and keeps
toolchain-specific binaries outside Git. The Mac path stays exactly as it is. The decision is the
owner's. Whichever option is chosen, cluster qualification begins from the **later
owner-approved canonical SHA** that contains the change (§10 permits this). It cannot begin from
`232565a…`.

**NumPy:** whether NumPy is needed only for configure/matplotlibcpp or also by the scientific
path or tests is audited in §8/§13. In either case, the bootstrap installs a pinned NumPy into a
user-local, hashed virtual environment (§4.2), not into any system location.

---

## 3. Mac reference platform (recorded, not re-measured)

From `docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:24-45`:

| Item | Value |
|---|---|
| OS | macOS 26.6.2 build 25G83; Darwin 25.6.0 |
| architecture | arm64 |
| compiler | Apple clang 21.0.0 (`clang-2100.3.34.2`), target `arm64-apple-darwin25.6.0` |
| CMake | 4.2.1 |
| GSL | 2.7.1 |
| build | Debug; effective `-g -std=c++17 -arch arm64 -fPIC -pthread -Wall -Wextra -Xclang -fopenmp` |
| FP flags | none unsafe (no `-ffast-math`, finite-math-only, unsafe-math, explicit FMA, or contraction-changing flag) |

The Mac is the development and numerical **reference** platform. The purpose of cluster
qualification is controlled **cross-platform numerical equivalence**, not identical binaries.

Known unavoidable sources of cross-platform difference:

1. Apple libm vs glibc 2.34 libm (`exp`, `log`, `pow`, and others are not correctly rounded, and
   their last bits differ).
2. AppleClang arm64 defaults to `-ffp-contract=on` (FMA contraction is possible even in Debug).
   GCC in ISO `-std=c++17` mode defaults to `-ffp-contract=off`, and baseline x86-64 has no FMA.
3. Possibly GSL 2.6 vs 2.7.1, which §4 removes.

None of these is a scientific change. Together they may change adaptive RKF45 step placement.

---

## 4. Toolchain recommendation

### 4.1 Primary qualification toolchain (one)

| Component | Selection | Reason |
|---|---|---|
| compiler | **Spack GCC 14.1.0** (`module load gcc/14.1.0-gcc-11.4.1-yd5t5aj`) | Modern C++17, bundled libgomp, and one set of NFS-shared `/opt` bytes identical on every node. Already proven in compute-node Slurm builds by PWNSpall. System GCC was rejected because the compute image (el9.4) and login image (el9.7) differ, and whether the compute image even has a compiler is UNKNOWN |
| C++ runtime | link with `-static-libstdc++ -static-libgcc` | Compute-node system libstdc++ is older than GCC 14's. Static linking puts the runtime inside the hashed executable. PWNSpall hit the related problem that Spack's GCC specs inject `-rpath`; the static runtime avoids it |
| CMake | Spack `cmake/3.29.4-gcc-11.4.1-gybrkgt` | only CMake available |
| generator | Spack `ninja/1.12.0-gcc-11.4.1-2stao5p` | deterministic parallel build |
| GSL | **user-local GSL 2.7.1** (option B below), linked **statically** by pointing FindGSL cache variables `GSL_LIBRARY`/`GSL_CBLAS_LIBRARY` at `libgsl.a`/`libgslcblas.a` | exact Mac version |
| Python | user-local venv from Spack `python/3.11.9`, with pinned `numpy`, `mpmath`, and `scipy` wheels (versions and wheel SHA-256 recorded). The 18 Python tests import all three (§13.1) | satisfies B3 |
| Zaki / Confind | built from owner-pinned source SHAs with the same GCC 14.1.0 and flags (§5), under a versioned prefix | satisfies B1/B2 per the owner's §2.1 decision |

### 4.2 GSL decision: option B (user-local GSL 2.7.1) recommended

Option A (qualify the Spack GSL 2.6 stack) is technically acceptable because cross-platform
equivalence is tested by tolerance anyway. **Option B is recommended** for three reasons:

1. A dependency-bootstrap task is already unavoidable (B1–B3). Adding GSL 2.7.1 costs almost
   nothing extra.
2. It removes a variable that is otherwise UNKNOWN: whether `ode-initval2`, `interpolation`,
   `roots`, and `integration` changed between 2.6 and 2.7.1.
3. It follows the minimum-difference principle.

The GSL 2.7.1 tarball SHA-256 must be recorded against the GNU release signature. Install under
`/mnt/sdd/zaki/CompactStar/toolchains/<toolchain-id>/install/gsl-2.7.1`, built with the same
GCC 14.1.0, Debug-neutral release flags `-O2 -g` with `-ffp-contract=off`, and **no**
`-march`. (GSL's own optimization level is recorded as a toolchain key. It is not CompactStar's
build type.) **Nothing is built or installed in this task.**

### 4.3 Toolchain bootstrap is a separate owner-authorized task (CQ-B)

The bootstrap covers the GSL 2.7.1, Zaki, Confind, and NumPy venv builds. It runs as an
owner-authorized task *before* CQ0, on a compute node, under `toolchains/`, with manifests
(source SHA/tarball hash, configure lines, `sha256` of every installed library).

---

## 5. Floating-point / compiler policy

**First cluster qualification build type: `Debug`**, the same as the Mac reference. It is set
explicitly with `-DCMAKE_BUILD_TYPE=Debug`, and GCC's `CMAKE_CXX_FLAGS_DEBUG` defaults to `-g`
(the code optimizes at `-O0`). Release/RelWithDebInfo are **not** used for the first equivalence
test. Optimization is a later, separate qualification (§20).

Required, and recorded in the qualification key:

- `-std=c++17` (from `CMAKE_CXX_EXTENSIONS OFF`); **explicit `-ffp-contract=off`** passed via
  `CMAKE_CXX_FLAGS`. It is inert at `-O0` on baseline x86-64, but pins contraction regardless of
  future optimization level;
- no `-march`/`-mtune`/`-mfma`/`-mavx*`. The login CPU (EPYC 7232P) and compute CPU (EPYC 7352)
  differ, and any host-tuned code is forbidden;
- SSE2 double arithmetic: CQ0 asserts `FLT_EVAL_METHOD == 0` and `__FLT_EVAL_METHOD__ == 0`.

**Forbidden** (CQ0 refuses if any appears in the effective compile or link line, the CMake cache,
or the GSL/Zaki/Confind build records): `-ffast-math`, `-Ofast`, `-funsafe-math-optimizations`,
`-freciprocal-math`, `-ffinite-math-only`, `-fassociative-math`, `-fno-signed-zeros`,
`-fno-trapping-math`, `-fcx-limited-range`, `-ffp-contract=fast`, `-fno-math-errno` added
explicitly, `-march=*`, `-mfma`, and `-mavx*`. Also forbidden are any `#pragma GCC optimize` or
`__attribute__((optimize))` enabling them. CQ0 verifies the effective GCC settings with
`g++ -std=c++17 -O0 -ffp-contract=off -Q --help=optimizers` and `--help=common`, and archives
the output.

**Environment:** `GLIBC_TUNABLES`, `LD_PRELOAD`, and `LD_AUDIT` must be unset in every job. This
mirrors PWNSpall's governed launch environment
(`/mnt/sdd/zaki/PWNSpall/releases/cluster-runtime-04/launch_environment.sh`). glibc selects libm
variants by hardware capability, so the qualification key includes the compute-node `libm.so.6`
hash and CPU family/model (§18).

---

## 6. OpenMP / thread policy

CompactStar links `OpenMP::OpenMP_CXX` publicly (`CMakeLists.txt`). Build parallelism
(`ninja -j N`) and **scientific runtime threading** are different things:

- ADR-0017 §6 states that "No scientific thread-level sharing or parallelism is introduced". The
  accepted Mac qualification ran "one process. No scientific local solve was multithreaded"
  (`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:458-459`).
- Source-level audit (§13.3): **no** OpenMP or thread use on any Phase-5/Phase-6 or test path.
  OpenMP exists only inside the vendored `libConfind.a`, which only the unused `TaskManager`
  calls. The canonical scientific paths are single-threaded.

**Policy:** every scientific Slurm task exports `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=false`,
`OMP_PROC_BIND=true`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and
`NUMEXPR_NUM_THREADS=1`, and runs with `--cpu-bind=cores`. Because `ConstrainCores=no`, an unset
`OMP_NUM_THREADS` would let libgomp start up to 96 threads on a shared node. The explicit setting
is therefore also an isolation requirement. CQ0 refuses a job script that does not set it.
Multithreaded star evolution is **not** pre-authorized. Thread scaling is a later, separate
qualification.

---

## 7. Process-level parallelism, repository layout, scratch

### 7.1 Process model

One immutable run card per Slurm task and process. For every concurrent task the following are
unique: scratch root, output root, stdout/stderr log, run-card identity, and provenance manifest.
Nothing is shared: no candidate JSON, no Phase-5D producer/comparator scratch, and no build tree
writes. All tasks read one immutable executable set (§7.4).

### 7.2 Observed PWNSpall structure (`/mnt/sdd/zaki/PWNSpall`, 283 G)

```text
PWNSpall/
  repo/PWNSpall/              primary clone (origin git@github.com:ZAKI1905/pwn-libeb-spallation.git, main)
  worktrees/PWNSpall-<task>/  ~26 task worktrees on qualification/*, analysis/*, implementation/* branches
  runtimes/<name>-<sha>/      7 detached-HEAD worktrees pinned to exact SHAs
  toolchains/<id>/{source,sources,build,install,manifests,sysroot-rocky9.4}
  releases/<id>/              immutable executables + launch_environment.sh
  qualification/<id>/         durable qualification evidence (slurm logs, probes, prelaunch)
  campaigns/<id>/             production campaign outputs
  scratch/<id>/               transient job data (some detached worktrees also live here)
  acquisitions/, diagnostics/ project-specific data acquisition / diagnostics
  notes/                      chats, prompts, markdown notes
  .agents/ .codex/ .git/      empty directories
```

**Mirror:** `repo/` plus `worktrees/` (one primary clone, and worktrees per task);
`toolchains/<id>/…/manifests` (versioned, hashed dependency prefixes); `releases/<id>/` (build
once, run many; a governed launch environment); `qualification/` versus `campaigns/`
separation; a top-level `scratch/`.

**Do not mirror:**
- `runtimes/` (fold into `worktrees/` with detached HEADs to keep one worktree namespace);
- scratch-hosted worktrees (the source never lives under scratch);
- `sysroot-rocky9.4` (PWNSpall built on the el9.7 login node for el9.4 compute nodes; CompactStar
  builds **on a compute node**, so no sysroot is needed);
- `acquisitions/`, `diagnostics/`, and `notes/` (notes belong in Git/Mac);
- the empty `.agents`/`.codex`/`.git` directories;
- PWNSpall's `origin` over SSH (the cluster has no GitHub SSH access to CompactStar; it clones
  read-only over HTTPS and **never pushes**).

### 7.3 Proposed CompactStar layout (only the top directory exists; nothing below is created)

```text
/mnt/sdd/zaki/CompactStar/                       <- created 2026-09-25 (empty)
  repo/CompactStar/                              primary FULL-history clone (HTTPS; never shallow: tests/rotochemical/manifest.py needs commit f7116c14…; never built in)
  worktrees/CompactStar-cq-<shortsha>/           detached HEAD at the qualified SHA; chmod -R a-w during runs
  toolchains/<toolchain-id>/{source,build,install,manifests}/    GSL 2.7.1, Zaki, Confind, venv+NumPy
  builds/<qualkey12>/                            CMake build tree (Debug); write-protected after CQ1
  releases/<qualkey12>/                          staged immutable executables + launch_environment.sh + SHA256SUMS
  inputs/<input-manifest12>/                     hash-verified EOS/profile/frozen-fixture inputs imported from Mac (outside Git)
  qualification/<cq-campaign-id>/{cq0..cq7}/     durable compact evidence (manifests, results, logs)
  campaigns/<campaign-id>/                       future BA12R durable evidence (after acceptance only)
  scratch/<campaign-or-cq-id>/<jobid>_<arrayidx|0>_<attempt>/    transient job roots
```

`<toolchain-id>` = `gcc14.1.0-gsl2.7.1-zaki<sha7>-confind<sha7>-py3.11.9-numpy<ver>`.
`<qualkey12>` = first 12 hex characters of the qualification key (§18).

Transient scratch: node-local storage is preferred *if* CQ0 shows node-local `/tmp` is disk-backed
with at least 20 GB free. Otherwise, and by default, `/mnt/sdd/zaki/CompactStar/scratch/…`, which
PWNSpall has used for this purpose. CompactStar I/O is small: Phase-6 tables are about 8193 rows.
Raw scratch is preserved until validation and hashing complete (§16).

**Qualification worktree path (proposed):**
`/mnt/sdd/zaki/CompactStar/worktrees/CompactStar-cq-<first 12 hex of the approved SHA>`.

### 7.4 Build once, run many

1. One CQ1 build job produces one immutable executable set.
2. It is staged into `releases/<qualkey12>/` together with `SHA256SUMS`, covering executables,
   `libCompactStar.a`, the linked static libraries, and the Python interpreter/`libpython`/NumPy
   files.
3. The release is made read-only.
4. Every scientific task re-verifies `SHA256SUMS` before it executes.
5. If node-local staging is used, the staged bytes are re-verified against the same sums.
6. There is no per-task rebuild.

---

## 8. Source and input authentication (CQ0 prerequisites)

Before any cluster build:

- cluster worktree `HEAD` = expected canonical SHA = `origin/master` (after a `git fetch` into
  `repo/`) = live `refs/heads/master` from `git ls-remote`; `git status --porcelain` empty. The
  expected SHA is the owner-approved successor of `232565a…` that resolves §2.1.
- **Source manifest:** SHA-256 of `'<sha256>  <path>\n'` records for every tracked file (from
  `git ls-files -s`, with blob contents re-hashed from the worktree) in bytewise path order.
  Exact equality is required for all of the following:
  - the 33 protected Phase-5D paths and the 11 governed baselines
    (`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:79-131`);
  - ADR-0015 `795c8bc8…`, ADR-0016 `b2b0bba0…`, ADR-0017 `9c05fac0…`;
  - `ScaledRKF45.hpp` `27bd510b…`;
  - the Phase-6 BNV production sources (`CompactStar/Physics/BNV/**`).
- **Scientific inputs:** exact byte identity with the Mac authority for the thermal source tree
  `1cdb9850…`, frozen certificate `9217e278…`, frozen coefficients `3efa060d…`, profile root tree
  `23311486…` (`profile.tsv` `e9cd03b0…`, `model.txt` `3ea70de7…`, `freegas.tsv` `7cd44c92…`),
  qualification certificate `7fc892b5…`, the passive schedule `43ec23ad…`, and the EOS data root
  (`COMPACTSTAR_EOS_DATA_ROOT`) tree manifest. These are imported from the Mac into
  `inputs/<id>/` and hash-verified. **Scientific input bytes must match exactly; numerical
  outputs need not.**
- **Literature:** `literature/SHA256SUMS.txt` is outside Git on the Mac
  (`/Users/keeper/Documents/CompactStar/literature`). The cluster does not need the PDFs. CQ0
  records the Mac manifest SHA-256 as authority and requires the Mac-side 22/22 PASS in the same
  Mac session that exports the inputs.

---

## 9. Qualification gates CQ0–CQ7

A failure at any gate stops all later gates. Diagnostic continuation is allowed only if the owner
authorizes it in writing for that gate, and never counts toward acceptance. All thresholds below
are fixed now; **no post-result change is allowed**.

### CQ-B — toolchain bootstrap (prerequisite; separate owner-authorized task)

§2.1 decision implemented; Linux Zaki/Confind; GSL 2.7.1; NumPy venv. Manifests only. No
CompactStar build.

### CQ0 — environment / source authentication (no scientific executable before PASS)

One `phy01` job, `-n1 -c1`, `--time=00:15:00`, capturing:

1. Live compute identity: `hostname`, `lscpu` (vendor, family 23, model 49, stepping, flags),
   `numactl -H` or `lscpu -e`, `free -m`, `df -hT` and `stat -f` of `/tmp`, `/dev/shm`,
   `/mnt/sdd`, `/opt`, and `mount` lines for `/mnt/sdd` and `/opt`. Also `/etc/os-release`,
   `uname -a`, `getconf GNU_LIBC_VERSION`, and `sha256sum /lib64/libm.so.6 /lib64/libc.so.6`.
2. Source SHA/manifest checks and input-hash checks (§8).
3. Toolchain capture (§17) and the forbidden-flag scan (§5); `FLT_EVAL_METHOD==0`.
4. `OMP_NUM_THREADS=1` and the other thread variables set; `GLIBC_TUNABLES`/`LD_PRELOAD`/`LD_AUDIT`
   unset.
5. Explicit Slurm allocation recorded (`scontrol show job $SLURM_JOB_ID` captured in-job).
6. Scratch/output roots confirmed unique and absent before creation.
7. Mac-side reference extraction, from Mac authority and not from any cluster run: minimum
   accepted main step excluding the first; maximum accepted steps per observation interval; and
   the full diagnostic packet. All come from the retained Arm-E/production artifacts
   (`checkpoints.tsv` `a182f07a…`, accepted-step history `fe9afd8c…`) and are frozen into the CQ
   predeclaration file before CQ1.

PASS requires every item. Any mismatch → STOP.

### CQ1 — clean build on a compute node

One `phy01` job: `-N1 -n1 -c16 --mem=32G --time=02:00:00`. It builds a fresh
`builds/<qualkey12>/` with:

```text
cmake -G Ninja -S <worktree> -B builds/<qualkey12> -DCMAKE_BUILD_TYPE=Debug
      -DCMAKE_C_COMPILER=<spack gcc>/bin/gcc -DCMAKE_CXX_COMPILER=<spack gcc>/bin/g++
      -DCMAKE_CXX_FLAGS="-ffp-contract=off" -DCMAKE_C_FLAGS="-ffp-contract=off"
      -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc"
      -DGSL_ROOT_DIR=<gsl-2.7.1> -DGSL_LIBRARY=<…>/libgsl.a -DGSL_CBLAS_LIBRARY=<…>/libgslcblas.a
      -DPython3_EXECUTABLE=<venv>/bin/python -DCOMPACTSTAR_EOS_DATA_ROOT=<inputs/…>
      [+ the §2.1-approved Zaki/Confind override]
ninja -j 16
```

PASS criteria:
- configure and build exit 0;
- `CMakeCache.txt` and `compile_commands.json` archived; zero forbidden flags;
- `readelf -d` shows no unexpected `RPATH`/`RUNPATH` outside `/opt/spack` Python and the venv;
- `ldd` resolves every library on the compute node;
- max RSS and wall time recorded with `/usr/bin/time -v`;
- the build tree and staged release are made read-only and `SHA256SUMS` is written.

Build parallelism is allowed. This is not scientific runtime parallelism.

### CQ2 — focused / data-free tests

See §13.1 for the exact list. It runs on the CQ1 release.

### CQ3 — Phase-5 cross-platform qualification

See §12.

### CQ4 — bounded ADR-0017 Phase-6 cross-platform qualification

See §10–§11.

### CQ5 — concurrency / isolation

See §14.

### CQ6 — reproducibility

See §15.

### CQ7 — platform ACCEPT/REFUSE

See §21.

---

## 10. CQ4 — the bounded Phase-6 reference fixture (exact)

The fixture is identical to the accepted Mac qualification
(`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:185-216`):

- `CPL-P2-LINEAR-QSS-v1`, source ON, P2, spin OFF, Me/Mmu ON, De/Dmu OFF;
- `t = 0` → `462269531250 s`; initial `(x, eta_e, eta_mu) = (0, 0, 0)`;
- one uninterrupted adaptive RKF45 main integration, `rtol=1e-11`, `atol=(1e-16,1e-22,1e-22)`,
  initial `h = 1 s`, exactly one positive GSL `t1 = 462269531250 s`;
- the passive 241-point schedule, SHA-256 `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e`;
- self-qualified two-level isolated rk8pd reconstruction: O1 `rtol=1e-12`,
  `atol=(1e-17,1e-23,1e-23)`; O2 `rtol=1e-13`, `atol=(1e-18,1e-24,1e-24)`;
- the frozen fixture inputs imported byte-identically from the Mac (§8), so only platform numerics
  differ.

Mac reference values:

| Quantity | Mac value |
|---|---|
| final `x` / `eta_e` / `eta_mu` | `0.49240008824076903` / `-2.5123474256442210e-7` / `-4.7906773046561003e-7` |
| accepted / rejected main steps | `232 / 60` (diagnostic only) |
| internal accepted-step SHA-256 | `fe9afd8c…` (diagnostic only) |
| strict-interior checkpoints | 239/239 self-qualified; knot observations 82, 117, 228 |
| max `d_O/D_O1`, max `U_O/(0.20 F_i)` | `0.5794839113173227`, `0.5793505315921852` (obs 117, `x`) |
| `R20_residual`, `N_R20`, `R20_normalized` | `5.7186274768377777e39 erg`, `1.0940924194731047e46 erg`, `5.2268230499136139e-7` |
| reconstruction uncertainty | `7.0988433612780846e31 erg` |
| Mac wall time | fixture construction 264.3 s; main 23.15 s; checkpoint output 1.81 s |

Mac-vs-cluster values are compared on parsed max-digits-10 binary64 values (`Q=0`).

**Implementation prerequisite:** the accepted verifier `tests/bnv/adr0017_production_verify.py`
requires bit identity with the Mac oracle and must remain unchanged. A **separate**
cluster-qualification comparator implementing §11 must be written, reviewed, and accepted in a
bounded implementation task before CQ4 executes (B6).

---

## 11. Cross-platform acceptance formulas (predeclared; derived from governed budgets)

### 11.1 Governed scales reused (no new constants minted)

- The BA12R/ADR-0017 state resolution scale
  `F_i = max(atol_ULTRA,i + rtol_ULTRA M_i, 64 ulp(M_i), Q_i)` with `rtol_ULTRA=1e-11` and
  `atol_ULTRA=(1e-16,1e-22,1e-22)` (ADR-0017 §5.3;
  `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:326-340`).
- BA12R's floor-limited allowance `10 F`, below which the campaign declares differences
  unresolvable and forms no contraction ratio
  (`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:308-319`).
- The refined-stability scale `D_R = atol_REFINED + rtol_REFINED M ≈ 100 D_U`.
- The adjacent-tier RKF4(5) global-error factor `0.01^(4/5) = 0.0251188643150958`
  (`…RECOVERY_PREFLIGHT.md:321-329`).
- The ledger power scale `G_P` and `F_P = max(1e-11 G_P, 64 ulp(M_P), Q_P)`, with the refined/ULTRA
  ledger stability ratio `1e-9` (`…RECOVERY_PREFLIGHT.md:333-354`;
  `…CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:356-373`).
- R20: governed gate `|R20|/N_R20 ≤ 2e-4`; reconstruction subsidiary bound `5e-6 N_R20`.

### 11.2 Why `K = 10`

A platform change acts like re-solving the same IVP with a perturbed accepted-step sequence. Two
ULTRA-tier solutions of the same problem should differ by at most about twice the ULTRA global
error. From the governed contraction logic, the ULTRA global error is about
`0.0251 × D_R ≈ 2.51 D_U`, so the expected cross-platform difference is `≲ 5.0 F`.

`K = 10` gives about a 2× margin over that estimate. It is **exactly the governed floor-limited
allowance `10F`**, the level BA12R already treats as numerically unresolvable and scientifically
meaningless. It is also about 10× below the refined stability scale `D_R`, which is the smallest
difference BA12R treats as a meaningful tier signal. A platform discrepancy that passes
therefore cannot masquerade as a BA12R convergence or BNV signal.

`K` is fixed from this reasoning and is not selected after cluster results.

### 11.3 State agreement (every observation `j` = 0…240, component `i ∈ {x, eta_e, eta_mu}`)

```text
Delta_ij  = |y_cl(i,j) - y_Mac(i,j)|
M_ij      = max(|y_cl(i,j)|, |y_Mac(i,j)|)
F^X_ij    = max(atol_ULTRA,i + rtol_ULTRA * M_ij, 64 ulp(M_ij), Q_ij)        Q_ij = 0
U_Mac,ij  = Mac O2 reconstruction uncertainty 2 max(d_O,F_O)  (0 at exact endpoints)
U_cl,ij   = cluster O2 reconstruction uncertainty            (0 at exact endpoints)

REQUIRE   Delta_ij + U_Mac,ij + U_cl,ij  <=  Tau_platform,ij = 10 * F^X_ij
```

Special cases:
- Observation 0 (the initial state) must be **bit-identical**.
- Observation 240 (the terminal main endpoint on both platforms) uses the main state with
  `U = 0`.
- `ulp(0)` is the binary64 spacing at zero.

This is the **individual-component** rule. There is no aggregation and no averaging.

### 11.4 Diagnostics agreement (each of the 241 checkpoints, from the governed Phase-6 packet)

| Class | Quantities | Agreement requirement (Mac vs cluster) |
|---|---|---|
| P — powers/luminosities | `P_dir_eq`, `P_dir_actual`, `L_H`, `DeltaLnu`, `DeltaPbeta`, `Lnu_eq`, `Lnu_full`, `L_out_fluid` (if emitted), `Lgamma`, `Lother`, `Pnet` | **`|O_cl - O_Mac| <= max(1e-9 * G_P,Mac, 64 ulp(max(|O_cl|,|O_Mac|)))`**. `G_P,Mac` is the Mac packet's governed gross power scale at that checkpoint. `1e-9 G_P` is BA12R's governed refined/ULTRA ledger constant, reused unchanged. The class-T `1e-9` is a ceiling, not the expected difference. Expected state agreement is `≲ 5F` (§11.2), about `5e-11` relative in `x`. Propagated through the steepest power law (`Lnu ∝ T^8`), that is `≲ 4e-10` relative, which stays within `1e-9 G_P` because each power is at most `G_P`. A T difference near the class-T ceiling could fail class P. That is intended: class P is then the binding, stricter test |
| S — state-identical diagnostics | `sigma` components where the owner defines them as the `eta` state components | the §11.3 state rule |
| T — scalar thermodynamic diagnostics | `Tinf`, `Tsurface_inf` (if emitted), `eta`/`xi` derived scalars, `mu_n_actual`, `mu_B`, `Echem`, `Uth` (if emitted), Cstar, baryon count, and frozen-validity numeric margins | `|O_cl - O_Mac| <= max(1e-9 * max(|O_cl|,|O_Mac|), 64 ulp(max(|O_cl|,|O_Mac|)))`. `1e-9` = 10 × (state relative scale `1e-11` at the rtol-dominated limit) × a condition-number allowance of 10. It equals BA12R's governed `1e-9` refined/ULTRA ledger ratio |
| E — exact semantic | source/domain/revision/partition/product-fate identities, currentness flags, frozen-validity **flags**, `MAIN_ENDPOINT`/`RK8PD_RECONSTRUCTED` classification, schedule size/hash | exact equality |

A quantity absent from the accepted Mac packet has no agreement requirement. Its individual-run
gate still applies.

### 11.5 Individual-run acceptance versus agreement (both required, neither implies the other)

| Quantity | Cluster individual gate | Mac-vs-cluster agreement |
|---|---|---|
| R20 | `|R20_cl|/N_R20,cl <= 2e-4` **and** cluster reconstruction uncertainty `<= 5e-6 N_R20,cl` | `|R20n_cl - R20n_Mac| <= 5e-6`, where `R20n = R20_residual/N_R20`. Also `|N_R20,cl - N_R20,Mac| <= 1e-9 N_R20,Mac` |
| R18 | the governed nonzero-eta/roundoff budget, applied **only where the harness emits it** (the bounded Mac record reports none) | none. R18 is a roundoff-level residual, and agreement between two roundoff realizations is not meaningful |
| ADR-0017 self-qualification | §11.6 | none (platform tolerance cannot substitute for it) |

### 11.6 ADR-0017 self-qualification must pass on the cluster

Every cluster strict-interior checkpoint must itself satisfy `d_O <= D_O1` and
`U_O <= 0.20 F_i` for every component, under the exact O1/O2 configurations. Every observation
whose bracketing cluster accepted step crosses a Cstar knot must pass under the same rule. The
cluster's knot set may differ from Mac's 82/117/228 if steps differ, and is recorded. No
cross-platform tolerance overrides a cluster self-qualification failure.

### 11.7 Main-integrator health (cross-platform; step history not required identical)

All of the following are required:

1. There is no GSL failure, nonfinite state, currentness refusal, frozen-validity refusal, or
   100000-step guard.
2. The terminal time is exactly `462269531250 s`, and there is exactly one positive GSL `t1`.
3. The schedule hash is `43ec23ad…`, there are 241 observations, and observation 240 coincides
   with the terminal accepted endpoint.
4. The accepted count satisfies `174 <= N_acc <= 290` (Mac `232` ± 25 %).
5. The rejected count satisfies `N_rej <= 120` (2 × Mac `60`), and the rejection fraction satisfies
   `N_rej/(N_acc+N_rej) <= 0.30` (Mac `0.2055` + 0.10, the same additive allowance as the
   BA12R §6.2 rejection threshold).
6. The minimum accepted step excluding the first satisfies
   `>= 1e-3 ×` the Mac minimum accepted step excluding the first (Mac value frozen in CQ0 step 7).
7. The maximum accepted steps within any observation interval is `<= 2 ×` the Mac maximum
   (frozen in CQ0 step 7).

If the cluster history happens to be identical to Mac's, that is recorded but not required.

---

## 12. CQ3 — Phase-5 cross-platform qualification

### 12.1 Scope

**Scope boundary (declared):** CQ3 qualifies the cluster as a structural/ODE solve platform that
consumes **Mac-generated, hash-identical** EOS/profile/thermal tables (§8). Cluster-side
*generation* of EOS tables is **not** qualified by this plan. It is also not part of the cluster
production scope: the cluster never generates scientific EOS/profile/frozen inputs for campaigns.
Any future cluster-side input generation requires its own qualification.

**Cluster-mode producers are required as well as comparators.** The governed producers cannot
run on the cluster:
- the 5B producer calls the hard-coded `/usr/bin/clang++`;
- the 5D `qualified_suite.py` hard-codes regenerated-table SHA-256 values and re-runs the 5B/5C
  governed regressions live.

Several affected files are among the 33 protected Phase-5D paths:
`tests/analysis/produce_particle_number_reference.py`,
`tests/analysis/chemical_production_evidence.py`,
`tests/analysis/chemical_production_fixture.cpp`, and `tests/analysis/chemical_trackr_budget.py`.
The cluster-qualification producers and comparators must therefore be **new, parallel files**
under a new cluster-qualification test directory. **No protected or governed file is modified.**
They are implemented and reviewed in a bounded, owner-accepted implementation task before CQ3.

The three governed regressions: Phase-5B structural response, Phase-5C Z/W coefficients, and
Phase-5D controlled evolution. The governed comparators
(`tests/analysis/particle_number_response_regression.py`,
`tests/analysis/chemical_coefficient_regression.py`,
`tests/rotochemical/phase5d1_controlled_evolution_regression.py`/`compare_artifacts.py`) and all
baselines stay **unchanged**. On the cluster they are expected to fail on platform fields
(`architecture`, `platform`, `GSL`) and on bit-exact doubles. That failure is **informational**
and is neither a scientific failure nor grounds to touch a baseline.

A **new, separate cluster-qualification comparator** is required. It is implemented and reviewed
in a bounded task before CQ3. It compares the cluster-produced artifacts against the governed
Mac baselines using this field classification:

| Class | Fields |
|---|---|
| BYTE_IDENTITY_REQUIRED | source hashes (5B: 19 files; 5C: 17; 5D: 91 production + 33 protected + 10 baseline SHAs); governed baseline files; imported scientific input bytes |
| EXACT_SEMANTIC_IDENTITY | schema/status/units/fixture/classification strings, `predeclaration_sha`, `goals`, `EOS_revision`, key sets, list lengths (for example the 5C certificate interval count), refusal lists, `solver.configuration` tolerances, and 5D output-cadence count (402) |
| PLATFORM_METADATA_DIFFERENCE_ALLOWED (recorded, non-empty, truthful) | 5B `build.compiler`/`architecture`; 5C `toolchain.{compiler,platform,architecture,GSL}`; 5D execution sidecar; generated-intermediate hashes (`EOS_table_sha256`, `profile_sha256`, `table_sha256`, `partition_sha256`, 5D logical qualification hashes). Their *content* is compared numerically instead |
| NUMERICAL_TOLERANCE | as below |

Numerical tolerances. Each cluster run must *also* pass its own governed per-run budgets unchanged:

- **5B** (`A`, `B`, `K`, `N`, `I_phys`, each with a declared `*_error`):
  `|X_cl - X_Mac| <= E_X,Mac` (the Mac baseline's declared numerical error for that coefficient),
  and `E_X,cl` satisfies the same governed goals. Fields with no declared error are compared under
  class T of §11.4 (relative `1e-9`, 64-ulp floor).
- **5C:**
  - `|Z_cl - Z_Mac| <= E_Z,Mac` componentwise (Mac `E_Z` =
    `[[4.499198710915918e-61, 4.285424700647402e-61],[4.285424702045283e-61, 2.908175030741565e-59]]`);
  - `|W_cl - W_Mac| <= E_W_numerical,Mac` = `[1.59100967128208e-13, 1.4843208824397632e-12]`;
  - G, Q, I by their declared `E_*`;
  - the cluster must satisfy `E_W_numerical <= G_W_numerical = [1e-12, 4e-12]` and
    `V_W <= G_W_validation = [1e-10, 3e-10]`
    (`docs/validation/PHASE5C2_PRODUCTION_ACCEPTANCE_PREDECLARATION.md:39-68`);
  - undeclared-error fields use class T.
- **5D:**
  - trajectory gating quantities `q ∈ {Tinf_K, eta_e_MeV, eta_mu_MeV}`, using the governed
    metrics (`|ΔT|/|T|`; `|Δeta|/max(|eta_Mac|, 1e-10 MeV)`), maximized over all 402 checkpoints:
    `Tau_5D,q = min(2e-4, 10 × e_BR,q)`.
    - `e_BR,q` is the **achieved** Mac baseline-vs-refined metric in the governed baseline
      `tests/baselines/phase5d1_controlled_evolution.json` → `convergence.ODE`: `Tinf_K`
      `1.7408584682611815e-05`, `eta_e_MeV` `1.5866651419012378e-07`, `eta_mu_MeV`
      `1.9479546821964215e-07`.
    - Resulting exact budgets: **Tinf `1.7408584682611815e-4`; eta_e `1.5866651419012378e-6`;
      eta_mu `1.9479546821964215e-6`**.
    - Reasoning (the same model as §11.2): the governed checkpoints are the baseline-tier
      (`rtol=1e-7`) run, whose global error is about `e_BR` because the refined error is
      smaller. Two baseline-tier solutions on different platforms therefore differ by `≲ 2 e_BR`,
      and `K = 10` leaves a 5× margin. The cap at the governed `2e-4` means platform agreement can
      never be looser than the governed accuracy gate.
    - All other 5D trajectory quantities (powers, `Pnet`, `x_dot`, and so on) are reported but not
      gating. The governed 5D gates are T and eta only;
  - ledger residuals and convergence summaries must pass the governed validator budgets on the
    cluster run itself;
  - step counts (Mac 8828/10415 accepted, 1618/2183 rejected) are diagnostic only, with the §11.7
    health rules scaled to the 5D run (±25 % accepted, rejection fraction ≤ Mac + 0.10);
  - the 5D producer's hard-coded intermediate hashes and its live 5B/5C gate make the *governed*
    producer unusable on the cluster. The cluster comparator consumes a cluster producer run that
    uses the **Mac-imported** frozen inputs (§8), so input bytes are exact.

---

## 13. CQ2 focused tests, full-suite policy, threading audit

### 13.1 Inventory facts (read-only audit of `tests/CMakeLists.txt` at `232565a…`)

- **76** registered ctest identities, not 77: 52 always registered (data-free) plus 24 registered
  only when `COMPACTSTAR_EOS_DATA_ROOT` exists (`tests/CMakeLists.txt:370,561,637`). The
  "full 77-test campaign" wording (`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:304,492`)
  refers to the unmerged Phase-6A1 controlled-BNV test set, which is absent from canonical. This is
  recorded as a documentation discrepancy, not reconciled.
- There is no `RESOURCE_LOCK`, `WORKING_DIRECTORY`, or `DEPENDS`. There are 13 `RUN_SERIAL` tests
  and one fixture pair (`chemical_trackr_budget` → `chemical_production_validation`).
- **There is no registered Phase-6 BNV source/projection or ADR-0016 tangent-adapter test at
  canonical.** `adr0017_production_qualification` is built but not `add_test`-ed; it takes 10
  external arguments. BNV source/projection is therefore exercised on the cluster only through the
  CQ4 actual-fixture run. This is recorded as a coverage gap for the owner.
- **Shared-root hazards (require serialization or isolation):**
  - all 12 `tov_surface_*` tests use the hard-coded `/tmp/compactstar-tov-surf-ir/`
    (`tests/core/tov_surface_contract.cpp:358,500`);
  - about 15 tests use fixed names under `std::filesystem::temp_directory_path()` with
    `remove_all`, for example `heat_capacity_v1.cpp:273`, `photon_cooling_conformance.cpp:159`,
    and `cache_thermal_contract.cpp:469`;
  - `chemical_production_validation.py:24-26` picks the newest `<bin>/phase5c-evidence/run-*`;
  - `phase5d_protected_manifest` runs `git worktree add/remove` **on the source repository**
    (`tests/rotochemical/protected_manifest.py:51-70`);
  - `manifest.py` needs full history (commit `f7116c14…`); a shallow clone fails;
  - the Phase-5D regression requires `git status --porcelain --untracked-files=all` to be empty,
    and it builds a nested tree.
- The 18 Python tests need `numpy`, `mpmath`, and `scipy`. The CQ-B venv therefore pins all three.
- The only portability-risk construct noted is `using enum StateTag;` (C++20) in `Tags.hpp:72`.
  AppleClang accepts it in C++17 mode. **Whether GCC 14 accepts it under `-std=c++17` is
  UNVERIFIED.** A CQ1 compile error here stops CQ1 and returns to the owner; it is not patched
  inside qualification.

### 13.2 CQ2 exact focused test set (data-free; `COMPACTSTAR_EOS_DATA_ROOT` not needed)

The set is chosen to prove each required property with no multi-hour stellar calculations:

| Property | ctest identities |
|---|---|
| binary loads / library links | `compactstar_library_smoke` |
| core unit semantics | `relativistic_unit_boundary`, `relativistic_unit_background` |
| local thermodynamics / EOS contracts | `rotochemical_local_thermodynamics`, `rotochemical_trackr_freegas_local`, `rotochemical_trackr_npe`, `rotochemical_trackr_pe`, `eos_derivative_contract` |
| thermal / cache contracts | `heat_capacity_v1`, `evolution_stepper_contract`, `photon_cooling_conformance`, `cache_contract`, `cache_thermal_contract` |
| Phase-5 analytic / oracle contracts | `particle_number_analytic`, `chemical_production_contract`, `chemical_exact_oracles`, `chemical_curved_gc9`, `phase5d_response`, `phase5d_independent_oracles`, `phase5d_component_tolerances` |
| provenance / currentness refusals | `phase5d_harness_controls` |
| ADR-0017 checkpoint machinery (P1–P10) | `passive_checkpoint_output_contract` |

That is 22 tests. On the Mac, all except `heat_capacity_v1` (10.28 s) ran in under about 9 s each
(`docs/validation/phase5b_resume_evidence.json:661-721`).

Execution:
- one job with `-n1 -c1`;
- serial `ctest -R '^(…22 names…)$' -j1 --output-on-failure`;
- `TMPDIR=<job scratch>/tmp`, which isolates the fixed `temp_directory_path()` names.

Process parallelism is **not** used in CQ2: the total is short, and several of these tests share
fixed temp names.

None of the 22 takes a baseline argument (`tests/CMakeLists.txt`: baseline paths appear only for
the three governed regressions, `passive_cooling_regression`, `hartle_monopole_regression`, and
`baryon_number_cmf`, all outside CQ2). A spot-check of the Python sources found `tests/baselines`
only as a refusal-string check in `harness_controls.py:67-69`. CQ2 is therefore expected to be
platform-independent.

PASS requires all 22 to pass. **Any failure stops CQ2 and returns to the owner for adjudication.**
No test is reclassified as "platform-expected" after the fact.

### 13.3 Threading audit result

There is no `#pragma omp`, `omp_*`, `<omp.h>`, `std::async`, TBB, or parallel-algorithm use
anywhere in `CompactStar/`, `tests/`, or `main/`. No code reads `OMP_NUM_THREADS`.

- The only `std::thread` is `Core::TaskManager` (`CompactStar/Core/src/TaskManager.cpp:47-53,112`,
  sized by `hardware_concurrency()`). No Phase-5, Phase-6, or test path references it.
- OpenMP appears only inside the vendored `libConfind.a` (`__kmpc_fork_call`), which only
  `TaskManager` calls.
- The Phase-6 BNV, Rotochemical (including `ScaledRKF45.hpp`), and Analysis sources contain no
  threading.

**Conclusion:** the canonical scientific paths are single-threaded. OpenMP is a link and load-time
dependency only, and `OMP_NUM_THREADS=1` matches the qualified Mac execution. The Linux
`libConfind.a` rebuild in CQ-B must keep its OpenMP semantics unchanged. They are unused by the
qualified paths.

### 13.4 Full-suite policy — recommendation **B**

**Recommendation:** focused qualification first; the full suite after CQ4 passes, as a
prerequisite of CQ7.

**Why not A (full suite as the first gate):**
- about 4717 s of Mac serial runtime (75/75,
  `docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md:169`), and more in
  cluster Debug;
- it overlaps heavily with CQ2/CQ3;
- many tests compare against Mac-generated artifacts bit-exactly, so they validate same-platform
  identity rather than cross-platform numerics.

**Why not C (skip the full suite):** a substitute cannot be defined without the per-test
comparator classification in the next paragraph.

**Classification before execution:** the full-suite task begins with a **static comparator
classification of all 76 tests**, frozen and committed before any execution. Each test is classed
from source as either:
- `ANALYTIC/INTERNAL` (checks analytic limits, identities, or internally generated oracles), or
- `MAC-ARTIFACT-EXACT` (bit-exact comparison with a Mac-produced baseline, such as `*_debug.tsv`
  or the three governed Phase-5 regressions).

`ANALYTIC/INTERNAL` tests must PASS on the cluster. `MAC-ARTIFACT-EXACT` failures are informational
only where the same quantity is covered by a CQ3/CQ4 tolerance comparator. Otherwise the gap is
reported to the owner before CQ7.

**Safe process-parallel execution** replaces the historical indiscriminate `ctest -j1`:
- **Lanes** use independent *copies* of the release build tree (`cp -a` to lane scratch, verified
  against release hashes), each with its own `TMPDIR`. Execution within a lane is serial.
- **Lane T:** all 12 `tov_surface_*`, because of the hard-coded `/tmp/compactstar-tov-surf-ir`.
  At most one lane T per node.
- **Lane 5B:** `phase5b_*` and `particle_number_*`.
- **Lane 5C:** `chemical_*` and `phase5c_*`, with the fixture pair kept in order.
- **Lane 5D:** `phase5d_*` except the protected-manifest test, plus the Phase-5D regression.
- **Lane D:** the remaining data-free tests.
- **Lane E:** the remaining EOS-data tests.
- **Excluded from the shared worktree:** `phase5d_protected_manifest` runs only against a
  disposable full clone in its own scratch, because it creates git worktrees.
- Lanes run as separate `-n1 -c1` jobs. Lanes 5B, 5C, and 5D are serialized relative to each other
  because the 5D regression re-runs the 5B and 5C regressions.

---

## 14. CQ5 — concurrency / isolation experiment

Everything runs on the **same** CQ1 release, on `phyn001` (`-w phyn001`), with `OMP_NUM_THREADS=1`:

1. **Serialized reference S:** one CQ4 fixture run alone (`-n1 -c1`).
2. **Concurrent pair C1, C2:** a later single job with `-N1 -n2 -c1 -w phyn001`, launched as two
   `srun --exact -n1 -c1 --cpu-bind=cores` steps in the background and joined with `wait`. They
   use separate scratch/output roots and identical run cards.

PASS requires S, C1, and C2 to have **bit-identical** parsed scientific values (all states,
O1/O2 values, diagnostics, R20) and identical accepted-step history hashes. The code is
single-threaded and deterministic, and the binary, libm, and node are the same. If exact identity
fails, CQ5 FAILS and the cause must be explained (for example nondeterministic library behaviour)
before any fallback to the §11 numerical rule, and that fallback needs owner authorization.

## 15. CQ6 — reproducibility

- **R1 (same node):** a second independent serialized CQ4 run on `phyn001`. It must be
  **bit-identical** to S.
- **R2 (cross-node, same node class):** one CQ4 run on `phyn002`, and on `phyn004` if production
  will use it. It must be bit-identical to S. That is expected because the CPU model, libm, and
  binary are the same. A difference means the node classes are not equivalent: that node is
  excluded pending explanation.
- Exact repeatability failing ⇒ CQ6 FAIL. Qualification must explain the cause before production.

---

## 16. Scratch / output / evidence policy

- Scratch root:
  `/mnt/sdd/zaki/CompactStar/scratch/<campaign-or-cq-id>/<SLURM_JOB_ID>_<SLURM_ARRAY_TASK_ID|0>_a<attempt>`.
  It is created by the job with `mkdir` (no `-p` beyond the campaign root), and the job refuses if
  the root already exists.
- All paths are absolute. No output depends on the CWD, and the job does `cd` into its own
  scratch root.
- On success, compact durable artifacts are copied to `qualification/<id>/cqN/…` or
  `campaigns/<id>/runs/<run-card-id>/a<attempt>/`: the manifests, result tables, checkpoint
  tables, logs, and SHA256SUMS. The copy is hash-verified. Raw scratch is kept until the owner
  accepts the gate.
- Durable evidence is copied back to the Mac and committed to Git, as compact documents and JSON
  results. Large raw outputs stay outside Git.
- No job writes to `repo/`, `worktrees/`, `builds/`, `releases/`, `inputs/`, or another job's
  root.

## 17. Provenance schema (every cluster job evidence packet)

`sacct` is unavailable (accounting disabled), so everything is captured **in-job** by a wrapper
into `provenance.json`, plus `start.json`/`end.json` sentinels:

- `canonical_sha`, `source_manifest_sha256`, `input_manifest_sha256`, `run_card_id`, `run_card_sha256`;
- `release_sha256sums_sha256` and per-binary SHA-256 (executables, `libCompactStar.a`, Zaki,
  Confind, `libgsl.a`, `libgslcblas.a`, Python, `libpython`, NumPy);
- `toolchain_manifest_sha256`: absolute paths and `--version` for `g++`/`gcc`, `cmake`, `ninja`,
  and `python`; the GSL version from `gsl_version.h` plus the library hash; the libgomp path and
  hash if dynamically linked; `module list` output; the CMake configure command;
  `CMakeCache.txt` SHA-256; `compile_commands.json` SHA-256; effective compile/link flags; build
  type;
- `env`: `OMP_*`, `*_NUM_THREADS`, `GLIBC_TUNABLES`, `LD_*`, `PATH`, `PYTHON*`, and
  `COMPACTSTAR_*` (full values);
- OS/kernel (`uname -a`, `/etc/os-release`), glibc version, `libm.so.6`/`libc.so.6` SHA-256, CPU
  vendor/family/model/stepping/flags, hostname;
- Slurm: `SLURM_JOB_ID`, `SLURM_ARRAY_JOB_ID`, `SLURM_ARRAY_TASK_ID`, `SLURM_JOB_PARTITION`,
  `SLURM_JOB_NODELIST`, `SLURM_CPUS_ON_NODE`, `SLURM_JOB_CPUS_PER_NODE`, and a captured
  `scontrol show job $SLURM_JOB_ID`;
- start and end UTC time, exit code, and the terminating signal if any (a `trap` on SIGTERM
  records `TERMINATED_BY_SIGNAL`); `/usr/bin/time -v` output (max RSS, wall, user, sys);
- output file SHA-256 list; validation disposition (`PASS`/`SCIENTIFIC_FAIL`/`EXTERNAL_INTERRUPTION`).

Shell history is never provenance.

## 18. Qualification version key and invalidation triggers

`CQKEY = SHA-256` of canonical JSON (sorted keys) containing:

| Key component | Invalidates on change? |
|---|---|
| canonical CompactStar SHA and source-manifest hash | yes |
| critical input manifest hash (§8) | yes |
| compiler identity: `g++ -v` full string + `g++`/`cc1plus` SHA-256 | yes (any version, including patch, or bytes) |
| C++ runtime linkage mode (static libstdc++/libgcc) | yes |
| GSL version + `libgsl.a`/`libgslcblas.a` SHA-256 + GSL build flags | yes |
| Zaki/Confind source SHAs + library SHA-256 | yes |
| Python version + interpreter/`libpython` SHA-256 + NumPy version | yes |
| CompactStar effective compile/link flags + build type | yes |
| compute OS image: `/etc/os-release` `VERSION_ID`, glibc version, `libm.so.6` and `libc.so.6` SHA-256 | yes (an OS update that changes libm/libc requires requalification) |
| CPU vendor/family/model (AMD 23/49, EPYC 7352) | yes (a new node class requires requalification) |
| ADR-0017 SHA-256 `9c05fac0…` + Phase-6 BNV production source hashes | yes |
| kernel version, hostname, CMake/Ninja version | **no** (recorded only). CMake/Ninja matter only through the resulting binary; a rebuild produces a new release whose binaries must be re-hashed |

The release binds `CQKEY` to its binary SHA-256 set. Any rebuild under the same `CQKEY` must reproduce
identical binary hashes. Otherwise it counts as a new release that requires CQ1, CQ2, CQ4, CQ5, and
CQ6 again. Hostname is never a key.

## 19. Slurm resource model (initial; `phy01` unless noted)

Allocation facts: `CR_CORE` means `-c1` allocates one physical core (two hardware threads).
Memory is not enforced (`ConstrainRAMSpace=no`), so `--mem` is a declared expectation used for
node selection, and packing is controlled by the wrapper. The Mac evidence records no RSS, so the
memory values below are **provisional resource requests, not tolerances**. CQ1/CQ2 measure max RSS
with `/usr/bin/time -v`. From then on the production request is `ceil(1.5 × measured max RSS)`, and
per-node concurrency is `min(48, floor(0.8 × 125000 MB / that request))`.

| Job | partition | nodes | ntasks | cpus-per-task | mem | time | array |
|---|---|---|---|---|---|---|---|
| CQ-B bootstrap build | phy01 | 1 | 1 | 16 | 32G | 04:00:00 | — |
| CQ0 authentication | phy01 | 1 | 1 | 1 | 2G | 00:15:00 | — |
| CQ1 build | phy01 | 1 | 1 | 16 | 32G | 02:00:00 | — |
| CQ2 focused tests | phy01 | 1 | 1 per test (process-parallel only where §13 allows) | 1 | 4G | 01:00:00 | — |
| CQ3 Phase-5 regressions | phy01 | 1 | 1 (serial 5B → 5C → 5D in one job; RUN_SERIAL semantics) | 1 | 8G | 12:00:00 (Mac 5D ≈ 2043 s; 5B/5C TIMEOUTs 900/3600 s) | — |
| CQ4 bounded Phase-6 | phy01 `-w phyn001` | 1 | 1 | 1 | 4G | 02:00:00 (Mac ≈ 290 s) | — |
| CQ5 concurrency | phy01 `-w phyn001` | 1 | 2 | 1 | 8G | 02:00:00 | — |
| full suite (§13.4), per lane | phy01 | 1 | 1 | 1 | 8G | 12:00:00 (Mac full suite ≈ 4717 s serial) | — |
| CQ6 repeat / cross-node | phy01 `-w phyn001` / `-w phyn002` | 1 | 1 | 1 | 4G | 02:00:00 | — |
| future six-run BA12R | phy01 | 1 per task | 1 | 1 | 4G (replaced by the measured rule) | 12:00:00 (Mac segmented ULTRA 344 s + 8191-checkpoint estimate 62 s + fixture 264 s) | `--array=0-5%6` |

Every scientific job also uses `--hint=nomultithread --cpu-bind=cores` and the `OMP_NUM_THREADS=1`
environment. `phyn003` is excluded while drained. `phy02/phyn004` is used only if CQ6 R2 includes
it. No whole-node requests are made except where the owner separately authorizes them.

**Safe process concurrency per node:** 48 (physical cores), capped further by the measured-RSS rule.

## 20. Future six-run BA12R job array (design only; nothing submitted)

The same release, canonical SHA, fixture family, passive 8193-point schedule, and ADR-0017
reconstruction are used for all six runs. The run cards are immutable, and their bytes and
SHA-256 are created and owner-reviewed in the campaign-preparation task. Common parameters come
from `docs/validation/PHASE6A1_CONTROLLED_BNV_BA12R.md:112-125`: `CPL-P2-LINEAR-QSS-v1`, P2, spin
OFF, Me/Mmu ON, De/Dmu OFF, duration `5.0e5 yr`, 8193 outputs, initial `Tinf=1e8 K`,
`eta_npe=eta_npmu=0`.

| `SLURM_ARRAY_TASK_ID` | Run-card identity | Tier `rtol` / `atol(x,eta_e,eta_mu)` | Source |
|---|---|---|---|
| 0 | `CS-P6-BA12R-CLEAN-v1-BASELINE-SOURCE` | `1e-7` / `(1e-12,1e-18,1e-18)` | ON, `Bdot = -2.4136520263641375e37 count/s` (`Bdot/B0=-1.0e-12 yr^-1`) |
| 1 | `CS-P6-BA12R-CLEAN-v1-BASELINE-CONTROL` | same | OFF (`Bdot = 0`) |
| 2 | `CS-P6-BA12R-CLEAN-v1-REFINED-SOURCE` | `1e-9` / `(1e-14,1e-20,1e-20)` | ON |
| 3 | `CS-P6-BA12R-CLEAN-v1-REFINED-CONTROL` | same | OFF |
| 4 | `CS-P6-BA12R-CLEAN-v1-ULTRA-SOURCE` | `1e-11` / `(1e-16,1e-22,1e-22)` | ON |
| 5 | `CS-P6-BA12R-CLEAN-v1-ULTRA-CONTROL` | same | OFF |

Each task maps its index to exactly one run-card file, verified by SHA-256 before execution. It
uses a unique scratch/output root and never modifies shared source or evidence.

## 21. Failure and resubmission policy

- **SCIENTIFIC FAILURE** means any of: a nonzero exit from the scientific executable with a
  scientific cause, solver refusal, nonfinite state, validity/currentness refusal, a
  self-qualification failure, a failed §11 comparison, or a failed health threshold. Any of these
  means **no retry and no settings change**. The evidence is retained with disposition
  `SCIENTIFIC_FAIL` and the stop is reported.
- **EXTERNAL INTERRUPTION** means any of: a node failure or reboot, the job disappearing without
  `end.json`, a SIGTERM/SIGKILL from the scheduler unrelated to scientific exit, a filesystem
  error (`EIO`/`ESTALE`/`ENOSPC`) recorded in logs, or a wall-time kill.
- **Identical resubmission** is allowed at most **2** times per run card. It requires the same
  source SHA, release hashes, run card, and environment, and a new attempt suffix and scratch
  root. After a wall-time kill, one resubmission may double `--time`. That is not a scientific
  setting, and it is recorded.
- Every attempt is recorded. Failed evidence is never deleted or silently replaced.

## 22. Stop conditions (future execution)

Any of the following STOPS cluster qualification:

- the cluster source SHA differs;
- an input hash differs;
- a forbidden FP flag is present;
- `OMP_NUM_THREADS` or the other thread variables are unset or not equal to 1;
- toolchain provenance is incomplete;
- the Phase-5 cluster comparator (§12) fails materially;
- ADR-0017 cluster self-qualification fails;
- any Cstar-knot checkpoint fails;
- the §11.3/§11.4/§11.5 Mac-comparison budget is exceeded;
- the cluster R20 or R18 individual gates fail;
- a §11.7 health threshold fails;
- CQ6 repeatability is not bit-identical;
- CQ5 concurrent and serialized runs differ;
- scratch/output isolation is violated;
- any job would need a scientific setting change;
- a production run would start before platform acceptance.

No post-result tolerance change is permitted.

## 23. Cluster acceptance gate for a future BA12R

Only when all of the following hold:

- CQ-B complete and accepted, and the §2.1 ADR accepted;
- CQ0, CQ1, CQ2, CQ3, CQ4, CQ5, and CQ6 each PASS;
- the full-suite step (§13.4) is complete: every `ANALYTIC/INTERNAL` test passes, and every
  `MAC-ARTIFACT-EXACT` failure is covered by a CQ3/CQ4 tolerance comparator or has been
  reported to and dispositioned by the owner;
- the toolchain and provenance are frozen under one `CQKEY` and release;
- `docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION.md` is complete (predeclaration; environment;
  toolchain; source/input hashes; Slurm configuration; CQ0–CQ7; cross-platform comparisons;
  repeatability; concurrency; resources and timing; `CQKEY`; final disposition);
- there is no unresolved MATERIAL numerical discrepancy.

Only then may CQ7 record:

**CLUSTER PLATFORM QUALIFIED — SIX-RUN CLEAN BA12R CAMPAIGN READY FOR OWNER AUTHORIZATION.**

Mac evidence is never rewritten. Mac and cluster evidence coexist: the Mac remains the
development and reference platform, and the cluster becomes the production campaign platform
only after CQ7.

### 23.1 Mac-side steps to place and commit this preflight (owner or Mac-hosted agent)

```text
cd /Users/keeper/Documents/CompactStar/repo/CompactStar
git status --porcelain   # must be empty
git fetch origin
# require: HEAD == master == origin/master == live master == 232565a32303a4953f3f516d1d5286b6663f8f99
git worktree add -b docs/phase6-cluster-qualification-preflight \
  /Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6-cluster-qualification-preflight \
  232565a32303a4953f3f516d1d5286b6663f8f99
cp <this file> <worktree>/docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md
cd <worktree>; git diff --check; git add docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md
git diff --cached --name-status   # exactly: A docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md
git commit -m "docs: preflight eku cluster qualification"
git push -u origin docs/phase6-cluster-qualification-preflight   # non-force; do not merge
```

## 23.2 Roadmap / current-architecture banner

No `docs/architecture/CURRENT_ARCHITECTURE.md` or roadmap change is made. A PROPOSED preflight
alters no current behavior, ownership, or component boundary, so governance requires no banner.
The cluster remains **NOT QUALIFIED** in every existing status record.

## 24. Execution accounting for this preflight and disposition

| Operation | Count / result |
|---|---|
| Slurm jobs submitted (`sbatch`/`srun`/`salloc`) | 0 |
| remote builds / compilations | 0 |
| remote package or module installs | 0 |
| remote source mutation / clone / copy | 0 / 0 / 0 |
| remote Git branches / worktrees created | 0 / 0 |
| remote project directories created | **1**: `/mnt/sdd/zaki/CompactStar` (empty, `drwxr-xr-x zaki:zaki`, created 2026-09-25 09:57:19 EDT) |
| other remote writes | agent-session scratch only, under `/tmp/claude-1038/…` (a Git path listing and this draft), not under any project tree |
| ODE / BA12R / six-trajectory runs | 0 / 0 / 0 |
| production, test, CMake, baseline, EOS/data, literature changes | none |

**Disposition:** EKU CLUSTER QUALIFICATION PREFLIGHT COMPLETE — ENVIRONMENT AUTHENTICATED —
CQ0–CQ7 PLAN READY FOR OWNER ACCEPTANCE. This disposition carries the §2 blockers: the canonical
SHA `232565a…` is **not buildable on Linux**. CQ-B and the §2.1 ADR decision, plus two
comparator implementations (§10, §12), must precede CQ0. GCC acceptance of `using enum` in
C++17 mode is unverified (§13.1). Compute-node facts marked SECONDARY must
be re-authenticated live in CQ0.
