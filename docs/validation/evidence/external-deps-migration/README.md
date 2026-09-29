# External dependency migration — early link-authority stop evidence

This is a stopped implementation pass, not a complete migration qualification.
See ../../PHASE6_EXTERNAL_DEPS_MIGRATION_STOP.md for disposition and scope.

The retained execution root is:
`/Users/keeper/Documents/CompactStar/external/qualification/compactstar-migration/e263a6e-b0cbd510`.
No build tree or binary is committed. Scripts here are exact archived execution
scripts, except this evidence README; their relative root is the execution root.
Do not execute the archived build scripts inside this repository. The driver is
copied unchanged from authorized Zaki candidate tests/fp_preservation/t2_driver.cpp.
The recorded build/test/log artifacts are outputs, never instructions.

Reproduction sequence in a fresh external execution root, with the same local
source identities and toolchain: configure/build canonical CompactStar Debug
and Release with the flags in old-*-configure.log and check_old_baselines.py;
run build_packages.py (source SHA and clean-tree checks precede each package);
export canonical CompactStar 812463a with git archive to provider-probe-source;
apply only the hash-guarded two-file plotting patch from the Zaki candidate's
remove_scratch_plots.py; run provider_probe.py, then audit_provider_probe.py.
The build scripts use local paths only and perform no source acquisition.
The old link is derived from the actual CMake smoke target link command.
The NEW diagnostic rebuilds exactly the 12 archive members selected by the
OLD TaskManager link, preserving original archive member ordering and consumer
numerical flags. Only the two separately manifested plotting deletions occur in
that disposable source export. This is not a full NEW CompactStar build or a
CMake package-discovery qualification. OpenMP remains during this staged probe.

*.map.gz are complete final executable link maps, losslessly compressed.
*-selected-disassembly.txt retains the full VecSaver function and all extant
Coord3D comparator bodies; *-selected-nm.txt records their linkage classes.
retained-large-artifact-hashes.json authenticates complete external disassembly
and nm inventories. provider-probe-audit.json records common-symbol counts,
all changed provider classes and their floating-point instruction inventories.
Counts include namespace-containing symbols, template instantiations, stubs,
GOT entries and switch tables; they are not unique function counts.

headers.sha256 and configs.sha256 use sorted prefix-relative paths, SHA256,
two spaces and LF, including the include/ or lib/cmake/ path component.
Their own hashes, installed archives, exact source SHAs and compiler are in
packages.json. Cache excerpts retain GSL and SDK zlib resolutions. Fresh Zaki
suites are 12/12 per mode, not the prior qualification's opt-in expanded 20/20.
All five OLD Debug hash-only baselines reproduce their committed bytes.
No T1/T2 execution or full CompactStar suite result is asserted by this bundle.
