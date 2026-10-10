# EKU candidate reproduction

These are cluster-specific recipes, not relocatable release packaging. The Linux
lock is an integrity record for the reviewed installation and remains
`CANDIDATE_UNQUALIFIED`. Preserve all historical logs before rerunning a recipe:
the scripts use fixed build/log locations. Never run their compilation or numerical
steps directly on the login node. Do not rerun the full scientific suite merely
to reproduce a single refusal; select the relevant CTest name.

## Existing authenticated installation

From a login shell:

```sh
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
sbatch --export=ALL --partition=phyall --nodes=1 --ntasks=1 \
  --cpus-per-task=4 --mem=12G --time=01:00:00 \
  --chdir="$TC" --output="$TC/logs/rebuild-%j.log" \
  /mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar/tools/linux/eku_build_compactstar.sh
```

The build recipe configures Debug/Release separately, supplies all dependency and
data paths explicitly, builds with four workers and runs the 13 readiness tests.
It regenerates the local candidate lock from exact clean dependency sources;
compare its bytes with the tracked lock rather than silently accepting a changed
installation. Preserve previous package provenance and locks before any rebuild.

The recorded broad test command is in `eku_scientific_tests.sh`; submit it with
the same exported environment, four CPUs, 12 GiB and a 90 minute Slurm limit.
CTest runs serially and excludes the entire 10 Gyr controlled-evolution regression
and the three long-horizon scalar thermal regressions named in the report.
Individual scientific test timeouts are retained. Inspect its `.exit` files and
JUnit output: the recipe preserves both build-mode results and does not treat
the shell's final exit code as an aggregate scientific verdict.

The relocation adapter is `run_bounded.py`. It authenticates all transferred
scientific identities and the original six-input authority before executing an
existing bounded mode. The recorded `eku_bounded_job.sh` attempts source, control
and events. They currently refuse at the owning-star profile guard; repeating
them cannot close the scientific authority gap. Do not change certificate bytes
or disable the guard to obtain a successful run.

## Rebuilding the native dependency installation

1. Verify the supplied ZIP against its external receipt and owner-provided digest.
   Follow the supplied README and use its exact `verify_transfer.py` and
   `import_sources.py` in a new destination. Retain the promisor sparse import.
2. Read `eku_environment.sh` for exact shared module identities. Load modules on
   the login node, before Slurm export: compute nodes have no Tcl module runtime.
3. `stage_sysroot.py --output NEW_DIRECTORY` copies installed-RPM-content-verified
   headers/libraries into `NEW_DIRECTORY/sysroot` and creates its manifest. It
   refuses an existing output. Check exact RPM versions against the retained
   runtime provenance. Place the resulting sysroot at the candidate prefix before
   using the compiler wrappers; a different prefix requires a new reviewed lock.
4. Install the supplied GCC/G++ wrapper scripts as executable `bin/gcc` and
   `bin/g++` at the toolchain prefix. Their explicit sysroot, startup-file search,
   loader and RPATH are needed because compute nodes lack development packages.
5. Compute nodes also lack Git. The recorded installation copies `/usr/bin/git`
   to `toolchain/libexec/git`, copies the resolved `libpcre2-8.so.0` and
   `libcrypto.so.3` host files into `sysroot/usr/lib64` preserving loader names,
   and invokes Git with the staged loader and `--library-path`. Versions, hashes
   and wrapper bytes are retained with the evidence. This helper is used only
   for local source provenance; remote push uses the login node's Git.
6. Submit `eku_build_dependencies.sh` in a Slurm allocation (four CPUs, 8–12 GiB,
   explicit wall limit). It builds authenticated GSL 2.7.1, then exact clean Zaki
   2.0.1 and CONFIND 2.0 production sources in Debug and Release. The script
   records failures instead of converting them to passes. Submit the separate packaging
   repair recipe with the shared pkgconf 2.2 module loaded; it supplies explicit
   GSL/zlib consumer prefixes. Preserve
   initial failures and rerun outputs separately.
7. Use shared Python 3.11.9 to create the isolated venv. Install pinned NumPy
   2.2.6, SciPy 1.15.3 and mpmath 1.3.0 from the retained wheels, after checking
   their SHA-256 values in runtime provenance. No system Python is modified.
8. `record_dependencies.py --toolchain ... --sources ... --output ...` requires
   clean exact production commits, records install/compile/cache hashes and
   writes the candidate Linux lock. Read its differences before use. The existing
   Mac lock must remain unchanged. The final CompactStar build recipe performs
   this recording on a compute node before configuration.

The source/build/test recipes preserve scientific calculations and historical
floating-point protections. Rebuilding creates a new candidate installation;
matching package versions alone is insufficient to authenticate its artifacts or
claim numerical acceptance. Review the qualification report and retained failure
evidence before making any further execution decision.
