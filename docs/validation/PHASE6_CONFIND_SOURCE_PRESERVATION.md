# Phase-6 CONFIND Recovered Source Preservation

**Status:** CONFIND RECOVERED SOURCE SNAPSHOT PRESERVED — CANDIDATE AUTHORITY READY FOR OWNER REVIEW

## 1. Scope and authority

This record preserves the exact non-Git CONFIND composite source snapshot recovered by
`PHASE6_CONFIND_SOURCE_RECOVERY_SHA`:

```text
0049b9ad828ac6e517a7314ae05b07e5361b6ea2
```

The preservation task was local-Mac-only. It did not access the EKU cluster, build or test
CONFIND, Zaki, or CompactStar, change CompactStar production source, modify either external Git
repository, assign a Git identity to the recovered snapshot, or ratify ADR-0018.

## 2. Original source authentication

The exact original roots recorded by the recovery evidence were used without substitution:

- implementation/build root:
  `/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/Confind`;
- exact consumed-header root:
  `/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/dependencies/include/Confind`.

All 18 files named by the recovery manifest remained available. Their byte sizes and SHA-256
values matched the recovery record before reconstruction. No missing file was synthesized and no
alternate copy was used.

## 3. Preserved source identity

The reconstructed composite was copied byte-for-byte to:

```text
/Users/keeper/Documents/CompactStar/external/recovered/CONFIND/ed76163c22e0a1f8/source/
```

The canonical manifest is stored at:

```text
/Users/keeper/Documents/CompactStar/external/recovered/CONFIND/ed76163c22e0a1f8/CONFIND_SOURCE_MANIFEST.sha256
```

The manifest contains logical path, byte size, and SHA-256 for every source member, sorted by
logical path and terminated by a newline. Its SHA-256 is:

```text
ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083
```

Regenerating the manifest from the preserved `source/` tree produced that exact identity and an
exact byte comparison with the recovery manifest. The source tree has 18 files. Its files are
stored with mode `0444` and directories with mode `0555`; these reversible permissions are a
handling convention, while the manifest hash is the scientific identity.

The human-readable authority metadata is:

```text
/Users/keeper/Documents/CompactStar/external/recovered/CONFIND/ed76163c22e0a1f8/CONFIND_RECOVERED_SOURCE_AUTHORITY.txt
```

It records the original roots, all original file hashes, and the classification:

```text
EXACT NON-GIT SOURCE SNAPSHOT — CANDIDATE ONLY
```

No Git SHA or Git tree was assigned: none is defensibly available. The snapshot remains **NOT
ACCEPTED SOURCE AUTHORITY**, **NOT BUILD-QUALIFIED**, and **NOT LINUX-QUALIFIED**. Historical
identity confidence remains **MODERATE**; the header contract remains **EXACT_HEADER_MATCH** and
the consumed symbol contract remains **COMPLETE_SYMBOL_CONTRACT_MATCH**.

## 4. Immutable transfer ZIP

The manual-transfer artifact is:

```text
/Users/keeper/Documents/CompactStar/external/recovered/CONFIND/ed76163c22e0a1f8/CONFIND_RECOVERED_SOURCE_ed76163c22e0a1f8.zip
```

Its SHA-256 is:

```text
be41ee1c71b88627f0b11703c9d732abd7378ebbc162ee10649af771dec9ae0f
```

The ZIP contains 27 entries: 21 files and six directory entries. The files are the 18 members
under `source/`, `CONFIND_SOURCE_MANIFEST.sha256`,
`CONFIND_RECOVERED_SOURCE_AUTHORITY.txt`, and `TRANSFER_MANIFEST.txt`. The transfer manifest
SHA-256 is:

```text
639a60475fff176847d7a3f68297e911f5c8b9124475c9b81091d7f9b8562a4f
```

It lists and hashes every non-self ZIP file, states the source/target intent and scope, and records
credential exclusion. The ZIP contains no `.git` directory, credentials, keys, tokens, passwords,
browser data, cloud-provider metadata, Finder metadata, `.DS_Store`, compiled library, build
output, or vendored Darwin `libConfind.a`.

## 5. Round-trip verification

The completed ZIP was extracted into a second temporary directory. A manifest regenerated from
the extracted `source/` tree had SHA-256:

```text
ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083
```

It was byte-identical to `CONFIND_SOURCE_MANIFEST.sha256`. Every extracted payload member matched
the hash and size recorded by `TRANSFER_MANIFEST.txt`, and recursive source comparison against
the preserved tree reported zero differences. Round-trip classification: **EXACT**.

## 6. Dependency and repository guardrails

The authenticated CompactStar vendored Darwin archive remains referenced by hash only:

```text
09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2
```

It was not copied into the preservation artifact. The current divergent CONFIND Git repository
remained at:

```text
89c5d9b731534e4289d9f686549d9f0ac178e567
```

with tracked worktree state unchanged. The Zaki candidate remained at:

```text
b9ddebaded24962468954846f47238aec2726fd4
```

with tracked state unchanged. No build or test ran in either repository. CompactStar production,
tests, baselines, EOS/data, and literature were unchanged.

ADR-0018 remains **PROPOSED — OWNER RATIFICATION REQUIRED**. This preservation record neither
accepts the recovered snapshot as dependency authority nor authorizes source-build equivalence.

## 7. Next gate

Return ADR-0018 and this preserved candidate to the human owner. The next decision must:

1. ratify or revise ADR-0018;
2. accept Zaki `b9ddebaded24962468954846f47238aec2726fd4` as a source-equivalence
   candidate only;
3. accept CONFIND manifest
   `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083` as a
   source-equivalence candidate only; and
4. explicitly authorize a Mac-only source-build equivalence experiment.

No build, test, source transfer, cluster action, or subsequent qualification step begins from
this preservation record alone.
