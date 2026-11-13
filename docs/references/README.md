# Tracked scientific references

Source PDFs that CompactStar documents cite, kept in Git so that every checkout
(including the EKU cluster) can read the exact bytes a citation refers to.

## What is here

| Short name | File | arXiv | Pages |
|---|---|---|---|
| FR2005 | `rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf` | astro-ph/0502116v1 | 34 |
| RFJ2006 | `rotochemical/2006-Reisenegger-Fernandez-Jofre-Internal-Heating-and-Thermal-Emission-from-Old-Neutron-Stars.pdf` | astro-ph/0610955v1 | 6 |
| BGZ2022 | `bnv/zakeri/2022-Berryman-Gardner-Zakeri-Neutron-Stars-with-Baryon-Number-Violation-Probing-Dark-Sectors.pdf` | 2201.02637v2 | 81 |
| BGZ2023 | `bnv/zakeri/2023-Berryman-Gardner-Zakeri-Macroscopic-Limits-on-Neutron-Star-Baryon-Loss.pdf` | 2305.13377v2 | 87 |
| MPR2020 | `bnv/2020-McKeen-Pospelov-Raj-Cosmological-and-Astrophysical-Probes-of-Dark-Baryons.pdf` | 2012.09865v1 | 12 |
| MPR2021 | `bnv/2021-McKeen-Pospelov-Raj-Smoke-and-Mirrors-Neutron-Star-Internal-Heating.pdf` | 2105.09951v2 | 6 |
| SNOplus2022 | `bnv/2022-SNOplus-Improved-Search-for-Invisible-Modes-of-Nucleon-Decay.pdf` | 2205.06400v2 | 7 |

[catalog.tsv](catalog.tsv) is the machine-readable index: path, SHA-256, size,
arXiv identifier and version, the file name as supplied, title, authors, and
notes. [SHA256SUMS](SHA256SUMS) lists the same hashes for `sha256sum -c`.

## Provenance

The owner supplied all seven files on 2026-11-12 in
`NSHeating_references_Nov_12_2026.zip` (SHA-256
`678ef2139779e44d345d26466a9309163ead1111e7b07c86e1815ca28ddc48fc`). They are
arXiv PDFs. The bytes are unchanged; only the file names were changed, and the
original names are in `catalog.tsv`. The archive's `__MACOSX` metadata entries
were not imported.

The owner supplied these openly available arXiv PDFs for the existing **public**
`ZAKI1905/CompactStar` repository and explicitly requested that the tracked
collection be preserved and synchronized on 2026-11-12. The earlier description
of this repository as private was incorrect. Availability in this collection
does not assign scientific authority or alter the papers' licenses.

## Relation to the shared source library

`AGENTS.md` names `/Users/keeper/Documents/CompactStar/literature` on the macOS
host as the shared source root, outside Git, with `catalog.tsv` role labels that
define source authority. That library is not available on EKU. This directory
does not replace it and does not assign authority roles.

- `FR2005` is byte-identical to the macOS-catalogue file of the same name. Its
  SHA-256 equals the value recorded in
  `docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:44`, where its
  role is primary non-superfluid formalism.
- `RFJ2006` is **not** the governed `R2006`. `R2006` is the electrostatic-
  potential formalism paper with SHA-256 `a286f15e…` (ADR-0014 §1), which
  supersedes several FR2005 equations. It is not in this directory.
- The other five files have no catalogue role recorded in this repository. A
  role is the owner's to assign.

Folder names follow the macOS library (`rotochemical/`, `bnv/zakeri/` for
owner-authored BNV papers). File names follow its pattern
`YYYY-Authors-Short-Title.pdf`.

## Using a reference

Cite the repository path and the PDF page, and name the short name, for example
`docs/references/rotochemical/2005-Fernandez-Reisenegger-Formalism-NonSuperfluid.pdf`,
PDF p. 15, eq. (68). Verify the file first:

```bash
cd docs/references && sha256sum -c SHA256SUMS
```

These PDFs are read-only references. Do not edit, re-save or re-compress them;
any change alters the hash.

## Adding a reference

1. Place the PDF in the matching topic folder with a `YYYY-Authors-Short-Title.pdf`
   name.
2. Add a row to `catalog.tsv` and a line to `SHA256SUMS`.
3. Record who supplied it and when. Leave the role unassigned unless the owner
   states it.

## Not yet here

Sources that repository documents rely on and that are still only on the macOS
host: `R2006`, `R1995`, `Y2020`, `JRF2006`, `R1997` (all listed with hashes in
`docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md:55`).

For **R2006**, the historical library filename is recorded as
`2006-Reisenegger-...-Electrostatic-Potential-Perturbations.pdf` and the expected
SHA-256 is `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8`
([preflight source table](../validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md#2-literature-authority)).
The ellipsis is part of that historical description, not an imported filename.
The 2026-11-12 housekeeping audit found no matching PDF among 43 accessible PDFs
under the EKU CompactStar root, `/home/zaki`, and `/tmp`, or inside the supplied
reference ZIP. The macOS library was not mounted; inaccessible system temporary
directories were outside the search. See the [audit](../maintenance/REPOSITORY_AUDIT_2026-11-12.md).

Import the original library bytes when available, verify the full digest, and
read their bibliographic metadata before adding a catalog row. No placeholder
PDF, invented metadata, or replacement by RFJ2006 has been added. `catalog.tsv`
and `SHA256SUMS` enumerate available files only.

For navigation and sparse-checkout instructions, see the [repository guide](../README.md).
