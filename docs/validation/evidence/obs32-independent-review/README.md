# Preserved obs-32 independent review

Imported from the actual Downloads files on 2026-10-03. The original tar.gz is
byte-unchanged, including invalid historical evidence. `validity.json` records
external hashes, reviewed commit, and explicit use/exclusion rules. All 76
internal SHA256SUMS entries were authenticated and the readable standalone
review is byte-identical to its bundled counterpart. Consult the bundled README
for detailed provenance and interpretation. No compiled executable is imported.

Preservation does not validate excluded data, ratify the amendment, or reopen the
campaign. BA12, BA12R and the obs-32 campaign failure remain unchanged.

Run1 is excluded. Run2 split-main M-* results are excluded; independent studies,
replays and unsplit results remain usable. Reject the entire run3 control
observations file as input: it is corrupted from data row 599, not just truncated.
Run3 in-memory gate/error summaries remain usable as described by its README.
Use run4 complete observations for corrected comparisons. Metadata control
prediction signs are wrong; the sign-correct analysis is run2 `analysis.json`.

The compressed bundle is the immutable source; extraction belongs in isolated
scratch output roots. Diagnostic code must be reviewed and promoted separately.
