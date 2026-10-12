# Final correction-pass review

Performed after the bounded compute qualification, 2026-10-12. This is the
implementation agent's final source/evidence review, not a claim of a new
independent reviewer or owner acceptance. The supplied independent review is
preserved verbatim and remains the external review authority for the starting
candidate.

Reviewed the production/test diff against f4e7c18, the seven protected-source
transitions and their history, pre-scoring ADR/width rules, native receipts,
dependency overlay scope, uncertainty comparisons and final report wording.

* The FMA expression and inverse relation preserve the physical formula and
  constants; the arbitrary-precision oracle is independent of native long-double
  width. PE-V6 retains its criterion. U11 changes only synthetic boundary inputs.
* Historical manifests and scientific baselines are unchanged. The proposed
  successor retains all 33 paths, ten baselines and three special artifacts;
  41 negative controls refuse mutations. Historical FAIL and proposed candidate
  PASS are reported separately. No owner-ratified status is inferred.
* Both dependency locks retain their authenticated hashes. CONFIND/Zaki diffs are
  confined to tests/linux; installed accepted libraries are not replaced.
* All scored source hashes match source-receipt.json. Build-binding.json in the
  evidence identifies executable/archive bytes and the production build source
  commit 81a5f3a; subsequent edits affect tooling/evidence/documentation. The
  coefficient comparator additionally binds its input to the exact repeated
  artifact and refuses nonfinite inputs. No acceptance bound was enlarged after
  scoring. The final ADR wording only clarifies which predeclared section governs.
* All 31 original-review manifest entries and all 147 final evidence manifest
  entries verify. The original report hash matches the supplied hash exactly.
  No compiled executable or build tree was imported into the evidence directory.
* Actual Slurm results, including failed/retired infrastructure attempts and
  original exact-gate failures, are retained. The completed static job has exit
  0:0 and wall time 24:41. Its selected CTest failures are explicitly recorded,
  not erased by the batch job's successful completion.
* The finest coefficient artifact is exactly repeatable and agrees within the
  predeclared propagated numerical uncertainties. Native profile and table hashes
  remain attached to it. The static tail certificates are not mislabeled as the
  absent Phase-6 achieved-B 21-star frozen sensitivity certificate.

No additional production-code correction was required by this final review.
The following scientific limitations remain material and explicit:

1. Missing original 21-star producer/source/input authority blocks native frozen
   certificate construction and the requested certificate-dependent pilots.
2. Actual AppleClang confirmation and owner acceptance of the proposed amendment
   remain pending. No cluster simulation substitutes for either.
3. CMF exact golden comparisons and complete M1 stellar contour equivalence are
   unqualified across platforms; selected independent accuracy checks do not
   prove those broader claims.
4. Legacy T2 lifetime sensitivity remains unresolved, with a specific conditioning
   follow-up and evidence that it does not feed the controlled Phase-6 call path.
5. Nonlinear stability-ceiling qualification and a full-campaign cost/readiness
   conclusion are not supplied by static tests or linear diagnostics.

The deliverable is a reviewable correction candidate with reproducible bounded
evidence and a specific missing-input request. No merge or baseline promotion.
