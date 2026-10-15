# Current status

Date: 15 October 2026. Scholarly creator: Anonymous. Status: unrefereed computer-assisted theorem candidate, version 1.0.0-candidate.

The current manuscript states the sharp moment theorem and full equality classification for p>1. This is one research output, not an aggregation requiring access to an earlier release. Current documentation supersedes overbroad descriptive claims in historical records without modifying those records.

## Completed current-source checks

Fresh same-code macOS replay passed all 11 production certificates in both normal and optimised Python. Both regression-suite runs passed; eight deliberately false bulk controls and eighteen malformed-input/overwrite requests were rejected. Python 3.11.14, python-flint 0.9.0, FLINT 3.6.0 (Arb integrated; separate version not exposed). This is not independent reimplementation or external peer review. Frozen-archive replay and Linux CI outcomes are recorded separately in the release checks; these local receipts do not assert them.

Exact endpoint export and containment regressions passed. Receipt/source bindings were read back with no source drift. See ASSURANCE.md for the exact receipt digests.

## Archive and release evidence

- The eleven conventional certificate files are regenerated exact normal-mode outputs, not the older rounded exports.
- The PDF was rebuilt from the integrated source and all 16 pages inspected. The archive manifest binds the final payload.
- Clean frozen-extraction and Linux CI outcomes are separate release gates, inspectable in the release and GitHub Actions records; a configured workflow alone is not evidence of success.
- Archive identity: DOI 10.5281/zenodo.23405173; repository ipitchford/sharp-logconcave-moment-comparison; version 1.0.0-candidate. A reserved DOI is not proof that the record is public; provider state and exact public asset checks establish publication separately.

Historical working logs and reviewer material are excluded from this curated package. Original input archives remain unchanged; their PASS labels are not reused as proof for revised bytes. Public portable replay summaries bind the current executed source and outputs.

## Claim ceiling

Exactly two extremal families at p*, modulo nonzero scaling and reflection. Every admissible law is extremal at p=2. No assertion that p* is irrational. No priority claim for the real-exponent reduction or crossing method. No inference of external peer review from internal model reports or of theorem correctness from a REF-style score.

The integrated corollary `cor:distributional-stability` states fixed-p qualitative total-variation and Wq stability for each finite q≥1, on mean-zero variance-one log-concave laws. It includes p=p* and excludes p=2. Compact subsets of a single open regime admit a common modulus; there is no explicit rate or uniformity through p* or 2. This analytic consequence still relies on the equality theorem and its computational obligations; it does not replace frozen-archive or platform verification.
