# Sharp Lp–L2 comparison for centred log-concave laws

Anonymous · version 1.0.0-candidate · 15 October 2026

Reserved DOI: [10.5281/zenodo.23405173](https://doi.org/10.5281/zenodo.23405173). Intended repository: [ipitchford/sharp-logconcave-moment-comparison](https://github.com/ipitchford/sharp-logconcave-moment-comparison). Reservation and intended URLs do not assert completed publication.

This standalone package presents one evolving research output: sharp variance-normalised moment comparison for nondegenerate mean-zero real log-concave random variables, including equality cases. It combines a self-contained exposition of an established reduction with the real-exponent optimisation and its computer-assisted certificates. Earlier input versions are historical sources, not separate results that must be combined by the reader.

## Result

Let E and E' be independent standard exponential variables, L=E−E', and F(p)=E|E−1|^p. The transition p* is the unique solution in (1,2) of Gamma(p+1)^(1/p)/sqrt(2)=F(p)^(1/p). The proof uses the rational bracket 1.682457208 < p* < 1.682457209; longer decimal displays are not an irrationality claim.

| Exponent | Sharp ratio ||X||p / ||X||2 | Equality, modulo scaling and reflection |
|---|---|---|
| 1 < p < p* | Minimum Gamma(p+1)^(1/p)/sqrt(2) | Symmetric Laplace L |
| p = p* | Minimum, the two formulas coincide | Exactly two families: L and E−1 |
| p* < p < 2 | Minimum F(p)^(1/p) | Centred exponential E−1 |
| p = 2 | Ratio exactly 1 | Every admissible law |
| p > 2 | Maximum F(p)^(1/p) | Centred exponential E−1 |

The full statement and proof are in [paper.tex](paper.tex). The accompanying PDF was generated from this source and visually inspected on all 16 pages.

## Proof and evidence

Eitan already established the reduction for real exponent pairs. The manuscript supplies a self-contained, equality-sensitive proof using crossing arguments and a Chebyshev system, then optimises over Zs=s(E−1)−(1−s)(E'−1). The p≥7 comparison is analytic; the remaining regions use rigorous interval arithmetic and special treatment of endpoints, p=2 and p*. Equality follows from strict separation and the reduction. See [CLAIMS.json](CLAIMS.json) and [AI_INDEX.md](AI_INDEX.md) for precise dependencies.

The conjecture and approximate transition were previously identified by Melbourne, Roysdon, Tang and Tkocz. No first-proof claim for Eitan's reduction or new crossing method is made.

## Reproduction

Run `bash code/reproduce.sh /path/to/python /path/to/new-output` from a clean extraction, adding `--optimized` for optimised Python. Use an existing interpreter with python-flint 0.9.0; the runner does not install dependencies or access the network. Output directories must be new. It checks mathematical statuses, exact bounds and coverage rather than byte equality of timing or diagnostic fields.

Fresh same-code macOS replay passed all 11 production certificates in both normal and optimised Python. Both regression-suite runs passed; eight deliberately false bulk controls and eighteen malformed-input/overwrite requests were rejected. Python 3.11.14, python-flint 0.9.0, FLINT 3.6.0 (Arb integrated; separate version not exposed). This is not independent reimplementation or external peer review. Frozen-archive replay and Linux CI outcomes are recorded separately in the release checks; these local receipts do not assert them. [ASSURANCE.md](ASSURANCE.md) distinguishes universal proofs, interval certificates and sampled diagnostics.

## Reading and reuse

- [STATUS.md](STATUS.md): current versus historical status.
- [PROVENANCE.md](PROVENANCE.md): inputs, revision and review limitations.
- [SOURCES.md](SOURCES.md) and [CITATION_AUDIT.md](CITATION_AUDIT.md): source roles and access boundaries.
- [LICENSES.md](LICENSES.md): original prose/data CC0, original code MIT; third-party material excluded.

The manuscript also proves qualitative stability: for each fixed p>1, p≠2 and finite q≥1, a vanishing moment deficit among mean-zero variance-one log-concave laws forces convergence in total variation and Wq to the equality set. This includes fixed p=p*, with both extremal families (three oriented variance-one laws). The modulus can be chosen uniformly on a compact subset of a single open exponent regime. No explicit rate or uniform modulus through p* or 2 is asserted. This is a standard compactness consequence of the equality theorem, not a new general stability method. The general two-exponent optimisation remains unsolved. No completed public release, external peer review or final archive verification is asserted.
