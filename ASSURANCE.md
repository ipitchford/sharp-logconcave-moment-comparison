# Assurance and trust boundary

## Mathematical argument

The manuscript's universal class is all nondegenerate mean-zero real log-concave densities, including bounded and half-line supports. Reduction to the two-exponential family and equality transfer are proved in the text, following established work by Eitan and crossing-method antecedents. The large-p argument covers every real p≥7 analytically, not by extrapolation from samples.

The remaining continuum regions require the production interval computations plus the analytic enclosure lemmas and regional assembly. A passing executable alone does not prove the correctness of those formulas, the domain decomposition, the library, or the full theorem.

## Current computational status

Fresh same-code macOS replay passed all 11 production certificates in both normal and optimised Python. Both regression-suite runs passed; eight deliberately false bulk controls and eighteen malformed-input/overwrite requests were rejected. Python 3.11.14, python-flint 0.9.0, FLINT 3.6.0 (Arb integrated; separate version not exposed). This is not independent reimplementation or external peer review. Frozen-archive replay and Linux CI outcomes are recorded separately in the release checks; these local receipts do not assert them. The supplied reviewer inspected the argument and ran non-production checks but reported that production Arb replay could not be completed in that review environment. A historical coordinator replay reports success for its recorded older hashes. Neither is being relabelled as a fresh independent replay of this candidate.

The reported export defect is substantive: rounding two refined root endpoints to the same short decimal loses the enclosure; similarly, a displayed rounded endpoint need not preserve its directed-bound property. The repaired exports retain exact rational/dyadic endpoints and include endpoint round-trip containment tests. The manuscript's coarse rational crossing bracket is separate.

## Proof obligations and evidence map

| Obligation | Manuscript label | Production code | Output family |
|---|---|---|---|
| Large p≥7 strict inequality | prop:largep | code/largep/largep_certify.py | certificates/largep/largep_cert.json; analytic proof is primary |
| Endpoint derivatives L1–L3 | prop:local | code/local/local_cert.py | certificates/local/cert_L1.json, cert_L2.json, cert_L3.json |
| Strip S1/S2 around p=2 | prop:strip | code/strip/strip_cert.py | certificates/strip/S1_cert.json, S2_cert.json |
| Unique crossing | prop:crossing | code/bulk/crossing_cert.py | certificates/bulk/cert_X.json |
| Bulk strict inequalities B1–B4 | prop:bulk | code/bulk/bulk_cert.py | certificates/bulk/cert_B1.json through cert_B4.json |
| Equality in the original law class | thm:U, lem:S, thm:R2 | No standalone numerical test proves this | Reduction + strict regional inequalities |

Exact domains and signs are in [CLAIMS.json](CLAIMS.json) and the manuscript. The supplements contain enclosure derivations. Historical alternative-parameter and higher-precision outputs are not current replay evidence. The current receipt binds the eleven production outputs.

## Qualitative stability

Section `sec:stability`, Corollary `cor:distributional-stability`, derives convergence in total variation and every fixed finite Wq (q≥1) from vanishing log moment deficit, for mean-zero variance-one log-concave laws and fixed p>1, p≠2. The proof uses a uniform exponential tail, weak compactness, moment convergence, the equality classification, and existing log-concave density convergence results. It supplies a finite nondecreasing modulus tending to zero, not an explicit rate. A common modulus exists on compact subsets of one of (1,p*), (p*,2), (2,infinity). Fixed p=p* is covered with the full equality set; uniformity through p* or 2 is explicitly excluded. The method is a standard compactness consequence, not claimed as a new general theorem about log-concave convergence. No new numerical certificate is needed, but the equality theorem's existing obligations remain load-bearing.

## Exact replay binding

Normal receipt SHA-256: `85ba0c527a7d8c23e9e7e28ca7e5a1d79f8557cac83224c1dd85c912d51037b5`.

Optimised receipt SHA-256: `f38690d3650bbfb31f69590d837a83a6802924aa7c7a028cf1ae83f1a776a5f0`.

Each receipt binds source, commands, process return codes, semantic results, certificate hashes and stdout/stderr hashes. Documentation preparation checked every bound source hash against the working files and found no drift. The retained coordinator receipts are evidence records, not private scientific dependencies; the curated release must carry a portable evidence summary and its adopted exact outputs. This does not pre-approve the eventual archive.

## Non-proof diagnostics

`checks/reduction_check.py`, `checks/tsystem_check.py`, formula checks, finite-difference checks, sampled boxes and numerical scans can reveal errors but cannot establish a universal inequality or exclude every tangency. Decimal root estimates do not prove uniqueness or irrationality. Byte-identical output is an integrity property; mathematical success must be checked from exact bounds and stated obligations, allowing irrelevant runtime formatting to differ.

The strip helper `G_direct` is a diagnostic approximation, not a certified full-function enclosure; production S1/S2 use the separately bounded evaluator. Large-p sampled Part B is diagnostic; the analytic proof and Part A constants establish the unbounded range. Optional long reduction/sign-system scans were not rerun in this repair.

No formal proof assistant verification, external human peer review, independent-model agreement guarantee, or exhaustive priority verification is claimed. Final archive acceptance remains separate from documentation consistency.
