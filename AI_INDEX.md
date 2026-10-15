# Agent reading index

## Identity and status

Anonymous candidate: sharp Lp–L2 moment comparison for centred log-concave random variables, with equality cases. Read [STATUS.md](STATUS.md) before interpreting evidence. Fresh same-code macOS replay passed all 11 production certificates in both normal and optimised Python. Both regression-suite runs passed; eight deliberately false bulk controls and eighteen malformed-input/overwrite requests were rejected. Python 3.11.14, python-flint 0.9.0, FLINT 3.6.0 (Arb integrated; separate version not exposed). This is not independent reimplementation or external peer review. Frozen-extraction replay and Linux outcomes are separately recorded release checks, not asserted by the local receipts. This is one standalone output, with no private-development dependency.

## Start here

1. [README.md](README.md): admissible laws and exponent regimes.
2. [paper.tex](paper.tex), labels `thm:main`, `thm:U`: exact theorem and equality cases.
3. [CLAIMS.json](CLAIMS.json): machine-readable claims, conditions, proof locators and limits.
4. [ASSURANCE.md](ASSURANCE.md): which obligations are analytic, interval-assisted or diagnostic.

## Dependency route

- `thm:R2` and `thm:R1`, with lettered lemmas A–G: equality-sensitive reduction and signed-moment bound. Existing Eitan methodology; no new-priority claim.
- `sec:reduction`: Zs and G(s,p), with reflection reducing s to [0,1/2].
- `prop:largep`: all real p≥7.
- `prop:local`: L1–L3 endpoint derivatives; [supplement/local.md](supplement/local.md), [code/local/local_cert.py](code/local/local_cert.py).
- `prop:strip`: S1/S2 near p=2; [supplement/strip.md](supplement/strip.md), [code/strip/strip_cert.py](code/strip/strip_cert.py).
- `prop:crossing`, `prop:bulk`: [supplement/bulk.md](supplement/bulk.md), [code/bulk/crossing_cert.py](code/bulk/crossing_cert.py), [code/bulk/bulk_cert.py](code/bulk/bulk_cert.py).
- `sec:assembly`, `lem:S`, `thm:U`: global coverage, strictness, and transfer back to arbitrary laws.
- `sec:stability`, `cor:distributional-stability`: fixed-p qualitative total-variation plus Wq convergence to the variance-one equality set, every finite q≥1. Uniformity only on compact subsets of one open exponent regime; fixed p=p* is included but no uniform modulus through p* or 2, and no explicit rate. This is a standard compactness consequence using established log-concave convergence results, not a novel general method.

[Supplement/largep.md](supplement/largep.md) accompanies the analytic tail argument. Main certificate locations are mapped in ASSURANCE. Read exact endpoint values as rationals/dyadics, not rounded diagnostic decimals.

## Replay route

Run [code/reproduce.sh](code/reproduce.sh) with an explicit Python executable and a new output directory; add `--optimized` for the second mode. Frozen-extraction replay is a separate release gate. Retain environment/build versions, source hashes, exit codes, exact bound comparisons and negative-control results. An optimised versus normal run checks assertion dependence, not independence of the mathematical method.

## Negative knowledge and limitations

- p=2: ratio one for every admissible law.
- p*: exactly two families modulo scaling/reflection, not three.
- No proof that p* is irrational.
- No solution of the general two-exponent optimisation.
- Qualitative distributional stability is stated with variance normalisation; no quantitative rate or cross-transition uniformity is claimed.
- Equality-sensitive exposition is not priority for Eitan's reduction.
- Samples and root decimal estimates are not universal proof.
- Historical review grades and internal agent agreement are not external peer review.

## Source and rights boundary

[CITATION_AUDIT.md](CITATION_AUDIT.md) separates inherited source inspection from current verification. [PROVENANCE.md](PROVENANCE.md) identifies superseded historical wording; [LICENSES.md](LICENSES.md) excludes third-party reviews, papers, software and fonts from blanket licensing. Do not cite private paths appearing in historical records as required resources.
