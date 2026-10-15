# Sources and their roles

The current manuscript contains its own proofs and local computational dependencies. The following external works provide attribution and context; no private development file is required.

| Source | Role and precise comparison |
|---|---|
| Melbourne, Roysdon, Tang, Tkocz, *From simplex slicing to sharp reverse Hölder inequalities*, JLMS 113 (2026), e70461; [arXiv v2](https://arxiv.org/abs/2505.00944v2) | Section 3.1 poses the variance-normalised problem and predicts the transition near 1.68. The p≤1 result is prior work, not a new theorem of this package. |
| Eitan, *The centered convex body whose marginals have the heaviest tails*, Studia Math. 274 (2024), 201–215; [DOI](https://doi.org/10.4064/sm211027-20-11), [arXiv](https://arxiv.org/abs/2110.14382) | Established real-exponent reduction; arXiv Theorems 12/14 and proof of Theorem 4 contain relevant equality reasoning. Even higher-moment variance-normalised case already known. |
| Eskenazis, Nayar, Tkocz, *Sharp comparison of moments and the log-concave moment problem*, Adv. Math. 334 (2018), 389–416; [arXiv](https://arxiv.org/abs/1801.07597) | Crossing and sharp moment-comparison methodology. |
| Murawski, *Comparing moments of real log-concave random variables*; [arXiv](https://arxiv.org/abs/2211.05210) | General centred p/q bound and non-centred sharp comparisons; distinct from the exact centred variance-normalised optimisation here. |
| Brazitikos, Tang, Tkocz, *Moments of sums of exponentials, beyond CHS*; [arXiv](https://arxiv.org/abs/2602.03058) | Related exponential sums with different optimisation direction/range. |
| Brazitikos, Pandis, *Sharp inequalities for symmetric polynomials, Hunter's conjecture, and moments of exponential random variables*; [arXiv](https://arxiv.org/abs/2512.12254) | Related symmetric-polynomial/even-degree moment work, not a substitute for the real-exponent theorem. |
| Pólya and Szegő, *Problems and Theorems in Analysis II* (1976), Part V, Ch.1, §6 | Classical generalised-polynomial sign rule; proof reproduced in the manuscript. |
| Johansson, *Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic* (2017); *Numerical integration in arbitrary-precision ball arithmetic* (2018) | Numerical trust base and rigorous integration methods. |
| [python-flint Arb documentation](https://python-flint.readthedocs.io/en/stable/arb.html) | Ball/endpoint representation; human-readable rounded strings are not automatically directed enclosures. |
| Saumard and Wellner, *Log-concavity and strong log-concavity: A review*, Statist. Surv. 8 (2014), 45–114; [DOI](https://doi.org/10.1214/14-SS107) | Stability proof: Theorem 3.3 marginalisation, Proposition 3.6 weak closure, Theorem 2.7 existence of a density for the nondegenerate limit. |
| Cule and Samworth, *Theoretical properties of the log-concave maximum likelihood estimator of a multidimensional density*, Electron. J. Stat. 4 (2010), 254–270; [DOI](https://doi.org/10.1214/09-EJS505), [arXiv](https://arxiv.org/abs/0908.4400) | Proposition 2 upgrades weak convergence to a density to total-variation convergence in the log-concave setting. The manuscript verifies the limit is nondegenerate before using it. |

See [CITATION_AUDIT.md](CITATION_AUDIT.md) for the inspection and search ceiling. Links identify sources; they do not imply that all full texts were freshly reviewed during documentation preparation.
