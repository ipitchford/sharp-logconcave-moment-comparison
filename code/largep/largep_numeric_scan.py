#!/usr/bin/env python3
"""Non-rigorous mpmath scan (exploration record, not part of the proof).

Run:
  uv run --no-project --with mpmath python3 largep_numeric_scan.py [--out scan.json]

For p in {7, 12, 20, 50, 100, 1000}: G(s,p) - G0(p) on s = i/800, i = 1..400,
plus (G - G0)/s^2 at s = 1e-2, 1e-3, 1e-4 against the predicted limit
(1/2) d2G/ds2(0,p) = ((p-1)/(2p)) (1/F(p) - 1).
Evaluator: closed forms (incomplete gamma, 1F1), cross-checked against the SPEC
reference values and an independent mpmath quadrature; raises on mismatch.
"""
import hashlib
import json
import sys

import mpmath as mp

mp.mp.dps = 50


def Epos(a, b, p):  # E (a + b R)^p, a >= 0, b > 0
    if a == 0:
        return b ** p * mp.gamma(p + 1)
    return b ** p * mp.e ** (a / b) * mp.gammainc(p + 1, a / b, mp.inf)


def J(x, p):  # int_0^x u^p e^u du
    if x == 0:
        return mp.mpf(0)
    return x ** (p + 1) / (p + 1) * mp.hyp1f1(p + 1, p + 2, x)


def h(x, p):  # E|x - R|^p
    return mp.e ** (-x) * (mp.gamma(p + 1) + J(x, p))


def M(s, p):
    s = mp.mpf(s); p = mp.mpf(p); c = 1 - s; d = 1 - 2 * s
    if s == 0:
        return h(mp.mpf(1), p)
    return s * Epos(d, s, p) + c * c ** p * h(d / c, p)


def Mquad(s, p):
    s = mp.mpf(s); p = mp.mpf(p); d = 1 - 2 * s; c = 1 - s
    A = mp.quad(lambda r: (d + s * r) ** p * mp.e ** (-r), [0, p, 4 * p, mp.inf]) if s > 0 else d ** p
    x = d / c
    B = mp.quad(lambda r: abs(d - c * r) ** p * mp.e ** (-r), [0, x, p, 4 * p, mp.inf])
    return s * A + c * B


def G(s, p):
    s = mp.mpf(s); p = mp.mpf(p)
    return mp.log(M(s, p)) / p - mp.log(s ** 2 + (1 - s) ** 2) / 2


def main():
    out_path = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    refs = [(0.1, 1.3, '-0.213049905296183'), (0.3, 1.9, '-0.0294149150315886'),
            (0.25, 2.5, '0.143192764694943'), (0.4, 5, '0.641602497665675'),
            (0.05, 12, '1.58099627165213'), (0.5, 1.2, '-0.265784034621107'),
            (0, 3.3, '0.375483768956575')]
    for s, p, ref in refs:
        if abs(G(mp.mpf(str(s)), mp.mpf(str(p))) - mp.mpf(ref)) > mp.mpf('1e-13'):
            raise RuntimeError(f"reference mismatch at {(s, p)}")
    for s, p in [(0.05, 12), (0.3, 20), (0.45, 50)]:
        if abs(M(s, p) / Mquad(s, p) - 1) > mp.mpf('1e-30'):
            raise RuntimeError(f"quadrature cross-check failed at {(s, p)}")
    mp.mp.dps = 40
    rows = []
    for pv in [7, 12, 20, 50, 100, 1000]:
        p = mp.mpf(pv); G0 = G(0, p); F = M(0, p)
        worst, ws, n = -mp.inf, None, 400
        for i in range(1, n + 1):
            s = mp.mpf(i) / (2 * n)
            g = G(s, p) - G0
            if g > worst:
                worst, ws = g, s
        if not worst < 0:
            raise RuntimeError(f"G - G0 >= 0 found at p={pv}")
        pred = (p - 1) / (2 * p) * (1 / F - 1)
        ratios = {str(sv): mp.nstr((G(mp.mpf(sv), p) - G0) / mp.mpf(sv) ** 2, 10)
                  for sv in ['1e-2', '1e-3', '1e-4']}
        rows.append({"p": pv, "max_grid(G-G0)": mp.nstr(worst, 8), "argmax_s": mp.nstr(ws, 6),
                     "predicted_half_G_ss(0)": mp.nstr(pred, 10), "(G-G0)/s^2": ratios})
        print(rows[-1])
    with open(__file__, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    rec = {"record": "MRTT cycle42 LARGE-P numerical scan (non-rigorous, mpmath)",
           "script": "largep_numeric_scan.py", "script_sha256": sha, "mpmath": mp.__version__,
           "s_grid": "i/800, i=1..400", "rows": rows}
    if out_path:
        with open(out_path, "w") as fh:
            fh.write(json.dumps(rec, indent=2) + "\n")


if __name__ == "__main__":
    main()
