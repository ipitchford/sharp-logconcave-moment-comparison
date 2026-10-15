#!/usr/bin/env python3
"""
check_formulas.py -- independent (non-rigorous) cross-check of the closed forms used by
bulk_cert.py against mpmath quadrature of the defining integrals.

  M(s,p)   = s int_0^inf (d+sr)^p e^{-r} dr + c int_0^inf |d-cr|^p e^{-r} dr
  M_p      = same with an extra factor log|.|            (quadrature)
  M_s      = mpmath numerical derivative (mp.diff) of the quadrature M in s
  SPEC M_s = p( s E[(d+sR)^{p-1}(R-1)] + c E[|d-cR|^{p-2}(d-cR)(R-1)] )   (quadrature)

Run: uv run --no-project --with python-flint==0.9.0 --with mpmath python3 check_formulas.py
"""
import json, random, sys
import mpmath as mp
from flint import arb, ctx, fmpq
import bulk_cert as bc

mp.mp.dps = 30


def quadM(s, p, logw=False):
    s = mp.mpf(s); p = mp.mpf(p); d = 1 - 2 * s; c = 1 - s
    w = (lambda y: mp.log(y)) if logw else (lambda y: 1)
    A = mp.quad(lambda r: (d + s * r) ** p * w(d + s * r) * mp.e ** (-r), [0, 1, 10, mp.inf]) if s > 0 else 0
    x = d / c
    pts = [0, x, x + 1, 10, mp.inf] if x > 0 else [0, 1, 10, mp.inf]
    B = mp.quad(lambda r: abs(d - c * r) ** p * w(abs(d - c * r)) * mp.e ** (-r), pts)
    return s * A + c * B


def spec_Ms(s, p):
    s = mp.mpf(s); p = mp.mpf(p); d = 1 - 2 * s; c = 1 - s; x = d / c
    A = mp.quad(lambda r: (d + s * r) ** (p - 1) * (r - 1) * mp.e ** (-r), [0, 1, 10, mp.inf]) if s > 0 else 0
    pts = [0, x, x + 1, 10, mp.inf] if x > 0 else [0, 1, 10, mp.inf]
    B = mp.quad(lambda r: abs(d - c * r) ** (p - 2) * (d - c * r) * (r - 1) * mp.e ** (-r), pts)
    return p * (s * A + c * B)


def main():
    ctx.prec = 128
    random.seed(20261014)
    pts = [(fmpq(1, 10), fmpq(13, 10)), (fmpq(1, 20), fmpq(39, 20)), (fmpq(9, 20), fmpq(41, 20)),
           (fmpq(1, 2), fmpq(17, 10)), (fmpq(1, 64), fmpq(3, 2)), (fmpq(3, 10), fmpq(7)), (fmpq(1, 5), fmpq(1))]
    for _ in range(8):
        pts.append((fmpq(random.randint(1, 500), 1000), fmpq(random.randint(1000, 7000), 1000)))
    rows = []; worst = 0.0
    for s, p in pts:
        sf = mp.mpf(int(s.p)) / int(s.q); pf = mp.mpf(int(p.p)) / int(p.q)
        S, P = arb(s), arb(p)
        M_arb = bc.M_pt(S, P)
        Ms_arb = bc.Ms_box(s, s, P)
        Mp_arb = bc.Mp_pt(s, p)
        M_q = quadM(sf, pf)
        Mp_q = quadM(sf, pf, True)
        Ms_q = mp.diff(lambda t: quadM(t, pf), sf) if sf < 0.5 else mp.mpf(0)
        Ms_spec = spec_Ms(sf, pf) if pf > 1 else mp.nan
        def rel(a, b):
            a = mp.mpf(a.mid().str(25, radius=False)) if hasattr(a, "mid") else a
            return float(abs(a - b) / max(abs(b), mp.mpf(1)))
        e = [rel(M_arb, M_q), rel(Ms_arb, Ms_q), rel(Mp_arb, Mp_q)]
        if pf > 1:
            e.append(rel(Ms_arb, Ms_spec))
        worst = max(worst, max(e))
        rows.append({"s": str(s), "p": str(p), "M": M_arb.str(15, radius=False), "M_quad": mp.nstr(M_q, 15),
                     "Ms": Ms_arb.str(15, radius=False), "Ms_numdiff": mp.nstr(Ms_q, 15),
                     "Ms_SPEC_formula": mp.nstr(Ms_spec, 15), "Mp": Mp_arb.str(15, radius=False),
                     "Mp_quad": mp.nstr(Mp_q, 15), "max_rel_err": max(e)})
    out = {"check": "closed forms vs mpmath quadrature / numerical differentiation", "points": rows,
           "worst_relative_discrepancy": worst, "ok": worst < 1e-12}
    print(json.dumps(out, indent=1))
    sys.exit(0 if out["ok"] else 1)


if __name__ == "__main__":
    main()
