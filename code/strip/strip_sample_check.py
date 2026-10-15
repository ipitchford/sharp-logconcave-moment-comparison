#!/usr/bin/env python3
"""
strip_sample_check.py -- independent (non-rigorous) sanity check of STRIP certificates.

For every leaf box of an S1/S2 certificate, sample points (s,p) inside the box, evaluate the
certified quantity with mpmath quadrature of the derivative formulas (code path independent of
the Arb certifier: no closed forms, no Taylor expansion in p), and check that the value lies in
the certified enclosure [lower, upper] of that leaf.

  S1: H(s,p)  = G_p(s,p) - G0'(p)
  S2: G_pss(s,p) = d/dp d^2/ds^2 G(s,p)
Usage: strip_sample_check.py CERT.json [--per-leaf 2] [--seed 1]
"""
import argparse
import json
import random
import sys
from fractions import Fraction as Fr

import mpmath as mp

mp.mp.dps = 25


def E(s, g):
    s = mp.mpf(s); d = 1 - 2*s; c = 1 - s
    pos = 0
    if s > 0:
        pos = s*mp.quad(lambda r: g(d + s*r, r)*mp.e**(-r), [0, 1, 5, 20, mp.inf])
    r0 = d/c
    pts = [0, r0, r0 + 1, 20, mp.inf] if r0 > 0 else [0, 1, 20, mp.inf]
    neg = c*mp.quad(lambda r: g(d - c*r, r)*mp.e**(-r), pts)
    return pos + neg


def quantities(s, p):
    s = mp.mpf(s); p = mp.mpf(p)
    a = abs
    M = E(s, lambda z, r: a(z)**p)
    Mp = E(s, lambda z, r: a(z)**p*mp.log(a(z)))
    out = dict(Gp=-mp.log(M)/p**2 + Mp/(p*M))
    return M, Mp, out


def H(s, p):
    _, _, q1 = quantities(s, p)
    _, _, q0 = quantities(0, p)
    return q1['Gp'] - q0['Gp']


def Gpss(s, p):
    s = mp.mpf(s); p = mp.mpf(p)
    a = abs
    M = E(s, lambda z, r: a(z)**p)
    Mp = E(s, lambda z, r: a(z)**p*mp.log(a(z)))
    A1 = E(s, lambda z, r: a(z)**(p-2)*z*(r-1))
    A1L = E(s, lambda z, r: a(z)**(p-2)*z*(r-1)*mp.log(a(z)))
    B1 = E(s, lambda z, r: a(z)**(p-2)*((r-1)**2 + 1))
    B1L = E(s, lambda z, r: a(z)**(p-2)*((r-1)**2 + 1)*mp.log(a(z)))
    Ms = p*A1; Msp = A1 + p*A1L
    Mss = p*(p-1)*B1; Mssp = (2*p-1)*B1 + p*(p-1)*B1L
    return (-(Mss/M - (Ms/M)**2)/p**2
            + (Mssp/M - Mss*Mp/M**2 - 2*(Ms/M)*(Msp/M - Ms*Mp/M**2))/p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cert')
    ap.add_argument('--per-leaf', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    cert = json.load(open(args.cert))
    lemma = cert['lemma']
    f = H if lemma == 'S1' else Gpss
    rng = random.Random(args.seed)
    bad = 0
    n = 0
    worst_rel_pos = 0.0
    for leaf in cert['leaves']:
        slo, shi = Fr(leaf['s'][0]), Fr(leaf['s'][1])
        qlo, qhi = Fr(leaf['q'][0]), Fr(leaf['q'][1])
        def exact_mpf(text):
            q = Fr(text)
            return mp.mpf(q.numerator) / q.denominator
        lo, hi = exact_mpf(leaf['lower']), exact_mpf(leaf['upper'])
        tol = mp.mpf('1e-11')*(abs(lo) + abs(hi)) + mp.mpf('1e-15')
        for _ in range(args.per_leaf):
            s = slo + (shi - slo)*Fr(rng.randrange(0, 10**6 + 1), 10**6)
            q = qlo + (qhi - qlo)*Fr(rng.randrange(0, 10**6 + 1), 10**6)
            v = f(mp.mpf(s.numerator)/s.denominator, 2 + mp.mpf(q.numerator)/q.denominator)
            n += 1
            if not (lo - tol <= v <= hi + tol):
                bad += 1
                print(f"OUTSIDE: s={float(s)} p={2+float(q)} value={mp.nstr(v, 15)} "
                      f"enclosure=[{leaf['lower']}, {leaf['upper']}]", flush=True)
            # how close is the sample to the upper bound, relative to the enclosure width
            if hi > lo:
                worst_rel_pos = max(worst_rel_pos, float((v - lo)/(hi - lo)))
    print(json.dumps(dict(cert=args.cert, lemma=lemma, samples=n, outside=bad,
                          max_relative_position_in_enclosure=round(worst_rel_pos, 6),
                          result="PASS" if bad == 0 else "FAIL")))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
