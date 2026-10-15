#!/usr/bin/env python3
"""LOCAL worker certifier (cycle 42, MRTT L_p-L_2 conjecture).

Certifies the sign of d2G/ds2 (with a quantitative margin) on the regions

  L1:  s in [9/20, 1/2],  p in (1, 7/4]       d2G/ds2 >  thr   (production thr = 6/125)
  L2:  s in [0, 1/20],    p in [8/5, 39/20]   d2G/ds2 >  thr   (production thr = 29/2000)
  L3:  s in [0, 1/20],    p in [41/20, 7]     d2G/ds2 < -thr   (production thr = 31/2000)

(each contains the SPEC region: p*+1/20 < 7/4, p*-1/20 > 8/5; L3 upper end 7 per coordinator update).

Formula (local.md, Lemma 2):  with c = 1-s, d = 1-2s, N = s^2+c^2, x = d/c, y = d/s,
   G_ss = (p-1) b / M - p A^2 + 2 (d/N) A - 1/N^2,   A = G_s = (s c / N) a / M,
where M, a, (p-1) b are explicit combinations of Gamma(p), Gamma(p+1), upper incomplete gamma
functions, 1F1 and Tricomi U functions (local.md, Lemma 3).

Rigour: python-flint 0.9.0 Arb ball arithmetic (outward rounding).  Region endpoints and box
endpoints are exact rationals (fractions.Fraction); every point value enters Arb through
arb(fmpq(.)).  The p-dependence is handled by a single Arb ball containing [p0,p1]; the
s-dependence by exact monotone hulls (endpoint evaluations) of the special-function pieces and by
independent tight interval enclosures of the rational quantities s, c, d, N, x, y.  A box is
accepted only if the Arb ball lies strictly beyond the threshold (sign*G_ss > thr); undecided comparisons are
failures.  Failures raise / exit non-zero.

Usage:
  uv run --no-project --with python-flint==0.9.0 python3 local_cert.py selftest
  uv run --no-project --with python-flint==0.9.0 python3 local_cert.py L1|L2|L3 [--threshold Q] [--out FILE]
  (the default threshold in REGIONS is a weaker sign certificate; the production certificates
   pass --threshold 6/125, 29/2000, 31/2000 respectively; exit status 0 iff PASS)
"""
import argparse
import functools
import random
import hashlib
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import flint
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from certificate_io import point, interval, positive_int, rational, runtime
from flint import arb, fmpq, ctx

PREC = 128
ctx.prec = PREC


class CertError(Exception):
    pass


# ----------------------------------------------------------------------------- helpers
def A(q):
    """exact rational -> Arb ball containing it (outward rounded)."""
    q = Fr(q)
    return arb(fmpq(q.numerator, q.denominator))


def hull(*xs):
    r = xs[0]
    for t in xs[1:]:
        r = r.union(t)
    return r


def fin(x, what):
    if not x.is_finite():
        raise CertError("non-finite enclosure for %s: %s" % (what, x))
    return x


def sqr(x):
    """enclosure of x^2 that is nonnegative."""
    if x.contains(0):
        m = abs(x).upper()
        return arb(0).union(m * m)
    return x * x


def x_of(s):
    return (1 - 2 * s) / (1 - s)


def y_of(s):
    return (1 - 2 * s) / s


def N_of(s):
    return s * s + (1 - s) ** 2


def J_pt(a, xF):
    """int_0^x u^a e^u du at an exact point x >= 0; a an Arb ball with a > -1 (here a >= p-1 > 0)."""
    if xF == 0:
        return arb(0)
    x = A(xF)
    return x ** (a + 1) / (a + 1) * x.hypgeom_1f1(a + 1, a + 2)


def J_hull(a, x0, x1):
    """J_a is increasing in x >= 0: hull of the endpoint values encloses J_a([x0,x1])."""
    return hull(J_pt(a, x0), J_pt(a, x1))


def E_pt(a, yF):
    """e^y Gamma(a, y) = int_0^inf (y+t)^(a-1) e^-t dt at exact point y >= 0; a ball, a > 0."""
    if yF == 0:
        return a.gamma()
    y = A(yF)
    return y.exp() * y.gamma_upper(a)


def E_hull(a, y0, y1):
    """e^y Gamma(a,y) is monotone in y (direction sign(a-1)) for each fixed a: hull of endpoints."""
    return hull(E_pt(a, y0), E_pt(a, y1))


def K_pt(b, zF):
    """1F1(1; b; z) = sum_k z^k / (b)_k at an exact point z with |z| <= 8, for an Arb ball b > 1/2.
    (The ball for b = p may stick out below 1 by a rounding radius, so only b > 1/2 is assumed.)
    For b > 1/2 and k >= 1, (b)_k = b (b+1)...(b+k-1) >= (k-1)!/2, so the tail after K terms is
    sum_{k>=K} |z|^k / (b)_k <= 2 sum_{k>=K} |z|^k/(k-1)! <= (2 |z|^K / (K-1)!) K / (K - |z|);
    it is added as an explicit error radius."""
    zF = Fr(zF)
    if abs(zF) > 8:
        raise CertError("K_pt needs |z| <= 8")
    if not (b > 0.5):
        raise CertError("K_pt needs b > 1/2")
    K = 48 + 8 * math.ceil(abs(zF))
    z = A(zF)
    t = arb(1)
    tot = arb(1)
    for k in range(1, K):
        t = t * z / (b + (k - 1))
        tot += t
    az = A(abs(zF))
    tail = 2 * az ** K / math.factorial(K - 1) * K / (K - az)
    return tot + arb(0, tail.upper())


def K_hull(b, z0, z1):
    """1F1(1; b; z) for z in [z0, z1] subset [-8, 8], b >= 1: increasing in z, because
    d/dz 1F1(1;b;z) = (1/b) 1F1(2;b+1;z) > 0 (positive Euler integral for b > 1, e^z for b = 1)."""
    return hull(K_pt(b, z0), K_pt(b, z1))


def powbox(b0, b1, e0, e1):
    """Enclosure of b^e for b in [b0,b1] in [0,1], e in [e0,e1] in [0,inf), e > 0 where b = 0.
    b^e is increasing in b and nonincreasing in e, so the range is [b0^e1, b1^e0]."""
    if not (0 <= b0 <= b1 <= 1 and 0 <= e0 <= e1):
        raise CertError("powbox domain %s %s %s %s" % (b0, b1, e0, e1))
    lo = arb(0) if b0 == 0 else A(b0) ** A(e1)
    hi = arb(1) if e0 == 0 else A(b1) ** A(e0)
    return lo.union(hi)


@functools.lru_cache(maxsize=400000)
def P_pt(p0, p1, shift, j, dF, sF):
    """E[(d+sR)^q R^j] for all q in [p0+shift, p1+shift], at the exact point (d,s), d > 0, s >= 0.
    Tricomi-U form E[(d+sR)^q R^j] = d^q y^(j+1) j! U(j+1, q+j+2, y), y = d/s; s = 0 gives d^q j!.
    (Memoised on exact rational arguments; the Arb ball for q is rebuilt identically each time.)"""
    q = hull(A(p0), A(p1)) + shift
    d = A(dF)
    if sF == 0:
        return d ** q * math.factorial(j)
    y = A(dF / sF)
    return d ** q * y ** (j + 1) * math.factorial(j) * y.hypgeom_u(j + 1, q + j + 2)


def P_hull(p0, p1, shift, j, d0, d1, s0, s1):
    """E[(d+sR)^q R^j] is monotone in d and in s separately (R >= 0, t -> t^q monotone),
    for each fixed q.  Hence over [d0,d1] x [s0,s1] it lies between the four corner values."""
    return hull(P_pt(p0, p1, shift, j, d0, s0), P_pt(p0, p1, shift, j, d0, s1),
                P_pt(p0, p1, shift, j, d1, s0), P_pt(p0, p1, shift, j, d1, s1))


# ----------------------------------------------------------------------------- core enclosure
def parts_box(s0, s1, p0, p1, mode):
    """Enclosures of M, a, (p-1)b, A, G_ss valid for all (s,p) in [s0,s1] x [p0,p1], p > 1.

    mode 'L1'      : plus side in Gamma(.,y) form with the d^(p-1) singular part cancelled
                     (requires s0 >= 1/4).
    mode 'small_s' : plus side in Tricomi-U form via monotone corner hulls (requires s1 <= 1/4).
    """
    s0, s1, p0, p1 = map(rational, (s0,s1,p0,p1))
    if not (0 <= s0 <= s1 <= Fr(1, 2)):
        raise CertError("bad s range")
    if not (1 <= p0 <= p1):
        raise CertError("bad p range")
    if mode == 'L1' and s0 < Fr(1, 4):
        raise CertError("L1 mode needs s0 >= 1/4")
    if mode == 'small_s' and s1 > Fr(1, 4):
        raise CertError("small_s mode needs s1 <= 1/4")
    if mode == 'small_s' and p0 <= 1:
        raise CertError("uncancelled small_s formula requires p > 1")
    P = hull(A(p0), A(p1))
    Pm1 = hull(A(p0 - 1), A(p1 - 1))
    S = hull(A(s0), A(s1))
    C = hull(A(1 - s1), A(1 - s0))
    Dd = hull(A(1 - 2 * s1), A(1 - 2 * s0))
    N = hull(A(N_of(s1)), A(N_of(s0)))          # N decreasing on [0, 1/2]
    x0, x1 = x_of(s1), x_of(s0)                   # x decreasing in s
    X = hull(A(x0), A(x1))
    ex = (-X).exp()
    gp = P.gamma()
    gp1 = (P + 1).gamma()
    cp1 = C ** Pm1                                # c^(p-1)
    N_C = N / C
    N_CC = N / (C * C)
    SS_CN = S * S / (C * N)
    SS_NN = S * S / (N * N)
    # ---- minus side (Z = d - cR), u = x - R inside, t = R - x outside
    Jp = J_hull(P, x0, x1)
    Jpm1 = J_hull(Pm1, x0, x1)
    M_minus = cp1 * C * C * ex * (Jp + gp1)
    a_minus = -cp1 * ex * (N_C * Jpm1 + Jp + gp1 - N_C * gp)
    b_minus = cp1 * ex * (Pm1 * 2 * SS_CN * Jpm1 + Pm1 * SS_NN * Jp + N_CC * gp
                          - Pm1 * 2 * SS_CN * gp + Pm1 * SS_NN * gp1)
    Dp = powbox(1 - 2 * s1, 1 - 2 * s0, p0 - 1, p1 - 1)   # d^(p-1)
    K1 = K_hull(P, -x1, -x0)                               # 1F1(1;p;-x)
    if mode == 'L1':
        y0, y1 = y_of(s1), y_of(s0)                        # y decreasing in s
        Kp = K_hull(P, y0, y1)                             # 1F1(1;p;y)
        sp1 = S ** Pm1
        Ep1 = E_hull(P + 1, y0, y1)
        Ep = E_hull(P, y0, y1)
        ey = hull(A(y0).exp(), A(y1).exp())
        N_S = N / S
        CC_NN = C * C / (N * N)
        CC_SN = C * C / (S * N)
        M_plus = sp1 * S * S * Ep1
        a_plus = sp1 * (Ep1 - N_S * Ep)
        b_plus = sp1 * (Pm1 * CC_NN * Ep1 - Pm1 * 2 * CC_SN * Ep + (N / (S * S)) * ey * gp)
        sing = Dp * (N_CC * K1 - (N / (S * S)) * Kp)
    elif mode == 'small_s':
        d0, d1 = 1 - 2 * s1, 1 - 2 * s0
        Pp_0 = P_hull(p0, p1, 0, 0, d0, d1, s0, s1)
        Pa_1 = P_hull(p0, p1, -1, 1, d0, d1, s0, s1)
        Pa_0 = P_hull(p0, p1, -1, 0, d0, d1, s0, s1)
        Pb_2 = P_hull(p0, p1, -2, 2, d0, d1, s0, s1)
        Pb_1 = P_hull(p0, p1, -2, 1, d0, d1, s0, s1)
        Pb_0 = P_hull(p0, p1, -2, 0, d0, d1, s0, s1)
        CC_NN = C * C / (N * N)
        M_plus = S * Pp_0
        a_plus = Pa_1 - 2 * S * Pa_0
        b_plus = Pm1 * S * (CC_NN * Pb_2 - 4 * S * CC_NN * Pb_1 + (4 * S * S * CC_NN + 1) * Pb_0)
        sing = Dp * N_CC * K1
    else:
        raise CertError("unknown mode")
    M = fin(M_plus + M_minus, "M")
    if not (M > 0):
        raise CertError("M not certified positive")
    a = fin(a_plus + a_minus, "a")
    pb = fin(b_plus + b_minus + sing, "(p-1)b")
    Aq = fin((S * C / N) * a / M, "A")
    G = fin(pb / M - P * sqr(Aq) + 2 * (Dd / N) * Aq - 1 / (N * N), "G_ss")
    return dict(M=M, a=a, pb=pb, A=Aq, Gss=G, N=N)


def gss_box(s0, s1, p0, p1, mode):
    return parts_box(s0, s1, p0, p1, mode)['Gss']


# ----------------------------------------------------------------------------- self tests
def selftest():
    rep = {}
    # (1) SPEC reference values of G (mpmath, 40 digits; quoted to 15 significant digits)
    refs = [("1/10", "13/10", "-0.213049905296183"), ("3/10", "19/10", "-0.0294149150315886"),
            ("1/4", "5/2", "0.143192764694943"), ("2/5", "5", "0.641602497665675"),
            ("1/20", "12", "1.58099627165213"), ("1/2", "6/5", "-0.265784034621107"),
            ("0", "33/10", "0.375483768956575")]
    tol = Fr(1, 10 ** 14)
    rows = []
    for s, p, g in refs:
        sF, pF = Fr(s), Fr(p)
        modes = ['small_s'] if sF <= Fr(1, 4) else ['L1']
        if Fr(1, 4) <= sF <= Fr(1, 4):
            modes = ['small_s', 'L1']
        for mode in modes:
            r = parts_box(sF, sF, pF, pF, mode)
            G = r['M'].log() / A(pF) - r['N'].log() / 2
            err = abs(G - arb(g))
            ok = err.upper() < A(tol)
            rows.append(dict(s=s, p=p, mode=mode, G=G.str(16), ref=g, ok=bool(ok)))
            if not ok:
                raise CertError("reference value mismatch at %s %s: %s vs %s" % (s, p, G, g))
    rep["G_reference_values"] = rows
    # G0(13/10), Gh(13/10)
    p = A(Fr(13, 10))
    F = (-arb(1)).exp() * ((p + 1).gamma() + J_pt(p, Fr(1)))
    G0 = F.log() / p
    Gh = (p + 1).gamma().log() / p - arb(2).log() / 2
    for val, g in ((G0, "-0.213229500079513"), (Gh, "-0.227966317234103")):
        if not (abs(val - arb(g)).upper() < A(tol)):
            raise CertError("G0/Gh reference mismatch")
    rep["G0_Gh_13_10"] = [G0.str(16), Gh.str(16)]
    # (2) s = 0: G_ss = ((p-1)/p)(1/F - 1);  (3) s = 1/2: G_ss = 4 (p-2)^2 / p
    rows = []
    for pq in ["8/5", "19/10", "39/20", "41/20", "5/2", "7", "12"]:
        pF = Fr(pq); p = A(pF)
        F = (-arb(1)).exp() * ((p + 1).gamma() + J_pt(p, Fr(1)))
        pred = (p - 1) / p * (1 / F - 1)
        r = parts_box(0, 0, pF, pF, 'small_s')
        ok = r['Gss'].overlaps(pred) and r['A'].contains(0) and abs(r['A']).upper() == 0
        rows.append(dict(s="0", p=pq, Gss=r['Gss'].str(16), pred=pred.str(16), A=r['A'].str(5), ok=bool(ok)))
        if not ok:
            raise CertError("s=0 identity failed at p=%s" % pq)
    for pq in ["1", "21/20", "3/2", "7/4", "5/2"]:
        pF = Fr(pq); p = A(pF)
        pred = 4 * (p - 2) ** 2 / p
        r = parts_box(Fr(1, 2), Fr(1, 2), pF, pF, 'L1')
        ok = r['Gss'].overlaps(pred) and r['A'].overlaps(arb(0))
        rows.append(dict(s="1/2", p=pq, Gss=r['Gss'].str(16), pred=pred.str(16), A=r['A'].str(5), ok=bool(ok)))
        if not ok:
            raise CertError("s=1/2 identity failed at p=%s" % pq)
    rep["endpoint_identities"] = rows
    # (4) the two plus-side forms agree where both are usable (s = 1/4)
    rows = []
    for pq in ["11/10", "19/10", "5/2", "6"]:
        pF = Fr(pq)
        g1 = parts_box(Fr(1, 4), Fr(1, 4), pF, pF, 'small_s')['Gss']
        g2 = parts_box(Fr(1, 4), Fr(1, 4), pF, pF, 'L1')['Gss']
        ok = g1.overlaps(g2) and (g1 - g2).rad() < 1e-20
        rows.append(dict(s="1/4", p=pq, U_form=g1.str(16), Gamma_form=g2.str(16), ok=bool(ok)))
        if not ok:
            raise CertError("U-form / Gamma-form mismatch at p=%s" % pq)
    rep["U_vs_Gamma_form"] = rows
    # (5) consistency of the box enclosures: point values at random interior points must lie in
    #     the box enclosure (tests the monotone-hull logic; seed fixed for reproducibility)
    rng = random.Random(20261014)
    rows = []
    for name in ("L1", "L2", "L3"):
        S0, S1, P0, P1, mode = REGIONS[name][:5]
        nin = 0
        for _ in range(40):
            ws = Fr(rng.choice([1, 2, 5, 10, 50]), 1000)
            wpp = Fr(rng.choice([1, 2, 5, 10, 50]), 1000)
            a0 = S0 + (S1 - S0 - ws) * Fr(rng.randrange(1001), 1000) if S1 - S0 > ws else S0
            b0 = P0 + (P1 - P0 - wpp) * Fr(rng.randrange(1001), 1000)
            a1, b1 = min(a0 + ws, S1), b0 + wpp
            if rng.random() < 0.2:       # also exercise the region edges s = s0, s1, p = p0
                a0, a1 = (S0, S0 + ws) if rng.random() < 0.5 else (S1 - ws, S1)
            box = parts_box(a0, a1, b0, b1, mode)
            for _k in range(3):
                sp = a0 + (a1 - a0) * Fr(rng.randrange(1001), 1000)
                pp = b0 + (b1 - b0) * Fr(rng.randrange(1001), 1000)
                if pp == 1:
                    pp = b0 + (b1 - b0) / 2
                pt = parts_box(sp, sp, pp, pp, mode)
                for key in ("M", "a", "pb", "A", "Gss"):
                    if not box[key].overlaps(pt[key]):
                        raise CertError("box enclosure misses point value (%s) at s=%s p=%s" % (key, sp, pp))
                    if box[key].contains(pt[key]):
                        nin += 1
        rows.append(dict(region=name, boxes=40, points=120, quantities_contained=nin, of=600))
        if nin < 600:
            raise CertError("box enclosure does not strictly contain a point enclosure in %s" % name)
    rep["random_box_containment"] = rows
    # (6) the four corners of every region (boxes of width 1/1000 touching the region corners,
    #     in particular the singular corner (s,p) = (1/2,1) of L1): finite enclosures that
    #     contain the point values at the box centres
    rows = []
    h = Fr(1, 1000)
    for name in ("L1", "L2", "L3"):
        S0, S1, P0, P1, mode = REGIONS[name][:5]
        for a0 in (S0, S1 - h):
            for b0 in (P0, P1 - h):
                box = parts_box(a0, a0 + h, b0, b0 + h, mode)
                pt = parts_box(a0 + h / 2, a0 + h / 2, b0 + h / 2, b0 + h / 2, mode)
                for key in ("M", "a", "pb", "A", "Gss"):
                    if not box[key].contains(pt[key]):
                        raise CertError("corner box check failed (%s) %s %s %s" % (name, key, a0, b0))
                rows.append(dict(region=name, s0=str(a0), p0=str(b0), Gss_box=box["Gss"].str(6),
                                 Gss_centre=pt["Gss"].str(10)))
    rep["region_corner_boxes"] = rows
    return rep


# ----------------------------------------------------------------------------- branch and bound
REGIONS = {
    # name: (s0, s1, p0, p1, mode, sign, threshold, initial s-cuts, initial p-cuts, weight_s, weight_p)
    # claim: sign * G_ss > threshold on the whole region
    'L1': (Fr(9, 20), Fr(1, 2), Fr(1), Fr(7, 4), 'L1', +1, Fr(1, 40), 5, 30, 2, 1),
    'L2': (Fr(0), Fr(1, 20), Fr(8, 5), Fr(39, 20), 'small_s', +1, Fr(1, 150), 5, 35, 4, 1),
    'L3': (Fr(0), Fr(1, 20), Fr(41, 20), Fr(7), 'small_s', -1, Fr(1, 150), 5, 100, 4, 1),
}
MIN_WIDTH = Fr(1, 2 ** 22)


def certify(name, max_boxes=10 ** 7, verbose=True, threshold=None):
    runtime()
    positive_int(max_boxes, "max_boxes")
    s0, s1, p0, p1, mode, sign, thr, ns, np_, ws, wp = REGIONS[name]
    if not (0 <= s0 < s1 <= Fr(1,2) and 1 <= p0 < p1) or isinstance(sign,bool) or sign not in (-1,1):
        raise CertError("invalid certification region")
    positive_int(ns, 's grid'); positive_int(np_, 'p grid')
    if threshold is not None:
        thr = rational(threshold)
    if thr < 0:
        raise CertError("threshold must be >= 0")
    THR = A(thr)
    t0 = time.time()
    stack = []
    for i in range(ns):
        for k in range(np_):
            stack.append((s0 + (s1 - s0) * i / ns, s0 + (s1 - s0) * (i + 1) / ns,
                          p0 + (p1 - p0) * k / np_, p0 + (p1 - p0) * (k + 1) / np_))
    leaves = []
    fails = []
    nev = 0
    best = None      # min certified |G_ss| lower bound over leaves
    worst_leaf = None
    while stack:
        a0, a1, b0, b1 = stack.pop()
        nev += 1
        if nev > max_boxes:
            raise CertError("box budget exceeded")
        try:
            G = gss_box(a0, a1, b0, b1, mode)
            okb = (G > THR) if sign > 0 else (G < -THR)
        except CertError:
            okb = False
            G = None
        if okb:
            bound = G.lower() if sign > 0 else -G.upper()   # arb lower bound of |G_ss| on the leaf
            bf = Fr(point(bound))
            if best is None or bf < best:
                best = bf
                worst_leaf = (str(a0), str(a1), str(b0), str(b1), interval(G))
            leaves.append((a0, a1, b0, b1))
            continue
        if (a1 - a0) <= MIN_WIDTH and (b1 - b0) <= MIN_WIDTH:
            fails.append((str(a0), str(a1), str(b0), str(b1), None if G is None else G.str(10)))
            continue
        if (a1 - a0) * ws >= (b1 - b0) * wp and (a1 - a0) > MIN_WIDTH or (b1 - b0) <= MIN_WIDTH:
            m = (a0 + a1) / 2
            stack += [(a0, m, b0, b1), (m, a1, b0, b1)]
        else:
            m = (b0 + b1) / 2
            stack += [(a0, a1, b0, m), (a0, a1, m, b1)]
        if verbose and nev % 5000 == 0:
            print("  %s: evaluated %d, leaves %d, stack %d, %.0fs" % (name, nev, len(leaves), len(stack),
                                                                   time.time() - t0), flush=True)
    # exact coverage check: leaves tile the region (areas add up, all inside, generated by bisection)
    area = sum((l[1] - l[0]) * (l[3] - l[2]) for l in leaves)
    full = (s1 - s0) * (p1 - p0)
    for l in leaves:
        if not (s0 <= l[0] < l[1] <= s1 and p0 <= l[2] < l[3] <= p1):
            raise CertError("leaf outside region")
    h = hashlib.sha256()
    for l in sorted(leaves):
        h.update(("%s,%s,%s,%s;" % l).encode())
    res = dict(region=name, s_range=[str(s0), str(s1)], p_range=[str(p0), str(p1)], mode=mode,
               claim=("d2G/ds2 > %s" % thr if sign > 0 else "d2G/ds2 < -%s" % thr), threshold=str(thr),
               boxes_evaluated=nev, leaves=len(leaves), failures=len(fails), failed_boxes=fails[:50],
               coverage_area_ok=bool(area == full and not fails),
               min_certified_abs_Gss=str(best) if best is not None else None, worst_leaf=worst_leaf,
               bound_encoding="exact rational strings",
               leaves_sha256=h.hexdigest(), runtime_s=round(time.time() - t0, 1))
    return res


def file_sha256(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["selftest", "L1", "L2", "L3"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--threshold", default=None, help="exact rational; claim sign*G_ss > threshold")
    args = ap.parse_args()
    here = os.path.abspath(__file__)
    meta = dict(script=os.path.basename(here), script_sha256=file_sha256(here),
                python_flint=flint.__version__, arb_prec_bits=PREC, python=sys.version.split()[0],
                min_box_width=str(MIN_WIDTH), started=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    try:
        st = selftest()
        if args.what == "selftest":
            out = dict(meta, selftest=st, status="PASS")
        else:
            s0, s1, p0, p1, mode, sign, thr, ns, np_, ws, wp = REGIONS[args.what]
            if args.threshold is not None:
                thr = Fr(args.threshold)
            meta["parameters"] = dict(s0=str(s0), s1=str(s1), p0=str(p0), p1=str(p1), mode=mode,
                                      sign=sign, threshold=str(thr), initial_grid=[ns, np_],
                                      split_weights=[ws, wp])
            res = certify(args.what, threshold=thr)
            status = "PASS" if (res["failures"] == 0 and res["coverage_area_ok"]) else "FAIL"
            out = dict(meta, selftest="PASS", result=res, status=status)
    except CertError as e:
        out = dict(meta, status="FAIL", error=str(e))
    js = json.dumps(out, indent=1)
    if args.out:
        with open(args.out, "x") as f:
            f.write(js + "\n")
    print(js if args.what == "selftest" or out["status"] != "PASS" else json.dumps(
        {k: v for k, v in out.get("result", {}).items() if k != "failed_boxes"}, indent=1))
    if out["status"] != "PASS":
        print("CERTIFICATION FAILED", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
