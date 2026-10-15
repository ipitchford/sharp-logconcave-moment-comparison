#!/usr/bin/env python3
"""
strip_cert.py -- rigorous Arb certificates for the STRIP lemmas (MRTT L_p-L_2 conjecture, cycle 42).

Lemma S1:  H(s,p) := d/dp [G(s,p) - G0(p)] < 0   on  s in [1/20, 1/2],  p in [2-1/20, 2+1/20].
Lemma S2:  d/dp d^2/ds^2 G(s,p) < 0              on  s in [0, 1/20],    p in [2-1/20, 2+1/20].
(--s-lo/--s-hi/--p-lo/--p-hi certify the same inequalities on other rectangles.)

Notation (proofs in strip.md):
  Z = Z_s = s(E-1) - (1-s)(E'-1),  d = 1-2s,  c = 1-s,  r0 = d/c,  q = p-2,  L = log|Z|.
  With probability s:  Z = d + sR;  with probability c:  Z = d - cR;   R ~ Exp(1).
  Phi_A(s,q) = E|Z|^{2+q} = M(s,p),
  Phi_B(s,q) = E[|Z|^q Z (R-1)]         (M_s  = p Phi_B),
  Phi_C(s,q) = E[|Z|^q ((R-1)^2 + 1)]   (M_ss = p(p-1) Phi_C).
  G    = log(Phi_A)/p - log(N)/2,      G_s = Phi_B/Phi_A - N'/(2N),
  G_ss = (p-1) Phi_C/Phi_A - p (Phi_B/Phi_A)^2 - (log N)''/2.
Every Phi_V(s, .) is a truncated power series in q (python-flint arb_series, orders 0..2):
  * the part with Z <= d ("negative side") is evaluated in closed form (Gamma functions and a
    rapidly convergent series in r0), exactly as an analytic function of q;
  * the part with Z > d ("positive side") is, on a q-box with centre qT and half-width a, the
    Taylor polynomial sum_{k<K} P_k (q-qT)^k/k! with P_k = E[V |Z|^qT L^k ; Z > d] a rigorous
    Arb integral (acb.integral + explicit tail), plus an explicit remainder for every order.
Box enclosures: S1 uses a mean-value form in s and in q; S2 uses direct interval evaluation over
the s-box with a mean-value form in q.  Undecided or non-finite comparisons raise CertError,
and the program then exits with status 2.
"""
import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from certificate_io import point, interval, positive_int, rational, runtime
import time
from fractions import Fraction as Fr

from flint import arb, acb, arb_series, ctx

PREC = 128
K_TAYLOR = 12        # Taylor order in q for the positive-side part
N_SER = 40           # truncation of the r0-series in the closed form for 0 < Z < d
L_GAMMA = 24         # length of the exact Gamma Taylor series at the centre of a q-box
CAP = 3              # series length used for all q-series (orders 0, 1, 2)
T_BOUND = 400        # integration cut-off inside the remainder moment bound
MAX_QW = Fr(1, 10)   # maximal q-box width (so |q - qT| <= 1/20)


class CertError(Exception):
    pass


def fail(msg):
    raise CertError(msg)


def qa(fr):
    """Arb ball containing the exact rational fr."""
    fr = Fr(fr)
    return arb(fr.numerator) / arb(fr.denominator)


def hull(lo, hi):
    """Arb ball containing the closed interval [lo, hi] (exact rationals)."""
    lo, hi = Fr(lo), Fr(hi)
    if lo > hi:
        fail("hull: empty interval")
    return qa(lo).union(qa(hi))


def finite(x, what):
    if not x.is_finite():
        fail(f"non-finite enclosure: {what}")
    return x


def ball0(r):
    """Ball [-r, r] for an arb r (uses an upper bound of |r|)."""
    return arb(0, r.abs_upper())


def t_cut(qT):
    """Truncation point of the positive-side r-integrals."""
    return 120 + 20 * max(0, math.ceil(qT))


# ----------------------------------------------------------------------------------------------
# Gamma(beta + q) as a series in t (q = q0 + t), uniformly for q0 in [qlo, qhi].
# Exact Taylor series at zc = beta + qc (length L_GAMMA), recentred on u = q0 - qc in [-rho, rho],
# plus a Cauchy-estimate tail: |Gamma^{(n)}(zc)/n!| <= M / r^n on the disc |w - zc| <= r with
# r = min(9/10, zc - 1/20) and M = max(Gamma(zc - r), Gamma(zc + r)) (|Gamma(w)| <= Gamma(Re w),
# and Gamma is log-convex on (0, inf), so its maximum on [zc-r, zc+r] is at an endpoint).
# ----------------------------------------------------------------------------------------------
_gamma_cache = {}


def gamma_series(beta, qlo, qhi):
    key = (beta, qlo, qhi)
    if key in _gamma_cache:
        return _gamma_cache[key]
    qc = (qlo + qhi) / 2
    rho = (qhi - qlo) / 2
    zc = beta + qc
    r = min(Fr(9, 10), zc - Fr(1, 20))
    if r < Fr(1, 2):
        fail("gamma_series: centre too close to the pole")
    theta = rho / r
    if theta > Fr(1, 10):
        fail("gamma_series: q-box too wide")
    old = ctx.cap
    ctx.cap = L_GAMMA
    g = arb_series([qa(zc), arb(1)]).gamma()
    gc = [finite(g[n], "gamma coeff") for n in range(L_GAMMA)]
    ctx.cap = old
    u = hull(-rho, rho)
    y = arb_series([u, arb(1)])
    acc = arb_series([gc[L_GAMMA - 1]])
    for n in range(L_GAMMA - 2, -1, -1):
        acc = acc * y + gc[n]
    M = qa(zc - r).gamma().union(qa(zc + r).gamma())
    M = arb(M.upper())
    coeffs = []
    for j in range(CAP):
        # tail_j <= (M / r^j) * binom(L,j) theta^(L-j) / (1 - 2 theta)
        tb = M / qa(r) ** j * arb(math.comb(L_GAMMA, j)) * qa(theta) ** (L_GAMMA - j) / (1 - 2 * qa(theta))
        coeffs.append(acc[j] + ball0(tb))
    out = arb_series(coeffs)
    _gamma_cache[key] = out
    return out


# ----------------------------------------------------------------------------------------------
# Negative side (Z <= d), closed form, exact in q.  sb is an arb ball for s.
#   0 < Z < d :  exp(-r0) d^{m+1+q} sum_n r0^n/n! sum_i w_i/(beta_i + q + n)
#   Z < 0     :  Gamma-function terms.
# For an s-ball reaching s = 1/2 (d = 0) the 0<Z<d piece is replaced by a rigorous bound.
# ----------------------------------------------------------------------------------------------

def inv_lin_series(z):
    """Series of 1/(z + t) (z an arb ball, z > 1/2)."""
    if not (z.lower() > arb('0.5')):
        fail("inv_lin_series: argument too small")
    iz = 1 / z
    return arb_series([iz, -iz * iz, iz * iz * iz][:CAP])


def pow_series(base, beta, xq):
    """base^(beta + q) as a series (base > 0)."""
    if not (base.lower() > 0):
        fail("pow_series: base not positive")
    return ((xq + beta) * base.log()).exp()


def apiece_bound(Dhi, beta, psimax):
    """Coefficient bounds for the 0<Z<d piece when d in [0, Dhi], Dhi <= 1/10:
    |d^j/dq^j piece| <= psimax * int_0^Dhi z^beta |log z|^j dz  (beta = m + qlo > -1);
    coefficient j = derivative / j!."""
    if Dhi > Fr(1, 10):
        fail("apiece_bound: D too large")
    D = qa(Dhi)
    beta1 = beta + 1
    if not (beta1.lower() > 0):
        fail("apiece_bound: exponent too small")
    lD = -D.log()                     # -log D >= 0
    coeffs = []
    for j in range(CAP):
        tot = arb(0)
        for i in range(j + 1):
            tot += arb(math.factorial(j) // math.factorial(j - i)) * lD ** (j - i) / beta1 ** (i + 1)
        b = arb(psimax) * (beta1 * D.log()).exp() * tot / math.factorial(j)
        coeffs.append(ball0(b))
    return arb_series(coeffs)


def neg_series(kind, sb, qlo, qhi, d_touches_zero=False, Dhi=None, d_zero_exact=False):
    """Series (in t, q = q0 + t, q0 in [qlo,qhi]) of the Z <= d contribution to Phi_kind."""
    q0 = hull(qlo, qhi)
    xq = arb_series([q0, arb(1)])
    d = 1 - 2 * sb
    c = 1 - sb
    r0 = d / c
    er0 = (-r0).exp()
    if kind == 'A':
        G3 = gamma_series(3, qlo, qhi)
        bpart = er0 * pow_series(c, 3, xq) * G3
        m, betas, ws, psimax = 2, [3], [arb(1)], 1
    elif kind == 'B':
        G2 = gamma_series(2, qlo, qhi)
        G3 = gamma_series(3, qlo, qhi)
        bpart = -er0 * (pow_series(c, 2, xq) * G3 - sb * pow_series(c, 1, xq) * G2)
        m, betas, ws, psimax = 1, [2, 3], [r0 - 1, -r0], 1
    elif kind == 'C':
        G1 = gamma_series(1, qlo, qhi)
        G2 = gamma_series(2, qlo, qhi)
        G3 = gamma_series(3, qlo, qhi)
        bpart = er0 / (c * c) * (pow_series(c, 3, xq) * G3 - 2 * sb * pow_series(c, 2, xq) * G2
                                 + (sb * sb + c * c) * pow_series(c, 1, xq) * G1)
        al = r0 - 1
        m, betas, ws, psimax = 0, [1, 2, 3], [al * al + 1, -2 * al * r0, r0 * r0], 2
    else:
        fail("unknown kind")
    if d_zero_exact:
        return bpart
    if d_touches_zero:
        if m == 0:
            fail("neg_series: kind C not supported for d touching 0")
        return bpart + apiece_bound(Dhi, m + qa(qlo), psimax)
    # 0 < Z < d piece
    tot = arb_series([arb(0)])
    fact = arb(1)
    r0n = arb(1)
    for n in range(N_SER + 1):
        term = arb_series([arb(0)])
        for b, w in zip(betas, ws):
            term = term + w * inv_lin_series(q0 + (b + n))
        tot = tot + (r0n / fact) * term
        r0n = r0n * r0
        fact = fact * (n + 1)
    # tail: sum_{n>N} r0^n/n! * sum|w_i| * |coeff_j of 1/(beta+q+n+t)| <= sum|w| * 2/(N+1)!
    # (r0 <= 1 and beta + q + n >= 1 for n >= 1)
    wsum = sum((abs(w) for w in ws), arb(0))
    tail = wsum * 2 / arb(math.factorial(N_SER + 1))
    tot = tot + arb_series([ball0(tail)] * CAP)
    apart = er0 * pow_series(d, m + 1, xq) * tot
    return bpart + apart


# ----------------------------------------------------------------------------------------------
# Positive side (Z > d):  P_{V,k}(s; qT) = s * int_0^inf v(r,x) x^qT log(x)^k e^{-r} dr,
#   x = 1 + s(r-2);  v = x^2 (A), x (r-1) (B), (r-1)^2 + 1 (C).
# Integral on [r_lo, T] by acb.integral; [T, inf) bounded by s*Gamma(k+3+max(qT,0), T); [0, r_lo]
# bounded explicitly when the s-ball reaches 1/2 (d = 0).
# ----------------------------------------------------------------------------------------------

def _integrand(kind, sacb, k, qT):
    qTa = acb(qa(qT))
    plain = (k == 0 and qT == 0)

    def f(r, analytic):
        x = 1 + sacb * (r - 2)
        if kind == 'A':
            v = x * x
        elif kind == 'B':
            v = x * (r - 1)
        else:
            v = (r - 1) * (r - 1) + 1
        w = v * (-r).exp()
        if plain:
            return w
        L = x.log(analytic=analytic)
        if qT != 0:
            w = w * (qTa * L).exp()
        if k > 0:
            w = w * L ** k
        return w
    return f


def g_sup(mu, k, X):
    """Upper bound for sup_{0 < x <= X} x^mu |log x|^k   (mu > 0 arb, X > 0 arb)."""
    if not (mu.lower() > 0):
        fail("g_sup: exponent must be positive")
    Xu = arb(X.upper())
    if k == 0:
        return arb((Xu.log() * mu).exp().upper()) if Xu.upper() > 0 else arb(0)
    thr = (arb(-k) / mu).exp()
    cand = []
    if Xu.upper() < thr.lower():
        cand.append((mu * Xu.log()).exp() * (-Xu.log()) ** k)
    else:
        cand.append((arb(k) / mu) ** k * arb(-k).exp())
        if Xu.upper() > 1:
            cand.append((mu * Xu.log()).exp() * Xu.log() ** k)
    best = cand[0]
    for cnd in cand[1:]:
        best = best.union(cnd)
    return arb(best.upper())


_pos_cache = {}
STATS = dict(integrals=0)


def pos_moments(kind, slo, shi, qT, point=False):
    """[P_0, ..., P_{K-1}] enclosing the positive-side moments (Taylor centre qT) for s in [slo,shi]."""
    key = (kind, slo, shi, qT, point)
    if key in _pos_cache:
        return _pos_cache[key]
    if point and slo != shi:
        fail("pos_moments: point requires slo == shi")
    sb = qa(slo) if point else hull(slo, shi)
    if shi == 0:
        out = [arb(0)] * K_TAYLOR
        _pos_cache[key] = out
        return out
    special = (1 - 2 * shi) <= 0
    if special:
        r_lo = max(4 * (shi - slo), Fr(1, 10 ** 6))
        if kind == 'C':
            fail("pos_moments: kind C with d touching 0")
    else:
        r_lo = Fr(0)
    sacb = acb(sb)
    if point:
        opts = dict(rel_tol=arb(2) ** -100, abs_tol=arb(2) ** -110, eval_limit=400000)
    else:
        opts = dict(rel_tol=arb(2) ** -40, abs_tol=arb(2) ** -60, eval_limit=400000)
    T = t_cut(qT)
    Ta = arb(T)
    shi_a = qa(shi)
    out = []
    for k in range(K_TAYLOR):
        I = acb.integral(_integrand(kind, sacb, k, qT), qa(r_lo), Ta, **opts)
        STATS['integrals'] += 1
        Ir = I.real
        if not Ir.is_finite() or not I.imag.is_finite():
            fail(f"positive-side integral not finite: kind={kind} k={k} s=[{slo},{shi}] qT={qT}")
        if not I.imag.contains(0):
            fail("positive-side integral: imaginary part excludes 0")
        tailb = shi_a * Ta.gamma_upper(arb(k + 3) + qa(max(qT, 0)))   # Gamma(k+3+qT+, T)
        val = sb * Ir + ball0(tailb)
        if special:
            # |s int_0^{r_lo} v x^qT log^k x e^{-r} dr| <= s_hi r_lo psimax sup x^{m+qT} |log x|^k,
            # x in (0, X], X = d_hi + s_hi r_lo, d_hi = 1 - 2 slo; psimax = 1 for A and B.
            m = 2 if kind == 'A' else 1
            X = (1 - 2 * qa(slo)) + shi_a * qa(r_lo)
            b = shi_a * qa(r_lo) * g_sup(arb(m) + qa(qT), k, X)
            val = val + ball0(b)
        out.append(finite(val, f"P_{kind},{k}"))
    _pos_cache[key] = out
    return out


# ----------------------------------------------------------------------------------------------
# Remainder bound for the positive-side Taylor polynomial (uniform in s in [0,1/2], |q-qT| <= a):
#   E_pos[|V| |Z|^qT |L|^K e^{a|L|}]
#     <= 2 C_psi e^{1/2} sum_{sigma=+-1} int_0^inf x^{m+qT+sigma a}|log x|^K e^{-x/2} dx.
# ----------------------------------------------------------------------------------------------
_B_cache = {}


def moment_bound(kind, qT, a):
    key = (kind, qT, a)
    if key in _B_cache:
        return _B_cache[key]
    K = K_TAYLOR
    m, Cpsi = {'A': (2, 1), 'B': (1, 1), 'C': (0, 2)}[kind]
    tot = arb(0)
    for sg in (1, -1):
        beta = arb(m) + qa(qT) + sg * qa(a)
        if not (beta.lower() > -1):
            fail("moment_bound: exponent <= -1")
        I0 = arb(math.factorial(K)) / (beta + 1) ** (K + 1)    # int_0^1 x^beta |log x|^K dx

        def f(x, analytic, beta=beta):
            return x.pow(acb(beta), analytic=analytic) * x.log(analytic=analytic) ** K * (-x / 2).exp()
        Tb = arb(T_BOUND)
        J = acb.integral(f, 1, Tb, rel_tol=arb(2) ** -40, abs_tol=arb(2) ** -60, eval_limit=400000)
        STATS['integrals'] += 1
        if not J.real.is_finite():
            fail("moment_bound integral not finite")
        # tail (log x <= x): int_Tb^inf x^{beta+K} e^{-x/2} dx = 2^{beta+K+1} Gamma(beta+K+1, Tb/2)
        tail = arb(2) ** (beta + K + 1) * (Tb / 2).gamma_upper(beta + K + 1)
        tot += I0 + J.real + tail
    B = arb((2 * Cpsi * arb('0.5').exp() * tot).upper())
    _B_cache[key] = B
    return B


def pos_series(kind, P, qlo, qhi, qT):
    """Series of the positive-side part of Phi_kind at q0 in [qlo,qhi] (Taylor centre qT)."""
    u0 = hull(qlo - qT, qhi - qT)
    xu = arb_series([u0, arb(1)])
    K = K_TAYLOR
    acc = arb_series([P[K - 1] / math.factorial(K - 1)])
    for k in range(K - 2, -1, -1):
        acc = acc * xu + P[k] / math.factorial(k)
    umax = max(abs(qlo - qT), abs(qhi - qT))
    if umax > MAX_QW / 2:
        fail("pos_series: |q - qT| too large")
    if umax == 0:
        return acc
    B = moment_bound(kind, qT, umax)
    rem = []
    for j in range(CAP):
        rb = qa(umax) ** (K - j) / math.factorial(K - j) * B / math.factorial(j)
        rem.append(ball0(rb))
    return acc + arb_series(rem)


def phi_series(kind, slo, shi, qlo, qhi, point_s=False):
    """Series in t of Phi_kind(s, q0 + t), for all s in [slo,shi] (or s = slo if point_s) and
    q0 in [qlo,qhi]; the positive side uses the Taylor centre qT = (qlo+qhi)/2."""
    qT = (qlo + qhi) / 2
    if point_s:
        sb = qa(slo)
        special = False
        Dhi = None
    else:
        sb = hull(slo, shi)
        special = (1 - 2 * shi) <= 0
        Dhi = 1 - 2 * slo
    neg = neg_series(kind, sb, qlo, qhi, d_touches_zero=special, Dhi=Dhi,
                     d_zero_exact=(point_s and slo == Fr(1, 2)))
    if slo == 0 and shi == 0:
        return neg      # s = 0: the positive side has probability 0
    P = pos_moments(kind, slo, shi, qT, point=point_s)
    return neg + pos_series(kind, P, qlo, qhi, qT)


# ----------------------------------------------------------------------------------------------
# Derived quantities
# ----------------------------------------------------------------------------------------------

def gp_series(PhiA, qlo, qhi):
    """Series of G(s, 2+q) + log(N)/2 = log(Phi_A)/(2+q)."""
    q0 = hull(qlo, qhi)
    xq = arb_series([q0, arb(1)])
    if not (PhiA[0].lower() > 0):
        fail("Phi_A not certified positive")
    return PhiA.log() / (xq + 2)


def H_center(sc, qc):
    """H(sc, qc) = G_p(sc, 2+qc) - G0'(2+qc) at an exact point."""
    gA = gp_series(phi_series('A', sc, sc, qc, qc, point_s=True), qc, qc)
    g0 = gp_series(phi_series('A', Fr(0), Fr(0), qc, qc, point_s=True), qc, qc)
    return gA[1] - g0[1]


def H_q(sc, qlo, qhi):
    """Enclosure of dH/dq (s = sc fixed) over q in [qlo,qhi]."""
    gA = gp_series(phi_series('A', sc, sc, qlo, qhi, point_s=True), qlo, qhi)
    g0 = gp_series(phi_series('A', Fr(0), Fr(0), qlo, qhi, point_s=True), qlo, qhi)
    return 2 * (gA[2] - g0[2])


def H_s(slo, shi, qlo, qhi):
    """Enclosure of dH/ds = d/dq (Phi_B/Phi_A) over the box."""
    A = phi_series('A', slo, shi, qlo, qhi)
    B = phi_series('B', slo, shi, qlo, qhi)
    if not (A[0].lower() > 0):
        fail("Phi_A not positive on box")
    return (B / A)[1]


def gss_series(slo, shi, qlo, qhi, point_s=False):
    """Series of G_ss(s, 2+q) + (log N)''/2 = (1+q) Phi_C/Phi_A - (2+q) (Phi_B/Phi_A)^2."""
    A = phi_series('A', slo, shi, qlo, qhi, point_s)
    B = phi_series('B', slo, shi, qlo, qhi, point_s)
    C = phi_series('C', slo, shi, qlo, qhi, point_s)
    q0 = hull(qlo, qhi)
    xq = arb_series([q0, arb(1)])
    if not (A[0].lower() > 0):
        fail("Phi_A not positive on box")
    rB = B / A
    return (xq + 1) * C / A - (xq + 2) * rB * rB


# ----------------------------------------------------------------------------------------------
# Box checks
# ----------------------------------------------------------------------------------------------

def check_S1_box(slo, shi, qlo, qhi):
    """Enclosure of H over the box: H(sc,qc) + H_q(sc,Q)(Q-qc) + H_s(S,Q)(S-sc)."""
    sc = (slo + shi) / 2
    qc = (qlo + qhi) / 2
    h0 = H_center(sc, qc)
    hq = H_q(sc, qlo, qhi)
    hs = H_s(slo, shi, qlo, qhi)
    tq = hq * hull(qlo - qc, qhi - qc)
    ts = hs * hull(slo - sc, shi - sc)
    enc = h0 + tq + ts
    finite(enc, "S1 enclosure")
    return enc, dict(c=h0, tq=tq, ts=ts)


def check_S2_box(slo, shi, qlo, qhi):
    """Enclosure of G_pss over the box: G_pss(S,qc) + G_ppss(S,Q)(Q-qc)."""
    qc = (qlo + qhi) / 2
    g_c = gss_series(slo, shi, qc, qc)          # G_pss at q = qc over the s-box
    g_Q = gss_series(slo, shi, qlo, qhi)        # G_ppss over the box
    tq = 2 * g_Q[2] * hull(qlo - qc, qhi - qc)
    enc = g_c[1] + tq
    finite(enc, "S2 enclosure")
    return enc, dict(c=g_c[1], tq=tq)


def adaptive(name, check, s_range, q_range, s_init, q_init, max_boxes=2000000, verbose=False):
    """Bisect until every box has a certified negative upper bound."""
    s0, s1 = s_range
    q0, q1 = q_range
    for value, label in [(s_init, 's_init'), (q_init, 'q_init'), (max_boxes, 'max_boxes')]:
        positive_int(value, label)
    if not (0 <= s0 < s1 <= Fr(1, 2) and Fr(-1, 5) <= q0 < q1 <= 10):
        fail("invalid adaptive rectangle")
    if name not in ('S1', 'S2') or (name == 'S1' and s0 <= 0) or (name == 'S2' and s1 > Fr(1, 10)):
        fail("invalid lemma/domain")
    if (q1 - q0) / q_init > MAX_QW:
        fail("initial q-boxes too wide")
    work = []
    for i in range(s_init):
        for j in range(q_init):
            work.append((s0 + (s1 - s0) * i / s_init, s0 + (s1 - s0) * (i + 1) / s_init,
                         q0 + (q1 - q0) * j / q_init, q0 + (q1 - q0) * (j + 1) / q_init))
    work.reverse()
    leaves = []
    nchecks = 0
    while work:
        slo, shi, qlo, qhi = work.pop()
        nchecks += 1
        if nchecks > max_boxes:
            fail(f"{name}: box budget exhausted")
        enc, parts = check(slo, shi, qlo, qhi)
        if not enc.is_finite():
            fail("nonfinite leaf enclosure")
        if enc.upper() < 0:
            lo, hi = interval(enc)
            leaves.append(dict(s=[str(slo), str(shi)], q=[str(qlo), str(qhi)],
                               lower=lo, upper=hi))
            continue
        if enc.lower() >= 0:
            fail(f"{name}: enclosure >= 0 on box s=[{slo},{shi}] p=[{2 + qlo},{2 + qhi}] "
                 f"(the inequality is FALSE there): {enc}")
        if (shi - slo) < Fr(1, 10 ** 7) or (qhi - qlo) < Fr(1, 10 ** 7):
            fail(f"{name}: cannot certify box s=[{slo},{shi}] q=[{qlo},{qhi}] enc={enc}")
        ws = parts['ts'].rad() if 'ts' in parts else parts['c'].rad()
        wq = parts['tq'].rad()
        split_s = (ws > wq)
        if split_s:
            sm = (slo + shi) / 2
            work.append((sm, shi, qlo, qhi))
            work.append((slo, sm, qlo, qhi))
        else:
            qm = (qlo + qhi) / 2
            work.append((slo, shi, qm, qhi))
            work.append((slo, shi, qlo, qm))
        if verbose:
            print(f"  split {'s' if split_s else 'q'}: s=[{float(slo):.6f},{float(shi):.6f}] "
                  f"p=[{2 + float(qlo):.5f},{2 + float(qhi):.5f}] enc={enc}", flush=True)
    return leaves, nchecks


def verify_cover(leaves, s_range, q_range):
    """Exact-rational check that the leaves tile the rectangle: containment, pairwise disjoint
    interiors, and total area equal to the area of the rectangle."""
    if not leaves or not (s_range[0] < s_range[1] and q_range[0] < q_range[1]):
        fail("empty/reversed coverage domain")
    boxes = []
    area = Fr(0)
    for b in leaves:
        slo, shi = Fr(b['s'][0]), Fr(b['s'][1])
        qlo, qhi = Fr(b['q'][0]), Fr(b['q'][1])
        if not (s_range[0] <= slo < shi <= s_range[1] and q_range[0] <= qlo < qhi <= q_range[1]):
            fail("leaf outside the region")
        boxes.append((slo, shi, qlo, qhi))
        area += (shi - slo) * (qhi - qlo)
    if area != (s_range[1] - s_range[0]) * (q_range[1] - q_range[0]):
        fail("leaves do not tile the region (area mismatch)")
    boxes.sort()
    # pairwise interior-disjointness (sweep over s)
    for i in range(len(boxes)):
        a = boxes[i]
        for j in range(i + 1, len(boxes)):
            b = boxes[j]
            if b[0] >= a[1]:
                break
            if b[2] < a[3] and a[2] < b[3]:
                fail("overlapping leaves")
    return True


# ----------------------------------------------------------------------------------------------
# Validation of the evaluator
# ----------------------------------------------------------------------------------------------

def G_direct(s, p):
    """Non-Taylor evaluation of G(s,p) at a point (validates closed forms and integrands)."""
    s = Fr(s)
    p = Fr(p)
    q = p - 2
    sb = qa(s)
    d = 1 - 2 * sb
    c = 1 - sb
    r0 = d / c
    er0 = (-r0).exp()
    qb = qa(q)
    bpart = er0 * c ** (3 + qb) * (3 + qb).gamma()
    tot = arb(0)
    fact = arb(1)
    r0n = arb(1)
    for n in range(80):
        tot += r0n / fact / (3 + qb + n)
        r0n *= r0
        fact *= (n + 1)
    apart = er0 * d ** (3 + qb) * tot if s < Fr(1, 2) else arb(0)
    pos = arb(0)
    if s > 0:
        sacb = acb(sb)
        pb = acb(qa(p))

        def f(r, analytic):
            x = 1 + sacb * (r - 2)
            return x.pow(pb, analytic=analytic) * (-r).exp()
        lo = 0 if s < Fr(1, 2) else Fr(1, 10 ** 30)
        I = acb.integral(f, qa(lo), 400, rel_tol=arb(2) ** -90, abs_tol=arb(2) ** -100)
        pos = sb * I.real
    M = bpart + apart + pos
    Nn = sb * sb + c * c
    return M.log() / qa(p) - Nn.log() / 2


def validate(report):
    t0 = time.time()
    # 1. SPEC reference values of G (mpmath, 40 digits)
    refs = [(Fr(1, 10), Fr(13, 10), '-0.213049905296183'), (Fr(3, 10), Fr(19, 10), '-0.0294149150315886'),
            (Fr(1, 4), Fr(5, 2), '0.143192764694943'), (Fr(2, 5), Fr(5), '0.641602497665675'),
            (Fr(1, 20), Fr(12), '1.58099627165213'), (Fr(1, 2), Fr(6, 5), '-0.265784034621107'),
            (Fr(0), Fr(33, 10), '0.375483768956575')]
    out = []
    for s, p, v in refs:
        g = G_direct(s, p)
        ok = abs(g - arb(v)).upper() < arb('1e-13')
        out.append(dict(s=str(s), p=str(p), G=g.str(16, radius=False), ref=v, ok=bool(ok)))
        if not ok:
            fail(f"reference value mismatch at s={s}, p={p}: {g} vs {v}")
    # 2. independent mpmath values (quadrature of the derivative formulas, checked against
    #    finite differences) of G_p and G_pss; the s = 1/2 values of G_pss are 4(2q/p - q^2/p^2).
    mp_refs = [(Fr(0), Fr(49, 25), '0.303066746695', '-0.285602440745'),
               (Fr(0), Fr(2), '0.302825116765', '-0.302825116765'),
               (Fr(1, 50), Fr(51, 25), '0.30248041025', '-0.32905848903'),
               (Fr(1, 5), Fr(2), '0.29703399374', '-0.102302989802'),
               (Fr(1, 5), Fr(49, 25), '0.297880638504', '-0.0520501911232'),
               (Fr(1, 2), Fr(49, 25), '0.291778657056', '-0.164931278634'),
               (Fr(1, 2), Fr(51, 25), '0.284526345309', '0.155324875048'),
               (Fr(9, 20), Fr(19, 10), '0.297035534249786', '-0.114298244072489'),
               (Fr(0), Fr(3), '0.278693363853009', '-0.377068248826325'),
               (Fr(3, 100), Fr(6), '0.156962562092662', '-0.0405492504447171'),
               (Fr(1, 10), Fr(9, 2), '0.207614739600153', '-0.264936077108685'),
               (Fr(3, 10), Fr(12), '0.0781360751243442', '-0.0420693239401909'),
               (Fr(1, 20), Fr(12), '0.0786346344575135', '-0.00850472825590062'),
               (Fr(1, 2), Fr(7), '0.113965897645385', '3.6734693877551')]
    out2 = []
    for s, p, gp_ref, gpss_ref in mp_refs:
        q = p - 2
        gp = gp_series(phi_series('A', s, s, q, q, point_s=True), q, q)[1]
        gpss = gss_series(s, s, q, q, point_s=True)[1] if s < Fr(1, 2) else None
        tol = arb('1e-10')
        ok1 = abs(gp - arb(gp_ref)).upper() < tol
        ok2 = True if gpss is None else bool(abs(gpss - arb(gpss_ref)).upper() < tol)
        out2.append(dict(s=str(s), p=str(p), Gp=gp.str(15, radius=False), Gp_ref=gp_ref,
                         Gpss=(gpss.str(15, radius=False) if gpss is not None else None), Gpss_ref=gpss_ref,
                         ok=bool(ok1 and ok2)))
        if not (ok1 and ok2):
            fail(f"mpmath cross-check failed at s={s} p={p}: Gp={gp} Gpss={gpss}")
    # 3. structural identities: G_ss(s,2) = 0, G_s(0,p) = 0, and Taylor-centre consistency
    for s in [Fr(0), Fr(1, 40), Fr(3, 10), Fr(9, 20)]:
        gss = gss_series(s, s, Fr(0), Fr(0), point_s=True)[0]
        Nn = qa(s) ** 2 + (1 - qa(s)) ** 2
        N1 = 4 * qa(s) - 2
        val = gss - (4 / Nn - (N1 / Nn) ** 2) / 2
        if not (abs(val).upper() < arb('1e-25')):
            fail(f"G_ss(s,2) != 0 at s={s}: {val}")
    for q in [Fr(-1, 20), Fr(1, 37), Fr(7)]:
        A = phi_series('A', Fr(0), Fr(0), q, q, point_s=True)
        B = phi_series('B', Fr(0), Fr(0), q, q, point_s=True)
        gs0 = (B / A)[0] + 1   # G_s(0,p) = Phi_B/Phi_A - N'(0)/(2N(0)) = Phi_B/Phi_A + 1
        if not (abs(gs0).upper() < arb('1e-25')):
            fail(f"G_s(0,p) != 0 at q={q}: {gs0}")
    # Taylor-centre consistency: Phi evaluated on a q-box (centre qT) must contain the point values
    for kind in 'ABC':
        s, qlo, qhi = Fr(7, 100), Fr(-1, 20), Fr(1, 20)
        box = phi_series(kind, s, s, qlo, qhi, point_s=True)
        for q in [qlo, Fr(-1, 50), qhi]:
            pt = phi_series(kind, s, s, q, q, point_s=True)
            for j in range(CAP):
                if not box[j].contains(pt[j]):
                    fail(f"Taylor-centre consistency failed: kind={kind} q={q} j={j}")
    report['validation'] = dict(reference_G=out, mpmath_Gp_Gpss=out2,
                                identities="G_ss(s,2)=0 (4 points), G_s(0,p)=0 (3 points), "
                                           "Taylor-centre consistency (3 kinds x 3 points): OK",
                                seconds=round(time.time() - t0, 2))


# ----------------------------------------------------------------------------------------------

def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        h.update(fh.read())
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lemma', choices=['S1', 'S2', 'validate'], required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--s-lo', default=None)
    ap.add_argument('--s-hi', default=None)
    ap.add_argument('--p-lo', default='39/20')
    ap.add_argument('--p-hi', default='41/20')
    ap.add_argument('--s-init', type=int, default=None)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('--prec', type=int, default=None, help='override PREC (bits)')
    ap.add_argument('--k-taylor', type=int, default=None, help='override K_TAYLOR')
    ap.add_argument('--n-ser', type=int, default=None, help='override N_SER')
    ap.add_argument('--l-gamma', type=int, default=None, help='override L_GAMMA')
    args = ap.parse_args()
    global PREC, K_TAYLOR, N_SER, L_GAMMA
    runtime()
    PREC = PREC if args.prec is None else args.prec
    K_TAYLOR = K_TAYLOR if args.k_taylor is None else args.k_taylor
    N_SER = N_SER if args.n_ser is None else args.n_ser
    L_GAMMA = L_GAMMA if args.l_gamma is None else args.l_gamma
    if args.s_init is not None:
        positive_int(args.s_init, 's_init')
    if K_TAYLOR < 4 or N_SER < 20 or L_GAMMA < 12 or PREC < 64:
        fail("parameter override out of the supported range")
    ctx.prec = PREC
    ctx.cap = CAP
    if os.path.exists(args.out):
        fail(f"refusing to overwrite existing {args.out}")
    t0 = time.time()
    p_lo, p_hi = Fr(args.p_lo), Fr(args.p_hi)
    q_range = (p_lo - 2, p_hi - 2)
    if not (Fr(-1, 5) <= q_range[0] < q_range[1] <= 10):
        fail("p-range outside the supported window [9/5, 12]")
    q_init = math.ceil((q_range[1] - q_range[0]) / MAX_QW)
    report = dict(lemma=args.lemma, script=os.path.basename(__file__), script_sha256=sha256_of(__file__),
                  flint_version=__import__('flint').__version__, python=sys.version.split()[0],
                  params=dict(PREC=PREC, K_TAYLOR=K_TAYLOR, N_SER=N_SER, L_GAMMA=L_GAMMA, CAP=CAP,
                              T_BOUND=T_BOUND, MAX_QW=str(MAX_QW), T_CUT="120 + 20*ceil(max(qT,0))"))
    validate(report)
    if args.lemma == 'S1':
        s_range = (Fr(args.s_lo or '1/20'), Fr(args.s_hi or '1/2'))
        if s_range[1] > Fr(1, 2) or s_range[0] <= 0:
            fail("S1 s-range must lie in (0, 1/2]")
        s_init = 9 if args.s_init is None else args.s_init
        leaves, n = adaptive('S1', check_S1_box, s_range, q_range, s_init, q_init, verbose=args.verbose)
        statement = (f"H(s,p) = d/dp[G(s,p)-G0(p)] < 0 for s in [{s_range[0]},{s_range[1]}], "
                     f"p in [{p_lo},{p_hi}]")
    elif args.lemma == 'S2':
        s_range = (Fr(args.s_lo or '0'), Fr(args.s_hi or '1/20'))
        if s_range[0] < 0 or s_range[1] > Fr(1, 10):
            fail("S2 s-range must lie in [0, 1/10]")
        s_init = 1 if args.s_init is None else args.s_init
        leaves, n = adaptive('S2', check_S2_box, s_range, q_range, s_init, q_init, verbose=args.verbose)
        statement = (f"d/dp d^2G/ds^2 (s,p) < 0 for s in [{s_range[0]},{s_range[1]}], "
                     f"p in [{p_lo},{p_hi}]")
    else:
        leaves, n, statement = [], 0, "validation only"
    if args.lemma in ('S1', 'S2'):
        verify_cover(leaves, s_range, q_range)
        worst = max(leaves, key=lambda b: Fr(b['upper']))
        report.update(statement=statement, s_range=[str(s_range[0]), str(s_range[1])],
                      p_range=[str(p_lo), str(p_hi)], n_leaves=len(leaves),
                      n_box_evaluations=n, n_integrals=STATS['integrals'],
                      max_upper_bound=worst['upper'], worst_box=worst,
                      status="CERTIFIED", leaves=leaves, bound_encoding="exact rational strings")
    report['seconds'] = round(time.time() - t0, 2)
    with open(args.out, 'x') as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps({k: v for k, v in report.items() if k not in ('leaves', 'validation')}, indent=1))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except CertError as e:
        print(f"CERTIFICATION FAILURE: {e}", file=sys.stderr)
        sys.exit(2)
