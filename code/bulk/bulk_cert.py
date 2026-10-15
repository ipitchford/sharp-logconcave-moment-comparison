#!/usr/bin/env python3
"""
bulk_cert.py -- rigorous Arb certificate for the bulk regions B1..B4 of the
MRTT L_p-L_2 programme (cycle 42).  See bulk.md for the mathematics.

Notation.  E, E', R i.i.d. Exp(1);  Z_s = s(E-1) - (1-s)(E'-1),  s in [0,1/2];
d = 1-2s, c = 1-s, x = d/c, z = d/s;  N(s) = s^2 + (1-s)^2;
M(s,p) = E|Z_s|^p;  G(s,p) = log(M)/p - log(N)/2;
G0(p) = G(0,p) = log F(p)/p, F(p) = E|E-1|^p;  Gh(p) = G(1/2,p) = log Gamma(p+1)/p - log(2)/2.

Closed forms (bulk.md, Lemmas B.1-B.3):
  M      = s*J_p + c^{p+1} phi(x,p),
  J_q    = E(d+sR)^q = s^q e^{z} Gamma(q+1,z)          (s>0),  J_q = d^q (s=0)
  phi    = E|x-R|^p = e^{-x} (Gamma(p+1) + I(x,p)),  I(x,p) = int_0^x u^p e^u du
         = x^{p+1}/(p+1) 1F1(p+1;p+2;x)
  M_s    = (p+1) J_p - p J_{p-1} - c^{p-1} phi(x,p) ((p+1)c-1) - d^p/c
  M_p    = log(s) s J_p + s^{p+1} e^z dGamma(a,z)/da|_{a=p+1}
           + log(c) c^{p+1} phi + c^{p+1} e^{-x} (Gamma(p+1) psi(p+1) + I_p(x,p))
  I_p    = log(x) I(x,p) - x^{p+1}/(p+1)^2 2F2(p+1,p+1;p+2,p+2;x)
  dGamma(a,z)/da = Gamma(a)psi(a) - [log(z) gamma(a,z) - z^a/a^2 2F2(a,a;a+1,a+1;-z)]

Box enclosure (bulk.md, Proposition B.6): for a box [s0,s1]x[p0,p1] with centre (sc,pc),
half-widths (hs,hp) and reference R in {G0, Gh}:
  G(s,p) - R(p)  in  D(sc,pc) + [D_p](sc,[p0,p1]) [-hp,hp] + [G_s](box) [-hs,hs].
p-intervals are never passed to Gamma(.,z) as wide balls on their own: every quantity is
enclosed over [p0,p1] from exact endpoint evaluations using the monotonicity facts of
Lemmas B.4-B.5 (M_p, F', Gamma(p+1), psi, log Gamma nondecreasing in p; I(x,p) nonincreasing
in p for x <= 1; J_q convex in q and J_q/d^q nondecreasing in q).  For sc < 1/64 the
p-direction uses the interval enclosure of D(sc,[p0,p1]) instead of the mean-value form.

Usage:  bulk_cert.py REGION [--out FILE]
REGION in {B1,B2,B3,B4} (B4ext = optional p in [7,12]; CTRL_* = negative controls that must fail).
Exit status: 0 = certified, 1 = not certified, 2 = error (self-test/precision/input).
"""
import sys, time, json, hashlib, os, platform, argparse

import flint
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from certificate_io import point, interval, positive_int, rational, runtime
from flint import arb, ctx, fmpq

FLINT_REQUIRED = "0.9.0"
BASE_PREC = 128

# ----------------------------------------------------------------------------
# exact parameters (SPEC.md, cycle 42) -- all rationals
# ----------------------------------------------------------------------------
ETA = fmpq(1, 20)
SIGMA = fmpq(1, 20)       # near s = 1/2
SIGMA0 = fmpq(1, 20)      # near s = 0
DELTA = fmpq(1, 20)       # around p = 2
PMAX_B4 = fmpq(7)         # coordinator update: large-p tail proved for p >= 7
PSTAR_LO = fmpq(1682457208, 10**9)   # certified by crossing_cert.py: p* in (PSTAR_LO, PSTAR_HI)
PSTAR_HI = fmpq(1682457209, 10**9)
HALF = fmpq(1, 2)

# region: (s0, s1, p0, p1, refs, sense)
#   sense +1: certify  G - R > 0 for at least one R in refs on every box
#   sense -1: certify  R - G > 0 (refs has one element)
REGIONS = {
    # p in [1, p*-eta] is covered since p* < PSTAR_HI
    "B1": (fmpq(0), HALF - SIGMA, fmpq(1), PSTAR_HI - ETA, ("h",), +1),
    # [p*-eta, p*+eta] subset of [PSTAR_LO-eta, PSTAR_HI+eta]; G > min(G0,Gh) iff G>G0 or G>Gh
    "B2": (SIGMA0, HALF - SIGMA, PSTAR_LO - ETA, PSTAR_HI + ETA, ("0", "h"), +1),
    "B3": (SIGMA0, HALF, PSTAR_LO + ETA, 2 - DELTA, ("0",), +1),
    "B4": (SIGMA0, HALF, 2 + DELTA, PMAX_B4, ("0",), -1),
    # optional, redundant with the analytic large-p lemma (p >= 7): B4 continued to p = 12
    "B4ext": (SIGMA0, HALF, PMAX_B4, fmpq(12), ("0",), -1),
}
# negative controls: false statements, the certifier MUST report failure on these
CONTROLS = {
    "CTRL_B3_wrong_side_of_2": (SIGMA0, fmpq(1, 10), 2 + DELTA, fmpq(21, 10), ("0",), +1),
    "CTRL_B4_wrong_side_of_2": (SIGMA0, fmpq(1, 10), fmpq(19, 10), 2 - DELTA, ("0",), -1),
    "CTRL_B1_past_pstar": (fmpq(0), fmpq(1, 10), fmpq(17, 10), fmpq(18, 10), ("h",), +1),
    "CTRL_B3_overshoot": (SIGMA0, HALF, fmpq(19, 10), fmpq(201, 100), ("0",), +1),
}
REGIONS.update(CONTROLS)

MIN_WIDTH = fmpq(1, 2**34)     # boxes narrower than this in both directions -> failure
MAX_BOXES = 5 * 10**6
S_MV_P = fmpq(1, 64)          # use the mean-value form in p only when sc >= S_MV_P
ASPECT = 2 ** 12              # maximal ratio of normalised box side lengths


class CertError(Exception):
    pass


def A(q):
    """arb ball containing the rational q (exact if q is dyadic)."""
    return arb(q)


def lower_decimal(x, digits):
    """decimal string <= the exact float x (x an exact arb), rounded down at 'digits' significant digits."""
    m, e = x.mid().man_exp()
    m, e = int(m), int(e)
    q = fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))
    if q <= 0:
        return "nonpositive"
    k = 0
    while q * fmpq(10) ** k < 10 ** (digits - 1):
        k += 1
    while q * fmpq(10) ** k >= 10 ** digits:
        k -= 1
    num = (q * fmpq(10) ** k).floor()
    return f"{int(num)}e{-k}"


def finite(x):
    return x.is_finite()


def hull(a, b):
    return a.union(b)


def nonneg_mul(a, b):
    """enclosure of {u*v : u in a, v in b} for quantities known to be >= 0."""
    alo = a.lower(); blo = b.lower()
    if alo < 0: alo = arb(0)
    if blo < 0: blo = arb(0)
    return hull(alo * blo, a.upper() * b.upper())


# ----------------------------------------------------------------------------
# special-function building blocks (valid for ball arguments; Arb is rigorous)
# ----------------------------------------------------------------------------
def J(q, s, d):
    """J_q(s,d) = E(d + s R)^q = int_0^inf (d+sr)^q e^{-r} dr,  s >= 0, d >= 0, q >= 0."""
    if s == 0:
        return d ** q if d != 0 else arb(0)
    if d == 0:
        return s ** q * (q + 1).gamma()
    z = d / s
    return s ** q * z.exp() * z.gamma_upper(q + 1)


def Iint(x, p):
    """I(x,p) = int_0^x u^p e^u du, x >= 0."""
    if x == 0:
        return arb(0)
    return x ** (p + 1) / (p + 1) * x.hypgeom_1f1(p + 1, p + 2)


def Ip(x, p):
    """d/dp I(x,p) = int_0^x u^p log(u) e^u du, x > 0 point/ball, p ball."""
    if x == 0:
        return arb(0)
    return x.log() * Iint(x, p) - x ** (p + 1) / (p + 1) ** 2 * x.hypgeom([p + 1, p + 1], [p + 2, p + 2])


def phi(x, p):
    """phi(x,p) = E|x - R|^p, x >= 0."""
    return (-x).exp() * ((p + 1).gamma() + Iint(x, p))


def dGamma_upper_da(a, z):
    """d/da Gamma(a,z) for z >= 0 (cancellation: call with enough precision)."""
    g = a.gamma()
    if z == 0:
        return g * a.digamma()
    glow = g - z.gamma_upper(a)
    dglow = z.log() * glow - z ** a / a ** 2 * (-z).hypgeom([a, a], [a + 1, a + 1])
    return g * a.digamma() - dglow


def M_pt(s, p):
    """M(s,p) at a (tiny-ball) point s in [0,1/2], p point or ball."""
    d = 1 - 2 * s; c = 1 - s; x = d / c
    TA = s * J(p, s, d) if s != 0 else arb(0)
    return TA + c ** (p + 1) * phi(x, p)


def Mp_pt(sq, pq):
    """dM/dp at exact rationals s = sq in [0,1/2], p = pq >= 1.  The evaluation of
    dGamma(a,z)/da loses about 2.9 z bits to cancellation (z = (1-2s)/s); the inputs are
    therefore converted to balls only after the precision has been raised."""
    sf = float(sq.p) / float(sq.q)
    z = (1 - 2 * sf) / sf if sf > 0 else 0.0
    prec = BASE_PREC + 64 + int(3 * z) + 1
    old = ctx.prec
    ctx.prec = prec
    try:
        s = arb(sq); p = arb(pq)
        d = 1 - 2 * s; c = 1 - s; x = d / c
        out = arb(0)
        if s != 0:
            TA = s * J(p, s, d)
            if d == 0:
                ka = dGamma_upper_da(p + 1, arb(0))
                out += s.log() * TA + s ** (p + 1) * ka
            else:
                zz = d / s
                out += s.log() * TA + s ** (p + 1) * zz.exp() * dGamma_upper_da(p + 1, zz)
        TB = c ** (p + 1) * phi(x, p)
        out += c.log() * TB + c ** (p + 1) * (-x).exp() * ((p + 1).gamma() * (p + 1).digamma() + Ip(x, p))
    finally:
        ctx.prec = old
    out = +out   # round to working precision (outward)
    if not finite(out) or out.rad() > 1e-20 * (abs(out.mid()) + 1):
        raise CertError(f"M_p evaluation lost accuracy at s={sq}, p={pq}")
    return out


def N(S):
    return S * S + (1 - S) ** 2


def F_pt(p):
    return phi(arb(1), p)


def Fp_pt(p):
    """F'(p) = E|E-1|^p log|E-1| = e^{-1}(Gamma(p+1)psi(p+1) + I_p(1,p))."""
    one = arb(1)
    return (-one).exp() * ((p + 1).gamma() * (p + 1).digamma() + Ip(one, p))


def G_pt(S, P):
    return M_pt(S, P).log() / P - N(S).log() / 2


def G0_pt(P):
    return F_pt(P).log() / P


def Gh_pt(P):
    return (P + 1).lgamma() / P - arb(2).log() / 2


# ----------------------------------------------------------------------------
# enclosures over p-intervals [p0,p1] (exact rationals) -- Lemma B.4 / B.5
# ----------------------------------------------------------------------------
def intersect(*encs):
    """intersection of several rigorous enclosures of the same real number."""
    lo = None; hi = None
    for e in encs:
        if e is None or not finite(e):
            continue
        l, u = e.lower(), e.upper()
        lo = l if lo is None or l > lo else lo
        hi = u if hi is None or u < hi else hi
    if lo is None:
        return arb.nan()
    if lo > hi:
        raise CertError("empty intersection of enclosures (inconsistent evaluation)")
    return hull(lo, hi)


class PInt:
    """data for an exact p-interval [p0,p1], 1 <= p0 <= p1."""
    def __init__(self, p0, p1):
        rational(p0); rational(p1)
        if not (1 <= p0 <= p1):
            raise CertError("invalid p interval")
        self.p0, self.p1 = p0, p1
        self.A0, self.A1 = A(p0), A(p1)
        self.P = hull(self.A0, self.A1)
        # Gamma(p+1), psi(p+1), log Gamma(p+1) are nondecreasing for p >= 1
        self.gam = hull((self.A0 + 1).gamma(), (self.A1 + 1).gamma())
        self.psi = hull((self.A0 + 1).digamma(), (self.A1 + 1).digamma())
        self.lgam = hull((self.A0 + 1).lgamma(), (self.A1 + 1).lgamma())

    def I(self, x):
        """I(x,p) over p in [p0,p1] for a point x in [0,1]: nonincreasing in p."""
        return hull(Iint(x, self.A1), Iint(x, self.A0))

    def I_box(self, x0, x1):
        """I(x,p) over x in [x0,x1] (0<=x0<=x1<=1), p in [p0,p1]: increasing in x, nonincreasing in p."""
        return hull(Iint(x0, self.A1), Iint(x1, self.A0))


def J_int(q0, q1, s, d):
    """enclosure of J_q(s,d) for q in [q0,q1] (exact rationals, 0 <= q0 <= q1), s,d point balls.
    Intersection of: (i) the convexity bound J <= max(J(q0),J(q1)) (Lemma B.4);
    (ii) J_q = d^q Phi_q(z), Phi_q(z) = E(1+R/z)^q nondecreasing in q (Lemma B.5);
    (iii) Arb ball evaluation of s^q e^z Gamma(q+1,z) on the ball q in [q0,q1]."""
    Q0, Q1 = A(q0), A(q1)
    if q0 == q1:
        return J(Q0, s, d)
    Q = hull(Q0, Q1)
    J0, J1 = J(Q0, s, d), J(Q1, s, d)
    upper = hull(arb(0), hull(J0, J1).upper())          # [0, max]  (J > 0)
    if s == 0:                                         # J = d^q, monotone in q
        return hull(J0, J1)
    if d == 0:                                         # J = s^q Gamma(q+1)
        return intersect(upper, s ** Q * (Q + 1).gamma())
    phi0, phi1 = J0 / d ** Q0, J1 / d ** Q1
    e_phi = d ** Q * hull(phi0, phi1)
    z = d / s
    e_ball = s ** Q * z.exp() * z.gamma_upper(Q + 1) if z < 1 else None
    return intersect(upper, e_phi, e_ball)


def M_int(s, d, c, x, pi):
    """M(s,p) for a point s (balls s,d,c,x) and p in [p0,p1]."""
    TA = s * J_int(pi.p0, pi.p1, s, d) if s != 0 else arb(0)
    TB = nonneg_mul(nonneg_mul(c ** (pi.P + 1), (-x).exp()), pi.gam + pi.I(x))
    return TA + TB


def nonneg_pow(b, P):
    """b^P for an exact/tiny ball b >= 0 and P ball > 0."""
    if b == 0:
        return arb(0)
    return b ** P


def s_box_data(s0, s1, pi):
    """enclosures of M and M_s over [s0,s1] x [p0,p1] (Lemmas B.1, B.2, B.5)."""
    S0, S1 = A(s0), A(s1)
    P = pi.P
    D0, D1 = 1 - 2 * S1, 1 - 2 * S0                  # d range
    # J_q(s,d) nondecreasing in s and d: corners (s0, 1-2 s1) and (s1, 1-2 s0)
    Jp = hull(J_int(pi.p0, pi.p1, S0, D0), J_int(pi.p0, pi.p1, S1, D1))
    Jq = hull(J_int(pi.p0 - 1, pi.p1 - 1, S0, D0), J_int(pi.p0 - 1, pi.p1 - 1, S1, D1))
    TA = hull(S0 * J_int(pi.p0, pi.p1, S0, D0), S1 * J_int(pi.p0, pi.p1, S1, D1)) if s0 != 0 else \
        hull(arb(0), S1 * J_int(pi.p0, pi.p1, S1, D1))
    x0, x1 = D0 / (1 - S1), D1 / (1 - S0)             # x = d/c decreasing in s
    ex = hull((-x1).exp(), (-x0).exp())
    ph = nonneg_mul(ex, pi.gam + pi.I_box(x0, x1))     # phi(x,p) = e^{-x}(Gamma(p+1)+I(x,p))
    cp1 = hull(nonneg_pow(1 - S1, P + 1), nonneg_pow(1 - S0, P + 1))
    cpm1 = hull(nonneg_pow(1 - S1, P - 1), nonneg_pow(1 - S0, P - 1))
    C = hull(1 - S1, 1 - S0)
    dp = hull(nonneg_pow(D0, P), nonneg_pow(D1, P))
    Mb = TA + nonneg_mul(cp1, ph)
    Ms = (P + 1) * Jp - P * Jq - cpm1 * ph * ((P + 1) * C - 1) - dp / C
    return Mb, Ms


def Gs_box(s0, s1, pi):
    """enclosure of dG/ds = M_s/(pM) - (2s-1)/N(s) over the box."""
    S0, S1 = A(s0), A(s1)
    Mb, Ms = s_box_data(s0, s1, pi)
    if not (Mb.lower() > 0):
        return arb.nan()
    # (2s-1)/N(s) is nondecreasing on [0,1/2]: derivative 4s(1-s)/N^2 >= 0
    q = hull((2 * S0 - 1) / N(S0), (2 * S1 - 1) / N(S1))
    return Ms / (pi.P * Mb) - q


def M_box(s0, s1, P):
    """(kept for the self-test interface) enclosure of M over a box with P = PInt."""
    return s_box_data(s0, s1, P)[0]


def Ms_box(s0, s1, P):
    """M_s at a point (s0 == s1) or over a box; P a PInt or an arb point."""
    if not isinstance(P, PInt):
        # point evaluation of Lemma B.2 at exact p (used by the self-tests)
        S = A(s0); d = 1 - 2 * S; c = 1 - S; x = d / c
        Jp = J(P, S, d) if s0 != 0 else arb(1)
        Jq = J(P - 1, S, d) if s0 != 0 else arb(1)
        dp = d ** P if d != 0 else arb(0)
        return (P + 1) * Jp - P * Jq - c ** (P - 1) * phi(x, P) * ((P + 1) * c - 1) - dp / c
    return s_box_data(s0, s1, P)[1]


def F_int(pi):
    one = arb(1)
    return (-one).exp() * (pi.gam + pi.I(one))


def dR_box(ref, pi):
    """enclosure of R'(p) over [p0,p1], R in {G0, Gh}."""
    P = pi.P
    if ref == "0":
        Fp = hull(Fp_pt(pi.A0), Fp_pt(pi.A1))     # F' nondecreasing in p (Lemma B.4)
        F = F_int(pi)
        return Fp / (P * F) - F.log() / P ** 2
    return pi.psi / P - pi.lgam / P ** 2


def R_pt(ref, P):
    return G0_pt(P) if ref == "0" else Gh_pt(P)


def R_int(ref, pi):
    if ref == "0":
        return F_int(pi).log() / pi.P
    return pi.lgam / pi.P - arb(2).log() / 2


def box_enclosures(s0, s1, p0, p1, refs):
    """returns {ref: (enclosure of G-R on box, p-width contribution, s-width contribution, mode, centre value)}"""
    for value in (s0,s1,p0,p1):
        rational(value)
    if not (0 <= s0 <= s1 <= HALF and 1 <= p0 <= p1) or not refs or any(r not in ('0','h') for r in refs):
        raise CertError("invalid box/reference domain")
    sc = (s0 + s1) / 2; pc = (p0 + p1) / 2
    hs = A((s1 - s0) / 2); hp = A((p1 - p0) / 2)
    Sc, Pc = A(sc), A(pc)
    pi = PInt(p0, p1)
    P = pi.P
    Gs = Gs_box(s0, s1, pi)
    out = {}
    use_mv_p = sc >= S_MV_P
    Gc = G_pt(Sc, Pc)
    dc = 1 - 2 * Sc; cc = 1 - Sc; xc = dc / cc
    Mc = M_int(Sc, dc, cc, xc, pi)                    # M(sc, p), p in [p0,p1]
    if use_mv_p:
        Mp = hull(Mp_pt(sc, p0), Mp_pt(sc, p1))       # M_p nondecreasing in p (Lemma B.4)
        dG = Mp / (P * Mc) - Mc.log() / P ** 2
    else:
        GcP = Mc.log() / P - N(Sc).log() / 2
    for r in refs:
        Dc = Gc - R_pt(r, Pc)
        if use_mv_p:
            Dp = dG - dR_box(r, pi)
            encP = Dc + Dp * arb(0, hp)
            wp = abs(Dp).upper() * hp if finite(Dp) else arb.pos_inf()
        else:
            encP = GcP - R_int(r, pi)
            wp = encP.rad() if finite(encP) else arb.pos_inf()
        enc = encP + Gs * arb(0, hs)
        ws = abs(Gs).upper() * hs if finite(Gs) else arb.pos_inf()
        out[r] = (enc, wp, ws, "mv-p" if use_mv_p else "ball-p", Dc)
    return out


# ----------------------------------------------------------------------------
# branch and bound
# ----------------------------------------------------------------------------
def certify(region):
    validate_parameters()
    s0, s1, p0, p1, refs, sense = REGIONS[region]
    t0 = time.time()
    stack = [(s0, s1, p0, p1, 0)]
    nboxes = 0; nleaves = 0; fails = []; ncounter = 0
    min_lb = None; min_box = None; maxdepth = 0
    modes = {"mv-p": 0, "ball-p": 0}
    which_ref = {r: 0 for r in refs}
    leafhash = hashlib.sha256()
    ws_tot = s1 - s0; wp_tot = p1 - p0
    while stack:
        a0, a1, b0, b1, dep = stack.pop()
        nboxes += 1
        maxdepth = max(maxdepth, dep)
        res = box_enclosures(a0, a1, b0, b1, refs)
        ok = False; best = None; counterexample = True
        for r in refs:
            enc, wp, ws, mode, Dc = res[r]
            if not ((Dc if sense > 0 else -Dc).upper() < 0):
                counterexample = False
            val = enc if sense > 0 else -enc
            lb = val.lower()
            if finite(val) and lb > 0:
                ok = True
                nleaves += 1; modes[mode] += 1; which_ref[r] += 1
                if min_lb is None or lb < min_lb:
                    min_lb = lb; min_box = (a0, a1, b0, b1)
                leafhash.update(f"{a0},{a1},{b0},{b1},{r}\n".encode())
                break
            if best is None or (finite(wp) and finite(ws) and (wp + ws) < (best[0] + best[1])):
                best = (wp, ws)
        if ok:
            continue
        if counterexample:
            # the claimed inequality is certified FALSE at the box centre: stop, do not split
            fails.append((a0, a1, b0, b1))
            ncounter += 1
            break           # a certified counterexample settles the outcome: stop
        if nboxes > MAX_BOXES or ((a1 - a0) < MIN_WIDTH and (b1 - b0) < MIN_WIDTH):
            fails.append((a0, a1, b0, b1))
            break           # stop at the first undecidable box: the region is not certified
        wp, ws = best
        if finite(wp) and finite(ws):
            split_s = ws > wp
        else:
            split_s = (a1 - a0) / ws_tot > (b1 - b0) / wp_tot
        ns_, np_ = (a1 - a0) / ws_tot, (b1 - b0) / wp_tot
        if ns_ * ASPECT < np_:
            split_s = False        # keep boxes from becoming needle-shaped
        elif np_ * ASPECT < ns_:
            split_s = True
        if split_s and (a1 - a0) < MIN_WIDTH:
            split_s = False
        if (not split_s) and (b1 - b0) < MIN_WIDTH:
            split_s = True
        if split_s:
            am = (a0 + a1) / 2
            stack.append((am, a1, b0, b1, dep + 1)); stack.append((a0, am, b0, b1, dep + 1))
        else:
            bm = (b0 + b1) / 2
            stack.append((a0, a1, bm, b1, dep + 1)); stack.append((a0, a1, b0, bm, dep + 1))
    runtime = time.time() - t0
    return {
        "boxes_evaluated": nboxes, "leaves_certified": nleaves, "failures": len(fails),
        "centre_counterexamples": ncounter,
        "failure_boxes": [[str(v) for v in f] for f in fails[:50]],
        "min_certified_margin_lower_bound": (point(min_lb) if min_lb is not None else None),
        "bound_encoding": "exact rational strings",
        "min_margin_box": [str(v) for v in min_box] if min_box else None,
        "max_depth": maxdepth, "leaf_modes": modes, "leaves_by_reference": which_ref,
        "leaves_sha256": leafhash.hexdigest(), "runtime_s": round(runtime, 2),
    }


# ----------------------------------------------------------------------------
# self-tests (reference values from SPEC.md; independent mpmath quadrature)
# ----------------------------------------------------------------------------
SPEC_REFS = [((fmpq(1, 10), fmpq(13, 10)), "-0.213049905296183"),
             ((fmpq(3, 10), fmpq(19, 10)), "-0.0294149150315886"),
             ((fmpq(1, 4), fmpq(5, 2)), "0.143192764694943"),
             ((fmpq(2, 5), fmpq(5)), "0.641602497665675"),
             ((fmpq(1, 20), fmpq(12)), "1.58099627165213"),
             ((fmpq(1, 2), fmpq(6, 5)), "-0.265784034621107"),
             ((fmpq(0), fmpq(33, 10)), "0.375483768956575")]


def selftest():
    rows = []
    for (s, p), ref in SPEC_REFS:
        g = G_pt(A(s), A(p))
        diff = g - arb(ref)
        ok = abs(diff).upper() < arb("1e-14")
        rows.append({"s": str(s), "p": str(p), "G_arb": g.str(18, radius=False), "spec": ref, "ok": bool(ok)})
        if not ok:
            raise CertError(f"self-test failed at s={s}, p={p}: {g} vs {ref}")
    for nm, val, ref in [("G0(1.3)", G0_pt(A(fmpq(13, 10))), "-0.213229500079513"),
                         ("Gh(1.3)", Gh_pt(A(fmpq(13, 10))), "-0.227966317234103")]:
        ok = abs(val - arb(ref)).upper() < arb("1e-14")
        rows.append({"quantity": nm, "arb": val.str(18, radius=False), "spec": ref, "ok": bool(ok)})
        if not ok:
            raise CertError(f"self-test failed for {nm}")
    # consistency: M_s, M_p closed forms against central differences of the closed form for M
    h = fmpq(1, 10**12)
    for (s, p) in [(fmpq(1, 10), fmpq(13, 10)), (fmpq(1, 20), fmpq(39, 20)), (fmpq(9, 20), fmpq(41, 20)), (fmpq(3, 10), fmpq(7))]:
        old = ctx.prec; ctx.prec = 256
        fd_s = (M_pt(A(s + h), A(p)) - M_pt(A(s - h), A(p))) / A(2 * h)
        fd_p = (M_pt(A(s), A(p + h)) - M_pt(A(s), A(p - h))) / A(2 * h)
        ctx.prec = old
        ms = Ms_box(s, s, A(p)); mp_ = Mp_pt(s, p)
        es = abs(ms - fd_s).upper(); ep = abs(mp_ - fd_p).upper()
        ok = es < arb("1e-15") and ep < arb("1e-15")
        rows.append({"fd_check_s": str(s), "p": str(p), "Ms": ms.str(15, radius=False), "fd_s": fd_s.str(15, radius=False),
                     "Mp": mp_.str(15, radius=False), "fd_p": fd_p.str(15, radius=False), "ok": bool(ok)})
        if not ok:
            raise CertError(f"derivative self-test failed at s={s}, p={p}")
    rows.append(containment_test())
    return rows


def containment_test(nbox=24, npts=4):
    """Containment test of the box enclosures: at exact random interior points of random
    boxes, the point values of M, M_s, dG/ds, dD/dp and D = G - R must lie in the box
    enclosures.  (A necessary condition for correctness of the monotonicity logic.)"""
    import random
    rng = random.Random(20261014)
    checked = 0
    for k in range(nbox):
        s0 = fmpq(rng.randint(0, 480), 1000); ws = fmpq(rng.choice([1, 5, 20, 80]), 1000)
        s1 = min(s0 + ws, HALF)
        p0 = fmpq(rng.randint(1000, 6900), 1000); wp = fmpq(rng.choice([1, 10, 50, 200]), 1000)
        p1 = p0 + wp
        pi = PInt(p0, p1)
        Mb, Msb = s_box_data(s0, s1, pi)
        Gsb = Gs_box(s0, s1, pi)
        refs = ("0", "h")
        enc = box_enclosures(s0, s1, p0, p1, refs)
        sc = (s0 + s1) / 2
        dc = 1 - 2 * A(sc); cc = 1 - A(sc)
        Mc = M_int(A(sc), dc, cc, dc / cc, pi)
        for j in range(npts):
            s = s0 + (s1 - s0) * fmpq(rng.randint(0, 1000), 1000)
            p = p0 + (p1 - p0) * fmpq(rng.randint(0, 1000), 1000)
            S, Pp = A(s), A(p)
            m = M_pt(S, Pp); ms = Ms_box(s, s, Pp)
            gs = ms / (Pp * m) - (2 * S - 1) / N(S)
            mcp = M_pt(A(sc), Pp)
            g = m.log() / Pp - N(S).log() / 2
            ok = Mb.contains(m) and Msb.contains(ms) and Gsb.contains(gs) and Mc.contains(mcp)
            for r in refs:
                ok = ok and enc[r][0].contains(g - R_pt(r, Pp))
            if not ok:
                raise CertError(f"containment test failed: box {[str(v) for v in (s0, s1, p0, p1)]} point {s},{p}")
            checked += 1
    return {"containment_test_points": checked, "ok": True}


def validate_parameters():
    if flint.__version__ != FLINT_REQUIRED:
        raise CertError(f"python-flint {FLINT_REQUIRED} required, found {flint.__version__}")
    if not (fmpq(1682457208, 10**9) == PSTAR_LO and PSTAR_HI == fmpq(1682457209, 10**9)):
        raise CertError("bad p* bracket")
    for name, (s0, s1, p0, p1, refs, sense) in REGIONS.items():
        if not (0 <= s0 < s1 <= HALF and 1 <= p0 < p1):
            raise CertError(f"bad region {name}")
        if isinstance(sense, bool) or sense not in (+1, -1) or not refs or any(r not in ('0','h') for r in refs):
            raise CertError(f"bad region spec {name}")
    if not (ETA == SIGMA == SIGMA0 == DELTA == fmpq(1, 20)):
        raise CertError("parameters differ from SPEC")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("region")
    ap.add_argument("--out", default=None, help="certificate JSON path")
    ap.add_argument("--prec", type=int, default=128, help="Arb working precision in bits (>= 128)")
    args = ap.parse_args()
    global BASE_PREC
    if args.prec < 128 or args.prec > 1024:
        raise CertError("--prec must be in [128, 1024]")
    BASE_PREC = args.prec
    if args.region not in REGIONS:
        raise CertError(f"unknown region {args.region!r}; expected one of {sorted(REGIONS)}")
    ctx.prec = BASE_PREC
    validate_parameters()
    tst = selftest()
    ctx.prec = BASE_PREC
    s0, s1, p0, p1, refs, sense = REGIONS[args.region]
    res = certify(args.region)
    stmt = {"+1": "G(s,p) - R(p) > 0 for some R in refs, on every box", "-1": "R(p) - G(s,p) > 0 on every box"}[f"{sense:+d}"]
    summary = {
        "region": args.region,
        "s_range": [str(s0), str(s1)], "p_range": [str(p0), str(p1)],
        "references": ["G0" if r == "0" else "Gh" for r in refs], "statement": stmt,
        "parameters": {"eta": str(ETA), "sigma": str(SIGMA), "sigma0": str(SIGMA0), "delta": str(DELTA),
                       "p_max_B4": str(PMAX_B4), "pstar_lo": str(PSTAR_LO), "pstar_hi": str(PSTAR_HI),
                       "min_width": str(MIN_WIDTH), "s_threshold_mv_p": str(S_MV_P)},
        "status": "CERTIFIED" if res["failures"] == 0 else "FAILED",
        **res,
        "flint_version": flint.__version__, "arb_precision_bits": BASE_PREC,
        "python": platform.python_version(),
        "script": os.path.basename(__file__), "script_sha256": sha256_file(__file__),
        "selftest": tst,
    }
    txt = json.dumps(summary, indent=1)
    print(txt)
    if args.out:
        with open(args.out, "x") as f:
            f.write(txt + "\n")
    if res["failures"] != 0:
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except CertError as e:
        print(json.dumps({"status": "ERROR", "error": str(e)}))
        sys.exit(2)
