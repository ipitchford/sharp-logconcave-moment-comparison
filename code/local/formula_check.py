#!/usr/bin/env python3
"""Independent (non-rigorous, high-precision) checks of the LOCAL-worker formulas.

Run:  uv run --no-project --with python-flint==0.9.0 --with mpmath python3 formula_check.py [out.json]

Checks performed (all in mpmath, 50 digits):
  C1  closed form of M(s,p) (Gamma / 1F1 / Tricomi U) vs direct quadrature of E|Z_s|^p,
      and the SPEC reference values of G.
  C2  SPEC derivative formulas M_s, M_ss (conditioning integrals, by quadrature)
      vs finite differences of M (quadrature-based, independent of the closed forms).
  C3  reformulated closed forms  G_s = A,
      G_ss = (p-1) B - p A^2 + 2 (d/N) A - 1/N^2
      vs finite differences of G (quadrature-based).
  C4  s = 0 (central FD; Z_s for s < 0 via the hypoexponential density):   G_s = 0, (1/2) G_ss = ((p-1)/(2p)) (1/F - 1);
      s = 1/2: G_s = 0, (1/2) G_ss = 2 (p-2)^2 / p.
  C5  Arb evaluator of local_cert.py at points vs mpmath closed forms.
"""
import json
import sys

import mpmath as mp

mp.mp.dps = 50


# ----------------------------------------------------------------- quadrature reference
def M_quad(s, p):
    """E|Z_s|^p by direct quadrature over the density of X = sE - (1-s)E' (valid 0<=s<=1)."""
    s = mp.mpf(s); p = mp.mpf(p); c = 1 - s; d = 1 - 2 * s
    if s == 0:
        return mp.quad(lambda r: abs(1 - r) ** p * mp.exp(-r), [0, 1, mp.inf])
    if s < 0:
        # X = s E - c E' = -(|s| E + c E'), |s|E + cE' has the hypoexponential density
        a = -s
        f = lambda t: (mp.exp(-t / c) - mp.exp(-t / a)) / (c - a)
        return mp.quad(lambda t: abs(d - t) ** p * f(t), [0, d, mp.inf])
    if c == 0:
        return mp.quad(lambda r: abs(-1 + r) ** p * mp.exp(-r), [0, 1, mp.inf])
    # Z = d + X, X has density e^{-x/s} (x>0), e^{x/c} (x<0)
    pos = mp.quad(lambda x: abs(d + x) ** p * mp.exp(-x / s), [0, max(0, -d), mp.inf] if d < 0 else [0, mp.inf])
    neg = mp.quad(lambda x: abs(d + x) ** p * mp.exp(x / c), [-mp.inf, -d, 0] if d > 0 else [-mp.inf, 0])
    return pos + neg


def G_quad(s, p):
    s = mp.mpf(s); p = mp.mpf(p)
    return mp.log(M_quad(s, p)) / p - mp.log(s * s + (1 - s) ** 2) / 2


def spec_derivs_quad(s, p):
    """SPEC formulas for M_s, M_ss via quadrature of the conditioning integrals (0<s<1/2 or s=0)."""
    s = mp.mpf(s); p = mp.mpf(p); c = 1 - s; d = 1 - 2 * s; x = d / c
    h = lambda r: (r - 1) ** 2 + 1
    if s > 0:
        Ps = mp.quad(lambda r: (d + s * r) ** (p - 1) * (r - 1) * mp.exp(-r), [0, mp.inf])
        Pss = mp.quad(lambda r: (d + s * r) ** (p - 2) * h(r) * mp.exp(-r), [0, mp.inf])
    else:
        Ps = mp.mpf(0); Pss = mp.mpf(0)
    Ns = mp.quad(lambda r: abs(d - c * r) ** (p - 2) * (d - c * r) * (r - 1) * mp.exp(-r), [0, x, mp.inf])
    # |d - cR|^(p-2) is singular at R = x; substitute u = |d - cR| = v^(1/(p-1)), so that
    # int_0^U u^(p-2) g du = (1/(p-1)) int_0^(U^(p-1)) g(v^(1/(p-1))) dv  (smooth integrand)
    e = p - 1
    g = lambda r: h(r) * mp.exp(-r)
    inner = mp.quad(lambda v: g((d - v ** (1 / e)) / c), [0, d ** e]) / (c * e)
    outer = mp.quad(lambda v: g((d + v ** (1 / e)) / c), [0, 1]) / (c * e) \
        + mp.quad(lambda u: u ** (p - 2) * g((d + u) / c), [1, mp.inf]) / c
    Nss = inner + outer
    Ms = p * (s * Ps + c * Ns)
    Mss = p * (p - 1) * (s * Pss + c * Nss)
    return Ms, Mss


# ----------------------------------------------------------------- closed forms (mpmath)
def J(a, x):
    """int_0^x u^a e^u du  (a > -1)."""
    if x == 0:
        return mp.mpf(0)
    return x ** (a + 1) / (a + 1) * mp.hyp1f1(a + 1, a + 2, x)


def P_U(q, j, d, s):
    """E[(d+sR)^q R^j], Tricomi-U form; s=0 gives d^q j!."""
    if s == 0:
        return d ** q * mp.factorial(j)
    y = d / s
    return d ** q * y ** (j + 1) * mp.factorial(j) * mp.hyperu(j + 1, q + j + 2, y)


def P_G(q, coeffs, d, s):
    """E[(d+sR)^q f(R)] with f(R) = sum_k coeffs[k] (R + y)^k  (v = y + R), Gamma form."""
    y = d / s
    return s ** q * mp.exp(y) * sum(ck * mp.gammainc(q + 1 + k, y) for k, ck in enumerate(coeffs))


def closed(s, p):
    """Return dict with M, a, b, A, B, Gs, Gss from the closed forms."""
    s = mp.mpf(s); p = mp.mpf(p)
    c = 1 - s; d = 1 - 2 * s; x = d / c; N = s * s + c * c
    ex = mp.exp(-x)
    g = mp.gamma
    # plus side
    if s == 0:
        Pp = Pa = Pb = mp.mpf(0)
        Pa_raw = d ** (p - 1) * (1 - 0)  # E[(d+sR)^{p-1}(R-2s)] at s=0
    else:
        Pp = P_U(p, 0, d, s)
        Pa_raw = P_U(p - 1, 1, d, s) - 2 * s * P_U(p - 1, 0, d, s)
        Pb = (c * c / N ** 2) * P_U(p - 2, 2, d, s) - (4 * s * c * c / N ** 2) * P_U(p - 2, 1, d, s) \
            + (4 * s * s * c * c / N ** 2 + 1) * P_U(p - 2, 0, d, s)
    M = s * Pp + c * (c ** p * ex * J(p, x) + c ** p * ex * g(p + 1))
    a = Pa_raw - c ** (p - 1) * ex * ((N / c) * J(p - 1, x) + J(p, x)) \
        - c ** (p - 1) * ex * (g(p + 1) - (N / c) * g(p))
    Qb = c ** (p - 2) * ex * ((N / c ** 2) * J(p - 2, x) + (2 * s * s / (c * N)) * J(p - 1, x) + (s * s / N ** 2) * J(p, x))
    Ob = c ** (p - 2) * ex * ((N / c ** 2) * g(p - 1) - (2 * s * s / (c * N)) * g(p) + (s * s / N ** 2) * g(p + 1))
    b = (s * Pb if s > 0 else 0) + c * (Qb + Ob)
    A = (s * c / N) * a / M
    B = b / M
    Gss = (p - 1) * B - p * A ** 2 + 2 * (d / N) * A - 1 / N ** 2
    return dict(M=M, a=a, b=b, A=A, B=B, Gs=A, Gss=Gss)


def closed_gamma_form(s, p):
    """Same quantities, plus side in Gamma(.,y) form with cancelled singular part (L1 form)."""
    s = mp.mpf(s); p = mp.mpf(p)
    c = 1 - s; d = 1 - 2 * s; x = d / c; y = d / s; N = s * s + c * c
    ex = mp.exp(-x); ey = mp.exp(y); g = mp.gamma; G = mp.gammainc
    M = s ** (p + 1) * ey * G(p + 1, y) + c ** (p + 1) * ex * (J(p, x) + g(p + 1))
    a = s ** (p - 1) * ey * (G(p + 1, y) - (N / s) * G(p, y)) \
        - c ** (p - 1) * ex * ((N / c) * J(p - 1, x) + J(p, x)) - c ** (p - 1) * ex * (g(p + 1) - (N / c) * g(p))
    # (p-1) b, cancelled
    D = d ** (p - 1)
    F1 = lambda z: mp.hyp1f1(p - 1, p, z)
    S = D * (ex * (N / c ** 2) * F1(x) - ey * (N / s ** 2) * F1(-y))
    plus = s ** (p - 1) * ey * ((p - 1) * (c * c / N ** 2) * G(p + 1, y) - (p - 1) * (2 * c * c / (s * N)) * G(p, y)
                                + (N / s ** 2) * g(p))
    inner = c ** (p - 1) * ex * ((p - 1) * (2 * s * s / (c * N)) * J(p - 1, x) + (p - 1) * (s * s / N ** 2) * J(p, x))
    outer = c ** (p - 1) * ex * ((N / c ** 2) * g(p) - (p - 1) * (2 * s * s / (c * N)) * g(p) + (p - 1) * (s * s / N ** 2) * g(p + 1))
    pm1b = S + plus + inner + outer
    A = (s * c / N) * a / M
    Gss = pm1b / M - p * A ** 2 + 2 * (d / N) * A - 1 / N ** 2
    return dict(M=M, a=a, pm1b=pm1b, A=A, Gss=Gss)


def F(p):
    p = mp.mpf(p)
    return mp.e ** -1 * (mp.gamma(p + 1) + mp.quad(lambda u: u ** p * mp.exp(u), [0, 1]))


def main():
    out = {"dps": mp.mp.dps, "checks": {}}
    worst = {}

    def rec(name, val):
        worst[name] = max(worst.get(name, 0), float(val))

    # C1: M closed form vs quadrature, and SPEC reference values
    ref = [(0.1, 1.3, '-0.213049905296183'), (0.3, 1.9, '-0.0294149150315886'), (0.25, 2.5, '0.143192764694943'),
           (0.4, 5, '0.641602497665675'), (0.05, 12, '1.58099627165213'), (0.5, 1.2, '-0.265784034621107'),
           (0, 3.3, '0.375483768956575')]
    c1 = []
    for s, p, gref in ref:
        s = mp.mpf(str(s)); p = mp.mpf(str(p))
        Mc = closed(s, p)['M']
        Mq = M_quad(s, p)
        Gc = mp.log(Mc) / p - mp.log(s * s + (1 - s) ** 2) / 2
        c1.append(dict(s=float(s), p=float(p), M_closed=mp.nstr(Mc, 20), M_quad=mp.nstr(Mq, 20),
                       relerr=mp.nstr(abs(Mc / Mq - 1), 3), G=mp.nstr(Gc, 16), G_ref=gref,
                       G_abs_diff=mp.nstr(abs(Gc - mp.mpf(gref)), 3)))
        rec("C1_M_relerr", abs(Mc / Mq - 1)); rec("C1_G_vs_ref", abs(Gc - mp.mpf(gref)))
    G0 = mp.log(F(1.3)) / mp.mpf(1.3); Gh = mp.log(mp.gamma(mp.mpf(2.3))) / mp.mpf(1.3) - mp.log(2) / 2
    c1.append(dict(G0_13=mp.nstr(G0, 16), G0_ref='-0.213229500079513', Gh_13=mp.nstr(Gh, 16), Gh_ref='-0.227966317234103'))
    rec("C1_G_vs_ref", abs(G0 - mp.mpf('-0.213229500079513'))); rec("C1_G_vs_ref", abs(Gh - mp.mpf('-0.227966317234103')))
    out["checks"]["C1_M_closed_vs_quad_and_refs"] = c1

    # C2 and C3 at interior points
    pts = [(s, p) for s in ['0.001', '0.02', '0.05', '0.3', '0.45', '0.49'] for p in ['1.05', '1.9', '2.5', '10']]
    c2 = []; c3 = []
    for s, p in pts:
        s = mp.mpf(s); p = mp.mpf(p)
        Ms_spec, Mss_spec = spec_derivs_quad(s, p)
        fM = lambda t: M_quad(t, p)
        Ms_fd = mp.diff(fM, s, 1)
        Mss_fd = mp.diff(fM, s, 2)
        c2.append(dict(s=float(s), p=float(p), Ms_spec=mp.nstr(Ms_spec, 15), Ms_fd=mp.nstr(Ms_fd, 15),
                       Mss_spec=mp.nstr(Mss_spec, 15), Mss_fd=mp.nstr(Mss_fd, 15),
                       Ms_relerr=mp.nstr(abs(Ms_spec - Ms_fd) / abs(Ms_fd), 3),
                       Mss_relerr=mp.nstr(abs(Mss_spec - Mss_fd) / abs(Mss_fd), 3)))
        rec("C2_Ms_relerr", abs(Ms_spec - Ms_fd) / abs(Ms_fd)); rec("C2_Mss_relerr", abs(Mss_spec - Mss_fd) / abs(Mss_fd))
        cl = closed(s, p)
        fG = lambda t: G_quad(t, p)
        Gs_fd = mp.diff(fG, s, 1); Gss_fd = mp.diff(fG, s, 2)
        # SPEC combination
        Mq = M_quad(s, p); N = s * s + (1 - s) ** 2
        Gss_spec = (Mss_spec * Mq - Ms_spec ** 2) / (p * Mq ** 2) - (4 * N - (4 * s - 2) ** 2) / (2 * N ** 2)
        row = dict(s=float(s), p=float(p), Gs_closed=mp.nstr(cl['Gs'], 15), Gs_fd=mp.nstr(Gs_fd, 15),
                   Gss_closed=mp.nstr(cl['Gss'], 15), Gss_fd=mp.nstr(Gss_fd, 15), Gss_spec=mp.nstr(Gss_spec, 15),
                   Gs_abserr=mp.nstr(abs(cl['Gs'] - Gs_fd), 3), Gss_abserr=mp.nstr(abs(cl['Gss'] - Gss_fd), 3),
                   Gss_spec_vs_closed=mp.nstr(abs(cl['Gss'] - Gss_spec), 3))
        rec("C3_Gs_abserr", abs(cl['Gs'] - Gs_fd)); rec("C3_Gss_abserr", abs(cl['Gss'] - Gss_fd))
        rec("C3_Gss_spec_vs_closed", abs(cl['Gss'] - Gss_spec))
        if s > mp.mpf('0.2'):
            cg = closed_gamma_form(s, p)
            row["Gss_gammaform"] = mp.nstr(cg['Gss'], 15)
            row["gammaform_vs_U"] = mp.nstr(abs(cg['Gss'] - cl['Gss']), 3)
            rec("C3_gammaform_vs_Uform", abs(cg['Gss'] - cl['Gss']))
        c3.append(row)
    out["checks"]["C2_spec_Ms_Mss_vs_fd"] = c2
    out["checks"]["C3_closed_Gs_Gss_vs_fd"] = c3

    # C4: s=0 (one-sided FD) and s=1/2 (central FD, G symmetric)
    c4 = []
    for p in ['1.05', '1.3', '1.6824572', '1.9', '2.5', '5', '10', '12']:
        p = mp.mpf(p)
        fG = lambda t: G_quad(t, p)
        Gs0 = mp.diff(fG, 0, 1)      # central differences; M_quad covers s < 0 as well
        Gss0 = mp.diff(fG, 0, 2)
        cl0 = closed(0, p)
        pred0 = (p - 1) / p * (1 / F(p) - 1)
        Gsh = mp.diff(fG, mp.mpf(1) / 2, 1)
        Gssh = mp.diff(fG, mp.mpf(1) / 2, 2)
        clh = closed_gamma_form(mp.mpf(1) / 2, p)
        predh = 4 * (p - 2) ** 2 / p
        c4.append(dict(p=float(p), Gs0_fd=mp.nstr(Gs0, 5), Gss0_fd=mp.nstr(Gss0, 15), Gss0_closed=mp.nstr(cl0['Gss'], 15),
                       Gss0_pred=mp.nstr(pred0, 15), Gsh_fd=mp.nstr(Gsh, 5), Gssh_fd=mp.nstr(Gssh, 15),
                       Gssh_closed=mp.nstr(clh['Gss'], 15), Gssh_pred=mp.nstr(predh, 15)))
        rec("C4_Gs0_fd", abs(Gs0)); rec("C4_Gsh_fd", abs(Gsh))
        rec("C4_Gss0_fd_vs_pred", abs(Gss0 - pred0)); rec("C4_Gss0_closed_vs_pred", abs(cl0['Gss'] - pred0))
        rec("C4_Gssh_fd_vs_pred", abs(Gssh - predh)); rec("C4_Gssh_closed_vs_pred", abs(clh['Gss'] - predh))
    out["checks"]["C4_endpoints"] = c4

    # C5: Arb evaluator (local_cert.py) at points vs mpmath closed forms
    try:
        import os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import local_cert as lc
        from fractions import Fraction as Fr
        c5 = []
        for s, p in [('1/1000', '19/10'), ('1/50', '33/20'), ('1/20', '39/20'), ('1/20', '41/20'), ('1/40', '6'),
                     ('1/20', '12'), ('0', '19/10'), ('0', '12'), ('9/20', '1'), ('9/20', '21/20'), ('19/40', '3/2'),
                     ('1/2', '7/4'), ('1/2', '21/20'), ('499/1000', '1001/1000'), ('9/20', '7/4')]:
            sF = Fr(s); pF = Fr(p)
            mode = 'L1' if sF >= Fr(2, 5) else 'small_s'
            v = lc.gss_box(sF, sF, pF, pF, mode)
            sm = mp.mpf(sF.numerator) / sF.denominator; pm = mp.mpf(pF.numerator) / pF.denominator
            if pF == 1:
                ref_v = closed_gamma_form(sm, pm + mp.mpf('1e-40'))['Gss'] if sm < mp.mpf(1) / 2 else mp.mpf(4)
            elif mode == 'L1':
                ref_v = closed_gamma_form(sm, pm)['Gss']
            else:
                ref_v = closed(sm, pm)['Gss']
            mid = mp.mpf(v.mid().str(30, radius=False))
            rad = mp.mpf(v.rad().str(5, radius=False))
            ok = abs(mid - ref_v) <= rad + mp.mpf('1e-25')
            c5.append(dict(s=s, p=p, arb=v.str(15), mpmath=mp.nstr(ref_v, 20), contained=bool(ok)))
            rec("C5_arb_not_containing", 0 if ok else 1)
        out["checks"]["C5_arb_point_vs_mpmath"] = c5
    except ImportError as e:
        out["checks"]["C5_arb_point_vs_mpmath"] = "skipped: %s" % e
    out["worst"] = {k: "%.3e" % v for k, v in worst.items()}
    js = json.dumps(out, indent=1)
    if len(sys.argv) > 1:
        open(sys.argv[1], 'w').write(js)
    print(json.dumps(out["worst"], indent=1))


if __name__ == "__main__":
    main()
