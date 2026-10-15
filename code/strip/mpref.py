"""Independent mpmath reference for the STRIP worker (non-rigorous; used for validation)."""
import mpmath as mp
mp.mp.dps = 40


def E(s, g):
    """E[g(Z, r)] where, with prob s, Z = d + s r and with prob c, Z = d - c r; r ~ Exp(1)."""
    s = mp.mpf(s); d = 1 - 2*s; c = 1 - s
    pos = 0
    if s > 0:
        pos = s*mp.quad(lambda r: g(d + s*r, r)*mp.e**(-r), [0, 1, 5, 20, mp.inf])
    r0 = d/c
    if r0 > 0:
        neg = c*mp.quad(lambda r: g(d - c*r, r)*mp.e**(-r), [0, r0, r0 + 1, 20, mp.inf])
    else:
        neg = c*mp.quad(lambda r: g(d - c*r, r)*mp.e**(-r), [0, 1, 20, mp.inf])
    return pos + neg


def ab(z):
    return abs(z)


def moments(s, p):
    p = mp.mpf(p)
    M = E(s, lambda z, r: ab(z)**p)
    Mp = E(s, lambda z, r: ab(z)**p*mp.log(ab(z)))
    A1 = E(s, lambda z, r: ab(z)**(p-2)*z*(r-1))
    A1L = E(s, lambda z, r: ab(z)**(p-2)*z*(r-1)*mp.log(ab(z)))
    B1 = E(s, lambda z, r: ab(z)**(p-2)*((r-1)**2 + 1))
    B1L = E(s, lambda z, r: ab(z)**(p-2)*((r-1)**2 + 1)*mp.log(ab(z)))
    Ms = p*A1
    Msp = A1 + p*A1L
    Mss = p*(p-1)*B1
    Mssp = (2*p-1)*B1 + p*(p-1)*B1L
    return dict(M=M, Mp=Mp, Ms=Ms, Msp=Msp, Mss=Mss, Mssp=Mssp)


def N(s):
    s = mp.mpf(s); return s**2 + (1-s)**2


def G(s, p):
    p = mp.mpf(p)
    return mp.log(E(s, lambda z, r: ab(z)**p))/p - mp.log(N(s))/2


def derived(s, p):
    s = mp.mpf(s); p = mp.mpf(p)
    m = moments(s, p)
    M, Mp, Ms, Msp, Mss, Mssp = (m[k] for k in ['M', 'Mp', 'Ms', 'Msp', 'Mss', 'Mssp'])
    Nn = N(s); N1 = 4*s - 2; N2 = 4
    Gp = -mp.log(M)/p**2 + Mp/(p*M)
    Gs = Ms/(p*M) - N1/(2*Nn)
    Gss = (Mss/M - (Ms/M)**2)/p - (N2/Nn - (N1/Nn)**2)/2
    Gpss = -(Mss/M - (Ms/M)**2)/p**2 + (Mssp/M - Mss*Mp/M**2 - 2*(Ms/M)*(Msp/M - Ms*Mp/M**2))/p
    return dict(m, Gp=Gp, Gs=Gs, Gss=Gss, Gpss=Gpss)


if __name__ == '__main__':
    import sys
    # reference values check
    for s, p, v in [(0.1, 1.3, '-0.213049905296183'), (0.3, 1.9, '-0.0294149150315886'),
                    (0.25, 2.5, '0.143192764694943'), (0.4, 5, '0.641602497665675'),
                    (0.05, 12, '1.58099627165213'), (0.5, 1.2, '-0.265784034621107'),
                    (0, 3.3, '0.375483768956575')]:
        g = G(mp.mpf(s), mp.mpf(p))
        print(s, p, mp.nstr(g, 16), v, mp.nstr(g - mp.mpf(v), 3))
