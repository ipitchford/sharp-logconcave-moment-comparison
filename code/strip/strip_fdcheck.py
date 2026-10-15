"""Finite-difference validation of the derivative formulas (mpmath, non-rigorous)."""
import mpmath as mp
from mpref import E, moments, derived, G, N
mp.mp.dps = 40


def Mdirect(s, p):
    return E(s, lambda z, r: abs(z)**mp.mpf(p))


def Mdouble(s, p):
    """Independent 2D quadrature of E|s(E-1)-(1-s)(E'-1)|^p (no conditioning)."""
    s = mp.mpf(s); p = mp.mpf(p); c = 1 - s
    def inner(y):
        if s == 0:
            return abs(c*(y-1))**p
        k = 1 + c*(y-1)/s
        pts = [0, k, mp.inf] if k > 0 else [0, mp.inf]
        return mp.quad(lambda x: abs(s*(x-1) - c*(y-1))**p*mp.e**(-x), pts)
    return mp.quad(lambda y: inner(y)*mp.e**(-y), [0, 1, 1/c if c else 1, 10, mp.inf])


def ds1(f, s, h):
    """first s-derivative, 4th-order; one-sided at the endpoints of [0,1/2]."""
    if s - 2*h < 0:
        return (-25*f(s) + 48*f(s+h) - 36*f(s+2*h) + 16*f(s+3*h) - 3*f(s+4*h))/(12*h)
    if s + 2*h > mp.mpf(1)/2:
        return -(-25*f(s) + 48*f(s-h) - 36*f(s-2*h) + 16*f(s-3*h) - 3*f(s-4*h))/(12*h)
    return (f(s-2*h) - 8*f(s-h) + 8*f(s+h) - f(s+2*h))/(12*h)


def ds2(f, s, h):
    if s - 2*h < 0:
        return (45*f(s) - 154*f(s+h) + 214*f(s+2*h) - 156*f(s+3*h) + 61*f(s+4*h) - 10*f(s+5*h))/(12*h**2)
    if s + 2*h > mp.mpf(1)/2:
        return (45*f(s) - 154*f(s-h) + 214*f(s-2*h) - 156*f(s-3*h) + 61*f(s-4*h) - 10*f(s-5*h))/(12*h**2)
    return (-f(s-2*h) + 16*f(s-h) - 30*f(s) + 16*f(s+h) - f(s+2*h))/(12*h**2)


def dp1(f, p, h):
    return (f(p-2*h) - 8*f(p-h) + 8*f(p+h) - f(p+2*h))/(12*h)


if __name__ == '__main__':
    worst = {}
    hs = mp.mpf('1e-5'); hp = mp.mpf('1e-5')
    for s in ['0', '0.02', '0.2', '0.5']:
        for p in ['1.96', '2', '2.04']:
            s_ = mp.mpf(s); p_ = mp.mpf(p)
            m = moments(s_, p_)
            dd = derived(s_, p_)
            fd = {}
            fd['Mp'] = dp1(lambda q: Mdirect(s_, q), p_, hp)
            fd['Ms'] = ds1(lambda t: Mdirect(t, p_), s_, hs)
            fd['Mss'] = ds2(lambda t: Mdirect(t, p_), s_, mp.mpf('1e-4'))
            fd['Msp'] = dp1(lambda q: moments(s_, q)['Ms'], p_, hp)
            fd['Mssp'] = dp1(lambda q: moments(s_, q)['Mss'], p_, hp)
            fd['Gp'] = dp1(lambda q: G(s_, q), p_, hp)
            fd['Gs'] = ds1(lambda t: G(t, p_), s_, hs)
            fd['Gss'] = ds2(lambda t: G(t, p_), s_, mp.mpf('1e-4'))
            # fully independent mixed third derivative of G (nested finite differences of G only)
            fd['Gpss'] = dp1(lambda q: ds2(lambda t: G(t, q), s_, mp.mpf('2e-3')), p_, mp.mpf('1e-3'))
            line = []
            for k, v in fd.items():
                err = abs(dd[k] - v)
                worst[k] = max(worst.get(k, 0), err)
                line.append(f"{k}={mp.nstr(dd[k], 12)} (fd err {mp.nstr(err, 2)})")
            print(f"s={s} p={p}: " + "; ".join(line), flush=True)
    print("worst abs discrepancies:", {k: mp.nstr(v, 3) for k, v in worst.items()})
    # independent 2D check of M itself
    mp.mp.dps = 20
    for s, p in [('0.02', '1.96'), ('0.2', '2.04'), ('0.5', '1.96'), ('0.37', '2.01')]:
        a = Mdouble(s, p); b = Mdirect(s, p)
        print('2D check', s, p, mp.nstr(a, 15), mp.nstr(b, 15), mp.nstr(a - b, 3))
