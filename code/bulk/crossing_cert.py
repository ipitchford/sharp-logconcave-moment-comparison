#!/usr/bin/env python3
"""
crossing_cert.py -- rigorous Arb certificate of the crossing lemma X (cycle 42, MRTT).

H(p) = G0(p) - Gh(p),  G0(p) = log F(p)/p,  F(p) = E|E-1|^p = e^{-1}(Gamma(p+1) + I(1,p)),
I(1,p) = int_0^1 u^p e^u du = 1F1(p+1;p+2;1)/(p+1),   Gh(p) = log Gamma(p+1)/p - log(2)/2.

Certifies, with PSTAR_LO = 1682457208/10^9 and PSTAR_HI = 1682457209/10^9:
  X1:  H(p) > 0  for p in [1, PSTAR_LO]
  X2:  H(p) < 0  for p in [PSTAR_HI, 2 - 1/20]
  X3:  H'(p) < 0 for p in [PSTAR_LO, PSTAR_HI]
Consequently H has exactly one zero p* in [1, 2-1/20], and PSTAR_LO < p* < PSTAR_HI;
moreover Gh - G0 = -H > 0 on [p* + 1/20, 2 - 1/20].

Enclosures: H(P) is enclosed by the mean-value form H(pc) + H'(P)[-h,h], where
H'(p) = F'(p)/(pF(p)) - log F(p)/p^2 - psi(p+1)/p + log Gamma(p+1)/p^2 and F'(P) is
enclosed by [F'(p0), F'(p1)] (F' = E|E-1|^p log|E-1| is nondecreasing in p since
F'' = E|E-1|^p log^2|E-1| >= 0).  F'(p) = e^{-1}(Gamma(p+1)psi(p+1) + I_p(1,p)),
I_p(1,p) = int_0^1 u^p log(u) e^u du = -2F2(p+1,p+1;p+2,p+2;1)/(p+1)^2.

Usage: crossing_cert.py [--out FILE].  Exit status 0 iff all three parts are certified.
"""
import sys, time, json, hashlib, os, platform, argparse
import flint
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from certificate_io import point, interval, positive_int, rational, runtime
from flint import arb, ctx, fmpq

FLINT_REQUIRED = "0.9.0"
PREC = 128
PSTAR_LO = fmpq(1682457208, 10**9)
PSTAR_HI = fmpq(1682457209, 10**9)
DELTA = fmpq(1, 20)
MIN_WIDTH = fmpq(1, 2**60)


class CertError(Exception):
    pass


def F(p):
    one = arb(1)
    return (-one).exp() * ((p + 1).gamma() + one.hypgeom_1f1(p + 1, p + 2) / (p + 1))


def Fp(p):
    one = arb(1)
    Ip = -one.hypgeom([p + 1, p + 1], [p + 2, p + 2]) / (p + 1) ** 2
    return (-one).exp() * ((p + 1).gamma() * (p + 1).digamma() + Ip)


def H(p):
    return F(p).log() / p - (p + 1).lgamma() / p + arb(2).log() / 2


def dH_box(p0, p1):
    P = arb(p0).union(arb(p1))
    Fpb = Fp(arb(p0)).union(Fp(arb(p1)))          # F' nondecreasing in p
    Fb = F(P)
    return Fpb / (P * Fb) - Fb.log() / P ** 2 - (P + 1).digamma() / P + (P + 1).lgamma() / P ** 2


def sign_cert(p0, p1, sign):
    """certify sign*H > 0 on [p0,p1] by adaptive bisection with the mean-value form."""
    rational(p0); rational(p1)
    if isinstance(sign, bool) or sign not in (-1, 1) or not (1 <= p0 < p1 < 2):
        raise CertError("invalid crossing sign/domain")
    stack = [(p0, p1)]; n = 0; fails = []; minlb = None; minbox = None; leaves = 0
    while stack:
        a, b = stack.pop(); n += 1
        c = (a + b) / 2; h = arb((b - a) / 2)
        enc = H(arb(c)) + dH_box(a, b) * arb(0, h)
        val = sign * enc
        if val.is_finite() and val.lower() > 0:
            leaves += 1
            if minlb is None or val.lower() < minlb:
                minlb = val.lower(); minbox = (a, b)
            continue
        if b - a < MIN_WIDTH:
            fails.append((a, b)); continue
        stack.append((c, b)); stack.append((a, c))
    return {"interval": [str(p0), str(p1)], "claim": ("H > 0" if sign > 0 else "H < 0"),
            "boxes_evaluated": n, "leaves": leaves, "failures": len(fails),
            "failure_boxes": [[str(a), str(b)] for a, b in fails[:20]],
            "min_certified_margin": point(minlb) if minlb is not None else None,
            "min_margin_box": [str(v) for v in minbox] if minbox else None}


def deriv_cert(p0, p1):
    rational(p0); rational(p1)
    if not (1 <= p0 < p1 < 2):
        raise CertError("invalid derivative domain")
    stack = [(p0, p1)]; n = 0; fails = []; maxub = None
    while stack:
        a, b = stack.pop(); n += 1
        d = dH_box(a, b)
        if d.is_finite() and d.upper() < 0:
            if maxub is None or d.upper() > maxub:
                maxub = d.upper()
            continue
        if b - a < MIN_WIDTH:
            fails.append((a, b)); continue
        c = (a + b) / 2
        stack.append((c, b)); stack.append((a, c))
    return {"interval": [str(p0), str(p1)], "claim": "H' < 0", "boxes_evaluated": n, "failures": len(fails),
            "max_certified_upper_bound_of_H'": point(maxub) if maxub is not None else None,
            "enclosure_on_whole_interval": dH_box(p0, p1).str(15)}


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if flint.__version__ != FLINT_REQUIRED:
        raise CertError(f"python-flint {FLINT_REQUIRED} required, found {flint.__version__}")
    if not (1 < PSTAR_LO < PSTAR_HI < 2 - DELTA):
        raise CertError("bad bracket")
    ctx.prec = PREC
    # self-test against SPEC values: G0(1.3), Gh(1.3)
    g0 = F(arb(fmpq(13, 10))).log() / arb(fmpq(13, 10))
    gh = (arb(fmpq(23, 10))).lgamma() / arb(fmpq(13, 10)) - arb(2).log() / 2
    if not (abs(g0 - arb("-0.213229500079513")).upper() < arb("1e-14") and abs(gh - arb("-0.227966317234103")).upper() < arb("1e-14")):
        raise CertError("self-test failed")
    # self-test of F' against a central difference of F
    hh = fmpq(1, 10**12); p = fmpq(17, 10)
    ctx.prec = 256
    fd = (F(arb(p + hh)) - F(arb(p - hh))) / arb(2 * hh)
    ctx.prec = PREC
    if not abs(Fp(arb(p)) - fd).upper() < arb("1e-15"):
        raise CertError("F' self-test failed")
    t0 = time.time()
    x1 = sign_cert(fmpq(1), PSTAR_LO, +1)
    x2 = sign_cert(PSTAR_HI, 2 - DELTA, -1)
    x3 = deriv_cert(PSTAR_LO, PSTAR_HI)
    ok = x1["failures"] == 0 and x2["failures"] == 0 and x3["failures"] == 0
    # informational: high-precision location of p* by Arb-checked bisection of the sign of H
    a, b = PSTAR_LO, PSTAR_HI
    for _ in range(40):
        c = (a + b) / 2
        hc = H(arb(c))
        if hc.lower() > 0: a = c
        elif hc.upper() < 0: b = c
        else: break
    out = {"lemma": "X (crossing)", "status": "CERTIFIED" if ok else "FAILED",
           "pstar_lo": str(PSTAR_LO), "pstar_hi": str(PSTAR_HI),
           "pstar_lo_decimal": "1.682457208", "pstar_hi_decimal": "1.682457209",
           "X1": x1, "X2": x2, "X3": x3,
           "pstar_refined_bracket": [str(a), str(b)],
           "bound_encoding": "exact rational strings",
           "consequence": "H has a unique zero p* in [1, 39/20]; PSTAR_LO < p* < PSTAR_HI; Gh - G0 > 0 on [p*+1/20, 39/20]",
           "runtime_s": round(time.time() - t0, 2), "flint_version": flint.__version__, "arb_precision_bits": PREC,
           "python": platform.python_version(), "script": os.path.basename(__file__), "script_sha256": sha256_file(__file__)}
    txt = json.dumps(out, indent=1)
    print(txt)
    if args.out:
        with open(args.out, "x") as f:
            f.write(txt + "\n")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except CertError as e:
        print(json.dumps({"status": "ERROR", "error": str(e)}))
        sys.exit(2)
