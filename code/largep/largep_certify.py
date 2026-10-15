#!/usr/bin/env python3
"""LARGE-P certifier for the MRTT L_p-L_2 programme (cycle 42), lemma LP in largep.md.

Run:
  uv run --no-project --with python-flint==0.9.0 python3 largep_certify.py [--out cert.json]

The proof of Lemma LP (largep.md) is analytic.  This script does two things.

PART A (part of the proof).  Verifies, with exact rationals and one Arb
enclosure of e, every numerical constant used in Step 7 of the proof:
    e < 87/32,   e - 2 < 23/32,   Gamma(7) = 720,
    11/2 - 2*(87/32) = 1/16 > 0,
    720/16 = 45,   K_lin := 45 - 87/64 > 0,
    7*K_lin - (87/32)*(23/32)/8 > 0.
These give K(p) > 0 for every p >= 7 (see largep.md, Step 7).

PART B (audit, NOT needed for the proof).  Rigorous Arb evaluation, at a grid
of exact dyadic points (tau, p), of
    S(tau)      exact slack  e(1+tau)c^{-p}(N^{p/2}F - M)   (tau-form, Step 3),
    S2(tau)     the same quantity computed directly from M(s,p) at s = tau/(1+tau),
    Delta(tau)  the Jensen-relaxed slack (Step 4),
    tau^2 K(p)  the final lower bound (Step 6),
and checks  S == S2 (to relative 1e-40),  S - Delta == e*tau*jg where
    jg = (1-tau) + tau*Gamma(p+1) - E(1-tau+tau R)^p is the Jensen gap,  jg > 0 (tau<1),
    jg == 0 (tau=1),
Delta > tau^2 K(p) > 0.  It also checks the closed-form evaluator against the
SPEC reference values.  Any failed or undecided ball comparison raises.

PART C (negative control).  The Part B audit must fail at p = 6 and p = 34/5
(below 4e-4, where K(p) < 0); K(34/5) < 0 < K(69/10) is certified.

Exit status is non-zero on any failure.
"""
import hashlib
import json
import math
import sys
from fractions import Fraction

import flint
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from certificate_io import point, interval, positive_int, rational, runtime
from flint import arb, fmpq, ctx

PREC = 320
ctx.prec = PREC


class CertificationError(Exception):
    pass


def require(cond, msg):
    # cond must be the Python bool True; an undecided arb comparison yields False.
    if cond is not True:
        raise CertificationError(msg)


def A_(q):
    """Exact rational -> arb ball."""
    q = Fraction(q)
    return arb(fmpq(q.numerator, q.denominator))


E = arb.const_e()

# ---------------------------------------------------------------------------
# PART A: constants used in the analytic proof (Step 7 of largep.md)
# ---------------------------------------------------------------------------
def part_a():
    out = {}
    e_hi = Fraction(87, 32)
    require(E < A_(e_hi), "e < 87/32 not certified")
    require(E - 2 < A_(Fraction(23, 32)), "e - 2 < 23/32 not certified")
    require(E > A_(Fraction(5, 2)), "e > 5/2 not certified")
    # Gamma(7) = 6! = 720 (exact integer identity)
    require(math.factorial(6) == 720, "6! != 720")
    g7 = arb(7).gamma()
    require(abs(g7 - 720) < A_(Fraction(1, 10**30)), "Arb Gamma(7) != 720")
    kappa = Fraction(11, 2) - 2 * e_hi           # lower bound for p/2 - 2e + 2 at p >= 7
    require(kappa == Fraction(1, 16), "kappa != 1/16")
    lin = kappa * 720                            # (1/16)*720 = 45
    require(lin == 45, "kappa*720 != 45")
    c_lin = e_hi / 2                             # e(p-2)/2 <= (87/64) p
    c_const = e_hi * Fraction(23, 32) / 8        # e(e-2)/(p+1) <= (87/32)(23/32)/8 for p >= 7
    K_lin = lin - c_lin                          # coefficient of p in the lower bound
    require(K_lin > 0, "45 - 87/64 <= 0")
    K7 = 7 * K_lin - c_const                     # lower bound for K(p) at p = 7; increasing in p
    require(K7 > 0, "lower bound for K(7) not positive")
    out["e_upper"] = str(e_hi)
    out["e_minus_2_upper"] = "23/32"
    out["kappa_lower(p>=7)"] = str(kappa)
    out["Gamma(p+1) >= 720 p for p >= 7"] = "Gamma increasing on [2,inf), Gamma(7)=720"
    out["K(p) >= (45 - 87/64) p - 2001/8192"] = {"K_lin": str(K_lin), "const": str(c_const)}
    out["K_lower_at_p=7"] = str(K7)
    out["K_lower_at_p=7_float"] = float(K7)
    return out


# ---------------------------------------------------------------------------
# Closed-form Arb evaluators
# ---------------------------------------------------------------------------
def Jfun(x, p):
    """J(x) = int_0^x u^p e^u du = x^{p+1}/(p+1) * 1F1(p+1; p+2; x), x >= 0."""
    if x == 0:
        return arb(0)
    return x ** (p + 1) / (p + 1) * x.hypgeom_1f1(p + 1, p + 2)


def Eplus(a, b, p):
    """E (a + b R)^p for a >= 0, b > 0:  b^p e^{a/b} Gamma(p+1, a/b)."""
    if a == 0:
        return b ** p * (p + 1).gamma()
    z = a / b
    return b ** p * z.exp() * z.gamma_upper(p + 1)


def h(x, p):
    """E|x - R|^p = e^{-x} (Gamma(p+1) + J(x)), x >= 0."""
    return (-x).exp() * ((p + 1).gamma() + Jfun(x, p))


def M_direct(s, p):
    """M(s,p) = s E(d+sR)^p + c E|d-cR|^p  (SPEC representation), s in (0,1/2]."""
    c = 1 - s
    d = 1 - 2 * s
    return s * Eplus(d, s, p) + c * c ** p * h(d / c, p)


def F(p):
    return h(arb(1), p)


def G(s, p):
    if s == 0:
        return F(p).log() / p
    N = s * s + (1 - s) ** 2
    return M_direct(s, p).log() / p - N.log() / 2


def slack_terms(t, p):
    """tau-form quantities of Steps 3-6."""
    Gm = (p + 1).gamma()
    J1 = Jfun(arb(1), p)
    A = (1 + t) * (1 + t * t) ** (p / 2)
    Jt = Jfun(1 - t, p) if t != 1 else arb(0)
    S = A * (Gm + J1) - E * t * Eplus(1 - t, t, p) - t.exp() * (Gm + Jt)
    Delta = (A - t.exp() - E * t * t) * Gm + A * J1 - t.exp() * Jt - E * t * (1 - t)
    K = (p / 2 - 2 * E + 2) * Gm - E * (p - 2) / 2 - E * (E - 2) / (p + 1)
    return S, Delta, t * t * K, K


def part_b(p_list, n_tau):
    positive_int(n_tau, "n_tau")
    require(bool(p_list), "empty audit grid")
    p_list = [rational(p) for p in p_list]
    require(all(p > 1 for p in p_list), "audit requires p > 1")
    refs = [  # SPEC reference values (mpmath, 40 digits): (s, p, G)
        (Fraction(1, 10), Fraction(13, 10), "-0.213049905296183"),
        (Fraction(3, 10), Fraction(19, 10), "-0.0294149150315886"),
        (Fraction(1, 4), Fraction(5, 2), "0.143192764694943"),
        (Fraction(2, 5), Fraction(5), "0.641602497665675"),
        (Fraction(1, 20), Fraction(12), "1.58099627165213"),
        (Fraction(1, 2), Fraction(6, 5), "-0.265784034621107"),
        (Fraction(0), Fraction(33, 10), "0.375483768956575"),
    ]
    tol = A_(Fraction(1, 10**13))
    for s, p, ref in refs:
        g = G(A_(s), A_(p))
        diff = g - arb(ref)
        require(abs(diff) < tol, f"reference mismatch at s={s}, p={p}: {g} vs {ref}")
    # G0 and Gh reference values at p = 1.3
    p13 = A_(Fraction(13, 10))
    require(abs(G(arb(0), p13) - arb("-0.213229500079513")) < tol, "G0(1.3) mismatch")
    require(abs(G(A_(Fraction(1, 2)), p13) - arb("-0.227966317234103")) < tol, "Gh(1.3) mismatch")

    rows = []
    rel = A_(Fraction(1, 10**40))
    for pq in p_list:
        p = A_(pq)
        F_p = F(p)
        min_ratio = None
        for k in range(1, n_tau + 1):
            tq = Fraction(k, n_tau)
            t = A_(tq)
            S, Delta, L, K = slack_terms(t, p)
            # identity check against the direct s-evaluator
            sq = tq / (1 + tq)
            s = A_(sq)
            c = 1 - s
            N = s * s + c * c
            S2 = E * (1 + t) * c ** (-p) * (N ** (p / 2) * F_p - M_direct(s, p))
            require(abs(S - S2) < rel * abs(S), f"identity S==S2 fails at tau={tq}, p={pq}")
            require(S > 0, f"S > 0 fails at tau={tq}, p={pq}")
            # Jensen step: S - Delta = e*tau*jg with jg = (1-tau) + tau*Gamma(p+1) - E(1-tau+tau R)^p
            jg = (1 - t) + t * (p + 1).gamma() - Eplus(1 - t, t, p)
            require(abs(S - Delta - E * t * jg) < rel * abs(S),
                    f"identity S - Delta == e*tau*jg fails at tau={tq}, p={pq}")
            if tq < 1:
                require(jg > 0, f"Jensen gap > 0 fails at tau={tq}, p={pq}")
            else:
                require(abs(jg) < rel * (p + 1).gamma(), f"Jensen gap == 0 at tau=1 fails, p={pq}")
            require(Delta > L, f"Delta > tau^2 K fails at tau={tq}, p={pq}")
            require(L > 0, f"tau^2 K > 0 fails at tau={tq}, p={pq}")
            # G-level statement: G(s,p) < G0(p)
            require(G(s, p) < G(arb(0), p), f"G < G0 fails at s={sq}, p={pq}")
            r = L / S
            rmid = float(r.mid())
            if min_ratio is None or rmid < min_ratio[0]:
                min_ratio = (rmid, str(tq))
        rows.append({"p": str(pq), "n_tau": n_tau, "K(p)_mid": K.mid().str(10, radius=False),
                     "min_over_grid(tau^2 K / S)": min_ratio[0], "argmin_tau": min_ratio[1]})
    return rows


def part_c():
    """Negative control: the audit must FAIL below the threshold 4e-4 of K(p) > 0."""
    out = {}
    for pq in [Fraction(6), Fraction(34, 5)]:
        try:
            part_b([pq], 64)
        except CertificationError as exc:
            out[str(pq)] = "failed as expected: " + str(exc)
        else:
            raise CertificationError(f"negative control unexpectedly passed at p={pq}")
    K = lambda p: (p / 2 - 2 * E + 2) * (p + 1).gamma() - E * (p - 2) / 2 - E * (E - 2) / (p + 1)
    k68, k69 = K(A_(Fraction(34, 5))), K(A_(Fraction(69, 10)))
    require(k68 < 0, "K(6.8) < 0 not certified")
    require(k69 > 0, "K(6.9) > 0 not certified")
    out["K(6.8)"] = k68.str(8, radius=False)
    out["K(6.9)"] = k69.str(8, radius=False)
    return out


def main():
    runtime()
    out_path = None
    if "--out" in sys.argv:
        out_path = sys.argv[sys.argv.index("--out") + 1]
    p_list = [Fraction(7), Fraction(15, 2), Fraction(8), Fraction(10), Fraction(12),
              Fraction(20), Fraction(50), Fraction(100), Fraction(1000)]
    n_tau = 64
    with open(__file__, "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    try:
        a = part_a()
        b = part_b(p_list, n_tau)
        c = part_c()
    except CertificationError as exc:
        print("CERTIFICATION FAILED:", exc, file=sys.stderr)
        sys.exit(2)
    cert = {
        "certificate": "MRTT cycle42 LARGE-P (Lemma LP, largep.md)",
        "claim": "For all p >= 7 and s in [0,1/2]: M(s,p) <= N(s)^{p/2} F(p), strict for s>0.",
        "proof_type": "analytic; PART A certifies the numerical constants of Step 7",
        "script": "largep_certify.py",
        "script_sha256": sha,
        "python_flint": flint.__version__,
        "arb_prec_bits": PREC,
        "part_A_constants": a,
        "part_B_audit": {
            "note": "audit only, not used by the proof",
            "reference_values_checked": 9,
            "grid": {"p": [str(x) for x in p_list], "tau": f"k/{n_tau}, k=1..{n_tau}"},
            "checks": ["S==S2 (rel 1e-40)", "S>0", "S-Delta==e*tau*jensen_gap (rel 1e-40)",
                       "jensen_gap>0 (tau<1)", "jensen_gap==0 (tau=1)",
                       "Delta>tau^2K", "tau^2K>0", "G(s,p)<G0(p)"],
            "rows": b,
        },
        "part_C_negative_control": c,
        "status": "PASS",
    }
    text = json.dumps(cert, indent=2)
    if out_path:
        with open(out_path, "x") as fh:
            fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
