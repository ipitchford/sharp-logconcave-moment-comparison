"""Exact certificate endpoints. Decimal displays are not proof bounds."""
from fractions import Fraction
import flint


def require(condition, message):
    if condition is not True:
        raise ValueError(message)


def positive_int(value, name):
    require(type(value) is int and value > 0, name + " must be a positive integer")


def rational(value):
    require(not isinstance(value, (bool, float)), "exact rational input required")
    return Fraction(str(value))


def point(value):
    """Serialize a finite, exact Arb dyadic without decimal rounding."""
    require(value.is_finite() and value.is_exact(), "finite exact Arb endpoint required")
    mantissa, exponent = value.man_exp()
    m, e = int(mantissa), int(exponent)
    q = Fraction(m * (2 ** max(e, 0)), 2 ** max(-e, 0))
    return str(q)


def interval(value):
    require(value.is_finite(), "nonfinite enclosure")
    lo, hi = point(value.lower()), point(value.upper())
    require(Fraction(lo) <= Fraction(hi), "reversed enclosure")
    # Exact dyadic round-trip, at the current precision, must contain the ball.
    from flint import arb, fmpq
    require(arb(fmpq(lo)) <= value.lower() and arb(fmpq(hi)) >= value.upper(),
            "serialized endpoints do not contain the enclosure")
    return [lo, hi]


def runtime():
    require(flint.__version__ == "0.9.0", "python-flint 0.9.0 required")
    return {"python_flint": flint.__version__, "FLINT": flint.__FLINT_VERSION__,
            "Arb": "integrated in FLINT; separate Arb version not exposed",
            "bound_encoding": "exact rational strings (dyadic Arb endpoints)"}
