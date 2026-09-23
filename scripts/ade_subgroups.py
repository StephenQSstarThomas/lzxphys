"""Exact label arithmetic and explicit SU(2) representatives for binary ADE groups.

The Dic_n labels are exact presentation elements; their quaternion coordinates are
only used after the group relation has already been checked at the label level.
The three binary polyhedral constructors use SymPy algebraic expressions so that
closure checks do not depend on floating point tolerances.
"""
from dataclasses import dataclass
from itertools import permutations, product
from math import cos, pi, sin, sqrt


@dataclass(frozen=True)
class DicElement:
    """a^r b^epsilon in Dic_n, where a^(2n)=1, b^2=a^n, bab^-1=a^-1."""

    n: int
    r: int
    epsilon: int

    def __post_init__(self):
        if self.n < 2 or self.epsilon not in (0, 1):
            raise ValueError("Dic_n requires n>=2 and epsilon in {0,1}")
        object.__setattr__(self, "r", self.r % (2 * self.n))

    def __mul__(self, other):
        if not isinstance(other, DicElement) or other.n != self.n:
            return NotImplemented
        r = self.r + (-1 if self.epsilon else 1) * other.r
        if self.epsilon and other.epsilon:
            r += self.n
        return DicElement(self.n, r, (self.epsilon + other.epsilon) % 2)

    def inverse(self):
        for candidate in dic_elements(self.n):
            if self * candidate == identity_dic(self.n) and candidate * self == identity_dic(self.n):
                return candidate
        raise ArithmeticError("presentation inverse search failed")


def dic_elements(n):
    return [DicElement(n, r, epsilon) for epsilon in (0, 1) for r in range(2 * n)]


def identity_dic(n):
    return DicElement(n, 0, 0)


def dic_generators(n):
    return DicElement(n, 1, 0), DicElement(n, 0, 1)


def dic_quaternion(n, element):
    """Return (q0,q1,q2,q3) for a^r b^epsilon, with b=i*sigma_2."""
    if element.n != n:
        raise ValueError("element belongs to a different Dic_n")
    angle = pi * element.r / n
    a = (cos(angle), sin(angle), 0.0, 0.0)
    if not element.epsilon:
        return a
    # (q0,q1,q2,q3) multiplication by b=(0,0,1,0) with the +i sigma convention.
    return (0.0, 0.0, cos(angle), -sin(angle))


def dic_quaternion_exact_q8(element):
    """Exact rational quaternion for Dic_2=Q8."""
    from fractions import Fraction
    if element.n != 2:
        raise ValueError("the exact helper is only for Dic_2")
    powers = (
        (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(-1), Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(-1), Fraction(0), Fraction(0)),
    )
    a = powers[element.r % 4]
    if not element.epsilon:
        return a
    b = (Fraction(0), Fraction(0), Fraction(1), Fraction(0))
    x0, x1, x2, x3 = a
    y0, y1, y2, y3 = b
    return (x0*y0-x1*y1-x2*y2-x3*y3,
            x0*y1+y0*x1-x2*y3+x3*y2,
            x0*y2+y0*x2-x3*y1+x1*y3,
            x0*y3+y0*x3-x1*y2+x2*y1)


def _unique(points):
    result = []
    seen = set()
    for point in points:
        key = tuple(str(x) for x in point)
        if key not in seen:
            seen.add(key)
            result.append(tuple(point))
    return tuple(result)


def binary_tetrahedral():
    """The 24 Hurwitz units (2T), as exact SymPy tuples."""
    import sympy as sp
    half = sp.Rational(1, 2)
    points = [(s, 0, 0, 0) for s in (1, -1)]
    points += [(0, s, 0, 0) for s in (1, -1)]
    points += [(0, 0, s, 0) for s in (1, -1)]
    points += [(0, 0, 0, s) for s in (1, -1)]
    points += [tuple(s * half for s in signs) for signs in product((1, -1), repeat=4)]
    return _unique(points)


def binary_octahedral():
    import sympy as sp
    root = sp.Symbol("s2")
    points = list(binary_tetrahedral())
    h = root / 2
    for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        for signs in product((1, -1), repeat=2):
            point = [sp.Integer(0)] * 4
            point[i], point[j] = signs[0] * h, signs[1] * h
            points.append(tuple(point))
    return _unique(points)


def binary_icosahedral():
    import sympy as sp
    root = sp.Symbol("s5")
    phi = (1 + root) / 2
    invphi = phi - 1
    half = sp.Rational(1, 2)
    points = [(s, 0, 0, 0) for s in (1, -1)]
    points += [(0, s, 0, 0) for s in (1, -1)]
    points += [(0, 0, s, 0) for s in (1, -1)]
    points += [(0, 0, 0, s) for s in (1, -1)]
    points += [tuple(s * half for s in signs) for signs in product((1, -1), repeat=4)]
    base = (sp.Integer(0), half, phi * half, invphi * half)
    # Even coordinate permutations, with all independent signs on nonzero entries.
    for perm in permutations(range(4)):
        inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4))
        if inversions % 2:
            continue
        permuted = tuple(base[perm[i]] for i in range(4))
        for signs in product((1, -1), repeat=3):
            it = iter(signs)
            points.append(tuple(0 if value == 0 else value * next(it) for value in permuted))
    return _unique(points)


def sympy_quaternion_multiply(left, right):
    """Exact +i sigma quaternion multiplication for SymPy tuples."""
    a, x, y, z = left
    b, u, v, w = right
    return (
        a * b - x * u - y * v - z * w,
        a * u + b * x - y * w + z * v,
        a * v + b * y - z * u + x * w,
        a * w + b * z - x * v + y * u,
    )


def sympy_key(point, modulus=None):
    import sympy as sp
    if modulus is None:
        return tuple(sp.expand(value) for value in point)
    root, square = modulus
    polynomial = sp.Poly(root**2 - square, root, domain=sp.QQ)
    keys = []
    for value in point:
        remainder = sp.rem(sp.Poly(sp.expand(value), root, domain=sp.QQ), polynomial).as_expr()
        keys.append(sp.expand(remainder))
    return tuple(keys)


def exact_closure(points, modulus=None):
    """Check multiplication closure and return the number of distinct elements."""
    lookup = {sympy_key(point, modulus) for point in points}
    for left in points:
        for right in points:
            if sympy_key(sympy_quaternion_multiply(left, right), modulus) not in lookup:
                raise AssertionError((left, right, sympy_quaternion_multiply(left, right)))
    return len(lookup)


def normalized_numeric(point, dps=80):
    import mpmath as mp
    import sympy as sp
    with mp.workdps(dps):
        values = []
        for value in point:
            value = sp.sympify(value)
            roots = {symbol: sp.sqrt(2 if symbol.name == 's2' else 5)
                     for symbol in value.free_symbols if symbol.name in ('s2', 's5')}
            values.append(mp.mpf(str(value.subs(roots).evalf(dps))))
        norm = mp.sqrt(mp.fsum(value * value for value in values))
        return tuple(value / norm for value in values)
