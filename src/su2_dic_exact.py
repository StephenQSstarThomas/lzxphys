"""All-n proof data for the binary dihedral (Dic_n = 2D_n) phase table.

The table of scripts/dic_formula.py is n-independent once written on the
normaliser K = N_SU(2)(U(1)) = {V^s exp(2 pi i x sigma_1)}: elements (x, s),
x in [0,1), s in {0,1}, product (x,s)(y,t) = ({(-1)^t x + y + st/2}, s xor t),
and Dic_n = {x in Z/2n}.  The phase is exp(-2 pi i k q) with q = Vol/(2 pi^2):

 000,100: -x fl(y+z)            001: (fr(z-x-y)-1/2) fl(x+y)      010: x z
 011: -(fr(y-x)-1/2) fl(x+fr(1/2-y+z))                             101: fr(1/2-x-y+z) y
 110: (x-1/2) fl(fr(1/2-x+y)+z)                                    111: -fr(1/2-x+y) fr(1/2-y+z)

Proved here (exactly, for every real x, hence every n):
 * closure: delta q is a Z-combination of products of floors (formal identity);
 * normalisation: q vanishes when any argument is the identity;
 * class: on every cyclic Z_M in U(1) the pairing with the lens-space fundamental
   cycle is 1/M, equal to the geometric representative; with Wigner's theorem and
   H^4(BK;Z) = Z c_2 this identifies the class on K, hence on every Dic_n.
 * finite cross-checks: explicit fundamental cycle of S^3/Dic_n (join triangulation,
   contracting-homotopy chain map) gives pairing exactly 1/(4n).
No floating point anywhere.  Proof text: docs/review-0923/SU2_POLAR_DUAL_EXACT.md.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product

import sympy as s

HALF = Fraction(1, 2)


def _floor(v):
    return v.numerator//v.denominator


def _frac(v):
    return v-_floor(v)


def k_product(g, h):
    (x, a), (y, b) = g, h
    return (_frac((-1)**b*x+y+HALF*a*b), a ^ b)


def k_inverse(g):
    x, a = g
    return (_frac(-x), 0) if a == 0 else (_frac(x+HALF), 1)


def table_exponent(g1, g2, g3):
    """q(g1, g2, g3) exactly (Fractions); the phase is exp(-2 pi i k q)."""
    (x, s1), (y, s2), (z, s3) = g1, g2, g3
    br = (s1, s2, s3)
    if br in ((0, 0, 0), (1, 0, 0)):
        return -x*_floor(y+z)
    if br == (0, 0, 1):
        return (_frac(z-x-y)-HALF)*_floor(x+y)
    if br == (0, 1, 0):
        return x*z
    if br == (0, 1, 1):
        return -(_frac(y-x)-HALF)*_floor(x+_frac(HALF-y+z))
    if br == (1, 0, 1):
        return _frac(HALF-x-y+z)*y
    if br == (1, 1, 0):
        return (x-HALF)*_floor(_frac(HALF-x+y)+z)
    return -_frac(HALF-x+y)*_frac(HALF-y+z)


def dic_element(n, r, epsilon):
    """Table coordinates of V^epsilon U^r in Dic_n."""
    return (Fraction(r % (2*n), 2*n), epsilon)


# ---------------------------------------------------------------------------
# Formal floor algebra on K (closure and normalisation for all real x)
# ---------------------------------------------------------------------------

class _FloorAlgebra:
    """Polynomials in x1..x4 and integer floor atoms F[L], L affine.

    Only two identities are used: fr(L) = L - F[L] and F[L + m] = F[L] + m for an
    integer combination m of atoms and integers.
    """

    def __init__(self):
        self.x = s.symbols('x1:5', real=True)
        self.atoms = {}
        self.affine = {}

    def floor(self, expr):
        expr = s.expand(expr)
        atom_part = s.S.Zero
        atoms = set(self.atoms.values())
        for term in s.Add.make_args(expr):
            if term.free_symbols & atoms:
                coefficient, monomial = term.as_coeff_Mul()
                # only integer multiples of products of atoms may leave the floor
                if not coefficient.is_integer or not monomial.free_symbols <= atoms:
                    raise ArithmeticError('non-integer term inside a floor')
                atom_part += term
        affine = s.expand(expr-atom_part)
        constant = affine.subs({v: 0 for v in self.x})
        shift = s.floor(constant)
        affine = s.expand(affine-shift)
        if not affine.free_symbols:
            return shift+atom_part
        key = s.srepr(affine)
        if key not in self.atoms:
            self.atoms[key] = s.Symbol('F[%s]' % affine, integer=True)
            self.affine[self.atoms[key]] = affine
        return self.atoms[key]+shift+atom_part

    def frac(self, expr):
        return s.expand(expr-self.floor(expr))

    def q(self, g1, g2, g3):
        (x, s1), (y, s2), (z, s3) = g1, g2, g3
        fl, fr, h = self.floor, self.frac, s.Rational(1, 2)
        br = (s1, s2, s3)
        if br in ((0, 0, 0), (1, 0, 0)):
            return -x*fl(y+z)
        if br == (0, 0, 1):
            return (fr(z-x-y)-h)*fl(x+y)
        if br == (0, 1, 0):
            return x*z
        if br == (0, 1, 1):
            return -(fr(y-x)-h)*fl(x+fr(h-y+z))
        if br == (1, 0, 1):
            return fr(h-x-y+z)*y
        if br == (1, 1, 0):
            return (x-h)*fl(fr(h-x+y)+z)
        return -fr(h-x+y)*fr(h-y+z)

    def mul(self, g, h):
        (x, a), (y, b) = g, h
        return (self.frac((-1)**b*x+y+s.Rational(a*b, 2)), a ^ b)

    def non_integer_part(self, expr):
        """Terms of expr that are not manifestly integers (integer coefficient times atoms)."""
        rest = s.S.Zero
        atoms = set(self.atoms.values())
        for term in s.Add.make_args(s.expand(expr)):
            coefficient, monomial = term.as_coeff_Mul()
            if coefficient.is_integer and monomial.free_symbols <= atoms:
                continue
            rest += term
        return rest


def formal_closure_certificate():
    """delta q for the 16 sign patterns; each residual must be exactly 0."""
    out = {}
    for signs in product((0, 1), repeat=4):
        A = _FloorAlgebra()
        g = [(A.x[i], signs[i]) for i in range(4)]
        g12, g23, g34 = A.mul(g[0], g[1]), A.mul(g[1], g[2]), A.mul(g[2], g[3])
        delta = (A.q(g[1], g[2], g[3])-A.q(g12, g[2], g[3])+A.q(g[0], g23, g[3])
                 - A.q(g[0], g[1], g34)+A.q(g[0], g[1], g[2]))
        out[signs] = A.non_integer_part(delta)
    return out


def formal_normalisation_certificate():
    """q with the identity in any slot, for all 8 branches; uses only x_i in [0,1)."""
    out = {}
    for signs in product((0, 1), repeat=2):
        for slot in range(3):
            A = _FloorAlgebra()
            e = (s.S.Zero, 0)
            others = [(A.x[0], signs[0]), (A.x[1], signs[1])]
            args = others[:slot]+[e]+others[slot:]
            value = s.expand(A.q(*args))
            # a single-variable atom F[x_i] vanishes because 0 <= x_i < 1
            zero = {atom: 0 for atom, affine in A.affine.items() if affine in A.x[:2]}
            out[(slot, signs)] = s.expand(value.subs(zero))
    return out


# ---------------------------------------------------------------------------
# Class identification data
# ---------------------------------------------------------------------------

def lens_pairing(M):
    """<q, [S^3/Z_M]> for Z_M = <u>, u = (1/M, 0), via the join triangulation.

    Vertices a_r = u^r (circle A) and b_s = u^{-s} b_0 (circle B) with the
    u-invariant order a < b; the vertex map is an equivariant chain map and the
    orbit representatives (a0, a1, b_s, b_{s+1}) give the bar cycle
    sum_j [u | u^j | u^{-1}].  The geometric representative pairs to +1/M.
    """
    u = (Fraction(1, M), 0)
    total = sum((table_exponent(u, (Fraction(j, M), 0), k_inverse(u)) for j in range(M)), Fraction(0))
    return total % 1


def _perm_sign(p):
    p, sign = list(p), 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            sign = -sign
    return sign


@lru_cache(maxsize=None)
def dic_fundamental_cycle(n):
    """Bar 3-cycle representing [S^3/Dic_n] (homogeneous tuples -> coefficient).

    Join triangulation of S^3 with vertex set Dic_n: a_r = U^r, b_s = V U^s, positive
    tetrahedra (a_r, a_{r+1}, b_s, b_{s+1}); orbit representatives s = 0..n-1.  The
    equivariant chain map is phi(rep) = h(phi(boundary rep)), h(g0..gk) = (e,g0..gk),
    extended by the (free) group action.
    """
    elements = [(Fraction(r, 2*n), t) for t in (0, 1) for r in range(2*n)]
    e = (Fraction(0), 0)
    reps, cache = {}, {}

    def act(g, simplex):
        return tuple(k_product(g, v) for v in simplex)

    def locate(simplex):
        for rep in reps.setdefault(len(simplex), []):
            for g in elements:
                image = act(g, rep)
                if set(image) == set(simplex):
                    return rep, g, _perm_sign([image.index(v) for v in simplex])
        reps[len(simplex)].append(tuple(simplex))
        return tuple(simplex), e, 1

    def phi(simplex):
        if len(simplex) == 1:
            return {simplex: 1}
        rep, g, sign = locate(simplex)
        if rep not in cache:
            out = {}
            for i in range(len(rep)):
                for t, c in phi(rep[:i]+rep[i+1:]).items():
                    out[(e,)+t] = out.get((e,)+t, 0)+(-1)**i*c
            cache[rep] = {t: c for t, c in out.items() if c}
        return {act(g, t): sign*c for t, c in cache[rep].items()}
    total = {}
    for sidx in range(n):
        simplex = ((Fraction(0), 0), (Fraction(1, 2*n), 0),
                   (Fraction(sidx, 2*n), 1), (Fraction((sidx+1) % (2*n), 2*n), 1))
        for t, c in phi(simplex).items():
            total[t] = total.get(t, 0)+c
    return tuple(sorted((t, c) for t, c in total.items() if c))


def is_cycle_in_coinvariants(chain):
    boundary = {}
    for t, c in chain:
        for i in range(4):
            face = t[:i]+t[i+1:]
            g = k_inverse(face[0])
            face = tuple(k_product(g, v) for v in face)
            boundary[face] = boundary.get(face, 0)+(-1)**i*c
    return all(v == 0 for v in boundary.values())


def inhomogeneous(t):
    h0, h1, h2, h3 = t
    return (k_product(k_inverse(h0), h1), k_product(k_inverse(h1), h2), k_product(k_inverse(h2), h3))


def dic_pairing(n):
    """<q, phi[S^3/Dic_n]> modulo 1 (expected 1/(4n))."""
    chain = dic_fundamental_cycle(n)
    if not is_cycle_in_coinvariants(chain):
        raise ArithmeticError('Constructed chain is not a cycle')
    return sum((c*table_exponent(*inhomogeneous(t)) for t, c in chain), Fraction(0)) % 1


# ---------------------------------------------------------------------------
# Nondegenerate Dic_n tetrahedra: exact double-arc volumes (Li2 identity proved)
# ---------------------------------------------------------------------------

def dic_point(g):
    """Unit quaternion (as exact sympy tuple) of the table element (x, s)."""
    x, t = g
    angle = 2*s.pi*s.Rational(x.numerator, x.denominator)
    c, v = s.cos(angle), s.sin(angle)
    return (c, v, 0, 0) if t == 0 else (0, 0, c, v)


def dic_nondegenerate_volume(n, g1, g2, g3):
    """Oriented volume of the geometric tetrahedron (1, g1, g1g2, g1g2g3), det != 0.

    Four linearly independent points of Dic_n lie two on each great circle
    A = span(1, i), B = span(j, k); the tetrahedron is the join of an A-arc and a
    B-arc, V = alpha beta / 2 (the double-arc identity), alpha, beta in (pi/n) Z.
    """
    p1 = g1
    p2 = k_product(p1, g2)
    p3 = k_product(p2, g3)
    cumulative = ((Fraction(0), 0), p1, p2, p3)
    A = [g for g in cumulative if g[1] == 0]
    B = [g for g in cumulative if g[1] == 1]
    if len(A) != 2 or len(B) != 2 or A[0] == A[1] or B[0] == B[1]:
        raise ValueError('Degenerate Dic_n input: use the global chain')
    arcs = []
    for P, Q in (A, B):
        d = _frac(Q[0]-P[0])
        arcs.append(min(d, 1-d))
    if HALF in arcs:
        raise ValueError('Antipodal pair: degenerate tetrahedron')
    matrix = s.Matrix.hstack(*(s.Matrix(dic_point(g)) for g in cumulative))
    orientation = int(s.sign(s.simplify(matrix.det())))
    volume = s.Rational(2)*s.pi**2*s.Rational(arcs[0].numerator, arcs[0].denominator) \
        * s.Rational(arcs[1].numerator, arcs[1].denominator)
    return {'n': n, 'cumulative': cumulative, 'arcs_over_2pi': tuple(arcs),
            'orientation': orientation, 'volume': volume,
            'oriented_volume': orientation*volume}
