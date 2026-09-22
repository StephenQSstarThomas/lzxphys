"""Finite ADE phase evaluator with exact algebraic branch selection.

This adapter closes the remaining gap between the exact SymPy representatives
of 2T/2O/2I and the mpmath Li_2 volume routine. Convex-hull, rank, antipodal,
and cone choices are made in exact algebraic arithmetic; only each nonzero
tetrahedron volume is converted to high-precision real numbers for Li_2.
"""
from itertools import combinations
import sympy as sp
import mpmath as mp

from scripts.ade_subgroups import (
    binary_tetrahedral, binary_octahedral, binary_icosahedral,
    sympy_quaternion_multiply,
)
from su2_omega import oriented_volume


I4 = (sp.Integer(1), sp.Integer(0), sp.Integer(0), sp.Integer(0))
J1 = (sp.Integer(0), sp.Integer(1), sp.Integer(0), sp.Integer(0))


def _root_data(expr):
    expr = sp.sympify(expr)
    for symbol in expr.free_symbols:
        if symbol.name == "s2":
            return symbol, 2
        if symbol.name == "s5":
            return symbol, 5
    return None, None


def _reduce_expr(expr):
    expr = sp.expand(sp.sympify(expr))
    root, square = _root_data(expr)
    if root is None:
        return expr
    poly = sp.Poly(root**2-square, root, domain=sp.QQ)
    return sp.expand(sp.rem(sp.Poly(expr, root, domain=sp.QQ), poly).as_expr())


def _reduce_point(point):
    return tuple(_reduce_expr(value) for value in point)


def _actual(expr):
    expr = _reduce_expr(expr)
    root, square = _root_data(expr)
    if root is None:
        return sp.simplify(expr)
    return sp.simplify(expr.subs(root, sp.sqrt(square)))


def _eq(a, b):
    return all(sp.simplify(_actual(x-y)) == 0 for x, y in zip(a, b))


def _rank(vertices):
    vertices = tuple(_reduce_point(v) for v in vertices)
    rows, columns = 4, len(vertices)
    for size in range(min(rows, columns), 0, -1):
        for row_indices in combinations(range(rows), size):
            for col_indices in combinations(range(columns), size):
                matrix = sp.Matrix([[ _actual(vertices[j][i]) for j in col_indices]
                                     for i in row_indices])
                if matrix.det() != 0:
                    return size
    return 0


def _sign(expr):
    value = sp.simplify(_actual(expr))
    sign = sp.sign(value)
    if sign in (-1, 0, 1):
        return int(sign)
    raise ArithmeticError(f"could not determine algebraic sign: {expr}")


def _safe(vertices):
    """Exact zero-outside-convex-hull test for at most four algebraic points."""
    vertices = tuple(_reduce_point(v) for v in vertices)
    if _rank(vertices) == len(vertices):
        return True
    for size in range(2, len(vertices)+1):
        for subset in combinations(vertices, size):
            matrix = sp.Matrix.hstack(*(sp.Matrix([_actual(x) for x in v]) for v in subset))
            if _rank(subset) != size-1:
                continue
            null = matrix.nullspace()
            if len(null) != 1:
                continue
            coeff = [_reduce_expr(x) for x in null[0]]
            signs = [_sign(x) for x in coeff if _sign(x) != 0]
            if signs and (all(x > 0 for x in signs) or all(x < 0 for x in signs)):
                return False
    return True


def _moment(t):
    # Only the positive ray matters. Avoid introducing sqrt(85) into K_d.
    return (sp.Integer(1), sp.Integer(t), sp.Integer(t*t), sp.Integer(t*t*t))


def _radial(vertices):
    if any(_eq(a, b) for a, b in zip(vertices, vertices[1:])):
        return []
    return [(1, tuple(vertices))]


def _edge(a, b):
    if _eq(a, b):
        return []
    if _safe((a, b)):
        return _radial((a, b))
    midpoint = _reduce_point(sympy_quaternion_multiply(a, J1))
    return _radial((a, midpoint)) + _radial((midpoint, b))


def _cone(anchor, chain):
    ranks = [_rank(vertices) for _, vertices in chain]
    for t in range(3*len(chain)+1):
        apex = _reduce_point(sympy_quaternion_multiply(anchor, _moment(t)))
        if all(_rank(vertices + (apex,)) > rank
               for (_, vertices), rank in zip(chain, ranks)):
            return [(sign, (apex,) + vertices) for sign, vertices in chain]
    raise ArithmeticError("exact ADE cone candidate bound was exceeded")


def _face(a, b, c):
    if _eq(a, b) or _eq(b, c):
        return []
    if _safe((a, b, c)):
        return _radial((a, b, c))
    boundary = (_edge(b, c) + [(-s, v) for s, v in _edge(a, c)]
                + _edge(a, b))
    return _cone(a, boundary)


def exact_tetrahedron_chain(vertices):
    vertices = tuple(_reduce_point(v) for v in vertices)
    a, b, c, d = vertices
    if _eq(a, b) or _eq(b, c) or _eq(c, d):
        return []
    if _safe(vertices):
        return _radial(vertices)
    boundary = (_face(b, c, d) + [(-s, v) for s, v in _face(a, c, d)]
                + _face(a, b, d) + [(-s, v) for s, v in _face(a, b, c)])
    return _cone(a, boundary)


def _numeric(point, dps):
    with mp.workdps(dps):
        converted = []
        for value in point:
            value = sp.sympify(value)
            substitutions = {}
            for symbol in value.free_symbols:
                if symbol.name == "s2":
                    substitutions[symbol] = sp.sqrt(2)
                elif symbol.name == "s5":
                    substitutions[symbol] = sp.sqrt(5)
            numeric = sp.N(value.subs(substitutions), dps)
            converted.append(mp.mpf(str(numeric)))
        return tuple(converted)


def phase_from_exact_quaternions(g1, g2, g3, k=1, dps=80):
    """Evaluate the source-note A convention for exact algebraic quaternions."""
    p0 = I4
    p1 = _reduce_point(g1)
    p2 = _reduce_point(sympy_quaternion_multiply(g1, g2))
    p3 = _reduce_point(sympy_quaternion_multiply(p2, g3))
    chain = exact_tetrahedron_chain((p0, p1, p2, p3))
    with mp.workdps(dps):
        volume = mp.mpf(0)
        for sign, vertices in chain:
            determinant = sp.Matrix.hstack(*(sp.Matrix([_actual(x) for x in v]) for v in vertices)).det()
            if sp.simplify(determinant) == 0:
                continue
            volume += sign * oriented_volume(tuple(_numeric(v, dps) for v in vertices), dps=dps-5)
        return mp.exp(-1j * int(k) * volume / mp.pi)


def ade_group(name):
    if name == "2T":
        return binary_tetrahedral()
    if name == "2O":
        return binary_octahedral()
    if name == "2I":
        return binary_icosahedral()
    raise ValueError("name must be 2T, 2O, or 2I")


def representative_phase(name, indices=(0, 1, 2), k=1, dps=80):
    points = ade_group(name)
    return phase_from_exact_quaternions(*(points[i] for i in indices), k=k, dps=dps)
