"""Exact root-of-unity representative for E6/2T, E7/2O and E8/2I.

This is COHOMOLOGOUS to, not pointwise equal to, ade_phase's geometric phase.
See docs/review-0922/SU2_E_ALGEBRAIC_RESEARCH.md, equations (1)--(9).
All branch and orbit decisions use exact arithmetic. No small float is used
for the formal probe (1, epsilon, epsilon**2, epsilon**3).
"""
from fractions import Fraction
from functools import lru_cache
from numbers import Integral

import mpmath as mp
import sympy as sp

from scripts.ade_phase import (
    I4, ade_group, exact_tetrahedron_chain, phase_details, _reduce_point, _safe,
)
from scripts.ade_subgroups import sympy_quaternion_multiply


def _pair(value, d):
    """Embed exact input in Q + Q sqrt(d); reject floats and other fields."""
    value = sp.sympify(value)
    if value.has(sp.Float):
        raise TypeError('Exact rational or quadratic coordinates required, not floats')
    symbols = value.free_symbols
    if symbols and (not d or any(symbol.name != 's'+str(d) for symbol in symbols)):
        raise ValueError('Coordinate is not in the selected group field')
    value = sp.expand(sp.radsimp(value.subs({s: sp.sqrt(d) for s in symbols})))
    b = value.coeff(sp.sqrt(d)) if d else sp.Integer(0)
    a = sp.expand(value-b*sp.sqrt(d)) if d else value
    if not (a.is_Rational and b.is_Rational):
        raise ValueError('Coordinate is not in the selected group field')
    return Fraction(int(a.p), int(a.q)), Fraction(int(b.p), int(b.q))


def _add(x, y):
    return x[0]+y[0], x[1]+y[1]


def _neg(x):
    return -x[0], -x[1]


def _mul(x, y, d):
    return x[0]*y[0]+d*x[1]*y[1], x[0]*y[1]+x[1]*y[0]


def _sign(x, d):
    a, b = x
    sa, sb = (a > 0)-(a < 0), (b > 0)-(b < 0)
    if not b:
        return sa
    if not a or sa == sb:
        return sb
    difference = a*a-d*b*b
    return sa*((difference > 0)-(difference < 0))


def _dot(x, y, d):
    result = (0, 0)
    for a, b in zip(x, y):
        result = _add(result, _mul(a, b, d))
    return result


def _det3(rows, d):
    result = (0, 0)
    for j in range(3):
        k, l = (j+1) % 3, (j+2) % 3
        minor = _add(_mul(rows[1][k], rows[2][l], d),
                     _neg(_mul(rows[1][l], rows[2][k], d)))
        result = _add(result, _mul(rows[0][j], minor, d))
    return result


def _cramer_rows(vertices, d):
    """Rows of adj(V), where vertices are V's columns; no field division."""
    rows = []
    for j in range(4):
        row = []
        for i in range(4):
            minor = [[vertices[c][r] for c in range(4) if c != j]
                     for r in range(4) if r != i]
            value = _det3(minor, d)
            row.append(_neg(value) if (i+j) % 2 else value)
        rows.append(tuple(row))
    return tuple(rows)


def _left_columns(h):
    """Columns h*e_s for the MINUS-cross quaternion convention."""
    a, x, y, z = h
    return ((a, x, y, z), (_neg(x), a, _neg(z), y),
            (_neg(y), z, a, _neg(x)), (_neg(z), _neg(y), x, a))


@lru_cache(maxsize=3)
def _group_data(name):
    points = tuple(map(_reduce_point, ade_group(name)))
    d = {'2T': 0, '2O': 2, '2I': 5}[name]
    pairs = tuple(tuple(_pair(x, d) for x in q) for q in points)
    lookup = {q: i for i, q in enumerate(pairs)}
    return points, d, lookup, tuple(_left_columns(h) for h in pairs)


def _input_indices(name, inputs):
    points, d, lookup, _ = _group_data(name)
    indices = []
    for q in inputs:
        if len(q) != 4:
            raise ValueError('A group element must have exactly four coordinates')
        key = tuple(_pair(x, d) for x in q)
        if key not in lookup:
            raise ValueError('Input is not an element of '+name)
        indices.append(lookup[key])
    return tuple(indices)


def _level(k):
    if not isinstance(k, (Integral, sp.Integer)):
        raise TypeError('The level k must be an integer')
    return int(k)


@lru_cache(maxsize=16384)
def _count(name, indices):
    """Cached immutable exact result, independent of level and precision."""
    points, d, _, left_matrices = _group_data(name)
    g1, g2, g3 = (points[i] for i in indices)
    p2 = _reduce_point(sympy_quaternion_multiply(g1, g2))
    p3 = _reduce_point(sympy_quaternion_multiply(p2, g3))
    cumulative = (I4, g1, p2, p3)
    chain = exact_tetrahedron_chain(cumulative)
    terms, total = [], 0
    for coefficient, vertices in chain:
        pairs = tuple(tuple(_pair(x, d) for x in v) for v in vertices)
        rows = _cramer_rows(pairs, d)
        orientation = _sign(_dot(rows[0], pairs[0], d), d)
        hits = []
        if orientation:
            for h, columns in enumerate(left_matrices):
                inside = True
                for row in rows:
                    first_sign = next((_sign(value, d) for column in columns
                                       if (value := _dot(row, column, d)) != (0, 0)), 0)
                    if first_sign != orientation:
                        inside = False
                        break
                if inside:
                    hits.append(h)
        total += coefficient*orientation*len(hits)
        terms.append((coefficient, vertices, orientation, tuple(hits)))
    branch = 'empty' if not chain else ('direct' if _safe(cumulative) else 'cone')
    return total, cumulative, branch, tuple(terms)


def count_details(name, g1, g2, g3, k=1):
    """JSON-ready signed count and exact root exponent, not a decimal phase.

    Inputs must be exact MEMBERS of the named group, not arbitrary unit rays.
    The result means exp(2*pi*i*exponent_mod_order/order). `count` is the
    unreduced signed count needed in the real coboundary comparison.
    """
    k = _level(k)
    indices = _input_indices(name, (g1, g2, g3))
    points = _group_data(name)[0]
    total, cumulative, branch, terms = _count(name, indices)
    return {
        'group': name, 'order': len(points), 'k': k,
        'representative': 'orbit-count; cohomologous, not equal, to source-A geometry',
        'probe': 'lexicographic positive (1, epsilon, epsilon^2, epsilon^3)',
        'indices': list(indices),
        'inputs': [[str(x) for x in points[i]] for i in indices],
        'cumulative_vertices': [[str(x) for x in v] for v in cumulative],
        'branch': branch, 'count': total,
        'exponent_mod_order': (-k*total) % len(points),
        'terms': [{'coefficient': c, 'orientation': orientation,
                   'vertices': [[str(x) for x in v] for v in vertices],
                   'hit_indices': list(hits),
                   'hit_elements': [[str(x) for x in points[h]] for h in hits],
                   'signed_count': c*orientation*len(hits)}
                  for c, vertices, orientation, hits in terms],
    }


def algebraic_phase(name, g1, g2, g3, k=1, dps=80):
    """Numerical value of the exact root; large levels reduced BEFORE conversion."""
    k = _level(k)
    indices = _input_indices(name, (g1, g2, g3))
    total = _count(name, indices)[0]
    order = len(_group_data(name)[0])
    with mp.workdps(int(dps)+15):
        return mp.exp(2j*mp.pi*((-k*total) % order)/order)


@lru_cache(maxsize=256)
def _gauge_average(name, indices, digits):
    """Unreduced real average B(g,h)=mean_x(A(x,g,h)/N-t(x,g,h))."""
    points = _group_data(name)[0]
    g, h = (points[i] for i in indices)
    order = len(points)
    terms, differences = [], []
    with mp.workdps(digits+15):
        for index, x in enumerate(points):
            count = _count(name, (index,)+indices)[0]
            geometry = phase_details(x, g, h, dps=digits)
            t = mp.mpf(geometry['volume_over_pi2'])/2
            difference = mp.mpf(count)/order-t
            differences.append(difference)
            terms.append((index, count, mp.nstr(t, digits),
                          mp.nstr(difference, digits)))
        average = mp.fsum(differences)/order
        return mp.nstr(average, digits), tuple(terms)


def gauge_details(name, g, h, k=1, dps=50):
    """Explicit beta with ALL N summands; w_count = w_geometry * delta(beta).

    Unlike the integer counting formula, this comparison uses numerical Li_2
    volumes. Reported errors are precision checks, not interval certificates.
    Counts MUST remain unreduced when computing this real cochain.
    """
    k = _level(k)
    indices = _input_indices(name, (g, h))
    points = _group_data(name)[0]
    level_digits = (abs(k).bit_length()*30103)//100000+1
    digits = int(dps)+level_digits+20
    average, terms = _gauge_average(name, indices, digits)
    with mp.workdps(digits+15):
        beta = mp.exp(-2j*mp.pi*k*mp.mpf(average))
        return {'group': name, 'order': len(points), 'k': k, 'dps': int(dps),
                'convention': 'beta=exp(-2*pi*i*k*B); count_phase=geometric_phase*delta(beta)',
                'inputs': [[str(x) for x in points[i]] for i in indices],
                'B': average,
                'terms': [{'x_index': index, 'x': [str(x) for x in points[index]],
                           'count': count, 'geometric_volume_over_2pi2': t,
                           'count_over_N_minus_t': difference}
                          for index, count, t, difference in terms],
                'beta': {'real': mp.nstr(mp.re(beta), int(dps)),
                         'imag': mp.nstr(mp.im(beta), int(dps))}}
