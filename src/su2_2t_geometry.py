"""Actual Hurwitz tetrahedron volumes from the F4 spherical chamber tiling.

This is a geometric subdivision, not group-cohomological orbit counting.
Every face is an F4 mirror, so each accepted probe represents one ENTIRE chamber
of volume pi**2/576. See the accompanying proof for completeness of the tiling.
"""
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, gcd

import sympy as s
from su2_symbolic import IDENTITY, exact_level, quaternion_product, unit_quaternion

PAIRS = tuple(combinations(range(4),2))
PAIR_INDEX = {pair:i for i,pair in enumerate(PAIRS)}
PERMUTED_EDGES = tuple(tuple(PAIR_INDEX[tuple(sorted((p[i],p[j])))] for i,j in PAIRS)
                       for p in permutations(range(4)))


def hurwitz_integer_vertices():
    """Twice the unit vertices, identity first; all entries are integers."""
    axes = tuple(tuple(sign*2 if i==j else 0 for i in range(4))
                 for j in range(4) for sign in (1,-1))
    return axes+tuple(product((1,-1), repeat=4))


def det3(rows):
    a,b,c = rows
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def cramer_rows(vertices):
    """Adjugate rows of the matrix whose columns are the vertices."""
    matrix = tuple(zip(*vertices))
    return tuple(tuple((-1)**(i+j)*det3(tuple(tuple(matrix[r][c] for c in range(4) if c!=i)
                                               for r in range(4) if r!=j))
                       for j in range(4)) for i in range(4))


def determinant(vertices):
    row = cramer_rows(vertices)[0]
    return sum(x*y for x,y in zip(row,vertices[0]))


def primitive_normal(normal):
    divisor = 0
    for x in normal:
        divisor = gcd(divisor,abs(x))
    if not divisor:
        raise ValueError('Zero normal')
    sign = 1 if next(x for x in normal if x) > 0 else -1
    return tuple(sign*x//divisor for x in normal)


@lru_cache(maxsize=1)
def hyperplane_normals():
    roots = {primitive_normal(v) for v in hurwitz_integer_vertices()}
    for i,j in combinations(range(4),2):
        for sign in (1,-1):
            roots.add(tuple(1 if k==i else sign if k==j else 0 for k in range(4)))
    return tuple(sorted(roots))


@lru_cache(maxsize=1)
def chamber_probes():
    return tuple(sorted(tuple(sign*x for sign,x in zip(signs,p))
                        for seed in ((8,3,2,1),(7,4,3,2),(6,5,4,1))
                        for p in permutations(seed) for signs in product((1,-1),repeat=4)))


def canonical_gram_key(entries):
    return min(tuple(entries[i] for i in indices) for indices in PERMUTED_EDGES)


def _geometry(vertices):
    rows = cramer_rows(vertices)
    det = sum(x*y for x,y in zip(rows[0],vertices[0]))
    if det == 0:
        raise ValueError('Singular input is not covered by the twelve-type table')
    roots = set(hyperplane_normals())
    if any(primitive_normal(row) not in roots for row in rows):
        raise ValueError('A face is not an F4 mirror; chamber counting is invalid')
    orientation = 1 if det>0 else -1
    count = sum(all(orientation*sum(x*y for x,y in zip(row,p))>0 for row in rows)
                for p in chamber_probes())
    return count,orientation


def tetrahedron_2t(vertices, k=1):
    k = exact_level(k)
    if len(vertices)!=4:
        raise ValueError('Four vertices required')
    vertices = tuple(unit_quaternion(v) for v in vertices)
    doubled = tuple(tuple(2*x for x in v) for v in vertices)
    members = set(hurwitz_integer_vertices())
    if any(v not in members for v in doubled):
        raise ValueError('All vertices must be exact Hurwitz units')
    doubled = tuple(tuple(int(x) for x in v) for v in doubled)
    count,orientation = _geometry(doubled)
    volume = Fraction(count,576)
    return {'chambers': count, 'orientation': orientation, 'volume_over_pi2': volume,
            'oriented_volume_over_pi2': orientation*volume,
            'phase': s.exp(-s.I*k*s.pi*s.Rational(orientation*count,576)),
            'method': 'actual F4 chamber subdivision; same geometric representative'}


def phase_2t(g1,g2,g3,k=1):
    g1,g2,g3 = map(unit_quaternion,(g1,g2,g3))
    p2 = quaternion_product(g1,g2)
    return tetrahedron_2t((IDENTITY,g1,p2,quaternion_product(p2,g3)),k=k)


@lru_cache(maxsize=1)
def _classification():
    vertices = hurwitz_integer_vertices()
    dot = tuple(tuple(sum(x*y for x,y in zip(a,b))//2 for b in vertices) for a in vertices)
    classes = {}
    singular = 0
    for indices in combinations(range(24),4):
        points = tuple(vertices[i] for i in indices)
        if not determinant(points):
            singular += 1
            continue
        key = canonical_gram_key(tuple(dot[indices[i]][indices[j]] for i,j in PAIRS))
        if key not in classes:
            count,orientation = _geometry(points)
            classes[key] = {'gram_key':key,'multiplicity':0,'representative_indices':indices,
                            'representative':tuple(tuple(Fraction(x,2) for x in v) for v in points),
                            'chambers':count,'volume_over_pi2':Fraction(count,576),
                            'representative_orientation':orientation}
        classes[key]['multiplicity'] += 1
    return {'total_subsets':comb(24,4),'degenerate_subsets':singular,
            'nondegenerate_subsets':sum(row['multiplicity'] for row in classes.values()),
            'types':tuple(classes[key] for key in sorted(classes))}


def classify_2t():
    """A caller cannot mutate the cached proof data."""
    return deepcopy(_classification())
