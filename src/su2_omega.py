"""SU(2) anomaly phases with exact geometric predicates and finite Li_2 volumes.

Quaternions use g=q0*I+i*q.sigma, hence the vector product has MINUS cross.
Each rational ray represents its unit normalization, never an antipodal quotient.
Numerical inputs are interpreted as their exact decimal/binary rational values,
then normalized. No determinant tolerance is used to select geometric branches.
"""
from fractions import Fraction
from itertools import combinations
from numbers import Integral
from math import isfinite

import mpmath as mp


def _fraction(value):
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError('Finite quaternion coordinates required')
        return Fraction(value)
    if isinstance(value, mp.mpf):
        if not mp.isfinite(value):
            raise ValueError('Finite quaternion coordinates required')
        sign, mantissa, exponent, _ = value._mpf_
        return Fraction((-1)**sign * int(mantissa)) * Fraction(2)**int(exponent)
    return Fraction(str(value))


def quaternion(values):
    """Exact positive-ray canonicalization; the represented element has norm one."""
    q = tuple(_fraction(x) for x in values)
    if len(q) != 4 or not any(q):
        raise ValueError('A quaternion needs four real coordinates and nonzero norm')
    scale = abs(next(x for x in q if x))
    return tuple(x / scale for x in q)


IDENTITY = quaternion((1, 0, 0, 0))


def quaternion_multiply(left, right):
    a, x, y, z = quaternion(left)
    b, u, v, w = quaternion(right)
    return quaternion((a*b-x*u-y*v-z*w, a*u+b*x-y*w+z*v,
                       a*v+b*y-z*u+x*w, a*w+b*z-x*v+y*u))


def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def _rref(rows):
    a = [list(map(Fraction, row)) for row in rows]
    pivots = []
    for col in range(len(a[0])):
        pivot = next((r for r in range(len(pivots), len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        row = len(pivots)
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [x/scale for x in a[row]]
        for r in range(len(a)):
            if r != row and a[r][col]:
                scale = a[r][col]
                a[r] = [x-scale*y for x, y in zip(a[r], a[row])]
        pivots.append(col)
        if len(pivots) == len(a):
            break
    return a, pivots


def _det(rows):
    a = [list(map(Fraction, row)) for row in rows]
    result = Fraction(1)
    for col in range(len(a)):
        pivot = next((r for r in range(col, len(a)) if a[r][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            result = -result
        scale = a[col][col]
        result *= scale
        for r in range(col + 1, len(a)):
            factor = a[r][col] / scale
            for j in range(col + 1, len(a)):
                a[r][j] -= factor * a[col][j]
    return result


def _inverse(rows):
    size = len(rows)
    reduced, pivots = _rref([list(row) + [Fraction(i == j) for j in range(size)]
                             for i, row in enumerate(rows)])
    if pivots[:size] != list(range(size)):
        raise ValueError('Singular matrix')
    return [row[size:] for row in reduced]


def _mp(value):
    value = _fraction(value)
    return mp.mpf(value.numerator) / value.denominator


def _rank(vertices):
    return len(_rref(list(zip(*vertices)))[1]) if vertices else 0


def hemisphere_safe(vertices):
    """Exact test: zero is outside the convex hull of the unit vertices.

    Positive dependencies are unchanged by independent positive ray scalings.
    A minimal positive dependence has a one-dimensional kernel on its support.
    """
    vertices = tuple(quaternion(v) for v in vertices)
    if _rank(vertices) == len(vertices):
        return True
    for count in range(2, len(vertices) + 1):
        for subset in combinations(vertices, count):
            reduced, pivots = _rref(list(zip(*subset)))
            if len(pivots) != count - 1:
                continue
            free = next(j for j in range(count) if j not in pivots)
            kernel = [Fraction(j == free) for j in range(count)]
            for row, col in enumerate(pivots):
                kernel[col] = -reduced[row][free]
            if all(value >= 0 for value in kernel) or all(value <= 0 for value in kernel):
                return False
    return True


def _radial(vertices):
    vertices = tuple(quaternion(v) for v in vertices)
    if any(a == b for a, b in zip(vertices, vertices[1:])):
        return []
    return [(1, vertices)]


def edge_chain(a, b):
    a, b = quaternion(a), quaternion(b)
    if hemisphere_safe((a, b)):
        return _radial((a, b))
    midpoint = quaternion_multiply(a, (0, 1, 0, 0))
    return _radial((a, midpoint)) + _radial((midpoint, b))


def _cone(anchor, chain):
    if not chain:
        return []
    ranks = [_rank(vertices) for _, vertices in chain]
    for t in range(3 * len(chain) + 1):
        apex = quaternion_multiply(anchor, (1, t, t*t, t*t*t))
        if all(_rank(vertices + (apex,)) > rank
               for (_, vertices), rank in zip(chain, ranks)):
            return [(sign, (apex,) + vertices) for sign, vertices in chain]
    raise ArithmeticError('Finite cone construction violated its rank bound')


def face_chain(a, b, c):
    a, b, c = map(quaternion, (a, b, c))
    if a == b or b == c:
        return []
    if hemisphere_safe((a, b, c)):
        return _radial((a, b, c))
    boundary = edge_chain(b, c) + [(-sign, v) for sign, v in edge_chain(a, c)] + edge_chain(a, b)
    return _cone(a, boundary)


def tetrahedron_chain(vertices):
    """Global left-equivariant filling, at most 24 radial tetrahedra.

    Vertices denote UNIT-normalized points for chain parameterizations. Raw rays
    are retained only for exact predicates and reparameterization-invariant volumes.
    Reversing vertex order is not simplified to a negative singular simplex.
    """
    a, b, c, d = map(quaternion, vertices)
    if a == b or b == c or c == d:
        return []
    if hemisphere_safe((a, b, c, d)):
        return _radial((a, b, c, d))
    boundary = (face_chain(b, c, d)
                + [(-sign, v) for sign, v in face_chain(a, c, d)]
                + face_chain(a, b, d)
                + [(-sign, v) for sign, v in face_chain(a, b, c)])
    return _cone(a, boundary)


# Murakami's face pairs: 01,02,12,23,13,03.
# Actual edges, respectively: 23,13,03,01,02,12.
FACE_PAIRS = ((0, 1), (0, 2), (1, 2), (2, 3), (1, 3), (0, 3))


def oriented_volume(vertices, dps=50):
    """Oriented convex spherical tetrahedron volume (hemisphere-safe input).

    Eight principal dilogarithms, Murakami arXiv:1011.2584v4 Thm. 1.1.
    Exact rank/sign checks; work precision increases with the angle Gram condition.
    """
    vertices = tuple(quaternion(q) for q in vertices)
    if len(vertices) != 4:
        raise ValueError('Four vertices required')
    determinant = _det(list(zip(*vertices)))
    if not determinant:
        if hemisphere_safe(vertices):
            return mp.mpf(0)
        raise ValueError('Radial simplex hits zero: use the global chain evaluator')
    gram = [[_dot(a, b) for b in vertices] for a in vertices]
    inverse = _inverse(gram)
    det_angle_gram = 1 / (_det(gram) * mp_product(inverse[i][i] for i in range(4)))
    # For unit vertices |Qt|<=1, so |V|>=|D_unit|/6. This lower bound
    # budgets cancellation in a tiny volume, separately from angular conditioning.
    det_unit_squared = determinant**2 / mp_product(gram[i][i] for i in range(4))
    with mp.workdps(30):
        extra = max(0, int(mp.ceil(-mp.log10(_mp(det_angle_gram)))),
                    int(mp.ceil(-mp.log10(mp.sqrt(_mp(det_unit_squared))/6))))
    with mp.workdps(int(dps) + extra + 25):
        angles, phases = [], []
        for i, j in FACE_PAIRS:
            diag = inverse[i][i] * inverse[j][j]
            cosine = -_mp(inverse[i][j]) / mp.sqrt(_mp(diag))
            sine = mp.sqrt(_mp((diag - inverse[i][j]**2) / diag))
            angles.append(mp.atan2(sine, cosine))
            phases.append(mp.mpc(cosine, sine))
        a, b, c, d, e, f = phases
        r0 = a*d+b*e+c*f+a*b*f+a*c*e+b*c*d+d*e*f+a*b*c*d*e*f
        r1 = 4 * sum(mp.sin(angles[j]) * mp.sin(angles[j+3]) for j in range(3))
        r2 = mp.conj(r0)
        discriminant_root = 4 * mp.sqrt(_mp(det_angle_gram))
        # Rationalized form of (-r1 + sqrt(discriminant))/(2*r2).
        z = -2*r0 / (r1 + discriminant_root)
        if not abs(z) < 1 or not abs(r2):
            raise ArithmeticError('Murakami root check failed; increase precision')
        positive = (z, z/(a*b*d*e), z/(a*c*d*f), z/(b*c*e*f))
        negative = (-z/(a*b*c), -z/(a*e*f), -z/(b*d*f), -z/(c*d*e))
        real_l = (sum(mp.re(mp.polylog(2, u)) for u in positive)
                  - sum(mp.re(mp.polylog(2, u)) for u in negative)
                  - sum(angles[j]*angles[j+3] for j in range(3))) / 2
        raw = -real_l + mp.pi*(mp.arg(-r2) + sum(angles)/2) - 3*mp.pi**2/2
        volume = raw % (2*mp.pi**2)
        if not 0 < volume < mp.pi**2:
            raise ArithmeticError('Volume branch outside (0,pi^2); increase precision')
        return +volume if determinant > 0 else -volume


def mp_product(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def matrix_quaternion(matrix, tolerance=1e-10):
    """Validate a numerical 2x2 SU(2) matrix, then extract its normalized ray.

    A tolerated numerical residual is projected to the SU(2) quaternion form.
    Use quaternion inputs for exact algebraic products on branch loci.
    """
    if hasattr(matrix, 'shape'):
        valid_shape = matrix.shape == (2, 2)
    elif hasattr(matrix, 'rows') and hasattr(matrix, 'cols'):
        valid_shape = (matrix.rows, matrix.cols) == (2, 2)
    else:
        try:
            valid_shape = len(matrix) == 2 and all(len(row) == 2 for row in matrix)
        except TypeError:
            valid_shape = False
    if not valid_shape:
        raise ValueError('A 2x2 SU(2) matrix is required')
    try:
        a, b = matrix[0][0], matrix[0][1]
        c, d = matrix[1][0], matrix[1][1]
    except (IndexError, TypeError):
        try:
            a, b, c, d = matrix[0, 0], matrix[0, 1], matrix[1, 0], matrix[1, 1]
        except (IndexError, TypeError) as exc:
            raise ValueError('A 2x2 SU(2) matrix is required') from exc
    # Preserve the stored components BEFORE any mpmath arithmetic can round them.
    (ar, ai), (br, bi), (cr, ci), (dr, di) = [
        (_fraction(z.real), _fraction(z.imag)) for z in (a, b, c, d)]
    tol = _fraction(tolerance)
    if tol < 0:
        raise ValueError('A nonnegative matrix tolerance is required')
    if ((dr-ar)**2+(di+ai)**2 > tol**2
            or (cr+br)**2+(ci-bi)**2 > tol**2
            or abs(ar**2+ai**2+br**2+bi**2-1) > tol):
        raise ValueError('Input is not an SU(2) matrix within the stated tolerance')
    return quaternion(((ar+dr)/2, (bi+ci)/2, (br-cr)/2, (ai-di)/2))


def omega_quaternions(g1, g2, g3, k=1, dps=50):
    if not isinstance(k, Integral):
        raise TypeError('The level k must be an integer')
    first = quaternion(g1)
    second = quaternion_multiply(first, g2)
    third = quaternion_multiply(second, g3)
    chain = tetrahedron_chain((IDENTITY, first, second, third))
    # Phase error grows with |k|. bit_length bounds decimal digits without the
    # interpreter's limit on converting huge integers to decimal strings.
    level_digits = (abs(int(k)).bit_length() * 30103) // 100000 + 1
    phase_dps = int(dps) + level_digits
    with mp.workdps(phase_dps + 15):
        volume = mp.fsum(sign * oriented_volume(vertices, dps=phase_dps+5) for sign, vertices in chain)
        return mp.exp(1j * int(k) * volume / mp.pi)


def omega_su2(g1, g2, g3, k=1, dps=50):
    """Evaluate the phase from three numerical 2x2 SU(2) matrices."""
    return omega_quaternions(*(matrix_quaternion(g) for g in (g1, g2, g3)), k=k, dps=dps)


def omega_quaternions_source(g1, g2, g3, k=1, dps=50):
    """The 0922 source-note convention, Eq. (7.1).

    The original repository helper uses +i*k*V/pi for its selected F-move.
    The source note reports Z_before/Z_after with -i*k*V/pi, so this is the
    inverse representative. The level and geometric chain are unchanged.
    """
    return 1 / omega_quaternions(g1, g2, g3, k=k, dps=dps)


def omega_su2_source(g1, g2, g3, k=1, dps=50):
    """Source-note phase for numerical 2x2 SU(2) matrices (Eq. (7.1))."""
    return 1 / omega_su2(g1, g2, g3, k=k, dps=dps)
