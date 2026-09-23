"""Exact Murakami expressions in the user's edge order (01,02,03,23,13,12).

No floating evaluation, fitting, or replacement of the geometric representative.
The scalar ``raw_volume`` is a volume modulo 2*pi**2, not a globally selected
positive volume. Its integer-level phase is unambiguous. Singular inputs belong
to the separate global-chain interface, not to this single-tetrahedron formula.
"""
from collections import defaultdict
from numbers import Integral

import sympy as s

EDGE_PAIRS = ((0,1), (0,2), (0,3), (2,3), (1,3), (1,2))
PLUS_SUPPORTS = ((0,1,3,4), (0,2,3,5), (1,2,4,5))
MINUS_SUPPORTS = ((3,4,5), (1,2,3), (0,2,4), (0,1,5))
IDENTITY = (s.S.One, s.S.Zero, s.S.Zero, s.S.Zero)
PAULI = (s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-s.I],[s.I,0]]),
         s.Matrix([[1,0],[0,-1]]))


def exact_level(k):
    if isinstance(k, bool) or not isinstance(k, (Integral, s.Integer)):
        raise TypeError('The level must be an exact integer')
    return s.Integer(k)


def _exact(value):
    value = s.sympify(value)
    if value.has(s.Float) or value.has(s.nan, s.oo, -s.oo, s.zoo):
        raise TypeError('Only finite exact expressions are accepted')
    return value


def unit_quaternion(value):
    """Validate, rather than approximate/project, an exact SU(2) element."""
    if isinstance(value, s.MatrixBase):
        if value.shape != (2,2):
            raise ValueError('An SU(2) matrix must be 2 by 2')
        for x in value:
            _exact(x)
        q = (s.trace(value)/2,) + tuple(s.trace(p*value)/(2*s.I) for p in PAULI)
    else:
        if len(value) != 4:
            raise ValueError('Four quaternion coordinates required')
        q = tuple(value)
    q = tuple(s.simplify(_exact(x)) for x in q)
    if any(x.is_real is not True for x in q):
        raise ValueError('Quaternion coordinates must be provably real')
    if s.simplify(sum(x*x for x in q)-1) != 0:
        raise ValueError('An exact unit quaternion is required')
    return q


def quaternion_product(left, right):
    """The negative-cross convention for q0 I + i q.sigma."""
    a,x,y,z = left
    b,u,v,w = right
    return tuple(s.expand(t) for t in
                 (a*b-x*u-y*v-z*w, a*u+b*x-y*w+z*v,
                  a*v+b*y-z*u+x*w, a*w+b*z-x*v+y*u))


def quaternion_matrix(q):
    q = unit_quaternion(q)
    return q[0]*s.eye(2) + s.I*sum((q[j+1]*PAULI[j] for j in range(3)), s.zeros(2))


def _inside_disk(z):
    return s.simplify(1-z*s.conjugate(z)).is_positive is True


def reduce_real_dilog(terms):
    """Conservative exact identities, with a record of every reduction.

    Input is (rational coefficient, argument). Unknown domains are left alone.
    This is not an assertion that these three identities generate every possible
    dilogarithm identity. Remaining terms are returned explicitly.
    """
    coefficients = defaultdict(lambda: s.S.Zero)
    steps = []

    def canonical(z):
        z = s.cancel(s.expand(z))
        if _inside_disk(z):
            bar = s.cancel(s.expand(s.conjugate(z)))
            if s.default_sort_key(bar) < s.default_sort_key(z):
                steps.append(('principal-conjugation', str(z)))
                return bar
        return z

    for coefficient, z in terms:
        coefficient = _exact(coefficient)
        if coefficient.is_Rational is not True:
            raise TypeError('Dilogarithm coefficients must be rational')
        coefficients[canonical(_exact(z))] += coefficient
    correction = s.S.Zero
    for z in tuple(coefficients):
        c = coefficients[z]
        if not c or not _inside_disk(z):
            continue
        partner = canonical(1-z)
        if partner != z and _inside_disk(partner) and coefficients.get(partner) == c:
            correction += c*(s.pi**2/6-s.re(s.log(z)*s.log(1-z)))
            coefficients[z] = coefficients[partner] = s.S.Zero
            steps.append(('Euler-reflection', str(z)))
    for z in tuple(coefficients):
        c = coefficients[z]
        partner = canonical(-z)
        if c and partner != z and _inside_disk(z) and coefficients.get(partner) == c:
            coefficients[z] = coefficients[partner] = s.S.Zero
            coefficients[canonical(z*z)] += c/2
            steps.append(('duplication', str(z)))
    remaining = tuple((c,z) for z,c in sorted(coefficients.items(), key=lambda kv:s.default_sort_key(kv[0])) if c)
    expression = correction + sum(c*s.re(s.polylog(2,z)) for c,z in remaining)
    return {'expression': expression, 'remaining': remaining, 'identities': tuple(steps)}


def tetrahedron_formula(vertices, k=1, *, reduce=True):
    """Complete single-tetrahedron Li2/log expression; det != 0 is required.

    Symbolic determinants whose zero status is undecided retain an explicit
    nondegeneracy condition. Matrix input is accepted for each unit vertex too.
    """
    k = exact_level(k)
    if len(vertices) != 4:
        raise ValueError('Four vertices required')
    vertices = tuple(unit_quaternion(v) for v in vertices)
    V = s.Matrix.hstack(*(s.Matrix(v) for v in vertices))
    det = s.factor(V.det())
    if det.is_zero is True:
        raise ValueError('Degenerate tetrahedron: use a global exact chain, not zero volume')
    gram = V.T*V
    cosines = tuple(s.simplify(gram[i,j]) for i,j in EDGE_PAIRS)
    lengths = tuple(s.acos(c) for c in cosines)
    b = tuple(c+s.I*s.sqrt(1-c*c) for c in cosines)
    a = tuple(-s.conjugate(b[(j+3)%6]) for j in range(6))
    q0 = s.expand(a[0]*a[3]+a[1]*a[4]+a[2]*a[5]
                  +a[0]*a[1]*a[5]+a[0]*a[2]*a[4]+a[1]*a[2]*a[3]
                  +a[3]*a[4]*a[5]+s.prod(a))
    q1 = s.simplify(4*sum(s.sqrt(1-cosines[j]**2)*s.sqrt(1-cosines[j+3]**2) for j in range(3)))
    q2 = s.conjugate(q0)
    discriminant = 16*det**2
    z = s.cancel(s.expand_complex(-2*q0/(q1+s.sqrt(discriminant))))
    positive = tuple(s.cancel(s.expand(z*s.prod(b[j] for j in support))) for support in PLUS_SUPPORTS)
    negative = tuple(s.cancel(s.expand(z*s.prod(b[j] for j in support))) for support in MINUS_SUPPORTS)
    Z = (z,)+positive+negative
    log_reduce = s.simplify if reduce else lambda x:x
    logarithm = lambda x:s.log(x,evaluate=reduce)
    Sigma = tuple(log_reduce(
        sum(logarithm(1-t) for support,t in zip(PLUS_SUPPORTS,positive) if j in support)
        -sum(logarithm(1-t) for support,t in zip(MINUS_SUPPORTS,negative) if j in support)) for j in range(6))
    terms = tuple((s.S.One if j<4 else -s.S.One,t) for j,t in enumerate(Z))
    unreduced = sum(c*s.polylog(2,t,evaluate=reduce) for c,t in terms)
    reduction = reduce_real_dilog(terms) if reduce else {
        'expression': s.re(unreduced,evaluate=False), 'remaining': terms, 'identities': ()}
    quadratic = -sum((s.pi-lengths[j])*(s.pi-lengths[j+3]) for j in range(3))/2
    logarithmic = -s.pi*s.arg(-q2,evaluate=reduce)-sum(lengths[j]*(s.pi-lengths[(j+3)%6]+s.im(Sigma[j],evaluate=reduce)) for j in range(6))/2-s.pi**2/2
    raw = reduction['expression']/2 + quadratic + logarithmic
    if not reduction['remaining']:
        raw = s.simplify(raw)
    orientation = s.sign(det)
    return {'vertices': vertices, 'edge_pairs': EDGE_PAIRS, 'cosines': cosines,
            'lengths': lengths, 'gram': s.ImmutableMatrix(gram), 'determinant': det,
            'orientation': orientation, 'condition': s.Ne(det,0), 'a': a,
            'q0': q0, 'q1': q1, 'q2': q2, 'discriminant': discriminant,
            'z': z, 'original_root': (-q1+s.sqrt(discriminant))/(2*q2),
            'Z': Z, 'Sigma': Sigma, 'li2_sum': unreduced,
            'dilog_reduction': reduction, 'raw_volume': raw,
            'volume_modulus': 2*s.pi**2,
            'physical_volume': s.Mod(raw, 2*s.pi**2,evaluate=False),
            'phase': s.exp(-s.I*k*orientation*raw/s.pi,evaluate=reduce)}


def group_formula(g1, g2, g3, k=1, *, reduce=True):
    """Three group increments -> six traces -> exact Li2 expression."""
    g1,g2,g3 = map(unit_quaternion, (g1,g2,g3))
    p2 = quaternion_product(g1,g2)
    p3 = quaternion_product(p2,g3)
    result = tetrahedron_formula((IDENTITY,g1,p2,p3), k=k, reduce=reduce)
    A,B,C = map(quaternion_matrix, (g1,p2,p3))
    result['orientation_trace'] = s.simplify(s.trace(A*(B*C-C*B))/4)
    result['trace_cosines'] = tuple(s.simplify(q[0]) for q in
        (g1,p2,p3,g3,quaternion_product(g2,g3),g2))
    return result


def orthogonal_join_certificate(alpha, beta):
    """Exact join identity, with its explicit domain (not an unchecked fit).

    The Li2 identity is proved by differentiating its four terms, using the
    principal logarithm factorizations given in the accompanying manuscript.
    Negative alpha reverses orientation, not the six edge lengths.
    """
    alpha,beta = _exact(alpha),_exact(beta)
    if alpha.is_real is not True or beta.is_real is not True:
        raise ValueError('Join angles must be provably real')
    condition = s.And(s.Abs(alpha)>0, s.Abs(alpha)<s.pi, beta>0, beta<s.pi)
    if condition is s.false:
        raise ValueError('Requires 0 < |alpha| < pi and 0 < beta < pi')
    rho = s.Symbol('rho', real=True)
    integral = s.integrate(s.cos(rho)*s.sin(rho), (rho,0,s.pi/2))
    a = s.Abs(alpha)
    p,q = s.exp(s.I*a),s.exp(s.I*beta)
    S = p*q+p+q-1
    eta = s.arg(S)-(a+beta)/2
    return {'condition': condition, 'rho_integral': integral,
            'oriented_volume': alpha*beta*integral,
            'lengths': (a,s.pi/2,s.pi/2,beta,s.pi/2,s.pi/2),
            'z': -2/S, 'li2_four_term_real': s.pi*eta-s.pi**2/4,
            'proof': 'principal-log derivative identity + exact orthant constant'}


def formula_template():
    """Executable straight-line Li2 formula shared by every catalogue row.

    Each pair binds its symbol to an expression involving only input symbols
    and earlier bindings. This compact exact representation avoids expanding
    hundreds of identical transcendental expression trees in the catalogue.
    Inputs: c1..c6 (cosines), detG>0, orientation in {-1,1}, integer k.
    """
    c=s.symbols('c1:7',real=True)
    l=s.symbols('l1:7',real=True)
    a=s.symbols('a1:7')
    q0,q1,q2,z=s.symbols('q0 q1 q2 z')
    Z=s.symbols('Z0:8')
    Sigma=s.symbols('Sigma1:7')
    detG=s.Symbol('detG',positive=True)
    orientation=s.Symbol('orientation',integer=True)
    k=s.Symbol('k',integer=True)
    bindings=[(l[j],s.acos(c[j])) for j in range(6)]
    bindings += [(a[j],-c[(j+3)%6]+s.I*s.sqrt(1-c[(j+3)%6]**2)) for j in range(6)]
    bindings += [(q0,a[0]*a[3]+a[1]*a[4]+a[2]*a[5]+a[0]*a[1]*a[5]
                  +a[0]*a[2]*a[4]+a[1]*a[2]*a[3]+a[3]*a[4]*a[5]+s.prod(a)),
                 (q1,4*sum(s.sqrt(1-c[j]**2)*s.sqrt(1-c[j+3]**2) for j in range(3))),
                 (q2,s.conjugate(q0)),(z,-2*q0/(q1+4*s.sqrt(detG))),(Z[0],z)]
    for j,support in enumerate(PLUS_SUPPORTS+MINUS_SUPPORTS,1):
        bindings.append((Z[j],z*s.prod(c[r]+s.I*s.sqrt(1-c[r]**2) for r in support)))
    for j in range(6):
        bindings.append((Sigma[j],sum(s.log(1-Z[r+1]) for r,S in enumerate(PLUS_SUPPORTS) if j in S)
                         -sum(s.log(1-Z[r+4]) for r,S in enumerate(MINUS_SUPPORTS) if j in S)))
    W=s.Symbol('W',real=True)
    bindings.append((W,s.re(sum(s.polylog(2,t) for t in Z[:4])-sum(s.polylog(2,t) for t in Z[4:]))/2
                     -sum((s.pi-l[j])*(s.pi-l[j+3]) for j in range(3))/2-s.pi*s.arg(-q2)
                     -sum(l[j]*(s.pi-l[(j+3)%6]+s.im(Sigma[j])) for j in range(6))/2-s.pi**2/2))
    bindings.append((s.Symbol('phase'),s.exp(-s.I*k*orientation*W/s.pi)))
    return tuple(bindings)
