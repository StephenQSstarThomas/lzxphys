"""Li2-free exact phase of the ORIGINAL E/F/T representative for every E-group input.

The chain T(1, g1, g1g2, g1g2g3) is exactly the one of su2_exact_ade (same edge
split J=e1, same face and tetrahedron cone points a(1,t,t^2,t^3)); nothing about
the representative is changed.  Only its volume modulo 2 pi^2 is rewritten:

1. Nondegenerate safe input: one group tetrahedron, su2_polar_dual.
2. Cone input T = C_a'(Z), Z = sum eps_F F closed.  For an unsafe face
   F = C_y(L) (y = the ORIGINAL face cone point) put F' = C_m(L) with a group
   point m and D_F = sum_e c_e [y, m, u_e, w_e].  Then F = F' - dD_F exactly,
   and for any admissible group point g
       vol T == vol C_g(Z') - sum eps_F vol D_F        (mod 2 pi^2),
   Z' = sum eps_F F'.  C_g(Z') consists of group tetrahedra (polar formula).
3. vol D_F is elementary:
   * L on one great circle (plane P): D_F is W copies of the lune between the
     half-spheres through y and m, vol = sigma W pi psi, psi = angle(y_perp, m_perp);
   * L = two half great circles from p to -p through m1, m2 (m = m1): D_F is the
     suspension of the triangle (y, m1, m2) between p and -p,
     vol = t (pi/2) Omega, Omega = area of its shadow on the equator p^perp.
Proof: docs/review-0923/SU2_POLAR_DUAL_EXACT.md, section on singular inputs.
"""
from itertools import combinations

import sympy as s
from su2_exact_ade import (_clean_point, _cone, _edge, _face, _rank, _safe,
                           exact_group)
from su2_polar_dual import edge_fraction, _surd_from_sympy, FIELD, tetrahedron_volume
from su2_symbolic import IDENTITY, exact_level, quaternion_product, unit_quaternion


def _vec(p):
    return s.Matrix(p)


def _det(points):
    return s.simplify(s.Matrix.hstack(*(_vec(p) for p in points)).det())


def _sgn(x):
    value = s.sign(s.simplify(x))
    if value not in (-1, 0, 1):
        raise ArithmeticError('Undecided exact sign')
    return int(value)


def _dot(u, v):
    return sum(a*b for a, b in zip(u, v))


def _loop(x, y, z):
    """The E-loop of an unsafe face, in exactly the order used by _face."""
    return _edge(y, z)+[(-w, v) for w, v in _edge(x, z)]+_edge(x, y)


def _admissible(apex, simplices):
    return all(_safe((apex,)+v) for v in simplices)


def _circle_correction(name, y, m, loop):
    """vol D for a loop on one great circle: sigma * W * pi * psi."""
    points = [p for _, e in loop for p in e]
    base = [points[0]]
    for p in points[1:]:
        if _rank(tuple(base)+(p,)) == 2:
            base.append(p)
            break
    U = s.Matrix.hstack(*(_vec(b) for b in base))
    proj = U*(U.T*U).inv()*U.T

    def perp(v):
        return [s.simplify(x) for x in (_vec(v)-proj*_vec(v))]
    yp, mp = perp(y), perp(m)
    c2 = s.radsimp(_dot(yp, mp)**2/(_dot(yp, yp)*_dot(mp, mp)))
    if c2 == 1 and _sgn(_dot(yp, mp)) < 0:
        return None, None               # psi = pi: segment [y, m] meets P; choose another m
    psi = s.acos(_sgn(_dot(yp, mp))*s.sqrt(c2))

    def plane_sign(u, w):             # orientation of (u, w) in P w.r.t. basis (base0, base1)
        cu = (U.T*U).inv()*U.T*_vec(u)
        cw = (U.T*U).inv()*U.T*_vec(w)
        return _sgn(cu[0]*cw[1]-cu[1]*cw[0])
    winding = s.S.Zero
    sigma = None
    for c, (u, w) in loop:
        q = edge_fraction(_surd_from_sympy(_dot(u, w), FIELD[name]))
        winding += c*plane_sign(u, w)*s.Rational(q.numerator, q.denominator)/2
        tetra = (y, m, u, w)
        if _rank(tetra) == 4:
            local = c*_sgn(_det(tetra))*c*plane_sign(u, w)   # = sigma_perp
            if sigma not in (None, local):
                raise ArithmeticError('Inconsistent lune orientation')
            sigma = local
    if not winding.is_integer:
        raise ArithmeticError('Loop does not close on the circle')
    if winding == 0 or sigma is None:
        return s.S.Zero, {'type': 'circle', 'winding': 0}
    return sigma*winding*s.pi*psi, {'type': 'circle', 'winding': int(winding), 'sigma': sigma, 'psi': psi}


def _digon_correction(p, y, m1, m2, loop):
    """vol D = t (pi/2) Omega for the suspension of triangle (y, m1, m2) over +-p."""
    t = None
    for c, (u, w) in loop:
        tetra = (y, m1, u, w)
        if _rank(tetra) == 4:
            local = c*_sgn(_det(tetra))
            if t not in (None, local):
                raise ArithmeticError('Inconsistent suspension orientation')
            t = local
    if t is None:
        return s.S.Zero, {'type': 'digon', 'degenerate': True}
    pp = _dot(p, p)

    def eq(v):
        k = _dot(v, p)/pp
        return [s.simplify(a-k*b) for a, b in zip(v, p)]
    Y, M1, M2 = eq(y), eq(m1), eq(m2)
    nY, n1, n2 = (s.sqrt(s.radsimp(_dot(v, v))) for v in (Y, M1, M2))
    det = s.Abs(_det((p, Y, M1, M2)))/s.sqrt(pp)
    den = nY*n1*n2+_dot(Y, M1)*n2+_dot(M1, M2)*nY+_dot(M2, Y)*n1
    omega = 2*s.atan2(det, den)
    return t*s.pi/2*omega, {'type': 'digon', 't': t, 'Omega': omega}


def _face_replacement(name, group, face):
    """Unsafe face -> (F' chain of group triangles, elementary vol D_F, record)."""
    x, y0, z = face
    loop = _loop(x, y0, z)
    cone = _cone(x, loop)
    if cone != _face(x, y0, z):
        raise AssertionError('Face cone does not reproduce the original chain')
    apex = cone[0][1][0]
    loop_points = []
    for _, e in loop:
        for p in e:
            if p not in loop_points:
                loop_points.append(p)
    rank = _rank(tuple(loop_points))
    if rank == 2:
        for m in group:
            if _rank(tuple(loop_points[:2])+(m,)) != 3 or m in loop_points:
                continue
            new = [(c, (m,)+e) for c, e in loop]
            D = [(c, (apex, m)+e) for c, e in loop]
            if all(_safe(v) for _, v in new) and all(_safe(v) for _, v in D):
                vol, record = _circle_correction(name, apex, m, loop)
                if vol is None:
                    continue
                return new, vol, dict(record, face_apex=apex, group_apex=m)
        raise ArithmeticError('No admissible group apex for a circular face loop')
    if rank == 3:
        antipodal = [(a, b) for a, b in combinations(loop_points, 2)
                     if all(s.simplify(u+v) == 0 for u, v in zip(a, b))]
        if len(antipodal) != 1:
            raise ArithmeticError('Unexpected digon structure')
        p = antipodal[0][0]
        mids = [v for v in loop_points if v not in antipodal[0]]
        if len(mids) != 2:
            raise ArithmeticError('Unexpected digon structure')
        for m1, m2 in (mids, mids[::-1]):
            new = [(c, (m1,)+e) for c, e in loop]
            D = [(c, (apex, m1)+e) for c, e in loop]
            if all(_safe(v) for _, v in new) and all(_safe(v) for _, v in D):
                vol, record = _digon_correction(p, apex, m1, m2, loop)
                return new, vol, dict(record, face_apex=apex, group_apex=m1)
        raise ArithmeticError('No admissible digon apex')
    raise ArithmeticError('Unexpected face loop rank')


def global_polar_formula(name, g1, g2, g3, k=1):
    """Exact elementary phase of the original global E/F/T representative."""
    k = exact_level(k)
    group = exact_group(name)
    members = set(group)
    gs = tuple(unit_quaternion(g) for g in (g1, g2, g3))
    if any(g not in members for g in gs):
        raise ValueError('An input is not a member of '+name)
    a, b, c = gs
    p2 = _clean_point(quaternion_product(a, b))
    vertices = (IDENTITY, _clean_point(a), p2, _clean_point(quaternion_product(p2, c)))
    result = {'group': name, 'inputs': gs, 'cumulative_vertices': vertices,
              'representative': 'original E/F/T geometric phase (unchanged)'}
    if any(u == v for u, v in zip(vertices, vertices[1:])):
        volume, branch, terms, corrections = s.S.Zero, 'empty', (), ()
    elif _safe(vertices):
        branch, corrections = 'direct', ()
        if _rank(vertices) == 4:
            r = tetrahedron_volume(name, vertices)
            volume, terms = r['oriented_volume'], ((1, vertices, r['oriented_volume']),)
        else:
            volume, terms = s.S.Zero, ()
    else:
        branch = 'cone'
        p0, p1, q2, q3 = vertices
        faces = ((1, (p1, q2, q3)), (-1, (p0, q2, q3)), (1, (p0, p1, q3)), (-1, (p0, p1, q2)))
        z_prime, corrections = [], []
        correction_total = s.S.Zero
        for eps, face in faces:
            if face[0] == face[1] or face[1] == face[2]:
                continue
            if _safe(face):
                z_prime.append((eps, face))
                continue
            new, vol, record = _face_replacement(name, group, face)
            z_prime.extend((eps*w, v) for w, v in new)
            correction_total += eps*vol
            corrections.append(dict(record, sign=eps, volume=vol, face=face))
        simplices = [v for _, v in z_prime]
        apex = next((g for g in group if _admissible(g, simplices)), None)
        if apex is None:
            raise ArithmeticError('No admissible group apex for the closed face cycle')
        terms, total = [], s.S.Zero
        for w, v in z_prime:
            tetra = (apex,)+v
            if _rank(tetra) == 4:
                r = tetrahedron_volume(name, tetra)
                terms.append((w, tetra, r['oriented_volume']))
                total += w*r['oriented_volume']
        volume = total-correction_total
        result['group_apex'] = apex
        terms, corrections = tuple(terms), tuple(corrections)
    result.update(branch=branch, terms=terms, corrections=corrections,
                  oriented_volume_mod=volume,
                  phase=s.exp(-s.I*k*volume/s.pi))
    return result
