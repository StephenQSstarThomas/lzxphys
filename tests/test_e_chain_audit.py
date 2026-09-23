"""Persistent coverage of the exact E-adapter requirements in the human review."""
from collections import Counter
import unittest

import sympy as sp

from scripts import ade_phase as ade
from scripts.ade_subgroups import sympy_quaternion_multiply as mul
from scripts.check_ade_phase import analytic_inputs


def canonical(v):
    values = tuple(map(ade._actual, v))
    first = next(x for x in values if x != 0)
    scale = first*ade._sign(first)
    return tuple(sp.simplify(x/scale) for x in values)


def collected(chain):
    result = Counter()
    for coefficient, vertices in chain:
        vertices = tuple(map(canonical, vertices))
        if any(a == b for a, b in zip(vertices, vertices[1:])):
            continue
        result[vertices] += coefficient
    return {k: v for k, v in result.items() if v}


def boundary(chain):
    return [(c*(-1)**j, v[:j]+v[j+1:])
            for c, v in chain for j in range(len(v))]


class ExactEChainTests(unittest.TestCase):
    def test_all_group_elements_have_unit_norm_and_inverse(self):
        for name, order in [('2T', 24), ('2O', 48), ('2I', 120)]:
            points = tuple(map(ade._reduce_point, ade.ade_group(name)))
            self.assertEqual(len(set(points)), order)
            self.assertIn(ade.I4, points)
            self.assertIn((-1, 0, 0, 0), points)
            for q in points:
                self.assertEqual(ade._sign(sum(x*x for x in q)-1), 0)
                inverse = (q[0], -q[1], -q[2], -q[3])
                self.assertIn(inverse, points)
                self.assertEqual(ade._reduce_point(mul(q, inverse)), ade.I4)

    def test_quadratic_relations_are_reduced_exactly(self):
        s2, s5 = sp.Symbol('s2'), sp.Symbol('s5')
        self.assertEqual(ade._reduce_expr(s2*s2-2), 0)
        self.assertEqual(ade._reduce_expr(s5*s5-5), 0)
        q = (s2/2, s2/2, 0, 0)
        self.assertEqual(ade._reduce_point(mul(q, q)), (0, 1, 0, 0))
        phi = (1+s5)/2
        self.assertEqual(ade._reduce_expr(phi*phi-phi-1), 0)
        self.assertEqual(ade._reduce_expr(phi*(phi-1)-1), 0)

    def test_exact_boundary_safety_and_equivariance_in_every_branch(self):
        # Dropping rank-deficient chains, using float hull tests, or moving the
        # cone anchor independently of left multiplication breaks these identities.
        half = sp.Rational(1, 2)
        r = (-half, half, half, half)
        r2 = (-half, -half, -half, -half)
        for name in ('2T', '2O', '2I'):
            a, b, x = analytic_inputs(name)
            p2 = ade._reduce_point(mul(a, b))
            p3 = ade._reduce_point(mul(p2, x))
            cases = ((ade.I4, a, p2, p3), (ade.I4, x, ade.I4, x),
                     (ade.I4, (-1, 0, 0, 0), x, ade.I4),
                     (ade.I4, (-1, 0, 0, 0), ade.I4, (-1, 0, 0, 0)),
                     (ade.I4, ade.I4, x, ade.I4), (ade.I4, r, r2, x))
            self.assertFalse(ade._safe((ade.I4, r, r2)))
            for index, vertices in enumerate(cases):
                with self.subTest(group=name, branch=index):
                    chain = ade.exact_tetrahedron_chain(vertices)
                    faces = [(c*(-1)**j, w) for j in range(4)
                             for c, w in ade._face(*(vertices[:j]+vertices[j+1:]))]
                    self.assertEqual(collected(boundary(chain)), collected(faces))
                    self.assertTrue(all(ade._safe(w) for _, w in chain))
                    self.assertLessEqual(len(chain), 24)
                    translate = ade.ade_group(name)[-1]
                    translated = ade.exact_tetrahedron_chain(tuple(
                        ade._reduce_point(mul(translate, v)) for v in vertices))
                    expected = [(c, tuple(ade._reduce_point(mul(translate, v)) for v in w))
                                for c, w in chain]
                    self.assertEqual(collected(translated), collected(expected))
                    for j in range(4):
                        face = vertices[:j]+vertices[j+1:]
                        edges = [(c*(-1)**i, w) for i in range(3)
                                 for c, w in ade._edge(*(face[:i]+face[i+1:]))]
                        self.assertEqual(collected(boundary(ade._face(*face))), collected(edges))

    def test_positive_ray_scaling_preserves_chain(self):
        for name in ('2T', '2O', '2I'):
            a, b, c = analytic_inputs(name)
            vertices = (ade.I4, a, ade._reduce_point(mul(a, b)),
                        ade._reduce_point(mul(mul(a, b), c)))
            cone = (ade.I4, (-1, 0, 0, 0), c, ade.I4)
            for case in (vertices, cone):
                scaled = tuple(tuple((i+2)*x for x in v) for i, v in enumerate(case))
                self.assertEqual(collected(ade.exact_tetrahedron_chain(case)),
                                 collected(ade.exact_tetrahedron_chain(scaled)))


if __name__ == '__main__':
    unittest.main()
