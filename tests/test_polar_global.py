import importlib
import importlib.util
import unittest
import sympy as s


class PolarGlobalTests(unittest.TestCase):
    """Exact tests of the elementary global phase (original E/F/T representative)."""

    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('su2_polar_global'),
                             'Elementary global phase module is missing')
        self.api = importlib.import_module('su2_polar_global')

    def assert_elementary(self, expr):
        self.assertFalse(expr.has(s.polylog))
        self.assertFalse(expr.has(s.Float))

    def test_normalisation(self):
        e, x = (1, 0, 0, 0), (0, 1, 0, 0)
        for name in ('2T', '2O', '2I'):
            for args in ((e, x, x), (x, e, x), (x, x, e)):
                r = self.api.global_polar_formula(name, *args)
                self.assertEqual(r['branch'], 'empty')
                self.assertEqual(r['phase'], 1)

    def test_center_triple(self):
        m = (-1, 0, 0, 0)
        for name in ('2T', '2O'):
            for k in (1, 2, 3):
                r = self.api.global_polar_formula(name, m, m, m, k=k)
                self.assertEqual(r['branch'], 'cone')
                self.assert_elementary(r['oriented_volume_mod'])
                self.assertEqual(s.simplify(r['phase']-(-1)**k), 0)

    def test_direct_branch_is_the_polar_formula(self):
        polar = importlib.import_module('su2_polar_dual')
        g = ((0, 1, 0, 0), (0, 0, 0, 1), (s.Rational(1, 2),)*4)
        r = self.api.global_polar_formula('2T', *g)
        self.assertEqual(r['branch'], 'direct')
        self.assertEqual(r['oriented_volume_mod'], polar.polar_formula('2T', *g)['oriented_volume'])

    def test_degenerate_corrections_are_elementary_and_keep_cone_points(self):
        exact = importlib.import_module('su2_exact_ade')
        seen = set()
        h = (s.Rational(1, 2),)*4
        for g in (((-1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)),
                  ((0, 1, 0, 0), (0, 1, 0, 0), (0, 1, 0, 0)),
                  (h, h, h), ((0, 0, 1, 0), (-1, 0, 0, 0), h)):
            r = self.api.global_polar_formula('2T', *g)
            self.assert_elementary(r['oriented_volume_mod'])
            chain = exact.exact_chain(r['cumulative_vertices'])
            original_apexes = {v[0] for _, v in chain}
            for c in r['corrections']:
                seen.add(c['type'])
                # the correction is built on the ORIGINAL face cone point
                self.assertTrue(any(c['face_apex'] in v for _, v in chain) or c['face_apex'] in original_apexes)
                self.assertIn(c['group_apex'], exact.exact_group('2T'))
            for _, tetra, _ in r['terms']:
                self.assertTrue(all(p in exact.exact_group('2T') for p in tetra))
        self.assertEqual(seen, {'circle', 'digon'})

    def test_sphere_case_when_cone_point_lies_on_the_digon_sphere(self):
        # Found by the exhaustive 2I check: the original face cone point y = (1,1,1,1)
        # lies on the great sphere of the digon, so F - F' = k[S] and vol D_F == k pi^2.
        group = importlib.import_module('su2_exact_ade').exact_group('2I')
        r = self.api.global_polar_formula('2I', group[5], group[44], group[1])
        kinds = {c['type']: c.get('multiplicity') for c in r['corrections']}
        self.assertEqual(kinds.get('sphere'), 1)
        self.assert_elementary(r['oriented_volume_mod'])

    def test_rejects_non_members_and_non_integer_level(self):
        with self.assertRaises(ValueError):
            self.api.global_polar_formula('2T', (1/s.sqrt(2), 1/s.sqrt(2), 0, 0), (1, 0, 0, 0), (1, 0, 0, 0))
        with self.assertRaises(TypeError):
            self.api.global_polar_formula('2T', (1, 0, 0, 0), (1, 0, 0, 0), (1, 0, 0, 0), k=0.5)


if __name__ == '__main__':
    unittest.main()
