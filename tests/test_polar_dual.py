import importlib
import importlib.util
import unittest
from fractions import Fraction
import sympy as s


class PolarDualTests(unittest.TestCase):
    """Exact tests only: no floating point enters any assertion."""

    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('su2_polar_dual'),
                             'Polar-dual exact volume module is missing')
        self.api = importlib.import_module('su2_polar_dual')

    def test_reflection_group_orders(self):
        for name, order in (('2T', 192), ('2O', 1152), ('2I', 14400)):
            self.assertEqual(self.api.weyl_chamber_masks(name)[0], order)

    def test_rational_angle_table_is_exact(self):
        for d in (0, 2, 5):
            for (a, b), r in self.api.rational_cos2_table(d).items():
                value = s.Rational(a.numerator, a.denominator)+s.Rational(b.numerator, b.denominator)*s.sqrt(d)
                exact = s.cos(s.pi*s.Rational(r.numerator, r.denominator))**2
                self.assertEqual(s.simplify(s.expand(exact)-value), 0)

    def test_group_inner_products_are_rational_angles(self):
        for name in ('2T', '2O', '2I'):
            points = self.api.group_vectors(name)
            values = {self.api.dot(points[0], p) for p in points}
            for c in values:
                q = self.api.edge_fraction(c)
                exact = s.cos(s.pi*s.Rational(q.numerator, q.denominator))
                self.assertEqual(s.simplify(exact-c.sympy()), 0)

    def _bracket(self, T, d):
        """Exact rational bounds lo < theta/pi < hi for theta = arctan(sqrt T)."""
        grid = [(Fraction(0), self.api.Surd(0, 0, d))]
        for (a, b), r in sorted(self.api.rational_cos2_table(d).items(), key=lambda kv: kv[1]):
            if 0 < r < Fraction(1, 2):
                c2 = self.api.Surd(a, b, d)
                grid.append((r, (1-c2)/c2))
        grid.append((Fraction(1, 2), None))
        for (r0, t0), (r1, t1) in zip(grid, grid[1:]):
            if T > t0 and (t1 is None or T < t1):
                return r0, r1
            if t1 is not None and T == t1:
                self.fail('certificate angle is rational')
        self.fail('no bracket')

    def test_angle_certificate_independently(self):
        """Re-verify every identity by exp(2 i theta) products and exact brackets."""
        x = s.Symbol('x')
        for name, rank in (('2O', 1), ('2I', 3)):
            d = self.api.FIELD[name]
            system = self.api.angle_system(name)
            self.assertEqual(len(system['basis']), rank)
            for _, rel, r in system['certificate']:
                product, lo, hi = s.S.One, Fraction(0), Fraction(0)
                for key, c in rel.items():
                    T = self.api.Surd(key[0], key[1], d)
                    t = s.sqrt(T.sympy())
                    product *= ((1+s.I*t)**2/(1+T.sympy()))**int(c)
                    b0, b1 = self._bracket(T, d)
                    lo += c*b0 if c > 0 else c*b1
                    hi += c*b1 if c > 0 else c*b0
                self.assertTrue(lo < r < hi and hi-lo < 1)
                target = s.exp(2*s.I*s.pi*s.Rational(r.numerator, r.denominator))
                self.assertEqual(s.minimal_polynomial(s.expand(product-target), x), x)

    def test_every_dihedral_angle_is_certified(self):
        for name in ('2O', '2I'):
            system = self.api.angle_system(name)
            values = self.api._all_irrational_tan2(name)
            self.assertEqual(set(values), set(system['angles']))

    def test_catalogue_and_rational_types(self):
        expected = {'2T': (12, 12), '2O': (84, 60), '2I': (563, 268)}
        for name, (types, rational) in expected.items():
            rows = self.api.polar_catalogue(name)['types']
            self.assertEqual(len(rows), types)
            self.assertEqual(sum(r['rational_volume'] for r in rows), rational)
            for r in rows:
                self.assertTrue(0 < r['rational_part'] <= 1 or any(r['basis_coefficients']))

    def test_2t_agrees_with_independent_f4_chambers(self):
        geometry = importlib.import_module('su2_2t_geometry')
        points = importlib.import_module('su2_exact_ade').exact_group('2T')
        for r in self.api.polar_catalogue('2T')['types']:
            f4 = geometry.tetrahedron_2t([points[i] for i in r['representative_indices']])
            self.assertEqual(f4['volume_over_pi2'], r['rational_part'])

    def test_orthant_and_orientation(self):
        e = [tuple(int(i == j) for i in range(4)) for j in range(4)]
        r = self.api.tetrahedron_volume('2T', e)
        self.assertEqual(r['volume'], s.pi**2/8)
        m = self.api.tetrahedron_volume('2T', (e[1], e[0], e[2], e[3]))
        self.assertEqual(m['oriented_volume'], -r['oriented_volume'])
        p = self.api.polar_formula('2O', (0, 1, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0), k=3)
        self.assertEqual(s.simplify(p['phase']-s.exp(-3*s.I*s.pi/8)), 0)

    def test_irrational_certificates_exactly(self):
        r2, phi = s.sqrt(2), (1+s.sqrt(5))/2
        half = s.Rational(1, 2)
        o = self.api.polar_formula('2O', (0, 1, 0, 0), (1/r2, 0, 0, 1/r2), (half, half, -half, half))
        beta = s.atan(1/r2)
        self.assertEqual(s.expand(o['volume']-s.pi/2*(beta-s.pi/8)), 0)
        self.assertEqual(s.simplify(s.expand_trig(s.tan(beta-s.pi/8))-1/(3+r2)), 0)
        i = self.api.polar_formula('2I', (0, 1, 0, 0), (0, 0, 0, 1), (phi/2, (phi-1)/2, 0, -half))
        A, B = s.atan(half), s.atan(1/s.sqrt(5))
        theta = s.pi/4-B-A/2                      # claimed = arctan(1/(4 phi+1))
        self.assertEqual(s.expand(i['volume']-s.pi/2*theta), 0)
        t = 1/(4*phi+1)                           # tan(2 theta) = cot(2B + A)
        self.assertEqual(s.simplify(s.expand_trig(1/s.tan(2*B+A))-2*t/(1-t*t)), 0)

    def test_rejections(self):
        with self.assertRaises(ValueError):
            self.api.polar_formula('2T', (1, 0, 0, 0), (0, 1, 0, 0), (0, 1, 0, 0))
        with self.assertRaises(ValueError):
            self.api.polar_formula('2T', (0, 1, 0, 0), (0, 0, 0, 1), (1/s.sqrt(2), 1/s.sqrt(2), 0, 0))
        with self.assertRaises(TypeError):
            self.api.polar_formula('2T', (0, 1, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0), k=s.Rational(1, 2))


if __name__ == '__main__':
    unittest.main()
