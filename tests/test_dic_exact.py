import importlib
import importlib.util
import itertools
import unittest
from fractions import Fraction
import sympy as s


class DicTableProofTests(unittest.TestCase):
    """Exact all-n proof data for the Dic_n table (no floating point)."""

    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('su2_dic_exact'),
                             'Dic_n proof module is missing')
        self.api = importlib.import_module('su2_dic_exact')

    def test_formal_closure_on_the_whole_normaliser(self):
        residuals = self.api.formal_closure_certificate()
        self.assertEqual(len(residuals), 16)
        self.assertTrue(all(r == 0 for r in residuals.values()))

    def test_formal_normalisation(self):
        values = self.api.formal_normalisation_certificate()
        self.assertEqual(len(values), 12)
        self.assertTrue(all(v == 0 for v in values.values()))

    def test_same_table_as_the_historical_implementation(self):
        formula = importlib.import_module('scripts.dic_formula')
        subgroups = importlib.import_module('scripts.ade_subgroups')
        for n in (2, 3, 4):
            elements = subgroups.dic_elements(n)
            coords = {g: self.api.dic_element(n, (-1 if g.epsilon else 1)*g.r, g.epsilon) for g in elements}
            for a, b, c in itertools.product(elements, repeat=3):
                self.assertEqual(Fraction(formula.dic_volume_units(n, a, b, c), 2),
                                 self.api.table_exponent(coords[a], coords[b], coords[c]))
            for a, b in itertools.product(elements, repeat=2):
                self.assertEqual(coords[a*b], self.api.k_product(coords[a], coords[b]))

    def test_lens_pairing_for_every_cyclic_subgroup(self):
        # q|U(1) = -x floor(y+z): only j >= 1 contribute, total -(M-1)/M = 1/M mod 1.
        for M in range(1, 41):
            self.assertEqual(self.api.lens_pairing(M), Fraction(1, M) % 1)

    def test_explicit_fundamental_cycle_pairing(self):
        for n in range(2, 11):
            self.assertTrue(self.api.is_cycle_in_coinvariants(self.api.dic_fundamental_cycle(n)))
            self.assertEqual(self.api.dic_pairing(n), Fraction(1, 4*n))

    def test_geometric_pairings_agree_exactly(self):
        geometric = importlib.import_module('su2_polar_global')
        # Dic_2 = Q8 inside 2T: fundamental cycle of S^3/Q8.
        total = s.S.Zero
        for t, c in self.api.dic_fundamental_cycle(2):
            g = [self.api.dic_point(x) for x in self.api.inhomogeneous(t)]
            total += c*geometric.global_polar_formula('2T', *g)['oriented_volume_mod']
        self.assertEqual(s.simplify(total/(2*s.pi**2)), s.Rational(1, 8))
        # Z_4 inside 2T: lens cycle sum_j [u | u^j | u^-1].
        u = lambda j: (s.cos(s.pi*j/2), s.sin(s.pi*j/2), 0, 0)
        total = sum((geometric.global_polar_formula('2T', u(1), u(j), u(-1))['oriented_volume_mod']
                     for j in range(4)), s.S.Zero)
        self.assertEqual(s.simplify(total/(2*s.pi**2)), s.Rational(1, 4))

    def test_nondegenerate_dicyclic_tetrahedra_are_double_arcs(self):
        n = 4
        U = lambda r: self.api.dic_element(n, r, 0)
        V = lambda r: self.api.dic_element(n, r, 1)
        r = self.api.dic_nondegenerate_volume(n, U(1), V(0), U(3))
        self.assertIn(r['volume'], (s.pi**2/32, 3*s.pi**2/32))
        polar = importlib.import_module('su2_polar_dual')
        points = [self.api.dic_point(x) for x in r['cumulative']]
        exact = polar.tetrahedron_volume('2O', points)
        self.assertEqual(exact['oriented_volume'], r['oriented_volume'])
        with self.assertRaises(ValueError):
            self.api.dic_nondegenerate_volume(n, U(1), U(1), U(1))


if __name__ == '__main__':
    unittest.main()
