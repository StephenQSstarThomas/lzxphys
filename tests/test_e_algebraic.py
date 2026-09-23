"""Exact counting representative: expectations from independent hand examples."""
import importlib
import importlib.util
import unittest

import mpmath as mp
import sympy as sp

from scripts.check_ade_phase import analytic_inputs


class AlgebraicPhaseTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('scripts.ade_count'),
                             'The exact E-type counting representative is missing')
        self.api = importlib.import_module('scripts.ade_count')

    def test_three_hand_counted_examples(self):
        # Detect wrong quaternion convention, probe ordering, or volume substitution.
        for name, count, hits in [('2T', 1, [0]), ('2O', 2, [0, 28]),
                                  ('2I', 4, [0, 64, 80, 112])]:
            with self.subTest(name=name):
                result = self.api.count_details(name, *analytic_inputs(name))
                self.assertEqual(result['count'], count)
                self.assertEqual(result['terms'][0]['hit_indices'], hits)
                self.assertEqual(result['terms'][0].get('hit_elements', [None])[0],
                                 ['1', '0', '0', '0'])
                self.assertEqual(result['exponent_mod_order'], -count % result['order'])

    def test_center_is_not_discarded_as_zero_determinant(self):
        for name, count in [('2T', 12), ('2O', 24), ('2I', 60)]:
            with self.subTest(name=name):
                result = self.api.count_details(name, *((-1, 0, 0, 0),)*3)
                self.assertEqual(result['count'] % result['order'], count)
                self.assertEqual(result['branch'], 'cone')

    def test_identity_normalization_in_each_position(self):
        for name in ('2T', '2O', '2I'):
            triple = analytic_inputs(name)
            for j in range(3):
                args = triple[:j] + ((1, 0, 0, 0),) + triple[j+1:]
                self.assertEqual(self.api.count_details(name, *args)['count'], 0)

    def test_level_is_exact_and_phase_keeps_requested_precision(self):
        triple = analytic_inputs('2I')
        for k in [0, -7, 10**120 + 1]:
            result = self.api.count_details('2I', *triple, k=k)
            self.assertEqual(result['exponent_mod_order'], (-4*k) % 120)
            with mp.workdps(15):
                phase = self.api.algebraic_phase('2I', *triple, k=k, dps=55)
            with mp.workdps(65):
                expected = mp.exp(2j*mp.pi*((-4*k) % 120)/120)
                self.assertLess(abs(phase-expected), mp.mpf('1e-53'))

    def test_reject_nonintegral_level_and_nonmember_inputs(self):
        triple = analytic_inputs('2T')
        with self.assertRaises(TypeError):
            self.api.count_details('2T', *triple, k=0.5)
        for bad in [(2, 0, 0, 0), (1, 0, 0),
                    (sp.Rational(3, 5), sp.Rational(4, 5), 0, 0),
                    (0.5, 0.5, 0.5, 0.5), (sp.Symbol('x'), 0, 0, 0),
                    (1+sp.Symbol('s0'), 0, 0, 0)]:
            with self.subTest(bad=bad), self.assertRaises((ValueError, TypeError)):
                self.api.count_details('2T', bad, *triple[1:])

    def test_positive_radicals_and_symbolic_fields_agree(self):
        for name, d in [('2O', 2), ('2I', 5)]:
            triple = analytic_inputs(name)
            actual = tuple(tuple(sp.sympify(x).subs(sp.Symbol('s'+str(d)), sp.sqrt(d))
                                 for x in g) for g in triple)
            self.assertEqual(self.api.count_details(name, *actual)['count'],
                             {'2O': 2, '2I': 4}[name])

    def test_gauge_coboundary_has_the_derived_sign(self):
        # A sign-reversed beta, or reducing counts before averaging, breaks this.
        self.assertTrue(callable(getattr(self.api, 'gauge_details', None)))
        from scripts.ade_subgroups import sympy_quaternion_multiply as mul
        a, b, c = analytic_inputs('2T')
        with mp.workdps(60):
            def beta(x, y):
                report = self.api.gauge_details('2T', x, y, dps=42)
                self.assertEqual(len(report['terms']), 24)
                return mp.mpc(report['beta']['real'], report['beta']['imag'])
            delta = beta(b, c)*beta(a, mul(b, c))/(beta(mul(a, b), c)*beta(a, b))
            self.assertLess(abs(delta-mp.exp(-5j*mp.pi/96)), mp.mpf('1e-38'))
            identity = self.api.gauge_details('2T', (1, 0, 0, 0), a, dps=42)
            self.assertEqual(mp.mpc(identity['beta']['real'], identity['beta']['imag']), 1)

    def test_c4_bar_cycle_calibrates_source_A_sign(self):
        x = (0, 1, 0, 0)
        powers = ((1, 0, 0, 0), x, (-1, 0, 0, 0), (0, -1, 0, 0))
        for name, order in [('2T', 24), ('2O', 48), ('2I', 120)]:
            exponent = sum(self.api.count_details(name, x, h, x)['exponent_mod_order']
                           for h in powers) % order
            self.assertEqual(exponent, order//4)  # exactly +i, not -i

    def test_exact_pentagon_exercises_non_Q8_inputs(self):
        from scripts.ade_subgroups import sympy_quaternion_multiply as mul
        from scripts.ade_phase import ade_group
        for name, order in [('2T', 24), ('2O', 48), ('2I', 120)]:
            a, b, c = analytic_inputs(name)
            d = ade_group(name)[-1]
            A = lambda *g: self.api.count_details(name, *g)['count']
            self.assertEqual((A(b, c, d)-A(mul(a, b), c, d)+A(a, mul(b, c), d)
                              - A(a, b, mul(c, d))+A(a, b, c)) % order, 0)


if __name__ == '__main__':
    unittest.main()
