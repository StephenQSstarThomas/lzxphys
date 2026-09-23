"""Human-review regressions: group inputs, edge formula, exact E examples."""
import unittest

import mpmath as mp
import sympy as sp

import su2_omega as su2
import scripts.ade_phase as ade
from scripts.ade_phase import phase_from_exact_quaternions
from scripts.ade_subgroups import DicElement, dic_quaternion, dic_quaternion_exact_q8, normalized_numeric
from scripts.dic_formula import dic_phase


class EdgeFormulaTests(unittest.TestCase):
    def test_orthant_volume_and_reflection(self):
        self.assertTrue(callable(getattr(su2, 'oriented_volume_from_edges', None)))
        vertices = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))
        with mp.workdps(65):
            self.assertLess(abs(su2.oriented_volume_from_edges(vertices, dps=55)
                                - mp.pi**2/8), mp.mpf('1e-53'))
            reflected = vertices[:3] + ((0, 0, 0, -1),)
            self.assertLess(abs(su2.oriented_volume_from_edges(reflected, dps=55)
                                + mp.pi**2/8), mp.mpf('1e-53'))

    def test_three_group_elements_determine_ordered_edges(self):
        self.assertTrue(callable(getattr(su2, 'edge_lengths_quaternions', None)))
        # p=(e, e1, e2, (e2+e3)/sqrt(2)); the 23 edge is pi/4.
        with mp.workdps(65):
            lengths = su2.edge_lengths_quaternions(
                (0, 1, 0, 0), (0, 0, 0, 1), (1, 1, 0, 0), dps=55)
            expected = (mp.pi/4,) + (mp.pi/2,)*5
            self.assertLess(max(abs(a-b) for a, b in zip(lengths, expected)),
                            mp.mpf('1e-53'))

    def test_global_edge_phase_keeps_center_branch(self):
        self.assertTrue(callable(getattr(su2, 'omega_quaternions_edge', None)))
        center = (-1, 0, 0, 0)
        with mp.workdps(55):
            self.assertLess(abs(su2.omega_quaternions_edge(center, center, center,
                                                         dps=50) + 1), mp.mpf('1e-48'))

    def test_source_phase_preserves_requested_precision_at_default_context(self):
        triple = ((0, 1, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0))
        with mp.workdps(15):
            value = su2.omega_quaternions_source(*triple, dps=55)
        with mp.workdps(65):
            expected = 1/su2.omega_quaternions(*triple, dps=55)
            self.assertLess(abs(value-expected), mp.mpf('1e-53'))


class ADEFormulaTests(unittest.TestCase):
    def test_symbolic_group_coordinates_embed_positive_square_roots(self):
        s2 = sp.Symbol('s2')
        point = normalized_numeric((s2/2, s2/2, 0, 0), dps=55)
        with mp.workdps(65):
            self.assertLess(abs(point[0]-1/mp.sqrt(2)), mp.mpf('1e-53'))
            self.assertLess(abs(sum(x*x for x in point)-1), mp.mpf('1e-53'))

    def test_dic_phase_rejects_fractional_level(self):
        identity = DicElement(2, 0, 0)
        with self.assertRaises(TypeError):
            dic_phase(2, identity, identity, identity, k=0.5)

    def test_audit_exposes_exact_input_and_nonzero_center_filling(self):
        self.assertTrue(callable(getattr(ade, 'phase_details', None)))
        center = (-1, 0, 0, 0)
        report = ade.phase_details(center, center, center, dps=50)
        self.assertEqual(report['branch'], 'cone')
        self.assertEqual(report['inputs'][0], ['-1', '0', '0', '0'])
        self.assertGreater(len(report['terms']), 1)
        self.assertTrue(any(term['orientation'] != 0 for term in report['terms']))
        with mp.workdps(55):
            phase = mp.mpc(report['phase']['real'], report['phase']['imag'])
            self.assertLess(abs(phase+1), mp.mpf('1e-48'))

    def test_dic_coordinates_represent_a_power_times_b(self):
        for r in range(4):
            element = DicElement(2, r, 1)
            self.assertLess(max(abs(float(x)-float(y)) for x, y in zip(
                dic_quaternion(2, element), dic_quaternion_exact_q8(element))), 1e-14)

    def test_e_phase_rejects_fractional_level(self):
        with self.assertRaises(TypeError):
            phase_from_exact_quaternions((1, 0, 0, 0), (1, 0, 0, 0),
                                         (1, 0, 0, 0), k=0.5)

    def test_e_group_examples_have_independent_geometric_values(self):
        # Orthant subdivision (2T), bisected orthant (2O), and spherical join (2I).
        x, z = (0, 1, 0, 0), (0, 0, 0, 1)
        half = sp.Rational(1, 2)
        s2, s5 = sp.Symbol('s2'), sp.Symbol('s5')
        examples = ((half, half, -half, -half), (s2/2, s2/2, 0, 0),
                    ((1+s5)/4, (s5-1)/4, 0, -half))
        with mp.workdps(65):
            phi = (1+mp.sqrt(5))/2
            phase_angles = (-mp.pi/32, -mp.pi/16, -mp.atan(1/(4*phi+1))/2)
            for third, angle in zip(examples, phase_angles):
                with self.subTest(third=third):
                    phase = phase_from_exact_quaternions(x, z, third, dps=55)
                    self.assertLess(abs(phase-mp.exp(1j*angle)), mp.mpf('1e-49'))


if __name__ == '__main__':
    unittest.main()
