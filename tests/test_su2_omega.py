"""Independent fixtures catch normalization, orientation, and branch errors."""
import unittest
from fractions import Fraction
from collections import Counter

import mpmath as mp
import numpy as np

from scripts.su2_audit_checks import PAULI, exp_vector
from su2_omega import omega_su2, oriented_volume, omega_quaternions, omega_quaternions_source, quaternion_multiply, quaternion
from su2_omega import face_chain, tetrahedron_chain, hemisphere_safe


class ClosedFormulaTests(unittest.TestCase):
    def test_high_precision_matrix_input_is_preserved(self):
        with mp.workdps(100):
            qs = [(mp.cos(a), *(mp.sin(a)*int(j == i) for j in range(3)))
                  for i, a in enumerate(map(mp.mpf, ('.7', '.8', '.9')))]
            matrices = [mp.matrix([[a+1j*z, y+1j*x], [-y+1j*x, a-1j*z]]) for a,x,y,z in qs]
            expected = omega_quaternions(*qs, dps=80)
            self.assertLess(abs(omega_su2(*matrices, dps=80)-expected), mp.mpf('1e-75'))

    def test_extremely_small_volume_keeps_relative_precision(self):
        eps = Fraction(1, 10**100)
        vertices = [(1,0,0,0)]
        for g in [(1,eps,0,0), (1,0,eps,0), (1,0,0,eps)]:
            vertices.append(quaternion_multiply(vertices[-1], g))
        with mp.workdps(80):
            expected_leading = mp.mpf('1e-300')/6
            self.assertLess(abs(oriented_volume(vertices, dps=50)/expected_leading-1), mp.mpf('1e-45'))

    def test_nonfinite_quaternions_are_rejected(self):
        for bad in (mp.inf, -mp.inf, mp.nan):
            with self.assertRaises(ValueError):
                quaternion((1, bad, 0, 0))

    def test_large_integer_level_keeps_exact_period(self):
        # This triple gives the positive orthant: omega_k=exp(i*k*pi/8).
        x, z = (0, 1, 0, 0), (0, 0, 0, 1)
        with mp.workdps(55):
            self.assertLess(abs(omega_quaternions(x, z, x, k=16*10**80, dps=50)-1), mp.mpf('1e-45'))

    def test_orthant_and_reversed_orientation(self):
        # An orthant is 1/16 of S^3; swapping vertices reverses orientation.
        vertices = np.eye(4).tolist()
        with mp.workdps(50):
            self.assertLess(abs(oriented_volume(vertices) - mp.pi**2 / 8), mp.mpf('1e-45'))
            vertices[1], vertices[2] = vertices[2], vertices[1]
            self.assertLess(abs(oriented_volume(vertices) + mp.pi**2 / 8), mp.mpf('1e-45'))

    def test_two_parameter_exact_volume(self):
        with mp.workdps(70):
            for alpha, beta in [(mp.mpf('1.1'), mp.mpf('2')), (mp.mpf('.2'), mp.mpf('2.9'))]:
                vertices = [[1, 0, 0, 0], [mp.cos(alpha), mp.sin(alpha), 0, 0],
                            [0, 0, 1, 0], [0, 0, mp.cos(beta), mp.sin(beta)]]
                self.assertLess(abs(oriented_volume(vertices, dps=55) - alpha * beta / 2), mp.mpf('1e-50'))

    def test_general_matrix_fixture(self):
        matrices = [exp_vector(axis * a) for axis, a in zip(np.eye(3), [.7, .8, .9])]
        result = complex(omega_su2(*matrices, 1))
        self.assertLess(abs(result - np.exp(1j * .022543184393801173)), 2e-15)

    def test_source_note_convention_is_inverse_representative(self):
        x, z = (0, 1, 0, 0), (0, 0, 0, 1)
        with mp.workdps(55):
            old = omega_quaternions(x, z, x, k=1)
            source = omega_quaternions_source(x, z, x, k=1)
            self.assertLess(abs(source * old - 1), mp.mpf('1e-45'))

    def test_invalid_matrix_and_level_are_rejected(self):
        with self.assertRaises(ValueError):
            omega_su2(np.eye(2) * 2, np.eye(2), np.eye(2), 1)
        with self.assertRaises((TypeError, ValueError)):
            omega_su2(np.eye(2), np.eye(2), np.eye(2), .5)

    def test_non_two_by_two_shapes_are_rejected_for_all_matrix_types(self):
        for bad in (mp.eye(3), np.eye(3), [[1,0,7],[0,1,8]], [[1,0],[0,1],[9,9]]):
            with self.subTest(matrix_type=type(bad).__name__):
                with self.assertRaises(ValueError):
                    omega_su2(bad,mp.eye(2),mp.eye(2))


class GlobalFormulaTests(unittest.TestCase):
    def test_positive_dependence_without_antipodes(self):
        # Four points surround zero in R^3, with no antipodal pair: a hemispherical
        # filling still has phase -1, despite the zero original determinant.
        triple = ((0,1,0,0), (0,0,0,1), (-1,0,1,1))
        with mp.workdps(50):
            self.assertLess(abs(omega_quaternions(*triple)+1), mp.mpf('1e-45'))

    def test_exact_chain_boundary_and_left_equivariance(self):
        def collect(terms):
            result = Counter()
            for sign, simplex in terms:
                if not any(a == b for a,b in zip(simplex,simplex[1:])):
                    result[simplex] += sign
            return {simplex: sign for simplex,sign in result.items() if sign}

        def boundary(chain):
            return collect((sign*(-1)**i, v[:i]+v[i+1:])
                           for sign,v in chain for i in range(len(v)))

        examples = [[(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)],
                    [(1,0,0,0), (-1,0,0,0), (1,0,0,0), (-1,0,0,0)],
                    [(1,0,0,0), (0,1,0,0), (0,0,1,0), (-1,-1,-1,0)]]
        for raw in examples:
            vertices = tuple(map(quaternion,raw))
            chain = tetrahedron_chain(vertices)
            expected = collect(((-1)**i*sign, v)
                               for i in range(4)
                               for sign,v in face_chain(*(vertices[:i]+vertices[i+1:])))
            self.assertEqual(boundary(chain), expected)
            self.assertLessEqual(len(chain),24)
            self.assertTrue(all(hemisphere_safe(v) for _,v in chain))
            left = (1,2,-1,3)
            translate = lambda v: quaternion_multiply(left,v)
            translated = tetrahedron_chain(tuple(map(translate,vertices)))
            self.assertEqual(collect(translated), collect((sign,tuple(map(translate,v))) for sign,v in chain))

    def test_center_is_nontrivial_at_odd_level(self):
        # Restriction to the central C2 detects the odd integral level.
        center = (-1, 0, 0, 0)
        with mp.workdps(50):
            for k in (-1, 0, 1, 2):
                self.assertLess(abs(omega_quaternions(center, center, center, k) - (-1)**k), mp.mpf('1e-45'))

    def test_normalization_including_antipodes(self):
        identity, center, element = (1, 0, 0, 0), (-1, 0, 0, 0), (1, 2, 3, 4)
        for triple in [(identity, center, element), (center, identity, element), (center, element, identity)]:
            self.assertEqual(omega_quaternions(*triple), 1)

    def test_safe_lower_dimensional_volume(self):
        self.assertEqual(oriented_volume([(1, 0, 0, 0), (1, 1, 0, 0),
                                         (1, 0, 1, 0), (1, 1, 1, 0)]), 0)

    def test_exact_pentagon_at_exceptional_inputs(self):
        # Exact products retain relations; floating matrix products cannot do so.
        c, x, y, h = (-1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (1, 2, -1, 3)
        mul = quaternion_multiply
        with mp.workdps(45):
            for a, b, c1, d in [(c, c, c, c), (x, x, y, h), (h, c, c, y), (x, y, x, y)]:
                w = lambda g, h1, j: omega_quaternions(g, h1, j, dps=45)
                lhs = w(b, c1, d) * w(a, mul(b, c1), d) * w(a, b, c1)
                rhs = w(mul(a, b), c1, d) * w(a, b, mul(c1, d))
                self.assertLess(abs(lhs-rhs), mp.mpf('1e-39'))


if __name__ == '__main__':
    unittest.main()
