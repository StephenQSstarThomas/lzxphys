import importlib
import importlib.util
import unittest
from fractions import Fraction
import sympy as s


class HurwitzGeometryTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('su2_2t_geometry'),
                             'Exact 2T geometric solver is missing')
        self.api = importlib.import_module('su2_2t_geometry')

    def test_hand_volumes_and_mirror(self):
        e = tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
        h = (s.Rational(1,2),)*4
        for vertices,count in ((e,72), ((e[0],e[1],e[2],h),18)):
            r = self.api.tetrahedron_2t(vertices)
            self.assertEqual(r['chambers'], count)
            self.assertEqual(r['volume_over_pi2'], Fraction(count,576))
            mirrored = self.api.tetrahedron_2t((vertices[1],vertices[0],*vertices[2:]))
            self.assertEqual(mirrored['oriented_volume_over_pi2'], -r['oriented_volume_over_pi2'])

    def test_complete_twelve_type_distribution(self):
        result = self.api.classify_2t()
        self.assertEqual(result['total_subsets'], 10626)
        self.assertEqual(result['degenerate_subsets'], 5586)
        self.assertEqual(result['nondegenerate_subsets'], 5040)
        self.assertEqual([row['chambers'] for row in result['types']],
                         [270,84,120,204,36,66,60,30,24,72,18,12])
        self.assertEqual([row['multiplicity'] for row in result['types']],
                         [192,576,576,288,576,576,576,576,576,48,192,288])

    def test_chambers_are_regular_and_have_distinct_signatures(self):
        points = self.api.chamber_probes()
        roots = self.api.hyperplane_normals()
        self.assertEqual(len(points),1152)
        self.assertEqual(len(roots),24)
        signatures = set()
        for p in points:
            products = tuple(sum(x*y for x,y in zip(p,n)) for n in roots)
            self.assertNotIn(0,products)
            self.assertEqual(sum(x*x for x in p),78)
            signatures.add(tuple(x>0 for x in products))
        self.assertEqual(len(signatures),1152)

    def test_actual_phase_is_not_old_orbit_count_representative(self):
        r = self.api.phase_2t((0,1,0,0),(0,0,0,1),(s.Rational(1,2),s.Rational(1,2),-s.Rational(1,2),-s.Rational(1,2)))
        self.assertEqual(r['oriented_volume_over_pi2'], Fraction(1,32))
        self.assertEqual(s.simplify(r['phase']/s.exp(-s.I*s.pi/32)),1)

    def test_nonmember_and_singular_not_silently_zero(self):
        e = tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
        with self.assertRaises(ValueError):
            self.api.tetrahedron_2t((e[0],e[0],e[2],e[3]))
        with self.assertRaises(ValueError):
            self.api.tetrahedron_2t((e[0],(s.Rational(3,5),s.Rational(4,5),0,0),e[2],e[3]))


if __name__ == '__main__':
    unittest.main()
