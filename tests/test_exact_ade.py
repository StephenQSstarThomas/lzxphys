import importlib
import importlib.util
import unittest
from collections import Counter
import sympy as s


class ExactADETests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('su2_exact_ade'),
                             'Exact E-group catalogue and global formula are missing')
        self.api = importlib.import_module('su2_exact_ade')

    def test_complete_2o_and_2i_catalogues(self):
        expected = {'2O':(84,10080,6135,{24:2,32:6,48:7,96:39,192:30}),
                    '2I':(563,212520,61299,{20:3,40:1,60:3,80:9,120:16,160:2,240:186,480:343})}
        for name,(types,full,singular,histogram) in expected.items():
            r = self.api.classify_group(name)
            self.assertEqual(len(r['types']),types)
            self.assertEqual(r['anchored_nondegenerate'],full)
            self.assertEqual(r['anchored_degenerate'],singular)
            self.assertEqual(dict(Counter(t['anchored_multiplicity'] for t in r['types'])),histogram)

    def test_catalogue_vertices_have_exact_unit_norm(self):
        for name,size in (('2T',24),('2O',48),('2I',120)):
            points = self.api.exact_group(name)
            self.assertEqual(len(set(points)),size)
            self.assertTrue(all(s.simplify(sum(x*x for x in p))==1 for p in points))

    def test_exact_irrational_geometric_examples(self):
        for name,det,tangent in (('2O',s.Rational(1,2),1/(3+s.sqrt(2))),
                                 ('2I',(s.sqrt(5)-1)/4,1/(3+2*s.sqrt(5)))):
            r = self.api.irrational_volume_certificate(name)
            self.assertEqual(s.simplify(r['determinant']-det),0)
            self.assertEqual(s.simplify(r['tan_half_area']-tangent),0)
            self.assertEqual(r['volume'],s.pi*s.atan(r['tan_half_area'])/2)
            self.assertEqual(r['root_of_unity_order_bound'],8 if name=='2O' else 4)

    def test_global_identity_and_nontrivial_degenerate_chain(self):
        e=(1,0,0,0); x=(0,1,0,0)
        for name in ('2T','2O','2I'):
            r=self.api.global_formula(name,e,x,x,expand=False)
            self.assertEqual(r['phase'],1)
        r=self.api.global_formula('2T',*((-1,0,0,0),)*3,expand=False)
        self.assertEqual(r['branch'],'cone')
        self.assertGreater(len(r['terms']),0)
        self.assertTrue(any(t['rank']==4 for t in r['terms']))
        self.assertLessEqual(len(r['terms']),24)

    def test_global_general_position_is_same_analytic_formula(self):
        r=self.api.global_formula('2O',(0,1,0,0),(0,0,0,1),(0,1,0,0))
        self.assertEqual(r['branch'],'direct')
        self.assertEqual(s.simplify(r['phase']/s.exp(-s.I*s.pi/8)),1)

    def test_expanded_center_keeps_exact_li2_without_forced_cas_reduction(self):
        r=self.api.global_formula('2T',*((-1,0,0,0),)*3,reduce=False)
        self.assertTrue(r['expanded'])
        self.assertIsNotNone(r['phase'])
        self.assertFalse(r['phase'].has(s.Float))
        self.assertTrue(all('formula' in t for t in r['terms'] if t['rank']==4))

    def test_membership_and_group_validation(self):
        with self.assertRaises(ValueError):
            self.api.exact_group('E9')
        with self.assertRaises(ValueError):
            self.api.global_formula('2T',(s.Rational(3,5),s.Rational(4,5),0,0),(0,1,0,0),(0,0,1,0))


if __name__=='__main__':
    unittest.main()
