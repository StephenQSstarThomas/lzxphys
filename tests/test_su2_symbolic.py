"""Exact expectations; no decimal approximation is evidence in these tests."""
import importlib
import importlib.util
import unittest

import sympy as s


class SymbolicFormulaTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('su2_symbolic'),
                             'Pure symbolic group-element formula is missing')
        self.api = importlib.import_module('su2_symbolic')

    def test_orthant_reduces_by_conjugation(self):
        r = self.api.group_formula((0,1,0,0), (0,0,0,1), (0,1,0,0))
        self.assertEqual(r['lengths'], (s.pi/2,)*6)
        self.assertEqual(r['determinant'], 1)
        self.assertEqual(r['q0'], -4-4*s.I)
        self.assertEqual(r['q1'], 12)
        self.assertEqual(r['z'], (1+s.I)/2)
        self.assertEqual(r['Sigma'], (-s.I*s.pi,)*6)
        self.assertEqual(r['raw_volume'], s.pi**2/8)
        self.assertEqual(s.simplify(r['phase']/s.exp(-s.I*s.pi/8)), 1)
        self.assertEqual(r['dilog_reduction']['remaining'], ())

    def test_matrix_trace_uses_cumulative_vertices(self):
        qs = ((0,1,0,0), (0,0,1,0), (0,s.Rational(3,5),0,s.Rational(4,5)))
        matrices = tuple(self.api.quaternion_matrix(q) for q in qs)
        r = self.api.group_formula(*matrices)
        self.assertEqual(r['determinant'], s.Rational(3,5))
        self.assertEqual(r['orientation_trace'], r['determinant'])
        self.assertEqual(r['cosines'], self.api.group_formula(*qs)['cosines'])
        A,B,C = matrices
        self.assertEqual(s.trace(A*(B*C-C*B))/4, s.Rational(4,5))

    def test_mirror_flips_only_oriented_phase(self):
        v = tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
        a = self.api.tetrahedron_formula(v)
        b = self.api.tetrahedron_formula((v[1],v[0],v[2],v[3]))
        self.assertEqual(a['raw_volume'], b['raw_volume'])
        self.assertEqual(s.simplify(a['phase']*b['phase']), 1)

    def test_rejects_floats_nonunit_and_degenerate(self):
        good = ((0,1,0,0), (0,0,0,1), (0,1,0,0))
        for bad in ((0.0,1,0,0), (0,2,0,0), (0,1,0), (s.I,0,0,0)):
            with self.subTest(bad=bad), self.assertRaises((TypeError,ValueError)):
                self.api.group_formula(bad, *good[1:])
        for k in (s.Rational(1,2), 1.0, True):
            with self.assertRaises(TypeError):
                self.api.group_formula(*good, k=k)
        with self.assertRaises(ValueError):
            self.api.group_formula((1,0,0,0), *good[1:])
        with self.assertRaises(ValueError):
            self.api.group_formula(*((-1,0,0,0),)*3)

    def test_exact_large_and_negative_levels(self):
        for k in (0,-1,10**80+1):
            r = self.api.group_formula((0,1,0,0), (0,0,0,1), (0,1,0,0), k=k)
            self.assertEqual(s.simplify(r['phase']/s.exp(-s.I*k*s.pi/8)), 1)
            self.assertFalse(r['raw_volume'].has(s.Float))

    def test_euler_duplication_and_unknown_domain(self):
        a = self.api.reduce_real_dilog(((1,s.Rational(1,3)), (1,s.Rational(2,3))))
        self.assertEqual(a['remaining'], ())
        self.assertEqual(s.simplify(a['expression']-s.pi**2/6+s.log(s.Rational(1,3))*s.log(s.Rational(2,3))), 0)
        b = self.api.reduce_real_dilog(((1,s.Rational(1,3)), (1,-s.Rational(1,3))))
        self.assertEqual(b['remaining'], ((s.Rational(1,2),s.Rational(1,9)),))
        z = s.Symbol('z')
        c = self.api.reduce_real_dilog(((1,z), (1,1-z)))
        self.assertEqual(len(c['remaining']), 2)

    def test_join_symbolic_integral_and_orientation(self):
        for a,b,w in ((s.pi/3,s.pi/4,s.pi**2/24),(-s.pi/3,s.pi/4,-s.pi**2/24)):
            r = self.api.orthogonal_join_certificate(a,b)
            self.assertEqual(r['oriented_volume'], w)
            self.assertEqual(r['rho_integral'], s.Rational(1,2))

    def test_join_requires_real_exact_angles(self):
        for a,b in ((s.I,s.pi/2),(s.pi/2,s.I),(s.Symbol('z'),s.pi/2)):
            with self.subTest(a=a,b=b), self.assertRaises(ValueError):
                self.api.orthogonal_join_certificate(a,b)
        alpha,beta=s.symbols('alpha beta',real=True)
        r=self.api.orthogonal_join_certificate(alpha,beta)
        self.assertNotEqual(r['condition'],s.true)
        self.assertEqual(r['oriented_volume'],alpha*beta/2)

    def test_explicit_straight_line_template_at_orthant(self):
        self.assertTrue(callable(getattr(self.api,'formula_template',None)))
        template=self.api.formula_template()
        env={s.Symbol('c'+str(j),real=True):s.S.Zero for j in range(1,7)}
        env.update({s.Symbol('detG',positive=True):s.S.One,
                    s.Symbol('orientation',integer=True):s.S.One,
                    s.Symbol('k',integer=True):s.S.One})
        for name,expression in template:
            env[name]=expression.xreplace(env)
        self.assertEqual(s.simplify(env[s.Symbol('q0')]+4+4*s.I),0)
        self.assertEqual(s.simplify(env[s.Symbol('z')]-(1+s.I)/2),0)
        self.assertEqual(len([name for name,_ in template if str(name).startswith('Z')]),8)


if __name__ == '__main__':
    unittest.main()
