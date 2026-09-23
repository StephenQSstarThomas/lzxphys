"""Independent finite certificates; integer arithmetic only (no approximation)."""
from itertools import combinations, permutations
import unittest

import sympy as s
from su2_exact_ade import integer_pair_group, classify_group
from su2_2t_geometry import chamber_probes, hyperplane_normals, hurwitz_integer_vertices


def plus(*pairs):
    return tuple(sum(p[i] for p in pairs) for i in (0,1))


def neg(p):
    return -p[0],-p[1]


def times(p,q,d):
    return p[0]*q[0]+d*p[1]*q[1],p[0]*q[1]+p[1]*q[0]


def determinant_by_permutations(matrix,d):
    value=(0,0)
    for p in permutations(range(len(matrix))):
        term=(1,0)
        for i,j in enumerate(p):
            term=times(term,matrix[i][j],d)
        if sum(p[i]>p[j] for i,j in combinations(range(len(p)),2))%2:
            term=neg(term)
        value=plus(value,term)
    return value


class IndependentExactCertificates(unittest.TestCase):
    def test_every_product_in_three_groups(self):
        for name,d in (('2T',0),('2O',2),('2I',5)):
            points=integer_pair_group(name)
            members=set(points)
            for a in points:
                for b in points:
                    t=[[times(x,y,d) for y in b] for x in a]
                    numerators=(plus(t[0][0],neg(t[1][1]),neg(t[2][2]),neg(t[3][3])),
                                plus(t[0][1],t[1][0],neg(t[2][3]),t[3][2]),
                                plus(t[0][2],t[2][0],neg(t[3][1]),t[1][3]),
                                plus(t[0][3],t[3][0],neg(t[1][2]),t[2][1]))
                    self.assertTrue(all(x%4==0 for p in numerators for x in p))
                    self.assertIn(tuple(tuple(x//4 for x in p) for p in numerators),members)

    def test_all_catalogue_gram_determinants(self):
        for name in ('2T','2O','2I'):
            catalogue=classify_group(name)
            d=catalogue['field_d']
            points=integer_pair_group(name)
            for row in catalogue['types']:
                matrix=[[(16,0) if i==j else (0,0) for j in range(4)] for i in range(4)]
                for (i,j),code in zip(combinations(range(4),2),row['gram_key']):
                    matrix[i][j]=matrix[j][i]=catalogue['gram_alphabet'][code]
                determinant=determinant_by_permutations(matrix,d)
                detv=row['determinant_pair_over_64']
                self.assertEqual(determinant,tuple(16*x for x in times(detv,detv,d)))
                reps=[points[j] for j in row['representative_indices']]
                original=[[plus(*(times(x,y,d) for x,y in zip(a,b))) for b in reps] for a in reps]
                original_key=min(tuple(original[p[i]][p[j]] for i,j in combinations(range(4),2))
                                 for p in permutations(range(4)))
                self.assertEqual(original_key,tuple(catalogue['gram_alphabet'][j] for j in row['gram_key']))

    def test_four_wall_reflections_preserve_every_probe(self):
        probes=set(chamber_probes())
        for a,b,c,d in probes:
            self.assertIn((a,c,b,d),probes)
            self.assertIn((a,b,d,c),probes)
            self.assertIn((a,b,c,-d),probes)
            self.assertEqual((a-b-c-d)%2,0)
            t=(a-b-c-d)//2
            self.assertIn((a-t,b+t,c+t,d+t),probes)

    def test_every_hurwitz_face_normal_is_a_wall(self):
        roots=hyperplane_normals()
        for face in combinations(hurwitz_integer_vertices(),3):
            normal=[]
            for j in range(4):
                matrix=[[(((-1)**j if i==0 else 1)*row[k],0)
                         for k in range(4) if k!=j] for i,row in enumerate(face)]
                normal.append(determinant_by_permutations(matrix,0)[0])
            if not any(normal):
                continue
            self.assertTrue(any(all(normal[i]*root[j]==normal[j]*root[i]
                                    for i,j in combinations(range(4),2)) for root in roots))

    def test_join_principal_log_sign_at_orthant(self):
        # This fixture caught the wrong sign in the first manuscript draft.
        p=q=s.I
        z=s.cancel(-2/(p*q+p+q-1))
        LU,LV,LX,LY=[s.expand_complex(s.log(s.expand(w)))
                      for w in (1-z,1+p*q*z,1+p*z,1+q*z)]
        self.assertEqual(s.simplify(LU+LV-LX-LY),-s.I*s.pi)
        self.assertEqual(s.simplify(LX-LV),s.I*s.pi/2)
        self.assertEqual(s.simplify(LY-LV),s.I*s.pi/2)


if __name__=='__main__':
    unittest.main()
