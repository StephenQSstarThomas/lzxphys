"""Exact E-group catalogues and the original global geometric representative.

All finite-group enumeration uses integer pairs; all analytic evaluation uses
SymPy expressions. No numerical evaluator, cohomology transfer, or fitting.
The E/F/T choices here reproduce the existing geometric rules with real radical
coordinates; singular-input cone vertices are NOT assumed to be group vertices.
"""
from copy import deepcopy
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb

import sympy as s
from su2_symbolic import (IDENTITY, exact_level, quaternion_product,
                         tetrahedron_formula, unit_quaternion)
from su2_2t_geometry import canonical_gram_key


def _add(a,b):
    return a[0]+b[0],a[1]+b[1]


def _sub(a,b):
    return a[0]-b[0],a[1]-b[1]


def _mul(a,b,d):
    return a[0]*b[0]+d*a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def _det3(a,b,c,d):
    return _add(_sub(_mul(a[0],_sub(_mul(b[1],c[2],d),_mul(b[2],c[1],d)),d),
                     _mul(a[1],_sub(_mul(b[0],c[2],d),_mul(b[2],c[0],d)),d)),
                _mul(a[2],_sub(_mul(b[0],c[1],d),_mul(b[1],c[0],d)),d))


@lru_cache(maxsize=3)
def integer_pair_group(name):
    """Coordinates mean (a+b sqrt(d))/4; Gram numerators have denominator 16."""
    if name not in ('2T','2O','2I'):
        raise ValueError('Use 2T, 2O, or 2I')
    points = []
    for j in range(4):
        for sign in (1,-1):
            q=[(0,0)]*4
            q[j]=(4*sign,0)
            points.append(tuple(q))
    points.extend(tuple((2*t,0) for t in signs) for signs in product((1,-1),repeat=4))
    if name=='2O':
        for i,j in combinations(range(4),2):
            for u,v in product((1,-1),repeat=2):
                q=[(0,0)]*4
                q[i],q[j]=(0,2*u),(0,2*v)
                points.append(tuple(q))
    if name=='2I':
        base=((0,0),(2,0),(1,1),(-1,1))
        for p in permutations(range(4)):
            if sum(p[i]>p[j] for i,j in combinations(range(4),2))%2:
                continue
            for signs in product((1,-1),repeat=3):
                it=iter(signs)
                q=[]
                for i in p:
                    a,b=base[i]
                    t=next(it) if (a,b)!=(0,0) else 1
                    q.append((t*a,t*b))
                points.append(tuple(q))
    return tuple(points)


@lru_cache(maxsize=3)
def exact_group(name):
    points=integer_pair_group(name)
    d={'2T':0,'2O':2,'2I':5}[name]
    return tuple(tuple(s.Rational(a,4)+s.Rational(b,4)*s.sqrt(d) for a,b in q) for q in points)


@lru_cache(maxsize=3)
def _classify(name):
    points=integer_pair_group(name)
    d={'2T':0,'2O':2,'2I':5}[name]
    dot=[[tuple(sum(_mul(a,b,d)[r] for a,b in zip(u,v)) for r in range(2))
          for v in points] for u in points]
    alphabet=tuple(sorted({z for row in dot for z in row}))
    code={z:i for i,z in enumerate(alphabet)}
    G=[[code[z] for z in row] for row in dot]
    classes={}
    singular=0
    for a,b,c in combinations(range(1,len(points)),3):
        det=_det3(points[a][1:],points[b][1:],points[c][1:],d)
        if det==(0,0):
            singular+=1
            continue
        key=canonical_gram_key((G[0][a],G[0][b],G[0][c],G[a][b],G[a][c],G[b][c]))
        if key not in classes:
            # Canonical key order is 01,02,03,12,13,23, NOT formula edge order.
            classes[key]={'gram_key':key,'anchored_multiplicity':0,
                          'representative_indices':(0,a,b,c),'determinant_pair_over_64':det,
                          'cosines_in_formula_order':tuple(alphabet[key[j]] for j in (0,1,2,5,4,3))}
        classes[key]['anchored_multiplicity']+=1
    count=sum(row['anchored_multiplicity'] for row in classes.values())
    for row in classes.values():
        row['all_subsets_multiplicity']=row['anchored_multiplicity']*len(points)//4
    return {'group':name,'order':len(points),'field_d':d,'gram_denominator':16,
            'gram_alphabet':alphabet,'total_subsets':comb(len(points),4),
            'anchored_nondegenerate':count,'anchored_degenerate':singular,
            'nondegenerate_subsets':count*len(points)//4,
            'nondegenerate_input_triples':count*6,
            'degenerate_input_triples':len(points)**3-count*6,
            'types':tuple(classes[key] for key in sorted(classes))}


def classify_group(name):
    return deepcopy(_classify(name))


def _clean_point(point):
    return tuple(s.simplify(x) for x in point)


@lru_cache(maxsize=32768)
def _rank(vertices):
    return s.Matrix.hstack(*(s.Matrix(v) for v in vertices)).rank()


@lru_cache(maxsize=32768)
def _safe(vertices):
    if _rank(vertices)==len(vertices):
        return True
    for size in range(2,len(vertices)+1):
        for subset in combinations(vertices,size):
            if _rank(subset)!=size-1:
                continue
            kernel=s.Matrix.hstack(*(s.Matrix(v) for v in subset)).nullspace()
            if len(kernel)!=1:
                continue
            signs=[s.sign(s.simplify(x)) for x in kernel[0] if x!=0]
            if any(x not in (-1,1) for x in signs):
                raise ArithmeticError('Undecidable exact convex-hull sign')
            if all(x==1 for x in signs) or all(x==-1 for x in signs):
                return False
    return True


def _radial(vertices):
    return [] if any(a==b for a,b in zip(vertices,vertices[1:])) else [(1,vertices)]


def _edge(a,b):
    if a==b:
        return []
    if _safe((a,b)):
        return [(1,(a,b))]
    midpoint=_clean_point(quaternion_product(a,(0,1,0,0)))
    return _radial((a,midpoint))+_radial((midpoint,b))


def _cone(anchor,chain):
    ranks=[_rank(v) for _,v in chain]
    for t in range(3*len(chain)+1):
        apex=_clean_point(quaternion_product(anchor,(1,t,t*t,t*t*t)))
        if all(_rank(v+(apex,))>rank for (_,v),rank in zip(chain,ranks)):
            return [(c,(apex,)+v) for c,v in chain]
    raise ArithmeticError('Exact cone bound exceeded')


def _face(a,b,c):
    if a==b or b==c:
        return []
    if _safe((a,b,c)):
        return [(1,(a,b,c))]
    return _cone(a,_edge(b,c)+[(-w,v) for w,v in _edge(a,c)]+_edge(a,b))


def exact_chain(vertices):
    """Same E/F/T convention as the previous geometric evaluator, without floats."""
    vertices=tuple(_clean_point(v) for v in vertices)
    a,b,c,d=vertices
    if any(x==y for x,y in zip(vertices,vertices[1:])):
        return []
    if _safe(vertices):
        return [(1,vertices)]
    return _cone(a,_face(b,c,d)+[(-w,v) for w,v in _face(a,c,d)]
                 +_face(a,b,d)+[(-w,v) for w,v in _face(a,b,c)])


def global_formula(name,g1,g2,g3,k=1,*,expand=True,reduce=False):
    """Exact original phase for every E-group triple, including exceptional ones.

    expand=False is a chain-inspection mode, explicitly NOT an expanded phase:
    it returns phase=None whenever a nonzero-volume term remains.
    """
    k=exact_level(k)
    members=set(exact_group(name))
    gs=tuple(unit_quaternion(g) for g in (g1,g2,g3))
    if any(g not in members for g in gs):
        raise ValueError('An input is not a member of '+name)
    a,b,c=gs
    p2=_clean_point(quaternion_product(a,b))
    vertices=(IDENTITY,a,p2,_clean_point(quaternion_product(p2,c)))
    chain=exact_chain(vertices)
    terms=[]
    volume=s.S.Zero
    pending=False
    for coefficient,points in chain:
        rank=_rank(points)
        record={'coefficient':coefficient,'vertices':points,'rank':rank,'safe':_safe(points)}
        if rank==4:
            unit=tuple(_clean_point(tuple(x/s.sqrt(sum(t*t for t in p)) for x in p)) for p in points)
            record['unit_vertices']=unit
            if expand:
                # Orthants have an immediate proved conjugate cancellation.
                orthant=all(sum(x*y for x,y in zip(unit[i],unit[j]))==0
                            for i,j in combinations(range(4),2))
                formula=tetrahedron_formula(unit,k=k,reduce=reduce or orthant)
                record['formula']=formula
                volume+=coefficient*formula['orientation']*formula['raw_volume']
            else:
                pending=True
        else:
            record['oriented_volume']=s.S.Zero
        terms.append(record)
    return {'group':name,'inputs':gs,'cumulative_vertices':vertices,
            'branch':'empty' if not chain else 'direct' if _safe(vertices) else 'cone',
            'terms':tuple(terms),'oriented_volume_mod':None if pending else volume,
            'phase':None if pending else s.exp(-s.I*k*volume/s.pi,evaluate=(volume==0)),
            'expanded':not pending,'representative':'original E/F/T geometric phase'}


def irrational_volume_certificate(name):
    """Exact spatial-triangle joins disproving all-rational 2O/2I volume tables."""
    e0,e1,e2=(IDENTITY,(0,1,0,0),(0,0,1,0))
    if name=='2O':
        vertices=(e0,e1,(0,1/s.sqrt(2),1/s.sqrt(2),0),(0,0,1/s.sqrt(2),1/s.sqrt(2)))
        root_order=8
    elif name=='2I':
        phi=(1+s.sqrt(5))/2
        vertices=(e0,e1,e2,(0,s.Rational(1,2),phi/2,(phi-1)/2))
        root_order=4
    else:
        raise ValueError('The irrational certificates are for 2O and 2I')
    A,B,C=(s.Matrix(v[1:]) for v in vertices[1:])
    det=s.simplify(s.Matrix.hstack(A,B,C).det())
    tangent=s.radsimp(det/(1+A.dot(B)+B.dot(C)+C.dot(A)))
    theta=s.atan(tangent)
    return {'group':name,'vertices':vertices,'determinant':det,
            'tan_half_area':tangent,'volume':s.pi*theta/2,
            'phase_level_one':s.exp(-s.I*theta/2),
            'root_of_unity_order_bound':root_order,
            'proof':'spherical triangle excess and join; cyclotomic field exclusion',
            'irrational_volume_over_pi2':True}
