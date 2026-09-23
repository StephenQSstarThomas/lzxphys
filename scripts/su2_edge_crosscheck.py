"""Independent edge-length formula (Murakami v4 Theorem 1.2).

No inverse Gram or production dihedral-angle helpers enter this calculation.
The partial derivatives in Theorem 1.2 are evaluated by explicit principal logs,
holding z fixed. This is a verification implementation for nondegenerate inputs.
"""
import json
from pathlib import Path

import mpmath as mp

from su2_omega import IDENTITY, oriented_volume, oriented_volume_from_edges, quaternion_multiply

RESULTS = Path(__file__).resolve().parents[1] / 'results'


def edge_volume(vertices, dps=85):
    with mp.workdps(dps):
        q = [[mp.mpf(v.numerator)/v.denominator if hasattr(v, 'denominator')
              else mp.mpf(v) for v in point] for point in vertices]
        q = [[v/mp.sqrt(mp.fsum(x*x for x in point)) for v in point] for point in q]
        edges = ((2,3), (1,3), (0,3), (0,1), (0,2), (1,2))
        lengths = [mp.acos(mp.fsum(x*y for x,y in zip(q[i],q[j]))) for i,j in edges]
        theta = [mp.pi-lengths[(j+3)%6] for j in range(6)]
        a = [mp.exp(1j*t) for t in theta]
        r0 = (a[0]*a[3]+a[1]*a[4]+a[2]*a[5]
              +a[0]*a[1]*a[5]+a[0]*a[2]*a[4]+a[1]*a[2]*a[3]
              +a[3]*a[4]*a[5]+mp.fprod(a))
        r1 = 4*sum(mp.sin(theta[j])*mp.sin(theta[j+3]) for j in range(3))
        discriminant = r1*r1-4*abs(r0)**2
        if discriminant <= 0:
            raise ValueError('Nondegenerate spherical tetrahedron required')
        z = -2*r0/(r1+mp.sqrt(discriminant))
        pos_sets = ((0,1,3,4), (0,2,3,5), (1,2,4,5))
        neg_sets = ((0,1,2), (0,4,5), (1,3,5), (2,3,4))
        positive = [z/mp.fprod(a[j] for j in s) for s in pos_sets]
        negative = [-z/mp.fprod(a[j] for j in s) for s in neg_sets]
        le = (mp.polylog(2,z)+sum(mp.polylog(2,v) for v in positive)
              -sum(mp.polylog(2,v) for v in negative)
              -sum(theta[j]*theta[j+3] for j in range(3)))/2
        derivatives = []
        for j in range(6):
            varied = (j+3)%6
            derivative = (sum(mp.im(mp.log(1-v)) for s,v in zip(pos_sets,positive) if varied in s)
                          -sum(mp.im(mp.log(1-v)) for s,v in zip(neg_sets,negative) if varied in s)
                          +theta[j])/2
            derivatives.append(derivative)
        raw = (mp.re(le)-mp.pi*mp.arg(-mp.conj(r0))
               -sum(length*derivative for length,derivative in zip(lengths,derivatives))
               -mp.pi**2/2)
        volume = raw % (2*mp.pi**2)
        if not 0 < volume < mp.pi**2:
            raise ArithmeticError('Edge formula volume outside geometric range')
        sign = mp.sign(mp.det(mp.matrix(q)))
        return sign*volume


def main():
    samples = json.loads((RESULTS / 'SU2_validation.json').read_text())['samples']
    result = {'route': 'Murakami v4 Theorem 1.2, explicit fixed-z log derivatives',
              'edge_dps': 85, 'angle_dps': 70, 'indices': list(range(len(samples))), 'samples': []}
    for index in result['indices']:
        vertices = [IDENTITY]
        for g in samples[index]['raw_quaternions']:
            vertices.append(quaternion_multiply(vertices[-1],g))
        with mp.workdps(85):
            edge = edge_volume(vertices)
            angle = oriented_volume(vertices,dps=70)
            production = oriented_volume_from_edges(vertices, dps=70)
            error = abs(edge-angle)
            production_error = abs(production-angle)
            result['samples'].append({'index':index, 'edge_volume':mp.nstr(edge,65),
                                      'angle_volume':mp.nstr(angle,65), 'error':float(error),
                                      'production_edge_volume':mp.nstr(production,65),
                                      'production_error':float(production_error)})
            print(index, 'independent edge/angle:', mp.nstr(error,8),
                  'production edge/angle:', mp.nstr(production_error,8))
    result['maximum_error'] = max(x['error'] for x in result['samples'])
    result['maximum_production_error'] = max(x['production_error'] for x in result['samples'])
    result['passed'] = max(result['maximum_error'], result['maximum_production_error']) < 1e-65
    (RESULTS / 'SU2_edge_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:
        raise SystemExit('Independent edge formula disagrees; inspect the saved report')


if __name__ == '__main__':
    main()
