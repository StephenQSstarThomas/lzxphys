"""Deterministic, independent integral and structural validation; writes JSON."""
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss

from scripts.su2_audit_checks import PAULI, geometry
from su2_omega import (IDENTITY, omega_quaternions, oriented_volume,
                       quaternion, quaternion_multiply, tetrahedron_chain)

RESULTS = Path(__file__).resolve().parents[1] / 'results'
REPORT_PATH = RESULTS / 'SU2_validation.json'


def matrix(q):
    v = np.array([float(x) for x in q])
    v /= np.linalg.norm(v)
    return v[0]*np.eye(2) + 1j*np.einsum('a,aij->ij', v[1:], PAULI)


def reference_integral(gs, order):
    """Original triple integral via independently multiplied 2x2 matrices.

    Duffy map followed by tensor Gauss-Legendre; O(order^2) memory.
    """
    determinant, gram = geometry([matrix(g) for g in gs])
    nodes, weights = leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    b, c = np.meshgrid(nodes, nodes, indexing='ij')
    w2 = weights[:, None]*weights[None, :]
    total = 0.
    for a, weight in zip(nodes, weights):
        t = [(1-a)*(1-b)*(1-c), np.full_like(b, a), (1-a)*b, (1-a)*(1-b)*c]
        norm = sum(gram[i, i]*t[i]**2 for i in range(4))
        norm += sum(2*gram[i, j]*t[i]*t[j] for i in range(4) for j in range(i+1, 4))
        total += weight*np.sum(w2*(1-a)**2*(1-b)/norm**2)
    return float(determinant*total)


def vertices(gs):
    result = [IDENTITY]
    for g in gs:
        result.append(quaternion_multiply(result[-1], g))
    return result


def regular_raw_volume(theta):
    """Independent symmetry reduction, used only to expose the Arg branch jump."""
    a = mp.exp(1j*theta)
    r0 = 3*a**2 + 4*a**3 + a**6
    r1 = 12*mp.sin(theta)**2
    z = -2*r0/(r1+mp.sqrt(r1**2-4*abs(r0)**2))
    ell = (mp.polylog(2,z)+3*mp.polylog(2,z/a**4)
           -4*mp.polylog(2,-z/a**3)-3*theta**2)/2
    return -mp.re(ell)+mp.pi*(mp.arg(-mp.conj(r0))+3*theta)-3*mp.pi**2/2


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    seed = 20260913
    rng = np.random.default_rng(seed)
    report = {'seed': seed, 'mpmath': mp.__version__, 'numpy': np.__version__,
              'formula_dps': 55, 'integral_dtype': 'float64', 'samples': [],
              'integral_tolerance': 2e-8,
              'note': 'Quadrature convergence estimates and discrepancies are empirical, not certified bounds.'}
    triples = []
    for index in range(24):
        raw = rng.integers(-12, 13, size=(3, 4)).tolist()
        gs = [quaternion(q) for q in raw]
        tetra = vertices(gs)
        value = float(oriented_volume(tetra, dps=55))
        history = []
        for order in (48, 96, 192, 384, 768):
            reference = reference_integral(gs, order)
            history.append({'order': order, 'volume': reference})
            if len(history)>1 and abs(history[-1]['volume']-history[-2]['volume']) < 2e-9:
                break
        error = abs(value-reference)
        report['samples'].append({'index': index, 'raw_quaternions': raw,
                                  'closed_volume': value, 'reference_history': history,
                                  'absolute_error': error,
                                  'converged': len(history)>1 and abs(history[-1]['volume']-history[-2]['volume'])<2e-9})
        triples.append(gs)
        print(f'integral {index+1}/24: order={order}, error={error:.3g}', flush=True)
        REPORT_PATH.write_text(json.dumps(report, indent=2)+'\n')

    with mp.workdps(60):
        conjugation_errors, level_errors, norm_errors, pentagon_errors = [], [], [], []
        for gs in triples[:8]:
            h = quaternion(rng.integers(-8, 9, size=4).tolist())
            inv = (h[0], -h[1], -h[2], -h[3])
            transformed = [quaternion_multiply(quaternion_multiply(h, g), inv) for g in gs]
            w1 = omega_quaternions(*gs, k=1, dps=55)
            w2 = omega_quaternions(*gs, k=2, dps=55)
            w3 = omega_quaternions(*gs, k=3, dps=55)
            conjugation_errors.append(float(abs(w1-omega_quaternions(*transformed, dps=55))))
            level_errors.append(float(abs(w1*w2-w3)))
            norm_errors.append(float(abs(abs(w1)-1)))
        for _ in range(8):
            a, b, c, d = [quaternion(rng.integers(-8, 9, size=4).tolist()) for _ in range(4)]
            mul = quaternion_multiply
            w = lambda x, y, z: omega_quaternions(x, y, z, dps=55)
            lhs = w(b,c,d)*w(a,mul(b,c),d)*w(a,b,c)
            rhs = w(mul(a,b),c,d)*w(a,b,mul(c,d))
            pentagon_errors.append(float(abs(lhs-rhs)))
        report['structure'] = {
            'conjugation_count': 8, 'conjugation_max_error': max(conjugation_errors),
            'level_addition_count': 8, 'level_addition_max_error': max(level_errors),
            'unit_modulus_max_error': max(norm_errors),
            'generic_pentagon_count': 8, 'generic_pentagon_max_error': max(pentagon_errors)}
        report['cyclic_subgroups'] = []
        for n, h in [(4, (0,1,0,0)), (8, (1,1,0,0))]:
            power, product = IDENTITY, mp.mpc(1)
            for _ in range(n):
                product *= omega_quaternions(h, power, h, dps=55)
                power = quaternion_multiply(power, h)
            expected = mp.exp(-2j*mp.pi/n)
            report['cyclic_subgroups'].append({'order': n, 'product': mp.nstr(product, 35),
                                                'expected': mp.nstr(expected,35), 'error': float(abs(product-expected))})
        report['near_degenerate'] = []
        for exponent in (6,12,24):
            eps = Fraction(1,10**exponent)
            for name, gs in [
                ('near_identity', [(1,eps,0,0), (1,0,eps,0), (1,0,0,eps)]),
                ('near_minus_identity_positive', [(-1,eps,0,0), (1,2,3,4), (2,-1,1,3)]),
                ('near_minus_identity_negative', [(-1,-eps,0,0), (1,2,3,4), (2,-1,1,3)]),
                ('near_degenerate_face', [(1,1,0,0), (1,2,eps,0), (1,0,0,1)])]:
                low = omega_quaternions(*gs, dps=50)
                high = omega_quaternions(*gs, dps=80)
                entry = {'case': name, 'epsilon': f'1e-{exponent}',
                         'phase': mp.nstr(high, 35), 'precision_difference': float(abs(low-high))}
                if name == 'near_identity':
                    vol = oriented_volume(vertices(gs), dps=80)
                    e = mp.mpf(1)/10**exponent
                    entry['volume_over_epsilon_cubed_sixth'] = mp.nstr(vol/(e**3/6), 35)
                report['near_degenerate'].append(entry)
        center_tetra = vertices([(-1,0,0,0)]*3)
        center_chain = tetrahedron_chain(center_tetra)
        center_volume = mp.fsum(s*oriented_volume(t, dps=60) for s,t in center_chain)
        report['center'] = {'tetrahedron_count': len(center_chain), 'volume': mp.nstr(center_volume,40),
                            'error_from_pi_squared': float(abs(center_volume-mp.pi**2))}

        theta_cut = mp.findroot(lambda t: mp.im(3*mp.exp(2j*t)+4*mp.exp(3j*t)+mp.exp(6j*t)), (2.2,2.4))
        regular = []
        for shift in (-mp.mpf('1e-12'), mp.mpf('1e-12')):
            theta = theta_cut+shift
            c = mp.cos(theta)/(1-2*mp.cos(theta))
            gram = mp.matrix([[1 if i==j else c for j in range(4)] for i in range(4)])
            q = mp.cholesky(gram).T
            vol = oriented_volume([[q[i,j] for i in range(4)] for j in range(4)], dps=55)
            raw = regular_raw_volume(theta)
            regular.append((raw, vol, mp.exp(1j*raw/mp.pi)))
        raw_jump = abs(regular[1][0]-regular[0][0])
        formula_match = max(abs((raw % (2*mp.pi**2))-vol) for raw,vol,_ in regular)
        report['branch_jump'] = {
            'theta_cut': mp.nstr(theta_cut,40), 'theta_offset': '1e-12',
            'raw_volume_jump': mp.nstr(raw_jump,35),
            'jump_difference_from_two_pi_squared': float(abs(raw_jump-2*mp.pi**2)),
            'phase_difference': float(abs(regular[1][2]-regular[0][2])),
            'specialized_general_formula_error': float(formula_match),
            'meaning': 'Arg changes the raw expression by 2*pi^2; the phase remains continuous.'}
        positive = next(x for x in report['near_degenerate'] if x['case']=='near_minus_identity_positive' and x['epsilon']=='1e-24')
        negative = next(x for x in report['near_degenerate'] if x['case']=='near_minus_identity_negative' and x['epsilon']=='1e-24')
        eps = Fraction(1,10**24)
        side1 = omega_quaternions((-1,eps,0,0),(1,2,3,4),(2,-1,1,3),dps=60)
        side2 = omega_quaternions((-1,-eps,0,0),(1,2,3,4),(2,-1,1,3),dps=60)
        exact = omega_quaternions((-1,0,0,0),(1,2,3,4),(2,-1,1,3),dps=60)
        report['geometric_discontinuity'] = {
            'positive_side': positive['phase'], 'negative_side': negative['phase'],
            'exact_assigned_phase': mp.nstr(exact,35), 'side_difference': float(abs(side1-side2)),
            'meaning': 'Different limiting face fillings can change the representative phase, beyond a 2*pi^2 lift.'}

    report['max_integral_error'] = max(x['absolute_error'] for x in report['samples'])
    report['integral_unconverged_indices'] = [x['index'] for x in report['samples'] if not x['converged']]
    report['positive_orientation_count'] = sum(x['closed_volume']>0 for x in report['samples'])
    report['negative_orientation_count'] = sum(x['closed_volume']<0 for x in report['samples'])
    checks = {
        'independent_integrals': report['max_integral_error']<2e-8 and not report['integral_unconverged_indices'],
        'both_orientations': report['positive_orientation_count']>0 and report['negative_orientation_count']>0,
        'conjugation': report['structure']['conjugation_max_error']<1e-48,
        'level_addition': report['structure']['level_addition_max_error']<1e-48,
        'unit_modulus': report['structure']['unit_modulus_max_error']<1e-48,
        'generic_pentagon': report['structure']['generic_pentagon_max_error']<1e-48,
        'cyclic_subgroups': max(x['error'] for x in report['cyclic_subgroups'])<1e-48,
        'center': report['center']['error_from_pi_squared']<1e-48,
        'near_degenerate_precision': max(x['precision_difference'] for x in report['near_degenerate'])<1e-45,
        'arg_branch_jump': report['branch_jump']['jump_difference_from_two_pi_squared']<1e-9
                           and report['branch_jump']['phase_difference']<1e-9
                           and report['branch_jump']['specialized_general_formula_error']<1e-48,
        'distinct_face_limits': report['geometric_discontinuity']['side_difference']>.1}
    report['checks'] = checks
    report['passed'] = all(checks.values())
    REPORT_PATH.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('samples','near_degenerate')}, indent=2))
    if not report['passed']:
        raise SystemExit('Validation has unresolved checks; see SU2_validation.json')


if __name__ == '__main__':
    main()
