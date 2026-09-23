"""Auditable E-type examples and pentagons, with exact inputs and all five terms.

These are selected implementation checks, not exhaustive E-type verification.
The independent analytic examples use orthant subdivisions and a spherical join.
"""
import json
from pathlib import Path

import mpmath as mp

from scripts.ade_phase import (
    ade_group, phase_from_exact_quaternions, phase_details, _moment, _safe, _reduce_point,
)
from scripts.ade_subgroups import sympy_quaternion_multiply
import sympy as sp


CASES = {
    "2T": ((1, 5, 9, 13), (2, 7, 14, 20)),
    "2O": ((1, 5, 9, 17), (3, 11, 21, 31)),
    "2I": ((1, 5, 17, 41), (3, 13, 29, 67)),
}


def analytic_inputs(name):
    half = sp.Rational(1, 2)
    s2, s5 = sp.Symbol('s2'), sp.Symbol('s5')
    third = {'2T': (half, half, -half, -half),
             '2O': (s2/2, s2/2, 0, 0),
             '2I': ((1+s5)/4, (s5-1)/4, 0, -half)}[name]
    return ((0, 1, 0, 0), (0, 0, 0, 1), third)


def complex_record(value, digits=65):
    return {'real': mp.nstr(mp.re(value), digits), 'imag': mp.nstr(mp.im(value), digits)}


def complex_value(record):
    return mp.mpc(record['real'], record['imag'])


def check_group(name, index_sets):
    points = ade_group(name)
    errors, checks = [], []
    cache = {}
    def detail(*triple):
        key = tuple(_reduce_point(q) for q in triple)
        if key not in cache:
            cache[key] = phase_details(*key, dps=65)
        return cache[key]
    mul = sympy_quaternion_multiply
    with mp.workdps(85):
        # Add a case containing the explicit non-Q8 example for each E group.
        lookup = {_reduce_point(q): i for i, q in enumerate(points)}
        example_indices = tuple(lookup[_reduce_point(q)] for q in analytic_inputs(name))
        index_sets = tuple(index_sets) + (example_indices + (8,),)
        for a_i, b_i, c_i, d_i in index_sets:
            a, b, c, d = (points[i] for i in (a_i, b_i, c_i, d_i))
            terms = {'b,c,d': detail(b, c, d), 'a,bc,d': detail(a, mul(b, c), d),
                     'a,b,c': detail(a, b, c), 'ab,c,d': detail(mul(a, b), c, d),
                     'a,b,cd': detail(a, b, mul(c, d))}
            values = {label: complex_value(item['phase']) for label, item in terms.items()}
            lhs = values['b,c,d']*values['a,bc,d']*values['a,b,c']
            rhs = values['ab,c,d']*values['a,b,cd']
            errors.append(float(abs(lhs-rhs)))
            checks.append({'indices': [a_i, b_i, c_i, d_i],
                           'quadruple': [[str(x) for x in q] for q in (a, b, c, d)],
                           'terms': terms, 'lhs': complex_record(lhs), 'rhs': complex_record(rhs),
                           'residual': mp.nstr(abs(lhs-rhs), 15)})
        example = detail(*analytic_inputs(name))
        phi = (1+mp.sqrt(5))/2
        angle = {'2T': -mp.pi/32, '2O': -mp.pi/16,
                 '2I': -mp.atan(1/(4*phi+1))/2}[name]
        expected = mp.exp(1j*angle)
        observed = complex_value(example['phase'])
        angle_route = phase_from_exact_quaternions(*analytic_inputs(name), dps=65, method='angle')
        analytic_error = abs(observed-expected)
        crosscheck_error = abs(observed-angle_route)
        center = (-1, 0, 0, 0)
        center_example = detail(center, center, center)
        center_error = abs(complex_value(center_example['phase'])+1)
        identity_example = detail((1, 0, 0, 0), *analytic_inputs(name)[1:])
        normalization_error = abs(complex_value(identity_example['phase'])-1)
        report = {'order': len(points), 'cases': len(checks),
                  'coverage': '3 selected quadruples; not exhaustive; k=1',
                  'max_pentagon_error': max(errors), 'pentagons': checks,
                  'analytic_example': example, 'analytic_expected': complex_record(expected),
                  'analytic_error': float(analytic_error),
                  'edge_angle_error': float(crosscheck_error),
                  'center_example': center_example, 'center_error': float(center_error),
                  'identity_example': identity_example, 'normalization_error': float(normalization_error)}
        report['passed'] = max(max(errors), analytic_error, crosscheck_error,
                               center_error, normalization_error) < mp.mpf('1e-55')
        return report


def main():
    s2 = sp.Symbol("s2")
    # Reviewer regression: this convex hull contains zero and must not take the
    # direct radial branch; the cone ray must stay in the quadratic field.
    assert not _safe(((1, 0, 0, 0), (0, 1, 0, 0), (-s2/2, -s2/2, 0, 0)))
    assert _moment(2) == (1, 2, 4, 8)
    result = {name: check_group(name, cases) for name, cases in CASES.items()}
    path = Path(__file__).resolve().parents[1] / "docs" / "review-0922" / "SU2_ADE_phase_validation.json"
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    for name, item in result.items():
        print(name, 'quadruples:', item['cases'], 'max pentagon:', item['max_pentagon_error'],
              'analytic:', item['analytic_error'], 'edge/angle:', item['edge_angle_error'],
              'passed:', item['passed'])
    print('Full inputs, five factors and chain terms:', path.relative_to(path.parents[2]))
    if not all(item['passed'] for item in result.values()):
        raise SystemExit('E-type validation failed; inspect the saved report')


if __name__ == "__main__":
    main()
