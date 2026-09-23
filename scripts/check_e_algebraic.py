"""Reproduce integer E pentagons and the explicit geometry/count gauge change.

The default checks 16 fixed quadruples per E group and one gauge-conversion
triple per group, with complete summands. This is NOT exhaustive E testing.
Use --skip-gauge for a faster exact-arithmetic-only check.
"""
import argparse
import json
from pathlib import Path
from random import Random

import mpmath as mp

from scripts.ade_count import count_details, algebraic_phase, gauge_details
from scripts.ade_phase import ade_group, phase_from_exact_quaternions, _reduce_point
from scripts.ade_subgroups import sympy_quaternion_multiply as mul
from scripts.check_ade_phase import analytic_inputs


def _complex(z, digits):
    return {'real': mp.nstr(mp.re(z), digits), 'imag': mp.nstr(mp.im(z), digits)}


def check_group(name, include_gauge=True, dps=45):
    points = tuple(map(_reduce_point, ade_group(name)))
    lookup = {q: i for i, q in enumerate(points)}
    order = len(points)
    example = analytic_inputs(name)
    examples = tuple(lookup[_reduce_point(q)] for q in example)
    cases = [(1, 1, 1, 1), examples+(order-1,), (0, 2, 3, order-1),
             (2, 2, 2, 2)]
    rng = Random(9222026+order)
    while len(cases) < 16:
        case = tuple(rng.randrange(order) for _ in range(4))
        if case not in cases:
            cases.append(case)
    pentagons = []
    for indices in cases:
        a, b, c, d = (points[i] for i in indices)
        triples = ((b, c, d), (mul(a, b), c, d), (a, mul(b, c), d),
                   (a, b, mul(c, d)), (a, b, c))
        factors = [count_details(name, *triple) for triple in triples]
        counts = [factor['count'] for factor in factors]
        delta = sum(sign*count for sign, count in zip((1, -1, 1, -1, 1), counts))
        pentagons.append({'indices': list(indices),
                         'inputs': [[str(x) for x in points[i]] for i in indices],
                         'term_order': ['b,c,d', 'ab,c,d', 'a,bc,d', 'a,b,cd', 'a,b,c'],
                         'signs': [1, -1, 1, -1, 1], 'counts': counts,
                         'delta_count': delta, 'delta_mod_order': delta % order,
                         'degree': delta//order if delta % order == 0 else None,
                         'factors': factors})
    x = (0, 1, 0, 0)
    powers = ((1, 0, 0, 0), x, (-1, 0, 0, 0), (0, -1, 0, 0))
    cycle = [count_details(name, x, h, x) for h in powers]
    cycle_exponent = sum(item['exponent_mod_order'] for item in cycle) % order
    report = {'order': order, 'coverage': '16 selected quadruples, exact integer arithmetic; NOT exhaustive',
              'pentagons': pentagons, 'analytic_example': count_details(name, *example),
              'C4_cycle': cycle, 'C4_exponent_mod_order': cycle_exponent,
              'C4_expected_exponent': order//4,
              'passed': all(item['delta_mod_order'] == 0 for item in pentagons)
                        and cycle_exponent == order//4}
    print(name, '16 exact pentagons:', report['passed'], flush=True)
    if include_gauge:
        a, b, c = example
        pairs = {'b,c': (b, c), 'ab,c': (mul(a, b), c),
                 'a,bc': (a, mul(b, c)), 'a,b': (a, b)}
        gauges = {}
        for label, pair in pairs.items():
            gauges[label] = gauge_details(name, *pair, dps=dps)
            print(name, 'gauge', label, 'summands:', len(gauges[label]['terms']), flush=True)
        with mp.workdps(dps+20):
            beta = {label: mp.mpc(item['beta']['real'], item['beta']['imag'])
                    for label, item in gauges.items()}
            delta_beta = beta['b,c']*beta['a,bc']/(beta['ab,c']*beta['a,b'])
            geometry = phase_from_exact_quaternions(a, b, c, dps=dps+10)
            algebraic = algebraic_phase(name, a, b, c, dps=dps+10)
            error = abs(algebraic-geometry*delta_beta)
            report['gauge_check'] = {'dps': dps, 'k': 1, 'gauges': gauges,
                                     'geometry': _complex(geometry, dps),
                                     'algebraic': _complex(algebraic, dps),
                                     'delta_beta': _complex(delta_beta, dps),
                                     'residual': mp.nstr(error, dps),
                                     'tolerance': '1e-'+str(dps-6),
                                     'passed': error < mp.mpf(10)**(-(dps-6))}
            report['passed'] &= report['gauge_check']['passed']
            print(name, 'gauge residual:', mp.nstr(error, 8), flush=True)
    else:
        report['gauge_check'] = {'status': 'not run (--skip-gauge)'}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-gauge', action='store_true')
    parser.add_argument('--dps', type=int, default=45)
    parser.add_argument('--output', type=Path, default=Path(
        'docs/review-0922/SU2_E_algebraic_validation.json'))
    args = parser.parse_args()
    if args.dps < 20:
        parser.error('--dps must be at least 20')
    report = {name: check_group(name, not args.skip_gauge, args.dps)
              for name in ('2T', '2O', '2I')}
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')
    print('Full report:', args.output)
    if not all(item['passed'] for item in report.values()):
        raise SystemExit('Exact E representative validation failed')


if __name__ == '__main__':
    main()
