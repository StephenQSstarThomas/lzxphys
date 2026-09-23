"""Exact polar-dual catalogue, angle certificates and the all-n Dic table proof data.

Everything written here is decided by exact rational/Q(sqrt d) arithmetic or
formal algebra; no floating point is used.  See scripts/audit_polar_numeric.py
for the separate, clearly labelled numerical audit.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path

import sympy as s
from su2_dic_exact import (dic_pairing, formal_closure_certificate,
                           formal_normalisation_certificate, lens_pairing)
from su2_polar_dual import (BASIS_TAN2, FIELD, Surd, angle_system, basis_symbols,
                            polar_catalogue, weyl_chamber_masks)


def frac(x):
    return str(Fraction(x))


def surd_text(key, d):
    return str(Surd(key[0], key[1], d).sympy())


def group_report(name):
    d = FIELD[name]
    catalogue = polar_catalogue(name)
    system = angle_system(name)
    rows = []
    for index, row in enumerate(catalogue['types']):
        rows.append({'type_index': index, 'gram_key': list(row['gram_key']),
                     'anchored_multiplicity': row['anchored_multiplicity'],
                     'representative_indices': list(row['representative_indices']),
                     'chambers': row['chambers'], 'rational_part_over_pi2': frac(row['rational_part']),
                     'basis_coefficients': [frac(b) for b in row['basis_coefficients']],
                     'rational_volume': row['rational_volume'], 'volume': str(row['volume'])})
    certificate = [{'identity': desc, 'relation': {surd_text(k, d): frac(c) for k, c in rel.items()},
                    'rhs_over_pi': frac(r)} for desc, rel, r in system['certificate']]
    angles = {surd_text(k, d): {'rational_over_pi': frac(r),
                                'basis': {surd_text(b, d): frac(c) for b, c in coeff.items()}}
              for k, (r, coeff) in system['angles'].items()}
    return {'weyl_order': weyl_chamber_masks(name)[0],
            'basis_tan2': [frac(T) for T in BASIS_TAN2[name]],
            'basis': [str(t) for t in basis_symbols(name)],
            'irrational_dihedral_tan2_count': len(system['angles']),
            'identity_candidates_checked': system['identity_count'],
            'certificate': certificate, 'angle_table': angles,
            'types': rows, 'type_count': len(rows),
            'rational_volume_types': sum(r['rational_volume'] for r in rows)}


def dic_report(max_n, max_m):
    closure = formal_closure_certificate()
    normal = formal_normalisation_certificate()
    return {'closure_residuals': {''.join(map(str, k)): str(v) for k, v in closure.items()},
            'normalisation_values': {f'slot{k[0]}_{k[1][0]}{k[1][1]}': str(v) for k, v in normal.items()},
            'lens_pairings': {str(M): frac(lens_pairing(M)) for M in range(1, max_m+1)},
            'dic_fundamental_pairings': {str(n): frac(dic_pairing(n)) for n in range(2, max_n+1)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/SU2_polar_dual_exact.json'))
    parser.add_argument('--max-n', type=int, default=16)
    parser.add_argument('--max-m', type=int, default=60)
    args = parser.parse_args()
    report = {'method': 'polar dual + Schlaefli + Coxeter chambers; certified angle identities',
              'numeric_evaluation_used': False, 'representative_changed': False,
              'groups': {name: group_report(name) for name in ('2T', '2O', '2I')},
              'dic_table': dic_report(args.max_n, args.max_m)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1, ensure_ascii=False)+'\n')
    summary = {n: (g['type_count'], g['rational_volume_types'], g['weyl_order'], len(g['certificate']))
               for n, g in report['groups'].items()}
    print('Wrote', args.output, summary)


if __name__ == '__main__':
    main()
