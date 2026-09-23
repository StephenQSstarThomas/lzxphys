"""Exhaustive exact check that the elementary decomposition exists for every input.

Every ordered triple (g1, g2, g3) is classified exactly: nondegenerate cumulative
vertices (det != 0, decided in Q(sqrt d)) are the direct branch; degenerate ones are
classified by the ORIGINAL safety test of su2_exact_ade; for every cone input the
certified face replacements and an admissible group apex of su2_polar_global must
exist (it raises ArithmeticError otherwise).  Direct inputs need no check: every
nondegenerate group tetrahedron is isometric to a catalogue representative, so its
dihedral angles lie in the certified table.  No floating point is used.
"""
import argparse
import json
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

from su2_exact_ade import _det3, _safe, exact_group, integer_pair_group
from su2_polar_dual import FIELD, _index, group_vectors, quaternion_product
from su2_polar_global import cone_structure

NAME = None
TABLE = None
PAIRS = None


def _init(name):
    global NAME, TABLE, PAIRS
    NAME = name
    PAIRS = integer_pair_group(name)
    points = group_vectors(name)
    index = _index(name)
    TABLE = [[index[tuple(x.key() for x in quaternion_product(p, q))] for q in points] for p in points]


def _row(i1):
    points, sym = group_vectors(NAME), exact_group(NAME)
    one = next(i for i, p in enumerate(points) if p[0] == 1)
    branches, kinds, failures = Counter(), Counter(), []
    for i2 in range(len(points)):
        j2 = TABLE[i1][i2]
        for i3 in range(len(points)):
            j3 = TABLE[j2][i3]
            idx = (one, i1, j2, j3)
            if any(a == b for a, b in zip(idx, idx[1:])):
                branches['empty'] += 1
                continue
            # det(e0, p1, p2, p3) = 3x3 det of the spatial parts; exact integer pairs
            if _det3(PAIRS[i1][1:], PAIRS[j2][1:], PAIRS[j3][1:], FIELD[NAME]) != (0, 0):
                branches['direct'] += 1
                continue
            vertices = tuple(sym[i] for i in idx)
            if _safe(vertices):
                branches['direct'] += 1
                continue
            branches['cone'] += 1
            try:
                _, corrections, _, _ = cone_structure(NAME, vertices)
                kinds['+'.join(sorted(c['type'] for c in corrections)) or 'none'] += 1
            except ArithmeticError as error:
                failures.append({'indices': (i1, i2, i3), 'error': repr(error)[:200]})
    return branches, kinds, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group', required=True, choices=('2T', '2O', '2I'))
    parser.add_argument('--jobs', type=int, default=6)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    size = len(exact_group(args.group))
    branches, kinds, failures = Counter(), Counter(), []
    with Pool(args.jobs, initializer=_init, initargs=(args.group,)) as pool:
        for b, k, f in pool.imap_unordered(_row, range(size)):
            branches.update(b)
            kinds.update(k)
            failures.extend(f)
    report = {'group': args.group, 'ordered_triples': size**3, 'branches': dict(branches),
              'cone_correction_patterns': dict(sorted(kinds.items())), 'failure_count': len(failures),
              'failures': failures[:50], 'numeric_evaluation_used': False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'failures'}))


if __name__ == '__main__':
    main()
