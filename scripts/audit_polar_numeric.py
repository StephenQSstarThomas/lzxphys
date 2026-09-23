"""NUMERICAL AUDIT ONLY -- not part of any proof.

Cross-checks the exact polar-dual results against independent historical
implementations with mpmath, to catch implementation slips:
  1. all 659 type volumes vs the dihedral-angle Murakami formula (su2_omega);
  2. elementary global phases vs the historical mpmath E/F/T chain (2T all, 2O sample);
  3. 2I global phases vs the exact Li2 global chain of su2_exact_ade (small sample).
The derivations and certificates never use these numbers.
"""
import argparse
import itertools
import json
import random
import time
from multiprocessing import Pool
from pathlib import Path

import mpmath as mp
import sympy as s
from su2_exact_ade import _clean_point, _safe, exact_group, global_formula
from su2_omega import omega_quaternions_edge, oriented_volume
from su2_polar_dual import polar_catalogue
from su2_polar_global import global_polar_formula
from su2_symbolic import IDENTITY, quaternion_product

NAME = None


def catalogue_audit():
    mp.mp.dps = 40
    worst = mp.mpf(0)
    for name in ('2T', '2O', '2I'):
        group = exact_group(name)
        for row in polar_catalogue(name)['types']:
            ref = abs(oriented_volume([tuple(mp.mpf(s.N(x, 50)) for x in group[i])
                                       for i in row['representative_indices']], dps=35))
            worst = max(worst, abs(mp.mpf(str(s.N(row['volume'], 40)))-ref))
    return float(worst)


def _ray(q):
    if NAME == '2O' and any(x.has(s.sqrt(2)) for x in q):
        q = tuple(s.simplify(x*s.sqrt(2)) for x in q)
    return tuple(s.Rational(x) for x in q)


def _history(indices):
    group = exact_group(NAME)
    g = [group[i] for i in indices]
    r = global_polar_formula(NAME, *g)
    ours = complex(s.N(r['phase'], 30))
    ref = complex(omega_quaternions_edge(*[_ray(x) for x in g], k=1, dps=30))
    return r['branch'], abs(ours-ref)


def _init(name):
    global NAME
    NAME = name


def history_audit(name, sample, jobs):
    size = len(exact_group(name))
    triples = list(itertools.product(range(size), repeat=3))
    if sample:
        random.seed(5)
        triples = random.sample(triples, sample)
    worst, counts = {}, {}
    with Pool(jobs, initializer=_init, initargs=(name,)) as pool:
        for branch, diff in pool.imap_unordered(_history, triples, chunksize=32):
            counts[branch] = counts.get(branch, 0)+1
            worst[branch] = max(worst.get(branch, 0.0), diff)
    return {'triples': len(triples), 'branches': counts, 'max_phase_difference': worst}


def _li2(indices):
    group = exact_group('2I')
    g = [group[i] for i in indices]
    ours = complex(s.N(global_polar_formula('2I', *g)['phase'], 25))
    ref = complex(s.N(global_formula('2I', *g, expand=True, reduce=False)['phase'], 25))
    return abs(ours-ref)


def li2_audit(sample, jobs):
    group = exact_group('2I')
    random.seed(11)
    chosen = []
    while len(chosen) < sample:
        idx = tuple(random.randrange(len(group)) for _ in range(3))
        a, b, c = (group[i] for i in idx)
        p2 = _clean_point(quaternion_product(a, b))
        v = (IDENTITY, _clean_point(a), p2, _clean_point(quaternion_product(p2, c)))
        if not any(x == y for x, y in zip(v, v[1:])) and not _safe(v):
            chosen.append(idx)
    with Pool(jobs) as pool:
        diffs = pool.map(_li2, chosen, chunksize=1)
    return {'cone_triples': len(chosen), 'max_phase_difference': max(diffs)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/SU2_polar_numeric_audit.json'))
    parser.add_argument('--jobs', type=int, default=6)
    parser.add_argument('--o-sample', type=int, default=8000)
    parser.add_argument('--i-sample', type=int, default=12)
    args = parser.parse_args()
    start = time.time()
    report = {'purpose': 'numerical audit only; not used in any derivation or certificate',
              'catalogue_vs_murakami_max_volume_difference': catalogue_audit(),
              'global_vs_historical_mpmath': {'2T': history_audit('2T', 0, args.jobs),
                                              '2O': history_audit('2O', args.o_sample, args.jobs)},
              'global_2I_vs_exact_li2_chain': li2_audit(args.i_sample, args.jobs)}
    report['seconds'] = round(time.time()-start)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
