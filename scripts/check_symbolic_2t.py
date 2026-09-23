"""Exact Li2/E-group report; the historical name includes the requested 2T entry.

Symbolic trial time limits bound CAS work, NOT mathematical validity. A timeout
means only that this attempt did not finish; the explicit master formula remains.
No approximate value is calculated or used to assign a volume or phase.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
from pathlib import Path
import signal
import sys

import sympy as s
from su2_symbolic import formula_template, tetrahedron_formula, reduce_real_dilog
from su2_2t_geometry import classify_2t
from su2_exact_ade import classify_group, exact_group, irrational_volume_certificate


def json_exact(value):
    if isinstance(value,dict):
        return {str(k):json_exact(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [json_exact(v) for v in value]
    if isinstance(value,Fraction):
        return {'numerator':value.numerator,'denominator':value.denominator}
    if isinstance(value,s.MatrixBase):
        return [[str(x) for x in row] for row in value.tolist()]
    if isinstance(value,s.Basic):
        if value.has(s.Float):
            raise ValueError('An approximate number entered the exact report')
        return str(value)
    return value


class SymbolicBudgetExceeded(Exception):
    pass


def _alarm(*_):
    raise SymbolicBudgetExceeded()


def symbolic_trial(vertices,seconds):
    old=signal.signal(signal.SIGALRM,_alarm)
    signal.setitimer(signal.ITIMER_REAL,seconds)
    stage='formula_construction'
    r=None
    try:
        r=tetrahedron_formula(vertices,reduce=False)
        stage='li2_reduction'
        old_reduction=r['dilog_reduction']
        reduction=reduce_real_dilog(old_reduction['remaining'])
        raw=r['raw_volume']+(reduction['expression']-old_reduction['expression'])/2
        return {'status':'li2_eliminated' if not reduction['remaining'] else 'li2_remaining',
                'identities':reduction['identities'],'remaining':reduction['remaining'],
                'raw_volume_mod_2pi2':raw,'Z':r['Z']}
    except SymbolicBudgetExceeded:
        return {'status':'budget_exhausted','stage':stage,
                'meaning':'CAS trial unfinished; not a proof of impossibility',
                **({'Z':r['Z']} if r is not None else {})}
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        signal.signal(signal.SIGALRM,old)


def _trial_job(job):
    name,index,row,seconds=job
    points=exact_group(name)
    trial=symbolic_trial(tuple(points[j] for j in row['representative_indices']),seconds)
    return {'type_index':index,'gram_key':row['gram_key'],**trial}


def make_report(trials='all',seconds=5,jobs=1):
    if trials not in ('all','2t','none') or seconds<=0 or jobs<1:
        raise ValueError('Invalid exact-report trial mode, budget, or worker count')
    catalogues={name:classify_group(name) for name in ('2T','2O','2I')}
    for name,catalogue in catalogues.items():
        catalogue['vertices']=exact_group(name)
        d=catalogue['field_d']
        for row in catalogue['types']:
            a,b=row['determinant_pair_over_64']
            det=(s.Integer(a)+b*s.sqrt(d))/64
            row['formula_inputs']={
                **{'c'+str(j+1):(s.Integer(a)+b*s.sqrt(d))/16
                   for j,(a,b) in enumerate(row['cosines_in_formula_order'])},
                'detG':s.expand(det**2),'orientation':1,'k':'any integer'}
            row['formula_scope']='volume modulo 2*pi**2 of canonical Gram; ordered input supplies orientation separately'
    result={'method':'exact principal Li2 followed by exact geometric reasoning',
            'numeric_evaluation_used':False,'transfer_used':False,
            'formula_template':formula_template(),'catalogues':catalogues,
            'hurwitz_actual_volumes':classify_2t(),
            'irrational_geometric_certificates':{name:irrational_volume_certificate(name) for name in ('2O','2I')},
            'symbolic_trial_budget_seconds':seconds,'symbolic_workers':jobs,
            'symbolic_trial_order':'construct exact expression, then Li2 identities; logs stay explicit',
            'symbolic_trials':{}}
    names=('2T','2O','2I') if trials=='all' else ('2T',) if trials=='2t' else ()
    for name in names:
        records=[]
        tasks=[(name,index,row,seconds) for index,row in enumerate(catalogues[name]['types'])]
        def collect(iterator):
            for record in iterator:
                records.append(record)
                if len(records)%25==0:
                    print(name,len(records),'/',len(tasks),dict(Counter(r['status'] for r in records)),file=sys.stderr,flush=True)
        if jobs==1:
            collect(map(_trial_job,tasks))
        else:
            with ProcessPoolExecutor(max_workers=jobs) as pool:
                collect(pool.map(_trial_job,tasks))
        result['symbolic_trials'][name]=records
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--symbolic-trials',choices=('none','2t','all'),default='all')
    parser.add_argument('--symbolic-budget-seconds',type=int,default=5)
    parser.add_argument('--jobs',type=int,default=1,help='Independent CAS worker processes; default 1')
    parser.add_argument('--group',choices=('2T','2O','2I'))
    parser.add_argument('--type-index',type=int)
    args=parser.parse_args()
    if args.symbolic_budget_seconds<=0:
        parser.error('symbolic budget must be positive')
    if args.jobs<1:
        parser.error('jobs must be positive')
    if (args.group is None)!=(args.type_index is None):
        parser.error('--group and --type-index must be used together')
    if args.group:
        catalogue=classify_group(args.group)
        if not 0<=args.type_index<len(catalogue['types']):
            parser.error('type index outside catalogue')
        row=catalogue['types'][args.type_index]
        points=exact_group(args.group)
        report={'group':args.group,'type_index':args.type_index,
                'formula':tetrahedron_formula(tuple(points[j] for j in row['representative_indices']),reduce=False)}
    else:
        report=make_report(args.symbolic_trials,args.symbolic_budget_seconds,args.jobs)
    encoded=json_exact(report)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(encoded,ensure_ascii=False,indent=2)+'\n')
        print('Wrote',args.output)
    if 'catalogues' in report:
        print(json.dumps({name:{'types':len(c['types']),'nondegenerate_subsets':c['nondegenerate_subsets']}
                          for name,c in report['catalogues'].items()}))
        print(json.dumps({name:dict(Counter(r['status'] for r in rows)) for name,rows in report['symbolic_trials'].items()}))
    elif not args.output:
        print(json.dumps(encoded,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
