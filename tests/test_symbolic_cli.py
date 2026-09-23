import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import sympy as s


class SymbolicCLITests(unittest.TestCase):
    def test_printable_complete_parameter_tables(self):
        self.assertIsNotNone(importlib.util.find_spec('scripts.render_exact_catalogue'))
        from scripts.render_exact_catalogue import render_catalogue
        text=render_catalogue()
        self.assertEqual(text.count('\\texttt{'),84+563)
        self.assertIn('\\sqrt2',text)
        self.assertIn('\\sqrt5',text)
        source=Path(__file__).resolve().parents[1]/'paper'/'su2_exact_catalogue.tex'
        self.assertEqual(source.read_text(),text)

    def test_li2_reduction_precedes_log_simplification(self):
        from scripts import check_symbolic_2t as cli
        vertices=tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
        original=cli.tetrahedron_formula
        with patch.object(cli,'tetrahedron_formula',wraps=original) as construct:
            result=cli.symbolic_trial(vertices,5)
        self.assertIs(construct.call_args.kwargs.get('reduce'),False)
        self.assertEqual(result['status'],'li2_eliminated')
        self.assertEqual(s.simplify(result['raw_volume_mod_2pi2']-s.pi**2/8),0)

    def test_timeout_records_actual_stage_and_keeps_arguments(self):
        from scripts import check_symbolic_2t as cli
        vertices=tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
        with patch.object(cli,'reduce_real_dilog',side_effect=cli.SymbolicBudgetExceeded):
            result=cli.symbolic_trial(vertices,5)
        self.assertEqual(result['status'],'budget_exhausted')
        self.assertEqual(result['stage'],'li2_reduction')
        self.assertEqual(len(result['Z']),8)

    def test_parallel_trials_keep_all_type_indices(self):
        from scripts import check_symbolic_2t as cli
        report=cli.make_report(trials='2t',seconds=1,jobs=2)
        rows=report['symbolic_trials']['2T']
        self.assertEqual([row['type_index'] for row in rows],list(range(12)))
        self.assertTrue(all(row['status'] in ('li2_eliminated','li2_remaining','budget_exhausted') for row in rows))

    def test_end_to_end_exact_catalogue_export(self):
        self.assertIsNotNone(importlib.util.find_spec('scripts.check_symbolic_2t'),
                             'Exact symbolic ADE CLI is missing')
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'report.json'
            run=subprocess.run([sys.executable,'-m','scripts.check_symbolic_2t',
                                '--symbolic-trials','none','--output',str(output)],
                               capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            report=json.loads(output.read_text())
            self.assertEqual([len(report['catalogues'][n]['types']) for n in ('2T','2O','2I')],[12,84,563])
            self.assertEqual(len(report['hurwitz_actual_volumes']['types']),12)
            self.assertEqual(report['numeric_evaluation_used'],False)
            self.assertIn('formula_template',report)
            self.assertEqual(report['symbolic_trials'],{})


if __name__=='__main__':
    unittest.main()
