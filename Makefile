PYTHON ?= python
TECTONIC ?= tectonic

.PHONY: install test baseline validate edge adecheck eexact symboliccheck symbolicreport polarcheck polarreport polarexistence polaraudit verify pdf adepdf completepdf all clean

install:
	$(PYTHON) -m pip install -e '.[validation]'

test:
	$(PYTHON) -m unittest discover -s tests -v

baseline:
	$(PYTHON) -m scripts.su2_audit_checks

validate:
	$(PYTHON) -m scripts.validate_su2

edge:
	$(PYTHON) -m scripts.su2_edge_crosscheck

adecheck:
	$(PYTHON) -m scripts.check_ade_phase

eexact:
	$(PYTHON) -m scripts.check_e_algebraic

# New exact route: no floating computation or change of representative.
symboliccheck:
	$(PYTHON) -m unittest tests.test_su2_symbolic tests.test_2t_geometry tests.test_exact_ade tests.test_exact_certificates tests.test_symbolic_cli -v

# CAS timeouts are recorded as unfinished attempts, not impossibility claims.
symbolicreport:
	$(PYTHON) -m scripts.check_symbolic_2t --symbolic-trials all --symbolic-budget-seconds 4 --jobs 4 --output results/SU2_symbolic_ADE_exact.json

# Second round: polar dual, exact only (no floating point).
polarcheck:
	$(PYTHON) -m unittest tests.test_polar_dual tests.test_polar_global tests.test_dic_exact -v

polarreport:
	PYTHONPATH=src $(PYTHON) -m scripts.check_polar_dual
	PYTHONPATH=src $(PYTHON) -m scripts.render_exact_catalogue

# Exhaustive per-input existence check of the elementary decomposition (2I takes hours).
polarexistence:
	for g in 2T 2O 2I; do PYTHONPATH=src $(PYTHON) -m scripts.check_polar_global_existence --group $$g --output build/polar_existence_$$g.json; done

# Numerical audit only; never part of a proof.
polaraudit:
	PYTHONPATH=src $(PYTHON) -m scripts.audit_polar_numeric

# Keep these sequential: the edge crosscheck consumes the main validation report.
verify:
	$(MAKE) test
	$(MAKE) baseline
	$(MAKE) validate
	$(MAKE) edge
	$(MAKE) adecheck
	$(MAKE) eexact

pdf:
	$(PYTHON) scripts/build_paper.py --engine $(TECTONIC)

adepdf:
	mkdir -p build/ade_review
	$(TECTONIC) --keep-logs --outdir build/ade_review paper/su2_ade_review.tex
	cp build/ade_review/su2_ade_review.pdf paper/su2_ade_review.pdf

completepdf:
	$(PYTHON) scripts/build_complete_report.py --engine $(TECTONIC)

all:
	$(MAKE) verify
	$(MAKE) completepdf

clean:
	$(PYTHON) -c "import shutil; shutil.rmtree('build', ignore_errors=True)"
