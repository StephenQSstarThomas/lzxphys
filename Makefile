PYTHON ?= python
TECTONIC ?= tectonic

.PHONY: install test baseline validate edge adecheck eexact verify pdf adepdf completepdf all clean

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
