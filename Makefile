PYTHON ?= python
TECTONIC ?= tectonic

.PHONY: install test baseline validate edge verify pdf all clean

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

# Keep these sequential: the edge crosscheck consumes the main validation report.
verify:
	$(MAKE) test
	$(MAKE) baseline
	$(MAKE) validate
	$(MAKE) edge

pdf:
	$(PYTHON) scripts/build_paper.py --engine $(TECTONIC)

all:
	$(MAKE) verify
	$(MAKE) pdf

clean:
	$(PYTHON) -c "import shutil; shutil.rmtree('build', ignore_errors=True)"
