PYTHON ?= python3

.PHONY: verify

verify:
	$(PYTHON) papers/01-directional-stability/analysis/verify_release.py
	$(PYTHON) papers/02-present-state-resolution-and-lineage-history/analysis/verify_release.py
