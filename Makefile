PYTHON ?= python3

.PHONY: verify

verify:
	$(PYTHON) papers/01-directional-stability/analysis/verify_release.py
