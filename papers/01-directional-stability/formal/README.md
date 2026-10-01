# Formal Publication Controls

The `standalone/` directory contains publication-control Lean files whose only external library dependency is Mathlib.

Additional publication checks import broader private RMMO research modules. Those dependency-backed sources are not mirrored here. Their public authority is represented by verification receipts and the exact completion commits in `../docs/METHODS_PROVENANCE.md`.

The frozen publication surface records zero `sorry`, zero `admit`, and zero project-specific scientific axiom declarations. Public finite-result verification is performed by `make verify`; it does not require the private theorem environment.
