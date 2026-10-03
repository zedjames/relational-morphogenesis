# Reproducibility

## Finite verification

From repository root:

    make verify

Paper 2 finite verification uses the Python standard library only. It checks the released manuscript-facing tables and decision files, including the canonical molecular result, representation sensitivity, target-PCA result, history exposure, lineage family calibration, complete discrete diagnostic family, and finite-state qualification surface.

## Raw replay

The released analysis code under `analysis/` is the publication-facing replay specification. For raw-data replay, obtain the three hash-identified MELA H5TD files from the official archive listed in `data/authority/source_authority.json`.

Replay dependencies:

    pip install -r papers/02-present-state-resolution-and-lineage-history/analysis/replay_requirements.txt

The full lineage-family null uses the prespecified 999 synchronized conditional assignments with nested retuning. The public finite result tables preserve the complete observed and null summary surfaces required to inspect the manuscript claims.

## Biological replication

The biological replicate count is three embryos. Cell pairs and permutation replicates are finite descriptive or conditional quantities and do not increase biological n.
