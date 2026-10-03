# Predictor/target configuration accounting

The previous bank contains **19 prespecified predictor/target configurations**, not “19 representations” and not a full state × target cross-product:

1. One canonical state hash64 / canonical target hash64 configuration.
2. Seven additional state hash64 seeds, holding the canonical target fixed.
3. Three additional canonical-seed state dimensions: 32, 128 and 256, holding the canonical target fixed.
4. One training-fitted state PCA32 configuration, still predicting the canonical **hashed** target.
5. Seven additional target hash64 seeds, holding the canonical state hash64 fixed.

Thus 1 + 7 + 3 + 1 + 7 = 19. Original and E1 grids each compare jointMetric to independently tuned molecular kNN and history-augmented ridge to independently tuned stage-matched molecular ridge: 19 × 2 × 2 = 76 three-embryo multiplicity units, or 228 outer-fold history comparisons.

The new canonical state hash64 / training-fitted **target** PCA32 experiment is a prespecified non-hash outcome check, not a twentieth lineage-family member: no new history model was tested against it. Original state PCA32 does not substitute for target PCA32.

Configurations using different target coordinate systems are evaluated against their own matched baselines. Absolute NSE values from different target coordinate systems are not ranked against one another. Target PCA improves a within-frame comparison; it does not license “representation independent.”
