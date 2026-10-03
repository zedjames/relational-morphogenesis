# Inner-fold eligibility and panel construction

For outer training on two embryos, each eligible gene must be detected (count > 0) in at least 1% of cells in **each** of those embryos. It must have positive mean count across fit-owned cells, and its source gene name must not begin with `mt-`, `Mt-`, `MT-`, `Rpl` or `Rps`. No held-out embryo contributes to eligibility.

Within nested tuning the fit set contains **one** embryo. Apply the same detection >= 0.01 rule to that one available inner-fit embryo, plus positive fit-set mean and the same exclusions. Recompute the eligible vocabulary, stable-ID hash partition, panel-specific library denominators and representations. The inner-validation embryo contributes neither its gene detection fractions nor target values to feature discovery or frame fitting. The third outer-test embryo contributes to none of these operations or tuning scores.

The state panel consists of eligible genes for which the first-eight-byte big-endian SHA256 integer of `RMMO12_STATE_REALIZATION_V1|geneId` is < 7 modulo 10; remaining genes form the target panel. Retain versioned deposited `gene_ids`, original source gene order and separate log1p(10000 × count / panel-total) normalizations. A zero panel total gives a zero normalized row.

`utils.Bank.matrices` implements exactly these rules and asserts exact ordered state/target gene indices against all six immutable fit vocabularies. `target_pca.py` fits each of the three single-embryo inner target frames and three two-embryo outer target frames separately. The lineage null reuses the already fitted, outcome-invariant molecular frames but repeats every history-bearing fit and tuning candidate under each assignment.
