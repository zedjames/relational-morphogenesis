# Historical witness controls and simultaneous correction

The pair universe consists of all 1,556,131 unordered same-S0 pairs, including both within- and cross-embryo pairs. The historical witness W requires both different source `cell_type` and different frozen L4 history. It has 780,783 pairs. Source cell_type is expression-derived; this is not an independent phenotype.

C0 includes all same-S0 pairs. C1 requires cell_type discordance without requiring history discordance. C2 requires L4 discordance without requiring annotation discordance. C3 is the **non-witness complement** matched by S0 and unordered embryo pair. C3A is the same complement matched additionally by unordered cell_type pair. Neither control is obtained by assuming pairs are independent observations.

For each matching stratum s, let w_s and c_s be its witness and complement counts. Restrict to common support w_s > 0 and c_s > 0. Each control pair receives weight w_s/c_s; each supported witness pair has unit weight. Witness survival is recomputed on the same supported witness subset as its particular control. C3A covers 98.6443% of W; unsupported witnesses are not silently assigned zero survival. C3 and C3A have different estimands/support.

Freeze the five full-bank state codes at K = 8,16,32,64,128. For each of 128 seed IDs, jointly reassign the complete five-code vector at the **cell** level within (S0, embryo), using one permutation index across all five K values. Keep annotations, histories, W membership, matching weights, embryo membership and shared-cell pair structure fixed. Never independently permute pair outcomes or code columns.

For each K and control, define delta_K = supported-witness survival minus weighted-control survival. The directional preference statistic is -delta_K. In each null realization record max_K(-delta_K). The K-specific simultaneous probability is `(1 + count(null maximum >= -observed delta_K))/129`, with inclusive ties and attainable minimum 1/129 ≈ 0.007752. C3's five-K maximum and C3A's five-K maximum are separately recorded controls; the reported C3A value is **not** a correction over both control families or all later analyses.

At K32, full W survival is 0.0313378 and C1 survival is 0.0313525, giving ratio ≈ 0.99953. C3A common-support witness survival is 0.0315789 and control survival 0.0415866, so the observed preference is 0.0100077. The null maximum reaches or exceeds this value in 67 of 128 realizations; `(1+67)/129 = 0.5271318`. The tighter historical comparison therefore has no corrected history-specific departure. Conditional permutation assignments do not increase the biological replicate count beyond three embryos.

This witness family is distinct from the immutable 100-contrast diagnostic maximum (family p = 0.31008) and from the new 76-member prospective predictive lineage maximum. Do not transfer significance from one family to another.
