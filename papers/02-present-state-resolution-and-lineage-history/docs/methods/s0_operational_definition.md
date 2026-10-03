# Historical and operational coarse state

The historical key remains `S0 = (germ_layer, phase, type)` exactly as used in the source join and prior computations. Source `type` denotes donor/host origin, not `cell_type`.

Every one of the selected 3,704 cells has type = donor. Therefore all operational S0 variation in this selected population is through **germ layer × cell-cycle phase**. Retain the historical three-field key for provenance and exact replay rather than retrospectively dropping its constant component.

Germ-layer and phase labels derive from source expression annotation/scoring. The S0 predictor is an expression-derived coarse comparator, not an independent phenotype or causal state variable. Its target prediction is the mean of fit-owned target vectors among matching keys, falling back to the fit-global target mean only if the key is absent.
