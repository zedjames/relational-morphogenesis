# Paper 3 — Reproducibility and data coverage

The authoritative experimental specification remains the manuscript's complete Methods section. It specifies:

1. Source carrier membership, 3,704 lineage-qualified cells and coarse metadata categories;
2. Fixed gene eligibility, disjoint G_S/G_R feature partition (14,660 / 6,339);
3. Own-lane normalization, nine hashed construction lanes, q50 pair-symmetric thresholds and seven-of-nine consensus;
4. q25 source domains, active landmarks, set-valued target fibers and source overlaps (with q10/q50 sensitivities);
5. Fixed-stratum spatial, matched-edge, matched-target, and degree-only simple-graph reference laws;
6. Independent B_cal=1,024 same-law calibration, B=1,024 synchronized inferential banks, plus-one empirical tails and family maxT / T_sum / T_min;
7. Competitive kNN/MNN baselines, single-lane PCA controls, R0–R3 retained-structure references, and exact routewise composition/existential-bridge tests.

## Verification now available

Run from repository root:

```sh
python3 papers/03-overlapping-local-relational-geometry/analysis/verify_release.py
```

The pre-DOI verifier checks the presence and internal consistency of manuscript, Methods, supplement, primary figure references, and declared MELA source records. It **does not independently recompute** the numerical experiment. Such a replay requires the verified binary banks and source/figure evidence described in [release boundary](RELEASE_BOUNDARY.md).

The manuscript reports 108/108 preserved component maxT decisions and 72/72 aggregate/minimum decisions after separate same-law calibration. Those are reported experimental findings awaiting transfer of the corresponding audit banks for independent public replay.

## Original source

MELA v1: https://doi.org/10.5281/zenodo.19892785. The official archive retains raw source data; the local metadata file documents original filenames, byte counts and SHA-256.

## Expected frozen publication receipt

After the DOI and PDF are available, the companion should record final version DOI and all-versions DOI, PDF SHA-256, manuscript source SHA-256, exact evidence manifest, passing full verifier output, public commit SHA, and immutable release tag. These values are not assigned in advance.
