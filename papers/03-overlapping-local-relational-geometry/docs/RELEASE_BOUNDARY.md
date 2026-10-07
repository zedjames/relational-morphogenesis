# Paper 3 — Public release boundary

This companion exports the frozen paper-specific Final08 reviewer package and approved public documents. Of 473 source-manifest files, 466 scientific files are released byte-exactly. Six earlier combined-manuscript TeX/PDF/log files and the source Git LFS configuration are explicitly excluded. The approved complete Methods and detached supplement remain intact; an older PDF is never published as the final scholarly PDF.

Released materials include ten vector figure PDFs with exact names, PNGs, source CSVs and generation provenance; derived carrier/geometry inputs; unchanged inferential banks; independent calibration banks/constants; original/calibrated comparison tables; and standalone replay. Source/export, binary and complete release SHA-256 manifests serve distinct checks. Pointer files fail the binary gate.

No private worktree, Git history, remotes, credentials, configuration, unrelated scientific program or raw MELA H5TD archive is exported. Source commit identifiers are provenance hashes only; public Git history is independent. The private source repository remains private. Historical paper-specific scripts are retained as frozen provenance; only documented public entry points run during release verification.

MELA v1 ([Zenodo 19892785](https://doi.org/10.5281/zenodo.19892785)) retains its original rights. Author manuscript, figures and original result tables follow CC BY-NC-ND 4.0; author code follows PolyForm Noncommercial 1.0.0. Derived third-party inputs are not blanket-relicensed.

Local public-clone full replay and document compilation passed. Hosted full replay also passed after pinning the Haswell OpenBLAS execution kernel, without changing any scientific authority, algorithm or comparison. Unsupported runtime portability is not claimed; failed diagnostics and the fixed execution environment are documented in `REPRODUCIBILITY.md`. Public integration still requires all final-head numerical, packaging, privacy and document checks before merge.

Publication fields also remain pending: final author-supplied scholarly PDF and SHA-256, version DOI, all-versions DOI, and associated immutable release/tag. Compilation PDFs remain validation outputs outside the release. No final Paper 3 tag is created here.
