# Paper 3 — Standalone reproducibility

The approved manuscript's complete Methods and detached supplement remain the scientific specification. Final08 is frozen: three E7.5 embryos / 3,704 cells, strict nine-lane 59/17/44 consensus, 36 independently calibrated families, 108 component and 72 aggregate/minimum decisions, unchanged original inferential banks.

## Numerical verification

From a clean repository clone, create a Python 3.9–3.12 environment and install `analysis/requirements.txt`. NumPy 2.0.2 and SciPy 1.13.1 are the original numerical pins. Run:

```sh
python papers/03-overlapping-local-relational-geometry/analysis/verify_numerical.py
```

The gate checks the frozen source mapping and strict binary signatures, byte counts and SHA-256 before and after replay. It invokes the unchanged Final08 `scalar_ranks`, `matched_replay` and `degree_replay` actions. These reproduce sample-SD calibration (divisor 1023), epsilon 1e-8, inclusive ties and plus-one ranks /1025; all 108 component and 72 aggregate/minimum ranks; five matched banks / 30,720 graphs and two degree-scale banks / 3,072 shared pair graphs; and the complete calibrated decision ledger. No scientific bank or table is rewritten.

The geometry gate reconstructs q10/q25/q50 domains, active anchors, fibers and overlaps from released frozen coordinates and cores; verifies nine-lane provenance and 59/17/44 edges; independently replays all three compositions and 746 sourcewise existential-bridge checks. The figure gate preserves all 954 original Final08 displayed-statistic comparisons, ten one-page vector PDFs, matching PNGs and source CSV authorities.

`make verify` also runs the existing Papers 1 and 2 gates and Paper 3 document/table structure checks. Replay takes a few minutes on a typical CPU. The wrapper sets single-threaded BLAS. Generated receipts go under ignored `reproducibility/rmmo_paper3_final08/completion/`; optional `--receipt FILE` records only a completed successful run.

## Document validation

```sh
python papers/03-overlapping-local-relational-geometry/analysis/verify_documents.py --engine latexmk --output /tmp/rmmo-paper3-documents
```

Install `latexmk` and TeX Live latex-extra, science and recommended fonts; alternatively pass a Tectonic executable. The gate compares complete approved document bodies and Methods with their hashes, then compiles both sources using six main and four supplementary figures. It rejects LaTeX errors, unresolved citations/references, overfull boxes and oversized floats. Only preamble figure paths, heading wrapping and URL wrapping differ from the approved sources. Compiled PDFs are validation builds outside the companion, never substitutes for the pending final scholarly PDF.

## Provenance and exclusions

`verification/source_export_manifest.json` maps 467 scientific destinations / 195,430,732 bytes to the 473-file / 196,288,214-byte Final08 source manifest. Six old combined-manuscript TeX/PDF/log files are explicitly excluded. All exported scientific content is byte-identical. Historical original verifier/layout scripts and receipts are provenance; use public `analysis/` entry points, which distinguish numerical authority from the revised document layout.

Raw third-party H5TD files are absent. Original MELA v1 integrity records, licensing and attribution remain at [Zenodo record 19892785](https://doi.org/10.5281/zenodo.19892785). Released derived inputs suffice for this gate. Raw-carrier re-extraction is a distinct historical operation and is not claimed here. No private checkout, Git, raw archive, credentials or Lean is needed for numerical verification.

Final author PDF/checksum, version DOI, all-versions DOI and immutable release/tag remain pending.
