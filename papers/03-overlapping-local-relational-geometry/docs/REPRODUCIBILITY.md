# Paper 3 — Standalone reproducibility

The approved manuscript's complete Methods and detached supplement remain the scientific specification. Final08 is frozen: three E7.5 embryos / 3,704 cells, strict nine-lane 59/17/44 consensus, 36 independently calibrated families, 108 component and 72 aggregate/minimum decisions, unchanged original inferential banks.

## Numerical verification

From a clean repository clone, install `analysis/requirements.txt` in an isolated environment. NumPy 2.0.2 and SciPy 1.13.1 are the original numerical pins. The validated hosted runtime is Ubuntu 24.04/x86_64 with AVX2/FMA, Python 3.12.15, PyPI numerical wheels and single-thread OpenBLAS Haswell kernels. The wrapper pins this kernel automatically on Linux x86_64. Independent local replay also passed on Apple Silicon/macOS 27.0.1, Python 3.9.6 and Accelerate. Run:

```sh
python papers/03-overlapping-local-relational-geometry/analysis/verify_numerical.py
```

The gate checks the frozen source mapping and strict binary signatures, byte counts and SHA-256 before and after replay. It invokes the unchanged Final08 `scalar_ranks`, `matched_replay` and `degree_replay` actions. These reproduce sample-SD calibration (divisor 1023), epsilon 1e-8, inclusive ties and plus-one ranks /1025; all 108 component and 72 aggregate/minimum ranks; five matched banks / 30,720 graphs and two degree-scale banks / 3,072 shared pair graphs; and the complete calibrated decision ledger. No scientific bank or table is rewritten.

The geometry gate reconstructs q10/q25/q50 domains, active anchors, fibers and overlaps from released frozen coordinates and cores; verifies nine-lane provenance and 59/17/44 edges; independently replays all three compositions and 746 sourcewise existential-bridge checks. The figure gate preserves all 954 original Final08 displayed-statistic comparisons, ten one-page vector PDFs, matching PNGs and source CSV authorities.

`make verify` also runs the existing Papers 1 and 2 gates and Paper 3 document/table structure checks. Replay takes a few minutes on a typical CPU. The wrapper sets single-threaded BLAS. Generated receipts go under ignored `reproducibility/rmmo_paper3_final08/completion/`; optional `--receipt FILE` records only a completed successful run.

### Exact-replay runtime boundary

The full hosted replay passed all seven banks / 33,792 computations in [Actions run 37695015369](https://github.com/zedjames/relational-morphogenesis/actions/runs/37695015369), using Ubuntu 24.04, Python 3.12.15, the original NumPy/SciPy pins and fixed single-thread `OPENBLAS_CORETYPE=HASWELL`. The public wrapper requires AVX2/FMA and sets that fixed kernel on Linux x86_64; it never searches kernels or selects banks dynamically. This is execution-environment pinning, not a scientific algorithm, null-model or tolerance change.

Unrestricted runtime portability is not claimed. An unpinned Linux run 37692093364 passed, while unpinned runs 37692384954 and 37692395029 failed the frozen hash-S0 q25 first-draw assertion. Hosted macOS 15/ARM64/Accelerate also failed in runs 37693535965 and 37693539295, despite local macOS 27.0.1/ARM64/Python 3.9.6/Accelerate full success. Diagnostics showed identical frozen state bytes but different reconstructed decile assignments and fiber memberships before sampling. On one Linux runner, changing only the OpenBLAS kernel changed these contexts: Haswell and Zen exactly matched all six frozen first draws; Sandybridge did not. The divergence is therefore sensitive to numerical execution kernels; exact internal instruction/rounding differences were not further characterized. Failed runs remain visible, and unsupported runtimes still fail every original assertion.

`analysis/diagnose_replay_runtime.py` reports sanitized backend/version identifiers, frozen state and sampler-context hashes, and expected/actual first-draw hashes for all six orientations. It is read-only, does not generate a calibration bank, and never substitutes for or weakens the full gate. Linux document/packaging validation passed in runs 37693535965 and 37693539295.

The diagnostic-only kernel comparison is recorded in [Actions run 37694743965](https://github.com/zedjames/relational-morphogenesis/actions/runs/37694743965). Kernel pinning uses the documented [OpenBLAS runtime control](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/). A first-draw diagnostic alone is never release approval: all original per-draw graph hashes, numerical scores, ranks, geometry and figure comparisons must pass the full gate. No frozen scientific source, bank, statistic, rounding rule or comparison is changed.

## Document validation

```sh
python papers/03-overlapping-local-relational-geometry/analysis/verify_documents.py --engine latexmk --output /tmp/rmmo-paper3-documents
```

Install `latexmk` and TeX Live latex-extra, science and recommended fonts; alternatively pass a Tectonic executable. The gate compares complete approved document bodies and Methods with their hashes, then compiles both sources using six main and four supplementary figures. It rejects LaTeX errors, unresolved citations/references, overfull boxes and oversized floats. Only preamble figure paths, heading wrapping and URL wrapping differ from the approved sources. Compiled PDFs are validation builds outside the companion, never substitutes for the pending final scholarly PDF.

## Provenance and exclusions

`verification/source_export_manifest.json` maps 466 scientific destinations / 195,430,565 bytes to the 473-file / 196,288,214-byte Final08 source manifest. Six old combined-manuscript TeX/PDF/log files and the source Git LFS configuration are explicitly excluded. All exported scientific content is byte-identical. Historical original verifier/layout scripts and receipts are provenance; use public `analysis/` entry points, which distinguish numerical authority from the revised document layout.

Raw third-party H5TD files are absent. Original MELA v1 integrity records, licensing and attribution remain at [Zenodo record 19892785](https://doi.org/10.5281/zenodo.19892785). Released derived inputs suffice for this gate. Raw-carrier re-extraction is a distinct historical operation and is not claimed here. No private checkout, Git, raw archive, credentials or Lean is needed for numerical verification.

Final author PDF/checksum, version DOI, all-versions DOI and immutable release/tag remain pending.

## Independent clean-room receipt

`verification/cleanroom_validation_receipt.json` records a fresh public GitHub clone, a new dependency environment, all 95 Paper 1 and 25 Paper 2 checks, full Paper 3 replay and document compilation, including the updated public wrapper. `analysis/verify_git_payloads.py` additionally checks actual committed Git blobs: all 466 payloads have frozen sizes/hashes and none is an LFS pointer. This packaging gate requires Git; numerical verification itself does not. GitHub Actions requires both the complete pinned-kernel Linux numerical gate and the separate Linux document/packaging gate before merge.
