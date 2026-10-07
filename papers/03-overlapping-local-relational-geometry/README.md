# Paper 3 — An Overlapping Local Relational Geometry of Early Mouse Embryogenesis

**Relational Morphogenesis and Multicellular Organization (RMMO), Paper 3**  
**Zed James | October 2026**

**Publication status: published Zenodo preprint; public scientific companion verified; Paper 3 release tag pending.** The version DOI is [10.5281/zenodo.23225035](https://doi.org/10.5281/zenodo.23225035). The author's final scholarly PDF is designated as the canonical Zenodo artifact; its submitted copy has SHA-256 `5d507094748517024a7f06d41716435bdbf9137601268f5593b585d88cb4a95f`.

The [frozen, previously verified manuscript source](manuscript/RMMO_Paper3_Manuscript_preDOI.tex) preserves its complete document body and full Methods. The [publication-facing LaTeX source](manuscript/RMMO_Paper3_Manuscript_publication_updated.tex) updates availability statements while retaining the scientific document. The [detached supplement](supplement/SupplementaryEvidence.tex) retains the former end appendix. All ten actual frozen scientific figure PDFs and source-data tables are in [figures/](figures/). Only preamble paths, heading wrapping and URL wrapping were adapted for public compilation.

## Scientific scope

The study describes three lineage-qualified E7.5 mouse embryo carriers (3,704 cells) using pair-local relations, q25 source domains, set-valued target fibers and spatial overlaps. A frozen G_S-only nine-lane construction (14,660 genes) precedes evaluation in disjoint G_R (6,339 genes). The strict consensus has 59/17/44 edges (P12/P13/P23). The principal results are construction-held-out molecular coherence, connectivity-controlled source fidelity, spatial overlap and a tested failure of exact global composition.

## Public evidence status

- Released: 466 byte-exact Final08 scientific files, 195,430,565 bytes, including derived inputs, original inferential banks, independent calibration banks, numerical replay code, ten vector figures and their source tables.
- Independent verification: seven authorized calibration banks; 33,792 exact graph computations; 108 component rankings; 72 aggregate/minimum rankings; all 36 families; zero changed support decisions or inferential banks; exact finite composition and all 746 cellwise bridge checks.
- Document validation: approved manuscript and detached supplement compile using the actual figures. Validation-build PDFs are not published as the final scholarly PDF.
- Publication metadata: Zenodo version DOI `10.5281/zenodo.23225035`; checksum of the author-submitted final PDF `5d507094748517024a7f06d41716435bdbf9137601268f5593b585d88cb4a95f` (715930 bytes). Zenodo all-versions DOI has not been independently resolved, and the immutable Paper 3 release tag remains to be created.

The full verifier uses only these public files and declared software dependencies. Raw MELA acquisition, a private checkout, Git and Lean are unnecessary for numerical verification. The manifest records six excluded older combined-manuscript files and one excluded source Git LFS configuration; no scientific dependency was omitted.

## Run verification

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r papers/03-overlapping-local-relational-geometry/analysis/requirements.txt
make verify PYTHON=.venv/bin/python
.venv/bin/python papers/03-overlapping-local-relational-geometry/analysis/verify_documents.py --engine latexmk --output /tmp/rmmo-paper3-documents
```

Numerical and document validation are separate gates. LaTeX needs `latexmk`, TeX Live latex-extra, science and recommended fonts, or Tectonic. See [reproducibility](docs/REPRODUCIBILITY.md) for coverage and the [source-to-destination manifest](verification/source_export_manifest.json).

Full hosted exact replay passes on Ubuntu 24.04/x86_64, Python 3.12.15 and fixed Haswell OpenBLAS kernels (AVX2/FMA required). The wrapper pins the execution kernel automatically on Linux x86_64. Independent local macOS replay also passes. Unrestricted backend portability is not claimed; see the documented [runtime boundary and validation](docs/REPRODUCIBILITY.md). Frozen scientific code, banks and strict comparisons remain unchanged.

## Source and licensing

Primary third-party resource: [MELA v1](https://doi.org/10.5281/zenodo.19892785). Raw third-party H5TD files are not mirrored. The public companion is a selective scientific release, not a copy of any research worktree or history. See the repository's [licensing policy](../../LICENSE.md).

## Release documentation

- [Reproducibility and data coverage](docs/REPRODUCIBILITY.md)
- [Result-to-artifact map](docs/RESULT_MAP.md)
- [Release boundary and outstanding items](docs/RELEASE_BOUNDARY.md)

## Published preprint and release handoff

- [Zenodo paper (version DOI)](https://doi.org/10.5281/zenodo.23225035)
- Author-submitted final PDF SHA-256: `5d507094748517024a7f06d41716435bdbf9137601268f5593b585d88cb4a95f`
- Publication-facing LaTeX source: [`manuscript/RMMO_Paper3_Manuscript_publication_updated.tex`](manuscript/RMMO_Paper3_Manuscript_publication_updated.tex)
- The scholarly PDF is cited at Zenodo, as in Paper 2's GitHub publication model. No repository copy of the PDF is included in this metadata commit.
- The scientific inputs, replay banks, exact statistics, and previously verified manuscript/supplement sources have not been altered.
- Pending final handoff: independently confirm the Zenodo all-versions DOI, pass CI at the publication commit, and create the immutable Paper 3 release tag.
