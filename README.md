# Relational Morphogenesis and Multicellular Organization (RMMO)

Public reproducibility materials for the **Relational Morphogenesis and Multicellular Organization (RMMO)** research series.

## Current public paper

**Paper 1 — Relational Organization and Directional Stability Across Mouse Organogenesis**

- Author: Zed James
- Preprint date: October 1, 2026
- Zenodo DOI: https://doi.org/10.5281/zenodo.23090479
- Public companion: `papers/01-directional-stability/`

The scientific manuscript is frozen on Zenodo. This repository contains the finite manuscript-facing result surfaces, provenance, sensitivity analyses, verification code, exact manuscript/figure release assets and checksums, and publication-facing formal controls for Paper 1.

## Purpose

RMMO studies multicellular development through relational organization: which biological states are comparable, what relational structure persists, how direction is expressed, and which levels of description remain stable across developmental change.

This repository is the public publication and reproducibility surface for the series. It is intentionally smaller than the private development environment. Each paper receives a self-contained folder only after its scientific object and publication boundary have been independently hardened.

## Publication principles

1. **Every paper must stand independently.**
2. **The series-level story is a hypothesis under attack.**

See `docs/PUBLICATION_PRINCIPLES.md`.

## Quick start

    make verify

## Structure

    docs/                              series-level publication and release policy
    papers/01-directional-stability/  self-contained Paper 1 public companion

## Source data

Third-party raw data are not mirrored. Paper 1 uses the Mouse Organogenesis Spatiotemporal Transcriptomic Atlas (MOSTA), dataset **STDS0000058**, raw project **CNP0001543**.

## Repository boundary

This repository releases paper-specific empirical results, provenance, validation code, figure inputs/release assets, manuscript metadata, and publication-facing formal controls. It does not release the broader private research environment, unrelated theorem libraries, generic orchestration systems, or research programs outside RMMO.

## Licensing

See `LICENSE.md`.

- manuscripts, figures, documentation, and original derived result tables: **CC BY-NC-ND 4.0**
- paper-specific executable code and build helpers: **PolyForm Noncommercial 1.0.0**
- third-party materials retain their original licenses and terms
