# Relational Morphogenesis and Multicellular Organization (RMMO)

Public reproducibility materials for the **Relational Morphogenesis and Multicellular Organization (RMMO)** research series.

## Current public release

The repository currently contains the prepublication computational companion for:

**Paper 1 — Relational Organization and Directional Stability Across Mouse Organogenesis**

The final manuscript source, PDF, DOI, release tag, and manuscript hash will be added when the paper is frozen for publication.

## Purpose

RMMO studies multicellular development through relational organization: which biological states are comparable, what relational structure persists, how direction is expressed, and which levels of description remain stable across developmental change.

This repository is the public publication and reproducibility surface for the series. It is intentionally smaller than the private development environment. Each paper receives a self-contained folder only after its scientific object and publication boundary have been independently hardened.

The public reproducibility direction is:

    public source authority
            ↓
    paper-specific finite carriers / released summaries
            ↓
    registered controls and sensitivity surfaces
            ↓
    manuscript claims, tables, and figures

## Publication principles

Two rules govern the series:

1. **Every paper must stand independently.**
2. **The series-level story is a hypothesis under attack.**

See `docs/PUBLICATION_PRINCIPLES.md`.

## Quick start

    make verify

## Structure

    docs/                              series-level publication and release policy
    papers/01-directional-stability/  self-contained Paper 1 public companion

Future paper folders will be added only when those papers reach their own publication-ready boundaries.

## Source data

Third-party raw data are not mirrored. Paper 1 uses the Mouse Organogenesis Spatiotemporal Transcriptomic Atlas (MOSTA), dataset **STDS0000058**, raw project **CNP0001543**.

## Repository boundary

This repository releases paper-specific empirical results, provenance, validation code, figure inputs, and publication-facing formal controls. It does not release the broader private research environment, unrelated theorem libraries, generic orchestration systems, or research programs outside RMMO.

## Licensing

See `LICENSE.md`.

- manuscripts, figures, documentation, and original derived result tables: **CC BY-NC-ND 4.0**
- paper-specific executable code and build helpers: **PolyForm Noncommercial 1.0.0**
- third-party materials retain their original licenses and terms
