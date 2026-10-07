#!/usr/bin/env python3
"""Verify the paper-specific pre-DOI PUBLIC STAGING surface.

This checks document identity, manuscript/supplement separation, source
metadata and main result references. It is intentionally not a numerical
replay of the scientific reference banks.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def main():
    manuscript = (ROOT / "manuscript/RMMO_Paper3_Manuscript_preDOI.tex").read_text()
    supplement = (ROOT / "supplement/SupplementaryEvidence.tex").read_text()
    source = json.loads((ROOT / "data/authority/source_normalization_provenance.json").read_text())

    require(r"\\section{Methods}" in manuscript, "Methods section missing")
    require(r"\\section{Results}" in manuscript, "Results section missing")
    require(r"\\section{Discussion}" in manuscript, "Discussion section missing")
    require(r"\\section{Data availability}" in manuscript, "Data availability missing")
    require(r"\\section{Code availability}" in manuscript, "Code availability missing")
    require(r"\\appendix" not in manuscript, "Detached appendix remains in manuscript")
    require(r"\\appendix" in supplement, "Supplement is missing appendix boundary")

    required_methods = (
        "Carrier construction and coarse state",
        "Prespecified gene partition",
        "Own-lane normalization and signed-hash projection",
        "Construction/evaluation firewall",
        "State distance, lane thresholds, and consensus edges",
        "Domains, active anchors, fibers, and overlaps",
        "Independent calibration of component and family scores",
        "Joint family-draw construction",
        "Retention references",
        "Exact relation composition and historical common-refinement criterion",
    )
    for name in required_methods:
        require(name in manuscript, "Missing Methods subsection: " + name)
    figures = re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", manuscript)
    require(len(figures) == 6, "Expected six primary figure references")
    for index in range(1, 7):
        require(any(f"Paper3_Figure{index}_" in f for f in figures),
                "Missing primary Figure " + str(index))

    expected_pairs = ("P12", "P13", "P23")
    for lane in range(9):
        for pair in expected_pairs:
            require(f"$\\ell_{lane}$&{pair}&" in supplement,
                    f"Missing supplementary lane {lane}/{pair}")
    for fig in ("FigureS1_", "FigureS2_", "FigureS3_", "FigureS4_"):
        require(fig in supplement, "Missing supplementary figure " + fig)

    raw = source.get("raw_sources", [])
    require(len(raw) == 3, "Expected three original source records")
    for row in raw:
        require(re.fullmatch(r"[a-f0-9]{64}", row["sha256"]) is not None,
                "Invalid source SHA-256")
        require(int(row["bytes"]) > 0, "Invalid source byte count")
    require(source.get("raw_H5TD_included_in_release") is False,
            "Source-data coverage statement changed")

    print("PASS — RMMO Paper 3 pre-DOI document and source-metadata checks")
    print("Methods retained; appendix detached; main figures referenced: 6")
    print("Supplementary lane/pair authorities: 27; raw source records: 3")
    print("LIMIT: no independent numerical replay or figure-binary verification")


if __name__ == "__main__":
    main()
