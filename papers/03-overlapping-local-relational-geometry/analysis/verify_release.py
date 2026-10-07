#!/usr/bin/env python3
"""Verify the paper-specific pre-DOI PUBLIC STAGING surface.

This checks document identity, manuscript/supplement separation, source
metadata and main result references. It is intentionally not a numerical
replay of the scientific reference banks.
"""
import csv
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

    require(r"\section{Methods}" in manuscript, "Methods section missing")
    require(r"\section{Results}" in manuscript, "Results section missing")
    require(r"\section{Discussion}" in manuscript, "Discussion section missing")
    require(r"\section{Data availability}" in manuscript, "Data availability missing")
    require(r"\section{Code availability}" in manuscript, "Code availability missing")
    require(r"\appendix" not in manuscript, "Detached appendix remains in manuscript")
    require(r"\appendix" in supplement, "Supplement is missing appendix boundary")

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
            require(f"$\ell_{lane}$&{pair}&" in supplement,
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

    expected = {
        "nine_lane_construction.csv": 27,
        "coarse_state_counts.csv": 11,
        "fibers_q25.csv": 6,
        "heldout_edge_agreement.csv": 3,
        "matched_target_q25.csv": 6,
        "degree_controlled_source_fidelity_q25.csv": 6,
        "spatial_overlap_q25.csv": 3,
        "construction_identity_sensitivity.csv": 5,
        "single_lane_representation_sensitivity.csv": 3,
        "domain_scale_sensitivity.csv": 6,
        "matched_edge_pool_diagnostics.csv": 3,
        "matched_target_pool_diagnostics.csv": 6,
        "exact_composition_counts.csv": 3,
        "existential_bridge_counts.csv": 3,
    }
    tables = ROOT / "results/manuscript_tables"
    imported = {}
    for filename, expected_rows in expected.items():
        with (tables / filename).open(newline="") as stream:
            data = list(csv.DictReader(stream))
        require(len(data) == expected_rows,
                f"{filename}: expected {expected_rows} rows; got {len(data)}")
        imported[filename] = data

    lanes = imported["nine_lane_construction.csv"]
    require({r["Lane"] for r in lanes} == {f"lane_{j}" for j in range(9)},
            "Nine-lane authority mismatch")
    require(all({r["Pair"] for r in lanes if r["Lane"] == f"lane_{j}"}
                == {"P12", "P13", "P23"} for j in range(9)),
            "Expected three pair thresholds per lane")

    coarse = imported["coarse_state_counts.csv"]
    for key, expected_total in (("E1", 1032), ("E2", 1298), ("E3", 1374),
                                ("Total", 3704)):
        require(sum(int(r[key]) for r in coarse) == expected_total,
                "Coarse-state carrier total mismatch: " + key)

    fidelity = imported["degree_controlled_source_fidelity_q25.csv"]
    o1 = [r for r in fidelity if r["Panel"] == "G_R"
          and r["Overlap"] == "O1"]
    require(len(o1) == 1 and abs(float(o1[0]["maxT-adjusted p"]) -
                0.122927) < 1e-7, "Held-out O1 fidelity authority changed")

    for row in imported["exact_composition_counts.csv"]:
        require(int(row["Direct"]) == int(row["Common"]) +
                int(row["Direct only"]), "Direct relation decomposition invalid")
        require(int(row["Composed"]) == int(row["Common"]) +
                int(row["Composed only"]), "Composed decomposition invalid")
        require(row["Equal"] == "No", "Exact global boundary changed")

    for row in imported["existential_bridge_counts.csv"]:
        require(int(row["Bridge present"]) + int(row["Failing"]) ==
                int(row["Tested"]), "Bridge counts do not close")
        require(row["Complete"] == "No", "Existential bridge scope changed")

    print("PASS — RMMO Paper 3 pre-DOI document and source-metadata checks")
    print("Methods retained; appendix detached; main figures referenced: 6")
    print("Supplementary lane/pair authorities: 27; raw source records: 3")
    print(f"Manuscript-derived public CSV tables checked: {len(expected)}")
    print("LIMIT: no independent numerical replay or figure-binary verification")


if __name__ == "__main__":
    main()
