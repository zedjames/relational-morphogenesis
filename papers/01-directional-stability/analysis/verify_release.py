#!/usr/bin/env python3
"""Verify the finite public result surfaces for RMMO Paper 1."""

from __future__ import annotations
import csv
import json
import math
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
checks = []

def check(condition, label):
    if not condition:
        raise SystemExit(f"FAIL: {label}")
    checks.append(label)

def read_json(relative):
    return json.loads((PAPER / relative).read_text())

def read_csv(relative):
    with (PAPER / relative).open(newline="") as handle:
        return list(csv.DictReader(handle))

def close(a, b, tol=1e-12):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)

source = read_json("data/authority/source_hashes.json")
check(source["officialDataset"] == "STDS0000058", "MOSTA dataset authority")
check(source["rawProject"] == "CNP0001543", "MOSTA raw-project authority")
check(source["primaryCarrierSectionCount"] == 53, "53 official sections")
check(source["primaryCarrierRootCount"] == 4_099_397, "4,099,397 spatial roots")

stages = read_csv("data/authority/stage_registry.csv")
check(len(stages) == 8, "eight developmental stages")
check(sum(int(r["sectionCount"]) for r in stages) == 53, "stage registry section total")
check(sum(int(r["biologicalEmbryoCount"]) for r in stages) == 37, "37 biological embryos")
check(sum(int(r["rootCount"]) for r in stages) == 4_099_397, "stage registry root total")

annotations = read_csv("data/authority/annotation_registry.csv")
check(len(annotations) == 49, "49 source annotation classes")
check(sum(int(r["rootCount"]) for r in annotations) == 4_099_397, "annotation registry root total")

carrier = read_csv("results/carrier/carrier_robustness_summary.csv")
by_variant = {r["variant"]: r for r in carrier}
check(set(by_variant) == {
    "historical_k32", "k16", "k64", "sentinel_k32",
    "densityAblated_k32", "rootAnnotationOnlyGate_k32"
}, "six released carrier variants")
check(close(by_variant["historical_k32"]["chronologicalChainScore"], -0.13977910974666596), "primary chain score")
check(float(by_variant["k16"]["chronologicalChainScore"]) < 0, "k16 observed chain score negative")
check(float(by_variant["k64"]["chronologicalChainScore"]) < 0, "k64 observed chain score negative")
check(float(by_variant["k16"]["bootstrapQ975"]) > 0, "k16 bootstrap interval crosses zero")
check(int(by_variant["historical_k32"]["localResolvedCount"]) == 0, "primary 0/7 local resolution")
check(int(by_variant["k16"]["localResolvedCount"]) == 0, "k16 0/7 local resolution")
check(int(by_variant["k64"]["localResolvedCount"]) == 0, "k64 0/7 local resolution")
check(int(by_variant["rootAnnotationOnlyGate_k32"]["localResolvedCount"]) == 2, "root-only gate 2/7 local resolution")
check(int(by_variant["rootAnnotationOnlyGate_k32"]["relationPairCountAdjacentTotal"]) == 144_921_780, "root-only relation expansion")

boot = read_csv("results/direction/carrier_scale_bootstrap_summary.csv")
boot_by = {r["variant"]: r for r in boot}
check(int(boot_by["historical_k32"]["bootstrapReplicates"]) == 4096, "primary 4096 bootstrap replicates")
check(close(boot_by["historical_k32"]["negativeSignFrequency"], 0.999755859375), "primary bootstrap sign frequency")
check(close(boot_by["k16"]["negativeSignFrequency"], 0.969482421875), "k16 bootstrap sign frequency")

tol = read_csv("results/sensitivity/tolerance_summary.csv")
check([r["epsilon"] for r in tol] == ["0.15", "0.25", "0.35"], "three registered tolerances")
check(all(float(r["chronologicalPathScore"]) < 0 for r in tol), "negative chain score at all tolerances")
check(all(int(r["edgeSpecificResolvedCount"]) == 0 for r in tol), "0/7 endpoint resolution at all tolerances")
check(all(close(r["mappedFullSignConcordance"], 1.0) for r in tol), "mapped/full sign preservation at all tolerances")

mapped = read_csv("results/sensitivity/sink_mapped_decomposition.csv")
check(len(mapped) == 28, "28 pairwise mapped/sink sensitivity rows")
check(sum(r["signMappedEqualsFull"] == "True" for r in mapped) == 28, "28/28 mapped/full sign concordance")

leave = read_csv("results/sensitivity/leave_one_edge.csv")
check(len(leave) == 7, "seven leave-one-edge analyses")
check(all(r["remainingSign"] == "negative" for r in leave), "7/7 leave-one-edge negative sign")

effects = read_json("results/direction/directional_error_summary.json")
check(effects["identityVerified"] is True, "directional error identity")
check(close(effects["meanOmega"], -0.019968444249523678), "mean Omega")

curve = read_csv("results/composition/curveball_summary.csv")
check(len(curve) == 6, "six composition triples")
check(close(max(abs(float(r["absoluteDifference"])) for r in curve), 0.0007732060236939065), "maximum Curveball residual")

multi = read_csv("results/composition/multiple_testing.csv")
check(len(multi) == 6, "six multiple-testing rows")
check(sum(float(r["holmAdjusted"]) < 0.05 for r in multi) == 5, "five Holm-adjusted composition tails below 0.05")

sign = read_json("results/calibration/exact_sign_null.json")
check(close(sign["twoSidedProbability"], 0.15625), "exact sign-flip reference")

centered = read_json("results/calibration/centered_bootstrap_null.json")
check(close(centered["variants"]["bootstrapMean"]["twoSidedEmpiricalP"], 0.000244081034903588), "bootstrap-mean centered reference")
check(close(centered["variants"]["observedVector"]["twoSidedEmpiricalP"], 0.000488162069807176), "observed-vector centered reference")

pair_ref = read_json("results/calibration/stage_balanced_accumulation_null.json")
check(close(pair_ref["twoSidedP"], 0.2317000133513895), "stage-balanced embryo-pair reference-tail fraction")

density = read_csv("results/sensitivity/density_bin_occupancy.csv")
check([int(r["rootCount"]) for r in density] == [453601, 0, 0, 3645796], "density-bin occupancy")

root_only = read_json("results/sensitivity/root_only_gate_summary.json")
check(root_only["negativeChainSignPreserved"] is True, "root-only negative aggregate sign")
check(root_only["endpointResolvedCount"] == 2, "root-only two endpoint-reference exceedances")

formal = read_json("verification/formal_verification.json")
check(formal["passed"] is True, "private formal verification receipt passed")
check(formal["sorryCount"] == 0 and formal["admitCount"] == 0, "zero sorry/admit")
check(formal["newScientificAxioms"] is False, "zero new scientific axioms")

print(f"PASS: {len(checks)} Paper 1 public verification checks")
for item in checks:
    print(f"  - {item}")
