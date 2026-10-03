#!/usr/bin/env python3
"""Verify the finite public result surfaces for RMMO Paper 2."""
from __future__ import annotations
import csv, json, math
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
checks=[]

def check(v,label):
    if not v:
        raise SystemExit(f"FAIL: {label}")
    checks.append(label)

def rows(rel):
    with (PAPER/rel).open(newline="") as f:
        return list(csv.DictReader(f))

def js(rel):
    return json.loads((PAPER/rel).read_text())

def close(a,b,tol=1e-12):
    return math.isclose(float(a),float(b),rel_tol=0.0,abs_tol=tol)

pub=js("verification/publication_authority.json")
check(pub["versionDOI"]=="10.5281/zenodo.23125763","version DOI")
check(pub["allVersionsDOI"]=="10.5281/zenodo.23125762","all-versions DOI")
check(pub["canonicalPDF"]["sha256"]=="216707a03c0abde7bf9e6f3cfd12598830b2d2d2fd0ec07d40c66a194bc33b92","canonical PDF SHA-256")
check(pub["canonicalPDF"]["bytes"]==629019,"canonical PDF byte count")

src=js("data/authority/source_authority.json")
check(src["biologicalReplicates"]==3,"three biological embryos")
check(len(src["rawFiles"])==3,"three frozen MELA source files")

w=rows("results/witness/relative_survival.csv")
k32c1=next(r for r in w if r["analysis"]=="historical_annotation" and r["K"]=="32" and r["control"]=="C1")
check(close(k32c1["witnessSurvival"],0.03133777246686979),"K32 witness survival")
check(close(k32c1["controlSurvival"],0.03135251274115473),"K32 annotation control survival")
check(close(k32c1["ratio"],0.9995298534951048),"K32 W/C1 ratio")

pred=rows("results/prediction/baseline_comparison.csv")
canonical=[r for r in pred if r["variant"]=="hash64_seed0"]
check(len(canonical)==3 and all(float(r["gainOverS0"])>0 for r in canonical),"canonical molecular-over-S0 in 3/3 embryos")
check(all(float(r["gainOverGlobalMean"])>0 for r in pred),"all 19 configurations beat matched global mean")
h32=[r for r in pred if r["variant"]=="hash32_seed0"]
check(sum(float(r["gainOverS0"])>0 for r in h32)==1,"hash32 S0 sensitivity retained")

tp=rows("results/target_pca/baseline_comparison.csv")
check(len(tp)==3 and all(float(r["gain_over_S0"])>0 for r in tp),"target PCA molecular-over-S0 in 3/3 embryos")

ex=rows("results/history/exposure_summary.csv")
l4=[r for r in ex if r["historyLevel"]=="L4"]
check(len(l4)==3 and all(float(r["fraction0"])>0.97 for r in l4),"L4 exact-key zero-match exposure above 97 percent")

lin=js("results/lineage_family/decision.json")
check(lin["memberCount"]==76 and lin["permutations"]==999,"76-member 999-permutation lineage family")
check(close(lin["combinedCorrectedP"],0.553),"lineage family p=0.553")
check(lin["status"]=="FINITE_POSITIVE_CONFIGURATION_WITHOUT_FAMILY_WISE_SUPPORT","bounded lineage decision")
check(lin["effect"]["positive_fold_count"]==3 and lin["winnerMaterial"] is False,"finite 3/3 lineage result remains sub-material")

diag=js("results/diagnostics/global_null_decision.json")
check(diag["contrasts"]==100 and diag["correctedDepartures"]==0,"100 contrasts and zero corrected departures")
check(close(diag["conditionalFamilyP"],0.31007751937984496),"discrete family p=0.31008")

qual=rows("results/diagnostics/qualification_sensitivity.csv")
check(len(qual)==405,"405 qualification combinations")
check(sum(r["passesDiagnosticGates"]=="True" for r in qual)==45,"45 numerical gate combinations pass")
check(all(r["certifiesBiologicalIdentity"]=="False" for r in qual),"zero qualification rows certify biological identity")

formal=js("formal/verification_receipt.json")
check(formal["status"]=="PASS","formal receipt passed")
check(formal["biologicalExchangeabilityCertified"] is False,"formal layer does not certify exchangeability")

print(f"PASS: {len(checks)} Paper 2 public verification checks")
for c in checks:
    print("  -",c)
