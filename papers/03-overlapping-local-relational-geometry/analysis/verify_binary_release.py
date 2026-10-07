#!/usr/bin/env python3
"""Strict ordinary-file validation of every released Final08 binary."""
import argparse
import hashlib
import json
import re
from pathlib import Path

FIGURES = (
    "RMMO_Paper3_Figure1_StrictConstructionConsensus.pdf",
    "RMMO_Paper3_Figure2_PrimaryMolecularResults.pdf",
    "RMMO_Paper3_Figure3_FiberSpatialBiologicalGeometry.pdf",
    "RMMO_Paper3_Figure4_CompetitiveBaselines.pdf",
    "RMMO_Paper3_Figure5_ConstructionRepresentationSensitivity.pdf",
    "RMMO_Paper3_Figure6_RetentionGlobalBoundary.pdf",
    "RMMO_Paper3_FigureS1_ConsensusVoteHistogram.pdf",
    "RMMO_Paper3_FigureS2_RetentionTrajectories.pdf",
    "RMMO_Paper3_FigureS3_DomainScaleSensitivity.pdf",
    "RMMO_Paper3_FigureS4_RetentionEffectSurface.pdf",
)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1048576), b""):
            h.update(part)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    root = ap.parse_args().root.resolve()
    mp = root / "verification" / "binary_authority.json"
    if not mp.is_file():
        raise SystemExit("INCOMPLETE: binary_authority.json is not published")
    authority = json.loads(mp.read_text())
    figures, banks = authority.get("figures", {}), authority.get("banks", {})
    if set(figures) != set(FIGURES) or not isinstance(banks, dict) or not banks:
        raise SystemExit("INCOMPLETE: ten figure PDFs and a nonempty bank manifest are required")
    records = authority.get("files", {})
    exports = json.loads((root / "verification/source_export_manifest.json").read_text())["files"]
    expected = {rel: rec for rel, rec in exports.items()
                if Path(rel).suffix in (".pdf", ".png", ".npz", ".gz")}
    if records != expected:
        raise SystemExit("Binary inventory differs from the frozen export authority")
    if figures != {Path(rel).name: rec for rel, rec in records.items()
                   if rel.startswith("figures/") and rel.endswith(".pdf")}:
        raise SystemExit("Figure inventory differs from the frozen export authority")
    for rel, rec in banks.items():
        if not rel.startswith("reproducibility/") or not rel.endswith(".npz") or rec != records.get(rel):
            raise SystemExit("Invalid bank path: " + rel)
    if banks != {rel: rec for rel, rec in records.items() if rel.endswith(".npz")}:
        raise SystemExit("Incomplete NPZ inventory")
    for rel, rec in records.items():
        path = root / rel
        if Path(rel).is_absolute() or ".." in Path(rel).parts or not path.resolve().is_relative_to(root) or not path.is_file() or path.is_symlink():
            raise SystemExit("Missing or unsafe artifact: " + rel)
        start = path.read_bytes()[:160]
        if start.startswith(b"version https://git-lfs.github.com/spec/v1"):
            raise SystemExit("Git LFS pointer instead of actual artifact: " + rel)
        magic = {".pdf": b"%PDF-", ".png": b"\x89PNG\r\n\x1a\n", ".npz": b"PK\x03\x04", ".gz": b"\x1f\x8b"}
        if path.suffix not in magic or not start.startswith(magic[path.suffix]):
            raise SystemExit("Invalid artifact format: " + rel)
        if not isinstance(rec.get("bytes"), int) or path.stat().st_size != rec["bytes"]:
            raise SystemExit("Byte count mismatch: " + rel)
        if not re.fullmatch(r"[a-f0-9]{64}", rec.get("sha256", "")):
            raise SystemExit("Expected SHA-256 missing: " + rel)
        if sha256(path) != rec["sha256"]:
            raise SystemExit("SHA-256 mismatch: " + rel)
    print(f"PASS: {len(FIGURES)} verified figure PDFs; {len(records)} binaries including {len(banks)} NPZ inputs/banks")
    print("NOTE: numerical replay remains a separate scientific gate")

if __name__ == "__main__":
    main()
