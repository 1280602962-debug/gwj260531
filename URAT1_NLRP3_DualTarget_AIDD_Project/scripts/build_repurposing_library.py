#!/usr/bin/env python3
"""
Build URAT1/NLRP3 repurposing library from ChEMBL Excel exports.

Expected raw inputs (user download from ChEMBL Explore Drugs):
  - Phase1_2_3.xls(x)   clinical-phase export
  - Level1 ATC.xls(x)   Level-1 ATC filtered export
  - Level 2 ATC.xls(x)  Level-2 ATC filtered export

Outputs:
  data/repurposing/repurposing_manifest.csv      # merged, deduplicated, tagged
  data/repurposing/repurposing_primary.csv         # recommended screening panel
  data/repurposing/repurposing_build_summary.json

Example:
  python3 scripts/build_repurposing_library.py \\
    --input-dir data/repurposing/raw \\
    --phase-file "Phase1_2_3.xls" \\
    --atc-l1-file "Level1 ATC.xls" \\
    --atc-l2-file "Level 2 ATC.xls"
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from utils_ml import canonicalize, murcko_scaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS_CSV = PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks.csv"
DEFAULT_OUT = PROJECT_ROOT / "data" / "repurposing"

# ChEMBL column aliases (case-insensitive match)
COL_ALIASES = {
    "chembl_id": [
        "molecule chembl id",
        "chembl id",
        "chembl_id",
        "molecule_chembl_id",
    ],
    "smiles": [
        "smiles",
        "canonical smiles",
        "canonical_smiles",
        "structure",
    ],
    "pref_name": [
        "pref name",
        "pref_name",
        "molecule name",
        "name",
    ],
    "max_phase": [
        "max phase",
        "max_phase",
        "maximum clinical phase",
    ],
    "mol_type": [
        "molecule type",
        "molecule_type",
        "type",
    ],
    "atc": [
        "atc classification",
        "atc codes",
        "atc",
        "atc code",
    ],
    "mw": [
        "mw",
        "molecular weight",
        "molecular_weight",
        "full mol wt",
    ],
}

MANIFEST_COLS = [
    "canonical_smiles",
    "scaffold",
    "inchikey_block1",
    "molecule_chembl_id",
    "pref_name",
    "max_phase",
    "mw",
    "mol_type",
    "atc_codes",
    "source_phase",
    "source_atc_l1",
    "source_atc_l2",
    "library_panel",
    "benchmark_ref",
    "include_screen",
]


def _norm_col(name: str) -> str:
    return re.sub(r"\s+", " ", str(name).strip().lower())


def _resolve_col(columns: list[str], key: str) -> str | None:
    norm_map = {_norm_col(c): c for c in columns}
    for alias in COL_ALIASES[key]:
        if alias in norm_map:
            return norm_map[alias]
    return None


def read_chembl_excel(path: Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.suffix.lower()
    if suffix in {".xls", ".xlsx", ".xlsm"}:
        try:
            return pd.read_excel(path, engine="xlrd" if suffix == ".xls" else "openpyxl")
        except ImportError as exc:
            hint = "pip install openpyxl" if suffix != ".xls" else "pip install xlrd"
            raise ImportError(f"Cannot read {path.name}: {exc}. Try: {hint}") from exc
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported format: {path}")


def desalt_largest_fragment(smiles: str) -> str | None:
    from rdkit import Chem
    from rdkit.Chem.SaltRemover import SaltRemover

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    remover = SaltRemover()
    mol = remover.StripMol(mol, dontRemoveEverything=True)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        return None
    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    try:
        Chem.SanitizeMol(largest)
        return Chem.MolToSmiles(largest)
    except Exception:
        return None


def inchikey_block1(smiles: str) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol).split("-")[0]
    except Exception:
        return None


def mol_weight(smiles: str) -> float | None:
    from rdkit.Chem import Descriptors

    mol = __import__("rdkit").Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return float(Descriptors.MolWt(mol))


def normalize_export(df: pd.DataFrame, source_tag: str) -> pd.DataFrame:
    cols = list(df.columns)
    smi_col = _resolve_col(cols, "smiles")
    if smi_col is None:
        raise ValueError(f"No SMILES column in export. Columns: {cols[:20]}")

    out = pd.DataFrame()
    out["smiles_raw"] = df[smi_col].astype(str)
    cid_col = _resolve_col(cols, "chembl_id")
    out["molecule_chembl_id"] = df[cid_col].astype(str) if cid_col else pd.NA
    pn_col = _resolve_col(cols, "pref_name")
    out["pref_name"] = df[pn_col].astype(str) if pn_col else pd.NA
    mp_col = _resolve_col(cols, "max_phase")
    out["max_phase"] = pd.to_numeric(df[mp_col], errors="coerce") if mp_col else pd.NA
    mt_col = _resolve_col(cols, "mol_type")
    out["mol_type"] = df[mt_col].astype(str) if mt_col else None
    atc_col = _resolve_col(cols, "atc")
    out["atc_codes"] = df[atc_col].astype(str) if atc_col else None
    mw_col = _resolve_col(cols, "mw")
    out["mw_raw"] = pd.to_numeric(df[mw_col], errors="coerce") if mw_col else pd.NA
    out["source_file"] = source_tag
    return out


def standardize_rows(df: pd.DataFrame, desalt: bool, mw_min: float, mw_max: float) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        raw = str(r["smiles_raw"]).strip()
        if not raw or raw.lower() in {"nan", "none"}:
            continue
        smi = desalt_largest_fragment(raw) if desalt else canonicalize(raw)
        if not smi:
            smi = canonicalize(raw)
        if not smi:
            continue
        mw = r["mw_raw"] if pd.notna(r["mw_raw"]) else mol_weight(smi)
        if mw is not None and (mw < mw_min or mw > mw_max):
            continue
        mt = str(r.get("mol_type") or "").lower()
        if mt and "small molecule" not in mt and mt not in {"nan", "none"}:
            if any(x in mt for x in ("protein", "peptide", "oligo", "antibody")):
                continue
        rows.append(
            {
                "canonical_smiles": smi,
                "scaffold": murcko_scaffold(smi),
                "inchikey_block1": inchikey_block1(smi),
                "molecule_chembl_id": r.get("molecule_chembl_id"),
                "pref_name": r.get("pref_name"),
                "max_phase": r.get("max_phase"),
                "mw": mw,
                "mol_type": r.get("mol_type"),
                "atc_codes": r.get("atc_codes"),
                "source_file": r["source_file"],
            }
        )
    return pd.DataFrame(rows)


def load_benchmarks() -> pd.DataFrame:
    if not BENCHMARKS_CSV.exists():
        return pd.DataFrame(columns=["canonical_smiles", "benchmark_ref", "compound_name"])
    bench = pd.read_csv(BENCHMARKS_CSV, low_memory=False)
    bench = bench.dropna(subset=["canonical_smiles"]).drop_duplicates("canonical_smiles")
    return bench[["canonical_smiles", "compound_id", "compound_name"]].rename(
        columns={"compound_id": "benchmark_ref", "compound_name": "pref_name"}
    )


def merge_sources(
    phase_df: pd.DataFrame,
    atc_l1_df: pd.DataFrame,
    atc_l2_df: pd.DataFrame,
    min_phase_primary: float,
) -> pd.DataFrame:
    combined = pd.concat([phase_df, atc_l1_df, atc_l2_df], ignore_index=True)
    if combined.empty:
        return combined

    combined = combined[combined["canonical_smiles"].notna()].copy()
    missing_ik = combined["inchikey_block1"].isna()
    if missing_ik.any():
        combined.loc[missing_ik, "inchikey_block1"] = combined.loc[missing_ik, "canonical_smiles"]

    def _agg_phase(s: pd.Series) -> bool:
        return (s == "phase").any()

    def _agg_l1(s: pd.Series) -> bool:
        return (s == "atc_l1").any()

    def _agg_l2(s: pd.Series) -> bool:
        return (s == "atc_l2").any()

    grouped = combined.groupby("inchikey_block1", as_index=False).agg(
        canonical_smiles=("canonical_smiles", "first"),
        scaffold=("scaffold", "first"),
        inchikey_block1=("inchikey_block1", "first"),
        molecule_chembl_id=("molecule_chembl_id", "first"),
        pref_name=("pref_name", "first"),
        max_phase=("max_phase", "max"),
        mw=("mw", "first"),
        mol_type=("mol_type", "first"),
        atc_codes=("atc_codes", lambda x: "; ".join(sorted({str(v) for v in x if str(v) not in {"nan", "None", ""}}))),
        source_phase=("source_file", _agg_phase),
        source_atc_l1=("source_file", _agg_l1),
        source_atc_l2=("source_file", _agg_l2),
    )

    def _panel(row) -> str:
        atc_hit = row["source_atc_l1"] or row["source_atc_l2"]
        phase_ok = pd.notna(row["max_phase"]) and float(row["max_phase"]) >= min_phase_primary
        if atc_hit and phase_ok:
            return "primary_atc_phase"
        if atc_hit:
            return "atc_only"
        if row["source_phase"]:
            return "phase_only"
        return "other"

    grouped["library_panel"] = grouped.apply(_panel, axis=1)
    grouped["benchmark_ref"] = None
    grouped["include_screen"] = grouped["library_panel"].isin(
        ["primary_atc_phase", "atc_only", "phase_only"]
    )
    return grouped[MANIFEST_COLS]


def force_add_benchmarks(manifest: pd.DataFrame, bench: pd.DataFrame) -> pd.DataFrame:
    if bench.empty:
        return manifest
    existing = set(manifest["inchikey_block1"].dropna())
    add_rows = []
    for _, b in bench.iterrows():
        smi = canonicalize(b["canonical_smiles"])
        if not smi:
            continue
        ik = inchikey_block1(smi)
        if ik in existing:
            idx = manifest.index[manifest["inchikey_block1"] == ik]
            if len(idx):
                manifest.loc[idx, "benchmark_ref"] = b["benchmark_ref"]
            continue
        add_rows.append(
            {
                "canonical_smiles": smi,
                "scaffold": murcko_scaffold(smi),
                "inchikey_block1": ik,
                "molecule_chembl_id": None,
                "pref_name": b.get("pref_name"),
                "max_phase": None,
                "mw": mol_weight(smi),
                "mol_type": "benchmark",
                "atc_codes": None,
                "source_phase": False,
                "source_atc_l1": False,
                "source_atc_l2": False,
                "library_panel": "benchmark_forced",
                "benchmark_ref": b["benchmark_ref"],
                "include_screen": True,
            }
        )
    if add_rows:
        manifest = pd.concat([manifest, pd.DataFrame(add_rows)], ignore_index=True)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge ChEMBL repurposing Excel exports")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_OUT / "raw")
    parser.add_argument("--phase-file", type=str, default="Phase1_2_3.xls")
    parser.add_argument("--atc-l1-file", type=str, default="Level1 ATC.xls")
    parser.add_argument("--atc-l2-file", type=str, default="Level 2 ATC.xls")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--desalt", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--mw-min", type=float, default=150.0)
    parser.add_argument("--mw-max", type=float, default=800.0)
    parser.add_argument(
        "--min-phase-primary",
        type=float,
        default=3.0,
        help="max_phase threshold for primary_atc_phase panel (3=Phase III+approved)",
    )
    parser.add_argument(
        "--primary-mode",
        choices=["atc_phase", "atc_only", "phase_only", "union"],
        default="atc_phase",
        help="Which panel to write repurposing_primary.csv",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def _load(tag: str, fname: str) -> pd.DataFrame:
        path = args.input_dir / fname
        if not path.exists():
            # try alternate extensions
            for ext in (".xlsx", ".csv"):
                alt = path.with_suffix(ext)
                if alt.exists():
                    path = alt
                    break
        raw = read_chembl_excel(path)
        norm = normalize_export(raw, tag)
        return standardize_rows(norm, desalt=args.desalt, mw_min=args.mw_min, mw_max=args.mw_max)

    phase_df = _load("phase", args.phase_file)
    atc_l1_df = _load("atc_l1", args.atc_l1_file)
    atc_l2_df = _load("atc_l2", args.atc_l2_file)

    manifest = merge_sources(phase_df, atc_l1_df, atc_l2_df, args.min_phase_primary)
    manifest = force_add_benchmarks(manifest, load_benchmarks())

    manifest_path = args.output_dir / "repurposing_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    panel_map = {
        "atc_phase": manifest["library_panel"] == "primary_atc_phase",
        "atc_only": manifest["source_atc_l1"] | manifest["source_atc_l2"],
        "phase_only": manifest["source_phase"],
        "union": manifest["include_screen"],
    }
    primary = manifest[panel_map[args.primary_mode]].copy()
    primary_path = args.output_dir / "repurposing_primary.csv"
    primary.to_csv(primary_path, index=False)

    summary = {
        "input_files": {
            "phase": str(args.input_dir / args.phase_file),
            "atc_l1": str(args.input_dir / args.atc_l1_file),
            "atc_l2": str(args.input_dir / args.atc_l2_file),
        },
        "n_raw_after_clean": {
            "phase": int(len(phase_df)),
            "atc_l1": int(len(atc_l1_df)),
            "atc_l2": int(len(atc_l2_df)),
        },
        "n_manifest_unique_inchikey": int(len(manifest)),
        "n_by_library_panel": manifest["library_panel"].value_counts().to_dict(),
        "n_primary_export": int(len(primary)),
        "primary_mode": args.primary_mode,
        "filters": {
            "desalt": args.desalt,
            "mw_min": args.mw_min,
            "mw_max": args.mw_max,
            "min_phase_primary": args.min_phase_primary,
        },
        "outputs": {
            "manifest": str(manifest_path),
            "primary": str(primary_path),
        },
    }
    summary_path = args.output_dir / "repurposing_build_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"\nWrote manifest ({len(manifest)} compounds) -> {manifest_path}")
    print(f"Wrote primary panel ({len(primary)} compounds, mode={args.primary_mode}) -> {primary_path}")


if __name__ == "__main__":
    main()
