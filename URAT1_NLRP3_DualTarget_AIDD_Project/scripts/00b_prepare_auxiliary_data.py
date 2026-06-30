#!/usr/bin/env python3
"""
Curate OAT1/OAT3 (transfer) and OCT1/OCT2 (detargeting) ChEMBL cf12 exports.

Input (default): data/raw/auxiliary/*_chembl_cf12.csv
Output:
  data/auxiliary/oat1_chembl_curated.csv
  data/auxiliary/oat3_chembl_curated.csv
  data/auxiliary/oct1_chembl_curated.csv
  data/auxiliary/oct2_chembl_curated.csv
  data/auxiliary/oat_combined_transfer.csv
  data/auxiliary/auxiliary_data_summary.json
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from utils_ml import canonicalize, curate_urat1_raw

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "targets.yaml"
DEFAULT_UPLOADS = Path("/home/ubuntu/.cursor/projects/workspace/uploads")

# ChEMBL may return legacy target IDs for OCT; both map to human SLC22A1/A2.
OCT_TARGET_ID_ALIASES = {
    "OCT1": {"CHEMBL2073664", "CHEMBL5685"},
    "OCT2": {"CHEMBL1770032", "CHEMBL1743122"},
}


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


UPLOAD_NAME_MAP = {
    "OAT1_chembl_cf12.csv": ["OAT1_chembl_cf12.csv", "OAT1_0f82.csv", "OAT1.csv"],
    "OAT3_chembl_cf12.csv": ["OAT3_chembl_cf12.csv", "OAT3_74e9.csv", "OAT3.csv"],
    "OCT1_chembl_cf12.csv": ["OCT1_chembl_cf12.csv", "OCT1_54fb.csv", "OCT1.csv"],
    "OCT2_chembl_cf12.csv": ["OCT2_chembl_cf12.csv", "OCT2_eea2.csv", "OCT2.csv"],
}


def resolve_raw(path: Path | None, default_name: str) -> Path:
    if path and path.exists():
        return path
    local = PROJECT_ROOT / "data" / "raw" / "auxiliary" / default_name
    if local.exists():
        return local
    for candidate in UPLOAD_NAME_MAP.get(default_name, [default_name]):
        upload = DEFAULT_UPLOADS / candidate
        if upload.exists():
            return upload
    raise FileNotFoundError(f"Missing {default_name}; pass explicit --csv path")


def filter_target(df: pd.DataFrame, target_ids: set[str]) -> pd.DataFrame:
    return df[df["Target ChEMBL ID"].astype(str).isin(target_ids)].copy()


def curate_auxiliary(
    path: Path,
    *,
    source_target: str,
    target_ids: set[str],
    pactivity_min: float,
    pactivity_max: float,
    max_std: float,
    max_range: float,
) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    raw = filter_target(raw, target_ids)
    tmp = path.parent / f"_tmp_{source_target}.csv"
    raw.to_csv(tmp, index=False)
    out = curate_urat1_raw(
        tmp,
        pactivity_min=pactivity_min,
        pactivity_max=pactivity_max,
        max_std=max_std,
        max_range=max_range,
    )
    tmp.unlink(missing_ok=True)
    out["source_target"] = source_target
    out["source_chembl_id"] = raw["Target ChEMBL ID"].astype(str).mode().iloc[0]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Curate SLC22 auxiliary ChEMBL libraries")
    parser.add_argument("--oat1-csv", type=Path, default=None)
    parser.add_argument("--oat3-csv", type=Path, default=None)
    parser.add_argument("--oct1-csv", type=Path, default=None)
    parser.add_argument("--oct2-csv", type=Path, default=None)
    parser.add_argument("--copy-raw", action="store_true")
    args = parser.parse_args()

    config = load_config()["data_curation"]
    out_dir = PROJECT_ROOT / "data" / "auxiliary"
    raw_dir = PROJECT_ROOT / "data" / "raw" / "auxiliary"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    specs = {
        "OAT1": (args.oat1_csv, "OAT1_chembl_cf12.csv", {"CHEMBL1641347"}),
        "OAT3": (args.oat3_csv, "OAT3_chembl_cf12.csv", {"CHEMBL1641348"}),
        "OCT1": (args.oct1_csv, "OCT1_chembl_cf12.csv", OCT_TARGET_ID_ALIASES["OCT1"]),
        "OCT2": (args.oct2_csv, "OCT2_chembl_cf12.csv", OCT_TARGET_ID_ALIASES["OCT2"]),
    }

    curated: dict[str, pd.DataFrame] = {}
    raw_stats: dict[str, dict] = {}

    for name, (arg_path, default_name, target_ids) in specs.items():
        src = resolve_raw(arg_path, default_name)
        if args.copy_raw:
            shutil.copy2(src, raw_dir / default_name)
        raw = pd.read_csv(src, low_memory=False)
        raw_stats[name] = {
            "source_file": str(src),
            "n_rows": int(len(raw)),
            "target_ids_in_file": sorted(raw["Target ChEMBL ID"].astype(str).unique().tolist()),
            "ic50_equals_smiles": int(
                raw[
                    (raw["Standard Type"] == "IC50")
                    & (raw["Standard Relation"].astype(str).str.strip("'\"") == "=")
                ]["Smiles"].nunique()
            ),
        }
        curated[name] = curate_auxiliary(
            src,
            source_target=name,
            target_ids=target_ids,
            pactivity_min=config["pactivity_range"][0],
            pactivity_max=config["pactivity_range"][1],
            max_std=config["conflict_std_threshold"],
            max_range=config["conflict_range_threshold"],
        )
        out_path = out_dir / f"{name.lower()}_chembl_curated.csv"
        curated[name].to_csv(out_path, index=False)
        print(f"{name}: {len(curated[name])} curated -> {out_path}")

    oat_combined = pd.concat([curated["OAT1"], curated["OAT3"]], ignore_index=True)
    oat_combined = oat_combined.drop_duplicates("canonical_smiles").reset_index(drop=True)
    oat_combined.to_csv(out_dir / "oat_combined_transfer.csv", index=False)
    print(f"OAT combined: {len(oat_combined)} unique SMILES")

    urat1_path = PROJECT_ROOT / "data" / "processed" / "urat1_curated.csv"
    urat1_smiles: set[str] = set()
    if urat1_path.exists():
        urat1_smiles = set(pd.read_csv(urat1_path)["canonical_smiles"])

    summary = {
        "curated": {
            k: {
                "n_compounds": int(len(v)),
                "n_scaffolds": int(v["scaffold"].nunique()),
                "pactivity_mean": float(v["pActivity"].mean()),
                "pactivity_std": float(v["pActivity"].std()),
                "n_pactivity_ge_6": int((v["pActivity"] >= 6).sum()),
                "source_chembl_id": v["source_chembl_id"].iloc[0],
                "overlap_urat1_smiles": int(len(set(v["canonical_smiles"]) & urat1_smiles)),
            }
            for k, v in curated.items()
        },
        "oat_combined_n": int(len(oat_combined)),
        "oat_combined_overlap_urat1": int(len(set(oat_combined["canonical_smiles"]) & urat1_smiles)),
        "raw": raw_stats,
        "notes": [
            "IC50-only ChEMBL exports are often <100 compounds after curation; "
            "for transfer learning re-export with IC50+Ki+EC50 and broader assay filters.",
            "OCT exports may use legacy ChEMBL target IDs CHEMBL5685 / CHEMBL1743122 "
            "(same SLC22A1/A2 as CHEMBL2073664 / CHEMBL1770032).",
        ],
    }
    summary_path = out_dir / "auxiliary_data_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
