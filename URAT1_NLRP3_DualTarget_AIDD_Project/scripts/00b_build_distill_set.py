#!/usr/bin/env python3
"""Build distill subsets A/B/C/E (URAT1 retrospective only; not the production funnel).

Subsets:
  A  URAT1 training actives with pIC50 labels (~822)
  B  Murcko scaffold cluster centers from URAT1
  C  ChEMBL SLC22 neighborhood (URAT1 + OAT1/OAT3)
  D  Unlabeled diversity negatives (pre-built via sample_distill_subset_d.py)
  E  Benchmark boundary analogs

Not the TrueDecoy protocol-selection pool and not the clinical-library funnel.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from utils_ml import canonicalize, featurize_smiles, murcko_scaffold, scaffold_cv_indices

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "data" / "processed"
AUXILIARY = PROJECT_ROOT / "data" / "auxiliary"
BENCHMARKS = PROJECT_ROOT / "data" / "benchmarks"
DISTILL_DIR = PROJECT_ROOT / "data" / "distill"
SPLITS_DIR = PROJECT_ROOT / "data" / "splits"

MANIFEST_COLS = [
    "canonical_smiles",
    "scaffold",
    "subset",
    "label_type",
    "source_library",
    "has_bioactivity_label",
    "pActivity",
    "target_arm",
    "benchmark_ref",
    "molecule_chembl_id",
]

SUBSET_PRIORITY = {"A": 0, "E": 1, "B": 2, "C": 3, "D": 4}


def _row(
    smi: str,
    subset: str,
    label_type: str,
    source: str,
    has_label: bool,
    pactivity: float | None = None,
    target_arm: str | None = None,
    benchmark_ref: str | None = None,
    chembl_id: str | None = None,
) -> dict:
    return {
        "canonical_smiles": smi,
        "scaffold": murcko_scaffold(smi),
        "subset": subset,
        "label_type": label_type,
        "source_library": source,
        "has_bioactivity_label": has_label,
        "pActivity": pactivity,
        "target_arm": target_arm,
        "benchmark_ref": benchmark_ref,
        "molecule_chembl_id": chembl_id,
    }


def build_subset_a(urat1: pd.DataFrame) -> pd.DataFrame:
    """URAT1 curated training actives with pIC50."""
    rows = []
    for _, r in urat1.iterrows():
        rows.append(
            _row(
                r["canonical_smiles"],
                "A",
                "urat1_pactivity",
                "urat1_curated.csv",
                True,
                float(r["pActivity"]),
                "URAT1",
                chembl_id=str(r.get("molecule_chembl_id", "")) or None,
            )
        )
    return pd.DataFrame(rows)


def build_subset_b(urat1: pd.DataFrame, n_target: int, seed: int) -> pd.DataFrame:
    """One Murcko scaffold representative per cluster (highest pActivity medoid)."""
    rows = []
    for scaf, g in urat1.groupby("scaffold"):
        best = g.loc[g["pActivity"].idxmax()]
        rows.append(
            _row(
                best["canonical_smiles"],
                "B",
                "scaffold_representative",
                "urat1_curated.csv",
                True,
                float(best["pActivity"]),
                "URAT1",
                chembl_id=str(best.get("molecule_chembl_id", "")) or None,
            )
        )
    df = pd.DataFrame(rows)
    if len(df) > n_target:
        df = df.sample(n=n_target, random_state=seed).reset_index(drop=True)
    return df


def _fps_pick(smiles: list[str], n_pick: int, seed: int) -> list[int]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    valid = [(i, m) for i, m in enumerate(mols) if m is not None]
    if not valid:
        return []
    idx_map, mols_v = zip(*valid)
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols_v]

    n_pick = min(n_pick, len(fps))

    def distfunc(i: int, j: int) -> float:
        return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

    picker = MaxMinPicker()
    picked = list(picker.LazyPick(distfunc, len(fps), n_pick, firstPicks=[], seed=seed))
    return [idx_map[i] for i in picked]


def build_subset_c(
    urat1: pd.DataFrame,
    oat: pd.DataFrame,
    n_target: int,
    seed: int,
) -> pd.DataFrame:
    """FPS diversity pick from SLC22 neighborhood (URAT1 + OAT1/OAT3)."""
    pool_rows: list[dict] = []
    seen: set[str] = set()

    for _, r in urat1.iterrows():
        smi = r["canonical_smiles"]
        if smi in seen:
            continue
        seen.add(smi)
        pool_rows.append(
            {
                "canonical_smiles": smi,
                "pActivity": float(r["pActivity"]),
                "target_arm": "URAT1",
                "source_library": "urat1_curated.csv",
                "molecule_chembl_id": str(r.get("molecule_chembl_id", "")) or None,
            }
        )

    for _, r in oat.iterrows():
        smi = r["canonical_smiles"]
        if smi in seen:
            continue
        seen.add(smi)
        pool_rows.append(
            {
                "canonical_smiles": smi,
                "pActivity": float(r["pActivity"]),
                "target_arm": r.get("source_target", "OAT"),
                "source_library": "oat_combined_transfer.csv",
                "molecule_chembl_id": str(r.get("molecule_chembl_id", "")) or None,
            }
        )

    pool = pd.DataFrame(pool_rows)
    smiles = pool["canonical_smiles"].tolist()
    picked_idx = _fps_pick(smiles, n_target, seed)

    rows = []
    for i in picked_idx:
        r = pool.iloc[i]
        rows.append(
            _row(
                r["canonical_smiles"],
                "C",
                "slc22_fps_neighbor",
                r["source_library"],
                True,
                float(r["pActivity"]),
                str(r["target_arm"]),
                chembl_id=r["molecule_chembl_id"],
            )
        )
    return pd.DataFrame(rows)


def load_subset_d(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Subset D not found: {path}. Run scripts/sample_distill_subset_d.py first."
        )
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        smi = r["canonical_smiles"]
        rows.append(
            _row(
                smi,
                "D",
                r.get("label_type", "unlabeled_negative"),
                str(r.get("source_library", "external_library")),
                bool(r.get("has_bioactivity_label", False)),
            )
        )
    return pd.DataFrame(rows)


def build_subset_e(
    benchmarks: pd.DataFrame,
    urat1: pd.DataFrame,
    nlrp3: pd.DataFrame,
    n_target: int,
    min_tanimoto: float,
    seed: int,
) -> pd.DataFrame:
    """Benchmark references + scaffold-novel analogs from project libraries."""
    from rdkit import DataStructs
    from rdkit.Chem import AllChem

    bench = benchmarks.drop_duplicates("canonical_smiles").copy()
    bench["canonical_smiles"] = bench["canonical_smiles"].map(
        lambda s: canonicalize(s) if pd.notna(s) else None
    )
    bench = bench[bench["canonical_smiles"].notna()]

    rows: list[dict] = []
    seen: set[str] = set()

    # Reference benchmark compounds (may be absent from training CSV)
    for _, r in bench.iterrows():
        smi = r["canonical_smiles"]
        if smi in seen:
            continue
        seen.add(smi)
        role = str(r.get("validation_role", ""))
        arm = "URAT1" if "URAT1" in role or r.get("target_gene") == "SLC22A12" else "NLRP3"
        pact = r.get("pactivity")
        pact_f = float(pact) if pd.notna(pact) else None
        rows.append(
            _row(
                smi,
                "E",
                "benchmark_reference",
                "literature_benchmarks.csv",
                pact_f is not None,
                pact_f,
                arm,
                str(r.get("compound_name", r.get("compound_id", ""))),
                str(r.get("chembl_id", "")) or None,
            )
        )

    # Analog neighbors from URAT1 + NLRP3 pools
    lib = pd.concat(
        [
            urat1[["canonical_smiles", "pActivity"]].assign(target_arm="URAT1"),
            nlrp3.groupby("canonical_smiles", as_index=False)
            .agg(pActivity=("pActivity", "median"), active=("active", "max"))
            .assign(
                pActivity=lambda d: d["pActivity"].where(d["pActivity"].notna(), np.where(d["active"] == 1, 6.5, 4.5)),
                target_arm="NLRP3",
            )[["canonical_smiles", "pActivity", "target_arm"]],
        ],
        ignore_index=True,
    ).drop_duplicates("canonical_smiles")

    lib_smiles = lib["canonical_smiles"].tolist()
    lib_fps = {
        s: AllChem.GetMorganFingerprintAsBitVect(AllChem.MolFromSmiles(s), 2, nBits=2048)
        for s in lib_smiles
        if AllChem.MolFromSmiles(s) is not None
    }

    analog_candidates: list[tuple[float, dict]] = []
    for _, br in bench.iterrows():
        bsmi = br["canonical_smiles"]
        bmol = AllChem.MolFromSmiles(bsmi)
        if bmol is None:
            continue
        bfp = AllChem.GetMorganFingerprintAsBitVect(bmol, 2, nBits=2048)
        bname = str(br.get("compound_name", br.get("compound_id", "")))
        for smi, lfp in lib_fps.items():
            if smi == bsmi or smi in seen:
                continue
            tc = DataStructs.TanimotoSimilarity(bfp, lfp)
            if tc < min_tanimoto:
                continue
            lr = lib[lib["canonical_smiles"] == smi].iloc[0]
            analog_candidates.append(
                (
                    tc,
                    _row(
                        smi,
                        "E",
                        "benchmark_analog",
                        "project_library_neighbor",
                        True,
                        float(lr["pActivity"]) if pd.notna(lr["pActivity"]) else None,
                        str(lr["target_arm"]),
                        bname,
                    ),
                )
            )

    analog_candidates.sort(key=lambda x: -x[0])
    for _, row in analog_candidates:
        if len(rows) >= n_target:
            break
        smi = row["canonical_smiles"]
        if smi in seen:
            continue
        seen.add(smi)
        rows.append(row)

    return pd.DataFrame(rows)


def merge_subsets(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """Deduplicate by canonical_smiles; higher-priority subset wins."""
    combined = pd.concat(parts, ignore_index=True)
    combined["_pri"] = combined["subset"].map(SUBSET_PRIORITY)
    combined = combined.sort_values("_pri").drop_duplicates("canonical_smiles", keep="first")
    return combined.drop(columns=["_pri"]).reset_index(drop=True)


def write_scaffold_folds(urat1: pd.DataFrame, n_splits: int, out_dir: Path) -> None:
    smiles = urat1["canonical_smiles"].tolist()
    out_dir.mkdir(parents=True, exist_ok=True)
    for fold_i, (tr_idx, te_idx) in enumerate(scaffold_cv_indices(smiles, n_splits=n_splits)):
        fold_df = pd.DataFrame(
            {
                "canonical_smiles": [smiles[i] for i in tr_idx] + [smiles[i] for i in te_idx],
                "fold": fold_i,
                "split": ["train"] * len(tr_idx) + ["test"] * len(te_idx),
            }
        )
        fold_df.to_csv(out_dir / f"scaffold_fold_{fold_i}.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="URAT1 distill set builder (A/B/C/E + merge D)")
    parser.add_argument("--output-dir", type=Path, default=DISTILL_DIR)
    parser.add_argument("--n-b", type=int, default=1000, help="Target size subset B (500–1000)")
    parser.add_argument("--n-c", type=int, default=5000, help="Target size subset C (2000–5000)")
    parser.add_argument("--n-e", type=int, default=200, help="Target size subset E (~200)")
    parser.add_argument("--min-tanimoto-e", type=float, default=0.35, help="Min Tanimoto for E analogs")
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-d", action="store_true", help="Do not reload subset D")
    parser.add_argument("--d-csv", type=Path, default=None, help="Override distill_subset_d.csv path")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    urat1 = pd.read_csv(PROCESSED / "urat1_curated.csv")
    nlrp3 = pd.read_csv(PROCESSED / "nlrp3_records.csv")
    oat = pd.read_csv(AUXILIARY / "oat_combined_transfer.csv")
    benchmarks = pd.read_csv(BENCHMARKS / "literature_benchmarks.csv")

    print("Building subset A (URAT1 training actives)...")
    sub_a = build_subset_a(urat1)
    sub_a.to_csv(args.output_dir / "distill_subset_a.csv", index=False)
    print(f"  A: {len(sub_a)} compounds")

    print("Building subset B (Murcko scaffold representatives)...")
    sub_b = build_subset_b(urat1, n_target=max(500, min(args.n_b, 1000)), seed=args.seed)
    sub_b.to_csv(args.output_dir / "distill_subset_b.csv", index=False)
    print(f"  B: {len(sub_b)} scaffolds (data has {urat1['scaffold'].nunique()} unique Murcko scaffolds)")

    print("Building subset C (SLC22 FPS neighborhood)...")
    sub_c = build_subset_c(urat1, oat, n_target=args.n_c, seed=args.seed)
    sub_c.to_csv(args.output_dir / "distill_subset_c.csv", index=False)
    pool_n = len(set(urat1["canonical_smiles"]) | set(oat["canonical_smiles"]))
    print(f"  C: {len(sub_c)} picked from SLC22 pool of {pool_n} unique compounds")

    print("Building subset E (benchmark boundary set)...")
    sub_e = build_subset_e(
        benchmarks, urat1, nlrp3, n_target=args.n_e, min_tanimoto=args.min_tanimoto_e, seed=args.seed
    )
    sub_e.to_csv(args.output_dir / "distill_subset_e.csv", index=False)
    print(f"  E: {len(sub_e)} compounds ({(sub_e['label_type']=='benchmark_reference').sum()} refs)")

    parts = [sub_a, sub_b, sub_c, sub_e]
    sub_d = None
    if not args.skip_d:
        d_path = args.d_csv or (args.output_dir / "distill_subset_d.csv")
        print(f"Loading subset D from {d_path}...")
        sub_d = load_subset_d(d_path)
        sub_d.to_csv(args.output_dir / "distill_subset_d.csv", index=False)
        parts.append(sub_d)
        print(f"  D: {len(sub_d)} unlabeled negatives")

    manifest = merge_subsets(parts)
    manifest = manifest[MANIFEST_COLS]
    manifest.to_csv(args.output_dir / "distill_manifest.csv", index=False)
    print(f"\nMerged manifest: {len(manifest)} unique compounds")

    counts = manifest["subset"].value_counts().to_dict()
    summary = {
        "framework": "MASFL_v3.1",
        "stage": "0.3_distill_set",
        "n_total_unique": int(len(manifest)),
        "n_by_subset_manifest": {k: int(counts.get(k, 0)) for k in ["A", "B", "C", "D", "E"]},
        "n_by_subset_file": {
            "A": int(len(sub_a)),
            "B": int(len(sub_b)),
            "C": int(len(sub_c)),
            "D": int(len(sub_d)) if not args.skip_d else 0,
            "E": int(len(sub_e)),
        },
        "n_labeled": int(manifest["has_bioactivity_label"].sum()),
        "n_unlabeled": int((~manifest["has_bioactivity_label"]).sum()),
        "targets_masfl": {
            "A": "~650/fold URAT1 train (all curated used)",
            "B": "500–1000 scaffold reps",
            "C": "2000–5000 SLC22 FPS",
            "D": "3000–10000 unlabeled",
            "E": "~200 benchmark boundary",
            "total": "5k–20k",
        },
        "data_limits": {
            "urat1_scaffolds": int(urat1["scaffold"].nunique()),
            "slc22_pool_unique": pool_n,
            "note": "Subset C below 2k target because ChEMBL SLC22 exports are small (~822 URAT1 + 73 OAT).",
        },
        "seed": args.seed,
        "files": {
            "manifest": str(args.output_dir / "distill_manifest.csv"),
            "subsets": [f"distill_subset_{s}.csv" for s in "abcde"],
        },
    }
    with open(args.output_dir / "distill_set_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Writing URAT1 scaffold CV fold maps...")
    write_scaffold_folds(urat1, args.n_splits, SPLITS_DIR)

    print(f"\n=== Distill set summary ===")
    for k in "ABCDE":
        file_n = summary["n_by_subset_file"].get(k, 0)
        man_n = summary["n_by_subset_manifest"].get(k, 0)
        print(f"  Subset {k}: {file_n} in file, {man_n} primary in manifest")
    print(f"  Total unique: {summary['n_total_unique']}")
    print(f"  Labeled / unlabeled: {summary['n_labeled']} / {summary['n_unlabeled']}")
    print(f"  Summary -> {args.output_dir / 'distill_set_summary.json'}")


if __name__ == "__main__":
    main()
