#!/usr/bin/env python3
"""
Collect per-molecule top poses/scores from Vina/gnina/RTM outputs.
Uses work/mol_index_map.csv for SMILES/role join keys (LigPrep order).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from rdkit import Chem


def parse_vina_log_or_pdbqt(pdbqt: Path, log: Path | None) -> list[dict]:
    scores = []
    if log and log.exists():
        for line in log.read_text(errors="ignore").splitlines():
            m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)", line)
            if m and int(m.group(1)) <= 20:
                scores.append({"mode": int(m.group(1)), "vina_affinity": float(m.group(2))})
        if scores:
            return scores
    if pdbqt.exists():
        mode = 0
        for line in pdbqt.read_text(errors="ignore").splitlines():
            if "REMARK VINA RESULT:" in line:
                mode += 1
                scores.append({"mode": mode, "vina_affinity": float(line.split()[3])})
    return scores


def parse_gnina_sdf(sdf: Path) -> list[dict]:
    out = []
    if not sdf.exists():
        return out
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    for i, mol in enumerate(suppl, start=1):
        if mol is None:
            continue
        row = {"mode": i}
        for k in ("minimizedAffinity", "CNNscore", "CNNaffinity"):
            if mol.HasProp(k):
                row[k] = float(mol.GetProp(k))
        out.append(row)
    return out


def load_rtm_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    col_score = "score" if "score" in df.columns else df.columns[-1]
    col_id = "id" if "id" in df.columns else df.columns[0]
    rows = []
    for _, r in df.iterrows():
        sid = str(r[col_id])
        m = re.search(r"(mol_\d+).*?-(\d+)$", sid)
        if not m:
            m = re.search(r"(mol_\d+)", sid)
            pose = 0
            mol_id = m.group(1) if m else sid
        else:
            mol_id, pose = m.group(1), int(m.group(2))
        rows.append({"mol_id": mol_id, "mode": pose + 1, "rtmscore": float(r[col_score])})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--pool", required=True, help="unused if mol_index_map exists; kept for CLI compat")
    ap.add_argument("--index-map", default="")
    args = ap.parse_args()
    work = Path(args.work)
    idx_path = Path(args.index_map) if args.index_map else work / "mol_index_map.csv"
    if not idx_path.exists():
        raise SystemExit(f"Missing {idx_path}. Run scripts/build_index_map.py first.")

    pool = pd.read_csv(idx_path)
    score_dir = work / "scores"
    score_dir.mkdir(parents=True, exist_ok=True)

    vina_rows = []
    for _, r in pool.iterrows():
        mid = r["mol_id"]
        modes = parse_vina_log_or_pdbqt(
            work / "vina" / f"{mid}_out.pdbqt",
            work / "logs" / "vina" / f"{mid}.log",
        )
        for m in modes:
            vina_rows.append({"mol_id": mid, "canonical_smiles": r.get("canonical_smiles", ""), **m})
    vina_df = pd.DataFrame(vina_rows)
    vina_df.to_csv(score_dir / "vina_modes.csv", index=False)

    gnina_rows = []
    for _, r in pool.iterrows():
        mid = r["mol_id"]
        modes = parse_gnina_sdf(work / "gnina" / f"{mid}_out.sdf")
        for m in modes:
            gnina_rows.append({"mol_id": mid, "canonical_smiles": r.get("canonical_smiles", ""), **m})
    gnina_df = pd.DataFrame(gnina_rows)
    gnina_df.to_csv(score_dir / "gnina_modes.csv", index=False)

    keep = ["mol_id", "canonical_smiles"]
    for c in ("role", "in_true", "in_random", "source_file_index"):
        if c in pool.columns:
            keep.append(c)
    out = pool[keep].copy()

    def add_top(df, score_col, out_col, higher_better: bool):
        if df.empty or score_col not in df.columns:
            out[out_col] = pd.NA
            return
        d = df.dropna(subset=[score_col]).copy()
        d = d.sort_values(["mol_id", score_col], ascending=[True, not higher_better])
        top = d.groupby("mol_id", as_index=False).first()[["mol_id", score_col]]
        out[out_col] = out["mol_id"].map(dict(zip(top["mol_id"], top[score_col])))

    add_top(vina_df, "vina_affinity", "P1_vina_affinity", higher_better=False)
    add_top(gnina_df, "CNNaffinity", "P2_CNNaffinity", higher_better=True)
    add_top(gnina_df, "minimizedAffinity", "P3_gnina_affinity", higher_better=False)
    add_top(gnina_df, "CNNscore", "P0_CNNscore", higher_better=True)

    for eng, col in (("vina", "P4_RTMScore"), ("gnina", "P5_RTMScore")):
        cand = list((work / f"rtmscore_{eng}").glob("*.csv"))
        cand = sorted(cand, key=lambda p: ("score" not in p.name, p.name))
        if not cand:
            out[col] = pd.NA
            continue
        rtm = load_rtm_csv(cand[0])
        add_top(rtm, "rtmscore", col, True)

    out.to_csv(score_dir / "mol_protocol_scores.csv", index=False)
    print(f"wrote {score_dir}/mol_protocol_scores.csv n={len(out)}")
    print(out.isna().mean().to_string())


if __name__ == "__main__":
    main()
