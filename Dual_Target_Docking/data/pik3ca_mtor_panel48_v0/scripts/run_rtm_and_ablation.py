#!/usr/bin/env python3
"""RTMScore rescoring + ablation tables for pik3ca_mtor_panel48_v0.

Pipeline:
  1) Convert all mode_*.pdbqt -> combined SDF per target (obabel via cadd_tools)
  2) Score with RTMScore model1 against frozen 10A pockets
  3) Write scores_rtm*.csv and ablation_* tables (EGFR-aligned arms)

Usage:
  /home/gwj/miniconda3/envs/rtmscore/bin/python scripts/run_rtm_and_ablation.py
  # or stepwise:
  ... --prep-only
  ... --score-only
  ... --ablation-only
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
POSES = ROOT / "poses"
LOGS = ROOT / "logs" / "rtmscore"
RECEPTORS = ROOT / "receptors"
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
TARGETS = ["4L23", "4JT6"]
ARMS = ["vina_mean", "vina_min", "rtm_mean", "rtm_min", "rtm_min_z"]
KEY_LIGS = [
    "PM48_01",  # PI-103 cognate
    "PM48_02",
    "PM48_04",
    "PM48_10",
    "PM48_19",
    "PM48_25",
    "PM48_34",
    "PM48_45",
]
CLASH_CUTOFF = 2.2
CLASH_FAIL_N = 3


def bash(cmd: str, log_path: Path | None = None) -> int:
    if log_path is None:
        proc = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + "\n" + proc.stderr + "\n")
        return proc.returncode
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as fh:
        proc = subprocess.run(["bash", "-lc", cmd], stdout=fh, stderr=subprocess.STDOUT)
    return proc.returncode


LIGANDS_SDF = ROOT / "ligands_sdf"


def _pdbqt_xyz(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return xyz


def _smiles_idx_pairs(path: Path):
    nums: list[int] = []
    for line in path.read_text().splitlines():
        if line.startswith("REMARK SMILES IDX"):
            nums.extend(int(x) for x in line.split()[3:])
    if len(nums) % 2:
        raise ValueError(f"odd SMILES IDX count in {path}")
    return list(zip(nums[0::2], nums[1::2]))


def _vina_from_pdbqt(path: Path):
    for line in path.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def _load_templates():
    from rdkit import Chem

    templates = {}
    for p in sorted(LIGANDS_SDF.glob("PM48_*.sdf")):
        mol = Chem.SDMolSupplier(str(p), removeHs=False)[0]
        if mol is None:
            raise RuntimeError(f"failed to read template {p}")
        templates[p.stem] = Chem.RemoveHs(mol)
    return templates


def prep_sdfs() -> Path:
    from rdkit import Chem

    LOGS.mkdir(parents=True, exist_ok=True)
    templates = _load_templates()
    index_rows = []
    for target in TARGETS:
        out_sdf = LOGS / f"{target}_all_poses.sdf"
        writer = Chem.SDWriter(str(out_sdf))
        n = 0
        lig_dirs = sorted((POSES / target).glob("PM48_*"))
        print(f"convert {target}: ...", flush=True)
        for lig_dir in lig_dirs:
            lig = lig_dir.name
            tmpl = templates[lig]
            for pose in sorted(lig_dir.glob("mode_*.pdbqt")):
                m = re.search(r"mode_(\d+)", pose.name)
                if not m:
                    continue
                mode = int(m.group(1))
                title = f"{lig}_mode{mode}"
                xyz = _pdbqt_xyz(pose)
                pairs = _smiles_idx_pairs(pose)
                if len(pairs) != tmpl.GetNumAtoms():
                    raise RuntimeError(
                        f"{pose}: SMILES IDX heavy={len(pairs)} != template {tmpl.GetNumAtoms()}"
                    )
                mol = Chem.Mol(tmpl)
                conf = Chem.Conformer(mol.GetNumAtoms())
                for s_idx, p_idx in pairs:
                    conf.SetAtomPosition(s_idx - 1, xyz[p_idx - 1])
                mol.RemoveAllConformers()
                mol.AddConformer(conf, assignId=True)
                mol.SetProp("_Name", title)
                writer.write(mol)
                n += 1
                index_rows.append(
                    {
                        "target": target,
                        "ligand": lig,
                        "mode": mode,
                        "vina_score": _vina_from_pdbqt(pose),
                        "sdf_name": title,
                        "pose_path": str(pose),
                    }
                )
        writer.close()
        suppl = Chem.SDMolSupplier(str(out_sdf), removeHs=False)
        ok = sum(1 for m in suppl if m is not None)
        print(f"wrote {out_sdf} n={n} rdkit_ok={ok}", flush=True)
        if ok != n:
            raise RuntimeError(f"SDF parse mismatch for {target}: ok={ok} n={n}")
    idx = TABLES / "poses_index_rtm.csv"
    with idx.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(index_rows[0].keys()))
        w.writeheader()
        w.writerows(index_rows)
    print(f"wrote {idx} n={len(index_rows)}", flush=True)
    return idx


def run_rtm() -> None:
    rtm_python = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
    for target in TARGETS:
        sdf = LOGS / f"{target}_all_poses.sdf"
        pocket = RECEPTORS / f"{target}_pocket_10.0.pdb"
        if not sdf.exists():
            raise FileNotFoundError(sdf)
        if not pocket.exists():
            raise FileNotFoundError(pocket)
        out_prefix = LOGS / f"{target}_rtmscore_all_poses"
        log_path = LOGS / f"{target}_rtmscore.log"
        print(f"RTM {target} ...", flush=True)
        with log_path.open("w") as fh:
            proc = subprocess.run(
                [
                    str(rtm_python),
                    str(RTM_PY),
                    "-p",
                    str(pocket),
                    "-l",
                    str(sdf),
                    "-m",
                    str(MODEL),
                    "-o",
                    str(out_prefix),
                ],
                cwd=str(RTM_ROOT / "example"),
                stdout=fh,
                stderr=subprocess.STDOUT,
            )
        if proc.returncode != 0:
            raise RuntimeError(f"RTM failed for {target}; see {log_path}")
        csv_path = Path(f"{out_prefix}.csv")
        if not csv_path.exists():
            alt = RTM_ROOT / "example" / f"{out_prefix.name}.csv"
            if alt.exists():
                alt.rename(csv_path)
        if not csv_path.exists():
            # sometimes written as out_prefix without checking path with spaces
            found = list(LOGS.glob(f"{target}_rtmscore_all_poses*.csv"))
            if found:
                csv_path = found[0]
            else:
                raise FileNotFoundError(f"missing RTM csv for {target}: {csv_path}")
        print(f"OK {csv_path}", flush=True)


def parse_rtm_id(s: str) -> tuple[str, int]:
    m = re.match(r"(PM48_\d+)_mode(\d+)", str(s))
    if not m:
        raise ValueError(f"bad rtm id: {s}")
    return m.group(1), int(m.group(2))


def roc_auc_score(y_true, y_score):
    try:
        from sklearn.metrics import roc_auc_score as sk_auc

        return float(sk_auc(y_true, y_score))
    except Exception:
        y_true = np.asarray(y_true)
        y_score = np.asarray(y_score)
        pos = y_score[y_true == 1]
        neg = y_score[y_true == 0]
        if len(pos) == 0 or len(neg) == 0:
            return float("nan")
        correct = 0.0
        for p in pos:
            correct += np.sum(p > neg) + 0.5 * np.sum(p == neg)
        return correct / (len(pos) * len(neg))


def parse_pdbqt_heavy(path: Path):
    coords = []
    if not path.exists():
        return np.zeros((0, 3))
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line.split()[-1] in ("H", "HD", "HS"):
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords) if coords else np.zeros((0, 3))


def parse_protein_heavy(pdb_path: Path):
    coords = []
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        elem = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        if elem.startswith("H"):
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords)


def clash_count(lig_xyz, prot_xyz, cutoff=CLASH_CUTOFF):
    if len(lig_xyz) == 0 or len(prot_xyz) == 0:
        return 0
    n = 0
    for i in range(0, len(lig_xyz), 20):
        L = lig_xyz[i : i + 20]
        d = np.linalg.norm(L[:, None, :] - prot_xyz[None, :, :], axis=2)
        n += int((d < cutoff).sum())
    return n


def build_tables() -> None:
    long_rows = []
    for target in TARGETS:
        csv_path = LOGS / f"{target}_rtmscore_all_poses.csv"
        df = pd.read_csv(csv_path)
        # columns typically: id, score
        id_col = "id" if "id" in df.columns else df.columns[0]
        score_col = "score" if "score" in df.columns else df.columns[1]
        for _, r in df.iterrows():
            lig, mode = parse_rtm_id(r[id_col])
            long_rows.append(
                {
                    "target": target,
                    "ligand": lig,
                    "mode": mode,
                    "rtmscore": float(r[score_col]),
                }
            )
    long_df = pd.DataFrame(long_rows)
    # attach vina from poses_index if present
    idx_path = TABLES / "poses_index_rtm.csv"
    if idx_path.exists():
        idx = pd.read_csv(idx_path)
        long_df = long_df.merge(
            idx[["target", "ligand", "mode", "vina_score"]],
            on=["target", "ligand", "mode"],
            how="left",
        )
    long_df.to_csv(TABLES / "scores_rtm_all_poses.csv", index=False)

    # per-target best
    best = (
        long_df.sort_values(["target", "ligand", "rtmscore"], ascending=[True, True, False])
        .groupby(["target", "ligand"], as_index=False)
        .first()
        .rename(columns={"mode": "best_rtm_mode", "rtmscore": "rtmscore"})
    )
    for target in TARGETS:
        sub = best[best.target == target][
            ["ligand", "rtmscore", "best_rtm_mode", "vina_score"]
        ].copy()
        sub.to_csv(TABLES / f"scores_rtm_{target}.csv", index=False)

    # wide scores_rtm.csv
    panel = pd.read_csv(TABLES / "panel_v0_48.csv")
    vina = pd.read_csv(TABLES / "scores_vina.csv")
    a = best[best.target == "4L23"][["ligand", "rtmscore", "best_rtm_mode"]].rename(
        columns={"rtmscore": "4L23_rtmscore", "best_rtm_mode": "4L23_best_rtm_mode"}
    )
    b = best[best.target == "4JT6"][["ligand", "rtmscore", "best_rtm_mode"]].rename(
        columns={"rtmscore": "4JT6_rtmscore", "best_rtm_mode": "4JT6_best_rtm_mode"}
    )
    wide = panel.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "pref_name", "molecule_chembl_id"]
    ].merge(a, on="ligand").merge(b, on="ligand")
    wide = wide.merge(
        vina.rename(columns={"panel_id": "ligand"})[
            ["ligand", "vina_4L23_mode1", "vina_4JT6_mode1", "vina_mean", "vina_min"]
        ],
        on="ligand",
        how="left",
    )
    # rename vina columns for clarity in combined file
    wide = wide.rename(
        columns={
            "vina_4L23_mode1": "4L23_affinity",
            "vina_4JT6_mode1": "4JT6_affinity",
            "vina_mean": "mean_affinity",
            "vina_min": "min_affinity_raw",  # note: scores_vina vina_min was max of affinities (worse)
        }
    )
    wide["mean_rtmscore"] = (wide["4L23_rtmscore"] + wide["4JT6_rtmscore"]) / 2
    wide["min_rtmscore"] = wide[["4L23_rtmscore", "4JT6_rtmscore"]].min(axis=1)
    wide = wide.sort_values("mean_rtmscore", ascending=False).reset_index(drop=True)
    wide.insert(0, "dual_rank_by_mean", range(1, len(wide) + 1))
    wide_min = wide.sort_values("min_rtmscore", ascending=False).reset_index(drop=True)
    rank_min = {r.ligand: i + 1 for i, r in wide_min.iterrows()}
    wide["dual_rank_by_min"] = wide["ligand"].map(rank_min)
    wide.to_csv(TABLES / "scores_rtm.csv", index=False)

    # --- ablation (higher-better scores) ---
    df = panel.rename(columns={"panel_id": "ligand"})[
        ["ligand", "class", "pref_name", "molecule_chembl_id"]
    ].copy()
    df = df.merge(
        vina.rename(columns={"panel_id": "ligand"})[
            ["ligand", "vina_4L23_mode1", "vina_4JT6_mode1"]
        ],
        on="ligand",
        how="left",
    )
    df = df.merge(
        a.rename(columns={"4L23_rtmscore": "rtm_4L23", "4L23_best_rtm_mode": "rtm_mode_4L23"}),
        on="ligand",
        how="left",
    )
    df = df.merge(
        b.rename(columns={"4JT6_rtmscore": "rtm_4JT6", "4JT6_best_rtm_mode": "rtm_mode_4JT6"}),
        on="ligand",
        how="left",
    )
    # higher-better vina
    df["vina_4L23_hb"] = -df["vina_4L23_mode1"]
    df["vina_4JT6_hb"] = -df["vina_4JT6_mode1"]
    df["vina_mean"] = (df["vina_4L23_hb"] + df["vina_4JT6_hb"]) / 2
    df["vina_min"] = df[["vina_4L23_hb", "vina_4JT6_hb"]].min(axis=1)
    df["rtm_mean"] = (df["rtm_4L23"] + df["rtm_4JT6"]) / 2
    df["rtm_min"] = df[["rtm_4L23", "rtm_4JT6"]].min(axis=1)
    for col, zcol in [("rtm_4L23", "rtm_4L23_z"), ("rtm_4JT6", "rtm_4JT6_z")]:
        mu, sd = df[col].mean(), df[col].std(ddof=0)
        df[zcol] = (df[col] - mu) / (sd if sd > 0 else 1.0)
    df["rtm_min_z"] = df[["rtm_4L23_z", "rtm_4JT6_z"]].min(axis=1)

    # optional clash gate on RTM-best poses
    prot = {
        "4L23": parse_protein_heavy(RECEPTORS / "4L23_protein.pdb"),
        "4JT6": parse_protein_heavy(RECEPTORS / "4JT6_protein.pdb"),
    }
    clash_rows = []
    for _, r in df.iterrows():
        lig = r["ligand"]
        rec = {"ligand": lig}
        fail = False
        for tgt, mode_col in [("4L23", "rtm_mode_4L23"), ("4JT6", "rtm_mode_4JT6")]:
            mode = int(r[mode_col]) if pd.notna(r[mode_col]) else 1
            pose = POSES / tgt / lig / f"mode_{mode:02d}.pdbqt"
            n = clash_count(parse_pdbqt_heavy(pose), prot[tgt])
            rec[f"clash_{tgt}"] = n
            if n >= CLASH_FAIL_N:
                fail = True
        rec["clash_fail"] = int(fail)
        clash_rows.append(rec)
    clash_df = pd.DataFrame(clash_rows)
    df = df.merge(clash_df, on="ligand", how="left")
    df["gated_rtm_min"] = df["rtm_min"].where(df["clash_fail"] == 0, other=-1e9)
    arms = ARMS + ["gated_rtm_min"]

    df.to_csv(TABLES / "ablation_ligand_scores.csv", index=False)

    # ranks (1 = best)
    ranks = df[["ligand", "class"]].copy()
    for arm in arms:
        ranks[arm] = df[arm].rank(ascending=False, method="min").astype(int)
    ranks.to_csv(TABLES / "ablation_ranks.csv", index=False)

    # metrics Dual vs A∪B∪neither (hardneg = not dual)
    y = (df["class"] == "dual").astype(int).values
    hard = df["class"] != "dual"
    metrics = []
    for arm in arms:
        score = df[arm].values
        auroc = roc_auc_score(y, score)
        order = np.argsort(-score)
        top10 = df.iloc[order[:10]]
        n_dual = int((top10["class"] == "dual").sum())
        n_hard = int((top10["class"] != "dual").sum())
        metrics.append(
            {
                "arm": arm,
                "auroc_dual_vs_rest": auroc,
                "top10_dual": n_dual,
                "top10_hardneg": n_hard,
                "top10_hardneg_frac": n_hard / 10.0,
                "n_dual": int(y.sum()),
                "n_rest": int((1 - y).sum()),
            }
        )
    met = pd.DataFrame(metrics)
    met.to_csv(TABLES / "ablation_metrics.csv", index=False)

    # key ligand ranks
    key = ranks[ranks.ligand.isin(KEY_LIGS)].copy()
    key.to_csv(TABLES / "ablation_key_ligand_ranks.csv", index=False)

    print("Wrote ablation tables:")
    print(met.to_string(index=False))
    print("\nKey ligand ranks:")
    print(key.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prep-only", action="store_true")
    ap.add_argument("--score-only", action="store_true")
    ap.add_argument("--ablation-only", action="store_true")
    args = ap.parse_args()
    do_all = not (args.prep_only or args.score_only or args.ablation_only)
    if do_all or args.prep_only:
        prep_sdfs()
    if do_all or args.score_only:
        run_rtm()
    if do_all or args.ablation_only:
        build_tables()


if __name__ == "__main__":
    main()
