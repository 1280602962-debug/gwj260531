#!/usr/bin/env python3
"""RTMScore best-of-9 on production Track B Vina poses.

Copies the K=4 panel_rtm_scores protocol (rtmscore_model1, SMILES-IDX SDF
rebuild). Does not dock. Does not replace Table 2.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
TABLES = ROOT / "tables"
CACHE_SDF = ROOT / "cache" / "track_b_ligands" / "sdf"

RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_PY = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"
RTM_PYTHON = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")

PAIRS = [
    ("F2/F10", ["4UDW", "2JKH"], TABLES / "track_b_panels" / "panel_F2_F10_v1.csv"),
    ("JAK1/TYK2", ["6N7A", "3LXP"], TABLES / "track_b_panels" / "panel_JAK1_TYK2_v1.csv"),
    ("JAK1/JAK2", ["6N7A", "8BXH"], TABLES / "track_b_panels" / "panel_JAK1_JAK2_v1.csv"),
    ("PPARG/PPARA", ["9V8H", "6LXA"], TABLES / "track_b_panels" / "panel_PPARG_PPARA_v1.csv"),
    ("PPARA/PPARD", ["6LXA", "5U3Q"], TABLES / "track_b_panels" / "panel_PPARA_PPARD_v1.csv"),
]


def pdbqt_xyz(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return xyz


def smiles_idx_pairs(path: Path):
    nums = []
    for line in path.read_text().splitlines():
        if line.startswith("REMARK SMILES IDX"):
            nums.extend(int(x) for x in line.split()[3:])
    if not nums:
        return None
    return list(zip(nums[0::2], nums[1::2]))


def write_sdfs(targets: list[str]) -> None:
    logs = LOCAL / "logs" / "rtmscore"
    logs.mkdir(parents=True, exist_ok=True)
    for target in targets:
        out_sdf = logs / f"{target}_poses.sdf"
        if out_sdf.exists() and out_sdf.stat().st_size > 1000:
            print("reuse", out_sdf, flush=True)
            continue
        pose_root = LOCAL / "poses" / target
        if not pose_root.exists():
            raise SystemExit(f"missing poses {pose_root} — regenerate production Vina first")
        w = Chem.SDWriter(str(out_sdf))
        n = 0
        for lig_dir in sorted(pose_root.iterdir()):
            if not lig_dir.is_dir():
                continue
            lig = lig_dir.name
            sdf_lig = CACHE_SDF / f"{lig}.sdf"
            if not sdf_lig.exists():
                sdf_lig = LOCAL / "ligands_sdf" / f"{lig}.sdf"
            modes = sorted(lig_dir.glob("mode_*.pdbqt"))
            if not modes or not sdf_lig.exists():
                continue
            tmpl = Chem.RemoveHs(Chem.SDMolSupplier(str(sdf_lig), removeHs=False)[0])
            if tmpl is None:
                continue
            for mp in modes:
                mode = int(re.search(r"mode_(\d+)", mp.name).group(1))
                xyz = pdbqt_xyz(mp)
                pairs = smiles_idx_pairs(mp)
                if not pairs or len(pairs) != tmpl.GetNumAtoms():
                    print(f"skip pairs {mp}", flush=True)
                    continue
                mol = Chem.Mol(tmpl)
                conf = Chem.Conformer(mol.GetNumAtoms())
                for s_idx, p_idx in pairs:
                    conf.SetAtomPosition(s_idx - 1, xyz[p_idx - 1])
                mol.RemoveAllConformers()
                mol.AddConformer(conf, assignId=True)
                mol.SetProp("_Name", f"{lig}_mode{mode}")
                w.write(mol)
                n += 1
        w.close()
        print("wrote", out_sdf, n, flush=True)


def run_rtm(target: str, protein: Path, reflig: Path | None, rtm_python: Path, rtm_py: Path, model: Path, rtm_root: Path):
    logs = LOCAL / "logs" / "rtmscore"
    sdf = logs / f"{target}_poses.sdf"
    out_prefix = logs / f"{target}_rtmscore"
    csv_path = Path(f"{out_prefix}.csv")
    if csv_path.exists() and csv_path.stat().st_size > 100:
        print("reuse", csv_path, flush=True)
        return csv_path
    pocket = LOCAL / "receptors" / f"{target}_pocket_10.0.pdb"
    log = logs / f"{target}_rtmscore.log"
    print("RTM", target, flush=True)
    if pocket.exists() and pocket.stat().st_size > 100:
        cmd = [str(rtm_python), str(rtm_py), "-p", str(pocket), "-l", str(sdf), "-m", str(model), "-o", str(out_prefix)]
    else:
        if reflig is None or not reflig.exists():
            raise SystemExit(f"need pocket or cognate reflig for {target}")
        cmd = [
            str(rtm_python),
            str(rtm_py),
            "-p",
            str(protein),
            "-l",
            str(sdf),
            "-m",
            str(model),
            "-o",
            str(out_prefix),
            "-gen_pocket",
            "-c",
            "10.0",
            "-rl",
            str(reflig),
        ]
    with log.open("w") as fh:
        proc = subprocess.run(cmd, cwd=str(rtm_root / "example"), stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise SystemExit(f"RTM failed {target}; see {log}")
    gen = protein.with_name(protein.stem + "_pocket_10.0.pdb")
    if gen.exists() and not pocket.exists():
        pocket.write_bytes(gen.read_bytes())
    if not csv_path.exists():
        alt = rtm_root / "example" / f"{out_prefix.name}.csv"
        if alt.exists():
            alt.rename(csv_path)
    print("OK", csv_path, flush=True)
    return csv_path


def parse_id(s):
    m = re.search(r"([A-Za-z0-9]+_\d+)_mode(\d+)", str(s))
    if not m:
        raise ValueError(s)
    return m.group(1), int(m.group(2))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rtm-python", type=Path, default=RTM_PYTHON)
    ap.add_argument("--rtm-root", type=Path, default=RTM_ROOT)
    ap.add_argument("--model", type=Path, default=MODEL)
    args = ap.parse_args()
    rtm_py = args.rtm_root / "example" / "rtmscore.py"
    if not args.rtm_python.exists():
        raise SystemExit(f"RTM python missing: {args.rtm_python}")
    if not (LOCAL / "poses").exists():
        raise SystemExit("local_track_b_v0/poses/ missing — regenerate production Vina first")

    targets = sorted({t for _, ts, _ in PAIRS for t in ts})
    write_sdfs(targets)
    for t in targets:
        protein = LOCAL / "receptors" / f"{t}_protein.pdb"
        if not protein.exists():
            protein = LOCAL / "receptors" / f"{t}_receptor.pdb"
        reflig = None
        cog = LOCAL / "cognates"
        if cog.exists():
            hits = list(cog.glob(f"{t}_*.sdf")) + list(cog.glob(f"{t}_*.pdb"))
            reflig = hits[0] if hits else None
        run_rtm(t, protein, reflig, args.rtm_python, rtm_py, args.model, args.rtm_root)

    rows = []
    for pair, ts, panel_csv in PAIRS:
        panel = list(csv.DictReader(panel_csv.open()))
        for t in ts:
            d = pd.read_csv(LOCAL / "logs" / "rtmscore" / f"{t}_rtmscore.csv")
            id_col = "id" if "id" in d.columns else d.columns[0]
            sc_col = "score" if "score" in d.columns else d.columns[1]
            by = {}
            for _, r in d.iterrows():
                lig, mode = parse_id(r[id_col])
                sc = float(r[sc_col])
                prev = by.get(lig)
                if prev is None or sc > prev[0]:
                    by[lig] = (sc, mode)
            for rec in panel:
                lig = rec["panel_id"]
                if lig in by:
                    rows.append(
                        {
                            "pair": pair,
                            "target": t,
                            "ligand": lig,
                            "rtm_best": by[lig][0],
                            "best_mode": by[lig][1],
                            "status": "success",
                        }
                    )
    out = LOCAL / "tables" / "scores_rtm_best9_v1.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair", "target", "ligand", "rtm_best", "best_mode", "status"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, "n=", len(rows), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
