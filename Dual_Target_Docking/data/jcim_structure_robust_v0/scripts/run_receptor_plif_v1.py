#!/usr/bin/env python3
"""Residue-level interaction snapshot for PIK3CA receptor swaps (4L23 / 4JPS / 5DXT).

Optional P1 support — docs/AGENT_COMMAND_RECEPTOR_PLIF_V1.md

Uses local Vina mode-1 poses. Redocks PM48 to alt receptors if missing.
Contact occupancy is a ProLIF-*equivalent* geometric fingerprint:
  heavy-atom proximity ≤ 4.5 Å to each of the 20 frozen pocket residues,
  with optional H-bond / hydrophobic / aromatic flags for the top-shift table.
(Full ProLIF/RDKit conversion segfaults on these all-atom PIK3CA PDBs.)
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from vina import Vina

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_structure_robust_v0"
PM48_ROOT = ROOT / "results/pik3ca_mtor_panel48_rdkit_v0"
LIG_DIR = PM48_ROOT / "ligands_pdbqt"
PANEL = PM48_ROOT / "tables/panel_v0_48.csv"
SEED = 20260727
EXHAUST = 16
OBABEL = Path("/mnt/d/CADD paper exercise/gnina/conda_env/bin/obabel")
CONTACT_CUTOFF = 4.5

POCKET_RESIDUES = [
    ("772", "Met"), ("780", "Trp"), ("800", "Ile"), ("802", "Lys"), ("807", "Leu"),
    ("810", "Asp"), ("814", "Leu"), ("836", "Tyr"), ("838", "Cys"), ("848", "Ile"),
    ("849", "Glu"), ("850", "Val"), ("851", "Val"), ("854", "Ser"), ("856", "Thr"),
    ("859", "Gln"), ("922", "Met"), ("930", "Phe"), ("932", "Ile"), ("933", "Asp"),
]
POCKET_LABELS = [f"{aa}{num}" for num, aa in POCKET_RESIDUES]
AROMATIC = {"Phe", "Tyr", "Trp", "His"}
HYDROPHOBIC = {"Ala", "Val", "Leu", "Ile", "Met", "Phe", "Pro", "Trp", "Tyr"}
HB_SIDECHAIN = {"Asp", "Glu", "Asn", "Gln", "Ser", "Thr", "Tyr", "Lys", "Arg", "His", "Cys"}

RECEPTORS = {
    "4L23": {
        "protein": PM48_ROOT / "receptors/4L23_protein.pdb",
        "receptor_pdbqt": PM48_ROOT / "receptors/4L23_receptor.pdbqt",
        "box": PM48_ROOT / "boxes/4L23_box.json",
        "pose_root": PM48_ROOT / "poses/4L23",
    },
    "4JPS": {
        "protein": OUT / "receptors/4JPS_protein.pdb",
        "receptor_pdbqt": OUT / "receptors/4JPS_receptor.pdbqt",
        "box": OUT / "receptors/4JPS_box.json",
        "pose_root": OUT / "poses/4JPS",
    },
    "5DXT": {
        "protein": OUT / "receptors/5DXT_protein.pdb",
        "receptor_pdbqt": OUT / "receptors/5DXT_receptor.pdbqt",
        "box": OUT / "receptors/5DXT_box.json",
        "pose_root": OUT / "poses/5DXT",
    },
}

REPRESENTATIVE = {
    "dual": ["PM48_01", "PM48_02", "PM48_05"],
    "A_only": ["PM48_22", "PM48_26", "PM48_28"],
    "B_only": ["PM48_33", "PM48_40", "PM48_44"],
}


def split_modes(out_pdbqt: Path, dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    text = out_pdbqt.read_text().splitlines()
    models, cur = [], []
    for line in text:
        if line.startswith("MODEL"):
            cur = [line]
        elif line.startswith("ENDMDL"):
            cur.append(line)
            models.append(cur)
            cur = []
        elif cur:
            cur.append(line)
    for i, m in enumerate(models, 1):
        (dest_dir / f"mode_{i:02d}.pdbqt").write_text("\n".join(m) + "\n")


def dock_one(args):
    receptor_pdbqt, box_path, lig_id, lig_pdbqt, pose_dir = args
    pose_dir = Path(pose_dir)
    pose_dir.mkdir(parents=True, exist_ok=True)
    mode1 = pose_dir / "mode_01.pdbqt"
    if mode1.exists():
        return {"ligand": lig_id, "status": "exists", "seconds": 0.0}
    out_pdbqt = pose_dir / "out.pdbqt"
    box = json.loads(Path(box_path).read_text())
    t0 = time.time()
    try:
        v = Vina(sf_name="vina", cpu=1, seed=SEED, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(lig_pdbqt))
        v.compute_vina_maps(
            center=[box["center_x"], box["center_y"], box["center_z"]],
            box_size=[box["size_x"], box["size_y"], box["size_z"]],
        )
        v.dock(exhaustiveness=EXHAUST, n_poses=9)
        v.write_poses(str(out_pdbqt), n_poses=9, overwrite=True, energy_range=3)
        split_modes(out_pdbqt, pose_dir)
        return {"ligand": lig_id, "status": "success", "seconds": round(time.time() - t0, 1)}
    except Exception as exc:  # noqa: BLE001
        return {
            "ligand": lig_id,
            "status": "fail",
            "reason": str(exc)[:300],
            "seconds": round(time.time() - t0, 1),
        }


def ensure_poses(receptor: str, ligands: list[str], workers: int) -> None:
    cfg = RECEPTORS[receptor]
    jobs = []
    for lig in ligands:
        pdbqt = LIG_DIR / f"{lig}.pdbqt"
        if not pdbqt.exists():
            continue
        jobs.append(
            (
                str(cfg["receptor_pdbqt"]),
                str(cfg["box"]),
                lig,
                str(pdbqt),
                str(cfg["pose_root"] / lig),
            )
        )
    print(f"{receptor}: docking {len(jobs)} ligands", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(dock_one, j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            if i % 10 == 0 or res["status"] not in ("success", "exists"):
                print(f"  [{i}/{len(jobs)}] {res}", flush=True)


def pdb_heavy_atoms(path: Path) -> dict[str, np.ndarray]:
    """resid label -> Nx3 heavy-atom coords. Label uses 3-letter capitalize + resid."""
    by_res: dict[str, list] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if name.startswith("H"):
            continue
        elem = (line[76:78].strip() if len(line) >= 78 else "") or name[:1]
        if elem.upper() == "H":
            continue
        resname = line[17:20].strip().capitalize()
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        label = f"{resname}{resid}"
        by_res.setdefault(label, []).append(xyz)
    return {k: np.asarray(v, dtype=float) for k, v in by_res.items()}


def pdbqt_heavy_xyz(path: Path) -> np.ndarray:
    xyz = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if name.startswith("H"):
            continue
        # AutoDock atom type is last token; skip hydrogens typed HD/H
        toks = line.split()
        if toks and toks[-1] in {"HD", "H"}:
            continue
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(xyz, dtype=float)


def residue_contacts(lig_xyz: np.ndarray, pocket_xyz: dict[str, np.ndarray]) -> dict[str, dict]:
    out = {}
    if lig_xyz.size == 0:
        return out
    for label in POCKET_LABELS:
        # pocket maps may use different capitalization / missing residues
        coords = pocket_xyz.get(label)
        if coords is None:
            # try uppercase 3-letter
            aa = "".join(c for c in label if c.isalpha())
            num = "".join(c for c in label if c.isdigit())
            for k, v in pocket_xyz.items():
                if k.endswith(num) and k.lower().startswith(aa[:3].lower()):
                    coords = v
                    break
        if coords is None or len(coords) == 0:
            continue
        d = np.linalg.norm(lig_xyz[:, None, :] - coords[None, :, :], axis=2)
        mind = float(d.min())
        if mind > CONTACT_CUTOFF:
            continue
        aa = "".join(c for c in label if c.isalpha())
        out[label] = {
            "occupancy": 1,
            "min_dist": round(mind, 3),
            "hydrophobic": int(aa in HYDROPHOBIC),
            "aromatic": int(aa in AROMATIC),
            "hbond_capable": int(aa in HB_SIDECHAIN),
        }
    return out


def occupancy_matrix(ligands: list[str]) -> pd.DataFrame:
    rows = []
    panel = pd.read_csv(PANEL)
    cls_map = {r["panel_id"]: r["class"] for _, r in panel.iterrows()}

    pocket_cache = {}
    for receptor, cfg in RECEPTORS.items():
        print(f"loading pocket atoms: {receptor}", flush=True)
        all_xyz = pdb_heavy_atoms(cfg["protein"])
        # keep only pocket resid numbers of interest (any aa name at that resid)
        pocket = {}
        for num, aa in POCKET_RESIDUES:
            label = f"{aa}{num}"
            if label in all_xyz:
                pocket[label] = all_xyz[label]
            else:
                hits = [k for k in all_xyz if k.endswith(num)]
                if hits:
                    pocket[label] = all_xyz[hits[0]]
                    if hits[0] != label:
                        print(f"  {receptor}: mapped {hits[0]} -> {label}", flush=True)
                else:
                    print(f"  WARN {receptor}: missing residue {label}", flush=True)
        pocket_cache[receptor] = pocket
        print(f"  {receptor}: {len(pocket)}/{len(POCKET_LABELS)} pocket residues found", flush=True)

    for receptor, cfg in RECEPTORS.items():
        print(f"contacts: {receptor}", flush=True)
        for i, lig in enumerate(ligands, 1):
            pose = cfg["pose_root"] / lig / "mode_01.pdbqt"
            if not pose.exists():
                print(f"MISSING pose {pose}", flush=True)
                continue
            lig_xyz = pdbqt_heavy_xyz(pose)
            contacts = residue_contacts(lig_xyz, pocket_cache[receptor])
            if i % 10 == 0 or i == 1:
                print(f"  [{i}/{len(ligands)}] {lig} n_contacts={len(contacts)}", flush=True)
            for res in POCKET_LABELS:
                info = contacts.get(res, {})
                rows.append(
                    {
                        "receptor": receptor,
                        "ligand": lig,
                        "class": cls_map.get(lig, ""),
                        "residue": res,
                        "occupancy": int(info.get("occupancy", 0)),
                        "min_dist": info.get("min_dist", ""),
                        "hydrophobic": int(info.get("hydrophobic", 0)),
                        "aromatic": int(info.get("aromatic", 0)),
                        "hbond_capable": int(info.get("hbond_capable", 0)),
                    }
                )
    return pd.DataFrame(rows)


def plot_heatmap(df: pd.DataFrame, ligands: list[str], out_png: Path) -> None:
    mean_by_rec = df.groupby(["receptor", "residue"])["occupancy"].mean().unstack("receptor")
    for col in ("4L23", "4JPS", "5DXT"):
        if col not in mean_by_rec.columns:
            mean_by_rec[col] = 0.0
    shift = (mean_by_rec["4JPS"] - mean_by_rec["4L23"]).abs()
    shift2 = (mean_by_rec["5DXT"] - mean_by_rec["4L23"]).abs()
    top = (shift + shift2).sort_values(ascending=False).head(8).index.tolist()

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    im = None
    for ax, rec in zip(axes, ["4L23", "4JPS", "5DXT"]):
        sub = df[(df["receptor"] == rec) & (df["ligand"].isin(ligands))]
        mat = sub.pivot_table(index="ligand", columns="residue", values="occupancy", fill_value=0)
        mat = mat.reindex(index=ligands, columns=top, fill_value=0)
        im = ax.imshow(mat.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
        ax.set_title(rec)
        ax.set_xticks(range(len(top)))
        ax.set_xticklabels(top, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(ligands)))
        ax.set_yticklabels(ligands, fontsize=7)
    if im is not None:
        fig.colorbar(im, ax=axes, fraction=0.02, label="contact (0/1)")
    fig.suptitle("PIK3CA pocket contact occupancy — representative ligands")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_verdict(top: pd.DataFrame, path: Path) -> None:
    residues = ", ".join(top.index.tolist()[:8])
    text = f"""# PLIF / pocket-contact verdict v1

Method: geometric heavy-atom proximity ≤ {CONTACT_CUTOFF} Å to the 20 frozen
PIK3CA pocket residues (ProLIF/RDKit conversion segfaults on these PDBs;
this is the SOP-allowed equivalent occupancy snapshot).

Top occupancy-shift residues (4JPS/5DXT vs 4L23): {residues}

Allowed claim:
> The performance shift coincided with altered interaction patterns at residues
> {residues}, providing a structural **hypothesis** for the receptor sensitivity.

Forbidden: residue X caused the AUROC change; PLIF explains the opposite
PIK3CA/PIK3CB shift.
"""
    path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--skip-dock", action="store_true")
    args = ap.parse_args()

    panel = pd.read_csv(PANEL)
    all_ligs = panel["panel_id"].tolist()
    rep_ligs = []
    for ids in REPRESENTATIVE.values():
        rep_ligs.extend(ids)

    if not args.skip_dock:
        for rec in ("4JPS", "5DXT"):
            ensure_poses(rec, all_ligs, args.workers)

    df = occupancy_matrix(all_ligs)
    tab_dir = OUT / "analysis/plif_v1"
    tab_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(tab_dir / "plif_occupancy_all_v1.csv", index=False)

    mean_occ = df.groupby(["receptor", "residue"])["occupancy"].mean().reset_index()
    wide = mean_occ.pivot(index="residue", columns="receptor", values="occupancy").fillna(0)
    for col in ("4L23", "4JPS", "5DXT"):
        if col not in wide.columns:
            wide[col] = 0.0
    wide["shift_4JPS_vs_4L23"] = wide["4JPS"] - wide["4L23"]
    wide["shift_5DXT_vs_4L23"] = wide["5DXT"] - wide["4L23"]
    wide["abs_shift_max"] = wide[["shift_4JPS_vs_4L23", "shift_5DXT_vs_4L23"]].abs().max(axis=1)
    top = wide.sort_values("abs_shift_max", ascending=False).head(10)
    top.to_csv(tab_dir / "plif_residue_shift_top10_v1.csv")
    print("Top shifted residues:\n", top[["4L23", "4JPS", "5DXT", "shift_4JPS_vs_4L23", "shift_5DXT_vs_4L23"]])

    rep_df = df[df["ligand"].isin(rep_ligs)]
    rep_df.to_csv(tab_dir / "plif_occupancy_representative_v1.csv", index=False)
    plot_heatmap(rep_df, rep_ligs, tab_dir / "plif_heatmap_representative_v1.png")
    write_verdict(top, tab_dir / "PLIF_VERDICT_V1.md")
    print(f"wrote {tab_dir}", flush=True)


if __name__ == "__main__":
    main()
