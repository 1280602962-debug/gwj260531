import csv
from pathlib import Path
from collections import defaultdict
import numpy as np

ROOT = Path(__file__).resolve().parents[3]

CONFIGS = {
    "HOAB": {
        "pair": "AChE/BChE",
        "pocket_A": ("4EY7", ROOT / "data/ache_bche_panel_v0/receptors/4EY7_receptor.pdbqt"),
        "pocket_B": ("4BDS", ROOT / "data/ache_bche_panel_v0/receptors/4BDS_receptor.pdbqt"),
        "panel": ROOT / "data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv",
        "pose_root": ROOT / "data/jcim_holdout_v0/HOAB/poses",
    },
    "HOPM": {
        "pair": "PIK3CA/mTOR",
        "pocket_A": ("4L23", ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_receptor.pdbqt"),
        "pocket_B": ("4JT6", ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_receptor.pdbqt"),
        "panel": ROOT / "data/jcim_holdout_v0/tables/holdout_panel_HOPM.csv",
        "pose_root": ROOT / "data/jcim_holdout_v0/HOPM/poses",
    },
    "HOAP": {
        "pair": "PIK3CA/PIK3CB",
        "pocket_A": ("4L23", ROOT / "data/pik3ca_pik3cb_panel_v0/receptors/4L23_receptor.pdbqt"),
        "pocket_B": ("2WXF", ROOT / "data/pik3ca_pik3cb_panel_v0/receptors/2WXF_receptor.pdbqt"),
        "panel": ROOT / "data/jcim_holdout_v0/tables/holdout_panel_HOAP.csv",
        "pose_root": ROOT / "data/jcim_holdout_v0/HOAP/poses",
    },
}


def read_pdbqt_heavy_atoms(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                elem = line[77:79].strip() or line[12:16].strip()[0]
                if elem.upper().startswith("H"):
                    continue
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(coords)


def read_first_model_heavy_atoms(path):
    coords, in_model = [], False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if in_model:
                    break
                in_model = True
                continue
            if line.startswith(("ATOM", "HETATM")):
                elem = line[77:79].strip() or line[12:16].strip()[0]
                if elem.upper().startswith("H"):
                    continue
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            if line.startswith("ENDMDL") and in_model:
                break
    return np.array(coords)


def auroc(pos, neg):
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p, n = np.asarray(pos, float), np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


CUTOFF = 4.0
all_summary = []

for prefix, cfg in CONFIGS.items():
    if not cfg["pose_root"].exists():
        print(f"{prefix}: pose root missing, skip")
        continue
    panel = {}
    with open(cfg["panel"]) as f:
        for r in csv.DictReader(f):
            panel[r["holdout_id"]] = r["class"]

    rec_coords = {}
    for key in ("pocket_A", "pocket_B"):
        pdb_id, path = cfg[key]
        if path.exists():
            rec_coords[pdb_id] = read_pdbqt_heavy_atoms(path)

    records = []
    for key, label in [("pocket_A", "A"), ("pocket_B", "B")]:
        pdb_id, _ = cfg[key]
        pose_dir = cfg["pose_root"] / pdb_id
        if not pose_dir.exists() or pdb_id not in rec_coords:
            continue
        rc = rec_coords[pdb_id]
        for lig_dir in sorted(pose_dir.iterdir()):
            lig_id = lig_dir.name
            out = lig_dir / "out.pdbqt"
            if not out.exists() or lig_id not in panel:
                continue
            lc = read_first_model_heavy_atoms(out)
            if len(lc) == 0:
                continue
            dmat = np.linalg.norm(lc[:, None, :] - rc[None, :, :], axis=2)
            contact = int((dmat.min(axis=1) <= CUTOFF).sum())
            records.append({"ligand": lig_id, "pocket": label, "cls": panel[lig_id],
                             "n_heavy": len(lc), "contact": contact})

    by = defaultdict(list)
    for r in records:
        by[(r["pocket"], r["cls"])].append(r["contact"])

    dA, AA = by[("A", "dual")], by[("A", "A_only")]
    dB, BB = by[("B", "dual")], by[("B", "B_only")]
    auroc_A = auroc(dA, AA)
    auroc_B = auroc(dB, BB)
    print(f"\n=== {prefix} ({cfg['pair']}) ===")
    print(f"  n records: {len(records)}")
    for (pocket, cls), vals in sorted(by.items()):
        if vals:
            print(f"  pocket={pocket} cls={cls:8} n={len(vals):3d} mean_contact={np.mean(vals):.2f}")
    print(f"  own-pocket contact-count AUROC: D vs A_only (pocket A) = {auroc_A:.3f}  |  D vs B_only (pocket B) = {auroc_B:.3f}")
    all_summary.append((prefix, cfg["pair"], auroc_A, auroc_B))

print("\n\n=== Summary: geometric own-pocket contact-count AUROC (mirrors wrong_pocket_control_vina) ===")
print("prefix\tpair\tauroc_D_vs_A_pocketA\tauroc_D_vs_B_pocketB")
for prefix, pair, a, b in all_summary:
    print(f"{prefix}\t{pair}\t{a:.3f}\t{b:.3f}")
