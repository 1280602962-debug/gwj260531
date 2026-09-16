#!/usr/bin/env python3
"""Adopt corrected EGFR/HER2 boxes and AChE/BChE panel as official canonical.

Pre-fix copies go to data/_legacy_archive/. Official files are overwritten
only after those copies exist. Does not delete history.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

ROOT = Path("/tmp/pr38_audit/Dual_Target_Docking")
REM = ROOT / "remediation_outputs"
LEGACY = ROOT / "data/_legacy_archive"
BOX40 = ROOT / "data/egfr_her2_panel40_v0/boxes"
P120 = ROOT / "data/egfr_her2_panel120_v0"
ACHE = ROOT / "data/ache_bche_panel_v0"


def copy_keep(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.exists() and not dest.exists():
        shutil.copy2(src, dest)


def write_box_txt(box: dict, path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def adopt_boxes() -> None:
    dest_dir = LEGACY / "egfr_her2_boxes_pre_fix"
    dest_dir.mkdir(parents=True, exist_ok=True)
    (dest_dir / "README.md").write_text(
        "historical pre-audit EGFR/HER2 boxes (hydrogen-inclusive AABB);\n"
        "not used in post-remediation primary analysis.\n",
        encoding="utf-8",
    )
    for name in ("3POZ_box.json", "3RCD_box.json", "3POZ_box.txt", "3RCD_box.txt", "all_boxes.json"):
        copy_keep(BOX40 / name, dest_dir / name)
        copy_keep(BOX40 / name, REM / "legacy_egfr_her2_boxes" / name)

    boxes = {}
    for pdb in ("3POZ", "3RCD"):
        box = json.loads((REM / "phase1_boxes" / f"{pdb}_box_corrected.json").read_text())
        box["status"] = "canonical_post_remediation"
        box["n_ligand_atoms"] = int(box.get("n_heavy_atoms", 38))
        boxes[pdb] = box
        (BOX40 / f"{pdb}_box.json").write_text(json.dumps(box, indent=2) + "\n", encoding="utf-8")
        write_box_txt(box, BOX40 / f"{pdb}_box.txt")
        p120_boxes = P120 / "boxes"
        p120_boxes.mkdir(parents=True, exist_ok=True)
        (p120_boxes / f"{pdb}_box.json").write_text(json.dumps(box, indent=2) + "\n", encoding="utf-8")
    (BOX40 / "all_boxes.json").write_text(json.dumps(boxes, indent=2) + "\n", encoding="utf-8")
    rec_src = ROOT / "data/egfr_her2_panel40_v0/receptors"
    rec_dst = P120 / "receptors"
    rec_dst.mkdir(parents=True, exist_ok=True)
    for pdb in ("3POZ", "3RCD"):
        src = rec_src / f"{pdb}_receptor.pdbqt"
        dst = rec_dst / f"{pdb}_receptor.pdbqt"
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)
        prot = rec_src / f"{pdb}_protein.pdb"
        if prot.exists() and not (rec_dst / f"{pdb}_protein.pdb").exists():
            shutil.copy2(prot, rec_dst / f"{pdb}_protein.pdb")
    lig_src = Path(
        "/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel120_v0/ligands_pdbqt"
    )
    lig_dst = P120 / "ligands_pdbqt"
    if lig_src.exists() and not lig_dst.exists():
        lig_dst.symlink_to(lig_src)
    print("adopted EGFR/HER2 canonical boxes")


def adopt_ache_panel() -> None:
    dest = LEGACY / "ache_bche_panel_pre_fix"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "README.md").write_text(
        "historical AChE/BChE panel built with ChEMBL-ID prefix cap;\n"
        "not used in post-remediation primary analysis.\n",
        encoding="utf-8",
    )
    for name in (
        "panel_v0_strict.csv",
        "panel_v0_strict_with_smiles.csv",
        "ablation_ligand_scores.csv",
    ):
        copy_keep(ACHE / "tables" / name, dest / name)

    new_panel = list(csv.DictReader((REM / "phase2_ache_bche/panel_v0_strict_no_id_prefix_cap.csv").open()))
    fields_strict = [
        "panel_id",
        "molecule_chembl_id",
        "class",
        "pchembl_ACHE",
        "pchembl_BCHE",
        "min_pchembl",
        "label_rule",
        "gray_excluded",
        "smiles",
    ]
    with (ACHE / "tables/panel_v0_strict.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields_strict, extrasaction="ignore")
        w.writeheader()
        for r in new_panel:
            w.writerow(r)
    with (ACHE / "tables/panel_v0_strict_with_smiles.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields_strict, extrasaction="ignore")
        w.writeheader()
        for r in new_panel:
            w.writerow(r)

    old_abl = {r["ligand"]: r for r in csv.DictReader((dest / "ablation_ligand_scores.csv").open())}
    new_sc = {}
    with (REM / "phase2_ache_bche/vina_new_ligand/scores_vina_mode1_new_ligand.csv").open() as fh:
        for r in csv.DictReader(fh):
            new_sc[r["target"]] = float(r["affinity"])
    rows = []
    for r in new_panel:
        pid = r["panel_id"]
        cid = r["molecule_chembl_id"]
        if cid == "CHEMBL3960861":
            ea, eb = new_sc["ACHE"], new_sc["BCHE"]
            rec = {
                "ligand": pid,
                "molecule_chembl_id": cid,
                "class": r["class"],
                "pchembl_ACHE": r["pchembl_ACHE"],
                "pchembl_BCHE": r["pchembl_BCHE"],
                "min_pchembl": r["min_pchembl"],
                "label_rule": r["label_rule"],
                "gray_excluded": r["gray_excluded"],
                "smiles": r["smiles"],
                "vina_ACHE": ea,
                "vina_BCHE": eb,
                "rtm_ACHE": "",
                "rtm_BCHE": "",
                "vina_ACHE_hb": -ea,
                "vina_BCHE_hb": -eb,
                "vina_mean": 0.5 * ((-ea) + (-eb)),
                "vina_min": min(-ea, -eb),
                "rtm_mean": "",
                "rtm_min": "",
                "rtm_ACHE_z": "",
                "rtm_BCHE_z": "",
                "rtm_min_z": "",
                "prep": "rdkit_meeko",
            }
            rows.append(rec)
            continue
        old = old_abl.get(pid)
        if old is None:
            # panel_id may have shifted; match by chembl
            old = next((v for v in old_abl.values() if v.get("molecule_chembl_id") == cid), None)
        if old is None:
            raise SystemExit(f"missing old ablation row for {pid} {cid}")
        rec = dict(old)
        rec["ligand"] = pid
        rec["molecule_chembl_id"] = cid
        rec["class"] = r["class"]
        rec["smiles"] = r["smiles"]
        rows.append(rec)
    fields = list(rows[0].keys())
    with (ACHE / "tables/ablation_ligand_scores.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"adopted AChE panel n={len(rows)}")


def adopt_ache_holdout() -> None:
    ho = ROOT / "data/jcim_holdout_v0/tables"
    dest = LEGACY / "ache_bche_holdout_pre_fix"
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("holdout_panel_HOAB.csv", "holdout_ligand_scores_v1.csv", "scores_vina_mode1_HOAB.csv"):
        copy_keep(ho / name, dest / name)
    new_ho = list(csv.DictReader((REM / "phase2_ache_bche/holdout_panel_HOAB_corrected.csv").open()))
    with (ho / "holdout_panel_HOAB.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(new_ho[0].keys()))
        w.writeheader()
        w.writerows(new_ho)

    old_scores = [r for r in csv.DictReader((dest / "holdout_ligand_scores_v1.csv").open())]
    old_by_chembl = {r["chembl"]: r for r in old_scores if r.get("pair") == "AChE/BChE"}
    new_aff = {}
    npath = REM / "phase2_ache_bche/vina_holdout_new/scores_vina_mode1_new_holdout.csv"
    if npath.exists():
        for r in csv.DictReader(npath.open()):
            new_aff.setdefault(r["molecule_chembl_id"], {})[r["target"]] = float(r["affinity"])
    kept = {r["ligand_id"] for r in csv.DictReader((REM / "phase2_ache_bche/holdout_old_vs_corrected.csv").open()) if r["old_holdout"] == "yes" and r["new_holdout"] == "yes"}
    out_ab = []
    for r in new_ho:
        cid = r["molecule_chembl_id"]
        hid = r["holdout_id"]
        if cid in kept and cid in old_by_chembl:
            rec = dict(old_by_chembl[cid])
            rec["ligand"] = hid
            rec["chembl"] = cid
            rec["cls"] = r["class"]
            rec["smiles"] = r["smiles"]
            out_ab.append(rec)
            continue
        sc = new_aff.get(cid, {})
        ea, eb = sc.get("4EY7"), sc.get("4BDS")
        if ea is None or eb is None:
            print("WARN holdout missing scores", cid, hid)
            continue
        out_ab.append(
            {
                "prefix": "HOAB",
                "pair": "AChE/BChE",
                "ligand": hid,
                "chembl": cid,
                "cls": r["class"],
                "smiles": r["smiles"],
                "vina_A_raw": ea,
                "vina_B_raw": eb,
                "vina_A": -ea,
                "vina_B": -eb,
                "vina_mean": 0.5 * ((-ea) + (-eb)),
                "heavy": "",
                "mw": "",
                "clogp": "",
                "tpsa": "",
            }
        )
    others = [r for r in old_scores if r.get("pair") != "AChE/BChE"]
    all_rows = others + out_ab
    with (ho / "holdout_ligand_scores_v1.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    print(f"adopted AChE holdout n={len(out_ab)}")


def update_egfr_primary_vina() -> None:
    dest = LEGACY / "egfr_her2_scores_pre_fix"
    dest.mkdir(parents=True, exist_ok=True)
    abl = P120 / "tables/ablation_ligand_scores.csv"
    copy_keep(abl, dest / "ablation_ligand_scores.csv")
    scores = {}
    with (REM / "phase1_vina/scores_vina_mode1_corrected_box.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["status"] in {"ok", "cached"}:
                scores.setdefault(r["ligand"], {})[r["target"]] = float(r["affinity"])
    rows = list(csv.DictReader(abl.open()))
    for r in rows:
        lig = r["ligand"]
        ea, eb = scores[lig]["3POZ"], scores[lig]["3RCD"]
        r["3POZ_affinity"] = ea
        r["3RCD_affinity"] = eb
        r["vina_3POZ_hb"] = -ea
        r["vina_3RCD_hb"] = -eb
        r["vina_mean"] = 0.5 * ((-ea) + (-eb))
        r["vina_min"] = min(-ea, -eb)
        # RTM/CNN leftover columns stay until pose rescoring overwrites them
    with abl.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"updated EGFR primary vina n={len(rows)}")


def main() -> int:
    adopt_boxes()
    adopt_ache_panel()
    adopt_ache_holdout()
    update_egfr_primary_vina()
    (LEGACY / "README.md").write_text(
        "Pre-fix official files archived before canonical adoption of\n"
        "corrected EGFR/HER2 heavy-atom boxes and AChE/BChE no-ID-prefix-cap panel.\n"
        "Do not use these files in current primary analysis.\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
