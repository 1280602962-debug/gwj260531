#!/usr/bin/env python3
"""Inventory cognate redock artifacts that can be recomputed from git.

EGFR/HER2 3POZ/3RCD ranked RMSD cannot be recomputed unless local pose files
are recovered. This script records presence/absence; it does not redock.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANALYSIS = ROOT / "data" / "jcim_novelty_v0" / "analysis"
TAB.mkdir(parents=True, exist_ok=True)
ANALYSIS.mkdir(parents=True, exist_ok=True)

SPECS = (
    {
        "pair": "AChE/BChE",
        "target": "AChE",
        "pdb": "4EY7",
        "ligand": "E20",
        "crystal_sdf": "data/ache_bche_panel_v0/cognate_qc/4EY7_E20_crystal.sdf",
        "pose_pdbqt": "data/ache_bche_panel_v0/cognate_qc/4EY7_cognate_out.pdbqt",
        "ranked_recomputable": True,
        "status_in_paper": "recomputed in Table S3 / S38 cognate re-audit",
    },
    {
        "pair": "AChE/BChE",
        "target": "BChE",
        "pdb": "4BDS",
        "ligand": "THA",
        "crystal_sdf": "data/ache_bche_panel_v0/cognate_qc/4BDS_THA_crystal.sdf",
        "pose_pdbqt": "data/ache_bche_panel_v0/cognate_qc/4BDS_cognate_out_E8.pdbqt",
        "ranked_recomputable": True,
        "status_in_paper": "fails top-1, passes top-3",
    },
    {
        "pair": "PIK3CA/PIK3CB",
        "target": "PIK3CB",
        "pdb": "2WXF",
        "ligand": "039",
        "crystal_sdf": "data/pik3ca_pik3cb_panel_v0/cognate_qc/2WXF_039_crystal.sdf",
        "pose_pdbqt": "data/pik3ca_pik3cb_panel_v0/cognate_qc/2WXF_cognate_out_E8.pdbqt",
        "ranked_recomputable": True,
        "status_in_paper": "passes top-1",
    },
    {
        "pair": "PIK3CA/mTOR",
        "target": "PIK3CA",
        "pdb": "4L23",
        "ligand": "X6K",
        "crystal_sdf": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/COGNATE_QC_VERDICT_E16.md",
        "pose_pdbqt": "",
        "ranked_recomputable": False,
        "status_in_paper": "best-of-9 verdict deposited; cognate SDF/PDBQT not in git",
    },
    {
        "pair": "PIK3CA/mTOR",
        "target": "mTOR",
        "pdb": "4JT6",
        "ligand": "PI-103",
        "crystal_sdf": "data/pik3ca_mtor_panel48_v0/analysis/cognate_redock_v0/COGNATE_QC_VERDICT_E16.md",
        "pose_pdbqt": "",
        "ranked_recomputable": False,
        "status_in_paper": "fails top-1, passes top-3 in historical ranked summary; poses not in git",
    },
    {
        "pair": "EGFR/HER2",
        "target": "EGFR",
        "pdb": "3POZ",
        "ligand": "03P",
        "crystal_sdf": "data/egfr_her2_panel40_v0/analysis/exhaustiveness_sensitivity_v1/tables/3POZ_cocrystal_03P.pdb",
        "pose_pdbqt": "",
        "ranked_recomputable": False,
        "status_in_paper": "crystal coordinates present; nine-mode cognate PDBQT absent from git",
    },
    {
        "pair": "EGFR/HER2",
        "target": "HER2",
        "pdb": "3RCD",
        "ligand": "TAK-285",
        "crystal_sdf": "data/egfr_her2_panel40_v0/analysis/exhaustiveness_sensitivity_v1/tables/3RCD_cocrystal_03P.pdb",
        "pose_pdbqt": "",
        "ranked_recomputable": False,
        "status_in_paper": "file name may not be the HER2 cognate; nine-mode PDBQT absent",
    },
)


def exists(relative: str) -> bool:
    return bool(relative) and (ROOT / relative).exists()


def main() -> None:
    extra_egfr = [
        "data/egfr_her2_panel40_v0/receptors/3POZ_receptor.pdbqt",
        "data/egfr_her2_panel40_v0/receptors/3RCD_receptor.pdbqt",
        "data/egfr_her2_panel40_v0/protocol/protocol.yaml",
    ]
    rows = []
    for spec in SPECS:
        crystal_ok = exists(spec["crystal_sdf"])
        pose_ok = exists(spec["pose_pdbqt"])
        action = "recompute_ranked_rmsd_from_git"
        if spec["pair"] == "EGFR/HER2":
            action = "local_recover_or_reconstructed_qc_redock"
        elif not pose_ok:
            action = "local_recover_pose_or_use_historical_verdict"
        rows.append(
            {
                "pair": spec["pair"],
                "target": spec["target"],
                "pdb": spec["pdb"],
                "cognate_ligand": spec["ligand"],
                "crystal_present": int(crystal_ok),
                "crystal_path": spec["crystal_sdf"],
                "pose_pdbqt_present": int(pose_ok),
                "pose_pdbqt_path": spec["pose_pdbqt"],
                "ranked_recomputable_from_git": int(bool(spec["ranked_recomputable"] and crystal_ok and pose_ok)),
                "cloud_complete": int(bool(spec["ranked_recomputable"] and crystal_ok and pose_ok)),
                "required_local_action": action,
                "status_in_paper": spec["status_in_paper"],
            }
        )
    for path in extra_egfr:
        rows.append(
            {
                "pair": "EGFR/HER2",
                "target": "support",
                "pdb": Path(path).name.split("_")[0],
                "cognate_ligand": "",
                "crystal_present": int(exists(path)),
                "crystal_path": path,
                "pose_pdbqt_present": 0,
                "pose_pdbqt_path": "",
                "ranked_recomputable_from_git": 0,
                "cloud_complete": int(exists(path)),
                "required_local_action": "keep; insufficient without nine-mode cognate poses",
                "status_in_paper": "receptor/box support file",
            }
        )
    out = TAB / "cognate_artifact_inventory_v1.csv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    missing = [row for row in rows if row["pair"] == "EGFR/HER2" and not row["ranked_recomputable_from_git"]]
    lines = [
        "# Cognate artifact inventory",
        "",
        "AChE/BChE and PIK3CB ranked RMSD can be recomputed from git.",
        "EGFR/HER2 nine-mode cognate poses are **not** in the repository.",
        "Cloud action stops at this inventory. Local recovery SOP:",
        "`docs/EGFR_HER2_COGNATE_RECOVERY_SOP.md`.",
        "",
        f"EGFR/HER2 missing ranked-recompute files: {len(missing)}.",
        "",
    ]
    (ANALYSIS / "COGNATE_ARTIFACT_INVENTORY.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
