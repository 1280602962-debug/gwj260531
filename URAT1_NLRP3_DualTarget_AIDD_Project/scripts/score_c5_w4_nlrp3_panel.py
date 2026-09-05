#!/usr/bin/env python3
"""C5 W4 scoring: NLRP3 structural gate on the rebuilt panel (0 new docks).

Thresholds are the frozen C1 floors (overlap>=0.50, IFP>=0.50, key>=5/7, no clash).
Do not retune after seeing names or p-values.

Empty / timed-out SDFs count as keep_loose=False and keep_structural=False
(no pose = did not enter the pocket). Also report a complete-case table that
drops the 5 decoys that failed all three seeds.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from c1_nlrp3_pose_metrics import (  # noqa: E402
    crystal_reference_ifp,
    evaluate_nlrp3_structural,
    load_key_map,
    load_receptor_heavy,
)

OUT = PROJECT_ROOT / "data/campaigns/c5/02_nlrp3_panel"
PANEL = PROJECT_ROOT / "data/campaigns/c1/05_metrics/nlrp3_structural_panel/panel_ligands.csv"
SEEDS = [42, 43, 44]
FAILED_DECOYS = {"C5W4D_002", "C5W4D_022", "C5W4D_026", "C5W4D_027", "C5W4D_029"}
POS_IDS = [
    "NP3-146",
    "MCC950",
    "CHEMBL4204644",
    "CHEMBL5219789",
    "CHEMBL4212407",
    "CHEMBL6143743",
    "CHEMBL4209503",
    "CHEMBL4216836",
    "CHEMBL6171925",
]
BG_ROOT = {
    42: PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual/nlrp3_7alv/seed42",
    43: PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_7alv/seed43",
    44: PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual_a2/nlrp3_7alv/seed44",
}
W4_BG_OVERRIDE = {
    (43, "REP_07837"): OUT / "background/seed43/REP_07837_out.sdf",
    (44, "REP_07837"): OUT / "background/seed44/REP_07837_out.sdf",
}


def empty_row(ligand_id: str, seed: int, role: str, sdf: Path, reason: str) -> dict:
    return {
        "ligand_id": ligand_id,
        "role": role,
        "seed": seed,
        "error": reason,
        "keep_nlrp3_pose": False,
        "keep_nlrp3_structural": False,
        "sdf": str(sdf),
    }


def fisher_block(pos: pd.DataFrame, bg: pd.DataFrame, col: str) -> dict:
    a = int(pos[col].sum())
    b = int((~pos[col]).sum())
    c = int(bg[col].sum())
    d = int((~bg[col]).sum())
    oddsr, p = fisher_exact([[a, b], [c, d]])
    or_out = None if oddsr == float("inf") else float(oddsr)
    return {
        "tp": a,
        "fn": b,
        "fp": c,
        "tn": d,
        "positive_n": int(len(pos)),
        "background_n": int(len(bg)),
        "positive_pass_rate": float(pos[col].mean()) if len(pos) else float("nan"),
        "background_pass_rate": float(bg[col].mean()) if len(bg) else float("nan"),
        "odds_ratio": or_out,
        "odds_ratio_note": "infinite because all positives passed (fn=0)" if or_out is None else None,
        "fisher_exact_p": float(p),
    }


def score_one(sdf: Path, ligand_id: str, seed: int, role: str, caches: dict) -> dict:
    if not sdf.exists() or sdf.stat().st_size == 0:
        return empty_row(ligand_id, seed, role, sdf, "empty_or_missing_sdf")
    ev = evaluate_nlrp3_structural(
        sdf,
        ligand_id,
        seed,
        key_map=caches["key_map"],
        ref_heavy=caches["ref_heavy"],
        ref_com=caches["ref_com"],
        receptor_heavy=caches["receptor_heavy"],
        ref_ifp=caches["ref_ifp"],
    )
    ev["role"] = role
    return ev


def main() -> None:
    panel = pd.read_csv(PANEL)
    bg_ids = [
        r
        for r in panel.loc[panel.role == "clinical_acid_background", "ligand_id"].tolist()
    ]
    decoy_ids = [f"C5W4D_{i:03d}" for i in range(1, 41)]

    caches = {
        "key_map": load_key_map(),
        "receptor_heavy": load_receptor_heavy(),
    }
    ref_ifp, ref_heavy, ref_com = crystal_reference_ifp()
    caches["ref_ifp"] = ref_ifp
    caches["ref_heavy"] = ref_heavy
    caches["ref_com"] = ref_com

    rows: list[dict] = []
    for seed in SEEDS:
        for lid in POS_IDS:
            sdf = OUT / "positives" / f"seed{seed}" / f"{lid}_out.sdf"
            rows.append(score_one(sdf, lid, seed, "positive", caches))
        for lid in decoy_ids:
            sdf = OUT / "decoys" / f"seed{seed}" / f"{lid}_out.sdf"
            rows.append(score_one(sdf, lid, seed, "decoy", caches))
        for lid in bg_ids:
            sdf = W4_BG_OVERRIDE.get((seed, lid), BG_ROOT[seed] / f"{lid}_out.sdf")
            rows.append(score_one(sdf, lid, seed, "clinical_acid_background", caches))

    df = pd.DataFrame(rows)
    for col in ("keep_nlrp3_pose", "keep_nlrp3_structural"):
        df[col] = df[col].fillna(False).astype(bool)

    per_path = OUT / "w4_panel_metrics_all_seeds.csv"
    df.to_csv(per_path, index=False)

    summary: dict = {
        "n_rows": int(len(df)),
        "failed_decoys_all_seeds": sorted(FAILED_DECOYS),
        "n_failed_decoy_jobs": 15,
        "thresholds_frozen": {
            "overlap_ge": 0.50,
            "ifp_jaccard_ge": 0.50,
            "key_contacts_ge": 5,
            "note": "copied from c1_nlrp3_pose_metrics; not retuned",
        },
        "empty_sdf_rule": "keep_loose=False and keep_structural=False",
        "by_seed": {},
        "pass_rule": {
            "structural_gate_specificity_gt_loose_gate": None,
            "fisher_p_lt_0.05": None,
        },
    }

    for seed in SEEDS:
        sub = df[df.seed == seed]
        pos = sub[sub.role == "positive"]
        decoy_all = sub[sub.role == "decoy"]
        decoy_cc = decoy_all[~decoy_all.ligand_id.isin(FAILED_DECOYS)]
        bg = sub[sub.role == "clinical_acid_background"]
        block = {
            "positives_vs_decoys_intent_to_score_n40": {
                "loose": fisher_block(pos, decoy_all, "keep_nlrp3_pose"),
                "structural": fisher_block(pos, decoy_all, "keep_nlrp3_structural"),
            },
            "positives_vs_decoys_complete_case_n35": {
                "loose": fisher_block(pos, decoy_cc, "keep_nlrp3_pose"),
                "structural": fisher_block(pos, decoy_cc, "keep_nlrp3_structural"),
            },
            "positives_vs_clinical_background_n20": {
                "loose": fisher_block(pos, bg, "keep_nlrp3_pose"),
                "structural": fisher_block(pos, bg, "keep_nlrp3_structural"),
            },
        }
        for name, pair in block.items():
            loose_fp = pair["loose"]["background_pass_rate"]
            struct_fp = pair["structural"]["background_pass_rate"]
            pair["structural_specificity_gt_loose"] = bool(struct_fp < loose_fp)
            pair["structural_and_loose_identical"] = (
                pair["loose"]["tp"] == pair["structural"]["tp"]
                and pair["loose"]["fp"] == pair["structural"]["fp"]
            )
        summary["by_seed"][str(seed)] = block

    # Primary decision uses seed 42, positives vs decoys, intent-to-score (fail=not-pass).
    primary = summary["by_seed"]["42"]["positives_vs_decoys_intent_to_score_n40"]
    summary["primary"] = {
        "seed": 42,
        "contrast": "positives_vs_decoys_intent_to_score_n40",
        **primary,
    }
    spec_gt = bool(primary["structural_specificity_gt_loose"])
    p_ok = bool(primary["structural"]["fisher_exact_p"] < 0.05)
    summary["pass_rule"]["structural_gate_specificity_gt_loose_gate"] = spec_gt
    summary["pass_rule"]["fisher_p_lt_0.05"] = p_ok
    summary["gate_pass"] = bool(spec_gt and p_ok)
    summary["on_fail"] = (
        None
        if summary["gate_pass"]
        else "report_as_pose_qc_not_structural_compatibility_gate"
    )

    out_json = OUT / "w4_structural_gate_summary.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["primary"], indent=2))
    print("gate_pass", summary["gate_pass"])
    print("wrote", per_path)
    print("wrote", out_json)


if __name__ == "__main__":
    main()
