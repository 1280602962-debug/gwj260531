#!/usr/bin/env python3
"""J0 — Expand strict hard-negative supply audit to ≥40–50 candidate pairs.

Uses cached mols_*.json under data/public_pair_selection/ (no docking).
If a target dictionary is missing, the pair is recorded in j0_fetch_queue.csv
instead of fabricating counts. ChEMBL Web API may be down; this script never
pretends missing targets were audited.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data" / "public_pair_selection"
OUT = Path(__file__).resolve().parents[1] / "tables"

THETA = 6.0
HI = 6.5
LO = 5.5
MIN_HARDNEG_STRICT = 50
MIN_HARDNEG_THIN = 20

# (pair_id, name_A, name_B, family, motivation, literature_anchor, metal_risk, ppi_risk, exclude_note)
# name_* must match mols_{name}.json stem when available.
CANDIDATE_PAIRS = [
    # --- original 12 (controls) ---
    ("P01", "EGFR", "HER2", "kinase_homolog", "classic dual TKI; pose-gold", "TAK-285 3POZ/3RCD", False, False, ""),
    ("P02", "PIK3CA", "MTOR", "pathway", "largest paired ChEMBL set", "PI-103 4L23/4JT6", False, False, ""),
    ("P03", "ACHE", "BCHE", "hydrolase_homolog", "cholinesterase dual", "literature dual AChE/BChE", False, False, ""),
    ("P04", "PIK3CA", "PIK3CB", "isozyme", "isoform control (near)", "PI3Kα/β duals", False, False, "too_close_for_primary"),
    ("P05", "MCL1", "BCL2L1_BclxL", "PPI_BH3", "linked pose-gold", "Tanaka LC6 3WIY/3WIZ", False, True, ""),
    ("P06", "MCL1", "BCL2", "PPI_BH3", "BH3 family", "Mcl-1/Bcl-2 duals", False, True, ""),
    ("P07", "AKT1", "RPS6KB1_p70S6K", "pathway", "AKT/S6K pathway", "M2698-class", False, False, ""),
    ("P08", "BRD4", "HDAC1", "epigenetic_hybrid", "literature hotspot", "JMC HDAC–BET hybrids", True, False, "metal_Zn_HDAC"),
    ("P09", "BRD4", "HDAC6", "epigenetic_hybrid", "literature hotspot", "HDAC6/BET hybrids", True, False, "metal_Zn_HDAC"),
    ("P10", "JAK2", "HDAC1", "epigenetic_hybrid", "JAK/HDAC duals", "JMC JAK/HDAC", True, False, "metal_Zn_HDAC"),
    ("P11", "PARP1", "MET", "synthetic_lethalish", "PARP1/c-Met duals", "JMC PARP1/c-Met", False, False, ""),
    ("P12", "CDK6", "BRD4", "kinase_epigenetic", "CDK6/BRD4 duals", "JMC CDK6/BRD4", False, False, ""),
    # --- expanded among cached targets (literature / pathway motivated) ---
    ("P13", "EGFR", "MET", "kinase_combo", "EGFR/MET resistance duals", "literature EGFR/MET", False, False, ""),
    ("P14", "EGFR", "PIK3CA", "pathway", "EGFR–PI3K axis", "pathway duals", False, False, ""),
    ("P15", "EGFR", "MTOR", "pathway", "EGFR–mTOR axis", "pathway duals", False, False, ""),
    ("P16", "EGFR", "BRD4", "kinase_epigenetic", "EGFR/BRD hybrids", "occasional hybrids", False, False, ""),
    ("P17", "EGFR", "CDK6", "kinase_combo", "cell-cycle + EGFR", "CDK/EGFR literature", False, False, ""),
    ("P18", "EGFR", "PARP1", "DDR_kinase", "DDR + EGFR", "occasional duals", False, False, ""),
    ("P19", "EGFR", "HDAC1", "kinase_epigenetic", "EGFR/HDAC hybrids", "JMC hybrids", True, False, "metal_Zn_HDAC"),
    ("P20", "EGFR", "JAK2", "kinase_combo", "EGFR/JAK axis", "inflammation/oncology", False, False, ""),
    ("P21", "HER2", "PIK3CA", "pathway", "HER2–PI3K", "clinical dual context", False, False, ""),
    ("P22", "HER2", "MTOR", "pathway", "HER2–mTOR", "pathway duals", False, False, ""),
    ("P23", "HER2", "CDK6", "kinase_combo", "HER2/CDK", "cell-cycle", False, False, ""),
    ("P24", "HER2", "HDAC6", "kinase_epigenetic", "HER2/HDAC", "hybrids", True, False, "metal_Zn_HDAC"),
    ("P25", "HER2", "PARP1", "DDR_kinase", "HER2/PARP", "DDR", False, False, ""),
    ("P26", "PIK3CA", "AKT1", "pathway", "PI3K–AKT", "pathway duals", False, False, ""),
    ("P27", "PIK3CB", "MTOR", "pathway", "PI3Kβ–mTOR", "pathway", False, False, ""),
    ("P28", "PIK3CB", "AKT1", "pathway", "PI3Kβ–AKT", "pathway", False, False, ""),
    ("P29", "MTOR", "AKT1", "pathway", "mTOR–AKT", "pathway duals", False, False, ""),
    ("P30", "MTOR", "RPS6KB1_p70S6K", "pathway", "mTOR–S6K", "pathway", False, False, ""),
    ("P31", "AKT1", "MTOR", "pathway", "AKT–mTOR (alias order)", "pathway", False, False, "alias_of_P29"),
    ("P32", "JAK2", "BRD4", "kinase_epigenetic", "JAK/BET", "literature", False, False, ""),
    ("P33", "JAK2", "CDK6", "kinase_combo", "JAK/CDK", "oncology", False, False, ""),
    ("P34", "JAK2", "PIK3CA", "kinase_pathway", "JAK/PI3K", "oncology", False, False, ""),
    ("P35", "CDK6", "PARP1", "DDR_cellcycle", "PARP/CDK6 duals", "JMC PARP/CDK6", False, False, ""),
    ("P36", "CDK6", "HDAC1", "kinase_epigenetic", "CDK/HDAC", "hybrids", True, False, "metal_Zn_HDAC"),
    ("P37", "CDK6", "PIK3CA", "kinase_pathway", "CDK/PI3K", "oncology", False, False, ""),
    ("P38", "BRD4", "PARP1", "epigenetic_DDR", "BET/PARP", "literature", False, False, ""),
    ("P39", "BRD4", "HDAC2", "epigenetic_hybrid", "BET/HDAC2", "HDAC family", True, False, "metal_Zn_HDAC"),
    ("P40", "BRD4", "PIK3CA", "epigenetic_pathway", "BET/PI3K", "hybrids", False, False, ""),
    ("P41", "PARP1", "HDAC1", "DDR_epigenetic", "PARP/HDAC", "hybrids", True, False, "metal_Zn_HDAC"),
    ("P42", "PARP1", "EGFR", "DDR_kinase", "alias order of P18", "DDR", False, False, "alias_of_P18"),
    ("P43", "MET", "EGFR", "kinase_combo", "alias order of P13", "EGFR/MET", False, False, "alias_of_P13"),
    ("P44", "MET", "PIK3CA", "kinase_pathway", "MET–PI3K", "resistance", False, False, ""),
    ("P45", "MET", "AKT1", "kinase_pathway", "MET–AKT", "pathway", False, False, ""),
    ("P46", "HDAC1", "HDAC6", "isozyme", "HDAC isoform control", "selectivity panels", True, False, "metal_both;isozyme"),
    ("P47", "HDAC1", "PIK3CA", "epigenetic_pathway", "HDAC/PI3K hybrids", "literature", True, False, "metal_Zn_HDAC"),
    ("P48", "HDAC6", "PIK3CA", "epigenetic_pathway", "HDAC6/PI3K", "literature", True, False, "metal_Zn_HDAC"),
    ("P49", "BCL2", "BCL2L1_BclxL", "PPI_BH3", "Bcl-2/Bcl-xL", "BH3 duals", False, True, ""),
    ("P50", "MCL1", "AKT1", "PPI_pathway", "apoptosis + AKT", "combo literature", False, True, ""),
    ("P51", "ACHE", "HDAC1", "hydrolase_epigenetic", "rare hybrid", "exploratory", True, False, "metal_Zn_HDAC;weak_lit"),
    ("P52", "BCHE", "EGFR", "cross_fold", "negative-control-ish pair", "exploratory", False, False, "weak_lit"),
    # --- explicit exclusions (still listed for transparency; not audited as main) ---
    ("X01", "NLRP3", "JNK1", "inflammasome", "private holdout", "project holdout", False, False, "EXCLUDED_private_holdout"),
]

# Targets not in cache that would unlock additional literature pairs (fetch queue).
EXTRA_TARGETS_FOR_QUEUE = [
    {"name": "VEGFR2_KDR", "chembl": "CHEMBL279", "uniprot": "P35968", "reason": "VEGFR2/HDAC duals (JMC)"},
    {"name": "AXL", "chembl": "CHEMBL4895", "uniprot": "P30530", "reason": "MER/AXL duals"},
    {"name": "MERTK", "chembl": "CHEMBL3983", "uniprot": "Q12866", "reason": "MER/AXL duals"},
    {"name": "SYK", "chembl": "CHEMBL2599", "uniprot": "P43405", "reason": "SYK/HDAC duals"},
    {"name": "HSP90AA1", "chembl": "CHEMBL3880", "uniprot": "P07900", "reason": "Hsp90/HDAC6 duals"},
    {"name": "WEE1", "chembl": "CHEMBL2421", "uniprot": "P30291", "reason": "Wee1/HDAC duals"},
    {"name": "TOP1", "chembl": "CHEMBL1781", "uniprot": "P11387", "reason": "Top/HDAC duals"},
    {"name": "ROCK1", "chembl": "CHEMBL3231", "uniprot": "Q13464", "reason": "ROCK/HDAC duals"},
    {"name": "PIM1", "chembl": "CHEMBL2147", "uniprot": "P11309", "reason": "PIM/HDAC duals"},
    {"name": "SERT_SLC6A4", "chembl": "CHEMBL228", "uniprot": "P31645", "reason": "SERT/ER duals"},
    {"name": "ESR1", "chembl": "CHEMBL206", "uniprot": "P03372", "reason": "SERT/ER duals"},
    {"name": "FGFR1", "chembl": "CHEMBL3650", "uniprot": "P11362", "reason": "FGFR dual TKIs"},
    {"name": "ALK", "chembl": "CHEMBL4247", "uniprot": "Q9UM73", "reason": "ALK dual TKIs"},
    {"name": "BRAF", "chembl": "CHEMBL5145", "uniprot": "P15056", "reason": "BRAF/MEK pathway pairs"},
    {"name": "MAP2K1_MEK1", "chembl": "CHEMBL3587", "uniprot": "Q02750", "reason": "BRAF/MEK"},
    {"name": "BTK", "chembl": "CHEMBL5251", "uniprot": "Q06187", "reason": "BTK duals"},
    {"name": "FLT3", "chembl": "CHEMBL1974", "uniprot": "P36888", "reason": "FLT3 dual TKIs"},
    {"name": "SRC", "chembl": "CHEMBL267", "uniprot": "P12931", "reason": "SRC dual TKIs"},
    {"name": "ABL1", "chembl": "CHEMBL1862", "uniprot": "P00519", "reason": "BCR-ABL dual context"},
    {"name": "BACE1", "chembl": "CHEMBL4822", "uniprot": "P56817", "reason": "CNS duals / transporter pairs"},
]


def load(target: str):
    path = SRC / f"mols_{target}.json"
    if not path.exists():
        return None
    with path.open() as fh:
        return {k: float(v) for k, v in json.load(fh).items()}


def audit(a_vals: dict, b_vals: dict) -> dict:
    both = set(a_vals) & set(b_vals)
    rec = {
        "n_both_measured": len(both),
        "theta_dual": 0,
        "theta_A_only": 0,
        "theta_B_only": 0,
        "strict_dual": 0,
        "strict_A_only": 0,
        "strict_B_only": 0,
        "strict_neither": 0,
        "gray": 0,
    }
    for mol in both:
        x, y = a_vals[mol], b_vals[mol]
        if x >= THETA and y >= THETA:
            rec["theta_dual"] += 1
        elif x >= THETA:
            rec["theta_A_only"] += 1
        elif y >= THETA:
            rec["theta_B_only"] += 1
        if x >= HI and y >= HI:
            rec["strict_dual"] += 1
        elif x >= HI and y <= LO:
            rec["strict_A_only"] += 1
        elif y >= HI and x <= LO:
            rec["strict_B_only"] += 1
        elif x <= LO and y <= LO:
            rec["strict_neither"] += 1
        else:
            rec["gray"] += 1
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    available = {p.stem.replace("mols_", "") for p in SRC.glob("mols_*.json")}

    cand_rows = []
    supply_rows = []
    fetch_targets = {}
    seen_unordered = set()

    for pair_id, a, b, family, motiv, lit, metal, ppi, note in CANDIDATE_PAIRS:
        key = tuple(sorted([a, b]))
        is_dup = key in seen_unordered and "alias" not in note
        if "alias" not in note:
            seen_unordered.add(key)
        cand_rows.append(
            {
                "pair_id": pair_id,
                "pair": f"{a}/{b}",
                "target_A": a,
                "target_B": b,
                "family": family,
                "motivation": motiv,
                "literature_anchor": lit,
                "metal_enzyme_risk": metal,
                "ppi_risk": ppi,
                "notes": note,
                "dict_A_present": a in available,
                "dict_B_present": b in available,
                "auditable_now": a in available and b in available,
                "unordered_dup": is_dup,
            }
        )
        if a not in available:
            fetch_targets[a] = {
                "name": a,
                "chembl": "",
                "uniprot": "",
                "reason": f"needed_for_pair:{pair_id}",
            }
        if b not in available:
            fetch_targets[b] = {
                "name": b,
                "chembl": "",
                "uniprot": "",
                "reason": f"needed_for_pair:{pair_id}",
            }
        if a not in available or b not in available:
            continue
        if note.startswith("EXCLUDED") or note.startswith("alias_of_"):
            continue
        rec = audit(load(a), load(b))
        rec.update(
            {
                "pair_id": pair_id,
                "pair": f"{a}/{b}",
                "family": family,
                "metal_enzyme_risk": metal,
                "ppi_risk": ppi,
                "notes": note,
            }
        )
        n = rec["n_both_measured"]
        rec["gray_frac"] = round(rec["gray"] / n, 3) if n else None
        rec["min_strict_hardneg"] = min(rec["strict_A_only"], rec["strict_B_only"])
        rec["supports_strict_panel"] = rec["min_strict_hardneg"] >= MIN_HARDNEG_STRICT
        rec["supports_thin_panel"] = rec["min_strict_hardneg"] >= MIN_HARDNEG_THIN
        supply_rows.append(rec)

    for t in EXTRA_TARGETS_FOR_QUEUE:
        if t["name"] not in available:
            fetch_targets[t["name"]] = t

    supply_rows.sort(key=lambda r: (-r["min_strict_hardneg"], -r["n_both_measured"]))

    fields_s = [
        "pair_id",
        "pair",
        "family",
        "n_both_measured",
        "theta_dual",
        "theta_A_only",
        "theta_B_only",
        "strict_dual",
        "strict_A_only",
        "strict_B_only",
        "strict_neither",
        "gray",
        "gray_frac",
        "min_strict_hardneg",
        "supports_strict_panel",
        "supports_thin_panel",
        "metal_enzyme_risk",
        "ppi_risk",
        "notes",
    ]
    path_s = OUT / "j0_strict_label_supply.csv"
    with path_s.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields_s)
        w.writeheader()
        w.writerows(supply_rows)

    path_c = OUT / "j0_candidate_pairs.csv"
    with path_c.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cand_rows[0].keys()))
        w.writeheader()
        w.writerows(cand_rows)

    fetch_rows = list(fetch_targets.values())
    path_f = OUT / "j0_fetch_queue.csv"
    with path_f.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["name", "chembl", "uniprot", "reason"], extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(fetch_rows)

    n_ok = sum(1 for r in supply_rows if r["supports_strict_panel"])
    n_thin = sum(
        1
        for r in supply_rows
        if r["supports_thin_panel"] and not r["supports_strict_panel"]
    )
    n_aud = len(supply_rows)
    print(f"wrote {path_s} audited={n_aud}")
    print(f"wrote {path_c} candidates={len(cand_rows)}")
    print(f"wrote {path_f} fetch_queue={len(fetch_rows)}")
    print(
        f"supports_strict_panel(Y, ≥{MIN_HARDNEG_STRICT} both sides): {n_ok}/{n_aud}"
    )
    print(f"thin_panel_only (≥{MIN_HARDNEG_THIN}, <{MIN_HARDNEG_STRICT}): {n_thin}")
    print("TOP by min_strict_hardneg:")
    for r in supply_rows[:15]:
        flag = (
            "Y"
            if r["supports_strict_panel"]
            else ("T" if r["supports_thin_panel"] else "-")
        )
        print(
            f"  {r['pair']:28s} both={r['n_both_measured']:5d} "
            f"strict A/B={r['strict_A_only']:4d}/{r['strict_B_only']:4d} "
            f"min={r['min_strict_hardneg']:3d} {flag}"
        )


if __name__ == "__main__":
    main()
