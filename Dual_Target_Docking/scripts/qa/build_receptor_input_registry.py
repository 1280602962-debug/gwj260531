#!/usr/bin/env python3
"""Inspect deposited 14-slot receptors. Do not re-prepare anything."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/provenance/receptor_input_registry.csv"
RMSD = ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv"

SLOTS = [
    {
        "pair": "EGFR/HER2",
        "pocket": "A",
        "receptor_pdb": "3POZ",
        "unique_14slot": 1,
        "cognate_ligand": "03P",
        "receptor_pdbqt": "data/egfr_her2_uniform_rdkit_v1/receptors/3POZ_receptor.pdbqt",
        "protein_pdb": "",
        "box_file": "data/egfr_her2_panel120_v0/boxes/3POZ_box_corrected.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "chain A kinase; git blob 2b398c62 from c213c485",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "EGFR/HER2",
        "pocket": "B",
        "receptor_pdb": "3RCD",
        "unique_14slot": 1,
        "cognate_ligand": "03P",
        "receptor_pdbqt": "data/egfr_her2_uniform_rdkit_v1/receptors/3RCD_receptor.pdbqt",
        "protein_pdb": "",
        "box_file": "data/egfr_her2_panel120_v0/boxes/3RCD_box_corrected.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "chain A kinase; git blob c338094b from c213c485",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "PIK3CA/mTOR",
        "pocket": "A",
        "receptor_pdb": "4L23",
        "unique_14slot": 1,
        "cognate_ligand": "X6K",
        "receptor_pdbqt": "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_receptor.pdbqt",
        "protein_pdb": "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_protein.pdb",
        "box_file": "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4L23_box.json",
        "box_definition": "deposited_box_json",
        "construct_note": "chain A; protein PDB deposited",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "PIK3CA/mTOR",
        "pocket": "B",
        "receptor_pdb": "4JT6",
        "unique_14slot": 1,
        "cognate_ligand": "X6K",
        "receptor_pdbqt": "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_receptor.pdbqt",
        "protein_pdb": "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_protein.pdb",
        "box_file": "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4JT6_box.json",
        "box_definition": "deposited_box_json",
        "construct_note": "chain A; protein PDB deposited",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "AChE/BChE",
        "pocket": "A",
        "receptor_pdb": "4EY7",
        "unique_14slot": 1,
        "cognate_ligand": "E20",
        "receptor_pdbqt": "data/ache_bche_panel_v0/receptors/ACHE_receptor.pdbqt",
        "protein_pdb": "data/ache_bche_panel_v0/receptors/ACHE_protein.pdb",
        "box_file": "data/ache_bche_panel_v0/boxes/4EY7_box.json",
        "box_definition": "deposited_box_json",
        "construct_note": "dimer AB; alias 4EY7_receptor.pdbqt",
        "prep_script": "data/ache_bche_panel_v0/scripts/freeze_receptors_cognate_qc_v2.py",
        "prep_log": "",
    },
    {
        "pair": "AChE/BChE",
        "pocket": "B",
        "receptor_pdb": "4BDS",
        "unique_14slot": 1,
        "cognate_ligand": "THA",
        "receptor_pdbqt": "data/ache_bche_panel_v0/receptors/BCHE_receptor.pdbqt",
        "protein_pdb": "data/ache_bche_panel_v0/receptors/BCHE_protein.pdb",
        "box_file": "data/ache_bche_panel_v0/boxes/4BDS_box.json",
        "box_definition": "deposited_box_json",
        "construct_note": "chain A; alias 4BDS_receptor.pdbqt",
        "prep_script": "data/ache_bche_panel_v0/scripts/freeze_receptors_cognate_qc_v2.py",
        "prep_log": "",
    },
    {
        "pair": "F2/F10",
        "pocket": "A",
        "receptor_pdb": "4UDW",
        "unique_14slot": 1,
        "cognate_ligand": "N6L",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/4UDW_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/4UDW_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/4UDW_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "thrombin chains H/I/L; yaml site thrombin S1 chain H; pocket_10.0.pdb is a pocket extract not a full protein PDB",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "F2/F10",
        "pocket": "B",
        "receptor_pdb": "2JKH",
        "unique_14slot": 1,
        "cognate_ligand": "BI7",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/2JKH_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/2JKH_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/2JKH_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "chains A/L; SI notes CCD BI7 after OpenBabel SDF issue (cognate QC, not receptor rebuild)",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "JAK1/JAK2",
        "pocket": "A",
        "receptor_pdb": "6N7A",
        "unique_14slot": 1,
        "cognate_ligand": "KEV",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6N7A_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6N7A_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6N7A_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "shared JAK1 with JAK1/TYK2; yaml JH1 ATP",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "JAK1/JAK2",
        "pocket": "B",
        "receptor_pdb": "8BXH",
        "unique_14slot": 1,
        "cognate_ligand": "C87",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/8BXH_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/8BXH_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/8BXH_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "yaml JH1 ATP",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "JAK1/TYK2",
        "pocket": "A",
        "receptor_pdb": "6N7A",
        "unique_14slot": 0,
        "cognate_ligand": "KEV",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6N7A_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6N7A_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6N7A_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "same file as JAK1/JAK2 pocket A",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "JAK1/TYK2",
        "pocket": "B",
        "receptor_pdb": "3LXP",
        "unique_14slot": 1,
        "cognate_ligand": "IZA",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/3LXP_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/3LXP_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/3LXP_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "yaml JH1 ATP not JH2",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "PPARG/PPARA",
        "pocket": "A",
        "receptor_pdb": "9V8H",
        "unique_14slot": 1,
        "cognate_ligand": "BRL",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/9V8H_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/9V8H_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/9V8H_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "yaml: LBD ternary keep PG08-NL peptide",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "PPARG/PPARA",
        "pocket": "B",
        "receptor_pdb": "6LXA",
        "unique_14slot": 1,
        "cognate_ligand": "EPA",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6LXA_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6LXA_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6LXA_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "shared PPARA with PPARA/PPARD",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "PPARA/PPARD",
        "pocket": "A",
        "receptor_pdb": "6LXA",
        "unique_14slot": 0,
        "cognate_ligand": "EPA",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6LXA_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/6LXA_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/6LXA_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "same file as PPARG/PPARA pocket B",
        "prep_script": "",
        "prep_log": "",
    },
    {
        "pair": "PPARA/PPARD",
        "pocket": "B",
        "receptor_pdb": "5U3Q",
        "unique_14slot": 1,
        "cognate_ligand": "7UJ",
        "receptor_pdbqt": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/5U3Q_receptor.pdbqt",
        "protein_pdb": "data/jcim_chembl_universe_v0/local_track_b_v0/receptors/5U3Q_pocket_10.0.pdb",
        "box_file": "data/jcim_chembl_universe_v0/local_track_b_v0/boxes/5U3Q_box.json",
        "box_definition": "cognate_heavy_atom_AABB+5A_min20",
        "construct_note": "yaml: LBD agonist 1 not PEG",
        "prep_script": "",
        "prep_log": "",
    },
]


def evidence(path: str) -> str:
    if not path:
        return "not_recoverable"
    p = ROOT / path
    if p.is_file() and p.stat().st_size > 0:
        return "confirmed"
    return "not_recoverable"


def inspect_pdbqt(path: Path) -> dict:
    chains = []
    resnames = []
    waters = 0
    n_atom = 0
    n_het = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        n_atom += 1
        chain = line[21].strip() if len(line) > 21 else ""
        resn = line[17:20].strip() if len(line) > 20 else ""
        if chain:
            chains.append(chain)
        if resn:
            resnames.append(resn)
        if resn in {"HOH", "WAT", "H2O", "DOD"}:
            waters += 1
        if line.startswith("HETATM"):
            n_het += 1
    chain_str = "".join(sorted(set(chains))) or "not_recoverable"
    counts = Counter(resnames)
    extra = [f"{k}:{v}" for k, v in counts.most_common() if k not in {"HOH", "WAT"} and v >= 5]
    peptide = "none_in_pdbqt"
    if path.name.startswith("9V8H"):
        peptide = "chain_B_peptide_present_in_pdbqt" if "B" in chain_str else "yaml_keep_peptide_not_confirmed_in_file"
    elif path.name.startswith("4UDW"):
        peptide = "thrombin_HIL_construct"
    elif path.name.startswith("ACHE") or "4EY7" in path.name:
        peptide = "dimer_AB"
    return {
        "chains": chain_str,
        "n_atom_records": n_atom,
        "n_hetatm_records": n_het,
        "waters_in_pdbqt": "none_in_pdbqt" if waters == 0 else f"{waters}_water_atoms",
        "residue_summary": ";".join(extra[:8]),
        "cofactor_peptide": peptide,
    }


def load_rmsd() -> dict:
    out = {}
    if not RMSD.is_file():
        return out
    with RMSD.open(encoding="utf-8-sig", newline="") as handle:
        for r in csv.DictReader(handle):
            out[r["pdb"]] = r
    return out


def status_three(ok: bool) -> str:
    return "confirmed" if ok else "not_recoverable"


def main() -> int:
    rmsd = load_rmsd()
    rows = []
    for spec in SLOTS:
        pdbqt = ROOT / spec["receptor_pdbqt"]
        protein = ROOT / spec["protein_pdb"] if spec["protein_pdb"] else None
        box = ROOT / spec["box_file"]
        rec = inspect_pdbqt(pdbqt) if pdbqt.is_file() else {
            "chains": "not_recoverable",
            "n_atom_records": 0,
            "n_hetatm_records": 0,
            "waters_in_pdbqt": "not_recoverable",
            "residue_summary": "",
            "cofactor_peptide": "not_recoverable",
        }
        box_meta = {}
        if box.is_file():
            try:
                box_meta = json.loads(box.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                box_meta = {}
        rms = rmsd.get(spec["receptor_pdb"], {})
        protein_kind = "not_recoverable"
        if protein is not None and protein.is_file():
            protein_kind = "confirmed" if protein.name.endswith("_protein.pdb") else "metadata_only"
        prep_method = "not_recoverable"
        prep_evidence = "not_recoverable"
        prep_script_status = evidence(spec["prep_script"]) if spec["prep_script"] else "not_recoverable"
        if spec["receptor_pdb"] == "4EY7" and prep_script_status == "confirmed":
            # Recovered script uses mk_prepare_receptor for 4EY7, but also lists BCHE as 6ZWI not 4BDS.
            # No log that the deposited ACHE_receptor.pdbqt is the output of this script.
            prep_method = "meeko_mk_prepare_receptor_recipe_in_recovered_script"
            prep_evidence = "metadata_only"
        elif spec["receptor_pdb"] == "4BDS":
            prep_method = "not_recoverable"
            prep_evidence = "not_recoverable"
            prep_script_status = "metadata_only"
        rows.append(
            {
                "pair": spec["pair"],
                "pocket": spec["pocket"],
                "receptor_pdb": spec["receptor_pdb"],
                "unique_14slot": spec["unique_14slot"],
                "receptor_chain_or_construct": rec["chains"],
                "construct_note": spec["construct_note"],
                "cognate_ligand": spec["cognate_ligand"],
                "cognate_evidence": "confirmed" if rms else "metadata_only",
                "receptor_pdbqt": spec["receptor_pdbqt"],
                "receptor_pdbqt_exists": int(pdbqt.is_file()),
                "receptor_file_status": status_three(pdbqt.is_file() and pdbqt.stat().st_size > 0),
                "protein_pdb": spec["protein_pdb"] or "not_deposited",
                "protein_pdb_status": protein_kind,
                "box_file": spec["box_file"],
                "box_file_status": status_three(box.is_file()),
                "box_definition": spec["box_definition"],
                "box_construction_field": box_meta.get("construction") or box_meta.get("status") or "metadata_only",
                "waters": rec["waters_in_pdbqt"],
                "waters_evidence": "confirmed" if pdbqt.is_file() else "not_recoverable",
                "cofactor_peptide": rec["cofactor_peptide"],
                "cofactor_peptide_evidence": "confirmed" if pdbqt.is_file() else "not_recoverable",
                "n_atom_records": rec["n_atom_records"],
                "residue_summary": rec["residue_summary"],
                "calcrrms_top1_A": rms.get("calcrrms_top1_A", ""),
                "calcrrms_best_A": rms.get("calcrrms_best_A", ""),
                "protein_preparation_method": prep_method,
                "protein_preparation_evidence": prep_evidence,
                "recovered_prep_script": spec["prep_script"] or "not_recoverable",
                "recovered_prep_script_status": prep_script_status,
                "prep_log": spec["prep_log"] or "not_recoverable",
                "prep_log_status": "not_recoverable",
                "protonation": "not_recoverable",
                "missing_loops": "not_recoverable",
                "notes": "deposited PDBQT is the receptor; receptors were not re-prepared",
            }
        )
    fields = list(rows[0].keys())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT, "rows", len(rows), "unique", sum(int(r["unique_14slot"]) for r in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
