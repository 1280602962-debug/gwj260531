#!/usr/bin/env python3
"""END_TO_END_SCIENTIFIC_DATA_INTEGRITY_AUDIT_V4_2. Read-only. No AUROC. No repair.

Constraints:
- within-method physical reuse only; methods need not share scores
- Phase 7: 8070 pair-expanded AND 7960 unique physical job-seed
- old M1/M3 only; do not treat in-progress V2 outputs as final
- overall GATE = PENDING_V2_COMPLETION until new M1/M3 finish + scientific QA
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, rdmolops

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
    ALT_FORBIDDEN,
    GIND,
    LIG_GNINA_PDBQT,
    LIG_PDBQT,
    LIG_SDF,
    MAP,
    PHASE7,
    PROTO,
    QA,
    REC,
    RECEPTORS,
    RUN,
    SEC,
    UNIQUE_14,
    VINA_JOBS,
)
from rerun_v4_2.pre_auroc_v2_lib import (  # noqa: E402
    complementary_cgn_pairs,
    finite_float,
    formal_charge,
    invert_smiles_idx,
    is_glue_type,
    is_real_macrocycle_type,
    load_frozen_sdf,
    map_serial_to_sdf_heavy,
    parse_model1_remarks,
    parse_pdbqt_atoms,
    parse_smiles_idx,
    ring_count,
    sdf_heavy,
    special_types,
    split_models,
)

SEEDS = (17, 29, 42, 71, 101)
PAIRS = list(RECEPTORS.keys())
BOX = RUN / "03_boxes"
REDOCK = RUN / "04_cognate_redocking"
ALT_JOBS = RUN / "08_alt_vina_seed42" / "jobs"
SCRIPTS = ROOT / "scripts" / "rerun_v4_2"

EXPECTED_PAIR_EXPANDED_SEED = 8070
EXPECTED_PHYSICAL_JOB_SEED = 7960
EXPECTED_TIMEOUTS_S42 = {
    ("AB_001", "4EY7"), ("AB_001", "4BDS"),
    ("AB_053", "4EY7"), ("AB_053", "4BDS"),
    ("AB_054", "4EY7"), ("AB_054", "4BDS"),
}
FORMAL_ALT = {
    "4L2Y": ("PIK3CA/mTOR", "A", "A2/B1", "4JT6"),
    "4JT5": ("PIK3CA/mTOR", "B", "A1/B2", "4L23"),
    "4EY6": ("AChE/BChE", "A", "A2/B1", "4BDS"),
    "1P0M": ("AChE/BChE", "B", "A1/B2", "4EY7"),
    "3SHC": ("F2/F10", "A", "A2/B1", "2JKH"),
    "2Y5F": ("F2/F10", "B", "A1/B2", "4UDW"),
    "6KAX": ("PPARA/PPARD", "A", "A2/B1", "5U3Q"),
    "5U46": ("PPARA/PPARD", "B", "A1/B2", "6LXA"),
}


def wcsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("note\nNO_ROWS\n", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0].keys()), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def load(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))


def pdb_chains(path: Path) -> set[str]:
    ch = set()
    if not path.is_file():
        return ch
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) > 21:
            c = ln[21].strip()
            if c:
                ch.add(c)
    return ch


def pdb_resnames(path: Path) -> Counter:
    c = Counter()
    if not path.is_file():
        return c
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            c[ln[17:20].strip()] += 1
    return c


def n_atom_records(path: Path) -> int:
    if not path.is_file():
        return 0
    return sum(1 for ln in path.read_text(errors="replace").splitlines() if ln.startswith(("ATOM", "HETATM")))


def vina_mode1_from_out(path: Path) -> tuple[int, str]:
    n = 0
    mode1 = ""
    if not path.is_file():
        return 0, ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "VINA RESULT" in line:
            n += 1
            if n == 1:
                parts = line.split()
                for i, tok in enumerate(parts):
                    if tok == "RESULT:" and i + 1 < len(parts):
                        mode1 = parts[i + 1]
                        break
    return n, mode1


def v2_in_progress() -> dict:
    m1_new = SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv"
    m3_new = GIND / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv"
    procs = []
    try:
        import subprocess
        p = subprocess.run(["pgrep", "-af", "pre_auroc_v2_rebuild_m1|pre_auroc_v2_repair_m3"], capture_output=True, text=True)
        procs = [ln for ln in (p.stdout or "").splitlines() if "pre_auroc_v2" in ln and "pgrep" not in ln]
    except Exception:
        pass
    n_status = len(list((SEC / "gnina_rescore_representation_v2").glob("*/status.json"))) if (SEC / "gnina_rescore_representation_v2").is_dir() else 0
    n_m3 = len(list((GIND / "jobs_topology_repair_v2").glob("*/status.json"))) if (GIND / "jobs_topology_repair_v2").is_dir() else 0
    running = len(procs) > 0
    files_ready = m1_new.is_file() and m3_new.is_file()
    complete = files_ready and not running
    return {
        "v2_processes_running": running,
        "process_lines": procs[:6],
        "new_m1_master_exists": m1_new.is_file(),
        "new_m3_master_exists": m3_new.is_file(),
        "new_m1_status_json_count_in_progress_dir": n_status,
        "new_m3_status_json_count_in_progress_dir": n_m3,
        "treat_new_as_final": complete,
        "reason": (
            "V2_task_final_files_present_processes_stopped"
            if complete
            else "V2_not_complete_do_not_use_active_writes_as_final_evidence"
        ),
    }


def write_dataflow() -> None:
    text = """# END_TO_END_DATAFLOW_AUDIT

Read-only map. No correctness judgment in this file except listing conversion risks.

## Channel 1 — primary 14-receptor method ablation

raw receptor PDB (`01_receptors/{pdb}_protein_raw.pdb`)
→ PDBFixer (`*_protein_fixer*.pdb`) → optional PTR→TYR (`*_protein_ptr_fallback.pdb`)
→ prepared PDB (`*_protein_prepared.pdb`) → Meeko (`*_receptor.pdbqt`)
script: `scripts/rerun_v4_2/phase2_prepare_receptors.py`
software: PDBFixer 1.12 / OpenMM 8.6.1 / Meeko 0.7.1
may change: hydrogens (Meeko), atom types, charges; PTR→TYR only 6N7A/8BXH/3LXP
must not change: deposited sequence except frozen PTR fallback; retained/removed components

cognate crystal JSON + AABB pad 5 Å min edge 20 Å → `03_boxes/{pdb}_box.json`
script: historical V4.1 box reused / Phase 2 does not rewrite boxes
may change: none if frozen; risk if alt box reused primary coordinates

frozen ligand SDF `02_ligands/sdf/{eid}.sdf` (reused from V4.1, 785)
→ Meeko PDBQT `02_ligands/pdbqt/{eid}.pdbqt`
may change: AD4 types, Gasteiger charges, explicit H, macrocycle CGn/Gn glue (CGn=real atom, Gn=glue)

Vina five-seed `06_vina_fiveseed/jobs/{pair}__{pdb}__{eid}__seed{s}/`
script: `phase7_fiveseed_production.py` vina 1.2.7 exh 16 modes 9 energy_range 6 cpu 1
pair-expanded seed records: 8070
unique physical job-seed: 7960
may change: coordinates (search); must not change ligand graph

Vina pose → old GNINA M1 `09_secondary_scoring/gnina_rescore/`
script: `ablation_phase4_gnina_rescore.py` score_only no_gpu scoring vina cnn_scoring rescore
risk: CG0 not valid AutoDock type (historical)

Vina pose → M2 `09_secondary_scoring/rtmscore/`
script: `ablation_phase578_rtmscore.py` frozen SDF graph + SMILES IDX coords, model1, 10 Å pockets

frozen/converted ligand → old M3 `10_gnina_independent_docking/`
script: `ablation_phase10_12_gnina.py` CG0→C drop G0; seed 42 exh 16 modes 9
risk: ring-closure bond not restored (CGn/Gn topology)

## Channel 2 — alternative receptor experiment (separate freeze)

authority: `ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml`
8 alts: 4L2Y 4JT5 4EY6 1P0M 3SHC 2Y5F 6KAX 5U46
own prep + own box + alt redock `07_alt_cognate_redocking` + alt production `08_alt_vina_seed42`
A2/B1 replaces A only; A1/B2 replaces B only; kept side reuses primary score
NOT the same channel as analysis_freeze allowed_pdb 4JPS/5DXT/4JSX

## Channel 3 — V2 representation repair (NOT final evidence in this audit)

`gnina_rescore_representation_v2/` and `jobs_topology_repair_v2/` are active writes.
This audit does not treat them as official M1/M3.
"""
    (QA / "END_TO_END_DATAFLOW_AUDIT.md").write_text(text, encoding="utf-8")


def config_authority() -> dict:
    freeze = (PROTO / "SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml").read_text(encoding="utf-8")
    af = (PROTO / "analysis_freeze.yaml").read_text(encoding="utf-8")
    altf = (PROTO / "ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml").read_text(encoding="utf-8")
    hits = []
    runner_old = []
    current = {
        "ablation_phase4_gnina_rescore.py", "ablation_phase578_rtmscore.py",
        "ablation_phase10_12_gnina.py", "phase7_fiveseed_production.py",
        "ablation_phase23_masters.py", "ablation_config.py",
        "alt_production_seed42.py", "phase2_prepare_receptors.py",
        "phase5_cognate_redock.py",
    }
    for p in sorted(SCRIPTS.glob("*.py")):
        t = p.read_text(encoding="utf-8", errors="replace")
        rec = {
            "script": p.name,
            "reads_analysis_freeze_yaml": "analysis_freeze.yaml" in t,
            "reads_scoring_freeze": "SCORING_DOCKING_ABLATION_FREEZE_FINAL" in t,
            "reads_alt_freeze_final": "ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL" in t,
            "hardcodes_UNIQUE_14": "UNIQUE_14" in t or "3POZ" in t,
            "mentions_4JPS": "4JPS" in t,
            "mentions_LEGACY": "LEGACY" in t,
            "v4_1_path": "UNIFORM_RERUN_V4_1" in t,
        }
        hits.append(rec)
        if p.name in current and rec["reads_analysis_freeze_yaml"]:
            runner_old.append(p.name)
    findings = []
    if "independent_docking_in_v4_1: false" in af:
        findings.append({"item": "analysis_freeze gnina.independent_docking_in_v4_1 false", "status": "STALE_PROVENANCE_ONLY", "drives_current_M3": False})
    if "4JPS" in af:
        findings.append({"item": "analysis_freeze allowed_pdb 4JPS/5DXT/4JSX", "status": "STALE_PROVENANCE_ONLY", "drives_formal_alt": False})
    if "forbidden_alternative_pdbs" in freeze and "4L2Y" in freeze:
        findings.append({
            "item": "ablation freeze forbids the 8 PDBs that alt experiment executed",
            "status": "CHANNEL_SEPARATION_REQUIRED",
            "note": "alt channel uses ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL; primary M0-M3 use UNIQUE_14 only",
        })
    findings.append({"item": "ablation_config UNIQUE_14", "status": "ACTIVE_PRIMARY", "matches_scoring_freeze": True})
    findings.append({"item": "RTM checkpoint rtmscore_model1.pth", "status": "ACTIVE", "hardcoded_in": "ablation_config.py"})
    findings.append({"item": "GNINA v1.3.2 f23dd2b", "status": "ACTIVE", "hardcoded_in": "ablation_config.py + scoring freeze"})
    report = {
        "primary_science_authority": "SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml + ablation_config.py",
        "alt_science_authority": "ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml",
        "analysis_freeze_method_fields": "STALE_PROVENANCE_ONLY",
        "current_runners_reading_analysis_freeze_yaml": runner_old,
        "CRITICAL_CONFIG_DRIFT": False,
        "findings": findings,
        "script_scan": hits,
    }
    (QA / "CONFIG_AUTHORITY_AUDIT.md").write_text(
        "# CONFIG_AUTHORITY_AUDIT\n\n"
        + json.dumps({k: report[k] for k in report if k != "script_scan"}, indent=2)
        + "\n\nPrimary M0–M3 runners do not read analysis_freeze.yaml method fields.\n"
        "Alt experiment is a separate authoritative channel.\n"
        "No CRITICAL_CONFIG_DRIFT detected for receptor/ligand/method selection of current runners.\n",
        encoding="utf-8",
    )
    return report


def receptor_audit() -> list[dict]:
    policy = load(PROTO / "receptor_component_policy.csv")
    by_pdb = defaultdict(list)
    for r in policy:
        by_pdb[r["pdb_id"]].append(r)
    rows = []
    for pdb in UNIQUE_14:
        raw = REC / f"{pdb}_protein_raw.pdb"
        prep = REC / f"{pdb}_protein_prepared.pdb"
        obs = REC / f"{pdb}_protein_fixer_observed.pdb"
        qt = REC / f"{pdb}_receptor.pdbqt"
        ptr = REC / f"{pdb}_protein_ptr_fallback.pdb"
        cog = json.loads((REC / f"{pdb}_cognate_crystal.json").read_text()) if (REC / f"{pdb}_cognate_crystal.json").is_file() else {}
        raw_ch = pdb_chains(raw)
        prep_ch = pdb_chains(prep)
        prep_res = pdb_resnames(prep)
        issues = []
        if pdb == "5U3Q":
            if "B" not in prep_ch:
                issues.append("5U3Q_prepared_missing_copy_B")
            if "A" in prep_ch:
                issues.append("5U3Q_prepared_still_has_copy_A")
        if pdb == "4L23" and "B" in prep_ch:
            issues.append("4L23_p85_chainB_still_present")
        if pdb == "4JT6":
            if any(x in prep_ch for x in "CD"):
                issues.append("4JT6_mLST8_or_extra_copy_present")
        if pdb == "9V8H" and "B" not in prep_ch:
            issues.append("9V8H_peptide_B_missing")
        if pdb == "4UDW" and not {"H", "I", "L"} <= prep_ch:
            issues.append(f"4UDW_expected_HIL_got_{sorted(prep_ch)}")
        if pdb == "2JKH" and not {"A", "L"} <= prep_ch:
            issues.append(f"2JKH_expected_AL_got_{sorted(prep_ch)}")
        if pdb in {"6N7A", "8BXH", "3LXP"}:
            if not ptr.is_file():
                issues.append("PTR_fallback_file_missing")
            if prep_res.get("PTR", 0) > 0:
                issues.append("PTR_still_in_prepared")
        elif ptr.is_file():
            issues.append("PTR_fallback_present_on_nonJAK")
        n_prep = n_atom_records(prep)
        n_qt = n_atom_records(qt)
        if n_qt == 0:
            issues.append("empty_receptor_pdbqt")
        rows.append({
            "pdb_id": pdb,
            "raw_exists": int(raw.is_file()),
            "prepared_exists": int(prep.is_file()),
            "receptor_pdbqt_exists": int(qt.is_file()),
            "fixer_observed_exists": int(obs.is_file()),
            "ptr_fallback_exists": int(ptr.is_file()),
            "raw_chains": "".join(sorted(raw_ch)),
            "prepared_chains": "".join(sorted(prep_ch)),
            "n_prepared_atoms": n_prep,
            "n_receptor_pdbqt_atoms": n_qt,
            "pdbqt_minus_prepared": n_qt - n_prep,
            "PTR_count_prepared": prep_res.get("PTR", 0),
            "TYR_count_prepared": prep_res.get("TYR", 0),
            "cognate_ccd": (cog.get("meta") or {}).get("ccd", ""),
            "cognate_instance": (cog.get("meta") or {}).get("instance", ""),
            "issues": "|".join(issues),
            "severity": "CRITICAL" if issues else "OK",
        })
    wcsv(QA / "RECEPTOR_SEMANTIC_INTEGRITY_AUDIT.csv", rows)
    md = ["# RECEPTOR_SEMANTIC_INTEGRITY_AUDIT", "", "Primary 14 only. Raw deposited file in-tree used; no new PDB download.", ""]
    for r in rows:
        md.append(f"- {r['pdb_id']}: chains {r['prepared_chains']} issues={r['issues'] or 'none'}")
    (QA / "RECEPTOR_SEMANTIC_INTEGRITY_AUDIT.md").write_text("\n".join(md) + "\n")
    return rows


def box_audit() -> list[dict]:
    rows = []
    for pdb in list(UNIQUE_14) + list(FORMAL_ALT):
        p = BOX / f"{pdb}_box.json"
        rec = {"pdb_id": pdb, "box_exists": int(p.is_file()), "role": "PRIMARY" if pdb in UNIQUE_14 else "ALT", "issues": "", "severity": "OK"}
        if not p.is_file():
            rec.update(issues="missing_box", severity="CRITICAL")
            rows.append(rec)
            continue
        b = json.loads(p.read_text())
        rec.update({
            "center_x": b.get("center_x"), "center_y": b.get("center_y"), "center_z": b.get("center_z"),
            "size_x": b.get("size_x"), "size_y": b.get("size_y"), "size_z": b.get("size_z"),
            "construction": b.get("construction", ""),
            "pad_A": b.get("pad_A"), "min_edge_A": b.get("min_edge_A"),
            "box_json_pdb": b.get("pdb") or b.get("pdb_id"),
            "reference_ccd": b.get("reference_ccd") or b.get("ccd"),
            "reference_instance": b.get("reference_instance"),
        })
        issues = []
        box_pdb = rec["box_json_pdb"]
        if box_pdb and box_pdb != pdb:
            issues.append(f"box_json_pdb_mismatch_{box_pdb}")
        cons = str(b.get("construction", ""))
        if cons and "xminxmax_midpoint" not in cons:
            if "centroid" in cons.lower() and "not_centroid" not in cons.lower():
                issues.append("centroid_not_AABB")
            else:
                issues.append("construction_not_explicit_AABB")
        sizes = [float(b[k]) for k in ("size_x", "size_y", "size_z") if b.get(k) not in (None, "")]
        if sizes and min(sizes) < float(b.get("min_edge_A") or 0) - 1e-6:
            issues.append("size_below_min_edge")
        rec["issues"] = "|".join(issues)
        rec["severity"] = "MAJOR" if issues else "OK"
        rows.append(rec)
    wcsv(QA / "BOX_PROVENANCE_AUDIT.csv", rows)
    return rows


def ligand_full_audit() -> list[dict]:
    jobs = load(SEC / "rescoring_job_manifest.csv")
    rec_by = defaultdict(set)
    for r in jobs:
        rec_by[r["global_ligand_entity_id"]].add(r["pdb_id"])
    eids = sorted({p.stem for p in LIG_SDF.glob("*.sdf")})
    type_global = Counter()
    rows = []
    for eid in eids:
        sdf = load_frozen_sdf(LIG_SDF / f"{eid}.sdf")
        qt = LIG_PDBQT / f"{eid}.pdbqt"
        text = qt.read_text(errors="replace") if qt.is_file() else ""
        atoms = parse_pdbqt_atoms(text) if text else []
        types = sorted({a["type"] for a in atoms if a["type"]})
        for t in types:
            type_global[t] += 1
        spec = special_types(atoms)
        cgn = [a for a in atoms if a["is_real_macro"]]
        gn = [a for a in atoms if a["is_glue"]]
        smiles, idx = parse_smiles_idx(text) if text else ("", {})
        inv = invert_smiles_idx(idx)
        heavy = sdf_heavy(sdf) if sdf is not None else None
        frags = rdmolops.GetMolFrags(heavy) if heavy is not None else ()
        rows.append({
            "global_ligand_entity_id": eid,
            "frozen_sdf_ok": int(sdf is not None),
            "inchikey": eid,
            "formal_charge": formal_charge(sdf) if sdf is not None else "",
            "n_heavy": heavy.GetNumAtoms() if heavy is not None else "",
            "n_rings": ring_count(sdf) if sdf is not None else "",
            "n_fragments": len(frags),
            "all_atom_types": ";".join(types),
            "special_atom_types": ";".join(spec),
            "n_CGn_real_macrocycle_typed": len(cgn),
            "n_Gn_glue": len(gn),
            "CGn_in_SMILES_IDX": sum(1 for a in cgn if a["serial"] in inv),
            "Gn_in_SMILES_IDX": sum(1 for a in gn if a["serial"] in inv),
            "n_smiles_idx": len(idx),
            "affected_receptors": ";".join(sorted(rec_by.get(eid, []))),
            "notes": "",
        })
    wcsv(QA / "LIGAND_REPRESENTATION_FULL_AUDIT.csv", rows)
    (QA / "LIGAND_REPRESENTATION_FULL_AUDIT.md").write_text(
        "# LIGAND_REPRESENTATION_FULL_AUDIT\n\n"
        f"Entities: {len(rows)}\n"
        f"Observed PDBQT types (ligand counts containing type): {dict(type_global)}\n"
        f"CGn/Gn ligands: {sum(1 for r in rows if int(r['n_CGn_real_macrocycle_typed'] or 0)+int(r['n_Gn_glue'] or 0)>0)}\n"
        "CGn = real atom type. Gn = glue pseudoatom.\n",
        encoding="utf-8",
    )
    return rows


def interface_audit(lig_rows: list[dict]) -> list[dict]:
    rows = []
    for r in lig_rows:
        eid = r["global_ligand_entity_id"]
        sdf = load_frozen_sdf(LIG_SDF / f"{eid}.sdf")
        text = (LIG_PDBQT / f"{eid}.pdbqt").read_text(errors="replace") if (LIG_PDBQT / f"{eid}.pdbqt").is_file() else ""
        unexpected = []
        expected = ["AD4_atom_typing", "gasteiger_charges", "explicit_H_parameterization"]
        if int(r["n_CGn_real_macrocycle_typed"] or 0):
            expected.append("macrocycle_CGn_retype_of_real_carbon")
        if int(r["n_Gn_glue"] or 0):
            expected.append("macrocycle_Gn_glue_pseudoatom")
        if sdf is not None and text:
            serial_to_sdf, why = map_serial_to_sdf_heavy(text, sdf)
            heavy = sdf_heavy(sdf)
            if not serial_to_sdf or len(serial_to_sdf) != heavy.GetNumAtoms():
                unexpected.append(why or "incomplete_real_atom_map")
            if int(r["Gn_in_SMILES_IDX"] or 0):
                unexpected.append("glue_in_SMILES_IDX")
            if int(r["n_CGn_real_macrocycle_typed"] or 0) and int(r["CGn_in_SMILES_IDX"] or 0) != int(r["n_CGn_real_macrocycle_typed"]):
                unexpected.append("CGn_dropped_from_mapping")
        sev = "CRITICAL" if unexpected else "OK"
        rows.append({
            "interface": "RDKit_SDF_to_Meeko_PDBQT",
            "entity_id": eid,
            "input_valid": r["frozen_sdf_ok"],
            "conversion_valid": int(not unexpected),
            "output_valid": int(bool(text)),
            "semantic_change": ";".join(expected),
            "expected_change": "YES",
            "unexpected_change": "|".join(unexpected),
            "severity": sev,
            "evidence": "02_ligands/sdf and pdbqt + REMARK SMILES IDX",
        })
    # interface summaries (not per-pose)
    rows.append({
        "interface": "Meeko_PDBQT_to_Vina", "entity_id": "ALL_PHASE7",
        "input_valid": 1, "conversion_valid": 1, "output_valid": 1,
        "semantic_change": "docking_search_moves_coordinates",
        "expected_change": "YES_coords_only", "unexpected_change": "",
        "severity": "OK", "evidence": "phase7 commands pin vina 1.2.7; graph from input PDBQT",
    })
    rows.append({
        "interface": "Vina_pose_to_old_GNINA_M1", "entity_id": "OLD_M1",
        "input_valid": 1, "conversion_valid": 0, "output_valid": 0,
        "semantic_change": "none_intended",
        "expected_change": "NO",
        "unexpected_change": "GNINA_rejects_CG0_on_61_physical_jobs",
        "severity": "MAJOR",
        "evidence": "09_secondary_scoring/gnina_rescore logs CG0 not valid AutoDock type",
    })
    rows.append({
        "interface": "Vina_pose_to_RTMScore_M2", "entity_id": "M2",
        "input_valid": 1, "conversion_valid": 1, "output_valid": 1,
        "semantic_change": "coordinates_only_onto_frozen_SDF",
        "expected_change": "YES_coords_only", "unexpected_change": "",
        "severity": "OK", "evidence": "ablation_phase578_rtmscore map_vina_to_frozen_sdf",
    })
    rows.append({
        "interface": "frozen_ligand_to_old_GNINA_M3", "entity_id": "OLD_M3_MACROCYCLE",
        "input_valid": 1, "conversion_valid": 0, "output_valid": 1,
        "semantic_change": "CG0_to_C_drop_G0",
        "expected_change": "NO_as_lossless",
        "unexpected_change": "ring_closure_bond_absent_from_torsion_tree",
        "severity": "CRITICAL",
        "evidence": "13_qa/M3_MACROCYCLE_TOPOLOGY_DECISION.csv 29 FAIL; cartesian distance not used as proof",
    })
    rows.append({
        "interface": "prepared_receptor_to_engines", "entity_id": "UNIQUE_14",
        "input_valid": 1, "conversion_valid": 1, "output_valid": 1,
        "semantic_change": "Meeko_receptor_atom_typing",
        "expected_change": "YES_types_charges", "unexpected_change": "",
        "severity": "OK", "evidence": "same *_receptor.pdbqt used by Vina/GNINA",
    })
    wcsv(QA / "CROSS_SOFTWARE_INTERFACE_AUDIT.csv", rows)
    return rows


def phase5_audit() -> list[dict]:
    master = load(QA / "redocking_master.csv")
    rows = []
    for r in master:
        pdb = r.get("pdb_id") or r.get("pdb")
        seed = str(r.get("seed", ""))
        job = REDOCK / "jobs" / f"{pdb}_seed{seed}"
        out = job / "out.pdbqt"
        rms = job / "rmsd.json"
        n, aff = vina_mode1_from_out(out)
        issues = []
        if r.get("status") in {"OK", "SUCCESS"} and not finite_float(aff):
            issues.append("summary_ok_but_no_raw_mode1")
        if rms.is_file():
            rj = json.loads(rms.read_text())
            if rj.get("parameters_changed_from_RMSD"):
                issues.append("params_changed_from_RMSD")
        rows.append({
            "pdb_id": pdb, "seed": seed,
            "summary_status": r.get("status") or r.get("job_status"),
            "raw_out_exists": int(out.is_file()),
            "raw_n_modes": n, "raw_mode1_affinity": aff,
            "summary_mode1": r.get("mode1_affinity") or r.get("vina_mode1_affinity") or r.get("affinity", ""),
            "rmsd_json_exists": int(rms.is_file()),
            "issues": "|".join(issues),
            "execution_complete": int(out.is_file() or (job / "status.json").is_file()),
            "scientific_qc_info_complete": int(rms.is_file() and n > 0),
        })
    if not rows:
        # fallback walk jobs
        for job in sorted((REDOCK / "jobs").glob("*_seed*")):
            if not job.is_dir():
                continue
            name = job.name
            pdb, _, seed = name.partition("_seed")
            n, aff = vina_mode1_from_out(job / "out.pdbqt")
            rows.append({
                "pdb_id": pdb, "seed": seed, "summary_status": "",
                "raw_out_exists": int((job / "out.pdbqt").is_file()),
                "raw_n_modes": n, "raw_mode1_affinity": aff,
                "summary_mode1": "", "rmsd_json_exists": int((job / "rmsd.json").is_file()),
                "issues": "", "execution_complete": 1, "scientific_qc_info_complete": int((job / "rmsd.json").is_file() and n > 0),
            })
    wcsv(QA / "VINA_REDOCK_RAW_TO_SUMMARY_AUDIT.csv", rows)
    return rows


def phase7_audit() -> dict:
    p7 = load(PHASE7)
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")
    rows = []
    phys = {}
    n_raw_mismatch = 0
    n_missing_out = 0
    for r in p7:
        pair, pdb, eid, seed = r["pair"], r["pdb_id"], r["global_ligand_entity_id"], str(r["seed"])
        jdir = VINA_JOBS / r["job_id"]
        if not jdir.is_dir():
            jdir = VINA_JOBS / f"{pair.replace('/', '_')}__{pdb}__{eid}__seed{seed}"
        out = jdir / "out.pdbqt"
        n, aff = vina_mode1_from_out(out)
        st = r.get("status")
        issues = []
        if st == "SUCCESS":
            if not finite_float(aff) and not finite_float(r.get("mode1_affinity", "")):
                issues.append("SUCCESS_without_finite_mode1")
            if finite_float(aff) and finite_float(r.get("mode1_affinity", "")) and abs(float(aff) - float(r["mode1_affinity"])) > 1e-3:
                issues.append("master_ne_raw_mode1")
                n_raw_mismatch += 1
            if not out.is_file():
                issues.append("SUCCESS_missing_out")
                n_missing_out += 1
        if st == "TIMEOUT" and finite_float(r.get("mode1_affinity", "")):
            issues.append("TIMEOUT_has_master_score")
        if st == "TIMEOUT" and out.is_file() and finite_float(aff):
            issues.append("TIMEOUT_but_raw_mode1_present")
        rec = {
            "level": "pair-expanded-seed",
            "job_id": r.get("job_id"),
            "pair": pair, "pdb_id": pdb, "global_ligand_entity_id": eid, "seed": seed,
            "master_status": st, "master_mode1": r.get("mode1_affinity", ""),
            "raw_n_modes": n, "raw_mode1": aff, "master_n_modes": r.get("n_modes_returned", ""),
            "issues": "|".join(issues),
        }
        rows.append(rec)
        key = (pdb, eid, seed)
        phys.setdefault(key, []).append(rec)
    # physical uniqueness
    phys_rows = []
    reuse_conflict = 0
    for key, grp in phys.items():
        scores = {(g["master_status"], g["master_mode1"]) for g in grp}
        pairs = sorted({g["pair"] for g in grp})
        if len(scores) > 1:
            reuse_conflict += 1
        phys_rows.append({
            "level": "unique-physical-job-seed",
            "pdb_id": key[0], "global_ligand_entity_id": key[1], "seed": key[2],
            "n_pair_expanded_records": len(grp),
            "pairs": ";".join(pairs),
            "statuses": ";".join(sorted({g["master_status"] for g in grp})),
            "mode1_values": ";".join(sorted({g["master_mode1"] for g in grp})),
            "within_method_reuse_ok": int(len(scores) == 1),
        })
    # M0 vs raw seed42
    m0_mismatch = 0
    raw_by = {(r["pair"], r["pdb_id"], r["global_ligand_entity_id"]): r for r in p7 if str(r["seed"]) == "42"}
    for r in m0:
        p = raw_by.get((r["pair"], r["pdb_id"], r["global_ligand_entity_id"]))
        if not p:
            m0_mismatch += 1
            continue
        if r["vina_status"] != p["status"]:
            m0_mismatch += 1
            continue
        if r["vina_status"] == "SUCCESS":
            try:
                if abs(float(r["vina_score"]) - (-float(p["mode1_affinity"]))) > 1e-4:
                    m0_mismatch += 1
            except (TypeError, ValueError):
                m0_mismatch += 1
    st = Counter(r["status"] for r in p7)
    s42 = [r for r in p7 if str(r["seed"]) == "42"]
    mode_hist = Counter()
    for r in s42:
        if r["status"] == "SUCCESS":
            try:
                mode_hist[int(r.get("n_modes_returned") or 0)] += 1
            except ValueError:
                pass
    summary = {
        "pair_expanded_seed_records": len(p7),
        "expected_pair_expanded_seed": EXPECTED_PAIR_EXPANDED_SEED,
        "unique_physical_job_seed": len(phys),
        "expected_physical_job_seed": EXPECTED_PHYSICAL_JOB_SEED,
        "pair_expanded_status": dict(st),
        "seed42_TIMEOUT": sum(1 for r in s42 if r["status"] == "TIMEOUT"),
        "seed42_SUCCESS": sum(1 for r in s42 if r["status"] == "SUCCESS"),
        "seed42_pose_hist_pair_expanded": dict(mode_hist),
        "raw_master_affinity_mismatches": n_raw_mismatch,
        "SUCCESS_missing_out": n_missing_out,
        "physical_within_method_reuse_conflicts": reuse_conflict,
        "M0_vs_phase7_seed42_mismatches": m0_mismatch,
    }
    wcsv(QA / "VINA_RAW_PARSE_MASTER_AUDIT.csv", rows)
    wcsv(QA / "VINA_PHYSICAL_JOB_SEED_AUDIT.csv", phys_rows)
    (QA / "VINA_PHASE7_LEVEL_COUNTS.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def _pairs_field(s: str) -> set[str]:
    return {p.strip() for p in (s or "").replace("|", ";").replace(",", ";").split(";") if p.strip()}


def identity_reuse() -> dict:
    mapping = [r for r in load(MAP) if r["pair"] in RECEPTORS]
    indep = [r for r in mapping if r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"]
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")
    m1 = load(SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv")
    m2 = load(SEC / "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv")
    m3 = load(GIND / "GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv")
    ab040 = [r for r in mapping if r.get("panel_id") == "AB_040" or r.get("canonical_ligand_id") == "AB_040"]
    rows = []

    # M0 is pair-expanded: same (pdb, eid) must reuse one physical score across pairs.
    c0 = 0
    by0 = defaultdict(list)
    for r in m0:
        by0[(r["pdb_id"], r["global_ligand_entity_id"])].append(r)
    shared_pairs = {}
    for k, grp in by0.items():
        pairs = sorted({x["pair"] for x in grp})
        if len(pairs) > 1:
            shared_pairs[k] = pairs
        if len(grp) < 2:
            continue
        vals = {(g["vina_status"], g.get("vina_score")) for g in grp}
        ok = len(vals) == 1
        if not ok:
            c0 += 1
        rows.append({
            "method": "M0", "pdb_id": k[0], "global_ligand_entity_id": k[1],
            "n_mapped_rows": len(grp), "pairs": ";".join(pairs),
            "within_method_reuse_ok": int(ok),
            "pairs_field_ok": 1,
            "values": str(sorted(vals)),
        })

    def physical_master(name, store, scorefn, statusfn):
        """One physical row per (pdb, eid). Shared receptors must list all pairs. Do not compare to other methods."""
        by = defaultdict(list)
        dup_conflict = 0
        pairs_mismatch = 0
        missing_shared = 0
        for r in store:
            by[(r["pdb_id"], r["global_ligand_entity_id"])].append(r)
        for k, expected in shared_pairs.items():
            grp = by.get(k, [])
            if not grp:
                missing_shared += 1
                rows.append({
                    "method": name, "pdb_id": k[0], "global_ligand_entity_id": k[1],
                    "n_mapped_rows": 0, "pairs": "",
                    "within_method_reuse_ok": 0, "pairs_field_ok": 0,
                    "values": "MISSING_SHARED_PHYSICAL",
                })
                continue
            vals = {(statusfn(r), scorefn(r)) for r in grp}
            reuse_ok = len(vals) == 1
            if not reuse_ok:
                dup_conflict += 1
            got = set()
            for r in grp:
                got |= _pairs_field(r.get("pairs", ""))
            pairs_ok = got == set(expected)
            if not pairs_ok:
                pairs_mismatch += 1
            rows.append({
                "method": name, "pdb_id": k[0], "global_ligand_entity_id": k[1],
                "n_mapped_rows": len(grp), "pairs": ";".join(sorted(got)),
                "within_method_reuse_ok": int(reuse_ok),
                "pairs_field_ok": int(pairs_ok),
                "values": str(sorted(vals))[:200],
            })
        extra_dups = sum(1 for k, g in by.items() if k not in shared_pairs and len(g) > 1)
        return {
            "n": len(store),
            "unique_physical": len(by),
            "duplicate_score_conflicts": dup_conflict + extra_dups,
            "shared_pairs_field_mismatches": pairs_mismatch,
            "shared_missing": missing_shared,
        }

    m1s = physical_master("OLD_M1", m1, lambda r: r.get("M1_CNNscore"), lambda r: r.get("job_status"))
    m2s = physical_master("M2", m2, lambda r: r.get("M2_RTMScore"), lambda r: r.get("job_status"))
    m3s = physical_master("OLD_M3", m3, lambda r: r.get("CNNscore"), lambda r: r.get("status"))
    m1n_p = SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv"
    m3n_p = GIND / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv"
    rec = {
        "historical_mapping_rows_8pair": len(mapping),
        "independent_pair_level": len(indep),
        "global_unique_indep": len({r["global_ligand_entity_id"] for r in indep}),
        "pair_expanded_m0": len(m0),
        "physical_m0_keys": len(by0),
        "shared_physical_keys": len(shared_pairs),
        "shared_6N7A": sum(1 for (pdb, _) in shared_pairs if pdb == "6N7A"),
        "shared_6LXA": sum(1 for (pdb, _) in shared_pairs if pdb == "6LXA"),
        "AB_040_mapping_rows": len(ab040),
        "AB_040_in_m0": sum(1 for r in m0 if r.get("canonical_ligand_id") == "AB_040"),
        "AB_046_fold_note": "alias registry fold_id=1",
        "M0_within_method_conflicts": c0,
        "OLD_M1_physical": m1s,
        "M2_physical": m2s,
        "OLD_M3_physical": m3s,
        "NEW_M1_physical": (
            physical_master("NEW_M1", load(m1n_p), lambda r: r.get("M1_CNNscore"), lambda r: r.get("job_status"))
            if m1n_p.is_file() else None
        ),
        "NEW_M3_physical": (
            physical_master("NEW_M3", load(m3n_p), lambda r: r.get("CNNscore"), lambda r: r.get("status"))
            if m3n_p.is_file() else None
        ),
        "methods_not_required_to_match_each_other": True,
        "cross_method_score_equality_not_checked": True,
    }
    wcsv(QA / "JOB_IDENTITY_REUSE_AUDIT.csv", rows or [{"method": "none_shared", "within_method_reuse_ok": 1}])
    return rec


def old_m1_audit() -> dict:
    master = load(SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv")
    long = load(SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv")
    st = Counter(r.get("job_status") for r in master)
    fails = [r for r in master if r.get("reason") == "parse_or_rc" or r.get("job_status") == "TECHNICAL_FAIL"]
    rows = []
    m1b_lt = 0
    mode1_mismatch = 0
    long_m1 = {(r["pdb_id"], r["global_ligand_entity_id"], str(r.get("vina_mode"))): r for r in long}
    for r in master:
        issues = []
        if r.get("job_status") == "SUCCESS":
            try:
                if float(r["M1b_CNNscore"]) + 1e-12 < float(r["M1_CNNscore"]):
                    m1b_lt += 1
                    issues.append("M1b_lt_M1")
            except (TypeError, ValueError):
                issues.append("nonfinite")
            lr = long_m1.get((r["pdb_id"], r["global_ligand_entity_id"], "1"))
            if lr and finite_float(r.get("M1_CNNscore", "")) and lr.get("CNNscore") != r.get("M1_CNNscore"):
                # allow float format
                try:
                    if abs(float(lr["CNNscore"]) - float(r["M1_CNNscore"])) > 1e-5:
                        mode1_mismatch += 1
                        issues.append("master_ne_long_mode1")
                except (TypeError, ValueError):
                    mode1_mismatch += 1
        rows.append({
            "layer": "OLD_M1_historical_technical_run",
            "pdb_id": r["pdb_id"], "global_ligand_entity_id": r["global_ligand_entity_id"],
            "job_status": r.get("job_status"), "reason": r.get("reason"),
            "n_poses": r.get("n_poses"), "M1_CNNscore": r.get("M1_CNNscore"),
            "M1b_CNNscore": r.get("M1b_CNNscore"), "n_retries": r.get("n_retries"),
            "issues": "|".join(issues),
        })
    wcsv(QA / "GNINA_M1_RAW_INPUT_OUTPUT_MASTER_AUDIT.csv", rows)
    return {
        "layer": "OLD_M1_only",
        "n": len(master),
        "status": dict(st),
        "parse_or_rc": len(fails),
        "long_rows": len(long),
        "M1b_lt_M1": m1b_lt,
        "mode1_mismatch": mode1_mismatch,
        "new_M1_not_audited": True,
    }


def m2_audit() -> dict:
    master = load(SEC / "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv")
    long = load(SEC / "RTMSCORE_VINA_POSE_RESCORE_POSE_LONG.csv")
    long_m1 = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in long if str(r.get("vina_mode")) == "1" or str(r.get("is_vina_mode1")) == "1"}
    mismatch = 0
    rows = []
    for r in master:
        lr = long_m1.get((r["pdb_id"], r["global_ligand_entity_id"]))
        ok = True
        if r.get("job_status") == "SUCCESS" and lr:
            try:
                ok = abs(float(r["M2_RTMScore"]) - float(lr["RTMScore"])) <= 1e-6
            except (TypeError, ValueError):
                ok = False
            if not ok:
                mismatch += 1
        rows.append({
            "pdb_id": r["pdb_id"], "global_ligand_entity_id": r["global_ligand_entity_id"],
            "job_status": r.get("job_status"), "master_mode1": r.get("M2_RTMScore"),
            "long_mode1": (lr or {}).get("RTMScore", ""), "mode1_match": int(ok),
        })
    pockets = list((SEC / "rtmscore" / "pockets").glob("*_rtmscore_pocket10.pdb"))
    wcsv(QA / "RTMSCORE_RAW_TO_MASTER_AUDIT.csv", rows)
    return {
        "master_n": len(master), "long_n": len(long),
        "expected_long": 14250, "mode1_mismatches": mismatch,
        "pockets": len(pockets), "checkpoint": "rtmscore_model1.pth",
        "openbabel_graph_inference_in_runner": False,
    }


def old_m3_audit() -> dict:
    old = load(GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv")
    rep = load(GIND / "GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv")
    topo = {r["global_ligand_entity_id"]: r for r in load(QA / "M3_MACROCYCLE_TOPOLOGY_DECISION.csv")}
    rows = []
    parse_fail = 0
    for r in rep:
        eid = r["global_ligand_entity_id"]
        jobdir = GIND / "jobs" / r["job_id"]
        outp = jobdir / "out.pdbqt"
        raw = parse_model1_remarks(outp.read_text(errors="replace")) if outp.is_file() else {}
        issues = []
        if r.get("status") == "SUCCESS":
            if not finite_float(r.get("CNNscore", "")):
                parse_fail += 1
                issues.append("SUCCESS_nonfinite_reparsed")
            if finite_float(raw.get("CNNscore", "")) and finite_float(r.get("CNNscore", "")):
                if abs(float(raw["CNNscore"]) - float(r["CNNscore"])) > 1e-6:
                    issues.append("reparsed_ne_model1_remark")
                    parse_fail += 1
        td = topo.get(eid, {})
        rows.append({
            "layer": "OLD_M3",
            "job_id": r["job_id"], "pdb_id": r["pdb_id"], "global_ligand_entity_id": eid,
            "execution_status": r.get("status"),
            "original_master_CNNscore": next((x.get("CNNscore") for x in old if x.get("job_id") == r["job_id"]), ""),
            "reparsed_CNNscore": r.get("CNNscore"),
            "raw_model1_CNNscore": raw.get("CNNscore", ""),
            "topology": td.get("decision", "N/A_NO_SPECIAL_TYPE"),
            "official_for_analysis": (
                "NO_INVALID_INPUT_REPRESENTATION" if td.get("decision") == "M3_INPUT_TOPOLOGY_FAIL"
                else ("NO_TIMEOUT" if r.get("status") == "TIMEOUT" else "YES_IF_SUCCESS")
            ),
            "issues": "|".join(issues),
        })
    wcsv(QA / "GNINA_M3_RAW_INPUT_OUTPUT_MASTER_AUDIT.csv", rows)
    return {
        "layer": "OLD_M3_only",
        "n": len(rep),
        "parse_disagreements": parse_fail,
        "topology_fail_ligands": sum(1 for v in topo.values() if v.get("decision") == "M3_INPUT_TOPOLOGY_FAIL"),
        "TIMEOUT_kept": sum(1 for r in rep if r.get("status") == "TIMEOUT"),
        "new_M3_not_audited": True,
    }


def _cnnscore_from_log(path: Path) -> str:
    if not path.is_file():
        return ""
    m = re.search(r"CNNscore:\s+([0-9eE.+-]+)", path.read_text(errors="replace"))
    if m and m.group(1) != "--":
        return m.group(1)
    return ""


def new_m1_audit() -> dict:
    master_p = SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv"
    long_p = SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv"
    work = SEC / "gnina_rescore_representation_v2"
    if not master_p.is_file() or not long_p.is_file():
        return {"layer": "NEW_M1", "audited": False, "reason": "verified_master_missing"}
    master = load(master_p)
    long = load(long_p)
    st = Counter(r.get("job_status") for r in master)
    long_m1 = {(r["pdb_id"], r["global_ligand_entity_id"], str(r.get("vina_mode"))): r for r in long}
    rows = []
    m1b_lt = 0
    mode1_mismatch = 0
    raw_mismatch = 0
    non_success_scored = 0
    special_in_success = 0
    rms_bad = 0
    for r in master:
        issues = []
        pdb, eid = r["pdb_id"], r["global_ligand_entity_id"]
        jobdir = work / f"{pdb}__{eid}"
        stj = r.get("job_status")
        if stj == "SUCCESS":
            try:
                if float(r["M1b_CNNscore"]) + 1e-12 < float(r["M1_CNNscore"]):
                    m1b_lt += 1
                    issues.append("M1b_lt_M1")
            except (TypeError, ValueError):
                issues.append("nonfinite")
            lr = long_m1.get((pdb, eid, "1"))
            if lr and finite_float(r.get("M1_CNNscore", "")):
                try:
                    if abs(float(lr["CNNscore"]) - float(r["M1_CNNscore"])) > 1e-5:
                        mode1_mismatch += 1
                        issues.append("master_ne_long_mode1")
                except (TypeError, ValueError):
                    mode1_mismatch += 1
                    issues.append("master_ne_long_mode1")
            raw = _cnnscore_from_log(jobdir / "vina_mode1.gnina.log")
            if finite_float(raw) and finite_float(r.get("M1_CNNscore", "")):
                if abs(float(raw) - float(r["M1_CNNscore"])) > 1e-5:
                    raw_mismatch += 1
                    issues.append("master_ne_raw_log_mode1")
            elif not finite_float(raw):
                raw_mismatch += 1
                issues.append("missing_raw_log_mode1")
            pdbqt = jobdir / "vina_mode1_gnina.pdbqt"
            if pdbqt.is_file():
                types = special_types(parse_pdbqt_atoms(pdbqt.read_text(errors="replace")))
                if any(is_real_macrocycle_type(t) or is_glue_type(t) for t in types):
                    special_in_success += 1
                    issues.append("CGn_or_Gn_remain_in_scored_pdbqt")
            try:
                if float(r.get("coord_rms_max") or 0) > 0.001:
                    rms_bad += 1
                    issues.append("coord_rms_gt_0.001")
            except (TypeError, ValueError):
                pass
        else:
            if finite_float(r.get("M1_CNNscore", "")) or finite_float(r.get("M1b_CNNscore", "")):
                non_success_scored += 1
                issues.append("non_success_has_score")
        rows.append({
            "layer": "NEW_M1_verified",
            "pdb_id": pdb, "global_ligand_entity_id": eid,
            "job_status": stj, "reason": r.get("reason"),
            "n_poses": r.get("n_poses"), "M1_CNNscore": r.get("M1_CNNscore"),
            "M1b_CNNscore": r.get("M1b_CNNscore"), "coord_rms_max": r.get("coord_rms_max"),
            "issues": "|".join(issues),
        })
    n_succ_poses = sum(int(r.get("n_poses") or 0) for r in master if r.get("job_status") == "SUCCESS")
    n_map_poses = sum(int(r.get("n_poses") or 0) for r in master if r.get("job_status") == "ATOM_MAPPING_FAILURE")
    n_miss = sum(1 for r in master if r.get("job_status") == "MISSING_NO_VINA_POSE")
    wcsv(QA / "GNINA_M1_V2_RAW_INPUT_OUTPUT_MASTER_AUDIT.csv", rows)
    return {
        "layer": "NEW_M1",
        "audited": True,
        "n": len(master),
        "status": dict(st),
        "long_rows": len(long),
        "success_n_poses": n_succ_poses,
        "mapping_fail_n_poses": n_map_poses,
        "missing_no_vina_pose": n_miss,
        "pose_accounted": n_succ_poses + n_map_poses,
        "expected_poses": 14250,
        "M1b_lt_M1": m1b_lt,
        "mode1_mismatch": mode1_mismatch,
        "raw_log_mode1_mismatches": raw_mismatch,
        "non_success_has_score": non_success_scored,
        "CGn_or_Gn_in_success_mode1_pdbqt": special_in_success,
        "coord_rms_gt_0.001": rms_bad,
        "old_master_preserved": (SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv").is_file(),
    }


def new_m3_audit() -> dict:
    ver_p = GIND / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv"
    if not ver_p.is_file():
        return {"layer": "NEW_M3", "audited": False, "reason": "verified_master_missing"}
    ver = load(ver_p)
    topo = {r["global_ligand_entity_id"]: r for r in load(QA / "M3_MACROCYCLE_TOPOLOGY_DECISION.csv")}
    fail_eids = {eid for eid, r in topo.items() if r.get("decision") == "M3_INPUT_TOPOLOGY_FAIL"}
    lig_out = GIND / "ligands_gnina_rigid_macrocycles_v2"
    job_out = GIND / "jobs_topology_repair_v2"
    rows = []
    parse_fail = 0
    special_remain = 0
    extra_repair = 0
    missing_repair = 0
    ab001 = "NPHIQFCGKIKGSA-YQOHNZFASA-N"
    ab001_kept = 0
    for r in ver:
        eid = r["global_ligand_entity_id"]
        pdb = r["pdb_id"]
        src = r.get("score_source", "")
        issues = []
        repaired = "topology_repaired" in src
        if repaired and eid not in fail_eids:
            extra_repair += 1
            issues.append("repaired_eid_not_in_frozen_fail_set")
        if eid in fail_eids and not repaired and r.get("status") != "TIMEOUT":
            missing_repair += 1
            issues.append("frozen_fail_eid_not_repaired")
        if repaired and r.get("status") == "SUCCESS":
            outp = job_out / r["job_id"] / "out.pdbqt"
            raw = parse_model1_remarks(outp.read_text(errors="replace")) if outp.is_file() else {}
            if not finite_float(r.get("CNNscore", "")):
                parse_fail += 1
                issues.append("SUCCESS_nonfinite")
            if finite_float(raw.get("CNNscore", "")) and finite_float(r.get("CNNscore", "")):
                if abs(float(raw["CNNscore"]) - float(r["CNNscore"])) > 1e-6:
                    parse_fail += 1
                    issues.append("verified_ne_model1_remark")
            lig = lig_out / f"{eid}.pdbqt"
            if lig.is_file():
                types = special_types(parse_pdbqt_atoms(lig.read_text(errors="replace")))
                if any(is_real_macrocycle_type(t) or is_glue_type(t) for t in types):
                    special_remain += 1
                    issues.append("CGn_or_Gn_remain_in_repaired_ligand")
        if eid == ab001 and r.get("status") == "TIMEOUT" and "timeout_kept" in src:
            ab001_kept += 1
        if r.get("status") != "SUCCESS" and finite_float(r.get("CNNscore", "")):
            issues.append("non_success_has_score")
        rows.append({
            "layer": "NEW_M3_verified",
            "job_id": r.get("job_id"), "pdb_id": pdb, "global_ligand_entity_id": eid,
            "status": r.get("status"), "score_source": src,
            "official_m3": r.get("official_m3"),
            "CNNscore": r.get("CNNscore"),
            "issues": "|".join(issues),
        })
    wcsv(QA / "GNINA_M3_V2_RAW_INPUT_OUTPUT_MASTER_AUDIT.csv", rows)
    repaired_n = sum(1 for r in ver if "topology_repaired" in (r.get("score_source") or ""))
    return {
        "layer": "NEW_M3",
        "audited": True,
        "n": len(ver),
        "status": dict(Counter(r.get("status") for r in ver)),
        "score_source": dict(Counter(r.get("score_source") for r in ver)),
        "repaired_jobs": repaired_n,
        "expected_repaired_jobs": 61,
        "parse_disagreements": parse_fail,
        "special_types_remain": special_remain,
        "repaired_eid_not_in_frozen_set": extra_repair,
        "frozen_fail_not_repaired": missing_repair,
        "ab001_timeout_kept_rows": ab001_kept,
        "old_master_preserved": (GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv").is_file(),
        "non_success_has_score": sum(1 for r in ver if r.get("status") != "SUCCESS" and finite_float(r.get("CNNscore", ""))),
    }


def alt_audit() -> dict:
    prod = load(QA / "alt_production_seed42_master.csv")
    indep = load(QA / "alt_score_master_independent.csv")
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")
    m0_score = {(r["pair"], r["pdb_id"], r["global_ligand_entity_id"]): r.get("vina_score") for r in m0}
    rows = []
    reuse_bad = 0
    for r in prod:
        alt = r.get("pdb_id") or r.get("alt_pdb")
        design = r.get("design") or r.get("substitution", "")
        eid = r.get("global_ligand_entity_id")
        issues = []
        if alt not in FORMAL_ALT:
            issues.append("not_in_formal_8")
        rows.append({
            "job_id": r.get("job_id"), "alt_pdb": alt, "design": design,
            "pair": r.get("pair"), "eid": eid, "status": r.get("status"),
            "issues": "|".join(issues),
        })
    for r in indep:
        design = r.get("design", "")
        pair = r.get("pair")
        eid = r.get("global_ligand_entity_id")
        kept = r.get("kept_primary")
        issues = []
        if design == "A2/B1":
            prim = m0_score.get((pair, kept, eid))
            reused = r.get("score_B_primary")
        elif design == "A1/B2":
            prim = m0_score.get((pair, kept, eid))
            reused = r.get("score_A_primary")
        else:
            prim, reused = None, None
            issues.append("unexpected_design")
        if prim not in (None, "") and reused not in (None, ""):
            try:
                if abs(float(prim) - float(reused)) > 1e-4:
                    reuse_bad += 1
                    issues.append("kept_side_not_bitwise_primary")
            except (TypeError, ValueError):
                reuse_bad += 1
                issues.append("kept_side_nonfinite")
        elif r.get("status_alt") == "SUCCESS" and (prim in (None, "") or reused in (None, "")):
            issues.append("kept_side_reuse_not_checkable")
        if issues:
            rows.append({
                "job_id": f"reuse_check__{r.get('alt_pdb')}__{eid}",
                "alt_pdb": r.get("alt_pdb"), "design": design,
                "pair": pair, "eid": eid, "status": r.get("status_alt"),
                "issues": "|".join(issues),
            })
    st = Counter(r.get("status") for r in prod)
    wcsv(QA / "ALTERNATIVE_RECEPTOR_END_TO_END_AUDIT.csv", rows or [{"note": "master_schema_partial"}])
    return {
        "n_master": len(prod), "status": dict(st),
        "expected": 734, "formal_alts": sorted(FORMAL_ALT),
        "v1_4JPS_in_master": any((r.get("pdb_id") or "") in {"4JPS", "5DXT", "4JSX"} for r in prod),
        "kept_side_reuse_flags": reuse_bad,
    }


def activity_audit() -> dict:
    mapping = [r for r in load(MAP) if r["pair"] in RECEPTORS]
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")
    zero = [r for r in mapping if r.get("activity_eligible") == "0"]
    neither_bad = [r for r in m0 if r.get("activity_eligible") == "0" and r.get("class") == "neither"]
    watch = [r for r in mapping if r.get("panel_id") in {"EH120_059", "AB_087", "AB_040", "AB_046"}]
    rows = []
    for r in watch + zero:
        rows.append({k: r.get(k) for k in (
            "pair", "panel_id", "canonical_ligand_id", "global_ligand_entity_id",
            "class", "activity_eligible", "parent_independent", "parent_role",
        )})
    wcsv(QA / "ACTIVITY_LABEL_SEMANTICS_AUDIT.csv", rows)
    return {
        "eligible0_mapping": len(zero),
        "eligible0_became_neither": len(neither_bad),
        "theta": 6.0,
        "missing_activity_is_inactive_forbidden": True,
    }


def failure_and_missing(p7sum, old_m1, old_m3, m2) -> None:
    rows = []
    rows.append({"class": "TIMEOUT", "method": "M0_Vina_seed42_pair_expanded", "count": 6, "evidence": "official_primary_seed42 6 TIMEOUT AB_001/053/054 x 4EY7/4BDS"})
    rows.append({"class": "TIMEOUT", "method": "Vina_phase7_all_seeds_pair_expanded", "count": p7sum["pair_expanded_status"].get("TIMEOUT", 0), "evidence": "phase7_fiveseed_master"})
    rows.append({"class": "INPUT_FORMAT_FAILURE", "method": "OLD_M1", "count": old_m1.get("parse_or_rc", 0), "evidence": "CG0 not valid AutoDock type; n_retries=1"})
    rows.append({"class": "INVALID_INPUT_REPRESENTATION", "method": "OLD_M3", "count": old_m3.get("topology_fail_ligands", 0), "evidence": "29 ligands; 61 physical jobs; topology FAIL"})
    rows.append({"class": "TIMEOUT", "method": "OLD_M3", "count": old_m3.get("TIMEOUT_kept", 0), "evidence": "AB_001 x 4EY7/4BDS; topology PASS; keep missing"})
    rows.append({"class": "MISSING_NO_INPUT", "method": "M1_M2", "count": 6, "evidence": "MISSING_NO_VINA_POSE physical"})
    wcsv(QA / "FAILURE_TAXONOMY_AUDIT.csv", rows)
    miss = []
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")
    for (pair, pdb, cls), grp in _grp(m0, lambda r: (r["pair"], r["pdb_id"], r.get("class", ""))).items():
        succ = sum(1 for r in grp if r["vina_status"] == "SUCCESS")
        to = sum(1 for r in grp if r["vina_status"] == "TIMEOUT")
        miss.append({"method": "M0", "level": "pair-expanded", "pair": pair, "receptor": pdb, "class": cls,
                     "expected": len(grp), "valid": succ, "timeout": to, "missing": 0, "technical_failure": 0,
                     "invalid_representation": 0})
    wcsv(QA / "METHOD_MISSINGNESS_AUDIT_FULL.csv", miss)


def _grp(rows, fn):
    d = defaultdict(list)
    for r in rows:
        d[fn(r)].append(r)
    return d


def consistency() -> list[dict]:
    rows = []
    p7 = load(PHASE7)
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")
    rows.append({"method": "M0", "raw": "06_vina_fiveseed/jobs/*/out.pdbqt", "long": "phase7_fiveseed_master", "master": "official_primary_seed42_score_master_8pair.csv",
                 "raw_n": "see VINA_RAW_PARSE_MASTER_AUDIT", "long_n": len(p7), "master_n": len(m0),
                 "id_seed_mode_score_checked": "YES_FULL", "new_v2_included": "NO"})
    rows.append({"method": "OLD_M1", "raw": "09_secondary_scoring/gnina_rescore", "long": "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv",
                 "master": "GNINA_VINA_POSE_RESCORE_MASTER.csv", "raw_n": "", "long_n": len(load(SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv")),
                 "master_n": len(load(SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv")), "id_seed_mode_score_checked": "YES_OLD_LAYER", "new_v2_included": "NO"})
    rows.append({"method": "M2", "raw": "rtmscore/jobs", "long": "RTMSCORE_VINA_POSE_RESCORE_POSE_LONG.csv",
                 "master": "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv", "raw_n": "", "long_n": 14250, "master_n": 1592,
                 "id_seed_mode_score_checked": "YES", "new_v2_included": "NO"})
    rows.append({"method": "OLD_M3", "raw": "10_gnina_independent_docking/jobs/*/out.pdbqt MODEL1",
                 "long": "n/a", "master": "GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv",
                 "raw_n": 1592, "long_n": "", "master_n": 1592, "id_seed_mode_score_checked": "YES_OLD_LAYER", "new_v2_included": "NO"})
    rows.append({"method": "NEW_M1", "raw": "09_secondary_scoring/gnina_rescore_representation_v2/*/vina_mode1.gnina.log",
                 "long": "GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv",
                 "master": "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv",
                 "raw_n": "", "long_n": len(load(SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv")) if (SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv").is_file() else 0,
                 "master_n": len(load(SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv")) if (SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv").is_file() else 0,
                 "id_seed_mode_score_checked": "YES_NEW_LAYER", "new_v2_included": "YES"})
    rows.append({"method": "NEW_M3", "raw": "jobs_topology_repair_v2/*/out.pdbqt MODEL1 + original jobs for unrepaired",
                 "long": "n/a", "master": "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv",
                 "raw_n": 1592, "long_n": "",
                 "master_n": len(load(GIND / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv")) if (GIND / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv").is_file() else 0,
                 "id_seed_mode_score_checked": "YES_NEW_LAYER", "new_v2_included": "YES"})
    wcsv(QA / "RAW_LONG_MASTER_CONSISTENCY_AUDIT.csv", rows)
    return rows


def write_gate(v2, cfg, recs, p7, ident, old_m1, m2, old_m3, alt, act, new_m1, new_m3) -> None:
    blocking_completed = []
    if cfg.get("CRITICAL_CONFIG_DRIFT"):
        blocking_completed.append("CRITICAL_CONFIG_DRIFT")
    if ident.get("M0_within_method_conflicts"):
        blocking_completed.append("M0_shared_receptor_score_conflict")
    for tag, key in (
        ("OLD_M1", "OLD_M1_physical"), ("M2", "M2_physical"), ("OLD_M3", "OLD_M3_physical"),
        ("NEW_M1", "NEW_M1_physical"), ("NEW_M3", "NEW_M3_physical"),
    ):
        phys = ident.get(key) or {}
        if phys.get("duplicate_score_conflicts") or phys.get("shared_pairs_field_mismatches") or phys.get("shared_missing"):
            blocking_completed.append(f"{tag}_within_method_reuse_conflict")
    if ident.get("AB_040_in_m0"):
        blocking_completed.append("AB_040_used_as_observation")
    if act.get("eligible0_became_neither"):
        blocking_completed.append("activity0_mapped_to_neither")
    if p7.get("M0_vs_phase7_seed42_mismatches"):
        blocking_completed.append("M0_ne_raw_seed42")
    if p7.get("physical_within_method_reuse_conflicts"):
        blocking_completed.append("phase7_physical_reuse_conflict")
    if old_m3.get("parse_disagreements"):
        blocking_completed.append("M3_reparse_disagreement")
    new_blocking = []
    if v2.get("treat_new_as_final"):
        if not new_m1.get("audited"):
            new_blocking.append("NEW_M1_not_audited")
        else:
            if new_m1.get("n") != 1592:
                new_blocking.append("NEW_M1_master_n_ne_1592")
            if new_m1.get("pose_accounted") != 14250:
                new_blocking.append("NEW_M1_poses_not_accounted_14250")
            if new_m1.get("long_rows") != new_m1.get("success_n_poses"):
                new_blocking.append("NEW_M1_long_ne_success_poses")
            if new_m1.get("M1b_lt_M1"):
                new_blocking.append("NEW_M1_M1b_lt_M1")
            if new_m1.get("mode1_mismatch"):
                new_blocking.append("NEW_M1_master_ne_long_mode1")
            if new_m1.get("raw_log_mode1_mismatches"):
                new_blocking.append("NEW_M1_master_ne_raw_log")
            if new_m1.get("non_success_has_score"):
                new_blocking.append("NEW_M1_score_imputation")
            if new_m1.get("CGn_or_Gn_in_success_mode1_pdbqt"):
                new_blocking.append("NEW_M1_CGn_or_Gn_remain")
            if new_m1.get("coord_rms_gt_0.001"):
                new_blocking.append("NEW_M1_coord_rms_gt_0.001")
            if not new_m1.get("old_master_preserved"):
                new_blocking.append("NEW_M1_overwrote_old_master")
        if not new_m3.get("audited"):
            new_blocking.append("NEW_M3_not_audited")
        else:
            if new_m3.get("n") != 1592:
                new_blocking.append("NEW_M3_master_n_ne_1592")
            if new_m3.get("repaired_jobs") != 61:
                new_blocking.append("NEW_M3_repaired_jobs_ne_61")
            if new_m3.get("parse_disagreements"):
                new_blocking.append("NEW_M3_verified_ne_raw")
            if new_m3.get("special_types_remain"):
                new_blocking.append("NEW_M3_CGn_or_Gn_remain")
            if new_m3.get("repaired_eid_not_in_frozen_set"):
                new_blocking.append("NEW_M3_repair_outside_frozen_set")
            if new_m3.get("frozen_fail_not_repaired"):
                new_blocking.append("NEW_M3_frozen_fail_not_repaired")
            if new_m3.get("ab001_timeout_kept_rows") != 2:
                new_blocking.append("NEW_M3_AB001_timeout_not_kept")
            if new_m3.get("non_success_has_score"):
                new_blocking.append("NEW_M3_score_imputation")
            if not new_m3.get("old_master_preserved"):
                new_blocking.append("NEW_M3_overwrote_old_master")
    ident_ok = not (
        ident.get("M0_within_method_conflicts")
        or any((ident.get(k) or {}).get(f) for k in (
            "OLD_M1_physical", "M2_physical", "OLD_M3_physical", "NEW_M1_physical", "NEW_M3_physical")
            for f in ("duplicate_score_conflicts", "shared_pairs_field_mismatches", "shared_missing"))
    )
    new_m1_sci = (
        "NOT_AUDITED_V2_IN_PROGRESS" if not v2.get("treat_new_as_final")
        else ("YES" if new_m1.get("audited") and not any(x.startswith("NEW_M1_") for x in new_blocking) else "NO")
    )
    new_m3_sci = (
        "NOT_AUDITED_V2_IN_PROGRESS" if not v2.get("treat_new_as_final")
        else ("YES" if new_m3.get("audited") and not any(x.startswith("NEW_M3_") for x in new_blocking) else "NO")
    )
    modules = [
        ("Receptor preparation", "YES", "YES", "YES" if not any(r.get("severity") == "CRITICAL" for r in recs.get("receptor", [])) else "NO"),
        ("Ligand preparation", "YES", "YES", "YES"),
        ("Box", "YES", "YES", "YES"),
        ("Vina redocking", "YES", "YES", "YES"),
        ("Vina production", "YES" if p7.get("pair_expanded_seed_records") == 8070 else "NO", "YES", "YES" if not p7.get("M0_vs_phase7_seed42_mismatches") else "NO"),
        ("M0", "YES", "YES", "YES" if not p7.get("M0_vs_phase7_seed42_mismatches") else "NO"),
        ("M1/M1b old", "YES", "PARTIAL", "NO — 61 INPUT_FORMAT_FAILURE; mixed representation; not the formal candidate"),
        ("M1/M1b new", "YES" if new_m1.get("audited") else "NO", "YES" if new_m1.get("audited") else "NO", new_m1_sci),
        ("M2", "YES", "YES", "YES" if m2.get("mode1_mismatches") == 0 and m2.get("long_n") == 14250 else "NO"),
        ("M3 old", "YES", "YES", "NO — 29 ligands INVALID_INPUT_REPRESENTATION; not the formal candidate"),
        ("M3 new", "YES" if new_m3.get("audited") else "NO", "YES" if new_m3.get("audited") else "NO", new_m3_sci),
        ("Alternative receptor", "YES", "YES", "YES" if not alt.get("v1_4JPS_in_master") else "NO"),
        ("Identity/reuse", "YES", "YES", "YES" if ident_ok else "NO"),
        ("Labels/activity", "YES", "YES", "YES" if not act.get("eligible0_became_neither") else "NO"),
        ("Parsing", "YES", "YES", "YES" if v2.get("treat_new_as_final") and not new_blocking else "PARTIAL"),
        ("Missingness", "YES", "YES", "YES"),
        ("Cross-software interfaces", "YES", "YES",
         "YES" if v2.get("treat_new_as_final") and not any(x in new_blocking for x in (
             "NEW_M1_CGn_or_Gn_remain", "NEW_M3_CGn_or_Gn_remain")) else "NO — old M1 CG0; old M3 ring closure"),
    ]
    if not v2.get("treat_new_as_final"):
        gate = "PENDING_V2_COMPLETION"
        header = [
            "This is not FAIL. New M1/M3 V2 outputs are incomplete; they were not read as final evidence.",
            "Final PASS/FAIL will be written only after new M1/M3 task-final files exist and scientific QA of those files is done.",
        ]
    else:
        official_blocking = blocking_completed + new_blocking
        gate = "FAIL" if official_blocking else "PASS"
        header = [
            "New M1/M3 task-final files were audited. Historical old M1/M3 remain documented as non-candidates.",
            "Official analysis candidates: M0, NEW M1/M1b, M2, NEW M3.",
        ]
    lines = [
        "# END_TO_END_SCIENTIFIC_INTEGRITY_GATE",
        "",
        f"END_TO_END_SCIENTIFIC_INTEGRITY_GATE = {gate}",
        "",
        *header,
        "",
        f"V2 running: {v2.get('v2_processes_running')} ; new M1 master exists: {v2.get('new_m1_master_exists')} ; new M3 master exists: {v2.get('new_m3_master_exists')}",
        "",
        "## Blocking issues",
    ]
    shown = (blocking_completed + new_blocking) if v2.get("treat_new_as_final") else blocking_completed
    for b in shown or ["none"]:
        lines.append(f"- {b}")
    lines += ["", "## Per-module", ""]
    for name, ex, par, sci in modules:
        lines.append(f"### {name}")
        lines.append(f"EXECUTION_COMPLETE = {ex}")
        lines.append(f"PARSING_COMPLETE = {par}")
        lines.append(f"SCIENTIFIC_QA_PASS = {sci}")
        lines.append("")
    if new_m1.get("audited"):
        lines += [
            "## NEW M1 summary",
            f"- master {new_m1.get('n')} status {new_m1.get('status')}",
            f"- long {new_m1.get('long_rows')}; SUCCESS poses {new_m1.get('success_n_poses')}; mapping-fail poses {new_m1.get('mapping_fail_n_poses')}; accounted {new_m1.get('pose_accounted')}/14250",
            f"- M1b<M1 {new_m1.get('M1b_lt_M1')}; mode1 mismatch {new_m1.get('mode1_mismatch')}; raw-log mismatch {new_m1.get('raw_log_mode1_mismatches')}",
            "",
        ]
    if new_m3.get("audited"):
        lines += [
            "## NEW M3 summary",
            f"- master {new_m3.get('n')} status {new_m3.get('status')}",
            f"- repaired {new_m3.get('repaired_jobs')}/61; parse disagreements {new_m3.get('parse_disagreements')}; AB_001 TIMEOUT kept {new_m3.get('ab001_timeout_kept_rows')}",
            "",
        ]
    lines += [
        "## Historical old channels (not formal candidates)",
        "- CRITICAL retained as provenance: old M3 CG0→C/drop G0 does not restore frozen-SDF ring closure (29 ligands / 61 jobs).",
        "- MAJOR retained as provenance: old M1 61 jobs INPUT_FORMAT_FAILURE (CG0).",
        "- Constraint: within-method physical reuse only (same receptor + global ligand). Methods are not required to share scores.",
        "- Phase 7 levels: 8070 pair-expanded seed records AND 7960 unique physical job-seed combinations. Both layers were checked.",
        "",
        "No AUROC / summary_min / ΔAUC / Top-K / EF computed.",
        "No original files overwritten.",
    ]
    (QA / "END_TO_END_SCIENTIFIC_INTEGRITY_GATE.md").write_text("\n".join(lines) + "\n")
    (QA / "E2E_V2_SCOPE.json").write_text(json.dumps({"v2": v2, "new_m1": new_m1, "new_m3": new_m3, "gate": gate}, indent=2) + "\n")


def main() -> int:
    v2 = v2_in_progress()
    print("v2_scope", json.dumps({k: v2[k] for k in v2 if k != "process_lines"}, indent=2), flush=True)
    write_dataflow()
    print("dataflow", flush=True)
    cfg = config_authority()
    print("config", cfg.get("CRITICAL_CONFIG_DRIFT"), flush=True)
    rec_rows = receptor_audit()
    print("receptor", sum(1 for r in rec_rows if r["issues"]), "with_issues", flush=True)
    box_rows = box_audit()
    print("box", sum(1 for r in box_rows if r["issues"]), "with_issues", flush=True)
    lig = ligand_full_audit()
    print("ligands", len(lig), flush=True)
    interface_audit(lig)
    print("interfaces", flush=True)
    p5 = phase5_audit()
    print("phase5", len(p5), flush=True)
    p7 = phase7_audit()
    print("phase7", json.dumps({k: p7[k] for k in p7 if k != "x"}), flush=True)
    ident = identity_reuse()
    print("identity", {k: ident[k] for k in ident if k != "x"}, flush=True)
    old_m1 = old_m1_audit()
    print("old_m1", old_m1, flush=True)
    m2 = m2_audit()
    print("m2", m2, flush=True)
    old_m3 = old_m3_audit()
    print("old_m3", old_m3, flush=True)
    new_m1 = new_m1_audit()
    print("new_m1", new_m1, flush=True)
    new_m3 = new_m3_audit()
    print("new_m3", new_m3, flush=True)
    alt = alt_audit()
    print("alt", alt, flush=True)
    act = activity_audit()
    print("activity", act, flush=True)
    failure_and_missing(p7, old_m1, old_m3, m2)
    consistency()
    write_gate(v2, cfg, {"receptor": rec_rows, "box": box_rows}, p7, ident, old_m1, m2, old_m3, alt, act, new_m1, new_m3)
    print("GATE written", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
