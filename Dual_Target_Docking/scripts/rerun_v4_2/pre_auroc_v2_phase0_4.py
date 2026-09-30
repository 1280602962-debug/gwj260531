#!/usr/bin/env python3
"""PRE_AUROC V2 phases 0-4: authority, counts, ligand representation, M3 topology, M1 61-class.
No AUROC. Does not overwrite original masters.
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import (  # noqa: E402
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
    has_branch_edge,
    invert_smiles_idx,
    load_frozen_sdf,
    map_serial_to_sdf_heavy,
    parse_pdbqt_atoms,
    parse_smiles_idx,
    parse_torsion_tree,
    ring_count,
    same_rigid_block,
    sdf_has_bond,
    sdf_heavy,
    special_types,
    split_models,
)

MEEKO_BIN = Path("/home/gwj/miniconda3/bin/mk_prepare_ligand.py")
SCRIPTS = ROOT / "scripts" / "rerun_v4_2"
M1_MASTER = SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv"
M1_WORK = SEC / "gnina_rescore"
MANIFEST = SEC / "rescoring_job_manifest.csv"
M3_MASTER = GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv"
M3_JOBS = GIND / "jobs"
FREEZE_FINAL = PROTO / "SCORING_DOCKING_ABLATION_FREEZE_FINAL.yaml"
ANALYSIS_FREEZE = PROTO / "analysis_freeze.yaml"


def grep_scripts() -> dict:
    hits = []
    for p in sorted(SCRIPTS.glob("*.py")):
        text = p.read_text(encoding="utf-8", errors="replace")
        if "analysis_freeze" in text or "independent_docking_in_v4_1" in text:
            hits.append({
                "script": str(p.relative_to(ROOT)),
                "reads_analysis_freeze_path": "analysis_freeze.yaml" in text or "analysis_freeze" in text,
                "mentions_old_independent_flag": "independent_docking_in_v4_1" in text,
                "mentions_allowed_pdb": "allowed_pdb" in text,
            })
    current_runners = {
        "ablation_phase4_gnina_rescore.py",
        "ablation_phase578_rtmscore.py",
        "ablation_phase10_12_gnina.py",
        "phase7_fiveseed_production.py",
        "ablation_phase23_masters.py",
        "ablation_config.py",
    }
    runner_reads_old = []
    for name in current_runners:
        p = SCRIPTS / name
        text = p.read_text(encoding="utf-8", errors="replace") if p.is_file() else ""
        if "analysis_freeze.yaml" in text or "independent_docking_in_v4_1" in text:
            runner_reads_old.append(name)
    return {
        "scripts_mentioning_analysis_freeze": hits,
        "current_m0_m3_runners_reading_old_fields": runner_reads_old,
        "authority_status": (
            "STALE_PROVENANCE_ONLY" if not runner_reads_old else "AUTHORITY_CONFLICT"
        ),
    }


def meeko_help() -> dict:
    proc = subprocess.run([str(MEEKO_BIN), "-h"], capture_output=True, text=True)
    help_text = (proc.stdout or "") + (proc.stderr or "")
    has_flag = "--rigid_macrocycles" in help_text
    ver_proc = subprocess.run(
        ["/home/gwj/miniconda3/bin/python", "-c", "import meeko; print(meeko.__version__)"],
        capture_output=True, text=True,
    )
    return {
        "mk_prepare_ligand": str(MEEKO_BIN),
        "meeko_version": (ver_proc.stdout or "").strip(),
        "has_rigid_macrocycles": has_flag,
        "help_excerpt": "\n".join(
            ln for ln in help_text.splitlines() if "macrocycle" in ln.lower() or "rigid_macro" in ln.lower()
        ),
    }


def physical_jobs() -> list[dict]:
    return list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig", newline="")))


def pair_map_independent() -> list[dict]:
    return [
        r for r in csv.DictReader(MAP.open(encoding="utf-8-sig", newline=""))
        if r["pair"] in RECEPTORS and r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"
    ]


def count_hierarchy() -> dict:
    jobs = physical_jobs()
    n_phys = len(jobs)
    n_elig = sum(1 for r in jobs if r.get("rescoring_eligible") == "1")
    n_no_pose = n_phys - n_elig
    pose_hist = Counter()
    for r in jobs:
        if r.get("rescoring_eligible") != "1" or not r.get("pose_path"):
            continue
        n = len(split_models(Path(r["pose_path"]).read_text(errors="replace")))
        pose_hist[n] += 1
    n9 = pose_hist.get(9, 0)
    n8 = pose_hist.get(8, 0)
    other = {k: v for k, v in sorted(pose_hist.items()) if k not in (8, 9)}
    expected_long = n9 * 9 + n8 * 8 + sum(k * v for k, v in other.items())
    m0 = list(csv.DictReader((QA / "official_primary_seed42_score_master_8pair.csv").open(encoding="utf-8-sig")))
    m0_status = Counter(r["vina_status"] for r in m0)
    timeouts = [
        (r["canonical_ligand_id"], r["pdb_id"])
        for r in m0 if r["vina_status"] == "TIMEOUT"
    ]
    return {
        "level": {
            "pair_expanded": {
                "rows": len(m0),
                "SUCCESS": m0_status.get("SUCCESS", 0),
                "TIMEOUT": m0_status.get("TIMEOUT", 0),
                "other": {k: v for k, v in m0_status.items() if k not in {"SUCCESS", "TIMEOUT"}},
            },
            "physical": {
                "jobs": n_phys,
                "with_vina_pose": n_elig,
                "MISSING_NO_VINA_POSE": n_no_pose,
                "n_jobs_9_poses": n9,
                "n_jobs_8_poses": n8,
                "other_pose_counts": other,
                "expected_physical_long_rows": expected_long,
                "formula": f"{n9}×9 + {n8}×8",
            },
        },
        "pair_expanded_timeouts": timeouts,
        "forbidden_expectation_1584x9_plus_24x8": "DO_NOT_USE_that_is_pair_expanded_layer",
    }


def ligand_receptors(jobs: list[dict]) -> dict[str, set[str]]:
    out = defaultdict(set)
    for r in jobs:
        out[r["global_ligand_entity_id"]].add(r["pdb_id"])
    return out


def ligand_physical_jobs(jobs: list[dict], eligible_only: bool = False) -> dict[str, list[str]]:
    out = defaultdict(list)
    for r in jobs:
        if eligible_only and r.get("rescoring_eligible") != "1":
            continue
        out[r["global_ligand_entity_id"]].append(f"{r['pdb_id']}__{r['global_ligand_entity_id']}")
    return out


def audit_one_ligand(eid: str, recs: set[str], m1_jobs: list[str], m3_jobs: list[str]) -> dict:
    sdf_p = LIG_SDF / f"{eid}.sdf"
    pdbqt_p = LIG_PDBQT / f"{eid}.pdbqt"
    m3_in_p = LIG_GNINA_PDBQT / f"{eid}.pdbqt"
    sdf = load_frozen_sdf(sdf_p)
    frozen_text = pdbqt_p.read_text(errors="replace") if pdbqt_p.is_file() else ""
    atoms = parse_pdbqt_atoms(frozen_text) if frozen_text else []
    spec = special_types(atoms)
    cgn = [a for a in atoms if a["is_real_macro"]]
    gn = [a for a in atoms if a["is_glue"]]
    smiles, idx = parse_smiles_idx(frozen_text) if frozen_text else ("", {})
    inv = invert_smiles_idx(idx) if idx else {}
    cgn_in_idx = sum(1 for a in cgn if a["serial"] in inv)
    gn_in_idx = sum(1 for a in gn if a["serial"] in inv)
    m3_atoms = parse_pdbqt_atoms(m3_in_p.read_text(errors="replace")) if m3_in_p.is_file() else []
    m3_spec = special_types(m3_atoms)
    has_special = "YES" if (cgn or gn or spec) else "NO"
    notes = []
    if cgn and cgn_in_idx != len(cgn):
        notes.append(f"CGn_not_all_in_SMILES_IDX:{cgn_in_idx}/{len(cgn)}")
    if gn and gn_in_idx:
        notes.append(f"Gn_incorrectly_in_SMILES_IDX:{gn_in_idx}")
    if has_special == "YES" and not cgn and gn:
        notes.append("glue_without_real_macrocycle_typed_atom")
    if has_special == "YES" and cgn and not gn:
        notes.append("real_macrocycle_typed_atom_without_glue")
    return {
        "global_ligand_entity_id": eid,
        "frozen_sdf_exists": int(sdf_p.is_file()),
        "frozen_pdbqt_exists": int(pdbqt_p.is_file()),
        "has_macrocycle_special_type": has_special,
        "special_atom_types": ";".join(spec),
        "real_macrocycle_typed_atoms": ";".join(sorted({a["type"] for a in cgn})) if cgn else "",
        "glue_pseudoatom_types": ";".join(sorted({a["type"] for a in gn})) if gn else "",
        "n_special_atoms": len(cgn) + len(gn),
        "n_real_macrocycle_typed_atoms": len(cgn),
        "n_glue_pseudoatoms": len(gn),
        "n_CGn_in_SMILES_IDX": cgn_in_idx,
        "n_Gn_in_SMILES_IDX": gn_in_idx,
        "affected_receptors": ";".join(sorted(recs)),
        "affected_physical_jobs_M1": ";".join(m1_jobs),
        "n_affected_physical_jobs_M1": len(m1_jobs),
        "affected_physical_jobs_M3": ";".join(m3_jobs),
        "n_affected_physical_jobs_M3": len(m3_jobs),
        "original_frozen_sdf_ring_count": ring_count(sdf),
        "frozen_sdf_formal_charge": formal_charge(sdf),
        "frozen_sdf_heavy_atoms": sdf_heavy(sdf).GetNumAtoms() if sdf is not None else "",
        "m3_converted_input_exists": int(m3_in_p.is_file()),
        "m3_converted_special_types": ";".join(m3_spec),
        "notes": "|".join(notes),
    }


def topology_one(eid: str) -> dict:
    sdf_p = LIG_SDF / f"{eid}.sdf"
    pdbqt_p = LIG_PDBQT / f"{eid}.pdbqt"
    m3_in_p = LIG_GNINA_PDBQT / f"{eid}.pdbqt"
    rec = {
        "global_ligand_entity_id": eid,
        "decision": "",
        "n_broken_ring_closure_pairs": 0,
        "broken_ring_closure_pairs": "",
        "frozen_sdf_has_closure_bonds": "",
        "frozen_pdbqt_tree_has_closure_edges": "",
        "m3_tree_has_closure_edges": "",
        "m3_real_heavy_count_match": "",
        "m3_formal_charge_note": "not_inferred_from_pdbqt",
        "cartesian_distances_A_note_only": "",
        "evidence": "",
    }
    if not pdbqt_p.is_file() or not sdf_p.is_file() or not m3_in_p.is_file():
        rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
        rec["evidence"] = "missing_frozen_or_m3_input"
        return rec
    frozen_text = pdbqt_p.read_text(errors="replace")
    m3_text = m3_in_p.read_text(errors="replace")
    sdf = load_frozen_sdf(sdf_p)
    if sdf is None:
        rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
        rec["evidence"] = "frozen_sdf_unreadable"
        return rec
    atoms = parse_pdbqt_atoms(frozen_text)
    pairs = complementary_cgn_pairs(atoms)
    rec["n_broken_ring_closure_pairs"] = len(pairs)
    serial_to_sdf, map_why = map_serial_to_sdf_heavy(frozen_text, sdf)
    heavy = sdf_heavy(sdf)
    frozen_tree = parse_torsion_tree(frozen_text)
    m3_tree = parse_torsion_tree(m3_text)
    m3_atoms = parse_pdbqt_atoms(m3_text)
    m3_glue = [a for a in m3_atoms if a["is_glue"]]
    m3_macro = [a for a in m3_atoms if a["is_real_macro"]]
    frozen_real = [a for a in atoms if not a["is_glue"]]
    # heavy-atom comparison: frozen real non-H vs M3 non-H
    def n_heavy(atom_list):
        return sum(1 for a in atom_list if a["type"] not in {"H", "HD", "HS"})

    rec["m3_real_heavy_count_match"] = int(
        n_heavy(frozen_real) == n_heavy(m3_atoms) and not m3_glue and not m3_macro
    )
    sdf_ok = []
    frozen_tree_hit = []
    m3_tree_hit = []
    pair_desc = []
    dists = []
    if not pairs:
        rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
        rec["evidence"] = "cannot_pair_CGn_via_Gn_provenance|" + map_why
        return rec
    if not serial_to_sdf:
        rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
        rec["evidence"] = "no_unique_atom_mapping|" + map_why
        return rec
    for p in pairs:
        a, b = p["serial_a"], p["serial_b"]
        ia, ib = serial_to_sdf.get(a), serial_to_sdf.get(b)
        if ia is None or ib is None:
            rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
            rec["evidence"] = f"CGn_serial_not_in_SMILES_IDX:{a},{b}|{map_why}"
            return rec
        has_sdf = sdf_has_bond(heavy, ia, ib)
        sdf_ok.append(has_sdf)
        ft = has_branch_edge(frozen_tree, a, b) or same_rigid_block(frozen_tree, a, b)
        mt = has_branch_edge(m3_tree, a, b) or same_rigid_block(m3_tree, a, b)
        frozen_tree_hit.append(ft)
        m3_tree_hit.append(mt)
        pair_desc.append(f"pdbqt:{a}-{b}->sdf:{ia}-{ib}")
        dists.append(str(p["cartesian_distance_A"]))
    rec["broken_ring_closure_pairs"] = ";".join(pair_desc)
    rec["frozen_sdf_has_closure_bonds"] = int(all(sdf_ok))
    rec["frozen_pdbqt_tree_has_closure_edges"] = int(any(frozen_tree_hit))
    rec["m3_tree_has_closure_edges"] = int(any(m3_tree_hit))
    rec["cartesian_distances_A_note_only"] = ";".join(dists)
    if all(sdf_ok) and all(m3_tree_hit):
        rec["decision"] = "M3_INPUT_TOPOLOGY_PASS"
        rec["evidence"] = (
            "frozen_SDF_closure_bond_present_AND_M3_torsion_tree_has_same_edge_or_same_rigid_block;"
            "cartesian_distance_not_used_as_proof"
        )
    elif all(sdf_ok) and not any(m3_tree_hit):
        rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
        rec["evidence"] = (
            "frozen_SDF_has_ring_closure_bond; Meeko CGn/Gn provenance identifies the broken pair; "
            "historical CG0->C drop-G0 left no BRANCH edge and the two real carbons are not in one rigid block; "
            "cartesian_distance_not_used_as_proof"
        )
    else:
        rec["decision"] = "M3_INPUT_TOPOLOGY_FAIL"
        rec["evidence"] = (
            f"cannot_prove_equivalence sdf_bonds={sdf_ok} m3_tree={m3_tree_hit} map={map_why}"
        )
    return rec


def classify_m1_fail(row: dict) -> dict:
    pdb, eid = row["pdb_id"], row["global_ligand_entity_id"]
    jobdir = M1_WORK / f"{pdb}__{eid}"
    status_p = jobdir / "status.json"
    status = {}
    if status_p.is_file():
        status = json.loads(status_p.read_text())
    logs = sorted(jobdir.glob("vina_mode*.gnina.log"))
    log_text = "\n".join(p.read_text(errors="replace") for p in logs)
    pose_files = sorted(jobdir.glob("vina_mode*.pdbqt"))
    n_pose_inputs = len(pose_files)
    type_err = bool(re.search(r"is not a valid AutoDock type", log_text))
    cg_err = bool(re.search(r'"(CG\d+|G\d+)" is not a valid AutoDock type', log_text))
    rc_nonzero = "Parse error" in log_text or re.search(r"\n[1-9]\d*\n*$", log_text) is not None
    cnn_vals = re.findall(r"CNNscore[:\s]+([0-9eE.+-]+)", log_text)
    finite_cnn = [v for v in cnn_vals if finite_float(v) and v != "--"]
    n_retries = status.get("n_retries", row.get("n_retries", ""))
    if type_err or cg_err:
        cat = "GNINA_INPUT_TYPE_FAILURE"
    elif finite_cnn and len(finite_cnn) >= n_pose_inputs:
        cat = "PARSER_RECOVERABLE" if (row.get("M1_CNNscore") in ("", None) or row.get("job_status") != "SUCCESS") else "OTHER"
    elif rc_nonzero and n_pose_inputs and all((jobdir / f"vina_mode{i}.pdbqt").is_file() for i in range(1, n_pose_inputs + 1)):
        # output pose files exist because they are inputs; check if gnina wrote scores
        cat = "TRUE_GNINA_FAIL" if not finite_cnn else "RC_NONZERO_BUT_OUTPUT_COMPLETE"
    elif not logs:
        cat = "OTHER"
    else:
        cat = "TRUE_GNINA_FAIL" if not finite_cnn else "OTHER"
    return {
        "pdb_id": pdb,
        "global_ligand_entity_id": eid,
        "pairs": row.get("pairs", ""),
        "vina_status": row.get("vina_status", ""),
        "n_poses_master": row.get("n_poses", ""),
        "n_pose_inputs": n_pose_inputs,
        "n_logs": len(logs),
        "job_status": row.get("job_status", status.get("job_status", "")),
        "reason": row.get("reason", status.get("reason", "")),
        "n_retries": n_retries,
        "status_json_exists": int(status_p.is_file()),
        "type_error": int(type_err),
        "cgn_or_gn_type_error": int(cg_err),
        "n_finite_cnn_in_logs": len(finite_cnn),
        "same_input_retry_forbidden": int(str(n_retries) == "1"),
        "category": cat,
        "log_excerpt": next(
            (ln.strip() for ln in log_text.splitlines() if "valid AutoDock type" in ln or "Parse error" in ln),
            "",
        ),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    QA.mkdir(parents=True, exist_ok=True)
    auth = grep_scripts()
    meeko = meeko_help()
    gate0 = {
        "freeze_final_exists": FREEZE_FINAL.is_file(),
        "analysis_freeze_exists": ANALYSIS_FREEZE.is_file(),
        "analysis_freeze_old_fields_status": "STALE_PROVENANCE_ONLY",
        "old_fields": [
            "gnina.independent_docking_in_v4_1: false",
            "alternative_receptors.allowed_pdb: [4JPS, 5DXT, 4JSX]",
        ],
        **auth,
        "meeko": meeko,
        "STOP": False,
        "STOP_reason": "",
    }
    if auth["authority_status"] == "AUTHORITY_CONFLICT":
        gate0["STOP"] = True
        gate0["STOP_reason"] = "current_execution_script_reads_analysis_freeze_old_fields"
    if not meeko["has_rigid_macrocycles"]:
        gate0["STOP"] = True
        gate0["STOP_reason"] = (gate0["STOP_reason"] + "|").lstrip("|") + "meeko_lacks_rigid_macrocycles"
    (QA / "PRE_AUROC_V2_PHASE0_AUTHORITY.json").write_text(json.dumps(gate0, indent=2) + "\n")
    print(json.dumps({"phase0": gate0["authority_status"], "meeko": meeko["has_rigid_macrocycles"], "STOP": gate0["STOP"]}, indent=2), flush=True)
    if gate0["STOP"]:
        print("STOP_AUTHORITY_OR_MEEKO", flush=True)
        return 2

    counts = count_hierarchy()
    (QA / "PRE_AUROC_V2_COUNT_HIERARCHY.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(json.dumps(counts["level"], indent=2), flush=True)

    jobs = physical_jobs()
    recs = ligand_receptors(jobs)
    m1_by = ligand_physical_jobs(jobs, eligible_only=True)
    m3_by = ligand_physical_jobs(jobs, eligible_only=False)
    eids = sorted({p.stem for p in LIG_PDBQT.glob("*.pdbqt")} | {p.stem for p in LIG_SDF.glob("*.sdf")})
    lig_rows = []
    for i, eid in enumerate(eids, 1):
        lig_rows.append(audit_one_ligand(eid, recs.get(eid, set()), m1_by.get(eid, []), m3_by.get(eid, [])))
        if i % 100 == 0:
            print(f"ligand_audit {i}/{len(eids)}", flush=True)
    write_csv(QA / "GNINA_LIGAND_REPRESENTATION_AUDIT.csv", lig_rows)
    special = [r for r in lig_rows if r["has_macrocycle_special_type"] == "YES"]
    type_counter = Counter()
    for r in special:
        for t in r["special_atom_types"].split(";"):
            if t:
                type_counter[t] += 1
    md = []
    md.append("# GNINA ligand representation audit")
    md.append("")
    md.append("Level: frozen ligand assets, not performance.")
    md.append("CGn = real carbon (or other real atom) with macrocycle-specific type. Gn = glue pseudoatom only.")
    md.append("")
    md.append(f"- Frozen ligand entities scanned: {len(lig_rows)}")
    md.append(f"- Missing frozen SDF: {sum(1 for r in lig_rows if r['frozen_sdf_exists'] != 1)}")
    md.append(f"- Missing frozen PDBQT: {sum(1 for r in lig_rows if r['frozen_pdbqt_exists'] != 1)}")
    md.append(f"- has_macrocycle_special_type=YES: {len(special)}")
    md.append(f"- Observed special types (ligand counts): {dict(type_counter)}")
    md.append(f"- Ligands where every CGn is in SMILES IDX: {sum(1 for r in special if r['n_CGn_in_SMILES_IDX'] == r['n_real_macrocycle_typed_atoms'])}")
    md.append(f"- Ligands where any Gn is in SMILES IDX (should be 0): {sum(1 for r in special if int(r['n_Gn_in_SMILES_IDX']) > 0)}")
    md.append("")
    md.append("No AUROC or other performance metric was read.")
    (QA / "GNINA_LIGAND_REPRESENTATION_AUDIT.md").write_text("\n".join(md) + "\n")
    print(f"special_ligands {len(special)} types {dict(type_counter)}", flush=True)

    topo_rows = []
    for r in special:
        topo_rows.append(topology_one(r["global_ligand_entity_id"]))
    if topo_rows:
        write_csv(QA / "M3_MACROCYCLE_TOPOLOGY_DECISION.csv", topo_rows)
        dec = Counter(r["decision"] for r in topo_rows)
        print("m3_topology", dict(dec), flush=True)
    else:
        print("m3_topology none_affected", flush=True)

    m1 = list(csv.DictReader(M1_MASTER.open(encoding="utf-8-sig")))
    fails = [r for r in m1 if r.get("reason") == "parse_or_rc" or r.get("job_status") == "TECHNICAL_FAIL"]
    m1_aud = [classify_m1_fail(r) for r in fails]
    if m1_aud:
        write_csv(QA / "M1_PARSE_OR_RC_AUDIT.csv", m1_aud)
        cats = Counter(r["category"] for r in m1_aud)
        md2 = [
            "# M1 parse_or_rc audit",
            "",
            "Classification only. Same-input retry is forbidden when n_retries=1.",
            "These rows are not the final M1 missing set.",
            "",
            f"- Rows classified: {len(m1_aud)}",
            f"- Categories: {dict(cats)}",
            f"- n_retries=1: {sum(1 for r in m1_aud if str(r['n_retries']) == '1')}",
            f"- GNINA_INPUT_TYPE_FAILURE: {cats.get('GNINA_INPUT_TYPE_FAILURE', 0)}",
            "",
            "No AUROC was computed.",
        ]
        (QA / "M1_PARSE_OR_RC_AUDIT.md").write_text("\n".join(md2) + "\n")
        print("m1_class", dict(cats), "n", len(m1_aud), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
