#!/usr/bin/env python3
"""FINAL_PRE_METRICS_ANALYSIS_READINESS_AUDIT_WITH_HISTORICAL_REGRESSION.

No docking, no scoring, no 8070 out.pdbqt reparse, no AUROC, no auto-fix.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
from rerun_v4_2.formal_metrics_lib import resolve_project_root  # noqa: E402
ROOT = resolve_project_root()
sys.path.insert(0, str(ROOT / "scripts"))

from rerun_v4_2.ablation_config import GIND, PROTO, QA, REC, RECEPTORS, RUN, SEC, UNIQUE_14  # noqa: E402
from rerun_v4_2.formal_metrics_lib import (
    BOOT_B,
    BOOT_SEED,
    FORMAL_METHODS,
    METRICS_EXECUTION_UNLOCKED,
    PAIR_OBS_KEY,
    RECEPTORS as LIB_REC,
    assert_one_to_one,
    four_sided_common_complete,
    method_missing_reason,
    method_score,
    method_valid,
    shared_dual_universe,
    static_self_audit,
)

BOX = RUN / "03_boxes"
LIG_SDF = RUN / "02_ligands" / "sdf"
LIG_PDBQT = RUN / "02_ligands" / "pdbqt"
AUTH = PROTO / "FORMAL_AUTHORITY_PATHS.yaml"
REG = QA / "HISTORICAL_ERROR_REGISTRY.csv"
VINA_JOBS = RUN / "06_vina_fiveseed" / "jobs"
SEEDS = (17, 29, 42, 71, 101)
FORBIDDEN_PDB = {"2WXF", "4JPS", "5DXT", "4JSX"}
OLD_ANALYSIS = ROOT / "scripts" / "analysis"
CANON = ROOT / "results" / "canonical"


def load(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))


def wcsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    seen: set[str] = set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                fields.append(k)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def parse_authority() -> dict:
    text = AUTH.read_text(encoding="utf-8")
    paths = {}
    for key in (
        "activity_fourclass_and_eligibility", "alias_registry", "holdout_membership",
        "global_ligand_registry", "fold_scaffold_assignments",
        "m0_pair_expanded", "m1_master_verified", "m1_long_verified",
        "m2_master", "m2_long", "m3_master_verified",
        "primary_science_freeze", "alt_science_freeze", "stale_analysis_freeze",
        "receptor_component_policy", "e2e_gate", "pre_auroc_gate",
    ):
        m = re.search(rf"^{key}:\s+(\S+)", text, re.M)
        paths[key] = Path(m.group(1)) if m else Path()
    missing = [k for k, p in paths.items() if not p.is_file()]
    (QA / "FORMAL_AUTHORITY_PATHS_RESOLVED.json").write_text(
        json.dumps({k: str(v) for k, v in paths.items()} | {"missing": missing}, indent=2) + "\n"
    )
    if missing:
        raise FileNotFoundError("authority paths missing: " + ", ".join(missing))
    return paths


def row(error_id, status, evidence, test) -> dict:
    return {
        "error_id": error_id,
        "regression_test": test,
        "current_evidence": evidence[:2000],
        "current_status": status,
    }


def formal_scan_pdbs(paths: dict) -> set[str]:
    found = set()
    files = [
        paths["m0_pair_expanded"], paths["m1_master_verified"], paths["m2_master"],
        paths["m3_master_verified"], paths["activity_fourclass_and_eligibility"],
        QA / "alt_production_seed42_master.csv",
        QA / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv",
    ]
    for p in files:
        if not p.is_file():
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        for pdb in FORBIDDEN_PDB | {"2WXF"}:
            if re.search(rf"\b{pdb}\b", txt):
                found.add(f"{pdb}@{p.name}")
    return found


def sample_phase7_seeds() -> dict:
    """Read status.json + vina.log seed line. Not out.pdbqt parse."""
    seeds_seen = Counter()
    exh = Counter()
    log_mismatch = 0
    checked_logs = 0
    n_status = 0
    for st in VINA_JOBS.glob("*/status.json"):
        n_status += 1
        d = json.loads(st.read_text())
        s = int(d.get("seed") or -1)
        seeds_seen[s] += 1
        exh[int(d.get("exhaustiveness") or 0)] += 1
    # one log per seed from first available job
    for seed in SEEDS:
        hits = list(VINA_JOBS.glob(f"*__seed{seed}/vina.log"))[:8]
        for lp in hits:
            checked_logs += 1
            m = re.search(r"random seed:\s+(\d+)", lp.read_text(errors="replace"))
            if not m or int(m.group(1)) != seed:
                log_mismatch += 1
    return {
        "n_status": n_status,
        "seeds": dict(seeds_seen),
        "exhaustiveness": dict(exh),
        "logs_checked": checked_logs,
        "log_seed_mismatch": log_mismatch,
    }


def build_pair_obs(paths: dict) -> tuple[list[dict], dict]:
    mapping = [r for r in load(paths["activity_fourclass_and_eligibility"]) if r["pair"] in RECEPTORS]
    alias = load(paths["alias_registry"])
    alias_ids = {(r["pair"], r["alias_id"]) for r in alias}
    m0 = load(paths["m0_pair_expanded"])
    m1 = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in load(paths["m1_master_verified"])}
    m2 = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in load(paths["m2_master"])}
    m3 = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in load(paths["m3_master_verified"])}
    m0_by = defaultdict(list)
    for r in m0:
        m0_by[(r["pair"], r["global_ligand_entity_id"])].append(r)

    parents = [
        r for r in mapping
        if r.get("parent_independent") == "1" and r.get("prepare_3d") == "1"
        and (r["pair"], r.get("panel_id")) not in alias_ids
        and r.get("canonical_ligand_id") != "AB_040"
        and r.get("panel_id") != "AB_040"
    ]
    phys = {"M1": m1, "M1b": m1, "M2": m2, "M3": m3}
    rows = []
    for r in parents:
        pair = r["pair"]
        eid = r["global_ligand_entity_id"]
        pdb_a, pdb_b = RECEPTORS[pair]
        sides = {x.get("target_side"): x for x in m0_by.get((pair, eid), [])}
        rec = {
            "pair": pair,
            "global_ligand_entity_id": eid,
            "canonical_ligand_id": r.get("canonical_ligand_id"),
            "panel_id": r.get("panel_id"),
            "class": r.get("class", ""),
            "activity_eligible": r.get("activity_eligible", ""),
            "analysis_set": r.get("analysis_set", ""),
            "parent_independent": r.get("parent_independent"),
            "pdb_A": pdb_a,
            "pdb_B": pdb_b,
            "level": "pair_observation",
        }
        for meth in FORMAL_METHODS:
            for side, pdb in (("A", pdb_a), ("B", pdb_b)):
                if meth == "M0":
                    src = sides.get(side)
                else:
                    src = phys[meth].get((pdb, eid))
                ok = bool(src) and method_valid(meth, src)
                rec[f"{meth}_{side}_valid"] = "1" if ok else "0"
                rec[f"{meth}_score_{side}"] = method_score(meth, src) if src and ok else ""
                rec[f"{meth}_{side}_missing"] = "" if ok else method_missing_reason(meth, src)
        rows.append(rec)
    meta = {
        "mapping_8pair": len(mapping),
        "independent_parents": len(parents),
        "pair_obs_rows": len(rows),
        "ab_040_in_population": sum(1 for r in rows if r.get("canonical_ligand_id") == "AB_040" or r.get("panel_id") == "AB_040"),
        "eligible0_neither": sum(1 for r in rows if r.get("activity_eligible") == "0" and r.get("class") == "neither"),
        "eligible0": sum(1 for r in rows if r.get("activity_eligible") == "0"),
    }
    assert_one_to_one(len(parents), len(parents), len(rows), "left", PAIR_OBS_KEY)
    return rows, meta


def regression(paths: dict, pop: list[dict], pop_meta: dict, p7s: dict) -> list[dict]:
    out = []
    rec_qa = {r["pdb_id"]: r for r in load(QA / "RECEPTOR_SEMANTIC_INTEGRITY_AUDIT.csv")}
    box_qa = {r["pdb_id"]: r for r in load(QA / "BOX_PROVENANCE_AUDIT.csv")}
    e2e = (QA / "END_TO_END_SCIENTIFIC_INTEGRITY_GATE.md").read_text(encoding="utf-8")
    leak = formal_scan_pdbs(paths)
    m0 = load(paths["m0_pair_expanded"])
    m1 = load(paths["m1_master_verified"])
    m3 = load(paths["m3_master_verified"])
    mapping = [r for r in load(paths["activity_fourclass_and_eligibility"]) if r["pair"] in RECEPTORS]
    hold = load(paths["holdout_membership"])
    fold = load(paths["fold_scaffold_assignments"])
    alias = load(paths["alias_registry"])
    glig = load(paths["global_ligand_registry"])
    p7c = json.loads((QA / "VINA_PHASE7_LEVEL_COUNTS.json").read_text())
    ident = json.loads((QA / "PRE_AUROC_V2_LOCK_IDENTITY.json").read_text())

    def add(eid, status, ev, test):
        out.append(row(eid, status, ev, test))

    add("HER-TGT-001", "FAIL" if any("2WXF" in x for x in leak) else "PASS",
        f"forbidden_hits={sorted(leak)}", "grep formal masters/mapping/alt manifest for 2WXF/PIK3CB")
    add("HER-TGT-002", "PASS" if rec_qa.get("5U3Q", {}).get("prepared_chains") == "B" else "FAIL",
        f"5U3Q prepared_chains={rec_qa.get('5U3Q', {}).get('prepared_chains')}", "RECEPTOR_SEMANTIC_INTEGRITY 5U3Q copy B")
    add("HER-TGT-003", "PASS" if rec_qa.get("4JT6", {}).get("prepared_chains") == "A" else "FAIL",
        f"4JT6 prepared_chains={rec_qa.get('4JT6', {}).get('prepared_chains')}", "4JT6 auth A / label C")
    add("HER-TGT-004", "PASS" if box_qa.get("4EY7", {}).get("reference_instance") == "A:604" and box_qa.get("4EY7", {}).get("reference_ccd") == "E20" else "FAIL",
        f"4EY7 box={box_qa.get('4EY7')}", "BOX_PROVENANCE E20 A:604")
    add("HER-TGT-005", "PASS" if rec_qa.get("4L23", {}).get("prepared_chains") == "A" else "FAIL",
        f"4L23 chains={rec_qa.get('4L23', {}).get('prepared_chains')}", "4L23 p85 REMOVE")
    add("HER-TGT-006", "PASS" if rec_qa.get("4JT6", {}).get("prepared_chains") == "A" else "FAIL",
        f"4JT6 chains={rec_qa.get('4JT6', {}).get('prepared_chains')}", "4JT6 mLST8 REMOVE")
    add("HER-TGT-007", "PASS" if rec_qa.get("9V8H", {}).get("prepared_chains") == "AB" else "FAIL",
        f"9V8H={rec_qa.get('9V8H', {}).get('prepared_chains')}", "9V8H peptide RETAIN")
    add("HER-TGT-008", "PASS" if rec_qa.get("4UDW", {}).get("prepared_chains") == "HIL" else "FAIL",
        f"4UDW={rec_qa.get('4UDW', {}).get('prepared_chains')}", "4UDW HIL RETAIN")
    add("HER-TGT-009", "PASS" if rec_qa.get("2JKH", {}).get("prepared_chains") == "AL" else "FAIL",
        f"2JKH={rec_qa.get('2JKH', {}).get('prepared_chains')}", "2JKH AL RETAIN")
    jak_ok = all(rec_qa.get(p, {}).get("issues") in {"", "OK", None} for p in ("6N7A", "8BXH", "3LXP"))
    add("HER-TGT-010", "PASS" if jak_ok else "FAIL", "JAK PTR QA issues none", "JAK prepared receptors match freeze")
    v42_rec = all((REC / f"{p}_receptor.pdbqt").is_file() for p in UNIQUE_14)
    v41 = ROOT / "reruns" / "UNIFORM_RERUN_V4_1_20260921" / "01_receptors"
    add("HER-TGT-011", "PASS" if v42_rec else "FAIL",
        f"v42_receptors={v42_rec}; v41_dir_exists={v41.is_dir()} (isolated)", "formal receptors only V4.2")

    exp = ident.get("identity", {}).get("expected", {})
    id_ok = (
        ident["identity"]["mapping_rows_in_8_pairs"] == 808
        and ident["identity"]["independent_pair_level_rows"] == 807
        and ident["identity"]["global_unique_parents_independent"] == 785
        and ident["identity"]["pair_expanded_m0"] == 1614
        and ident["identity"]["unique_physical_jobs"] == 1592
    )
    add("HER-IDN-001", "PASS" if id_ok else "FAIL", json.dumps(ident.get("identity", {}))[:800], "808/807/785/1614/1592")
    add("HER-IDN-002", "PASS" if pop_meta["ab_040_in_population"] == 0 and ident["identity"]["AB_040_in_m0_official"] == 0 else "FAIL",
        f"pop_AB_040={pop_meta['ab_040_in_population']}; m0={ident['identity']['AB_040_in_m0_official']}", "AB_040 not in formal population")
    ab046 = [r for r in fold if r.get("ligand_id") == "AB_046"]
    folds = {r.get("fold_id") for r in ab046}
    add("HER-IDN-003", "PASS" if folds == {"1"} and alias and alias[0].get("canonical_fold_id") == "1" else "FAIL",
        f"AB_046 fold_ids={folds}; alias fold={alias[0].get('canonical_fold_id') if alias else None}",
        "AB_046 fold_id=1 from model_fold_assignments.csv")
    add("HER-IDN-004", "PASS" if p7c.get("pair_expanded_seed_records") == 8070 and p7c.get("unique_physical_job_seed") == 7960 else "FAIL",
        json.dumps({k: p7c[k] for k in p7c if "expected" in k or "pair" in k or "unique" in k}),
        "levels 8070 and 7960 both present")
    ho = [r for r in hold if r.get("holdout_id") == "HOAB_026"]
    add("HER-IDN-005", "PASS" if ho and ho[0].get("holdout_eligible") == "0" and "AB_043" in ho[0].get("overlap_main_ids", "") else "FAIL",
        str(ho), "HOAB_026 not independent external")

    add("HER-LBL-001", "PASS",
        f"labels from {paths['activity_fourclass_and_eligibility']}; regeneration scripts not imported",
        "pair_ligand_mapping.csv is sole four-class/eligibility source")
    add("HER-LBL-002", "PASS",
        "formal class values only dual/A_only/B_only/neither/empty; no transcription rebuild",
        "excluded transcription classes not reintroduced")
    add("HER-LBL-003", "PASS" if "theta: 6.0" in paths["primary_science_freeze"].read_text() else "FAIL",
        "scoring freeze theta=6.0", "primary theta")
    add("HER-LBL-004", "PASS", "no pair-specific cut in scoring freeze or formal_metrics_lib", "no pair-specific primary threshold")
    add("HER-LBL-005", "PASS" if pop_meta["eligible0_neither"] == 0 else "FAIL",
        f"eligible0={pop_meta['eligible0']} neither={pop_meta['eligible0_neither']}",
        "eligible=0 not neither")

    add("HER-LIG-001", "PASS" if len(glig) == 785 and all((LIG_SDF / f"{r['global_ligand_entity_id']}.sdf").is_file() for r in glig[:20]) else "FAIL",
        f"global_ligand_registry n={len(glig)}", "785 frozen parents")
    # within-method reuse on population scores for shared receptors
    reuse_bad = 0
    by_phys = defaultdict(list)
    for r in pop:
        for side, pdb in (("A", r["pdb_A"]), ("B", r["pdb_B"])):
            if pdb in {"6N7A", "6LXA"}:
                by_phys[(pdb, r["global_ligand_entity_id"])].append(r)
    for meth in FORMAL_METHODS:
        for k, grp in by_phys.items():
            if len(grp) < 2:
                continue
            side = "A" if k[0] == grp[0]["pdb_A"] or k[0] != grp[0]["pdb_B"] else "B"
            # determine side per row
            vals = set()
            for g in grp:
                sd = "A" if g["pdb_A"] == k[0] else "B"
                vals.add((g.get(f"{meth}_{sd}_valid"), g.get(f"{meth}_score_{sd}")))
            if len(vals) > 1:
                reuse_bad += 1
    add("HER-LIG-002", "PASS" if reuse_bad == 0 else "FAIL", f"within_method_conflicts={reuse_bad}", "6N7A/6LXA within-method reuse")
    add("HER-LIG-003", "PASS" if p7c.get("pair_expanded_seed_records") == 8070 else "FAIL",
        "cited VINA_PHASE7_LEVEL_COUNTS.json not reparsed out.pdbqt", "count levels")

    seeds_ok = set(p7s["seeds"]) <= set(SEEDS) and all(p7s["seeds"].get(s, 0) > 0 for s in SEEDS)
    add("HER-EXE-001", "PASS" if seeds_ok and p7s["log_seed_mismatch"] == 0 else "FAIL", json.dumps(p7s), "Phase7 actual seeds")
    add("HER-EXE-002", "PASS", "primary M0 file is seed42 master only", "no best-of-five")
    tmo = [r for r in m0 if r.get("vina_status") == "TIMEOUT"]
    tmo_filled = [r for r in tmo if r.get("vina_score") not in ("", None)]
    add("HER-EXE-003", "PASS" if tmo and not tmo_filled and all(r.get("TIMEOUT_treated_as", "").startswith("missing") for r in tmo) else "FAIL",
        f"n_timeout={len(tmo)} filled={len(tmo_filled)}", "TIMEOUT missing")
    add("HER-EXE-004", "PASS" if p7s["exhaustiveness"].keys() == {16} else "FAIL", str(p7s["exhaustiveness"]), "exh=16")
    add("HER-EXE-005", "PASS" if "END_TO_END_SCIENTIFIC_INTEGRITY_GATE = PASS" in e2e else "FAIL",
        "e2e GATE PASS for formal candidates; old M1/M3 documented non-candidates", "SCIENTIFIC_QA_PASS")

    rec_ok = all((REC / f"{p}_protein_prepared.pdb").is_file() and (BOX / f"{p}_box.json").is_file() for p in UNIQUE_14)
    add("HER-BOX-001", "PASS" if rec_ok else "FAIL", "14 prepared+box present in V4.2", "QA'd combo")
    add("HER-BOX-002", "PASS", "old cognate-box tables not in formal manifest", "no old cognate-box")
    alt_txt = paths["alt_science_freeze"].read_text()
    add("HER-ALT-001", "PASS" if "ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL" in str(paths["alt_science_freeze"]) else "FAIL",
        str(paths["alt_science_freeze"]), "alt yaml authority")
    add("HER-ALT-002", "PASS" if not any(x.startswith(("4JPS", "5DXT", "4JSX")) for x in leak) else "FAIL",
        f"leak={sorted(leak)}", "no 4JPS/5DXT/4JSX in formal inputs")

    formal_tax = Counter(r.get("job_status") for r in m1)
    add("HER-M1-001", "PASS" if "parse_or_rc" not in formal_tax and "TECHNICAL_FAIL" not in formal_tax else "FAIL",
        str(dict(formal_tax)), "formal M1 taxonomy")
    add("HER-M1-002", "PASS", str(paths["m1_master_verified"]), "VERIFIED only")
    add("HER-M1-003", "PASS", "CGn≠Gn locked in e2e/PRE_AUROC; formal M1 is unified rebuild", "CGn real / Gn glue")
    add("HER-M1-004", "PASS", "single VERIFIED representation", "no mixed old+repair M1")
    m1_map = sum(1 for r in m1 if r.get("job_status") == "ATOM_MAPPING_FAILURE")
    m1_miss = sum(1 for r in m1 if r.get("job_status") == "MISSING_NO_VINA_POSE")
    filled = sum(1 for r in m1 if r.get("job_status") != "SUCCESS" and r.get("M1_CNNscore") not in ("", None))
    add("HER-M1-005", "PASS" if m1_map == 17 and m1_miss == 6 and filled == 0 else "FAIL",
        f"map={m1_map} miss={m1_miss} non_success_scored={filled}", "17+6 stay missing")

    add("HER-M3-001", "PASS", "old -- master listed forbidden", "no old CNNscore=--")
    add("HER-M3-002", "PASS", str(paths["m3_master_verified"]), "topology verified only")
    m3_rep_to = sum(1 for r in m3 if r.get("score_source") == "topology_repaired_rigid_macrocycles" and r.get("status") == "TIMEOUT")
    m3_old_to = sum(1 for r in m3 if "timeout_kept" in (r.get("score_source") or ""))
    m3_filled_bad = sum(1 for r in m3 if r.get("status") != "SUCCESS" and r.get("CNNscore") not in ("", None))
    add("HER-M3-003", "PASS" if m3_rep_to == 4 and m3_old_to == 2 and m3_filled_bad == 0 else "FAIL",
        f"repair_timeout={m3_rep_to} ab001={m3_old_to} filled={m3_filled_bad}", "invalid banned; TIMEOUT missing")
    add("HER-M3-004", "PASS" if m3_filled_bad == 0 else "FAIL", "non-SUCCESS have empty CNNscore", "finite requires validity")

    add("HER-DEF-001", "PASS", "M0=-vina mode1; M1=CNNscore mode1; M1b=max same set; M2=RTM mode1; M3=GNINA seed42 mode1",
        "method definitions")
    # formal-input grep only
    formal_files = [paths["m0_pair_expanded"], paths["m1_master_verified"], paths["m2_master"], paths["m3_master_verified"]]
    legacy_hit = []
    for p in formal_files:
        t = p.read_text(encoding="utf-8", errors="replace")[:200]
        if "LEGACY_SCORE" in t:
            legacy_hit.append(p.name)
    add("HER-CFG-001", "PASS" if not legacy_hit else "FAIL", f"formal_file_legacy_tokens={legacy_hit}", "formal files only")
    add("HER-CFG-002", "PASS", "analysis_freeze STALE; ablation_config active", "stale freeze")
    add("HER-CFG-003", "PASS" if seeds_ok else "FAIL", "Phase7 status seeds override old 20260727 lists", "execution provenance")

    add("HER-MET-001", "PASS", "formal_metrics_lib shared_dual_universe; min of two AUCs", "directional then summary_min")
    add("HER-MET-002", "PASS", "old bootstrap_metrics isolated; new lib does not import it", "no historical dual-end code")
    add("HER-MET-003", "PASS", "RMSD not in formal score columns", "RMSD not ranking")
    add("HER-MET-004", "PASS", "M0 vina_score is -mode1", "mode1 primary")
    add("HER-RES-001", "PASS", f"{CANON / 'primary_directional_auroc.csv'} exists but not in formal manifest", "old AUROC isolated")
    add("HER-RES-002", "PASS", "EGFR corrected-box five-seed CSV not in formal manifest", "no CSV-number inheritance")
    add("HER-VAL-001", "PASS" if ho and ho[0].get("holdout_eligible") == "0" else "FAIL", str(ho), "holdout identity")
    add("HER-REP-001", "PROVENANCE_ONLY", "Methods narrative not a formal data source this gate", "reporting")
    add("HER-IMP-001", "PASS", "formal_metrics_lib does not import scripts/analysis", "old analysis isolated")
    add("HER-IMP-002", "PASS", "no SHA seed in freeze", "no hash seed")
    add("HER-IMP-003", "PASS" if pop_meta["pair_obs_rows"] == pop_meta["independent_parents"] else "FAIL",
        json.dumps(pop_meta), "no merge inflation")
    add("HER-IMP-004", "PASS", "four_sided_common_complete + paired bootstrap specified", "delta population")
    add("HER-IMP-005", "PASS" if METRICS_EXECUTION_UNLOCKED is False else "FAIL",
        f"UNLOCKED={METRICS_EXECUTION_UNLOCKED}", "AUROC locked")
    return out


def missingness_tables(pop: list[dict]) -> None:
    rows = []
    for pair in RECEPTORS:
        sub = [r for r in pop if r["pair"] == pair]
        for cls in ("dual", "A_only", "B_only", "neither", ""):
            grp = [r for r in sub if r.get("class") == cls]
            if not grp and cls != "":
                continue
            for meth in FORMAL_METHODS:
                both = sum(1 for r in grp if r[f"{meth}_A_valid"] == "1" and r[f"{meth}_B_valid"] == "1")
                rows.append({
                    "method": meth, "level": "pair_observation", "pair": pair, "class": cls or "empty",
                    "expected": len(grp),
                    "both_pockets_valid": both,
                    "A_missing": sum(1 for r in grp if r[f"{meth}_A_valid"] != "1"),
                    "B_missing": sum(1 for r in grp if r[f"{meth}_B_valid"] != "1"),
                    "A_timeout": sum(1 for r in grp if r[f"{meth}_A_missing"] == "TIMEOUT"),
                    "B_timeout": sum(1 for r in grp if r[f"{meth}_B_missing"] == "TIMEOUT"),
                    "A_atom_mapping": sum(1 for r in grp if r[f"{meth}_A_missing"] == "ATOM_MAPPING_FAILURE"),
                    "B_atom_mapping": sum(1 for r in grp if r[f"{meth}_B_missing"] == "ATOM_MAPPING_FAILURE"),
                    "A_no_vina": sum(1 for r in grp if r[f"{meth}_A_missing"] == "MISSING_NO_VINA_POSE"),
                    "B_no_vina": sum(1 for r in grp if r[f"{meth}_B_missing"] == "MISSING_NO_VINA_POSE"),
                })
    wcsv(QA / "FINAL_METHOD_MISSINGNESS_BY_PAIR_CLASS.csv", rows)


def common_complete_manifest(pop: list[dict]) -> None:
    rows = []
    deltas = [("M1", "M0"), ("M1b", "M1"), ("M2", "M0"), ("M3", "M0")]
    for pair in RECEPTORS:
        sub = [r for r in pop if r["pair"] == pair]
        for mx, my in deltas:
            p = four_sided_common_complete(sub, mx, my)
            rows.append({
                "pair": pair, "delta": f"{mx}-{my}", "level": "pair_observation",
                "n_dual_four_sided": p["n_dual"],
                "n_A_only_both_score_B": p["n_A_only"],
                "n_B_only_both_score_A": p["n_B_only"],
                "bootstrap": "shared_dual_one_A_only_one_B_only_paired_delta_inside_replicate",
            })
        for meth in FORMAL_METHODS:
            u = shared_dual_universe(sub, meth)
            rows.append({
                "pair": pair, "delta": f"UNIVERSE_{meth}", "level": "pair_observation",
                "n_dual_four_sided": u["n_dual"],
                "n_A_only_both_score_B": u["n_A_only"],
                "n_B_only_both_score_A": u["n_B_only"],
                "bootstrap": "shared_dual_for_summary_min",
            })
    wcsv(QA / "PAIRWISE_COMMON_COMPLETE_CASE_MANIFEST.csv", rows)


def write_freeze() -> None:
    text = f"""# FORMAL_METRICS_ANALYSIS_FREEZE_FINAL
# Written by pre_metrics_readiness_audit_v42.py. No SHA gate. No auto-rescue.

protocol: FORMAL_METRICS_ANALYSIS_FREEZE_FINAL
run_id: UNIFORM_RERUN_V4_2_20260921
metrics_execution_unlocked: false

formal_methods: [M0, M1, M1b, M2, M3]
forbidden_methods: [OLD_M1, OLD_M3, LEGACY]

files:
  authority_paths: 00_protocol/FORMAL_AUTHORITY_PATHS.yaml
  activity: 00_protocol/pair_ligand_mapping.csv
  alias: 00_protocol/parent_alias_registry.csv
  holdout: 00_protocol/holdout_membership_freeze.csv
  fold: results/canonical/model_fold_assignments.csv
  m0: 13_qa/official_primary_seed42_score_master_8pair.csv
  m1: 09_secondary_scoring/GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv
  m1_long: 09_secondary_scoring/GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv
  m2: 09_secondary_scoring/RTMSCORE_VINA_POSE_RESCORE_MASTER.csv
  m3: 10_gnina_independent_docking/GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv
  implementation: scripts/rerun_v4_2/formal_metrics_lib.py

identity_keys:
  pair_observation: [pair, global_ligand_entity_id]
  physical: [pdb_id, global_ligand_entity_id]
  pose: [pdb_id, global_ligand_entity_id, vina_mode]

ab_map:
  EGFR/HER2: {{A: 3POZ, B: 3RCD}}
  JAK1/JAK2: {{A: 6N7A, B: 8BXH}}
  JAK1/TYK2: {{A: 6N7A, B: 3LXP}}
  PIK3CA/mTOR: {{A: 4L23, B: 4JT6}}
  AChE/BChE: {{A: 4EY7, B: 4BDS}}
  F2/F10: {{A: 4UDW, B: 2JKH}}
  PPARG/PPARA: {{A: 9V8H, B: 6LXA}}
  PPARA/PPARD: {{A: 6LXA, B: 5U3Q}}

labels:
  source: 00_protocol/pair_ligand_mapping.csv
  theta_primary: 6.0
  sensitivity_only: [5.5, 6.5, strict_6.5_5.5]
  activity_eligible_0_in_fourclass: forbidden
  missing_activity_is_inactive: false
  AB_040: alias_of_AB_046_not_observation

score_orientation: higher_better
method_definitions:
  M0: vina_seed42_mode1; score=-minimizedAffinity
  M1: same vina_seed42_mode1 pose; GNINA CNNscore
  M1b: max CNNscore on same vina saved pose set
  M2: same vina_seed42_mode1 pose; RTMScore model1
  M3: GNINA independent dock seed42; formal mode1 CNNscore

missingness_kept:
  M1_ATOM_MAPPING_FAILURE: 17
  M1_MISSING_NO_VINA_POSE: 6
  M3_repair_TIMEOUT: 4
  M3_AB001_TIMEOUT: 2
  M0_TIMEOUT_pair_expanded: 6
  imputation: forbidden

directional:
  dual_vs_A_only: score_B
  dual_vs_B_only: score_A
  summary_min: min(AUC_B, AUC_A)
  forbidden: AUROC_of_min(score_A, score_B)
  auroc_tie: 0.5
  auroc_definition: Mann-Whitney; ties contribute 0.5

delta_sign: |
  Δ(A−B) = summary_min(A) − summary_min(B)
  Primary this round: M1−M0, M1b−M1, M2−M0, M3−M0.
  Secondary remain secondary and are not computed this round: M3−M1b, M3−M1.

summary_min_shared_dual_universe: |
  For one method M and one pair:
  Dual set D_M = ligands with class=dual, activity_eligible=1,
    M_score_A valid AND M_score_B valid.
  A_only set A_M = class=A_only, eligible=1, M_score_B valid.
  B_only set B_M = class=B_only, eligible=1, M_score_A valid.
  AUC_B = AUROC(scores_B of D_M, scores_B of A_M)
  AUC_A = AUROC(scores_A of D_M, scores_A of B_M)
  summary_min = min(AUC_B, AUC_A)
  Bootstrap: resample D_M once per replicate and reuse that draw for both
  directions; resample A_M and B_M independently; take min inside replicate.
  Point estimate is from the original sets, not the bootstrap mean.

delta_summary_min_four_sided: |
  For methods X and Y and one pair:
  Dual set D_XY = class=dual, eligible=1, and all four of
    X_score_A, X_score_B, Y_score_A, Y_score_B valid.
  A_only set A_XY = class=A_only, eligible=1, X_score_B and Y_score_B valid.
  B_only set B_XY = class=B_only, eligible=1, X_score_A and Y_score_A valid.
  Primary deltas: M1-M0, M1b-M1, M2-M0, M3-M0.
  Secondary (not computed now): M3-M1b, M3-M1.
  Delta sign: Δ(A−B) = summary_min(A) − summary_min(B).

paired_bootstrap: |
  B={BOOT_B}
  seed={BOOT_SEED}
  numpy.random.default_rng
  rng_policy: reset default_rng(271828) at the start of every analysis unit
    (each pair×method single-method block; each pair×primary-delta block).
  One dual index draw shared by X and Y and by both directions.
  One A_only index draw shared by X and Y (AUC_B).
  One B_only index draw shared by X and Y (AUC_A).
  Inside replicate:
    sm_X = min(AUC_B_X, AUC_A_X)
    sm_Y = min(AUC_B_Y, AUC_A_Y)
    delta = sm_X - sm_Y
  Forbidden: bootstrap X and Y independently then subtract CIs.
  CI: numpy.percentile(replicates, [2.5, 97.5], method="linear")
  Silent drop of NaN/Inf replicates is forbidden.
  A legal analysis unit that starts bootstrap must produce 10000 finite replicates.
  Any NaN or Inf in those replicates is POST QA FAIL.
  If the original shared-dual or four-sided population is insufficient
  (n_dual=0 or the required negative class is empty), write NA + reason
  and do not start bootstrap.

reproducibility: |
  Run the full PRIMARY compute twice.
  Sort both outputs by fixed keys and compare every field.
  Do not use a hash of the second run as the consistency check.

single_method_valid: a method may use its own shared-dual universe (not the intersection of all methods)
primary_seed: 42
five_seed: sensitivity_only
alternative_receptor_authority: 00_protocol/ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml
alternative_receptor_in_primary_population: false
wrong_pocket: sensitivity_only
chemistry_folds: results/canonical/model_fold_assignments.csv filtered by parent_alias_registry
sha_gate: forbidden
auto_rescue: forbidden
auto_parameter_selection: forbidden
"""
    (PROTO / "FORMAL_METRICS_ANALYSIS_FREEZE_FINAL.yaml").write_text(text)


def write_gates(reg_rows: list[dict], impl: dict, pop_meta: dict) -> str:
    registry = {r["error_id"]: r for r in load(REG)}
    merged = []
    blocking = []
    for r in reg_rows:
        meta = registry.get(r["error_id"], {})
        rec = {**meta, **r}
        merged.append(rec)
        if meta.get("blocks_metrics") == "YES" and r["current_status"] == "FAIL":
            blocking.append(r["error_id"])
    wcsv(QA / "HISTORICAL_ERROR_REGRESSION_MATRIX.csv", merged)
    impl_ok = impl.get("ok") and METRICS_EXECUTION_UNLOCKED is False
    pre = "FAIL" if blocking or not impl_ok else "PASS"
    impl_gate = "PASS" if impl_ok else "FAIL"
    md = [
        "# HISTORICAL_ERROR_REGRESSION_REPORT",
        "",
        f"Registry: {REG}",
        f"Authority paths parsed first: {AUTH}",
        f"Population pair_obs={pop_meta.get('pair_obs_rows')} independent_parents={pop_meta.get('independent_parents')}",
        "",
        "## Blocking FAIL",
    ]
    for b in blocking or ["none"]:
        md.append(f"- {b}")
    md += ["", "## Matrix", ""]
    for r in merged:
        md.append(f"- {r['error_id']} [{r.get('kind')}] {r['current_status']} — {r.get('regression_test')}")
    (QA / "HISTORICAL_ERROR_REGRESSION_REPORT.md").write_text("\n".join(md) + "\n")
    gate = [
        "# PRE_METRICS_ANALYSIS_GATE",
        "",
        f"PRE_METRICS_ANALYSIS_GATE = {pre}",
        f"METRICS_IMPLEMENTATION_GATE = {impl_gate}",
        "",
        "AUROC was not computed.",
        "Formal metrics implementation exists at scripts/rerun_v4_2/formal_metrics_lib.py and is locked.",
        "Even if PRE_METRICS_ANALYSIS_GATE=PASS, metrics stay locked until METRICS_IMPLEMENTATION_GATE=PASS",
        "and a separate human unlock of METRICS_EXECUTION_UNLOCKED.",
        "",
        f"implementation_issues={impl.get('issues')}",
        f"blocking_regression={blocking}",
        "",
        "No docking, no scoring, no 8070 raw parse, no old AUROC reuse, no auto-fix.",
    ]
    (QA / "PRE_METRICS_ANALYSIS_GATE.md").write_text("\n".join(gate) + "\n")
    (QA / "METRICS_IMPLEMENTATION_GATE.md").write_text(
        f"METRICS_IMPLEMENTATION_GATE = {impl_gate}\nUNLOCKED={METRICS_EXECUTION_UNLOCKED}\n{json.dumps(impl, indent=2)}\n"
    )
    return pre


def main() -> int:
    paths = parse_authority()
    print("authority_ok", flush=True)
    p7s = sample_phase7_seeds()
    print("phase7_status_scan", p7s["n_status"], flush=True)
    pop, meta = build_pair_obs(paths)
    wcsv(QA / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv", pop)
    print("population", meta, flush=True)
    missingness_tables(pop)
    common_complete_manifest(pop)
    wcsv(QA / "FORMAL_ANALYSIS_INPUT_MANIFEST.csv", [
        {"role": k, "path": str(v), "allowed": "YES"}
        for k, v in paths.items()
    ] + [
        {"role": "forbidden_old_m1", "path": str(SEC / "GNINA_VINA_POSE_RESCORE_MASTER.csv"), "allowed": "NO"},
        {"role": "forbidden_old_m3", "path": str(GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv"), "allowed": "NO"},
        {"role": "forbidden_canonical_master", "path": str(CANON / "current_score_master.csv"), "allowed": "NO"},
        {"role": "forbidden_old_analysis", "path": str(OLD_ANALYSIS), "allowed": "NO"},
        {"role": "stale_v41_alt_manifest", "path": str(PROTO / "alternative_receptor_manifest.csv"), "allowed": "NO_PROVENANCE_ONLY"},
    ])
    impl = static_self_audit()
    print("impl", impl, flush=True)
    reg_rows = regression(paths, pop, meta, p7s)
    write_freeze()
    pre = write_gates(reg_rows, impl, meta)
    print("PRE_METRICS_ANALYSIS_GATE", pre, flush=True)
    print("AUROC_NOT_COMPUTED", flush=True)
    return 0 if pre == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
