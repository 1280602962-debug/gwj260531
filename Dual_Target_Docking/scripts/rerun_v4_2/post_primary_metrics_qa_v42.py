#!/usr/bin/env python3
"""POST_PRIMARY_METRICS_QA. Binary PASS/FAIL. Does not interpret AUROC values."""
from __future__ import annotations

import ast
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path("/tmp/pr39_fiveseed/Dual_Target_Docking")
RUN = ROOT / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
QA = RUN / "13_qa"
PROTO = RUN / "00_protocol"
OUT = ROOT / "results" / "formal_metrics"
SCRIPTS = ROOT / "scripts" / "rerun_v4_2"
SEC = RUN / "09_secondary_scoring"
M3P = RUN / "10_gnina_independent_docking"

sys.path.insert(0, str(SCRIPTS))
import formal_metrics_lib as F  # noqa: E402
from run_primary_metrics_v42 import (  # noqa: E402
    DELTA_FIELDS, DELTA_KEYS, DIR_FIELDS, DIR_KEYS, MEM_FIELDS,
    compute, load_csv, rows_equal, sort_rows,
)

ISSUES: list[str] = []

M1_ATOM = [
    ("4BDS", "CFFDQXIMLVONET-UHFFFAOYSA-N"),
    ("4EY7", "CFFDQXIMLVONET-UHFFFAOYSA-N"),
    ("4UDW", "SCZBXBBBMCBNSK-TVGQLCNQSA-N"),
    ("5U3Q", "GUFHDBFHNDUDJA-GXJYFLSQSA-N"),
    ("5U3Q", "KGCCMHIRXXLQQG-UHFFFAOYSA-N"),
    ("5U3Q", "XWDWDQGJOOSFCU-CFVCEQMKSA-N"),
    ("5U3Q", "XWDWDQGJOOSFCU-UHFFFAOYSA-N"),
    ("5U3Q", "ZYTAJJSXWKDKGN-UHFFFAOYSA-N"),
    ("6LXA", "CADIUVOXJVKRSR-UHFFFAOYSA-N"),
    ("6LXA", "GUFHDBFHNDUDJA-GXJYFLSQSA-N"),
    ("6LXA", "KGCCMHIRXXLQQG-UHFFFAOYSA-N"),
    ("6LXA", "KSIQZOUSAASEIJ-UHFFFAOYSA-N"),
    ("6LXA", "XWDWDQGJOOSFCU-CFVCEQMKSA-N"),
    ("6LXA", "XWDWDQGJOOSFCU-UHFFFAOYSA-N"),
    ("6LXA", "ZYTAJJSXWKDKGN-UHFFFAOYSA-N"),
    ("9V8H", "CADIUVOXJVKRSR-UHFFFAOYSA-N"),
    ("9V8H", "KSIQZOUSAASEIJ-UHFFFAOYSA-N"),
]
M1_NOVINA = [
    ("4BDS", "HDAQKLAKMZNBQB-UHFFFAOYSA-N"),
    ("4BDS", "KHBDHGUQOPEDSD-UHFFFAOYSA-N"),
    ("4BDS", "NPHIQFCGKIKGSA-YQOHNZFASA-N"),
    ("4EY7", "HDAQKLAKMZNBQB-UHFFFAOYSA-N"),
    ("4EY7", "KHBDHGUQOPEDSD-UHFFFAOYSA-N"),
    ("4EY7", "NPHIQFCGKIKGSA-YQOHNZFASA-N"),
]
M3_TIMEOUT = [
    ("2JKH", "DWNHMQOMYAIZCG-HHHXNRCGSA-N"),
    ("2JKH", "FNERQDASDNXFMT-IRXDYDNUSA-N"),
    ("4BDS", "NPHIQFCGKIKGSA-YQOHNZFASA-N"),
    ("4EY7", "NPHIQFCGKIKGSA-YQOHNZFASA-N"),
    ("4UDW", "DWNHMQOMYAIZCG-HHHXNRCGSA-N"),
    ("4UDW", "FNERQDASDNXFMT-IRXDYDNUSA-N"),
]
M0_TIMEOUT_LIGANDS = {"AB_001", "AB_053", "AB_054"}
ELIGIBLE0 = {("EGFR/HER2", "EH120_059"), ("AChE/BChE", "AB_087")}
ALT_PRIMARY_FORBIDDEN = {
    "4L2Y", "4JT5", "4EY6", "1P0M", "3SHC", "2Y5F", "6KAX", "5U46",
    "4JPS", "5DXT", "4JSX",
}
FORBIDDEN_NAME = {
    "GNINA_VINA_POSE_RESCORE_MASTER.csv",
    "GNINA_VINA_POSE_RESCORE_POSE_LONG.csv",
    "GNINA_SEED42_PRODUCTION_MASTER.csv",
    "current_score_master.csv",
}


def fail(msg: str) -> None:
    ISSUES.append(msg)


def load(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def pop_side(row: dict, pdb: str) -> str | None:
    if row.get("pdb_A") == pdb:
        return "A"
    if row.get("pdb_B") == pdb:
        return "B"
    return None


def check_missing_physical(pop, jobs, meth, label):
    seen = 0
    for pdb, eid in jobs:
        hits = [r for r in pop if r["global_ligand_entity_id"] == eid and pop_side(r, pdb)]
        if not hits:
            fail(f"{label}_not_in_population {pdb} {eid}")
            continue
        for r in hits:
            seen += 1
            side = pop_side(r, pdb)
            if r.get(f"{meth}_{side}_valid") != "0" or r.get(f"{meth}_score_{side}", "") != "":
                fail(f"{label}_not_missing {r['pair']} {pdb} {eid} {meth}")
    if label.startswith("M1_ATOM") and seen < 17:
        fail(f"{label}_hit_rows={seen}<17")
    if label.startswith("M1_NOVINA") and seen < 6:
        fail(f"{label}_hit_rows={seen}<6")
    if label.startswith("M3") and seen < 6:
        fail(f"{label}_hit_rows={seen}<6")


def ast_no_old_analysis(*paths: Path) -> None:
    for p in paths:
        tree = ast.parse(p.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name == "analysis" or a.name.startswith("analysis."):
                        fail(f"import_old_analysis {p.name}")
            if isinstance(node, ast.ImportFrom) and node.module:
                if node.module == "analysis" or node.module.startswith("analysis."):
                    fail(f"import_old_analysis {p.name}")
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name == "sklearn" or a.name.startswith("sklearn."):
                        fail(f"import_sklearn {p.name}")
            if isinstance(node, ast.ImportFrom) and node.module:
                if node.module == "sklearn" or node.module.startswith("sklearn."):
                    fail(f"import_sklearn {p.name}")


def main() -> int:
    dir_path = OUT / "PRIMARY_DIRECTIONAL_METRICS.csv"
    delta_path = OUT / "PRIMARY_METHOD_DELTA_METRICS.csv"
    mem_path = QA / "PRIMARY_METRICS_UNIVERSE_MEMBERSHIP.csv"
    man_path = QA / "PRIMARY_METRICS_EXECUTION_MANIFEST.csv"
    pop_path = QA / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv"
    cc_path = QA / "PAIRWISE_COMMON_COMPLETE_CASE_MANIFEST.csv"

    for p in (dir_path, delta_path, mem_path, man_path, pop_path, cc_path):
        if not p.is_file():
            fail(f"missing_file {p}")
    if ISSUES:
        return write_gate()

    directional = load(dir_path)
    deltas = load(delta_path)
    membership = load(mem_path)
    manifest = {r["key"]: r["value"] for r in load(man_path)}
    pop = load(pop_path)
    cc = load(cc_path)
    m1 = load(SEC / "GNINA_VINA_POSE_RESCORE_MASTER_VERIFIED.csv")
    m1_long = load(SEC / "GNINA_VINA_POSE_RESCORE_POSE_LONG_VERIFIED.csv")
    m2 = load(SEC / "RTMSCORE_VINA_POSE_RESCORE_MASTER.csv")
    m3 = load(M3P / "GNINA_SEED42_PRODUCTION_MASTER_TOPOLOGY_VERIFIED.csv")
    m0 = load(QA / "official_primary_seed42_score_master_8pair.csv")

    if len(pop) != 807:
        fail(f"population_n={len(pop)}")
    if any(r.get("canonical_ligand_id") == "AB_040" or r.get("panel_id") == "AB_040" for r in pop):
        fail("AB_040_observation")
    if any(r.get("canonical_ligand_id") == "AB_040" for r in membership):
        fail("AB_040_in_membership")

    elig0 = [(r["pair"], r["canonical_ligand_id"]) for r in pop if r.get("activity_eligible") != "1"]
    if set(elig0) != ELIGIBLE0:
        fail(f"eligible0_ids={elig0}")
    if any((r["pair"], r["canonical_ligand_id"]) in ELIGIBLE0 for r in membership):
        fail("eligible0_entered_metrics")
    if any(r.get("class") == "neither" for r in membership):
        fail("neither_entered_metrics")
    n_neither = sum(1 for r in pop if r.get("class") == "neither" and r.get("activity_eligible") == "1")
    if n_neither != 101:
        fail(f"neither_count={n_neither}")

    pdbs = {r["pdb_A"] for r in pop} | {r["pdb_B"] for r in pop}
    leak = sorted(pdbs & ALT_PRIMARY_FORBIDDEN)
    if leak:
        fail(f"alt_or_withdrawn_pdb={leak}")

    if len(directional) != 8 * 5 * 3:
        fail(f"directional_nrow={len(directional)}")
    if len(deltas) != 8 * 4:
        fail(f"delta_nrow={len(deltas)}")
    if any(r.get("delta") in {"M3-M1b", "M3-M1"} for r in deltas):
        fail("secondary_delta_computed")
    if {r["delta"] for r in deltas} != {"M1-M0", "M1b-M1", "M2-M0", "M3-M0"}:
        fail(f"delta_set={ {r['delta'] for r in deltas} }")
    if any(r.get("delta_class") != "primary" for r in deltas):
        fail("delta_class_not_primary")

    exp_cc = {(r["pair"], r["delta"]): r for r in cc}
    by_pair = defaultdict(list)
    for r in pop:
        by_pair[r["pair"]].append(r)

    mem_sets = defaultdict(set)
    for r in membership:
        mem_sets[(r["pair"], r["universe_type"], r["method_or_delta"], r["class"])].add(
            r["global_ligand_entity_id"]
        )

    for pair in F.RECEPTORS:
        sub = by_pair[pair]
        n_duals = {}
        for meth in F.FORMAL_METHODS:
            u = F.shared_dual_universe(sub, meth)
            exp = exp_cc[(pair, f"UNIVERSE_{meth}")]
            got = (str(u["n_dual"]), str(u["n_A_only"]), str(u["n_B_only"]))
            want = (exp["n_dual_four_sided"], exp["n_A_only_both_score_B"], exp["n_B_only_both_score_A"])
            if got != want:
                fail(f"universe_n {pair} {meth} {got}!={want}")
            for cls, key in (("dual", "dual_shared"), ("A_only", "A_only_for_AUC_B"), ("B_only", "B_only_for_AUC_A")):
                ids = {x["global_ligand_entity_id"] for x in u[key]}
                if ids != mem_sets[(pair, "UNIVERSE", meth, cls)]:
                    fail(f"membership_set {pair} UNIVERSE {meth} {cls}")
            n_duals[meth] = u["n_dual"]
        drows = [r for r in directional if r["pair"] == pair]
        for meth in F.FORMAL_METHODS:
            sides = [r for r in drows if r["method"] == meth]
            nd = {r["contrast"]: r["n_dual"] for r in sides}
            if nd.get("dual_vs_A_only") != nd.get("dual_vs_B_only") or nd.get("dual_vs_A_only") != nd.get("summary_min"):
                fail(f"n_dual_not_shared {pair} {meth} {nd}")
            if str(n_duals[meth]) != str(nd.get("dual_vs_A_only")):
                fail(f"n_dual_ne_universe {pair} {meth}")
        for mx, my in F.PRIMARY_DELTAS:
            u = F.four_sided_common_complete(sub, mx, my)
            exp = exp_cc[(pair, f"{mx}-{my}")]
            got = (str(u["n_dual"]), str(u["n_A_only"]), str(u["n_B_only"]))
            want = (exp["n_dual_four_sided"], exp["n_A_only_both_score_B"], exp["n_B_only_both_score_A"])
            if got != want:
                fail(f"delta_n {pair} {mx}-{my} {got}!={want}")
            drow = next(r for r in deltas if r["pair"] == pair and r["delta"] == f"{mx}-{my}")
            if (str(drow["n_dual_four_sided"]), str(drow["n_A_only_both_score_B"]), str(drow["n_B_only_both_score_A"])) != want:
                fail(f"delta_row_n {pair} {mx}-{my}")
            for cls, key in (("dual", "dual_four_sided"), ("A_only", "A_only_both_score_B"), ("B_only", "B_only_both_score_A")):
                ids = {x["global_ligand_entity_id"] for x in u[key]}
                if ids != mem_sets[(pair, "DELTA", f"{mx}-{my}", cls)]:
                    fail(f"membership_set {pair} DELTA {mx}-{my} {cls}")

    for r in directional:
        if r["bootstrap_started"] == "YES":
            if r["n_replicates"] != "10000" or r["n_finite"] != "10000":
                fail(f"boot_count {r['pair']} {r['method']} {r['contrast']} {r['n_replicates']}/{r['n_finite']}")
            if r["auroc"] == "" or r["reason"] != "":
                fail(f"boot_without_value {r['pair']} {r['method']} {r['contrast']}")
            try:
                lo, hi, val = float(r["ci_lo"]), float(r["ci_hi"]), float(r["auroc"])
                if not all(map(np_finite := (lambda x: x == x and abs(x) != float("inf")), (lo, hi, val))):
                    fail(f"nonfinite_metric {r['pair']} {r['method']} {r['contrast']}")
            except ValueError:
                fail(f"ci_parse {r['pair']} {r['method']} {r['contrast']}")
        else:
            if r["n_replicates"] != "0" or r["ci_lo"] != "" or r["auroc"] != "" or not r["reason"]:
                fail(f"na_shape {r['pair']} {r['method']} {r['contrast']}")
        if r["boot_B"] != "10000" or r["boot_seed"] != "271828":
            fail("boot_constants_directional")
        if r["ci_method"] != "linear" or r["rng_policy"] != "reset_default_rng_per_analysis_unit":
            fail("ci_or_rng_policy_directional")
        if r["contrast"] == "dual_vs_A_only" and r["score_used"] != "score_B":
            fail("wrong_score_dual_vs_A_only")
        if r["contrast"] == "dual_vs_B_only" and r["score_used"] != "score_A":
            fail("wrong_score_dual_vs_B_only")
        keys_l = {k.lower() for k in r}
        if keys_l & {"p_value", "pvalue", "auprc", "ef", "top_k", "topk"}:
            fail("forbidden_metric_column")

    for r in directional:
        if r["contrast"] != "summary_min":
            continue
        b = next(x for x in directional if x["pair"] == r["pair"] and x["method"] == r["method"] and x["contrast"] == "dual_vs_A_only")
        a = next(x for x in directional if x["pair"] == r["pair"] and x["method"] == r["method"] and x["contrast"] == "dual_vs_B_only")
        if r["bootstrap_started"] == "YES":
            if b["auroc"] == "" or a["auroc"] == "":
                fail(f"summary_min_without_directions {r['pair']} {r['method']}")
            else:
                expect = min(float(b["auroc"]), float(a["auroc"]))
                if abs(float(r["auroc"]) - expect) > 1e-12:
                    fail(f"summary_min_ne_min_directions {r['pair']} {r['method']}")
        else:
            if b["bootstrap_started"] == "YES" and a["bootstrap_started"] == "YES":
                fail(f"summary_min_na_but_directions_ok {r['pair']} {r['method']}")

    for r in deltas:
        if r["bootstrap_started"] == "YES":
            if r["n_replicates"] != "10000" or r["n_finite"] != "10000":
                fail(f"delta_boot {r['pair']} {r['delta']}")
            try:
                dlt = float(r["delta_summary_min"])
                sx, sy = float(r["summary_min_x"]), float(r["summary_min_y"])
                if abs(dlt - (sx - sy)) > 5e-10:
                    fail(f"delta_sign {r['pair']} {r['delta']}")
            except ValueError:
                fail(f"delta_parse {r['pair']} {r['delta']}")
        elif r["reason"] == "" or r["delta_summary_min"] != "":
            fail(f"delta_na_shape {r['pair']} {r['delta']}")
        if r["boot_B"] != "10000" or r["boot_seed"] != "271828" or r["paired"] != "YES":
            fail("delta_constants")
        if r["delta_sign"] != "summary_min(A)-summary_min(B)":
            fail("delta_sign_label")

    check_missing_physical(pop, M1_ATOM, "M1", "M1_ATOM")
    check_missing_physical(pop, M1_ATOM, "M1b", "M1b_ATOM")
    check_missing_physical(pop, M1_NOVINA, "M1", "M1_NOVINA")
    check_missing_physical(pop, M3_TIMEOUT, "M3", "M3_TIMEOUT")
    for r in pop:
        if r["pair"] == "AChE/BChE" and r.get("canonical_ligand_id") in M0_TIMEOUT_LIGANDS:
            for side in ("A", "B"):
                if r.get(f"M0_{side}_valid") != "0" or r.get(f"M0_score_{side}", "") != "":
                    fail(f"M0_TIMEOUT_filled {r['canonical_ligand_id']} {side}")

    m1_tax = {}
    for r in m1:
        m1_tax[r["job_status"]] = m1_tax.get(r["job_status"], 0) + 1
    if m1_tax.get("ATOM_MAPPING_FAILURE") != 17 or m1_tax.get("MISSING_NO_VINA_POSE") != 6:
        fail(f"m1_master_tax={m1_tax}")
    if any(r["job_status"] != "SUCCESS" and r.get("M1_CNNscore") not in ("", None) for r in m1):
        fail("m1_nonsuccess_scored")
    if any(r["status"] != "SUCCESS" and r.get("CNNscore") not in ("", None) for r in m3):
        fail("m3_nonsuccess_scored")
    if sum(1 for r in m3 if r["status"] == "TIMEOUT") != 6:
        fail("m3_timeout_count")
    if {r["seed"] for r in m3} != {"42"}:
        fail("m3_seed_not_42")
    if any(r["status"] == "TIMEOUT" and r.get("official_m3") == "YES" for r in m3):
        fail("m3_timeout_official")

    sign_bad = 0
    for r in m0:
        if r["vina_status"] == "SUCCESS":
            if abs(float(r["vina_score"]) + float(r["vina_mode1_affinity"])) > 1e-6:
                sign_bad += 1
        elif r["vina_status"] == "TIMEOUT" and r.get("vina_score") not in ("", None):
            fail("m0_timeout_filled")
    if sign_bad:
        fail(f"m0_sign_bad={sign_bad}")

    m1b_lt = 0
    for r in pop:
        for side in ("A", "B"):
            if r.get(f"M1_{side}_valid") == "1" and r.get(f"M1b_{side}_valid") == "1":
                if float(r[f"M1b_score_{side}"]) + 1e-12 < float(r[f"M1_score_{side}"]):
                    m1b_lt += 1
    if m1b_lt:
        fail(f"M1b_lt_M1={m1b_lt}")

    long_max = defaultdict(lambda: -1.0)
    for r in m1_long:
        try:
            s = float(r["CNNscore"])
        except (TypeError, ValueError):
            continue
        k = (r["pdb_id"], r["global_ligand_entity_id"])
        if s > long_max[k]:
            long_max[k] = s
    m1b_mis = 0
    for r in m1:
        if r["job_status"] != "SUCCESS":
            continue
        k = (r["pdb_id"], r["global_ligand_entity_id"])
        if abs(float(r["M1b_CNNscore"]) - long_max[k]) > 1e-8:
            m1b_mis += 1
    if m1b_mis:
        fail(f"M1b_ne_long_max={m1b_mis}")

    m2_by = {(r["pdb_id"], r["global_ligand_entity_id"]): r for r in m2}
    m2_wrong = 0
    for r in pop:
        for side, pdb in (("A", r["pdb_A"]), ("B", r["pdb_B"])):
            if r.get(f"M2_{side}_valid") != "1":
                continue
            src = m2_by.get((pdb, r["global_ligand_entity_id"]))
            if not src or abs(float(r[f"M2_score_{side}"]) - float(src["M2_RTMScore"])) > 1e-6:
                m2_wrong += 1
    if m2_wrong:
        fail(f"M2_not_mode1={m2_wrong}")

    if manifest.get("two_run_compare") != "PASS":
        fail("runner_two_run_compare_not_PASS")
    if manifest.get("secondary_deltas_computed") != "NO":
        fail("manifest_secondary_computed")
    if manifest.get("source_METRICS_EXECUTION_UNLOCKED") != "False":
        fail("manifest_source_unlock")
    opened = [r["value"] for r in load(man_path) if r["key"] == "opened"]
    for p in opened:
        name = Path(p).name
        if name in FORBIDDEN_NAME:
            fail(f"opened_forbidden {p}")
        if "/scripts/analysis/" in p.replace("\\", "/"):
            fail(f"opened_old_analysis {p}")

    ast_no_old_analysis(
        SCRIPTS / "formal_metrics_lib.py",
        SCRIPTS / "run_primary_metrics_v42.py",
        SCRIPTS / "post_primary_metrics_qa_v42.py",
    )
    if F.METRICS_EXECUTION_UNLOCKED is not False:
        fail("qa_entered_with_unlock_true")
    src_assign = None
    tree = ast.parse((SCRIPTS / "formal_metrics_lib.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if "METRICS_EXECUTION_UNLOCKED" in names:
                src_assign = ast.literal_eval(node.value)
    if src_assign is not False:
        fail("source_unlock_not_false")

    keys = [("pair", r["pair"]) for r in pop]
    if len(keys) != len({(r["pair"], r["global_ligand_entity_id"]) for r in pop}):
        fail("population_key_not_unique")
    if len(membership) != len({(
        r["pair"], r["universe_type"], r["method_or_delta"], r["class"], r["global_ligand_entity_id"]
    ) for r in membership}):
        fail("membership_duplicated")

    # Independent recompute: same seed, field compare to written tables.
    F.METRICS_EXECUTION_UNLOCKED = True
    d3, g3, m3r = compute(pop)
    F.METRICS_EXECUTION_UNLOCKED = False
    c1 = rows_equal(directional, d3, DIR_KEYS, DIR_FIELDS)
    c2 = rows_equal(deltas, g3, DELTA_KEYS, DELTA_FIELDS)
    c3 = rows_equal(membership, m3r, ("pair", "universe_type", "method_or_delta", "class", "global_ligand_entity_id"), MEM_FIELDS)
    for msg in c1 + c2 + c3:
        fail("qa_recompute_mismatch:" + msg)

    return write_gate()


def write_gate() -> int:
    status = "FAIL" if ISSUES else "PASS"
    lines = [
        "# POST_PRIMARY_METRICS_QA",
        "",
        f"POST_PRIMARY_METRICS_QA = {status}",
        "",
        "No result interpretation.",
        "No sensitivity.",
        "No population or parameter change.",
        "",
    ]
    if ISSUES:
        lines.append("issues:")
        lines.extend(f"- {x}" for x in ISSUES)
    else:
        lines.append("issues: none")
    (QA / "POST_PRIMARY_METRICS_QA.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"POST_PRIMARY_METRICS_QA = {status}")
    if ISSUES:
        for x in ISSUES:
            print("ISSUE", x)
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
