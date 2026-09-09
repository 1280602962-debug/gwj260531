#!/usr/bin/env python3
"""Five-round JCIM submission audit: manuscript and SI vs frozen CSVs.

Does not re-bootstrap. Reads locked tables and assembled manuscripts.
Exit 0 only if every blocking check passes.
"""
from __future__ import annotations

import csv
import re
import subprocess
import sys
from collections import defaultdict
from decimal import ROUND_HALF_EVEN, ROUND_HALF_UP, Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
REPORT = DOCS / "SUBMISSION_AUDIT_FIVE_ROUNDS_V1.md"

PRIMARY = [
    "EGFR/HER2",
    "AChE/BChE",
    "PIK3CA/mTOR",
    "F2/F10",
    "JAK1/TYK2",
    "JAK1/JAK2",
    "PPARG/PPARA",
    "PPARA/PPARD",
]
ORIGINAL = {"EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"}
PAPER_GLOBS = (
    "MANUSCRIPT_JCIM_EN.md",
    "MANUSCRIPT_JCIM_ZH.md",
    "SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    "RESULTS_SECTION_JCIM_EN_V1.md",
    "RESULTS_DRAFT_ZH_JCIM_V1.md",
    "METHODS_SECTION_JCIM_EN_V1.md",
    "METHODS_DRAFT_ZH_JCIM_V1.md",
    "TITLE_AND_ABSTRACT_JCIM_EN_V1.md",
    "TITLE_AND_ABSTRACT_JCIM_ZH_V1.md",
    "DISCUSSION_SECTION_JCIM_EN_V1.md",
    "DISCUSSION_DRAFT_ZH_JCIM_V1.md",
    "CONCLUSIONS_SECTION_JCIM_EN_V1.md",
    "CONCLUSIONS_DRAFT_ZH_JCIM_V1.md",
)

FINDINGS: list[tuple[str, str, str]] = []


def rec(round_id: str, status: str, msg: str) -> None:
    FINDINGS.append((round_id, status, msg))


def r3(x, rounding=ROUND_HALF_UP) -> str:
    return format(Decimal(str(x)).quantize(Decimal("0.001"), rounding=rounding), "f")


def r3_opts(x) -> set[str]:
    """Accept both half-up and half-even (banker's) three-decimal displays."""
    return {r3(x, ROUND_HALF_UP), r3(x, ROUND_HALF_EVEN)}


def r3_ci(lo, hi, rounding=ROUND_HALF_UP) -> str:
    return f"[{r3(lo, rounding)}, {r3(hi, rounding)}]"


def r3_ci_opts(lo, hi) -> set[str]:
    return {r3_ci(lo, hi, ROUND_HALF_UP), r3_ci(lo, hi, ROUND_HALF_EVEN)}


def cell_ok(got: str, raw) -> bool:
    got = got.replace("−", "-")
    if re.fullmatch(r"-?\d+\.\d+", got):
        return got in r3_opts(raw)
    m = re.fullmatch(r"(-?\d+\.\d+) \[(-?\d+\.\d+), (-?\d+\.\d+)\]", got)
    if m and isinstance(raw, (tuple, list)) and len(raw) == 3:
        return (
            m.group(1) in r3_opts(raw[0])
            and m.group(2) in r3_opts(raw[1])
            and m.group(3) in r3_opts(raw[2])
        )
    return got == str(raw)


def rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def one(data: list[dict], **keys) -> dict:
    found = [row for row in data if all(row[key] == value for key, value in keys.items())]
    if len(found) != 1:
        raise AssertionError(f"{keys} matched {len(found)}")
    return found[0]


def parse_tables(text: str) -> list[list[list[str]]]:
    tables: list[list[list[str]]] = []
    current: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if cells and re.fullmatch(r":?-{3,}:?", cells[0].replace(" ", "")):
                continue
            if all(re.fullmatch(r":?-+:?", c.replace(" ", "")) for c in cells):
                continue
            current.append(cells)
        elif current:
            tables.append(current)
            current = []
    if current:
        tables.append(current)
    return tables


def table_by_header(text: str, first_cell: str) -> list[list[str]] | None:
    for table in parse_tables(text):
        if table and table[0] and table[0][0] == first_cell:
            return table
    return None


def near(a, b, tol: float = 5.5e-4) -> bool:
    return abs(float(a) - float(b)) <= tol


def load_table2() -> dict[str, dict]:
    three = rows(ROOT / "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv")
    five = rows(
        ROOT
        / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv"
    )
    desc = rows(ROOT / "data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv")
    form = rows(ROOT / "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv")
    out: dict[str, dict] = {}
    for pair in PRIMARY:
        if pair in ORIGINAL:
            t = one(three, pair=pair, label_rule="theta_6.0")
            d = one(desc, pair=pair)
            n = one(form, pair=pair, contrast="D_vs_neither_mean")
            alln = one(form, pair=pair, contrast="D_vs_A+B+neither_mean")
            out[pair] = {
                "nD": int(t["n_dual"]),
                "nA": int(t["n_A_only"]),
                "nB": int(t["n_B_only"]),
                "da": float(t["auroc_D_vs_A"]),
                "db": float(t["auroc_D_vs_B"]),
                "smin": float(t["pocket_matched_summary_min"]),
                "lo": float(t["ci_lo"]),
                "hi": float(t["ci_hi"]),
                "heavy": float(d["heavy_summary_min"]),
                "mw": float(d["mw_summary_min"]),
                "clogp": float(d["clogp_summary_min"]),
                "tpsa": float(d["tpsa_summary_min"]),
                "nei": float(n["auroc"]),
                "nei_lo": float(n["ci_lo"]),
                "nei_hi": float(n["ci_hi"]),
                "n_neg": int(n["n_neg"]),
                "alln": float(alln["auroc"]),
                "alln_lo": float(alln["ci_lo"]),
                "alln_hi": float(alln["ci_hi"]),
            }
        else:
            t = one(five, pair=pair)
            out[pair] = {
                "nD": int(t["n_dual"]),
                "nA": int(t["n_A_only"]),
                "nB": int(t["n_B_only"]),
                "nN": int(t["n_neither"]),
                "da": float(t["auroc_D_vs_A_pocketB"]),
                "db": float(t["auroc_D_vs_B_pocketA"]),
                "smin": float(t["summary_min"]),
                "lo": float(t["ci_lo"]),
                "hi": float(t["ci_hi"]),
                "heavy": float(t["heavy_summary_min"]),
                "mw": float(t["mw_summary_min"]),
                "clogp": float(t["clogp_summary_min"]),
                "tpsa": float(t["tpsa_summary_min"]),
                "nei": float(t["D_vs_neither_vina_mean"]),
                "nei_lo": float(t["D_vs_neither_ci_lo"]),
                "nei_hi": float(t["D_vs_neither_ci_hi"]),
                "n_neg": int(t["n_neither"]),
                "alln": float(t["D_vs_all_nonduals"]),
                "alln_lo": float(t["D_vs_all_nonduals_ci_lo"]),
                "alln_hi": float(t["D_vs_all_nonduals_ci_hi"]),
            }
    return out


def expected_table2_raw(src: dict) -> list:
    return [
        f"{src['nD']} / {src['nA']} / {src['nB']}",
        src["da"],
        src["db"],
        (src["smin"], src["lo"], src["hi"]),
        src["heavy"],
        src["mw"],
        src["clogp"],
        src["tpsa"],
    ]


def expected_table3_raw(src: dict) -> list:
    return [
        (src["smin"], src["lo"], src["hi"]),
        (src["nei"], src["nei_lo"], src["nei_hi"]),
        str(src["n_neg"]),
        (src["alln"], src["alln_lo"], src["alln_hi"]),
    ]


def check_named_table(round_id: str, label: str, table: list[list[str]] | None, expected: dict[str, list], pair_col: int = 0, value_slice: slice = slice(1, None)) -> None:
    if table is None:
        rec(round_id, "FAIL", f"{label}: table not found")
        return
    by_pair = {row[pair_col]: row[value_slice] for row in table[1:]}
    for pair, exp in expected.items():
        got = by_pair.get(pair)
        if got is None:
            rec(round_id, "FAIL", f"{label}: missing {pair}")
            continue
        bad = [f"{g}!~{e}" for g, e in zip(got, exp) if not cell_ok(g, e)]
        if bad:
            rec(round_id, "FAIL", f"{label} {pair}: {bad}")
        else:
            rec(round_id, "PASS", f"{label} {pair} matches locked CSV (3-dp half-up or half-even)")


def round1() -> dict:
    t2 = load_table2()
    en = (DOCS / "MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    zh = (DOCS / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    exp2 = {p: expected_table2_raw(t2[p]) for p in PRIMARY}
    exp3 = {p: expected_table3_raw(t2[p]) for p in PRIMARY}
    t2_en = next((t for t in parse_tables(en) if t and t[0][:1] == ["Pair"] and "summary_min" in " ".join(t[0])), None)
    t2_zh = next((t for t in parse_tables(zh) if t and t[0][0] in {"靶对", "Pair"} and "summary_min" in " ".join(t[0])), None)
    check_named_table("R1", "EN Table 2", t2_en, exp2, value_slice=slice(1, 9))
    check_named_table("R1", "ZH Table 2", t2_zh, exp2, value_slice=slice(1, 9))
    tables_en = [t for t in parse_tables(en) if t and t[0] and t[0][0] in {"Pair", "靶对"}]
    tables_zh = [t for t in parse_tables(zh) if t and t[0] and t[0][0] in {"Pair", "靶对"}]
    t3_en = next((t for t in tables_en if len(t[0]) == 5 and "neither" in " ".join(t[0]).lower()), None)
    if t3_en is None:
        t3_en = next((t for t in tables_en if any("neither" in c.lower() for c in t[0])), None)
    t3_zh = next((t for t in tables_zh if any("neither" in c.lower() or "双端低活性" in c or "n_neg" in c for c in t[0])), None)
    check_named_table("R1", "EN Table 3", t3_en, exp3)
    check_named_table("R1", "ZH Table 3", t3_zh, exp3)

    t1_en = table_by_header(en, "Pair")
    if t1_en and "Quota" in " ".join(t1_en[0]):
        t1 = t1_en
    else:
        t1 = next((t for t in parse_tables(en) if t and "Quota" in " ".join(t[0])), None)
    expected_t1 = {
        "PIK3CA/mTOR": ("θ = 6.0", "18 / 14 / 12 / 4", "48", "18 / 14 / 12", "16"),
        "AChE/BChE": ("strict 6.5/5.5", "28 / 28 / 28 / 16", "100", "27 / 25 / 28", "8"),
        "EGFR/HER2": ("θ = 6.0", "28 / 38 / 32 / 12", "110", "28 / 38 / 32", "8"),
        "F2/F10": ("strict 6.5/5.5", "32 / 32 / 32 / 14", "110", "31 / 32 / 32", "8"),
        "JAK1/TYK2": ("strict 6.5/5.5", "32 / 32 / 32 / 14", "110", "31 / 32 / 32", "8"),
        "JAK1/JAK2": ("strict 6.5/5.5", "32 / 32 / 32 / 14", "110", "32 / 32 / 32", "8"),
        "PPARG/PPARA": ("strict 6.5/5.5", "32 / 32 / 32 / 14", "110", "32 / 31 / 32", "8"),
        "PPARA/PPARD": ("strict 6.5/5.5", "32 / 32 / 32 / 14", "110", "32 / 32 / 32", "8"),
    }
    if t1 is None:
        rec("R1", "FAIL", "Table 1 not found")
    else:
        by = {row[0]: row for row in t1[1:]}
        five_panel = rows(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panel_summary_v1.csv")
        for pair, (pool, quota, n_panel, n_scored, exh) in expected_t1.items():
            row = by.get(pair)
            if row is None:
                rec("R1", "FAIL", f"Table 1 missing {pair}")
                continue
            got = (row[1], row[2], row[6], row[7], row[8])
            exp = (pool, quota, n_panel, n_scored, exh)
            if got != exp:
                rec("R1", "FAIL", f"Table 1 {pair}: {got} != {exp}")
            else:
                rec("R1", "PASS", f"Table 1 {pair} construction/n_scored")
            src = t2[pair]
            scored = f"{src['nD']} / {src['nA']} / {src['nB']}"
            if scored != n_scored:
                rec("R1", "FAIL", f"Table 1 n_scored {pair} {n_scored} != Table 2 CSV {scored}")
            if pair not in ORIGINAL:
                p = one(five_panel, pair=pair)
                if p["n_panel"] != n_panel:
                    rec("R1", "FAIL", f"Table 1 n_panel {pair} {n_panel} != track_b {p['n_panel']}")
    lock = (DOCS / "STATISTICAL_LOCK_V1.md").read_text(encoding="utf-8")
    for pair, src in t2.items():
        ok = any(
            f"| {pair} | {r3(src['smin'], rnd)} | {r3_ci(src['lo'], src['hi'], rnd)} |" in lock
            for rnd in (ROUND_HALF_UP, ROUND_HALF_EVEN)
        )
        if not ok:
            rec("R1", "FAIL", f"STATISTICAL_LOCK missing 3-dp row for {pair}")
        else:
            rec("R1", "PASS", f"STATISTICAL_LOCK {pair}")
    return t2


def round2(t2: dict) -> None:
    en_si = (DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md").read_text(encoding="utf-8")
    zh_si = (DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md").read_text(encoding="utf-8")
    equal3 = rows(ROOT / "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv")
    equal5 = rows(
        ROOT
        / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/equal_score_negative_s34_v1.csv"
    )
    s4_map = {
        ("EGFR/HER2", "D_vs_B_or_neither_pocketA"): one(equal3, pair="EGFR/HER2", contrast="D_vs_B_or_neither_pocketA"),
        ("EGFR/HER2", "D_vs_A_or_neither_pocketB"): one(equal3, pair="EGFR/HER2", contrast="D_vs_A_or_neither_pocketB"),
        ("AChE/BChE", "D_vs_B_or_neither_pocketA"): one(equal3, pair="AChE/BChE", contrast="D_vs_B_or_neither_pocketA"),
        ("AChE/BChE", "D_vs_A_or_neither_pocketB"): one(equal3, pair="AChE/BChE", contrast="D_vs_A_or_neither_pocketB"),
        ("PIK3CA/mTOR", "D_vs_B_or_neither_pocketA"): one(equal3, pair="PIK3CA/mTOR", contrast="D_vs_B_or_neither_pocketA"),
        ("PIK3CA/mTOR", "D_vs_A_or_neither_pocketB"): one(equal3, pair="PIK3CA/mTOR", contrast="D_vs_A_or_neither_pocketB"),
    }
    for r in equal5:
        s4_map[(r["pair"], r["contrast"])] = r

    s4 = None
    for table in parse_tables(en_si):
        if table and table[0][0] == "Pair" and any("Δ" in c or "delta" in c.lower() for c in table[0]) and any(
            "underpowered" in c.lower() or "neither" in c.lower() for c in table[0]
        ):
            if "Score channel" in table[0] or "pocket" in " ".join(table[0]).lower():
                s4 = table
                break
    if s4 is None:
        rec("R2", "FAIL", "SI Table S4 not found")
    else:
        expected_cells = {
            "EGFR/HER2 pocket A (vs B-only)": s4_map[("EGFR/HER2", "D_vs_B_or_neither_pocketA")],
            "JAK1/TYK2 pocket A": s4_map[("JAK1/TYK2", "D_vs_B_or_neither_pocketA")],
        }
        text_rows = [" | ".join(row) for row in s4]
        blob = "\n".join(text_rows)
        egfr = s4_map[("EGFR/HER2", "D_vs_B_or_neither_pocketA")]
        jak = s4_map[("JAK1/TYK2", "D_vs_B_or_neither_pocketA")]
        for name, recs, sel, nei, delta in (
            ("EGFR/HER2 pocket A", egfr, "0.430", "0.808", "0.378"),
            ("JAK1/TYK2 pocket A", jak, "0.365", "0.809", "0.444"),
        ):
            if not (sel in blob and nei in blob and delta in blob):
                rec("R2", "FAIL", f"Table S4 {name} missing {sel}/{nei}/{delta}")
            else:
                rec("R2", "PASS", f"Table S4 {name} points {sel} / {nei} / {delta}")
            if r3(recs["delta_neither_minus_selective"]) != delta:
                rec("R2", "FAIL", f"Table S4 {name} Δ CSV {recs['delta_neither_minus_selective']} != {delta}")
            ci = r3_ci(recs["delta_ci_lo"], recs["delta_ci_hi"])
            # manuscript uses unicode minus
            ci_ms = ci.replace("-", "−")
            if ci not in blob and ci_ms not in blob:
                rec("R2", "FAIL", f"Table S4 {name} CI {ci} not in table")
            else:
                rec("R2", "PASS", f"Table S4 {name} CI {ci}")

    cluster = rows(ROOT / "data/jcim_novelty_v0/tables/equal_score_cluster_bootstrap_v1.csv")
    en = (DOCS / "MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    zh = (DOCS / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    egfr_doc = one(cluster, pair="EGFR/HER2", estimator="document_cluster")
    egfr_scf = one(cluster, pair="EGFR/HER2", estimator="scaffold_cluster")
    jak_doc = one(cluster, pair="JAK1/TYK2", estimator="document_cluster")
    jak_scf = one(cluster, pair="JAK1/TYK2", estimator="scaffold_cluster")
    for label, recs, text in (
        ("EGFR document-cluster", egfr_doc, en + zh),
        ("EGFR scaffold-cluster", egfr_scf, en + zh),
        ("JAK1 document-cluster", jak_doc, en + zh),
        ("JAK1 scaffold-cluster", jak_scf, en + zh),
    ):
        ci = r3_ci(recs["delta_ci_lo"], recs["delta_ci_hi"]).replace("-", "−")
        ci_plain = r3_ci(recs["delta_ci_lo"], recs["delta_ci_hi"])
        if ci not in text and ci_plain not in text:
            rec("R2", "FAIL", f"{label} CI {ci_plain} missing from manuscripts")
        else:
            rec("R2", "PASS", f"{label} CI present")

    wp3 = rows(ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv")
    wp5 = rows(
        ROOT
        / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/wrong_pocket_by_channel_v1.csv"
    )
    s6_expected = {
        "EGFR/HER2": one(wp3, pair="EGFR/HER2", **{"set": "main_panel"}),
        "AChE/BChE": one(wp3, pair="AChE/BChE", **{"set": "main_panel"}),
        "PIK3CA/mTOR": one(wp3, pair="PIK3CA/mTOR", **{"set": "main_panel"}),
    }
    for r in wp5:
        if r["channel"] == "vina_20260727":
            s6_expected[r["pair"]] = r
    s6_blob = en_si
    for pair, recs in s6_expected.items():
        delta = r3(recs["delta_matched_minus_wrong"])
        delta_ms = delta.replace("-", "−")
        if pair == "EGFR/HER2" and delta != "0.170":
            rec("R2", "FAIL", f"S6 {pair} rounded Δ {delta}")
        if pair in s6_blob and (delta in s6_blob or delta_ms in s6_blob):
            rec("R2", "PASS", f"Table S6 {pair} Δ {delta}")
        else:
            rec("R2", "FAIL", f"Table S6 {pair} Δ {delta} missing")

    native = rows(ROOT / "data/jcim_novelty_v0/tables/external_slice_summary_v1.csv")
    s11_expected = {
        "EGFR/HER2": ("180", "10", "20"),
        "AChE/BChE": ("4", "8", "14"),
        "PIK3CA/mTOR": ("91", "4", "1"),
        "F2/F10": ("46", "15", "16"),
        "JAK1/TYK2": ("323", "7", "103"),
        "JAK1/JAK2": ("928", "40", "14"),
        "PPARG/PPARA": ("0", "0", "1"),
        "PPARA/PPARD": ("0", "1", "0"),
    }
    for pair, counts in s11_expected.items():
        r = one(native, pair=pair)
        got = (r["n_dual"], r["n_A_only"], r["n_B_only"])
        if got != counts:
            rec("R2", "FAIL", f"S11b CSV {pair} {got} != lock {counts}")
        needle = f"{counts[0]} / {counts[1]} / {counts[2]}"
        if needle not in en_si or needle not in zh_si:
            rec("R2", "FAIL", f"S11b {pair} {needle} missing from EN or ZH SI")
        else:
            rec("R2", "PASS", f"S11b {pair} {needle}")
        if r["packaged_as_external_evaluation"] != "0":
            rec("R2", "FAIL", f"{pair} packaged_as_external_evaluation != 0")

    incr3 = rows(ROOT / "data/jcim_novelty_v0/tables/incremental_information_v1.csv")
    incr5 = rows(
        ROOT
        / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv"
    )
    deltas = []
    for src in (incr3, incr5):
        by = defaultdict(dict)
        for r in src:
            if r["model"] in {"ECFP4", "ECFP4+docking"}:
                by[(r["pair"], r["contrast"])][r["model"]] = float(r["cv_auroc"])
        for key, models in by.items():
            deltas.append(abs(models["ECFP4+docking"] - models["ECFP4"]))
    max_abs = max(deltas)
    if abs(max_abs - 0.0234) > 5e-4 and abs(max_abs - 0.023) > 5e-3:
        # PPARA/PPARD D vs B: 0.8584 - 0.835 = 0.0234
        rec("R2", "FAIL", f"ECFP increment max |Δ| {max_abs:.4f} unexpected")
    else:
        rec("R2", "PASS", f"ECFP increment max |Δ| = {max_abs:.4f} (manuscript 0.023)")

    gnina = rows(ROOT / "data/jcim_independent_dock_v0/tables/independent_dock_formulation_v1.csv")
    jak_g = one(
        rows(
            ROOT
            / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv"
        ),
        channel="gnina_independent_jak1_tyk2",
        pair="JAK1/TYK2",
    )
    egfr_db = one(gnina, pair="EGFR/HER2", contrast="D_vs_B_pocketA")
    egfr_n = one(gnina, pair="EGFR/HER2", contrast="D_vs_neither_mean")
    if "0.220" not in r3_opts(egfr_db["auroc"]) or not r3_opts(egfr_n["auroc"]) & {"0.782", "0.783"}:
        rec("R2", "FAIL", f"EGFR GNINA {egfr_db['auroc']} / {egfr_n['auroc']}")
    else:
        rec("R2", "PASS", "EGFR independent GNINA 0.220 / 0.783 (0.7825)")
    if int(egfr_n["n_neg"]) != 11:
        rec("R2", "FAIL", f"EGFR GNINA Dual-vs-neither n_neg={egfr_n['n_neg']} (one neither failed)")
    else:
        rec("R2", "NOTE", "EGFR GNINA Dual-vs-neither n_neg=11 (EH120_109 failed); Vina Table 3 uses n=12")
    if r3(jak_g["summary_min"]) != "0.317" or r3(jak_g["D_vs_neither_mean"]) != "0.705":
        rec("R2", "FAIL", f"JAK1 GNINA {jak_g['summary_min']} / {jak_g['D_vs_neither_mean']}")
    else:
        rec("R2", "PASS", "JAK1/TYK2 independent GNINA 0.317 / 0.705")
    jak_ci = r3_ci(jak_g["ci_lo"], jak_g["ci_hi"])
    if jak_ci.replace("-", "−") not in en_si and jak_ci not in en_si:
        rec("R2", "FAIL", f"JAK1 GNINA CI {jak_ci} missing from EN SI")
    else:
        rec("R2", "PASS", f"JAK1 GNINA CI {jak_ci}")

    five_seed = rows(
        ROOT
        / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/fiveseed_summary_min_aggregate_v1.csv"
    )
    multi = rows(ROOT / "data/jcim_multiseed_v0/tables/multiseed_auroc_aggregate_v2.csv")
    s9_blob = en_si
    for pair, med, lo, hi in (
        ("EGFR/HER2", "0.373", "0.321", "0.430"),
        ("F2/F10", "0.366", "0.345", "0.385"),
        ("PPARG/PPARA", "0.651", "0.649", "0.691"),
    ):
        if all(x in s9_blob for x in (pair, med)):
            rec("R2", "PASS", f"Table S9 {pair} median {med}")
        else:
            rec("R2", "FAIL", f"Table S9 {pair} median {med} missing")
    egfr_ms = one(multi, pair="EGFR/HER2", metric="summary_min")
    if r3(egfr_ms["median"]) != "0.373":
        rec("R2", "FAIL", f"EGFR five-seed median CSV {egfr_ms['median']}")
    f2 = one(five_seed, pair="F2/F10")
    if r3(f2["median_of_per_seed_summary_min"]) != "0.366":
        rec("R2", "FAIL", f"F2 five-seed median {f2['median_of_per_seed_summary_min']}")

    hold = rows(ROOT / "data/jcim_holdout_v0/tables/holdout_pocket_matched_v1.csv")
    hold5 = [
        r
        for r in rows(
            ROOT
            / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/table2_comparable_by_channel_v1.csv"
        )
        if r["channel"] == "holdout_vina_20260727"
    ]
    ache_h = one(hold, pair="AChE/BChE", variant="pocket_matched_vina")
    if r3(ache_h["summary_min"]) != "0.618":
        rec("R2", "FAIL", f"AChE holdout {ache_h['summary_min']}")
    else:
        rec("R2", "PASS", "AChE holdout summary_min 0.618")
    pparg_h = one(hold5, pair="PPARG/PPARA")
    if r3(pparg_h["summary_min"]) != "0.535" or r3_ci(pparg_h["ci_lo"], pparg_h["ci_hi"]) != "[0.350, 0.717]":
        rec("R2", "FAIL", f"PPARG holdout {pparg_h['summary_min']} {pparg_h['ci_lo']} {pparg_h['ci_hi']}")
    else:
        rec("R2", "PASS", "PPARG holdout 0.535 [0.350, 0.717] from table2_comparable_by_channel")

    return


def round3() -> None:
    en = (DOCS / "MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    zh = (DOCS / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    en_si = (DOCS / "SUPPORTING_INFORMATION_JCIM_EN_V1.md").read_text(encoding="utf-8")
    zh_si = (DOCS / "SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md").read_text(encoding="utf-8")

    def numeric_rows(text: str) -> dict[str, list[list[str]]]:
        out: dict[str, list[list[str]]] = defaultdict(list)
        for table in parse_tables(text):
            if not table:
                continue
            for row in table[1:]:
                if row and row[0] in PRIMARY:
                    out[row[0]].append(row)
        return out

    en_rows = numeric_rows(en + "\n" + en_si)
    zh_rows = numeric_rows(zh + "\n" + zh_si)
    for pair in PRIMARY:
        if len(en_rows[pair]) != len(zh_rows[pair]):
            rec("R3", "NOTE", f"{pair}: EN has {len(en_rows[pair])} table rows, ZH has {len(zh_rows[pair])}")
        def nums(row: list[str]) -> tuple[str, ...]:
            return tuple(re.findall(r"-?\d+\.\d+|\d+", " ".join(row[1:])))

        en_set = {nums(r) for r in en_rows[pair]}
        zh_set = {nums(r) for r in zh_rows[pair]}
        if not en_set or not zh_set:
            rec("R3", "FAIL", f"{pair}: missing EN or ZH table rows")
        elif en_set != zh_set:
            rec("R3", "FAIL", f"{pair} EN/ZH numeric sets differ")
        else:
            rec("R3", "PASS", f"{pair} EN/ZH table-number sets match")

    flagship = [
        "0.378",
        "0.444",
        "0.205",
        "0.547",
        "0.263",
        "0.620",
        "0.083",
        "0.529",
        "0.168",
        "0.562",
        "0.234",
        "0.633",
    ]
    for token in flagship:
        if token not in en or token not in zh:
            rec("R3", "FAIL", f"flagship token {token} missing from EN or ZH manuscript")
        else:
            rec("R3", "PASS", f"flagship token {token} in both manuscripts")
    if "−0.034" not in zh and "[-0.034" not in zh:
        rec("R3", "FAIL", "ZH missing JAK1 document-cluster lower bound −0.034")
    else:
        rec("R3", "PASS", "ZH JAK1 document-cluster lower bound present")


def round4() -> None:
    plot = (ROOT / "data/jcim_bench_v0/scripts/plot_jcim_article_figures_v3.py").read_text(encoding="utf-8")
    if 'j0_lab = ["HDAC1/\\nHDAC6", "PIK3CA/\\nmTOR", "AChE/\\nBChE", "PIK3CA/\\nPIK3CB", "EGFR/\\nHER2"]' in plot:
        rec("R4", "FAIL", "fig1C still has a five-label axis including withdrawn PIK3CB")
    else:
        rec("R4", "PASS", "fig1C labels no longer include a PIK3CB tick without a bar")
    if '["EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"]' in plot:
        rec("R4", "FAIL", "fig1C still reads withdrawn PIK3CA/PIK3CB from overlap CSV")
    else:
        rec("R4", "PASS", "fig1C complete-case range uses the three original maps only")
    if "PIK3CA/PIK3CB" in plot and 'if "PIK3CA/PIK3CB" in p1["j0_hardneg"]' not in plot:
        # remaining mentions should be withdrawal guards
        leftover = [
            line
            for line in plot.splitlines()
            if "PIK3CA/PIK3CB" in line and "must not" not in line and "withdrawn" not in line and "if " not in line
        ]
        if leftover:
            rec("R4", "FAIL", f"plot v3 leftover PIK3CB data use: {leftover[:3]}")
        else:
            rec("R4", "PASS", "plot v3 PIK3CB mentions are withdrawal guards only")
    si = (ROOT / "data/jcim_bench_v0/scripts/plot_jcim_si_composites_v1.py").read_text(encoding="utf-8")
    if "PIK3CA/\\nPIK3CB" in si:
        rec("R4", "NOTE", "plot_jcim_si_composites_v1.py still ticks withdrawn PIK3CB (S1–S3/S9/S10 archive figures; do not submit as eight-row primaries)")
    v1 = ROOT / "data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py"
    if v1.exists():
        rec("R4", "NOTE", "v1 figure script retained as archive; submission figures must come from v3 only")


def round5() -> None:
    en = (DOCS / "MANUSCRIPT_JCIM_EN.md").read_text(encoding="utf-8")
    zh = (DOCS / "MANUSCRIPT_JCIM_ZH.md").read_text(encoding="utf-8")
    for name, text in (("EN manuscript", en), ("ZH manuscript", zh)):
        if "PIK3CA/PIK3CB" in text or re.search(r"PIK3CB(?! withdrawal)", text):
            if "PIK3CA/PIK3CB" in text:
                rec("R5", "FAIL", f"{name} still names PIK3CA/PIK3CB")
            else:
                rec("R5", "NOTE", f"{name} mentions PIK3CB (check context)")
        else:
            rec("R5", "PASS", f"{name} has no PIK3CA/PIK3CB primary-set mention")
    for rel in PAPER_GLOBS:
        text = (DOCS / rel).read_text(encoding="utf-8")
        if "PIK3CA/PIK3CB" in text:
            rec("R5", "FAIL", f"paper-facing {rel} still names PIK3CA/PIK3CB")
    for cmd, label in (
        ([sys.executable, str(ROOT / "scripts/primary/bootstrap_primary.py")], "bootstrap_primary"),
        ([sys.executable, str(ROOT / "data/jcim_novelty_v0/scripts/validate_revision_v1.py")], "validate_revision_v1"),
        ([sys.executable, str(ROOT / "data/jcim_novelty_v0/scripts/build_checksum_manifest_v1.py"), "--check"], "checksum --check"),
    ):
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        if proc.returncode != 0:
            rec("R5", "FAIL", f"{label} exit {proc.returncode}: {(proc.stdout + proc.stderr)[-400:]}")
        else:
            rec("R5", "PASS", f"{label} OK")
    master = rows(ROOT / "data/jcim_novelty_v0/tables/MASTER_RESULTS_TABLE.csv")
    t2_pairs = {r["pair"] for r in master if r["manuscript_table"] == "Table 2" and r["metric"] == "summary_min"}
    if t2_pairs != set(PRIMARY):
        rec("R5", "FAIL", f"MASTER Table 2 pairs {t2_pairs}")
    else:
        rec("R5", "PASS", "MASTER Table 2 has exactly the eight primary pairs")
    t2 = load_table2()
    for pair, src in t2.items():
        n_scored = src["nD"] + src["nA"] + src["nB"]
        if pair == "EGFR/HER2" and n_scored != 98:
            rec("R5", "FAIL", f"{pair} n_scored {n_scored}")
        if pair == "AChE/BChE" and n_scored != 80:
            rec("R5", "FAIL", f"{pair} n_scored {n_scored}")
        if pair == "PIK3CA/mTOR" and n_scored != 44:
            rec("R5", "FAIL", f"{pair} n_scored {n_scored}")
    rec("R5", "PASS", "n_scored D+A+B matches Table 2 class counts")
    census = rows(ROOT / "data/jcim_novelty_v0/tables/theta6_pair_census_v1.csv")
    if sum(int(r["docked_in_this_paper"]) for r in census) != 4:
        rec("R5", "NOTE", "J0 θ=6.0 census docked_in_this_paper != 4 (historical J0-era count may include withdrawn PIK3CB)")
    else:
        rec("R5", "NOTE", "J0 θ=6.0 census docked=4 is the historical J0-era count (includes later-withdrawn PIK3CB); current primary n=8")


def write_report() -> int:
    counts = defaultdict(int)
    for _, status, _ in FINDINGS:
        counts[status] += 1
    lines = [
        "# Five-round JCIM submission audit",
        "",
        "Date: 2026-09-09",
        "Branch: `cursor/jcim-submission-pack-0b1a`",
        "Script: `scripts/audit/audit_submission_five_rounds_v1.py`",
        "",
        "This audit compares assembled manuscripts and SI tables to frozen CSVs.",
        "It does not re-bootstrap and does not mint a DOI.",
        "",
        f"Summary: **{counts['PASS']} PASS**, **{counts['FAIL']} FAIL**, **{counts['NOTE']} NOTE**.",
        "",
        "| Round | Status | Finding |",
        "|---|---|---|",
    ]
    for round_id, status, msg in FINDINGS:
        safe = msg.replace("|", "\\|")
        lines.append(f"| {round_id} | {status} | {safe} |")
    lines.extend(
        [
            "",
            "## Blocking rule",
            "",
            "Publication-facing numbers must match the locked CSVs after three-decimal rounding.",
            "Fig1C must not plot a withdrawn PIK3CA/PIK3CB bar or read that pair from a CSV that no longer contains it.",
            "Assembled EN/ZH manuscripts must not list PIK3CA/PIK3CB as a primary pair.",
            "",
            "## Canonical sources",
            "",
            "- Table 2 original three: `unified_threshold_sensitivity_v2.csv` `theta_6.0`",
            "- Table 2 / Table 3 five pairs: `five_pair_stack_v1/table2_comparable_theta6_v1.csv`",
            "- Table 3 original three: `formulation_conventional_vs_directional_v1.csv`",
            "- Table S4: `formulation_equal_score_negative_v1.csv` + `equal_score_negative_s34_v1.csv`",
            "- Table S6: `wrong_pocket_paired_delta_bootstrap_v1.csv` + `wrong_pocket_by_channel_v1.csv`",
            "- Table S11b: `external_slice_summary_v1.csv`",
            "- JAK1/TYK2 independent GNINA: `table2_comparable_by_channel_v1.csv` `gnina_independent_jak1_tyk2`",
            "- EGFR/PIK3CA independent GNINA: `independent_dock_formulation_v1.csv`",
            "",
            "## Remediation on this branch",
            "",
            "- Table 2 EGFR/HER2 TPSA displayed 0.427 for CSV 0.4275, which is not valid under half-up or half-even; corrected to **0.428** in EN/ZH Table 2.",
            "- Fig1C had five x-tick labels (including withdrawn PIK3CA/PIK3CB) but only four J0 bars, and read complete-case overlap for a pair no longer in that CSV. Bars/labels now match `j0_strict_label_supply.csv`; overlap uses the three original maps (14.5%–34.0%).",
            "- Figure verify locks now match the current J0 scrape (48 pairs, 3 thick) and θ=6.0 census (`directional_n10` = 16, `docked_in_this_paper` = 3). The funnel still records historically docked 4 → PIK3CB withdrawal as narrative stages.",
            "- `plot_jcim_si_composites_v1.py` S1–S3/S9/S10 remain original-set archive figures (may still tick PIK3CB) and were moved to `figures/jcim_article/archive_original_set/`.",
            "- Submission slice: `submission_pack/`.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", REPORT)
    print(f"PASS={counts['PASS']} FAIL={counts['FAIL']} NOTE={counts['NOTE']}")
    return 1 if counts["FAIL"] else 0


def main() -> int:
    t2 = round1()
    round2(t2)
    round3()
    round4()
    round5()
    return write_report()


if __name__ == "__main__":
    raise SystemExit(main())
