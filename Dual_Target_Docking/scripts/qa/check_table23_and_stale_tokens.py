#!/usr/bin/env python3
"""Lock Table 2/3 manuscript cells to canonical CSVs and scan stale tokens.

Does not rewrite canonical results. Exit 1 if any Table cell mismatches
or any STALE_ERROR remains in current-facing publication files.
"""
from __future__ import annotations

import csv
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CANON = ROOT / "results" / "canonical"
QA = ROOT / "results" / "qa"
PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
CURRENT_FILES = [
    ROOT / "docs/MANUSCRIPT_JCIM_EN.md",
    ROOT / "docs/MANUSCRIPT_JCIM_ZH.md",
    ROOT / "docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    ROOT / "docs/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    ROOT / "docs/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
    ROOT / "docs/WRITING_INDEX_FREEZE.md",
    ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_EN.md",
    ROOT / "submission_pack/manuscript/MANUSCRIPT_JCIM_ZH.md",
    ROOT / "submission_pack/SI/SUPPORTING_INFORMATION_JCIM_EN_V1.md",
    ROOT / "submission_pack/SI/SUPPORTING_INFORMATION_DRAFT_ZH_JCIM_V1.md",
    ROOT / "submission_pack/manuscript/FIGURE_TABLE_LOCK_POSTFIX_V4.md",
]
TOKENS = [
    ("0.3237", "historical EGFR ablation summary_min"),
    ("0.554", "stale EGFR D-vs-neither CI lo"),
    ("0.354", "stale EGFR EF"),
    ("0.4275", "stale EGFR TPSA"),
    ("[0.235, 0.665]", "stale EGFR scaffold cluster CI"),
    ("[0.125, 0.644]", "stale EGFR document cluster CI"),
    ("1 / 5 / 5 / 0", "stale EGFR Top-10 composition"),
    ("1/5/5/0", "stale EGFR Top-10 composition compact"),
]


def r3(x) -> str:
    d = Decimal(str(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)
    return f"{d:.3f}"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_table(text: str, after_marker: str, header_token: str) -> dict[str, list[str]]:
    lines = text.splitlines()
    mark = None
    for i, line in enumerate(lines):
        if after_marker in line:
            mark = i
            break
    if mark is None:
        raise SystemExit(f"marker not found: {after_marker}")
    start = None
    for i, line in enumerate(lines[mark:], start=mark):
        if header_token in line and line.strip().startswith("|"):
            start = i
            break
    if start is None:
        raise SystemExit(f"table header not found after {after_marker}: {header_token}")
    rows = {}
    for line in lines[start + 2 :]:
        if not line.startswith("|"):
            break
        if set(line.replace("|", "").strip()) <= set("-: "):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells and cells[0] in PAIRS:
            rows[cells[0]] = cells
    if set(rows) != set(PAIRS):
        raise SystemExit(f"table pair set incomplete for {header_token}: {sorted(rows)}")
    return rows


def classify_token(ctx: str) -> str:
    low = ctx.lower()
    if any(
        m in low
        for m in (
            "historical",
            "superseded",
            "not current",
            "do not",
            "archived",
            "不得",
            "不是当前",
            "历史",
        )
    ):
        return "HISTORICAL_WITH_EXPLICIT_CONTEXT"
    return "STALE_ERROR"


def main() -> int:
    dau = {(r["pair"], r["estimand"]): r for r in read_csv(CANON / "primary_directional_auroc.csv")}
    smin = {r["pair"]: r for r in read_csv(CANON / "primary_summary_min.csv")}
    rank = {r["pair"]: r for r in read_csv(CANON / "two_pocket_mean_ranking.csv")}
    top = {r["pair"]: r for r in read_csv(CANON / "top10_operating_points.csv")}
    lock_rows = []
    stale = 0
    for path, label in (
        (ROOT / "docs/MANUSCRIPT_JCIM_EN.md", "EN"),
        (ROOT / "docs/MANUSCRIPT_JCIM_ZH.md", "ZH"),
    ):
        text = path.read_text(encoding="utf-8")
        t2 = parse_table(text, "**Table 2.**", "n_scored (dual / A-only / B-only)")
        t3 = parse_table(text, "**Table 3.**", "Top 10% D/A/B/N")
        for pair in PAIRS:
            da = dau[(pair, "AUROC_D_vs_A_pocketB")]
            db = dau[(pair, "AUROC_D_vs_B_pocketA")]
            sm = smin[pair]
            want_n = f"{int(sm['n_dual'])} / {int(sm['n_A_only'])} / {int(sm['n_B_only'])}"
            want_da = f"{r3(da['point'])} [{r3(da['ci_lo'])}, {r3(da['ci_hi'])}]"
            want_db = f"{r3(db['point'])} [{r3(db['ci_lo'])}, {r3(db['ci_hi'])}]"
            want_sm = f"{r3(sm['summary_min'])} [{r3(sm['ci_lo'])}, {r3(sm['ci_hi'])}]"
            got = t2[pair]
            for metric, want, gotv in (
                ("n_scored", want_n, got[1]),
                ("D_vs_A", want_da, got[2]),
                ("D_vs_B", want_db, got[3]),
                ("summary_min", want_sm, got[4]),
            ):
                status = "MATCH" if want == gotv else "MISMATCH"
                lock_rows.append(
                    {
                        "manuscript": label,
                        "table": "Table 2",
                        "pair": pair,
                        "metric": metric,
                        "canonical": want,
                        "typeset": gotv,
                        "status": status,
                    }
                )
                if status != "MATCH":
                    stale += 1
            rk = rank[pair]
            tp = top[pair]
            want_auc = (
                f"{r3(rk['two_pocket_mean_D_vs_neither'])} "
                f"[{r3(rk['ci_lo'])}, {r3(rk['ci_hi'])}]"
            )
            want_comp = (
                f"{int(tp['top_dual'])} / {int(tp['top_A_only'])} / "
                f"{int(tp['top_B_only'])} / {int(tp['top_neither'])}"
            )
            want_ef = r3(rk["ef_dual_10pct"])
            got3 = t3[pair]
            for metric, want, gotv in (
                ("D_vs_neither", want_auc, got3[1]),
                ("n_neither", str(int(float(tp["n_neither"]))), got3[2]),
                ("top10", want_comp, got3[3]),
                ("EF", want_ef, got3[4]),
            ):
                status = "MATCH" if want == gotv else "MISMATCH"
                lock_rows.append(
                    {
                        "manuscript": label,
                        "table": "Table 3",
                        "pair": pair,
                        "metric": metric,
                        "canonical": want,
                        "typeset": gotv,
                        "status": status,
                    }
                )
                if status != "MATCH":
                    stale += 1

    token_rows = []
    stale_err = 0
    for path in CURRENT_FILES:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        rel = str(path.relative_to(ROOT))
        for token, claim in TOKENS:
            start = 0
            while True:
                i = text.find(token, start)
                if i < 0:
                    break
                ctx = text[max(0, i - 70) : i + len(token) + 70].replace("\n", " ")
                status = classify_token(ctx)
                if status == "STALE_ERROR":
                    stale_err += 1
                token_rows.append(
                    {
                        "file": rel,
                        "token": token,
                        "claim": claim,
                        "status": status,
                        "context": ctx,
                    }
                )
                start = i + len(token)

    QA.mkdir(parents=True, exist_ok=True)
    with (QA / "table23_lock.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(lock_rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(lock_rows)
    with (QA / "stale_token_scan.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["file", "token", "claim", "status", "context"]
        w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(token_rows)

    n_lock = len(lock_rows)
    n_ok = sum(1 for r in lock_rows if r["status"] == "MATCH")
    print(f"Table 2/3 lock: {n_ok}/{n_lock} MATCH")
    for r in lock_rows:
        if r["status"] != "MATCH":
            print("MISMATCH", r)
    print(f"STALE_ERROR tokens: {stale_err}")
    for r in token_rows:
        if r["status"] == "STALE_ERROR":
            print("STALE_ERROR", r["file"], r["token"])
    if n_ok != n_lock or stale_err:
        return 1
    print("PASS: Table 2 MATCH; Table 3 MATCH; STALE_ERROR = 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
