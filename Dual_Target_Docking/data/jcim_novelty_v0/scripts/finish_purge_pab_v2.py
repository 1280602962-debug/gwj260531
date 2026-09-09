#!/usr/bin/env python3
"""Second-pass purge of the withdrawn PIK3CA/PIK3CB pair.

Removes leftover mixed CSV tokens, HOAP holdout artifacts, and exact
pair-list entries from Python sources. Does not touch PIK3CA/mTOR.
Does not delete PIK3CB/MTOR or PIK3CB/AKT1 universe-census rows.
"""
from __future__ import annotations

import ast
import csv
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

PAIR_EXACT = {"PIK3CA/PIK3CB", "PIK3CA_PIK3CB", "PIK3CA-PIK3CB"}
PAIR_TOKEN_RE = re.compile(r"^PIK3CA\s*/\s*PIK3CB$")
SKIP_DIRS = {".git", "node_modules", "__pycache__", ".venv"}

DELETE_PATHS = [
    "data/jcim_holdout_v0/HOAP",
    "data/jcim_holdout_v0/tables/scores_vina_mode1_HOAP.csv",
    "data/jcim_holdout_v0/tables/holdout_panel_HOAP.csv",
    "data/jcim_holdout_v0/tables/strict_pool_full_HOAP.csv",
    "data/jcim_holdout_v0/logs/dock_HOAP_local.log",
]

SKIP_PY = {
    Path("data/jcim_novelty_v0/scripts/purge_pik3ca_pik3cb_v1.py"),
    Path("data/jcim_novelty_v0/scripts/finish_purge_pab_v2.py"),
    # Superseded plotters; v3 is the article figure script.
    Path("data/jcim_bench_v0/scripts/plot_jcim_article_figures_v1.py"),
    Path("data/jcim_bench_v0/scripts/plot_jcim_article_figures_v2.py"),
}


def split_multi(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"[;|]", value or "") if part.strip()]


def is_withdrawn_pair_token(text: str) -> bool:
    return text in PAIR_EXACT or bool(PAIR_TOKEN_RE.match(text))


def row_is_withdrawn_only(row: dict) -> bool:
    identity_keys = ("pair", "system", "panel", "name", "target_pair", "pair_label")
    for key in identity_keys:
        if key in row and is_withdrawn_pair_token(str(row.get(key) or "").strip()):
            return True
    pdb = str(row.get("pdb") or row.get("pdb_id") or "").upper()
    if pdb == "2WXF":
        return True
    pairs = split_multi(str(row.get("pairs") or ""))
    if pairs and all(is_withdrawn_pair_token(part) for part in pairs):
        return True
    ligands = split_multi(str(row.get("ligands") or row.get("ligand_ids") or ""))
    if ligands and all(part.startswith(("PAB_", "HOAP_")) for part in ligands):
        other_pairs = [part for part in pairs if not is_withdrawn_pair_token(part)]
        if not other_pairs:
            return True
    return False


def clean_row_fields(row: dict) -> dict:
    for key, value in list(row.items()):
        text = str(value or "")
        if not text:
            continue
        if ";" not in text and key not in {"pairs", "targets", "ligands", "ligand_ids", "panel_ids"}:
            continue
        parts = split_multi(text)
        if len(parts) <= 1 and key not in {"pairs", "targets", "ligands"}:
            continue
        kept = []
        for part in parts:
            if is_withdrawn_pair_token(part):
                continue
            if part.startswith(("PAB_", "HOAP_")):
                continue
            if part == "PIK3CB" and key in {"targets", "genes", "target"}:
                continue
            kept.append(part)
        if kept != parts:
            row[key] = ";".join(kept)
    return row


def filter_csv(path: Path) -> str | None:
    try:
        raw = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return None
    if "PIK3CA/PIK3CB" not in raw and "PAB_" not in raw and "HOAP_" not in raw and "2WXF" not in raw:
        return None
    try:
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                return None
            fields = reader.fieldnames
            kept = []
            for row in reader:
                if row_is_withdrawn_only(row):
                    continue
                kept.append(clean_row_fields(row))
    except Exception:
        return None
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(kept)
    return f"csv cleaned {path.relative_to(ROOT)} kept={len(kept)}"


def skip_ws_and_comma(text: str, i: int) -> int:
    n = len(text)
    while i < n and text[i] in " \t\r\n":
        i += 1
    if i < n and text[i] == ",":
        i += 1
        while i < n and text[i] in " \t":
            i += 1
        if i < n and text[i] == "\n":
            i += 1
    return i


def consume_value(text: str, i: int) -> int:
    n = len(text)
    while i < n and text[i] in " \t\r\n":
        i += 1
    if i >= n:
        return i
    if text.startswith("dict(", i):
        i += 5
        depth = 1
        while i < n and depth:
            ch = text[i]
            if ch in "\"'":
                quote = ch
                i += 1
                while i < n and text[i] != quote:
                    if text[i] == "\\":
                        i += 2
                        continue
                    i += 1
                i += 1
                continue
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            i += 1
        return i
    if text[i] in "{[(":
        open_ch = text[i]
        close_ch = {"{": "}", "[": "]", "(": ")"}[open_ch]
        depth = 1
        i += 1
        while i < n and depth:
            ch = text[i]
            if ch in "\"'":
                quote = ch
                i += 1
                while i < n and text[i] != quote:
                    if text[i] == "\\":
                        i += 2
                        continue
                    i += 1
                i += 1
                continue
            if ch == open_ch:
                depth += 1
            elif ch == close_ch:
                depth -= 1
            i += 1
        return i
    while i < n and text[i] not in ",\n":
        i += 1
    return i


def remove_keyed_entries(text: str, keys: tuple[str, ...]) -> str:
    changed = True
    while changed:
        changed = False
        for key in keys:
            pattern = re.compile(rf'(?m)^(?P<indent>[ \t]*)(?P<q>["\']){re.escape(key)}(?P=q)\s*:')
            match = pattern.search(text)
            if not match:
                continue
            start = match.start()
            value_start = match.end()
            end = consume_value(text, value_start)
            end = skip_ws_and_comma(text, end)
            text = text[:start] + text[end:]
            changed = True
    return text


def remove_list_tokens(text: str) -> str:
    text = re.sub(r',\s*["\']PIK3CA/PIK3CB["\']', "", text)
    text = re.sub(r'["\']PIK3CA/PIK3CB["\']\s*,\s*', "", text)
    text = re.sub(r',\s*["\']HOAP["\']', "", text)
    text = re.sub(r'["\']HOAP["\']\s*,\s*', "", text)
    return text


def strip_python(path: Path) -> str | None:
    rel = path.relative_to(ROOT)
    if rel in SKIP_PY:
        return None
    try:
        original = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return None
    if "PIK3CA/PIK3CB" not in original and '"HOAP"' not in original and "'HOAP'" not in original:
        return None
    updated = remove_keyed_entries(original, ("PIK3CA/PIK3CB", "HOAP"))
    updated = remove_list_tokens(updated)
    if updated == original:
        return None
    try:
        ast.parse(updated)
    except SyntaxError as exc:
        return f"SKIP syntax {rel}: {exc}"
    path.write_text(updated, encoding="utf-8")
    return f"py stripped {rel}"


def walk_files(suffixes: set[str]):
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix.lower() in suffixes:
            yield path


def main() -> int:
    reports = []
    for rel in DELETE_PATHS:
        path = ROOT / rel
        if path.is_dir():
            shutil.rmtree(path)
            reports.append(f"deleted dir {rel}")
        elif path.is_file():
            path.unlink()
            reports.append(f"deleted file {rel}")

    hoap_sh = ROOT / "data/jcim_holdout_v0/logs/run_hoap_then_hopm.sh"
    if hoap_sh.exists():
        text = hoap_sh.read_text(encoding="utf-8")
        text = re.sub(r"echo \"=== HOAP resume.*?ec_hoap=\$\{PIPESTATUS\[0\]\}\n", "", text, flags=re.S)
        text = text.replace("ALL_DONE hoap=$ec_hoap hopm=$ec_hopm", "ALL_DONE hopm=$ec_hopm")
        hoap_sh.write_text(text, encoding="utf-8")
        reports.append("rewrote run_hoap_then_hopm.sh")

    for path in walk_files({".csv"}):
        note = filter_csv(path)
        if note:
            reports.append(note)

    for path in walk_files({".py"}):
        note = strip_python(path)
        if note:
            reports.append(note)

    plotted = ROOT / "figures/jcim_article/plotted_values.json"
    if plotted.exists():
        data = json.loads(plotted.read_text(encoding="utf-8"))
        fig2 = data.get("fig2A") or data.get("plotted", {}).get("fig2A")
        # stages live under plotted.fig1C
        plotted_block = data.get("plotted") or data
        fig1c = None
        if isinstance(plotted_block, dict):
            fig1c = plotted_block.get("fig1C")
        if isinstance(fig1c, dict) and "stages" in fig1c:
            fig1c["stages"] = {
                "candidate scrape": fig1c.get("n_pairs", 48),
                "min hard-neg ≥50": fig1c.get("n_thick", 3),
                "after HDAC metal exclusion": 3,
                "primary rows": 8,
            }
            fig1c.pop("after PIK3CB withdrawal", None)
            plotted.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
            reports.append("scrubbed plotted_values.json stages")

    log = ROOT / "data/jcim_novelty_v0/analysis/PURGE_PIK3CA_PIK3CB_V2.md"
    log.write_text("# PIK3CA/PIK3CB purge pass 2\n\n" + "\n".join(f"- {r}" for r in reports) + "\n", encoding="utf-8")
    print(f"wrote {log} ({len(reports)} actions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
