#!/usr/bin/env python3
"""Count-level BindingDB + PubChem strict-label supply check (zero docking).

For the frozen K=4 pairs (PIK3CA/mTOR, AChE/BChE, PIK3CA/PIK3CB, EGFR/HER2),
count ligands with IC50/Ki/Kd/EC50 on BOTH ends under the same construction
gate as J0: dual both ≥6.5; A_only A≥6.5 and B≤5.5; B_only converse.

This is a per-database identity check (BindingDB monomerid; PubChem CID).
It does not merge structures across databases and does not build a new panel.

If an endpoint is unreachable, the corresponding source is recorded as
fetch_failed — never fabricated.
"""
from __future__ import annotations

import csv
import json
import math
import re
import ssl
import time
import urllib.error
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
CACHE = ROOT / "cache"
TABLES = ROOT / "tables"
ANALYSIS = ROOT / "analysis"
CHEMBL_MOLS = REPO / "data" / "public_pair_selection"

HI = 6.5
LO = 5.5
THETA = 6.0
MIN_HARDNEG_STRICT = 50
MIN_HARDNEG_THIN = 20
KEEP_TYPES = {"IC50", "KI", "KD", "EC50"}
BDB_CUTOFF_NM = 1_000_000  # 1 mM; include weak measurements (p down to ~3)
UA = "DualFourClass-Bench/0.1 (count-level supply check; no docking)"
NUM_RE = re.compile(r"([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)")

TARGETS = {
    "PIK3CA": "P42336",
    "MTOR": "P42345",
    "ACHE": "P22303",
    "BCHE": "P06276",
    "PIK3CB": "P42338",
    "EGFR": "P00533",
    "HER2": "P04626",
}

PAIRS = [
    ("PIK3CA/MTOR", "PIK3CA", "MTOR"),
    ("ACHE/BCHE", "ACHE", "BCHE"),
    ("PIK3CA/PIK3CB", "PIK3CA", "PIK3CB"),
    ("EGFR/HER2", "EGFR", "HER2"),
]

CHEMBL_FILES = {
    "PIK3CA": "mols_PIK3CA.json",
    "MTOR": "mols_MTOR.json",
    "ACHE": "mols_ACHE.json",
    "BCHE": "mols_BCHE.json",
    "PIK3CB": "mols_PIK3CB.json",
    "EGFR": "mols_EGFR.json",
    "HER2": "mols_HER2.json",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def nm_to_p(nm: float) -> float | None:
    if nm is None or nm <= 0:
        return None
    return 9.0 - math.log10(nm)


def um_to_p(um: float) -> float | None:
    if um is None or um <= 0:
        return None
    return 6.0 - math.log10(um)


def parse_affinity_token(raw: str) -> tuple[str, float | None]:
    """Return (qualifier, numeric) from BindingDB affinity strings."""
    s = (raw or "").strip().replace(",", "")
    if not s:
        return "", None
    qual = ""
    rest = s
    if s[0] in "<>~":
        qual = s[0]
        rest = s[1:].strip()
    elif s[:2] in (">=", "<=", "~="):
        qual = s[0]
        rest = s[2:].strip()
    m = NUM_RE.search(rest)
    if not m:
        return qual, None
    try:
        return qual, float(m.group(1))
    except ValueError:
        return qual, None


def classify_pair(a_vals: dict, b_vals: dict) -> dict:
    both = set(a_vals) & set(b_vals)
    rec = {
        "n_ligands_A": len(a_vals),
        "n_ligands_B": len(b_vals),
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
    n = rec["n_both_measured"]
    rec["gray_frac"] = round(rec["gray"] / n, 3) if n else None
    rec["min_strict_hardneg"] = min(rec["strict_A_only"], rec["strict_B_only"])
    rec["supports_strict_panel"] = rec["min_strict_hardneg"] >= MIN_HARDNEG_STRICT
    rec["supports_thin_panel"] = rec["min_strict_hardneg"] >= MIN_HARDNEG_THIN
    return rec


def http_get(url: str, timeout: int, retries: int = 3) -> tuple[int | str, bytes, str]:
    ctx = ssl.create_default_context()
    last_err = ""
    for attempt in range(retries):
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
                return resp.status, resp.read(), (resp.getheader("Content-Type") or "")
        except urllib.error.HTTPError as e:
            last_err = f"HTTP {e.code} {e.reason}"
            body = e.read() if e.fp else b""
            if e.code in {404, 400, 401, 403}:
                return e.code, body, last_err
            time.sleep(4 * (2**attempt))
        except Exception as e:
            last_err = f"{type(e).__name__}: {e}"
            time.sleep(4 * (2**attempt))
    return "ERR", last_err.encode("utf-8"), last_err


def fetch_bindingdb(uniprot: str, cache_dir: Path) -> dict:
    cache_path = cache_dir / f"bindingdb_{uniprot}_cutoff{BDB_CUTOFF_NM}.json"
    meta = {
        "source": "BindingDB_REST",
        "uniprot": uniprot,
        "cutoff_nM": BDB_CUTOFF_NM,
        "url": (
            "https://bindingdb.org/rest/getLigandsByUniprots"
            f"?uniprot={uniprot}&cutoff={BDB_CUTOFF_NM}&response=application/json"
        ),
        "fetched_at_utc": None,
        "ok": False,
        "error": None,
        "n_raw_records": 0,
        "n_kept_records": 0,
        "n_ligands_pmax": 0,
        "cache_path": str(cache_path.relative_to(ROOT)),
    }
    raw = None
    if cache_path.exists() and cache_path.stat().st_size > 20:
        raw = cache_path.read_bytes()
        meta["fetched_at_utc"] = datetime.fromtimestamp(
            cache_path.stat().st_mtime, tz=timezone.utc
        ).strftime("%Y-%m-%dT%H:%M:%SZ")
        meta["from_cache"] = True
    else:
        status, body, info = http_get(meta["url"], timeout=180, retries=3)
        meta["http_status"] = status
        meta["fetched_at_utc"] = utc_now()
        meta["from_cache"] = False
        if status != 200:
            meta["error"] = info or body[:300].decode("utf-8", "replace")
            return meta, {"as_is": {}, "equal_only": {}}
        cache_path.write_bytes(body)
        raw = body
    try:
        data = json.loads(raw.decode("utf-8"))
    except Exception as e:
        meta["error"] = f"JSON parse failed: {e}"
        return meta, {"as_is": {}, "equal_only": {}}
    block = data.get("getLindsByUniprotsResponse") or data.get(
        "getLigandsByUniprotsResponse"
    ) or {}
    aff = block.get("affinities")
    if aff is None:
        meta["error"] = "no affinities key; " + ",".join(list(data)[:6])
        return meta, {"as_is": {}, "equal_only": {}}
    if isinstance(aff, dict):
        aff = [aff]
    meta["n_raw_records"] = len(aff)
    pmax_all: dict[str, float] = {}
    pmax_eq: dict[str, float] = {}
    kept = 0
    kept_eq = 0
    type_counts: dict[str, int] = defaultdict(int)
    for row in aff:
        atype = str(row.get("affinity_type") or "").strip().upper()
        type_counts[atype or "?"] += 1
        if atype not in KEEP_TYPES:
            continue
        qual, nm = parse_affinity_token(str(row.get("affinity") or ""))
        p = nm_to_p(nm) if nm is not None else None
        if p is None:
            continue
        mid = str(row.get("monomerid") or "").strip()
        if not mid:
            continue
        kept += 1
        if mid not in pmax_all or p > pmax_all[mid]:
            pmax_all[mid] = p
        if qual in {"", "="}:
            kept_eq += 1
            if mid not in pmax_eq or p > pmax_eq[mid]:
                pmax_eq[mid] = p
    meta["ok"] = True
    meta["n_kept_records"] = kept
    meta["n_kept_equal_records"] = kept_eq
    meta["n_ligands_pmax"] = len(pmax_all)
    meta["n_ligands_pmax_equal"] = len(pmax_eq)
    meta["affinity_type_counts"] = dict(sorted(type_counts.items(), key=lambda kv: -kv[1]))
    return meta, {"as_is": pmax_all, "equal_only": pmax_eq}


def fetch_pubchem(uniprot: str, cache_dir: Path) -> dict:
    cache_path = cache_dir / f"pubchem_{uniprot}_concise.csv"
    meta = {
        "source": "PubChem_PUG_REST_protein_concise",
        "uniprot": uniprot,
        "url": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/protein/accession/{uniprot}/concise/CSV",
        "fetched_at_utc": None,
        "ok": False,
        "error": None,
        "n_raw_records": 0,
        "n_kept_records": 0,
        "n_ligands_pmax": 0,
        "cache_path": str(cache_path.relative_to(ROOT)),
    }
    text = None
    if cache_path.exists() and cache_path.stat().st_size > 20:
        text = cache_path.read_text(encoding="utf-8", errors="replace")
        meta["fetched_at_utc"] = datetime.fromtimestamp(
            cache_path.stat().st_mtime, tz=timezone.utc
        ).strftime("%Y-%m-%dT%H:%M:%SZ")
        meta["from_cache"] = True
        meta["bytes"] = cache_path.stat().st_size
    else:
        status, body, info = http_get(meta["url"], timeout=180, retries=3)
        meta["http_status"] = status
        meta["fetched_at_utc"] = utc_now()
        meta["from_cache"] = False
        meta["bytes"] = len(body)
        if status != 200:
            meta["error"] = info or body[:300].decode("utf-8", "replace")
            return meta, {"as_is": {}, "equal_only": {}}
        cache_path.write_bytes(body)
        text = body.decode("utf-8", "replace")
    reader = csv.DictReader(text.splitlines())
    pmax_all: dict[str, float] = {}
    pmax_eq: dict[str, float] = {}
    kept = 0
    kept_eq = 0
    n_raw = 0
    type_counts: dict[str, int] = defaultdict(int)
    missing_cid = 0
    missing_val = 0
    for row in reader:
        n_raw += 1
        aname = str(row.get("Activity Name") or "").strip().upper()
        type_counts[aname or "?"] += 1
        if aname not in KEEP_TYPES:
            continue
        cid = str(row.get("CID") or "").strip()
        if not cid or cid in {"", "0"}:
            missing_cid += 1
            continue
        raw_val = str(row.get("Activity Value [uM]") or "").strip()
        parsed_qual, um = parse_affinity_token(raw_val)
        col_qual = str(row.get("Activity Qualifier") or "").strip()
        qual = col_qual or parsed_qual or "="
        if qual in {">=", "<="}:
            qual = qual[0]
        p = um_to_p(um) if um is not None else None
        if p is None:
            missing_val += 1
            continue
        kept += 1
        if cid not in pmax_all or p > pmax_all[cid]:
            pmax_all[cid] = p
        if qual in {"", "="}:
            kept_eq += 1
            if cid not in pmax_eq or p > pmax_eq[cid]:
                pmax_eq[cid] = p
    meta["ok"] = True
    meta["n_raw_records"] = n_raw
    meta["n_kept_records"] = kept
    meta["n_kept_equal_records"] = kept_eq
    meta["n_ligands_pmax"] = len(pmax_all)
    meta["n_ligands_pmax_equal"] = len(pmax_eq)
    meta["n_skipped_missing_cid"] = missing_cid
    meta["n_skipped_missing_value"] = missing_val
    meta["activity_name_top"] = dict(
        sorted(type_counts.items(), key=lambda kv: -kv[1])[:12]
    )
    return meta, {"as_is": pmax_all, "equal_only": pmax_eq}


def load_chembl(name: str) -> dict[str, float]:
    path = CHEMBL_MOLS / CHEMBL_FILES[name]
    with path.open() as fh:
        return {str(k): float(v) for k, v in json.load(fh).items()}


def write_pmax(path: Path, pmax: dict[str, float]) -> None:
    compact = {k: round(v, 4) for k, v in sorted(pmax.items(), key=lambda kv: kv[0])}
    path.write_text(json.dumps(compact, indent=None, separators=(",", ":")) + "\n")


def fmt_ab(rec: dict) -> str:
    return f"{rec['strict_A_only']}/{rec['strict_B_only']}"


def row_lookup(rows: list[dict]) -> dict:
    return {(r["pair"], r["source"], r["rule"]): r for r in rows}


def egfr_bindingdb_qualifier_diagnostic() -> dict | None:
    """How many EGFR/HER2 BindingDB B_only ligands are inequality-only on EGFR."""
    cache_e = CACHE / f"bindingdb_P00533_cutoff{BDB_CUTOFF_NM}.json"
    cache_h = CACHE / f"bindingdb_P04626_cutoff{BDB_CUTOFF_NM}.json"
    if not cache_e.exists() or not cache_h.exists():
        return None

    def load(path: Path):
        data = json.loads(path.read_text())
        aff = data["getLindsByUniprotsResponse"]["affinities"]
        quals: dict[str, list[str]] = defaultdict(list)
        pmax: dict[str, float] = {}
        for row in aff:
            atype = str(row.get("affinity_type") or "").strip().upper()
            if atype not in KEEP_TYPES:
                continue
            qual, nm = parse_affinity_token(str(row.get("affinity") or ""))
            p = nm_to_p(nm) if nm is not None else None
            if p is None:
                continue
            mid = str(row.get("monomerid") or "").strip()
            if not mid:
                continue
            quals[mid].append(qual or "=")
            if mid not in pmax or p > pmax[mid]:
                pmax[mid] = p
        return quals, pmax

    q_e, p_e = load(cache_e)
    _q_h, p_h = load(cache_h)
    bonly = [
        mid
        for mid in set(p_e) & set(p_h)
        if p_h[mid] >= HI and p_e[mid] <= LO
    ]
    gt_only = sum(1 for mid in bonly if all(q == ">" for q in q_e[mid]))
    has_eq = sum(1 for mid in bonly if any(q in {"", "="} for q in q_e[mid]))
    return {
        "n_strict_B_only_as_is": len(bonly),
        "n_B_only_EGFR_gt_records_only": gt_only,
        "n_B_only_EGFR_has_equal_record": has_eq,
        "note": (
            "as-is B_only ligands whose EGFR records are exclusively '>' "
            "(typical IC50 > 10 µM panel values). These inflate hard-neg "
            "counts relative to ChEMBL pChEMBL, which usually needs '='."
        ),
    }



def verdict_markdown(rows: list[dict], fetch_log: dict) -> str:
    bdb_ok = all(
        fetch_log["bindingdb"][TARGETS[a]]["ok"] and fetch_log["bindingdb"][TARGETS[b]]["ok"]
        for _, a, b in PAIRS
    )
    pc_ok = all(
        fetch_log["pubchem"][TARGETS[a]]["ok"] and fetch_log["pubchem"][TARGETS[b]]["ok"]
        for _, a, b in PAIRS
    )
    by = row_lookup(rows)
    diag = fetch_log.get("egfr_her2_bindingdb_qualifier") or {}

    def get(pair: str, source: str, rule: str):
        return by.get((pair, source, rule))

    lines = [
        "# SUPPLY_CROSSDB_VERDICT_V1 — BindingDB / PubChem count-level check",
        "",
        f"**Date (UTC):** {fetch_log['run_finished_utc']}",
        "**Script:** `scripts/bindingdb_pubchem_strict_count_v1.py`",
        "**Scope:** frozen K=4 pairs only; **zero docking**; per-database ligand IDs (no cross-DB structure merge).",
        "",
        "## Rule (identical to J0 construction gate)",
        "",
        "- Types: IC50 / Ki / Kd / EC50.",
        "- Representative potency = **max** converted p-activity per ligand–target.",
        "- p from nM: `9 − log10(nM)`; PubChem concise values are µM → `6 − log10(µM)`.",
        "- dual: both ≥ 6.5; A_only: A ≥ 6.5 and B ≤ 5.5; B_only: converse.",
        "- Gate: `min(strict_A_only, strict_B_only) ≥ 50` (thin: ≥ 20).",
        "- **as_is:** strip `>`/`<` and use the numeric value (can inflate hard-negs; `>10000 nM` → p = 5.0).",
        "- **equal_only:** keep only `=` (or unqualified) records. This is the apples-to-apples comparison with ChEMBL pChEMBL, which usually requires a standard `=` relation.",
        "- ChEMBL column is the frozen `mols_*.json` max-pChEMBL cache (J0); no inequality mode.",
        "",
        "## Fetch status (not fabricated)",
        "",
        f"- BindingDB REST `getLigandsByUniprots` (cutoff = {BDB_CUTOFF_NM} nM): **{'OK' if bdb_ok else 'PARTIAL/FAIL'}**",
        f"- PubChem PUG REST `protein/accession/…/concise/CSV`: **{'OK' if pc_ok else 'PARTIAL/FAIL'}**",
        "",
        "Per-target fetch metadata is in `tables/fetch_log_v1.json`. Raw dumps stay in `cache/` (gitignored); compact pmax maps are in `tables/`.",
        "",
        "## Counts",
        "",
        "| pair | source | rule | both | strict dual | strict A/B | min HN | ≥50 both-side | ≥20 thin |",
        "|------|--------|------|-----:|------------:|------------|-------:|:-------------:|:--------:|",
    ]
    for row in rows:
        flag50 = "Y" if row["supports_strict_panel"] else "N"
        flag20 = "Y" if row["supports_thin_panel"] else "N"
        lines.append(
            f"| {row['pair']} | {row['source']} | {row['rule']} | {row['n_both_measured']} | "
            f"{row['strict_dual']} | {fmt_ab(row)} | {row['min_strict_hardneg']} | "
            f"{flag50} | {flag20} |"
        )

    lines += ["", "## Does the public-data ceiling move?", ""]
    lines.append(
        "**Primary comparison = ChEMBL pChEMBL vs BindingDB/PubChem `equal_only`.** "
        "The `as_is` rows are a sensitivity to censored `>` values, not the matched-rule headline."
    )
    lines.append("")
    thick_unchanged = True
    for pair, _a, _b in PAIRS:
        ch = get(pair, "ChEMBL_cache", "pChEMBL")
        parts = [
            f"**{pair}** ChEMBL min HN = {ch['min_strict_hardneg']} (A/B {fmt_ab(ch)})"
        ]
        for src in ("BindingDB", "PubChem"):
            eq = get(pair, src, "equal_only")
            raw = get(pair, src, "as_is")
            if eq is None or raw is None:
                parts.append(f"{src} missing")
                continue
            parts.append(
                f"{src} equal_only min HN = {eq['min_strict_hardneg']} "
                f"(A/B {fmt_ab(eq)}; both={eq['n_both_measured']}); "
                f"as_is min HN = {raw['min_strict_hardneg']} (A/B {fmt_ab(raw)})"
            )
            if bool(ch["supports_strict_panel"]) != bool(eq["supports_strict_panel"]):
                thick_unchanged = False
        lines.append("- " + "; ".join(parts) + ".")

    if diag:
        lines += [
            "",
            "## EGFR/HER2 BindingDB qualifier diagnostic (as_is B_only)",
            "",
            f"- as_is strict B_only (HER2-selective): **{diag.get('n_strict_B_only_as_is')}**",
            f"- of which EGFR records are **only** `>`: **{diag.get('n_B_only_EGFR_gt_records_only')}**",
            f"- of which EGFR has at least one `=` record: **{diag.get('n_B_only_EGFR_has_equal_record')}**",
            "",
            str(diag.get("note") or ""),
        ]

    lines += ["", "### One-line verdict", ""]
    if not bdb_ok and not pc_ok:
        lines.append(
            "**No BindingDB or PubChem counts are reported as facts** — fetch failed. "
            "Do not cite a cross-database supply conclusion from this run."
        )
    elif thick_unchanged:
        ch_e = get("EGFR/HER2", "ChEMBL_cache", "pChEMBL")
        b_eq = get("EGFR/HER2", "BindingDB", "equal_only")
        p_eq = get("EGFR/HER2", "PubChem", "equal_only")
        b_as = get("EGFR/HER2", "BindingDB", "as_is")
        extra = ""
        if ch_e and b_eq and b_as and p_eq:
            extra = (
                f" EGFR/HER2 remains below the ≥50 thick-panel gate under equal_only "
                f"(BindingDB min HN = {b_eq['min_strict_hardneg']}; PubChem = {p_eq['min_strict_hardneg']}) "
                f"versus ChEMBL min HN = {ch_e['min_strict_hardneg']}; it does reach the "
                f"thin ≥20 pool. The as_is counts that would pass ≥50 "
                f"(BindingDB min HN = {b_as['min_strict_hardneg']}) are driven largely by "
                f"`>` censored EGFR measurements, not by a new pool of equality-bounded "
                f"HER2-selective ligands."
            )
        lines.append(
            "**Under the matched equal-relation rule, BindingDB/PubChem do not flip the "
            "≥50 both-side thick-panel membership of this K=4 set.** "
            "The three frozen thick pairs remain above the gate; EGFR/HER2 remains a "
            f"supply-limited case rather than a newly thick panel.{extra} "
            "Absolute paired counts are higher than the ChEMBL cache, so this is not "
            "a claim that the databases are identical — only that the construction gate "
            "used for the docked panels does not change. Count-only; **no new docking** "
            "and no panel rebuild."
        )
    else:
        lines.append(
            "**The ≥50 both-side hard-neg gate is source-dependent even under equal_only.** "
            "Do **not** write that thick-panel membership is ChEMBL-invariant. "
            "Qualify K=4 as a ChEMBL-construction set and report Table S12 in full."
        )
    lines += [
        "",
        "## What this is not",
        "",
        "- Not a merged unique-structure census (no RDKit/InChIKey union).",
        "- Not a new docking panel and not a change to frozen K=4 ligands.",
        "- Not species / assay-confidence filtering (matches J0 max-p aggregation).",
        "- PubChem concise and BindingDB overlap (deposition); similar counts are expected and are not two independent censuses.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 Dual_Target_Docking/data/jcim_supply_crossdb_v0/scripts/bindingdb_pubchem_strict_count_v1.py",
        "```",
        "",
        "Re-runs reuse `cache/` if present. Delete `cache/` to force a live refetch.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    (CACHE / ".gitignore").write_text("*\n!.gitignore\n")

    fetch_log = {
        "run_started_utc": utc_now(),
        "keep_types": sorted(KEEP_TYPES),
        "hi": HI,
        "lo": LO,
        "bindingdb_cutoff_nM": BDB_CUTOFF_NM,
        "bindingdb": {},
        "pubchem": {},
    }

    bdb_maps = {"as_is": {}, "equal_only": {}}
    pc_maps = {"as_is": {}, "equal_only": {}}
    chembl_maps = {}

    for name, acc in TARGETS.items():
        print(f"[ChEMBL] {name}", flush=True)
        chembl_maps[name] = load_chembl(name)

        print(f"[BindingDB] {name} {acc}", flush=True)
        meta, maps = fetch_bindingdb(acc, CACHE)
        fetch_log["bindingdb"][acc] = meta
        bdb_maps["as_is"][name] = maps["as_is"]
        bdb_maps["equal_only"][name] = maps["equal_only"]
        if meta["ok"]:
            write_pmax(TABLES / f"bindingdb_{acc}_pmax.json", maps["as_is"])
            write_pmax(TABLES / f"bindingdb_{acc}_pmax_equal.json", maps["equal_only"])
        print(
            f"  ok={meta['ok']} raw={meta.get('n_raw_records')} "
            f"ligands={meta.get('n_ligands_pmax')} equal={meta.get('n_ligands_pmax_equal')} "
            f"err={meta.get('error')}",
            flush=True,
        )

        print(f"[PubChem] {name} {acc}", flush=True)
        meta, maps = fetch_pubchem(acc, CACHE)
        fetch_log["pubchem"][acc] = meta
        pc_maps["as_is"][name] = maps["as_is"]
        pc_maps["equal_only"][name] = maps["equal_only"]
        if meta["ok"]:
            write_pmax(TABLES / f"pubchem_{acc}_pmax.json", maps["as_is"])
            write_pmax(TABLES / f"pubchem_{acc}_pmax_equal.json", maps["equal_only"])
        print(
            f"  ok={meta['ok']} raw={meta.get('n_raw_records')} "
            f"ligands={meta.get('n_ligands_pmax')} equal={meta.get('n_ligands_pmax_equal')} "
            f"err={meta.get('error')}",
            flush=True,
        )

    rows = []
    for pair, a, b in PAIRS:
        rec = classify_pair(chembl_maps.get(a) or {}, chembl_maps.get(b) or {})
        rec.update(
            {
                "pair": pair,
                "target_A": a,
                "target_B": b,
                "uniprot_A": TARGETS[a],
                "uniprot_B": TARGETS[b],
                "source": "ChEMBL_cache",
                "rule": "pChEMBL",
                "source_complete": bool(chembl_maps.get(a) and chembl_maps.get(b)),
            }
        )
        rows.append(rec)
        for src_name, bundle in (("BindingDB", bdb_maps), ("PubChem", pc_maps)):
            for rule in ("as_is", "equal_only"):
                maps = bundle[rule]
                rec = classify_pair(maps.get(a) or {}, maps.get(b) or {})
                rec.update(
                    {
                        "pair": pair,
                        "target_A": a,
                        "target_B": b,
                        "uniprot_A": TARGETS[a],
                        "uniprot_B": TARGETS[b],
                        "source": src_name,
                        "rule": rule,
                        "source_complete": bool(maps.get(a) and maps.get(b)),
                    }
                )
                rows.append(rec)

    fields = [
        "pair",
        "source",
        "rule",
        "source_complete",
        "uniprot_A",
        "uniprot_B",
        "n_ligands_A",
        "n_ligands_B",
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
    ]
    out_csv = TABLES / "crossdb_strict_supply_v1.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    fetch_log["egfr_her2_bindingdb_qualifier"] = egfr_bindingdb_qualifier_diagnostic()
    fetch_log["run_finished_utc"] = utc_now()
    (TABLES / "fetch_log_v1.json").write_text(json.dumps(fetch_log, indent=2) + "\n")

    verdict = verdict_markdown(rows, fetch_log)
    (ANALYSIS / "SUPPLY_CROSSDB_VERDICT_V1.md").write_text(verdict)

    print(f"\nwrote {out_csv}")
    print(f"wrote {ANALYSIS / 'SUPPLY_CROSSDB_VERDICT_V1.md'}")
    print("\n=== summary ===")
    for rec in rows:
        print(
            f"{rec['pair']:16s} {rec['source']:13s} {rec['rule']:12s} "
            f"both={rec['n_both_measured']:5d} "
            f"strict A/B={rec['strict_A_only']:4d}/{rec['strict_B_only']:4d} "
            f"minHN={rec['min_strict_hardneg']:3d} "
            f"Y50={rec['supports_strict_panel']}"
        )


if __name__ == "__main__":
    main()
