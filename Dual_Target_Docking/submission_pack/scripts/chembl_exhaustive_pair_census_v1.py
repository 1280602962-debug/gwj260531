#!/usr/bin/env python3
"""Exhaustive ChEMBL human SINGLE PROTEIN dual-target pair census.

Uses a local ChEMBL SQLite dump (default: ChEMBL 37). Does not call the
Web API. Does not dock. Does not rewrite Table 2 / Table S44 / K = 4.

Label rules match Dual_Target_Docking/data/jcim_j0j1_v0/scripts/run_j0_supply_audit.py.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

THETA = 6.0
HI = 6.5
LO = 5.5
MIN_HARDNEG_STRICT = 50
MIN_HARDNEG_THIN = 20
MIN_BOTH_WRITE = 10
STANDARD_OK = ("IC50", "Ki", "Kd", "EC50", "Potency", "IC50app", "Ki app")
METAL_RE = re.compile(
    r"histone deacetylase|\bhdac\b|carbonic anhydrase|matrix metallo",
    re.I,
)
PPI_RE = re.compile(r"\bbcl-2\b|\bbcl2\b|bh3|bromodomain", re.I)

# UniProt for fetch-queue names that are not in chembl_target_ids.json
NAMED_UNIPROT = {
    "VEGFR2_KDR": "P35968",
    "AXL": "P30530",
    "MERTK": "Q12866",
    "SYK": "P43405",
    "HSP90AA1": "P07900",
    "WEE1": "P30291",
    "TOP1": "P11387",
    "ROCK1": "Q13464",
    "PIM1": "P11309",
    "SERT_SLC6A4": "P31645",
    "ESR1": "P03372",
    "FGFR1": "P11362",
    "ALK": "Q9UM73",
    "BRAF": "P15056",
    "MAP2K1_MEK1": "Q02750",
    "BTK": "Q06187",
    "FLT3": "P36888",
    "SRC": "P12931",
    "ABL1": "P00519",
    "BACE1": "P56817",
}

# Literature pairs the fetch queue was meant to unlock (J0 EXTRA_TARGETS_FOR_QUEUE).
FETCH_QUEUE_INTENDED_PAIRS = [
    ("VEGFR2_KDR", "HDAC1", "VEGFR2/HDAC duals (JMC)"),
    ("VEGFR2_KDR", "HDAC6", "VEGFR2/HDAC duals"),
    ("AXL", "MERTK", "MER/AXL duals"),
    ("SYK", "HDAC1", "SYK/HDAC duals"),
    ("HSP90AA1", "HDAC6", "Hsp90/HDAC6 duals"),
    ("WEE1", "HDAC1", "Wee1/HDAC duals"),
    ("TOP1", "HDAC1", "Top/HDAC duals"),
    ("ROCK1", "HDAC1", "ROCK/HDAC duals"),
    ("PIM1", "HDAC1", "PIM/HDAC duals"),
    ("SERT_SLC6A4", "ESR1", "SERT/ER duals"),
    ("BRAF", "MAP2K1_MEK1", "BRAF/MEK"),
    ("FGFR1", "VEGFR2_KDR", "FGFR dual TKIs"),
    ("ALK", "EGFR", "ALK dual TKIs"),
    ("BTK", "EGFR", "BTK duals"),
    ("FLT3", "VEGFR2_KDR", "FLT3 dual TKIs"),
    ("SRC", "ABL1", "BCR-ABL dual context"),
    ("BACE1", "ACHE", "CNS duals"),
]

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = Path(__file__).resolve().parents[1]
J0_SCRIPT_HINT = ROOT / "data" / "jcim_j0j1_v0" / "scripts" / "run_j0_supply_audit.py"


def sha256_file(path: Path, buf: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(buf)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def connect(sqlite_path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(str(sqlite_path))
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    con.execute("PRAGMA temp_store = MEMORY")
    con.execute("PRAGMA cache_size = -250000")
    return con


def load_single_component_human_targets(con: sqlite3.Connection) -> list[dict]:
    sql = """
    SELECT
      td.tid,
      td.chembl_id AS target_chembl,
      td.pref_name,
      td.organism,
      cs.accession AS uniprot,
      cs.component_id,
      cs.description AS component_desc
    FROM target_dictionary td
    JOIN target_components tc ON tc.tid = td.tid
    JOIN component_sequences cs ON cs.component_id = tc.component_id
    WHERE td.target_type = 'SINGLE PROTEIN'
      AND td.organism = 'Homo sapiens'
      AND td.tid IN (
        SELECT tid FROM target_components GROUP BY tid HAVING COUNT(*) = 1
      )
    """
    rows = [dict(r) for r in con.execute(sql)]
    symbols: dict[int, str] = {}
    for r in con.execute(
        """
        SELECT component_id, component_synonym
        FROM component_synonyms
        WHERE syn_type = 'GENE_SYMBOL'
        """
    ):
        cid, syn = r[0], r[1]
        if cid not in symbols and syn:
            symbols[cid] = syn
    classes: dict[int, tuple[int, str, str]] = {}
    for r in con.execute(
        """
        SELECT cc.component_id, pc.class_level, pc.pref_name, pc.protein_class_desc
        FROM component_class cc
        JOIN protein_classification pc ON pc.protein_class_id = cc.protein_class_id
        """
    ):
        cid, level, pref, desc = r
        prev = classes.get(cid)
        if prev is None or (level or 0) > prev[0]:
            classes[cid] = (level or 0, pref or "", desc or "")
    out = []
    for rec in rows:
        cid = rec["component_id"]
        rec["gene_symbol"] = symbols.get(cid, "")
        lvl, pref, desc = classes.get(cid, (0, "", ""))
        rec["protein_class"] = pref
        rec["protein_class_desc"] = desc
        blob = " ".join(
            str(rec.get(k) or "")
            for k in ("pref_name", "component_desc", "protein_class", "protein_class_desc", "gene_symbol")
        )
        rec["metal_enzyme_risk"] = bool(METAL_RE.search(blob))
        rec["ppi_like"] = bool(PPI_RE.search(blob))
        out.append(rec)
    return out


def harvest_maps(
    con: sqlite3.Connection, tids: set[int], confidence_min: int | None
) -> dict[int, dict[int, float]]:
    placeholders = ",".join("?" * len(STANDARD_OK))
    conf_clause = ""
    params: list = list(STANDARD_OK)
    if confidence_min is not None:
        conf_clause = "AND ass.confidence_score >= ?"
        params.append(confidence_min)
    sql = f"""
    SELECT ass.tid, act.molregno, MAX(act.pchembl_value) AS pchembl
    FROM activities act
    JOIN assays ass ON ass.assay_id = act.assay_id
    JOIN target_dictionary td ON td.tid = ass.tid
    WHERE td.target_type = 'SINGLE PROTEIN'
      AND td.organism = 'Homo sapiens'
      AND act.pchembl_value IS NOT NULL
      AND act.standard_type IN ({placeholders})
      {conf_clause}
    GROUP BY ass.tid, act.molregno
    """
    maps: dict[int, dict[int, float]] = defaultdict(dict)
    n = 0
    t0 = time.time()
    for tid, mol, pchembl in con.execute(sql, params):
        if tid not in tids:
            continue
        maps[tid][int(mol)] = float(pchembl)
        n += 1
        if n % 500_000 == 0:
            print(f"  harvested {n:,} target-mol rows in {time.time()-t0:.1f}s", flush=True)
    return dict(maps)


def tri_index(i: int, j: int, n: int) -> int:
    if i > j:
        i, j = j, i
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def census_pairs(tids: list[int], maps: dict[int, dict[int, float]]):
    n = len(tids)
    n_tri = n * (n - 1) // 2
    n_both = [0] * n_tri
    th_d = [0] * n_tri
    th_a = [0] * n_tri
    th_b = [0] * n_tri
    st_d = [0] * n_tri
    st_a = [0] * n_tri
    st_b = [0] * n_tri
    st_n = [0] * n_tri
    st_g = [0] * n_tri
    inv: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for tidx, tid in enumerate(tids):
        for mol, pv in maps.get(tid, {}).items():
            inv[mol].append((tidx, pv))
    n_mols = 0
    t0 = time.time()
    for mol, hits in inv.items():
        n_mols += 1
        if n_mols % 200_000 == 0:
            print(f"  pair-accumulate mols={n_mols:,} elapsed={time.time()-t0:.1f}s", flush=True)
        if len(hits) < 2:
            continue
        hits.sort(key=lambda x: x[0])
        for a in range(len(hits)):
            ia, xa = hits[a]
            for b in range(a + 1, len(hits)):
                ib, xb = hits[b]
                k = tri_index(ia, ib, n)
                n_both[k] += 1
                if xa >= THETA and xb >= THETA:
                    th_d[k] += 1
                elif xa >= THETA:
                    th_a[k] += 1
                elif xb >= THETA:
                    th_b[k] += 1
                if xa >= HI and xb >= HI:
                    st_d[k] += 1
                elif xa >= HI and xb <= LO:
                    st_a[k] += 1
                elif xb >= HI and xa <= LO:
                    st_b[k] += 1
                elif xa <= LO and xb <= LO:
                    st_n[k] += 1
                else:
                    st_g[k] += 1
    return {
        "n_both": n_both,
        "theta_dual": th_d,
        "theta_A_only": th_a,
        "theta_B_only": th_b,
        "strict_dual": st_d,
        "strict_A_only": st_a,
        "strict_B_only": st_b,
        "strict_neither": st_n,
        "strict_gray": st_g,
    }


def pair_record(i: int, j: int, tids: list[int], meta: dict[int, dict], acc: dict, k: int) -> dict:
    ia, ib = (i, j) if i < j else (j, i)
    ta, tb = tids[ia], tids[ib]
    ma, mb = meta[ta], meta[tb]
    n_both = acc["n_both"][k]
    n_d = acc["theta_dual"][k]
    n_a = acc["theta_A_only"][k]
    n_b = acc["theta_B_only"][k]
    n_neither = n_both - n_d - n_a - n_b
    st_a = acc["strict_A_only"][k]
    st_b = acc["strict_B_only"][k]
    min_hn = min(st_a, st_b)
    gene_a = ma["gene_symbol"] or ma["pref_name"]
    gene_b = mb["gene_symbol"] or mb["pref_name"]
    return {
        "target_A_chembl": ma["target_chembl"],
        "target_B_chembl": mb["target_chembl"],
        "uniprot_A": ma["uniprot"] or "",
        "uniprot_B": mb["uniprot"] or "",
        "gene_A": gene_a,
        "gene_B": gene_b,
        "pref_A": ma["pref_name"],
        "pref_B": mb["pref_name"],
        "class_A": ma["protein_class"],
        "class_B": mb["protein_class"],
        "n_A": ma["n_mols"],
        "n_B": mb["n_mols"],
        "n_both_measured": n_both,
        "theta_dual": n_d,
        "theta_A_only": n_a,
        "theta_B_only": n_b,
        "theta_neither": n_neither,
        "strict_dual": acc["strict_dual"][k],
        "strict_A_only": st_a,
        "strict_B_only": st_b,
        "strict_neither": acc["strict_neither"][k],
        "strict_gray": acc["strict_gray"][k],
        "min_strict_hardneg": min_hn,
        "directional_n10": int(min(n_d, n_a, n_b) >= 10),
        "formulation_n10": int(min(n_d, n_a, n_b) >= 10 and n_neither >= 10),
        "supports_strict_panel": int(min_hn >= MIN_HARDNEG_STRICT),
        "supports_thin_panel": int(min_hn >= MIN_HARDNEG_THIN),
        "metal_either": int(ma["metal_enzyme_risk"] or mb["metal_enzyme_risk"]),
        "metal_both": int(ma["metal_enzyme_risk"] and mb["metal_enzyme_risk"]),
        "same_class": int(bool(ma["protein_class"]) and ma["protein_class"] == mb["protein_class"]),
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for rec in rows:
            w.writerow({k: rec.get(k, "") for k in fields})


def load_j0_named(root: Path) -> tuple[list[dict], list[dict]]:
    ids_path = root / "data" / "public_pair_selection" / "chembl_target_ids.json"
    named = json.loads(ids_path.read_text())
    cand_path = root / "data" / "jcim_j0j1_v0" / "tables" / "j0_candidate_pairs.csv"
    cands = list(csv.DictReader(cand_path.open()))
    return named, cands


def resolve_named(named_json: list[dict], meta_by_chembl: dict, meta_by_uniprot: dict) -> dict[str, dict | None]:
    out: dict[str, dict | None] = {}
    for rec in named_json:
        hit = meta_by_chembl.get(rec["chembl"]) or meta_by_uniprot.get(rec["uniprot"])
        out[rec["name"]] = hit
    for name, acc in NAMED_UNIPROT.items():
        if name not in out:
            out[name] = meta_by_uniprot.get(acc)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqlite", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--archive-sha256", default="")
    args = ap.parse_args()
    sqlite_path = args.sqlite.resolve()
    if not sqlite_path.exists():
        print(f"missing sqlite: {sqlite_path}", file=sys.stderr)
        return 2
    print(f"sqlite={sqlite_path}", flush=True)
    print(f"sqlite_sha256={sha256_file(sqlite_path)}", flush=True)
    if args.archive_sha256:
        print(f"archive_sha256_claimed={args.archive_sha256}", flush=True)

    con = connect(sqlite_path)
    try:
        version_row = con.execute("SELECT name, creation_date FROM version").fetchone()
        print(f"version={dict(version_row) if version_row else 'unknown'}", flush=True)
    except sqlite3.Error as exc:
        print(f"version lookup skipped: {exc}", flush=True)

    targets = load_single_component_human_targets(con)
    print(f"human SINGLE PROTEIN (1 component): {len(targets)}", flush=True)
    tids = {t["tid"] for t in targets}
    print("harvesting max-pChEMBL maps (no confidence cut)...", flush=True)
    maps_all = harvest_maps(con, tids, None)
    print("harvesting max-pChEMBL maps (confidence ≥ 8)...", flush=True)
    maps_hc = harvest_maps(con, tids, 8)

    def attach_n(maps):
        meta = {}
        kept = []
        for t in targets:
            n = len(maps.get(t["tid"], {}))
            rec = dict(t)
            rec["n_mols"] = n
            meta[t["tid"]] = rec
            if n > 0:
                kept.append(t["tid"])
        kept.sort()
        return meta, kept

    def run_slice(tag: str, maps: dict[int, dict[int, float]], out_dir: Path) -> dict:
        meta, kept = attach_n(maps)
        print(f"[{tag}] targets with ≥1 mol: {len(kept)}", flush=True)
        print(f"[{tag}] accumulating unordered pairs...", flush=True)
        acc = census_pairs(kept, maps)
        n = len(kept)
        pair_fields = [
            "target_A_chembl",
            "target_B_chembl",
            "uniprot_A",
            "uniprot_B",
            "gene_A",
            "gene_B",
            "pref_A",
            "pref_B",
            "class_A",
            "class_B",
            "n_A",
            "n_B",
            "n_both_measured",
            "theta_dual",
            "theta_A_only",
            "theta_B_only",
            "theta_neither",
            "strict_dual",
            "strict_A_only",
            "strict_B_only",
            "strict_neither",
            "strict_gray",
            "min_strict_hardneg",
            "directional_n10",
            "formulation_n10",
            "supports_strict_panel",
            "supports_thin_panel",
            "metal_either",
            "metal_both",
            "same_class",
        ]
        ge10 = []
        directional = []
        strict_y = []
        n_both_pos = 0
        n_dir = n_form = n_y = n_y_nometal = n_thin = 0
        for i in range(n):
            for j in range(i + 1, n):
                k = tri_index(i, j, n)
                if acc["n_both"][k] <= 0:
                    continue
                n_both_pos += 1
                if acc["n_both"][k] < MIN_BOTH_WRITE:
                    continue
                rec = pair_record(i, j, kept, meta, acc, k)
                ge10.append(rec)
                n_dir += rec["directional_n10"]
                n_form += rec["formulation_n10"]
                n_y += rec["supports_strict_panel"]
                n_thin += rec["supports_thin_panel"]
                if rec["supports_strict_panel"] and not rec["metal_either"]:
                    n_y_nometal += 1
                if rec["directional_n10"]:
                    directional.append(rec)
                if rec["supports_strict_panel"]:
                    strict_y.append(rec)
        ge10.sort(key=lambda r: (-r["min_strict_hardneg"], -r["n_both_measured"]))
        directional.sort(key=lambda r: (-r["min_strict_hardneg"], -r["n_both_measured"]))
        strict_y.sort(key=lambda r: (-r["min_strict_hardneg"], -r["n_both_measured"]))
        write_csv(out_dir / f"universe_pairs_n_both_ge10_{tag}.csv", ge10, pair_fields)
        write_csv(out_dir / f"universe_pairs_directional_n10_{tag}.csv", directional, pair_fields)
        write_csv(out_dir / f"universe_pairs_strict_thick_{tag}.csv", strict_y, pair_fields)
        tgt_rows = []
        for tid in kept:
            m = meta[tid]
            tgt_rows.append(
                {
                    "target_chembl": m["target_chembl"],
                    "uniprot": m["uniprot"] or "",
                    "gene_symbol": m["gene_symbol"],
                    "pref_name": m["pref_name"],
                    "protein_class": m["protein_class"],
                    "n_mols": m["n_mols"],
                    "metal_enzyme_risk": int(m["metal_enzyme_risk"]),
                }
            )
        tgt_rows.sort(key=lambda r: -r["n_mols"])
        write_csv(
            out_dir / f"universe_targets_{tag}.csv",
            tgt_rows,
            [
                "target_chembl",
                "uniprot",
                "gene_symbol",
                "pref_name",
                "protein_class",
                "n_mols",
                "metal_enzyme_risk",
            ],
        )
        summary = {
            "slice": tag,
            "n_human_single_protein_targets": len(targets),
            "n_targets_with_map": len(kept),
            "n_unordered_pairs_possible": n * (n - 1) // 2,
            "n_pairs_n_both_ge_1": n_both_pos,
            "n_pairs_n_both_ge_10": len(ge10),
            "n_directional_n10": n_dir,
            "n_formulation_n10": n_form,
            "n_strict_thick": n_y,
            "n_strict_thick_nonmetal": n_y_nometal,
            "n_thin_or_thick": n_thin,
        }
        return {"summary": summary, "meta": meta, "kept": kept, "acc": acc}

    tables = args.out / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    primary = run_slice("all", maps_all, tables)
    hc = run_slice("conf8", maps_hc, tables)
    write_csv(
        tables / "universe_census_summary_v1.csv",
        [primary["summary"], hc["summary"]],
        list(primary["summary"].keys()),
    )

    named_json, j0_cands = load_j0_named(args.root)
    meta_by_chembl = {m["target_chembl"]: m for m in primary["meta"].values()}
    meta_by_uniprot = {m["uniprot"]: m for m in primary["meta"].values() if m.get("uniprot")}
    resolved = resolve_named(named_json, meta_by_chembl, meta_by_uniprot)
    tid_index = {tid: i for i, tid in enumerate(primary["kept"])}

    def lookup_pair(name_a: str, name_b: str) -> dict:
        ma, mb = resolved.get(name_a), resolved.get(name_b)
        rec = {
            "name_A": name_a,
            "name_B": name_b,
            "resolved_A": int(ma is not None),
            "resolved_B": int(mb is not None),
            "chembl_A": (ma or {}).get("target_chembl", ""),
            "chembl_B": (mb or {}).get("target_chembl", ""),
            "uniprot_A": (ma or {}).get("uniprot", ""),
            "uniprot_B": (mb or {}).get("uniprot", ""),
            "n_A": (ma or {}).get("n_mols", 0),
            "n_B": (mb or {}).get("n_mols", 0),
        }
        if ma is None or mb is None:
            rec["in_universe"] = 0
            return rec
        ia = tid_index.get(ma["tid"])
        ib = tid_index.get(mb["tid"])
        if ia is None or ib is None or ia == ib:
            rec["in_universe"] = 0
            return rec
        k = tri_index(ia, ib, len(primary["kept"]))
        full = pair_record(ia, ib, primary["kept"], primary["meta"], primary["acc"], k)
        rec["in_universe"] = 1
        rec.update({f"u_{k}": full[k] for k in (
            "n_both_measured",
            "theta_dual",
            "theta_A_only",
            "theta_B_only",
            "theta_neither",
            "strict_A_only",
            "strict_B_only",
            "min_strict_hardneg",
            "directional_n10",
            "formulation_n10",
            "supports_strict_panel",
            "metal_either",
        )})
        return rec

    xwalk = []
    for cand in j0_cands:
        rec = lookup_pair(cand["target_A"], cand["target_B"])
        rec["pair_id"] = cand["pair_id"]
        rec["j0_notes"] = cand.get("notes", "")
        rec["j0_auditable_now"] = cand.get("auditable_now", "")
        xwalk.append(rec)
    write_csv(
        tables / "j0_universe_crosswalk_v1.csv",
        xwalk,
        [
            "pair_id",
            "name_A",
            "name_B",
            "j0_notes",
            "j0_auditable_now",
            "resolved_A",
            "resolved_B",
            "in_universe",
            "chembl_A",
            "chembl_B",
            "uniprot_A",
            "uniprot_B",
            "n_A",
            "n_B",
            "u_n_both_measured",
            "u_theta_dual",
            "u_theta_A_only",
            "u_theta_B_only",
            "u_theta_neither",
            "u_strict_A_only",
            "u_strict_B_only",
            "u_min_strict_hardneg",
            "u_directional_n10",
            "u_formulation_n10",
            "u_supports_strict_panel",
            "u_metal_either",
        ],
    )
    fetch_rows = []
    for name, acc in NAMED_UNIPROT.items():
        hit = resolved.get(name)
        fetch_rows.append(
            {
                "name": name,
                "uniprot": acc,
                "resolved": int(hit is not None),
                "target_chembl": (hit or {}).get("target_chembl", ""),
                "n_mols": (hit or {}).get("n_mols", 0),
                "pref_name": (hit or {}).get("pref_name", ""),
            }
        )
    write_csv(
        tables / "fetch_queue_universe_targets_v1.csv",
        fetch_rows,
        ["name", "uniprot", "resolved", "target_chembl", "n_mols", "pref_name"],
    )

    intended_rows = []
    for name_a, name_b, reason in FETCH_QUEUE_INTENDED_PAIRS:
        rec = lookup_pair(name_a, name_b)
        ma, mb = resolved.get(name_a), resolved.get(name_b)
        n_both = int(rec.get("u_n_both_measured") or 0)
        intended_rows.append(
            {
                "name_A": name_a,
                "name_B": name_b,
                "gene_A": (ma or {}).get("gene_symbol", ""),
                "gene_B": (mb or {}).get("gene_symbol", ""),
                "reason": reason,
                "found_n_both_ge10": int(n_both >= 10),
                "n_both_measured": n_both,
                "theta_dual": rec.get("u_theta_dual", 0),
                "theta_A_only": rec.get("u_theta_A_only", 0),
                "theta_B_only": rec.get("u_theta_B_only", 0),
                "theta_neither": rec.get("u_theta_neither", 0),
                "min_strict_hardneg": rec.get("u_min_strict_hardneg", 0),
                "directional_n10": rec.get("u_directional_n10", 0),
                "supports_strict_panel": rec.get("u_supports_strict_panel", 0),
            }
        )
    write_csv(
        tables / "fetch_queue_intended_pairs_v1.csv",
        intended_rows,
        [
            "name_A",
            "name_B",
            "gene_A",
            "gene_B",
            "reason",
            "found_n_both_ge10",
            "n_both_measured",
            "theta_dual",
            "theta_A_only",
            "theta_B_only",
            "theta_neither",
            "min_strict_hardneg",
            "directional_n10",
            "supports_strict_panel",
        ],
    )

    print("SUMMARY primary (no confidence cut):", primary["summary"], flush=True)
    print("SUMMARY conf8:", hc["summary"], flush=True)
    print(f"wrote tables under {tables}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
