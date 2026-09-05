#!/usr/bin/env python3
"""Count-only BindingDB + PubChem supply for the five post-census pairs.

Reuses the frozen K=4 fetch/classify rules from
jcim_supply_crossdb_v0/scripts/bindingdb_pubchem_strict_count_v1.py:
IC50/Ki/Kd/EC50, max converted p, as_is vs equal_only, no structure merge,
no docking, no panel rebuild.

ChEMBL comparison column is a fresh harvest from the same ChEMBL 37 dump
(STANDARD_OK + max pChEMBL), not the K=4 mols_*.json cache.

Does not overwrite jcim_supply_crossdb_v0 tables. Does not replace Table 2.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
OUT = ROOT / "local_track_b_v0" / "tables" / "five_pair_crossdb_v1"
AN = ROOT / "local_track_b_v0" / "analysis"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(ROOT.parent / "jcim_supply_crossdb_v0" / "scripts"))
from bindingdb_pubchem_strict_count_v1 import (  # noqa: E402
    CACHE as XDB_CACHE,
    classify_pair,
    fetch_bindingdb,
    fetch_pubchem,
    utc_now,
)
from pair_ligand_identity_qc_v1 import connect, harvest, resolve_targets  # noqa: E402

TARGETS = {
    "F2": "P00734",
    "F10": "P00742",
    "JAK1": "P23458",
    "TYK2": "P29597",
    "JAK2": "O60674",
    "PPARG": "P37231",
    "PPARA": "Q07869",
    "PPARD": "Q03181",
}
PAIRS = [
    ("F2/F10", "F2", "F10"),
    ("JAK1/TYK2", "JAK1", "TYK2"),
    ("JAK1/JAK2", "JAK1", "JAK2"),
    ("PPARG/PPARA", "PPARG", "PPARA"),
    ("PPARA/PPARD", "PPARA", "PPARD"),
]


def harvest_chembl_maps(sqlite: Path) -> dict[str, dict[int, float]]:
    con = connect(sqlite)
    meta = resolve_targets(con, set(TARGETS))
    missing = sorted(set(TARGETS) - set(meta))
    if missing:
        raise SystemExit(f"unresolved genes in dump: {missing}")
    maps = harvest(con, {m["tid"] for m in meta.values()})
    out = {}
    for gene, acc in TARGETS.items():
        if meta[gene]["uniprot"] != acc:
            raise SystemExit(
                f"UniProt mismatch {gene}: dump {meta[gene]['uniprot']} != lock {acc}"
            )
        out[gene] = maps.get(meta[gene]["tid"], {})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqlite", type=Path, required=True)
    args = ap.parse_args()
    if not args.sqlite.exists():
        raise SystemExit(f"missing sqlite: {args.sqlite}")
    OUT.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)
    XDB_CACHE.mkdir(parents=True, exist_ok=True)

    print("harvesting ChEMBL 37 maps…", flush=True)
    chembl_maps = harvest_chembl_maps(args.sqlite)

    fetch_log = {
        "run_started_utc": utc_now(),
        "scope": "five_post_census_pairs",
        "does_not_replace_table2": True,
        "bindingdb": {},
        "pubchem": {},
        "uniprot_check": TARGETS,
    }
    bdb_maps = {"as_is": {}, "equal_only": {}}
    pc_maps = {"as_is": {}, "equal_only": {}}

    for name, acc in TARGETS.items():
        print(f"[BindingDB] {name} {acc}", flush=True)
        meta, maps = fetch_bindingdb(acc, XDB_CACHE)
        fetch_log["bindingdb"][acc] = {k: v for k, v in meta.items() if k != "url"} | {
            "url": meta["url"]
        }
        # rewrite cache_path relative to this cache
        bdb_maps["as_is"][name] = maps["as_is"]
        bdb_maps["equal_only"][name] = maps["equal_only"]
        print(
            f"  ok={meta['ok']} ligands={meta.get('n_ligands_pmax')} "
            f"equal={meta.get('n_ligands_pmax_equal')} err={meta.get('error')}",
            flush=True,
        )
        print(f"[PubChem] {name} {acc}", flush=True)
        meta, maps = fetch_pubchem(acc, XDB_CACHE)
        fetch_log["pubchem"][acc] = meta
        pc_maps["as_is"][name] = maps["as_is"]
        pc_maps["equal_only"][name] = maps["equal_only"]
        print(
            f"  ok={meta['ok']} ligands={meta.get('n_ligands_pmax')} "
            f"equal={meta.get('n_ligands_pmax_equal')} err={meta.get('error')}",
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
                "source": "ChEMBL37_dump",
                "rule": "pChEMBL_STANDARD_OK",
                "source_complete": True,
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
    out_csv = OUT / "crossdb_strict_supply_v1.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    fetch_log["run_finished_utc"] = utc_now()
    (OUT / "fetch_log_v1.json").write_text(json.dumps(fetch_log, indent=2) + "\n")

    lines = [
        "# Five-pair BindingDB / PubChem count-only\n\n",
        "Same fetch rules as K=4 (`bindingdb_pubchem_strict_count_v1.py`). ",
        "Per-database IDs; no InChIKey merge; **no docking**; not Table 2.\n\n",
        "Primary comparison is **equal_only**. `as_is` can inflate hard-negatives.\n\n",
        "| pair | ChEMBL both / min HN | BindingDB equal both / min HN | PubChem equal both / min HN |\n",
        "|---|---:|---:|---:|\n",
    ]
    by = {(r["pair"], r["source"], r["rule"]): r for r in rows}
    for pair, _a, _b in PAIRS:
        ch = by[(pair, "ChEMBL37_dump", "pChEMBL_STANDARD_OK")]
        bd = by.get((pair, "BindingDB", "equal_only"), {})
        pc = by.get((pair, "PubChem", "equal_only"), {})
        lines.append(
            f"| {pair} | {ch['n_both_measured']} / {ch['min_strict_hardneg']} | "
            f"{bd.get('n_both_measured', '')} / {bd.get('min_strict_hardneg', '')} | "
            f"{pc.get('n_both_measured', '')} / {pc.get('min_strict_hardneg', '')} |\n"
        )
    lines.append(
        "\nDo not hard-dock BindingDB as external validation. "
        "Do not treat a thicker BindingDB `as_is` hard-neg count as a new panel.\n"
    )
    (AN / "FIVE_PAIR_CROSSDB_V1.md").write_text("".join(lines), encoding="utf-8")
    print("wrote", out_csv)
    print("".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
