#!/usr/bin/env python3
"""Extract co-tested four-state SMILES panels for the five Track B pairs.

Does not dock. Does not include CTSK/CTSS. Does not rebuild frozen K=4 panels.

Rules match DOCKING_PLAN_V1.md and pair_ligand_identity_qc_v1.py:
  - co-tested only (set(map_A) ∩ set(map_B))
  - max pChEMBL, same endpoints as the census
  - main labels: strict 6.5 / 5.5 (gray excluded)
  - small-molecule filter (MOL, Small molecule, MW 150–750, heavy 10–60, no metal)
  - class-quota + deterministic shuffle (seed 20260729)
  - neither arm only if strict_neither after the small-molecule filter ≥ 20
  - target depth 110 ligands
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
TABLES = Path(__file__).resolve().parents[1] / "tables"
sys.path.insert(0, str(SCRIPTS))

from pair_ligand_identity_qc_v1 import (  # noqa: E402
    HI,
    LO,
    THETA,
    classify,
    connect,
    harvest,
    mol_properties,
    resolve_targets,
)

SEED = 20260729
NEITHER_MIN = 20
QUOTA_WITH_NEITHER = {"dual": 32, "A_only": 32, "B_only": 32, "neither": 14}
QUOTA_NO_NEITHER = {"dual": 37, "A_only": 37, "B_only": 36}

PAIRS = [
    ("F2", "F10", "F2F10"),
    ("JAK1", "TYK2", "J1TYK2"),
    ("JAK1", "JAK2", "J1J2"),
    ("PPARG", "PPARA", "PGPA"),
    ("PPARA", "PPARD", "PAPD"),
]


def theta6_class(pa: float, pb: float) -> str:
    if pa >= THETA and pb >= THETA:
        return "dual"
    if pa >= THETA and pb < THETA:
        return "A_only"
    if pb >= THETA and pa < THETA:
        return "B_only"
    return "neither"


def strict_class(pa: float, pb: float) -> str | None:
    if pa >= HI and pb >= HI:
        return "dual"
    if pa >= HI and pb <= LO:
        return "A_only"
    if pb >= HI and pa <= LO:
        return "B_only"
    if pa <= LO and pb <= LO:
        return "neither"
    return None  # gray


def sample_class(pool: list[dict], need: int, rng: random.Random) -> list[dict]:
    shuffled = list(pool)
    rng.shuffle(shuffled)
    return shuffled[:need]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqlite", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=TABLES / "track_b_panels")
    args = ap.parse_args()
    if not args.sqlite.exists():
        print(f"missing sqlite: {args.sqlite}", file=sys.stderr)
        return 2

    genes = {g for a, b, _ in PAIRS for g in (a, b)}
    con = connect(args.sqlite)
    meta = resolve_targets(con, genes)
    missing = sorted(genes - set(meta))
    if missing:
        print(f"unresolved genes: {missing}", file=sys.stderr)
        return 2
    print(f"resolved {len(meta)} targets; harvesting maps...", flush=True)
    maps = harvest(con, {m["tid"] for m in meta.values()})

    pair_rows: dict[str, dict[str, list[dict]]] = {}
    all_mols: set[int] = set()
    theta_counts: dict[str, dict[str, int]] = {}
    for a, b, prefix in PAIRS:
        pair = f"{a}/{b}"
        ma = maps.get(meta[a]["tid"], {})
        mb = maps.get(meta[b]["tid"], {})
        both = set(ma) & set(mb)
        buckets: dict[str, list[int]] = defaultdict(list)
        t6 = defaultdict(int)
        for mol in both:
            pa, pb = ma[mol], mb[mol]
            t6[theta6_class(pa, pb)] += 1
            cls = strict_class(pa, pb)
            if cls:
                buckets[cls].append(mol)
                all_mols.add(mol)
        pair_rows[pair] = buckets
        theta_counts[pair] = dict(t6)
        print(
            f"{pair:16s} n_both={len(both):5d}  strict "
            f"D/A/B/N={len(buckets['dual'])}/{len(buckets['A_only'])}/"
            f"{len(buckets['B_only'])}/{len(buckets['neither'])}  "
            f"gray={len(both) - sum(len(v) for v in buckets.values())}",
            flush=True,
        )

    print(f"loading properties for {len(all_mols):,} co-tested molecules...", flush=True)
    props = mol_properties(con, all_mols)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    panel_fields = [
        "panel_id",
        "pair",
        "class",
        "molecule_chembl_id",
        "canonical_smiles",
        "pchembl_A",
        "pchembl_B",
        "theta6_class",
        "label_rule",
        "mw_freebase",
        "heavy_atoms",
        "gene_A",
        "gene_B",
        "uniprot_A",
        "uniprot_B",
    ]

    for a, b, prefix in PAIRS:
        pair = f"{a}/{b}"
        ma = maps[meta[a]["tid"]]
        mb = maps[meta[b]["tid"]]
        pools: dict[str, list[dict]] = {k: [] for k in ("dual", "A_only", "B_only", "neither")}
        reject = defaultdict(int)
        for cls, mols in pair_rows[pair].items():
            for mol in mols:
                p = props.get(mol) or {}
                good, why = classify(p)
                if not good:
                    reject[why] += 1
                    continue
                pa, pb = ma[mol], mb[mol]
                pools[cls].append(
                    {
                        "molecule_chembl_id": p["chembl_id"],
                        "canonical_smiles": p["canonical_smiles"],
                        "pchembl_A": pa,
                        "pchembl_B": pb,
                        "theta6_class": theta6_class(pa, pb),
                        "mw_freebase": p.get("mw_freebase"),
                        "heavy_atoms": p.get("heavy_atoms"),
                        "class": cls,
                    }
                )

        n_neither_sm = len(pools["neither"])
        include_neither = n_neither_sm >= NEITHER_MIN
        quota = QUOTA_WITH_NEITHER if include_neither else QUOTA_NO_NEITHER
        rng = random.Random(SEED)
        picked: list[dict] = []
        got = {}
        for cls, need in quota.items():
            take = sample_class(pools[cls], need, rng)
            got[cls] = len(take)
            if len(take) < need:
                print(
                    f"WARN {pair} class {cls}: got {len(take)}/{need} after small-mol filter",
                    file=sys.stderr,
                )
            picked.extend(take)

        picked.sort(key=lambda r: (r["class"], r["molecule_chembl_id"]))
        out_rows = []
        for i, r in enumerate(picked, 1):
            out_rows.append(
                {
                    "panel_id": f"{prefix}_{i:03d}",
                    "pair": pair,
                    "class": r["class"],
                    "molecule_chembl_id": r["molecule_chembl_id"],
                    "canonical_smiles": r["canonical_smiles"],
                    "pchembl_A": f"{r['pchembl_A']:.2f}",
                    "pchembl_B": f"{r['pchembl_B']:.2f}",
                    "theta6_class": r["theta6_class"],
                    "label_rule": "strict_6.5_5.5",
                    "mw_freebase": r["mw_freebase"],
                    "heavy_atoms": r["heavy_atoms"],
                    "gene_A": a,
                    "gene_B": b,
                    "uniprot_A": meta[a]["uniprot"],
                    "uniprot_B": meta[b]["uniprot"],
                }
            )
        out_csv = args.out_dir / f"panel_{a}_{b}_v1.csv"
        with out_csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=panel_fields)
            w.writeheader()
            w.writerows(out_rows)

        min_hn = min(len(pools["A_only"]), len(pools["B_only"]))
        rec = {
            "pair": pair,
            "n_both_theta_dual": theta_counts[pair].get("dual", 0),
            "n_both_theta_A_only": theta_counts[pair].get("A_only", 0),
            "n_both_theta_B_only": theta_counts[pair].get("B_only", 0),
            "n_both_theta_neither": theta_counts[pair].get("neither", 0),
            "n_strict_dual_smallmol": len(pools["dual"]),
            "n_strict_A_only_smallmol": len(pools["A_only"]),
            "n_strict_B_only_smallmol": len(pools["B_only"]),
            "n_strict_neither_smallmol": n_neither_sm,
            "min_hardneg_smallmol": min_hn,
            "include_neither": int(include_neither),
            "n_panel": len(out_rows),
            "n_panel_dual": got.get("dual", 0),
            "n_panel_A_only": got.get("A_only", 0),
            "n_panel_B_only": got.get("B_only", 0),
            "n_panel_neither": got.get("neither", 0),
            "seed": SEED,
            "quota": "+".join(f"{k}:{v}" for k, v in quota.items()),
            "csv": str(out_csv.name),
        }
        summary.append(rec)
        print(
            f"wrote {out_csv.name} n={len(out_rows)} "
            f"D/A/B/N={got.get('dual', 0)}/{got.get('A_only', 0)}/"
            f"{got.get('B_only', 0)}/{got.get('neither', 0)} "
            f"pool_neither_sm={n_neither_sm} include_neither={include_neither}",
            flush=True,
        )

    sum_csv = TABLES / "track_b_panel_summary_v1.csv"
    fields = list(summary[0].keys())
    with sum_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(summary)
    print(f"wrote {sum_csv}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
