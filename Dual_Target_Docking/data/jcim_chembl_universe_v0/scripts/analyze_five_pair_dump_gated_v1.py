#!/usr/bin/env python3
"""Dump-gated analyses for the five post-census pairs.

Requires a local ChEMBL 37 sqlite. Cloud VMs without the dump cannot run this.
Does not dock. Does not restock Table 2.

Items:
  - max vs median pChEMBL relabel (θ = 6.0) on frozen Vina scores
  - document year split (AUROC only if dual/A/B each n≥10 after 2018)
  - document-cluster bootstrap and document-blocked GroupKFold
  - leftover unused-pool holdout IDs (20/20/20, seed 20260731, Murcko cap 3)

JAK1/JAK2 leftover B-only is thin (21) but eligible. Do not dock holdout
until these member lists are written.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import random
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
LOCAL = ROOT / "local_track_b_v0"
TAB = ROOT / "tables"
OUT = LOCAL / "tables" / "five_pair_dump_gated_v1"
AN = LOCAL / "analysis"
sys.path.insert(0, str(SCRIPTS))
from pair_ligand_identity_qc_v1 import (  # noqa: E402
    HI,
    LO,
    STANDARD_OK,
    classify,
    connect,
    harvest,
    mol_properties,
    resolve_targets,
)

SEED = 20260729
HOLDOUT_SEED = 20260731
MURCKO_CAP = 3
QUOTA = {"dual": 20, "A_only": 20, "B_only": 20}
N_BOOT = 2000
THETA = 6.0
YEAR_CUTS = (2015, 2018, 2020)
PRIMARY_YEAR_CUT = 2018
EXPECTED_LEFTOVER = {
    # n_strict_smallmol − n_panel from frozen track_b_panel_summary_v1.csv
    "F2/F10": {"dual": 312, "A_only": 76, "B_only": 245},
    "JAK1/TYK2": {"dual": 1874, "A_only": 59, "B_only": 80},
    "JAK1/JAK2": {"dual": 5953, "A_only": 76, "B_only": 21},
    "PPARG/PPARA": {"dual": 408, "A_only": 50, "B_only": 59},
    "PPARA/PPARD": {"dual": 187, "A_only": 50, "B_only": 68},
}

PAIRS = [
    {
        "pair": "F2/F10",
        "genes": ("F2", "F10"),
        "prefix": "HOF2F10",
        "panel": TAB / "track_b_panels" / "panel_F2_F10_v1.csv",
        "target_a": "4UDW",
        "target_b": "2JKH",
    },
    {
        "pair": "JAK1/TYK2",
        "genes": ("JAK1", "TYK2"),
        "prefix": "HOJ1TYK2",
        "panel": TAB / "track_b_panels" / "panel_JAK1_TYK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "3LXP",
    },
    {
        "pair": "JAK1/JAK2",
        "genes": ("JAK1", "JAK2"),
        "prefix": "HOJ1J2",
        "panel": TAB / "track_b_panels" / "panel_JAK1_JAK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "8BXH",
    },
    {
        "pair": "PPARG/PPARA",
        "genes": ("PPARG", "PPARA"),
        "prefix": "HOPGPA",
        "panel": TAB / "track_b_panels" / "panel_PPARG_PPARA_v1.csv",
        "target_a": "9V8H",
        "target_b": "6LXA",
    },
    {
        "pair": "PPARA/PPARD",
        "genes": ("PPARA", "PPARD"),
        "prefix": "HOPAPD",
        "panel": TAB / "track_b_panels" / "panel_PPARA_PPARD_v1.csv",
        "target_a": "6LXA",
        "target_b": "5U3Q",
    },
]


def stable_offset(*parts, modulus=99991):
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus


def r4(x):
    if x is None or (isinstance(x, float) and x != x):
        return ""
    return round(float(x), 4)


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def assign_theta(pA, pB, cut=THETA):
    if pA is None or pB is None:
        return None
    a, b = pA >= cut, pB >= cut
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def load_vina():
    rows = list(csv.DictReader((LOCAL / "tables" / "scores_vina_mode1_v1.csv").open()))
    out = {}
    for r in rows:
        out.setdefault(r["pair"], {}).setdefault(r["ligand"], {})[r["target"]] = float(r["score_S"])
    return out


def harvest_activities(con, tids, chembl_ids):
    """Per-activity pChEMBL + document for the listed molecules/targets."""
    if not chembl_ids:
        return []
    ph_t = ",".join("?" * len(tids))
    ph_s = ",".join("?" * len(STANDARD_OK))
    out = []
    ids = list(chembl_ids)
    for i in range(0, len(ids), 400):
        chunk = ids[i : i + 400]
        ph_m = ",".join("?" * len(chunk))
        sql = f"""
        SELECT md.chembl_id AS molecule_chembl_id, ass.tid, act.pchembl_value,
               docs.year, docs.chembl_id AS document_chembl_id
        FROM activities act
        JOIN assays ass ON ass.assay_id = act.assay_id
        JOIN docs ON docs.doc_id = ass.doc_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        WHERE ass.tid IN ({ph_t})
          AND md.chembl_id IN ({ph_m})
          AND act.pchembl_value IS NOT NULL
          AND act.standard_type IN ({ph_s})
        """
        for r in con.execute(sql, list(tids) + chunk + list(STANDARD_OK)):
            out.append(dict(r))
    return out


def union_find(items):
    parent = {item: item for item in items}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    return find, union


def directional(recs):
    D = [r for r in recs if r["cls"] == "dual"]
    A = [r for r in recs if r["cls"] == "A_only"]
    B = [r for r in recs if r["cls"] == "B_only"]
    da = auroc([r["vina_B"] for r in D], [r["vina_B"] for r in A])
    db = auroc([r["vina_A"] for r in D], [r["vina_A"] for r in B])
    return da, db, min(da, db) if da == da and db == db else float("nan"), len(D), len(A), len(B)


def boot_pm_ci(recs, seed):
    usable = [r for r in recs if r["cls"] in ("dual", "A_only", "B_only")]
    if len(usable) < 8:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    idx = np.arange(len(usable))
    mins = []
    for _ in range(N_BOOT):
        sub = [usable[int(i)] for i in rng.choice(idx, size=len(idx), replace=True)]
        _, _, mn, *_ = directional(sub)
        if mn == mn:
            mins.append(mn)
    if len(mins) < N_BOOT // 2:
        return float("nan"), float("nan")
    lo, hi = np.percentile(mins, [2.5, 97.5])
    return float(lo), float(hi)


def murcko(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except Exception:
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqlite", type=Path, required=True)
    args = ap.parse_args()
    if not args.sqlite.exists():
        print(f"missing sqlite: {args.sqlite}", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    AN.mkdir(parents=True, exist_ok=True)

    con = connect(args.sqlite)
    genes = {g for spec in PAIRS for g in spec["genes"]}
    meta = resolve_targets(con, genes)
    missing = sorted(genes - set(meta))
    if missing:
        print(f"unresolved genes: {missing}", file=sys.stderr)
        return 2
    vina = load_vina()
    maps = harvest(con, {m["tid"] for m in meta.values()})

    maxmed_rows = []
    maxmed_auroc = []
    year_rows = []
    cluster_rows = []
    cv_rows = []
    holdout_all = []
    leftover_counts = []
    leftover_ok_all = True

    for spec in PAIRS:
        pair = spec["pair"]
        ga, gb = spec["genes"]
        tid_a, tid_b = meta[ga]["tid"], meta[gb]["tid"]
        panel = list(csv.DictReader(spec["panel"].open()))
        ids = [r["molecule_chembl_id"] for r in panel]
        acts = harvest_activities(con, {tid_a, tid_b}, ids)
        by = defaultdict(lambda: {"A": [], "B": [], "docs": set(), "years": []})
        for a in acts:
            end = "A" if int(a["tid"]) == tid_a else ("B" if int(a["tid"]) == tid_b else None)
            if end is None:
                continue
            by[a["molecule_chembl_id"]][end].append(float(a["pchembl_value"]))
            if a.get("document_chembl_id"):
                by[a["molecule_chembl_id"]]["docs"].add(a["document_chembl_id"])
            if a.get("year") is not None:
                by[a["molecule_chembl_id"]]["years"].append(int(a["year"]))

        recs_max, recs_med = [], []
        pack = []
        for row in panel:
            cid = row["molecule_chembl_id"]
            lig = row["panel_id"]
            sc = vina.get(pair, {}).get(lig, {})
            sa, sb = sc.get(spec["target_a"]), sc.get(spec["target_b"])
            if sa is None or sb is None:
                continue
            vals = by[cid]
            max_a = max(vals["A"]) if vals["A"] else float(row["pchembl_A"])
            max_b = max(vals["B"]) if vals["B"] else float(row["pchembl_B"])
            med_a = statistics.median(vals["A"]) if vals["A"] else None
            med_b = statistics.median(vals["B"]) if vals["B"] else None
            cls_max = assign_theta(max_a, max_b)
            cls_med = assign_theta(med_a, med_b) if med_a is not None and med_b is not None else None
            recs_max.append({"cls": cls_max, "vina_A": sa, "vina_B": sb})
            if cls_med:
                recs_med.append({"cls": cls_med, "vina_A": sa, "vina_B": sb})
            years = [y for y in vals["years"] if y is not None]
            first_year = min(years) if years else None
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol is not None else None
            pack.append(
                {
                    "ligand": lig,
                    "cls": row["theta6_class"],
                    "vina_A": sa,
                    "vina_B": sb,
                    "documents": sorted(vals["docs"]),
                    "first_year": first_year,
                    "fp": fp,
                    "flip_max_to_median": int(cls_max != cls_med) if cls_med else "",
                    "max_ne_median_A": int(med_a is not None and abs(max_a - med_a) > 1e-6),
                    "max_ne_median_B": int(med_b is not None and abs(max_b - med_b) > 1e-6),
                }
            )
            maxmed_rows.append(
                {
                    "pair": pair,
                    "ligand": lig,
                    "molecule_chembl_id": cid,
                    "n_act_A": len(vals["A"]),
                    "n_act_B": len(vals["B"]),
                    "max_A": r4(max_a),
                    "median_A": r4(med_a),
                    "max_B": r4(max_b),
                    "median_B": r4(med_b),
                    "class_max": cls_max,
                    "class_median": cls_med or "",
                    "flip_max_to_median": int(cls_max != cls_med) if cls_med else "",
                    "n_documents": len(vals["docs"]),
                    "first_year": first_year if first_year is not None else "",
                    "panel_pA": row["pchembl_A"],
                    "panel_pB": row["pchembl_B"],
                    "dump_max_matches_panel_A": int(abs(max_a - float(row["pchembl_A"])) < 0.015),
                    "dump_max_matches_panel_B": int(abs(max_b - float(row["pchembl_B"])) < 0.015),
                }
            )

        for label, recs in (("max_pchembl", recs_max), ("median_pchembl", recs_med)):
            da, db, sm, nD, nA, nB = directional(recs)
            lo, hi = boot_pm_ci(recs, SEED + stable_offset(pair, label))
            maxmed_auroc.append(
                {
                    "pair": pair,
                    "aggregation": label,
                    "n_dual": nD,
                    "n_A_only": nA,
                    "n_B_only": nB,
                    "auroc_D_vs_A": r4(da),
                    "auroc_D_vs_B": r4(db),
                    "summary_min": r4(sm),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "note": "same frozen Vina; labels from dump max vs median at θ=6.0",
                }
            )
        for cut in YEAR_CUTS:
            for split, pred in (
                ("train_first_year_lt", lambda y: y is not None and y < cut),
                ("test_first_year_ge", lambda y: y is not None and y >= cut),
            ):
                sub = [
                    {"cls": r["cls"], "vina_A": r["vina_A"], "vina_B": r["vina_B"]}
                    for r in pack
                    if pred(r.get("first_year"))
                ]
                da, db, sm, nD, nA, nB = directional(sub)
                nN = sum(r["cls"] == "neither" for r in pack if pred(r.get("first_year")))
                n_undated = sum(r.get("first_year") is None for r in pack)
                powered = split.startswith("test") and nD >= 10 and nA >= 10 and nB >= 10
                if powered:
                    lo, hi = boot_pm_ci(sub, SEED + stable_offset(pair, "year", cut, split))
                else:
                    lo = hi = float("nan")
                year_rows.append(
                    {
                        "pair": pair,
                        "cutoff_year": cut,
                        "split": split,
                        "year_definition": "min_document_year_STANDARD_OK_pchembl",
                        "n_dual": nD,
                        "n_A_only": nA,
                        "n_B_only": nB,
                        "n_neither": nN,
                        "n_undated_panel_scored": n_undated,
                        "auroc_reportable": int(powered),
                        "auroc_D_vs_A": r4(da) if powered else "",
                        "auroc_D_vs_B": r4(db) if powered else "",
                        "summary_min": r4(sm) if powered else "",
                        "ci_lo": r4(lo) if powered else "",
                        "ci_hi": r4(hi) if powered else "",
                        "note": (
                            "TIME_SPLIT_PROTOCOL_FREEZE: earliest document.year; "
                            "test AUROC only if dual/A/B each n≥10; do not shop cutoffs; "
                            "train is counts only"
                        ),
                    }
                )

        nodes = []
        for rec in pack:
            nodes.append(("lig", rec["ligand"]))
            nodes.extend(("doc", d) for d in rec["documents"])
        if nodes:
            find, union = union_find(nodes)
            for rec in pack:
                if rec["documents"]:
                    first = rec["documents"][0]
                    union(("lig", rec["ligand"]), ("doc", first))
                    for d in rec["documents"][1:]:
                        union(("doc", first), ("doc", d))
            roots, next_id = {}, 1
            for rec in pack:
                root = find(("lig", rec["ligand"]))
                if root not in roots:
                    roots[root] = f"G{next_id:03d}"
                    next_id += 1
                rec["group_id"] = roots[root]
        else:
            for rec in pack:
                rec["group_id"] = rec["ligand"]

        for contrast, pos_cls, neg_cls, key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            kept = [r for r in pack if r["cls"] in (pos_cls, neg_cls)]
            groups = defaultdict(list)
            for r in kept:
                groups[r["group_id"]].append(r)
            names = list(groups)
            point = auroc(
                [r[key] for r in kept if r["cls"] == pos_cls],
                [r[key] for r in kept if r["cls"] == neg_cls],
            )
            rng = np.random.default_rng(SEED + stable_offset(pair, contrast, "doc"))
            vals = []
            for _ in range(N_BOOT):
                if not names:
                    break
                chosen = rng.choice(names, size=len(names), replace=True)
                sub = [x for g in chosen for x in groups[g]]
                pos = [r[key] for r in sub if r["cls"] == pos_cls]
                neg = [r[key] for r in sub if r["cls"] == neg_cls]
                if len(pos) < 2 or len(neg) < 2:
                    continue
                vals.append(auroc(pos, neg))
            if len(vals) < N_BOOT // 2:
                lo = hi = float("nan")
            else:
                lo, hi = np.percentile(vals, [2.5, 97.5])
            cluster_rows.append(
                {
                    "pair": pair,
                    "contrast": contrast,
                    "n_ligands": len(kept),
                    "n_document_groups": len(names),
                    "n_boot_ok": len(vals),
                    "auroc": r4(point),
                    "ci_lo": r4(lo),
                    "ci_hi": r4(hi),
                    "note": "resample document-connected groups",
                }
            )
            y = np.array([1 if r["cls"] == pos_cls else 0 for r in kept], dtype=int)
            gids = np.array([r["group_id"] for r in kept])
            dock = np.array([[r[key]] for r in kept], dtype=float)
            fps = [r["fp"] for r in kept]
            if any(f is None for f in fps) or len(kept) < 8:
                cv_rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "model": "docking",
                        "cv_auroc": "",
                        "note": "too few ligands or missing fingerprints",
                    }
                )
                continue
            fp = np.vstack([np.asarray(f) for f in fps])
            n_splits = min(5, len(set(gids)), int(y.sum()), int((1 - y).sum()))
            for name, X in (("docking", dock), ("ECFP4", fp), ("ECFP4+docking", np.hstack([fp, dock]))):
                auc = float("nan")
                if n_splits >= 2:
                    try:
                        cv = GroupKFold(n_splits=n_splits)
                        lr = LogisticRegression(max_iter=4000, C=1.0)
                        from sklearn.model_selection import cross_val_predict

                        prob = cross_val_predict(lr, X, y, cv=cv, groups=gids, method="predict_proba")[:, 1]
                        auc = float(roc_auc_score(y, prob))
                    except Exception as exc:
                        auc = float("nan")
                        note_exc = str(exc)[:80]
                    else:
                        note_exc = ""
                cv_rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "model": name,
                        "n": len(kept),
                        "n_groups": len(set(gids)),
                        "n_splits": n_splits if n_splits >= 2 else 0,
                        "cv_auroc": r4(auc),
                        "note": note_exc or "document-blocked GroupKFold",
                    }
                )

        # leftover holdout IDs
        ma, mb = maps[tid_a], maps[tid_b]
        both = sorted(set(ma) & set(mb))
        used = {r["molecule_chembl_id"] for r in panel}
        buckets = defaultdict(list)
        for mol in both:
            pa, pb = ma[mol], mb[mol]
            if pa >= HI and pb >= HI:
                cls = "dual"
            elif pa >= HI and pb <= LO:
                cls = "A_only"
            elif pb >= HI and pa <= LO:
                cls = "B_only"
            else:
                continue
            buckets[cls].append(mol)
        all_mols = {m for v in buckets.values() for m in v}
        props = mol_properties(con, all_mols)
        leftover_pools = {k: [] for k in QUOTA}
        n_strict_sm = {k: 0 for k in QUOTA}
        for cls, mols in buckets.items():
            for mol in mols:
                p = props.get(mol) or {}
                good, _ = classify(p)
                if not good:
                    continue
                n_strict_sm[cls] += 1
                cid = p["chembl_id"]
                if cid in used:
                    continue
                leftover_pools[cls].append(
                    {
                        "molecule_chembl_id": cid,
                        "canonical_smiles": p["canonical_smiles"],
                        "pchembl_A": ma[mol],
                        "pchembl_B": mb[mol],
                        "class": cls,
                    }
                )
        for cls in QUOTA:
            leftover_pools[cls].sort(key=lambda r: r["molecule_chembl_id"])
        exp = EXPECTED_LEFTOVER[pair]
        leftover_ok = all(len(leftover_pools[c]) == exp[c] for c in QUOTA)
        leftover_ok_all = leftover_ok_all and leftover_ok
        leftover_counts.append(
            {
                "pair": pair,
                "n_strict_smallmol_dual": n_strict_sm["dual"],
                "n_strict_smallmol_A_only": n_strict_sm["A_only"],
                "n_strict_smallmol_B_only": n_strict_sm["B_only"],
                "leftover_dual": len(leftover_pools["dual"]),
                "leftover_A_only": len(leftover_pools["A_only"]),
                "leftover_B_only": len(leftover_pools["B_only"]),
                "expected_leftover_dual": exp["dual"],
                "expected_leftover_A_only": exp["A_only"],
                "expected_leftover_B_only": exp["B_only"],
                "matches_frozen_summary": int(leftover_ok),
                "holdout_20_20_20_eligible": int(all(len(leftover_pools[c]) >= 20 for c in QUOTA)),
                "holdout_thin_margin": int(
                    all(len(leftover_pools[c]) >= 20 for c in QUOTA)
                    and min(len(leftover_pools[c]) for c in QUOTA) < 25
                ),
                "holdout_drawn": 0,
                "note": (
                    "JAK1/JAK2 leftover B-only is thin (margin=1); still eligible"
                    if pair == "JAK1/JAK2"
                    else "exclude main-panel IDs; seed 20260731; Murcko cap 3"
                ),
            }
        )
        if not leftover_ok:
            print(
                f"ERROR {pair}: leftover "
                f"{len(leftover_pools['dual'])}/{len(leftover_pools['A_only'])}/{len(leftover_pools['B_only'])} "
                f"!= expected {exp['dual']}/{exp['A_only']}/{exp['B_only']}; "
                f"holdout IDs not drawn",
                flush=True,
            )
            continue
        rng = random.Random(HOLDOUT_SEED)
        picked = []
        scaffold_caps: dict[tuple[str, str], int] = {}
        for cls, need in QUOTA.items():
            pool = list(leftover_pools[cls])
            rng.shuffle(pool)
            got = 0
            for rec in pool:
                if got >= need:
                    break
                scaf = murcko(rec["canonical_smiles"]) or rec["molecule_chembl_id"]
                key = (cls, scaf)
                if scaffold_caps.get(key, 0) >= MURCKO_CAP:
                    continue
                scaffold_caps[key] = scaffold_caps.get(key, 0) + 1
                rec = dict(rec)
                rec["murcko_scaffold"] = scaf
                picked.append(rec)
                got += 1
            if got < need:
                print(f"WARN {pair} {cls}: got {got}/{need}", flush=True)
        picked.sort(key=lambda r: (r["class"], r["molecule_chembl_id"]))
        for i, rec in enumerate(picked, 1):
            holdout_all.append(
                {
                    "holdout_id": f"{spec['prefix']}_{i:03d}",
                    "pair": pair,
                    "class": rec["class"],
                    "molecule_chembl_id": rec["molecule_chembl_id"],
                    "canonical_smiles": rec["canonical_smiles"],
                    "pchembl_A": f"{rec['pchembl_A']:.2f}",
                    "pchembl_B": f"{rec['pchembl_B']:.2f}",
                    "murcko_scaffold": rec["murcko_scaffold"],
                    "label_rule": "strict_6.5_5.5",
                    "holdout_seed": HOLDOUT_SEED,
                    "murcko_cap": MURCKO_CAP,
                    "source": "chembl37_unused_pool_post_panel_freeze",
                }
            )
        leftover_counts[-1]["holdout_drawn"] = len(picked)
        write_csv(
            OUT / f"holdout_panel_{spec['prefix']}_v1.csv",
            [r for r in holdout_all if r["pair"] == pair],
        )
        print(
            f"{pair}: leftover "
            f"{len(leftover_pools['dual'])}/{len(leftover_pools['A_only'])}/{len(leftover_pools['B_only'])} "
            f"match={leftover_ok} holdout={len(picked)}",
            flush=True,
        )

    write_csv(OUT / "max_vs_median_ligand_v1.csv", maxmed_rows)
    write_csv(OUT / "max_vs_median_auroc_v1.csv", maxmed_auroc)
    write_csv(OUT / "time_split_v1.csv", year_rows)
    write_csv(OUT / "document_cluster_bootstrap_v1.csv", cluster_rows)
    write_csv(OUT / "document_blocked_cv_v1.csv", cv_rows)
    write_csv(OUT / "holdout_leftover_counts_v1.csv", leftover_counts)
    write_csv(OUT / "holdout_panels_all_v1.csv", holdout_all)
    n_flip = sum(int(r["flip_max_to_median"] or 0) for r in maxmed_rows)
    n_mismatch = sum(
        1
        for r in maxmed_rows
        if not (int(r["dump_max_matches_panel_A"]) and int(r["dump_max_matches_panel_B"]))
    )
    n_2018_ok = sum(
        1
        for r in year_rows
        if r["cutoff_year"] == PRIMARY_YEAR_CUT and r["split"].startswith("test") and int(r["auroc_reportable"])
    )
    lines = [
        "# Five-pair dump-gated stack\n\n",
        f"sqlite: `{args.sqlite}` (ChEMBL 37; tarball SHA-256 ",
        "`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`).\n\n",
        "Document year = **earliest** `docs.year` among STANDARD_OK pChEMBL rows ",
        "(same endpoints as panel harvest). This is the dump analogue of ",
        "`TIME_SPLIT_PROTOCOL_FREEZE.md`; it is not the later K=4 API ",
        "high-confidence audit.\n\n",
        f"- dump max vs panel pChEMBL mismatches (tol 0.015): **{n_mismatch}** / {len(maxmed_rows)}\n",
        f"- θ=6.0 class flips max→median (scored ligands): **{n_flip}** / {len(maxmed_rows)}\n",
        f"- leftover vs frozen `track_b_panel_summary_v1.csv`: **{'MATCH' if leftover_ok_all else 'MISMATCH'}**\n",
        f"- 2018 test pairs with dual/A/B each n≥10: **{n_2018_ok}** / 5 ",
        "(package as external validation only if ≥2)\n",
        f"- holdout ligands drawn: **{len(holdout_all)}** (seed {HOLDOUT_SEED}, Murcko cap {MURCKO_CAP})\n\n",
        "Holdout IDs are frozen here. Do not re-draw after seeing scores. ",
        "Do not dock until these CSVs exist. JAK1/JAK2 leftover B-only is ",
        "eligible but thin. Does **not** replace Table 2. BindingDB is count-only.\n",
    ]
    (AN / "FIVE_PAIR_DUMP_GATED_V1.md").write_text("".join(lines), encoding="utf-8")
    print("wrote", OUT)
    print("".join(lines))
    return 0 if leftover_ok_all and n_mismatch == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
