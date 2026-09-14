#!/usr/bin/env python3
"""Document- and scaffold-cluster bootstrap for fixed-score negative-class Δ.

Primary contrasts:
  EGFR/HER2 pocket A (dual vs B-only vs neither)
  JAK1/TYK2 pocket A (dual vs B-only vs neither)

Δ = AUROC(dual, neither) − AUROC(dual, selective) with shared dual draws.
Resamples connected groups, not ligands. Does not replace the ligand-level
Table S4 interval.
"""
from __future__ import annotations

import csv
import hashlib
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
SQLITE = Path("/tmp/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
N_BOOT = 2000
STANDARD_OK = ("IC50", "Ki", "Kd", "EC50", "Potency", "IC50app", "Ki app")


def stable_offset(*parts, modulus=99991):
    payload = "|".join(map(str, parts)).encode()
    return int(hashlib.sha256(payload).hexdigest()[:16], 16) % modulus


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, float)
    n = np.asarray(neg, float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def scaffold_key(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return "unparsed"
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or "acyclic"
    except Exception:
        return "unparsed"


def union_find(items):
    parent = {x: x for x in items}

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


def assign_doc_groups(recs: list[dict]) -> None:
    nodes = []
    for r in recs:
        nodes.append(("lig", r["ligand"]))
        nodes.extend(("doc", d) for d in r["documents"])
    find, union = union_find(nodes)
    for r in recs:
        if r["documents"]:
            first = r["documents"][0]
            union(("lig", r["ligand"]), ("doc", first))
            for d in r["documents"][1:]:
                union(("doc", first), ("doc", d))
    roots, n = {}, 1
    for r in recs:
        root = find(("lig", r["ligand"]))
        if root not in roots:
            roots[root] = f"G{n:03d}"
            n += 1
        r["doc_group"] = roots[root]


def harvest_docs(con, tids, chembl_ids) -> dict[str, set[str]]:
    out = defaultdict(set)
    ids = list(chembl_ids)
    ph_t = ",".join("?" * len(tids))
    ph_s = ",".join("?" * len(STANDARD_OK))
    for i in range(0, len(ids), 400):
        chunk = ids[i : i + 400]
        ph_m = ",".join("?" * len(chunk))
        sql = f"""
        SELECT md.chembl_id AS molecule_chembl_id, docs.chembl_id AS document_chembl_id
        FROM activities act
        JOIN assays ass ON ass.assay_id = act.assay_id
        JOIN docs ON docs.doc_id = ass.doc_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        WHERE ass.tid IN ({ph_t})
          AND md.chembl_id IN ({ph_m})
          AND act.pchembl_value IS NOT NULL
          AND act.standard_type IN ({ph_s})
        """
        for mol, doc in con.execute(sql, list(tids) + chunk + list(STANDARD_OK)):
            if doc:
                out[mol].add(doc)
    return out


def cluster_delta(recs, estimator: str, pair: str, contrast: str) -> dict:
    key = "doc_group" if estimator == "document_cluster" else "scaffold"
    groups = defaultdict(list)
    for r in recs:
        groups[r[key]].append(r)
    gids = list(groups)
    point_dual = [r["score"] for r in recs if r["cls"] == "dual"]
    point_sel = [r["score"] for r in recs if r["cls"] == "selective"]
    point_nei = [r["score"] for r in recs if r["cls"] == "neither"]
    point = auroc(point_dual, point_nei) - auroc(point_dual, point_sel)
    seed_offset = stable_offset(pair, contrast, estimator)
    rng = np.random.default_rng(seed_offset)
    deltas = []
    for _ in range(N_BOOT):
        draw = rng.choice(gids, size=len(gids), replace=True)
        bag = [r for gid in draw for r in groups[gid]]
        d = [r["score"] for r in bag if r["cls"] == "dual"]
        s = [r["score"] for r in bag if r["cls"] == "selective"]
        n = [r["score"] for r in bag if r["cls"] == "neither"]
        if min(len(d), len(s), len(n)) < 2:
            continue
        deltas.append(auroc(d, n) - auroc(d, s))
    if len(deltas) < N_BOOT // 2:
        lo = hi = float("nan")
        note = "too few valid cluster resamples"
    else:
        lo, hi = (float(x) for x in np.percentile(deltas, [2.5, 97.5]))
        note = "resamples connected groups, not ligands; shared dual within each resample"
    return {
        "pair": pair,
        "contrast": contrast,
        "estimator": estimator,
        "n_groups": len(gids),
        "n_dual": len(point_dual),
        "n_selective": len(point_sel),
        "n_neither": len(point_nei),
        "delta_point": round(float(point), 4),
        "delta_ci_lo": None if lo != lo else round(lo, 4),
        "delta_ci_hi": None if hi != hi else round(hi, 4),
        "n_valid_boot": len(deltas),
        "excludes_zero": bool(lo == lo and (lo > 0 or hi < 0)),
        "note": note,
        "seed_offset": seed_offset,
    }


def load_egfr() -> list[dict]:
    scores = {r["ligand"]: r for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")}
    smiles = {r["panel_id"]: r["smiles"] for r in read_csv(ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv")}
    groups = {
        r["ligand"]: r["group_id"]
        for r in read_csv(TAB / "document_blocked_ligand_groups_v1.csv")
        if r["pair"] == "EGFR/HER2"
    }
    recs = []
    for lig, row in scores.items():
        cls = row["class"]
        if cls == "B_only":
            use = "selective"
        elif cls == "neither":
            use = "neither"
        elif cls == "dual":
            use = "dual"
        else:
            continue
        aff = row.get("3POZ_affinity")
        if aff in ("", None):
            continue
        recs.append(
            {
                "ligand": lig,
                "cls": use,
                "score": -float(aff),
                "doc_group": groups.get(lig, lig),
                "scaffold": scaffold_key(smiles.get(lig, "")),
                "documents": [],
            }
        )
    return recs


def load_jak1() -> list[dict]:
    panel = read_csv(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_TYK2_v1.csv")
    vina_rows = read_csv(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv")
    score = {}
    for r in vina_rows:
        if r["pair"] != "JAK1/TYK2" or r["status"] != "success":
            continue
        score.setdefault(r["ligand"], {})[r["target"]] = float(r["score_S"])
    if not SQLITE.exists():
        raise SystemExit(f"missing {SQLITE}")
    con = sqlite3.connect(SQLITE)
    con.row_factory = sqlite3.Row
    docs = harvest_docs(con, jak_tids(con), [r["molecule_chembl_id"] for r in panel])
    con.close()
    recs = []
    for r in panel:
        lig = r["panel_id"]
        cls = r["theta6_class"]
        if cls == "B_only":
            use = "selective"
        elif cls == "neither":
            use = "neither"
        elif cls == "dual":
            use = "dual"
        else:
            continue
        ends = score.get(lig, {})
        if "6N7A" not in ends or "3LXP" not in ends:
            continue
        recs.append(
            {
                "ligand": lig,
                "cls": use,
                "score": float(ends["6N7A"]),
                "documents": sorted(docs.get(r["molecule_chembl_id"], set())),
                "scaffold": scaffold_key(r["canonical_smiles"]),
            }
        )
    assign_doc_groups(recs)
    return recs


def jak_tids(con) -> set[int]:
    rows = list(
        con.execute(
            """
            SELECT td.tid, c.accession
            FROM target_dictionary td
            JOIN target_type tt ON tt.target_type = td.target_type
            JOIN component_sequences c ON c.component_id = (
                SELECT tc.component_id FROM target_components tc
                WHERE tc.tid = td.tid LIMIT 1
            )
            WHERE td.organism = 'Homo sapiens'
              AND td.target_type = 'SINGLE PROTEIN'
              AND c.accession IN ('P23458', 'P29597')
            """
        )
    )
    tids = {int(r[0]) for r in rows}
    if len(tids) < 2:
        raise SystemExit(f"JAK1/TYK2 tids not resolved: {rows}")
    return tids


def main() -> int:
    TAB.mkdir(parents=True, exist_ok=True)
    rows = []
    packs = {
        ("EGFR/HER2", "D_vs_B_or_neither_pocketA"): load_egfr(),
        ("JAK1/TYK2", "D_vs_B_or_neither_pocketA"): load_jak1(),
    }
    for (pair, contrast), recs in packs.items():
        n = Counter(r["cls"] for r in recs)
        print(pair, dict(n), "n", len(recs), file=sys.stderr)
        for estimator in ("document_cluster", "scaffold_cluster"):
            row = cluster_delta(recs, estimator, pair, contrast)
            rows.append(row)
            print(row, file=sys.stderr)
    out = TAB / "equal_score_cluster_bootstrap_v1.csv"
    write_csv(out, rows)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
