#!/usr/bin/env python3
"""Verify that four-state ligands of candidate pairs are real binding small molecules.

For each pair, pulls the strict 6.5/5.5 dual / A-only / B-only sets from a local
ChEMBL SQLite dump and checks ChEMBL-side molecular identity:

  - molecule_type == 'Small molecule' and structure_type == 'MOL'
  - drug-like size window (MW 150-750, heavy atoms 10-60)
  - no peptide / oligosaccharide / oligonucleotide / antibody rows
  - metal-containing formula flag
  - Bemis-Murcko scaffold diversity of the dual set (RDKit), so a "thick" pair
    made of one congeneric series is visible

Does not dock. Does not change Table 2 / K = 4.
"""
from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

TABLES = Path(__file__).resolve().parents[1] / "tables"

THETA = 6.0
HI = 6.5
LO = 5.5
STANDARD_OK = ("IC50", "Ki", "Kd", "EC50", "Potency", "IC50app", "Ki app")

MW_MIN, MW_MAX = 150.0, 750.0
HEAVY_MIN, HEAVY_MAX = 10, 60

METAL_RE = re.compile(
    r"\b(Pt|Ru|Au|Ag|Hg|Cd|As|Sb|Bi|Tc|Re|Gd|Fe|Cu|Zn|Mn|Co|Ni|V|Mo|W|Ti|Sn|Pb|Cr|Pd|Rh|Ir|Os)\b"
)

# Pairs to audit: frozen K=4 first, then the H3-pass fresh-roster candidates.
PAIRS = [
    ("EGFR", "ERBB2", "K4_case"),
    ("ACHE", "BCHE", "K4_development"),
    ("PIK3CA", "PIK3CB", "K4_isoform_control"),
    ("PIK3CA", "MTOR", "K4_development"),
    ("CREBBP", "BRD4", "H3_candidate"),
    ("CNR1", "CNR2", "H3_candidate"),
    ("HCRTR1", "HCRTR2", "H3_candidate"),
    ("SLC6A4", "SLC6A3", "H3_candidate"),
    ("SLC6A2", "SLC6A4", "H3_candidate"),
    ("F2", "F10", "H3_candidate"),
    ("F2", "PRSS1", "H3_candidate"),
    ("JAK1", "TYK2", "H3_candidate"),
    ("JAK1", "JAK2", "H3_candidate"),
    ("JAK3", "TYK2", "H3_candidate"),
    ("PPARG", "PPARA", "H3_candidate"),
    ("PPARA", "PPARD", "H3_candidate"),
    ("CTSK", "CTSS", "H3_candidate"),
    ("OPRM1", "OPRD1", "H3_candidate"),
    ("OPRD1", "OPRK1", "H3_candidate"),
    ("OPRM1", "OPRK1", "H3_candidate"),
    ("S1PR3", "S1PR1", "H3_candidate"),
]


def connect(path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(str(path))
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    con.execute("PRAGMA temp_store = MEMORY")
    con.execute("PRAGMA cache_size = -250000")
    return con


def resolve_targets(con: sqlite3.Connection, genes: set[str]) -> dict[str, dict]:
    sql = """
    SELECT td.tid, td.chembl_id, td.pref_name, cs.accession, csyn.component_synonym AS gene
    FROM target_dictionary td
    JOIN target_components tc ON tc.tid = td.tid
    JOIN component_sequences cs ON cs.component_id = tc.component_id
    JOIN component_synonyms csyn ON csyn.component_id = cs.component_id
    WHERE td.target_type = 'SINGLE PROTEIN'
      AND td.organism = 'Homo sapiens'
      AND csyn.syn_type = 'GENE_SYMBOL'
      AND td.tid IN (SELECT tid FROM target_components GROUP BY tid HAVING COUNT(*) = 1)
    """
    out: dict[str, dict] = {}
    for r in con.execute(sql):
        g = (r["gene"] or "").upper()
        if g in genes and g not in out:
            out[g] = {
                "tid": r["tid"],
                "chembl": r["chembl_id"],
                "uniprot": r["accession"],
                "pref_name": r["pref_name"],
            }
    return out


def harvest(con: sqlite3.Connection, tids: set[int]) -> dict[int, dict[int, float]]:
    ph = ",".join("?" * len(STANDARD_OK))
    tph = ",".join("?" * len(tids))
    sql = f"""
    SELECT ass.tid, act.molregno, MAX(act.pchembl_value) AS pchembl
    FROM activities act
    JOIN assays ass ON ass.assay_id = act.assay_id
    WHERE ass.tid IN ({tph})
      AND act.pchembl_value IS NOT NULL
      AND act.standard_type IN ({ph})
    GROUP BY ass.tid, act.molregno
    """
    maps: dict[int, dict[int, float]] = defaultdict(dict)
    for tid, mol, pv in con.execute(sql, list(tids) + list(STANDARD_OK)):
        maps[int(tid)][int(mol)] = float(pv)
    return dict(maps)


def mol_properties(con: sqlite3.Connection, molregnos: set[int]) -> dict[int, dict]:
    out: dict[int, dict] = {}
    lst = list(molregnos)
    for i in range(0, len(lst), 900):
        chunk = lst[i : i + 900]
        ph = ",".join("?" * len(chunk))
        sql = f"""
        SELECT md.molregno, md.chembl_id, md.molecule_type, md.structure_type,
               md.therapeutic_flag, cp.mw_freebase, cp.heavy_atoms,
               cp.full_molformula, cs.canonical_smiles
        FROM molecule_dictionary md
        LEFT JOIN compound_properties cp ON cp.molregno = md.molregno
        LEFT JOIN compound_structures cs ON cs.molregno = md.molregno
        WHERE md.molregno IN ({ph})
        """
        for r in con.execute(sql, chunk):
            out[int(r["molregno"])] = dict(r)
    return out


def classify(props: dict) -> tuple[bool, str]:
    mt = props.get("molecule_type") or ""
    st = props.get("structure_type") or ""
    mw = props.get("mw_freebase")
    heavy = props.get("heavy_atoms")
    formula = props.get("full_molformula") or ""
    smiles = props.get("canonical_smiles") or ""
    if not smiles:
        return False, "no_structure"
    if st != "MOL":
        return False, f"structure_type_{st or 'NONE'}"
    if mt and mt != "Small molecule":
        return False, f"molecule_type_{mt.replace(' ', '_')}"
    if METAL_RE.search(formula):
        return False, "metal_containing"
    if mw is None or heavy is None:
        return False, "missing_props"
    if not (MW_MIN <= float(mw) <= MW_MAX):
        return False, "mw_out_of_window"
    if not (HEAVY_MIN <= int(heavy) <= HEAVY_MAX):
        return False, "heavy_out_of_window"
    return True, "drug_like_small_molecule"


def scaffold_stats(smiles_list: list[str]) -> dict:
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:
        return {"n_scaffold": "", "singleton_frac": "", "top_scaffold_frac": ""}
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    scaffs = Counter()
    n_ok = 0
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            key = Chem.MolToSmiles(core) if core is not None else ""
        except Exception:
            continue
        if not key:
            key = "(acyclic)"
        scaffs[key] += 1
        n_ok += 1
    if not n_ok:
        return {"n_scaffold": 0, "singleton_frac": "", "top_scaffold_frac": ""}
    singles = sum(1 for v in scaffs.values() if v == 1)
    return {
        "n_scaffold": len(scaffs),
        "singleton_frac": round(singles / len(scaffs), 3),
        "top_scaffold_frac": round(max(scaffs.values()) / n_ok, 3),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqlite", type=Path, required=True)
    args = ap.parse_args()
    if not args.sqlite.exists():
        print(f"missing sqlite: {args.sqlite}", file=sys.stderr)
        return 2
    con = connect(args.sqlite)

    genes = {g for a, b, _ in PAIRS for g in (a, b)}
    meta = resolve_targets(con, genes)
    missing = sorted(genes - set(meta))
    if missing:
        print(f"unresolved genes: {missing}", file=sys.stderr)
    tids = {m["tid"] for m in meta.values()}
    print(f"resolved {len(meta)} targets; harvesting maps...", flush=True)
    maps = harvest(con, tids)

    class_sets: dict[tuple[str, str], dict[str, list[int]]] = {}
    all_mols: set[int] = set()
    for a, b, _ in PAIRS:
        if a not in meta or b not in meta:
            continue
        ma = maps.get(meta[a]["tid"], {})
        mb = maps.get(meta[b]["tid"], {})
        sets = {"strict_dual": [], "strict_A_only": [], "strict_B_only": []}
        for mol in set(ma) & set(mb):
            x, y = ma[mol], mb[mol]
            if x >= HI and y >= HI:
                sets["strict_dual"].append(mol)
            elif x >= HI and y <= LO:
                sets["strict_A_only"].append(mol)
            elif y >= HI and x <= LO:
                sets["strict_B_only"].append(mol)
        class_sets[(a, b)] = sets
        for v in sets.values():
            all_mols.update(v)
    print(f"loading properties for {len(all_mols):,} molecules...", flush=True)
    props = mol_properties(con, all_mols)

    rows = []
    reason_rows = []
    for a, b, role in PAIRS:
        if (a, b) not in class_sets:
            continue
        sets = class_sets[(a, b)]
        rec = {
            "pair": f"{a}/{b}",
            "gene_A": a,
            "gene_B": b,
            "role": role,
            "uniprot_A": meta[a]["uniprot"],
            "uniprot_B": meta[b]["uniprot"],
        }
        for cls, mols in sets.items():
            ok = 0
            reasons = Counter()
            smis = []
            for m in mols:
                p = props.get(m) or {}
                good, why = classify(p)
                reasons[why] += 1
                if good:
                    ok += 1
                    smis.append(p.get("canonical_smiles") or "")
            n = len(mols)
            rec[f"n_{cls}"] = n
            rec[f"n_{cls}_smallmol"] = ok
            rec[f"frac_{cls}_smallmol"] = round(ok / n, 4) if n else ""
            if cls == "strict_dual":
                rec.update({f"dual_{k}": v for k, v in scaffold_stats(smis).items()})
            for why, cnt in reasons.items():
                if why == "drug_like_small_molecule":
                    continue
                reason_rows.append(
                    {"pair": f"{a}/{b}", "class": cls, "reject_reason": why, "n": cnt}
                )
        hn = min(rec.get("n_strict_A_only", 0), rec.get("n_strict_B_only", 0))
        hn_sm = min(rec.get("n_strict_A_only_smallmol", 0), rec.get("n_strict_B_only_smallmol", 0))
        rec["min_hardneg_raw"] = hn
        rec["min_hardneg_smallmol_only"] = hn_sm
        rec["thick_after_smallmol_filter"] = int(hn_sm >= 50)
        rec["hardneg_lost_to_filter"] = hn - hn_sm
        rows.append(rec)

    fields = [
        "pair", "gene_A", "gene_B", "role", "uniprot_A", "uniprot_B",
        "n_strict_dual", "n_strict_dual_smallmol", "frac_strict_dual_smallmol",
        "n_strict_A_only", "n_strict_A_only_smallmol", "frac_strict_A_only_smallmol",
        "n_strict_B_only", "n_strict_B_only_smallmol", "frac_strict_B_only_smallmol",
        "min_hardneg_raw", "min_hardneg_smallmol_only", "hardneg_lost_to_filter",
        "thick_after_smallmol_filter",
        "dual_n_scaffold", "dual_singleton_frac", "dual_top_scaffold_frac",
    ]
    rows.sort(key=lambda r: (-r["min_hardneg_smallmol_only"], r["pair"]))
    with (TABLES / "pair_ligand_identity_qc_v1.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    agg = defaultdict(int)
    for r in reason_rows:
        agg[(r["pair"], r["class"], r["reject_reason"])] += r["n"]
    with (TABLES / "pair_ligand_reject_reasons_v1.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair", "class", "reject_reason", "n"])
        w.writeheader()
        for (pair, cls, why), n in sorted(agg.items(), key=lambda x: (-x[1], x[0])):
            w.writerow({"pair": pair, "class": cls, "reject_reason": why, "n": n})

    print("\npair                 minHN raw->smallmol  thick  dual_scaffolds")
    for r in rows:
        print(
            f"{r['pair']:20s} {r['min_hardneg_raw']:>4} -> {r['min_hardneg_smallmol_only']:<4} "
            f"thick={r['thick_after_smallmol_filter']}  "
            f"dual_n={r['n_strict_dual']:>4} sm={r['n_strict_dual_smallmol']:>4} "
            f"scaff={r['dual_n_scaffold']} top={r['dual_top_scaffold_frac']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
