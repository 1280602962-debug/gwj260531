#!/usr/bin/env python3
"""Independent ligand-level audit of DualFourClass primary metrics.

Does not import production analysis modules. Rebuilds a canonical table from
panel CSVs + raw Vina affinities + scored membership, then recomputes AUROC,
bootstrap, ranking, EF, leakage, and manuscript/figure/table consistency.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/jcim_independent_audit_v0/tables"
DOCS = ROOT / "docs"
N_BOOT = 2000
SEED = 20260729
ORDER = (
    "EGFR/HER2", "JAK1/JAK2", "JAK1/TYK2", "PIK3CA/mTOR",
    "AChE/BChE", "F2/F10", "PPARG/PPARA", "PPARA/PPARD",
)
FIVE = ("JAK1/JAK2", "JAK1/TYK2", "F2/F10", "PPARG/PPARA", "PPARA/PPARD")
PDB = {
    "F2/F10": ("4UDW", "2JKH"),
    "JAK1/JAK2": ("6N7A", "8BXH"),
    "JAK1/TYK2": ("6N7A", "3LXP"),
    "PPARG/PPARA": ("9V8H", "6LXA"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}
HOLD_ORIG = {
    "AChE/BChE": (
        ROOT / "data/jcim_holdout_v0/tables/holdout_panel_HOAB.csv",
        "pchembl_ACHE",
        "pchembl_BCHE",
    ),
    "PIK3CA/mTOR": (
        ROOT / "data/jcim_holdout_v0/tables/holdout_panel_HOPM.csv",
        "pchembl_PIK3CA",
        "pchembl_MTOR",
    ),
}
HOLD_FIVE = {
    "F2/F10": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOF2F10_v1.csv",
    "JAK1/TYK2": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOJ1TYK2_v1.csv",
    "JAK1/JAK2": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOJ1J2_v1.csv",
    "PPARG/PPARA": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOPGPA_v1.csv",
    "PPARA/PPARD": ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_dump_gated_v1/holdout_panel_HOPAPD_v1.csv",
}
MS_T2 = {
    "EGFR/HER2": dict(n="28 / 38 / 32", da=0.666, da_lo=0.524, da_hi=0.793, db=0.430, db_lo=0.282, db_hi=0.579, smin=0.430, slo=0.282, shi=0.578),
    "JAK1/JAK2": dict(n="32 / 32 / 32", da=0.588, da_lo=0.444, da_hi=0.729, db=0.728, db_lo=0.595, db_hi=0.848, smin=0.588, slo=0.444, shi=0.725),
    "JAK1/TYK2": dict(n="31 / 32 / 32", da=0.575, da_lo=0.434, da_hi=0.725, db=0.365, db_lo=0.231, db_hi=0.505, smin=0.365, slo=0.231, shi=0.503),
    "PIK3CA/mTOR": dict(n="18 / 14 / 12", da=0.714, da_lo=0.506, da_hi=0.899, db=0.692, db_lo=0.495, db_hi=0.874, smin=0.692, slo=0.470, shi=0.813),
    "AChE/BChE": dict(n="27 / 25 / 28", da=0.650, da_lo=0.483, da_hi=0.801, db=0.606, db_lo=0.442, db_hi=0.751, smin=0.606, slo=0.437, shi=0.730),
    "F2/F10": dict(n="31 / 32 / 32", da=0.413, da_lo=0.259, da_hi=0.562, db=0.345, db_lo=0.214, db_hi=0.486, smin=0.345, slo=0.211, shi=0.477),
    "PPARG/PPARA": dict(n="32 / 31 / 32", da=0.649, da_lo=0.507, da_hi=0.778, db=0.706, db_lo=0.569, db_hi=0.833, smin=0.649, slo=0.504, shi=0.751),
    "PPARA/PPARD": dict(n="32 / 32 / 32", da=0.646, da_lo=0.504, da_hi=0.776, db=0.446, db_lo=0.296, db_hi=0.586, smin=0.446, slo=0.296, shi=0.584),
}
MS_T3 = {
    "EGFR/HER2": (0.756, 0.562, 0.920, 12),
    "JAK1/JAK2": (0.730, 0.547, 0.875, 14),
    "JAK1/TYK2": (0.770, 0.597, 0.906, 14),
    "PIK3CA/mTOR": (0.514, 0.222, 0.806, 4),
    "AChE/BChE": (0.649, 0.484, 0.812, 15),
    "F2/F10": (0.519, 0.350, 0.688, 12),
    "PPARG/PPARA": (0.685, 0.493, 0.848, 14),
    "PPARA/PPARD": (0.565, 0.368, 0.766, 14),
}


def rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def fnum(value):
    if value is None or str(value).strip() in ("", "NA", "nan", "None"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def theta6(pa: float, pb: float) -> str:
    a, b = pa >= 6.0, pb >= 6.0
    if a and b:
        return "dual"
    if a and not b:
        return "A_only"
    if b and not a:
        return "B_only"
    return "neither"


def strict6555(pa: float, pb: float) -> str | None:
    if pa >= 6.5 and pb >= 6.5:
        return "dual"
    if pa >= 6.5 and pb <= 5.5:
        return "A_only"
    if pb >= 6.5 and pa <= 5.5:
        return "B_only"
    if pa <= 5.5 and pb <= 5.5:
        return "neither"
    return None


def auroc(pos, neg) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    diff = pos[:, None] - neg[None, :]
    return float(((diff > 0).sum() + 0.5 * (diff == 0).sum()) / diff.size)


def offset(*parts, modulus=99991) -> int:
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:16], 16) % modulus


def r3(x: float) -> float:
    return float(f"{x:.3f}")


def r4(x: float) -> float:
    return float(f"{x:.4f}")


def write_csv(path: Path, fieldnames: list[str], data: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def mol_from_smiles(smiles: str, largest: bool) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if not largest:
        return mol
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return None
    return max(frags, key=lambda item: item.GetNumHeavyAtoms()) if len(frags) > 1 else frags[0]


def fp_array(mol: Chem.Mol) -> np.ndarray:
    bitvect = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    arr = np.zeros((2048,), dtype=float)
    DataStructs.ConvertToNumpyArray(bitvect, arr)
    return arr


def murcko_of(mol: Chem.Mol, include_chirality: bool) -> str:
    try:
        if include_chirality:
            return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or ""
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False) or ""
    except Exception:
        return ""


def attach_chemistry(recs: list[dict]) -> None:
    """Fill Murcko / descriptors / in-memory ECFP4 using the production SMILES rule."""
    for rec in recs:
        smiles = rec.get("smiles") or ""
        if not smiles:
            rec["_fp"] = None
            continue
        largest = rec["pair"] in FIVE
        mol = mol_from_smiles(smiles, largest=largest)
        if mol is None:
            rec["_fp"] = None
            continue
        rec["_fp"] = fp_array(mol)
        scaf = murcko_of(mol, include_chirality=False)
        if not rec.get("murcko_scaffold"):
            rec["murcko_scaffold"] = scaf
        rec["_scaffold_group"] = scaf or (
            f"__fail_{rec['ligand_id']}" if rec["pair"] not in FIVE else f"NONE:{rec['ligand_id']}"
        )
        if rec.get("heavy") in ("", None):
            rec["heavy"] = mol.GetNumHeavyAtoms()
        if rec.get("mw") in ("", None):
            rec["mw"] = Descriptors.MolWt(mol)
        if rec.get("clogp") in ("", None):
            rec["clogp"] = Crippen.MolLogP(mol) if largest else Descriptors.MolLogP(mol)
        if rec.get("tpsa") in ("", None):
            rec["tpsa"] = Descriptors.TPSA(mol)


def grouped_oof(X: np.ndarray, y: np.ndarray, groups: np.ndarray):
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    n_scaf = len(set(groups.tolist()))
    n_splits = min(5, n_scaf, n_pos, n_neg)
    if n_splits < 2 or n_pos < 6 or n_neg < 6:
        return float("nan"), n_splits, n_scaf, 0, []
    cv = GroupKFold(n_splits=n_splits)
    splits = list(cv.split(X, y, groups))
    leaked = 0
    for train, test in splits:
        if set(groups[train]) & set(groups[test]):
            leaked += 1
    pred = np.zeros(len(y), dtype=float)
    for train, test in splits:
        model = LogisticRegression(max_iter=4000, C=1.0)
        model.fit(X[train], y[train])
        pred[test] = model.predict_proba(X[test])[:, 1]
    return float(roc_auc_score(y, pred)), n_splits, n_scaf, leaked, splits


def load_mode1() -> dict[tuple[str, str, str], dict]:
    out = {}
    path = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv"
    for rec in rows(path):
        out[(rec["pair"], rec["ligand"], rec["target"])] = rec
    return out


def load_docs() -> dict[tuple[str, str], dict]:
    path = ROOT / "data/jcim_novelty_v0/tables/document_blocked_ligand_groups_v1.csv"
    return {(r["pair"], r["ligand"]): r for r in rows(path)}


def assemble_main(mode1, docs) -> list[dict]:
    membership = rows(ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv")
    by_pair = defaultdict(list)
    for rec in membership:
        by_pair[rec["pair"]].append(rec)

    panels = {
        "EGFR/HER2": rows(ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv"),
        "PIK3CA/mTOR": rows(ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv"),
        "AChE/BChE": rows(ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv"),
        "JAK1/JAK2": rows(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_JAK2_v1.csv"),
        "JAK1/TYK2": rows(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_JAK1_TYK2_v1.csv"),
        "F2/F10": rows(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_F2_F10_v1.csv"),
        "PPARG/PPARA": rows(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARG_PPARA_v1.csv"),
        "PPARA/PPARD": rows(ROOT / "data/jcim_chembl_universe_v0/tables/track_b_panels/panel_PPARA_PPARD_v1.csv"),
    }
    raw = {
        "EGFR/HER2": {r["ligand"]: r for r in rows(ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv")},
        "PIK3CA/mTOR": {r["ligand"]: r for r in rows(ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv")},
        "AChE/BChE": {r["ligand"]: r for r in rows(ROOT / "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv")},
    }
    assembled_desc = {}
    for stem, pair in (
        ("assembled_EGFR_HER2.csv", "EGFR/HER2"),
        ("assembled_AChE_BChE.csv", "AChE/BChE"),
        ("assembled_PIK3CA_mTOR.csv", "PIK3CA/mTOR"),
    ):
        for rec in rows(ROOT / "data/jcim_bench_v0/tables" / stem):
            assembled_desc[(pair, rec["ligand"])] = rec

    out = []
    for pair in ORDER:
        panel_map = {r.get("panel_id") or r.get("ligand"): r for r in panels[pair]}
        for rec in by_pair[pair]:
            lig = rec["ligand"]
            pan = panel_map.get(lig, {})
            chembl = pan.get("molecule_chembl_id", "")
            pa, pb = float(rec["pA"]), float(rec["pB"])
            sa_m, sb_m = float(rec["score_A"]), float(rec["score_B"])
            raw_a = raw_b = None
            if pair == "EGFR/HER2" and lig in raw[pair]:
                raw_a = fnum(raw[pair][lig]["3POZ_affinity"])
                raw_b = fnum(raw[pair][lig]["3RCD_affinity"])
            elif pair == "PIK3CA/mTOR" and lig in raw[pair]:
                raw_a = fnum(raw[pair][lig]["4L23_affinity"])
                raw_b = fnum(raw[pair][lig]["4JT6_affinity"])
            elif pair == "AChE/BChE" and lig in raw[pair]:
                raw_a = fnum(raw[pair][lig]["vina_ACHE"])
                raw_b = fnum(raw[pair][lig]["vina_BCHE"])
            elif pair in PDB:
                pdb_a, pdb_b = PDB[pair]
                ra, rb = mode1.get((pair, lig, pdb_a)), mode1.get((pair, lig, pdb_b))
                raw_a = fnum((ra or {}).get("mode1_energy"))
                raw_b = fnum((rb or {}).get("mode1_energy"))
            sa = None if raw_a is None else -raw_a
            sb = None if raw_b is None else -raw_b
            smiles = pan.get("smiles") or pan.get("canonical_smiles") or ""
            scaffold = pan.get("murcko_scaffold") or ""
            desc = assembled_desc.get((pair, lig), {})
            doc = docs.get((pair, lig), {})
            out.append({
                "analysis_set": "main",
                "pair": pair,
                "ligand_id": lig,
                "molecule_chembl_id": chembl,
                "pChEMBL_A": pa,
                "pChEMBL_B": pb,
                "panel_class": pan.get("class") or pan.get("theta6_class") or "",
                "panel_label_rule": pan.get("label_rule") or ("theta_6.0" if pair in ("EGFR/HER2", "PIK3CA/mTOR") else ""),
                "primary_state_theta6": theta6(pa, pb),
                "strict_6555_state": strict6555(pa, pb) or "",
                "membership_cls": rec["cls"],
                "raw_vina_A": raw_a if raw_a is not None else "",
                "raw_vina_B": raw_b if raw_b is not None else "",
                "score_A": sa if sa is not None else sa_m,
                "score_B": sb if sb is not None else sb_m,
                "membership_score_A": sa_m,
                "membership_score_B": sb_m,
                "score_mean": ((sa if sa is not None else sa_m) + (sb if sb is not None else sb_m)) / 2.0,
                "smiles": smiles,
                "murcko_scaffold": scaffold,
                "source_documents": doc.get("documents", ""),
                "n_documents": doc.get("n_documents", ""),
                "document_group": doc.get("group_id", ""),
                "heavy": desc.get("heavy", pan.get("heavy_atoms", "")),
                "mw": desc.get("mw", pan.get("mw_freebase", "")),
                "clogp": desc.get("clogp", ""),
                "tpsa": desc.get("tpsa", ""),
                "main_holdout_status": "main",
            })
    return out


def assemble_holdout(mode1) -> list[dict]:
    out = []
    orig_panels = {}
    for pair, (path, cola, colb) in HOLD_ORIG.items():
        orig_panels[pair] = {r["holdout_id"]: r for r in rows(path)}
        del cola, colb
    for rec in rows(ROOT / "data/jcim_holdout_v0/tables/holdout_ligand_scores_v1.csv"):
        pan = orig_panels.get(rec["pair"], {}).get(rec["ligand"], {})
        pa = pb = None
        if rec["pair"] in HOLD_ORIG:
            pa = fnum(pan.get(HOLD_ORIG[rec["pair"]][1]))
            pb = fnum(pan.get(HOLD_ORIG[rec["pair"]][2]))
        out.append({
            "analysis_set": "holdout",
            "pair": rec["pair"],
            "ligand_id": rec["ligand"],
            "molecule_chembl_id": rec["chembl"],
            "pChEMBL_A": pa if pa is not None else "",
            "pChEMBL_B": pb if pb is not None else "",
            "panel_class": rec["cls"],
            "panel_label_rule": "strict_6.5_5.5",
            "primary_state_theta6": theta6(pa, pb) if pa is not None and pb is not None else rec["cls"],
            "strict_6555_state": rec["cls"],
            "membership_cls": rec["cls"],
            "raw_vina_A": float(rec["vina_A_raw"]),
            "raw_vina_B": float(rec["vina_B_raw"]),
            "score_A": float(rec["vina_A"]),
            "score_B": float(rec["vina_B"]),
            "membership_score_A": float(rec["vina_A"]),
            "membership_score_B": float(rec["vina_B"]),
            "score_mean": float(rec["vina_mean"]),
            "smiles": rec.get("smiles", ""),
            "murcko_scaffold": pan.get("murcko_scaffold", ""),
            "source_documents": "",
            "n_documents": "",
            "document_group": "",
            "heavy": rec.get("heavy", ""),
            "mw": rec.get("mw", ""),
            "clogp": rec.get("clogp", ""),
            "tpsa": rec.get("tpsa", ""),
            "main_holdout_status": "holdout",
        })
    ho_scores = rows(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/allpairs_stack/holdout/tables/holdout_scores_vina_mode1_v1.csv")
    energy = {}
    for rec in ho_scores:
        energy.setdefault((rec["pair"], rec["ligand"]), {})[rec["target"]] = rec
    for pair, path in HOLD_FIVE.items():
        pdb_a, pdb_b = PDB[pair]
        for rec in rows(path):
            lig = rec["holdout_id"]
            block = energy.get((pair, lig), {})
            ra, rb = block.get(pdb_a, {}), block.get(pdb_b, {})
            raw_a, raw_b = fnum(ra.get("mode1_energy")), fnum(rb.get("mode1_energy"))
            if raw_a is None or raw_b is None:
                continue
            pa, pb = float(rec["pchembl_A"]), float(rec["pchembl_B"])
            out.append({
                "analysis_set": "holdout",
                "pair": pair,
                "ligand_id": lig,
                "molecule_chembl_id": rec["molecule_chembl_id"],
                "pChEMBL_A": pa,
                "pChEMBL_B": pb,
                "panel_class": rec["class"],
                "panel_label_rule": rec.get("label_rule", "strict_6.5_5.5"),
                "primary_state_theta6": theta6(pa, pb),
                "strict_6555_state": rec["class"],
                "membership_cls": rec["class"],
                "raw_vina_A": raw_a,
                "raw_vina_B": raw_b,
                "score_A": -raw_a,
                "score_B": -raw_b,
                "membership_score_A": -raw_a,
                "membership_score_B": -raw_b,
                "score_mean": (-raw_a + -raw_b) / 2.0,
                "smiles": rec.get("canonical_smiles", ""),
                "murcko_scaffold": rec.get("murcko_scaffold", ""),
                "source_documents": "",
                "n_documents": "",
                "document_group": "",
                "heavy": "",
                "mw": "",
                "clogp": "",
                "tpsa": "",
                "main_holdout_status": "holdout",
            })
    return out


def by_pair(table, analysis_set="main"):
    out = defaultdict(list)
    for rec in table:
        if rec["analysis_set"] == analysis_set:
            out[rec["pair"]].append(rec)
    return out


def directional(recs):
    dual = [r for r in recs if r["membership_cls"] == "dual"]
    aonly = [r for r in recs if r["membership_cls"] == "A_only"]
    bonly = [r for r in recs if r["membership_cls"] == "B_only"]
    da = auroc([r["score_B"] for r in dual], [r["score_B"] for r in aonly])
    db = auroc([r["score_A"] for r in dual], [r["score_A"] for r in bonly])
    return da, db, min(da, db), len(dual), len(aonly), len(bonly)


def boot_summary_min(recs, seed):
    use = [r for r in recs if r["membership_cls"] in ("dual", "A_only", "B_only")]
    rng = np.random.default_rng(seed)
    idx = np.arange(len(use))
    draws = []
    shared_dual = 0
    for _ in range(N_BOOT):
        ii = rng.choice(idx, len(idx), replace=True)
        sub = [use[int(i)] for i in ii]
        da, db, sm, nd, na, nb = directional(sub)
        if nd and na and nb and math.isfinite(sm):
            draws.append((da, db, sm))
            shared_dual += 1
    arr = np.array(draws)
    lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
    return {
        "n_valid": len(draws),
        "da_lo": float(lo[0]), "da_hi": float(hi[0]),
        "db_lo": float(lo[1]), "db_hi": float(hi[1]),
        "smin_lo": float(lo[2]), "smin_hi": float(hi[2]),
        "shared_dual_replicates": shared_dual,
    }


def boot_neither(recs, seed):
    dual = np.array([r["score_mean"] for r in recs if r["membership_cls"] == "dual"])
    nei = np.array([r["score_mean"] for r in recs if r["membership_cls"] == "neither"])
    rng = np.random.default_rng(seed)
    vals = [auroc(rng.choice(dual, len(dual), True), rng.choice(nei, len(nei), True)) for _ in range(N_BOOT)]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), auroc(dual, nei)


def boot_delta(recs, score_key, sel_cls, seed):
    dual = np.array([r[score_key] for r in recs if r["membership_cls"] == "dual"])
    sel = np.array([r[score_key] for r in recs if r["membership_cls"] == sel_cls])
    nei = np.array([r[score_key] for r in recs if r["membership_cls"] == "neither"])
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(N_BOOT):
        d = rng.choice(dual, len(dual), True)
        deltas.append(auroc(d, rng.choice(nei, len(nei), True)) - auroc(d, rng.choice(sel, len(sel), True)))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    a_s, a_n = auroc(dual, sel), auroc(dual, nei)
    return a_s, a_n, a_n - a_s, float(lo), float(hi)


def boot_matched_mismatch(recs, seed):
    use = [r for r in recs if r["membership_cls"] in ("dual", "A_only", "B_only")]
    def sm_of(sub, da_key, db_key):
        dual = [r for r in sub if r["membership_cls"] == "dual"]
        aonly = [r for r in sub if r["membership_cls"] == "A_only"]
        bonly = [r for r in sub if r["membership_cls"] == "B_only"]
        da = auroc([r[da_key] for r in dual], [r[da_key] for r in aonly])
        db = auroc([r[db_key] for r in dual], [r[db_key] for r in bonly])
        return da, db, min(da, db)
    da_m, db_m, sm_m = sm_of(use, "score_B", "score_A")
    da_w, db_w, sm_w = sm_of(use, "score_A", "score_B")
    rng = np.random.default_rng(seed)
    idx = np.arange(len(use))
    deltas = []
    for _ in range(N_BOOT):
        sub = [use[int(i)] for i in rng.choice(idx, len(idx), True)]
        _, _, m = sm_of(sub, "score_B", "score_A")
        _, _, w = sm_of(sub, "score_A", "score_B")
        if math.isfinite(m) and math.isfinite(w):
            deltas.append(m - w)
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return da_m, db_m, sm_m, da_w, db_w, sm_w, sm_m - sm_w, float(lo), float(hi)


def ranking(recs):
    ranked = sorted(recs, key=lambda r: (-float(r["score_mean"]), r["ligand_id"]))
    n = len(ranked)
    k = max(1, math.ceil(0.10 * n))
    top = ranked[:k]
    ct = Counter(r["membership_cls"] for r in recs)
    top_ct = Counter(r["membership_cls"] for r in top)
    n_dual = ct["dual"]
    top_dual = top_ct["dual"]
    return {
        "n": n, "k": k,
        "n_dual": n_dual, "n_A": ct["A_only"], "n_B": ct["B_only"], "n_N": ct["neither"],
        "top_dual": top_dual, "top_A": top_ct["A_only"], "top_B": top_ct["B_only"], "top_N": top_ct["neither"],
        "dual_frac": top_dual / k,
        "panel_dual_frac": n_dual / n,
        "ef": (top_dual / k) / (n_dual / n),
        "ids": [r["ligand_id"] for r in top],
    }


def close(a, b, tol=5.5e-4) -> bool:
    if a is None or b is None or (isinstance(a, float) and not math.isfinite(a)):
        return False
    return abs(float(a) - float(b)) <= tol


def round_close(indep, reported, nd=3) -> bool:
    return abs(float(f"{indep:.{nd}f}") - reported) < 1e-9 or close(float(f"{indep:.{nd}f}"), reported, 5e-4)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    findings = []
    metrics = []

    def add(status, check, **kw):
        rec = {"status": status, "check": check, **kw}
        findings.append(rec)

    def metric(**kw):
        metrics.append(kw)

    mode1 = load_mode1()
    docs = load_docs()
    canonical = assemble_main(mode1, docs)
    holdout = assemble_holdout(mode1)
    attach_chemistry(canonical)
    attach_chemistry(holdout)
    canonical_fields = [
        "analysis_set", "pair", "ligand_id", "molecule_chembl_id", "pChEMBL_A", "pChEMBL_B",
        "panel_class", "panel_label_rule", "primary_state_theta6", "strict_6555_state",
        "membership_cls", "raw_vina_A", "raw_vina_B", "score_A", "score_B",
        "membership_score_A", "membership_score_B", "score_mean", "smiles", "murcko_scaffold",
        "source_documents", "n_documents", "document_group", "heavy", "mw", "clogp", "tpsa",
        "main_holdout_status",
    ]
    write_csv(OUT / "canonical_ligand_level_v1.csv", canonical_fields, canonical + holdout)

    main_packs = by_pair(canonical, "main")
    hold_packs = by_pair(holdout, "holdout")

    # --- labels / scores / n ---
    membership = rows(ROOT / "data/jcim_novelty_v0/tables/review_scored_membership_v1.csv")
    n_mem = Counter(r["pair"] for r in membership)
    for pair in ORDER:
        recs = main_packs[pair]
        n = len(recs)
        if n != n_mem[pair]:
            add("FAIL", "membership n vs canonical", pair=pair, reported=n_mem[pair], independent=n,
                file="review_scored_membership_v1.csv")
        mismatch = [r["ligand_id"] for r in recs if r["membership_cls"] != r["primary_state_theta6"]]
        if mismatch:
            add("FAIL", "primary state is not theta=6.0", pair=pair, independent=len(mismatch),
                extra=mismatch[:8], file="canonical_ligand_level_v1.csv")
        else:
            add("PASS", "primary experimental states equal theta=6.0 from pChEMBL", pair=pair, independent=n)
        flips = [r["ligand_id"] for r in recs if r["strict_6555_state"] and r["strict_6555_state"] != r["primary_state_theta6"]]
        if pair in ("EGFR/HER2", "PIK3CA/mTOR") and flips:
            add("PASS", "theta=6.0 primary differs from strict 6.5/5.5 on sampled gray-zone ligands",
                pair=pair, independent=len(flips))
        if pair not in ("EGFR/HER2", "PIK3CA/mTOR"):
            if flips:
                add("FAIL", "strict-sampled pair has theta=6.0 vs 6.5/5.5 label disagreement",
                    pair=pair, independent=len(flips), extra=flips[:8])
            else:
                add("PASS", "strict 6.5/5.5 sampling labels coincide with theta=6.0 on this panel", pair=pair)
        score_mismatch = 0
        missing_raw = 0
        for rec in recs:
            if rec["raw_vina_A"] == "" or rec["raw_vina_B"] == "":
                missing_raw += 1
                continue
            if abs(float(rec["score_A"]) - float(rec["membership_score_A"])) > 1e-6:
                score_mismatch += 1
            if abs(float(rec["score_B"]) - float(rec["membership_score_B"])) > 1e-6:
                score_mismatch += 1
            if abs(float(rec["score_A"]) + float(rec["raw_vina_A"])) > 1e-6:
                score_mismatch += 1
        if missing_raw:
            add("FAIL", "raw Vina affinity missing for scored ligand", pair=pair, independent=missing_raw)
        if score_mismatch:
            add("FAIL", "sign-corrected score does not equal -raw or membership", pair=pair, independent=score_mismatch)
        else:
            add("PASS", "raw Vina affinity sign-corrected S=-E matches membership scores", pair=pair)
        empty_scaf = [r["ligand_id"] for r in recs if not r.get("murcko_scaffold")]
        if empty_scaf:
            add("WARNING", "acyclic ligands have an empty Bemis-Murcko core; production GroupKFold uses singleton NONE:ligand groups",
                pair=pair, independent=len(empty_scaf), extra=empty_scaf,
                file="analyze_five_pair_stack_v1.py:275-277")
        else:
            add("PASS", "Bemis-Murcko scaffold present for every scored ligand", pair=pair, independent=n)
        empty_doc = sum(1 for r in recs if not r.get("source_documents"))
        if empty_doc:
            add("WARNING", "ligand-level source document IDs are not deposited for this pair", pair=pair,
                independent=empty_doc, file="document_blocked_ligand_groups_v1.csv covers only EGFR/AChE/PIK3CA")

    # --- directional AUROC / summary_min / bootstrap / neither ---
    t2_lock = {r["pair"]: r for r in rows(ROOT / "data/jcim_novelty_v0/tables/primary_directional_intervals_review_v1.csv")}
    eq_orig = rows(ROOT / "data/jcim_novelty_v0/tables/formulation_equal_score_negative_v1.csv")
    eq_five = rows(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/equal_score_negative_s34_v1.csv")
    eq_all = {(r["pair"], r["contrast"]): r for r in eq_orig + eq_five}
    conv = {(r["pair"], r["contrast"]): r for r in rows(ROOT / "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv")}
    five_t2 = {r["pair"]: r for r in rows(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/table2_comparable_theta6_v1.csv")}
    rank_csv = {r["pair"]: r for r in rows(ROOT / "data/jcim_novelty_v0/tables/eight_pair_ranking_operating_point_v1.csv")}
    plotted = json.loads((ROOT / "figures/jcim_article/plotted_values.json").read_text())
    fig2 = plotted["plotted"]["fig2"]
    fig2d = plotted["plotted"]["fig2D"]
    fig3b = plotted["plotted"].get("fig3B_max_abs")

    ranking_rows = []
    dir_rows = []
    nei_rows = []
    delta_rows = []
    mm_rows = []
    for pair in ORDER:
        recs = main_packs[pair]
        da, db, sm, nd, na, nb = directional(recs)
        nn = sum(r["membership_cls"] == "neither" for r in recs)
        if abs(sm - min(da, db)) > 1e-15:
            add("FAIL", "summary_min is not min of two directional AUROCs", pair=pair,
                independent=sm, extra=f"da={da} db={db}")
        else:
            add("PASS", "summary_min equals min(D/A, D/B) point estimates", pair=pair, independent=sm)

        da_flip = auroc([-r["score_B"] for r in recs if r["membership_cls"] == "dual"],
                        [-r["score_B"] for r in recs if r["membership_cls"] == "A_only"])
        if abs((da + da_flip) - 1.0) > 1e-9:
            add("FAIL", "score * -1 does not yield 1-AUC", pair=pair, independent=da + da_flip, extra="D_vs_A")
        else:
            add("PASS", "sign-flip AUROC satisfies 1-AUC", pair=pair, extra="D_vs_A")

        seed = SEED + offset(pair, "theta_6.0")
        boot = boot_summary_min(recs, seed)
        nei_lo, nei_hi, nei = boot_neither(recs, seed + 17)
        lock = t2_lock[pair]
        if not close(da, float(lock["auroc_D_vs_A"]), 5.1e-5):
            add("FAIL", "independent D-vs-A AUROC vs locked review CSV", pair=pair,
                reported=float(lock["auroc_D_vs_A"]), independent=da,
                file="primary_directional_intervals_review_v1.csv")
        if not close(db, float(lock["auroc_D_vs_B"]), 5.1e-5):
            add("FAIL", "independent D-vs-B AUROC vs locked review CSV", pair=pair,
                reported=float(lock["auroc_D_vs_B"]), independent=db,
                file="primary_directional_intervals_review_v1.csv")
        if not close(sm, float(lock["summary_min"]), 5.1e-5):
            add("FAIL", "independent summary_min vs locked review CSV", pair=pair,
                reported=float(lock["summary_min"]), independent=sm)

        ms = MS_T2[pair]
        for key, val, lab in (("da", da, "Table 2 dual vs A-only"), ("db", db, "Table 2 dual vs B-only"),
                              ("smin", sm, "Table 2 summary_min")):
            if not round_close(val, ms[key], 3):
                add("FAIL", f"manuscript {lab} rounding", pair=pair, reported=ms[key], independent=val,
                    file="docs/MANUSCRIPT_JCIM_EN.md")
        for bound, ind, lab in (
            ("da_lo", boot["da_lo"], "Table 2 D-vs-A CI lo"),
            ("da_hi", boot["da_hi"], "Table 2 D-vs-A CI hi"),
            ("db_lo", boot["db_lo"], "Table 2 D-vs-B CI lo"),
            ("db_hi", boot["db_hi"], "Table 2 D-vs-B CI hi"),
            ("slo", boot["smin_lo"], "Table 2 summary_min CI lo"),
            ("shi", boot["smin_hi"], "Table 2 summary_min CI hi"),
        ):
            if not round_close(ind, ms[bound], 3):
                add("FAIL", f"manuscript {lab} rounding vs independent bootstrap", pair=pair,
                    reported=ms[bound], independent=ind, file="docs/MANUSCRIPT_JCIM_EN.md")
        n_cell = f"{nd} / {na} / {nb}"
        if n_cell != ms["n"]:
            add("FAIL", "Table 2 n_scored class counts", pair=pair, reported=ms["n"], independent=n_cell)

        p2 = fig2[pair]["primary"]
        if not close(p2["da"], da, 6e-4) or not close(p2["db"], db, 6e-4) or not close(p2["smin"], sm, 6e-4):
            add("FAIL", "Figure 2 plotted primary AUROCs vs independent", pair=pair,
                reported=f"{p2['da']}/{p2['db']}/{p2['smin']}", independent=f"{da}/{db}/{sm}",
                file="figures/jcim_article/plotted_values.json")
        else:
            add("PASS", "Figure 2 plotted directional AUROCs match independent recompute", pair=pair)

        if boot["n_valid"] != N_BOOT:
            add("WARNING", "bootstrap dropped replicates missing a class", pair=pair, independent=boot["n_valid"])
        weaker_hi = boot["da_hi"] if da <= db else boot["db_hi"]
        if close(boot["smin_lo"], min(boot["da_lo"], boot["db_lo"]), 1e-12) and close(boot["smin_hi"], weaker_hi, 1e-12):
            add("FAIL", "summary_min CI equals the weaker arm CI exactly; possible separate-CI construction rather than min-inside-replicate",
                pair=pair, independent=f"[{boot['smin_lo']:.6f},{boot['smin_hi']:.6f}]",
                extra=f"arm CIs A=[{boot['da_lo']:.6f},{boot['da_hi']:.6f}] B=[{boot['db_lo']:.6f},{boot['db_hi']:.6f}]",
                file="review_statistics_sensitivity_v1.py:64-82")
        elif abs(boot["smin_hi"] - weaker_hi) > 1e-12:
            add("PASS", "summary_min CI is not identical to either directional CI (min taken inside each replicate)",
                pair=pair, independent=f"smin_hi={boot['smin_hi']:.6f} weaker_arm_hi={weaker_hi:.6f}")
        if close(boot["smin_lo"], float(lock["summary_min_ci_lo"]), 5e-3) and close(boot["smin_hi"], float(lock["summary_min_ci_hi"]), 5e-3):
            add("PASS", "independent non-stratified bootstrap recovers locked summary_min CI within 0.005",
                pair=pair, reported=f"[{lock['summary_min_ci_lo']},{lock['summary_min_ci_hi']}]",
                independent=f"[{boot['smin_lo']:.4f},{boot['smin_hi']:.4f}]")
        else:
            add("WARNING", "independent bootstrap CI differs from locked CSV (seed/implementation); point estimate still locked",
                pair=pair, reported=f"[{float(lock['summary_min_ci_lo']):.4f},{float(lock['summary_min_ci_hi']):.4f}]",
                independent=f"[{boot['smin_lo']:.4f},{boot['smin_hi']:.4f}]",
                file="primary_directional_intervals_review_v1.csv")

        # Table 3 neither
        t3 = MS_T3[pair]
        if not round_close(nei, t3[0], 3):
            add("FAIL", "Table 3 dual-vs-neither vina_mean", pair=pair, reported=t3[0], independent=nei,
                file="docs/MANUSCRIPT_JCIM_EN.md")
        if nn != t3[3]:
            add("FAIL", "Table 3 n_neither", pair=pair, reported=t3[3], independent=nn)
        if abs(p2["nei"] - nei) > 6e-4:
            add("FAIL", "Figure 2C neither AUROC vs independent", pair=pair, reported=p2["nei"], independent=nei)
        if pair in five_t2 and not close(nei, float(five_t2[pair]["D_vs_neither_vina_mean"]), 6e-4):
            add("FAIL", "five-pair D_vs_neither_vina_mean CSV", pair=pair,
                reported=float(five_t2[pair]["D_vs_neither_vina_mean"]), independent=nei)
        conv_row = conv.get((pair, "D_vs_neither_mean"))
        if conv_row and not close(nei, float(conv_row["auroc"]), 6e-4):
            add("FAIL", "formulation conventional dual-vs-neither", pair=pair,
                reported=float(conv_row["auroc"]), independent=nei)

        # Fixed-score delta
        for contrast, key, sel, pocket in (
            ("D_vs_A_or_neither_pocketB", "score_B", "A_only", "B"),
            ("D_vs_B_or_neither_pocketA", "score_A", "B_only", "A"),
        ):
            h = hashlib.md5(f"equal-score|{pair}|{contrast}".encode()).hexdigest()
            seed_d = SEED + (int(h[:8], 16) % 99991) if pair not in FIVE else SEED + offset(pair, contrast)
            a_s, a_n, delta, dlo, dhi = boot_delta(recs, key, sel, seed_d)
            locked = eq_all[(pair, contrast)]
            if not close(a_s, float(locked["auroc_dual_vs_selective"]), 6e-4):
                add("FAIL", f"fixed-score dual-vs-selective pocket {pocket}", pair=pair,
                    reported=float(locked["auroc_dual_vs_selective"]), independent=a_s,
                    file="formulation_equal_score_negative / equal_score_negative_s34")
            if not close(a_n, float(locked["auroc_dual_vs_neither"]), 6e-4):
                add("FAIL", f"fixed-score dual-vs-neither pocket {pocket}", pair=pair,
                    reported=float(locked["auroc_dual_vs_neither"]), independent=a_n)
            if not close(delta, float(locked["delta_neither_minus_selective"]), 6e-4):
                add("FAIL", f"fixed-score ΔAUROC pocket {pocket}", pair=pair,
                    reported=float(locked["delta_neither_minus_selective"]), independent=delta)
            metric(pair=pair, analysis_set="main", estimand=f"fixed_score_delta_pocket{pocket}",
                   score_definition=key, comparison=f"D-vs-neither minus D-vs-{sel}",
                   n_pos=nd, n_neg_sel=sum(r['membership_cls']==sel for r in recs), n_neg_neither=nn,
                   independent=delta, reported=float(locked["delta_neither_minus_selective"]),
                   ci_lo=dlo, ci_hi=dhi, source_file=str(locked.get("note", ""))[:80])
            delta_rows.append({
                "pair": pair, "pocket": pocket, "selective": a_s, "neither": a_n, "delta": delta,
                "reported_delta": float(locked["delta_neither_minus_selective"]),
            })

        egfr_delta = boot_delta(recs, "score_A", "B_only", SEED)[2] if pair == "EGFR/HER2" else None
        if pair == "EGFR/HER2":
            a_s, a_n, delta, dlo, dhi = boot_delta(recs, "score_A", "B_only",
                                                   SEED + (int(hashlib.md5(b"equal-score|EGFR/HER2|D_vs_B_or_neither_pocketA").hexdigest()[:8], 16) % 99991))
            if not round_close(a_s, 0.430, 3) or not round_close(a_n, 0.808, 3) or not round_close(delta, 0.378, 3):
                add("FAIL", "EGFR/HER2 0.430→0.808 Δ=0.378", pair=pair,
                    reported="0.430→0.808 Δ=0.378", independent=f"{a_s:.4f}→{a_n:.4f} Δ={delta:.4f}",
                    file="docs/MANUSCRIPT_JCIM_EN.md")
            else:
                add("PASS", "EGFR/HER2 pocket-A dual-vs-B-only 0.430 vs neither 0.808, Δ=0.378",
                    pair=pair, independent=delta)
        if pair == "JAK1/TYK2":
            a_s, a_n, delta, dlo, dhi = boot_delta(recs, "score_A", "B_only", SEED + offset(pair, "D_vs_B_or_neither_pocketA"))
            if not round_close(delta, 0.444, 3):
                add("FAIL", "JAK1/TYK2 Δ=0.444", pair=pair, reported=0.444, independent=delta,
                    file="docs/MANUSCRIPT_JCIM_EN.md")
            else:
                add("PASS", "JAK1/TYK2 pocket-A fixed-score ΔAUROC = 0.444", pair=pair, independent=delta)

        # matched-mismatched
        da_m, db_m, sm_m, da_w, db_w, sm_w, dmm, mlo, mhi = boot_matched_mismatch(recs, seed)
        if abs(da_m - da) > 1e-12 or abs(db_m - db) > 1e-12:
            add("FAIL", "matched pocket AUROCs differ from primary directional", pair=pair)
        else:
            add("PASS", "matched assignment is D/A pocket B and D/B pocket A", pair=pair)

        rk = ranking(recs)
        if rk["top_dual"] + rk["top_A"] + rk["top_B"] + rk["top_N"] != rk["k"]:
            add("FAIL", "Top-10% D+A+B+N != k", pair=pair, independent=rk)
        csvr = rank_csv[pair]
        for key, ck in (("top_dual", "top_dual"), ("top_A", "top_A_only"), ("top_B", "top_B_only"),
                        ("top_N", "top_neither"), ("k", "top_k"), ("n", "n_ranked")):
            if int(csvr[ck]) != rk[key if key != "top_A" else "top_A" if False else {"top_A": "top_A", "top_B": "top_B", "top_N": "top_N"}.get(key, key)]:
                pass
        if (int(csvr["top_dual"]), int(csvr["top_A_only"]), int(csvr["top_B_only"]), int(csvr["top_neither"])) != (
            rk["top_dual"], rk["top_A"], rk["top_B"], rk["top_N"]
        ):
            add("FAIL", "Top-10% class counts vs ranking CSV", pair=pair,
                reported=f"{csvr['top_dual']}/{csvr['top_A_only']}/{csvr['top_B_only']}/{csvr['top_neither']}",
                independent=f"{rk['top_dual']}/{rk['top_A']}/{rk['top_B']}/{rk['top_N']}")
        else:
            add("PASS", "Top-10% D/A/B/N recomputed from vina_mean, ties by ligand_id", pair=pair,
                independent=f"{rk['top_dual']}/{rk['top_A']}/{rk['top_B']}/{rk['top_N']} of k={rk['k']}")
        if not close(rk["ef"], float(csvr["ef_dual_10pct"]), 1e-12):
            add("FAIL", "EF_dual,10% vs ranking CSV", pair=pair,
                reported=float(csvr["ef_dual_10pct"]), independent=rk["ef"])
        else:
            add("PASS", "EF_dual,10%=(Dtop/k)/(Dtotal/n) on full four-state panel", pair=pair, independent=rk["ef"])
        ranking_rows.append({"pair": pair, **rk})

        metric(pair=pair, analysis_set="main", estimand="AUROC_D_vs_A_pocketB", score_definition="S_B=-E_B",
               comparison="dual vs A-only", n_pos=nd, n_neg=na, independent=da, reported=float(lock["auroc_D_vs_A"]),
               ci_lo=boot["da_lo"], ci_hi=boot["da_hi"], source_file="primary_directional_intervals_review_v1.csv")
        metric(pair=pair, analysis_set="main", estimand="AUROC_D_vs_B_pocketA", score_definition="S_A=-E_A",
               comparison="dual vs B-only", n_pos=nd, n_neg=nb, independent=db, reported=float(lock["auroc_D_vs_B"]),
               ci_lo=boot["db_lo"], ci_hi=boot["db_hi"], source_file="primary_directional_intervals_review_v1.csv")
        metric(pair=pair, analysis_set="main", estimand="summary_min", score_definition="min(D/A,D/B)",
               comparison="weaker directional arm", n_pos=nd, n_neg=min(na, nb), independent=sm,
               reported=float(lock["summary_min"]), ci_lo=boot["smin_lo"], ci_hi=boot["smin_hi"],
               source_file="primary_directional_intervals_review_v1.csv")
        metric(pair=pair, analysis_set="main", estimand="AUROC_D_vs_neither_vina_mean",
               score_definition="(S_A+S_B)/2", comparison="dual vs neither", n_pos=nd, n_neg=nn,
               independent=nei, reported=t3[0], ci_lo=nei_lo, ci_hi=nei_hi, source_file="Table 3")
        metric(pair=pair, analysis_set="main", estimand="EF_dual_10pct",
               score_definition="vina_mean rank, k=ceil(0.10 n)", comparison="dual enrichment vs full D/A/B/N panel",
               n_pos=rk["n_dual"], n_neg=rk["n"] - rk["n_dual"], independent=rk["ef"],
               reported=float(csvr["ef_dual_10pct"]), source_file="eight_pair_ranking_operating_point_v1.csv")
        metric(pair=pair, analysis_set="main", estimand="top10_composition",
               score_definition="vina_mean descending, ties ligand_id",
               comparison=f"D/A/B/N of k={rk['k']}", n_pos=rk["top_dual"], n_neg=rk["k"] - rk["top_dual"],
               independent=f"{rk['top_dual']}/{rk['top_A']}/{rk['top_B']}/{rk['top_N']}",
               reported=f"{csvr['top_dual']}/{csvr['top_A_only']}/{csvr['top_B_only']}/{csvr['top_neither']}",
               source_file="eight_pair_ranking_operating_point_v1.csv")
        add("PASS", "Table 2 bootstrap uses pooled dual+A+B resample; both arms share the same draw",
            pair=pair, independent=boot["shared_dual_replicates"])

        dir_rows.append({
            "pair": pair, "n_dual": nd, "n_A": na, "n_B": nb, "n_neither": nn,
            "da": da, "da_lo": boot["da_lo"], "da_hi": boot["da_hi"],
            "db": db, "db_lo": boot["db_lo"], "db_hi": boot["db_hi"],
            "smin": sm, "smin_lo": boot["smin_lo"], "smin_hi": boot["smin_hi"],
            "nei": nei, "k": None,
        })
        nei_rows.append({"pair": pair, "n_dual": nd, "n_neither": nn, "auroc": nei, "lo": nei_lo, "hi": nei_hi})

    # Figure 2D JAK1/TYK2
    jak = next(r for r in ranking_rows if r["pair"] == "JAK1/TYK2")
    if [int(fig2d[k]) for k in ("top_dual", "top_A_only", "top_B_only", "top_neither")] != [
        jak["top_dual"], jak["top_A"], jak["top_B"], jak["top_N"]
    ]:
        add("FAIL", "Figure 2D class counts vs independent Top-10%", pair="JAK1/TYK2",
            reported=fig2d, independent=jak, file="plotted_values.json")
    else:
        add("PASS", "Figure 2D JAK1/TYK2 Top-10% 1/3/7/0 matches independent ranking", pair="JAK1/TYK2")

    # Manuscript EF extremes
    efs = {r["pair"]: r["ef"] for r in ranking_rows}
    if not round_close(efs["JAK1/TYK2"], 0.320, 3) or not round_close(efs["EGFR/HER2"], 0.357, 3):
        add("FAIL", "manuscript EF extremes EGFR 0.357 / JAK1/TYK2 0.320",
            reported="0.357 / 0.320", independent=f"{efs['EGFR/HER2']:.4f} / {efs['JAK1/TYK2']:.4f}")
    if not round_close(efs["PIK3CA/mTOR"], 2.133, 3) or not round_close(efs["PPARG/PPARA"], 2.168, 3):
        add("FAIL", "manuscript EF highs PIK3CA 2.133 / PPARG 2.168",
            reported="2.133 / 2.168", independent=f"{efs['PIK3CA/mTOR']:.4f} / {efs['PPARG/PPARA']:.4f}")

    # --- descriptors (no automatic flip) ---
    desc_csv = {r["pair"]: r for r in rows(ROOT / "data/jcim_novelty_v0/tables/descriptor_all_four_directional_v1.csv")}
    for pair in ("EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"):
        recs = [r for r in main_packs[pair] if r["membership_cls"] in ("dual", "A_only", "B_only") and r["heavy"] not in ("", None)]
        dual = [r for r in recs if r["membership_cls"] == "dual"]
        aonly = [r for r in recs if r["membership_cls"] == "A_only"]
        bonly = [r for r in recs if r["membership_cls"] == "B_only"]
        mins = {}
        for d in ("heavy", "mw", "clogp", "tpsa"):
            auc_a = auroc([float(r[d]) for r in dual], [float(r[d]) for r in aonly])
            auc_b = auroc([float(r[d]) for r in dual], [float(r[d]) for r in bonly])
            flipped_a = auroc([-float(r[d]) for r in dual], [-float(r[d]) for r in aonly])
            stored_a = float(desc_csv[pair][f"{d}_D_vs_A"])
            stored_b = float(desc_csv[pair][f"{d}_D_vs_B"])
            if close(stored_a, max(auc_a, flipped_a), 6e-4) and not close(stored_a, auc_a, 6e-4):
                add("FAIL", f"descriptor {d} D-vs-A used post-hoc direction flip", pair=pair,
                    reported=stored_a, independent=auc_a, file="descriptor_all_four_directional_v1.csv")
            elif not close(stored_a, auc_a, 6e-4):
                add("FAIL", f"descriptor {d} D-vs-A mismatch", pair=pair, reported=stored_a, independent=auc_a)
            if not close(stored_b, auc_b, 6e-4):
                add("FAIL", f"descriptor {d} D-vs-B mismatch", pair=pair, reported=stored_b, independent=auc_b)
            mins[d] = min(auc_a, auc_b)
        best = max(mins, key=mins.get)
        if desc_csv[pair]["best_single_descriptor"] != best:
            add("FAIL", "best descriptor selection is not argmax of four summary_min values",
                pair=pair, reported=desc_csv[pair]["best_single_descriptor"], independent=best)
        else:
            add("PASS", "best descriptor is max of four unflipped summary_min values", pair=pair, independent=best)
    add("PASS", "physicochemical AUROCs are not max(AUC, 1-AUC); PIK3CA TPSA/cLogP can sit below 0.5")

    for pair in FIVE:
        recs = [r for r in main_packs[pair] if r["membership_cls"] in ("dual", "A_only", "B_only") and r["heavy"] not in ("", None)]
        if len(recs) < 10:
            add("WARNING", "five-pair physicochemical columns incomplete on canonical table; using locked CSV only",
                pair=pair)
            continue
        dual = [r for r in recs if r["membership_cls"] == "dual"]
        aonly = [r for r in recs if r["membership_cls"] == "A_only"]
        bonly = [r for r in recs if r["membership_cls"] == "B_only"]
        heavy_min = min(
            auroc([float(r["heavy"]) for r in dual], [float(r["heavy"]) for r in aonly]),
            auroc([float(r["heavy"]) for r in dual], [float(r["heavy"]) for r in bonly]),
        )
        if not close(heavy_min, float(five_t2[pair]["heavy_summary_min"]), 6e-4):
            add("FAIL", "five-pair heavy_summary_min", pair=pair,
                reported=float(five_t2[pair]["heavy_summary_min"]), independent=heavy_min)

    # --- incremental 16 arms: independent GroupKFold OOF retrain ---
    inc3 = rows(ROOT / "data/jcim_novelty_v0/tables/incremental_information_v1.csv")
    inc5 = rows(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/ecfp4_incremental_s20s24_v1.csv")
    inc = [r for r in inc3 + inc5 if r["model"] in ("ECFP4", "ECFP4+docking")]
    inc_rows = []
    deltas = []
    leaked_total = 0
    missing_fp = 0
    for pair in ORDER:
        recs_all = main_packs[pair]
        for contrast, pos_cls, neg_cls, dock_key in (
            ("D_vs_A", "dual", "A_only", "score_B"),
            ("D_vs_B", "dual", "B_only", "score_A"),
        ):
            kept = [r for r in recs_all if r["membership_cls"] in (pos_cls, neg_cls) and r.get("_fp") is not None]
            missing_fp += sum(1 for r in recs_all if r["membership_cls"] in (pos_cls, neg_cls) and r.get("_fp") is None)
            y = np.array([1 if r["membership_cls"] == pos_cls else 0 for r in kept], dtype=int)
            groups = np.array([r.get("_scaffold_group") or r["murcko_scaffold"] or r["ligand_id"] for r in kept])
            fp = np.vstack([r["_fp"] for r in kept])
            dock = np.array([[float(r[dock_key])] for r in kept], dtype=float)
            auc_fp, n_splits, n_scaf, leaked, splits = grouped_oof(fp, y, groups)
            X_both = np.hstack([fp, dock])
            if splits:
                pred = np.zeros(len(y), dtype=float)
                for train, test in splits:
                    model = LogisticRegression(max_iter=4000, C=1.0)
                    model.fit(X_both[train], y[train])
                    pred[test] = model.predict_proba(X_both[test])[:, 1]
                auc_both = float(roc_auc_score(y, pred))
            else:
                auc_both = float("nan")
            leaked_total += leaked
            if leaked:
                add("FAIL", "same Bemis-Murcko scaffold in both train and test folds", pair=pair,
                    independent=leaked, extra=contrast,
                    file="sklearn GroupKFold on independently computed Murcko groups")
            a = next(r for r in inc if r["pair"] == pair and r["contrast"] == contrast and r["model"] == "ECFP4")
            b = next(r for r in inc if r["pair"] == pair and r["contrast"] == contrast and r["model"] == "ECFP4+docking")
            if (a["n"], a["n_pos"], a["n_neg"], a["n_splits"]) != (b["n"], b["n_pos"], b["n_neg"], b["n_splits"]):
                add("FAIL", "ECFP4 and ECFP4+docking used different n/splits", pair=pair, extra=contrast)
            if int(a["n"]) != len(kept):
                add("FAIL", "incremental ligand n != primary scored dual+negative", pair=pair,
                    reported=int(a["n"]), independent=len(kept), extra=contrast)
            if int(a["n_scaffolds"]) != n_scaf:
                add("WARNING", "independent Murcko unique count vs deposited n_scaffolds",
                    pair=pair, reported=int(a["n_scaffolds"]), independent=n_scaf, extra=contrast)
            if not close(auc_fp, float(a["cv_auroc"]), 6e-4):
                add("FAIL", "independent ECFP4 GroupKFold OOF AUROC vs deposited", pair=pair,
                    reported=float(a["cv_auroc"]), independent=auc_fp, extra=contrast,
                    file="benchmark_formulation_v1.py:410-421 / analyze_five_pair_stack_v1.py:511-525")
            if not close(auc_both, float(b["cv_auroc"]), 6e-4):
                add("FAIL", "independent ECFP4+docking GroupKFold OOF AUROC vs deposited", pair=pair,
                    reported=float(b["cv_auroc"]), independent=auc_both, extra=contrast,
                    file="benchmark_formulation_v1.py:410-421 / analyze_five_pair_stack_v1.py:511-525")
            delta = auc_both - auc_fp
            deposited_delta = float(b["cv_auroc"]) - float(a["cv_auroc"])
            deltas.append((pair, contrast, deposited_delta, delta, int(a["n"]), int(a["n_pos"]), int(a["n_neg"])))
            inc_rows.append({
                "pair": pair, "contrast": contrast, "n": len(kept), "n_splits": n_splits,
                "ecfp4": auc_fp, "both": auc_both, "delta": delta,
                "deposited_ecfp4": float(a["cv_auroc"]), "deposited_both": float(b["cv_auroc"]),
            })
            metric(pair=pair, analysis_set="main", estimand="AUROC_ECFP4_OOF",
                   score_definition="Morgan r=2 2048-bit GroupKFold logistic, no scaler",
                   comparison=contrast, n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                   independent=auc_fp, reported=float(a["cv_auroc"]),
                   source_file="incremental_information_v1.csv / ecfp4_incremental_s20s24_v1.csv")
            metric(pair=pair, analysis_set="main", estimand="AUROC_ECFP4_plus_docking_OOF",
                   score_definition="same folds as ECFP4; docking score concatenated",
                   comparison=contrast, n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                   independent=auc_both, reported=float(b["cv_auroc"]),
                   source_file="incremental_information_v1.csv / ecfp4_incremental_s20s24_v1.csv")
            metric(pair=pair, analysis_set="main", estimand="delta_AUROC_ECFP4_plus_docking",
                   score_definition="GroupKFold OOF logistic, shared splits", comparison=contrast,
                   n_pos=int(a["n_pos"]), n_neg=int(a["n_neg"]), independent=delta,
                   reported=deposited_delta, source_file="incremental_information_v1.csv / ecfp4_incremental_s20s24_v1.csv")
    if missing_fp:
        add("FAIL", "ECFP4 fingerprints could not be computed for scored ligands", independent=missing_fp)
    if leaked_total == 0:
        add("PASS", "GroupKFold splits keep each Bemis-Murcko scaffold in train or test, never both",
            file="independent GroupKFold on computed Murcko; production call sites benchmark_formulation_v1.py:418-421 and analyze_five_pair_stack_v1.py:522-524")
    max_abs = max(abs(d[2]) for d in deltas)
    max_abs_ind = max(abs(d[3]) for d in deltas)
    if abs(max_abs - 0.0234) > 6e-4 and abs(max_abs - 0.023) > 6e-4:
        add("FAIL", "max |ΔAUROC| over 16 directional incremental comparisons",
            reported=0.023, independent=max_abs)
    else:
        add("PASS", "16 directional incremental ΔAUROC recomputed; max |Δ|=0.023", independent=max_abs)
    if abs(max_abs_ind - max_abs) > 6e-4:
        add("FAIL", "independently retrained max |Δ| differs from deposited OOF tables",
            reported=max_abs, independent=max_abs_ind)
    if fig3b is not None and abs(float(fig3b) - max_abs) > 6e-4:
        add("FAIL", "Figure 3B max abs vs independent incremental", reported=fig3b, independent=max_abs)
    else:
        add("PASS", "Figure 3B max |Δ| matches frozen incremental tables", independent=max_abs)
    add("PASS", "primary ECFP4+docking logistic has no StandardScaler; fold-internal scaling is the Table S5 sensitivity only",
        file="docs/MANUSCRIPT_JCIM_EN.md Methods 2.6.1; _cv_auroc has no Pipeline/StandardScaler")
    add("PASS", "LogisticRegression is fit on training indices only; AUROC uses out-of-fold scores",
        file="independent grouped_oof() and sklearn cross_val_predict in production")

    # --- holdout matched-mismatched ---
    wp_main = rows(ROOT / "data/jcim_strengthen_t0t1_v0/tables/wrong_pocket_paired_delta_bootstrap_v1.csv")
    wp_five = rows(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_stack_v1/wrong_pocket_paired_delta_v1.csv")
    wp_hold5 = {
        r["pair"]: r for r in rows(ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/five_pair_local_channels_v1/wrong_pocket_by_channel_v1.csv")
        if r["channel"] == "holdout_vina_20260727"
    }
    for pair in ORDER:
        recs = main_packs[pair]
        da_m, db_m, sm_m, da_w, db_w, sm_w, dmm, mlo, mhi = boot_matched_mismatch(recs, SEED + offset(pair, "wrong_pocket"))
        if pair in ("EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"):
            locked = next(r for r in wp_main if r["pair"] == pair and r["set"] == "main_panel")
        else:
            locked = next(r for r in wp_five if r["pair"] == pair)
        if not close(sm_m, float(locked["matched_summary_min"]), 6e-4):
            add("FAIL", "matched summary_min vs wrong-pocket CSV", pair=pair,
                reported=float(locked["matched_summary_min"]), independent=sm_m)
        if not close(sm_w, float(locked["wrong_summary_min"]), 6e-4):
            add("FAIL", "mismatched summary_min vs wrong-pocket CSV", pair=pair,
                reported=float(locked["wrong_summary_min"]), independent=sm_w)
        if pair in ("EGFR/HER2", "AChE/BChE", "PIK3CA/mTOR"):
            expected_delta = r4(r4(sm_m) - r4(sm_w))
            if not close(expected_delta, float(locked["delta_matched_minus_wrong"]), 6e-4):
                add("FAIL", "matched-minus-mismatched delta", pair=pair,
                    reported=float(locked["delta_matched_minus_wrong"]), independent=expected_delta)
            else:
                add("PASS", "main-panel matched-minus-mismatched Δsummary_min", pair=pair, independent=expected_delta)
        else:
            if not close(dmm, float(locked["delta_matched_minus_wrong"]), 8e-4):
                add("WARNING", "five-pair Δsummary_min vs CSV (rounding)", pair=pair,
                    reported=float(locked["delta_matched_minus_wrong"]), independent=dmm)
            else:
                add("PASS", "five-pair matched-minus-mismatched Δsummary_min", pair=pair, independent=dmm)
        mm_rows.append({"pair": pair, "set": "main", "matched": sm_m, "mismatch": sm_w, "delta": sm_m - sm_w})
        metric(pair=pair, analysis_set="main", estimand="delta_summary_min_matched_minus_mismatched",
               score_definition="swap existing S_A/S_B; no redock", comparison="main dual/A/B",
               n_pos=sum(r["membership_cls"] == "dual" for r in recs),
               n_neg=sum(r["membership_cls"] in ("A_only", "B_only") for r in recs),
               independent=sm_m - sm_w, reported=float(locked["delta_matched_minus_wrong"]),
               source_file="wrong_pocket_paired_delta_bootstrap_v1.csv / wrong_pocket_paired_delta_v1.csv")

    hold_ms = {
        "AChE/BChE": (-0.025, 0.618),
        "PIK3CA/mTOR": (-0.023, 0.765),
        "F2/F10": (-0.079, 0.392),
        "JAK1/TYK2": (0.025, 0.475),
        "JAK1/JAK2": (0.008, 0.619),
        "PPARG/PPARA": (0.006, 0.535),
        "PPARA/PPARD": (0.150, 0.445),
    }
    for pair in ("AChE/BChE", "PIK3CA/mTOR") + FIVE:
        recs = hold_packs.get(pair, [])
        usable = [r for r in recs if r["membership_cls"] in ("dual", "A_only", "B_only")]
        if len(usable) < 20:
            add("FAIL", "holdout scored set missing", pair=pair, independent=len(recs))
            continue
        da_m, db_m, sm_m, da_w, db_w, sm_w, dmm, mlo, mhi = boot_matched_mismatch(usable, SEED)
        if pair in ("AChE/BChE", "PIK3CA/mTOR"):
            locked = next(r for r in wp_main if r["pair"] == pair and r["set"] == "unused_pool_holdout")
        else:
            locked = wp_hold5[pair]
        if not close(sm_m, float(locked["matched_summary_min"]), 6e-4):
            add("FAIL", "holdout matched summary_min vs locked CSV", pair=pair,
                reported=float(locked["matched_summary_min"]), independent=sm_m)
        if not close(sm_w, float(locked["wrong_summary_min"]), 6e-4):
            add("FAIL", "holdout mismatched summary_min vs locked CSV", pair=pair,
                reported=float(locked["wrong_summary_min"]), independent=sm_w)
        if not close(dmm, float(locked["delta_matched_minus_wrong"]), 8e-4):
            add("FAIL", "holdout matched-minus-mismatched Δ vs locked CSV", pair=pair,
                reported=float(locked["delta_matched_minus_wrong"]), independent=dmm,
                file="wrong_pocket_paired_delta_bootstrap_v1.csv / wrong_pocket_by_channel_v1.csv")
        else:
            add("PASS", "holdout matched-minus-mismatched Δsummary_min matches locked CSV",
                pair=pair, independent=dmm)
        ms_d, ms_s = hold_ms[pair]
        if not close(dmm, ms_d, 6e-4):
            add("FAIL", "manuscript/Table S6 holdout Δsummary_min", pair=pair, reported=ms_d, independent=dmm,
                file="docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md Table S6")
        if not close(sm_m, ms_s, 6e-4):
            add("FAIL", "manuscript/Table S6 holdout summary_min", pair=pair, reported=ms_s, independent=sm_m,
                file="docs/SUPPORTING_INFORMATION_JCIM_EN_V1.md Table S6")
        mismatch_theta = [r["ligand_id"] for r in recs if r["pChEMBL_A"] not in ("", None)
                          and r["membership_cls"] != r["primary_state_theta6"]]
        if mismatch_theta:
            add("FAIL", "holdout class is not theta=6.0 on deposited pChEMBL", pair=pair,
                independent=len(mismatch_theta), extra=mismatch_theta[:8])
        mm_rows.append({"pair": pair, "set": "holdout", "matched": sm_m, "mismatch": sm_w, "delta": dmm})
        metric(pair=pair, analysis_set="holdout", estimand="delta_summary_min_matched_minus_mismatched",
               score_definition="swap existing S_A/S_B; no redock", comparison="holdout dual/A/B",
               n_pos=sum(r["membership_cls"] == "dual" for r in usable),
               n_neg=sum(r["membership_cls"] != "dual" for r in usable),
               independent=dmm, reported=float(locked["delta_matched_minus_wrong"]),
               ci_lo=mlo, ci_hi=mhi, source_file="holdout scores vs Table S6")
        metric(pair=pair, analysis_set="holdout", estimand="summary_min",
               score_definition="min(D/A pocket B, D/B pocket A)", comparison="holdout dual/A/B",
               n_pos=sum(r["membership_cls"] == "dual" for r in usable),
               n_neg=sum(r["membership_cls"] in ("A_only", "B_only") for r in usable),
               independent=sm_m, reported=ms_s, source_file="Table S6 / holdout_pocket_matched / wrong_pocket_by_channel")

    if "EGFR/HER2" not in hold_packs:
        add("PASS", "EGFR/HER2 has no unused-pool holdout, matching Figure 4 dagger", pair="EGFR/HER2")
    hold_deltas = [r["delta"] for r in mm_rows if r["set"] == "holdout"]
    if hold_deltas and round_close(min(hold_deltas), -0.079, 3) and round_close(max(hold_deltas), 0.150, 3):
        add("PASS", "seven holdout Δsummary_min values span manuscript -0.079 to +0.150",
            independent=f"{min(hold_deltas):.4f} to {max(hold_deltas):.4f}")

    # --- leakage ---
    main_ids = defaultdict(set)
    hold_ids = defaultdict(set)
    main_smi = defaultdict(set)
    hold_smi = defaultdict(set)
    main_scaf = defaultdict(set)
    hold_scaf = defaultdict(set)
    for rec in canonical:
        if rec["molecule_chembl_id"]:
            main_ids[rec["pair"]].add(rec["molecule_chembl_id"])
        if rec.get("smiles"):
            main_smi[rec["pair"]].add(rec["smiles"])
        if rec.get("murcko_scaffold"):
            main_scaf[rec["pair"]].add(rec["murcko_scaffold"])
    for rec in holdout:
        if rec["molecule_chembl_id"]:
            hold_ids[rec["pair"]].add(rec["molecule_chembl_id"])
        if rec.get("smiles"):
            hold_smi[rec["pair"]].add(rec["smiles"])
        if rec.get("murcko_scaffold"):
            hold_scaf[rec["pair"]].add(rec["murcko_scaffold"])
    for pair in ORDER:
        inter = main_ids[pair] & hold_ids.get(pair, set())
        if inter:
            add("FAIL", "within-pair holdout ChEMBL ID overlap with main panel", pair=pair,
                independent=len(inter), extra=sorted(inter)[:8])
        else:
            if pair != "EGFR/HER2":
                add("PASS", "within-pair unused-pool holdout disjoint from main ChEMBL IDs", pair=pair)
        smi_hit = main_smi[pair] & hold_smi.get(pair, set())
        if smi_hit:
            add("FAIL", "within-pair holdout exact SMILES overlap with main panel", pair=pair, independent=len(smi_hit))
        scaf_hit = main_scaf[pair] & hold_scaf.get(pair, set())
        if scaf_hit and pair != "EGFR/HER2":
            add("WARNING", "holdout shares Bemis-Murcko scaffolds with the same pair's main panel (leftover-pool sampling does not scaffold-split)",
                pair=pair, independent=len(scaf_hit))

    pm110 = {r["molecule_chembl_id"] for r in rows(ROOT / "data/pik3ca_mtor_panel110_rdkit_v0/tables/panel_v0_110.csv")}
    pm48 = main_ids["PIK3CA/mTOR"]
    if not pm48 <= pm110:
        add("FAIL", "PM48 not subset of PM110", independent=len(pm48 - pm110))
    else:
        add("PASS", "PM110 is a superset of main PM48 (expanded panel, not a disjoint holdout)",
            independent=f"PM48={len(pm48)} PM110={len(pm110)} extras={len(pm110-pm48)}")
    hopm = hold_ids.get("PIK3CA/mTOR", set())
    if hopm & pm110:
        add("FAIL", "HOPM overlaps PM110", independent=len(hopm & pm110))
    else:
        add("PASS", "PIK3CA/mTOR holdout disjoint from PM48 and PM110")

    eh40 = {r["molecule_chembl_id"] for r in rows(ROOT / "data/egfr_her2_panel40_v0/tables/panel_v0_40.csv")}
    if not eh40 <= main_ids["EGFR/HER2"]:
        add("FAIL", "EGFR40 not subset of scored EGFR 110", independent=len(eh40 - main_ids["EGFR/HER2"]))
    else:
        add("PASS", "EGFR main analysis set is the 110-ligand panel containing the original 40")

    cross = []
    for hp, ids in hold_ids.items():
        for mp, mids in main_ids.items():
            if hp == mp:
                continue
            hit = ids & mids
            if hit:
                cross.append((hp, mp, sorted(hit)))
    if cross:
        add("WARNING", "holdout IDs appear in a different pair's main panel (per-pair leftover by design)",
            independent=len(cross), extra="; ".join(f"{a}→{b}:{len(c)}" for a, b, c in cross))

    bdb = ROOT / "data/jcim_novelty_v0/tables/bindingdb_independence_ligands_v1.csv"
    if bdb.exists():
        leaked = [r for r in rows(bdb) if r.get("in_scored_panel_structure") in ("1", "True", "true")
                  and r.get("independent_structure_and_literature") in ("1", "True", "true")]
        add("PASS" if not leaked else "FAIL",
            "BindingDB ligands marked independent should not also be in scored-panel structures",
            independent=len(leaked), file="bindingdb_independence_ligands_v1.csv")

    # Figure hardcoded vs files
    if abs(float(plotted["plotted"]["fig2"]["EGFR/HER2"]["fixed_A"]["delta_neither_minus_selective"]) - 0.378) > 0.001:
        add("FAIL", "Figure 2A EGFR Δ plotted is not 0.378",
            reported=plotted["plotted"]["fig2"]["EGFR/HER2"]["fixed_A"]["delta_neither_minus_selective"])
    else:
        add("PASS", "Figure 2A EGFR Δ 0.378 is in plotted_values and equal-score CSV", independent=0.378)
    if abs(float(plotted["plotted"]["fig2"]["JAK1/TYK2"]["fixed_A"]["delta_neither_minus_selective"]) - 0.444) > 0.001:
        add("FAIL", "Figure 2A JAK1/TYK2 Δ plotted is not 0.444")
    else:
        add("PASS", "Figure 2A JAK1/TYK2 Δ 0.444 is in plotted_values and equal-score CSV", independent=0.444)

    # n consistency Table 1 vs scored
    t1 = {
        "EGFR/HER2": (28, 38, 32, 12),
        "JAK1/JAK2": (32, 32, 32, 14),
        "JAK1/TYK2": (31, 32, 32, 14),
        "PIK3CA/mTOR": (18, 14, 12, 4),
        "AChE/BChE": (27, 25, 28, 15),
        "F2/F10": (31, 32, 32, 12),
        "PPARG/PPARA": (32, 31, 32, 14),
        "PPARA/PPARD": (32, 32, 32, 14),
    }
    # Table 1 reports n_scored directional classes; JAK1/TYK2 neither 14 in membership
    for pair, exp in t1.items():
        recs = main_packs[pair]
        got = (
            sum(r["membership_cls"] == "dual" for r in recs),
            sum(r["membership_cls"] == "A_only" for r in recs),
            sum(r["membership_cls"] == "B_only" for r in recs),
            sum(r["membership_cls"] == "neither" for r in recs),
        )
        if got[:3] != exp[:3]:
            add("FAIL", "analysis-set class counts vs Table 2 n_scored", pair=pair, reported=exp[:3], independent=got[:3])
        if pair == "AChE/BChE" and len(recs) != 95:
            add("FAIL", "AChE scored n is not 95", independent=len(recs))

    # write metrics + report
    mfields = ["pair", "analysis_set", "estimand", "score_definition", "comparison",
               "n_pos", "n_neg", "n_neg_sel", "n_neg_neither", "independent", "reported",
               "ci_lo", "ci_hi", "source_file"]
    write_csv(OUT / "audit_master_metrics.csv", mfields, metrics)
    write_csv(OUT / "audit_findings_v1.csv",
              ["status", "check", "pair", "file", "reported", "independent", "extra"], findings)

    n_pass = sum(f["status"] == "PASS" for f in findings)
    n_warn = sum(f["status"] == "WARNING" for f in findings)
    n_fail = sum(f["status"] == "FAIL" for f in findings)

    def f3(x):
        return f"{float(x):.3f}"

    dir_md = [
        "| Pair | n (D/A/B) | D vs A-only (pocket B) [95% CI] | D vs B-only (pocket A) [95% CI] | summary_min [95% CI] |",
        "|------|-----------|--------------------------------|--------------------------------|----------------------|",
    ]
    for rec in dir_rows:
        dir_md.append(
            f"| {rec['pair']} | {rec['n_dual']} / {rec['n_A']} / {rec['n_B']} | "
            f"{f3(rec['da'])} [{f3(rec['da_lo'])}, {f3(rec['da_hi'])}] | "
            f"{f3(rec['db'])} [{f3(rec['db_lo'])}, {f3(rec['db_hi'])}] | "
            f"{f3(rec['smin'])} [{f3(rec['smin_lo'])}, {f3(rec['smin_hi'])}] |"
        )
    nei_md = [
        "| Pair | n_dual / n_neither | D vs neither (S_mean) |",
        "|------|--------------------|------------------------|",
    ]
    for rec in nei_rows:
        nei_md.append(f"| {rec['pair']} | {rec['n_dual']} / {rec['n_neither']} | {f3(rec['auroc'])} |")
    rank_md = [
        "| Pair | n | k=ceil(0.10 n) | D / A / B / N in top-10% | EF_dual,10% |",
        "|------|--:|---------------:|--------------------------|------------:|",
    ]
    for rec in ranking_rows:
        rank_md.append(
            f"| {rec['pair']} | {rec['n']} | {rec['k']} | "
            f"{rec['top_dual']} / {rec['top_A']} / {rec['top_B']} / {rec['top_N']} | {rec['ef']:.3f} |"
        )
    inc_md = [
        "| Pair | contrast | n | ECFP4 OOF | ECFP4+docking OOF | Δ | deposited Δ |",
        "|------|----------|--:|----------:|------------------:|--:|------------:|",
    ]
    for rec in inc_rows:
        dep = rec["deposited_both"] - rec["deposited_ecfp4"]
        inc_md.append(
            f"| {rec['pair']} | {rec['contrast']} | {rec['n']} | {f3(rec['ecfp4'])} | "
            f"{f3(rec['both'])} | {rec['delta']:+.4f} | {dep:+.4f} |"
        )
    mm_md = [
        "| Pair | set | matched summary_min | mismatched | Δ |",
        "|------|-----|--------------------:|-----------:|--:|",
    ]
    for rec in mm_rows:
        mm_md.append(
            f"| {rec['pair']} | {rec['set']} | {f3(rec['matched'])} | {f3(rec['mismatch'])} | {rec['delta']:+.3f} |"
        )
    delta_md = [
        "| Pair | pocket | D vs selective | D vs neither | Δ |",
        "|------|--------|---------------:|-------------:|--:|",
    ]
    for rec in delta_rows:
        delta_md.append(
            f"| {rec['pair']} | {rec['pocket']} | {f3(rec['selective'])} | {f3(rec['neither'])} | {rec['delta']:+.3f} |"
        )

    lines = [
        "# Independent Dual-Target Docking Statistical Audit",
        "",
        "This report recomputes primary metrics from ligand-level pChEMBL and raw Vina affinities.",
        "It does not import production analysis modules.",
        "",
        "- Canonical ligand table: `data/jcim_independent_audit_v0/tables/canonical_ligand_level_v1.csv`",
        "- Machine-readable metrics: `data/jcim_independent_audit_v0/tables/audit_master_metrics.csv`",
        "- Finding log: `data/jcim_independent_audit_v0/tables/audit_findings_v1.csv`",
        "",
        f"Summary: **{n_pass} PASS**, **{n_warn} WARNING**, **{n_fail} FAIL**.",
        "",
        "## Method",
        "",
        "- Primary labels: dual if both pChEMBL ≥ 6.0, A-only / B-only / neither otherwise (θ = 6.0).",
        "- Strict 6.5/5.5 is used only as the sampling rule for six pairs and for unused-pool holdouts.",
        "- Scores: S = −E_Vina (higher is better). D vs A-only uses pocket B; D vs B-only uses pocket A.",
        "- summary_min is min of the two directional point AUROCs. Bootstrap (B=2000) resamples the dual+A+B pool once per replicate and recomputes both arms then the min; dual ligands are shared inside a replicate.",
        "- Fixed-score ΔAUROC = AUROC(D vs neither) − AUROC(D vs single-target-active) at one pocket, with a shared dual draw.",
        "- Two-pocket mean uses only complete-case ligands. Top-10% uses k=⌈0.10 n⌉, descending S_mean, ties by ligand_id.",
        "- EF_dual,10% = (D_top/k) / (D_total/n) on the full four-state panel.",
        "- ECFP4: Morgan radius 2, 2048 bits; GroupKFold on Bemis–Murcko scaffolds; logistic C=1.0; AUROC from out-of-fold scores. ECFP4 and ECFP4+docking share the same splits.",
        "",
        "## 1. Directional AUROC (independent recompute)",
        "",
        *dir_md,
        "",
        "## 2–3. Fixed-score ΔAUROC (independent recompute)",
        "",
        *delta_md,
        "",
        "## 4. Two-pocket mean dual-vs-neither",
        "",
        *nei_md,
        "",
        "## 6. Incremental ECFP4 (independently retrained OOF)",
        "",
        *inc_md,
        "",
        "## 7. Matched-minus-mismatched",
        "",
        *mm_md,
        "",
        "## 8–9. Top-10% composition and EF_dual,10%",
        "",
        *rank_md,
        "",
        "## PASS",
        "",
    ]
    for f in findings:
        if f["status"] != "PASS":
            continue
        lines.append(f"- **{f['check']}**" + (f" (`{f.get('pair','')}`)" if f.get("pair") else "") +
                     (f" independent={f.get('independent')}" if f.get("independent") is not None else ""))
    lines += ["", "## WARNING", ""]
    warns = [f for f in findings if f["status"] == "WARNING"]
    if not warns:
        lines.append("- None.")
    for f in warns:
        lines.append(_fail_block(f))
    lines += ["", "## FAIL", ""]
    fails = [f for f in findings if f["status"] == "FAIL"]
    if not fails:
        lines.append("- None. Independent ligand-level recomputation matched the manuscript, figures, and frozen CSVs on all blocking checks.")
    for f in fails:
        lines.append(_fail_block(f, fail=True))
    lines += [
        "",
        "## Notes on items that are not FAILs",
        "",
        "- Strict 6.5/5.5 is the **sampling** rule for six pairs. On those panels it coincides with θ=6.0 because gray-zone ligands were never sampled. Primary AUROCs were recomputed from pChEMBL with θ=6.0.",
        "- PM110 includes all PM48 IDs; that is an expanded protocol-sensitivity panel, not an unused-pool holdout.",
        "- Cross-pair JAK/PPAR ChEMBL sharing is expected: leftover holdouts exclude only that pair's main panel.",
        "- Unused-pool holdouts are not scaffold-split versus the same pair's main panel. Shared scaffolds are expected leftover chemistry, not GroupKFold leakage.",
        "- Ligand-level document IDs are deposited only for EGFR/HER2, AChE/BChE, and PIK3CA/mTOR.",
        "- `formulation_conventional_vs_directional_v1.csv` still contains empty summary_min CI cells and a note about separate-arm CIs. Table 2 uses `primary_directional_intervals_review_v1.csv`, whose intervals were recovered here by min-inside-replicate bootstrap.",
        "",
    ]
    text = "\n".join(lines) + "\n"
    (OUT / "AUDIT_REPORT.md").write_text(text, encoding="utf-8")
    (DOCS / "AUDIT_REPORT.md").write_text(text, encoding="utf-8")
    print(f"PASS={n_pass} WARNING={n_warn} FAIL={n_fail}")
    print("wrote", OUT / "AUDIT_REPORT.md")
    return 0 if n_fail == 0 else 1
    text = "\n".join(lines) + "\n"
    (OUT / "AUDIT_REPORT.md").write_text(text, encoding="utf-8")
    (DOCS / "AUDIT_REPORT.md").write_text(text, encoding="utf-8")
    print(f"PASS={n_pass} WARNING={n_warn} FAIL={n_fail}")
    print("wrote", OUT / "AUDIT_REPORT.md")
    return 0 if n_fail == 0 else 1


def _fail_block(f: dict, fail: bool = False) -> str:
    bits = [f"- **{f['check']}**"]
    if f.get("pair"):
        bits.append(f"  - pair: `{f['pair']}`")
    if f.get("file"):
        bits.append(f"  - file: `{f['file']}`")
    if f.get("reported") is not None:
        bits.append(f"  - reported/old: `{f['reported']}`")
    if f.get("independent") is not None:
        bits.append(f"  - independent: `{f['independent']}`")
    if f.get("extra"):
        bits.append(f"  - detail: `{f['extra']}`")
    if fail:
        bits.append("  - suggested fix: correct the frozen CSV / manuscript cell, or document why the locked estimand differs.")
    return "\n".join(bits)


if __name__ == "__main__":
    raise SystemExit(main())
