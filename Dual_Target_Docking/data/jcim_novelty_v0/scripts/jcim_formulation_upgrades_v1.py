#!/usr/bin/env python3
"""Zero-docking JCIM upgrades on frozen maps and scores.

1. θ = 6.0 four-state census of the J0 candidate pairs (cached mols_*.json).
2. Multivariate property-caliper matching of Dual vs selectives on K=4 scores.
3. AND-score operating-point diagnostic (vina_worst / vina_mean).
4. Ligand-only ECFP4/property AUROC on the *full* ChEMBL maps of the four
   frozen pairs (subsampled; SMILES from ChEMBL, cached). Not a docking
   scale-up and not external validation.

No new docking. Seed 20260729.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "data" / "public_pair_selection"
J0 = ROOT / "data" / "jcim_j0j1_v0" / "tables"
TAB = ROOT / "data" / "jcim_novelty_v0" / "tables"
ANL = ROOT / "data" / "jcim_novelty_v0" / "analysis"
FIG = ROOT / "figures" / "jcim_article"
TAB.mkdir(parents=True, exist_ok=True)
ANL.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

THETA = 6.0
SEED = 20260729
N_BOOT = 2000
N_PER_CLASS = 120
CALIPERS = (0.5, 1.0)
SMILES_CACHE = TAB / "chembl_smiles_cache_theta6_v1.csv"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
UA = "DualFourClass-audit/1.0 (zero-docking formulation upgrade; academic)"

FROZEN_MAPS = {
    "EGFR/HER2": ("EGFR", "HER2"),
    "AChE/BChE": ("ACHE", "BCHE"),
    "PIK3CA/PIK3CB": ("PIK3CA", "PIK3CB"),
    "PIK3CA/mTOR": ("PIK3CA", "MTOR"),
}

SCORE_SPEC = {
    "EGFR/HER2": dict(
        scores="data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        vina_a="3POZ_affinity",
        vina_b="3RCD_affinity",
        panel="data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        panel_key="panel_id",
    ),
    "AChE/BChE": dict(
        scores="data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_ACHE",
        vina_b="vina_BCHE",
        panel=None,
        panel_key=None,
    ),
    "PIK3CA/PIK3CB": dict(
        scores="data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        vina_a="vina_PIK3CA",
        vina_b="vina_PIK3CB",
        panel=None,
        panel_key=None,
    ),
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        panel="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        panel_key="panel_id",
    ),
}


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def load_csv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0].keys()), extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_mols(target: str) -> dict[str, float] | None:
    path = SRC / f"mols_{target}.json"
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        return {k: float(v) for k, v in json.load(handle).items()}


def fourclass(a_map: dict[str, float], b_map: dict[str, float]) -> dict[str, list[str]]:
    buckets = {"dual": [], "A_only": [], "B_only": [], "neither": []}
    for mol in set(a_map) & set(b_map):
        x, y = a_map[mol], b_map[mol]
        if x >= THETA and y >= THETA:
            buckets["dual"].append(mol)
        elif x >= THETA:
            buckets["A_only"].append(mol)
        elif y >= THETA:
            buckets["B_only"].append(mol)
        else:
            buckets["neither"].append(mol)
    for key in buckets:
        buckets[key].sort()
    return buckets


def auroc(pos, neg) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    p = np.asarray(pos, dtype=float)
    n = np.asarray(neg, dtype=float)
    d = p[:, None] - n[None, :]
    return float(((d > 0).sum() + 0.5 * (d == 0).sum()) / (len(p) * len(n)))


def boot_auroc(pos, neg, n_boot=N_BOOT, seed=SEED):
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = [auroc(rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)) for _ in range(n_boot)]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(auroc(pos, neg)), float(lo), float(hi)


def census_rows() -> list[dict]:
    candidates = load_csv(J0 / "j0_candidate_pairs.csv")
    rows = []
    for rec in candidates:
        if rec.get("unordered_dup") == "True":
            continue
        if str(rec.get("notes", "")).startswith("alias_of"):
            continue
        if rec.get("auditable_now") != "True":
            continue
        pair = rec["pair"]
        a_map = load_mols(rec["target_A"])
        b_map = load_mols(rec["target_B"])
        if a_map is None or b_map is None:
            continue
        buckets = fourclass(a_map, b_map)
        n_d, n_a, n_b, n_n = (len(buckets[k]) for k in ("dual", "A_only", "B_only", "neither"))
        n_both = n_d + n_a + n_b + n_n
        directional = int(min(n_d, n_a, n_b) >= 10)
        formulation = int(directional and n_n >= 10)
        rows.append(
            {
                "pair_id": rec["pair_id"],
                "pair": pair.replace("ACHE/BCHE", "AChE/BChE").replace("PIK3CA/MTOR", "PIK3CA/mTOR"),
                "family": rec["family"],
                "metal_enzyme_risk": rec["metal_enzyme_risk"],
                "ppi_risk": rec["ppi_risk"],
                "n_both_measured": n_both,
                "n_dual": n_d,
                "n_A_only": n_a,
                "n_B_only": n_b,
                "n_neither": n_n,
                "directional_n10": directional,
                "formulation_n10": formulation,
                "docked_in_this_paper": int(
                    pair.replace("ACHE/BCHE", "AChE/BChE").replace("PIK3CA/MTOR", "PIK3CA/mTOR")
                    in FROZEN_MAPS
                ),
                "notes": rec.get("notes", ""),
            }
        )
    return rows


def load_scored(pair: str, cfg: dict) -> list[dict]:
    rows = load_csv(ROOT / cfg["scores"])
    smimap = {}
    if cfg["panel"]:
        for rec in load_csv(ROOT / cfg["panel"]):
            smimap[rec[cfg["panel_key"]]] = rec.get("smiles")
    out = []
    for rec in rows:
        a, b = fnum(rec.get(cfg["vina_a"])), fnum(rec.get(cfg["vina_b"]))
        if a is None or b is None:
            continue
        lig = rec.get("ligand") or rec.get("panel_id")
        smi = rec.get("smiles") or smimap.get(lig)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            continue
        out.append(
            {
                "pair": pair,
                "ligand": lig,
                "cls": rec.get("class"),
                "mw": float(Descriptors.MolWt(mol)),
                "clogp": float(Descriptors.MolLogP(mol)),
                "tpsa": float(Descriptors.TPSA(mol)),
                "heavy": float(mol.GetNumHeavyAtoms()),
                "vina_A": -a,
                "vina_B": -b,
                "vina_mean": -(a + b) / 2.0,
                "vina_worst": min(-a, -b),
            }
        )
    return out


def caliper_match(duals: list[dict], negs: list[dict], caliper: float, seed: int):
    feats = ("mw", "clogp", "tpsa", "heavy")
    pool = duals + negs
    mu = np.array([np.mean([r[f] for r in pool]) for f in feats])
    sd = np.array([np.std([r[f] for r in pool], ddof=0) for f in feats])
    sd = np.where(sd < 1e-12, 1.0, sd)

    def z(row):
        return (np.array([row[f] for f in feats]) - mu) / sd

    zd = [z(r) for r in duals]
    zn = [z(r) for r in negs]
    order = np.argsort([r["ligand"] for r in duals])
    used = set()
    kept_d, kept_n = [], []
    for i in order:
        best = None
        best_j = None
        for j, vec in enumerate(zn):
            if j in used:
                continue
            dist = float(np.linalg.norm(zd[i] - vec))
            if dist <= caliper and (best is None or dist < best):
                best = dist
                best_j = j
        if best_j is not None:
            used.add(best_j)
            kept_d.append(duals[i])
            kept_n.append(negs[best_j])
    return kept_d, kept_n


def caliper_table(packs: dict[str, list[dict]]) -> list[dict]:
    rows = []
    for pair, recs in packs.items():
        duals = [r for r in recs if r["cls"] == "dual"]
        for contrast, negs, score_key in (
            ("D_vs_A_pocketB", [r for r in recs if r["cls"] == "A_only"], "vina_B"),
            ("D_vs_B_pocketA", [r for r in recs if r["cls"] == "B_only"], "vina_A"),
        ):
            for caliper in CALIPERS:
                md, mn = caliper_match(duals, negs, caliper, SEED)
                n = min(len(md), len(mn))
                under = int(n < 8)
                if n == 0:
                    pt = lo = hi = float("nan")
                else:
                    pt, lo, hi = boot_auroc(
                        [r[score_key] for r in md],
                        [r[score_key] for r in mn],
                        seed=SEED + int(caliper * 10) + (0 if "A" in contrast else 1),
                    )
                full_pt = auroc(
                    [r[score_key] for r in duals],
                    [r[score_key] for r in negs],
                )
                rows.append(
                    {
                        "pair": pair,
                        "contrast": contrast,
                        "caliper_sd": caliper,
                        "n_dual_matched": len(md),
                        "n_neg_matched": len(mn),
                        "n_dual_full": len(duals),
                        "n_neg_full": len(negs),
                        "auroc_matched": "" if pt != pt else round(pt, 4),
                        "ci_lo": "" if lo != lo else round(lo, 4),
                        "ci_hi": "" if hi != hi else round(hi, 4),
                        "auroc_full": round(full_pt, 4),
                        "delta_matched_minus_full": "" if pt != pt else round(pt - full_pt, 4),
                        "underpowered": under,
                        "note": "1:1 greedy match on z-scored MW/cLogP/TPSA/heavy; Euclidean caliper in SD units",
                    }
                )
    return rows


def and_filter_table(packs: dict[str, list[dict]]) -> list[dict]:
    rows = []
    for pair, recs in packs.items():
        duals = [r for r in recs if r["cls"] == "dual"]
        aonly = [r for r in recs if r["cls"] == "A_only"]
        bonly = [r for r in recs if r["cls"] == "B_only"]
        lib = duals + aonly + bonly
        for score in ("vina_worst", "vina_mean"):
            dual_scores = np.array([r[score] for r in duals], dtype=float)
            percentiles = (10, 25, 50, 75, 90)
            cuts = [float(np.percentile(dual_scores, p)) for p in percentiles]
            for pct, cut in zip(percentiles, cuts):
                passed = [r for r in lib if r[score] >= cut]
                n_pass = len(passed)
                n_d = sum(r["cls"] == "dual" for r in passed)
                n_a = sum(r["cls"] == "A_only" for r in passed)
                n_b = sum(r["cls"] == "B_only" for r in passed)
                precision = n_d / n_pass if n_pass else float("nan")
                recall = n_d / len(duals) if duals else float("nan")
                hardneg = (n_a + n_b) / n_pass if n_pass else float("nan")
                rows.append(
                    {
                        "pair": pair,
                        "score": score,
                        "dual_percentile": pct,
                        "threshold": round(cut, 4),
                        "n_library": len(lib),
                        "n_dual_library": len(duals),
                        "n_pass": n_pass,
                        "n_dual_pass": n_d,
                        "n_A_only_pass": n_a,
                        "n_B_only_pass": n_b,
                        "precision_dual": "" if precision != precision else round(precision, 4),
                        "recall_dual": "" if recall != recall else round(recall, 4),
                        "hardneg_fraction_pass": "" if hardneg != hardneg else round(hardneg, 4),
                        "note": "AND-like dual filter on Dual+A-only+B-only; neither excluded. Not a docking scale-up.",
                    }
                )
    return rows


def load_smiles_cache() -> dict[str, str]:
    if not SMILES_CACHE.exists() or SMILES_CACHE.stat().st_size == 0:
        return {}
    return {r["chembl_id"]: r["smiles"] for r in load_csv(SMILES_CACHE) if r.get("smiles")}


def save_smiles_cache(cache: dict[str, str]) -> None:
    rows = [{"chembl_id": k, "smiles": cache[k]} for k in sorted(cache)]
    write_csv(SMILES_CACHE, rows)


def fetch_smiles(ids: list[str], cache: dict[str, str], batch=50) -> dict[str, str]:
    missing = [i for i in ids if i not in cache]
    if not missing:
        return cache
    for start in range(0, len(missing), batch):
        chunk = missing[start : start + batch]
        query = urllib.parse.urlencode(
            {"molecule_chembl_id__in": ",".join(chunk), "limit": batch}
        )
        url = f"{CHEMBL}?{query}"
        req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=45) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            print("SMILES fetch failed for batch", start, exc, file=sys.stderr)
            time.sleep(1.0)
            continue
        for mol in payload.get("molecules") or []:
            cid = mol.get("molecule_chembl_id")
            smi = (mol.get("molecule_structures") or {}).get("canonical_smiles") or ""
            if cid and smi:
                cache[cid] = smi
        time.sleep(0.1)
        if start == 0 or (start // batch) % 5 == 0:
            print(f"SMILES fetch {min(start + batch, len(missing))}/{len(missing)} missing", flush=True)
    return cache


def subsample(ids: list[str], k: int, seed: int) -> list[str]:
    if len(ids) <= k:
        return list(ids)
    rng = np.random.default_rng(seed)
    return sorted(rng.choice(ids, size=k, replace=False).tolist())


def morgan(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    try:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        scaffold = None
    return mol, fp, scaffold or f"__fail_{smi[:12]}"


def oof_auroc(X, y, groups) -> float:
    n_groups = len(set(groups))
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    splits = min(5, n_pos, n_neg, n_groups)
    if splits < 2:
        return float("nan")
    model = LogisticRegression(C=1.0, max_iter=2000, solver="liblinear")
    gkf = GroupKFold(n_splits=splits)
    try:
        proba = cross_val_predict(model, X, y, groups=groups, cv=gkf, method="predict_proba")[:, 1]
    except ValueError:
        return float("nan")
    if len(set(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, proba))


def ligand_only_table(cache: dict[str, str]) -> tuple[list[dict], dict[str, str]]:
    rows = []
    needed = []
    labeled = {}
    for pair, (ta, tb) in FROZEN_MAPS.items():
        a_map, b_map = load_mols(ta), load_mols(tb)
        buckets = fourclass(a_map, b_map)
        sampled = {
            cls: subsample(
                ids,
                N_PER_CLASS,
                SEED + (int(hashlib.sha256(f"{pair}|{cls}".encode()).hexdigest()[:8], 16) % 99991),
            )
            for cls, ids in buckets.items()
        }
        labeled[pair] = sampled
        for ids in sampled.values():
            needed.extend(ids)
    cache = fetch_smiles(sorted(set(needed)), cache)
    save_smiles_cache(cache)

    for pair, sampled in labeled.items():
        recs = {cls: [] for cls in sampled}
        for cls, ids in sampled.items():
            for cid in ids:
                smi = cache.get(cid)
                if not smi:
                    continue
                mol, fp, scaffold = morgan(smi)
                if mol is None:
                    continue
                recs[cls].append(
                    {
                        "id": cid,
                        "fp": fp,
                        "scaffold": scaffold,
                        "mw": float(Descriptors.MolWt(mol)),
                        "clogp": float(Descriptors.MolLogP(mol)),
                        "tpsa": float(Descriptors.TPSA(mol)),
                        "heavy": float(mol.GetNumHeavyAtoms()),
                    }
                )

        def pack(pos, neg):
            X_fp = np.array([r["fp"] for r in pos + neg])
            X_prop = np.array([[r["mw"], r["clogp"], r["tpsa"], r["heavy"]] for r in pos + neg], dtype=float)
            y = np.array([1] * len(pos) + [0] * len(neg), dtype=int)
            groups = np.array([r["scaffold"] for r in pos + neg])
            return oof_auroc(X_fp, y, groups), oof_auroc(X_prop, y, groups)

        contrasts = (
            ("D_vs_A", recs["dual"], recs["A_only"]),
            ("D_vs_B", recs["dual"], recs["B_only"]),
            ("D_vs_neither", recs["dual"], recs["neither"]),
            ("D_vs_all_nondual", recs["dual"], recs["A_only"] + recs["B_only"] + recs["neither"]),
        )
        da_fp = db_fp = float("nan")
        for name, pos, neg in contrasts:
            fp_auc, prop_auc = pack(pos, neg) if pos and neg else (float("nan"), float("nan"))
            if name == "D_vs_A":
                da_fp = fp_auc
            if name == "D_vs_B":
                db_fp = fp_auc
            rows.append(
                {
                    "pair": pair,
                    "contrast": name,
                    "n_pos_sampled": len(pos),
                    "n_neg_sampled": len(neg),
                    "n_pos_with_smiles": len(pos),
                    "n_neg_with_smiles": len(neg),
                    "ecfp4_groupkfold_auroc": "" if fp_auc != fp_auc else round(fp_auc, 4),
                    "property_groupkfold_auroc": "" if prop_auc != prop_auc else round(prop_auc, 4),
                    "cap_per_class": N_PER_CLASS,
                    "note": "Full ChEMBL θ=6.0 maps, not the docked panel. Ligand-only; no docking.",
                }
            )
        summary = min(da_fp, db_fp) if da_fp == da_fp and db_fp == db_fp else float("nan")
        rows.append(
            {
                "pair": pair,
                "contrast": "summary_min_ecfp4",
                "n_pos_sampled": len(recs["dual"]),
                "n_neg_sampled": min(len(recs["A_only"]), len(recs["B_only"])),
                "n_pos_with_smiles": len(recs["dual"]),
                "n_neg_with_smiles": min(len(recs["A_only"]), len(recs["B_only"])),
                "ecfp4_groupkfold_auroc": "" if summary != summary else round(summary, 4),
                "property_groupkfold_auroc": "",
                "cap_per_class": N_PER_CLASS,
                "note": "min(D_vs_A, D_vs_B) on ligand-only ECFP4; not a docking result.",
            }
        )
    return rows, cache


def plot_upgrades(census, caliper, and_rows, ligand_rows):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))

    n_dir = sum(int(r["directional_n10"]) for r in census)
    n_form = sum(int(r["formulation_n10"]) for r in census)
    n_dock = sum(int(r["docked_in_this_paper"]) for r in census)
    axes[0].bar(
        ["audited\npairs", "directional\nn≥10", "formulation\nn≥10", "docked\nhere"],
        [len(census), n_dir, n_form, n_dock],
        color=["#4C78A8", "#F58518", "#54A24B", "#E45756"],
    )
    axes[0].set_ylabel("count")
    axes[0].set_title("A  θ=6.0 label supply")
    axes[0].set_ylim(0, max(len(census), 1) + 5)

    for pair, color in zip(
        ("EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"),
        ("#4C78A8", "#F58518", "#54A24B", "#E45756"),
    ):
        sub = [r for r in and_rows if r["pair"] == pair and r["score"] == "vina_worst"]
        if not sub:
            continue
        rec = [float(r["recall_dual"]) for r in sub if r["recall_dual"] != ""]
        prec = [float(r["precision_dual"]) for r in sub if r["precision_dual"] != ""]
        axes[1].plot(rec, prec, marker="o", color=color, label=pair, linewidth=1.2)
    axes[1].set_xlabel("dual recall")
    axes[1].set_ylabel("dual precision")
    axes[1].set_title("B  AND filter (vina_worst)")
    axes[1].set_xlim(0, 1.02)
    axes[1].set_ylim(0, 1.02)
    axes[1].legend(fontsize=7, frameon=False)

    labels, neither, directional = [], [], []
    for pair in ("EGFR/HER2", "AChE/BChE", "PIK3CA/PIK3CB", "PIK3CA/mTOR"):
        n = next((r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "D_vs_neither"), None)
        s = next((r for r in ligand_rows if r["pair"] == pair and r["contrast"] == "summary_min_ecfp4"), None)
        if n and s and n["ecfp4_groupkfold_auroc"] != "" and s["ecfp4_groupkfold_auroc"] != "":
            labels.append(pair.replace("/", "/\n"))
            neither.append(float(n["ecfp4_groupkfold_auroc"]))
            directional.append(float(s["ecfp4_groupkfold_auroc"]))
    x = np.arange(len(labels))
    if len(labels):
        axes[2].bar(x - 0.18, neither, 0.36, label="Dual vs neither", color="#4C78A8")
        axes[2].bar(x + 0.18, directional, 0.36, label="ECFP summary_min", color="#E45756")
        axes[2].set_xticks(x, labels, fontsize=7)
        axes[2].axhline(0.5, color="0.5", linestyle="--", linewidth=0.8)
        axes[2].set_ylim(0.4, 1.02)
        axes[2].set_title("C  ligand-only full maps")
        axes[2].legend(fontsize=7, frameon=False)
    else:
        axes[2].text(0.5, 0.5, "SMILES fetch incomplete", ha="center", va="center")
        axes[2].set_axis_off()

    fig.tight_layout()
    out = FIG / "FigS_formulation_upgrades_v1.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print("wrote", out)


def write_verdict(census, caliper, and_rows, ligand_rows) -> None:
    n_dir = sum(int(r["directional_n10"]) for r in census)
    n_form = sum(int(r["formulation_n10"]) for r in census)
    metal_dir = sum(int(r["directional_n10"]) and r["metal_enzyme_risk"] == "True" for r in census)
    lines = [
        "# Zero-docking JCIM formulation upgrades (v1)",
        "",
        "Not a docking scale-up. Not wet-lab. Not BindingDB docking. Not external validation.",
        "",
        f"- Audited pairs with both cached maps (aliases dropped): **{len(census)}**.",
        f"- θ = 6.0 directional gate (dual/A/B each n≥10): **{n_dir}** pairs.",
        f"- θ = 6.0 formulation gate (also neither n≥10): **{n_form}** pairs.",
        f"- Of directional pairs, metal-enzyme-risk: **{metal_dir}**.",
        "- Docked evaluation remains K = 4.",
        "",
        "Property-caliper matching uses z-scored MW/cLogP/TPSA/heavy on the frozen scored panels.",
        "AND-filter tables score Dual+A-only+B-only libraries at Dual-percentile cuts of `vina_worst`/`vina_mean`.",
        "Ligand-only ECFP uses the full ChEMBL maps of the four frozen pairs, capped at 250/class, scaffold GroupKFold.",
        "",
        "Do not write that docking was evaluated on the census pairs.",
        "Do not write that ligand-only full-map AUROC replaces Table 2.",
    ]
    path = ANL / "FORMULATION_UPGRADES_VERDICT_V1.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", path)


def write_contract() -> None:
    contract = {
        "name": "DualFourClass-Bench",
        "version": "1.0",
        "claim": "four-pair formulation audit of docking-based dual-target recognition",
        "not_claimed": [
            "general dual-target docking benchmark",
            "assay-harmonized ground truth",
            "external validation",
            "wet-lab dual activity",
        ],
        "theta_pchembl": 6.0,
        "primary_score": "AutoDock Vina 1.2.7 mode-1, S = -E",
        "primary_estimands": [
            "AUROC(dual, A-only; S_B)",
            "AUROC(dual, B-only; S_A)",
            "summary_min = min of the two arms",
        ],
        "pairs": list(FROZEN_MAPS),
        "seed": SEED,
        "bootstrap": {"B": N_BOOT, "seed": SEED},
        "tables": {
            "primary_directional": "data/jcim_strengthen_t0t1_v0/tables/unified_threshold_sensitivity_v2.csv",
            "formulation": "data/jcim_novelty_v0/tables/formulation_conventional_vs_directional_v1.csv",
            "theta6_census": "data/jcim_novelty_v0/tables/theta6_pair_census_v1.csv",
            "property_caliper": "data/jcim_novelty_v0/tables/property_caliper_match_v1.csv",
            "and_filter": "data/jcim_novelty_v0/tables/and_filter_operating_point_v1.csv",
            "ligand_only_fullmap": "data/jcim_novelty_v0/tables/ligand_only_fullmap_auroc_v1.csv",
        },
    }
    path = TAB / "DUALFOURCLASS_EVALUATION_CONTRACT_v1.json"
    path.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    print("wrote", path)


def main() -> int:
    census = census_rows()
    write_csv(TAB / "theta6_pair_census_v1.csv", census)

    packs = {pair: load_scored(pair, cfg) for pair, cfg in SCORE_SPEC.items()}
    caliper = caliper_table(packs)
    write_csv(TAB / "property_caliper_match_v1.csv", caliper)
    and_rows = and_filter_table(packs)
    write_csv(TAB / "and_filter_operating_point_v1.csv", and_rows)

    cache = load_smiles_cache()
    ligand_rows, cache = ligand_only_table(cache)
    write_csv(TAB / "ligand_only_fullmap_auroc_v1.csv", ligand_rows)
    save_smiles_cache(cache)

    write_verdict(census, caliper, and_rows, ligand_rows)
    write_contract()
    plot_upgrades(census, caliper, and_rows, ligand_rows)

    n_dir = sum(int(r["directional_n10"]) for r in census)
    n_form = sum(int(r["formulation_n10"]) for r in census)
    print(f"census pairs={len(census)} directional_n10={n_dir} formulation_n10={n_form}")
    print(f"caliper rows={len(caliper)} and_rows={len(and_rows)} ligand_rows={len(ligand_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
