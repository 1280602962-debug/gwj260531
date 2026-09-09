#!/usr/bin/env python3
"""StandardScaler sensitivity for the ECFP4 vs ECFP4+docking increment.

Reuses the published GroupKFold rule, ECFP4 (radius 2, 2048 bits), and
LogisticRegression(C=1.0, max_iter=4000). The unscaled arm must reproduce
incremental_information_v1.csv and ecfp4_incremental_s20s24_v1.csv.

The scaled arm fits StandardScaler on each training fold only
(Pipeline + GroupKFold). It does not replace the unscaled Table S5 primary.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_novelty_v0" / "tables"
UNIV = ROOT / "data" / "jcim_chembl_universe_v0"
LOCAL = UNIV / "local_track_b_v0"
PUB3 = ROOT / "data" / "jcim_novelty_v0" / "tables" / "incremental_information_v1.csv"
PUB5 = LOCAL / "tables" / "five_pair_stack_v1" / "ecfp4_incremental_s20s24_v1.csv"

SPEC3 = {
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
    "PIK3CA/mTOR": dict(
        scores="data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        vina_a="4L23_affinity",
        vina_b="4JT6_affinity",
        panel="data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        panel_key="panel_id",
    ),
}

PAIRS5 = [
    {
        "pair": "F2/F10",
        "panel": UNIV / "tables" / "track_b_panels" / "panel_F2_F10_v1.csv",
        "target_a": "4UDW",
        "target_b": "2JKH",
    },
    {
        "pair": "JAK1/TYK2",
        "panel": UNIV / "tables" / "track_b_panels" / "panel_JAK1_TYK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "3LXP",
    },
    {
        "pair": "JAK1/JAK2",
        "panel": UNIV / "tables" / "track_b_panels" / "panel_JAK1_JAK2_v1.csv",
        "target_a": "6N7A",
        "target_b": "8BXH",
    },
    {
        "pair": "PPARG/PPARA",
        "panel": UNIV / "tables" / "track_b_panels" / "panel_PPARG_PPARA_v1.csv",
        "target_a": "9V8H",
        "target_b": "6LXA",
    },
    {
        "pair": "PPARA/PPARD",
        "panel": UNIV / "tables" / "track_b_panels" / "panel_PPARA_PPARD_v1.csv",
        "target_a": "6LXA",
        "target_b": "5U3Q",
    },
]


def fnum(v):
    try:
        return float(v) if v not in ("", None) else None
    except (TypeError, ValueError):
        return None


def load_csv(path: Path):
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
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


def largest_fragment(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def morgan(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    return m, fp


def murcko(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return None


def assemble3(pair: str, cfg: dict) -> list[dict]:
    rows = load_csv(ROOT / cfg["scores"])
    smimap = {}
    if cfg["panel"]:
        for r in load_csv(ROOT / cfg["panel"]):
            smimap[r[cfg["panel_key"]]] = r.get("smiles")
    out = []
    for r in rows:
        a, b = fnum(r.get(cfg["vina_a"])), fnum(r.get(cfg["vina_b"]))
        if a is None or b is None:
            continue
        lig = r.get("ligand") or r.get("panel_id")
        smi = r.get("smiles") or smimap.get(lig)
        if not smi:
            continue
        mol, fp = morgan(smi)
        if mol is None:
            continue
        out.append(
            {
                "pair": pair,
                "ligand": lig,
                "cls": r.get("class"),
                "vina_A": -a,
                "vina_B": -b,
                "fp": fp,
                "scaffold": murcko(smi) or f"__fail_{lig}",
            }
        )
    return out


def load_scores5():
    rows = load_csv(LOCAL / "tables" / "scores_vina_mode1_v1.csv")
    out = {}
    for r in rows:
        out.setdefault(r["pair"], {}).setdefault(r["ligand"], {})[r["target"]] = float(r["score_S"])
    return out


def assemble5(spec, scores) -> list[dict]:
    panel = load_csv(spec["panel"])
    sc = scores.get(spec["pair"], {})
    recs = []
    for row in panel:
        lig = row["panel_id"]
        mol = largest_fragment(row["canonical_smiles"])
        if mol is None:
            continue
        m = sc.get(lig, {})
        sa = m.get(spec["target_a"])
        sb = m.get(spec["target_b"])
        if sa is None or sb is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        except Exception:
            scaf = ""
        recs.append(
            {
                "pair": spec["pair"],
                "ligand": lig,
                "cls": row["theta6_class"],
                "vina_A": sa,
                "vina_B": sb,
                "fp": fp,
                "scaffold": scaf or f"NONE:{lig}",
            }
        )
    return recs


def _estimator(scaled: bool):
    lr = LogisticRegression(max_iter=4000, C=1.0)
    if not scaled:
        return lr
    return Pipeline([("scaler", StandardScaler()), ("lr", lr)])


def _cv_auroc(X, y, groups, scaled: bool):
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    if n_pos < 6 or n_neg < 6 or len(y) < 16:
        return float("nan"), 0
    n_scaf = len(set(groups))
    n_splits = min(5, n_scaf, n_pos, n_neg)
    if n_splits < 2:
        return float("nan"), 0
    cv = GroupKFold(n_splits=n_splits)
    try:
        prob = cross_val_predict(
            _estimator(scaled), X, y, cv=cv, groups=groups, method="predict_proba"
        )[:, 1]
        return float(roc_auc_score(y, prob)), n_splits
    except Exception:
        return float("nan"), n_splits


def incremental_rows(packs, stack: str, scaled: bool):
    rows = []
    scaling = "standardscaler_train_fold" if scaled else "none"
    for pair, recs in packs.items():
        for contrast, pos_cls, neg_cls, dock_key in (
            ("D_vs_A", "dual", "A_only", "vina_B"),
            ("D_vs_B", "dual", "B_only", "vina_A"),
        ):
            kept = [r for r in recs if r["cls"] in (pos_cls, neg_cls)]
            if len(kept) < 16:
                continue
            y = np.array([1 if r["cls"] == pos_cls else 0 for r in kept], dtype=int)
            groups = np.array([r["scaffold"] for r in kept])
            dock = np.array([[r[dock_key]] for r in kept], dtype=float)
            fp = np.vstack([np.asarray(r["fp"], dtype=float) for r in kept])
            rank_dock = auroc(
                [r[dock_key] for r in kept if r["cls"] == pos_cls],
                [r[dock_key] for r in kept if r["cls"] == neg_cls],
            )
            for name, X in (
                ("docking", dock),
                ("ECFP4", fp),
                ("ECFP4+docking", np.hstack([fp, dock])),
            ):
                auc, n_splits = _cv_auroc(X, y, groups, scaled)
                rows.append(
                    {
                        "stack": stack,
                        "pair": pair,
                        "contrast": contrast,
                        "model": name,
                        "scaling": scaling,
                        "n": len(kept),
                        "n_pos": int(y.sum()),
                        "n_neg": int((1 - y).sum()),
                        "n_scaffolds": len(set(groups)),
                        "n_splits": n_splits,
                        "cv_auroc": "" if auc != auc else round(auc, 4),
                        "rank_auroc_docking": round(rank_dock, 4),
                        "note": (
                            "same GroupKFold as published increment; "
                            "StandardScaler fitted on training folds only"
                            if scaled
                            else "unscaled logistic; must match published increment"
                        ),
                    }
                )
    return rows


def published_lookup(path: Path):
    out = {}
    for r in load_csv(path):
        if r["model"] not in {"docking", "ECFP4", "ECFP4+docking"}:
            continue
        key = (r["pair"], r["contrast"], r["model"])
        out[key] = float(r["cv_auroc"])
    return out


def main() -> int:
    packs3 = {pair: assemble3(pair, cfg) for pair, cfg in SPEC3.items()}
    scores5 = load_scores5()
    packs5 = {spec["pair"]: assemble5(spec, scores5) for spec in PAIRS5}

    rows = []
    rows.extend(incremental_rows(packs3, "three_pair", scaled=False))
    rows.extend(incremental_rows(packs3, "three_pair", scaled=True))
    rows.extend(incremental_rows(packs5, "five_pair", scaled=False))
    rows.extend(incremental_rows(packs5, "five_pair", scaled=True))

    pub = {**published_lookup(PUB3), **published_lookup(PUB5)}
    mismatches = []
    for r in rows:
        if r["scaling"] != "none" or r["cv_auroc"] == "":
            continue
        key = (r["pair"], r["contrast"], r["model"])
        got = float(r["cv_auroc"])
        exp = pub[key]
        if abs(got - exp) > 5e-4:
            mismatches.append((key, got, exp))
    if mismatches:
        print("UNSCALED MISMATCH vs published increment:", file=sys.stderr)
        for item in mismatches:
            print(item, file=sys.stderr)
        raise SystemExit(1)

    by = {}
    for r in rows:
        if r["model"] not in {"ECFP4", "ECFP4+docking"} or r["cv_auroc"] == "":
            continue
        by.setdefault((r["scaling"], r["pair"], r["contrast"]), {})[r["model"]] = float(r["cv_auroc"])
    deltas = {}
    for (scaling, pair, contrast), models in by.items():
        deltas.setdefault(scaling, []).append(
            (abs(models["ECFP4+docking"] - models["ECFP4"]), pair, contrast, models)
        )
    max_unscaled = max(deltas["none"])[0]
    max_scaled = max(deltas["standardscaler_train_fold"])[0]
    print(f"unscaled max |Δ| = {max_unscaled:.4f} (published 0.0234)")
    print(f"scaled   max |Δ| = {max_scaled:.4f}")
    for scaling, items in deltas.items():
        worst = max(items, key=lambda t: t[0])
        print(f"  {scaling} worst {worst[1]} {worst[2]} |Δ|={worst[0]:.4f}")

    out_path = OUT / "ecfp4_docking_scaler_sensitivity_v1.csv"
    write_csv(out_path, rows)
    print(f"wrote {out_path} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
