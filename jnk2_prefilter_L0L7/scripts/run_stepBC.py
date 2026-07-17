#!/usr/bin/env python3
"""
Step B (tighten watch, rebuild L7b) + Step C (Glide/AF3 handoff export).
Reads L6 tracks + re-derives Step A heuristic keep for the full novel pool.
Does NOT run docking.

Run: python3 run_stepBC.py
"""
import csv, json
from pathlib import Path
from collections import Counter

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, FilterCatalog, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/mnt/d/CADD paper exercise/JNK2/chembl_amine_pipeline/prefilter_L0L7")
L7 = ROOT / "L7"
L7B = ROOT / "L7b"
STEPC = ROOT / "handoff_glide_af3"
L7B.mkdir(exist_ok=True)
STEPC.mkdir(exist_ok=True)

# ---------- reactivity split (Step B.1/2) ----------
HARD_BAD = {
    "haloacetamide": "[NX3]C(=O)C[Cl,Br,I]",
    "vinyl_sulfone": "C=CS(=O)(=O)",
    "maleimide": "O=C1C=CC(=O)N1",
    "isothiocyanate": "N=C=S",
    "alpha_halo_carbonyl": "[Cl,Br,I][CX4]C(=O)",
    "michael_vinylketone": "[CX3]=[CX3]C(=O)[!N;!O]",
    "beta_lactam": "O=C1CCN1",
    "epoxide": "C1OC1",
    "aldehyde": "[CX3H1](=O)[#6]",
    "acyl_halide": "C(=O)[Cl,Br,I]",
}
SOFT_WATCH = {
    "aromatic_acrylamide": "[c][NX3]C(=O)[C]=[C]",
    "extended_conjugation": "[c]C=CC(=O)[NX3]",
}
HARD_BAD_M = {k: Chem.MolFromSmarts(v) for k, v in HARD_BAD.items()}
SOFT_WATCH_M = {k: Chem.MolFromSmarts(v) for k, v in SOFT_WATCH.items()}


def classify(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return "hard_bad", ["parse_fail"]
    hb = [k for k, p in HARD_BAD_M.items() if p and m.HasSubstructMatch(p)]
    if hb:
        return "hard_bad", hb
    sw = [k for k, p in SOFT_WATCH_M.items() if p and m.HasSubstructMatch(p)]
    if sw:
        return "soft_watch", sw
    return "ok", []


# ---------- Step A heuristic (for full novel pool) ----------
ACRYLAMIDE = Chem.MolFromSmarts("[NX3]C(=O)[C]=[C]")
HINGE_LIKE = [Chem.MolFromSmarts(s) for s in [
    "c1nccc(n1)[NH,NH2]", "c1nccnc1[NH,NH2]", "n1ncc([NH,NH2])c1",
    "c1ccc2[nH]ccc2c1", "c1ccc2ncccc2c1", "c1nc2ccccc2[nH]1",
    "[nH]1ccc2ncccc12", "c1ncnc2[nH]ccc12", "c1ccncc1", "c1nc[nH]c1", "c1n[nH]c2ccccc12",
]]
HINGE_LIKE = [x for x in HINGE_LIKE if x]
AROM_HETERO = Chem.MolFromSmarts("[a;!#6]")
BAD_A = {
    "azo": Chem.MolFromSmarts("N=N"),
    "catechol": Chem.MolFromSmarts("c1ccc(O)c(O)c1"),
    "peptide_like": Chem.MolFromSmarts("[NX3]C(=O)[NX3]C(=O)[NX3]"),
}
_par = FilterCatalog.FilterCatalogParams()
_par.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
PAINS = FilterCatalog.FilterCatalog(_par)


def _phenol_oh(mol):
    p = Chem.MolFromSmarts("[c]O[H]")
    return len(mol.GetSubstructMatches(p)) if p else 0


def _acryl_dist(mol):
    ms = mol.GetSubstructMatches(ACRYLAMIDE)
    if not ms:
        return None, 0
    from rdkit.Chem import rdmolops
    dm = rdmolops.GetDistanceMatrix(mol)
    n_idx = ms[0][0]
    het = [i for i, a in enumerate(mol.GetAtoms()) if a.GetIsAromatic() and a.GetAtomicNum() not in (6, 1)]
    arom = [i for i, a in enumerate(mol.GetAtoms()) if a.GetIsAromatic()]
    tgt = het if het else arom
    return (min((int(dm[n_idx, t]) for t in tgt), default=99)), len(ms)


def stepA_keep(row):
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        return "no"
    hinge_hits = int(row.get("hinge_hits", 0))
    hinge_s = sum(1 for p in HINGE_LIKE if mol.HasSubstructMatch(p))
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    het = len(mol.GetSubstructMatches(AROM_HETERO)) if AROM_HETERO else 0
    r1 = (hinge_hits >= 1) or (hinge_s >= 1) or (arom >= 2 and het >= 1)
    min_d, n_acryl = _acryl_dist(mol)
    rb = Descriptors.NumRotatableBonds(mol)
    r2 = (n_acryl == 1) and (min_d is not None) and (min_d <= 8) and (rb <= 12)
    pains = PAINS.HasMatch(mol) or int(row.get("pains_flag", 0)) == 1
    bad = [k for k, p in BAD_A.items() if p and mol.HasSubstructMatch(p)]
    if _phenol_oh(mol) >= 3:
        bad.append("polyphenol")
    r3 = (not pains) and (not bad)
    tc = float(row["max_tc_core"]); erg = float(row["erg_max"])
    r4 = not (tc < 0.15 and erg >= 0.78 and not r1)
    risky = float(row["sa"]) > 5
    if not r3 or n_acryl != 1:
        return "no"
    if not r1 and not r4:
        return "no"
    if not r1 and hinge_hits == 0 and tc < 0.12:
        return "no"
    if r1 and r2 and r3 and r4 and not risky:
        return "yes"
    if r1 and r2 and r3 and r4 and risky:
        return "unsure"
    return "unsure"


def reclassify(df):
    cats, reasons = [], []
    for smi in df["smiles"]:
        c, r = classify(smi)
        cats.append(c); reasons.append("|".join(r))
    df = df.copy()
    df["react2"] = cats
    df["react2_reason"] = reasons
    return df


TARGET = {"sim_yl": 3000, "sim_56d": 1200, "novel": None, "pan": 300}
CAP = {"sim_yl": 4200, "sim_56d": 1800, "novel": 3500, "pan": 300}


def take(df, floor, cap):
    """ok first, then soft_watch; keep original order within group. Drop hard_bad."""
    df = df[df["react2"] != "hard_bad"]
    ok = df[df["react2"] == "ok"]
    soft = df[df["react2"] == "soft_watch"]
    picked = ok.copy()
    if cap is not None and len(picked) > cap:
        picked = picked.iloc[:cap]
    if cap is None or len(picked) < cap:
        need = (cap - len(picked)) if cap else len(soft)
        if need > 0:
            picked = pd.concat([picked, soft.iloc[:need]])
    return picked


def main():
    dfs = {n: reclassify(pd.read_csv(L7 / f"L7_dock_ready_{n}.csv")) for n in ("sim_yl", "sim_56d", "novel", "pan")}

    # audit
    reason_counts = Counter(); cat_counts = Counter(); per_track = {}
    for n, df in dfs.items():
        per_track[n] = df["react2"].value_counts().to_dict()
        for _, row in df.iterrows():
            cat_counts[row["react2"]] += 1
            for r in (row["react2_reason"].split("|") if row["react2_reason"] else []):
                reason_counts[r] += 1
    (ROOT / "watch_reason_counts.json").write_text(
        json.dumps({"by_reason": dict(reason_counts), "by_category": dict(cat_counts), "per_track": per_track}, indent=2),
        encoding="utf-8")

    dropped_rows = []
    delivered = {}

    # sim/pan: tighten watch
    for n in ("sim_yl", "sim_56d", "pan"):
        df = dfs[n]
        hard = df[df["react2"] == "hard_bad"]
        for _, r in hard.iterrows():
            dropped_rows.append({"id": r["id"], "track": n, "reason": "hard_bad:" + r["react2_reason"]})
        picked = take(df, TARGET[n], CAP[n])
        picked_ids = set(picked["id"])
        # record soft dropped beyond cap
        for _, r in df.iterrows():
            if r["id"] not in picked_ids and r["react2"] != "hard_bad":
                dropped_rows.append({"id": r["id"], "track": n, "reason": "over_cap_or_softwatch_excess:" + r["react2"]})
        delivered[n] = picked

    # novel: keep only Step A yes AND not hard_bad (novel already ok reactivity)
    ndf = dfs["novel"]
    ndf = ndf.assign(stepA_keep=[stepA_keep(r) for _, r in ndf.iterrows()])
    novel_keep = ndf[(ndf["stepA_keep"] == "yes") & (ndf["react2"] != "hard_bad")]
    if len(novel_keep) > CAP["novel"]:
        novel_keep = novel_keep.iloc[:CAP["novel"]]
    delivered["novel"] = novel_keep
    for _, r in ndf.iterrows():
        if r["id"] not in set(novel_keep["id"]):
            reason = "stepA_" + r["stepA_keep"] if r["react2"] != "hard_bad" else "hard_bad"
            dropped_rows.append({"id": r["id"], "track": "novel", "reason": reason})
    # hold unsure separately
    ndf[ndf["stepA_keep"] == "unsure"].to_csv(ROOT / "novel_unsure_hold.csv", index=False)

    # write L7b
    for n, df in delivered.items():
        df.to_csv(L7B / f"L7b_dock_ready_{n}.csv", index=False)
    pd.DataFrame(dropped_rows).to_csv(L7B / "L7b_dropped.csv", index=False)

    # summary
    def ratio(df):
        vc = df["react2"].value_counts().to_dict()
        return {"n": len(df), "ok": int(vc.get("ok", 0)), "soft_watch": int(vc.get("soft_watch", 0))}
    l7_sizes = {n: len(dfs[n]) for n in dfs}
    id_sets = {n: set(delivered[n]["id"]) for n in delivered}
    inter = {}
    ks = list(id_sets)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            inter[f"{ks[i]}∩{ks[j]}"] = len(id_sets[ks[i]] & id_sets[ks[j]])
    summary = {
        "L7_sizes": l7_sizes,
        "L7b_delivered": {n: ratio(delivered[n]) for n in delivered},
        "L7b_total": int(sum(len(delivered[n]) for n in delivered)),
        "retention_vs_L7": {n: round(len(delivered[n]) / l7_sizes[n], 3) for n in delivered},
        "targets": TARGET, "caps": CAP,
        "novel_rule": "Step A heuristic keep==yes only; unsure held in novel_unsure_hold.csv; no discard backfill",
        "track_id_intersections": inter,
        "note": "chemotype triage output; NOT activity hits",
    }
    (L7B / "L7b_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    # ---------- Step C: handoff export ----------
    for n, df in delivered.items():
        smi_path = STEPC / f"L7b_{n}.smi"
        with open(smi_path, "w", encoding="utf-8") as f:
            for _, r in df.iterrows():
                f.write(f"{r['smiles']}\t{r['id']}\n")
        man_path = STEPC / f"af3_manifest_{n}.csv"
        with open(man_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["id", "smiles", "warhead_atom_hint", "target", "cys", "template", "track"])
            for _, r in df.iterrows():
                w.writerow([r["id"], r["smiles"], "acrylamide_beta_C(C=C terminal CH2)", "JNK2", 116, "8ELC", n])
    # combined manifest
    allman = STEPC / "af3_manifest_ALL.csv"
    with open(allman, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "smiles", "warhead_atom_hint", "target", "cys", "template", "track"])
        for n, df in delivered.items():
            for _, r in df.iterrows():
                w.writerow([r["id"], r["smiles"], "acrylamide_beta_C(C=C terminal CH2)", "JNK2", 116, "8ELC", n])

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
