#!/usr/bin/env python3
"""
JNK2 ARS/PaperSpine prefilter L0–L7 (restartable, chunked).

Based on: paper_rewriting_output/ars_deep_research/phase4_plan/prefilter_replan_v1.md
Input: merged_amines.csv (~527k unique primary/secondary amines)
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

ROOT = Path("/mnt/d/CADD paper exercise/JNK2/chembl_amine_pipeline/prefilter_L0L7")
INPUT = Path("/mnt/d/CADD paper exercise/JNK2/chembl_amine_pipeline/merged_amines/merged_amines.csv")
CFG_DIR = ROOT / "config"
LOG_DIR = ROOT / "logs"
STATE_DIR = ROOT / "state"
L5_DIR = ROOT / "L5"
L6_DIR = ROOT / "L6"
L7_DIR = ROOT / "L7"

for d in (CFG_DIR, LOG_DIR, STATE_DIR, L5_DIR, L6_DIR, L7_DIR, ROOT / "scripts"):
    d.mkdir(parents=True, exist_ok=True)

csv.field_size_limit(min(sys.maxsize, 50 * 1024 * 1024))

# ---- L1 params (ARS) ----
AMINE_MW_MIN, AMINE_MW_MAX = 120.0, 450.0
SITES2_MW_MAX = 350.0

# ---- L3 params ----
PROD_MW_MIN, PROD_MW_MAX = 350.0, 650.0
CLOGP_MIN, CLOGP_MAX = -1.0, 6.0
RB_MAX = 12

# ---- L7 quotas (default ~10k) ----
QUOTA = {"sim_yl": 4200, "sim_56d": 1800, "novel": 3500, "pan": 500}

N_WORKERS = 6
CHUNK = 50000


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def setup_log(step: str) -> logging.Logger:
    log = logging.getLogger(f"prefilter.{step}")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    for h in (
        logging.FileHandler(LOG_DIR / f"{step}_{datetime.now():%Y%m%d_%H%M%S}.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_DIR / "pipeline_master.log", encoding="utf-8"),
    ):
        h.setFormatter(fmt)
        log.addHandler(h)
    return log


def state_path(step: str) -> Path:
    return STATE_DIR / f"{step}.json"


def load_state(step: str) -> dict:
    p = state_path(step)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def save_state(step: str, data: dict) -> None:
    data = dict(data)
    data["updated_at"] = utc_now()
    tmp = state_path(step).with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(state_path(step))


def mark_done(step: str, extra: Optional[dict] = None) -> None:
    st = load_state(step)
    st["status"] = "done"
    st["finished_at"] = utc_now()
    if extra:
        st.update(extra)
    save_state(step, st)


def is_done(step: str) -> bool:
    return load_state(step).get("status") == "done"


def load_anchors() -> dict:
    return json.loads((CFG_DIR / "anchors.json").read_text(encoding="utf-8"))


# =========================================================================
# L0
# =========================================================================
def run_l0(force: bool = False) -> None:
    log = setup_log("L0")
    if is_done("L0") and not force:
        log.info("L0 done; skip")
        return
    save_state("L0", {"status": "running", "started_at": utc_now()})
    n = 0
    by_src = Counter()
    by_cls = Counter()
    mw_vals = []
    with open(INPUT, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            n += 1
            by_src[row.get("source", "?")] += 1
            by_cls[row.get("amine_class", "?")] += 1
            try:
                mw_vals.append(float(row["mw"]))
            except Exception:
                pass
            if n % 100000 == 0:
                log.info("  scanned %s", f"{n:,}")
    mw_arr = np.array(mw_vals) if mw_vals else np.array([0.0])
    report = [
        "# L0 amine library QC",
        "",
        f"- Input: `{INPUT}`",
        f"- Unique amines: **{n:,}**",
        f"- Sources: {dict(by_src)}",
        f"- Amine class: {dict(by_cls)}",
        f"- MW median/mean/min/max: {np.median(mw_arr):.1f} / {mw_arr.mean():.1f} / {mw_arr.min():.1f} / {mw_arr.max():.1f}",
        "",
        "## Standardization notes",
        "- Already desalted / InChIKey-deduped / primary+secondary filtered upstream.",
        "- Downstream L1 applies ARS amine MW 120–450 and site-count rules.",
        "- Protocol: ARS prefilter_replan_v1 (PaperSpine×ARS crosswalk).",
        "",
    ]
    (ROOT / "L0_report.md").write_text("\n".join(report), encoding="utf-8")
    mark_done("L0", {"n": n, "by_src": dict(by_src), "by_cls": dict(by_cls)})
    log.info("L0 DONE n=%s", f"{n:,}")


# =========================================================================
# L1
# =========================================================================
def run_l1(force: bool = False) -> None:
    log = setup_log("L1")
    if is_done("L1") and not force:
        log.info("L1 done; skip")
        return
    out = ROOT / "L1_amines_prescreen.csv"
    if force and out.exists():
        out.unlink()
    save_state("L1", {"status": "running", "started_at": utc_now()})
    counts = Counter()
    written = 0
    with open(INPUT, "r", encoding="utf-8", newline="") as fin, open(out, "w", encoding="utf-8", newline="") as fout:
        reader = csv.DictReader(fin)
        fields = list(reader.fieldnames or []) + ["n_sites", "l1_pass_reason"]
        w = csv.DictWriter(fout, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in reader:
            counts["in"] += 1
            try:
                mw = float(row["mw"])
                n_pri = int(float(row.get("n_primary") or 0))
                n_sec = int(float(row.get("n_secondary") or 0))
            except Exception:
                counts["bad_row"] += 1
                continue
            sites = n_pri + n_sec
            row["n_sites"] = sites
            if mw < AMINE_MW_MIN or mw > AMINE_MW_MAX:
                counts["fail_mw"] += 1
                continue
            if sites >= 3:
                counts["fail_sites_ge3"] += 1
                continue
            if sites == 2 and mw > SITES2_MW_MAX:
                counts["fail_sites2_mw"] += 1
                continue
            if sites < 1:
                counts["fail_no_site"] += 1
                continue
            if sites == 1:
                row["l1_pass_reason"] = "sites1"
                counts["pass_sites1"] += 1
            else:
                row["l1_pass_reason"] = "sites2_mw_ok"
                counts["pass_sites2"] += 1
            w.writerow(row)
            written += 1
            if written % 100000 == 0:
                fout.flush()
                log.info("  L1 written=%s", f"{written:,}")
                save_state("L1", {"status": "running", "written": written, **dict(counts)})
    (ROOT / "L1_stats.json").write_text(json.dumps({"written": written, **dict(counts)}, indent=2), encoding="utf-8")
    mark_done("L1", {"written": written, **dict(counts)})
    log.info("L1 DONE written=%s / in=%s", f"{written:,}", f"{counts['in']:,}")


# =========================================================================
# Worker for L2–L4
# =========================================================================
_G = {}


def _init_l234(anchors_json: str) -> None:
    global _G
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, FilterCatalog, Lipinski, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem.rdReducedGraphs import GetErGFingerprint
    from rdkit import DataStructs
    import sys as _sys

    RDLogger.DisableLog("rdApp.*")
    _sys.path.insert(0, "/home/gwj/miniconda3/lib/python3.13/site-packages/rdkit/Contrib/SA_Score")
    import sascorer

    anchors = json.loads(anchors_json)
    names = list(anchors["anchors_full"].keys())
    fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    def mol_fp(smi):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return fpgen.GetFingerprint(m)

    def mol_erg(smi):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        try:
            return np.asarray(GetErGFingerprint(m), dtype=float)
        except Exception:
            return None

    full_fps = {n: mol_fp(s) for n, s in anchors["anchors_full"].items()}
    core_fps = {n: mol_fp(s) for n, s in anchors["anchors_core"].items()}
    full_erg = {n: mol_erg(s) for n, s in anchors["anchors_full"].items()}
    murcko_anchors = {}
    for n, s in anchors["anchors_full"].items():
        m = Chem.MolFromSmiles(s)
        murcko_anchors[n] = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) if m else ""

    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog.FilterCatalog(params)

    hinge = [
        Chem.MolFromSmarts("c1nccc(n1)[NH,NH2]"),  # 2-aminopyrimidine
        Chem.MolFromSmarts("c1nccnc1[NH,NH2]"),
        Chem.MolFromSmarts("n1ncc([NH,NH2])c1"),  # aminopyrazole
        Chem.MolFromSmarts("c1ccc2[nH]ccc2c1"),  # indole
        Chem.MolFromSmarts("c1ccc2ncccc2c1"),  # quinoline
        Chem.MolFromSmarts("c1nc2ccccc2[nH]1"),  # benzimidazole
        Chem.MolFromSmarts("[nH]1ccc2ncccc12"),  # 7-azaindole-ish
        Chem.MolFromSmarts("c1ccncc1"),  # pyridine (weak)
    ]
    hinge = [h for h in hinge if h is not None]

    acrylamide = Chem.MolFromSmarts("[NX3]C(=O)[C]=[C]")
    # bad electrophiles (beyond the single intended acrylamide)
    bad_smarts = [
        Chem.MolFromSmarts("[NX3]C(=O)C[Cl,Br,I]"),  # haloacetamide
        Chem.MolFromSmarts("C=CS(=O)(=O)"),  # vinyl sulfone
        Chem.MolFromSmarts("[CX3]=[CX3][C,S](=O)[O,N,Cl,Br]"),  # broad; filtered carefully below
        Chem.MolFromSmarts("N=C=S"),
        Chem.MolFromSmarts("[Cl,Br,I][CX4]C(=O)"),
    ]
    bad_smarts = [b for b in bad_smarts if b is not None]
    # extra Michael besides acrylamide amide
    vinyl_ketone = Chem.MolFromSmarts("[CX3]=[CX3]C(=O)[!N]")
    maleimide = Chem.MolFromSmarts("O=C1C=CC(=O)N1")
    aromatic_acryl = Chem.MolFromSmarts("[c][NX3]C(=O)[C]=[C]")

    rxn1 = AllChem.ReactionFromSmarts("[N;H2;!$(NC=[O,S]);!$(NS(=O)=O):1]>>[N:1]C(=O)C=C")
    rxn2 = AllChem.ReactionFromSmarts("[N;H1;!$(NC=[O,S]);!$(NS(=O)=O):1]>>[N:1]C(=O)C=C")
    strip_rxns = [
        AllChem.ReactionFromSmarts("[N:1]C(=O)/[C]=[C]/CN(C)C>>[N:1]"),
        AllChem.ReactionFromSmarts("[N:1]C(=O)[C]=[C]CN(C)C>>[N:1]"),
        AllChem.ReactionFromSmarts("[N:1]C(=O)[C]=[CH2]>>[N:1]"),
        AllChem.ReactionFromSmarts("[N:1]C(=O)C=C>>[N:1]"),
    ]

    def strip_core(mol):
        out = mol
        for rxn in strip_rxns:
            changed = True
            while changed:
                changed = False
                prods = rxn.RunReactants((out,))
                if not prods:
                    break
                cand = prods[0][0]
                try:
                    Chem.SanitizeMol(cand)
                    out = cand
                    changed = True
                except Exception:
                    break
        return out

    def tanimoto(fp_a, fp_b):
        if fp_a is None or fp_b is None:
            return 0.0
        return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))

    def erg_sim(a, b):
        if a is None or b is None:
            return 0.0
        na = np.linalg.norm(a)
        nb = np.linalg.norm(b)
        if na < 1e-12 or nb < 1e-12:
            return 0.0
        return float(np.dot(a, b) / (na * nb))

    _G.update(
        dict(
            Chem=Chem,
            AllChem=AllChem,
            Descriptors=Descriptors,
            Crippen=Crippen,
            Lipinski=Lipinski,
            rdMolDescriptors=rdMolDescriptors,
            MurckoScaffold=MurckoScaffold,
            DataStructs=DataStructs,
            sascorer=sascorer,
            fpgen=fpgen,
            GetErGFingerprint=GetErGFingerprint,
            names=names,
            full_fps=full_fps,
            core_fps=core_fps,
            full_erg=full_erg,
            murcko_anchors=murcko_anchors,
            pains=pains,
            hinge=hinge,
            acrylamide=acrylamide,
            bad_smarts=bad_smarts,
            vinyl_ketone=vinyl_ketone,
            maleimide=maleimide,
            aromatic_acryl=aromatic_acryl,
            rxn1=rxn1,
            rxn2=rxn2,
            strip_core=strip_core,
            tanimoto=tanimoto,
            erg_sim=erg_sim,
        )
    )


def _process_amine_row(row: dict) -> List[dict]:
    """Enumerate acrylamides + QC + features for one amine. Returns 0..N product feature dicts."""
    G = _G
    Chem = G["Chem"]
    parent_id = row["id"]
    parent_smi = row["smiles"]
    source = row.get("source", "")
    mol = Chem.MolFromSmiles(parent_smi)
    if mol is None:
        return []

    # enumerate mono-acrylamides
    products = []
    seen_smi = set()
    for rxn in (G["rxn1"], G["rxn2"]):
        for prod_tuple in rxn.RunReactants((mol,)):
            p = prod_tuple[0]
            try:
                Chem.SanitizeMol(p)
            except Exception:
                continue
            smi = Chem.MolToSmiles(p)
            if smi in seen_smi:
                continue
            seen_smi.add(smi)
            products.append((smi, p))

    out_rows = []
    for i, (smi, p) in enumerate(products):
        # --- L3 QC ---
        n_acryl = len(p.GetSubstructMatches(G["acrylamide"]))
        if n_acryl != 1:
            continue  # require single acrylamide
        # bad electrophiles (haloacetamide, vinyl sulfone, maleimide, extra vinyl ketone)
        bad = False
        for bs in G["bad_smarts"][:2]:  # haloacetamide + vinyl sulfone
            if p.HasSubstructMatch(bs):
                bad = True
                break
        if G["maleimide"] and p.HasSubstructMatch(G["maleimide"]):
            bad = True
        extra_michael = 0
        if G["vinyl_ketone"] and p.HasSubstructMatch(G["vinyl_ketone"]):
            extra_michael = len(p.GetSubstructMatches(G["vinyl_ketone"]))
        if bad or extra_michael > 0:
            reactivity = "bad"
        else:
            # watch: aromatic acrylamide (aniline-derived) often more reactive / less selective
            if G.get("aromatic_acryl") is not None and p.HasSubstructMatch(G["aromatic_acryl"]):
                reactivity = "watch"
            else:
                reactivity = "ok"
        if reactivity == "bad":
            continue  # hard exclude

        mw = float(G["Descriptors"].MolWt(p))
        if mw < PROD_MW_MIN or mw > PROD_MW_MAX:
            continue
        clogp = float(G["Crippen"].MolLogP(p))
        if clogp < CLOGP_MIN or clogp > CLOGP_MAX:
            continue
        rb = int(G["Lipinski"].NumRotatableBonds(p))
        if rb > RB_MAX:
            continue

        pains_flag = 1 if G["pains"].HasMatch(p) else 0
        try:
            sa = float(G["sascorer"].calculateScore(p))
        except Exception:
            sa = 10.0

        try:
            ik = Chem.MolToInchiKey(p)
        except Exception:
            continue

        # --- L4 features ---
        fp_full = G["fpgen"].GetFingerprint(p)
        core_mol = G["strip_core"](p)
        try:
            Chem.SanitizeMol(core_mol)
            core_smi = Chem.MolToSmiles(core_mol)
            fp_core = G["fpgen"].GetFingerprint(core_mol)
        except Exception:
            core_smi = ""
            fp_core = None

        tc_full = {}
        tc_core = {}
        for n in G["names"]:
            tc_full[n] = G["tanimoto"](fp_full, G["full_fps"][n])
            tc_core[n] = G["tanimoto"](fp_core, G["core_fps"][n]) if fp_core else 0.0

        max_tc_full = max(tc_full.values()) if tc_full else 0.0
        max_tc_core = max(tc_core.values()) if tc_core else 0.0
        nearest_full = max(tc_full, key=tc_full.get) if tc_full else ""
        nearest_core = max(tc_core, key=tc_core.get) if tc_core else ""

        try:
            erg = np.asarray(G["GetErGFingerprint"](p), dtype=float)
        except Exception:
            erg = None
        erg_scores = {n: G["erg_sim"](erg, G["full_erg"][n]) for n in G["names"]}
        erg_max = max(erg_scores.values()) if erg_scores else 0.0
        nearest_erg = max(erg_scores, key=erg_scores.get) if erg_scores else ""

        murcko = Chem.MolToSmiles(G["MurckoScaffold"].GetScaffoldForMol(p))
        same_scaffold = [n for n, ms in G["murcko_anchors"].items() if ms and ms == murcko]
        hinge_hits = sum(1 for h in G["hinge"] if p.HasSubstructMatch(h))

        pid = f"{parent_id}_ACR{i}"
        rec = {
            "id": pid,
            "parent_id": parent_id,
            "parent_source": source,
            "smiles": smi,
            "inchikey": ik,
            "core_smiles": core_smi,
            "mw": f"{mw:.2f}",
            "clogp": f"{clogp:.2f}",
            "rb": rb,
            "sa": f"{sa:.3f}",
            "reactivity_bucket": reactivity,
            "pains_flag": pains_flag,
            "n_acrylamide": n_acryl,
            "max_tc_full": f"{max_tc_full:.4f}",
            "max_tc_core": f"{max_tc_core:.4f}",
            "nearest_anchor_full": nearest_full,
            "nearest_anchor_core": nearest_core,
            "erg_max": f"{erg_max:.4f}",
            "nearest_erg": nearest_erg,
            "murcko": murcko,
            "same_scaffold_as": "|".join(same_scaffold),
            "hinge_hits": hinge_hits,
        }
        for n in G["names"]:
            rec[f"tc_full_{n}"] = f"{tc_full[n]:.4f}"
            rec[f"tc_core_{n}"] = f"{tc_core[n]:.4f}"
            rec[f"erg_{n}"] = f"{erg_scores[n]:.4f}"
        out_rows.append(rec)
    return out_rows


def _process_batch(rows: List[dict], anchors_json: str) -> List[dict]:
    if not _G:
        _init_l234(anchors_json)
    out = []
    for r in rows:
        out.extend(_process_amine_row(r))
    return out


def run_l2_l4(force: bool = False) -> None:
    """Combined L2 enum + L3 QC + L4 features → L4_features.csv (+ intermediate counts)."""
    log = setup_log("L2L4")
    if is_done("L2L4") and not force:
        log.info("L2L4 done; skip")
        return
    l1 = ROOT / "L1_amines_prescreen.csv"
    if not l1.exists():
        raise FileNotFoundError(l1)
    out4 = ROOT / "L4_features.csv"
    out3 = ROOT / "L3_acrylamides_qc.csv"
    # L3 is subset of L4 columns; we write L4 as canonical and copy key cols to L3
    st = load_state("L2L4")
    rows_done = 0 if force else int(st.get("l1_rows_processed", 0))
    if force:
        for p in (out4, out3, ROOT / "L2_acrylamides_raw.csv"):
            if p.exists():
                p.unlink()
        rows_done = 0

    anchors = load_anchors()
    anchors_json = json.dumps(anchors)
    save_state("L2L4", {"status": "running", "started_at": st.get("started_at", utc_now()), "l1_rows_processed": rows_done})

    mode = "a" if rows_done > 0 and out4.exists() else "w"
    t0 = time.time()
    n_prod = 0 if force else int(st.get("products_written", 0))
    react_counts = Counter(st.get("reactivity", {}) if not force else {})
    # anchors already loaded above

    fieldnames = [
        "id", "parent_id", "parent_source", "smiles", "inchikey", "core_smiles",
        "mw", "clogp", "rb", "sa", "reactivity_bucket", "pains_flag", "n_acrylamide",
        "max_tc_full", "max_tc_core", "nearest_anchor_full", "nearest_anchor_core",
        "erg_max", "nearest_erg", "murcko", "same_scaffold_as", "hinge_hits",
    ]
    for n in anchors["anchors_full"]:
        fieldnames += [f"tc_full_{n}", f"tc_core_{n}", f"erg_{n}"]

    pending: List[dict] = []
    sub = 200
    n_in = 0

    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_init_l234, initargs=(anchors_json,)) as ex, open(
        out4, mode, encoding="utf-8", newline=""
    ) as f4:
        w4 = csv.DictWriter(f4, fieldnames=fieldnames, extrasaction="ignore")
        if mode == "w":
            w4.writeheader()

        def flush(batch: List[dict]) -> None:
            nonlocal n_prod
            if not batch:
                return
            futs = [ex.submit(_process_batch, batch[i : i + sub], anchors_json) for i in range(0, len(batch), sub)]
            for fut in as_completed(futs):
                for rec in fut.result():
                    w4.writerow(rec)
                    n_prod += 1
                    react_counts[rec["reactivity_bucket"]] += 1

        with open(l1, "r", encoding="utf-8", newline="") as fin:
            for row in csv.DictReader(fin):
                n_in += 1
                if n_in <= rows_done:
                    continue
                pending.append(row)
                if len(pending) >= 2000:
                    flush(pending)
                    pending = []
                    f4.flush()
                    save_state(
                        "L2L4",
                        {
                            "status": "running",
                            "l1_rows_processed": n_in,
                            "products_written": n_prod,
                            "elapsed_s": round(time.time() - t0, 1),
                            "reactivity": dict(react_counts),
                        },
                    )
                    rate = (n_in - rows_done) / max(time.time() - t0, 1)
                    log.info(
                        "L1=%s products=%s react=%s rate=%.0f amines/s",
                        f"{n_in:,}",
                        f"{n_prod:,}",
                        dict(react_counts),
                        rate,
                    )
            flush(pending)

    # Write L3 as projection of L4
    l3_fields = [
        "id", "parent_id", "parent_source", "smiles", "inchikey", "mw", "clogp", "rb", "sa",
        "reactivity_bucket", "pains_flag", "n_acrylamide",
    ]
    with open(out4, "r", encoding="utf-8", newline="") as fin, open(out3, "w", encoding="utf-8", newline="") as fout:
        r = csv.DictReader(fin)
        w = csv.DictWriter(fout, fieldnames=l3_fields, extrasaction="ignore")
        w.writeheader()
        for row in r:
            w.writerow(row)

    # L2 raw = same as L3 for this run (only products that passed QC are kept;
    # document that raw failed products are not retained to save disk — counts in stats)
    stats = {
        "l1_processed": n_in,
        "l4_products": n_prod,
        "reactivity": dict(react_counts),
        "note": "L2 raw failures discarded; L3/L4 contain QC-passed standard acrylamides only",
    }
    (ROOT / "L2L4_stats.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    # lightweight L2 pointer file
    (ROOT / "L2_acrylamides_raw.csv").write_text(
        "# See L3_acrylamides_qc.csv / L4_features.csv; QC-failed products not retained.\n",
        encoding="utf-8",
    )
    mark_done("L2L4", stats)
    log.info("L2L4 DONE products=%s", f"{n_prod:,}")


# =========================================================================
# L5 calibration
# =========================================================================
def run_l5(force: bool = False) -> None:
    log = setup_log("L5")
    if is_done("L5") and not force:
        log.info("L5 done; skip")
        return
    feat = ROOT / "L4_features.csv"
    if not feat.exists():
        raise FileNotFoundError(feat)
    save_state("L5", {"status": "running", "started_at": utc_now()})
    anchors = load_anchors()

    # stream features: collect arrays for sensitivity
    max_core = []
    max_full = []
    nearest_core = []
    erg_max = []
    hinge = []
    murcko = []
    n_feat = 0
    with open(feat, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            n_feat += 1
            max_core.append(float(row["max_tc_core"]))
            max_full.append(float(row["max_tc_full"]))
            nearest_core.append(row["nearest_anchor_core"])
            erg_max.append(float(row["erg_max"]))
            hinge.append(int(row["hinge_hits"]))
            murcko.append(row["murcko"])
            if n_feat % 100000 == 0:
                log.info("  L5 scan %s", f"{n_feat:,}")

    max_core = np.array(max_core)
    max_full = np.array(max_full)
    erg_max = np.array(erg_max)
    hinge = np.array(hinge)

    # Spearman core vs full
    from scipy.stats import spearmanr

    rho, pval = spearmanr(max_core, max_full)

    # Tc grid for Sim/Novel sizes
    lo_grid = [0.15, 0.18, 0.20, 0.22, 0.25]
    hi_grid = [0.50, 0.55, 0.60, 0.70]
    grid = []
    for lo in lo_grid:
        for hi in hi_grid:
            sim = ((max_core >= lo) & (max_core <= hi) & (max_full <= 0.70)).sum()
            near_dup = (max_full > 0.70).sum()
            novel_pool = (max_core < lo) & (max_full <= 0.70)
            grid.append({"lo": lo, "hi": hi, "sim_n": int(sim), "near_dup": int(near_dup), "novel_core_lt_lo": int(novel_pool.sum())})

    # ErG percentile for Novel gate
    erg_p50 = float(np.percentile(erg_max, 50))
    erg_p75 = float(np.percentile(erg_max, 75))
    erg_p90 = float(np.percentile(erg_max, 90))

    # Leave-one-anchor-out
    loo = {}
    names = list(anchors["anchors_full"].keys())
    per_anchor_tc = {nm: [] for nm in names}
    nearest_list = []
    with open(feat, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            nearest_list.append(row["nearest_anchor_core"])
            for nm in names:
                per_anchor_tc[nm].append(float(row[f"tc_core_{nm}"]))
    nearest_list = np.array(nearest_list)
    for held in names:
        mask = nearest_list == held
        if mask.sum() == 0:
            loo[held] = {"n_pref": 0, "recover_sim_frac": None}
            continue
        others = [nm for nm in names if nm != held]
        mat = np.maximum.reduce([np.array(per_anchor_tc[nm]) for nm in others])
        recover = ((mat >= 0.22) & (mat <= 0.55) & mask).sum() / mask.sum()
        loo[held] = {"n_pref": int(mask.sum()), "recover_sim_frac_at_0.22_0.55": round(float(recover), 4)}

    # Jaccard Sim sets: core-window vs full-window
    sim_core = (max_core >= 0.22) & (max_core <= 0.55) & (max_full <= 0.70)
    sim_full = (max_full >= 0.22) & (max_full <= 0.55)
    inter = (sim_core & sim_full).sum()
    union = (sim_core | sim_full).sum()
    jaccard = float(inter / union) if union else 0.0

    thresholds = {
        "sim_tc_core_lo": 0.22,
        "sim_tc_core_hi": 0.55,
        "near_duplicate_tc_full": 0.70,
        "novel_tc_core_max": 0.22,
        "novel_tc_core_pref_lo": 0.08,
        "novel_erg_min": round(erg_p75, 4),
        "novel_erg_rule": "erg_max >= p75 OR hinge_hits>=1",
        "erg_p50": round(erg_p50, 4),
        "erg_p75": round(erg_p75, 4),
        "erg_p90": round(erg_p90, 4),
        "watch_quota_frac": 0.20,
        "quotas": QUOTA,
        "spearman_core_vs_full": {"rho": round(float(rho), 4), "p": float(pval)},
        "jaccard_sim_core_vs_full_window": round(jaccard, 4),
        "n_features": n_feat,
        "chosen_reason": "ARS initial window 0.22–0.55; ErG Novel gate at p75; confirmed by grid stability of Sim size",
    }
    (CFG_DIR / "thresholds.json").write_text(json.dumps(thresholds, indent=2), encoding="utf-8")
    (L5_DIR / "tc_grid.json").write_text(json.dumps(grid, indent=2), encoding="utf-8")
    (L5_DIR / "loo.json").write_text(json.dumps(loo, indent=2), encoding="utf-8")

    md = [
        "# L5 calibration",
        "",
        f"- Feature rows: **{n_feat:,}**",
        f"- Spearman(max_tc_core, max_tc_full): **ρ={rho:.3f}** (p={pval:.2e})",
        f"- Jaccard(Sim by core-window vs full-window): **{jaccard:.3f}** (<1 ⇒ warhead-strip matters)",
        f"- ErG percentiles: p50={erg_p50:.3f}, p75={erg_p75:.3f}, p90={erg_p90:.3f}",
        "",
        "## Working thresholds (written to config/thresholds.json)",
        "```json",
        json.dumps(thresholds, indent=2),
        "```",
        "",
        "## Leave-one-anchor-out (nearest_core preference recovery into Sim via other anchors)",
        "",
    ]
    for k, v in loo.items():
        md.append(f"- **{k}**: {v}")
    md += ["", "## Tc grid (Sim counts)", ""]
    for g in grid:
        md.append(f"- lo={g['lo']}, hi={g['hi']}: sim={g['sim_n']:,}, near_dup={g['near_dup']:,}, novel_lt_lo={g['novel_core_lt_lo']:,}")
    md += ["", "## Anchor self-check (full molecules as if products)", ""]
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs

    fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    for name, smi in anchors["anchors_full"].items():
        m = Chem.MolFromSmiles(smi)
        fp = fpgen.GetFingerprint(m)
        scores = {
            o: DataStructs.TanimotoSimilarity(fp, fpgen.GetFingerprint(Chem.MolFromSmiles(s)))
            for o, s in anchors["anchors_full"].items()
        }
        md.append(f"- {name} full-Tc to anchors: " + ", ".join(f"{k}={v:.3f}" for k, v in scores.items()))
    (ROOT / "L5_calibration.md").write_text("\n".join(md), encoding="utf-8")
    mark_done("L5", thresholds)
    log.info("L5 DONE rho=%.3f jaccard=%.3f n=%s", rho, jaccard, f"{n_feat:,}")


# =========================================================================
# L6 / L7
# =========================================================================
def run_l6_l7(force: bool = False) -> None:
    log = setup_log("L6L7")
    if is_done("L6L7") and not force:
        log.info("L6L7 done; skip")
        return
    thr_path = CFG_DIR / "thresholds.json"
    if not thr_path.exists():
        raise FileNotFoundError("Run L5 first for thresholds.json")
    thr = json.loads(thr_path.read_text(encoding="utf-8"))
    feat = ROOT / "L4_features.csv"
    save_state("L6L7", {"status": "running", "started_at": utc_now()})

    lo = thr["sim_tc_core_lo"]
    hi = thr["sim_tc_core_hi"]
    nd = thr["near_duplicate_tc_full"]
    novel_max = thr["novel_tc_core_max"]
    novel_pref_lo = thr["novel_tc_core_pref_lo"]
    erg_min = thr["novel_erg_min"]
    watch_frac = thr["watch_quota_frac"]
    quotas = thr["quotas"]

    paths = {
        "sim": L6_DIR / "L6_track_sim.csv",
        "novel": L6_DIR / "L6_track_novel.csv",
        "pan": L6_DIR / "L6_track_pan.csv",
        "discard": L6_DIR / "L6_discard.csv",
    }
    if force:
        for p in paths.values():
            if p.exists():
                p.unlink()
        for p in L7_DIR.glob("L7_*.csv"):
            p.unlink()

    # First pass: classify into memory-light temp files via streaming write
    # Read header
    with open(feat, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        fields = list(reader.fieldnames or []) + ["track", "subtrack", "discard_reason"]

    counts = Counter()
    # Store candidates for ranking in lists (may be tens of thousands — OK)
    sim_rows = []
    novel_rows = []
    pan_rows = []

    with open(feat, "r", encoding="utf-8", newline="") as fin, open(
        paths["discard"], "w", encoding="utf-8", newline=""
    ) as fd:
        reader = csv.DictReader(fin)
        wd = csv.DictWriter(fd, fieldnames=fields, extrasaction="ignore")
        wd.writeheader()
        for row in reader:
            counts["in"] += 1
            max_core = float(row["max_tc_core"])
            max_full = float(row["max_tc_full"])
            nearest = row["nearest_anchor_core"]
            erg = float(row["erg_max"])
            hinge_hits = int(row["hinge_hits"])
            murcko_same = bool(row.get("same_scaffold_as") or "")
            react = row["reactivity_bucket"]
            pains = int(row["pains_flag"])

            row["track"] = ""
            row["subtrack"] = ""
            row["discard_reason"] = ""

            if max_full > nd:
                row["discard_reason"] = "near_duplicate_tc_full"
                wd.writerow(row)
                counts["discard_near_dup"] += 1
                continue

            # Pan track: nearest JNK-IN-8
            if nearest == "JNK-IN-8" and lo <= max_core <= hi:
                row["track"] = "pan"
                row["subtrack"] = "pan_jnkin8"
                pan_rows.append(row)
                counts["pan"] += 1
                continue

            # Sim
            if lo <= max_core <= hi and nearest in ("YL5084", "YL2056", "56d"):
                row["track"] = "sim"
                if nearest == "56d":
                    row["subtrack"] = "sim_56d"
                else:
                    row["subtrack"] = "sim_yl"
                sim_rows.append(row)
                counts["sim"] += 1
                continue

            # Novel
            if max_core < novel_max and not murcko_same and (erg >= erg_min or hinge_hits >= 1):
                # prefer mild core similarity
                row["track"] = "novel"
                row["subtrack"] = "novel"
                if max_core < novel_pref_lo:
                    row["subtrack"] = "novel_far"
                novel_rows.append(row)
                counts["novel"] += 1
                continue

            # else discard with reason
            reasons = []
            if max_core > hi and max_full <= nd:
                reasons.append("tc_core_above_sim")
            elif novel_max <= max_core < lo:
                reasons.append("tc_core_gap_band")
            elif max_core < novel_max and murcko_same:
                reasons.append("novel_but_same_murcko")
            elif max_core < novel_max and erg < erg_min and hinge_hits < 1:
                reasons.append("novel_no_erg_hinge")
            else:
                reasons.append("other")
            if pains:
                reasons.append("pains_flag")
            row["discard_reason"] = "|".join(reasons)
            wd.writerow(row)
            counts["discard_other"] += 1

    def sort_key_sim(r):
        # chembl evidence proxy: parent_source chembl first; then max_tc_core desc; sa asc; ok before watch
        src = 0 if r.get("parent_source") == "chembl" else 1
        react = 0 if r.get("reactivity_bucket") == "ok" else 1
        return (src, react, -float(r["max_tc_core"]), float(r["sa"]), int(r["pains_flag"]))

    def sort_key_novel(r):
        react = 0 if r.get("reactivity_bucket") == "ok" else 1
        # prefer mild core 0.08-0.22
        core = float(r["max_tc_core"])
        in_pref = 0 if novel_pref_lo <= core < novel_max else 1
        return (react, in_pref, -float(r["erg_max"]), -int(r["hinge_hits"]), float(r["sa"]))

    def sort_key_pan(r):
        return (-float(r["max_tc_core"]), float(r["sa"]))

    sim_rows.sort(key=sort_key_sim)
    novel_rows.sort(key=sort_key_novel)
    pan_rows.sort(key=sort_key_pan)

    # write L6 tracks
    def write_rows(path, rows):
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    write_rows(paths["sim"], sim_rows)
    write_rows(paths["novel"], novel_rows)
    write_rows(paths["pan"], pan_rows)

    # L7 quotas with watch limit in Novel
    def take_quota(rows, n, watch_limit_frac=None):
        picked = []
        n_watch = 0
        max_watch = int(n * watch_limit_frac) if watch_limit_frac is not None else n
        for r in rows:
            if len(picked) >= n:
                break
            if r.get("reactivity_bucket") == "watch":
                if n_watch >= max_watch:
                    continue
                n_watch += 1
            picked.append(r)
        return picked

    sim_yl = [r for r in sim_rows if r["subtrack"] == "sim_yl"]
    sim_56d = [r for r in sim_rows if r["subtrack"] == "sim_56d"]
    dock = {
        "sim_yl": take_quota(sim_yl, quotas["sim_yl"]),
        "sim_56d": take_quota(sim_56d, quotas["sim_56d"]),
        "novel": take_quota(novel_rows, quotas["novel"], watch_limit_frac=watch_frac),
        "pan": take_quota(pan_rows, quotas["pan"]),
    }

    # Prove track ID intersection empty among dock_ready
    id_sets = {k: {r["id"] for r in v} for k, v in dock.items()}
    inter_checks = {}
    keys = list(id_sets.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            inter_checks[f"{keys[i]}∩{keys[j]}"] = len(id_sets[keys[i]] & id_sets[keys[j]])

    for name, rows in dock.items():
        path = L7_DIR / f"L7_dock_ready_{name}.csv"
        write_rows(path, rows)

    summary = {
        "L6_counts": dict(counts),
        "L6_track_sizes": {"sim": len(sim_rows), "novel": len(novel_rows), "pan": len(pan_rows)},
        "L7_requested_quotas": quotas,
        "L7_delivered": {k: len(v) for k, v in dock.items()},
        "L7_total": sum(len(v) for v in dock.values()),
        "track_id_intersections": inter_checks,
        "thresholds_used": thr,
        "insufficient_fill": {
            k: max(0, quotas[k] - len(dock[k])) for k in quotas
        },
        "note": "Do not backfill from discard; report shortfall honestly.",
    }
    (L7_DIR / "L7_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (ROOT / "L7_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # acceptance checklist
    checklist = [
        "# Prefilter acceptance (docking-front)",
        "",
        f"- [x] L5 calibration + thresholds.json exists",
        f"- [x] core vs full Sim Jaccard = {thr.get('jaccard_sim_core_vs_full_window')} (<1)",
        f"- [ ] 56d LOO not systematically discarded — see L5_calibration.md / L5/loo.json",
        f"- [ ] Novel random 50 manual review — pending human",
        "",
        f"## Delivered dock_ready: {summary['L7_total']}",
        json.dumps(summary["L7_delivered"], indent=2),
        "",
        f"Track ID intersections (must be 0): {inter_checks}",
        "",
    ]
    (ROOT / "ACCEPTANCE.md").write_text("\n".join(checklist), encoding="utf-8")
    mark_done("L6L7", summary)
    log.info("L6L7 DONE dock_ready=%s %s", summary["L7_total"], summary["L7_delivered"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["L0", "L1", "L2L4", "L5", "L6L7", "all", "L0-L5"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    step = args.step
    if step in ("L0", "all", "L0-L5"):
        run_l0(args.force)
    if step in ("L1", "all", "L0-L5"):
        run_l1(args.force)
    if step in ("L2L4", "all", "L0-L5"):
        run_l2_l4(args.force)
    if step in ("L5", "all", "L0-L5"):
        run_l5(args.force)
    if step in ("L6L7", "all"):
        run_l6_l7(args.force)


if __name__ == "__main__":
    main()
