#!/usr/bin/env python3
"""C5 pre-submission audit: leakage/novelty, gate intervals, Arg sensitivity,
limited receptor-geometry check. Zero new docking.

Verify-then-compute. Does not retune gates or reshuffle the frozen 12.
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import inchi, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import fisher_exact

MFP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
FP_CACHE: dict[str, object] = {}

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from parse_c1_sdf_readouts import _fprop, carboxylate_oxygens, load_poses, min_acid_arg_dist  # noqa: E402
from run_c5_w2_urat1_ifp_gate import kabsch_R_t, transform_xyz  # noqa: E402

OUT = PROJECT_ROOT / "data/campaigns/c5/05_presubmission_audit"
C1 = PROJECT_ROOT / "data/campaigns/c1"
C5 = PROJECT_ROOT / "data/campaigns/c5"
STRUCT = PROJECT_ROOT / "data/structures/pdb"

CRYSTAL_ARG_A = 6.7026966215098716
LOCKED_ARG_A = 7.7027
ARG_DELTAS = (0.5, 1.0, 1.5)
PHE_CAGE_9DKB = (241, 360, 364, 365, 449)
NLRP3_KEYS = (227, 228, 351, 408, 443, 575, 578)
CONTACT_A = 4.5
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# Vendor-reported SMILES for published dual-target comparators (not ChEMBL train).
# HNW005: DC Chemicals page citing Sun et al., Eur. J. Med. Chem. 2025, 117644.
# Compound 32: TargetMol NLRP3/URAT1-IN-1 = Nat. Commun. 2025, 16:7430 compound 32.
NAMED_DUAL = [
    {
        "name": "HNW005",
        "smiles": "O=C(Nc1ccccc1C(=O)O)c1cc(-c2ccc(C(F)(F)F)cc2)[nH]n1",
        "source": "vendor_DCChemicals_citing_EJMC_2025_117644",
        "note": "NLRP3 KD 204.6 nM / IC50 1.7 uM; URAT1 IC50 6.4 uM",
    },
    {
        "name": "NatCommun_compound_32",
        "smiles": "O=C(NS(=O)(=O)c1cccs1)c1cc2ccccc2n1Cc1ccc(Br)c2ccccc12",
        "source": "vendor_TargetMol_NLRP3_URAT1_IN_1_citing_NatCommun_2025_7430",
        "note": "URAT1 IC50 ~3.8 uM; NLRP3 SPR KD 27.8 uM",
    },
]


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    margin = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / den
    return max(0.0, center - margin), min(1.0, center + margin)


def haldane_or(tp: int, fn: int, fp: int, tn: int) -> float:
    return ((tp + 0.5) * (tn + 0.5)) / ((fn + 0.5) * (fp + 0.5))


def bootstrap_lr_or(y_true, y_pred, n_boot: int = 2000, seed: int = 42) -> dict:
    rng = random.Random(seed)
    n = len(y_true)
    ors, lrs = [], []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        yt = [y_true[i] for i in idx]
        yp = [y_pred[i] for i in idx]
        a = sum(1 for t, p in zip(yt, yp) if t == 1 and p == 1)
        b = sum(1 for t, p in zip(yt, yp) if t == 1 and p == 0)
        c = sum(1 for t, p in zip(yt, yp) if t == 0 and p == 1)
        d = sum(1 for t, p in zip(yt, yp) if t == 0 and p == 0)
        ors.append(haldane_or(a, b, c, d))
        sens = a / (a + b) if (a + b) else 0.0
        fpr = c / (c + d) if (c + d) else 0.0
        lrs.append(sens / fpr if fpr > 0 else (float("inf") if sens > 0 else 0.0))
    finite_or = [x for x in ors if math.isfinite(x)]
    finite_lr = [x for x in lrs if math.isfinite(x)]
    def pct(xs):
        if not xs:
            return float("nan"), float("nan")
        lo, hi = np.percentile(xs, [2.5, 97.5])
        return float(lo), float(hi)
    or_lo, or_hi = pct(finite_or)
    lr_lo, lr_hi = pct(finite_lr)
    return {
        "or_boot_lo": or_lo,
        "or_boot_hi": or_hi,
        "lr_boot_lo": lr_lo,
        "lr_boot_hi": lr_hi,
        "n_boot_finite_or": len(finite_or),
        "n_boot_finite_lr": len(finite_lr),
    }


def metrics_with_intervals(labels: list[int], passes: list[bool]) -> dict:
    y = np.array(labels, dtype=int)
    p = np.array(passes, dtype=bool)
    tp = int(((y == 1) & p).sum())
    fn = int(((y == 1) & ~p).sum())
    fp = int(((y == 0) & p).sum())
    tn = int(((y == 0) & ~p).sum())
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) else float("nan")
    fpr = 1.0 - spec if math.isfinite(spec) else float("nan")
    lr_plus = sens / fpr if (math.isfinite(sens) and fpr and fpr > 0) else (
        float("inf") if sens == 1.0 and fpr == 0.0 else float("nan")
    )
    table = [[tp, fn], [fp, tn]]
    oddsr, pval = fisher_exact(table)
    boot = bootstrap_lr_or(labels, passes)
    s_lo, s_hi = wilson_ci(tp, tp + fn)
    sp_lo, sp_hi = wilson_ci(tn, tn + fp)
    return {
        "n": int(len(labels)),
        "n_active": int((y == 1).sum()),
        "n_decoy": int((y == 0).sum()),
        "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "sensitivity": sens,
        "sensitivity_wilson95_lo": s_lo,
        "sensitivity_wilson95_hi": s_hi,
        "specificity": spec,
        "specificity_wilson95_lo": sp_lo,
        "specificity_wilson95_hi": sp_hi,
        "ppv": ppv,
        "lr_plus": float(lr_plus) if math.isfinite(lr_plus) else None,
        "odds_ratio_mle": float(oddsr) if math.isfinite(oddsr) else None,
        "odds_ratio_haldane": haldane_or(tp, fn, fp, tn),
        "fisher_exact_p": float(pval),
        **boot,
    }


def canon(smi: str) -> str | None:
    m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def inchikey(smi: str) -> str | None:
    m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if m is None:
        return None
    try:
        return inchi.MolToInchiKey(m)
    except Exception:
        return None


def scaffold(smi: str) -> str | None:
    if not smi or (isinstance(smi, float) and math.isnan(smi)):
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(str(smi))
    except Exception:
        return None


def fp_of(smi: str):
    if smi in FP_CACHE:
        return FP_CACHE[smi]
    m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if m is None:
        FP_CACHE[smi] = None
        return None
    fp = MFP.GetFingerprint(m)
    FP_CACHE[smi] = fp
    return fp


def nearest(query_smi: str, library: list[tuple[str, str, str]]) -> dict:
    q = fp_of(query_smi)
    if q is None or not library:
        return {"max_tc": None, "nn_id": None, "nn_name": None, "nn_smiles": None}
    best = -1.0
    hit = (None, None, None)
    for lid, name, smi in library:
        f = fp_of(smi)
        if f is None:
            continue
        tc = DataStructs.TanimotoSimilarity(q, f)
        if tc > best:
            best = tc
            hit = (lid, name, smi)
    return {"max_tc": round(float(best), 4), "nn_id": hit[0], "nn_name": hit[1], "nn_smiles": hit[2]}


def pubmed_count(term: str, retries: int = 3) -> dict:
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmode": "json", "retmax": 5})
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{q}"
    last_err = None
    for i in range(retries):
        try:
            time.sleep(0.35)
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.loads(r.read().decode())
            ids = data.get("esearchresult", {}).get("idlist", [])
            count = int(data.get("esearchresult", {}).get("count", 0))
            return {"term": term, "count": count, "pmids": ids, "ok": True}
        except Exception as e:
            last_err = str(e)
            time.sleep(1.0 * (i + 1))
    return {"term": term, "count": None, "pmids": [], "ok": False, "error": last_err}


def chain_seq_ca(path: Path, chain: str = "A"):
    st = gemmi.read_structure(str(path))
    model = st[0]
    ch = model[chain] if chain in [c.name for c in model] else model[0]
    seq, ca, atoms = [], {}, {}
    for res in ch:
        if res.name not in AA3:
            continue
        try:
            n = int(res.seqid.num)
        except Exception:
            continue
        seq.append((n, res.name, AA3[res.name]))
        xyz = {}
        for atom in res:
            name = atom.name.strip()
            xyz[name] = np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
        atoms[n] = xyz
        if "CA" in xyz:
            ca[n] = xyz["CA"]
    return seq, ca, atoms


def needleman(a: str, b: str):
    n, m = len(a), len(b)
    f = np.zeros((n + 1, m + 1))
    ptr = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        f[i, 0] = -i
        ptr[i, 0] = 1
    for j in range(1, m + 1):
        f[0, j] = -j
        ptr[0, j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = f[i - 1, j - 1] + (1 if a[i - 1] == b[j - 1] else -1)
            up = f[i - 1, j] - 1
            left = f[i, j - 1] - 1
            best = max(diag, up, left)
            f[i, j] = best
            ptr[i, j] = 0 if best == diag else (1 if best == up else 2)
    i, j = n, m
    pairs = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and ptr[i, j] == 0:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or ptr[i, j] == 1):
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    return pairs[::-1]


def align_map(src: Path, dst: Path):
    s1, c1, a1 = chain_seq_ca(src)
    s2, c2, a2 = chain_seq_ca(dst)
    pairs = needleman("".join(x[2] for x in s1), "".join(x[2] for x in s2))
    mapped = {}
    common = []
    n_id = 0
    n_paired = 0
    for i, j in pairs:
        if i is None or j is None:
            continue
        n_paired += 1
        if s1[i][2] == s2[j][2]:
            n_id += 1
        n1, n2 = s1[i][0], s2[j][0]
        mapped[n1] = (n2, s2[j][1])
        if n1 in c1 and n2 in c2:
            common.append((n1, n2))
    P = np.vstack([c1[n] for n, _ in common])
    Q = np.vstack([c2[m] for _, m in common])
    R, t = kabsch_R_t(P, Q)
    rms = float(np.sqrt(((transform_xyz(P, R, t) - Q) ** 2).sum(1).mean()))
    return {
        "R": R, "t": t, "mapped": mapped, "ca_rms": rms, "n_ca": len(common),
        "identity": n_id / n_paired if n_paired else None,
        "dst_atoms": a2, "src_atoms": a1,
    }


def cnn_top1(sdf: Path):
    poses = load_poses(sdf)
    if not poses:
        return None
    return max(poses, key=lambda m: _fprop(m, "CNNscore") or -1.0)


def mol_xyz(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    pts = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(atom.GetIdx())
        pts.append([p.x, p.y, p.z])
    return np.array(pts, dtype=float)


def arg_atoms_from_map(aln: dict, src_resi: int = 477) -> dict | None:
    hit = aln["mapped"].get(src_resi)
    if not hit:
        return None
    dst_n, dst_name = hit
    if dst_name != "ARG":
        return None
    xyz = aln["dst_atoms"].get(dst_n, {})
    out = {k: xyz[k].tolist() for k in ("NE", "NH1", "NH2") if k in xyz}
    return out or None


def phe_contacts(lig_xyz: np.ndarray, aln: dict, src_resis=PHE_CAGE_9DKB) -> dict:
    n_hit = 0
    mins = {}
    for r in src_resis:
        hit = aln["mapped"].get(r)
        if not hit:
            mins[f"F{r}"] = None
            continue
        dst_n, dst_name = hit
        ring = []
        for name, xyz in aln["dst_atoms"].get(dst_n, {}).items():
            if name in {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}:
                ring.append(xyz)
        if not ring:
            mins[f"F{r}"] = None
            continue
        ring = np.vstack(ring)
        d = float(np.linalg.norm(lig_xyz[:, None, :] - ring[None, :, :], axis=2).min())
        mins[f"d_F{r}_mapped{dst_n}_{dst_name}"] = round(d, 3)
        if dst_name == "PHE" and d <= CONTACT_A:
            n_hit += 1
    mins["n_phe_cage_le_4.5"] = n_hit
    mins["phe_cage_ok"] = n_hit >= 2
    return mins


def verify_inputs() -> dict:
    required = {
        "primary": C5 / "04_shortlist_frozen/primary_candidates.csv",
        "audited": C1 / "08_nomination/acid_a2_eligible_audited.csv",
        "a1_keep": C1 / "07_clinical_dock/acid_dual_a1_frozen/acid_dual_keep_seed42.csv",
        "a1_bench": C1 / "05_metrics/acid_gate_retrospective_benchmark/acid_gate_benchmark_per_mol.csv",
        "w2_sum": C5 / "02_urat1_ifp/w2_ifp_gate_summary.json",
        "w4_sum": C5 / "02_nlrp3_panel/w4_structural_gate_summary.json",
        "w4_metrics": C5 / "02_nlrp3_panel/w4_panel_metrics_all_seeds.csv",
        "nlrp3_train": PROJECT_ROOT / "data/processed/nlrp3_records.csv",
        "urat1": PROJECT_ROOT / "data/processed/urat1_curated.csv",
        "library": PROJECT_ROOT / "data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv",
        "w4_pos": C1 / "05_metrics/nlrp3_structural_panel/panel_ligands.csv",
        "9dkb": STRUCT / "9DKB.cif",
        "9dka": STRUCT / "9DKA.cif",
        "9dkc": STRUCT / "9DKC.cif",
        "9b1i": STRUCT / "9B1I.cif",
        "7alv": STRUCT / "7ALV.pdb",
        "8etr": STRUCT / "8ETR.cif",
    }
    missing = [k for k, p in required.items() if not p.exists()]
    return {"paths": {k: str(p) for k, p in required.items()}, "missing": missing}


def load_primary() -> pd.DataFrame:
    return pd.read_csv(C5 / "04_shortlist_frozen/primary_candidates.csv")


def build_ref_sets():
    nl = pd.read_csv(PROJECT_ROOT / "data/processed/nlrp3_records.csv")
    ur = pd.read_csv(PROJECT_ROOT / "data/processed/urat1_curated.csv")
    lib = pd.read_csv(
        PROJECT_ROOT / "data/repurposing/screening/nlrp3_ml_scores_clinical_all.csv",
        low_memory=False,
    )
    panel = pd.read_csv(C1 / "05_metrics/nlrp3_structural_panel/panel_ligands.csv")
    w4_ids = {
        "NP3-146", "MCC950", "CHEMBL4204644", "CHEMBL5219789", "CHEMBL4212407",
        "CHEMBL6143743", "CHEMBL4209503", "CHEMBL4216836", "CHEMBL6171925",
    }
    w4 = panel[panel["ligand_id"].isin(w4_ids)].copy()

    def pack(df, id_col, name_col, smi_col):
        rows = []
        seen = set()
        for _, r in df.iterrows():
            smi = canon(r.get(smi_col))
            if not smi or smi in seen:
                continue
            seen.add(smi)
            rows.append((str(r.get(id_col)), str(r.get(name_col) or r.get(id_col)), smi))
        return rows

    train_all = pack(nl, "Molecule ChEMBL ID", "Molecule Name", "canonical_smiles")
    train_act = pack(nl[nl["active"] == 1], "Molecule ChEMBL ID", "Molecule Name", "canonical_smiles")
    urat1 = pack(ur, "molecule_chembl_id", "molecule_name", "canonical_smiles")
    direct = pack(w4, "ligand_id", "name", "smiles")
    for item in NAMED_DUAL:
        smi = canon(item["smiles"])
        if smi:
            direct.append((item["name"], item["name"], smi))
            urat1.append((item["name"], item["name"], smi))
    return {
        "nlrp3_records": nl,
        "urat1_curated": ur,
        "library": lib,
        "w4": w4,
        "train_all": train_all,
        "train_act": train_act,
        "urat1": urat1,
        "direct_nlrp3": direct,
    }


def leakage_block(primary: pd.DataFrame, refs: dict) -> tuple[pd.DataFrame, dict]:
    lib = refs["library"]
    nl = refs["nlrp3_records"]
    train_ids = set(nl["Molecule ChEMBL ID"].astype(str))
    lib_ids = set(lib["chembl_id"].astype(str))
    train_in_lib = sorted(train_ids & lib_ids)
    train_in_lib_names = (
        lib.loc[lib["chembl_id"].astype(str).isin(train_in_lib), ["chembl_id", "name"]]
        .drop_duplicates()
        .to_dict("records")
    )
    train_keys = {inchikey(s) for s in nl["canonical_smiles"] if canon(s)}
    lib_keys = {k for k in lib["inchi_key"].astype(str) if k and k != "nan"}
    rows = []
    for _, r in primary.iterrows():
        smi = canon(r["canonical_smiles"])
        key = inchikey(smi) if smi else None
        cid = str(r["chembl_id"])
        rec = {
            "ligand_id": r["ligand_id"],
            "name": r["name"],
            "chembl_id": cid,
            "inchi_key": key,
            "canonical_smiles": smi,
            "scaffold": scaffold(smi),
            "in_nlrp3_train_chembl_id": cid in train_ids,
            "in_nlrp3_train_inchikey": key in train_keys if key else False,
            "in_clinical_library_8319": cid in lib_ids,
            "in_urat1_curated_chembl_id": cid in set(refs["urat1_curated"]["molecule_chembl_id"].astype(str)),
        }
        for label, pool in (
            ("nlrp3_train_all", refs["train_all"]),
            ("nlrp3_train_active", refs["train_act"]),
            ("nlrp3_direct_or_published_dual", refs["direct_nlrp3"]),
            ("urat1_known", refs["urat1"]),
        ):
            nn = nearest(smi, pool)
            rec[f"max_tc_{label}"] = nn["max_tc"]
            rec[f"nn_{label}_id"] = nn["nn_id"]
            rec[f"nn_{label}_name"] = nn["nn_name"]
        rows.append(rec)

    lit_rows = []
    for _, r in primary.iterrows():
        name = str(r["name"]).replace(" FREE ACID", "")
        cid = str(r["chembl_id"])
        for term in (
            f"{name} AND (URAT1 OR SLC22A12)",
            f"{name} AND NLRP3",
            f"{cid} AND (URAT1 OR NLRP3 OR inflammasome)",
        ):
            hit = pubmed_count(term)
            hit["ligand_id"] = r["ligand_id"]
            hit["name"] = r["name"]
            lit_rows.append(hit)

    extra_terms = [
        "PF-03882845 AND (URAT1 OR SLC22A12 OR NLRP3 OR gout)",
        "PF-03882845 AND mineralocorticoid",
        "HNW005 AND (URAT1 OR NLRP3)",
    ]
    extra_lit = [pubmed_count(t) for t in extra_terms]

    summary = {
        "n_primary": int(len(primary)),
        "n_nlrp3_train_mols": int(nl["Molecule ChEMBL ID"].nunique()),
        "n_nlrp3_assays": int(nl["Assay ChEMBL ID"].nunique()),
        "n_library": int(len(lib)),
        "n_train_ids_in_8319": int(len(train_in_lib)),
        "train_ids_in_8319": train_in_lib,
        "train_names_in_8319": train_in_lib_names,
        "n_train_inchikey_in_8319": int(len({k for k in train_keys if k} & lib_keys)),
        "any_primary_exact_train": bool(any(x["in_nlrp3_train_chembl_id"] or x["in_nlrp3_train_inchikey"] for x in rows)),
        "any_primary_in_urat1_curated": bool(any(x["in_urat1_curated_chembl_id"] for x in rows)),
        "literature_extra": extra_lit,
        "named_dual_smiles_source": NAMED_DUAL,
    }
    return pd.DataFrame(rows), pd.DataFrame(lit_rows), summary


def gate_intervals() -> tuple[pd.DataFrame, dict]:
    bench = pd.read_csv(C1 / "05_metrics/acid_gate_retrospective_benchmark/acid_gate_benchmark_per_mol.csv")
    a1 = bench[bench["pose_selection_rule"] == "a1"].copy()
    a2 = bench[bench["pose_selection_rule"] == "a2"].copy()
    w2 = json.loads((C5 / "02_urat1_ifp/w2_ifp_gate_summary.json").read_text())
    w4m = pd.read_csv(C5 / "02_nlrp3_panel/w4_panel_metrics_all_seeds.csv")
    w4s = json.loads((C5 / "02_nlrp3_panel/w4_structural_gate_summary.json").read_text())

    rows = []
    for name, df, pass_col in (
        ("A1_arg_locked_7.7027", a1, "pass_arg_A1"),
        ("A2_geometry_then_cnn", a2, "keep_urat1_acid"),
    ):
        labels = [int(x) for x in df["label"]]
        passes = [bool(x) for x in df[pass_col]]
        m = metrics_with_intervals(labels, passes)
        m["gate"] = name
        rows.append(m)

    # Recompute IFP CNN Top-1 from W2 stored 2x2
    ifp = w2["ifp_cnn_top1"]
    labels = [1] * ifp["n_active"] + [0] * ifp["n_decoy"]
    passes = [True] * ifp["tp"] + [False] * ifp["fn"] + [True] * ifp["fp"] + [False] * ifp["tn"]
    m = metrics_with_intervals(labels, passes)
    m["gate"] = "W2_IFP_cnn_top1"
    m["note"] = "reconstructed from stored 2x2; bootstrap is on reconstructed labels"
    rows.append(m)

    w4_seed42 = w4m[(w4m["seed"] == 42) & (w4m["role"].isin(["positive", "decoy"]))]
    labels = [1 if r == "positive" else 0 for r in w4_seed42["role"]]
    for col, gname in (("keep_nlrp3_pose", "W4_loose_seed42"), ("keep_nlrp3_structural", "W4_structural_seed42")):
        passes = [bool(x) for x in w4_seed42[col]]
        m = metrics_with_intervals(labels, passes)
        m["gate"] = gname
        rows.append(m)

    pos = w4m[(w4m["seed"] == 42) & (w4m["role"] == "positive")]
    scaf_rows = []
    panel = pd.read_csv(C1 / "05_metrics/nlrp3_structural_panel/panel_ligands.csv")
    smi_map = dict(zip(panel["ligand_id"], panel["smiles"]))
    for lid in pos["ligand_id"]:
        smi = smi_map.get(lid)
        scaf_rows.append({
            "ligand_id": lid,
            "smiles": smi,
            "scaffold": scaffold(smi) if smi else None,
            "keep_structural": bool(pos.loc[pos["ligand_id"] == lid, "keep_nlrp3_structural"].iloc[0]),
        })
    scaf = pd.DataFrame(scaf_rows)
    summary = {
        "w4_positive_n": int(len(scaf)),
        "w4_unique_murcko": int(scaf["scaffold"].nunique(dropna=True)),
        "w4_scaffolds": scaf["scaffold"].value_counts().to_dict(),
        "w4_chemotype_note": (
            "W4 positives are sulfonylurea / NP3-146–MCC950 family. "
            "9/9 sensitivity is within this chemotype, not a general NACHT-pocket classifier."
        ),
        "w2_or_ci_already_in_repo": {
            "a1": w2["a1"],
            "a2": w2["a2"],
            "ifp_cnn_top1": w2["ifp_cnn_top1"],
        },
        "w4_primary_seed42_structural": w4s["primary"]["structural"],
    }
    return pd.DataFrame(rows), scaf, summary


def threshold_sensitivity(primary: pd.DataFrame, refs: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    bench = pd.read_csv(C1 / "05_metrics/acid_gate_retrospective_benchmark/acid_gate_benchmark_per_mol.csv")
    a1 = bench[bench["pose_selection_rule"] == "a1"].copy()
    audited = pd.read_csv(C1 / "08_nomination/acid_a2_eligible_audited.csv")
    a1_keep = pd.read_csv(C1 / "07_clinical_dock/acid_dual_a1_frozen/acid_dual_keep_seed42.csv")
    # eligible = A2 structural >=2/3 and chemistry-audited (same as freeze)
    cnt = {}
    for s in (42, 43, 44):
        df = pd.read_csv(C1 / f"07_clinical_dock/acid_dual_a2/acid_dual_keep_structural_seed{s}.csv")
        for lid in df.loc[df["keep_dual_acid_structural"] == True, "ligand_id"]:  # noqa: E712
            cnt[lid] = cnt.get(lid, 0) + 1
    a2_ge2 = {lid for lid, c in cnt.items() if c >= 2}
    eligible = set(audited["ligand_id"]) & a2_ge2
    name_of = dict(zip(audited["ligand_id"], audited["name"]))
    frozen12 = set(primary["ligand_id"])
    arg_of = dict(zip(a1_keep["ligand_id"], a1_keep["acid_arg477_min_A"]))

    bench_rows = []
    name_rows = []
    for delta in ARG_DELTAS:
        thresh = CRYSTAL_ARG_A + delta
        labels = [int(x) for x in a1["label"]]
        passes = [
            bool(pd.notna(d) and float(d) <= thresh)
            for d in a1["acid_arg477_min_A"]
        ]
        m = metrics_with_intervals(labels, passes)
        m["delta_A"] = delta
        m["arg_thresh_A"] = thresh
        m["is_locked_plus_1A"] = abs(thresh - LOCKED_ARG_A) < 1e-4
        bench_rows.append(m)

        passing = {
            lid for lid in eligible
            if pd.notna(arg_of.get(lid)) and float(arg_of[lid]) <= thresh
        }
        name_rows.append({
            "delta_A": delta,
            "arg_thresh_A": thresh,
            "n_eligible_pass_A1": len(passing),
            "n_frozen12_retained": len(passing & frozen12),
            "n_frozen12_lost": len(frozen12 - passing),
            "lost_names": sorted(name_of.get(x, x) for x in frozen12 - passing),
            "gained_vs_locked12": sorted(name_of.get(x, x) for x in passing - frozen12),
            "pass_names": sorted(name_of.get(x, x) for x in passing),
        })

    # Per-molecule Arg vs each threshold
    per = []
    for lid in sorted(frozen12):
        d = arg_of.get(lid)
        rec = {
            "ligand_id": lid,
            "name": name_of.get(lid, lid),
            "acid_arg477_min_A_seed42_A1": d,
        }
        for delta in ARG_DELTAS:
            rec[f"pass_plus_{delta}"] = bool(pd.notna(d) and float(d) <= CRYSTAL_ARG_A + delta)
        per.append(rec)

    summary = {
        "crystal_arg_A": CRYSTAL_ARG_A,
        "locked_arg_A": LOCKED_ARG_A,
        "n_eligible_chemistry_and_A2ge2": len(eligible),
        "note": (
            "Name table holds A2-structural>=2/3 ∩ chemistry-audited fixed, "
            "and only varies the A1 Arg cutoff. No retuning."
        ),
    }
    return pd.DataFrame(bench_rows), pd.DataFrame(name_rows), pd.DataFrame(per), summary


def urat1_sdf(lid: str, seed: int) -> Path:
    if seed == 42:
        return C1 / f"07_clinical_dock/acid_dual/urat1_9dkb/seed42/{lid}_out.sdf"
    return C1 / f"07_clinical_dock/acid_dual_a2/urat1_9dkb/seed{seed}/{lid}_out.sdf"


def nlrp3_sdf(lid: str, seed: int) -> Path:
    if seed == 42:
        return C1 / f"07_clinical_dock/acid_dual/nlrp3_7alv/seed42/{lid}_out.sdf"
    return C1 / f"07_clinical_dock/acid_dual_a2/nlrp3_7alv/seed{seed}/{lid}_out.sdf"


def receptor_sensitivity(primary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    src_u = STRUCT / "9DKB.cif"
    targets = {
        "9DKB_native": STRUCT / "9DKB.cif",
        "9DKA_benzbromarone": STRUCT / "9DKA.cif",
        "9DKC_TD3": STRUCT / "9DKC.cif",
        "9B1I_verinurad": STRUCT / "9B1I.cif",
    }
    aln_u = {k: align_map(src_u, p) for k, p in targets.items()}
    src_n = STRUCT / "7ALV.pdb"
    aln_n = {
        "7ALV_native": align_map(src_n, STRUCT / "7ALV.pdb"),
        "8ETR_GDC2394": align_map(src_n, STRUCT / "8ETR.cif"),
    }

    u_rows = []
    for _, r in primary.iterrows():
        lid = r["ligand_id"]
        sdf = urat1_sdf(lid, 42)
        pose = cnn_top1(sdf) if sdf.exists() else None
        rec = {
            "ligand_id": lid,
            "name": r["name"],
            "seed": 42,
            "sdf_exists": sdf.exists(),
            "method": "pose_transfer_CA_Kabsch_not_redock",
        }
        if pose is None:
            rec["error"] = "missing_or_empty_sdf"
            u_rows.append(rec)
            continue
        xyz0 = mol_xyz(pose)
        for tname, aln in aln_u.items():
            xyz = transform_xyz(xyz0, aln["R"], aln["t"])
            arg = arg_atoms_from_map(aln, 477)
            # rebuild a temporary mol? min_acid_arg_dist needs mol coords.
            # Measure using transferred carboxylate oxygens.
            oxy = np.array(carboxylate_oxygens(pose), dtype=float) if carboxylate_oxygens(pose) else None
            if oxy is None or arg is None:
                rec[f"{tname}_arg_A"] = None
                rec[f"{tname}_arg_pass_7.7027"] = False
            else:
                oxy_t = transform_xyz(oxy, aln["R"], aln["t"])
                nit = np.array([arg[k] for k in arg], dtype=float)
                d = float(np.linalg.norm(oxy_t[:, None, :] - nit[None, :, :], axis=2).min())
                rec[f"{tname}_arg_A"] = round(d, 3)
                rec[f"{tname}_arg_pass_7.7027"] = d <= LOCKED_ARG_A
            phe = phe_contacts(xyz, aln)
            rec[f"{tname}_n_phe_cage"] = phe["n_phe_cage_le_4.5"]
            rec[f"{tname}_phe_ok"] = phe["phe_cage_ok"]
        rec["n_receptors_arg_pass"] = sum(
            bool(rec.get(f"{k}_arg_pass_7.7027")) for k in targets
        )
        rec["n_receptors_phe_ok"] = sum(bool(rec.get(f"{k}_phe_ok")) for k in targets)
        u_rows.append(rec)

    n_rows = []
    key_json = json.loads((C1 / "01_ligand_prep/selfdock_refs/nlrp3_key_residues.json").read_text())
    key_names = list(key_json["residues"].keys())
    for _, r in primary.iterrows():
        lid = r["ligand_id"]
        sdf = nlrp3_sdf(lid, 42)
        pose = cnn_top1(sdf) if sdf.exists() else None
        rec = {
            "ligand_id": lid,
            "name": r["name"],
            "seed": 42,
            "sdf_exists": sdf.exists(),
            "method": "pose_transfer_CA_Kabsch_not_redock",
        }
        if pose is None:
            rec["error"] = "missing_or_empty_sdf"
            n_rows.append(rec)
            continue
        xyz0 = mol_xyz(pose)
        for tname, aln in aln_n.items():
            xyz = transform_xyz(xyz0, aln["R"], aln["t"])
            hits = 0
            for src_n in NLRP3_KEYS:
                mapped = aln["mapped"].get(src_n)
                if not mapped:
                    continue
                dst_n, dst_name = mapped
                heavy = np.vstack(list(aln["dst_atoms"][dst_n].values())) if dst_n in aln["dst_atoms"] else None
                if heavy is None:
                    continue
                d = float(np.linalg.norm(xyz[:, None, :] - heavy[None, :, :], axis=2).min())
                rec[f"{tname}_d_{dst_name}{dst_n}"] = round(d, 3)
                if d <= CONTACT_A:
                    hits += 1
            rec[f"{tname}_n_key_le_4.5"] = hits
            rec[f"{tname}_key_ge5"] = hits >= 5
            ds = [rec[k] for k in rec if k.startswith(f"{tname}_d_") and rec[k] is not None]
            rec[f"{tname}_min_key_d"] = min(ds) if ds else None
            rec[f"{tname}_n_clash_lt_2.2"] = sum(1 for d in ds if d < 2.2)
            rec[f"{tname}_key_ge5_no_clash"] = bool(hits >= 5 and all(d >= 2.2 for d in ds))
        rec["both_nlrp3_receptors_key_ge5"] = bool(
            rec.get("7ALV_native_key_ge5") and rec.get("8ETR_GDC2394_key_ge5")
        )
        rec["both_nlrp3_key_ge5_no_clash"] = bool(
            rec.get("7ALV_native_key_ge5_no_clash") and rec.get("8ETR_GDC2394_key_ge5_no_clash")
        )
        n_rows.append(rec)

    summary = {
        "method": "pose_transfer_not_redock",
        "gnina_available": False,
        "urat1_align": {
            k: {"ca_rms": v["ca_rms"], "n_ca": v["n_ca"], "identity": v["identity"],
                "arg477_maps_to": v["mapped"].get(477)}
            for k, v in aln_u.items()
        },
        "nlrp3_align": {
            k: {"ca_rms": v["ca_rms"], "n_ca": v["n_ca"], "identity": v["identity"]}
            for k, v in aln_n.items()
        },
        "interpretation": (
            "This asks whether the frozen 9DKB/7ALV CNNscore Top-1 pose remains "
            "geometrically compatible after rigid CA superposition onto another "
            "experimental conformation. It does not ask whether GNINA would find "
            "a new pose on that receptor."
        ),
        "nlrp3_key_residue_file": key_names,
    }
    return pd.DataFrame(u_rows), pd.DataFrame(n_rows), summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inv = verify_inputs()
    (OUT / "inventory.json").write_text(json.dumps(inv, indent=2))
    if inv["missing"]:
        raise SystemExit(f"missing inputs: {inv['missing']}")

    primary = load_primary()
    if len(primary) != 12:
        raise SystemExit(f"expected 12 primary, got {len(primary)}")

    refs = build_ref_sets()
    leak, lit, leak_sum = leakage_block(primary, refs)
    leak.to_csv(OUT / "novelty_leakage_primary.csv", index=False)
    lit.to_csv(OUT / "literature_pubmed.csv", index=False)
    (OUT / "novelty_leakage_summary.json").write_text(json.dumps(leak_sum, indent=2, default=str))

    gates, w4scaf, gate_sum = gate_intervals()
    gates.to_csv(OUT / "gate_interval_metrics.csv", index=False)
    w4scaf.to_csv(OUT / "w4_positive_scaffolds.csv", index=False)
    (OUT / "gate_interval_summary.json").write_text(json.dumps(gate_sum, indent=2, default=str))

    bench_s, name_s, per_s, thr_sum = threshold_sensitivity(primary, refs)
    bench_s.to_csv(OUT / "arg_threshold_sensitivity_benchmark.csv", index=False)
    name_s.to_csv(OUT / "arg_threshold_sensitivity_names.csv", index=False)
    per_s.to_csv(OUT / "arg_threshold_sensitivity_primary.csv", index=False)
    (OUT / "arg_threshold_sensitivity_summary.json").write_text(json.dumps(thr_sum, indent=2, default=str))

    u_df, n_df, rec_sum = receptor_sensitivity(primary)
    u_df.to_csv(OUT / "receptor_sensitivity_urat1.csv", index=False)
    n_df.to_csv(OUT / "receptor_sensitivity_nlrp3.csv", index=False)
    (OUT / "receptor_sensitivity_summary.json").write_text(json.dumps(rec_sum, indent=2, default=str))

    print("wrote", OUT)
    print("leakage exact-train", leak_sum["any_primary_exact_train"],
          "urat1-curated", leak_sum["any_primary_in_urat1_curated"])
    print("train-in-8319", leak_sum["n_train_ids_in_8319"])
    print(gates[["gate", "tp", "fn", "fp", "tn", "sensitivity", "lr_plus", "odds_ratio_haldane"]].to_string(index=False))
    print(name_s.to_string(index=False))
    print(u_df[["ligand_id", "name", "n_receptors_arg_pass", "n_receptors_phe_ok"]].to_string(index=False))
    cols = [c for c in (
        "ligand_id", "name", "both_nlrp3_receptors_key_ge5", "both_nlrp3_key_ge5_no_clash",
        "8ETR_GDC2394_min_key_d", "8ETR_GDC2394_n_clash_lt_2.2",
    ) if c in n_df.columns]
    print(n_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
