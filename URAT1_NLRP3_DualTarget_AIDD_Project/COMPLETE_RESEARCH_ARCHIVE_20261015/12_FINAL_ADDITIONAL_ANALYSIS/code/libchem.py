"""Read-only chemistry helpers for the additional analysis. Does not alter source data."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

BETA_LACTAM_SMARTS = "[#6]1[#6][#7][#6]1=O"
COOH_SMARTS = "[CX3](=O)[OX2H1,OX1-]"
SOFT_EXCLUDE_SUBSTR = (
    "FLOXACIN",
    "PIROMIDIC",
    "TERBOGREL",
    "KIO-100",
    "MK-5108",
    "CEFCANEL",
    "CEFAZAFLUR",
    "CEFOVECIN",
    "MARBOFLOXACIN",
)
STRUCTURAL_CONTROL_NAMES = {"GSK-3008348 FREE BASE"}
AD_THRESHOLD_DEPOSITED = 0.5672473867595819


def mol_from_smiles(smi: str):
    if not isinstance(smi, str) or not smi.strip():
        return None
    return Chem.MolFromSmiles(smi)


def canonicalize(smi: str):
    mol = mol_from_smiles(smi)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def morgan_fp(smi: str):
    mol = mol_from_smiles(smi) if not isinstance(smi, Chem.Mol) else smi
    if mol is None:
        return None
    try:
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    except Exception:
        return None


def largest_fragment(mol: Chem.Mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not frags:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def standardize_record(smi: str) -> dict:
    """Standardize without dropping failures. Descriptors use the largest fragment."""
    rec = {
        "original_smiles": smi if isinstance(smi, str) else None,
        "parse_ok": False,
        "n_components": None,
        "standardized_smiles": None,
        "largest_fragment_smiles": None,
        "murcko_scaffold": None,
        "generic_murcko_scaffold": None,
        "murcko_empty": None,
        "mw": None,
        "clogp": None,
        "tpsa": None,
        "hba": None,
        "hbd": None,
        "rotatable_bonds": None,
        "aromatic_rings": None,
        "hetero_rings": None,
        "n_carboxylate": None,
        "formal_charge_full": None,
        "formal_charge_parent": None,
        "heavy_atoms": None,
        "fp_onbits": None,
        "unspecified_stereo": None,
    }
    mol = mol_from_smiles(smi) if isinstance(smi, str) else None
    if mol is None:
        return rec
    rec["parse_ok"] = True
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        rec["parse_ok"] = False
        return rec
    rec["n_components"] = len(Chem.GetMolFrags(mol))
    rec["standardized_smiles"] = Chem.MolToSmiles(mol, isomericSmiles=True)
    rec["formal_charge_full"] = int(Chem.GetFormalCharge(mol))
    parent = largest_fragment(mol)
    try:
        Chem.SanitizeMol(parent)
    except Exception:
        parent = mol
    rec["largest_fragment_smiles"] = Chem.MolToSmiles(parent, isomericSmiles=True)
    rec["formal_charge_parent"] = int(Chem.GetFormalCharge(parent))
    try:
        core = MurckoScaffold.GetScaffoldForMol(parent)
        rec["murcko_scaffold"] = Chem.MolToSmiles(core) if core is not None else ""
        rec["murcko_empty"] = rec["murcko_scaffold"] in ("", None)
        if core is not None and core.GetNumAtoms() > 0:
            generic = MurckoScaffold.MakeScaffoldGeneric(core)
            rec["generic_murcko_scaffold"] = Chem.MolToSmiles(generic)
        else:
            rec["generic_murcko_scaffold"] = ""
    except Exception:
        rec["murcko_scaffold"] = None
        rec["generic_murcko_scaffold"] = None
        rec["murcko_empty"] = None
    rec["mw"] = float(Descriptors.MolWt(parent))
    rec["clogp"] = float(Crippen.MolLogP(parent))
    rec["tpsa"] = float(Descriptors.TPSA(parent))
    rec["hba"] = int(Lipinski.NumHAcceptors(parent))
    rec["hbd"] = int(Lipinski.NumHDonors(parent))
    rec["rotatable_bonds"] = int(Lipinski.NumRotatableBonds(parent))
    rec["aromatic_rings"] = int(Lipinski.NumAromaticRings(parent))
    rec["hetero_rings"] = int(rdMolDescriptors.CalcNumHeterocycles(parent))
    patt = Chem.MolFromSmarts(COOH_SMARTS)
    rec["n_carboxylate"] = int(len(parent.GetSubstructMatches(patt))) if patt is not None else None
    rec["heavy_atoms"] = int(parent.GetNumHeavyAtoms())
    fp = morgan_fp(parent)
    rec["fp_onbits"] = int(fp.GetNumOnBits()) if fp is not None else None
    centers = Chem.FindMolChiralCenters(parent, includeUnassigned=True, useLegacyImplementation=False)
    rec["unspecified_stereo"] = int(sum(1 for _, lab in centers if lab == "?"))
    return rec


def build_catalogs():
    families = {
        "PAINS_A": FilterCatalogParams.FilterCatalogs.PAINS_A,
        "PAINS_B": FilterCatalogParams.FilterCatalogs.PAINS_B,
        "PAINS_C": FilterCatalogParams.FilterCatalogs.PAINS_C,
        "BRENK": FilterCatalogParams.FilterCatalogs.BRENK,
        "NIH": FilterCatalogParams.FilterCatalogs.NIH,
    }
    catalogs = {}
    for label, cat in families.items():
        params = FilterCatalogParams()
        params.AddCatalog(cat)
        catalogs[label] = FilterCatalog(params)
    return catalogs


def alert_record(smi: str, catalogs, name: str = "") -> dict:
    mol = mol_from_smiles(smi) if isinstance(smi, str) else None
    beta = Chem.MolFromSmarts(BETA_LACTAM_SMARTS)
    upper = str(name).upper()
    soft = any(x in upper for x in SOFT_EXCLUDE_SUBSTR)
    control = str(name) in STRUCTURAL_CONTROL_NAMES or "GSK-3008348" in upper
    if mol is None:
        return {
            "chem_parse_ok": False,
            "pains_any": None,
            "pains_families": None,
            "brenk": None,
            "nih": None,
            "alert_descriptions": None,
            "beta_lactam": None,
            "soft_name_excluded": soft,
            "structural_control_name": control,
            "alert_clean_pains_brenk": None,
            "method_comparable": False,
        }
    pains, desc = [], []
    for label in ("PAINS_A", "PAINS_B", "PAINS_C"):
        match = catalogs[label].GetFirstMatch(mol)
        if match is not None:
            pains.append(label)
            desc.append(f"{label}:{match.GetDescription()}")
    brenk = catalogs["BRENK"].GetFirstMatch(mol)
    nih = catalogs["NIH"].GetFirstMatch(mol)
    if brenk is not None:
        desc.append(f"BRENK:{brenk.GetDescription()}")
    if nih is not None:
        desc.append(f"NIH:{nih.GetDescription()}")
    beta_hit = bool(beta and mol.HasSubstructMatch(beta))
    pains_any = len(pains) > 0
    brenk_flag = brenk is not None
    return {
        "chem_parse_ok": True,
        "pains_any": pains_any,
        "pains_families": ";".join(pains),
        "brenk": brenk_flag,
        "nih": nih is not None,
        "alert_descriptions": " | ".join(desc),
        "beta_lactam": beta_hit,
        "soft_name_excluded": soft,
        "structural_control_name": control,
        "alert_clean_pains_brenk": (not pains_any) and (not brenk_flag),
        "method_comparable": True,
    }


def nearest_tanimoto(query_fps, ref_fps, exclude_self: bool = False) -> np.ndarray:
    """Same rule as scripts/nlrp3_v2_stage2.py:nearest_tanimoto."""
    ref = [f for f in ref_fps if f is not None]
    out = np.full(len(query_fps), np.nan)
    for i, q in enumerate(query_fps):
        if q is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(q, ref)
        if not sims:
            continue
        if exclude_self:
            sims = sorted(sims, reverse=True)
            out[i] = sims[1] if len(sims) > 1 else sims[0]
        else:
            out[i] = max(sims)
    return out


def quantiles(values, prefix: str = "") -> dict:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        keys = ["n", "min", "p05", "p25", "p50", "p75", "p95", "max", "mean"]
        return {f"{prefix}{k}": None for k in keys}
    qs = np.quantile(arr, [0.05, 0.25, 0.50, 0.75, 0.95])
    return {
        f"{prefix}n": int(len(arr)),
        f"{prefix}min": float(arr.min()),
        f"{prefix}p05": float(qs[0]),
        f"{prefix}p25": float(qs[1]),
        f"{prefix}p50": float(qs[2]),
        f"{prefix}p75": float(qs[3]),
        f"{prefix}p95": float(qs[4]),
        f"{prefix}max": float(arr.max()),
        f"{prefix}mean": float(arr.mean()),
    }


def dist(a, b) -> float:
    return math.sqrt(sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)))
