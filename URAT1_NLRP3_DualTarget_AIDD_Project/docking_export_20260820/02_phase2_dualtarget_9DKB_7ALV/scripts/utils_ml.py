"""Minimal SMILES helpers used by merge_docking_pareto (omitted from server upload pack)."""
from __future__ import annotations

from typing import Optional

from rdkit import Chem


def canonicalize(smiles: str) -> Optional[str]:
    """Return RDKit canonical SMILES, or None if parsing fails."""
    if smiles is None:
        return None
    s = str(smiles).strip()
    if not s or s.lower() == "nan":
        return None
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)
