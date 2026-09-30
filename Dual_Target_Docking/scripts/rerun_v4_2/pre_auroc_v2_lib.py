#!/usr/bin/env python3
"""Shared parsers for PRE_AUROC V2. No AUROC. CGn is a real atom; Gn is glue only."""
from __future__ import annotations

import math
import re
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, rdmolops

RDLogger.DisableLog("rdApp.*")

STANDARD_AD4 = {
    "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F",
    "Mg", "MG", "P", "SA", "S", "Cl", "CL", "Ca", "CA",
    "Mn", "MN", "Fe", "FE", "Zn", "ZN", "Br", "BR", "I", "Z",
}

GLUE_RE = re.compile(r"^G\d+$")
REAL_MACRO_RE = re.compile(r"^[A-Z]G\d+$")  # CGn; also NGn etc. if present
CGN_RE = re.compile(r"^CG\d+$")


def atom_type(line: str) -> str:
    if len(line) >= 78:
        return line[77:].strip()
    parts = line.split()
    return parts[-1] if parts else ""


def is_glue_type(typ: str) -> bool:
    return typ == "G" or bool(GLUE_RE.match(typ))


def is_real_macrocycle_type(typ: str) -> bool:
    if is_glue_type(typ):
        return False
    return bool(REAL_MACRO_RE.match(typ))


def parse_pdbqt_atoms(text: str) -> list[dict]:
    atoms = []
    for ln in text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        serial = int(ln[6:11])
        name = ln[12:16].strip()
        typ = atom_type(ln)
        x = float(ln[30:38])
        y = float(ln[38:46])
        z = float(ln[46:54])
        atoms.append({
            "ln": ln,
            "serial": serial,
            "name": name,
            "type": typ,
            "x": x, "y": y, "z": z,
            "is_glue": is_glue_type(typ),
            "is_real_macro": is_real_macrocycle_type(typ),
        })
    return atoms


def all_atom_types(atoms: list[dict]) -> list[str]:
    return sorted({a["type"] for a in atoms if a["type"]})


def special_types(atoms: list[dict]) -> list[str]:
    return sorted({a["type"] for a in atoms if a["type"] not in STANDARD_AD4})


def parse_smiles_idx(text: str) -> tuple[str, dict[int, int]]:
    smiles = ""
    pairs: list[int] = []
    for ln in text.splitlines():
        if ln.startswith("REMARK SMILES IDX"):
            pairs.extend(int(x) for x in ln.split()[3:])
        elif ln.startswith("REMARK SMILES") and "IDX" not in ln:
            smiles = ln[len("REMARK SMILES "):].strip()
    if len(pairs) % 2:
        return smiles, {}
    mapping = {}
    for i in range(0, len(pairs), 2):
        mapping[pairs[i]] = pairs[i + 1]  # smiles_idx 1-based -> pdbqt serial
    return smiles, mapping


def invert_smiles_idx(idx: dict[int, int]) -> dict[int, int]:
    inv = {}
    for smi, serial in idx.items():
        if serial in inv:
            return {}
        inv[serial] = smi
    return inv


def parse_h_parent(text: str) -> list[tuple[int, int]]:
    pairs: list[int] = []
    for ln in text.splitlines():
        if ln.startswith("REMARK H PARENT"):
            pairs.extend(int(x) for x in ln.split()[3:])
    if len(pairs) % 2:
        return []
    return [(pairs[i], pairs[i + 1]) for i in range(0, len(pairs), 2)]


def split_models(pdbqt_text: str) -> list[str]:
    chunks, cur = [], []
    for ln in pdbqt_text.splitlines(keepends=True):
        if ln.startswith("MODEL") and cur:
            chunks.append("".join(cur))
            cur = [ln]
        else:
            cur.append(ln)
    if cur:
        chunks.append("".join(cur))
    cleaned = []
    for chunk in chunks:
        keep = []
        for ln in chunk.splitlines(keepends=True):
            if ln.startswith(("MODEL", "ENDMDL")):
                continue
            keep.append(ln)
        text = "".join(keep)
        if "ATOM" in text or "HETATM" in text:
            cleaned.append(text)
    return cleaned or ([pdbqt_text] if ("ATOM" in pdbqt_text or "HETATM" in pdbqt_text) else [])


def parse_torsion_tree(text: str) -> dict:
    """Rigid blocks + explicit BRANCH edges. Does not invent bonds from coordinates."""
    blocks: list[list[int]] = []
    edges: set[tuple[int, int]] = set()
    current: list[int] = []
    in_block = False
    for ln in text.splitlines():
        if ln.startswith("ROOT") or ln.startswith("BRANCH"):
            if current:
                blocks.append(current)
                current = []
            in_block = True
            if ln.startswith("BRANCH"):
                parts = ln.split()
                if len(parts) >= 3:
                    a, b = int(parts[1]), int(parts[2])
                    edges.add((min(a, b), max(a, b)))
        elif ln.startswith(("ENDROOT", "ENDBRANCH", "TORSDOF")):
            if current:
                blocks.append(current)
                current = []
            in_block = False
        elif in_block and ln.startswith(("ATOM", "HETATM")):
            current.append(int(ln[6:11]))
    if current:
        blocks.append(current)
    serial_to_block = {}
    for i, block in enumerate(blocks):
        for s in block:
            serial_to_block[s] = i
    return {"blocks": blocks, "branch_edges": edges, "serial_to_block": serial_to_block}


def same_rigid_block(tree: dict, a: int, b: int) -> bool:
    ia = tree["serial_to_block"].get(a)
    ib = tree["serial_to_block"].get(b)
    return ia is not None and ia == ib


def has_branch_edge(tree: dict, a: int, b: int) -> bool:
    return (min(a, b), max(a, b)) in tree["branch_edges"]


def coord_key(atom: dict, nd: int = 3) -> tuple[float, float, float]:
    return (round(atom["x"], nd), round(atom["y"], nd), round(atom["z"], nd))


def complementary_cgn_pairs(atoms: list[dict]) -> list[dict]:
    """Pair real CGn atoms using Gn glue sitting on the complementary real atom."""
    cgn = [a for a in atoms if a["is_real_macro"]]
    gn = [a for a in atoms if a["is_glue"]]
    if not cgn:
        return []
    by_coord: dict[tuple[float, float, float], list[dict]] = {}
    for a in cgn:
        by_coord.setdefault(coord_key(a), []).append(a)
    pairs = []
    used = set()
    for g in gn:
        gkey = coord_key(g)
        hosts = by_coord.get(gkey, [])
        if len(hosts) != 1:
            continue
        other = hosts[0]
        # find the CGn that owns this G: typically same rigid neighborhood / nearest CGn
        owner = min(
            (c for c in cgn if c["serial"] != other["serial"]),
            key=lambda c: math.dist((c["x"], c["y"], c["z"]), (g["x"], g["y"], g["z"])),
            default=None,
        )
        if owner is None:
            continue
        key = tuple(sorted((owner["serial"], other["serial"])))
        if key in used or key[0] == key[1]:
            continue
        used.add(key)
        d = math.dist((owner["x"], owner["y"], owner["z"]), (other["x"], other["y"], other["z"]))
        pairs.append({
            "serial_a": key[0],
            "serial_b": key[1],
            "glue_serials": sorted({g["serial"] for g in gn if coord_key(g) in {coord_key(owner), coord_key(other)}}),
            "cartesian_distance_A": round(d, 4),
            "pairing_method": "complementary_Gn_on_other_CGn_coordinates",
        })
    return pairs


def load_frozen_sdf(sdf_path: Path) -> Chem.Mol | None:
    if not sdf_path.is_file():
        return None
    return Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=True)


def sdf_heavy(mol: Chem.Mol) -> Chem.Mol:
    return Chem.Mol(rdmolops.RemoveHs(mol))


def sdf_heavy_to_tmpl_match(heavy: Chem.Mol, smiles: str) -> tuple[tuple[int, ...] | None, str]:
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        return None, "ATOM_MAPPING_FAILURE_remark_smiles"
    if heavy.GetNumAtoms() != tmpl.GetNumAtoms():
        return None, "ATOM_MAPPING_FAILURE_atom_count"
    order_ok = all(
        heavy.GetAtomWithIdx(i).GetAtomicNum() == tmpl.GetAtomWithIdx(i).GetAtomicNum()
        for i in range(tmpl.GetNumAtoms())
    )
    if order_ok:
        return tuple(range(tmpl.GetNumAtoms())), "ok_sequential"
    matches = heavy.GetSubstructMatches(tmpl, uniquify=True)
    if len(matches) != 1:
        return None, "ATOM_MAPPING_FAILURE_ambiguous_isomorphism"
    return matches[0], "ok_isomorphism"


def map_serial_to_sdf_heavy(
    frozen_pdbqt_text: str,
    sdf: Chem.Mol,
) -> tuple[dict[int, int], str]:
    """pdbqt serial -> frozen-SDF heavy-atom index. Gn excluded because they are not in SMILES IDX."""
    smiles, idx = parse_smiles_idx(frozen_pdbqt_text)
    if not smiles or not idx:
        return {}, "ATOM_MAPPING_FAILURE_no_smiles_idx"
    heavy = sdf_heavy(sdf)
    match, why = sdf_heavy_to_tmpl_match(heavy, smiles)
    if match is None:
        return {}, why
    serial_to_sdf: dict[int, int] = {}
    for smi_1, serial in idx.items():
        if smi_1 < 1 or smi_1 > len(match):
            return {}, "ATOM_MAPPING_FAILURE_idx_range"
        serial_to_sdf[serial] = match[smi_1 - 1]
    return serial_to_sdf, why


def sdf_has_bond(mol: Chem.Mol, i: int, j: int) -> bool:
    return mol.GetBondBetweenAtoms(i, j) is not None


def ring_count(mol: Chem.Mol | None) -> int | str:
    if mol is None:
        return ""
    try:
        return int(rdMolDescriptors.CalcNumRings(mol))
    except Exception:
        return ""


def formal_charge(mol: Chem.Mol | None) -> int | str:
    if mol is None:
        return ""
    return int(rdmolops.GetFormalCharge(mol))


def first_model_text(pdbqt_text: str) -> str:
    models = split_models(pdbqt_text)
    return models[0] if models else ""


def parse_model1_remarks(out_pdbqt_text: str) -> dict[str, str]:
    """Parse REMARK scores from the first MODEL only. Do not search the whole file/log."""
    text = first_model_text(out_pdbqt_text)
    out = {"CNNscore": "", "CNNaffinity": "", "CNN_VS": "", "minimizedAffinity": ""}
    if not text:
        return out
    for ln in text.splitlines():
        if not ln.startswith("REMARK"):
            continue
        if re.match(r"REMARK\s+CNNscore\b", ln):
            m = re.search(r"CNNscore\s+([0-9eE.+-]+)", ln)
            if m:
                out["CNNscore"] = m.group(1)
        elif re.match(r"REMARK\s+CNNaffinity\b", ln):
            m = re.search(r"CNNaffinity\s+([0-9eE.+-]+)", ln)
            if m:
                out["CNNaffinity"] = m.group(1)
        elif re.match(r"REMARK\s+CNN_?VS\b", ln):
            m = re.search(r"CNN_?VS\s+([0-9eE.+-]+)", ln)
            if m:
                out["CNN_VS"] = m.group(1)
        elif re.match(r"REMARK\s+minimizedAffinity\b", ln):
            m = re.search(r"minimizedAffinity\s+([0-9eE.+-]+)", ln)
            if m:
                out["minimizedAffinity"] = m.group(1)
    return out


def finite_float(s: str) -> bool:
    if s in ("", None, "--"):
        return False
    try:
        v = float(s)
    except (TypeError, ValueError):
        return False
    return math.isfinite(v)
