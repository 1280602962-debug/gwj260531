#!/usr/bin/env python3
"""Stage-compare V4.2 HOLD residues. Read-only except writing QA diagnostics."""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_1.phase2_prepare_receptors import CIF, AA, parse_mmcif, coherent_altloc  # noqa: E402
from analysis.uniform_protocol_config import PRIMARY_PROTEIN_AUTH  # noqa: E402
from meeko.utils.rdkitutils import covalent_radius
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

RUN = ROOT / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
REC = RUN / "01_receptors"
QA = RUN / "13_qa"
ALLOWANCE = 1.2
ELEM_Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15, "SE": 34}


def parse_pdb(path: Path) -> list[dict]:
    atoms = []
    if not path.is_file():
        return atoms
    for ln in path.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        atoms.append({
            "name": ln[12:16].strip(),
            "alt": ln[16].strip() or ".",
            "comp": ln[17:20].strip(),
            "auth_asym": ln[21].strip(),
            "auth_seq": int(ln[22:26]),
            "ins": ln[26].strip() if len(ln) > 26 else "",
            "x": float(ln[30:38]),
            "y": float(ln[38:46]),
            "z": float(ln[46:54]),
            "occ": float(ln[54:60]) if ln[54:60].strip() else 1.0,
            "elem": (ln[76:78].strip() if len(ln) >= 78 else ln[12:16].strip()[:1]) or "C",
            "line": ln,
        })
    return atoms


def dist(a, b):
    return math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))


def atom_rec(a):
    return {
        "chain": a["auth_asym"],
        "seq": a["auth_seq"],
        "ins": a.get("ins", ""),
        "comp": a["comp"],
        "name": a["name"],
        "alt": a.get("alt", "."),
        "occ": a.get("occ", 1.0),
        "elem": a["elem"],
        "xyz": [round(a["x"], 3), round(a["y"], 3), round(a["z"], 3)],
    }


def meeko_thresh(e1, e2):
    z1, z2 = ELEM_Z.get(e1.upper(), 6), ELEM_Z.get(e2.upper(), 6)
    r1, r2 = covalent_radius.get(z1, 0.76), covalent_radius.get(z2, 0.76)
    return ALLOWANCE * (r1 + r2)


def residue(atoms, ch, seq):
    return [a for a in atoms if a["auth_asym"] == ch and a["auth_seq"] == seq]


def closest_pair(a_atoms, b_atoms):
    best = None
    for i in a_atoms:
        for j in b_atoms:
            d = dist(i, j)
            if best is None or d < best[0]:
                best = (d, i, j)
    return best


def neighbors_of(atoms, target, cutoff=2.5):
    out = []
    for a in atoms:
        if a is target:
            continue
        d = dist(a, target)
        if d <= cutoff:
            out.append((round(d, 3), atom_rec(a)))
    return sorted(out, key=lambda x: x[0])


def chain_continuity(atoms, ch, seq):
    by_seq = defaultdict(dict)
    for a in atoms:
        if a["auth_asym"] == ch:
            by_seq[a["auth_seq"]][a["name"]] = a
    seqs = sorted(by_seq)
    prevs = [s for s in seqs if s < seq]
    nexts = [s for s in seqs if s > seq]
    prev, nxt = (prevs[-1] if prevs else None), (nexts[0] if nexts else None)
    rec = {"chain": ch, "seq": seq, "prev_seq": prev, "next_seq": nxt}
    if prev is not None and "C" in by_seq[prev] and "N" in by_seq[seq]:
        rec["prev_CN"] = round(dist(by_seq[prev]["C"], by_seq[seq]["N"]), 3)
    if nxt is not None and "C" in by_seq[seq] and "N" in by_seq[nxt]:
        rec["next_CN"] = round(dist(by_seq[seq]["C"], by_seq[nxt]["N"]), 3)
    rec["is_internal_gap_prev"] = bool(prev is not None and seq - prev > 1)
    rec["is_internal_gap_next"] = bool(nxt is not None and nxt - seq > 1)
    rec["present_atoms"] = [a["name"] for a in residue(atoms, ch, seq)]
    return rec


def rdkit_residue_valence(atoms, ch, seq):
    recs = residue(atoms, ch, seq)
    if not recs:
        return {"error": "missing_residue"}
    block = "".join(
        (
            f"ATOM  {i+1:5d} {a['name']:>4s} {a['comp']:>3s} {a['auth_asym']}"
            f"{a['auth_seq']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            f"{a['occ']:6.2f}  1.00          {a['elem']:>2s}\n"
        )
        for i, a in enumerate(recs)
    ) + "END\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    if mol is None:
        return {"error": "MolFromPDBBlock_failed", "atoms": [atom_rec(a) for a in recs]}
    try:
        rdDetermineBonds.DetermineConnectivity(mol)
    except Exception as e:
        return {"error": f"DetermineConnectivity:{e}", "atoms": [atom_rec(a) for a in recs]}
    bonds = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        ai, aj = recs[i], recs[j]
        bonds.append({
            "a": f"{ai['name']}",
            "b": f"{aj['name']}",
            "d": round(dist(ai, aj), 3),
        })
    valences = []
    for i, atom in enumerate(mol.GetAtoms()):
        valences.append({
            "idx": i,
            "name": recs[i]["name"],
            "elem": recs[i]["elem"],
            "degree": atom.GetDegree(),
            "neighbors": [recs[n.GetIdx()]["name"] for n in atom.GetNeighbors()],
        })
    sanitize_error = None
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        sanitize_error = str(e)
    return {
        "n_atoms": len(recs),
        "atoms": [atom_rec(a) for a in recs],
        "determine_bonds": bonds,
        "valences": valences,
        "sanitize_error": sanitize_error,
    }


def stage_paths(pdb):
    return {
        "protein_raw": REC / f"{pdb}_protein_raw.pdb",
        "fixer_openmm": REC / f"{pdb}_protein_fixer_openmm.pdb",
        "fixer_cleaned": REC / f"{pdb}_protein_fixer.pdb",
        "ptr_fallback": REC / f"{pdb}_protein_ptr_fallback.pdb",
        "prepared": REC / f"{pdb}_protein_prepared.pdb",
    }


def cif_atoms(pdb, auths):
    atoms = parse_mmcif((CIF / f"{pdb}.cif").read_text(errors="replace"))
    chosen, kept = coherent_altloc(atoms, set(auths))
    return atoms, chosen, kept


def pair_stage_report(atoms, ch1, seq1, n1, ch2, seq2, n2):
    a = next((x for x in residue(atoms, ch1, seq1) if x["name"] == n1), None)
    b = next((x for x in residue(atoms, ch2, seq2) if x["name"] == n2), None)
    pair = closest_pair(residue(atoms, ch1, seq1), residue(atoms, ch2, seq2))
    rec = {
        "res1": chain_continuity(atoms, ch1, seq1),
        "res2": chain_continuity(atoms, ch2, seq2),
        "named_pair": None,
        "closest_pair": None,
        "meeko_would_bond_named": None,
        "meeko_would_bond_closest": None,
    }
    if a and b:
        d = dist(a, b)
        thr = meeko_thresh(a["elem"], b["elem"])
        rec["named_pair"] = {
            "a": atom_rec(a),
            "b": atom_rec(b),
            "distance": round(d, 3),
            "meeko_threshold": round(thr, 3),
            "meeko_covalent_hit": d < thr,
        }
        rec["meeko_would_bond_named"] = d < thr
    if pair:
        d, i, j = pair
        thr = meeko_thresh(i["elem"], j["elem"])
        rec["closest_pair"] = {
            "a": atom_rec(i),
            "b": atom_rec(j),
            "distance": round(d, 3),
            "meeko_threshold": round(thr, 3),
            "meeko_covalent_hit": d < thr,
        }
        rec["meeko_would_bond_closest"] = d < thr
    return rec


def cif_pair_all_altlocs(raw_cif, ch1, seq1, n1, ch2, seq2, n2):
    a_all = [a for a in raw_cif if a["auth_asym"] == ch1 and a["auth_seq"] == seq1]
    b_all = [a for a in raw_cif if a["auth_asym"] == ch2 and a["auth_seq"] == seq2]
    named = []
    for a in a_all:
        if a["name"] != n1:
            continue
        for b in b_all:
            if b["name"] != n2:
                continue
            named.append({
                "a": atom_rec(a),
                "b": atom_rec(b),
                "distance": round(dist(a, b), 3),
                "meeko_threshold": round(meeko_thresh(a["elem"], b["elem"]), 3),
                "meeko_covalent_hit": dist(a, b) < meeko_thresh(a["elem"], b["elem"]),
            })
    closest = []
    for a in a_all:
        for b in b_all:
            d = dist(a, b)
            closest.append((d, a, b))
    closest.sort(key=lambda x: x[0])
    top = []
    for d, a, b in closest[:8]:
        top.append({
            "a": atom_rec(a),
            "b": atom_rec(b),
            "distance": round(d, 3),
            "meeko_covalent_hit": d < meeko_thresh(a["elem"], b["elem"]),
        })
    comps = sorted({a["comp"] for a in a_all} | {b["comp"] for b in b_all})
    return {
        "residue_comps": comps,
        "named_pair_all_altlocs": named,
        "closest_contacts_any_atom": top,
        "res1_altlocs": sorted({a.get("alt", ".") for a in a_all}),
        "res2_altlocs": sorted({a.get("alt", ".") for a in b_all}),
        "res1_atoms": [atom_rec(a) for a in a_all],
        "res2_atoms": [atom_rec(a) for a in b_all],
    }


def coord_delta(a, b):
    if a is None or b is None:
        return None
    return {
        "dx": round(b["x"] - a["x"], 3),
        "dy": round(b["y"] - a["y"], 3),
        "dz": round(b["z"] - a["z"], 3),
        "d": round(dist(a, b), 3),
    }


def main():
    report = {
        "RUN_ID": "UNIFORM_RERUN_V4_2_20260921",
        "meeko_allowance": ALLOWANCE,
        "meeko_covalent_radius_C_N_O_S": {
            "C": covalent_radius[6],
            "N": covalent_radius[7],
            "O": covalent_radius[8],
            "S": covalent_radius[16],
        },
        "meeko_thresholds": {
            "C-N": round(meeko_thresh("C", "N"), 3),
            "C-O": round(meeko_thresh("C", "O"), 3),
            "N-O": round(meeko_thresh("N", "O"), 3),
        },
        "holds": {},
    }

    cases = {
        "6N7A": {
            "kind": "extra_interresidue_bond",
            "pair": ("A", 1024, "CD1", "A", 1041, "NH2"),
            "ptr": True,
        },
        "8BXH": {
            "kind": "extra_interresidue_bond",
            "pair": ("A", 839, "CE", "A", 918, "OH"),
            "ptr": True,
        },
        "4UDW": {
            "kind": "intraresidue_valence",
            "res": ("I", 565),
            "ptr": False,
        },
    }

    for pdb, spec in cases.items():
        auths = PRIMARY_PROTEIN_AUTH[pdb]
        raw_cif, chosen_alt, _ = cif_alt = cif_atoms(pdb, auths)
        stages = {k: parse_pdb(v) for k, v in stage_paths(pdb).items()}
        rec = {
            "primary_auth": auths,
            "altloc_choices_for_query_residues": {},
            "stages": {},
            "origin": {},
        }
        if spec["kind"] == "extra_interresidue_bond":
            ch1, s1, n1, ch2, s2, n2 = spec["pair"]
            rec["query"] = {"res1": f"{ch1}:{s1}:{n1}", "res2": f"{ch2}:{s2}:{n2}"}
            rec["cif_all_altlocs"] = cif_pair_all_altlocs(raw_cif, ch1, s1, n1, ch2, s2, n2)
            for key, alts in chosen_alt.items():
                if key[0] in {ch1, ch2} and key[1] in {s1, s2}:
                    rec["altloc_choices_for_query_residues"][f"{key[0]}:{key[1]}:{key[3]}"] = alts
            for name, atoms in [("frozen_cif_all_altlocs", raw_cif)] + list(stages.items()):
                if not atoms:
                    rec["stages"][name] = {"missing": True}
                    continue
                rec["stages"][name] = pair_stage_report(atoms, ch1, s1, n1, ch2, s2, n2)
            # coordinate drift of named atoms raw -> prepared
            raw_a = next((x for x in residue(stages["protein_raw"], ch1, s1) if x["name"] == n1), None)
            raw_b = next((x for x in residue(stages["protein_raw"], ch2, s2) if x["name"] == n2), None)
            prep_a = next((x for x in residue(stages["prepared"], ch1, s1) if x["name"] == n1), None)
            prep_b = next((x for x in residue(stages["prepared"], ch2, s2) if x["name"] == n2), None)
            rec["coord_drift_raw_to_prepared"] = {
                f"{ch1}:{s1}:{n1}": coord_delta(raw_a, prep_a),
                f"{ch2}:{s2}:{n2}": coord_delta(raw_b, prep_b),
            }
            cif_hits = rec["cif_all_altlocs"]["named_pair_all_altlocs"]
            cif_min = min((x["distance"] for x in cif_hits), default=None)
            prep_d = rec["stages"]["prepared"].get("named_pair", {})
            rec["origin"] = {
                "contact_in_frozen_cif": bool(cif_hits and any(x["meeko_covalent_hit"] for x in cif_hits)),
                "min_named_distance_cif": cif_min,
                "named_distance_prepared": None if not prep_d else prep_d.get("distance"),
                "introduced_by_preparation": bool(
                    cif_min is not None
                    and prep_d
                    and (not any(x["meeko_covalent_hit"] for x in cif_hits))
                    and prep_d.get("meeko_covalent_hit")
                ),
                "is_disulfide": False,
                "is_polymer_break": False,
                "is_real_covalent_chemistry": False,
                "template_or_boundary_can_represent_true_chemistry": False,
                "reason": (
                    "Meeko 0.7.1 find_inter_mols_bonds treats this non-canonical close contact "
                    "as a covalent inter-residue bond (allowance 1.2 * covalent radii). "
                    "It is not a disulfide and not a polymer break. "
                    "CYX / blunt-end / standard residue templates cannot represent a false C–N/C–O sidechain bond."
                ),
            }
        else:
            ch, seq = spec["res"]
            rec["query"] = {"res": f"{ch}:{seq}"}
            cif_res = [a for a in raw_cif if a["auth_asym"] == ch and a["auth_seq"] == seq]
            rec["cif_residue_all_altlocs"] = [atom_rec(a) for a in cif_res]
            rec["cif_altlocs"] = sorted({a.get("alt", ".") for a in cif_res})
            rec["cif_comps"] = sorted({a["comp"] for a in cif_res})
            for key, alts in chosen_alt.items():
                if key[0] == ch and key[1] == seq:
                    rec["altloc_choices_for_query_residues"][f"{key[0]}:{key[1]}:{key[3]}"] = alts
            intra = []
            for i, a in enumerate(cif_res):
                for b in cif_res[i + 1 :]:
                    d = dist(a, b)
                    if d <= 2.4:
                        intra.append({
                            "a": atom_rec(a),
                            "b": atom_rec(b),
                            "distance": round(d, 3),
                        })
            rec["cif_intraresidue_le_2.4A"] = sorted(intra, key=lambda x: x["distance"])
            for name, atoms in [("frozen_cif_all_altlocs", raw_cif)] + list(stages.items()):
                if not atoms:
                    rec["stages"][name] = {"missing": True}
                    continue
                rec["stages"][name] = {
                    "continuity": chain_continuity(atoms, ch, seq),
                    "rdkit_meeko_like": rdkit_residue_valence(atoms, ch, seq),
                }
            rec["nearby_to_OE1_prepared"] = []
            prep = stages["prepared"]
            oe1 = next((x for x in residue(prep, ch, seq) if x["name"] == "OE1"), None)
            if oe1:
                rec["nearby_to_OE1_prepared"] = neighbors_of(prep, oe1, 2.5)
            cif_val = rec["stages"]["frozen_cif_all_altlocs"]["rdkit_meeko_like"]
            prep_val = rec["stages"]["prepared"]["rdkit_meeko_like"]
            rec["origin"] = {
                "valence_error_in_frozen_cif_residue": bool(cif_val.get("sanitize_error")),
                "valence_error_in_prepared": bool(prep_val.get("sanitize_error")),
                "introduced_by_preparation": bool(
                    (not cif_val.get("sanitize_error")) and prep_val.get("sanitize_error")
                ),
                "is_disulfide": False,
                "is_polymer_break": rec["stages"]["prepared"]["continuity"].get("is_internal_gap_prev")
                or rec["stages"]["prepared"]["continuity"].get("is_internal_gap_next"),
                "component": "I hirudin variant-2 RETAIN",
                "template_or_boundary_can_represent_true_chemistry": False,
                "reason": (
                    "RDKit DetermineConnectivity + SanitizeMol on GLN I:565 (Meeko per-residue path). "
                    "Not a disulfide. Not a polymer-boundary fact that a CYX/blunt-end template can correct. "
                    "Deleting RETAIN I is forbidden."
                ),
            }
        report["holds"][pdb] = rec

    out = QA / "HOLD_STAGE_COMPARISON.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        pdb: {
            "query": rec["query"],
            "origin": rec["origin"],
            "altloc": rec.get("altloc_choices_for_query_residues"),
        }
        for pdb, rec in report["holds"].items()
    }, indent=2))
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
