#!/usr/bin/env python3
"""Experimental-structure pocket comparison. No docking scores or AUROC."""
from __future__ import annotations

import csv
import json
import math
import urllib.request
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np

CIFDIR = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa/alt_cif")
QA = Path("/tmp/pr39_fiveseed/Dual_Target_Docking/reruns/UNIFORM_RERUN_V4_2_20260921/13_qa")
CACHE = Path("/tmp/input_identity_rcsb")

# Primary copy used for pocket definition (V4.2 auth policy).
PRIMARY = {
    "EGFR": {"pdb": "3POZ", "auth": "A", "lig_auth": "A", "lig_ccd": "03P"},
    "HER2": {"pdb": "3RCD", "auth": "A", "lig_auth": "A", "lig_ccd": "03P"},
    "JAK1": {"pdb": "6N7A", "auth": "A", "lig_auth": "A", "lig_ccd": "KEV"},
    "JAK2": {"pdb": "8BXH", "auth": "A", "lig_auth": "A", "lig_ccd": "C87"},
    "TYK2": {"pdb": "3LXP", "auth": "A", "lig_auth": "A", "lig_ccd": "IZA"},
    "PIK3CA": {"pdb": "4L23", "auth": "A", "lig_auth": "A", "lig_ccd": "X6K"},
    "mTOR": {"pdb": "4JT6", "auth": "A", "lig_auth": "A", "lig_ccd": "X6K"},
    "AChE": {"pdb": "4EY7", "auth": "A", "lig_auth": "A", "lig_ccd": "E20"},
    "BChE": {"pdb": "4BDS", "auth": "A", "lig_auth": "A", "lig_ccd": "THA"},
    "F2": {"pdb": "4UDW", "auth": "H", "lig_auth": "H", "lig_ccd": "N6L"},
    "F10": {"pdb": "2JKH", "auth": "A", "lig_auth": "A", "lig_ccd": "BI7"},
    "PPARG": {"pdb": "9V8H", "auth": "A", "lig_auth": "A", "lig_ccd": "BRL"},
    "PPARA": {"pdb": "6LXA", "auth": "A", "lig_auth": "A", "lig_ccd": "EPA"},
    "PPARD": {"pdb": "5U3Q", "auth": "B", "lig_auth": "B", "lig_ccd": "7UJ"},
}

# Representative alternatives for RMSD (construct/site first; not docking-ranked).
# Includes all scarce HER2 kinase holos and the 4JT6 crystallization series.
ALTS = {
    "EGFR": ["1M17", "4HJO", "2ITY", "1XKK", "3W2S", "5UGC"],
    "HER2": ["3PP0", "7JXH", "7PCD", "8VB5", "8U8X"],
    "JAK1": ["4EI4", "4IVC", "5E1E", "6DBN"],
    "JAK2": ["3FUP", "4C61", "5CF4", "6BBV"],
    "TYK2": ["3NZ0", "3LXN", "6AAM", "4GVJ"],
    "PIK3CA": ["4JPS", "5DXT", "4J6I", "5XGH", "7K6M"],
    "mTOR": ["4JSX", "4JSP", "4JT5", "4JSV"],
    "AChE": ["4EY6", "4EY5", "1B41", "4M0E", "6O4W"],
    "BChE": ["1P0M", "6EMI", "3DJY", "4TPK"],
    "F2": ["1PPB", "2ZFF", "1H8D", "3SHC", "5AFY"],
    "F10": ["2P16", "1FAX", "2J94", "2W26"],
    "PPARG": ["7AWC", "2PRG", "1FM6", "3DZY"],
    "PPARA": ["1I7G", "2P54", "3KDU", "6K1O"],
    "PPARD": ["3GZ9", "5U46", "2ZNP", "3TKM"],
}

POCKET_CUTOFF = 4.5


def download_cif(pdb: str) -> Path:
    dest = CIFDIR / f"{pdb.upper()}.cif"
    if dest.is_file() and dest.stat().st_size > 1000:
        return dest
    src = CACHE / f"{pdb.upper()}.cif"
    if src.is_file() and src.stat().st_size > 1000:
        dest.write_bytes(src.read_bytes())
        return dest
    url = f"https://files.rcsb.org/download/{pdb.upper()}.cif"
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def load_doc(pdb: str) -> gemmi.cif.Document:
    return gemmi.cif.read(str(download_cif(pdb)))


def seq_difs(doc: gemmi.cif.Document) -> list[str]:
    block = doc.sole_block()
    tab = block.find("_struct_ref_seq_dif.", ["pdbx_pdb_strand_id", "mon_id", "db_mon_id", "db_align_seq_id", "details"])
    out = []
    if not tab:
        return out
    for row in tab:
        ch, mon, dbmon, seq, det = [gemmi.cif.as_string(x) if x else "" for x in row]
        if mon and dbmon and mon != dbmon:
            out.append(f"{ch}:{dbmon}{seq}{mon}")
        elif det and det not in (".", "?"):
            out.append(f"{ch}:{det}:{seq}")
    return out


def missing_res(doc: gemmi.cif.Document) -> list[str]:
    block = doc.sole_block()
    tab = block.find("_pdbx_unobs_or_zero_occ_residues.", ["auth_asym_id", "auth_comp_id", "auth_seq_id", "occupancy_flag"])
    out = []
    if not tab:
        return out
    for row in tab:
        ch, aa, seq, flag = [gemmi.cif.as_string(x) if x else "" for x in row]
        out.append(f"{ch}:{aa}{seq}")
    return out


def ligand_coords(st: gemmi.Structure, auth_chain: str, ccd: str) -> np.ndarray:
    xs = []
    model = st[0]
    for ch in model:
        if ch.name != auth_chain:
            continue
        for res in ch:
            if res.name == ccd:
                for atom in res:
                    if atom.element.name != "H":
                        xs.append(atom.pos.tolist())
    return np.array(xs, dtype=float) if xs else np.zeros((0, 3))


def residue_bb(res: gemmi.Residue) -> dict[str, np.ndarray]:
    out = {}
    for atom in res:
        if atom.name in {"N", "CA", "C", "O"} and atom.element.name != "H":
            out[atom.name] = np.array(atom.pos.tolist(), dtype=float)
    return out


def residue_heavy(res: gemmi.Residue) -> dict[str, np.ndarray]:
    out = {}
    for atom in res:
        if atom.element.name != "H":
            out[atom.name] = np.array(atom.pos.tolist(), dtype=float)
    return out


def chain_residues(st: gemmi.Structure, auth: str) -> dict[int, gemmi.Residue]:
    rec = {}
    for ch in st[0]:
        if ch.name != auth:
            continue
        for res in ch:
            if res.entity_type != gemmi.EntityType.Polymer:
                continue
            rec[res.seqid.num] = res
    return rec


def pocket_seqids(st: gemmi.Structure, auth: str, lig_xyz: np.ndarray) -> list[int]:
    if lig_xyz.size == 0:
        return []
    hits = []
    for seq, res in chain_residues(st, auth).items():
        for atom in res:
            if atom.element.name == "H":
                continue
            d = np.linalg.norm(np.array(atom.pos.tolist()) - lig_xyz, axis=1).min()
            if d <= POCKET_CUTOFF:
                hits.append(seq)
                break
    return sorted(hits)


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[float, np.ndarray]:
    """Return RMSD after aligning Q onto P."""
    if len(P) < 3:
        return float("nan"), Q
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Qc.T @ Pc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    Qr = Qc @ R
    rmsd = math.sqrt(((Pc - Qr) ** 2).sum() / len(P))
    return rmsd, Qr + P.mean(0)


def compare(target: str, alt_pdb: str, alt_auth: str | None, pocket: list[int], prim_st, prim_auth: str) -> dict:
    st = gemmi.read_structure(str(download_cif(alt_pdb)))
    st.setup_entities()
    # choose polymer chain with most pocket residues present
    best = None
    best_n = -1
    cand_chains = [alt_auth] if alt_auth else []
    if not cand_chains:
        cand_chains = [ch.name for ch in st[0]]
    for chn in cand_chains:
        resmap = chain_residues(st, chn)
        n = sum(1 for i in pocket if i in resmap)
        if n > best_n:
            best_n, best = n, chn
    res_p = chain_residues(prim_st, prim_auth)
    res_a = chain_residues(st, best)
    common = [i for i in pocket if i in res_p and i in res_a]
    missing_in_alt = [i for i in pocket if i not in res_a]
    # backbone
    P, Q, names = [], [], []
    for i in common:
        bp, ba = residue_bb(res_p[i]), residue_bb(res_a[i])
        for at in ("N", "CA", "C"):
            if at in bp and at in ba:
                P.append(bp[at]); Q.append(ba[at]); names.append((i, at))
    bb_rmsd, Qr = kabsch(np.array(P), np.array(Q)) if P else (float("nan"), None)
    # sidechain heavy after same rotation applied to all mapped atoms
    sc_d = []
    rotamer = []
    if Qr is not None and P:
        # rebuild transform from kabsch internals using CA only for residue-level
        Pc = np.array(P)
        for i in common:
            hp, ha = residue_heavy(res_p[i]), residue_heavy(res_a[i])
            # chi1-ish: CB if present
            if "CB" in hp and "CB" in ha and "CA" in hp and "CA" in ha:
                # after global bb alignment, compare CB direction
                pass
            for at, xyzp in hp.items():
                if at in {"N", "CA", "C", "O"}:
                    continue
                if at in ha:
                    sc_d.append((at, xyzp, ha[at]))
        # apply same R,t as backbone
        Pm = np.array(P)
        Qm = np.array(Q)
        Pc = Pm - Pm.mean(0)
        Qc = Qm - Qm.mean(0)
        H = Qc.T @ Pc
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1] *= -1
            R = Vt.T @ U.T
        t = Pm.mean(0) - Qm.mean(0) @ R
        diffs = []
        for at, xyzp, xyza in sc_d:
            aligned = xyza @ R + t
            diffs.append(np.linalg.norm(xyzp - aligned))
        sc_rmsd = math.sqrt(float(np.mean(np.square(diffs)))) if diffs else float("nan")
        # rotamer flags: CB deviation > 1.2 Å
        for i in common:
            hp, ha = residue_heavy(res_p[i]), residue_heavy(res_a[i])
            if "CB" in hp and "CB" in ha:
                aligned = ha["CB"] @ R + t
                d = float(np.linalg.norm(hp["CB"] - aligned))
                if d > 1.2:
                    rotamer.append(f"{res_p[i].name}{i}:{d:.2f}")
    else:
        sc_rmsd = float("nan")
    doc = load_doc(alt_pdb)
    return {
        "target": target,
        "primary_pdb": PRIMARY[target]["pdb"],
        "alt_pdb": alt_pdb,
        "alt_chain_used": best,
        "n_pocket_primary": len(pocket),
        "n_pocket_common": len(common),
        "n_pocket_missing_in_alt": len(missing_in_alt),
        "pocket_completeness": round(len(common) / len(pocket), 3) if pocket else "",
        "pocket_bb_rmsd_A": round(bb_rmsd, 3) if bb_rmsd == bb_rmsd else "",
        "pocket_sidechain_heavy_rmsd_A": round(sc_rmsd, 3) if sc_rmsd == sc_rmsd else "",
        "rotamer_changes_CB_gt_1.2A": ";".join(rotamer[:20]),
        "n_rotamer_changes": len(rotamer),
        "seq_dif": ";".join(seq_difs(doc)[:25]),
        "n_seq_dif": len(seq_difs(doc)),
        "n_missing_res_records": len(missing_res(doc)),
        "missing_pocket_seqids": ";".join(map(str, missing_in_alt[:20])),
    }


def infer_alt_chain(target: str, pdb: str) -> str | None:
    # defaults: first polymer often A; PPARD 5U3Q-like use B when present
    if target == "PPARD" and pdb == "5U3Q":
        return "B"
    if target == "F2":
        return "H" if pdb == "4UDW" else None
    if target == "mTOR":
        return "A"
    return None


def main() -> None:
    CIFDIR.mkdir(parents=True, exist_ok=True)
    pockets = {}
    prim_st = {}
    for target, spec in PRIMARY.items():
        print("primary", target, spec["pdb"], flush=True)
        path = download_cif(spec["pdb"])
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        prim_st[target] = st
        lig = ligand_coords(st, spec["lig_auth"], spec["lig_ccd"])
        pk = pocket_seqids(st, spec["auth"], lig)
        pockets[target] = pk
        print(" ", spec["pdb"], "lig_atoms", len(lig), "pocket_res", len(pk), flush=True)

    rows = []
    for target, alts in ALTS.items():
        for pdb in alts:
            print("compare", target, pdb, flush=True)
            try:
                rec = compare(target, pdb, infer_alt_chain(target, pdb), pockets[target], prim_st[target], PRIMARY[target]["auth"])
                rows.append(rec)
            except Exception as exc:
                rows.append({
                    "target": target, "primary_pdb": PRIMARY[target]["pdb"], "alt_pdb": pdb,
                    "error": str(exc),
                })
    out = QA / "ALTERNATIVE_RECEPTOR_STRUCTURAL_COMPARISON.csv"
    fields = sorted({k for r in rows for k in r})
    with out.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    (QA / "alt_receptor_pockets.json").write_text(json.dumps({k: v for k, v in pockets.items()}, indent=2), encoding="utf-8")
    print("wrote", out, "n", len(rows))


if __name__ == "__main__":
    main()
