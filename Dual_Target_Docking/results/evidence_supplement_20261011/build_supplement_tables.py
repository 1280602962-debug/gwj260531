#!/usr/bin/env python3
"""Read-only supplement: stereo identity of 167 ligands and macrocycle atom maps.

Does not rewrite formal structures, labels, or scores.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "reruns/UNIFORM_RERUN_V4_2_20260921"
REG = RUN / "00_protocol/global_ligand_registry.csv"
PREV = ROOT / "results/methods_provenance_recovery_20261011/07_ligand_preparation_provenance_785.csv"
SDF_DIR = RUN / "02_ligands/sdf"
PDBQT_DIR = RUN / "02_ligands/pdbqt"
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.pre_auroc_v2_lib import (  # noqa: E402
    complementary_cgn_pairs,
    map_serial_to_sdf_heavy,
    parse_pdbqt_atoms,
    parse_smiles_idx,
    parse_torsion_tree,
    same_rigid_block,
    has_branch_edge,
    sdf_has_bond,
    sdf_heavy,
)

AD4_Z = {
    "C": 6, "A": 6, "CG0": 6,
    "N": 7, "NA": 7, "NS": 7,
    "O": 8, "OA": 8, "OS": 8,
    "S": 16, "SA": 16,
    "P": 15, "F": 9,
    "Cl": 17, "CL": 17,
    "Br": 35, "BR": 35,
    "I": 53,
    "H": 1, "HD": 1, "HS": 1,
}


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def heavy_from_3d(path: Path):
    mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)
    if mol is None:
        return None
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    return Chem.RemoveHs(mol)


def stereo_flags(mol: Chem.Mol) -> dict:
    Chem.AssignStereochemistry(mol, force=True, cleanIt=False)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assigned = [(i, tag) for i, tag in centers if tag in {"R", "S"}]
    unassigned = [(i, tag) for i, tag in centers if tag not in {"R", "S"}]
    bonds = []
    for b in mol.GetBonds():
        st = b.GetStereo()
        if st not in {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}:
            bonds.append(f"{b.GetBeginAtomIdx()}-{b.GetEndAtomIdx()}:{st}")
    return {
        "n_assigned_centers": len(assigned),
        "n_unassigned_centers": len(unassigned),
        "assigned": ";".join(f"{i}:{tag}" for i, tag in assigned),
        "n_bond_stereo": len(bonds),
        "bond_stereo": ";".join(bonds),
    }


def cip_conflict(sdf_mol: Chem.Mol, ref: Chem.Mol) -> dict:
    matches = sdf_mol.GetSubstructMatches(ref, uniquify=True, useChirality=False)
    if not matches:
        return {"status": "NO_MATCH", "detail": "子结构没有匹配"}
    ref_centers = {i: tag for i, tag in Chem.FindMolChiralCenters(ref, includeUnassigned=True)}
    sdf_centers = {i: tag for i, tag in Chem.FindMolChiralCenters(sdf_mol, includeUnassigned=True)}
    specified = {i: tag for i, tag in ref_centers.items() if tag in {"R", "S"}}
    best = None
    for match in matches:
        bad = []
        ok = []
        for ref_i, tag in specified.items():
            sdf_i = match[ref_i]
            sdf_tag = sdf_centers.get(sdf_i, "?")
            item = f"登记原子{ref_i}:{tag}/SDF原子{sdf_i}:{sdf_tag}"
            if sdf_tag == tag:
                ok.append(item)
            else:
                bad.append(item)
        extra = []
        for sdf_i, sdf_tag in sdf_centers.items():
            if sdf_tag not in {"R", "S"}:
                continue
            if sdf_i not in match:
                continue
            ref_i = match.index(sdf_i)
            if ref_centers.get(ref_i) not in {"R", "S"}:
                extra.append(f"SDF原子{sdf_i}:{sdf_tag}/登记原子{ref_i}:未指定")
        cand = {"bad": bad, "ok": ok, "extra": extra}
        if best is None or (len(bad), -len(ok)) < (len(best["bad"]), -len(best["ok"])):
            best = cand
    if best["bad"]:
        return {"status": "CIP_CONFLICT", "detail": "；".join(best["bad"]), "extra": "；".join(best["extra"])}
    return {"status": "SPECIFIED_OK", "detail": "；".join(best["ok"]) or "登记没有已赋值手性中心", "extra": "；".join(best["extra"])}


def inchi_layers(text: str) -> dict[str, str]:
    out = {}
    for part in text.split("/")[1:]:
        if part:
            out[part[0]] = part[1:]
    return out


def parse_t_layer(body: str) -> dict[int, str]:
    out = {}
    if not body:
        return out
    for tok in body.split(","):
        num = "".join(ch for ch in tok if ch.isdigit())
        if not num:
            continue
        out[int(num)] = tok[len(num):]
    return out


def classify_one(reg: dict) -> dict:
    eid = reg["global_ligand_entity_id"]
    sdf = heavy_from_3d(SDF_DIR / f"{eid}.sdf")
    iso = Chem.MolFromSmiles(reg["isomeric_smiles"])
    can = Chem.MolFromSmiles(reg["canonical_smiles"])
    row = {
        "global_ligand_entity_id": eid,
        "canonical_ligand_id": reg["canonical_ligand_id"],
        "pairs": reg["pairs"],
        "specified_stereo": reg["specified_stereo"],
        "class": "",
        "evidence": "",
        "specified_center_status": "",
        "specified_center_detail": "",
        "extra_sdf_stereo": "",
        "sdf_inchi_t": "",
        "registry_inchi_t": "",
        "sdf_inchi_b": "",
        "registry_inchi_b": "",
        "inchi_m_normalized": "",
    }
    if sdf is None or iso is None or can is None:
        row["class"] = "其他结构差异"
        row["evidence"] = "SDF 或登记 SMILES 无法解析"
        return row
    k_sdf = Chem.MolToInchiKey(sdf)
    k_iso = Chem.MolToInchiKey(iso)
    k_can = Chem.MolToInchiKey(can)
    sdf_flat = Chem.Mol(sdf)
    iso_flat = Chem.Mol(iso)
    Chem.RemoveStereochemistry(sdf_flat)
    Chem.RemoveStereochemistry(iso_flat)
    sdf_flat.RemoveAllConformers()
    flat_same = Chem.MolToInchiKey(sdf_flat) == Chem.MolToInchiKey(iso_flat)
    conn_same = k_sdf.split("-")[0] == k_can.split("-")[0]
    sf = stereo_flags(sdf)
    rf = stereo_flags(iso)
    row.update({
        "sdf_inchikey": k_sdf,
        "isomeric_inchikey": k_iso,
        "canonical_inchikey": k_can,
        "connectivity_same": "YES" if conn_same else "NO",
        "nonstereo_inchikey_same": "YES" if flat_same else "NO",
        "sdf_assigned_centers": sf["n_assigned_centers"],
        "sdf_unassigned_centers": sf["n_unassigned_centers"],
        "sdf_bond_stereo": sf["n_bond_stereo"],
        "registry_assigned_centers": rf["n_assigned_centers"],
        "registry_bond_stereo": rf["n_bond_stereo"],
        "isomeric_equals_canonical": "YES" if reg["isomeric_smiles"] == reg["canonical_smiles"] else "NO",
    })
    tagged = Chem.Mol(sdf)
    tagged.RemoveAllConformers()
    sdf_inchi = Chem.MolToInchi(tagged)
    iso_inchi = Chem.MolToInchi(iso)
    ls = inchi_layers(sdf_inchi)
    li = inchi_layers(iso_inchi)
    row["sdf_inchi_t"] = ls.get("t", "")
    row["registry_inchi_t"] = li.get("t", "")
    row["sdf_inchi_b"] = ls.get("b", "")
    row["registry_inchi_b"] = li.get("b", "")
    if not conn_same or not flat_same or ls.get("c") != li.get("c"):
        row["class"] = "其他结构差异"
        row["evidence"] = "连接层或去掉立体后的 InChIKey 与登记 isomeric SMILES 不同"
        return row
    if k_sdf == k_iso and ls.get("t", "") == li.get("t", "") and ls.get("b", "") == li.get("b", ""):
        row["class"] = "完整标识一致"
        row["evidence"] = "SDF 与 isomeric SMILES 的 InChI 立体层相同"
        return row
    reg_t = parse_t_layer(li.get("t", ""))
    sdf_t_raw = parse_t_layer(ls.get("t", ""))
    invert = li.get("m", "0") != ls.get("m", "0")
    sdf_t = {
        a: ({"+": "-", "-": "+"}.get(p, p) if invert else p)
        for a, p in sdf_t_raw.items()
    }
    specified = {a: p for a, p in reg_t.items() if p in {"+", "-"}}
    conflicts = [f"InChI原子{a}:登记{p}/SDF校正后{sdf_t.get(a, '缺失')}" for a, p in specified.items() if sdf_t.get(a) != p]
    extras = [f"InChI原子{a}:{p}" for a, p in sdf_t.items() if p in {"+", "-"} and reg_t.get(a) not in {"+", "-"}]
    row["inchi_m_normalized"] = "YES" if invert else "NO"
    bond_same = li.get("b", "") == ls.get("b", "")
    bond_note = ""
    if li.get("b") and not bond_same:
        conflicts.append(f"双键层登记/{li.get('b')} SDF/{ls.get('b', '缺失')}")
    elif ls.get("b") and not li.get("b"):
        bond_note = f"SDF 双键层 /{ls['b']} 在登记 InChI 中没有对应项"
    row["specified_center_status"] = "CIP_CONFLICT" if conflicts else "SPECIFIED_OK"
    row["specified_center_detail"] = "；".join(
        f"InChI原子{a}:登记{p}/SDF{sdf_t.get(a, p)}" for a, p in specified.items()
    )
    row["extra_sdf_stereo"] = "；".join(extras)
    if conflicts:
        chiral_ok = sdf.HasSubstructMatch(iso, useChirality=True)
        if chiral_ok:
            row["class"] = "已指定构型与 SDF 一致；未指定中心由三维坐标赋值"
            row["evidence"] = (
                "isomeric SMILES 的已写立体可以手性子结构匹配 SDF。"
                "InChI /m 不同造成校正后 parity 看起来相反，不能单独当成构型冲突。"
                "登记层 t=" + li.get("t", "") + "；SDF 层 t=" + ls.get("t", "") + "。"
            )
            return row
        row["class"] = "实际构型不一致"
        row["evidence"] = "登记 InChI 已指定的立体与 SDF 不同，且手性子结构匹配失败。" + "；".join(conflicts)
        return row
    if specified and (extras or bond_note):
        row["class"] = "已指定构型与 SDF 一致；未指定中心由三维坐标赋值"
        row["evidence"] = (
            "登记 InChI 中已用 +/− 写出的原子在 SDF 中保持同一parity。"
            + row["specified_center_detail"]
            + "。SDF 另外赋值："
            + (row["extra_sdf_stereo"] or bond_note)
        )
        return row
    if not conflicts and li.get("b") and bond_same and (extras or ls.get("t")):
        row["class"] = "登记双键构型与 SDF 一致；未指定手性中心由三维坐标赋值"
        row["evidence"] = (
            f"登记双键层 /{li['b']} 与 SDF 相同。"
            f"SDF 额外手性层 /t {ls.get('t','')} 对应登记里没有用 +/− 指定的原子。"
        )
        return row
    if not specified and not li.get("b") and (extras or ls.get("b") or ls.get("t")):
        row["class"] = "登记未指定立体化学；SDF 由三维坐标赋值"
        row["evidence"] = (
            f"登记指定立体中心数为 {reg['specified_stereo']}，登记 InChI 没有已确定的 /t 或 /b。"
            f"SDF 的立体层为 t={ls.get('t', '') or '无'} b={ls.get('b', '') or '无'}。"
            "这是三维嵌入给出的构型。"
        )
        return row
    if specified and not extras and bond_same:
        row["class"] = "其他结构差异"
        row["evidence"] = "已指定立体层一致，完整 InChI 仍有其他层不同。"
        return row
    row["class"] = "其他结构差异"
    row["evidence"] = f"立体层未能归入上述三类。登记 t={li.get('t','')} b={li.get('b','')}；SDF t={ls.get('t','')} b={ls.get('b','')}"
    return row


def coordinate_map(heavy: Chem.Mol, smiles: str, idx: dict[int, int], atoms: list[dict]) -> tuple[dict[int, int], str]:
    """Map PDBQT serials to SDF heavy atoms by remark SMILES, choosing the isomorphism with the smallest coordinate error."""
    tmpl = Chem.MolFromSmiles(smiles) if smiles else None
    if tmpl is None or not idx:
        return {}, "NO_REMARK_SMILES"
    matches = heavy.GetSubstructMatches(tmpl, uniquify=True)
    if not matches:
        return {}, "NO_ISOMORPHISM"
    by_serial = {a["serial"]: a for a in atoms}
    best = None
    for match in matches[:32]:
        serial_to_sdf = {}
        dists = []
        ok = True
        for smi_1, serial in idx.items():
            if smi_1 < 1 or smi_1 > len(match) or serial not in by_serial:
                ok = False
                break
            sdf_i = match[smi_1 - 1]
            serial_to_sdf[serial] = sdf_i
            a = by_serial[serial]
            p = heavy.GetConformer().GetAtomPosition(sdf_i)
            dists.append(((p.x - a["x"]) ** 2 + (p.y - a["y"]) ** 2 + (p.z - a["z"]) ** 2) ** 0.5)
        if not ok or not dists:
            continue
        score = max(dists)
        if best is None or score < best[0]:
            best = (score, serial_to_sdf)
    if best is None:
        return {}, "MAP_FAILED"
    return best[1], f"best_isomorphism_max_dev_{best[0]:.4f}A_nmatch_{min(len(matches), 32)}"


def macro_one(eid: str) -> dict:
    sdf_path = SDF_DIR / f"{eid}.sdf"
    pdbqt_path = PDBQT_DIR / f"{eid}.pdbqt"
    text = pdbqt_path.read_text(errors="replace")
    sdf = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=True)
    atoms = parse_pdbqt_atoms(text)
    smiles, idx = parse_smiles_idx(text)
    heavy = sdf_heavy(sdf) if sdf is not None else None
    serial_to_sdf, why = coordinate_map(heavy, smiles, idx, atoms) if heavy is not None else ({}, "no_sdf")
    real = [a for a in atoms if not a["is_glue"] and a["type"] not in {"H", "HD", "HS"}]
    glue = [a for a in atoms if a["is_glue"]]
    cgn = [a for a in atoms if a["is_real_macro"]]
    inv = {}
    dup = False
    for smi, serial in idx.items():
        if serial in inv:
            dup = True
        inv[serial] = smi
    missing = [a["serial"] for a in real if a["serial"] not in serial_to_sdf]
    extra_glue_in_idx = [a["serial"] for a in glue if a["serial"] in inv]
    elem_bad = []
    coord_bad = []
    for a in real:
        sdf_i = serial_to_sdf.get(a["serial"])
        if sdf_i is None or heavy is None:
            continue
        z = AD4_Z.get(a["type"])
        atom = heavy.GetAtomWithIdx(sdf_i)
        if z is None or z != atom.GetAtomicNum():
            elem_bad.append(f"{a['serial']}:{a['type']}->sdf{sdf_i}:Z{atom.GetAtomicNum()}")
        conf = heavy.GetConformer()
        p = conf.GetAtomPosition(sdf_i)
        dist = ((p.x - a["x"]) ** 2 + (p.y - a["y"]) ** 2 + (p.z - a["z"]) ** 2) ** 0.5
        if dist > 0.05:
            coord_bad.append(f"{a['serial']}:{dist:.3f}")
    pairs = complementary_cgn_pairs(atoms)
    tree = parse_torsion_tree(text)
    bond_ok = []
    tree_direct = []
    glue_on_partner = []
    for p in pairs:
        ia = serial_to_sdf.get(p["serial_a"])
        ib = serial_to_sdf.get(p["serial_b"])
        bond_ok.append(ia is not None and ib is not None and heavy is not None and sdf_has_bond(heavy, ia, ib))
        tree_direct.append(
            has_branch_edge(tree, p["serial_a"], p["serial_b"])
            or same_rigid_block(tree, p["serial_a"], p["serial_b"])
        )
        # each glue should sit on one of the two closure carbons
        cgn_coords = []
        for serial in (p["serial_a"], p["serial_b"]):
            atom = next(x for x in cgn if x["serial"] == serial)
            cgn_coords.append((atom["x"], atom["y"], atom["z"]))
        for g in glue:
            dists = [((g["x"] - c[0]) ** 2 + (g["y"] - c[1]) ** 2 + (g["z"] - c[2]) ** 2) ** 0.5 for c in cgn_coords]
            glue_on_partner.append(min(dists) < 0.01)
    smiles_mol = Chem.MolFromSmiles(smiles) if smiles else None
    smiles_conn = "NO"
    if smiles_mol is not None and heavy is not None:
        smiles_conn = "YES" if Chem.MolToInchiKey(smiles_mol).split("-")[0] == Chem.MolToInchiKey(heavy).split("-")[0] else "NO"
    problems = []
    if missing:
        problems.append("真实原子未映射")
    if extra_glue_in_idx or dup:
        problems.append("伪原子进入 SMILES IDX 或序号重复")
    if elem_bad:
        problems.append("原子类型与 SDF 元素不一致")
    if coord_bad:
        problems.append("映射原子坐标偏离 SDF")
    if not pairs or not all(bond_ok):
        problems.append("闭环两端在 SDF 中没有键")
    if smiles_conn != "YES":
        problems.append("PDBQT 备注 SMILES 与 SDF 连接图不同")
    if not glue_on_partner or not all(glue_on_partner):
        problems.append("G0 没有落在互补 CG0 坐标上")
    return {
        "global_ligand_entity_id": eid,
        "n_real_heavy": len(real),
        "n_cg0": len(cgn),
        "n_g0": len(glue),
        "map_status": why,
        "n_unmapped_real": len(missing),
        "n_glue_in_smiles_idx": len(extra_glue_in_idx),
        "n_element_mismatch": len(elem_bad),
        "element_mismatch": ";".join(elem_bad),
        "n_coord_mismatch_gt_0.05A": len(coord_bad),
        "coord_mismatch": ";".join(coord_bad[:8]),
        "sdf_closure_bond": "YES" if pairs and all(bond_ok) else "NO",
        "direct_branch_or_same_block": "YES" if any(tree_direct) else "NO",
        "glue_on_complementary_cg0": "YES" if glue_on_partner and all(glue_on_partner) else "NO",
        "remark_smiles_connectivity": smiles_conn,
        "closure_pairs": ";".join(
            f"{p['serial_a']}-{p['serial_b']}" for p in pairs
        ),
        "problem": "；".join(problems) if problems else "无原子或拓扑不一致",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prev = {r["global_ligand_entity_id"]: r for r in read_csv(PREV)}
    reg = read_csv(REG)
    mismatch_ids = {
        eid for eid, r in prev.items()
        if r["structure_identity"].startswith("连接图一致；")
    }
    rows = []
    for rec in reg:
        if rec["global_ligand_entity_id"] in mismatch_ids:
            rows.append(classify_one(rec))
    write_csv(OUT / "01_stereo_identity_167.csv", rows)
    print("stereo", Counter(r["class"] for r in rows))
    print("specified>0", Counter(r["class"] for r in rows if r["specified_stereo"] not in {"0", ""}))

    macro_ids = [
        eid for eid, r in prev.items() if r["macrocycle_meeko_type"] == "CG0+G0"
    ]
    macro_rows = [macro_one(eid) for eid in macro_ids]
    write_csv(OUT / "03_macrocycle_atom_map_29.csv", macro_rows)
    print("macro", len(macro_rows), Counter(r["problem"] for r in macro_rows))
    print("direct tree", Counter(r["direct_branch_or_same_block"] for r in macro_rows))
    print("element", sum(int(r["n_element_mismatch"]) for r in macro_rows))
    print("coord", sum(int(r["n_coord_mismatch_gt_0.05A"]) for r in macro_rows))


if __name__ == "__main__":
    main()
