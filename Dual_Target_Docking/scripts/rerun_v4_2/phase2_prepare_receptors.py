#!/usr/bin/env python3
"""V4.2 receptor preparation: frozen extract → PDBFixer heavy → PTR fallback → Meeko.

No Reduce2. No chem_data/geostd patches. Does not overwrite V4.1.
Does not start Phase 5. Does not rewrite boxes. Does not search new PDBs.
"""
from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from meeko.utils.rdkitutils import covalent_radius

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from analysis.uniform_protocol_config import (  # noqa: E402
    BOX_DEFINING_LIGAND,
    COMPONENT_POLICY,
    PRIMARY_PROTEIN_AUTH,
    UNIQUE_14,
)
from rerun_v4_1.phase2_prepare_receptors import (  # noqa: E402
    AA,
    CIF,
    MEEKO,
    PY,
    REMOVE_HET,
    coherent_altloc,
    parse_mmcif,
    write_pdb,
    write_xyz,
)

RUN_ID = "UNIFORM_RERUN_V4_2_20260921"
RUN = ROOT / "reruns" / RUN_ID
REC = RUN / "01_receptors"
QA = RUN / "13_qa"
LOG = RUN / "14_logs"
PTR_FALLBACK = {"6N7A", "8BXH", "3LXP"}
TYR_HEAVY = {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "OXT"}
PTR_PHOS = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O4P", "OHP"}
SG_MIN, SG_MAX = 1.80, 2.30
PEPTIDE_CN_BREAK = 2.00
MEEKO_BOND_ALLOWANCE = 1.2
ELEM_Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15, "SE": 34}
OBS_TEMPLATE_PATH = RUN / "00_protocol" / "observed_fragment_templates.json"
OBS_TEMPLATE_BY_ATOMS = {
    frozenset({"N", "CA", "C", "O", "CB"}): "OBS_BB_CB",
    frozenset({"N", "CA", "C", "O", "CB", "CG"}): "OBS_BB_CB_CG",
}


def assert_v42(path: Path) -> None:
    if RUN_ID not in str(path.resolve()):
        raise SystemExit(f"refusing write outside V4.2: {path}")


def parse_pdb_atoms(text: str) -> list[dict]:
    atoms = []
    for ln in text.splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        try:
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
            })
        except ValueError:
            continue
    return atoms


def count_mod(path: Path, resn: str) -> dict:
    res = set()
    p_atoms = 0
    if not path.is_file():
        return {"n": 0, "n_P": 0, "ids": ""}
    for a in parse_pdb_atoms(path.read_text(errors="replace")):
        if a["comp"] != resn:
            continue
        res.add((a["auth_asym"], a["auth_seq"]))
        if a["elem"] == "P" or a["name"] == "P":
            p_atoms += 1
    return {"n": len(res), "n_P": p_atoms, "ids": ";".join(f"{c}{s}" for c, s in sorted(res))}


def pdbfixer_heavy(src: Path, dest: Path) -> dict:
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=str(src))
    fixer.findMissingResidues()
    missing = {str(k): v for k, v in fixer.missingResidues.items()}
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstd = [str(x) for x in fixer.nonstandardResidues]
    fixer.findMissingAtoms()
    miss_at = {str(k): (len(v) if hasattr(v, "__len__") else str(v)) for k, v in fixer.missingAtoms.items()}
    fixer.addMissingAtoms()
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return {
        "missing_residues_found_not_built": missing,
        "nonstandard_found_not_replaced": nonstd,
        "missing_atoms": miss_at,
        "replaceNonstandardResidues": False,
        "addMissingHydrogens": False,
        "missing_loop_reconstruction": False,
    }


def drop_fixer_added_oxt(raw: Path, fixed: Path, out: Path) -> dict:
    raw_keys = {
        (a["auth_asym"], a["auth_seq"], a.get("ins", ""))
        for a in parse_pdb_atoms(raw.read_text(errors="replace"))
        if a["name"] == "OXT"
    }
    kept, dropped = [], []
    for a in parse_pdb_atoms(fixed.read_text(errors="replace")):
        key = (a["auth_asym"], a["auth_seq"], a.get("ins", ""))
        if a["name"] == "OXT" and key not in raw_keys:
            dropped.append(f"{a['comp']} {a['auth_asym']}:{a['auth_seq']}")
            continue
        kept.append(a)
    write_pdb(out, kept)
    return {"dropped_fixer_added_OXT_not_in_raw": dropped}


def _atom_key(a: dict) -> tuple:
    return (a["auth_asym"], a["auth_seq"], a.get("ins", ""), a["name"])


def _res_key(a: dict) -> tuple:
    return (a["auth_asym"], a["auth_seq"], a.get("ins", ""))


def _meeko_thresh(e1: str, e2: str) -> float:
    z1 = ELEM_Z.get((e1 or "C").upper(), 6)
    z2 = ELEM_Z.get((e2 or "C").upper(), 6)
    return MEEKO_BOND_ALLOWANCE * (covalent_radius.get(z1, 0.76) + covalent_radius.get(z2, 0.76))


def _residue_sanitize_error(atoms: list[dict]) -> str | None:
    if not atoms:
        return "empty_residue"
    block = "".join(
        (
            f"ATOM  {i + 1:5d} {a['name']:>4s} {a['comp']:>3s} {a['auth_asym']}"
            f"{a['auth_seq']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            f"{a.get('occ', 1.0):6.2f}  1.00          {a['elem']:>2s}\n"
        )
        for i, a in enumerate(atoms)
    ) + "END\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    if mol is None:
        return "MolFromPDBBlock_failed"
    try:
        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return str(exc)
    return None


def drop_fixer_invented_unphysical(raw: Path, fixed: Path, out: Path) -> dict:
    """Drop fixer-added atoms that create new covalent-range contacts or valence failure.

    Crystal residues are not deleted. Only software-invented atoms absent from frozen raw
    are removed, matching the existing OXT-not-in-raw rule. Coordinates are not moved.
    """
    raw_atoms = parse_pdb_atoms(raw.read_text(errors="replace"))
    fixed_atoms = parse_pdb_atoms(fixed.read_text(errors="replace"))
    raw_keys = {_atom_key(a) for a in raw_atoms}
    added = [a for a in fixed_atoms if _atom_key(a) not in raw_keys]
    added_by_res: dict[tuple, list[dict]] = defaultdict(list)
    atoms_by_res: dict[tuple, list[dict]] = defaultdict(list)
    for a in added:
        added_by_res[_res_key(a)].append(a)
    for a in fixed_atoms:
        atoms_by_res[_res_key(a)].append(a)

    revert: dict[tuple, dict] = {}
    for res, added_atoms in added_by_res.items():
        reasons = []
        for a in added_atoms:
            for b in fixed_atoms:
                if _res_key(b) == res:
                    continue
                d = math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))
                if d < _meeko_thresh(a["elem"], b["elem"]):
                    reasons.append(
                        f"new_interresidue_meeko_covalent {a['name']}-{b['auth_asym']}:{b['auth_seq']}:{b['name']}={d:.3f}A"
                    )
        err = _residue_sanitize_error(atoms_by_res[res])
        if err:
            reasons.append(f"intraresidue_valence:{err}")
        if reasons:
            observed = sorted({x["name"] for x in raw_atoms if _res_key(x) == res})
            template = OBS_TEMPLATE_BY_ATOMS.get(frozenset(observed), "")
            revert[res] = {
                "res": f"{res[0]}:{res[1]}{res[2]}",
                "comp": next((x["comp"] for x in atoms_by_res[res]), ""),
                "dropped_fixer_atoms": [x["name"] for x in added_atoms],
                "observed_heavy": observed,
                "template": template,
                "reasons": reasons,
            }

    kept = []
    for a in fixed_atoms:
        if _res_key(a) in revert and _atom_key(a) not in raw_keys:
            continue
        kept.append(a)
    write_pdb(out, kept)

    missing_template = [v for v in revert.values() if not v["template"]]
    if missing_template:
        raise RuntimeError(f"no observed-atom template for {missing_template}")
    return {
        "n_fixer_added_atoms_total": len(added),
        "n_residues_reverted_to_raw": len(revert),
        "reverted": list(revert.values()),
        "moved_atoms": False,
        "deleted_residues": False,
        "deleted_retain_components": False,
    }


def ptr_to_tyr(src: Path, dest: Path) -> dict:
    out, converted, removed = [], [], []
    for a in parse_pdb_atoms(src.read_text(errors="replace")):
        if a["comp"] != "PTR":
            out.append(a)
            continue
        if a["name"] in PTR_PHOS or a["elem"] == "P":
            removed.append(f"{a['auth_asym']}:{a['auth_seq']}:{a['name']}")
            continue
        if a["name"] not in TYR_HEAVY:
            removed.append(f"{a['auth_asym']}:{a['auth_seq']}:{a['name']}")
            continue
        b = dict(a)
        b["comp"] = "TYR"
        out.append(b)
        converted.append(f"{a['auth_asym']}:{a['auth_seq']}")
    write_pdb(dest, out)
    return {
        "transformation": "PTR_to_TYR_drop_phosphate_keep_TYR_heavy",
        "converted_residues": sorted(set(converted)),
        "removed_atoms": removed,
        "n_PTR_in": count_mod(src, "PTR"),
        "n_PTR_out": count_mod(dest, "PTR"),
    }


def residue_atom_index(atoms: list[dict], chain: str, seq: int, name: str) -> int | None:
    recs = [a for a in atoms if a["auth_asym"] == chain and a["auth_seq"] == seq]
    for i, a in enumerate(recs):
        if a["name"] == name:
            return i
    return None


def detect_disulfides(atoms: list[dict]) -> list[dict]:
    sgs = [a for a in atoms if a["comp"] in {"CYS", "CYX"} and a["name"] == "SG"]
    pairs = []
    used = set()
    for i, a in enumerate(sgs):
        key_a = (a["auth_asym"], a["auth_seq"])
        if key_a in used:
            continue
        for b in sgs[i + 1 :]:
            key_b = (b["auth_asym"], b["auth_seq"])
            if key_b in used:
                continue
            d = math.dist((a["x"], a["y"], a["z"]), (b["x"], b["y"], b["z"]))
            if SG_MIN <= d <= SG_MAX:
                pairs.append({
                    "a": f"{a['auth_asym']}:{a['auth_seq']}",
                    "b": f"{b['auth_asym']}:{b['auth_seq']}",
                    "sg_sg_A": round(d, 3),
                    "fact": "crystal_SG-SG_distance_in_disulfide_range",
                    "template": "CYX",
                })
                used.add(key_a)
                used.add(key_b)
                break
    return pairs


def detect_polymer_breaks(atoms: list[dict]) -> list[dict]:
    by_chain: dict[str, dict[int, dict[str, dict]]] = defaultdict(lambda: defaultdict(dict))
    for a in atoms:
        by_chain[a["auth_asym"]][a["auth_seq"]][a["name"]] = a
    breaks = []
    for chain, seqs in by_chain.items():
        ordered = sorted(seqs)
        if not ordered:
            continue
        # Natural chain termini are left to Meeko auto-blunt.
        # Only crystal-verifiable internal polymer breaks get explicit blunt-ends.
        for i in range(len(ordered) - 1):
            r1, r2 = ordered[i], ordered[i + 1]
            c = seqs[r1].get("C")
            n = seqs[r2].get("N")
            seq_gap = r2 - r1 > 1
            dist = None
            if c and n:
                dist = math.dist((c["x"], c["y"], c["z"]), (n["x"], n["y"], n["z"]))
            peptide_break = dist is None or dist > PEPTIDE_CN_BREAK
            if seq_gap or peptide_break:
                c_idx = residue_atom_index(atoms, chain, r1, "C")
                n_idx = residue_atom_index(atoms, chain, r2, "N")
                fact = []
                if seq_gap:
                    fact.append(f"auth_seq_gap_{r1}_to_{r2}")
                if dist is None:
                    fact.append("missing_C_or_N")
                elif dist > PEPTIDE_CN_BREAK:
                    fact.append(f"C-N_{dist:.2f}A")
                if c_idx is not None:
                    breaks.append({"res": f"{chain}:{r1}", "atom_index": c_idx, "atom": "C", "fact": "+".join(fact)})
                if n_idx is not None:
                    breaks.append({"res": f"{chain}:{r2}", "atom_index": n_idx, "atom": "N", "fact": "+".join(fact)})
    return breaks


def format_set_template(pairs: list[dict]) -> str:
    ids = sorted({p["a"] for p in pairs} | {p["b"] for p in pairs})
    if not ids:
        return ""
    return ",".join(f"{rid}=CYX" for rid in ids)


def format_blunt_ends(breaks: list[dict]) -> str:
    # group by atom_index so A:123,200=2,A:1=0
    by_idx: dict[int, list[str]] = defaultdict(list)
    seen = set()
    for b in breaks:
        key = (b["res"], b["atom_index"])
        if key in seen:
            continue
        seen.add(key)
        by_idx[int(b["atom_index"])].append(b["res"])
    parts = []
    for idx in sorted(by_idx):
        res = ",".join(by_idx[idx])
        # parse_cmdline expects chain:res on first then more resnums can share chain
        # safest: emit each residue with chain:res=idx
        for rid in by_idx[idx]:
            parts.append(f"{rid}={idx}")
    return ",".join(parts)


def run_meeko(
    pdb: Path,
    out_base: Path,
    templates: str,
    blunt: str,
    log_stem: str,
    add_templates: str = "",
) -> dict:
    pdbqt = out_base.with_suffix(".pdbqt")
    cmd = [str(PY), str(MEEKO), "--read_pdb", str(pdb), "-o", str(out_base), "-p", str(pdbqt)]
    if templates:
        cmd += ["--set_template", templates]
    if blunt:
        cmd += ["--blunt_ends", blunt]
    if add_templates:
        cmd += ["--add_templates", add_templates]
    # forbidden flags are never added: --allow_bad_res, -d/--delete_residues, --forgive_extra_bonds
    proc = subprocess.run(cmd, capture_output=True, text=True)
    stdout_p = LOG / f"{log_stem}_meeko.stdout.txt"
    stderr_p = LOG / f"{log_stem}_meeko.stderr.txt"
    assert_v42(stdout_p)
    stdout_p.write_text(proc.stdout or "", encoding="utf-8")
    stderr_p.write_text(proc.stderr or "", encoding="utf-8")
    return {
        "rc": proc.returncode,
        "exists": int(pdbqt.is_file() and pdbqt.stat().st_size > 0),
        "cmd": " ".join(cmd),
        "stdout_bytes": len(proc.stdout or ""),
        "stderr_bytes": len(proc.stderr or ""),
        "stdout_path": str(stdout_p.relative_to(RUN)),
        "stderr_path": str(stderr_p.relative_to(RUN)),
        "stderr_tail": (proc.stderr or "")[-1200:],
        "used_forgive_extra_bonds": False,
        "used_delete_bad_res": False,
        "used_allow_bad_res": False,
    }


def extract_protein(pdb: str) -> tuple[Path, dict]:
    cif = CIF / f"{pdb}.cif"
    atoms = parse_mmcif(cif.read_text(errors="replace"))
    auths = set(PRIMARY_PROTEIN_AUTH[pdb])
    box_def = BOX_DEFINING_LIGAND[pdb]
    cognate = [
        a for a in atoms
        if a["comp"] == box_def["ccd"] and a["auth_asym"] == box_def["auth"]
        and a["auth_seq"] == box_def["res"] and a["elem"] != "H"
    ]
    cog_path = REC / f"{pdb}_cognate_crystal.json"
    assert_v42(cog_path)
    write_xyz(cog_path, cognate, {"pdb": pdb, **box_def, "source": "raw_holo_before_deletion"})
    _, atoms_alt = coherent_altloc(atoms, auths)
    remove_auth = {item[0] for item in COMPONENT_POLICY.get(pdb, ()) if item[2] == "REMOVE"}
    retain_auth = {item[0] for item in COMPONENT_POLICY.get(pdb, ()) if item[2] == "RETAIN"}
    prot = []
    for a in atoms_alt:
        if a["auth_asym"] in remove_auth:
            continue
        if a["auth_asym"] not in auths:
            continue
        if a["elem"] == "H":
            continue
        if a["comp"] == box_def["ccd"]:
            continue
        if a["comp"] in REMOVE_HET and a["comp"] not in AA:
            continue
        if a["comp"] not in AA:
            continue
        prot.append(a)
    raw = REC / f"{pdb}_protein_raw.pdb"
    assert_v42(raw)
    write_pdb(raw, prot)
    return raw, {
        "n_cognate_saved": len(cognate),
        "n_protein_raw": len(prot),
        "removed_auth": sorted(remove_auth),
        "retain_auth_kept": sorted(retain_auth),
        "cognate_deleted_from_receptor": True,
        "waters_deleted": True,
    }


def prepare_one(pdb: str) -> dict:
    rec = {
        "pdb_id": pdb,
        "protein_auth": "".join(PRIMARY_PROTEIN_AUTH[pdb]),
        "box_instance": BOX_DEFINING_LIGAND[pdb]["instance"],
        "box_reused_not_rewritten": 1,
        "reduce2_used": 0,
        "status": "",
        "redocking_eligible": 0,
        "reason": "",
        "meeko_rc": "",
        "ptr_status": "NA",
        "cyx_templates": "",
        "observed_templates": "",
        "blunt_ends": "",
    }
    reuse_fixer = "--reuse-fixer" in sys.argv and (REC / f"{pdb}_protein_raw.pdb").is_file() and (REC / f"{pdb}_protein_fixer.pdb").is_file()
    if reuse_fixer:
        raw = REC / f"{pdb}_protein_raw.pdb"
        cleaned = REC / f"{pdb}_protein_fixer.pdb"
        rec.update({
            "n_cognate_saved": "reused",
            "n_protein_raw": len(parse_pdb_atoms(raw.read_text(errors="replace"))),
            "cognate_deleted_from_receptor": True,
            "waters_deleted": True,
            "fixer_reused_not_rerun": 1,
        })
        oxt_log = LOG / f"{pdb}_oxt_cleanup.json"
        oxt = json.loads(oxt_log.read_text()) if oxt_log.is_file() else {"dropped_fixer_added_OXT_not_in_raw": []}
        rec["fixer_added_OXT_dropped"] = ";".join(oxt.get("dropped_fixer_added_OXT_not_in_raw", []))
    else:
        raw, meta = extract_protein(pdb)
        rec.update(meta)
        rec["fixer_reused_not_rerun"] = 0
        stages_raw = count_mod(raw, "PTR")
        fixed_openmm = REC / f"{pdb}_protein_fixer_openmm.pdb"
        flog = pdbfixer_heavy(raw, fixed_openmm)
        (LOG / f"{pdb}_pdbfixer.json").write_text(json.dumps(flog, indent=2) + "\n", encoding="utf-8")
        cleaned = REC / f"{pdb}_protein_fixer.pdb"
        oxt = drop_fixer_added_oxt(raw, fixed_openmm, cleaned)
        (LOG / f"{pdb}_oxt_cleanup.json").write_text(json.dumps(oxt, indent=2) + "\n", encoding="utf-8")
        rec["fixer_added_OXT_dropped"] = ";".join(oxt["dropped_fixer_added_OXT_not_in_raw"])
    stages = {"raw": count_mod(raw, "PTR")}

    observed = REC / f"{pdb}_protein_fixer_observed.pdb"
    clash = drop_fixer_invented_unphysical(raw, cleaned, observed)
    (LOG / f"{pdb}_fixer_unphysical.json").write_text(json.dumps(clash, indent=2) + "\n", encoding="utf-8")
    rec["fixer_residues_reverted_to_raw"] = ";".join(v["res"] for v in clash["reverted"])
    rec["observed_templates"] = ",".join(f"{v['res']}={v['template']}" for v in clash["reverted"])
    stages["after_pdbfixer"] = count_mod(observed, "PTR")

    work = observed
    if pdb in PTR_FALLBACK:
        fb = REC / f"{pdb}_protein_ptr_fallback.pdb"
        tlog = ptr_to_tyr(work, fb)
        (LOG / f"{pdb}_ptr_fallback.json").write_text(json.dumps(tlog, indent=2) + "\n", encoding="utf-8")
        rec["ptr_status"] = "DOCUMENTED_UNIFORM_PTR_TO_TYR_FALLBACK"
        stages["after_ptr_fallback"] = count_mod(fb, "PTR")
        work = fb
    else:
        stages["after_ptr_fallback"] = stages["after_pdbfixer"]

    atoms = parse_pdb_atoms(work.read_text(errors="replace"))
    disulf = detect_disulfides(atoms)
    breaks = detect_polymer_breaks(atoms)
    templates = format_set_template(disulf)
    observed_tmpl = rec.get("observed_templates") or ""
    if templates and observed_tmpl:
        templates = f"{templates},{observed_tmpl}"
    elif observed_tmpl:
        templates = observed_tmpl
    blunt = format_blunt_ends(breaks)
    rec["cyx_templates"] = format_set_template(disulf)
    rec["blunt_ends"] = blunt
    add_tmpl = str(OBS_TEMPLATE_PATH) if observed_tmpl else ""
    (LOG / f"{pdb}_structure_facts.json").write_text(
        json.dumps({
            "disulfides": disulf,
            "polymer_breaks": breaks,
            "templates": templates,
            "observed_templates": observed_tmpl,
            "add_templates": add_tmpl,
            "blunt": blunt,
        }, indent=2) + "\n",
        encoding="utf-8",
    )

    prepared = REC / f"{pdb}_protein_prepared.pdb"
    shutil.copy2(work, prepared)
    out_base = REC / f"{pdb}_receptor"
    assert_v42(out_base)
    meeko = run_meeko(prepared, out_base, templates, blunt, pdb, add_templates=add_tmpl)
    rec["meeko_rc"] = meeko["rc"]
    (LOG / f"{pdb}_meeko.json").write_text(json.dumps(meeko, indent=2) + "\n", encoding="utf-8")
    stages["after_meeko"] = count_mod(out_base.with_suffix(".pdbqt"), "PTR")
    rec["ptr_stage_raw_n"] = stages["raw"]["n"]
    rec["ptr_stage_fixer_n"] = stages["after_pdbfixer"]["n"]
    rec["ptr_stage_fallback_n"] = stages["after_ptr_fallback"]["n"]
    rec["ptr_stage_meeko_n"] = stages["after_meeko"]["n"]
    rec["ptr_lost_at"] = (
        "meeko_software_loss_then_uniform_fallback"
        if pdb in PTR_FALLBACK and stages["raw"]["n"]
        else ("none" if stages["raw"]["n"] == 0 else "see_stage_table")
    )
    if meeko["rc"] != 0 or not meeko["exists"]:
        rec.update(status="HOLD", redocking_eligible=0, reason=f"meeko_fail; {meeko['stderr_tail'][-240:]}")
    else:
        rec.update(status="PREPARED", redocking_eligible=1, reason="v4_2_meeko_ok")
    rec["_stages"] = stages
    return rec


def main() -> int:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    REC.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    rows = []
    stage_rows = []
    by_pdb = {}
    with ProcessPoolExecutor(max_workers=6) as pool:
        futs = {pool.submit(prepare_one, pdb): pdb for pdb in UNIQUE_14}
        for fut in as_completed(futs):
            rec = fut.result()
            by_pdb[rec["pdb_id"]] = rec
            print(rec["pdb_id"], rec["status"], str(rec["reason"])[:180], flush=True)
    for pdb in UNIQUE_14:
        rec = by_pdb[pdb]
        stages = rec.pop("_stages", {})
        stage_rows.append({"pdb_id": pdb, **{f"{k}_{sk}": v.get(sk, "") for k, v in stages.items() for sk in ("n", "n_P", "ids") if isinstance(v, dict)}})
        rows.append(rec)
    fields = sorted({k for r in rows for k in r})
    with (QA / "receptor_preparation_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with (QA / "ptr_loss_stage.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=sorted({k for r in stage_rows for k in r}), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(stage_rows)
    n_ok = sum(1 for r in rows if r["status"] == "PREPARED")
    print(f"PREPARED {n_ok}/14 Reduce2=NOT_USED PHASE5_NOT_STARTED")
    return 0 if n_ok == 14 else 1


if __name__ == "__main__":
    raise SystemExit(main())
