#!/usr/bin/env python3
"""Identity/component audit, save 8 cognate crystals, V4.2 prep, cognate AABB boxes.

Does not touch the 14 primary receptors. Does not search new PDBs.
"""
from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_1.phase2_prepare_receptors import (  # noqa: E402
    AA,
    CIF,
    REMOVE_HET,
    box_from,
    coherent_altloc,
    parse_mmcif,
    write_pdb,
    write_xyz,
)
from rerun_v4_2.alt_receptor_config import ALTS, CIF_ALT, CIF_CACHE, RUN  # noqa: E402
from rerun_v4_2.phase2_prepare_receptors import (  # noqa: E402
    OBS_TEMPLATE_PATH,
    assert_v42,
    count_mod,
    detect_disulfides,
    detect_polymer_breaks,
    drop_fixer_added_oxt,
    drop_fixer_invented_unphysical,
    format_blunt_ends,
    format_set_template,
    parse_pdb_atoms,
    pdbfixer_heavy,
    run_meeko,
)

REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
QA = RUN / "13_qa"
LOG = RUN / "14_logs"
AUDIT = QA / "alt_identity_component_audit.csv"
PREP_MANIFEST = QA / "alt_receptor_preparation_manifest.csv"
BOX_MANIFEST = QA / "alt_box_manifest.csv"


def ensure_cif(pdb: str) -> Path:
    dest = CIF_CACHE / f"{pdb}.cif"
    src = CIF_ALT / f"{pdb}.cif"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_file():
        if not dest.is_file() or dest.stat().st_size != src.stat().st_size:
            shutil.copy2(src, dest)
        return dest
    if dest.is_file():
        return dest
    raise FileNotFoundError(f"missing CIF for {pdb}: {src}")


def extract_alt(pdb: str) -> tuple[Path, dict]:
    cfg = ALTS[pdb]
    box_def = cfg["box"]
    cif = ensure_cif(pdb)
    atoms = parse_mmcif(cif.read_text(errors="replace"))
    auths = set(cfg["protein_auth"])
    cognate = [
        a for a in atoms
        if a["comp"] == box_def["ccd"] and a["auth_asym"] == box_def["auth"]
        and a["auth_seq"] == box_def["res"] and a["elem"] != "H"
    ]
    cog_path = REC / f"{pdb}_cognate_crystal.json"
    assert_v42(cog_path)
    write_xyz(cog_path, cognate, {"pdb": pdb, **box_def, "source": "raw_holo_before_deletion"})
    _, atoms_alt = coherent_altloc(atoms, auths)
    remove_auth = set(cfg["remove_auth"])
    retain_auth = set(cfg["retain_auth"])
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
    found_auths = sorted({a["auth_asym"] for a in prot})
    other_ccd = [
        a for a in atoms
        if a["comp"] == box_def["ccd"] and a["elem"] != "H"
        and not (a["auth_asym"] == box_def["auth"] and a["auth_seq"] == box_def["res"])
    ]
    return raw, {
        "n_cognate_saved": len(cognate),
        "n_protein_raw": len(prot),
        "removed_auth": sorted(remove_auth),
        "retain_auth_kept": sorted(retain_auth),
        "protein_auth_kept": found_auths,
        "nonprimary_cognate_atoms_not_used_for_box": len(other_ccd),
        "cognate_deleted_from_receptor": True,
        "waters_deleted": True,
        "cif_path": str(cif),
    }


def audit_one(pdb: str) -> dict:
    cfg = ALTS[pdb]
    box_def = cfg["box"]
    cif = ensure_cif(pdb)
    atoms = parse_mmcif(cif.read_text(errors="replace"))
    inst = [
        (a["auth_asym"], a["auth_seq"], a.get("label_asym", ""), a["comp"])
        for a in atoms if a["comp"] == box_def["ccd"] and a["elem"] != "H"
    ]
    inst = sorted(set(inst))
    match = any(a == box_def["auth"] and s == box_def["res"] and c == box_def["ccd"] for a, s, _, c in inst)
    chains = sorted({a["auth_asym"] for a in atoms if a["comp"] in AA})
    rec = {
        "pdb_id": pdb,
        "alt_id": cfg["id"],
        "pair": cfg["pair"],
        "side": cfg["side"],
        "design": cfg["design"],
        "kept_primary": cfg["kept_primary"],
        "protein_auth_declared": "".join(cfg["protein_auth"]),
        "remove_auth": ",".join(cfg["remove_auth"]),
        "box_ccd": box_def["ccd"],
        "box_instance_declared": box_def["instance"],
        "box_instance_found": int(match),
        "all_ccd_instances": ";".join(f"{a}:{s}" for a, s, _, _ in inst),
        "polymer_auth_observed": ",".join(chains),
        "equiv_nonprimary": cfg.get("equiv_nonprimary", ""),
        "note": cfg.get("note", ""),
        "status": "OK" if match else "HOLD",
    }
    return rec


def prepare_one(pdb: str) -> dict:
    cfg = ALTS[pdb]
    rec = {
        "pdb_id": pdb,
        "alt_id": cfg["id"],
        "protein_auth": "".join(cfg["protein_auth"]),
        "box_instance": cfg["box"]["instance"],
        "reduce2_used": 0,
        "status": "",
        "redocking_eligible": 0,
        "reason": "",
        "meeko_rc": "",
        "ptr_status": "NA",
        "cyx_templates": "",
        "observed_templates": "",
        "blunt_ends": "",
        "primary_receptors_untouched": 1,
    }
    raw, meta = extract_alt(pdb)
    rec.update(meta)
    rec["fixer_reused_not_rerun"] = 0
    fixed_openmm = REC / f"{pdb}_protein_fixer_openmm.pdb"
    flog = pdbfixer_heavy(raw, fixed_openmm)
    (LOG / f"{pdb}_pdbfixer.json").write_text(json.dumps(flog, indent=2) + "\n", encoding="utf-8")
    cleaned = REC / f"{pdb}_protein_fixer.pdb"
    oxt = drop_fixer_added_oxt(raw, fixed_openmm, cleaned)
    (LOG / f"{pdb}_oxt_cleanup.json").write_text(json.dumps(oxt, indent=2) + "\n", encoding="utf-8")
    rec["fixer_added_OXT_dropped"] = ";".join(oxt["dropped_fixer_added_OXT_not_in_raw"])

    observed = REC / f"{pdb}_protein_fixer_observed.pdb"
    clash = drop_fixer_invented_unphysical(raw, cleaned, observed)
    (LOG / f"{pdb}_fixer_unphysical.json").write_text(json.dumps(clash, indent=2) + "\n", encoding="utf-8")
    rec["fixer_residues_reverted_to_raw"] = ";".join(v["res"] for v in clash["reverted"])
    rec["observed_templates"] = ",".join(f"{v['res']}={v['template']}" for v in clash["reverted"])

    work = observed
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
    rec["meeko_exists"] = int(meeko["exists"])
    if meeko["rc"] != 0 or not meeko["exists"]:
        rec.update(status="HOLD", redocking_eligible=0, reason=f"meeko_fail; {meeko['stderr_tail'][-240:]}")
    else:
        rec.update(status="PREPARED", redocking_eligible=1, reason="v4_2_meeko_ok")
    rec["ptr_stage_raw_n"] = count_mod(raw, "PTR")["n"]
    rec["n_cognate_saved"] = meta["n_cognate_saved"]
    return rec


def write_box(pdb: str) -> dict:
    cfg = ALTS[pdb]
    box_def = cfg["box"]
    data = json.loads((REC / f"{pdb}_cognate_crystal.json").read_text())
    coords = [(a["x"], a["y"], a["z"]) for a in data["atoms"] if a.get("elem") != "H"]
    box = box_from(coords)
    box.update({
        "pdb_id": pdb,
        "ccd": box_def["ccd"],
        "reference_instance": box_def["instance"],
        "auth": box_def["auth"],
        "label": box_def["label"],
        "res": box_def["res"],
        "source": "box_defining_cognate_only",
        "equiv_nonprimary_not_used": cfg.get("equiv_nonprimary", ""),
    })
    dest = BOX / f"{pdb}_box.json"
    assert_v42(dest)
    dest.write_text(json.dumps(box, indent=2) + "\n", encoding="utf-8")
    return box


def main() -> int:
    REC.mkdir(parents=True, exist_ok=True)
    BOX.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)

    audits = [audit_one(pdb) for pdb in ALTS]
    with AUDIT.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(audits[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(audits)
    if any(a["status"] != "OK" for a in audits):
        print("AUDIT_HOLD")
        for a in audits:
            print(a["pdb_id"], a["status"], a["all_ccd_instances"])
        return 1

    rows = []
    boxes = []
    for pdb in ALTS:
        print("PREP", pdb, flush=True)
        rec = prepare_one(pdb)
        rows.append(rec)
        print(pdb, rec["status"], rec.get("n_cognate_saved"), rec.get("n_protein_raw"), rec["reason"][:160], flush=True)
        if rec["status"] == "PREPARED":
            boxes.append(write_box(pdb))

    fields = sorted({k for r in rows for k in r})
    with PREP_MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    if boxes:
        with BOX_MANIFEST.open("w", encoding="utf-8", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=sorted({k for b in boxes for k in b}), extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows(boxes)
    n_ok = sum(1 for r in rows if r["status"] == "PREPARED")
    print(f"ALT_PREPARED {n_ok}/8 PRIMARY_UNTOUCHED=YES")
    return 0 if n_ok == 8 else 1


if __name__ == "__main__":
    raise SystemExit(main())
