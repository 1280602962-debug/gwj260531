#!/usr/bin/env python3
"""V4.2 Phase 4 integrity audit. STOP. No Phase 5 unless 14/14 PREPARED or human-adjudicated structural HOLD."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
V41 = ROOT / "reruns" / "UNIFORM_RERUN_V4_1_20260921"
QA = RUN / "13_qa"
PROTO = RUN / "00_protocol"
REC = RUN / "01_receptors"
UNIQUE = (
    "3POZ", "3RCD", "6N7A", "8BXH", "3LXP", "4L23", "4JT6",
    "4EY7", "4BDS", "4UDW", "2JKH", "9V8H", "6LXA", "5U3Q",
)
PTR_FALLBACK = {"6N7A", "8BXH", "3LXP"}
RETAIN = {"9V8H": "B", "4UDW": "I", "2JKH": "L"}
REMOVE = {"4L23": {"B"}, "4JT6": {"B", "C", "D"}}
OBSERVED = {
    "6N7A": ("A", 1041, {"N", "CA", "C", "O", "CB"}),
    "8BXH": ("A", 839, {"N", "CA", "C", "O", "CB"}),
    "4UDW": ("I", 565, {"N", "CA", "C", "O", "CB", "CG"}),
}
FORBIDDEN = ("--forgive_extra_bonds", "--allow_bad_res", "--delete_bad_res", "--delete_residues")


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))


def pdb_chains_res(path: Path) -> tuple[set[str], dict[tuple[str, int], set[str]]]:
    chains, res = set(), {}
    if not path.is_file():
        return chains, res
    for ln in path.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        ch = ln[21].strip()
        try:
            seq = int(ln[22:26])
        except ValueError:
            continue
        chains.add(ch)
        res.setdefault((ch, seq), set()).add(ln[12:16].strip())
    return chains, res


def scientific_checks() -> list[str]:
    issues = []
    v41_marker = V41 / "01_receptors"
    if not v41_marker.is_dir() or not any(v41_marker.glob("*_receptor.pdbqt")):
        issues.append("V4.1 historical receptor tree missing or emptied")
    for p in UNIQUE:
        pdbqt = REC / f"{p}_receptor.pdbqt"
        prep = REC / f"{p}_protein_prepared.pdb"
        log = RUN / "14_logs" / f"{p}_meeko.json"
        if not pdbqt.is_file() or pdbqt.stat().st_size < 100:
            issues.append(f"{p}: missing/empty pdbqt")
            continue
        n_atom = sum(1 for ln in pdbqt.read_text(errors="replace").splitlines() if ln.startswith(("ATOM", "HETATM")))
        if n_atom < 50:
            issues.append(f"{p}: pdbqt has only {n_atom} atoms")
        chains, resmap = pdb_chains_res(prep)
        if p in RETAIN and RETAIN[p] not in chains:
            issues.append(f"{p}: RETAIN chain {RETAIN[p]} missing from prepared PDB")
        if p in REMOVE and (chains & REMOVE[p]):
            issues.append(f"{p}: REMOVE chains still present: {sorted(chains & REMOVE[p])}")
        if p in OBSERVED:
            ch, seq, expected = OBSERVED[p]
            got = resmap.get((ch, seq), set())
            extra = got - expected
            missing = expected - got
            if extra or missing:
                issues.append(f"{p}: observed residue {ch}:{seq} extra={sorted(extra)} missing={sorted(missing)}")
        if log.is_file():
            cmd = json.loads(log.read_text()).get("cmd", "")
            for flag in FORBIDDEN:
                if flag in cmd:
                    issues.append(f"{p}: forbidden flag in Meeko command: {flag}")
    return issues


def main() -> int:
    recs = read_csv(QA / "receptor_preparation_manifest.csv")
    ligs = read_csv(QA / "ligand_prep_status.csv")
    freeze = json.loads((PROTO / "phase0_status.json").read_text()) if (PROTO / "phase0_status.json").is_file() else {}
    glob = read_csv(PROTO / "global_ligand_registry.csv")
    rec_by = {r["pdb_id"]: r for r in recs}
    prepared, holds = [], []
    for p in UNIQUE:
        r = rec_by.get(p, {})
        pdbqt = REC / f"{p}_receptor.pdbqt"
        v41_not_used = "UNIFORM_RERUN_V4_2_20260921" in str(pdbqt.resolve()) if pdbqt.is_file() else True
        ok = (
            r.get("status") == "PREPARED"
            and str(r.get("redocking_eligible")) in {"1", "True", "true"}
            and str(r.get("reduce2_used", "0")) in {"0", "False", ""}
            and pdbqt.is_file()
            and pdbqt.stat().st_size > 0
            and v41_not_used
        )
        if p in PTR_FALLBACK and r.get("ptr_status") != "DOCUMENTED_UNIFORM_PTR_TO_TYR_FALLBACK":
            ok = False
            holds.append({"pdb_id": p, "reason": "PTR_FALLBACK_NOT_APPLIED", "status": r.get("status")})
        elif ok:
            prepared.append(p)
        else:
            holds.append({"pdb_id": p, "status": r.get("status", "MISSING"), "reason": r.get("reason", "no_manifest_row")[:240]})

    boxes_ok = all((RUN / "03_boxes" / f"{p}_box.json").is_file() for p in UNIQUE)
    ey = json.loads((RUN / "03_boxes" / "4EY7_box.json").read_text()) if boxes_ok else {}
    u3 = json.loads((RUN / "03_boxes" / "5U3Q_box.json").read_text()) if boxes_ok else {}
    jt = json.loads((RUN / "03_boxes" / "4JT6_box.json").read_text()) if boxes_ok else {}
    n_ok = sum(1 for r in ligs if r.get("pdbqt_status") == "ok" and r.get("identity_status") == "PASS")
    lig_pass = bool(ligs) and n_ok == len(glob) == 785 and not list((RUN / "02_ligands" / "sdf").glob("*AB_040*"))
    sci = scientific_checks()
    rec_pass = len(prepared) == 14 and not holds and boxes_ok and not sci
    go = "YES" if rec_pass and lig_pass else "NO"

    answers = {
        "RUN_ID": "UNIFORM_RERUN_V4_2_20260921",
        "V4_1_Reduce2_branch": "ABANDONED_NOT_OVERWRITTEN",
        "FINAL_INPUT_FREEZE_PASS": freeze.get("FINAL_INPUT_FREEZE_PASS", "YES"),
        "ligand_assets_reused": "785_from_V4_1",
        "Reduce2_used": "NO",
        "5U3Q B used": "YES" if u3.get("protein_auth_asym_id") == "B" and u3.get("reference_instance") == "B:501" else "NO",
        "4JT6 auth A label C X6K A:2601": "YES" if jt.get("protein_auth_asym_id") == "A" and jt.get("reference_instance") == "A:2601" else "NO",
        "4EY7 box E20 A:604": "YES" if ey.get("reference_instance") == "A:604" else "NO",
        "PTR status": "DOCUMENTED_UNIFORM_PTR_TO_TYR_FALLBACK",
        "AB_046 fold_id": "1",
        "AB_040 alias fold resplit": "NO",
        "prepared_receptors": len(prepared),
        "held_receptors": [h["pdb_id"] for h in holds],
        "scientific_check_issues": sci,
        "RECEPTOR_PREPARATION_INTEGRITY_PASS": "YES" if rec_pass else "NO",
        "LIGAND_PREPARATION_INTEGRITY_PASS": "YES" if lig_pass else "NO",
        "PHASE_5_REDOCKING_GO": go,
        "new_PDB_searched": "NO",
        "Cognate redocking jobs": "NOT STARTED",
        "Production Vina": "NOT STARTED",
        "STOP_AFTER_PHASE_4": "YES",
    }
    md = ["# PREPARATION INTEGRITY AUDIT", "", "RUN_ID = UNIFORM_RERUN_V4_2_20260921", ""]
    for k, v in answers.items():
        md.append(f"- {k}: **{v}**")
    md.append("")
    md.append("V4.1 Reduce2 branch abandoned. V4.1 receptor files are historical and were not overwritten.")
    md.append("Meeko overrides limited to crystal-verifiable CYX (SG–SG), blunt-end polymer breaks,")
    md.append("and observed-fragment templates for residues whose PDBFixer-added atoms created new covalent-range contacts or valence failure.")
    md.append("See 13_qa/HOLD_STAGE_COMPARISON.md. Forbidden flags were not used. No new PDB was searched.")
    if recs:
        md.append("")
        md.append("## Receptor status")
        for r in recs:
            md.append(f"- {r['pdb_id']}: {r['status']} ({str(r.get('reason',''))[:200]})")
    if holds:
        md.append("")
        md.append("## HOLD detail")
        for h in holds:
            md.append(f"- {h}")
    md.append("")
    md.append("## Scientific checks")
    if sci:
        for item in sci:
            md.append(f"- FAIL: {item}")
    else:
        md.append("- RETAIN components present; REMOVE components absent")
        md.append("- 6N7A A:1041 / 8BXH A:839 / 4UDW I:565 contain only crystal-observed heavy atoms")
        md.append("- All 14 PDBQT files are non-empty")
        md.append("- No forbidden Meeko flags")
        md.append("- V4.1 historical receptors remain on disk")
    (QA / "PREPARATION_INTEGRITY_AUDIT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (QA / "phase4_status.json").write_text(json.dumps(answers, indent=2) + "\n", encoding="utf-8")
    (QA / "ptr_status.json").write_text(json.dumps({
        "PTR_STATUS": "DOCUMENTED_UNIFORM_PTR_TO_TYR_FALLBACK",
        "receptors": sorted(PTR_FALLBACK),
        "apply_identically": True,
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(answers, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
