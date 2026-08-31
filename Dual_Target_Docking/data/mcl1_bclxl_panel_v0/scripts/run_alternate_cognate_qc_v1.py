#!/usr/bin/env python3
"""MCL1/Bcl-xL alternate-receptor cognate QC (6UDV / 3SP7) + failure sensitivity.

Preselected alternate receptors (before AUROC inspection):
  MCL1 6UDV cognate Q51; Bcl-xL 3SP7 cognate 03B.
Gate language matches primary stress-test: report top-1 / best-of-top3 / best-of-9
element-Hungarian RMSD as coordinate diagnostic only (not topology-aware gold).
Do not retune boxes after seeing AUROC.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/mcl1_bclxl_panel_v0"
PDB_CACHE = OUT / "cache" / "pdb"
REC = OUT / "receptors"
BOX = OUT / "boxes"
QC = OUT / "cognate_qc"
TAB = OUT / "tables"
AN = OUT / "analysis"
for d in (PDB_CACHE, REC, BOX, QC, TAB, AN):
    d.mkdir(parents=True, exist_ok=True)

VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
MK_REC = "/home/gwj/miniconda3/bin/mk_prepare_receptor.py"
MK_LIG = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
SEED = 20260727
EXHAUST = 8
N_MODES = 9
PAD = 5.0
MIN_EDGE = 22.0

ALTS = [
    {
        "target": "MCL1",
        "pdb": "6UDV",
        "chain": "A",
        "lig_resname": "Q51",
        "role": "alternate",
        "note": "preselected resolution/holo WT; cognate compound 3",
    },
    {
        "target": "BCL2L1",
        "pdb": "3SP7",
        "chain": "A",
        "lig_resname": "03B",
        "role": "alternate",
        "note": "preselected resolution/holo WT; cognate BM903",
    },
]


def utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def download_pdb(pdb: str) -> Path:
    dest = PDB_CACHE / f"{pdb}.pdb"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    url = f"https://files.rcsb.org/download/{pdb}.pdb"
    with urllib.request.urlopen(url, timeout=60) as r:
        dest.write_bytes(r.read())
    return dest


def extract_chain_protein(src: Path, chain: str, dst: Path) -> None:
    lines = []
    for ln in src.read_text().splitlines():
        if ln.startswith("ATOM") and ln[21] == chain:
            lines.append(ln)
        elif ln.startswith("TER") and len(ln) > 21 and ln[21] == chain:
            lines.append(ln)
    if not lines:
        raise RuntimeError(f"no ATOM for chain {chain} in {src}")
    dst.write_text("\n".join(lines) + "\nEND\n")


def extract_ligand(src: Path, chain: str, resname: str, dst: Path) -> None:
    lines = []
    for ln in src.read_text().splitlines():
        if not ln.startswith("HETATM"):
            continue
        if ln[21] != chain:
            continue
        if ln[17:20].strip() != resname:
            continue
        # skip waters etc already filtered by resname
        lines.append(ln)
    if not lines:
        raise RuntimeError(f"no ligand {resname} chain {chain}")
    dst.write_text("\n".join(lines) + "\nEND\n")


def pdb_hetatm_coords(path: Path) -> tuple[list[str], np.ndarray]:
    elems, coords = [], []
    for ln in path.read_text().splitlines():
        if not (ln.startswith("HETATM") or ln.startswith("ATOM")):
            continue
        el = ln[76:78].strip() or ln[12:14].strip()[0]
        if el.upper().startswith("H"):
            continue
        elems.append(el.upper())
        coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return elems, np.asarray(coords, float)


def pdbqt_heavy_coords(path: Path) -> tuple[list[str], np.ndarray]:
    elems, coords = [], []
    for ln in path.read_text().splitlines():
        if not ln.startswith("ATOM"):
            continue
        name = ln[12:16].strip()
        if name.startswith("H"):
            continue
        el = "".join([c for c in name if c.isalpha()])
        el = el[0].upper() + el[1:].lower() if el else "C"
        # AutoDock types at end
        parts = ln.split()
        adt = parts[-1] if parts else el
        el = adt[0].upper()
        if el == "H":
            continue
        elems.append(el)
        coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return elems, np.asarray(coords, float)


def hungarian_rmsd(ref_e, ref_xyz, mob_e, mob_xyz) -> float:
    # element-constrained assignment (legacy diagnostic; not topology-aware)
    n = len(ref_e)
    if n == 0 or len(mob_e) == 0:
        return float("nan")
    cost = np.full((n, len(mob_e)), 1e6)
    for i, e in enumerate(ref_e):
        for j, f in enumerate(mob_e):
            if e == f:
                d = ref_xyz[i] - mob_xyz[j]
                cost[i, j] = float(np.dot(d, d))
    ri, cj = linear_sum_assignment(cost)
    if not np.isfinite(cost[ri, cj]).all() or cost[ri, cj].max() >= 1e6:
        return float("nan")
    return float(math.sqrt(cost[ri, cj].mean()))


def box_from_ligand(xyz: np.ndarray) -> dict:
    mn = xyz.min(axis=0)
    mx = xyz.max(axis=0)
    center = 0.5 * (mn + mx)
    size = np.maximum(mx - mn + 2 * PAD, MIN_EDGE)
    return {
        "center_x": float(center[0]),
        "center_y": float(center[1]),
        "center_z": float(center[2]),
        "size_x": float(size[0]),
        "size_y": float(size[1]),
        "size_z": float(size[2]),
        "n_heavy_atoms": int(len(xyz)),
    }


def prepare_and_redock(spec: dict) -> dict:
    pdb = spec["pdb"]
    src = download_pdb(pdb)
    prot = QC / f"{pdb}_protein_chain{spec['chain']}.pdb"
    lig_pdb = QC / f"{pdb}_{spec['lig_resname']}_crystal.pdb"
    extract_chain_protein(src, spec["chain"], prot)
    extract_ligand(src, spec["chain"], spec["lig_resname"], lig_pdb)

    rec_pdbqt = REC / f"{pdb}_receptor.pdbqt"
    # Always rebuild: a prior bad call produced ligand-style ROOT/BRANCH PDBQT.
    base = rec_pdbqt.with_suffix("")
    proc = subprocess.run(
        [
            PY,
            MK_REC,
            "--read_pdb",
            str(prot),
            "-o",
            str(base),
            "-p",
            "-a",
            "--default_altloc",
            "A",
        ],
        capture_output=True,
        text=True,
    )
    produced = base.with_suffix(".pdbqt")
    if not produced.exists():
        cands = list(base.parent.glob(base.name + "*.pdbqt"))
        if not cands:
            raise RuntimeError(
                f"receptor prep failed for {pdb}: {(proc.stderr or proc.stdout)[-800:]}"
            )
        produced = cands[0]
    if produced.resolve() != rec_pdbqt.resolve():
        rec_pdbqt.write_bytes(produced.read_bytes())
    # Sanity: rigid receptor must not contain ROOT/BRANCH
    head = rec_pdbqt.read_text(errors="ignore")[:2000]
    if "ROOT" in head or "BRANCH" in head:
        raise RuntimeError(f"{rec_pdbqt} looks like a flexible ligand PDBQT")

    # Ligand: PDB -> SDF (+H) -> meeko; fallback obabel pdbqt
    lig_pdbqt = QC / f"{pdb}_{spec['lig_resname']}_crystal.pdbqt"
    obabel = "/mnt/d/CADD paper exercise/gnina/conda_env/bin/obabel"
    lig_sdf = QC / f"{pdb}_{spec['lig_resname']}_crystal.sdf"
    subprocess.run([obabel, str(lig_pdb), "-O", str(lig_sdf), "-h"], check=True, capture_output=True)
    p2 = subprocess.run(
        [PY, MK_LIG, "-i", str(lig_sdf), "-o", str(lig_pdbqt)],
        capture_output=True,
        text=True,
    )
    if p2.returncode != 0 or not lig_pdbqt.exists():
        subprocess.run([obabel, str(lig_sdf), "-O", str(lig_pdbqt)], check=True, capture_output=True)

    ref_e, ref_xyz = pdb_hetatm_coords(lig_pdb)
    box = box_from_ligand(ref_xyz)
    box_path = BOX / f"{pdb}_box.json"
    box_path.write_text(json.dumps(box, indent=2) + "\n")

    out_pdbqt = QC / f"{pdb}_cognate_out.pdbqt"
    conf = QC / f"{pdb}_cognate.conf"
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec_pdbqt}",
                f"ligand = {lig_pdbqt}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                f"energy_range = 3",
                f"seed = {SEED}",
            ]
        )
        + "\n"
    )
    log = QC / f"{pdb}_cognate_vina.log"
    proc = subprocess.run(
        [VINA, "--config", str(conf), "--out", str(out_pdbqt)],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    log.write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0:
        return {
            **spec,
            "status": "vina_fail",
            "stderr": proc.stderr[-400:],
        }

    # split modes roughly by MODEL
    text = out_pdbqt.read_text()
    models = text.split("MODEL")
    rmsds = []
    for block in models[1:]:
        tmp = QC / f"{pdb}_mode_tmp.pdbqt"
        tmp.write_text("MODEL" + block)
        mob_e, mob_xyz = pdbqt_heavy_coords(tmp)
        rmsds.append(hungarian_rmsd(ref_e, ref_xyz, mob_e, mob_xyz))
    rmsds = [r for r in rmsds if r == r]
    top1 = rmsds[0] if rmsds else float("nan")
    best3 = min(rmsds[:3]) if rmsds else float("nan")
    best9 = min(rmsds) if rmsds else float("nan")
    gate_top3 = bool(best3 < 2.0) if rmsds else False
    return {
        "target": spec["target"],
        "pdb": pdb,
        "role": spec["role"],
        "cognate": spec["lig_resname"],
        "chain": spec["chain"],
        "status": "ok",
        "rmsd_top1": round(top1, 3),
        "rmsd_best_top3": round(best3, 3),
        "rmsd_best9": round(best9, 3),
        "n_modes_scored": len(rmsds),
        "gate_best_top3_lt_2A": int(gate_top3),
        "matching": "element_hungarian_legacy_diagnostic",
        "topology_aware": 0,
        "note": spec["note"],
        "utc": utc(),
    }


def failure_sensitivity() -> None:
    panel = list(csv.DictReader((TAB / "mcl1_bclxl_chembl_panel96_v1.csv").open()))
    skips = list(csv.DictReader((TAB / "vina_skips_MBX_v1.csv").open()))
    scores = list(csv.DictReader((TAB / "vina_scores_MBX_v1.csv").open()))
    failed_ids = sorted({r["panel_id"] for r in skips})
    rows = []
    for pid in failed_ids:
        p = next(r for r in panel if r["panel_id"] == pid)
        smi = p.get("smiles") or p.get("canonical_smiles") or ""
        # try common smile columns
        for k in p:
            if "smile" in k.lower() and p[k]:
                smi = p[k]
                break
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is not None:
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        rows.append(
            {
                "panel_id": pid,
                "class": p.get("class", ""),
                "molecule_chembl_id": p.get("molecule_chembl_id", ""),
                "reason": next(r["reason"] for r in skips if r["panel_id"] == pid),
                "heavy_atoms": mol.GetNumHeavyAtoms() if mol else "",
                "mw": round(Descriptors.MolWt(mol), 2) if mol else "",
                "rotatable_bonds": Lipinski.NumRotatableBonds(mol) if mol else "",
                "clogp": round(Descriptors.MolLogP(mol), 3) if mol else "",
            }
        )
    with (TAB / "mcl1_bclxl_failure_properties_v1.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["panel_id"])
        w.writeheader()
        w.writerows(rows)

    # Rank-extreme bounds on directional AUROCs using complete-case + failed as best/worst
    # Build complete ligands from scores
    by = {}
    for r in scores:
        if r.get("status") not in ("ok", "success", "") and r.get("vina_mode1") in ("", None):
            continue
        # flexible schema
        pid = r.get("panel_id") or r.get("ligand")
        if not pid:
            continue
        by.setdefault(pid, {"class": r.get("class", "")})
        pocket = r.get("target") or r.get("pdb") or r.get("pocket")
        val = r.get("vina_mode1") or r.get("vina_score")
        if val in ("", None):
            continue
        if pocket in ("MCL1", "3WIY", "A"):
            by[pid]["A"] = -float(val)
        elif pocket in ("BCL2L1", "3WIZ", "B"):
            by[pid]["B"] = -float(val)

    # reload scores with known schema
    if scores:
        sample = scores[0]
        (AN / "mcl1_failure_score_schema_sample.json").write_text(json.dumps(sample, indent=2))

    # simpler: use analyze script table if present
    # Class distribution of failures
    from collections import Counter

    cls_counts = Counter(r["class"] for r in rows)
    with (TAB / "mcl1_bclxl_failure_class_counts_v1.csv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["class", "n_failed_ligands"])
        w.writeheader()
        for c, n in sorted(cls_counts.items()):
            w.writerow({"class": c, "n_failed_ligands": n})


def main():
    results = []
    for spec in ALTS:
        print(f"[{utc()}] cognate QC {spec['pdb']}...", flush=True)
        results.append(prepare_and_redock(spec))
        print(results[-1], flush=True)

    out = TAB / "mcl1_bclxl_alternate_cognate_qc_v1.csv"
    fields = [
        "target",
        "pdb",
        "role",
        "cognate",
        "chain",
        "status",
        "rmsd_top1",
        "rmsd_best_top3",
        "rmsd_best9",
        "n_modes_scored",
        "gate_best_top3_lt_2A",
        "matching",
        "topology_aware",
        "note",
        "utc",
    ]
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    failure_sensitivity()

    lines = [
        "# MCL1/Bcl-xL alternate-receptor cognate QC + failure sensitivity",
        "",
        "Role: applicability stress-test support. Not a fifth primary pair.",
        "",
        "## Alternate cognate QC (element-Hungarian diagnostic)",
        "",
        "| target | PDB | cognate | top-1 | best-top3 | best9 | gate_top3<2Å |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for r in results:
        if r.get("status") != "ok":
            lines.append(f"| {r.get('target')} | {r.get('pdb')} | {r.get('lig_resname', r.get('cognate'))} | FAIL | FAIL | FAIL | 0 |")
            continue
        lines.append(
            f"| {r['target']} | {r['pdb']} | {r['cognate']} | {r['rmsd_top1']} | {r['rmsd_best_top3']} | "
            f"{r['rmsd_best9']} | {r['gate_best_top3_lt_2A']} |"
        )
    sp7 = next((r for r in results if r.get("pdb") == "3SP7"), {})
    udv = next((r for r in results if r.get("pdb") == "6UDV"), {})
    lines += [
        "",
        "## Interpretation",
        "",
        f"- 6UDV gate_top3: {udv.get('gate_best_top3_lt_2A')}",
        f"- 3SP7 gate_top3: {sp7.get('gate_best_top3_lt_2A')}",
        "- Matching is legacy element-Hungarian (same class as primary LC6 diagnostic); not topology-aware.",
        "- If an alternate fails the <2 Å diagnostic, do **not** retune the box to rescue it.",
        "- Full-panel alternate-receptor AUROC is run only if the corresponding cognate diagnostic passes; otherwise stop at QC.",
        "",
        "## Failure sensitivity (primary 3WIY/3WIZ panel)",
        "",
        "See `tables/mcl1_bclxl_failure_properties_v1.csv` and `tables/mcl1_bclxl_failure_class_counts_v1.csv`.",
        "Failures were ligand-prep/embed failures, not Vina timeouts.",
        "",
    ]
    (AN / "MCL1_BCLXL_ALT_RECEPTOR_QC_V1.md").write_text("\n".join(lines) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
