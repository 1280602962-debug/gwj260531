#!/usr/bin/env python3
"""Prepare MCL1/Bcl-xL primary receptors and run LC6 pose-gold cognate gate.

Frozen receptors: 3WIY (MCL1) / 3WIZ (BCL2L1), cognate LC6 (Tanaka compound 10).
Chains from RCSB polymer entity 1 primary chain A (see mcl1_bclxl_receptor_freeze_v1.csv).

Protocol: Vina 1.2.x, seed 20260727, exhaustiveness 8, num_modes 9, energy_range 3.
Gate (JCIM_NO_WETLAB_DEEP_PLAN_V2): both ends best-of-top3 RMSD < 2.0 Å.
If gate fails, panel docking may still run as a predeclared applicability stress-test.
"""
from __future__ import annotations

import json
import math
import subprocess
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data/mcl1_bclxl_panel_v0"
PDB_CACHE = OUT / "cache" / "pdb"
REC = OUT / "receptors"
BOX = OUT / "boxes"
QC = OUT / "cognate_qc"
TAB = OUT / "tables"
ANALYSIS = OUT / "analysis"
for d in (PDB_CACHE, REC, BOX, QC, TAB, ANALYSIS):
    d.mkdir(parents=True, exist_ok=True)

VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
MK_REC = "/home/gwj/miniconda3/bin/mk_prepare_receptor.py"
MK_LIG = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
OBABEL = "/mnt/d/CADD paper exercise/gnina/conda_env/bin/obabel"
SEED = 20260727
EXHAUST = 8
N_MODES = 9
PAD = 5.0
MIN_EDGE = 22.0  # LC6 is large; keep box ≥ cognate extent + pad

SPECS = [
    {
        "target": "MCL1",
        "pdb": "3WIY",
        "chain": "A",
        "lig_resname": "LC6",
        "note": "primary; entity1 chain A; Tanaka compound 10",
    },
    {
        "target": "BCL2L1",
        "pdb": "3WIZ",
        "chain": "A",
        "lig_resname": "LC6",
        "note": "primary; entity1 chain A; Tanaka compound 10",
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


def extract_ligand(src: Path, chain: str, resname: str, dst: Path) -> str:
    resid = None
    lines = []
    for ln in src.read_text().splitlines():
        if not ln.startswith("HETATM"):
            continue
        if ln[21] != chain:
            continue
        if ln[17:20].strip() != resname:
            continue
        rid = ln[22:26]
        if resid is None:
            resid = rid
        if rid == resid:
            lines.append(ln)
    if not lines:
        raise RuntimeError(f"no {resname} on chain {chain} in {src.name}")
    dst.write_text("\n".join(lines) + "\nEND\n")
    return resid.strip()


def ligand_heavy(path: Path):
    """Return list of (element, xyz) for heavy atoms."""
    out = []
    for ln in path.read_text().splitlines():
        if not ln.startswith(("HETATM", "ATOM")):
            continue
        name = ln[12:16].strip()
        el = (ln[76:78].strip() if len(ln) >= 78 else "") or "".join(
            c for c in name if c.isalpha()
        )
        el = el[:1].upper() + el[1:].lower() if el else "C"
        if el.upper() == "H":
            continue
        # AD4 aromatic carbon marker
        if el.upper() == "A":
            el = "C"
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        out.append((el.upper(), xyz))
    return out


def box_from_heavy(heavy):
    xs = [p[1][0] for p in heavy]
    ys = [p[1][1] for p in heavy]
    zs = [p[1][2] for p in heavy]
    return {
        "center_x": round((min(xs) + max(xs)) / 2, 3),
        "center_y": round((min(ys) + max(ys)) / 2, 3),
        "center_z": round((min(zs) + max(zs)) / 2, 3),
        "size_x": round(max(max(xs) - min(xs) + 2 * PAD, MIN_EDGE), 3),
        "size_y": round(max(max(ys) - min(ys) + 2 * PAD, MIN_EDGE), 3),
        "size_z": round(max(max(zs) - min(zs) + 2 * PAD, MIN_EDGE), 3),
        "n_heavy_atoms": len(heavy),
        "pad": PAD,
        "min_edge": MIN_EDGE,
    }


def prepare_receptor(protein_pdb: Path, out_pdbqt: Path) -> None:
    base = out_pdbqt.with_suffix("")
    proc = subprocess.run(
        [
            PY,
            MK_REC,
            "--read_pdb",
            str(protein_pdb),
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
                f"receptor prep failed: {(proc.stderr or proc.stdout)[-800:]}"
            )
        produced = cands[0]
    out_pdbqt.write_bytes(produced.read_bytes())


def prepare_ligand(lig_pdb: Path, lig_pdbqt: Path) -> str:
    sdf = lig_pdb.with_suffix(".sdf")
    p1 = subprocess.run(
        [OBABEL, str(lig_pdb), "-O", str(sdf), "-h"], capture_output=True, text=True
    )
    if p1.returncode != 0 or not sdf.exists():
        raise RuntimeError(f"obabel pdb->sdf: {p1.stderr[-400:]}")
    p2 = subprocess.run(
        [PY, MK_LIG, "-i", str(sdf), "-o", str(lig_pdbqt)],
        capture_output=True,
        text=True,
    )
    if p2.returncode == 0 and lig_pdbqt.exists():
        return "meeko"
    p3 = subprocess.run(
        [OBABEL, str(sdf), "-O", str(lig_pdbqt)], capture_output=True, text=True
    )
    if p3.returncode != 0 or not lig_pdbqt.exists():
        raise RuntimeError(
            f"ligand prep failed meeko={p2.stderr[-200:]} obabel={p3.stderr[-200:]}"
        )
    return "obabel"


def vina_dock(rec, lig, box, out_pdbqt, log_path) -> None:
    conf = out_pdbqt.with_suffix(".txt")
    conf.write_text(
        "\n".join(
            [
                f"receptor = {rec}",
                f"ligand = {lig}",
                f"center_x = {box['center_x']}",
                f"center_y = {box['center_y']}",
                f"center_z = {box['center_z']}",
                f"size_x = {box['size_x']}",
                f"size_y = {box['size_y']}",
                f"size_z = {box['size_z']}",
                f"exhaustiveness = {EXHAUST}",
                f"num_modes = {N_MODES}",
                "energy_range = 3",
                "cpu = 2",
                f"seed = {SEED}",
                f"out = {out_pdbqt}",
            ]
        )
        + "\n"
    )
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log_path.write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(proc.stderr[-400:] or proc.stdout[-400:])


def pdbqt_models_heavy(path: Path):
    text = path.read_text().splitlines()
    models, cur = [], []
    for line in text:
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = []
        elif line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            el = (line[77:79].strip() if len(line) > 77 else "") or name[0]
            el = el.upper()
            if el in {"H", "HD", "HS"}:
                continue
            if el == "A":
                el = "C"
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            cur.append((el, xyz))
    if not models and cur:
        models = [cur]
    if not models:
        # single-model file without MODEL tags
        cur = []
        for line in text:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                el = (line[77:79].strip() if len(line) > 77 else "") or name[0]
                el = el.upper()
                if el in {"H", "HD", "HS"}:
                    continue
                if el == "A":
                    el = "C"
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                cur.append((el, xyz))
        models = [cur]
    return models


def hungarian_rmsd(ref, mob) -> float:
    """Legacy coordinate diagnostic; not topology-aware/symmetry-corrected RMSD.

    The global assignment can match same-element atoms across inequivalent graph
    positions.  Keep it only to reproduce the historical v1 table; a formal gate
    requires a graph-isomorphism-constrained v2 implementation.
    """
    by_el = {}
    for el, xyz in ref:
        by_el.setdefault(el, []).append(xyz)
    used = {el: 0 for el in by_el}
    pairs_r, pairs_m = [], []
    for el, xyz in mob:
        bucket = by_el.get(el)
        if not bucket:
            continue
        i = used[el]
        if i >= len(bucket):
            continue
        pairs_r.append(bucket[i])
        pairs_m.append(xyz)
        used[el] = i + 1
    if len(pairs_r) < 3:
        # fallback: ignore element, match by count
        R = np.asarray([p[1] for p in ref], float)
        M = np.asarray([p[1] for p in mob], float)
        n = min(len(R), len(M))
        R, M = R[:n], M[:n]
    else:
        R = np.asarray(pairs_r, float)
        M = np.asarray(pairs_m, float)
        n = len(R)
    cost = ((R[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    ri, ci = linear_sum_assignment(cost)
    return float(math.sqrt(cost[ri, ci].sum() / len(ri)))


def process(spec: dict) -> dict:
    pdb = spec["pdb"]
    print(f"download {pdb}", flush=True)
    raw = download_pdb(pdb)
    prot = QC / f"{pdb}_chain{spec['chain']}_protein.pdb"
    lig_pdb = QC / f"{pdb}_{spec['lig_resname']}_crystal.pdb"
    extract_chain_protein(raw, spec["chain"], prot)
    resid = extract_ligand(raw, spec["chain"], spec["lig_resname"], lig_pdb)
    heavy = ligand_heavy(lig_pdb)
    box = box_from_heavy(heavy)
    box.update(
        {
            "pdb": pdb,
            "target": spec["target"],
            "ligand": spec["lig_resname"],
            "chain": spec["chain"],
            "resid": resid,
            "exhaustiveness": EXHAUST,
            "seed": SEED,
            "num_modes": N_MODES,
        }
    )
    box_path = BOX / f"{pdb}_box.json"
    box_path.write_text(json.dumps(box, indent=2) + "\n")
    # also write target-named alias
    (BOX / f"{spec['target']}_box.json").write_text(json.dumps(box, indent=2) + "\n")

    rec_pdbqt = REC / f"{pdb}_receptor.pdbqt"
    print(f"prepare receptor {pdb}", flush=True)
    prepare_receptor(prot, rec_pdbqt)
    (REC / f"{spec['target']}_receptor.pdbqt").write_bytes(rec_pdbqt.read_bytes())

    lig_pdbqt = QC / f"{pdb}_{spec['lig_resname']}.pdbqt"
    print(f"prepare ligand {pdb}/{spec['lig_resname']}", flush=True)
    lig_method = prepare_ligand(lig_pdb, lig_pdbqt)

    out_pose = QC / f"{pdb}_cognate_out.pdbqt"
    log_path = QC / f"{pdb}_cognate_vina.log"
    print(f"cognate redock {pdb}", flush=True)
    vina_dock(rec_pdbqt, lig_pdbqt, box, out_pose, log_path)

    models = pdbqt_models_heavy(out_pose)
    rmsds = [hungarian_rmsd(heavy, m) for m in models]
    top1 = rmsds[0] if rmsds else float("nan")
    top3 = min(rmsds[:3]) if rmsds else float("nan")
    best9 = min(rmsds) if rmsds else float("nan")
    row = {
        "target": spec["target"],
        "pdb": pdb,
        "chain": spec["chain"],
        "ligand_resname": spec["lig_resname"],
        "ligand_resid": resid,
        "n_crystal_heavy": len(heavy),
        "n_modes": len(rmsds),
        "rmsd_top1": round(top1, 3),
        "rmsd_best_of_top3": round(top3, 3),
        "rmsd_best_of_9": round(best9, 3),
        "top1_pass_lt_2": int(top1 < 2.0),
        "best_top3_pass_lt_2": int(top3 < 2.0),
        "best9_pass_lt_2": int(best9 < 2.0),
        "ligand_prep": lig_method,
        "exhaustiveness": EXHAUST,
        "seed": SEED,
        "rmsd_method": "legacy_hungarian_coordinate_diagnostic_not_topology_aware",
        "note": spec["note"],
        "prepared_at_utc": utc(),
    }
    print(
        f"{pdb}: top1={top1:.3f} best3={top3:.3f} best9={best9:.3f} "
        f"gate_top3={row['best_top3_pass_lt_2']}",
        flush=True,
    )
    return row


def main():
    rows = [process(s) for s in SPECS]
    (TAB / "cognate_qc_lc6_v1.json").write_text(json.dumps(rows, indent=2) + "\n")
    # CSV
    keys = list(rows[0].keys())
    lines = [",".join(keys)]
    for r in rows:
        lines.append(",".join(str(r[k]) for k in keys))
    (TAB / "cognate_qc_lc6_v1.csv").write_text("\n".join(lines) + "\n")

    verdict = f"""# MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1

Updated: `{utc()}`

Primary receptors frozen before panel docking: **3WIY** (MCL1) / **3WIZ** (Bcl-xL), cognate **LC6** (Tanaka compound 10).
Protocol: Vina seed `{SEED}`, exhaustiveness `{EXHAUST}`, num_modes `{N_MODES}`, energy_range 3.
Coordinate diagnostic: Hungarian element-matched heavy-atom absolute displacement (AD4 `A`→C). This assignment is not molecular-graph-constrained and is not topology-aware symmetry-corrected RMSD.

| target | PDB | top1 Å | best-of-top3 Å | best-of-9 Å | top3 gate (<2Å) |
|--------|-----|-------:|---------------:|------------:|:---------------:|
"""
    for r in rows:
        verdict += (
            f"| {r['target']} | {r['pdb']} | {r['rmsd_top1']} | "
            f"{r['rmsd_best_of_top3']} | {r['rmsd_best_of_9']} | "
            f"{r['best_top3_pass_lt_2']} |\n"
        )
    verdict += (
        "\n**Formal gate: UNMET / not validly completed.** The legacy coordinate "
        "diagnostic cannot establish a pose-gold pass. A formal gate requires "
        "graph-isomorphism-constrained RMSD, physical-validity checks, interaction "
        "recovery, and the prespecified second seed. Independently, the 3WIZ point "
        "value is above 2.0 Å. Do **not** package this pair as standard "
        "screening-performance evidence; panel docking is an "
        "**applicability stress-test** only.\n"
    )
    role = "applicability_stress_test"
    verdict += f"\n`panel_role={role}`\n"
    (ANALYSIS / "MCL1_BCLXL_LC6_POSE_GOLD_GATE_V1.md").write_text(verdict)
    (TAB / "lc6_gate_summary_v1.json").write_text(
        json.dumps({"gate_ok": gate_ok, "panel_role": role, "rows": rows}, indent=2)
        + "\n"
    )
    print(verdict)
    print("GATE_OK" if gate_ok else "GATE_FAIL", flush=True)


if __name__ == "__main__":
    main()
