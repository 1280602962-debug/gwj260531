#!/usr/bin/env python3
"""Layer-3 cognate QC for Track B eight receptors (DOCKING_PLAN_V1).

Gate: best-of-9 heavy-atom RMSD < 2.0 Å at exhaustiveness 8;
      one fallback at exhaustiveness 16 if E=8 fails.
Does not start production docking.
"""
from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "local_track_b_v0"
PDB_CACHE = Path("/tmp/track_b_pdb")
SEED = 20260727
N_MODES = 9
ENERGY_RANGE = 3
PAD = 5.0
MIN_EDGE = 20.0
CPU = 8
VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
MK_REC = "/home/gwj/miniconda3/bin/mk_prepare_receptor.py"
MK_LIG = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
OBABEL = shutil.which("obabel") or "/home/gwj/miniconda3/envs/cadd_tools/bin/obabel"
if not Path(OBABEL).exists():
    raise SystemExit(f"obabel not found at {OBABEL}")

# Eight unique Track B receptors (CTSK/CTSS excluded).
SPECS = [
    {
        "pdb": "4UDW",
        "cognate": "N6L",
        "protein": "F2",
        "note": "thrombin S1 heavy chain H; keep H/L/I; all HIS set to HIE for meeko",
        "set_template": "H:57=HIE,H:71=HIE,H:91=HIE,H:119=HIE,H:230=HIE",
    },
    {"pdb": "2JKH", "cognate": "BI7", "protein": "F10", "note": "FXa S1/S4; keep A/L"},
    {"pdb": "6N7A", "cognate": "KEV", "protein": "JAK1", "note": "JH1 ATP; ASU has two copies, first cognate on A"},
    {"pdb": "3LXP", "cognate": "IZA", "protein": "TYK2", "note": "JH1 ATP not JH2"},
    {"pdb": "8BXH", "cognate": "C87", "protein": "JAK2", "note": "JH1 ATP momelotinib"},
    {
        "pdb": "9V8H",
        "cognate": "BRL",
        "protein": "PPARG",
        "note": "LBD ternary; keep PG08-NL chain B; remove BRL/water/PEG only",
        "keep_chains": {"A", "B"},
    },
    {"pdb": "6LXA", "cognate": "EPA", "protein": "PPARA", "note": "LBD"},
    {
        "pdb": "5U3Q",
        "cognate": "7UJ",
        "protein": "PPARD",
        "note": "LBD agonist 1; keep chain A only (B:233 has altloc B/C, no A)",
        "keep_chains": {"A"},
    },
]


def download_pdb(pdb_id: str) -> Path:
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    dst = PDB_CACHE / f"{pdb_id}.pdb"
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    subprocess.run(["curl", "-fsSL", url, "-o", str(dst)], check=True)
    return dst


def extract_protein(src: Path, dst: Path, keep_chains: set[str] | None = None) -> None:
    lines = []
    for ln in src.read_text().splitlines():
        if not ln.startswith(("ATOM", "TER")):
            continue
        if keep_chains is not None and ln.startswith("ATOM"):
            if ln[21:22] not in keep_chains:
                continue
        if keep_chains is not None and ln.startswith("TER") and len(ln) > 21:
            if ln[21:22] not in keep_chains:
                continue
        lines.append(ln)
    if not lines:
        raise RuntimeError(f"no protein atoms in {src}")
    dst.write_text("\n".join(lines) + "\nEND\n")


def extract_first_cognate(src: Path, resname: str, dst: Path) -> str:
    """Extract first cognate instance; keep altloc A (or blank), drop other altlocs."""
    resid = None
    chain = None
    lines = []
    for line in src.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != resname:
            continue
        alt = line[16:17]
        if alt not in (" ", "A"):
            continue
        rid = line[22:26]
        ch = line[21:22]
        if resid is None:
            resid = rid
            chain = ch
        if rid == resid and ch == chain:
            # Normalize altloc to blank for downstream tools.
            if len(line) >= 17 and line[16] == "A":
                line = line[:16] + " " + line[17:]
            lines.append(line)
    if not lines:
        # Fallback: first altloc letter present for that residue.
        resid = chain = None
        alt_keep = None
        for line in src.read_text().splitlines():
            if not line.startswith("HETATM"):
                continue
            if line[17:20].strip() != resname:
                continue
            rid = line[22:26]
            ch = line[21:22]
            alt = line[16:17]
            if resid is None:
                resid, chain, alt_keep = rid, ch, alt
            if rid == resid and ch == chain and alt == alt_keep:
                if len(line) >= 17 and alt not in (" ",):
                    line = line[:16] + " " + line[17:]
                lines.append(line)
    if not lines:
        raise RuntimeError(f"no cognate {resname} in {src}")
    dst.write_text("\n".join(lines) + "\nEND\n")
    return f"{chain}:{resid.strip()}"


def ligand_heavy_xyz(path: Path):
    xyz = []
    for line in path.read_text().splitlines():
        if not line.startswith(("HETATM", "ATOM")):
            continue
        name = line[12:16].strip()
        el = (line[76:78].strip() if len(line) >= 78 else "") or name[0]
        if name.upper().startswith("H") or el.upper() == "H":
            continue
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not xyz:
        raise RuntimeError(f"no heavy atoms in {path}")
    return xyz


def box_from_xyz(xyz):
    xs, ys, zs = zip(*xyz)
    return {
        "center_x": round((min(xs) + max(xs)) / 2, 3),
        "center_y": round((min(ys) + max(ys)) / 2, 3),
        "center_z": round((min(zs) + max(zs)) / 2, 3),
        "size_x": round(max(max(xs) - min(xs) + 2 * PAD, MIN_EDGE), 3),
        "size_y": round(max(max(ys) - min(ys) + 2 * PAD, MIN_EDGE), 3),
        "size_z": round(max(max(zs) - min(zs) + 2 * PAD, MIN_EDGE), 3),
        "n_heavy_atoms": len(xyz),
        "pad_A": PAD,
        "min_edge_A": MIN_EDGE,
    }


def prepare_receptor(protein_pdb: Path, out_pdbqt: Path, set_template: str | None = None) -> None:
    base = out_pdbqt.with_suffix("")
    cmd = [
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
    ]
    if set_template:
        cmd.extend(["--set_template", set_template])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    produced = base.with_suffix(".pdbqt")
    if not produced.exists():
        cands = list(base.parent.glob(base.name + "*.pdbqt"))
        if not cands:
            raise RuntimeError(proc.stderr[-800:] or proc.stdout[-800:] or "mk_prepare_receptor failed")
        produced = cands[0]
    if produced.resolve() != out_pdbqt.resolve():
        out_pdbqt.write_bytes(produced.read_bytes())


def prepare_ligand(lig_pdb: Path, lig_pdbqt: Path) -> str:
    sdf = lig_pdb.with_suffix(".sdf")
    p1 = subprocess.run([OBABEL, str(lig_pdb), "-O", str(sdf), "-h"], capture_output=True, text=True)
    if p1.returncode != 0 or not sdf.exists():
        raise RuntimeError(f"obabel pdb->sdf failed: {p1.stderr[-300:]}")
    p2 = subprocess.run([PY, MK_LIG, "-i", str(sdf), "-o", str(lig_pdbqt)], capture_output=True, text=True)
    if p2.returncode == 0 and lig_pdbqt.exists():
        return "meeko"
    p3 = subprocess.run([OBABEL, str(sdf), "-O", str(lig_pdbqt)], capture_output=True, text=True)
    if p3.returncode != 0 or not lig_pdbqt.exists():
        raise RuntimeError(f"ligand prep failed meeko={p2.stderr[-200:]} obabel={p3.stderr[-200:]}")
    return "obabel"


def vina_dock(rec: Path, lig: Path, box: dict, out_pdbqt: Path, log: Path, exhaustiveness: int) -> None:
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
                f"exhaustiveness = {exhaustiveness}",
                f"num_modes = {N_MODES}",
                f"energy_range = {ENERGY_RANGE}",
                f"cpu = {CPU}",
                f"seed = {SEED}",
                f"out = {out_pdbqt}",
            ]
        )
        + "\n"
    )
    proc = subprocess.run([VINA, "--config", str(conf)], capture_output=True, text=True)
    log.write_text(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(proc.stderr[-400:] or proc.stdout[-400:] or "vina failed")


def pdbqt_models_xyz(path: Path):
    text = path.read_text().splitlines()
    models, cur = [], None
    for line in text:
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            if cur is not None:
                models.append(cur)
            cur = None
        elif cur is not None and line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            if name.upper().startswith("H"):
                continue
            cur.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not models:
        xyz = []
        for line in text:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                if name.upper().startswith("H"):
                    continue
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        models = [xyz]
    return models


def rmsd_assign(ref, mob) -> float:
    R = np.asarray(ref, float)
    M = np.asarray(mob, float)
    n = min(len(R), len(M))
    if n == 0:
        return float("nan")
    R, M = R[:n], M[:n]
    cost = ((R[:, None, :] - M[None, :, :]) ** 2).sum(axis=2)
    ri, ci = linear_sum_assignment(cost)
    return float(math.sqrt(cost[ri, ci].sum() / len(ri)))


def score_modes(ref_xyz, out_pdbqt: Path):
    models = pdbqt_models_xyz(out_pdbqt)
    rows = []
    for i, mob in enumerate(models, 1):
        if not mob:
            continue
        rows.append((i, rmsd_assign(ref_xyz, mob)))
    if not rows:
        return None, None, []
    top1 = rows[0][1]
    best_mode, best = min(rows, key=lambda t: t[1])
    return top1, (best_mode, best), rows


def process_one(spec: dict) -> dict:
    pdb_id = spec["pdb"]
    ccd = spec["cognate"]
    rec_dir = OUT / "receptors"
    box_dir = OUT / "boxes"
    cog_dir = OUT / "cognates"
    qc_dir = OUT / "cognate_qc"
    for d in (rec_dir, box_dir, cog_dir, qc_dir, OUT / "tables", OUT / "logs"):
        d.mkdir(parents=True, exist_ok=True)

    src = download_pdb(pdb_id)
    protein = rec_dir / f"{pdb_id}_protein.pdb"
    extract_protein(src, protein, spec.get("keep_chains"))
    lig_pdb = cog_dir / f"{pdb_id}_{ccd}.pdb"
    instance = extract_first_cognate(src, ccd, lig_pdb)
    xyz = ligand_heavy_xyz(lig_pdb)
    box = box_from_xyz(xyz)
    box.update({"pdb": pdb_id, "cognate": ccd, "instance": instance, "protein": spec["protein"]})
    (box_dir / f"{pdb_id}_box.json").write_text(json.dumps(box, indent=2) + "\n")

    rec_pdbqt = rec_dir / f"{pdb_id}_receptor.pdbqt"
    print(f"  prepare receptor {pdb_id} ...", flush=True)
    prepare_receptor(protein, rec_pdbqt, spec.get("set_template"))
    lig_pdbqt = cog_dir / f"{pdb_id}_{ccd}.pdbqt"
    how = prepare_ligand(lig_pdb, lig_pdbqt)

    result = {
        "protein": spec["protein"],
        "pdb": pdb_id,
        "cognate": ccd,
        "instance": instance,
        "ligand_prep": how,
        "n_heavy_ref": len(xyz),
        "box_size_x": box["size_x"],
        "box_size_y": box["size_y"],
        "box_size_z": box["size_z"],
        "note": spec["note"],
    }

    used_e = None
    for exhaust in (8, 16):
        out_dock = qc_dir / f"{pdb_id}_{ccd}_out_E{exhaust}.pdbqt"
        log = qc_dir / f"{pdb_id}_{ccd}_vina_E{exhaust}.log"
        print(f"  vina {pdb_id}/{ccd} E={exhaust} ...", flush=True)
        vina_dock(rec_pdbqt, lig_pdbqt, box, out_dock, log, exhaust)
        top1, best, per_mode = score_modes(xyz, out_dock)
        (qc_dir / f"{pdb_id}_{ccd}_rmsd_E{exhaust}.csv").write_text(
            "mode,rmsd_A\n" + "\n".join(f"{m},{r:.4f}" for m, r in per_mode) + "\n"
        )
        result[f"top1_rmsd_E{exhaust}"] = None if top1 is None else round(top1, 3)
        result[f"best_mode_E{exhaust}"] = None if best is None else best[0]
        result[f"best_of_9_rmsd_E{exhaust}"] = None if best is None else round(best[1], 3)
        if best is not None and best[1] < 2.0:
            used_e = exhaust
            result["exhaustiveness_pass"] = exhaust
            result["top1_rmsd"] = result[f"top1_rmsd_E{exhaust}"]
            result["best_mode"] = best[0]
            result["best_of_9_rmsd"] = round(best[1], 3)
            result["pass_rmsd_lt_2"] = True
            result["status"] = "PASS"
            break
        if exhaust == 8:
            print(f"  E=8 failed gate (best={None if best is None else round(best[1],3)}); trying E=16", flush=True)

    if used_e is None:
        result["exhaustiveness_pass"] = ""
        result["top1_rmsd"] = result.get("top1_rmsd_E16")
        result["best_mode"] = result.get("best_mode_E16")
        result["best_of_9_rmsd"] = result.get("best_of_9_rmsd_E16")
        result["pass_rmsd_lt_2"] = False
        result["status"] = "FAIL_RMSD"
    return result


def main():
    import sys

    only = {x.upper() for x in sys.argv[1:]} if len(sys.argv) > 1 else None
    csv_path = OUT / "tables" / "layer3_cognate_rmsd_v1.csv"
    by_pdb: dict[str, dict] = {}
    if only and csv_path.exists():
        with csv_path.open() as f:
            for r in csv.DictReader(f):
                by_pdb[r["pdb"]] = r

    for spec in SPECS:
        if only and spec["pdb"] not in only:
            continue
        print(f"=== {spec['protein']} {spec['pdb']}/{spec['cognate']} ===", flush=True)
        try:
            row = process_one(spec)
        except Exception as e:
            row = {
                "protein": spec["protein"],
                "pdb": spec["pdb"],
                "cognate": spec["cognate"],
                "status": "ERROR",
                "error": str(e),
                "pass_rmsd_lt_2": False,
                "note": spec["note"],
            }
        by_pdb[spec["pdb"]] = row
        print(json.dumps(row, indent=2), flush=True)

    rows = [by_pdb[s["pdb"]] for s in SPECS if s["pdb"] in by_pdb]
    fields = [
        "protein",
        "pdb",
        "cognate",
        "instance",
        "ligand_prep",
        "n_heavy_ref",
        "box_size_x",
        "box_size_y",
        "box_size_z",
        "exhaustiveness_pass",
        "top1_rmsd",
        "best_mode",
        "best_of_9_rmsd",
        "top1_rmsd_E8",
        "best_of_9_rmsd_E8",
        "top1_rmsd_E16",
        "best_of_9_rmsd_E16",
        "pass_rmsd_lt_2",
        "status",
        "note",
        "error",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    n_pass = sum(1 for r in rows if r.get("status") == "PASS")
    md = [
        "# Layer-3 cognate QC — Track B eight receptors\n\n",
        f"Protocol: Vina 1.2.7, seed {SEED}, num_modes {N_MODES}, energy_range {ENERGY_RANGE}; "
        f"box = cognate heavy-atom AABB + {PAD:g} Å/axis (min edge {MIN_EDGE:g} Å).\n",
        "Gate: best-of-9 heavy-atom RMSD < 2.0 Å at E=8; fallback E=16 once.\n",
        "9V8H keeps PG08-NL peptide (chain B). Production docking not started.\n\n",
        f"**Result: {n_pass}/{len(rows)} PASS**\n\n",
        "| protein | PDB | cognate | E_pass | top-1 RMSD | best-of-9 RMSD | status |\n",
        "|---------|-----|---------|-------:|-----------:|---------------:|--------|\n",
    ]
    for r in rows:
        md.append(
            f"| {r.get('protein','')} | {r.get('pdb','')} | {r.get('cognate','')} | "
            f"{r.get('exhaustiveness_pass','')} | {r.get('top1_rmsd','')} | "
            f"{r.get('best_of_9_rmsd','')} | {r.get('status','')} |\n"
        )
    (OUT / "cognate_qc" / "LAYER3_COGNATE_QC.md").write_text("".join(md))
    (OUT / "cognate_qc" / "layer3_cognate_qc_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(f"\nWrote {csv_path}")
    print(f"PASS {n_pass}/{len(rows)}")
    if n_pass < len(rows):
        print("GATE FAILED for one or more receptors — do not start production docking.")
    else:
        print("All eight passed Layer-3. Await human go-ahead before production Vina.")


if __name__ == "__main__":
    main()
