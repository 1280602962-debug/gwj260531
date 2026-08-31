#!/usr/bin/env python3
"""Four-pair multi-seed Vina sensitivity (prespecified).

Frozen seeds: see protocol/SEED_LIST_FREEZE_v1.yaml
- 20260727: reuse existing primary scores (no redock)
- 20260811–20260814: redock with identical receptors/boxes/E/modes/energy_range

Ligand coordinates are prepared once with ETKDG seed 20260727 and reused.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "jcim_multiseed_v0"
TAB = OUT / "tables"
LIG = OUT / "ligands_pdbqt"
LOG = OUT / "logs"
POSE = OUT / "poses"
for d in (TAB, LIG, LOG, POSE):
    d.mkdir(parents=True, exist_ok=True)

VINA = "/home/gwj/miniconda3/bin/vina"
PY = "/home/gwj/miniconda3/bin/python"
MK_LIG = "/home/gwj/miniconda3/bin/mk_prepare_ligand.py"
PREP_SEED = 20260727
PRIMARY_SEED = 20260727
NEW_SEEDS = [20260811, 20260812, 20260813, 20260814]
ALL_SEEDS = [PRIMARY_SEED] + NEW_SEEDS
N_MODES = 9
ENERGY_RANGE = 3
WORKERS = 6
# Fail fast: do not let one hard ligand block a worker for long.
TIMEOUT = 480

PAIRS = {
    "EGFR/HER2": {
        "panel": ROOT / "data/egfr_her2_panel120_v0/tables/panel_v0_120.csv",
        "id": "panel_id",
        "class": "class",
        "smiles": "smiles",
        "E": 8,
        "primary_scores": ROOT / "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv",
        "primary_id": "ligand",
        "primary_A": "3POZ_affinity",
        "primary_B": "3RCD_affinity",
        "pockets": [
            {
                "name": "A",
                "pdb": "3POZ",
                "rec": ROOT / "data/egfr_her2_panel40_v0/receptors/3POZ_receptor.pdbqt",
                "box": ROOT / "data/egfr_her2_panel40_v0/boxes/3POZ_box.json",
            },
            {
                "name": "B",
                "pdb": "3RCD",
                "rec": ROOT / "data/egfr_her2_panel40_v0/receptors/3RCD_receptor.pdbqt",
                "box": ROOT / "data/egfr_her2_panel40_v0/boxes/3RCD_box.json",
            },
        ],
    },
    "AChE/BChE": {
        "panel": ROOT / "data/ache_bche_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "id": "panel_id",
        "class": "class",
        "smiles": "smiles",
        "E": 8,
        "primary_scores": ROOT / "data/ache_bche_panel_v0/tables/scores_vina_long.csv",
        "primary_long": True,
        "pockets": [
            {
                "name": "A",
                "pdb": "4EY7",
                "target_key": "ACHE",
                "rec": ROOT / "data/ache_bche_panel_v0/receptors/ACHE_receptor.pdbqt",
                "box": ROOT / "data/ache_bche_panel_v0/boxes/ACHE_box.json",
            },
            {
                "name": "B",
                "pdb": "4BDS",
                "target_key": "BCHE",
                "rec": ROOT / "data/ache_bche_panel_v0/receptors/BCHE_receptor.pdbqt",
                "box": ROOT / "data/ache_bche_panel_v0/boxes/BCHE_box.json",
            },
        ],
    },
    "PIK3CA/PIK3CB": {
        "panel": ROOT / "data/pik3ca_pik3cb_panel_v0/tables/panel_v0_strict_with_smiles.csv",
        "id": "panel_id",
        "class": "class",
        "smiles": "smiles",
        "E": 8,
        "primary_scores": ROOT / "data/pik3ca_pik3cb_panel_v0/tables/ablation_ligand_scores.csv",
        "primary_id": "ligand",
        "primary_A": "vina_PIK3CA",
        "primary_B": "vina_PIK3CB",
        "pockets": [
            {
                "name": "A",
                "pdb": "4L23",
                "rec": ROOT / "data/pik3ca_pik3cb_panel_v0/receptors/PIK3CA_receptor.pdbqt",
                "box": ROOT / "data/pik3ca_pik3cb_panel_v0/boxes/PIK3CA_box.json",
            },
            {
                "name": "B",
                "pdb": "2WXF",
                "rec": ROOT / "data/pik3ca_pik3cb_panel_v0/receptors/PIK3CB_receptor.pdbqt",
                "box": ROOT / "data/pik3ca_pik3cb_panel_v0/boxes/PIK3CB_box.json",
            },
        ],
    },
    "PIK3CA/mTOR": {
        "panel": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/panel_v0_48.csv",
        "id": "panel_id",
        "class": "class",
        "smiles": "smiles",
        "E": 16,
        "primary_scores": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv",
        "primary_id": "ligand",
        "primary_A": "4L23_affinity",
        "primary_B": "4JT6_affinity",
        "pockets": [
            {
                "name": "A",
                "pdb": "4L23",
                "rec": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4L23_receptor.pdbqt",
                "box": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4L23_box.json",
            },
            {
                "name": "B",
                "pdb": "4JT6",
                "rec": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/receptors/4JT6_receptor.pdbqt",
                "box": ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/boxes/4JT6_box.json",
            },
        ],
    },
}


def utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def prep_ligand(panel_id: str, smiles: str) -> tuple[Path | None, str]:
    pdbqt = LIG / f"{panel_id}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 0:
        return pdbqt, "exists"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "bad_smiles"
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = PREP_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        if AllChem.EmbedMolecule(mol, randomSeed=PREP_SEED) != 0:
            return None, "embed_fail"
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    sdf = OUT / "ligands_sdf" / f"{panel_id}.sdf"
    sdf.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(sdf))
    mol.SetProp("_Name", panel_id)
    w.write(mol)
    w.close()
    proc = subprocess.run(
        [PY, MK_LIG, "-i", str(sdf), "-o", str(pdbqt)],
        capture_output=True,
        text=True,
        timeout=300,
    )
    if proc.returncode != 0 or not pdbqt.exists():
        return None, f"meeko_fail:{proc.stderr[-200:]}"
    return pdbqt, "prepared"


def parse_mode1(out_pdbqt: Path) -> float | None:
    text = out_pdbqt.read_text(errors="ignore")
    m = re.search(r"REMARK VINA RESULT:\s+([-\d.]+)", text)
    return float(m.group(1)) if m else None


def _fail_row(job: dict, pocket: dict, status: str, reason: str) -> dict:
    return {
        "pair": job["pair"],
        "seed": job["seed"],
        "ligand": job["ligand"],
        "class": job["class"],
        "pocket": pocket["name"],
        "pdb": pocket["pdb"],
        "status": status,
        "vina_mode1": "",
        "reason": (reason or "")[:300],
    }


def dock_one(job: dict) -> dict:
    """Never raise: failed ligands are recorded and skipped so the pool keeps moving."""
    pair = job["pair"]
    seed = job["seed"]
    lig = job["ligand"]
    pocket = job["pocket"]
    try:
        pdbqt = Path(job["ligand_pdbqt"])
        if not pdbqt.exists() or pdbqt.stat().st_size <= 0:
            return _fail_row(job, pocket, "skip_missing_ligand", "ligand_pdbqt_missing")
        E = job["E"]
        box = json.loads(Path(pocket["box"]).read_text())
        out_dir = POSE / pair.replace("/", "_") / f"seed_{seed}" / pocket["pdb"] / lig
        out_dir.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_dir / "out.pdbqt"
        conf = out_dir / "vina.conf"
        log = out_dir / "vina.log"
        skip_marker = out_dir / "SKIPPED.fail"
        if skip_marker.exists():
            return _fail_row(job, pocket, "skip_cached_fail", skip_marker.read_text()[:200])
        if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            score = parse_mode1(out_pdbqt)
            return {
                "pair": pair,
                "seed": seed,
                "ligand": lig,
                "class": job["class"],
                "pocket": pocket["name"],
                "pdb": pocket["pdb"],
                "status": "exists" if score is not None else "parse_fail",
                "vina_mode1": score if score is not None else "",
                "reason": "",
            }
        conf.write_text(
            "\n".join(
                [
                    f"receptor = {pocket['rec']}",
                    f"ligand = {pdbqt}",
                    f"center_x = {box['center_x']}",
                    f"center_y = {box['center_y']}",
                    f"center_z = {box['center_z']}",
                    f"size_x = {box['size_x']}",
                    f"size_y = {box['size_y']}",
                    f"size_z = {box['size_z']}",
                    f"exhaustiveness = {E}",
                    f"num_modes = {N_MODES}",
                    f"energy_range = {ENERGY_RANGE}",
                    f"seed = {seed}",
                    f"cpu = 1",
                ]
            )
            + "\n"
        )
        t0 = time.time()
        try:
            proc = subprocess.run(
                [VINA, "--config", str(conf), "--out", str(out_pdbqt)],
                capture_output=True,
                text=True,
                timeout=TIMEOUT,
            )
        except subprocess.TimeoutExpired:
            skip_marker.write_text(f"timeout>{TIMEOUT}s")
            # remove partial out if any
            if out_pdbqt.exists() and out_pdbqt.stat().st_size == 0:
                out_pdbqt.unlink(missing_ok=True)
            return _fail_row(job, pocket, "timeout_skipped", f">{TIMEOUT}")

        log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
        if proc.returncode != 0 or not out_pdbqt.exists():
            reason = (proc.stderr or proc.stdout or "vina_fail")[-300:]
            skip_marker.write_text(reason)
            return _fail_row(job, pocket, "vina_fail_skipped", reason)
        score = parse_mode1(out_pdbqt)
        if score is None:
            skip_marker.write_text("parse_fail")
            return _fail_row(job, pocket, "parse_fail_skipped", "no_REMARK_VINA_RESULT")
        return {
            "pair": pair,
            "seed": seed,
            "ligand": lig,
            "class": job["class"],
            "pocket": pocket["name"],
            "pdb": pocket["pdb"],
            "status": "ok",
            "vina_mode1": score,
            "reason": f"elapsed_s={time.time()-t0:.1f}",
        }
    except Exception as exc:  # noqa: BLE001 — never abort the panel
        try:
            out_dir = POSE / pair.replace("/", "_") / f"seed_{seed}" / pocket["pdb"] / lig
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "SKIPPED.fail").write_text(repr(exc)[:300])
        except Exception:
            pass
        return _fail_row(job, pocket, "exception_skipped", repr(exc))


def load_primary_scores(cfg: dict) -> dict[str, dict[str, float]]:
    """ligand -> {A: score, B: score} using Vina affinity (more negative = better)."""
    out: dict[str, dict[str, float]] = {}
    path = cfg["primary_scores"]
    if cfg.get("primary_long"):
        for r in read_csv(path):
            if str(r.get("vina_mode", "1")) not in ("1", "01", "mode_01"):
                continue
            lig = r["ligand"]
            tgt = r["target"]
            score = float(r["vina_score"])
            out.setdefault(lig, {})
            for p in cfg["pockets"]:
                if tgt == p.get("target_key") or tgt == p["pdb"]:
                    out[lig][p["name"]] = score
        return out
    for r in read_csv(path):
        lig = r[cfg["primary_id"]]
        a, b = r.get(cfg["primary_A"], ""), r.get(cfg["primary_B"], "")
        if a in ("", None) or b in ("", None):
            continue
        out[lig] = {"A": float(a), "B": float(b)}
    return out


def write_rows(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def build_jobs(pair_filter: set[str] | None, seed_filter: set[int] | None) -> list[dict]:
    jobs = []
    prep_rows = []
    for pair, cfg in PAIRS.items():
        if pair_filter and pair not in pair_filter:
            continue
        panel = read_csv(cfg["panel"])
        for row in panel:
            lig = row[cfg["id"]]
            smi = row[cfg["smiles"]]
            cls = row[cfg["class"]]
            pdbqt, reason = prep_ligand(lig, smi)
            prep_rows.append(
                {"pair": pair, "ligand": lig, "class": cls, "status": "ok" if pdbqt else "prep_fail", "reason": reason}
            )
            if pdbqt is None:
                continue
            for seed in NEW_SEEDS:
                if seed_filter and seed not in seed_filter:
                    continue
                for pocket in cfg["pockets"]:
                    jobs.append(
                        {
                            "pair": pair,
                            "seed": seed,
                            "ligand": lig,
                            "class": cls,
                            "ligand_pdbqt": str(pdbqt),
                            "pocket": pocket,
                            "E": cfg["E"],
                        }
                    )
    write_rows(
        TAB / "multiseed_ligand_prep_v1.csv",
        prep_rows,
        ["pair", "ligand", "class", "status", "reason"],
    )
    return jobs


def export_primary_seed_rows() -> list[dict]:
    rows = []
    for pair, cfg in PAIRS.items():
        panel = {r[cfg["id"]]: r for r in read_csv(cfg["panel"])}
        scores = load_primary_scores(cfg)
        for lig, sc in scores.items():
            if lig not in panel:
                continue
            cls = panel[lig][cfg["class"]]
            for pocket in cfg["pockets"]:
                if pocket["name"] not in sc:
                    continue
                rows.append(
                    {
                        "pair": pair,
                        "seed": PRIMARY_SEED,
                        "ligand": lig,
                        "class": cls,
                        "pocket": pocket["name"],
                        "pdb": pocket["pdb"],
                        "status": "primary_reused",
                        "vina_mode1": sc[pocket["name"]],
                        "reason": "from_production_tables",
                    }
                )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", nargs="*", default=None, help="Subset of pairs")
    ap.add_argument("--seeds", nargs="*", type=int, default=None, help="Subset of NEW seeds only")
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    pair_filter = set(args.pairs) if args.pairs else None
    seed_filter = set(args.seeds) if args.seeds else None

    meta = {
        "started_utc": utc(),
        "primary_seed": PRIMARY_SEED,
        "new_seeds": NEW_SEEDS,
        "all_seeds": ALL_SEEDS,
        "prep_seed": PREP_SEED,
        "workers": args.workers,
        "note": "primary seed reused from production tables; new seeds redocked",
    }
    (TAB / "multiseed_run_meta_v1.json").write_text(json.dumps(meta, indent=2) + "\n")

    primary_rows = export_primary_seed_rows()
    if pair_filter:
        primary_rows = [r for r in primary_rows if r["pair"] in pair_filter]
    write_rows(
        TAB / "multiseed_scores_long_partial_v1.csv",
        primary_rows,
        ["pair", "seed", "ligand", "class", "pocket", "pdb", "status", "vina_mode1", "reason"],
    )
    print(f"[{utc()}] wrote {len(primary_rows)} primary-seed rows", flush=True)

    jobs = build_jobs(pair_filter, seed_filter)
    print(f"[{utc()}] jobs to dock: {len(jobs)}", flush=True)
    if args.dry_run:
        print("dry-run exit")
        return

    results = list(primary_rows)
    done = 0
    n_fail = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(dock_one, j): j for j in jobs}
        for fut in as_completed(futs):
            try:
                row = fut.result()
            except Exception as exc:  # noqa: BLE001
                j = futs[fut]
                row = _fail_row(j, j["pocket"], "future_exception_skipped", repr(exc))
            results.append(row)
            done += 1
            if str(row.get("status", "")).endswith("skipped") or row.get("status") in {
                "vina_fail",
                "timeout",
                "parse_fail",
            }:
                n_fail += 1
            if done % 25 == 0 or done == len(jobs):
                write_rows(
                    TAB / "multiseed_scores_long_partial_v1.csv",
                    results,
                    ["pair", "seed", "ligand", "class", "pocket", "pdb", "status", "vina_mode1", "reason"],
                )
                print(
                    f"[{utc()}] docked {done}/{len(jobs)} (skipped_or_failed={n_fail})",
                    flush=True,
                )

    write_rows(
        TAB / "multiseed_scores_long_v1.csv",
        results,
        ["pair", "seed", "ligand", "class", "pocket", "pdb", "status", "vina_mode1", "reason"],
    )
    meta["finished_utc"] = utc()
    meta["n_dock_jobs"] = len(jobs)
    meta["n_score_rows"] = len(results)
    (TAB / "multiseed_run_meta_v1.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[{utc()}] DONE rows={len(results)}", flush=True)


if __name__ == "__main__":
    main()
