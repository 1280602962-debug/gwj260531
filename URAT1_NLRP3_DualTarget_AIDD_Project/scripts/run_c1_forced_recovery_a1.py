#!/usr/bin/env python3
"""Forced-recovery URAT1 A1 gate table (8 known ligands).

Uses C1 free-dock settings (exhaustiveness 32, num_modes 9, CNNscore-selected).
Does NOT write new scientific locks — audit table only.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from parse_c1_sdf_readouts import (  # noqa: E402
    _fprop,
    carboxylate_oxygens,
    heavy_centroid,
    load_poses,
    load_ref_centroid,
    min_acid_arg_dist,
)

ARG_THRESH = 7.7027
CENTROID_MAX = 6.0
OUT = PROJECT_ROOT / "data/campaigns/c1/03_forced_recovery"
REFS = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"


def phenolate_or_acid_oxygens(mol: Chem.Mol) -> list[list[float]]:
    """Broader O set for audit only (carboxylate + phenolate/phenol)."""
    coords = list(carboxylate_oxygens(mol))
    pats = [
        Chem.MolFromSmarts("[OX1-]-[c,C]"),
        Chem.MolFromSmarts("[OH]-[c]"),
    ]
    conf = mol.GetConformer()
    seen = {tuple(c) for c in coords}
    for pat in pats:
        if pat is None:
            continue
        for match in mol.GetSubstructMatches(pat):
            idx = match[0]
            p = conf.GetAtomPosition(idx)
            t = (p.x, p.y, p.z)
            if t in seen:
                continue
            seen.add(t)
            coords.append([p.x, p.y, p.z])
    return coords


def min_oxy_arg(mol: Chem.Mol, arg_atoms: dict, oxy: list[list[float]]) -> float | None:
    if not oxy:
        return None
    nitrogens = [arg_atoms[k] for k in ("NE", "NH1", "NH2") if k in arg_atoms]
    best = float("inf")
    for o in oxy:
        for n in nitrogens:
            d = math.sqrt(sum((o[i] - n[i]) ** 2 for i in range(3)))
            best = min(best, d)
    return None if not math.isfinite(best) else float(best)


def run_gnina(gnina: Path, receptor: Path, ligand: Path, center, size, out_sdf: Path, seed: int, cpu: int) -> None:
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        print(f"SKIP {out_sdf.name}", flush=True)
        return
    if out_sdf.exists() and out_sdf.stat().st_size == 0:
        out_sdf.unlink()
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(gnina), "-r", str(receptor), "-l", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", "32", "--num_modes", "9", "--cpu", str(cpu),
        "--cnn_scoring", "rescore", "--seed", str(seed),
        "-o", str(out_sdf), "--log", str(log), "--no_gpu",
    ]
    print("RUN", out_sdf.name, flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    (out_sdf.parent / (out_sdf.stem + "_stdout.txt")).write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0 or not (out_sdf.exists() and out_sdf.stat().st_size > 0):
        if out_sdf.exists() and out_sdf.stat().st_size == 0:
            out_sdf.unlink(missing_ok=True)
        raise RuntimeError(f"gnina failed {out_sdf}: {(proc.stderr or '')[:400]}")


def evaluate(name: str, sdf: Path, arg: dict, ref_com, meta: dict) -> dict:
    poses = load_poses(sdf)
    if not poses:
        return {"ligand": name, "error": "no_poses", "pass_A1_carboxylate_gate": False}
    i_star = max(range(len(poses)), key=lambda j: _fprop(poses[j], "CNNscore") or -1.0)
    pose = poses[i_star]
    d_co2 = min_acid_arg_dist(pose, arg["atoms"])
    d_broad = min_oxy_arg(pose, arg["atoms"], phenolate_or_acid_oxygens(pose))
    # best carboxylate-Arg among all modes (audit; keep uses CNNscore pose)
    best_co2 = None
    for p in poses:
        d = min_acid_arg_dist(p, arg["atoms"])
        if d is None:
            continue
        best_co2 = d if best_co2 is None else min(best_co2, d)
    com = heavy_centroid(pose)
    d_com = float(math.sqrt(sum((com[i] - ref_com[i]) ** 2 for i in range(3))))
    has_co2 = bool(carboxylate_oxygens(pose))
    pass_arg = d_co2 is not None and d_co2 <= ARG_THRESH
    pass_pocket = d_com <= CENTROID_MAX
    # Official C1 Acid A1-style keep for carboxylate path
    pass_a1 = bool(has_co2 and pass_arg and pass_pocket)
    return {
        "ligand": name,
        "has_carboxylate_in_prep": bool(meta.get("has_carboxylate")),
        "prep_formal_charge": meta.get("formal_charge"),
        "n_poses": len(poses),
        "selected_mode_cnnscore": i_star + 1,
        "CNNscore": _fprop(pose, "CNNscore"),
        "CNNaffinity": _fprop(pose, "CNNaffinity"),
        "has_carboxylate_oxy_in_pose": has_co2,
        "acid_arg477_min_A_cnnscore_pose": d_co2,
        "acid_or_phenolate_arg477_min_A": d_broad,
        "best_carboxylate_arg477_among_modes_A": best_co2,
        "pass_arg_A1_carboxylate": pass_arg,
        "centroid_to_lesinurad_A": d_com,
        "pass_pocket_centroid_6A": pass_pocket,
        "pass_A1_carboxylate_gate": pass_a1,
        "arg_threshold_A": ARG_THRESH,
        "note": (
            "non-carboxylate (phenol etc.): A1 carboxylate gate N/A — see phenolate distance"
            if not meta.get("has_carboxylate")
            else "carboxylate acid; A1 = CO2 + Arg<=7.7027 + COM<=6"
        ),
        "sdf": str(sdf),
    }


def main() -> None:
    eng = yaml.safe_load((PROJECT_ROOT / "config/docking_c1_cpu.yaml").read_text())
    man = pd.read_csv(PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/forced_recovery/ligand_manifest.csv")
    assert set(man["status"]) == {"prepared"}
    tcfg = eng["targets"]["urat1_9dkb"]
    receptor = PROJECT_ROOT / tcfg["prepared_receptor"]
    gnina = PROJECT_ROOT / "tools" / "gnina"
    arg = json.loads((REFS / "arg477_coords.json").read_text())
    ref_com = load_ref_centroid(REFS / "lesinurad_crystal_ref.sdf")
    seed = 42
    cpu = int(eng["gnina"].get("cpu", 4))

    out_dir = OUT / "urat1_9dkb" / f"seed{seed}"
    rows = []
    for _, r in man.iterrows():
        name = str(r["repurposing_id"])
        lig = PROJECT_ROOT / r["pdbqt"]
        sdf = out_dir / f"{name}_out.sdf"
        run_gnina(gnina, receptor, lig, tcfg["center"], tcfg["size"], sdf, seed, cpu)
        rows.append(evaluate(name, sdf, arg, ref_com, r.to_dict()))

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "forced_recovery_A1_gate_table_seed42.csv"
    df.to_csv(csv_path, index=False)

    # markdown summary
    lines = [
        "# Forced-recovery URAT1 A1 gate table (seed 42, free dock)",
        "",
        f"- Protocol: gnina exhaustiveness=32, num_modes=9, CNNscore-selected pose, Arg threshold **{ARG_THRESH} Å** (A1), pocket COM ≤ {CENTROID_MAX} Å vs lesinurad crystal.",
        "- **pass_A1_carboxylate_gate** = carboxylate oxygens present AND Arg≤7.7027 AND COM≤6.",
        "- Phenol-class ligands (benzbromarone, dotinurad): carboxylate gate reported separately; see `acid_or_phenolate_arg477_min_A`.",
        "- This table is an **audit**, not a new scientific lock.",
        "",
        "| ligand | prep CO2 | Arg CO2 (Å) | pass Arg A1 | COM (Å) | pass pocket | **pass A1 gate** | phenolate/acid Arg (Å) | best CO2 Arg any mode (Å) |",
        "|---|:---:|---:|:---:|---:|:---:|:---:|---:|---:|",
    ]
    for _, r in df.iterrows():
        def fmt(x, nd=2):
            if x is None or (isinstance(x, float) and (math.isnan(x) if isinstance(x, float) else False)):
                return "—"
            try:
                if pd.isna(x):
                    return "—"
            except Exception:
                pass
            if isinstance(x, bool):
                return "yes" if x else "no"
            if isinstance(x, float):
                return f"{x:.{nd}f}"
            return str(x)

        lines.append(
            f"| {r['ligand']} | {fmt(r['has_carboxylate_in_prep'])} | {fmt(r['acid_arg477_min_A_cnnscore_pose'])} | "
            f"{fmt(r['pass_arg_A1_carboxylate'])} | {fmt(r['centroid_to_lesinurad_A'])} | {fmt(r['pass_pocket_centroid_6A'])} | "
            f"**{fmt(r['pass_A1_carboxylate_gate'])}** | {fmt(r['acid_or_phenolate_arg477_min_A'])} | "
            f"{fmt(r['best_carboxylate_arg477_among_modes_A'])} |"
        )
    n_pass = int(df["pass_A1_carboxylate_gate"].fillna(False).sum())
    n_co2 = int(df["has_carboxylate_in_prep"].fillna(False).sum())
    lines += [
        "",
        f"**Summary:** {n_pass} / {len(df)} pass A1 carboxylate gate; {n_co2} / {len(df)} prepared with carboxylate.",
        "",
        f"CSV: `{csv_path.relative_to(PROJECT_ROOT)}`",
    ]
    md_path = OUT / "forced_recovery_A1_gate_table_seed42.md"
    md_path.write_text("\n".join(lines) + "\n")
    print(df[["ligand", "pass_A1_carboxylate_gate", "acid_arg477_min_A_cnnscore_pose", "centroid_to_lesinurad_A"]].to_string(index=False))
    print("WROTE", csv_path)
    print("WROTE", md_path)


if __name__ == "__main__":
    main()
