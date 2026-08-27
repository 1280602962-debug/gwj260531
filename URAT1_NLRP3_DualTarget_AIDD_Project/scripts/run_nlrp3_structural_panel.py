#!/usr/bin/env python3
"""NLRP3 structural positive panel vs clinical Acid background.

Small clean experiment (not a giant decoy benchmark):
  positives: NP3-146 crystal + prepared known NACHT-pocket / sulfonylurea-class ligands
  background: random sample of Acid-track clinical chemistry-pass pool

Compares loose keep_nlrp3_pose vs structural keep (overlap + IFP + key recovery).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem
from scipy.stats import fisher_exact

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from gnina_dock import DEFAULT_DOCK_TIMEOUT_S, run_gnina_dock  # noqa: E402
from c1_nlrp3_pose_metrics import (  # noqa: E402
    crystal_reference_ifp,
    evaluate_nlrp3_structural,
    load_key_map,
    load_receptor_heavy,
)

OUT = PROJECT_ROOT / "data/campaigns/c1/05_metrics/nlrp3_structural_panel"
PANEL_CSV = OUT / "panel_ligands.csv"

# Explicit structural / tool compounds (SMILES literature / ChEMBL)
# NP3-146 = 7ALV RM5 co-crystal (Dekker 2021)
# MCC950 = tool sulfonylurea (same class; NOT the co-crystal ligand)
# GDC-2394, DFV890, selnoflast = clinical NACHT-pocket inhibitors
HARDCODED_POSITIVES = [
    {
        "ligand_id": "NP3-146",
        "name": "NP3-146",
        "role": "crystal_positive",
        "smiles": "CC(C)c1cc(Cl)cc(C(C)C)c1NC(=O)NS(=O)(=O)c1cc(C(C)(C)O)co1",
        "note": "7ALV co-crystal RM5",
    },
    {
        "ligand_id": "MCC950",
        "name": "MCC950",
        "role": "tool_positive",
        "smiles": "CC(C)(O)c1coc(S(=O)(=O)NC(=O)Nc2c3c(cc4c2CCC4)CCC3)c1",
        "note": "sulfonylurea tool; analog pocket, not co-crystal",
    },
]


def find_chembl_sulfonylurea_actives(n: int = 8) -> list[dict]:
    nl = pd.read_csv(PROJECT_ROOT / "data/processed/nlrp3_records.csv")
    rel = nl["Standard Relation"].astype(str).str.replace("'", "", regex=False).str.strip()
    nl = nl[rel.eq("=")].copy()
    nl = nl[nl["pChEMBL Value"].notna()]
    nl = nl[nl["pChEMBL Value"] >= 6.0]
    su = Chem.MolFromSmarts("NS(=O)(=O)")
    rows = []
    seen = set()
    for _, r in nl.sort_values("pChEMBL Value", ascending=False).iterrows():
        smi = str(r.get("canonical_smiles") or r.get("Smiles") or "")
        m = Chem.MolFromSmiles(smi)
        if m is None or not m.HasSubstructMatch(su):
            continue
        key = Chem.MolToInchiKey(m)
        if key in seen:
            continue
        seen.add(key)
        name = str(r.get("Molecule Name") or r.get("Molecule ChEMBL ID") or f"CHEMBL_{len(rows)}")
        if name == "nan":
            name = str(r.get("Molecule ChEMBL ID"))
        lid = str(r.get("Molecule ChEMBL ID") or f"NL_{len(rows)}")
        rows.append(
            {
                "ligand_id": lid,
                "name": name,
                "role": "chembl_sulfonylurea_active",
                "smiles": smi,
                "note": f"pChEMBL={r['pChEMBL Value']}",
            }
        )
        if len(rows) >= n:
            break
    return rows


def sample_background(n: int = 20, seed: int = 42) -> list[dict]:
    man = pd.read_csv(
        PROJECT_ROOT
        / "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv"
    )
    pool = pd.read_csv(PROJECT_ROOT / "data/repurposing/screening/docking_pool_p05.csv")
    m = man.merge(
        pool[["repurposing_id", "name"]],
        on="repurposing_id",
        how="left",
        suffixes=("", "_pool"),
    )
    m = m.sample(n=min(n, len(m)), random_state=seed)
    rows = []
    for _, r in m.iterrows():
        rows.append(
            {
                "ligand_id": r["repurposing_id"],
                "name": r.get("name", r["repurposing_id"]),
                "role": "clinical_acid_background",
                "smiles": r["canonical_smiles"],
                "note": "Acid-track chemistry-pass sample",
                "pdbqt": r.get("pdbqt"),
            }
        )
    return rows


def write_panel(path: Path) -> pd.DataFrame:
    rows = HARDCODED_POSITIVES + find_chembl_sulfonylurea_actives(8) + sample_background(20)
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


def prepare_missing_pdbqt(df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Minimal Meeko prep for panel ligands lacking pdbqt."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise SystemExit(f"meeko/rdkit required: {e}")

    pdbqt_dir = out_dir / "pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    prep = MoleculePreparation()
    paths = []
    for _, r in df.iterrows():
        if pd.notna(r.get("pdbqt")) and Path(str(r["pdbqt"])).exists():
            paths.append(r["pdbqt"])
            continue
        outp = pdbqt_dir / f"{r['ligand_id']}.pdbqt"
        if outp.exists() and outp.stat().st_size > 0:
            paths.append(str(outp))
            continue
        mol = Chem.MolFromSmiles(str(r["smiles"]))
        if mol is None:
            paths.append("")
            continue
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
        AllChem.MMFFOptimizeMolecule(mol)
        setup = prep.prepare(mol)[0]
        pdbqt_str, success, err = PDBQTWriterLegacy.write_string(setup)
        if not success:
            print("prep fail", r["ligand_id"], err)
            paths.append("")
            continue
        outp.write_text(pdbqt_str)
        paths.append(str(outp))
    df = df.copy()
    df["pdbqt"] = paths
    return df


def dock_one(
    gnina: Path,
    receptor: Path,
    ligand: Path,
    center,
    size,
    out_sdf: Path,
    seed: int,
    exh: int,
    cpu: int,
    timeout_s: int = DEFAULT_DOCK_TIMEOUT_S,
) -> str | None:
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        return None
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(gnina), "-r", str(receptor), "-l", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exh), "--num_modes", "9", "--cpu", str(cpu),
        "--cnn_scoring", "rescore", "--seed", str(seed), "--no_gpu",
        "-o", str(out_sdf),
    ]
    print(f"RUN panel {out_sdf.name}", flush=True)
    err = run_gnina_dock(cmd, out_sdf, timeout_s=timeout_s)
    if err:
        print(f"FAIL panel {out_sdf.name}: {err[:200]}", flush=True)
    return err


def metrics_table(df: pd.DataFrame) -> dict:
    def rate(sub, col):
        return float(sub[col].mean()) if len(sub) else float("nan")

    pos = df[df.role != "clinical_acid_background"]
    bg = df[df.role == "clinical_acid_background"]
    out = {}
    for col in ("keep_nlrp3_pose", "keep_nlrp3_structural"):
        a = int(pos[col].sum())
        b = int((~pos[col]).sum())
        c = int(bg[col].sum())
        d = int((~bg[col]).sum())
        oddsr, p = fisher_exact([[a, b], [c, d]])
        out[col] = {
            "positive_pass_rate": rate(pos, col),
            "background_pass_rate": rate(bg, col),
            "positive_n": int(len(pos)),
            "background_n": int(len(bg)),
            "tp": a,
            "fn": b,
            "fp": c,
            "tn": d,
            "odds_ratio": float(oddsr),
            "fisher_exact_p": float(p),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit-background", type=int, default=20)
    ap.add_argument("--dock", action="store_true", help="run gnina for missing panel SDFs")
    ap.add_argument("--metrics-only", action="store_true")
    ap.add_argument(
        "--dock-timeout",
        type=int,
        default=DEFAULT_DOCK_TIMEOUT_S,
        help="seconds before killing gnina and skipping ligand (default 1800)",
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    if PANEL_CSV.exists() and args.metrics_only:
        panel = pd.read_csv(PANEL_CSV)
    else:
        panel = write_panel(PANEL_CSV)
        # refresh background size
        if args.limit_background != 40:
            pos = panel[panel.role != "clinical_acid_background"]
            bg = sample_background(args.limit_background)
            panel = pd.concat([pos, pd.DataFrame(bg)], ignore_index=True)
            panel.to_csv(PANEL_CSV, index=False)

    # attach existing clinical pdbqt paths
    man = pd.read_csv(
        PROJECT_ROOT
        / "data/campaigns/c1/01_ligand_prep/acid_clinical_chemistry_pass/ligand_manifest.csv"
    )
    pdbqt_map = dict(zip(man.repurposing_id, man.pdbqt))
    panel["pdbqt"] = panel.apply(
        lambda r: r["pdbqt"] if pd.notna(r.get("pdbqt")) else pdbqt_map.get(r["ligand_id"], None),
        axis=1,
    )
    panel = prepare_missing_pdbqt(panel, OUT)

    eng = yaml.safe_load((PROJECT_ROOT / "config/docking_c1_cpu.yaml").read_text())
    tcfg = eng["targets"]["nlrp3_7alv"]
    receptor = PROJECT_ROOT / tcfg["prepared_receptor"]
    gnina = PROJECT_ROOT / "tools/gnina"
    if not gnina.exists():
        gnina = Path(eng["gnina"].get("binary", "gnina"))

    sdf_dir = OUT / f"poses_seed{args.seed}"
    fail_log = OUT / "panel_dock_failures.jsonl"
    if args.dock and not args.metrics_only:
        for _, r in panel.iterrows():
            if not r.get("pdbqt") or not Path(str(r["pdbqt"])).exists():
                continue
            out_sdf = sdf_dir / f"{r['ligand_id']}_out.sdf"
            # NP3-146: reuse self-dock if present
            if r["ligand_id"] == "NP3-146":
                src = PROJECT_ROOT / f"data/campaigns/c1/02_selfdock/nlrp3_7alv/seed{args.seed}/NP3-146_out.sdf"
                if src.exists():
                    out_sdf.parent.mkdir(parents=True, exist_ok=True)
                    out_sdf.write_bytes(src.read_bytes())
                    continue
            err = dock_one(
                gnina,
                receptor,
                Path(str(r["pdbqt"])),
                tcfg["center"],
                tcfg["size"],
                out_sdf,
                args.seed,
                8,  # panel exploratory exhaustiveness; clinical Acid remains 32
                int(eng["gnina"].get("cpu", 4)),
                timeout_s=args.dock_timeout,
            )
            if err:
                with fail_log.open("a") as f:
                    f.write(
                        json.dumps(
                            {
                                "ligand_id": r["ligand_id"],
                                "seed": args.seed,
                                "error": err[:500],
                            },
                            ensure_ascii=False,
                        )
                        + "\n"
                    )
                print(f"CONTINUE panel after fail {r['ligand_id']}", flush=True)

    key_map = load_key_map()
    ref_ifp, ref_heavy, ref_com = crystal_reference_ifp()
    receptor_heavy = load_receptor_heavy()
    rows = []
    for _, r in panel.iterrows():
        sdf = sdf_dir / f"{r['ligand_id']}_out.sdf"
        if r["ligand_id"] == "NP3-146" and not sdf.exists():
            src = PROJECT_ROOT / f"data/campaigns/c1/02_selfdock/nlrp3_7alv/seed{args.seed}/NP3-146_out.sdf"
            if src.exists():
                sdf.parent.mkdir(parents=True, exist_ok=True)
                sdf.write_bytes(src.read_bytes())
        if not sdf.exists():
            rows.append(
                {
                    "ligand_id": r["ligand_id"],
                    "name": r["name"],
                    "role": r["role"],
                    "error": "missing_sdf",
                    "keep_nlrp3_pose": False,
                    "keep_nlrp3_structural": False,
                }
            )
            continue
        ev = evaluate_nlrp3_structural(
            sdf,
            r["ligand_id"],
            args.seed,
            key_map=key_map,
            ref_heavy=ref_heavy,
            ref_com=ref_com,
            receptor_heavy=receptor_heavy,
            ref_ifp=ref_ifp,
        )
        ev["name"] = r["name"]
        ev["role"] = r["role"]
        ev["note"] = r.get("note")
        rows.append(ev)

    res = pd.DataFrame(rows)
    res.to_csv(OUT / f"nlrp3_panel_metrics_seed{args.seed}.csv", index=False)
    # drop missing for stats
    ok = res[res["error"].isna()] if "error" in res.columns else res
    if "CNNscore" in ok.columns:
        ok = ok[ok["CNNscore"].notna()]
    summary = {
        "seed": args.seed,
        "gates": {
            "loose": "COM<=6A and CNNscore>=0.5",
            "structural": "loose + overlap>=0.50 + IFP Jaccard>=0.50 + key_contacts>=5/7 + no clash",
        },
        "metrics": metrics_table(ok) if len(ok) else {},
        "n_rows": int(len(res)),
        "n_scored": int(len(ok)),
        "claim": "structural recovery enrichment only — not binding affinity proof",
    }
    (OUT / f"nlrp3_panel_summary_seed{args.seed}.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
