#!/usr/bin/env python3
"""C1 Acid-track dual docking (9DKB + 7ALV) with preflight checks.

Scientific lock: config/campaign_c1.yaml (Amendment A1)
Engine config: config/docking_c1_cpu.yaml (local) or docking_c1.yaml (GPU)

Hard rules enforced here:
  - NEVER pass campaign_c1.yaml to gnina
  - NEVER use docking percentiles for keep/drop
  - Arg gate = crystal_min + 1.0 = 7.7027 Å
  - Rank / L3 full decoy not opened
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from parse_c1_sdf_readouts import (  # noqa: E402
    load_poses,
    min_acid_arg_dist,
    _fprop,
    carboxylate_oxygens,
    heavy_centroid,
    load_ref_centroid,
)


ARG_THRESH = 7.7027
CRYSTAL_MIN = 6.7027
CENTROID_MAX_A = 6.0  # Phase I pose-QC style pocket proxy vs crystal ligand COM
REFS = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs"
CRYSTAL_COM = {
    "urat1_9dkb": REFS / "lesinurad_crystal_ref.sdf",
    "nlrp3_7alv": REFS / "NP3-146_RM5_crystal_ref.sdf",
}


def preflight(engine_yaml: Path, manifest: Path) -> dict:
    camp = yaml.safe_load((PROJECT_ROOT / "config/campaign_c1.yaml").read_text())
    assert camp.get("active_track") == "Acid", "campaign not on Acid track"
    assert "closed" in str(camp.get("rank_track_status", "")), "Rank must stay closed"
    acid = camp["acid_track"]
    thr = float(acid["urat1_acid_arg477_le_A"])
    assert abs(thr - ARG_THRESH) < 1e-6, f"unexpected Arg threshold {thr}"
    assert abs(float(acid["crystal_min_O_Arg477_N_A"]) - CRYSTAL_MIN) < 1e-3

    # engine yaml must be docking schema, not campaign lock
    eng = yaml.safe_load(engine_yaml.read_text())
    assert "targets" in eng and "gnina" in eng, "engine yaml missing targets/gnina"
    assert "acid_track" not in eng, "refusing to use campaign_c1.yaml as engine config"
    assert int(eng["gnina"].get("num_modes", 0)) == 9

    man = pd.read_csv(manifest)
    n_ok = int((man["status"] == "prepared").sum())
    assert n_ok > 0, "no prepared ligands"
    for col in ("repurposing_id", "pdbqt", "canonical_smiles"):
        assert col in man.columns, f"manifest missing {col}"

    for key in ("urat1_9dkb", "nlrp3_7alv"):
        rec = PROJECT_ROOT / eng["targets"][key]["prepared_receptor"]
        assert rec.exists(), f"missing receptor {rec}"

    arg_json = PROJECT_ROOT / "data/campaigns/c1/01_ligand_prep/selfdock_refs/arg477_coords.json"
    assert arg_json.exists(), "missing Arg477 coords; run extract_c1_crystal_refs.py"

    gnina = PROJECT_ROOT / "tools" / "gnina"
    assert gnina.exists() or True  # may be on PATH via wrapper

    info = {
        "active_track": camp["active_track"],
        "arg_threshold_A": thr,
        "crystal_min_A": CRYSTAL_MIN,
        "n_prepared": n_ok,
        "engine_yaml": str(engine_yaml),
        "num_modes": eng["gnina"]["num_modes"],
        "no_gpu": bool(eng["gnina"].get("no_gpu", True)),
        "rank_closed": True,
        "percentile_ranking": False,
    }
    print("PREFLIGHT_OK", json.dumps(info, indent=2), flush=True)
    return {"camp": camp, "eng": eng, "manifest": man, "arg_json": arg_json}


def run_one(gnina: Path, receptor: Path, ligand: Path, center, size, out_sdf: Path, seed: int, exh: int, modes: int, cpu: int, no_gpu: bool) -> str | None:
    """Run gnina. Returns None on success, error string on failure (does not abort batch)."""
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        print(f"SKIP {out_sdf}", flush=True)
        return None
    if out_sdf.exists() and out_sdf.stat().st_size == 0:
        out_sdf.unlink()
    log = out_sdf.with_suffix(".log")
    cmd = [
        str(gnina), "-r", str(receptor), "-l", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exh), "--num_modes", str(modes), "--cpu", str(cpu),
        "--cnn_scoring", "rescore", "--seed", str(seed),
        "-o", str(out_sdf), "--log", str(log),
    ]
    if no_gpu:
        cmd.append("--no_gpu")
    print("RUN", out_sdf.name, flush=True)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired as e:
        if out_sdf.exists() and out_sdf.stat().st_size == 0:
            out_sdf.unlink(missing_ok=True)
        err = f"timeout_7200s: {out_sdf.name}"
        print(f"FAIL {out_sdf.name}: {err}", flush=True)
        return err
    (out_sdf.parent / (out_sdf.stem + "_stdout.txt")).write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0 or not (out_sdf.exists() and out_sdf.stat().st_size > 0):
        if out_sdf.exists() and out_sdf.stat().st_size == 0:
            out_sdf.unlink(missing_ok=True)
        err = (proc.stderr or proc.stdout or "gnina_failed")[:400]
        print(f"FAIL {out_sdf.name}: {err[:200]}", flush=True)
        return err
    return None


def evaluate_acid_pose(sdf: Path, arg_json: Path, ligand_id: str, target: str, seed: int, pose_selection: str = "a1") -> dict:
    """URAT1 (A1/A2) or NLRP3 pose gates; no percentiles."""
    from c1_acid_pose_selection import evaluate_nlrp3_pose_sdf, evaluate_urat1_acid_sdf, load_ref_centroid

    if target == "urat1_9dkb":
        ref_com = load_ref_centroid(CRYSTAL_COM[target])
        row = evaluate_urat1_acid_sdf(
            sdf, arg_json, ref_com, ligand_id, seed, rule=pose_selection  # type: ignore[arg-type]
        )
        row["keep"] = row.get("keep_urat1_acid", False)
        return row
    ref_com = load_ref_centroid(CRYSTAL_COM[target])
    row = evaluate_nlrp3_pose_sdf(sdf, ref_com, ligand_id, seed)
    row["keep"] = row.get("keep_nlrp3_pose", False)
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--engine-config", type=Path, default=PROJECT_ROOT / "config/docking_c1_cpu.yaml")
    ap.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "data/campaigns/c1/07_clinical_dock/acid_dual")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--ids", nargs="*", default=None, help="optional repurposing_id subset")
    ap.add_argument("--targets", nargs="+", default=["urat1_9dkb", "nlrp3_7alv"])
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--metrics-only", action="store_true", help="re-score existing SDFs; do not call gnina")
    ap.add_argument(
        "--pose-selection",
        choices=["a1", "a2"],
        default="a1",
        help="URAT1 pose rule; a2 = geometry-first (Amendment A2)",
    )
    args = ap.parse_args()

    pf = preflight(args.engine_config, args.manifest)
    if args.preflight_only:
        return

    eng = pf["eng"]
    man = pf["manifest"]
    man = man[man["status"] == "prepared"].copy()
    if args.ids:
        # preserve caller ID order (isin alone follows manifest order)
        man = man.set_index("repurposing_id").loc[[i for i in args.ids if i in set(man["repurposing_id"])]].reset_index()
    if args.limit:
        man = man.head(args.limit)

    for ref in CRYSTAL_COM.values():
        assert ref.exists(), f"missing crystal ref {ref}"

    gnina = PROJECT_ROOT / "tools" / "gnina"
    if not gnina.exists():
        gnina = Path(eng["gnina"].get("binary", "gnina"))

    rows = []
    fail_log = args.output_dir / "dock_failures.jsonl"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def record_fail(rid: str, tkey: str, err: str) -> None:
        row = {
            "ligand_id": rid, "target": tkey, "seed": args.seed,
            "error": err, "keep": False, "keep_urat1_acid": False, "keep_nlrp3_pose": False,
        }
        rows.append(row)
        with fail_log.open("a") as f:
            f.write(json.dumps({"ligand_id": rid, "target": tkey, "seed": args.seed, "error": err[:500]}, ensure_ascii=False) + "\n")
        print(f"CONTINUE after fail {rid}/{tkey}", flush=True)

    for _, r in man.iterrows():
        rid = str(r["repurposing_id"])
        lig = Path(r["pdbqt"])
        if not lig.is_absolute():
            lig = PROJECT_ROOT / lig
        for tkey in args.targets:
            tcfg = eng["targets"][tkey]
            receptor = PROJECT_ROOT / tcfg["prepared_receptor"]
            out_sdf = args.output_dir / tkey / f"seed{args.seed}" / f"{rid}_out.sdf"
            if not args.metrics_only:
                err = run_one(
                    gnina, receptor, lig, tcfg["center"], tcfg["size"], out_sdf,
                    args.seed,
                    int(eng["gnina"]["exhaustiveness"]),
                    int(eng["gnina"]["num_modes"]),
                    int(eng["gnina"].get("cpu", 4)),
                    bool(eng["gnina"].get("no_gpu", True)),
                )
                if err:
                    record_fail(rid, tkey, err)
                    continue
            elif not (out_sdf.exists() and out_sdf.stat().st_size > 0):
                record_fail(rid, tkey, "missing_sdf")
                continue
            try:
                rows.append(evaluate_acid_pose(out_sdf, pf["arg_json"], rid, tkey, args.seed, args.pose_selection))
            except Exception as e:
                record_fail(rid, tkey, f"evaluate_error: {type(e).__name__}: {e}")

    df = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_dir / f"acid_pose_metrics_seed{args.seed}.csv", index=False)

    # join dual keep without percentiles
    if set(args.targets) >= {"urat1_9dkb", "nlrp3_7alv"}:
        u = df[df.target == "urat1_9dkb"][["ligand_id", "keep_urat1_acid", "acid_arg477_min_A", "CNNscore", "CNNaffinity"]].rename(
            columns={"CNNscore": "u_CNNscore", "CNNaffinity": "u_CNNaffinity"}
        )
        n = df[df.target == "nlrp3_7alv"][["ligand_id", "keep_nlrp3_pose", "CNNscore", "CNNaffinity"]].rename(
            columns={"CNNscore": "n_CNNscore", "CNNaffinity": "n_CNNaffinity"}
        )
        dual = u.merge(n, on="ligand_id", how="outer")
        dual["keep_dual_acid_geometry"] = dual["keep_urat1_acid"].fillna(False) & dual["keep_nlrp3_pose"].fillna(False)
        dual.to_csv(args.output_dir / f"acid_dual_keep_seed{args.seed}.csv", index=False)
        summary = {
            "seed": args.seed,
            "n_ligands": int(dual.shape[0]),
            "n_keep_urat1_arg": int(dual["keep_urat1_acid"].fillna(False).sum()),
            "n_keep_nlrp3_pose": int(dual["keep_nlrp3_pose"].fillna(False).sum()),
            "n_keep_dual": int(dual["keep_dual_acid_geometry"].sum()),
            "arg_threshold_A": ARG_THRESH,
            "percentile_used": False,
            "note": "keep flags are geometry/pose only; chemistry/qN already applied upstream in acid pool",
        }
        (args.output_dir / f"acid_dual_summary_seed{args.seed}.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
