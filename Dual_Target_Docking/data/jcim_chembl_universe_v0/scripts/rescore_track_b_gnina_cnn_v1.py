#!/usr/bin/env python3
"""GNINA CNN best-of-9 rescore of production Track B Vina poses.

Same flags as K=4: --cnn_scoring rescore --minimize --no_gpu.
Does not run independent GNINA search. Does not replace Table 2.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "local_track_b_v0"
BENCH = ROOT.parent / "jcim_bench_v0" / "scripts"
sys.path.insert(0, str(BENCH))
import gnina_rescore_panel as gn  # noqa: E402

TARGETS = ["4UDW", "2JKH", "6N7A", "3LXP", "8BXH", "9V8H", "6LXA", "5U3Q"]
PAIR_OF = {
    "4UDW": "F2/F10",
    "2JKH": "F2/F10",
    "6N7A": "JAK1/TYK2+JAK1/JAK2",
    "3LXP": "JAK1/TYK2",
    "8BXH": "JAK1/JAK2",
    "9V8H": "PPARG/PPARA",
    "6LXA": "PPARG/PPARA+PPARA/PPARD",
    "5U3Q": "PPARA/PPARD",
}


def receptor_for(target: str) -> Path:
    for name in (f"{target}_protein.pdb", f"{target}_receptor.pdb", f"{target}_receptor.pdbqt"):
        p = LOCAL / "receptors" / name
        if p.exists():
            return p
    raise SystemExit(f"missing receptor for {target} under {LOCAL / 'receptors'}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--gnina", type=Path, default=None)
    ap.add_argument("--obabel", type=Path, default=None)
    args = ap.parse_args()
    if not (LOCAL / "poses").exists():
        raise SystemExit("local_track_b_v0/poses/ missing — regenerate production Vina first")
    if args.gnina:
        gn.GNINA_BIN = args.gnina
        gn.GNINA_ROOT = args.gnina.parent.parent
    if args.obabel:
        gn.OBABEL = args.obabel

    rec_map = [f"{t}={receptor_for(t)}" for t in TARGETS]
    sys.argv = [
        "gnina_rescore_panel.py",
        "--root",
        str(LOCAL),
        "--targets",
        *TARGETS,
        "--receptor-map",
        *rec_map,
        "--workers",
        str(args.workers),
        "--timeout",
        str(args.timeout),
        "--modes",
        "all",
    ]
    rc = gn.main()
    long_csv = LOCAL / "tables" / "scores_gnina_long.csv"
    best_csv = LOCAL / "tables" / "scores_gnina_best.csv"
    if long_csv.exists():
        rows = list(csv.DictReader(long_csv.open()))
        tagged = LOCAL / "tables" / "scores_gnina_cnn_best9_v1.csv"
        # keep the K=4-style best file; also write a pair-tagged long copy
        with tagged.open("w", newline="") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["pair_hint", "target", "ligand", "cnn_score", "cnn_affinity"],
            )
            w.writeheader()
            best = {}
            for r in rows:
                if r.get("status") not in ("success", "exists"):
                    continue
                key = (r["ligand"], r["target"])
                score = r.get("cnn_score")
                aff = r.get("cnn_affinity")
                rank = float(score) if score not in ("", None) else (float(aff) if aff not in ("", None) else None)
                if rank is None:
                    continue
                if key not in best or rank > best[key][0]:
                    best[key] = (rank, score, aff)
            for (lig, t), (_, score, aff) in sorted(best.items()):
                w.writerow(
                    {
                        "pair_hint": PAIR_OF.get(t, ""),
                        "target": t,
                        "ligand": lig,
                        "cnn_score": score,
                        "cnn_affinity": aff,
                    }
                )
        print("wrote", tagged, flush=True)
    print("also", best_csv, flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
