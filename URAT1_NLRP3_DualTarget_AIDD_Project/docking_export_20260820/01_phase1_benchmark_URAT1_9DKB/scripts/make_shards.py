#!/usr/bin/env python3
"""Split ligand index range into shard files for array jobs."""
from __future__ import annotations

import argparse
import math
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mols", type=int, default=0, help="0 = count ligands_sdf")
    ap.add_argument("--ligands-dir", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--shard-size", type=int, default=50)
    args = ap.parse_args()

    n = args.n_mols
    if n <= 0:
        d = Path(args.ligands_dir)
        n = len(list(d.glob("mol_*.sdf"))) if d.is_dir() else 0
    if n <= 0:
        raise SystemExit("Provide --n-mols or --ligands-dir with mol_*.sdf")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    n_shards = int(math.ceil(n / args.shard_size))
    for s in range(n_shards):
        start = s * args.shard_size
        end = min(n, start + args.shard_size)
        (out / f"shard_{s:04d}.txt").write_text(f"{start} {end}\n")
    (out / "shards.tsv").write_text(
        "shard_id\tstart\tend\n"
        + "\n".join(
            f"{s}\t{s * args.shard_size}\t{min(n, (s + 1) * args.shard_size)}"
            for s in range(n_shards)
        )
        + "\n"
    )
    print(f"n_mols={n} shard_size={args.shard_size} n_shards={n_shards}")
    print(f"SLURM hint: --array=0-{n_shards - 1}")


if __name__ == "__main__":
    main()
