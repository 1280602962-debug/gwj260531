#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/mnt/d/CADD paper exercise/dual target docking/results/egfr_her2_panel40_v0")
OUT = ROOT / "analysis" / "exhaustiveness_sensitivity_v1"
TABLES = OUT / "tables"
POSES = OUT / "poses"
LOGS = OUT / "logs"
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_EXAMPLE = RTM_ROOT / "example" / "rtmscore.py"
MODEL = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"


def subset_ligands() -> list[str]:
    with (TABLES / "subset_ligands.csv").open() as fh:
        return [r["ligand_id"] for r in csv.DictReader(fh)]


def convert_pose_to_sdf_block(pose_path: Path, title: str) -> str:
    cmd = (
        "source /home/gwj/miniconda3/etc/profile.d/conda.sh && "
        "conda activate cadd_tools && "
        f"obabel -ipdbqt \"{pose_path}\" -osdf"
    )
    proc = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True, check=True)
    text = proc.stdout
    parts = text.split("$$$$")
    block = parts[0].rstrip() + "\n$$$$\n"
    lines = block.splitlines()
    if lines:
        lines[0] = title
    return "\n".join(lines) + ("\n" if not block.endswith("\n") else "")


def main() -> int:
    subset = subset_ligands()
    rows = []
    for exhaustiveness in [8, 16, 32]:
        for target in ["3POZ", "3RCD"]:
            combined_sdf = OUT / "logs" / f"rtm_input_E{exhaustiveness}_{target}.sdf"
            score_prefix = OUT / "logs" / f"rtm_E{exhaustiveness}_{target}"
            blocks = []
            for lig in subset:
                for mode in range(1, 10):
                    pose = POSES / f"E{exhaustiveness}_seed20260727" / target / lig / f"mode_{mode:02d}.pdbqt"
                    blocks.append(convert_pose_to_sdf_block(pose, f"{lig}_mode{mode}"))
            combined_sdf.write_text("".join(blocks))
            pocket = ROOT / "receptors" / f"{target}_pocket_10.0.pdb"
            cmd = (
                "source /home/gwj/miniconda3/etc/profile.d/conda.sh && "
                "conda activate rtmscore && "
                f"cd \"{RTM_ROOT / 'example'}\" && "
                f"python \"{RTM_EXAMPLE}\" "
                f"-p \"{pocket}\" "
                f"-l \"{combined_sdf}\" "
                f"-m \"{MODEL}\" "
                f"-o \"{score_prefix}\""
            )
            log_path = LOGS / f"rtm_E{exhaustiveness}_{target}.log"
            with log_path.open("w") as log_fh:
                proc = subprocess.run(["bash", "-lc", cmd], stdout=log_fh, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                print(f"RTM failed for E{exhaustiveness} {target}; see {log_path}")
                return 1
            csv_path = Path(f"{score_prefix}.csv")
            df = pd.read_csv(csv_path)
            for _, row in df.iterrows():
                lig, mode = str(row["id"]).split("_mode")
                rows.append(
                    {
                        "target": target,
                        "ligand_id": lig,
                        "exhaustiveness": exhaustiveness,
                        "seed": 20260727,
                        "mode": int(mode),
                        "rtmscore": float(row["score"]),
                    }
                )
    long_df = pd.DataFrame(rows)
    best = (
        long_df.sort_values(["target", "ligand_id", "exhaustiveness", "rtmscore"], ascending=[True, True, True, False])
        .groupby(["target", "ligand_id", "exhaustiveness", "seed"], as_index=False)
        .first()
        .rename(columns={"mode": "best_rtm_mode"})
    )
    best.to_csv(TABLES / "scores_rtm_experimentA.csv", index=False)
    long_df.to_csv(TABLES / "scores_rtm_experimentA_long.csv", index=False)
    print("Wrote", TABLES / "scores_rtm_experimentA.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
