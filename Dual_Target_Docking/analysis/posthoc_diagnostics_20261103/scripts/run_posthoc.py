#!/usr/bin/env python3
"""Run the post-hoc diagnostics twice and write deterministic numeric outputs."""
from __future__ import annotations

import json
import platform
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy
import pandas
import pytest
import rdkit
import scipy
import sklearn

from analysis_a import build_analysis_a
from analysis_b import build_analysis_b
from analysis_c import build_analysis_c
from figures import make_all_figures
from posthoc_common import (
    DIAG_ROOT,
    INPUTS,
    OUTPUT_DIR,
    PROJECT,
    QA_DIR,
    hashes_of_frames,
    input_manifest,
    load_config,
    sha256_file,
    verify_schemes,
    write_frame,
)
from reports import write_reports


OUTPUT_NAMES = [
    "A_FOLD_OVERLAP_SUMMARY",
    "A_FOLD_OVERLAP_PER_FOLD",
    "A_TEST_NEAREST_TRAIN",
    "A_CHEMICAL_ENTITY_DUPLICATES",
    "A_SCAFFOLD_FIELD_AUDIT",
    "A_TRIGGER_DECISION",
    "B_EGFR_SCORE_DISTRIBUTIONS",
    "B_EGFR_PAIRWISE_RANKING",
    "B_EGFR_DESCRIPTOR_COMPARISONS",
    "B_EGFR_SCORE_DESCRIPTOR_ASSOCIATIONS",
    "B_EGFR_ACTIVITY_DISTRIBUTIONS",
    "B_EGFR_ACTIVITY_MOLECULES",
    "B_EGFR_PROVENANCE_FLAGS",
    "B_SCORE_POINTS",
    "B_DESCRIPTOR_POINTS",
    "C_TOPK_SENSITIVITY_ALL",
    "C_TOPK_METHOD_COMPARISONS",
    "C_TOPK_CONCORDANCE_COUNTS",
    "C_TOPK_BASELINE_REPRODUCTION",
]


def build_once() -> dict:
    frames = {}
    frames.update(build_analysis_a())
    frames.update(build_analysis_b())
    frames.update(build_analysis_c())
    if "C_TOPK_METHOD_COMPARISONS" not in frames:
        raise SystemExit("Top-10% baseline was not reproduced; sensitivity results were not accepted")
    return frames


def numeric_subset(frames: dict) -> dict:
    return {name: frames[name] for name in OUTPUT_NAMES}


def assert_baseline_counts() -> None:
    cfg = load_config()
    pop = pandas.read_csv(PROJECT / INPUTS["population"])
    eligible = pop[pop["activity_eligible"] == 1]
    expected = cfg["expected_counts"]
    counts = eligible["class"].value_counts().to_dict()
    if len(eligible) != expected["eligible_observations"]:
        raise SystemExit("eligible observation count changed")
    if eligible["global_ligand_entity_id"].nunique() != expected["global_entities"]:
        raise SystemExit("global entity count changed")
    for label in ("dual", "A_only", "B_only", "neither"):
        if int(counts[label]) != expected[label]:
            raise SystemExit(f"class count changed: {label}")
    units = pandas.read_csv(PROJECT / INPUTS["chemistry_units"])
    if int(units["n_labeled_population"].sum()) != expected["labeled_directional"]:
        raise SystemExit("labeled directional count changed")
    if int(units["n_m0_directional_population"].sum()) != expected["m0_directional"]:
        raise SystemExit("M0 directional count changed")
    if int(units["n_final_chemistry_oof"].sum()) != expected["chemistry_oof"]:
        raise SystemExit("chemistry OOF count changed")


def main() -> int:
    verify_schemes()
    assert_baseline_counts()
    QA_DIR.mkdir(parents=True, exist_ok=True)
    before = input_manifest()
    before.to_csv(QA_DIR / "INPUT_SHA256_BEFORE.csv", index=False)
    first = build_once()
    second = build_once()
    if hashes_of_frames(numeric_subset(first)) != hashes_of_frames(numeric_subset(second)):
        raise SystemExit("numeric outputs were not deterministic across two runs")
    for name in OUTPUT_NAMES:
        write_frame(first[name], f"{name}.csv")
    pandas.concat([first["A_JOIN_AUDIT"], first["B_JOIN_AUDIT"]], ignore_index=True).to_csv(
        QA_DIR / "JOIN_AUDIT.csv", index=False
    )
    first["A_TRIGGER_DECISION"].to_csv(QA_DIR / "A_TRIGGER_DECISION.csv", index=False)
    pandas.DataFrame([
        {"file": f"{name}.csv", "sha256": digest}
        for name, digest in hashes_of_frames(numeric_subset(first)).items()
    ]).to_csv(QA_DIR / "NUMERIC_SHA256.csv", index=False)
    make_all_figures()
    write_reports(first)
    after = input_manifest()
    after.to_csv(QA_DIR / "INPUT_SHA256_AFTER.csv", index=False)
    if before["sha256"].tolist() != after["sha256"].tolist():
        raise SystemExit("input SHA-256 changed during the analysis")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT).decode().strip()
    branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=PROJECT).decode().strip()
    env = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "pytest": pytest.__version__,
        "rdkit": rdkit.__version__,
        "git_head_at_runtime": head,
        "git_branch": branch,
        "baseline_commit": load_config()["baseline_commit"],
        "random_seed": 271828,
        "determinism": "PASS",
    }
    (QA_DIR / "ENVIRONMENT.json").write_text(json.dumps(env, indent=2) + "\n")
    script_rows = []
    for path in sorted((DIAG_ROOT / "scripts").glob("*.py")):
        script_rows.append({"relative_path": str(path.relative_to(PROJECT)), "sha256": sha256_file(path)})
    pandas.DataFrame(script_rows).to_csv(QA_DIR / "SCRIPT_SHA256.csv", index=False)
    (QA_DIR / "TEST_STATUS.txt").write_text("numeric_determinism=PASS\ninput_hash_unchanged=PASS\n")
    print(f"wrote {len(list(OUTPUT_DIR.glob('*.csv')))} CSV files to {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
