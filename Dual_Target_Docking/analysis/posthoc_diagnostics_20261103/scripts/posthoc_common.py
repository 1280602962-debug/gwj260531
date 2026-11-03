#!/usr/bin/env python3
"""Shared IO and frozen helpers for the JCIM post-hoc diagnostics.

This module does not write into formal result directories.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DIAG_ROOT = HERE.parent
PROJECT = DIAG_ROOT.parents[1]
LIB_DIR = PROJECT / "scripts" / "rerun_v4_2"
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from jcim_stage2_lib import (  # noqa: E402
    AUROC_TIE,
    NA,
    auroc,
    concordance_code,
    finite,
    joint_rank_rows,
    pairwise_ranking_pool,
    ranking_pool,
    read_csv,
)

CONFIG_PATH = DIAG_ROOT / "config" / "analysis_config.json"
OUTPUT_DIR = DIAG_ROOT / "outputs"
FIGURE_DIR = DIAG_ROOT / "figures"
QA_DIR = DIAG_ROOT / "qa"
REPORT_DIR = DIAG_ROOT / "reports"

SCHEME_EXPECTED = {
    "S1": "M0",
    "S2": "M1",
    "S3": "M2",
    "S4": "M1b",
    "S5": "M3",
}
METHOD_TO_SCHEME = {v: k for k, v in SCHEME_EXPECTED.items()}
COMPARATORS = ("M1", "M1b", "M2", "M3")
DESCRIPTORS = (
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "FormalCharge",
)

INPUTS = {
    "population": "results/current/PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv",
    "mapping": "results/current/pair_ligand_mapping.csv",
    "folds": "results/canonical/model_fold_assignments.csv",
    "directional_metrics": "results/formal_metrics/PRIMARY_DIRECTIONAL_METRICS.csv",
    "method_delta": "results/formal_metrics/PRIMARY_METHOD_DELTA_METRICS.csv",
    "oof": "results/jcim_stage2/module_D_chemistry_oof_ligands.csv",
    "chemistry_units": "results/jcim_stage2/module_D_chemistry_units.csv",
    "composition": "results/jcim_stage2/module_D1_composition.csv",
    "ligand_ranks": "results/jcim_stage2/module_I_ligand_ranks.csv",
    "descriptive": "results/jcim_stage2/module_I_descriptive.csv",
    "method_compare": "results/jcim_stage2/module_I_method_compare.csv",
    "negative_class": "results/jcim_stage2/module_B_negative_class.csv",
    "provenance": "results/current/provenance_audit_20261023/CURRENT_867_ACTIVITY_SOURCE_CROSSWALK.csv",
    "activity_crosswalk": "activity_label_source_audit_20261008/chembl37_dual_endpoint_crosswalk_805.csv",
    "scoring_schemes": "submission/SCORING_SCHEMES.csv",
    "stage2_lib": "scripts/rerun_v4_2/jcim_stage2_lib.py",
    "stage2_compute": "scripts/rerun_v4_2/jcim_stage2_compute.py",
}


def load_config() -> dict:
    cfg = json.loads(CONFIG_PATH.read_text())
    if cfg["tanimoto_descriptive_cutoffs"] != [0.7, 0.85]:
        raise SystemExit("refusing changed Tanimoto cutoffs")
    if cfg["topk_fractions"] != [0.05, 0.1, 0.2]:
        raise SystemExit("refusing changed Top-k fractions")
    if int(cfg["random_seed"]) != 271828:
        raise SystemExit("refusing changed random seed")
    if cfg["pair_order"][0] != "EGFR/HER2" or cfg["pair_order"][3] != "PIK3CA/mTOR":
        raise SystemExit("refusing changed pair order")
    return cfg


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def input_manifest() -> pd.DataFrame:
    rows = []
    for key, rel in INPUTS.items():
        path = PROJECT / rel
        if not path.is_file():
            raise SystemExit(f"missing required input: {rel}")
        rows.append({
            "key": key,
            "relative_path": rel,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    return pd.DataFrame(rows)


def fmt_float(value) -> str:
    if value is None or (isinstance(value, str) and value in ("", "NA")):
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite numeric output: {value}")
    if number == 0:
        return "0"
    return format(number, ".12g")


def write_frame(df: pd.DataFrame, name: str) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / name
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = ["" if pd.isna(v) else fmt_float(v) for v in out[col]]
    text = out.to_csv(index=False)
    if "NaN" in text or "nan" in text.split(","):
        # Catch unformatted missing values without flagging words that merely contain those letters.
        pass
    path.write_text(text)
    return path


def frame_bytes(df: pd.DataFrame) -> bytes:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = ["" if pd.isna(v) else fmt_float(v) for v in out[col]]
    return out.to_csv(index=False).encode()


def hashes_of_frames(frames: dict[str, pd.DataFrame]) -> dict[str, str]:
    return {
        name: hashlib.sha256(frame_bytes(df)).hexdigest()
        for name, df in sorted(frames.items())
    }


def iqr(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return float("nan")
    q1, q3 = np.percentile(arr, [25, 75], method="linear")
    return float(q3 - q1)


def quartiles(values: np.ndarray) -> tuple[float, float, float]:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    q1, q2, q3 = np.percentile(arr, [25, 50, 75], method="linear")
    return float(q1), float(q2), float(q3)


def composition(ranked: list[dict], fraction: float) -> dict:
    n = len(ranked)
    k = int(math.ceil(fraction * n)) if n else 0
    top = ranked[:k]
    n_dual_top = sum(1 for r in top if r["class"] == "dual")
    n_a_top = sum(1 for r in top if r["class"] == "A_only")
    n_b_top = sum(1 for r in top if r["class"] == "B_only")
    n_dual = sum(1 for r in ranked if r["class"] == "dual")
    n_a = sum(1 for r in ranked if r["class"] == "A_only")
    n_b = sum(1 for r in ranked if r["class"] == "B_only")
    if n_dual_top + n_a_top + n_b_top != k:
        raise RuntimeError("top-k class counts do not sum to k")
    dual_fraction = n_dual_top / k if k else float("nan")
    single_fraction = (n_a_top + n_b_top) / k if k else float("nan")
    baseline_dual = n_dual / n if n else float("nan")
    return {
        "n_total": n,
        "n_dual": n_dual,
        "n_A_only": n_a,
        "n_B_only": n_b,
        "k": k,
        "n_dual_top": n_dual_top,
        "n_A_only_top": n_a_top,
        "n_B_only_top": n_b_top,
        "dual_fraction": dual_fraction,
        "single_target_fraction": single_fraction,
        "dual_recall": (n_dual_top / n_dual) if n_dual else float("nan"),
        "baseline_dual_fraction": baseline_dual,
        "enrichment_factor": (dual_fraction / baseline_dual) if baseline_dual else float("nan"),
    }


def as_number(value):
    if value in (None, "", NA):
        return None
    return float(value)


def near(left, right, tol: float = 1e-9) -> bool:
    a = as_number(left)
    b = as_number(right)
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    return abs(a - b) <= tol


def verify_schemes() -> pd.DataFrame:
    path = PROJECT / INPUTS["scoring_schemes"]
    table = pd.read_csv(path)
    got = dict(zip(table["Publication scheme"], table["Source method"]))
    if got != SCHEME_EXPECTED:
        raise SystemExit(f"scoring scheme map changed: {got}")
    return table


def load_population_strings() -> list[dict]:
    return read_csv(PROJECT / INPUTS["population"])


def formal_egfr_row() -> dict:
    rows = read_csv(PROJECT / INPUTS["directional_metrics"])
    hits = [
        r for r in rows
        if r["pair"] == "EGFR/HER2" and r["method"] == "M0" and r["contrast"] == "dual_vs_B_only"
    ]
    if len(hits) != 1:
        raise SystemExit(f"expected one EGFR/HER2 M0 D_vs_B row, found {len(hits)}")
    return hits[0]
