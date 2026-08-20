"""Shared docking score column aliases (gnina / Vina). Do not ingest Glide columns."""
from __future__ import annotations

import re
from typing import Iterable

import numpy as np
import pandas as pd

# Canonical column written by new pipelines
CANONICAL_SCORE_COL = "dock_score"

# Accepted input aliases (case-insensitive match on normalized header)
DOCK_SCORE_ALIASES = [
    CANONICAL_SCORE_COL,
    "vina_score",
    "vina affinity",
    "minimizedaffinity",
    "minimized_affinity",
    "binding_affinity",
    "affinity",
    "cnnaffinity",
    "cnn_affinity",
    "docking score",
    "score",
]

SMILES_ALIASES = [
    "canonical_smiles",
    "smiles",
]

STATUS_ALIASES = [
    "docking_status",
    "pose",
    "pose_status",
    "status",
]

NAME_ALIASES = ["name", "pref_name", "title", "compound_name", "ligand"]

REPURPOSING_ALIASES = ["repurposing_id"]

ENGINE_ALIASES = ["docking_engine", "engine", "tool"]


def _norm_header(s: str) -> str:
    return re.sub(r"\s+", " ", str(s).strip().lower())


def pick_col(columns: Iterable[str], aliases: list[str]) -> str | None:
    norm = {_norm_header(c): c for c in columns}
    for a in aliases:
        key = _norm_header(a)
        if key in norm:
            return norm[key]
    cols = list(columns)
    for c in cols:
        cn = _norm_header(c).replace(" ", "").replace("_", "")
        for a in aliases:
            an = _norm_header(a).replace(" ", "").replace("_", "")
            if an and an in cn:
                return c
    return None


def coerce_dock_score(series: pd.Series) -> pd.Series:
    """Lower (more negative) = better for Vina/gnina affinity-style scores."""
    return pd.to_numeric(series, errors="coerce")


def best_pose_per_compound(
    df: pd.DataFrame,
    smiles_col: str = "canonical_smiles",
    score_col: str = CANONICAL_SCORE_COL,
) -> pd.DataFrame:
    out = df.sort_values(score_col, ascending=True, na_position="last")
    return out.groupby(smiles_col, as_index=False).first()


def percentile_rank(series: pd.Series, higher_is_better: bool = False) -> pd.Series:
    if higher_is_better:
        return series.rank(method="average", pct=True, na_option="bottom") * 100.0
    return (1.0 - series.rank(method="average", pct=True, na_option="bottom")) * 100.0


def ensure_dock_score_column(df: pd.DataFrame, score_col: str | None = None) -> pd.DataFrame:
    """Add canonical dock_score (gnina/Vina)."""
    out = df.copy()
    if score_col is None:
        score_col = pick_col(out.columns, DOCK_SCORE_ALIASES)
    if score_col is None:
        raise ValueError(f"No docking score column found (tried {DOCK_SCORE_ALIASES})")
    out[CANONICAL_SCORE_COL] = coerce_dock_score(out[score_col])
    return out


def docking_status_from_score(scores: pd.Series, status: pd.Series | None = None) -> pd.Series:
    if status is not None:
        return status.astype(str)
    return np.where(scores.notna(), "docked", "missing")
