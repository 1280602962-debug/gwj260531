#!/usr/bin/env python3
"""C4 — Pre-registered JNK1/2/3 IC50 / selectivity-index analysis.

LOCKED BEFORE UNBLINDING.
Fill `results/assay/ic50_raw.csv` (template created if missing), then re-run.

Primary endpoint (RQ-A):
  any_active = IC50_uM <= PRIMARY_IC50_UM on ANY of JNK1/JNK2/JNK3
  for at least one of {690, 2157}.

Secondary endpoint (RQ-B; exploratory, underpowered):
  JNK1 preference only if SI_J2 = IC50_JNK2/IC50_JNK1 >= SI_THRESHOLD
  AND SI_J3 = IC50_JNK3/IC50_JNK1 >= SI_THRESHOLD
  (missing off-isoform IC50 → preference NOT claimed).

Forbidden post-hoc:
  Do not redefine SI using MD hinge occupancy.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ASSAY_DIR = ROOT / "results" / "assay"
OUT_DIR = ROOT / "results" / "assay_analysis"
ASSAY_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

RAW = ASSAY_DIR / "ic50_raw.csv"
TEMPLATE = ASSAY_DIR / "ic50_raw_TEMPLATE.csv"

# --- LOCKED constants (do not edit after unblinding without version bump) ---
ANALYSIS_VERSION = "c4_v2_locked_2026-07-16_purchase_690_2231"
PRIMARY_IC50_UM = 10.0
SI_THRESHOLD = 3.0
NEW_COMPOUNDS = ("690", "2231")
CONTROLS = ("E1", "CC-90001")
ISOFORMS = ("JNK1", "JNK2", "JNK3")
# Purchase change note: 2157 dropped from wet-lab panel; 2231 added as MD-bias hypothesis test.
PREVIOUS_NEW_COMPOUNDS = ("690", "2157")


def write_template():
    rows = []
    for cid in list(NEW_COMPOUNDS) + list(CONTROLS):
        for iso in ISOFORMS:
            rows.append(
                {
                    "compound_id": cid,
                    "isoform": iso,
                    "ic50_uM": "",
                    "ic50_nM": "",
                    "pct_inh_10uM": "",
                    "n_replicates": "",
                    "assay_date": "",
                    "notes": "",
                }
            )
    pd.DataFrame(rows).to_csv(TEMPLATE, index=False)
    # Refresh RAW layout when purchase set changes and file is still empty of numeric data
    if not RAW.exists():
        pd.DataFrame(rows).to_csv(RAW, index=False)
    else:
        existing = pd.read_csv(RAW)
        has_numeric = False
        for col in ("ic50_uM", "ic50_nM"):
            if col in existing.columns:
                vals = existing[col].astype(str).str.strip()
                has_numeric = has_numeric or bool((vals != "").any() and vals.str.lower().ne("nan").any())
        ids = set(existing.get("compound_id", pd.Series(dtype=str)).astype(str))
        if (not has_numeric) and (ids != set(list(NEW_COMPOUNDS) + list(CONTROLS))):
            pd.DataFrame(rows).to_csv(RAW, index=False)


def _to_float(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "nd", ">"}:
        return np.nan
    # support censored like >50
    if s.startswith(">"):
        return np.nan  # treat as inactive / non-finite for SI; handled via flag
    return float(s)


def load_raw() -> pd.DataFrame:
    write_template()
    df = pd.read_csv(RAW)
    # Prefer explicit uM; else convert nM
    ic50 = []
    for _, r in df.iterrows():
        u = _to_float(r.get("ic50_uM"))
        n = _to_float(r.get("ic50_nM"))
        if not np.isnan(u):
            ic50.append(u)
        elif not np.isnan(n):
            ic50.append(n / 1000.0)
        else:
            ic50.append(np.nan)
    df["ic50_uM_resolved"] = ic50
    return df


def wide_table(df: pd.DataFrame) -> pd.DataFrame:
    ids = list(NEW_COMPOUNDS) + list(CONTROLS)
    # Always emit one row per locked compound, even before assay data arrive.
    base = pd.DataFrame({"compound_id": ids})
    if df.empty or df["ic50_uM_resolved"].notna().sum() == 0:
        for iso in ISOFORMS:
            base[iso] = np.nan
        return base

    w = df.pivot_table(
        index="compound_id",
        columns="isoform",
        values="ic50_uM_resolved",
        aggfunc="mean",
    )
    if isinstance(w.columns, pd.MultiIndex):
        w.columns = w.columns.get_level_values(-1)
    w = w.reset_index()
    for iso in ISOFORMS:
        if iso not in w.columns:
            w[iso] = np.nan
    w["compound_id"] = w["compound_id"].astype(str)
    merged = base.merge(w[["compound_id", *ISOFORMS]], on="compound_id", how="left")
    return merged


def analyze(w: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in w.iterrows():
        cid = str(r["compound_id"])
        j1, j2, j3 = r["JNK1"], r["JNK2"], r["JNK3"]
        any_active = any(
            (not np.isnan(x)) and x <= PRIMARY_IC50_UM for x in (j1, j2, j3)
        )
        si_j2 = (j2 / j1) if (not np.isnan(j1) and not np.isnan(j2) and j1 > 0) else np.nan
        si_j3 = (j3 / j1) if (not np.isnan(j1) and not np.isnan(j3) and j1 > 0) else np.nan
        preference = (
            (not np.isnan(si_j2))
            and (not np.isnan(si_j3))
            and si_j2 >= SI_THRESHOLD
            and si_j3 >= SI_THRESHOLD
        )
        rows.append(
            {
                "compound_id": cid,
                "role": "new" if cid in NEW_COMPOUNDS else "control",
                "IC50_JNK1_uM": j1,
                "IC50_JNK2_uM": j2,
                "IC50_JNK3_uM": j3,
                "pIC50_JNK1": -np.log10(j1 * 1e-6) if not np.isnan(j1) and j1 > 0 else np.nan,
                "SI_J2_over_J1": si_j2,
                "SI_J3_over_J1": si_j3,
                "primary_any_active_le_10uM": any_active,
                "secondary_jnk1_preference_SI_ge_3": preference,
            }
        )
    return pd.DataFrame(rows)


def main():
    df = load_raw()
    filled = df["ic50_uM_resolved"].notna().sum()
    w = wide_table(df)
    out = analyze(w)

    out.to_csv(OUT_DIR / "c4_ic50_si_table.csv", index=False)
    meta = {
        "analysis_version": ANALYSIS_VERSION,
        "primary_ic50_uM": PRIMARY_IC50_UM,
        "si_threshold": SI_THRESHOLD,
        "n_numeric_ic50_cells": int(filled),
        "status": "WAITING_FOR_ASSAY_DATA" if filled == 0 else "COMPUTED",
        "rq_a_rule": "≥1 of {690,2231} with any isoform IC50 ≤ 10 µM",
        "rq_b_rule": "SI_J2≥3 AND SI_J3≥3 using IC50 ratios; preference hypothesis prioritized on 2231; else no preference claim",
        "purchase_set": list(NEW_COMPOUNDS),
        "supersedes": "c4_v1_locked_2026-07-15 (had 690+2157)",
        "forbidden": ["MD hinge as selectivity proof", "post-hoc SI threshold changes", "claiming kinome selectivity"],
    }
    (OUT_DIR / "c4_analysis_lock.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # RQ-A aggregate (only meaningful when data present)
    new = out[out["compound_id"].isin(NEW_COMPOUNDS)]
    rq_a = bool(new["primary_any_active_le_10uM"].any()) if filled else None
    verdict = {
        "RQ_A_family_enrichment_pass": rq_a,
        "RQ_B_preference_compounds": list(
            new.loc[new["secondary_jnk1_preference_SI_ge_3"], "compound_id"]
        )
        if filled
        else [],
    }
    (OUT_DIR / "c4_endpoint_verdict.json").write_text(json.dumps(verdict, indent=2), encoding="utf-8")

    md = [
        f"# C4 Pre-registered IC50 / SI Analysis (`{ANALYSIS_VERSION}`)",
        "",
        "## Locked rules",
        f"- Primary (RQ-A): any isoform IC50 ≤ **{PRIMARY_IC50_UM} µM** for ≥1 of {{690, 2157}}",
        f"- Secondary (RQ-B): JNK1 preference only if SI_J2 ≥ **{SI_THRESHOLD}** **and** SI_J3 ≥ **{SI_THRESHOLD}**",
        "- Controls: E1 (expect JNK1-leaning direction), CC-90001 (multi-isoform activity)",
        "",
        f"**Data status:** {meta['status']} (numeric IC50 cells = {filled})",
        "",
        "Fill `results/assay/ic50_raw.csv` then re-run this script.",
        "",
        "## Current table",
        "",
        out.to_markdown(index=False),
        "",
    ]
    (OUT_DIR / "C4_PREREGISTERED_ANALYSIS.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(meta, indent=2))
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
