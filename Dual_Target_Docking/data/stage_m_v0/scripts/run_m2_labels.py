#!/usr/bin/env python3
"""PART B — M2 label validity (margin, continuous, threshold, noise ceiling)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    ANALYSIS,
    SEED,
    TABLES,
    assign_fourclass,
    assign_margin_label,
    auroc,
    directional_metrics,
    load_merged,
    map_strict_to_class,
)

NOISE_SIGMAS = (0.3, 0.5, 0.7)
N_BOOT = 500
CUTOFFS = (5.5, 6.0, 6.5)
ARMS_SENS = ("vina_mean", "rtm_min_z", "heavy_atoms")


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    margin_counts = []
    margin_dir = []
    spearman_rows = []
    thresh_rows = []
    noise_rows = []

    for pair in ("EGFR_HER2", "PIK3CA_mTOR"):
        df = load_merged(pair)
        df["margin_label"] = [
            assign_margin_label(a, b) for a, b in zip(df["pA"], df["pB"])
        ]

        # B1 counts
        vc = df["margin_label"].value_counts().to_dict()
        n_both = int(((df["pA"].notna()) & (df["pB"].notna())).sum())
        for lab in (
            "dual_strict",
            "A_only_strict",
            "B_only_strict",
            "neither_strict",
            "gray",
            "incomplete",
        ):
            margin_counts.append(
                {
                    "pair": pair,
                    "label": lab,
                    "n": int(vc.get(lab, 0)),
                    "n_both_measured": n_both,
                    "frac_of_panel": round(int(vc.get(lab, 0)) / len(df), 3),
                }
            )

        # directional on strict (exclude gray/incomplete)
        strict = df[df["margin_label"].isin(
            ["dual_strict", "A_only_strict", "B_only_strict", "neither_strict"]
        )].copy()
        strict["class"] = strict["margin_label"].map(map_strict_to_class)
        under = {
            c: int((strict["class"] == c).sum()) < 8
            for c in ("dual", "A_only", "B_only")
        }
        for arm in ARMS_SENS:
            m = directional_metrics(strict, arm)
            m.update(
                {
                    "pair": pair,
                    "panel": "margin_strict",
                    "underpowered": any(under.values()),
                    "underpowered_detail": ";".join(
                        f"{k}<8" for k, v in under.items() if v
                    ),
                }
            )
            margin_dir.append(m)

        # B2 continuous Spearman
        for arm in list(ARMS_SENS) + ["MW", "cLogP", "TPSA"]:
            if arm not in df.columns:
                continue
            # overall: both ends measured
            mask = df["pA"].notna() & df["pB"].notna()
            y = df.loc[mask, "min_pchembl"].astype(float)
            # ensure min_pchembl
            if y.isna().all():
                y = pd.Series(
                    [
                        np.nanmin([a, b]) if a is not None and b is not None else np.nan
                        for a, b in zip(df.loc[mask, "pA"], df.loc[mask, "pB"])
                    ],
                    index=y.index,
                )
            s = df.loc[mask, arm].astype(float)
            ok = y.notna() & s.notna()
            rho, pval = spearmanr(s[ok], y[ok]) if ok.sum() >= 5 else (np.nan, np.nan)
            spearman_rows.append(
                {
                    "pair": pair,
                    "scope": "all_both_measured",
                    "arm": arm,
                    "n": int(ok.sum()),
                    "spearman_rho": round(float(rho), 4) if rho == rho else np.nan,
                    "pvalue": float(pval) if pval == pval else np.nan,
                }
            )
            for cls in ("dual", "A_only", "B_only", "neither"):
                sub = df[(df["class"] == cls) & mask]
                if len(sub) < 5:
                    continue
                yy = sub["min_pchembl"].astype(float)
                ss = sub[arm].astype(float)
                ok2 = yy.notna() & ss.notna()
                if ok2.sum() < 5:
                    continue
                rho2, p2 = spearmanr(ss[ok2], yy[ok2])
                spearman_rows.append(
                    {
                        "pair": pair,
                        "scope": f"class_{cls}",
                        "arm": arm,
                        "n": int(ok2.sum()),
                        "spearman_rho": round(float(rho2), 4),
                        "pvalue": float(p2),
                    }
                )

        # B3 threshold sensitivity
        for cut in CUTOFFS:
            lab = [
                assign_fourclass(a, b, cut) for a, b in zip(df["pA"], df["pB"])
            ]
            tmp = df.copy()
            tmp["class"] = lab
            tmp = tmp[tmp["class"].isin(["dual", "A_only", "B_only", "neither"])]
            for arm in ARMS_SENS:
                m = directional_metrics(tmp, arm)
                m.update({"pair": pair, "cutoff": cut})
                thresh_rows.append(m)

        # B4 noise ceiling
        base = df[df["pA"].notna() & df["pB"].notna()].copy()
        pA0 = base["pA"].astype(float).values
        pB0 = base["pB"].astype(float).values
        vina = base["vina_mean"].astype(float).values
        oracle = np.minimum(pA0, pB0)
        for sigma in NOISE_SIGMAS:
            da_v, db_v, da_o, db_o = [], [], [], []
            for _ in range(N_BOOT):
                pA = pA0 + rng.normal(0, sigma, size=len(pA0))
                pB = pB0 + rng.normal(0, sigma, size=len(pB0))
                labels = [assign_fourclass(a, b, 6.0) for a, b in zip(pA, pB)]
                labels = np.asarray(labels)
                # skip incomplete (none here)
                d_idx = np.where(labels == "dual")[0]
                a_idx = np.where(labels == "A_only")[0]
                b_idx = np.where(labels == "B_only")[0]
                if len(d_idx) == 0 or len(a_idx) == 0 or len(b_idx) == 0:
                    continue
                da_v.append(auroc(vina[d_idx], vina[a_idx]))
                db_v.append(auroc(vina[d_idx], vina[b_idx]))
                da_o.append(auroc(oracle[d_idx], oracle[a_idx]))
                db_o.append(auroc(oracle[d_idx], oracle[b_idx]))
            for score_name, da, db in (
                ("vina_mean_true_scores", da_v, db_v),
                ("oracle_min_pchembl", da_o, db_o),
            ):
                da, db = np.asarray(da), np.asarray(db)
                noise_rows.append(
                    {
                        "pair": pair,
                        "sigma": sigma,
                        "score": score_name,
                        "n_boot": len(da),
                        "median_D_vs_A": float(np.median(da)),
                        "p25_D_vs_A": float(np.percentile(da, 25)),
                        "p75_D_vs_A": float(np.percentile(da, 75)),
                        "median_D_vs_B": float(np.median(db)),
                        "p25_D_vs_B": float(np.percentile(db, 25)),
                        "p75_D_vs_B": float(np.percentile(db, 75)),
                        "median_min_DA_DB": float(
                            np.median(np.minimum(da, db))
                        ),
                    }
                )

    pd.DataFrame(margin_counts).to_csv(TABLES / "m2_margin_panel_counts.csv", index=False)
    md_df = pd.DataFrame(margin_dir)
    for c in ("auroc_D_vs_A", "auroc_D_vs_B", "auroc_pooled", "summary_min", "summary_mean"):
        if c in md_df.columns:
            md_df[c] = md_df[c].round(4)
    md_df.to_csv(TABLES / "m2_directional_on_margin.csv", index=False)
    pd.DataFrame(spearman_rows).to_csv(TABLES / "m2_continuous_spearman.csv", index=False)
    th = pd.DataFrame(thresh_rows)
    for c in ("auroc_D_vs_A", "auroc_D_vs_B", "auroc_pooled", "summary_min", "summary_mean"):
        if c in th.columns:
            th[c] = th[c].round(4)
    th.to_csv(TABLES / "m2_threshold_sensitivity.csv", index=False)
    nz = pd.DataFrame(noise_rows)
    for c in nz.columns:
        if c.startswith(("median", "p25", "p75")):
            nz[c] = nz[c].round(4)
    nz.to_csv(TABLES / "m2_noise_ceiling.csv", index=False)

    # Gate logic for doc
    # flip check: D/B sign (relative to 0.5) flip across cutoffs for vina_mean
    flip_notes = []
    for pair in ("EGFR_HER2", "PIK3CA_mTOR"):
        sub = th[(th.pair == pair) & (th.arm == "vina_mean")]
        signs = [1 if r.auroc_D_vs_B > 0.5 else -1 for _, r in sub.iterrows()]
        flipped = len(set(signs)) > 1
        vs_vol = []
        for cut in CUTOFFS:
            v = th[(th.pair == pair) & (th.arm == "vina_mean") & (th.cutoff == cut)].iloc[0]
            h = th[(th.pair == pair) & (th.arm == "heavy_atoms") & (th.cutoff == cut)].iloc[0]
            vs_vol.append(v.summary_min > h.summary_min)
        flip_notes.append(
            {
                "pair": pair,
                "DB_sign_flips_across_cutoff": flipped,
                "vina_beats_heavy_at_any_cutoff": any(vs_vol),
                "vina_beats_heavy_all_cutoffs": all(vs_vol),
            }
        )

    oracle_ok = True
    oracle_detail = []
    for pair in ("EGFR_HER2", "PIK3CA_mTOR"):
        row = nz[
            (nz.pair == pair)
            & (nz.sigma == 0.5)
            & (nz.score == "oracle_min_pchembl")
        ].iloc[0]
        ok_dir = (row.median_D_vs_A >= 0.65) or (row.median_D_vs_B >= 0.65)
        oracle_detail.append(
            f"{pair}: oracle@σ0.5 D/A={row.median_D_vs_A:.3f} D/B={row.median_D_vs_B:.3f}"
        )
        if not ok_dir:
            oracle_ok = False

    # margin same-sign as main for vina
    m1 = pd.read_csv(TABLES / "m1_directional_auroc.csv")
    same_sign = True
    for pair in ("EGFR_HER2", "PIK3CA_mTOR"):
        main = m1[(m1.pair == pair) & (m1.subset == "all") & (m1.arm == "vina_mean")].iloc[0]
        mar = md_df[(md_df.pair == pair) & (md_df.arm == "vina_mean")].iloc[0]
        for col in ("auroc_D_vs_A", "auroc_D_vs_B"):
            if np.sign(main[col] - 0.5) != np.sign(mar[col] - 0.5) and not (
                abs(main[col] - 0.5) < 0.05 and abs(mar[col] - 0.5) < 0.05
            ):
                # allow near-chance wobble
                if abs(main[col] - 0.5) >= 0.05 or abs(mar[col] - 0.5) >= 0.05:
                    if np.sign(main[col] - 0.5) != np.sign(mar[col] - 0.5):
                        same_sign = False

    # M2 gate
    gray_frac_eh = [
        r for r in margin_counts if r["pair"] == "EGFR_HER2" and r["label"] == "gray"
    ][0]["frac_of_panel"]
    any_flip = any(x["DB_sign_flips_across_cutoff"] for x in flip_notes)
    underpowered_margin = bool(md_df["underpowered"].any()) if "underpowered" in md_df else False

    if any_flip or not oracle_ok:
        m2_gate = "No-Go"
    elif underpowered_margin or not same_sign or gray_frac_eh >= 0.25:
        m2_gate = "Weak"
    else:
        m2_gate = "Go"
    # refine: if oracle ok and same sign and no flip → Go even with some gray;
    # Weak if margin n small or edge ceiling
    if m2_gate == "No-Go" and oracle_ok and not any_flip:
        m2_gate = "Weak"
    if underpowered_margin and m2_gate == "Go":
        m2_gate = "Weak"

    # write gate sidecar for verdict
    pd.DataFrame(
        [{"m2_gate": m2_gate, "same_sign_margin": same_sign, "oracle_ok_sigma0.5": oracle_ok}]
        + flip_notes
    ).to_csv(TABLES / "m2_gate_summary.csv", index=False)

    md = []
    md.append("# M2 — Label validity\n\n")
    md.append("## Rules (frozen)\n\n")
    md.append(
        "- `dual_strict`: both pChEMBL ≥ 6.5\n"
        "- `A_only_strict`: A≥6.5 and B≤5.5 (measured)\n"
        "- `B_only_strict`: B≥6.5 and A≤5.5\n"
        "- `neither_strict`: both ≤5.5\n"
        "- `gray`: both measured, else → excluded from strict directional analysis\n"
        "- Unmeasured end ≠ negative (incomplete excluded from re-binarization)\n\n"
    )
    md.append("## Margin counts\n\n")
    md.append(pd.DataFrame(margin_counts).to_markdown(index=False) + "\n\n")
    md.append("## Directional on strict margin\n\n")
    show = md_df[["pair", "arm", "n_dual", "n_A_only", "n_B_only",
                  "auroc_D_vs_A", "auroc_D_vs_B", "underpowered"]]
    md.append(show.to_markdown(index=False) + "\n\n")
    md.append("## Threshold sensitivity (vina_mean)\n\n")
    show2 = th[th.arm == "vina_mean"][
        ["pair", "cutoff", "n_dual", "n_A_only", "n_B_only", "auroc_D_vs_A", "auroc_D_vs_B"]
    ]
    md.append(show2.to_markdown(index=False) + "\n\n")
    md.append("## Noise / oracle ceiling (σ=0.5 highlight)\n\n")
    show3 = nz[nz.sigma == 0.5][
        ["pair", "score", "median_D_vs_A", "median_D_vs_B", "median_min_DA_DB"]
    ]
    md.append(show3.to_markdown(index=False) + "\n\n")
    md.append("## Gate\n\n")
    md.append(f"**M2 = {m2_gate}**\n\n")
    md.append(
        f"- Margin same-sign vs main panel (vina): {same_sign}\n"
        f"- Oracle distinguishable @σ=0.5 (≥1 direction median≥0.65): {oracle_ok}\n"
        f"- Details: {'; '.join(oracle_detail)}\n"
        f"- Cutoff flips D/B sign: {any_flip}\n"
        f"- Margin underpowered flag: {underpowered_margin}\n"
    )
    (ANALYSIS / "M2_LABELS.md").write_text("".join(md))
    print("M2 gate", m2_gate)


if __name__ == "__main__":
    main()
