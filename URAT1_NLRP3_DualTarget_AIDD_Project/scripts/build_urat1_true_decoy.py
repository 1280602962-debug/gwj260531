#!/usr/bin/env python3
"""
Build a URAT1 TrueDecoy enrichment benchmark (Gu et al. Nat Mach Intell 2025 framing).

TrueDecoy = property-matched putative inactives (harder enrichment test).
RandomDecoy companion = same-size random sample from the same pool (easier / VS-like).

Default actives: curated URAT1 with pActivity >= 6 (potent subset).
Default decoy pool: distill subset D (unlabeled diversity negatives).

Matching descriptors (DUD-E-inspired windows):
  MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds

Outputs under data/benchmarks/urat1_true_decoy/:
  actives.csv, true_decoys.csv, true_decoy_benchmark.csv,
  random_decoys.csv, random_decoy_benchmark.csv (optional),
  matching_assignments.csv, summary.json, README.md

Example:
  python3 scripts/build_urat1_true_decoy.py
  python3 scripts/build_urat1_true_decoy.py --pactivity-min 6 --ratio 30 --seed 42
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from utils_ml import canonicalize, murcko_scaffold

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ACTIVES = PROJECT_ROOT / "data" / "processed" / "urat1_curated.csv"
DEFAULT_POOL = PROJECT_ROOT / "data" / "distill" / "distill_subset_d.csv"
DEFAULT_OUT = PROJECT_ROOT / "data" / "benchmarks" / "urat1_true_decoy"

EXCLUDE_EXTRA = [
    PROJECT_ROOT / "data" / "distill" / "distill_subset_a.csv",
    PROJECT_ROOT / "data" / "benchmarks" / "literature_benchmarks.csv",
]

# DUD-E-inspired absolute windows (see Huang/Shoichet DUD family; Gu TrueDecoy spirit)
PROP_WINDOWS = {
    "MolWt": 40.0,
    "MolLogP": 1.0,
    "TPSA": 25.0,
    "NumHDonors": 1.5,  # allow ±1 integer donors
    "NumHAcceptors": 2.5,
    "NumRotatableBonds": 2.5,
}

PROP_COLS = list(PROP_WINDOWS.keys())


def _mol_props(smiles: str) -> dict | None:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return {
        "MolWt": float(Descriptors.MolWt(mol)),
        "MolLogP": float(Descriptors.MolLogP(mol)),
        "TPSA": float(Descriptors.TPSA(mol)),
        "NumHDonors": float(Descriptors.NumHDonors(mol)),
        "NumHAcceptors": float(Descriptors.NumHAcceptors(mol)),
        "NumRotatableBonds": float(Descriptors.NumRotatableBonds(mol)),
    }


def _morgan_fp(smiles: str, n_bits: int = 2048):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    return gen.GetFingerprint(mol)


def annotate_props(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    rows = []
    for smi in df[smiles_col].astype(str):
        canon = canonicalize(smi)
        if not canon:
            continue
        props = _mol_props(canon)
        if props is None:
            continue
        rows.append({"canonical_smiles": canon, **props})
    out = pd.DataFrame(rows).drop_duplicates("canonical_smiles").reset_index(drop=True)
    return out


def load_exclude_smiles(extra_paths: list[Path], active_smiles: set[str]) -> set[str]:
    excluded = set(active_smiles)
    for path in extra_paths:
        if not path.exists():
            continue
        df = pd.read_csv(path, low_memory=False)
        for col in ("canonical_smiles", "smiles", "SMILES", "Smiles"):
            if col not in df.columns:
                continue
            for smi in df[col].dropna().astype(str):
                c = canonicalize(smi)
                if c:
                    excluded.add(c)
            break
    return excluded


def max_tanimoto_to_actives(decoy_smi: str, active_fps: list) -> float:
    from rdkit import DataStructs

    fp = _morgan_fp(decoy_smi)
    if fp is None or not active_fps:
        return 0.0
    return float(max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)))


def property_distance(active_row: pd.Series, decoy_row: pd.Series) -> float:
    """Normalized L2 distance across matched properties."""
    diffs = []
    for col, window in PROP_WINDOWS.items():
        diffs.append((float(active_row[col]) - float(decoy_row[col])) / window)
    return float(np.sqrt(np.mean(np.square(diffs))))


def within_windows(
    active_row: pd.Series,
    decoy_row: pd.Series,
    windows: dict[str, float] | None = None,
) -> bool:
    wins = windows or PROP_WINDOWS
    for col, window in wins.items():
        if abs(float(active_row[col]) - float(decoy_row[col])) > window:
            return False
    return True


def _scaled_windows(scale: float) -> dict[str, float]:
    return {k: v * scale for k, v in PROP_WINDOWS.items()}


def match_true_decoys(
    actives: pd.DataFrame,
    pool: pd.DataFrame,
    *,
    ratio: int,
    seed: int,
    max_sim_to_active: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Round-robin property matching without decoy reuse.

    Round-robin avoids early actives monopolizing the pool (greedy-per-active
    left many potent ligands with zero decoys). A second pass with 1.5× windows
    fills remaining quota for hard-to-match actives.
    """
    rng = np.random.default_rng(seed)
    active_order = actives.sample(frac=1.0, random_state=seed).reset_index(drop=True)

    active_fps = []
    for smi in actives["canonical_smiles"]:
        fp = _morgan_fp(smi)
        if fp is not None:
            active_fps.append(fp)

    available = pool.copy().reset_index(drop=True)
    available["_used"] = False

    # Pre-filter by chemical similarity to actives (avoid near-analogs / latent actives)
    if max_sim_to_active < 1.0 and active_fps:
        sims = [max_tanimoto_to_actives(smi, active_fps) for smi in available["canonical_smiles"]]
        available["max_tc_active"] = sims
        available = available[available["max_tc_active"] <= max_sim_to_active].reset_index(drop=True)
    else:
        available["max_tc_active"] = np.nan

    # Precompute property matrix for speed
    decoy_props = available[PROP_COLS].to_numpy(dtype=float)
    used = np.zeros(len(available), dtype=bool)
    counts = {smi: 0 for smi in active_order["canonical_smiles"]}
    assignments: list[dict] = []

    def _best_match(active_row: pd.Series, windows: dict[str, float]) -> int | None:
        a = active_row[PROP_COLS].to_numpy(dtype=float)
        wins = np.array([windows[c] for c in PROP_COLS], dtype=float)
        free = ~used
        if not free.any():
            return None
        delta = np.abs(decoy_props[free] - a)
        ok = np.all(delta <= wins, axis=1)
        if not ok.any():
            return None
        free_idx = np.flatnonzero(free)
        cand = free_idx[ok]
        # normalized L2; break ties with RNG
        norm = delta[ok] / wins
        dist = np.sqrt(np.mean(norm * norm, axis=1))
        # among top-20 nearest, pick randomly for diversity
        order = np.argsort(dist)
        head = order[: min(20, len(order))]
        pick = int(rng.choice(head))
        return int(cand[pick])

    def _assign(active_row: pd.Series, idx: int, windows: dict[str, float], pass_name: str) -> None:
        used[idx] = True
        smi = active_row["canonical_smiles"]
        counts[smi] += 1
        a = active_row[PROP_COLS].to_numpy(dtype=float)
        d = decoy_props[idx]
        wins = np.array([windows[c] for c in PROP_COLS], dtype=float)
        norm = (a - d) / wins
        dist = float(np.sqrt(np.mean(norm * norm)))
        tc = available.at[idx, "max_tc_active"]
        assignments.append(
            {
                "active_smiles": smi,
                "decoy_smiles": available.at[idx, "canonical_smiles"],
                "property_distance": dist,
                "max_tc_active": float(tc) if pd.notna(tc) else None,
                "window_pass": pass_name,
            }
        )

    # Pass 1: strict windows, round-robin up to ratio
    for _round in range(ratio):
        made = 0
        for _, active in active_order.iterrows():
            smi = active["canonical_smiles"]
            if counts[smi] >= ratio:
                continue
            idx = _best_match(active, PROP_WINDOWS)
            if idx is None:
                continue
            _assign(active, idx, PROP_WINDOWS, "strict")
            made += 1
        if made == 0:
            break

    # Pass 2: 1.5× windows — top-up actives still below target (esp. zero-match)
    relaxed = _scaled_windows(1.5)
    for _round in range(ratio):
        made = 0
        for _, active in active_order.iterrows():
            smi = active["canonical_smiles"]
            if counts[smi] >= ratio:
                continue
            idx = _best_match(active, relaxed)
            if idx is None:
                continue
            _assign(active, idx, relaxed, "relaxed_1.5x")
            made += 1
        if made == 0:
            break

    selected_idx = [i for i, u in enumerate(used) if u]
    decoys = available.loc[selected_idx].drop(columns=["_used"]).reset_index(drop=True)
    assign_df = pd.DataFrame(assignments)
    return decoys, assign_df


def prop_stats(df: pd.DataFrame) -> dict:
    out = {}
    for col in PROP_COLS:
        s = df[col]
        out[col] = {
            "mean": float(s.mean()),
            "std": float(s.std(ddof=0)),
            "median": float(s.median()),
            "min": float(s.min()),
            "max": float(s.max()),
        }
    return out


def write_readme(
    path: Path,
    *,
    n_actives: int,
    n_true: int,
    n_random: int,
    pmin: float,
    ratio_target: int,
    ratio_achieved: float,
    max_sim: float,
    seed: int,
    n_covered: int,
    decoys_per_active: dict,
) -> None:
    text = f"""# URAT1 TrueDecoy Benchmark

Property-matched enrichment set for URAT1 docking / scoring validation, framed after
**Gu et al., *Nat. Mach. Intell.* 2025** (TrueDecoy vs RandomDecoy).

## Design

| Item | Choice |
|------|--------|
| Actives | Curated URAT1, `pActivity >= {pmin}` (n = {n_actives}) |
| Decoy pool | Distill subset D (unlabeled diversity negatives) |
| Matching | Round-robin MW / logP / TPSA / HBD / HBA / rotatable bonds (DUD-E-inspired windows); 1.5× relaxed top-up |
| Near-analog filter | Max Morgan TC to any active ≤ {max_sim} |
| Target ratio | 1 : {ratio_target} (active : decoy) |
| Achieved TrueDecoy ratio | 1 : {ratio_achieved:.2f} (n_decoy = {n_true}) |
| Actives with ≥1 decoy | {n_covered} / {n_actives} |
| Decoys / covered active | min {decoys_per_active.get('min', 0)}, median {decoys_per_active.get('median', 0):.1f}, mean {decoys_per_active.get('mean', 0):.1f}, max {decoys_per_active.get('max', 0)} |
| RandomDecoy companion | Same size random sample from pool (n = {n_random}) |
| Seed | {seed} |

### Property windows (strict)

| Descriptor | Window |
|------------|--------|
| MolWt | ±40 Da |
| MolLogP | ±1.0 |
| TPSA | ±25 Å² |
| NumHDonors | ±1 |
| NumHAcceptors | ±2 |
| NumRotatableBonds | ±2 |

## Files

| File | Content |
|------|---------|
| `actives.csv` | Potent URAT1 actives + descriptors |
| `true_decoys.csv` | Property-matched decoys |
| `true_decoy_benchmark.csv` | Combined set with `label` (1=active, 0=decoy) and `decoy_class=true` |
| `random_decoys.csv` | Random unmatched decoys (same n when possible) |
| `random_decoy_benchmark.csv` | Combined RandomDecoy set |
| `matching_assignments.csv` | Active→decoy pairs with property distance / window pass |
| `summary.json` | Counts, windows, property stats |

## Usage notes

- **TrueDecoy** is the harder test: physics-based docking/scoring should be judged here.
- **RandomDecoy** mirrors easier / library-like enrichment (subset A vs raw D is closer to this).
- Do **not** train ML models on these decoy labels; they are putative inactives for enrichment only.
- Pool size (~8k subset D) limits the achievable ratio below the 1:30–50 ideal; rebuild with a larger `--pool` if needed.

## Rebuild

```bash
python3 scripts/build_urat1_true_decoy.py \\
  --pactivity-min {pmin} \\
  --ratio {ratio_target} \\
  --max-sim-to-active {max_sim} \\
  --seed {seed}
```
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build URAT1 TrueDecoy / RandomDecoy benchmarks")
    parser.add_argument("--actives", type=Path, default=DEFAULT_ACTIVES)
    parser.add_argument("--pool", type=Path, default=DEFAULT_POOL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--pactivity-min", type=float, default=6.0)
    parser.add_argument("--ratio", type=int, default=30, help="Target decoys per active")
    parser.add_argument("--max-sim-to-active", type=float, default=0.5, help="Max Morgan TC to any active")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-random", action="store_true", help="Skip RandomDecoy companion")
    args = parser.parse_args()

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Actives ---
    raw_act = pd.read_csv(args.actives, low_memory=False)
    smi_col = "canonical_smiles" if "canonical_smiles" in raw_act.columns else "Smiles"
    raw_act = raw_act[pd.to_numeric(raw_act["pActivity"], errors="coerce") >= args.pactivity_min].copy()
    act_props = annotate_props(raw_act.rename(columns={smi_col: "canonical_smiles"}))
    # keep pActivity / ids
    meta = raw_act.copy()
    meta["canonical_smiles"] = meta[smi_col].map(canonicalize)
    keep_cols = [c for c in ("pActivity", "molecule_chembl_id") if c in meta.columns]
    meta = meta[["canonical_smiles"] + keep_cols].drop_duplicates("canonical_smiles")
    actives = act_props.merge(meta, on="canonical_smiles", how="left")
    actives["scaffold"] = actives["canonical_smiles"].map(murcko_scaffold)
    actives["label"] = 1
    actives["set_role"] = "active"

    active_set = set(actives["canonical_smiles"])
    excluded = load_exclude_smiles(EXCLUDE_EXTRA, active_set)

    # --- Pool ---
    raw_pool = pd.read_csv(args.pool, low_memory=False)
    pool_smi_col = "canonical_smiles" if "canonical_smiles" in raw_pool.columns else "Smiles"
    pool_props = annotate_props(raw_pool.rename(columns={pool_smi_col: "canonical_smiles"}))
    pool_props = pool_props[~pool_props["canonical_smiles"].isin(excluded)].reset_index(drop=True)

    print(f"Actives (pActivity>={args.pactivity_min}): {len(actives)}")
    print(f"Decoy pool after exclude: {len(pool_props)}")

    true_decoys, assignments = match_true_decoys(
        actives,
        pool_props,
        ratio=args.ratio,
        seed=args.seed,
        max_sim_to_active=args.max_sim_to_active,
    )
    true_decoys = true_decoys.copy()
    true_decoys["scaffold"] = true_decoys["canonical_smiles"].map(murcko_scaffold)
    true_decoys["label"] = 0
    true_decoys["set_role"] = "true_decoy"
    true_decoys["decoy_class"] = "true"

    ratio_achieved = len(true_decoys) / max(len(actives), 1)
    print(f"TrueDecoys selected: {len(true_decoys)} (ratio 1:{ratio_achieved:.2f})")
    print(f"Assignments (active-decoy pairs): {len(assignments)}")

    # Combined TrueDecoy benchmark
    true_bench = pd.concat(
        [
            actives.assign(decoy_class="true"),
            true_decoys,
        ],
        ignore_index=True,
    )

    # RandomDecoy companion: same n from unused pool members preferred
    n_random = 0
    if not args.skip_random:
        used = set(true_decoys["canonical_smiles"])
        unused = pool_props[~pool_props["canonical_smiles"].isin(used)].copy()
        n_want = len(true_decoys)
        if len(unused) < n_want:
            unused = pool_props.copy()
        random_decoys = unused.sample(n=min(n_want, len(unused)), random_state=args.seed).reset_index(drop=True)
        random_decoys["scaffold"] = random_decoys["canonical_smiles"].map(murcko_scaffold)
        random_decoys["label"] = 0
        random_decoys["set_role"] = "random_decoy"
        random_decoys["decoy_class"] = "random"
        n_random = len(random_decoys)
        random_bench = pd.concat(
            [actives.assign(decoy_class="random"), random_decoys],
            ignore_index=True,
        )
        random_decoys.to_csv(out_dir / "random_decoys.csv", index=False)
        random_bench.to_csv(out_dir / "random_decoy_benchmark.csv", index=False)
        print(f"RandomDecoys selected: {n_random}")

    # Write outputs
    actives.to_csv(out_dir / "actives.csv", index=False)
    true_decoys.to_csv(out_dir / "true_decoys.csv", index=False)
    true_bench.to_csv(out_dir / "true_decoy_benchmark.csv", index=False)
    assignments.to_csv(out_dir / "matching_assignments.csv", index=False)

    summary = {
        "target": "URAT1",
        "framing": "Gu et al. Nat Mach Intell 2025 TrueDecoy vs RandomDecoy",
        "pactivity_min": args.pactivity_min,
        "n_actives": int(len(actives)),
        "n_true_decoys": int(len(true_decoys)),
        "n_random_decoys": int(n_random),
        "target_ratio": args.ratio,
        "achieved_true_ratio": float(ratio_achieved),
        "max_sim_to_active": args.max_sim_to_active,
        "seed": args.seed,
        "property_windows": PROP_WINDOWS,
        "relaxed_window_scale": 1.5,
        "matching_algorithm": "round_robin_property_match_then_relaxed_1.5x",
        "pool_source": str(args.pool.relative_to(PROJECT_ROOT)),
        "actives_source": str(args.actives.relative_to(PROJECT_ROOT)),
        "n_pool_after_exclude": int(len(pool_props)),
        "n_assignment_pairs": int(len(assignments)),
        "property_stats": {
            "actives": prop_stats(actives),
            "true_decoys": prop_stats(true_decoys) if len(true_decoys) else {},
        },
        "mean_property_distance": float(assignments["property_distance"].mean())
        if len(assignments)
        else None,
        "actives_with_ge1_decoy": int(assignments["active_smiles"].nunique()) if len(assignments) else 0,
        "decoys_per_active": {
            "min": int(assignments.groupby("active_smiles").size().min()) if len(assignments) else 0,
            "median": float(assignments.groupby("active_smiles").size().median()) if len(assignments) else 0,
            "mean": float(assignments.groupby("active_smiles").size().mean()) if len(assignments) else 0,
            "max": int(assignments.groupby("active_smiles").size().max()) if len(assignments) else 0,
        },
        "window_pass_counts": (
            assignments["window_pass"].value_counts().to_dict() if "window_pass" in assignments.columns else {}
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    write_readme(
        out_dir / "README.md",
        n_actives=len(actives),
        n_true=len(true_decoys),
        n_random=n_random,
        pmin=args.pactivity_min,
        ratio_target=args.ratio,
        ratio_achieved=ratio_achieved,
        max_sim=args.max_sim_to_active,
        seed=args.seed,
        n_covered=int(summary["actives_with_ge1_decoy"]),
        decoys_per_active=summary["decoys_per_active"],
    )

    print(f"Wrote benchmark to {out_dir}")


if __name__ == "__main__":
    main()
