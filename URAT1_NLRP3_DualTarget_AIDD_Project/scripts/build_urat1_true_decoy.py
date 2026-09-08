#!/usr/bin/env python3
"""
Build a URAT1 TrueDecoy enrichment benchmark (Gu et al. Nat Mach Intell 2025 framing).

TrueDecoy negatives =
  (1) experimental weak/inactives from curated URAT1 (pActivity < inactive_max), plus
  (2) property-matched unlabeled molecules from a large commercial-like pool
  to reach target active:decoy ratio (default 1:10).

RandomDecoy companion (Gu-style library random draw) =
  same target count as TrueDecoy negatives, sampled ONLY from pool molecules
  NOT used in TrueDecoy (zero SMILES overlap). No fallback into TrueDecoy set.

Default actives: curated URAT1 with pActivity >= 6.
Default --pool: distill subset D (small local fallback only).
  Production protocol-selection set uses a taosu-prefiltered ~60k pool
  (see data/benchmarks/urat1_true_decoy/README.md); do not treat subset D
  as the official VS benchmark pool.

Matching descriptors (DUD-E-inspired windows):
  MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds

Example:
  python3 scripts/build_urat1_true_decoy.py --pool taosu_pool_prefiltered.csv \\
    --ratio 10 --inactive-pactivity-max 5 --seed 42
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def canonicalize(smiles: str) -> str | None:
    if smiles is None or (isinstance(smiles, float) and np.isnan(smiles)):
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def murcko_scaffold(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return None


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
    n_true_exp: int,
    n_true_matched: int,
    n_random: int,
    pmin: float,
    inactive_max: float,
    ratio_target: int,
    ratio_achieved: float,
    ratio_random: float,
    max_sim: float,
    seed: int,
    n_covered: int,
    decoys_per_active: dict,
    overlap: int,
) -> None:
    text = f"""# URAT1 TrueDecoy / RandomDecoy Benchmark

Framed after **Gu et al., *Nat. Mach. Intell.* 2025** (hard negatives for protocol
selection + library-random companion). This is a **single-target URAT1** set — not a
replica of Gu's multi-target BindingDB TrueDecoy.

## Design

| Item | Choice |
|------|--------|
| Actives | Curated URAT1, `pActivity >= {pmin}` (n = {n_actives}) |
| TrueDecoy negatives | (1) experimental weak/inactives `pActivity < {inactive_max}` (n = {n_true_exp}); (2) property-matched unlabeled from commercial/library pool (n = {n_true_matched}) |
| Matching | Round-robin MW / logP / TPSA / HBD / HBA / rotatable bonds; 1.5× relaxed top-up if needed |
| Near-analog filter (matched only) | Max Morgan TC to any active ≤ {max_sim} |
| Target ratio | 1 : {ratio_target} (active : TrueDecoy negative) |
| Achieved TrueDecoy ratio | 1 : {ratio_achieved:.2f} (n_decoy = {n_true}) |
| Actives with ≥1 matched decoy | {n_covered} / {n_actives} |
| Matched decoys / covered active | min {decoys_per_active.get('min', 0)}, median {decoys_per_active.get('median', 0):.1f}, mean {decoys_per_active.get('mean', 0):.1f}, max {decoys_per_active.get('max', 0)} |
| RandomDecoy | Gu-style random draw from **remaining** pool only; target n = TrueDecoy negatives; **zero SMILES overlap** with TrueDecoy |
| Achieved RandomDecoy ratio | 1 : {ratio_random:.2f} (n_decoy = {n_random}); True∩Random decoy overlap = {overlap} |
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
| `experimental_inactives.csv` | Curated weak/inactives used as TrueDecoy negatives |
| `true_decoys.csv` | All TrueDecoy negatives (`decoy_source` = experimental_inactive \\| property_matched) |
| `true_decoy_benchmark.csv` | Actives + TrueDecoy negatives (`label` 1/0) |
| `random_decoys.csv` | Non-overlapping random library decoys |
| `random_decoy_benchmark.csv` | Actives + RandomDecoy negatives |
| `matching_assignments.csv` | Active→property-matched decoy pairs |
| `summary.json` | Counts, windows, overlap check |
| `unique_docking_pool.csv` | Unique SMILES across both benchmarks (dock once) |

## Usage notes

- **TrueDecoy** is the harder test (experimental weak + property-matched).
- **RandomDecoy** is the easier / VS-like control; must not overlap TrueDecoy negatives.
- Dock **`unique_docking_pool.csv` once**, then join scores into each benchmark for EF/AUC.
- Do **not** train ML models on these decoy labels.

## Rebuild

```bash
python3 scripts/build_urat1_true_decoy.py \\
  --pactivity-min {pmin} \\
  --inactive-pactivity-max {inactive_max} \\
  --ratio {ratio_target} \\
  --max-sim-to-active {max_sim} \\
  --seed {seed}
```
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build URAT1 TrueDecoy / RandomDecoy benchmarks")
    parser.add_argument("--actives", type=Path, default=DEFAULT_ACTIVES, help="Curated URAT1 CSV")
    parser.add_argument("--pool", type=Path, default=DEFAULT_POOL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--pactivity-min", type=float, default=6.0, help="Active threshold")
    parser.add_argument(
        "--inactive-pactivity-max",
        type=float,
        default=5.0,
        help="Experimental inactives: pActivity < this value (default 5 → ~80 molecules)",
    )
    parser.add_argument("--ratio", type=int, default=10, help="Target decoys per active (TrueDecoy)")
    parser.add_argument("--max-sim-to-active", type=float, default=0.5, help="Max Morgan TC to any active")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-random", action="store_true", help="Skip RandomDecoy companion")
    args = parser.parse_args()

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.actives, low_memory=False)
    smi_col = "canonical_smiles" if "canonical_smiles" in raw.columns else "Smiles"
    raw = raw.copy()
    raw["canonical_smiles"] = raw[smi_col].map(canonicalize)
    raw["pActivity"] = pd.to_numeric(raw["pActivity"], errors="coerce")
    raw = raw.dropna(subset=["canonical_smiles", "pActivity"]).drop_duplicates("canonical_smiles")

    # --- Actives ---
    raw_act = raw[raw["pActivity"] >= args.pactivity_min].copy()
    act_props = annotate_props(raw_act)
    keep_cols = [c for c in ("pActivity", "molecule_chembl_id") if c in raw_act.columns]
    meta = raw_act[["canonical_smiles"] + keep_cols].drop_duplicates("canonical_smiles")
    actives = act_props.merge(meta, on="canonical_smiles", how="left")
    actives["scaffold"] = actives["canonical_smiles"].map(murcko_scaffold)
    actives["label"] = 1
    actives["set_role"] = "active"

    # --- Experimental inactives (TrueDecoy seed negatives) ---
    raw_inact = raw[raw["pActivity"] < args.inactive_pactivity_max].copy()
    # never allow overlap with actives
    raw_inact = raw_inact[~raw_inact["canonical_smiles"].isin(set(actives["canonical_smiles"]))]
    exp_props = annotate_props(raw_inact)
    exp_meta = raw_inact[["canonical_smiles"] + keep_cols].drop_duplicates("canonical_smiles")
    experimental = exp_props.merge(exp_meta, on="canonical_smiles", how="left")
    experimental["scaffold"] = experimental["canonical_smiles"].map(murcko_scaffold)
    experimental["label"] = 0
    experimental["set_role"] = "true_decoy"
    experimental["decoy_class"] = "true"
    experimental["decoy_source"] = "experimental_inactive"

    active_set = set(actives["canonical_smiles"])
    curated_all = set(raw["canonical_smiles"])
    excluded = load_exclude_smiles(EXCLUDE_EXTRA, active_set | curated_all | set(experimental["canonical_smiles"]))

    # --- Pool (unlabeled) ---
    raw_pool = pd.read_csv(args.pool, low_memory=False)
    pool_smi_col = "canonical_smiles" if "canonical_smiles" in raw_pool.columns else "Smiles"
    pool_props = annotate_props(raw_pool.rename(columns={pool_smi_col: "canonical_smiles"}))
    pool_props = pool_props[~pool_props["canonical_smiles"].isin(excluded)].reset_index(drop=True)

    n_target_decoys = int(len(actives) * args.ratio)
    n_exp = len(experimental)
    n_match_want = max(0, n_target_decoys - n_exp)

    print(f"Actives (pActivity>={args.pactivity_min}): {len(actives)}")
    print(f"Experimental inactives (pActivity<{args.inactive_pactivity_max}): {n_exp}")
    print(f"Target TrueDecoy negatives (1:{args.ratio}): {n_target_decoys}")
    print(f"Property-matched still needed: {n_match_want}")
    print(f"Decoy pool after exclude: {len(pool_props)}")

    # Match only as many as needed to fill the 1:ratio quota
    matched, assignments = match_true_decoys(
        actives,
        pool_props,
        ratio=args.ratio,  # per-active cap; we trim below to exact remaining quota
        seed=args.seed,
        max_sim_to_active=args.max_sim_to_active,
    )
    if len(matched) > n_match_want:
        # Keep earliest assignments (round-robin order) up to quota
        keep_smiles = set(assignments.head(n_match_want)["decoy_smiles"])
        matched = matched[matched["canonical_smiles"].isin(keep_smiles)].reset_index(drop=True)
        assignments = assignments.head(n_match_want).reset_index(drop=True)

    matched = matched.copy()
    matched["scaffold"] = matched["canonical_smiles"].map(murcko_scaffold)
    matched["label"] = 0
    matched["set_role"] = "true_decoy"
    matched["decoy_class"] = "true"
    matched["decoy_source"] = "property_matched"
    # matched rows may lack pActivity / chembl id
    if "pActivity" not in matched.columns:
        matched["pActivity"] = np.nan
    if "molecule_chembl_id" not in matched.columns:
        matched["molecule_chembl_id"] = np.nan

    # Align columns for concat
    decoy_cols = [
        "canonical_smiles",
        *PROP_COLS,
        "scaffold",
        "label",
        "set_role",
        "decoy_class",
        "decoy_source",
        "pActivity",
        "molecule_chembl_id",
    ]
    if "max_tc_active" in matched.columns:
        decoy_cols.insert(-2, "max_tc_active")
    for c in decoy_cols:
        if c not in experimental.columns:
            experimental[c] = np.nan if c not in ("label",) else experimental.get(c, 0)
        if c not in matched.columns:
            matched[c] = np.nan

    true_decoys = pd.concat([experimental[decoy_cols], matched[decoy_cols]], ignore_index=True)
    true_decoys = true_decoys.drop_duplicates("canonical_smiles").reset_index(drop=True)

    ratio_achieved = len(true_decoys) / max(len(actives), 1)
    print(f"TrueDecoys total: {len(true_decoys)} (exp={n_exp}, matched={len(matched)}; ratio 1:{ratio_achieved:.2f})")
    print(f"Matched assignments kept: {len(assignments)}")

    true_bench = pd.concat(
        [actives.assign(decoy_class="true", decoy_source="active"), true_decoys],
        ignore_index=True,
    )

    # --- RandomDecoy: Gu-style random from remaining pool ONLY (no overlap) ---
    n_random = 0
    overlap = 0
    random_decoys = pd.DataFrame()
    if not args.skip_random:
        used_true = set(true_decoys["canonical_smiles"])
        unused = pool_props[~pool_props["canonical_smiles"].isin(used_true)].copy()
        n_want = len(true_decoys)  # same count target as TrueDecoy
        if len(unused) < n_want:
            print(
                f"WARNING: only {len(unused)} unused pool molecules for RandomDecoy "
                f"(want {n_want}). Taking all remaining; ratio will be < 1:{args.ratio}."
            )
        n_take = min(n_want, len(unused))
        random_decoys = unused.sample(n=n_take, random_state=args.seed).reset_index(drop=True)
        random_decoys["scaffold"] = random_decoys["canonical_smiles"].map(murcko_scaffold)
        random_decoys["label"] = 0
        random_decoys["set_role"] = "random_decoy"
        random_decoys["decoy_class"] = "random"
        random_decoys["decoy_source"] = "library_random"
        random_decoys["pActivity"] = np.nan
        random_decoys["molecule_chembl_id"] = np.nan
        n_random = len(random_decoys)
        overlap = len(set(random_decoys["canonical_smiles"]) & used_true)
        assert overlap == 0, f"True/Random decoy overlap must be 0, got {overlap}"

        random_bench = pd.concat(
            [actives.assign(decoy_class="random", decoy_source="active"), random_decoys],
            ignore_index=True,
        )
        random_decoys.to_csv(out_dir / "random_decoys.csv", index=False)
        random_bench.to_csv(out_dir / "random_decoy_benchmark.csv", index=False)
        print(f"RandomDecoys selected: {n_random} (overlap with TrueDecoy={overlap})")

    # Unique docking pool (dock once)
    parts = [actives[["canonical_smiles"]].assign(in_true=1, in_random=1, role="active")]
    parts.append(
        true_decoys[["canonical_smiles"]].assign(in_true=1, in_random=0, role="true_decoy")
    )
    if n_random:
        parts.append(
            random_decoys[["canonical_smiles"]].assign(in_true=0, in_random=1, role="random_decoy")
        )
    uniq = pd.concat(parts, ignore_index=True)
    uniq = (
        uniq.groupby("canonical_smiles", as_index=False)
        .agg(in_true=("in_true", "max"), in_random=("in_random", "max"), role=("role", "first"))
        .reset_index(drop=True)
    )
    # fix role if somehow both — prefer active
    # already first within group; recompute for decoys that appear only once
    uniq.to_csv(out_dir / "unique_docking_pool.csv", index=False)

    # Write outputs
    actives.to_csv(out_dir / "actives.csv", index=False)
    experimental.to_csv(out_dir / "experimental_inactives.csv", index=False)
    true_decoys.to_csv(out_dir / "true_decoys.csv", index=False)
    true_bench.to_csv(out_dir / "true_decoy_benchmark.csv", index=False)
    assignments.to_csv(out_dir / "matching_assignments.csv", index=False)

    ratio_random = n_random / max(len(actives), 1)
    summary = {
        "target": "URAT1",
        "framing": "Gu et al. Nat Mach Intell 2025 TrueDecoy vs RandomDecoy (single-target adaptation)",
        "pactivity_min": args.pactivity_min,
        "inactive_pactivity_max": args.inactive_pactivity_max,
        "n_actives": int(len(actives)),
        "n_experimental_inactives": int(n_exp),
        "n_property_matched_decoys": int(len(matched)),
        "n_true_decoys": int(len(true_decoys)),
        "n_random_decoys": int(n_random),
        "n_unique_docking_pool": int(len(uniq)),
        "target_ratio": args.ratio,
        "achieved_true_ratio": float(ratio_achieved),
        "achieved_random_ratio": float(ratio_random),
        "true_random_decoy_smiles_overlap": int(overlap),
        "random_non_overlapping": True,
        "max_sim_to_active": args.max_sim_to_active,
        "seed": args.seed,
        "property_windows": PROP_WINDOWS,
        "relaxed_window_scale": 1.5,
        "matching_algorithm": "experimental_inactives_plus_round_robin_property_match",
        "random_sampling": "gu_style_random_from_remaining_pool_only_no_fallback",
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
        "actives_with_ge1_matched_decoy": int(assignments["active_smiles"].nunique())
        if len(assignments)
        else 0,
        "decoys_per_active": {
            "min": int(assignments.groupby("active_smiles").size().min()) if len(assignments) else 0,
            "median": float(assignments.groupby("active_smiles").size().median())
            if len(assignments)
            else 0,
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
        n_true_exp=n_exp,
        n_true_matched=len(matched),
        n_random=n_random,
        pmin=args.pactivity_min,
        inactive_max=args.inactive_pactivity_max,
        ratio_target=args.ratio,
        ratio_achieved=ratio_achieved,
        ratio_random=ratio_random,
        max_sim=args.max_sim_to_active,
        seed=args.seed,
        n_covered=int(summary["actives_with_ge1_matched_decoy"]),
        decoys_per_active=summary["decoys_per_active"],
        overlap=overlap,
    )

    print(f"Wrote benchmark to {out_dir}")
    print(f"Unique docking pool: {len(uniq)} molecules")


if __name__ == "__main__":
    main()
