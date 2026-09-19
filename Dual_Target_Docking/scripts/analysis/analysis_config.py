"""Single analysis configuration actually imported by compute scripts.

Do not duplicate these constants in publication prose. Changing a value here
must change the numbers that rebuild writes.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CANON = ROOT / "results" / "canonical"
DEFAULT_MASTER_NAME = "current_score_master.csv"

PRIMARY_PAIRS = (
    "EGFR/HER2",
    "JAK1/JAK2",
    "JAK1/TYK2",
    "PIK3CA/mTOR",
    "AChE/BChE",
    "F2/F10",
    "PPARG/PPARA",
    "PPARA/PPARD",
)
RECEPTORS = {
    "EGFR/HER2": ("3POZ", "3RCD"),
    "JAK1/JAK2": ("6N7A", "8BXH"),
    "JAK1/TYK2": ("6N7A", "3LXP"),
    "PIK3CA/mTOR": ("4L23", "4JT6"),
    "AChE/BChE": ("4EY7", "4BDS"),
    "F2/F10": ("4UDW", "2JKH"),
    "PPARG/PPARA": ("9V8H", "6LXA"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}

# Labels
THETA_PRIMARY = 6.0
STRICT_HI = 6.5
STRICT_LO = 5.5

# Scores: Vina/GNINA mode-1 energy is lower-better; stored score_S = -energy is higher-better.
SCORE_DIRECTION = "higher_better"
D_VS_A_SCORE = "score_B"  # dual vs A-only uses pocket B
D_VS_B_SCORE = "score_A"  # dual vs B-only uses pocket A
SUMMARY_MIN_RULE = "weaker_arm"  # min(D-vs-A pocket B, D-vs-B pocket A)
TWO_POCKET_MEAN_USES = ("dual_vs_neither", "ranking")

# Sample gate for primary estimands
REQUIRE_ANALYSIS_SET = "main"
REQUIRE_COMPLETE_CASE = True
REQUIRE_ACTIVITY_ELIGIBLE = True

# Statistics
N_BOOT = 2000
SEED = 20260729
CI_PERCENTILES = (2.5, 97.5)
BOOTSTRAP = "class_stratified_shared_dual"
TOP_FRACTION = 0.10
RANKING_SCORE = "score_mean"

# Detectable-effect simulation (not observed power)
N_MC_DETECTABLE = 1000
DETECTABLE_TRUE_AUCS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75)

# Ligand models
ECFP_RADIUS = 2
ECFP_NBITS = 2048
LOGREG_C = 1.0
LOGREG_MAX_ITER = 4000
GROUPKFOLD_MAX_SPLITS = 5
DESCRIPTOR_NAMES = ("tpsa", "clogp", "heavy", "mw")

# Deposited robustness inputs (zero-dock; scores already on disk)
GNINA_SOURCES = {
    "EGFR/HER2": (
        ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_EGFR_HER2.csv",
        "3POZ",
        "3RCD",
    ),
    "PIK3CA/mTOR": (
        ROOT / "data/jcim_independent_dock_v0/tables/gnina_dock_scores_PIK3CA_mTOR.csv",
        "4L23",
        "4JT6",
    ),
    "JAK1/TYK2": (
        ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_gnina_independent_jak1_tyk2_v1.csv",
        "6N7A",
        "3LXP",
    ),
}
RECEPTOR_SUB_SPECS = (
    ("4JPS", ROOT / "data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt4JPS.csv", "A"),
    ("5DXT", ROOT / "data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt5DXT.csv", "A"),
    ("4JSX", ROOT / "data/jcim_structure_robust_v0/tables/scores_vina_mode1_PM48_alt4JSX.csv", "B"),
)
COGNATE_RMSD_SOURCE = ROOT / "data/jcim_novelty_v0/tables/all14_cognate_rmsd_calcrrms_v1.csv"

EGFR_SCORE_SOURCE = "data/egfr_her2_panel120_v0/tables/ablation_ligand_scores.csv"
EGFR_UNIFORM_VINA_REL = "data/egfr_her2_uniform_rdkit_v1/tables/scores_vina_mode1_fiveseed.csv"
EGFR_UNIFORM_GNINA_REL = "data/egfr_her2_uniform_rdkit_v1/tables/gnina_dock_scores_EGFR_HER2.csv"
EGFR_UNIFORM_VINA_CSV = ROOT / EGFR_UNIFORM_VINA_REL
EGFR_UNIFORM_GNINA_CSV = ROOT / EGFR_UNIFORM_GNINA_REL
EGFR_PRODUCTION_SEED = 20260727
ACHE_SCORE_SOURCE = "data/ache_bche_panel_v0/tables/ablation_ligand_scores.csv"
PIK3CA_SCORE_SOURCE = "data/pik3ca_mtor_panel48_rdkit_v0/tables/ablation_ligand_scores.csv"
TRACK_B_SCORE_SOURCE = "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_vina_mode1_v1.csv"

FIVE_SEEDS = (20260727, 20260811, 20260812, 20260813, 20260814)
FIVE_SEED_LONG_SCORES = ROOT / "data/jcim_multiseed_v0/tables/multiseed_scores_long_v1.csv"
FIVE_SEED_TRACKB_TEMPLATE = (
    ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/multiseed/scores_vina_mode1_seed{seed}.csv"
)
FIVE_SEED_EGFR_CORRECTED_BOX_SCORES = (
    ROOT / "data/jcim_multiseed_v0/tables/scores_vina_mode1_EGFR_corrected_box_fiveseed.csv"
)
PM110_SCORE_SOURCE = ROOT / "data/pik3ca_mtor_panel110_rdkit_v0/tables/ablation_ligand_scores.csv"
E8_SCORE_SOURCE = ROOT / "data/pik3ca_mtor_panel48_rdkit_v0/tables/scores_vina_E8_best.csv"
EXTERNAL_ELIGIBILITY_SOURCE = ROOT / "data/jcim_novelty_v0/tables/external_slice_summary_v1.csv"
UNIVERSE_CENSUS_SOURCE = ROOT / "data/jcim_chembl_universe_v0/tables/universe_census_summary_v1.csv"
PAIR_ELIGIBILITY_SOURCE = ROOT / "data/jcim_chembl_universe_v0/tables/pair_eligibility_audit_s14_v1.csv"
RTM_BEST9_SOURCE = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_rtm_best9_v1.csv"
CNN_BEST9_SOURCE = ROOT / "data/jcim_chembl_universe_v0/local_track_b_v0/tables/scores_gnina_cnn_best9_v1.csv"

# Scientific CSVs written by the full zero-dock freeze chain (filename set after promotion).
CANONICAL_CSV_NAMES = (
    "current_score_master.csv",
    "primary_directional_auroc.csv",
    "primary_summary_min.csv",
    "fixed_score_negative_class_delta.csv",
    "two_pocket_mean_ranking.csv",
    "top10_operating_points.csv",
    "and_filter_operating_points.csv",
    "matched_mismatched_pocket.csv",
    "matched_minus_mismatched.csv",
    "holdout_metrics.csv",
    "class_counts.csv",
    "label_aggregation_sensitivity.csv",
    "max_vs_median_sensitivity.csv",
    "cluster_bootstrap_sensitivity.csv",
    "non_stratified_bootstrap_sensitivity.csv",
    "computational_robustness.csv",
    "receptor_substitution.csv",
    "cognate_rmsd.csv",
    "five_seed_summary_min.csv",
    "five_seed_fixed_membership_sensitivity.csv",
    "protocol_sensitivity.csv",
    "external_eligibility.csv",
    "descriptor_baselines.csv",
    "descriptor_nested_scaffold_cv.csv",
    "ecfp4_incremental_information.csv",
    "ecfp4_oof_predictions.csv",
    "ecfp4_scaler_sensitivity.csv",
    "model_fold_assignments.csv",
    "class_chemistry_summary.csv",
    "leave_one_document_delta.csv",
    "detectable_effect_simulation.csv",
)


def parse_finite(raw) -> float | None:
    """Return a finite float, else None. Adopted from PR #35 membership hardening."""
    if raw is None:
        return None
    text = str(raw).strip()
    if text == "":
        return None
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def is_true(value) -> bool:
    return value in (1, "1", True, "True")


def is_primary_row(row: dict, analysis_set: str = REQUIRE_ANALYSIS_SET) -> bool:
    if row.get("pair") not in PRIMARY_PAIRS:
        return False
    if row.get("analysis_set") != analysis_set:
        return False
    if REQUIRE_COMPLETE_CASE and not is_true(row.get("complete_case")):
        return False
    if REQUIRE_ACTIVITY_ELIGIBLE and not is_true(row.get("activity_eligible", "1")):
        return False
    return True


def io_paths(outdir=None, master=None) -> tuple[Path, Path]:
    out = Path(outdir) if outdir else DEFAULT_CANON
    if not out.is_absolute():
        out = (ROOT / out).resolve()
    else:
        out = out.resolve()
    if master:
        m = Path(master)
        m = m if m.is_absolute() else (ROOT / m).resolve()
    else:
        m = out / DEFAULT_MASTER_NAME
    return out, m


def egfr_uniform_ready() -> bool:
    """True when the CASE 2 uniform Vina table exists. That table is then current truth."""
    return EGFR_UNIFORM_VINA_CSV.is_file() and EGFR_UNIFORM_VINA_CSV.stat().st_size > 0


def current_egfr_score_source() -> str:
    if egfr_uniform_ready():
        return EGFR_UNIFORM_VINA_REL
    return EGFR_SCORE_SOURCE


def current_egfr_five_seed_path() -> Path:
    if egfr_uniform_ready():
        return EGFR_UNIFORM_VINA_CSV
    return FIVE_SEED_EGFR_CORRECTED_BOX_SCORES


def current_gnina_sources() -> dict:
    src = dict(GNINA_SOURCES)
    if EGFR_UNIFORM_GNINA_CSV.is_file() and EGFR_UNIFORM_GNINA_CSV.stat().st_size > 0:
        src["EGFR/HER2"] = (EGFR_UNIFORM_GNINA_CSV, "3POZ", "3RCD")
    return src


def add_io_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "--outdir",
        default=None,
        help="CSV output directory (default: results/canonical). Freeze rebuild uses a new directory.",
    )
    parser.add_argument(
        "--master",
        default=None,
        help="Score-master CSV (default: <outdir>/current_score_master.csv).",
    )
    return parser
