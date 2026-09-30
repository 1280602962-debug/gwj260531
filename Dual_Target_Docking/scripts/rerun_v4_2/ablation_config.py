"""Frozen identities for SCORING_AND_DOCKING_METHOD_ABLATION_V2."""
from __future__ import annotations

from pathlib import Path

ROOT = Path("/tmp/pr39_fiveseed/Dual_Target_Docking")
RUN_ID = "UNIFORM_RERUN_V4_2_20260921"
RUN = ROOT / "reruns" / RUN_ID

PAIRS = (
    "EGFR/HER2", "JAK1/JAK2", "JAK1/TYK2", "PIK3CA/mTOR",
    "AChE/BChE", "F2/F10", "PPARG/PPARA", "PPARA/PPARD",
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
UNIQUE_14 = (
    "3POZ", "3RCD", "6N7A", "8BXH", "3LXP", "4L23", "4JT6",
    "4EY7", "4BDS", "4UDW", "2JKH", "9V8H", "6LXA", "5U3Q",
)
ALT_FORBIDDEN = {"4L2Y", "4JT5", "4EY6", "1P0M", "3SHC", "2Y5F", "6KAX", "5U46"}

GNINA_ROOT = Path("/mnt/d/CADD paper exercise/gnina")
GNINA_BIN = GNINA_ROOT / "bin" / "gnina"
GNINA_LIB = GNINA_ROOT / "conda_env" / "lib"
GNINA_ACTIVATE = GNINA_ROOT / "activate.sh"
OBABEL = GNINA_ROOT / "conda_env" / "bin" / "obabel"

RTM_PY = Path("/home/gwj/miniconda3/envs/rtmscore/bin/python")
RTM_ROOT = Path("/home/gwj/software/RTMScore")
RTM_SCRIPT = RTM_ROOT / "example" / "rtmscore.py"
RTM_CKPT = RTM_ROOT / "trained_models" / "rtmscore_model1.pth"

REC = RUN / "01_receptors"
BOX = RUN / "03_boxes"
LIG_PDBQT = RUN / "02_ligands" / "pdbqt"
LIG_SDF = RUN / "02_ligands" / "sdf"
VINA_JOBS = RUN / "06_vina_fiveseed" / "jobs"
MAP = RUN / "00_protocol" / "pair_ligand_mapping.csv"
PHASE7 = RUN / "13_qa" / "phase7_fiveseed_master.csv"
QA = RUN / "13_qa"
SEC = RUN / "09_secondary_scoring"
GIND = RUN / "10_gnina_independent_docking"
# GNINA 1.3.2 cannot parse Meeko 0.7.1 extra types CG0/G0. Copies only; frozen ligands untouched.
LIG_GNINA_PDBQT = GIND / "ligands_gnina_ad4_types"
ANALYSIS = RUN / "08_analysis"
PROTO = RUN / "00_protocol"

BOOT_B = 10000
BOOT_SEED = 271828
TOP_FRACTION = 0.10
THETA = 6.0
POCKET_CUTOFF = 10.0
RESCORING_WORKERS = 2  # keep alt Vina at 6
