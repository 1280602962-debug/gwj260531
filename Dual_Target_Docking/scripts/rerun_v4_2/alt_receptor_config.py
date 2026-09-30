"""Frozen alternative-receptor identities from ALTERNATIVE_RECEPTOR_EXPERIMENT_FREEZE_FINAL.yaml."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "UNIFORM_RERUN_V4_2_20260921"
RUN = ROOT / "reruns" / RUN_ID
CIF_ALT = RUN / "13_qa" / "alt_cif"
CIF_CACHE = Path("/tmp/input_identity_rcsb")

FORMAL_PAIRS = ("PIK3CA/mTOR", "AChE/BChE", "F2/F10", "PPARA/PPARD")
PRIMARY = {
    "PIK3CA/mTOR": ("4L23", "4JT6"),
    "AChE/BChE": ("4EY7", "4BDS"),
    "F2/F10": ("4UDW", "2JKH"),
    "PPARA/PPARD": ("6LXA", "5U3Q"),
}

# Eight one-at-a-time substitutions. Do not add PDBs.
ALTS = {
    "4L2Y": {
        "id": "ALT_4L2Y", "pair": "PIK3CA/mTOR", "side": "A", "design": "A2/B1",
        "kept_primary": "4JT6",
        "protein_auth": ("A",),
        "box": {"ccd": "XXK", "auth": "A", "label": "C", "res": 1101, "instance": "A:1101"},
        "remove_auth": ("B",),  # p85
        "retain_auth": (),
    },
    "4JT5": {
        "id": "ALT_4JT5", "pair": "PIK3CA/mTOR", "side": "B", "design": "A1/B2",
        "kept_primary": "4L23",
        "protein_auth": ("A",),  # label C
        "box": {"ccd": "P2X", "auth": "A", "label": "F", "res": 2601, "instance": "A:2601"},
        "remove_auth": ("B", "C", "D"),  # second mTOR + both mLST8
        "retain_auth": ("A",),
        "equiv_nonprimary": "B:2601",
    },
    "4EY6": {
        "id": "ALT_4EY6", "pair": "AChE/BChE", "side": "A", "design": "A2/B1",
        "kept_primary": "4BDS",
        "protein_auth": ("A", "B"),
        "box": {"ccd": "GNT", "auth": "A", "label": "E", "res": 604, "instance": "A:604"},
        "remove_auth": (),
        "retain_auth": ("A", "B"),
        "equiv_nonprimary": "B:605",
    },
    "1P0M": {
        "id": "ALT_1P0M", "pair": "AChE/BChE", "side": "B", "design": "A1/B2",
        "kept_primary": "4EY7",
        "protein_auth": ("A",),
        "box": {"ccd": "CHT", "auth": "A", "label": "K", "res": 607, "instance": "A:607"},
        "remove_auth": (),
        "retain_auth": (),
    },
    "3SHC": {
        "id": "ALT_3SHC", "pair": "F2/F10", "side": "A", "design": "A2/B1",
        "kept_primary": "2JKH",
        "protein_auth": ("H", "I", "L"),
        "box": {"ccd": "B01", "auth": "H", "label": "G", "res": 3, "instance": "H:3"},
        "remove_auth": (),
        "retain_auth": ("I",),  # hirudin
    },
    "2Y5F": {
        "id": "ALT_2Y5F", "pair": "F2/F10", "side": "B", "design": "A1/B2",
        "kept_primary": "4UDW",
        "protein_auth": ("A", "L"),
        "box": {"ccd": "XWG", "auth": "A", "label": "C", "res": 1244, "instance": "A:1244"},
        "remove_auth": (),
        "retain_auth": ("L",),
    },
    "6KAX": {
        "id": "ALT_6KAX", "pair": "PPARA/PPARD", "side": "A", "design": "A2/B1",
        "kept_primary": "5U3Q",
        "protein_auth": ("A",),
        "box": {"ccd": "PLM", "auth": "A", "label": "C", "res": 502, "instance": "A:502"},
        "remove_auth": (),
        "retain_auth": (),
    },
    "5U46": {
        "id": "ALT_5U46", "pair": "PPARA/PPARD", "side": "B", "design": "A1/B2",
        "kept_primary": "6LXA",
        "protein_auth": ("A",),
        "box": {"ccd": "7T1", "auth": "A", "label": "G", "res": 505, "instance": "A:505"},
        "remove_auth": (),
        "retain_auth": (),
        "equiv_nonprimary": "B:502",
        "note": "7T1_not_B7G",
    },
}

BOOT_B = 10000
BOOT_SEED = 271828
THETA = 6.0
TOP_FRACTION = 0.10
JOB_TIMEOUT_S = 1586
WORKERS = 6
VINA = Path("/home/gwj/miniconda3/bin/vina")
SEEDS_REDOCK = (17, 29, 42, 71, 101)
PRODUCTION_SEED = 42
