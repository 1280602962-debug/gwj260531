#!/usr/bin/env python3
"""Offline checks for the D1 capability inventory and search/selection labels."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from urat1_capability_lib import (  # noqa: E402
    REUSE_SOURCES,
    classify_search_selection,
    d1_jobs,
    d1_new_jobs,
    evaluate_sdf,
    inventory_rows,
)


class CapabilityDockTests(unittest.TestCase):
    def test_inventory_counts(self) -> None:
        jobs = d1_jobs()
        self.assertEqual(len(jobs), 27)
        self.assertEqual(len(REUSE_SOURCES), 7)
        self.assertEqual(len(d1_new_jobs()), 20)
        self.assertTrue(all(j not in REUSE_SOURCES for j in d1_new_jobs()))
        self.assertNotIn("urat1_9dk9", {tgt for _lig, tgt, _seed in jobs})

    def test_inventory_table(self) -> None:
        rows = inventory_rows()
        self.assertEqual(len(rows), 27)
        self.assertEqual(sum(r["action"] == "reuse" for r in rows), 7)
        self.assertEqual(sum(r["action"] == "new_gnina" for r in rows), 20)
        self.assertTrue(all(r["pdb"] in {"9DKA", "9DKB", "9DKC"} for r in rows))

    def test_selection_failure_requires_native_like(self) -> None:
        self.assertEqual(
            classify_search_selection(3.59, 1.11), "search_ok_selection_fail"
        )
        self.assertEqual(classify_search_selection(0.8, 0.7), "selection_ok")
        self.assertEqual(classify_search_selection(4.31, 4.31), "search_incomplete")
        self.assertEqual(classify_search_selection(4.88, 3.87), "search_incomplete")
        self.assertEqual(classify_search_selection(None, 1.1), "unevaluable")

    def test_reuse_sources_exist(self) -> None:
        missing = [str(p) for p in REUSE_SOURCES.values() if not p.exists()]
        self.assertFalse(missing, missing)

    def test_score_existing_selfdocks(self) -> None:
        bbr = evaluate_sdf(
            "benzbromarone",
            "urat1_9dka",
            42,
            REUSE_SOURCES[("benzbromarone", "urat1_9dka", 42)],
        )
        self.assertEqual(bbr["capability_class"], "search_ok_selection_fail")
        self.assertGreater(bbr["top1_pose_rmsd_A"], 2.0)
        self.assertLess(bbr["best_of_9_pose_rmsd_A"], 2.0)
        self.assertGreaterEqual(bbr["best_rmsd_cnnscore_rank"], 2)
        self.assertFalse(bbr["success_at_1"])
        self.assertTrue(bbr["success_at_3"])

        les = evaluate_sdf(
            "lesinurad",
            "urat1_9dkb",
            43,
            REUSE_SOURCES[("lesinurad", "urat1_9dkb", 43)],
        )
        self.assertEqual(les["capability_class"], "search_incomplete")
        self.assertGreater(les["best_of_9_pose_rmsd_A"], 2.0)


if __name__ == "__main__":
    unittest.main()
