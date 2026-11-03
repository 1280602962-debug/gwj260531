#!/usr/bin/env python3
"""Acceptance tests for the post-hoc diagnostics."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analysis_a import build_analysis_a  # noqa: E402
from analysis_b import build_analysis_b  # noqa: E402
from analysis_c import build_analysis_c  # noqa: E402
from posthoc_common import (  # noqa: E402
    INPUTS,
    PROJECT,
    frame_bytes,
    input_manifest,
    load_config,
)


@pytest.fixture(scope="module")
def frames():
    before = input_manifest()
    out = {}
    out.update(build_analysis_a())
    out.update(build_analysis_b())
    out.update(build_analysis_c())
    after = input_manifest()
    assert before["sha256"].tolist() == after["sha256"].tolist()
    return out


def test_formal_population_counts():
    cfg = load_config()["expected_counts"]
    pop = pd.read_csv(PROJECT / INPUTS["population"])
    eligible = pop[pop["activity_eligible"] == 1]
    assert len(eligible) == cfg["eligible_observations"]
    assert eligible["global_ligand_entity_id"].nunique() == cfg["global_entities"]
    counts = eligible["class"].value_counts().to_dict()
    assert counts == {"dual": 248, "A_only": 256, "B_only": 252, "neither": 111}


def test_chemistry_oof_and_folds(frames):
    nearest = frames["A_TEST_NEAREST_TRAIN"]
    assert len(nearest) == 994
    assert set(nearest["fold_id"]).issubset({0, 1, 2, 3, 4})
    units = pd.read_csv(PROJECT / INPUTS["chemistry_units"])
    assert int(units["n_final_chemistry_oof"].sum()) == 994
    assert int(units["n_labeled_population"].sum()) == 1004
    assert int(units["n_m0_directional_population"].sum()) == 1000


def test_scaffold_audit_and_separation(frames):
    audit = frames["A_SCAFFOLD_FIELD_AUDIT"]
    assert set(audit["scaffold_relation"]).issubset({
        "STRING_MATCH", "CANONICAL_EQUIVALENT", "EMPTY_LABEL_MATCH"
    })
    summary = frames["A_FOLD_OVERLAP_SUMMARY"]
    assert int(summary["n_nonempty_scaffolds_in_multiple_folds"].sum()) == 0
    assert int(summary["n_global_entities_in_multiple_folds"].sum()) == 0
    assert int(summary["n_canonical_smiles_in_multiple_folds"].sum()) == 0
    per_fold = frames["A_FOLD_OVERLAP_PER_FOLD"]
    assert int(per_fold["n_shared_nonempty_scaffold"].sum()) == 0
    assert int(per_fold["n_test_global_entity_also_in_train"].sum()) == 0
    assert int(per_fold["n_test_canonical_smiles_also_in_train"].sum()) == 0
    nearest = frames["A_TEST_NEAREST_TRAIN"]
    assert not (nearest["canonical_ligand_id"] == nearest["nearest_train_ligand_id"]).any()
    assert frames["A_TRIGGER_DECISION"].iloc[0]["conditional_scaffold_group_cv"] == "NOT_RUN"
    assert not (ROOT / "outputs" / "A_SCAFFOLD_GROUP_CV_SUMMARY.csv").exists()


def test_scaler_uses_training_fold_only():
    train = np.array([[1.0], [2.0], [3.0]])
    test = np.array([[100.0]])
    scaler = StandardScaler().fit(train)
    assert scaler.mean_[0] == pytest.approx(2.0)
    scaler.transform(test)
    assert scaler.mean_[0] == pytest.approx(2.0)
    leaky = StandardScaler().fit(np.vstack([train, test]))
    assert leaky.mean_[0] != pytest.approx(scaler.mean_[0])


def test_egfr_auroc_reproduction(frames):
    row = frames["B_EGFR_PAIRWISE_RANKING"].iloc[0]
    assert int(row["n_dual"]) == 28
    assert int(row["n_B_only"]) == 32
    assert int(row["n_tie"]) == 0
    exact = 221 / 896
    assert float(row["auroc_ties_0_5"]) == pytest.approx(exact, abs=1e-12)
    assert round(float(row["auroc_ties_0_5"]) + 1e-15, 10) == 0.2466517857
    assert float(row["recomputed_auc_neither"]) == pytest.approx(0.7470238095238095, abs=1e-12)
    assert row["paired_delta_recomputed_here"] == "NO"


def test_topk_baseline_and_sensitivity(frames):
    repro = frames["C_TOPK_BASELINE_REPRODUCTION"]
    assert (repro["status"] == "PASS").all()
    counts = frames["C_TOPK_CONCORDANCE_COUNTS"]
    official = counts[(counts["cutoff"] == "top_10pct") & (counts["scope"] == "all_32")]
    got = dict(zip(official["concordance"], official["n"]))
    assert got == {
        "CONCORDANT_IMPROVE": 14,
        "CONCORDANT_WORSEN": 10,
        "DISCORDANT": 6,
        "NO_CHANGE": 2,
    }
    cmp = frames["C_TOPK_METHOD_COMPARISONS"]
    assert (cmp["member_ids_identical"] == "YES").all()
    assert len(cmp) == 96
    for _, row in cmp.iterrows():
        assert int(row["k"]) == math.ceil(float(row["fraction"]) * int(row["N"]))
        assert int(row["S1_n_dual_top"]) + int(row["S1_n_A_only_top"]) + int(row["S1_n_B_only_top"]) == int(row["k"])
        assert int(row["new_n_dual_top"]) + int(row["new_n_A_only_top"]) + int(row["new_n_B_only_top"]) == int(row["k"])
    f2 = cmp[(cmp["pair"] == "F2/F10") & (cmp["comparison"] == "S3_vs_S1")]
    assert set(f2["concordance"]) == {"CONCORDANT_IMPROVE"}
    pik = cmp[(cmp["pair"] == "PIK3CA/mTOR") & (cmp["comparison"] == "S2_vs_S1")].set_index("cutoff")
    assert pik.loc["top_10pct", "concordance"] == "DISCORDANT"
    assert pik.loc["top_20pct", "concordance"] == "CONCORDANT_WORSEN"
    ppar = cmp[(cmp["pair"] == "PPARG/PPARA") & (cmp["comparison"] == "S3_vs_S1")].set_index("cutoff")
    assert ppar.loc["top_5pct", "concordance"] == "NO_CHANGE"
    assert ppar.loc["top_10pct", "concordance"] == "CONCORDANT_WORSEN"


def test_outputs_have_no_unexplained_nonfinite_tokens(frames):
    for name, frame in frames.items():
        if name.endswith("JOIN_AUDIT"):
            continue
        text = frame_bytes(frame).decode()
        assert "NaN" not in text
        assert ",nan," not in text.lower()
        assert ",inf," not in text.lower()
        assert "Infinity" not in text


def test_numeric_frames_are_repeatable(frames):
    again = {}
    again.update(build_analysis_a())
    again.update(build_analysis_b())
    again.update(build_analysis_c())
    for name, frame in frames.items():
        if name.endswith("JOIN_AUDIT"):
            continue
        assert frame_bytes(frame) == frame_bytes(again[name])
