#!/usr/bin/env python3
"""Temporary-directory checks for replay_track_b_vina_mode1_v1."""
from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("replay_mod", HERE / "replay_track_b_vina_mode1_v1.py")
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def pose_text(energy: float) -> str:
    return f"MODEL 1\nREMARK VINA RESULT: {energy:.3f}  0  0\nENDMDL\n"


def run_case(tmp: Path, jobs, scores, poses: dict[tuple, str | None], extra_argv=None) -> int:
    MOD.SCORE = tmp / "scores.csv"
    MOD.STATUS = tmp / "jobs.csv"
    MOD.DEFAULT_OUT = tmp / "committed_replay.csv"
    write_csv(MOD.SCORE, scores)
    write_csv(MOD.STATUS, jobs)
    pose_root = tmp / "poses"
    for (target, ligand), body in poses.items():
        d = pose_root / target / ligand
        d.mkdir(parents=True, exist_ok=True)
        if body is not None:
            (d / "mode_01.pdbqt").write_text(body, encoding="utf-8")
    out = tmp / "out.csv"
    argv = ["--pose-root", str(pose_root), "--output", str(out)]
    if extra_argv:
        argv.extend(extra_argv)
    return MOD.main(argv), out


def test_match_and_skip():
    tmp = Path(tempfile.mkdtemp())
    jobs = [
        {"pair": "P", "target": "T", "ligand": "L1", "status": "success"},
        {"pair": "P", "target": "T", "ligand": "L2", "status": "timeout"},
    ]
    scores = [
        {"pair": "P", "target": "T", "ligand": "L1", "mode1_energy": "-9.5"},
        {"pair": "P", "target": "T", "ligand": "L2", "mode1_energy": ""},
    ]
    code, out = run_case(
        tmp,
        jobs,
        scores,
        {("T", "L1"): pose_text(-9.5)},
    )
    rows = list(csv.DictReader(out.open()))
    assert code == 0
    assert rows[0]["replay"] == "match"
    assert rows[0]["pose_location"] == "local_disk_gitignored"
    assert rows[1]["replay"] == "not_applicable_skipped"
    print("PASS match_and_skip")


def test_missing_is_nonzero():
    tmp = Path(tempfile.mkdtemp())
    jobs = [{"pair": "P", "target": "T", "ligand": "L1", "status": "success"}]
    scores = [{"pair": "P", "target": "T", "ligand": "L1", "mode1_energy": "-9.5"}]
    code, out = run_case(tmp, jobs, scores, {})
    rows = list(csv.DictReader(out.open()))
    assert code == 2
    assert rows[0]["replay"] == "pose_missing"
    print("PASS missing_is_nonzero")


def test_parse_fail_and_mismatch():
    tmp = Path(tempfile.mkdtemp())
    jobs = [
        {"pair": "P", "target": "T", "ligand": "L1", "status": "success"},
        {"pair": "P", "target": "T", "ligand": "L2", "status": "success"},
    ]
    scores = [
        {"pair": "P", "target": "T", "ligand": "L1", "mode1_energy": "-9.5"},
        {"pair": "P", "target": "T", "ligand": "L2", "mode1_energy": "-8.0"},
    ]
    code, out = run_case(
        tmp,
        jobs,
        scores,
        {("T", "L1"): "MODEL 1\nENDMDL\n", ("T", "L2"): pose_text(-1.0)},
    )
    rows = {r["ligand"]: r for r in csv.DictReader(out.open())}
    assert code == 2
    assert rows["L1"]["replay"] == "parse_fail"
    assert rows["L2"]["replay"] == "mismatch"
    print("PASS parse_fail_and_mismatch")


def test_refuse_committed_overwrite():
    tmp = Path(tempfile.mkdtemp())
    MOD.SCORE = tmp / "scores.csv"
    MOD.STATUS = tmp / "jobs.csv"
    committed = tmp / "committed_replay.csv"
    MOD.DEFAULT_OUT = committed
    write_csv(MOD.SCORE, [{"pair": "P", "target": "T", "ligand": "L1", "mode1_energy": "-1"}])
    write_csv(MOD.STATUS, [{"pair": "P", "target": "T", "ligand": "L1", "status": "success"}])
    committed.write_text("keep\n", encoding="utf-8")
    code = MOD.main(["--pose-root", str(tmp / "poses")])
    assert code == 2
    assert committed.read_text(encoding="utf-8") == "keep\n"
    print("PASS refuse_committed_overwrite")


def test_duplicate_score_key():
    tmp = Path(tempfile.mkdtemp())
    MOD.SCORE = tmp / "scores.csv"
    MOD.STATUS = tmp / "jobs.csv"
    MOD.DEFAULT_OUT = tmp / "committed_replay.csv"
    write_csv(
        MOD.SCORE,
        [
            {"pair": "P", "target": "T", "ligand": "L1", "mode1_energy": "-1"},
            {"pair": "P", "target": "T", "ligand": "L1", "mode1_energy": "-2"},
        ],
    )
    write_csv(MOD.STATUS, [{"pair": "P", "target": "T", "ligand": "L1", "status": "success"}])
    try:
        MOD.main(["--pose-root", str(tmp / "poses"), "--output", str(tmp / "out.csv")])
    except SystemExit as exc:
        assert "duplicate score" in str(exc)
        print("PASS duplicate_score_key")
        return
    raise AssertionError("duplicate score keys should abort")


if __name__ == "__main__":
    test_match_and_skip()
    test_missing_is_nonzero()
    test_parse_fail_and_mismatch()
    test_refuse_committed_overwrite()
    test_duplicate_score_key()
    print("all replay tests passed")
