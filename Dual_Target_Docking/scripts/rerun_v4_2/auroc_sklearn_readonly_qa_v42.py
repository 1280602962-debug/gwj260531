#!/usr/bin/env python3
"""Read-only sklearn cross-check of published PRIMARY directional rows.

80 directional AUROC rows: 8 pairs x 5 methods x 2 directions.
40 summary_min rows: min(AUC_B, AUC_A) only. Never AUROC(min(score_A, score_B)).

Does not write PRIMARY tables. Does not run bootstrap.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from sklearn.metrics import roc_auc_score

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
import formal_metrics_lib as F  # noqa: E402


def load(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    F.add_project_root_arg(p)
    args = p.parse_args(argv)
    root = F.resolve_project_root(args.project_root)
    run = root / "reruns" / "UNIFORM_RERUN_V4_2_20260921"
    qa = run / "13_qa"
    published = load(root / "results" / "formal_metrics" / "PRIMARY_DIRECTIONAL_METRICS.csv")
    pop = load(qa / "PRIMARY_DIRECTIONAL_ANALYSIS_POPULATION.csv")
    by_pair: dict[str, list[dict]] = defaultdict(list)
    for r in pop:
        by_pair[r["pair"]].append(r)
    pub = {(r["pair"], r["method"], r["contrast"]): r for r in published}

    n_dir = n_dir_mis = n_min = n_min_mis = 0
    issues = []
    for pair in F.RECEPTORS:
        sub = by_pair[pair]
        for meth in F.FORMAL_METHODS:
            u = F.shared_dual_universe(sub, meth)
            dual, ao, bo = u["dual_shared"], u["A_only_for_AUC_B"], u["B_only_for_AUC_A"]
            yb = [1] * len(dual) + [0] * len(ao)
            xb = [float(r[f"{meth}_score_B"]) for r in dual] + [float(r[f"{meth}_score_B"]) for r in ao]
            ya = [1] * len(dual) + [0] * len(bo)
            xa = [float(r[f"{meth}_score_A"]) for r in dual] + [float(r[f"{meth}_score_A"]) for r in bo]
            auc_b = float(roc_auc_score(yb, xb))
            auc_a = float(roc_auc_score(ya, xa))
            sm = min(auc_b, auc_a)
            for contrast, val in (
                ("dual_vs_A_only", auc_b),
                ("dual_vs_B_only", auc_a),
            ):
                n_dir += 1
                got = float(pub[(pair, meth, contrast)]["auroc"])
                if abs(got - val) > 5e-10:
                    n_dir_mis += 1
                    issues.append(f"{pair} {meth} {contrast} sklearn={val} published={got}")
            n_min += 1
            got_sm = float(pub[(pair, meth, "summary_min")]["auroc"])
            if abs(got_sm - sm) > 5e-10:
                n_min_mis += 1
                issues.append(f"{pair} {meth} summary_min derived={sm} published={got_sm}")

    report = {
        "n_directional_auroc_checked": n_dir,
        "n_directional_mismatch": n_dir_mis,
        "n_summary_min_checked": n_min,
        "n_summary_min_mismatch": n_min_mis,
        "n_checked_rows": n_dir + n_min,
        "forbidden_auroc_of_min_score": "not_computed",
        "issues": issues[:30],
        "PASS": n_dir == 80 and n_min == 40 and n_dir_mis == 0 and n_min_mis == 0,
    }
    out = qa / "AUROC_SKLEARN_READONLY_QA.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["PASS"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
