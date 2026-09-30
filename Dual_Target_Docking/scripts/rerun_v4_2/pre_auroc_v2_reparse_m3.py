#!/usr/bin/env python3
"""Reparse M3 CNNscore from first MODEL of each existing SUCCESS out.pdbqt. No redock. No AUROC."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from rerun_v4_2.ablation_config import GIND, QA, UNIQUE_14  # noqa: E402
from rerun_v4_2.pre_auroc_v2_lib import finite_float, first_model_text, parse_model1_remarks  # noqa: E402

OLD = GIND / "GNINA_SEED42_PRODUCTION_MASTER.csv"
OUT = GIND / "GNINA_SEED42_PRODUCTION_MASTER_REPARSED.csv"
JOBS = GIND / "jobs"


def main() -> int:
    rows = list(csv.DictReader(OLD.open(encoding="utf-8-sig", newline="")))
    out_rows = []
    n_ok = n_fail_parse = n_timeout = n_other = 0
    spot = {}
    for r in rows:
        rec = dict(r)
        rec["CNNscore_source"] = ""
        rec["reparse_status"] = ""
        rec["original_master_CNNscore"] = r.get("CNNscore", "")
        jobdir = JOBS / r["job_id"]
        outp = jobdir / "out.pdbqt"
        status = r.get("status", "")
        if status == "TIMEOUT":
            rec["reparse_status"] = "KEEP_TIMEOUT"
            rec["CNNscore_source"] = "none_timeout"
            n_timeout += 1
            out_rows.append(rec)
            continue
        if status != "SUCCESS":
            rec["reparse_status"] = f"KEEP_{status or 'OTHER'}"
            rec["CNNscore_source"] = "none_non_success"
            n_other += 1
            out_rows.append(rec)
            continue
        if not outp.is_file():
            rec["reparse_status"] = "PARSE_FAIL_MISSING_OUT"
            rec["CNNscore"] = ""
            rec["CNNscore_source"] = "missing_out_pdbqt"
            n_fail_parse += 1
            out_rows.append(rec)
            continue
        text = outp.read_text(errors="replace")
        scores = parse_model1_remarks(text)
        rec["CNNscore"] = scores["CNNscore"]
        rec["CNNaffinity"] = scores["CNNaffinity"] or rec.get("CNNaffinity", "")
        rec["CNN_VS"] = scores["CNN_VS"]
        if scores["minimizedAffinity"]:
            rec["empirical_affinity"] = scores["minimizedAffinity"]
            try:
                rec["M3_empirical_higher_better"] = f"{-float(scores['minimizedAffinity']):.6f}"
            except ValueError:
                pass
        rec["CNNscore_source"] = "out.pdbqt_MODEL1_REMARK"
        if finite_float(rec["CNNscore"]):
            rec["reparse_status"] = "SUCCESS_FINITE_CNNSCORE"
            n_ok += 1
        else:
            rec["reparse_status"] = "PARSE_FAIL_NO_FINITE_CNNSCORE"
            n_fail_parse += 1
        pdb = r["pdb_id"]
        if pdb in UNIQUE_14 and pdb not in spot and rec["reparse_status"] == "SUCCESS_FINITE_CNNSCORE":
            model1 = first_model_text(text)
            remarks = [ln for ln in model1.splitlines() if ln.startswith("REMARK") and ("CNN" in ln or "minimizedAffinity" in ln)]
            spot[pdb] = {
                "job_id": r["job_id"],
                "reparsed": {k: rec[k] for k in ("CNNscore", "CNNaffinity", "CNN_VS", "empirical_affinity")},
                "model1_remarks": remarks,
            }
        out_rows.append(rec)
    with OUT.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(out_rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    report = {
        "n_input": len(rows),
        "SUCCESS_FINITE_CNNSCORE": n_ok,
        "PARSE_FAIL": n_fail_parse,
        "KEEP_TIMEOUT": n_timeout,
        "KEEP_OTHER": n_other,
        "spot_check_first_success_per_receptor": {k: spot[k] for k in UNIQUE_14 if k in spot},
        "receptors_missing_success_spot": [k for k in UNIQUE_14 if k not in spot],
        "AUROC": "NOT_COMPUTED",
        "original_master_preserved": str(OLD),
        "reparsed_path": str(OUT),
    }
    (QA / "M3_MODEL1_REPARSE_SPOTCHECK.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in report if k != "spot_check_first_success_per_receptor"}, indent=2), flush=True)
    print("spot_receptors", len(spot), "/", 14, flush=True)
    return 0 if n_fail_parse == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
